!-------------------------------------- LICENCE BEGIN ------------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------

program midas_obsSelection
  !
  ! :Purpose: Main program for O-F computations, background check, and thinning
  !
  !           (O-F => Observation minus Forecast, i.e. y-H(x))
  !
  use version_mod
  use codePrecision_mod
  use ramDisk_mod
  use mpi_mod
  use utilities_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use gridStateVector_mod
  use backgroundCheck_mod
  use multi_ir_bgck_mod
  use innovation_mod
  use obsSpaceData_mod
  use columnData_mod
  use obsFiles_mod
  use obsErrors_mod
  use biasCorrectionSat_mod
  use biasCorrectionConv_mod
  use thinning_mod
  use bgckmicrowave_mod
  use bgckssmis_mod
  use bgckcsr_mod
  use bgckOcean_mod 
   
  implicit none

  integer :: datestamp, headerIndex, ierr, nulnam
  type(struct_columnData),target :: columnTrlOnAnlIncLev
  type(struct_columnData),target :: columnTrlOnTrlLev
  type(struct_obs),       target :: obsSpaceData
  type(struct_gsv)               :: stateVectorTrialHighRes
  type(struct_hco),      pointer :: hco_anl => null()
  type(struct_vco),      pointer :: vco_anl => null()
  type(struct_hco),      pointer :: hco_trl => null()
  type(struct_vco),      pointer :: vco_trl => null()
  type(struct_hco),      pointer :: hco_core => null()

  logical :: allocHeightSfc
  integer :: get_max_rss, fnom, fclos

  ! Namelist variables
  logical :: doThinning ! Control whether or not thinning is done

  namelist /namObsSelection/ doThinning

  call ver_printNameAndVersion('obsSelection','Obs Quality Control and Thinning')

  !- 1.0 mpi
  call mpi_initialize

  !- 1.1 timings
  call tmg_init(mpi_myid, 'TMG_INFO')
  call utl_tmg_start(0,'Main')

  !- 1.2 Read the namelist for obsSelection program (if it exists)
  doThinning = .false.
  if (utl_isNamelistPresent('namObsSelection','./flnml')) then
    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    if (ierr /= 0) call utl_abort('midas-obsSelection: Error opening file flnml')
    read(nulnam,nml=namObsSelection,iostat=ierr)
    if (ierr /= 0) call utl_abort('midas-obsSelection: Error reading namelist')
    if (mpi_myid == 0) write(*,nml=namObsSelection)
    ierr = fclos(nulnam)
  else
    write(*,*)
    write(*,*) 'midas-obsSelection: Namelist block namObsSelection is missing in the namelist.'
    write(*,*) '                    The default value will be taken.'
    if (mpi_myid == 0) write(*,nml=namObsSelection)
  end if

  !
  !- Initialize the ram disk
  !
  call ram_setup

  !     
  !- Initialize observation file names, but don't use datestamp
  !
  call obsf_setup( dateStamp, 'bgck' )

  !
  !- Initialize the Temporal grid
  !
  call tim_setup( fileNameForDate_opt='./trlm_01' )

  !
  !- Initialize constants
  !
  if ( mpi_myid == 0 ) then
    call mpc_printConstants(6)
    call pre_printPrecisions
  end if

  !
  !- Initialize the Analysis grid
  !
  if( mpi_myid == 0 ) write(*,*)
  if( mpi_myid == 0 ) write(*,*) 'var_setup: Set hco parameters for analysis grid'
  call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

  if ( hco_anl % global ) then
    hco_core => hco_anl
  else
    !- Initialize the core (Non-Extended) analysis grid
    if( mpi_myid == 0) write(*,*)'var_setup: Set hco parameters for core grid'
    call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore' ) ! IN
  end if

  !     
  !- Initialisation of the analysis grid vertical coordinate from analysisgrid file
  !
  call vco_SetupFromFile(vco_anl,'./analysisgrid')

  call col_setVco(columnTrlOnAnlIncLev,vco_anl)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- Setup and read observations
  !
  call inn_setupObs(obsSpaceData, hco_anl, 'ALL', 'LIKESPLITFILES', 'bgck')
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! if ssmis, compute the surface type ele and update obspacedata

  call ssbg_computeSsmisSurfaceType(obsSpaceData)

  !
  !- Basic setup of columnData module
  !
  call col_setup
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- Memory allocation for background column data
  !
  call col_allocate(columnTrlOnAnlIncLev,obs_numheader(obsSpaceData),mpiLocal_opt=.true.)

  !
  !- Initialize the observation error covariances
  !
  call oer_setObsErrors(obsSpaceData, 'bgck')

  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  ! Initialize list of analyzed variables.
  !
  call gsv_setup

  ! Apply optional bias corrections
  call bcc_applyAIBcor(obsSpaceData)    
  call bcc_applyGPBcor(obsSpaceData)
    
  ! Reading trials
  call inn_getHcoVcoFromTrlmFile( hco_trl, vco_trl )
  allocHeightSfc = ( vco_trl%Vcode /= 0 )

  call gsv_allocate( stateVectorTrialHighRes, tim_nstepobs, hco_trl, vco_trl,  &
                     dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                     mpi_distribution_opt='Tiles', dataKind_opt=4,  &
                     allocHeightSfc_opt=allocHeightSfc, hInterpolateDegree_opt='LINEAR', &
                     beSilent_opt=.false. )
  call gsv_zero( stateVectorTrialHighRes )
  call gsv_readTrials( stateVectorTrialHighRes )
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Horizontally interpolate trials to trial columns
  call inn_setupColumnsOnTrlLev( columnTrlOnTrlLev, obsSpaceData, hco_core, &
                                   stateVectorTrialHighRes )
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Interpolate trial columns to analysis levels and setup for linearized H
  call inn_setupColumnsOnAnlIncLev( columnTrlOnTrlLev, columnTrlOnAnlIncLev )

  ! Compute observation innovations and prepare obsSpaceData for minimization
  call inn_computeInnovation(columnTrlOnTrlLev,obsSpaceData)

  ! 2.2 Perform the background check
  !     The routine also calls compute_HBHT and writes to listings & obsSpaceData

  ! Do the conventional data background check
  call bgck_bgCheck_conv(columnTrlOnAnlIncLev, columnTrlOnTrlLev, hco_anl, obsSpaceData)

  if (obs_famExist(obsSpaceData,'TO')) then

    ! Satellite radiance bias correction
    call bcs_calcBias(obsSpaceData,columnTrlOnTrlLev) ! Fill in OBS_BCOR obsSpaceData column with computed bias correction
    call bcs_applyBiasCorrection(obsSpaceData,OBS_VAR,'TO') ! Apply bias correction to OBS
    call bcs_applyBiasCorrection(obsSpaceData,OBS_OMP,'TO') ! Apply bias correction to O-F

    ! Do the TO background check
    call irbg_bgCheckIR(columnTrlOnTrlLev,obsSpaceData)
    call mwbg_bgCheckMW(obsSpaceData)
    call csrbg_bgCheckCSR(obsSpaceData)
    call ssbg_bgCheckSsmis(obsSpaceData)

  end if

  ! Do the ocean data background check
  if ( obs_famExist ( obsSpaceData, 'TM' )) call ocebg_bgCheckSST( obsSpaceData, columnTrlOnTrlLev, hco_trl )

  if (doThinning) then

    ! Copy original obs files into another directory
    call obsf_copyObsDirectory('./obsOriginal',direction='TO')

    ! 2.3 Write obs files after background check, but before thinning
    call obsf_writeFiles(obsSpaceData, writeDiagFiles_opt=.false.)
    
    ! Add cloud parameter data to burp files (AIRS,IASI,CrIS,ATMS,AMSUA,...)
    if (obs_famExist(obsSpaceData,'TO')) then
      call obsf_updateMissingObsFlags(obsSpaceData)
      call obsf_addCloudParametersAndEmissivity(obsSpaceData)
    end if

    ! Copy the pre-thinning files into another directory
    call obsf_copyObsDirectory('./obsBeforeThinning',direction='TO')

    ! Copy original obs files back into usual directory
    call obsf_copyObsDirectory('./obsOriginal',direction='FROM')

    ! 2.4 Thinning:  Set bit 11 of flag, one observation type at a time
    call thn_thinHyper(obsSpaceData)
    call thn_thinTovs(obsSpaceData)
    call thn_thinCSR(obsSpaceData)
    call thn_thinRaobs(obsSpaceData)
    call thn_thinScat(obsSpaceData)
    call thn_thinSatWinds(obsSpaceData)
    call thn_thinAircraft(obsSpaceData)
    call thn_thinSurface(obsSpaceData)
    call thn_thinGbGps(obsSpaceData)
    call thn_thinGpsRo(obsSpaceData)
    call thn_thinAladin(obsSpaceData)
    ! if requested, dump the thinned predictors and coefficients to sqlite
    call bcs_dumpBiasToSqliteAfterThinning(obsSpaceData)

  end if

  
  ! 3 Write the final results

  ! 3.1 Into the listings
  write(*,*)
  write(*,*) '> midas-obsSelection: printing the FIRST header and body'
  do headerIndex = 1, min(1,obs_numHeader(obsSpaceData))
    call obs_prnthdr(obsSpaceData,headerIndex)
    call obs_prntbdy(obsSpaceData,headerIndex)
  end do

  ! 3.2 Into the observation files
  write(*,*)
  write(*,*) '> midas-obsSelection: writing to file'
  call obsf_writeFiles(obsSpaceData)
  
  ! Add cloud parameter data to burp files (AIRS,IASI,CrIS,...)
  if (obs_famExist(obsSpaceData,'TO')) then
    call obsf_updateMissingObsFlags(obsSpaceData)
    call obsf_addCloudParametersAndEmissivity(obsSpaceData)
  end if

  ! cleaning the observation files
  if ( doThinning ) call obsf_cleanObsFiles()

  !
  ! 4.  Ending
  !
  write(*,*)
  write(*,*) '> midas-obsSelection: Ending'
  call obs_finalize(obsSpaceData) ! deallocate obsSpaceData

  call rpn_comm_finalize(ierr)

  call utl_tmg_stop(0)
  call tmg_terminate(mpi_myid, 'TMG_INFO')

end program midas_obsSelection

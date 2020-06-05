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
  !           (O-F => Observation minus Forecast)
  !
  use ramDisk_mod
  use mpi_mod
  use utilities_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use analysisGrid_mod
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
  implicit none

  integer :: datestamp, headerIndex, ierr
  type(struct_columnData),target :: trlColumnOnAnlLev
  type(struct_columnData),target :: trlColumnOnTrlLev
  type(struct_obs),       target :: obsSpaceData
  type(struct_hco), pointer      :: hco_anl => null()
  type(struct_vco), pointer      :: vco_anl => null()
  type(struct_hco), pointer      :: hco_core => null()

  integer :: get_max_rss

  write(*,*) ' -------------------------------------------------'
  write(*,*) ' ---  START OF MAIN PROGRAM midas-obsSelection ---'
  write(*,*) ' ---  Computation of the innovation            ---'
  write(*,*) ' -------------------------------------------------'

  !- 1.0 mpi
  call mpi_initialize

  !- 1.1 timings
  call tmg_init(mpi_myid, 'TMG_OBSSELECTION' )
  call tmg_start(1,'MAIN')

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
  if(mpi_myid.eq.0) call mpc_printConstants(6)

  !
  !- Initialize the Analysis grid
  !
  if(mpi_myid.eq.0) write(*,*)
  if(mpi_myid.eq.0) write(*,*) 'var_setup: Set hco parameters for analysis grid'
  call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

  if ( hco_anl % global ) then
    call agd_SetupFromHCO( hco_anl ) ! IN
  else
    !- Initialize the core (Non-Extended) analysis grid
    if(mpi_myid.eq.0) write(*,*)'var_setup: Set hco parameters for core grid'
    call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore' ) ! IN
    !- Setup the LAM analysis grid metrics
    call agd_SetupFromHCO(hco_anl, hco_core)
  end if

  !     
  !- Initialisation of the analysis grid vertical coordinate from analysisgrid file
  !
  call vco_SetupFromFile(vco_anl,'./analysisgrid')

  call col_setVco(trlColumnOnAnlLev,vco_anl)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- Setup and read observations
  !
  call inn_setupObs(obsSpaceData, 'ALL', 'LIKESPLITFILES', 'bgck')
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- Basic setup of columnData module
  !
  call col_setup
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- Memory allocation for background column data
  !
  call col_allocate(trlColumnOnAnlLev,obs_numheader(obsSpaceData),mpiLocal_opt=.true.)

  !
  !- Initialize the observation error covariances
  !
  call oer_setObsErrors(obsSpaceData, 'bgck')
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Apply optional bias corrections
  call bcc_applyAIBcor(obsSpaceData)    
  call bcc_applyGPBcor(obsSpaceData)
    
  ! Reading, horizontal interpolation and unit conversions of the 3D trial fields
  call inn_setupBackgroundColumns( trlColumnOnTrlLev, obsSpaceData )

  ! Interpolate trial columns to analysis levels and setup for linearized H
  call inn_setupBackgroundColumnsAnl(trlColumnOnTrlLev,trlColumnOnAnlLev,obsSpaceData)

  ! Compute observation innovations and prepare obsSpaceData for minimization
  call inn_computeInnovation(trlColumnOnTrlLev,obsSpaceData)

  ! 2.2 Perform the background check
  !     The routine also calls compute_HBHT and writes to listings & obsSpaceData

  ! Do the conventional data background check
  call bgck_bgcheck_conv(trlColumnOnAnlLev, trlColumnOnTrlLev, obsSpaceData)

  if (obs_famExist(obsSpaceData,'TO')) then

    ! Satellite radiance bias correction
    call bcs_calcBias(obsSpaceData,trlColumnOnTrlLev) ! Fill in OBS_BCOR obsSpaceData column with computed bias correction
    call bcs_applyBiasCorrection(obsSpaceData,OBS_VAR,'TO') ! Apply bias correction to OBS
    call bcs_applyBiasCorrection(obsSpaceData,OBS_OMP,'TO') ! Apply bias correction to O-F

    ! Do the IR background check
    call irbg_bgCheckIR(trlColumnOnTrlLev,obsSpaceData)

  end if

  ! 2.3 Thinning1:  Set bit 11 of flag, one observation type at a time
  call thn_thinAladin(obsSpaceData)

  ! 3 Write the results

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

  !  Add cloud parameter data to burp files (AIRS,IASI,CrIS,...)
  if (obs_famExist(obsSpaceData,'TO')) then
    call obsf_addCloudParametersAndEmissivity(obsSpaceData)
  end if

  ! Delete the flagged observations, and make the files smaller
  call obsf_thinFiles(obsSpaceData)

  !
  ! 4.  Ending
  !
  write(*,*)
  write(*,*) '> midas-obsSelection: Ending'
  call obs_finalize(obsSpaceData) ! deallocate obsSpaceData

  call rpn_comm_finalize(ierr)

  call tmg_stop(1)
  call tmg_terminate(mpi_myid, 'TMG_OMINUSF' )

end program midas_obsSelection

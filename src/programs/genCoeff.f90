!--------------------------------------- LICENCE BEGIN -----------------------------------
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

program midas_gencoeff
  !
  ! :Purpose: Main program to compute radiance bias correction coefficients by linear regression.
  !
  use version_mod
  use codePrecision_mod
  use ramDisk_mod
  use utilities_mod
  use midasMpi_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use obsSpaceData_mod
  use columnData_mod  
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use obsFiles_mod
  use innovation_mod
  use obsErrors_mod
  use biasCorrectionSat_mod

  implicit none

  integer, external :: exdb,exfin,fnom, fclos, get_max_rss
  integer :: ierr,istamp

  type(struct_obs),         target :: obsSpaceData
  type(struct_columnData),  target :: columnTrlOnAnlIncLev
  type(struct_gsv)                 :: stateVectorTrialHighRes
  type(struct_hco),        pointer :: hco_anl => null()
  type(struct_hco),        pointer :: hco_core => null()
  type(struct_vco),        pointer :: vco_anl => null()
  type(struct_hco),        pointer :: hco_trl => null()
  type(struct_vco),        pointer :: vco_trl => null()

  logical :: allocHeightSfc

  character(len=48), parameter :: obsMpiStrategy = 'LIKESPLITFILES', &
                                  varMode        = 'analysis'

  istamp = exdb('GENCOEFF','DEBUT','NON')

  call ver_printNameAndVersion('genCoeff','Bias Correction Coefficient Computation')

  ! MPI initialization
  call mmpi_initialize

  call tmg_init(mmpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')

 
  ! 1. Top level setup
  call ram_setup()
 
  ! Do initial set up
  call gencoeff_setup('VAR') ! obsColumnMode

  ! Reading trials
  call inn_getHcoVcoFromTrlmFile( hco_trl, vco_trl )
  allocHeightSfc = ( vco_trl%Vcode /= 0 )

  call gsv_allocate( stateVectorTrialHighRes, tim_nstepobs, hco_trl, vco_trl,  &
                     dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                     mpi_distribution_opt='Tiles', dataKind_opt=4,  &
                     allocHeightSfc_opt=allocHeightSfc, hInterpolateDegree_opt='LINEAR', &
                     beSilent_opt=.false. )
  call gsv_zero( stateVectorTrialHighRes )
  call gio_readTrials( stateVectorTrialHighRes )
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Horizontally interpolate trials to trial columns
  call inn_setupColumnsOnTrlLev( columnTrlOnAnlIncLev, obsSpaceData, hco_core, &
                                   stateVectorTrialHighRes )
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  call utl_tmg_start(110,'--BiasCorrection')

  ! Remove bias correction if requested
  call bcs_removeBiasCorrection(obsSpaceData,"TO")

  call bcs_removeOutliers(obsSpaceData)

  call utl_tmg_stop(110)

  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Compute observation innovations
  call inn_computeInnovation(columnTrlOnAnlIncLev,obsSpaceData)
  
  call utl_tmg_start(110,'--BiasCorrection')

  ! Refresh bias correction if requested
  call bcs_refreshBiasCorrection(obsSpaceData,columnTrlOnAnlIncLev)

  call utl_tmg_start(111,'----Regression')
  call bcs_do_regression(columnTrlOnAnlIncLev,obsSpaceData)
  call utl_tmg_stop(111)

  ! Write coefficients to file
  call bcs_writebias()

  ! output O-F statistics befor bias coorection
  call bcs_computeResidualsStatistics(obsSpaceData,"_raw")

  ! fill OBS_BCOR with computed bias correction
  call bcs_calcBias(obsSpaceData,columnTrlOnAnlIncLev)

  ! output  O-F statistics after bias coorection
  call bcs_computeResidualsStatistics(obsSpaceData,"_corrected")

  ! Deallocate internal bias correction structures 
  call bcs_finalize()

  call utl_tmg_stop(110)

  ! Deallocate copied obsSpaceData
  call obs_finalize(obsSpaceData)
  
  ! 3. Job termination

  istamp = exfin('GENCOEFF','FIN','NON')

  call utl_tmg_stop(0)

  call tmg_terminate(mmpi_myid, 'TMG_INFO')

  call rpn_comm_finalize(ierr) 

contains

  subroutine gencoeff_setup(obsColumnMode)
    !
    ! :Purpose:  Control of the preprocessing of bais correction coefficient computation
    !

    implicit none
    !Arguments:
    character(len=*), intent(in) :: obsColumnMode
    !Locals:	
    integer :: datestamp

    write(*,*) ''
    write(*,*) '----------------------------------------'
    write(*,*) '-- Starting subroutine gencoeff_setup --'
    write(*,*) '----------------------------------------'

    !
    !- Initialize the Temporal grid
    !
    call tim_setup

    !     
    !- Initialize observation file names and set datestamp
    !
    call obsf_setup( dateStamp, varMode )

    dateStamp = tim_getDatestampFromFile("./trlm_01")

    if ( dateStamp > 0 ) then
      call tim_setDatestamp(datestamp)     ! IN
    else
      call utl_abort('var_setup: Problem getting dateStamp from observation first trial field')
    end if

    !
    !- Initialize constants
    !
    if ( mmpi_myid == 0 ) then
      call mpc_printConstants(6)
      call pre_printPrecisions
    end if

    !
    !- Initialize variables of the model states
    !
    call gsv_setup
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Initialize the Analysis grid
    !
    if(mmpi_myid == 0) write(*,*)''
    if(mmpi_myid == 0) write(*,*)'gencoeff_setup: Set hco parameters for analysis grid'
    call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

    if ( hco_anl % global ) then
      hco_core => hco_anl
    else
      !- Initialize the core (Non-Extended) analysis grid
      if(mmpi_myid == 0) write(*,*)'gencoeff_setup: Set hco parameters for core grid'
      call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore' ) ! IN
    end if

    !     
    !- Initialisation of the analysis grid vertical coordinate from analysisgrid file
    !
    call vco_SetupFromFile( vco_anl,        & ! OUT
                            './analysisgrid') ! IN

    call col_setVco(columnTrlOnAnlIncLev,vco_anl)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Setup and read observations
    !
    call inn_setupObs(obsSpaceData, hco_anl, obsColumnMode, obsMpiStrategy, varMode) ! IN
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Basic setup of columnData module
    !
    call col_setup
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Initialize the observation error covariances
    !
    if (.not. bcs_mimicSatbcor) call oer_setObsErrors(obsSpaceData, varMode) ! IN
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

 end subroutine gencoeff_setup



end program midas_gencoeff

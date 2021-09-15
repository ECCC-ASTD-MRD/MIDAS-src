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
  use mpi_mod
  use mpivar_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use obsSpaceData_mod
  use columnData_mod  
  use gridStateVector_mod
  use obsSpaceDiag_mod
  use obsFiles_mod
  use obsFilter_mod  
  use innovation_mod
  use tovs_nl_mod
  use obsErrors_mod
  use statetocolumn_mod
  use biasCorrectionSat_mod
  use increment_mod
  use stateToColumn_mod
  use backgroundCheck_mod

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
  call mpi_initialize

  call tmg_init(mpi_myid, 'TMG_GENCOEFF' )

  call tmg_start(1,'MAIN')

 
  ! 1. Top level setup
  call ram_setup()
 
  ! Do initial set up
  call tmg_start(2,'SETUP')
  call gencoeff_setup('VAR') ! obsColumnMode
  call tmg_stop(2)

  call tmg_start(3,'TRIALS')

  ! Reading trials
  call gsv_getHcoVcoFromTrlmFile( hco_trl, vco_trl )
  allocHeightSfc = ( vco_trl%Vcode /= 0 )

  call gsv_allocate( stateVectorTrialHighRes, tim_nstepobs, hco_trl, vco_trl,  &
                     dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                     mpi_distribution_opt='Tiles', dataKind_opt=4,  &
                     allocHeightSfc_opt=allocHeightSfc, hInterpolateDegree_opt='LINEAR', &
                     beSilent_opt=.false. )
  call gsv_zero( stateVectorTrialHighRes )
  call gsv_readTrials( stateVectorTrialHighRes )
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Horizontally interpolate 15-min trials to trial columns
  call inn_setupColumnsOnTrialLev( columnTrlOnAnlIncLev, obsSpaceData, hco_core, &
                                   stateVectorTrialHighRes )
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  call tmg_stop(3)

  !
  ! Remove bias correction if requested
  !
  call tmg_start(4,'REMOVE_BCOR')
  call bcs_removeBiasCorrection(obsSpaceData,"TO")
  call tmg_stop(4)

  call tmg_start(5,'REMOVE_OUTLIERS')
  call bcs_removeOutliers(obsSpaceData)
  call tmg_stop(5)
   
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Compute observation innovations
  call tmg_start(6,'COMP_INOV')
  call inn_computeInnovation(columnTrlOnAnlIncLev,obsSpaceData)
  call tmg_stop(6)
  
  !
  ! Refresh bias correction if requested
  !
  call tmg_start(8,'REFRESH_BCOR')
  call bcs_refreshBiasCorrection(obsSpaceData,columnTrlOnAnlIncLev)
  call tmg_stop(8)

  call tmg_start(9,'REGRESSION')
  call bcs_do_regression(columnTrlOnAnlIncLev,obsSpaceData)
  call tmg_stop(9)

  ! Write coefficients to file
  call tmg_start(12,'WRITECOEFFS')
  call bcs_writebias()
  call tmg_stop(12)

  !
  ! output O-F statistics befor bias coorection
  !
  call tmg_start(13,'STATS')
  call bcs_computeResidualsStatistics(obsSpaceData,"_raw")
  call tmg_stop(13)
  !
  ! fill OBS_BCOR with computed bias correction
  !
  call tmg_start(15,'COMPBIAS')
  call bcs_calcBias(obsSpaceData,columnTrlOnAnlIncLev)
  call tmg_stop(15)

  !
  ! output  O-F statistics after bias coorection
  !
  call tmg_start(13,'STATS')
  call  bcs_computeResidualsStatistics(obsSpaceData,"_corrected")
  call tmg_stop(13)

  ! Deallocate internal bias correction structures 
  call bcs_finalize()

  ! Deallocate copied obsSpaceData
  call obs_finalize(obsSpaceData)
  
  ! 3. Job termination

  istamp = exfin('GENCOEFF','FIN','NON')

  call tmg_stop(1)

  call tmg_terminate(mpi_myid, 'TMG_GENCOEFF' )

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
    if ( mpi_myid == 0 ) then
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
    if(mpi_myid == 0) write(*,*)''
    if(mpi_myid == 0) write(*,*)'gencoeff_setup: Set hco parameters for analysis grid'
    call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

    if ( hco_anl % global ) then
      hco_core => hco_anl
    else
      !- Initialize the core (Non-Extended) analysis grid
      if(mpi_myid == 0) write(*,*)'gencoeff_setup: Set hco parameters for core grid'
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

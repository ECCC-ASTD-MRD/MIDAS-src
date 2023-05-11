program midas_gencoeff
  !
  !:Purpose: Main program to compute radiance bias correction coefficients by linear regression.
  !
  !          ---
  !
  !:Algorithm: O-A are computed from input analysis fields resulting from a stand alone
  !            3DVar analysis assimilating only anchoring observations (considered as "trials")
  !            and bgckalt files. 
  !         
  !            --
  !
  !:File I/O: The required input files and produced output files are listed as follows.
  !
  !============================================== ==============================================================
  ! Input and Output Files (NWP application)        Description of file
  !============================================== ==============================================================
  ! ``flnml``                                      In - Main namelist file with parameters user may modify
  ! ``flnml_static``                               In - The "static" namelist that should not be modified
  ! ``trlm_$NN`` (e.g. ``trlm_01``)                In - Background state (a.k.a. trial) files for each timestep
  ! ``analysisgrid``                               In - File defining grid for computing the analysis increment
  ! ``obsfiles_$FAM/obs$FAM_$NNNN_$NNNN``          In - Observation file for each "family" and MPI task
  ! Remainder are files related to radiance obs:
  ! ``stats_tovs``                                 In - Observation error file for radiances
  ! ``stats_tovs_symmetricObsErr``                 In - user-defined symmetric TOVS errors for all sky
  ! ``Cmat_$PLATFORM_$SENSOR.dat``                 In - Inter-channel observation-error correlations
  ! ``rtcoef_$PLATFORM_$SENSOR.H5``                In - RTTOV coefficient file HDF-5 format 
  ! ``rtcoef_$PLATFORM_$SENSOR.dat``               In - RTTOV coefficient file ASCII format 
  ! ``ozoneclim98``                                In - ozone climatology standard file (Fortuin and Kelder)
  !============================================== ==============================================================
  !
  !           --
  !
  !:Synopsis: Below is a summary of the ``genCoeff`` program calling sequence:
  !
  !             - **Initial setups:**
  !
  !               - Setup time grid using initial time from trlm_01
  !
  !               - Setup horizontal and vertical grid objects for "analysis
  !                 grid" from ``analysisgrid`` file.
  !
  !               - Setup ``obsSpaceData`` object and read observations from
  !                 files: ``inn_setupObs``.
  !
  !               - Setup ``columnData`` and ``gridStateVector`` modules (read
  !                 list of analysis variables from namelist) and allocate column
  !                 object for storing trial on analysis levels.
  !
  !               - Setup the observation error statistics in ``obsSpaceData``
  !                 object: ``oer_setObsErrors``.
  !
  !               - Allocate a stateVector object on the trial grid and then
  !                 read the trials: ``gio_readTrials``.
  !
  !             - **Coefficients computation**
  !
  !               - Horizontally Interpolate "trial" fields to trial columns
  !
  !               - if needed substract bias correction from ObsSpaceData to get raw O-F and obs.
  !
  !               - Remove outliers
  !
  !               - Compute innovation from "trial" fields
  !
  !               - Refresh Bias Correction (probably useless)
  !
  !               - Perform linear regression
  !
  !               - write coefficients to output file.
  !
  !               - if requested compute and output to file raw (.i.e without bias correction) O-F statistics
  !
  !               - compute and apply bias coorection to obsSpace Data
  !
  !               - if requested compute and output to file bias corrected O-F statistics
  !
  !             - **Final step**
  !
  !               - deallocate memory (radiance bias correction module and obsSpaceData)
  !
  !           --
  !
  !:Options: `List of namelist blocks <../namelists_in_each_program.html#genCoeff>`_
  !          that can affect the ``genCoeff`` program.
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

  ! output O-F statistics befor bias correction
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

    ! Arguments:
    character(len=*), intent(in) :: obsColumnMode

    ! Locals:	
    integer :: dateStamp, dateStampFromObs

    write(*,*) ''
    write(*,*) '----------------------------------------'
    write(*,*) '-- Starting subroutine gencoeff_setup --'
    write(*,*) '----------------------------------------'

    !
    !- Initialize the Temporal grid
    !
    call tim_setup

    !     
    !- Initialize observation file names and set datestamp from trial file
    !
    call obsf_setup(dateStampFromObs, varMode)

    dateStamp = tim_getDatestampFromFile('./trlm_01')

    if (dateStamp > 0) then
      call tim_setDatestamp(datestamp)
    else
      call utl_abort('var_setup: Problem getting dateStamp from first trial field')
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

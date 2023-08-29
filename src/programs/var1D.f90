program midas_var1D
  !
  !:Purpose: Main program for one dimensional variational minimization
  !
  !          ---
  !
  !:Algorithm: This program performs a similar data assimilation procedure as the var program,
  !            except without taking into account the horizontal or temporal dimensions.
  !            The assimilation is performed separately at each horizontal location and time where
  !            observations are present. The B matrix can be either an explicit representation of
  !            the covariances produced by the program extractBmatrixFor1Dvar or from an ensemble
  !            (controlled by the namelists). The resulting analysis and analysis increment at the
  !            observation locations/times are output in standard files on a "Y" grid. There is an 
  !            option to conduct idealized simulation by simulating the background and observations
  !            based on a prescribed true state, the background and observation error covariance 
  !            matrices.
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
  ! ``Bmatrix_sea.bin``                            In - 1DVar Bmatrix file over sea (output from extractBmatrixFor1DVar) 
  ! ``Bmatrix_land.bin``                           In - 1DVar Bmatrix file over land (output from extractBmatrixFor1DVar)
  ! ``obsfiles_$FAM/obs$FAM_$NNNN_$NNNN``          In - Observation file for each "family" and MPI task
  ! ``obserr``                                     In - Observation error statistics
  ! ``rttov_h2o_limits.dat``                       In - minimum and maximum humidity profile to clip analysis
  ! ``pm1q``                                       In/Out - Preconditioning file (Hessian of the cost function)
  ! ``rebm_$MMMm`` (e.g. ``rebm_180m``)            Out - Analysis increment on Y grid
  ! ``anlm_$MMMm``                                 Out - Analysis on Y grid
  ! ``obsfiles_$FAM.updated/obs$FAM_$NNNN_$NNNN``  Out - Updated obs file for each "family" and MPI task
  ! ``SimTrialOnTrlLev_$MMMm``                     Out - Simulated background on trial vert. level (Idealized Simulation)
  ! ``TruthOnTrlLev_$MMMm``                        Out - True state on trial vert. level (Idealized Simulation)
  ! ``TruthOnAnlLev_$MMMm``                        Out - True state on analysis vert. level grid (Idealized Simulation)
  ! ``PertOnTrlLev_$MMMm``                         Out - Simulated perturbation on trial level grid (Idealized Simulation)
  ! ``PertOnAnlLev_$MMMm``                         Out - Simulated perturbation on analysis level grid (Idealized Simulation)
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
  !:Synopsis: Below is a summary of the ``var1D`` program calling sequence:
  !
  !             - **Initial setups:**
  !
  !               - Initialize temporal grid
  !
  !               - Initialize observation file names and set date stamp
  !
  !               - Setup horizontal and vertical grid objects for "analysis
  !                 grid" from ``analysisgrid`` file
  !
  !               - Setup ``obsSpaceData`` object and read observations from
  !                 files: ``inn_setupObs``
  !
  !               - Setup ``columnData`` module
  !
  !               - Setup the observation error statistics in ``obsSpaceData``
  !                 object: ``oer_setObsErrors``
  !
  !               - Setup the gridStateVector module (initialize list of analyzed variables)
  !
  !               - Get horizontal and vertical grid descriptors from trial fields: ``inn_getHcoVcoFromTrlmFile``,
  !                 and allocate a gridStateVector objects 
  !
  !               - Read the trials: ``gio_readTrials``
  !
  !               - Setup the 1DVar B matrix: ``bmat1D_bsetup``
  !
  !               - Setup the ``gridVariableTransforms`` and ``minimization`` modules
  !
  !               - Horizontally interpolate high-resolution stateVectorUpdate to trial columns: ``inn_setupColumnsOnTrlLev``
  !
  !               - Simulate background if specified in NAM1DVAR ``var1DIdealize_simulateBackgroundState``
  !
  !               - Interpolate trial columns to analysis levels and setup for linearized H: ``inn_setupColumnsOnAnlIncLev``
  !
  !               - Compute observation innovations and prepare obsSpaceData for minimization: ``inn_computeInnovation``
  !
  !             - **Minimization:**
  !
  !               - Do the minimization ``min_minimize`` to obtain ``controlVectorIncr``
  !
  !               - Get 1DVar increment from ``controlVectorIncr`` to ``columnAnlInc``: ``bmat1D_get1DvarIncrement``
  !               - Transfer increment to ``stateVectorIncr``: ``var1D_transferColumnToYGrid``
  !
  !               - Write increment to file: ``inc_writeIncrement``
  !
  !             - **Final step:**
  !
  !               - Release resources
  !
  !               - Write hessian if requested to: ``min_writeHessian``
  !
  !               - Write updated observation files: ``obsf_writeFiles``
  !
  !           --
  !
  !:Options: `List of namelist blocks <../namelists_in_each_program.html#var1d>`_
  !          that can affect the ``var1D`` program.
  !
  !          - The B matrix used by the ``var1D`` program is controlled by the namelist block
  !            ``&nambmat1D`` of module bmatrix1DVar_mod
  !              * scaleFactorHI scaling factors for HI variances
  !              * scaleFactorHIHumidity scaling factors for HI humidity variances
  !              * scaleFactorENs scaling factors for Ens variances
  !              * scaleFactorEnsHumidity scaling factors for Ens humidity variances
  !              * nEns number of ensemble members for the ensemble part of  B matrix (negative means no ensemble)
  !              * vLocalize vertical localization length scale for the ensemble part of B matrix
  !              * includeAnlVar list of variable for the B matrix
  !              * numIncludeAnlVar number of variables in  includeAnlVar (to be removed later)
  !
  !           - The option to simulate the background and observations in idealized simulation controlled 
  !             by namelist NAM1DVAR
  !
  !           --
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
  use controlVector_mod
  use obsFiles_mod
  use minimization_mod
  use innovation_mod
  use obsErrors_mod
  use gridVariableTransforms_mod
  use increment_mod
  use biasCorrectionSat_mod
  use var1D_mod
  use bMatrix1Dvar_mod
  use var1DIdealize_mod
 
  implicit none

  integer :: istamp, exdb, exfin
  integer :: ierr, dateStampFromObs
  integer :: get_max_rss
  character(len=48) :: obsMpiStrategy, varMode
  real(8), allocatable :: controlVectorIncr(:)
  real(8), allocatable :: controlVectorIncrSum(:)
  type(struct_obs),        target :: obsSpaceData
  type(struct_columnData), target :: columnTrlOnAnlIncLev
  type(struct_columnData), target :: columnTrlOnTrlLev
  type(struct_columnData), target :: columnTrlOnTrlLevTruth
  type(struct_columnData), target :: columnAnlInc
  type(struct_gsv)                :: stateVectorIncr
  type(struct_gsv)                :: stateVectorTrialHighRes
  type(struct_gsv)                :: stateVectorAnalysis
  type(struct_hco),       pointer :: hco_anl => null()
  type(struct_vco),       pointer :: vco_anl => null()
  type(struct_hco),       pointer :: hco_core => null()
  type(struct_hco),       pointer :: hco_trl => null()
  type(struct_vco),       pointer :: vco_trl => null()

  integer :: outerLoopIndex, numIterMaxInnerLoop
  logical :: allocHeightSfc
  integer :: nulnam, fclos, fnom

  ! Namelist variables:
  logical :: simBgAndObs  ! Simulate Background and Observation
  integer :: simBgSeed    ! Random seed used to generate perturbation sampling 

  NAMELIST /NAM1DVAR/ simBgAndObs, simBgSeed

  istamp = exdb('VAR1D', 'DEBUT', 'NON')

  call ver_printNameAndVersion('var1D', '1D Variational Assimilation')

  ! MPI initialization
  call mmpi_initialize

  call tmg_init(mmpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')

  write(*,*)
  write(*,*) 'Real Kind used for computing the increment =', pre_incrReal
  write(*,*)

  varMode='analysis'

  ! Setup the ram disk
  call ram_setup

  ! Do initial set up

  !  Set/Read values for the namelist NAM1DVAR
  ! setting default values
  simBgAndObs = .False.
  simBgSeed = 0

  ! Check if NAM1DVAR exist
  if (.not. utl_isNamelistPresent('NAM1DVAR','./flnml')) then
    write(*,*)
    write(*,*) 'var1D: Namelist block NAM1DVAR is missing in the namelist.'
    write(*,*) '       the default values will be taken.'
  else
    ! Read the namelist
    nulnam=0
    ierr=fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    read(nulnam, nml=NAM1DVAR, iostat=ierr)
    if(ierr.ne.0) call utl_abort('midas-var1D: Error reading namelist')
    ierr=fclos(nulnam)
  end if

  if( mmpi_myid == 0 ) write(*,nml=NAM1DVAR)

  obsMpiStrategy = 'LIKESPLITFILES'

  ! Initialize the Temporal grid and set dateStamp from env variable
  call tim_setup

  ! Initialize observation file names and set datestamp if not already
  call obsf_setup(dateStampFromObs, varMode)
  if (tim_getDateStamp() == 0) then
    if (dateStampFromObs > 0) then
      call tim_setDatestamp(dateStampFromObs)
    else
      call utl_abort('midas-var1D: DateStamp was not set')
    end if
  end if

  ! Initialize constants
  if (mmpi_myid == 0) call mpc_printConstants(6)

  ! Initialize the Analysis grid
  if (mmpi_myid == 0) write(*,*)
  if (mmpi_myid == 0) write(*,*) 'midas-var1D: Set hco parameters for analysis grid'
  call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

  if ( hco_anl % global ) then
    hco_core => hco_anl
  else
    !- Initialize the core (Non-Extended) analysis grid
    if (mmpi_myid == 0) write(*,*) 'midas-var1D: Set hco parameters for core grid'
    call hco_SetupFromFile(hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore' ) ! IN
  end if

  ! Initialisation of the analysis grid vertical coordinate from analysisgrid file
  call vco_SetupFromFile(vco_anl,        & ! OUT
                         './analysisgrid') ! IN

  call col_setVco(columnTrlOnAnlIncLev, vco_anl)
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  ! Setup and read observations
  call inn_setupObs(obsSpaceData, hco_anl, 'VAR', obsMpiStrategy, varMode) ! IN
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  ! Basic setup of columnData module
  call col_setup
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Memory allocation for background column data
  call col_allocate(columnTrlOnAnlIncLev, obs_numheader(obsSpaceData))

  ! Initialize the observation error covariances
  call oer_setObsErrors(obsSpaceData, varMode) ! IN
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  ! Initialize list of analyzed variables.
  call gsv_setup
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Reading trials
  call inn_getHcoVcoFromTrlmFile(hco_trl, vco_trl)
  allocHeightSfc = ( vco_trl%Vcode /= 0 )

  call gsv_allocate(stateVectorTrialHighRes, tim_nstepobs, hco_trl, vco_trl,            &
                    dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true.,             &
                    mpi_distribution_opt='Tiles', dataKind_opt=pre_incrReal,            &
                    allocHeightSfc_opt=allocHeightSfc, hInterpolateDegree_opt='LINEAR', &
                    beSilent_opt=.false.)
  call gsv_zero(stateVectorTrialHighRes)
  call gio_readTrials(stateVectorTrialHighRes)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Initialize the background-error covariance, also sets up control vector module (cvm)
  call bmat1D_bsetup(vco_anl, hco_anl, obsSpaceData)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Initialize the gridded variable transform module
  call gvt_setup(hco_anl, hco_core, vco_anl)

  ! Set up the minimization module, now that the required parameters are known
  ! NOTE: some global variables remain in minimization_mod that must be initialized before
  !       inn_setupBackgroundColumns
  call min_setup(cvm_nvadim, hco_anl, oneDVarMode_opt=.true.) ! IN
  allocate(controlVectorIncr(cvm_nvadim),stat=ierr)
  if (ierr /= 0) then
    write(*,*) 'midas-var1D: Problem allocating memory for ''controlVectorIncr''',ierr
    call utl_abort('aborting in VAR1D')
  end if
  call utl_reallocate(controlVectorIncrSum,cvm_nvadim)
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  ! Horizontally interpolate high-resolution stateVectorUpdate to trial columns
  call inn_setupColumnsOnTrlLev(columnTrlOnTrlLev, obsSpaceData, hco_core, &
                                stateVectorTrialHighRes )

  ! Simulate the Background and Observation in idealized experiements
  if (simBgAndObs) then
    call col_setVco(columnTrlOnTrlLevTruth, col_getVco(columnTrlOnTrlLev))
    call col_allocate(columnTrlOnTrlLevTruth, col_getNumCol(columnTrlOnTrlLev), &
                      setToZero_opt=.true.)
    call col_copy(columnTrlOnTrlLev, columnTrlOnTrlLevTruth)
    call col_deallocate(columnTrlOnTrlLev)
    
    ! Simulate Background state columnTrlOnTrlLev
    call var1DIdealize_simulateBackgroundState(columnTrlOnTrlLevTruth, columnTrlOnTrlLev, &
                                                obsSpaceData, vco_anl, simBgSeed)
  end if  

  ! Interpolate trial columns to analysis levels and setup for linearized H
  call inn_setupColumnsOnAnlIncLev( columnTrlOnTrlLev,columnTrlOnAnlIncLev )

  ! Compute observation innovations and prepare obsSpaceData for minimization
  call inn_computeInnovation(columnTrlOnTrlLev, obsSpaceData)

  ! Do minimization of cost function
  outerLoopIndex = 1
  numIterMaxInnerLoop = 0
  controlVectorIncr(:) = 0.0d0
  call min_minimize( outerLoopIndex, columnTrlOnAnlIncLev, obsSpaceData, controlVectorIncrSum, &
                     controlVectorIncr, numIterMaxInnerLoop )
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Compute satellite bias correction increment and write to file
  ! Is is still necessary ? (will do nothing, but does it make sense in 1DVar mode ?)
  call bcs_writebias(controlVectorIncr)

  call col_setVco(columnAnlInc, col_getVco(columnTrlOnAnlIncLev))
  call col_allocate(columnAnlInc, col_getNumCol(columnTrlOnAnlIncLev))
  call col_zero(columnAnlInc)
  ! get final increment
  call bmat1D_get1DvarIncrement(controlVectorIncr, columnAnlInc, columnTrlOnAnlIncLev, obsSpaceData, cvm_nvadim)
  call var1D_transferColumnToYGrid(stateVectorIncr, obsSpaceData, columnAnlInc, bmat1D_includeAnlVar)

  ! output the analysis increment
  call inc_writeIncrement(stateVectorIncr)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! compute and write the analysis (as well as the increment on the trial grid)
  call var1d_transferColumnToYGrid(stateVectorAnalysis, obsSpaceData, columnTrlOnAnlIncLev, bmat1D_includeAnlVar)

  if (mmpi_myId == 0) call gsv_add(statevectorIncr, statevectorAnalysis)

  call inc_writeAnalysis(stateVectorAnalysis)

  ! Deallocate memory related to B matrices
  call var1D_finalize()
  if (mmpi_myId == 0) then
    call gsv_deallocate(stateVectorIncr)
    call gsv_deallocate(stateVectorAnalysis)
  end if

  ! write the Hessian
  call min_writeHessian(controlVectorIncr)
  deallocate(controlVectorIncrSum)
  deallocate(controlVectorIncr)

  ! Deallocate memory related to variational bias correction
  call bcs_finalize()
  ! Now write out the observation data files
  if (min_niter > 0) then
    if ( .not. obsf_filesSplit() ) then 
      write(*,*) 'We read/write global observation files'
      call obs_expandToMpiGlobal(obsSpaceData)
      if (mmpi_myid == 0) call obsf_writeFiles(obsSpaceData)
    else
      ! redistribute obs data to how it was just after reading the files
      call obs_MpiRedistribute(obsSpaceData, OBS_IPF)
      call obsf_writeFiles(obsSpaceData)
    end if
  end if

  ! Deallocate copied obsSpaceData
  call obs_finalize(obsSpaceData)

  ! Job termination
  istamp = exfin('VAR1D','FIN','NON')

  call utl_tmg_stop(0)

  call tmg_terminate(mmpi_myid, 'TMG_INFO')

  call rpn_comm_finalize(ierr) 

end program midas_var1D

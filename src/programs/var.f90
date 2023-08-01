program midas_var
  !
  !:Purpose: Main program for performing data assimilation using one of the
  !          following algorithms based on the incremental variational
  !          approach:
  !
  !            - 3D-Var
  !            - 3D- or 4D-EnVar
  !
  !          ---
  !
  !:Algorithm: The incremental variational data assimilation approach uses
  !            observations to compute a correction to the background state
  !            (i.e. the analysis increment) while taking into account the
  !            specified uncertainties for both the observations and the
  !            background state (i.e. the R and B covariance matrices,
  !            respectively). The analysis increment is computed by finding the
  !            minimum value of a cost function. This cost function is a scalar
  !            measure of the weighted departure of the corrected state from
  !            both the background state and the observations. The minimization
  !            is perfomed with a standard minimization algorithm (currently
  !            the quasi-Newton algorithm) that employs the gradient of the
  !            cost function to enable the minimum to be found (or at least
  !            approximated with sufficient accuracy) with relatively few
  !            iterations. The analysis increment is not assumed to be on the
  !            same horizontal and vertical grid as the background state. In
  !            addition, the temporal resolution of the background state and
  !            analysis increment are not assumed to be equal.
  !
  !            --
  !
  !            In practical terms, the background state is used to compute the
  !            "innovation" vector, that is, the difference between the
  !            observations and the background state in observation space:
  !            ``y-H(xb)``. The function ``H()`` represents the non-linear
  !            observation operators that map a gridded state vector into
  !            observation space. The innovation vector is used in the cost
  !            function calculation. Versions of the observation operators that
  !            are linearized with respect to the background state (or updated
  !            background state when using the outer loop) are also used in the
  !            cost function calculation to transform the analysis increment
  !            into observation space. The gradient of the cost function is
  !            computed using the adjoint of these linearized operators.
  !
  !            --
  !
  !            After the minimization, the analysis increment is spatially
  !            interpolated to the same grid as the background state and then
  !            added to this state to obtain the analysis. For some variables a
  !            physical constraint is imposed on the analysis
  !            (e.g. non-negative humidity or sea-ice concentration between 0
  !            and 1) and then the analysis increment is recomputed.
  !
  !            --
  !
  !:File I/O: The required input files and produced output files can vary
  !           according to the application. Below are tables of files for
  !           typical NWP 4D-EnVar (e.g. GDPS) and sea ice or SST 3D-Var
  !           applications.
  !
  !           --
  !
  !============================================== ==============================================================
  ! Input and Output Files (NWP applicaton)        Description of file
  !============================================== ==============================================================
  ! ``flnml``                                      In - Main namelist file with parameters user may modify
  ! ``flnml_static``                               In - The "static" namelist that should not be modified
  ! ``trlm_$NN`` (e.g. ``trlm_01``)                In - Background state (a.k.a. trial) files for each timestep
  ! ``analysisgrid``                               In - File defining grid for computing the analysis increment
  ! ``bgcov``                                      In - Static (i.e. NMC) B matrix file for NWP fields
  ! ``bgchemcov``                                  In - Static B matrix file for chemistry fields
  ! ``ensemble/$YYYYMMDDHH_006_$NNNN``             In - Ensemble member files defining ensemble B matrix
  ! ``obsfiles_$FAM/obs$FAM_$NNNN_$NNNN``          In - Observation file for each "family" and MPI task
  ! ``obserr``                                     In - Observation error statistics
  ! ``obsinfo_chm``                                In - Something needed for chemistry assimilation?
  ! ``preconin``                                   In - Preconditioning file (Hessian of the cost function)
  ! ``rebm_$MMMm`` (e.g. ``rebm_180m``)            Out - Analysis increment on the (low-res) analysis grid
  ! ``rehm_$MMMm``                                 Out - Analysis increment on the (high-res) trial grid
  ! ``anlm_$MMMm``                                 Out - Analysis on the trial grid
  ! ``obsfiles_$FAM.updated/obs$FAM_$NNNN_$NNNN``  Out - Updated obs file for each "family" and MPI task
  ! Remainder are files related to radiance obs:
  ! ``stats_$SENSOR_assim``                        In - Satellite radiance observation errors of difference sensors
  ! ``stats_tovs``                                 In - Satellite radiance observation errors
  ! ``stats_tovs_symmetricObsErr``                 In - User-defined symmetric TOVS errors for all sky
  ! ``Cmat_$PLATFORM_$SENSOR.dat``                 In - Inter-channel observation-error correlations
  ! ``ceres_global.std``                           In - High-res surface type and water fraction for radiance obs
  ! ``rtcoef_$PLATFORM_$SENSOR.dat``               In - RTTOV coefficient files
  ! ``rttov_h2o_limits.dat``                       In - Min/max humidity limits applied to analysis
  ! ``ozoneclim98``                                In - Ozone climatology
  !============================================== ==============================================================
  !
  !           --
  !
  !============================================== ==============================================================
  ! Input and Output Files (SST or Sea ice)        Description of file
  !============================================== ==============================================================
  ! ``flnml``                                      In - Main namelist file with parameters user may modify
  ! ``trlm_$NN`` (e.g. ``trlm_01``)                In - Background state (a.k.a. trial) files for each timestep
  ! ``analysisgrid``                               In - File defining grid for computing the analysis increment
  ! ``sea_ice_obs-err``                            In - Observation error statistics (sea ice only)
  ! ``bgstddev``                                   In - Background error stddev
  ! ``diffusmod.std``                              In - Normalization factor for diffusion operator correlations
  ! ``obsfiles_$FAM/obs$FAM_$NNNN_$NNNN``          In - Observation file for each "family" and MPI task
  ! ``rebm_$MMMm`` (e.g. ``rebm_180m``)            Out - Analysis increment on the (low-res) analysis grid
  ! ``rehm_$MMMm``                                 Out - Analysis increment on the (high-res) trial grid
  ! ``anlm_$MMMm``                                 Out - Analysis on the trial grid
  ! ``obsfiles_$FAM.updated/obs$FAM_$NNNN_$NNNN``  Out - Updated obs file for each "family" and MPI task
  !============================================== ==============================================================
  !
  !           --
  !
  !:Synopsis: Below is a summary of the ``var`` program calling sequence:
  !
  !             - **Initial setups:**
  !
  !               - Setup horizontal and vertical grid objects for "analysis
  !                 grid" from ``analysisgrid`` file and for "trial grid" from
  !                 first trial file: ``trlm_01``.
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
  !               - Setup the B matrices: ``bmat_setup``.
  !
  !               - Setup the ``gridVariableTransforms`` and ``minimization``
  !                 modules.
  !
  !             - **Outer loop** (only 1 iteration when no outer loop is used)
  !               during which the analysis increment is computed to correct the
  !               current "updated state" which is used to compute the
  !               innovation and serve as the starting point for the
  !               minimization of the following outer loop iteration:
  !
  !               - Impose RTTOV humidity limits on updated state (initially the
  !                 trial) on trial grid: ``qlim_rttovLimit``.
  !
  !               - Use the updated state on trial grid to setup the reference
  !                 state used by ``gridVariableTransforms`` module for height
  !                 calculations.
  !
  !               - Compute ``columnTrlOnTrlLev`` and ``columnTrlOnAnlIncLev`` from
  !                 updated state: ``inn_setupColumnsOnTrlLev``,
  !                 ``inn_setupColumnsOnAnlIncLev``
  !
  !               - Compute innovation from updated state:
  !                 ``inn_computeInnovation``.
  !
  !               - Use the updated state on trial grid to setup the reference
  !                 state used by ``gridVariableTransforms`` module for ``LQ``
  !                 to ``HU`` calculations.
  !
  !               - Do the minimization for this outer loop iteration:
  !                 ``min_minimize`` to obtain ``controlVectorIncr``.
  !
  !               - Update sum of all computed increment control vectors:
  !                 ``controlVectorIncrSum(:)`` which is needed in the
  !                 background cost function for outer loop.
  !
  !               - Compute ``stateVectorIncr`` (on analysis grid) from
  !                 ``controlVectorIncr``: ``inc_getIncrement``.
  !
  !               - Interpolate and add ``stateVectorIncr`` to the updated
  !                 state: ``inc_computeHighResAnalysis``.
  !
  !               - If requested, impose saturation and RTTOV humidity limits on
  !                 the updated state.
  !
  !               - Write increment (or sum of increments when outer loop used)
  !                 to file: ``inc_writeIncrement``.
  !
  !               - End of outer loop.
  !
  !             - **Final steps**, after the outer loop:
  !
  !               - If requested, compute final cost function value using
  !                 non-linear observation operators.
  !
  !               - Write the final analysis and recomputed complete increment on
  !                 the trial grid: ``inc_writeincAndAnalHighRes``.
  !
  !               - Various final steps, including: write the Hessian to binary
  !                 file (``min_writeHessian``), update the observation files
  !                 (``obsf_writeFiles``).
  !
  !           --
  !
  !:Options: `List of namelist blocks <../namelists_in_each_program.html#var>`_
  !          that can affect the ``var`` program.
  !
  !          * The use of an outer loop is controlled by the namelist block
  !            ``&NAMVAR`` read by the ``var`` program.
  !
  !          * The choice of 3D-Var vs. EnVar is controlled by the weights given
  !            to the climatological (i.e. static) and ensemble-based background
  !            error covariance (i.e. B matrix) components. The weights for all
  !            B matrix components are zero be default and can be set to a
  !            nonzero value through the namelist variable ``SCALEFACTOR`` in
  !            the namelist block for each corresponding fortran module.
  !
  !          * Either algorithm can be used in combination with "First guess at
  !            appropriate time" (i.e. FGAT) which is controlled by the namelist
  !            variable ``DSTEPOBS`` (in ``NAMTIME``) that specifies the length
  !            of time (in hours) between times when the background state is
  !            used to compute the innovation (i.e. O-B).
  !
  !          * Similarly, the choice between 3D-EnVar and 4D-EnVar is controlled
  !            by the namelist variable ``DSTEPOBSINC`` (also in ``&NAMTIME``)
  !            which specifies the length of time (in hours) between times when
  !            the analysis increment is computed. In the context of EnVar, if
  !            this is less than the assimilation time window, then the
  !            ensembles used to construct the ensemble-based B matrix component
  !            will be used at multiple times within the assimilation window to
  !            obtain 4D covariances (i.e. 4D-EnVar).
  !
  !          * Some of the other relevant namelist blocks used to configure the
  !            variational analysis are listed in the following table:
  ! 
  !======================== ============ ==============================================================
  ! Module                   Namelist     Description of what is controlled
  !======================== ============ ==============================================================
  ! ``minimization_mod``     ``NAMMIN``   maximum number of iterations, convergence criteria and many
  !                                       additional parameters for controlling the minimization
  ! ``timeCoord_mod``        ``NAMTIME``  assimilation time window length, temporal resolution of
  !                                       the background state and increment
  ! ``bMatrixEnsemble_mod``  ``NAMBEN``   weight and other parameters for ensemble-based B matrix
  !                                       component
  ! ``bMatrixHI_mod``        ``NAMBHI``   weight and other parameters for the climatological B matrix
  !                                       component based on homogeneous-isotropic covariances
  !                                       represented in spectral space
  ! Other B matrix modules   various      weight and other parameters for each type of B matrix
  !======================== ============ ==============================================================
  !
  use version_mod
  use codePrecision_mod
  use ramDisk_mod
  use utilities_mod
  use midasMpi_mod
  use message_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use humidityLimits_mod
  use timeCoord_mod
  use obsSpaceData_mod
  use columnData_mod  
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use obsSpaceDiag_mod
  use controlVector_mod
  use obsFiles_mod
  use minimization_mod
  use innovation_mod
  use bMatrix_mod
  use rMatrix_mod
  use obsErrors_mod
  use gridVariableTransforms_mod
  use increment_mod
  use biasCorrectionSat_mod
  use varQC_mod
  use tovsNL_mod
  use stateToColumn_mod

  implicit none

  integer :: istamp, exdb, exfin
  integer :: ierr, dateStampFromObs, nulnam
  integer :: fclos, fnom
  character(len=9)  :: clmsg
  character(len=48) :: obsMpiStrategy, varMode
  real(8), allocatable :: controlVectorIncr(:)
  real(8), allocatable :: controlVectorIncrSum(:)

  type(struct_obs)       , target :: obsSpaceData
  type(struct_columnData), target :: columnTrlOnAnlIncLev
  type(struct_columnData), target :: columnTrlOnTrlLev
  type(struct_gsv)                :: stateVectorIncr
  type(struct_gsv)                :: stateVectorIncrSum
  type(struct_gsv)                :: stateVectorUpdateHighRes
  type(struct_gsv)                :: stateVectorTrial
  type(struct_gsv)                :: stateVectorPsfcHighRes
  type(struct_gsv)                :: stateVectorPsfc
  type(struct_gsv)                :: stateVectorAnal
  type(struct_hco)      , pointer :: hco_anl => null()
  type(struct_vco)      , pointer :: vco_anl => null()
  type(struct_hco)      , pointer :: hco_trl => null()
  type(struct_vco)      , pointer :: vco_trl => null()
  type(struct_hco)      , pointer :: hco_core => null()

  integer :: outerLoopIndex, numIterMaxInnerLoopUsed
  integer :: numIterWithoutVarqc, numInnerLoopIterDone

  logical :: allocHeightSfc, applyLimitOnHU
  logical :: deallocHessian, isMinimizationFinalCall
  logical :: varqcActive, applyVarqcOnNlJo
  logical :: filterObsAndInitOer
  logical :: deallocInterpInfoNL

  integer, parameter :: maxNumOuterLoopIter = 15

  ! namelist variables
  integer :: numOuterLoopIterations                    ! number of outer loop iterations (default=1)
  integer :: numIterMaxInnerLoop(maxNumOuterLoopIter)  ! number of each inner loop iterations
  logical :: limitHuInOuterLoop                        ! impose humidity limits on each outer loop iteration
  logical :: computeFinalNlJo                          ! compute final cost function using non-linear H()
  logical :: useTovsUtil                               ! do channel filtering based on UTIL column ot the stats_tovs file 
  NAMELIST /NAMVAR/ numOuterLoopIterations, numIterMaxInnerLoop, limitHuInOuterLoop
  NAMELIST /NAMVAR/ computeFinalNlJo, useTovsUtil

  istamp = exdb('VAR','DEBUT','NON')

  call ver_printNameAndVersion('var','Variational Assimilation')

  ! MPI initialization
  call mmpi_initialize

  call tmg_init(mmpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')

  if (mmpi_myid == 0) then
    clmsg = 'VAR3D_BEG'
    call utl_writeStatus(clmsg)
  end if 

  varMode='analysis'

  ! Setup the ram disk
  call ram_setup

  ! Do initial set up

  ! Set/Read values for the namelist NAMVAR
  ! Setting default namelist variable values
  numOuterLoopIterations = 1
  limitHuInOuterLoop = .false.
  numIterMaxInnerLoop(:) = 0
  computeFinalNlJo = .false.
  useTovsUtil = .false.

  if ( .not. utl_isNamelistPresent('NAMVAR','./flnml') ) then
  call msg('midas-var','namvar is missing in the namelist. '&
       //'The default values will be taken.', mpiAll_opt=.false.)

  else
    ! read in the namelist NAMVAR
    nulnam = 0
    ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    read(nulnam, nml=namvar, iostat=ierr)
    if( ierr /= 0) call utl_abort('midas-var: Error reading namelist')
    ierr = fclos(nulnam)
  end if
  if ( mmpi_myid == 0 ) write(*,nml=namvar)

  if ( numOuterLoopIterations > maxNumOuterLoopIter ) then
    call utl_abort('midas-var: numOuterLoopIterations is greater than max value')
  end if

  if ( numOuterLoopIterations > 1 .and. &
       .not. all(numIterMaxInnerLoop(1:numOuterLoopIterations) > 0) ) then
    call utl_abort('midas-var: some numIterMaxInnerLoop(:) in namelist are negative or zero')
  end if

  obsMpiStrategy = 'LIKESPLITFILES'

  ! Initialize the Temporal grid and set dateStamp from env variable
  call tim_setup()

  ! Initialize observation file names and set datestamp if not already
  call obsf_setup(dateStampFromObs, varMode)
  if (tim_getDateStamp() == 0) then
    if (dateStampFromObs > 0) then
      call tim_setDatestamp(datestampFromObs)
    else
      call utl_abort('var_setup: DateStamp was not set')
    end if
  end if

  ! Initialize constants
  if ( mmpi_myid == 0 ) then
    call mpc_printConstants(6)
    call pre_printPrecisions
  end if

  ! Initialize the Analysis grid
  call msg('midas-var','Set hco parameters for analysis grid', mpiAll_opt=.false.)
  call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

  if ( hco_anl % global ) then
    hco_core => hco_anl
  else
    ! Initialize the core (Non-Extended) analysis grid
    call msg('midas-var','Set hco parameters for core grid', mpiAll_opt=.false.)
    call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore' ) ! IN
  end if

  ! Initialisation of the analysis grid vertical coordinate from analysisgrid file
  call vco_SetupFromFile( vco_anl,        & ! OUT
                          './analysisgrid') ! IN

  call col_setVco(columnTrlOnAnlIncLev,vco_anl)
  call msg_memUsage('var')

  ! Setup and read observations
  call inn_setupObs(obsSpaceData, hco_anl, 'VAR', obsMpiStrategy, varMode) ! IN
  call msg_memUsage('var')

  ! Basic setup of columnData module
  call col_setup
  call msg_memUsage('var')

  !- Memory allocation for background column data
  call col_allocate(columnTrlOnAnlIncLev,obs_numheader(obsSpaceData),mpiLocal_opt=.true.)

  ! Initialize the observation error covariances
  call oer_setObsErrors(obsSpaceData, varMode, useTovsUtil_opt=useTovsUtil) ! IN

  
   ! Call suprep again to filter out channels according to 'util' column of stats_tovs
  if (useTovsUtil) call filt_suprep(obsSpaceData)
  call msg_memUsage('var')

  ! Initialize list of analyzed variables.
  call gsv_setup
  call msg_memUsage('var')

  ! Reading trials
  call inn_getHcoVcoFromTrlmFile( hco_trl, vco_trl )
  allocHeightSfc = ( vco_trl%Vcode /= 0 )

  call gsv_allocate( stateVectorUpdateHighRes, tim_nstepobs, hco_trl, vco_trl,  &
                     dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                     mpi_distribution_opt='Tiles', dataKind_opt=pre_incrReal,  &
                     allocHeightSfc_opt=allocHeightSfc, hInterpolateDegree_opt='LINEAR')
  call gsv_zero( stateVectorUpdateHighRes )
  call gio_readTrials( stateVectorUpdateHighRes )
  call msg_memUsage('var')

  ! Initialize the background-error covariance, also sets up control vector module (cvm)
  call bmat_setup(hco_anl,hco_core,vco_anl)
  call msg_memUsage('var')

  ! Initialize the gridded variable transform module
  call gvt_setup(hco_anl,hco_core,vco_anl)

  ! Set up the minimization module, now that the required parameters are known
  ! NOTE: some global variables remain in minimization_mod that must be initialized before
  !       inn_setupColumnsOnTrlLev
  call min_setup( cvm_nvadim, hco_anl,                                   & ! IN
                  varqc_opt=varqcActive, nwoqcv_opt=numIterWithoutVarqc )  ! OUT
  allocate(controlVectorIncr(cvm_nvadim),stat=ierr)
  if (ierr /= 0) then
    call msg('var','Problem allocating memory for controlVectorIncr'//str(ierr))
    call utl_abort('aborting in VAR')
  end if
  call utl_reallocate(controlVectorIncrSum,cvm_nvadim)
  call msg_memUsage('var')

  numInnerLoopIterDone = 0

  ! Enter outer-loop
  outer_loop: do outerLoopIndex = 1, numOuterLoopIterations
    call msg('var','start of outer-loop index='//str(outerLoopIndex))

    ! Impose limits on ALL cloud variables
    call qlim_rttovLimit(stateVectorUpdateHighRes, applyLimitToCloud_opt=.true.)

    ! Initialize stateVectorRefHeight for transforming TT/HU/P0 increments to
    ! height/pressure increments.
    if ( (gsv_varExist(stateVectorUpdateHighRes,'P_T') .and. &
          gsv_varExist(stateVectorUpdateHighRes,'P_M')) .or. &
         (gsv_varExist(stateVectorUpdateHighRes,'Z_T') .and. &
          gsv_varExist(stateVectorUpdateHighRes,'Z_M')) ) then

      call gvt_setupRefFromStateVector( stateVectorUpdateHighRes, 'height' )

      call msg_memUsage('var')
    end if

    ! Horizontally interpolate high-resolution stateVectorUpdate to trial columns
    deallocInterpInfoNL = ( numOuterLoopIterations <= 1 )
    call inn_setupColumnsOnTrlLev( columnTrlOnTrlLev, obsSpaceData, hco_core, &
                                   stateVectorUpdateHighRes, &
                                   deallocInterpInfoNL_opt=deallocInterpInfoNL )
    call msg_memUsage('var')

    ! Interpolate trial columns to analysis levels and setup for linearized H
    call inn_setupColumnsOnAnlIncLev( columnTrlOnTrlLev, columnTrlOnAnlIncLev )
    call msg_memUsage('var')

    ! Determine if to apply varqc to Jo of non-linear operator
    applyVarqcOnNlJo = ( varqcActive .and. numInnerLoopIterDone > numIterWithoutVarqc )
    if ( applyVarqcOnNlJo ) call msg('var','applying varqc to non-linear Jo', mpiAll_opt=.false.)

    ! Compute observation innovations and prepare obsSpaceData for minimization
    filterObsAndInitOer = ( outerLoopIndex == 1 )
    call inn_computeInnovation( columnTrlOnTrlLev, obsSpaceData, &
                                filterObsAndInitOer_opt=filterObsAndInitOer, &
                                applyVarqcOnNlJo_opt=applyVarqcOnNlJo, &
                                callSetErrGpsgb_opt=filterObsAndInitOer )
    call msg_memUsage('var')

    ! Initialize stateVectorRefHU for doing variable transformation of the increments.
    if ( gsv_varExist(stateVectorUpdateHighRes,'HU') ) then
      applyLimitOnHU = ( limitHuInOuterLoop .and. outerLoopIndex > 1 )

      call gvt_setupRefFromStateVector( stateVectorUpdateHighRes, 'HU', &
                                        applyLimitOnHU_opt=applyLimitOnHU )

      call msg_memUsage('var')
    end if

    ! Do minimization of cost function. Use numIterMaxInnerLoop from NAMVAR, instead of
    ! nitermax from NAMMIN, when numOuterLoopIterations > 1
    controlVectorIncr(:) = 0.0d0
    deallocHessian = ( numOuterLoopIterations == 1 )
    isMinimizationFinalCall = ( outerLoopIndex == numOuterLoopIterations )
    call min_minimize( outerLoopIndex, columnTrlOnAnlIncLev, obsSpaceData, controlVectorIncrSum, &
                       controlVectorIncr, numIterMaxInnerLoop(outerLoopIndex), &
                       deallocHessian_opt=deallocHessian, &
                       isMinimizationFinalCall_opt=isMinimizationFinalCall, &
                       numIterMaxInnerLoopUsed_opt=numIterMaxInnerLoopUsed )
    numInnerLoopIterDone = numInnerLoopIterDone + numIterMaxInnerLoopUsed
    call msg_memUsage('var')

    ! Accumulate control vector increments of all the previous iterations
    controlVectorIncrSum(:) = controlVectorIncrSum(:) + controlVectorIncr(:)

    ! Compute satellite bias correction increment and write to file on last outer-loop 
    ! iteration
    if ( outerLoopIndex == numOuterLoopIterations ) then
      call bcs_writebias(controlVectorIncr)
    end if

    call tvs_deallocateProfilesNlTlAd

    call gsv_allocate(stateVectorIncr, tim_nstepobsinc, hco_anl, vco_anl, &
         datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true., &
         dataKind_opt=pre_incrReal, allocHeight_opt=.false., allocPressure_opt=.false.)

    ! get final increment with mask if it exists
    call inc_getIncrement( controlVectorIncr, stateVectorIncr, cvm_nvadim )
    call gio_readMaskFromFile( stateVectorIncr, './analysisgrid' )
    call msg_memUsage('var')

    ! Compute high-resolution analysis on trial grid
    call inc_computeHighResAnalysis( stateVectorIncr,                                  & ! IN
                                     stateVectorUpdateHighRes, stateVectorPsfcHighRes )  ! OUT
    call msg_memUsage('var')

    ! Impose limits on stateVectorUpdateHighRes only when outer loop is used.
    if ( limitHuInOuterLoop ) then
      call msg('var','impose limits on stateVectorUpdateHighRes')
      call qlim_saturationLimit( stateVectorUpdateHighRes )
      call qlim_rttovLimit( stateVectorUpdateHighRes )
    end if

    ! prepare to write incremnt when no outer-loop, or sum of increments at last
    ! outer-loop iteration.
    if ( numOuterLoopIterations > 1 .and. &
         outerLoopIndex == numOuterLoopIterations ) then

      call gsv_allocate(stateVectorIncrSum, tim_nstepobsinc, hco_anl, vco_anl, &
           datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true., &
           dataKind_opt=pre_incrReal, allocHeight_opt=.false., allocPressure_opt=.false.)
      call inc_getIncrement( controlVectorIncrSum, stateVectorIncrSum, cvm_nvadim )
      call gio_readMaskFromFile( stateVectorIncrSum, './analysisgrid' )
      call msg_memUsage('var')

      call inc_writeIncrement( stateVectorIncrSum, &     ! IN
                               ip3ForWriteToFile_opt=0 ) ! IN
      call gsv_deallocate( stateVectorIncrSum )
    else if ( numOuterLoopIterations == 1 ) then
      call inc_writeIncrement( stateVectorIncr, &        ! IN
                               ip3ForWriteToFile_opt=0 ) ! IN
    end if
    call msg_memUsage('var')

    call gsv_deallocate( stateVectorIncr )

    call msg('var','end of outer-loop index='//str(outerLoopIndex))
  end do outer_loop

  ! Set the QC flags to be consistent with VAR-QC if control analysis
  if ( varqcActive ) call vqc_listrej(obsSpaceData)

  if ( computeFinalNlJo ) then
    ! Horizontally interpolate high-resolution stateVectorUpdate to trial columns
    call inn_setupColumnsOnTrlLev( columnTrlOnTrlLev, obsSpaceData, hco_core, &
                                   stateVectorUpdateHighRes )
    call msg_memUsage('var')

    ! Determine if to apply varqc to Jo of non-linear operator
    applyVarqcOnNlJo = ( varqcActive .and. numInnerLoopIterDone > numIterWithoutVarqc )
    if ( applyVarqcOnNlJo ) call msg('var','applying varqc to non-linear Jo', mpiAll_opt=.false.)

    ! Compute observation innovations and prepare obsSpaceData for minimization
    filterObsAndInitOer = .false.
    call inn_computeInnovation( columnTrlOnTrlLev, obsSpaceData, &
                                destObsColumn_opt=OBS_OMA, &
                                filterObsAndInitOer_opt=filterObsAndInitOer, &
                                applyVarqcOnNlJo_opt=applyVarqcOnNlJo , &
                                callSetErrGpsgb_opt=filterObsAndInitOer )
    call msg_memUsage('var')
  end if

  ! Memory deallocations for non diagonal R matrices for radiances
  call rmat_cleanup()

  ! Conduct obs-space post-processing diagnostic tasks (some diagnostic
  ! computations controlled by NAMOSD namelist in flnml)
  call osd_ObsSpaceDiag( obsSpaceData, columnTrlOnAnlIncLev, hco_anl )

  ! Deallocate memory related to B matrices and update stateVector
  call bmat_finalize()

  ! Deallocate structures needed for interpolation
  call s2c_deallocInterpInfo( inputStateVectorType='nl' )
  call s2c_deallocInterpInfo( inputStateVectorType='tlad' )

  ! Post processing of analyis before writing (variable transform+humidity clipping)
  call inc_analPostProcessing( stateVectorPsfcHighRes, stateVectorUpdateHighRes, &  ! IN 
                               stateVectorTrial, stateVectorPsfc, stateVectorAnal ) ! OUT
  call gsv_deallocate( stateVectorUpdateHighRes )

  ! compute and write the analysis (as well as the increment on the trial grid)
  call inc_writeIncAndAnalHighRes( stateVectorTrial, stateVectorPsfc, &
                                   stateVectorAnal )
  call msg_memUsage('var')

  if (mmpi_myid == 0) then
    clmsg = 'REBM_DONE'
    call utl_writeStatus(clmsg)
  end if

  ! write the Hessian
  call min_writeHessian(controlVectorIncr)
  deallocate(controlVectorIncr)

  ! Deallocate memory related to variational bias correction
  call bcs_finalize()

  ! Now write out the observation data files
  if (min_niter > 0) then
    if ( .not. obsf_filesSplit() ) then 
      call msg('var','reading/writing global observation files')
      call obs_expandToMpiGlobal(obsSpaceData)
      if (mmpi_myid == 0) call obsf_writeFiles(obsSpaceData)
    else
      ! redistribute obs data to how it was just after reading the files
      call obs_MpiRedistribute(obsSpaceData,OBS_IPF)
      call obsf_writeFiles(obsSpaceData)
    end if
  end if

  ! Deallocate copied obsSpaceData
  call obs_finalize(obsSpaceData)

  ! Job termination
  istamp = exfin('VAR','FIN','NON')

  if (mmpi_myid == 0) then
    clmsg = 'VAR3D_END'
    call utl_writeStatus(clmsg)
  end if

  call utl_tmg_stop(0)

  call tmg_terminate(mmpi_myid, 'TMG_INFO')

  call rpn_comm_finalize(ierr) 

end program midas_var

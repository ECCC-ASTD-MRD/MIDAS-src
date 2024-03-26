program midas_dfs
  !
  !:Purpose: Main program for computing a single random realization of
  !          background error in observation space, stored in ``obs_hbht``. This
  !          can then be used by python or other scripts to compute the
  !          background error variance (consistent with the specified B matrix)
  !          in observation space for comparison with the innovation variance
  !          and observation error variance.
  !
  !          --
  !
  !:Algorithm: The random realization of background error in observation space
  !            is computed following these steps:
  !
  !            1. Compute random values for the control vector with each element
  !               drawn from independent Gaussian distribution with variance of one
  !               and bias of zero.
  !            
  !            2. Multiply random vector by sqrt of B matrix.
  !
  !            3. Apply observation operator to obtain random perturbation in
  !               observation space.
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
  ! Input and Output Files (for NWP application)   Description of file
  !============================================== ==============================================================
  ! ``flnml``                                      In - Main namelist file with parameters user may modify
  ! ``flnml_static``                               In - The "static" namelist that should not be modified
  ! ``trlm_$NN`` (e.g. ``trlm_01``)                In - Background state (a.k.a. trial) files for each timestep
  ! ``analysisgrid``                               In - File defining grid for computing the random gridded perturbation
  ! ``bgcov``                                      In - Static (i.e. NMC) B matrix file for NWP fields
  ! ``bgchemcov``                                  In - Static B matrix file for chemistry fields
  ! ``ensemble/$YYYYMMDDHH_006_$NNNN``             In - Ensemble member files defining ensemble B matrix
  ! ``obsfiles_$FAM/obs$FAM_$NNNN_$NNNN``          In - Observation file for each "family" and MPI task
  ! ``obserr``                                     In - Observation error statistics
  ! ``obsinfo_chm``                                In - Something needed for chemistry assimilation?
  ! ``obsfiles_$FAM.updated/obs$FAM_$NNNN_$NNNN``  Out - Updated obs file for each "family" and MPI task
  ! Remainder are files related to radiance obs:
  ! ``stats_$SENSOR_assim``                        In - Satellite radiance observation errors of different sensors
  ! ``stats_tovs``                                 In - Satellite radiance observation errors
  ! ``stats_tovs_symmetricObsErr``                 In - User-defined symmetric TOVS errors for all sky
  ! ``ceres_global.std``                           In - High-res surface type and water fraction for radiance obs
  ! ``rtcoef_$PLATFORM_$SENSOR.dat``               In - RTTOV coefficient files
  ! ``ozoneclim98``                                In - Ozone climatology
  !============================================== ==============================================================
  !
  !           --
  !
  !:Synopsis: Below is a summary of the ``diagHBHt`` program calling sequence:
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
  !               - Setup the ``gridVariableTransforms``.
  !
  !             - **Calculation**
  !
  !               - Compute ``columnTrlOnTrlLev`` and ``columnTrlOnAnlIncLev``
  !                 from trials: ``inn_setupColumnsOnTrlLev``,
  !                 ``inn_setupColumnsOnAnlIncLev``
  !
  !               - Compute innovation from updated state:
  !                 ``inn_computeInnovation``.
  !
  !               - Compute an MPI global random vector, then extract only
  !                 portion needed for this MPI task (to reduce sensitivity of
  !                 results to MPI topology).
  !
  !               - Multiply random vector by sqrt of B matrix with resulting
  !                 gridded state random perturbation in ``statevector``.
  !
  !               - Apply linearized observation operators to the random gridded
  !                 state: ``s2c_tl`` and ``oop_Htl`` with final result in
  !                 observation space: ``obs_work`` column of ``obsSpaceData``.
  !
  !               - Copy result from ``obs_work`` to ``obs_hbht`` column.
  !
  !             - **Final steps**, after the outer loop:
  !
  !               - Various final steps, including: update the observation files
  !                 (``obsf_writeFiles``).
  !
  !           --
  !
  !:Options: `List of namelist blocks <../namelists_in_each_program.html#diaghbht>`_
  !          that can affect the ``diagHBHt`` program.
  !
  !          * The choice of what B matrix is used for the calculation is
  !            controlled for each individual B matrix component through it's
  !            own namelist block. The weights for all B matrix components are
  !            zero be default and can be set to a nonzero value through the
  !            namelist variable ``SCALEFACTOR`` in the namelist block for each
  !            corresponding fortran module.
  !
  !          * All other namelist blocks related to observations are relevant
  !            for the diagHBHt calculation, including ``NAMFILT`` and
  !            ``NAMTOV``.
  !
  !          * Some of the other relevant namelist blocks used to configure the
  !            diagHBHt calculation are listed in the following table:
  ! 
  !======================== ============ ==============================================================
  ! Module                   Namelist     Description of what is controlled
  !======================== ============ ==============================================================
  ! ``timeCoord_mod``        ``NAMTIME``  assimilation time window length, temporal resolution of
  !                                       the background state and increment (i.e. perturbation)
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
  use MathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use obsSpaceData_mod
  use columnData_mod  
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use controlVector_mod
  use obsFiles_mod
  use obsTimeInterp_mod
  use stateToColumn_mod
  use innovation_mod
  use bMatrix_mod
  use obsErrors_mod
  use gridVariableTransforms_mod
  use obsOperators_mod
  use tovsNL_mod
  use rMatrix_mod
  
  implicit none

  integer, external :: exdb, exfin
  integer :: ierr, istamp

  type(struct_obs),        target :: obsSpaceData
  type(struct_columnData), target :: columnTrlOnAnlIncLev
  type(struct_columnData), target :: columnTrlOnTrlLev
  type(struct_gsv)                :: stateVectorTrialHighRes
  type(struct_hco),       pointer :: hco_trl => null()
  type(struct_vco),       pointer :: vco_trl => null()

  logical           :: allocHeightSfc

  character(len=48) :: obsMpiStrategy, varMode

  type(struct_hco), pointer :: hco_anl => null()
  type(struct_hco), pointer :: hco_core => null()

  istamp = exdb('dfs', 'DEBUT', 'NON')

  call ver_printNameAndVersion('dfs', 'Dfs computation')

  ! MPI initilization
  call mmpi_initialize 

  call tmg_init(mmpi_myid, 'TMG_INFO')

  call utl_tmg_start(0, 'Main')
  call utl_printTime()

  varMode='analysis'

  ! Read the namelists
  call utl_readNml()

  call ram_setup

  ! Do initial set up

  call var_setup('VAR') ! obsColumnMode

  ! Reading trials
  call inn_getHcoVcoFromTrlmFile(hco_trl, vco_trl)
  allocHeightSfc = ( vco_trl%Vcode /= 0 )

  call gsv_allocate(stateVectorTrialHighRes, tim_nstepobs, hco_trl, vco_trl,  &
                    dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                    mpi_distribution_opt='Tiles', dataKind_opt=4,  &
                    allocHeightSfc_opt=allocHeightSfc, hInterpolateDegree_opt='LINEAR', &
                    beSilent_opt=.false.)
  call gsv_zero(stateVectorTrialHighRes)
  call gio_readTrials(stateVectorTrialHighRes)

  ! Horizontally interpolate trials to trial columns
  call inn_setupColumnsOnTrlLev(columnTrlOnTrlLev, obsSpaceData, hco_core, &
                                 stateVectorTrialHighRes )
  
  ! Interpolate trial columns to analysis levels and setup for linearized H
  call inn_setupColumnsOnAnlIncLev(columnTrlOnTrlLev, columnTrlOnAnlIncLev)

  ! Compute observation innovations and prepare obsSpaceData for minimization
  call inn_computeInnovation(columnTrlOnTrlLev, obsSpaceData)
  
  ! Compute perturbed
  call diagDFS(columnTrlOnAnlIncLev, obsSpaceData)

  ! Deallocate memory related to B matrices
  call bmat_finalize()

  ! Now write out the observation data files
  if ( .not. obsf_filesSplit() ) then 
    write(*,*) 'We read/write global observation files'
    call obs_expandToMpiGlobal(obsSpaceData)
    if (mmpi_myid == 0) call obsf_writeFiles(obsSpaceData)
  else
    ! redistribute obs data to how it was just after reading the files
    call obs_MpiRedistribute(obsSpaceData, OBS_IPF)
    call obsf_writeFiles(obsSpaceData)
  end if

  ! Deallocate copied obsSpaceData
  call obs_finalize(obsSpaceData)

  ! 3. Job termination

  istamp = exfin('dfs', 'FIN', 'NON')

  call utl_tmg_stop(0)
  call utl_printTime()

  call tmg_terminate(mmpi_myid, 'TMG_INFO')

  call rpn_comm_finalize(ierr) 

contains

  !--------------------------------------------------------------------------
  ! var_setup
  !--------------------------------------------------------------------------
  subroutine var_setup(obsColumnMode)
    !
    !:Purpose: Control of the preprocessing of the variational assimilation
    !
    implicit none

    ! Arguments:
    character (len=*), intent(in) :: obsColumnMode

    ! Locals:
    integer :: dateStampFromObs
    type(struct_vco), pointer :: vco_anl => null()
    integer, external :: get_max_rss

    write(*,*)
    write(*,*) '-----------------------------------'
    write(*,*) '-- Starting subroutine var_setup --'
    write(*,*) '-----------------------------------'

    !
    !- Initialize the Temporal grid and set dateStamp from env variable
    !
    call tim_setup()

    !     
    !- Initialize burp file names and set datestamp if not already
    !
    call obsf_setup(dateStampFromObs, 'analysis')
    if (tim_getDateStamp() == 0) then
      if (dateStampFromObs > 0) then
        call tim_setDatestamp(dateStampFromObs)
      else
        call utl_abort('var_setup: dateStamp was not set')
      end if
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
    write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

    !
    !- Initialize the Analysis grid
    !
    if (mmpi_myid == 0) write(*,*) ''
    if (mmpi_myid == 0) write(*,*) ' preproc: Set hco parameters for analysis grid'
    call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

    if (hco_anl % global) then
      hco_core => hco_anl
    else
      !- Initialized the core (Non-Exteded) analysis grid
      call hco_SetupFromFile(hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore') ! IN
    end if

    !     
    !- Initialisation of the analysis grid vertical coordinate from analysisgrid file
    !
    call vco_SetupFromFile(vco_anl, './analysisgrid') ! IN

    call col_setVco(columnTrlOnAnlIncLev, vco_anl)
    write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

    !
    !- Setup and read observations
    !
    obsMpiStrategy = 'LIKESPLITFILES'
    call inn_setupObs(obsSpaceData, hco_anl, obsColumnMode, obsMpiStrategy, varMode) ! IN
    write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

    !
    !- Basic setup of columnData module
    !
    call col_setup
    write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

    !
    !- Memory allocation for background column data
    !
    call col_allocate(columnTrlOnAnlIncLev, obs_numheader(obsSpaceData))

    !
    !- Initialize the observation error covariances
    !
    call oer_setObsErrors(obsSpaceData, varMode) ! IN
    write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

    !
    !- Initialize the background-error covariance, also sets up control vector module (cvm)
    !
    call bmat_setup(hco_anl, hco_core, vco_anl)
    write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

    !
    ! - Initialize the gridded variable transform module
    !
   
    call gvt_setup(hco_anl, hco_core, vco_anl)
    call gvt_setupRefFromTrialFiles('HU')
    call gvt_setupRefFromTrialFiles('height')
    
    write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  end subroutine var_setup

  !--------------------------------------------------------------------------
  ! diagDFS
  !--------------------------------------------------------------------------
  subroutine diagDFS(columnTrlOnAnlIncLev, obsSpaceData)
    !
    !:Purpose: Main DFS computation
    !
    implicit none

    ! Arguments:
    type(struct_obs),        intent(inout) :: obsSpaceData         ! Observation-related data
    type(struct_columnData), intent(inout) :: columnTrlOnAnlIncLev ! Background columns interpolated to anl levels and obs horiz locations

    ! Locals:
    type(struct_columnData) :: columnAnlInc
    type(struct_gsv)        :: stateVector
    type(struct_vco), pointer :: vco_anl
    integer :: headerIndex, bodyIndex1, bodyIndex2, obsIndex, taskIndex
    integer :: istart, iend, ierr
    integer :: procIndex
    integer, allocatable :: headerIndexList(:), channelList(:,:)
    integer, allocatable :: headerIndexListMpi(:,:), channelListMpi(:,:,:)
    integer, allocatable :: bodyIndexList(:,:), bodyIndexListMpi(:,:,:)
    integer, allocatable :: mpiTaskList(:), mpiTaskListMpi(:,:)
    real(8), allocatable :: stdDevList(:,:), stdDevListMpi(:,:,:)
    real(8), allocatable :: HBHtMatrix(:,:,:)
    real(8), pointer :: Rsub(:,:)
    real(8), pointer :: all_dfs(:)
    integer, pointer :: order(:)
    real(8) :: dfs
    integer :: localDimension
    logical :: first, llok
    integer :: channelNumber1, channelNumber2, channelIndex1, channelIndex2
    integer :: numHeader, numHeaderMaxMpi
    integer :: countObs, sumCountObsMpi, maxCountObsMpi, countChannel, maxCountChannelMpi
    integer :: sensorIndex
    integer, external :: get_max_rss
    real(8), allocatable :: perturbationVector(:)
    character(len=2) :: familyType
    integer :: nChannelsDfs = 103
    integer :: nObsMax = 2, nTaskMax = 2

    !
    !- 1.  Initialization

    write(*,*)
    write(*,*) 'Computing HBHT from selected observations start'

    
    vco_anl => col_getVco(columnTrlOnAnlIncLev)
    !- 1.3 Create a gridstateVector to store the perturbations
    call gsv_allocate(stateVector, tim_nstepobsinc, hco_anl, vco_anl, &
                      dataKind_opt=pre_incrReal, mpi_local_opt=.true.)

    !- 1.4 Create column vectors to store the perturbation interpolated to obs horizontal locations
    call col_setVco(columnAnlInc, vco_anl)
    call col_allocate(columnAnlInc, col_getNumCol(columnTrlOnAnlIncLev))

    !- 1.6
    call oti_timeBinning(obsSpaceData, tim_nstepobsinc)
    
    !numHeader = obs_numHeader(obsSpaceData)
    ! this is a temporary trick to be able to finish the run in reasonable time as it is SSSSLLLLOOOWWWWW
    ! some optimization or MPI topology change may help later
    numHeader = nObsMax
    call rpn_comm_allReduce(numHeader, numHeaderMaxMpi, 1, 'mpi_integer', 'mpi_max', 'grid', ierr)

    allocate(headerIndexList(numHeaderMaxMpi))
    allocate(mpiTaskList(numHeaderMaxMpi))
    allocate(channelList(numHeaderMaxMpi,nChannelsDfs))
    allocate(bodyIndexList(numHeaderMaxMpi,nChannelsDfs))
    allocate(stdDevList(numHeaderMaxMpi,nChannelsDfs))

    headerIndexList(:) = MPC_missingValue_INT
    mpiTaskList(:) = MPC_missingValue_INT
    channelList(:,:) = MPC_missingValue_INT
    bodyIndexList(:,:) = MPC_missingValue_INT
    stdDevList(:,:) = MPC_missingValue_R8
    
    !First step count the number of selected observation for each MPI task
    countObs = 0
    familyType = 'TO'
    first = .true.
    call obs_set_current_header_list(obsSpaceData,trim(familyType))
    HEADER1: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER1
      if (countObs == nObsMax) exit HEADER1
      istart = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
      iend = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) + istart - 1
      countChannel = 0
      BODY1:do bodyIndex1 = istart, iend
        llok = (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex1) == obs_assimilated)
        if (llok) then
          countChannel = countChannel + 1
        end if
      end do BODY1
      if (countChannel == nChannelsDfs) then
        write(*,*) "found one observation with ", nchannelsDfs, "levels/channels available ", headerIndex
        countObs = countObs + 1
        headerIndexList(countObs) = headerIndex
        mpiTaskList(countObs) = mmpi_myid
        countChannel = 0
        BODY2:do bodyIndex1 = istart, iend
          llok = (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex1) == obs_assimilated)
          if (llok) then
            call tvs_getChannelNumIndexFromPPP(obsSpaceData, headerIndex, bodyIndex1, &
                channelNumber1, channelIndex1)
            countChannel = countChannel + 1
            channelList(countObs, countChannel) = channelNumber1
            bodyIndexList(countObs, countChannel) = bodyIndex1
            stdDevList(countObs, countChannel) = obs_bodyElem_r(obsSpaceData, OBS_OER, bodyIndex1)
          end if
        end do BODY2
      end if
    end do HEADER1

    call rpn_comm_barrier('GRID', ierr)
    if (ierr /= 0) then
      call utl_abort('diagDFS: Error in call to rpn_comm_barrier 1')
    end if
    
    call rpn_comm_allReduce(countObs, sumCountObsMpi, 1, 'mpi_integer', 'mpi_sum', 'grid', ierr)
    if (ierr /= 0) then
      call utl_abort('diagDFS: Error in call to rpn_comm_allReduce 1')
    end if
    call rpn_comm_allReduce(countObs, maxCountObsMpi, 1, 'mpi_integer', 'mpi_max', 'grid', ierr)
    if (ierr /= 0) then
      call utl_abort('diagDFS: Error in call to rpn_comm_allReduce 2')
    end if
    call rpn_comm_allReduce(countChannel, maxCountChannelMpi, 1, 'mpi_integer', 'mpi_max', 'grid', ierr)
    if (ierr /= 0) then
      call utl_abort('diagDFS: Error in call to rpn_comm_allReduce 3')
    end if

    allocate(mpiTaskListMpi(maxCountObsMpi, mmpi_nprocs))
    allocate(headerIndexListMpi(maxCountObsMpi, mmpi_nprocs))
    allocate(channelListMpi(maxCountObsMpi, maxCountChannelMpi, mmpi_nprocs))
    allocate(bodyIndexListMpi(maxCountObsMpi, maxCountChannelMpi, mmpi_nprocs))
    allocate(stdDevListMpi(maxCountObsMpi, maxCountChannelMpi, mmpi_nprocs))
    
    mpiTaskListMpi(:,:) = MPC_missingValue_INT
    headerIndexListMpi(:,:) = MPC_missingValue_INT
    channelListMpi(:,:,:) = MPC_missingValue_INT
    bodyIndexListMpi(:,:,:) = MPC_missingValue_INT
    stdDevListMpi(:,:,:) = MPC_missingValue_R8
    
    call rpn_comm_allgather(mpiTaskList(1:maxCountObsMpi), maxCountObsMpi, 'mpi_integer', &
                            mpiTaskListMpi, maxCountObsMpi, 'mpi_integer', 'grid', ierr)
    if (ierr /= 0) then
      call utl_abort('diagDFS: Error in call to rpn_comm_allGather for mpiTaskList')
    end if
    call rpn_comm_allgather(headerIndexList(1:maxCountObsMpi), maxCountObsMpi, 'mpi_integer', &
                            headerIndexListMpi, maxCountObsMpi, 'mpi_integer', 'grid', ierr)
    if (ierr /= 0) then
      call utl_abort('diagDFS: Error in call to rpn_comm_allGather for headerIndexList')
    end if
    call rpn_comm_allgather(channelList(1:maxCountObsMpi,1:maxCountChannelMpi), &
        maxCountObsMpi * maxCountChannelMpi, 'mpi_integer',  &
        channelListMpi, maxCountObsMpi * maxCountChannelMpi, 'mpi_integer', 'grid', ierr)
    if (ierr /= 0) then
      call utl_abort('diagDFS: Error in call to rpn_comm_allGather for channelList')
    end if
    call rpn_comm_allgather(bodyIndexList(1:maxCountObsMpi,1:maxCountChannelMpi), &
        maxCountObsMpi * maxCountChannelMpi, 'mpi_integer',  &
        bodyIndexListMpi, maxCountObsMpi * maxCountChannelMpi, 'mpi_integer', 'grid', ierr)
    if (ierr /= 0) then
      call utl_abort('diagDFS: Error in call to rpn_comm_allGather for bodyIndexList')
    end if
    call rpn_comm_allgather(stdDevList(1:maxCountObsMpi,1:maxCountChannelMpi), &
        maxCountObsMpi * maxCountChannelMpi, 'mpi_real8',  &
        stdDevListMpi, maxCountObsMpi * maxCountChannelMpi, 'mpi_real8', 'grid', ierr)
    if (ierr /= 0) then
      call utl_abort('diagDFS: Error in call to rpn_comm_allGather for stdDevList')
    end if
    call rpn_comm_barrier('GRID', ierr)
    if (ierr /= 0) then
      call utl_abort('diagDFS: Error in call to rpn_comm_barrier 2')
    end if

    deallocate(bodyIndexList)
    deallocate(channelList)
    deallocate(mpiTaskList)
    deallocate(headerIndexList)
    deallocate(stdDevList)
    
    localDimension = cvm_nvadim
    allocate(perturbationVector(localDimension))
    allocate(HBHtMatrix(maxCountObsMpi,nChannelsDfs,nChannelsDfs))

    !do procIndex = 1, mmpi_nprocs
    do procIndex = 1, nTaskMax 
      do obsIndex = 1, maxCountObsMpi
        headerIndex = headerIndexListMpi(obsIndex,procIndex)
        taskIndex = mpiTaskListMpi(obsIndex,procIndex)
        write(*,*) "ICI HEADER", procIndex, obsIndex, headerIndex, taskIndex
        if (headerIndex /= MPC_missingValue_INT .and. taskIndex /= MPC_missingValue_INT) then
          do channelIndex1 = 1, maxCountChannelMpi
            bodyIndex1 = bodyIndexListMpi(obsIndex,channelIndex1,procIndex)
            channelNumber1 = channelListMpi(obsIndex,channelIndex1,procIndex)
            if (bodyIndex1 /= MPC_missingValue_INT) then
              write(*,*) "ICI body", procIndex, obsIndex, headerIndex, taskIndex, BODYiNDEX1, cHANNELnUMBER1
              ! We need to initialize the full OBS_WORK column to zero 
              do bodyIndex2 = 1, obs_numBody(obsSpaceData)
                call obs_bodySet_r(obsSpaceData, OBS_WORK, bodyIndex2, 0.d0)
              end do
              if (mmpi_myId == taskIndex) call obs_bodySet_r(obsSpaceData, OBS_WORK, bodyIndex1, 1.d0)
              write(*,*) "before oop_Had"
              write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
              call col_zero(columnAnlInc)
              call oop_Had(columnAnlInc, & !output
                  columnTrlOnAnlIncLev,  &
                  obsSpaceData,          & ! input
                  initializeLinearization_opt=first)
              first = .false.
              write(*,*) "before s2c_ad"
              write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
              call gsv_zero(stateVector)
              call s2c_ad(stateVector,  & ! output
                  columnAnlInc,         & ! input
                  columnTrlOnAnlIncLev, &
                  obsSpaceData)
              write(*,*) "before bmat_sqrtBT"
              write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
              perturbationVector(:) = 0.d0
              call bmat_sqrtBT(perturbationVector, & ! output
                  localDimension,                  &  
                  stateVector)                        ! input
              write(*,*) "before bmat_sqrtB"
              write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
              call gsv_zero(stateVector)
              call bmat_sqrtB(perturbationVector, & !input
                  localDimension,                 &
                  stateVector)                       ! output
              write(*,*) "before s2c_tl"
              write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
              call s2c_tl(stateVector,  & !input
                  columnAnlInc,         & ! output
                  columnTrlOnAnlIncLev, & 
                  obsSpaceData)
              write(*,*) "before oop_Htl"
              write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
              call oop_Htl(columnAnlInc, & ! input
                  columnTrlOnAnlIncLev,  &
                  obsSpaceData,          & !output
                  min_nsim=1, initializeLinearization_opt=.false.)
              
              do channelIndex2 = 1, maxCountChannelMpi
                bodyIndex2 = bodyIndexListMpi(obsIndex,channelIndex2,procIndex)
                channelNumber2 = channelListMpi(obsIndex,channelIndex2,procIndex)
                if (mmpi_myId == taskIndex .and. bodyIndex2 /= MPC_missingValue_INT) then
                  HBHtMatrix(obsIndex,channelIndex1,channelIndex2) = obs_bodyElem_r(obsSpaceData, OBS_WORK, bodyIndex2)
                  write(*,'(A4,1x,6(i12,1x),e14.6)') 'HBHt', obsIndex, headerIndex, channelNumber1, &
                      channelNumber2, bodyIndex1, bodyIndex2, obs_bodyElem_r(obsSpaceData, OBS_WORK, bodyIndex2)
                end if
              end do
            end if
          end do
          if (mmpi_myId == taskIndex) then
            sensorIndex = tvs_lsensor( tvs_tovsIndex(headerIndex) )
            call rmat_getRmatrix(sensorIndex,         &
                channelListMpi(obsIndex,:,procIndex), &
                stdDevListMpi(obsIndex,:,procIndex),  &
                Rsub)
            dfs = computeDfs(HBHtMatrix(obsIndex,:,:), Rsub)
            write(*,*) "dfs = ", dfs
            call selectChannels(HBHtMatrix(obsIndex,:,:), Rsub, all_dfs, order)
            write(*,*) 'all_dfs= ',all_dfs(:)
            write(*,*) 'order= ',order(:)
          end if
        end if
      end do
    end do

    deallocate(bodyIndexListMpi)
    deallocate(channelListMpi)
    deallocate(headerIndexListMpi)
    deallocate(mpiTaskListMpi)
    deallocate(stdDevListMpi)

    if (associated(Rsub)) deallocate(Rsub)
    if (allocated(HBHtMatrix)) deallocate(HBHtMatrix)
    
    deallocate(perturbationVector)
    call col_deallocate(columnAnlInc)
    call gsv_deallocate(stateVector)
    
    write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
    write(*,*)
    write(*,*) 'Computing DFS from selected observations end'

  end subroutine diagDFS

  !--------------------------------------------------------------------------
  ! subsetMatrix
  !--------------------------------------------------------------------------
  subroutine subsetMatrix(matrixInput, order, matrixOutput)
    !
    !:Purpose: extract sub-Matrix for a subset of channels/levels
    !
    implicit none
    
    ! Arguments:
    real(8), intent(in)  :: matrixInput(:,:)
    integer, intent(in)  :: order(:)
    real(8), intent(out) :: matrixOutput(size(order), size(order))

    matrixOutput(:, :) = matrixInput(order(:), order(:))
 
  end subroutine subsetMatrix
  
  !--------------------------------------------------------------------------
  ! computeDfs
  !--------------------------------------------------------------------------
  function computeDfs(HBHt, R) result(dfs)
    !
    !:Purpose: compute DFS from HBHt and R matrices
    !
    implicit none
    
    ! Arguments
    real(8), intent(in) :: HBHt(:,:)
    real(8), intent(in) :: R(:,:)
    real(8)             :: dfs
    
    ! Local variables
    integer :: nbChannels, channelIndex
    real(8), allocatable :: dMatrix(:,:), inverse(:,:), hk(:,:)

    nbChannels = size(R, dim=1)
    allocate(dMatrix(nbChannels,nbChannels))
    allocate(inverse(nbChannels,nbChannels))
    allocate(hk(nbChannels,nbChannels))
    dMatrix(:,:) =  HBHt(:,:) + R(:,:)
    call utl_pseudo_inverse(dMatrix, inverse)
    ! call utl_matInverse(matrix, rank, inverseSqrt_opt, printInformation_opt, &
    !                  eigenValueRelThreshold_opt)
    
    hk = matmul(HBHt, inverse)
    dfs = 0.d0
    do channelIndex = 1, nbChannels
       dfs = dfs + hk(channelIndex,channelIndex)
    end do
    deallocate(hk, inverse, dMatrix)
  end function computeDfs
  
  !--------------------------------------------------------------------------
  ! selectChannels
  !--------------------------------------------------------------------------
  subroutine selectChannels(HBHt, R, orderedDfs, orderedChannelIndexes, nChannelsOut_opt)
    !
    !:Purpose: perform DFS-based channel selection
    !
    
    implicit none
    ! Arguments
    real(8), intent(in) :: HBHt(:,:)
    real(8), intent(in) :: R(:,:)
    real(8), pointer, intent(inout):: orderedDfs(:)
    integer, pointer, intent(inout):: orderedChannelIndexes(:)
    integer, optional,  intent(in) :: nChannelsOut_opt
    ! Locals
    integer, allocatable :: tmpOrder(:), freeIndexList(:), tmpFree(:)
    real(8) :: optimalDfs, dfsTest
    integer :: channelIndex1, channelIndex2, optimalIndex, numberOfFreeChannels
    real(8), allocatable :: RSubset(:,:), HBHtSubset(:,:)
    integer :: nChannelsIn, nChannelsOut

    nChannelsIn = size(R, dim=1)
   
    if (present (nChannelsOut_opt)) then
      if (nChannelsOut_opt > nChannelsIn) then
        write(*,*) 'selectChannels: nChannelsIn, nChannelsOut_opt', nChannelsIn, nChannelsOut_opt
        call utl_abort('selectChannels: nChannelsOut_opt should be lower or equal than the dimension of input matrices.')
      end if
      nChannelsOut = nChannelsOut_opt
    else
      nChannelsOut = nChannelsIn
    end if
    
    allocate(freeIndexList(nChannelsIn))
    do channelIndex1 = 1, nChannelsIn
      freeIndexList(channelIndex1) = channelIndex1
    end do
    
    allocate(orderedChannelIndexes(nChannelsOut))
    orderedChannelIndexes(:) = -1
    allocate(tmpOrder(nChannelsOut))
    allocate(orderedDfs(nChannelsOut))
    do channelIndex1 = 1, nChannelsOut
      optimalIndex = 0
      optimalDfs = 0.0
      allocate(RSubset(channelIndex1,channelIndex1), HBHtSubset(channelIndex1,channelIndex1))
      numberOfFreeChannels = nChannelsIn - channelIndex1 + 1
      do channelIndex2 = 1, numberOfFreeChannels
        tmpOrder(1:channelIndex1-1) = orderedChannelIndexes(1:channelIndex1-1)
        tmpOrder(channelIndex1) = freeIndexList(channelIndex2)
        call subsetMatrix(R, tmpOrder(1:channelIndex1), RSubset)
        call subsetMatrix(HBHt, tmpOrder(1:channelIndex1), HBHtSubset)
        dfsTest = computeDfs(HBHtSubset, RSubset)
        if (dfsTest > optimalDfs) then
          optimalIndex =freeIndexList(channelIndex2)
          optimalDfs = dfsTest
        end if
      end do
      deallocate(RSubset, HBHtSubset)
      orderedDfs(channelIndex1) = optimalDfs
      orderedChannelIndexes(channelIndex1) = optimalIndex
      allocate(tmpFree(numberOfFreeChannels))
      !https://stackoverflow.com/questions/42140832/automatic-array-allocation-upon-assignment-in-fortran could help to simplify a bit this 
      !but it is needed to get rif of the (:)
      tmpFree(:) = freeIndexList(:)
      call utl_reAllocate(freeIndexList, numberOfFreeChannels - 1)
      freeIndexList(:) = pack(tmpFree, tmpFree /= optimalIndex)
      deallocate(tmpFree)
    end do
    
    deallocate(freeIndexList)
    deallocate(tmpOrder)

  end subroutine selectChannels
    
end program midas_dfs

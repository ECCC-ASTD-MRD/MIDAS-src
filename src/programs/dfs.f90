program midas_dfs
  !
  !:Purpose: Main program for computing DFS
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

  integer, external :: exdb, exfin,fnom, fclos
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

  integer :: nLevelsDfs, levelIndex
  !Namelist variables:
  character(len=2) :: familyType                 ! familyType to consider (TO, UA, AI, RO, etc... one at a time)
  logical :: doChannelSelection                  ! flag to perform DFS-based channel selection (TO only)
  integer :: maxSelect                           ! max number of channels to select (negative or zero to do all channels)
  logical :: outputHBHt                          ! flag to output HBHt
  integer :: nDfsMax                             ! maximum number of DFS computations
  integer :: vCoordList(tvs_maxNumberOfChannels) ! list of channels or levels (depending on FamilyType)
                                                 ! Dfs will be computed only for observation locations for which these levels are available
  NAMELIST /NAMDFS/ familyType, doChannelSelection, maxSelect, outputHBHt, nDfsMax, vCoordList
  
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
  
  ! Set/Read values for the namelist NAMDFS
  ! setting default values
  familyType = 'TO'
  doChannelSelection = .true.
  maxSelect = MPC_missingValue_INT
  outputHBHt = .true.
  nDfsMax = 4
  vCoordList(:) = MPC_missingValue_INT
  ! Check if NAMDFS exist
  if (.not. utl_isNamelistPresent('NAMDFS','./flnml')) then
    write(*,*)
    write(*,*) 'dfs: Namelist block NAMDFS is missing in the namelist.'
    write(*,*) '       the default values will be taken.'
  else
    ! Read the namelist
    call utl_tmg_start(181,'low-level--readNML')
    read(utl_flnml, nml=NAMDFS, iostat=ierr)
    call utl_tmg_stop(181)
    if (ierr /= 0) call utl_abort('midas-dfs: Error reading namelist')
  end if
    
  if (mmpi_myid == 0) write(*, nml=NAMDFS)

  nLevelsDfs = 0
  do levelIndex = 1, tvs_maxNumberOfChannels
    if (vCoordList(levelIndex) == MPC_missingValue_INT) exit
    nLevelsDfs = nLevelsDfs + 1
  end do

  if (nLevelsDfs == 0) call utl_abort('midas-dfs: empty vertical coordinate list vCoordList !') 
 
  if (doChannelSelection .and. familyType /= 'TO') then
    call utl_abort('midas-dfs: DFS-based channel selection does not make sense for families other than TO !')
  end if
  
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
  
  ! Compute HBHt, dfs and perform channel selection
  call diagDFS(columnTrlOnAnlIncLev, obsSpaceData)

  ! Deallocate memory related to B matrices
  call bmat_finalize()

  ! The following lines of code are no longer necessary unless it is decided to put some of the computed outputs in BURP or SQLITE observation files
  ! to be discussed with Mark
  
  ! Now write out the observation data files
  !if ( .not. obsf_filesSplit() ) then 
  !  write(*,*) 'We read/write global observation files'
  !  call obs_expandToMpiGlobal(obsSpaceData)
  !  if (mmpi_myid == 0) call obsf_writeFiles(obsSpaceData)
  !else
  !  ! redistribute obs data to how it was just after reading the files
  !  call obs_MpiRedistribute(obsSpaceData, OBS_IPF)
  !  call obsf_writeFiles(obsSpaceData)
  !end if

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
    call utl_printTime()
    
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
    write(*,*) 'var_setup: exiting...'
    call utl_printTime()
    
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
    integer, allocatable :: headerIndexList(:), levelList(:,:)
    integer, allocatable :: headerIndexListMpi(:,:), levelListMpi(:,:,:)
    integer, allocatable :: bodyIndexList(:,:), bodyIndexListMpi(:,:,:)
    integer, allocatable :: mpiTaskList(:), mpiTaskListMpi(:,:)
    real(8), allocatable :: stdDevList(:,:), stdDevListMpi(:,:,:)
    real(8), allocatable :: HBHtMatrix(:,:)
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
    integer :: dfsCount
    integer :: nulDfs, nulSelec, nulHBHt
    character(len=128) :: headerObs1
    character(len=16) :: headerObs2
    !
    !- 1.  Initialization

    write(*,*)
    write(*,*) 'Computing HBHT from selected observations start'
    call utl_printTime()
    
    vco_anl => col_getVco(columnTrlOnAnlIncLev)
    !- 1.3 Create a gridstateVector to store the perturbations
    call gsv_allocate(stateVector, tim_nstepobsinc, hco_anl, vco_anl, &
                      dataKind_opt=pre_incrReal, mpi_local_opt=.true.)

    !- 1.4 Create column vectors to store the perturbation interpolated to obs horizontal locations
    call col_setVco(columnAnlInc, vco_anl)
    call col_allocate(columnAnlInc, col_getNumCol(columnTrlOnAnlIncLev))

    !- 1.6
    call oti_timeBinning(obsSpaceData, tim_nstepobsinc)
    
    numHeader = obs_numHeader(obsSpaceData)
    call rpn_comm_allReduce(numHeader, numHeaderMaxMpi, 1, 'mpi_integer', 'mpi_max', 'grid', ierr)

    allocate(headerIndexList(numHeaderMaxMpi))
    allocate(mpiTaskList(numHeaderMaxMpi))
    allocate(levelList(numHeaderMaxMpi,nLevelsDfs))
    allocate(bodyIndexList(numHeaderMaxMpi,nLevelsDfs))
    allocate(stdDevList(numHeaderMaxMpi,nLevelsDfs))

    headerIndexList(:) = MPC_missingValue_INT
    mpiTaskList(:) = MPC_missingValue_INT
    levelList(:,:) = MPC_missingValue_INT
    bodyIndexList(:,:) = MPC_missingValue_INT
    stdDevList(:,:) = MPC_missingValue_R8
    
    !First step count the number of selected observation for each MPI task
    countObs = 0
    first = .true.
    countChannel = 0 ! necessary in the case where no obs in the file
    call obs_set_current_header_list(obsSpaceData,trim(familyType))
    HEADER1: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER1
      istart = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
      iend = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) + istart - 1
      countChannel = 0
      BODY1:do bodyIndex1 = istart, iend
        llok = (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex1) == obs_assimilated)
        if (llok) then
          call tvs_getChannelNumIndexFromPPP(obsSpaceData, headerIndex, bodyIndex1, &
              channelNumber1, channelIndex1)
          if (utl_findloc(vCoordList(1:nLevelsDfs), channelNumber1) > 0) then
            countChannel = countChannel + 1
          end if
        end if
      end do BODY1
      if (countChannel == nLevelsDfs) then
        write(*,*) "found one observation with all requested ", nlevelsDfs, "levels/channels available ", headerIndex
        countObs = countObs + 1
        headerIndexList(countObs) = headerIndex
        mpiTaskList(countObs) = mmpi_myid
        levelList(countObs,:) = vCoordList(1:nLevelsDfs)
        countChannel = 0
        BODY2:do bodyIndex1 = istart, iend
          llok = (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex1) == obs_assimilated)
          if (llok) then
            countChannel = countChannel + 1
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
    allocate(levelListMpi(maxCountObsMpi, maxCountChannelMpi, mmpi_nprocs))
    allocate(bodyIndexListMpi(maxCountObsMpi, maxCountChannelMpi, mmpi_nprocs))
    allocate(stdDevListMpi(maxCountObsMpi, maxCountChannelMpi, mmpi_nprocs))
    
    mpiTaskListMpi(:,:) = MPC_missingValue_INT
    headerIndexListMpi(:,:) = MPC_missingValue_INT
    levelListMpi(:,:,:) = MPC_missingValue_INT
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
    call rpn_comm_allgather(levelList(1:maxCountObsMpi,1:maxCountChannelMpi), &
        maxCountObsMpi * maxCountChannelMpi, 'mpi_integer',  &
        levelListMpi, maxCountObsMpi * maxCountChannelMpi, 'mpi_integer', 'grid', ierr)
    if (ierr /= 0) then
      call utl_abort('diagDFS: Error in call to rpn_comm_allGather for levelList')
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
    deallocate(levelList)
    deallocate(mpiTaskList)
    deallocate(headerIndexList)
    deallocate(stdDevList)
    
    if (mmpi_myId == 0) then
      if (outputHbHt) then
        nulhbht = 0
        ierr = fnom(nulhbht, './HBHt.dat', 'FTN+FMT+R/W', 0)
      end if
      if (doChannelSelection) then
        nulSelec = 0
        ierr = fnom(nulSelec, './selection.dat', 'FTN+FMT+R/W', 0)
      end if
      nulDfs = 0
      ierr = fnom(nulDfs, './dfs.dat', 'FTN+FMT+R/W', 0)
    end if
    
    localDimension = cvm_nvadim
    allocate(perturbationVector(localDimension))
    allocate(HBHtMatrix(nLevelsDfs,nLevelsDfs))
    dfsCount = 0
    mpiTaskLoop:do procIndex = 1, mmpi_nprocs
      do obsIndex = 1, maxCountObsMpi
        headerIndex = headerIndexListMpi(obsIndex,procIndex)
        taskIndex = mpiTaskListMpi(obsIndex,procIndex)
        if (headerIndex /= MPC_missingValue_INT .and. taskIndex /= MPC_missingValue_INT) then
          do channelIndex1 = 1, maxCountChannelMpi
            bodyIndex1 = bodyIndexListMpi(obsIndex,channelIndex1,procIndex)
            channelNumber1 = levelListMpi(obsIndex,channelIndex1,procIndex)
            if (bodyIndex1 /= MPC_missingValue_INT) then
              ! We need to initialize the full OBS_WORK column to zero 
              do bodyIndex2 = 1, obs_numBody(obsSpaceData)
                call obs_bodySet_r(obsSpaceData, OBS_WORK, bodyIndex2, 0.d0)
              end do
              if (mmpi_myId == taskIndex) call obs_bodySet_r(obsSpaceData, OBS_WORK, bodyIndex1, 1.d0)
              write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
              call col_zero(columnAnlInc)
              call oop_Had(columnAnlInc, & !output
                  columnTrlOnAnlIncLev,  &
                  obsSpaceData,          & !input
                  initializeLinearization_opt=first)
              first = .false.
              write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
              call gsv_zero(stateVector)
              call s2c_ad(stateVector,  & ! output
                  columnAnlInc,         & ! input
                  columnTrlOnAnlIncLev, &
                  obsSpaceData)
              write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
              perturbationVector(:) = 0.d0
              call bmat_sqrtBT(perturbationVector, & !output
                  localDimension,                  &  
                  stateVector)                       !input
              write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
              call gsv_zero(stateVector)
              call bmat_sqrtB(perturbationVector, & !input
                  localDimension,                 &
                  stateVector)                      !output
              write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
              call s2c_tl(stateVector,  & !input
                  columnAnlInc,         & !output
                  columnTrlOnAnlIncLev, & 
                  obsSpaceData)
              write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
              call oop_Htl(columnAnlInc, & !input
                  columnTrlOnAnlIncLev,  &
                  obsSpaceData,          & !output
                  min_nsim=1, initializeLinearization_opt=.false.)
              
              do channelIndex2 = 1, maxCountChannelMpi
                bodyIndex2 = bodyIndexListMpi(obsIndex,channelIndex2,procIndex)
                channelNumber2 = levelListMpi(obsIndex,channelIndex2,procIndex)
                if (mmpi_myId == taskIndex .and. bodyIndex2 /= MPC_missingValue_INT) then
                  HBHtMatrix(channelIndex1,channelIndex2) = obs_bodyElem_r(obsSpaceData, OBS_WORK, bodyIndex2)
                end if
              end do
            end if
            write(*,*) 'diagDfs: computed column ', channelIndex1, &
                'of HBHt'
            call utl_printTime()
          end do
          dfsCount = dfsCount + 1
          if (mmpi_myId == taskIndex) then
            write(headerObs1,"('# ',A12,1x,2e14.6,1x,i8.8,1x,i4.4)") &
                obs_elem_c(obsSpaceData, 'STID' ,headerIndex), &
                obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex), &
                obs_headElem_r(obsSpaceData, OBS_LON, headerIndex), &
                obs_headElem_i(obsSpaceData, OBS_DAT, headerIndex), &
                obs_headElem_i(obsSpaceData, OBS_ETM, headerIndex)
            headerObs2 = ''
            if (familyType == 'TO')  then
              write(headerObs2,'(1x,e14.6)') obs_headElem_r(obsSpaceData, OBS_SZA, headerIndex)
            end if
            
            sensorIndex = tvs_lsensor( tvs_tovsIndex(headerIndex) )
            
            if (familyType == 'TO') then
              call rmat_getRmatrix(sensorIndex,        &
                  levelListMpi(obsIndex,:,procIndex),  &
                  stdDevListMpi(obsIndex,:,procIndex), &
                  Rsub)
            else
              allocate(Rsub(nLevelsDfs,nLevelsDfs))
              Rsub(:,:) = 0.d0
              do channelIndex2 = 1, nLevelsDfs
                Rsub(channelIndex2,channelIndex2) = stdDevListMpi(obsIndex,channelIndex2,procIndex) ** 2
              end do
            end if
            dfs = computeDfs(HBHtMatrix, Rsub)
            if (doChannelSelection) then
              call selectChannels(HBHtMatrix, Rsub, all_dfs, order, maxSelect)
            end if
            deallocate(Rsub)
          else
            if (doChannelSelection) then
              allocate(all_dfs(nLevelsDfs))
              allocate(order(nLevelsDfs))
            end if
          end if
          
          call rpn_comm_barrier('GRID', ierr)
          if (ierr /= 0) then
            call utl_abort('diagDFS: Error in call to rpn_comm_barrier 3')
          end if

          call rpn_comm_bcastc(headerObs1, len(headerObs1), 'MPI_CHARACTER', taskIndex, 'GRID', ierr)
          call rpn_comm_bcastc(headerObs2, len(headerObs2), 'MPI_CHARACTER', taskIndex, 'GRID', ierr)
          headerObs1 = trim(headerObs1) // trim(headerObs2)
          if (outputHBHt) then
            call rpn_comm_bcast(HBHtMatrix, nLevelsDfs*nLevelsDfs, 'MPI_REAL8', taskIndex, 'GRID', ierr)
            if (mmpi_myId == 0) then
              write(nulHBHt,'(A)') trim(headerObs1)
              do channelIndex2 = 1, nLevelsDfs
                channelNumber2 = levelListMpi(obsIndex,channelIndex2,procIndex)
                do channelIndex1 = 1, nLevelsDfs
                  channelNumber1 = levelListMpi(obsIndex,channelIndex1,procIndex)
                  write(nulHBHt,'(A4,1x,2(i12,1x),e14.6)') 'HBHt', channelNumber1, &
                      channelNumber2, HBHtMatrix(channelIndex1,channelIndex2)
                end do
              end do
              write(nulHBHt,*)
            end if
          end if
          call rpn_comm_bcast(dfs, 1, 'MPI_REAL8', taskIndex, 'GRID', ierr)
          if (mmpi_myId == 0) write(nulDfs,'(A,1x,e14.6)') trim(headerObs1), dfs
          if (doChannelSelection) then
            call rpn_comm_bcast(all_dfs, nLevelsDfs, 'MPI_REAL8', taskIndex, 'GRID', ierr)
            call rpn_comm_bcast(order, nLevelsDfs, 'MPI_INTEGER', taskIndex, 'GRID', ierr)
            if (mmpi_myId == 0) then
              write(nulSelec,'(A)') trim(headerObs1)
              do channelIndex1 = 1, nLevelsDfs
                write(nulSelec,'(3(i5,1x),e14.6)')  channelIndex1, order(channelIndex1), &
                    levelListMpi(obsIndex,order(channelIndex1),procIndex), all_dfs(channelIndex1)
              end do
              write(nulSelec,*)
            end if
            deallocate(order)
            deallocate(all_dfs)
          end if
          if (dfsCount == nDfsMax) exit mpiTaskLoop
        end if
      end do
    end do mpiTaskLoop

    if (mmpi_myId == 0) then
      if (outputHbHt) then
        ierr = fclos(nulhbht)
      end if
      if (doChannelSelection) then
        ierr = fclos(nulSelec)
      end if
      ierr = fclos(nuldfs)
    end if
    
    deallocate(bodyIndexListMpi)
    deallocate(levelListMpi)
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
    !:Purpose: compute DFS from given HBHt and R matrices of the appropriate size
    !
    implicit none
    
    ! Arguments
    real(8), intent(in) :: HBHt(:,:)
    real(8), intent(in) :: R(:,:)
    real(8)             :: dfs
    
    ! Local variables
    integer :: nbLevels, levelIndex
    real(8), allocatable :: dMatrix(:,:), inverse(:,:), hk(:,:)

    nbLevels = size(R, dim=1)
    allocate(dMatrix(nbLevels,nbLevels))
    allocate(inverse(nbLevels,nbLevels))
    allocate(hk(nbLevels,nbLevels))
    
    dMatrix(:,:) =  HBHt(:,:) + R(:,:)
    call utl_pseudo_inverse(dMatrix, inverse)
    hk = matmul(HBHt, inverse)
    dfs = 0.d0
    do levelIndex = 1, nbLevels
       dfs = dfs + hk(levelIndex,levelIndex)
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
   
    if (present(nChannelsOut_opt)) then
      if (nChannelsOut_opt > nChannelsIn) then
        write(*,*) 'selectChannels: nChannelsIn, nChannelsOut_opt', nChannelsIn, nChannelsOut_opt
        call utl_abort('selectChannels: nChannelsOut_opt should be lower or equal than the dimension of input R and HBHt matrices.')
      end if
      if (nChannelsOut_opt < 1) then
        nChannelsOut = nChannelsIn
      else
        nChannelsOut = nChannelsOut_opt
      end if
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
      !but it is needed to get rid of the (:)
      tmpFree(:) = freeIndexList(:)
      call utl_reAllocate(freeIndexList, numberOfFreeChannels - 1)
      freeIndexList(:) = pack(tmpFree, tmpFree /= optimalIndex)
      deallocate(tmpFree)
    end do
    
    deallocate(freeIndexList)
    deallocate(tmpOrder)

  end subroutine selectChannels
    
end program midas_dfs

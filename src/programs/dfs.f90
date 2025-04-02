program midas_dfs
  !
  !:Purpose: Main program for computing degrees of freedom for signal (DFS)
  !          total and a iterative channel selection (optionally) for an infrared instrument.
  !            
  !
  !          --
  !
  !:Algorithm: The channel selection with the DFS is computed following these steps:            
  !
  !            1. HBHt is extracted column by column for the specified instrument for all the channels
  !
  !            2. The total DFS is calculated for all channels
  !                           
  !            3.1 Channel selection: DFS is calculated (DFS= tr(HBHt (HBHt+R)-1)) for each channel of
  !               the instrument individually and the channel with the largest DFS is selected
  !
  !            3.2 DFS is calculated to find a second channel after the first selected 
  !               channel with the largest total DFS combined
  !
  !            3.4 Channels are selected iteratively until all the channels are ordered
  !               or until a certain specified number of channels is reached
  !
  !            --
  !
  !:File I/O: The required input files and produced output files can vary
  !           according to the application. 
  !
  !           --
  !
  !============================================== ==============================================================
  ! Input and Output Files (for NWP application)   Description of file
  !============================================== ==============================================================
  ! ``flnml``                                      In - Main namelist file with parameters user may modify
  ! ``flnml_static``                               In - The "static" namelist that should not be modified
  ! ``trlm_$NN`` (e.g. ``trlm_01``)                In - Background state (a.k.a. trial) files for each timestep
  ! ``analysisgrid``                               In - File defining grid for computing HBHt
  ! ``bgcov``                                      In - Static (i.e. NMC) B matrix file for NWP fields
  ! ``ensemble/$YYYYMMDDHH_006_$NNNN``             In - Ensemble member files defining ensemble B matrix
  ! ``obsfiles_$FAM/obs$FAM_$NNNN_$NNNN``          In - Observation file for each "family" and MPI task
  ! ``obserr``                                     In - Conventionnal Observation error statistics
  ! Remainder are files related to radiance obs:
  ! ``stats_tovs``                                 In - Satellite radiance observation errors
  ! ``rtcoef_$PLATFORM_$SENSOR.dat``               In - RTTOV coefficient files
  ! ``ozoneclim98``                                In - Ozone climatology
  ! ``Cmat_$PLATFORM_$SENSOR.dat``                 In - Correlations structure of the R matrix
  ! ``dfs.dat``                                    Out - Total DFS of all channels available
  ! ``HBHt.dat``                                   Out - HBHt matrix for the specified instrument
  ! ``selection.dat``                              Out - Selected channels ordered by decreasing contribution to the DFS 
  !============================================== ==============================================================
  !
  !           --
  !
  !:Synopsis: Below is a summary of the ``dfs`` program calling sequence:
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
  !               - Create a MPI local list of observations corresponding to criteria
  !                 from the namelist
  !
  !               - Gather all selected observations
  !
  !               - Loop on the selected observations and a loop on channels on the
  !                 following steps
  !
  !               - Set the corresponding entry of ``obs_work`` column of
  !                   ``obsSpaceData`` to one
  !
  !               - Apply the adjoint of the observation operator ``oop_Had``
  !
  !               - Apply the adjoint of the horizontal interpolation ``s2c_ad``
  !
  !               - Apply the transpose of the sqrt of B matrix ``bmat_sqrtBT``
  !
  !               - Apply the sqrt of B matrix ``bmat_sqrtB``
  !
  !               - Apply the tangent linear of the horizontal interpolation
  !                  ``s2c_tl``
  !
  !               - Apply the tangent linear of the observation operator ``oop_Htl``
  !
  !               - At the end of the loop on the channels, we get a column of HBHt matrix
  !
  !               - Get the R matrix
  !
  !               - Compute the total DFS
  !
  !               - If requested, compute the channel selection 
  !
  !
  !
  !           --
  !
  !:Options: `List of namelist blocks <../namelists_in_each_program.html#dfs>`_
  !          that can affect the ``dfs`` program.
  !
  !          * The NAMDFS section is read by the DFS program and controls which
  !            observations and channels are considered.
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
  !            dfs calculation are listed in the following table:
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
  use message_mod
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
  use tovs_mod
  use rMatrix_mod
  use thinning_mod
  
  implicit none

  integer, external :: exdb, exfin, fnom, fclos
  integer :: ierr, istamp, obsIndex

  type(struct_obs),        target :: obsSpaceData
  type(struct_columnData), target :: columnTrlOnAnlIncLev
  type(struct_columnData), target :: columnTrlOnTrlLev
  type(struct_gsv)                :: stateVectorTrialHighRes
  type(struct_hco),       pointer :: hco_trl => null()
  type(struct_vco),       pointer :: vco_trl => null()

  logical           :: allocHeightSfc
  logical           :: selectSpecificObservationsFromList
  character(len=48) :: obsMpiStrategy, varMode

  type(struct_hco), pointer :: hco_anl => null()
  type(struct_hco), pointer :: hco_core => null()
  integer, parameter :: nObsMax=10

  integer :: nLevelsDfs, levelIndex
  
  ! Namelist variables:
  character(len=2) :: familyType                 ! familyType to consider (TO, UA, AI, RO, etc... one at a time)
  logical :: doChannelSelection                  ! flag to perform DFS-based channel selection (TO only)
  integer :: maxSelect                           ! max number of channels to select (negative or zero to do all channels)
  logical :: outputHBHt                          ! flag to output HBHt
  logical :: doThinning                          ! flag to perform thinning on the observations, if .true. thinning is done 
  integer :: nDfsMax                             ! maximum number of DFS computations
  integer :: vCoordList(tvs_maxNumberOfChannels) ! list of channels or levels (depending on FamilyType)
                                                 ! Dfs will be computed only for observation locations for which these levels are available
  real(8) :: latList(nObsMax)                    ! list of latitudes to select specific observations
  real(8) :: lonList(nObsMax)                    ! list of longitudes to select specific observations
  integer :: dayList(nObsMax)                    ! list of dates (yyyymmdd) to select specific observations
  integer :: timeList(nObsMax)                   ! list of hours (HHMM) to select specific observations
  real(8) :: satZenList(nObsMax)                 ! list of satellite zenith angles to select specific observations
  logical :: computeInParallel                   ! if .true. computation performed in parallel
                                                 ! if .false. (default) observation processed one at a time (slow)
  real(8) :: highResWaterFractionThreshold       ! value for the highResWaterFraction threshold
  real(8) :: lowResWaterFractionThreshold        ! value for the lowResWaterFraction threshold
  
  NAMELIST /NAMDFS/ familyType, doChannelSelection, maxSelect, outputHBHt, nDfsMax, vCoordList, latList, lonList, dayList, timeList, satZenList, computeInParallel, doThinning, highResWaterFractionThreshold, lowResWaterFractionThreshold
  
  istamp = exdb('dfs', 'DEBUT', 'NON')

  call ver_printNameAndVersion('dfs', 'Dfs computation')

  ! MPI initilization
  call mmpi_initialize 

  call tmg_init(mmpi_myid, 'TMG_INFO')

  call utl_tmg_start(0, 'Main')
  call utl_printTime()

  varMode = 'analysis'

  ! Read the namelists
  call utl_readNml()

  call ram_setup

  ! Setup the format of the output RPN standard files to 'XDF' or 'RSF'
  call gio_setup

  ! Do initial set up
  
  ! Set/Read values for the namelist NAMDFS
  ! setting default values
  familyType = 'TO'
  doChannelSelection = .true.
  maxSelect = MPC_missingValue_INT
  outputHBHt = .true.
  nDfsMax = 4
  vCoordList(:) = MPC_missingValue_INT
  latList(:) = MPC_missingValue_R8
  lonList(:) = MPC_missingValue_R8
  dayList(:) = MPC_missingValue_INT
  timeList(:) = MPC_missingValue_INT
  satZenList(:) = MPC_missingValue_R8
  computeInParallel = .false.
  doThinning = .false.
  highResWaterFractionThreshold=0.99d0
  lowResWaterFractionThreshold=0.97d0

  
  ! Check if NAMDFS exist
  if (.not. utl_isNamelistPresent('NAMDFS','./flnml')) then
    write(*,*)
    write(*,*) 'dfs: Namelist block NAMDFS is missing in the namelist.'
    write(*,*) '       the default values will be taken.'
  else
    ! Read the namelist
    call utl_tmg_start(181,'low-level--readNML')
    read(utl_flnml, nml = NAMDFS, iostat = ierr)
    call utl_tmg_stop(181)
    if (ierr /= 0) call utl_abort('midas-dfs: Error reading namelist')
  end if

  ! Longitudes in ObsSpaceData are positives
  do obsIndex = 1, nObsMax
    if (lonList(obsindex) == MPC_missingValue_R8) cycle
    if (lonList(obsIndex) < 0. ) lonList(obsIndex) = lonList(obsIndex) + 360.d0
  end do
  
  if (mmpi_myid == 0) write(*, nml = NAMDFS)

  nLevelsDfs = 0
  do levelIndex = 1, tvs_maxNumberOfChannels
    if (vCoordList(levelIndex) == MPC_missingValue_INT) exit
    nLevelsDfs = nLevelsDfs + 1
  end do

  if (nLevelsDfs == 0) call utl_abort('midas-dfs: empty vertical coordinate list vCoordList !') 
 
  if (doChannelSelection .and. familyType /= 'TO') then
    call utl_abort('midas-dfs: DFS-based channel selection does not make sense for families other than TO !')
  end if

  selectSpecificObservationsFromList = .not. all( latList(:) == MPC_missingValue_R8   .and. &
                                                  lonList(:) == MPC_missingValue_R8   .and. &
                                                  dayList(:) == MPC_missingValue_INT  .and. &
                                                  timeList(:) == MPC_missingValue_INT .and. &
                                                  satZenList(:) == MPC_missingValue_R8 )
  call dfs_setup('VAR') ! obsColumnMode

  ! Reading trials
  call inn_getHcoVcoFromTrlmFile(hco_trl, vco_trl, './trlm_01')
  allocHeightSfc = ( vco_trl%Vcode /= 0 )

  call gsv_allocate(stateVectorTrialHighRes, tim_nstepobs, hco_trl, vco_trl,  &
                    dateStamp_opt = tim_getDateStamp(), mpi_local_opt = .true., &
                    mpi_distribution_opt = 'Tiles', dataKind_opt=4,  &
                    allocHeightSfc_opt = allocHeightSfc, hInterpolateDegree_opt = 'LINEAR', &
                    beSilent_opt = .false.)
  call gsv_zero(stateVectorTrialHighRes)
  call gio_readTrials(stateVectorTrialHighRes)

  ! Horizontally interpolate trials to trial columns
  call inn_setupColumnsOnTrlLev(columnTrlOnTrlLev, obsSpaceData, hco_core, &
                                 stateVectorTrialHighRes )
  
  ! Interpolate trial columns to analysis levels and setup for linearized H
  call inn_setupColumnsOnAnlIncLev(columnTrlOnTrlLev, columnTrlOnAnlIncLev)

  ! Compute observation innovations and prepare obsSpaceData for minimization
  call inn_computeInnovation(columnTrlOnTrlLev, obsSpaceData)

  ! Set up parameters for water fraction
  call tvs_emis_read_climatology ()
  call tvs_allocateSurfaceParameters ()
  
  ! Compute HBHt, dfs and perform channel selection
  call diagDFS(columnTrlOnAnlIncLev, obsSpaceData)

  ! Deallocate memory related to B matrices
  call bmat_finalize()

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
  ! dfs_setup
  !--------------------------------------------------------------------------
  subroutine dfs_setup(obsColumnMode)
    !
    !:Purpose: Control of the preprocessing of the dfs computation
    !
    implicit none

    ! Arguments:
    character (len=*), intent(in) :: obsColumnMode ! obsSpaceData initialisation mode

    ! Locals:
    integer :: dateStampFromObs
    type(struct_vco), pointer :: vco_anl => null()

    write(*,*)
    write(*,*) '-----------------------------------'
    write(*,*) '-- Starting subroutine dfs_setup --'
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
        call utl_abort('dfs_setup: dateStamp was not set')
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
    call msg_memUsage('midas-dfs')
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
    call msg_memUsage('midas-dfs')

    !
    !- Setup and read observations
    !
    obsMpiStrategy = 'LIKESPLITFILES'
    call inn_setupObs(obsSpaceData, hco_anl, obsColumnMode, obsMpiStrategy, varMode) ! IN
    call msg_memUsage('midas-dfs')

    !
    !- Basic setup of columnData module
    !
    call col_setup
    call msg_memUsage('midas-dfs')

    !
    !- Memory allocation for background column data
    !
    call col_allocate(columnTrlOnAnlIncLev, obs_numheader(obsSpaceData))

    !
    !- Initialize the observation error covariances
    !
    call oer_setObsErrors(obsSpaceData, varMode) ! IN
    call msg_memUsage('midas-dfs')

    !
    !- Initialize the background-error covariance, also sets up control vector module (cvm)
    !
    call bmat_setup(hco_anl, hco_core, vco_anl)
    call msg_memUsage('midas-dfs')

    !
    ! - Initialize the gridded variable transform module
    !   
    call gvt_setup(hco_anl, hco_core, vco_anl)
    call gvt_setupRefFromTrialFiles('HU')
    call gvt_setupRefFromTrialFiles('height')
    
    call msg_memUsage('midas-dfs')
    write(*,*) 'dfs_setup: exiting...'
    call utl_printTime()
    
  end subroutine dfs_setup

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
    integer :: bodyIndexBeg, bodyIndexEnd, ierr, procIndex, outTaskIndex
    integer, allocatable :: headerIndexList(:), levelList(:,:)
    integer, allocatable :: headerIndexListMpi(:,:), levelListMpi(:,:,:)
    integer, allocatable :: bodyIndexList(:,:), bodyIndexListMpi(:,:,:)
    real(8), allocatable :: stdDevList(:,:), stdDevListMpi(:,:,:)
    real(8), allocatable :: HBHtMatrix(:,:), Rsub(:,:), dfsIncremental(:)
    real(8), allocatable :: HBHtMatrixForOutput(:,:,:), dfsIncrementalForOutput(:,:)
    integer, allocatable :: order(:),orderForOutput(:,:)
    real(8) :: dfs
    real(8), allocatable :: dfsForOutput(:)
    logical :: llok
    integer :: channelNumber1, channelNumber2, channelIndex1, channelIndex2,stringIndex
    integer :: numHeader, numHeaderMaxMpi, sensorIndex
    integer :: countObs, sumCountObsMpi, maxCountObsMpi, countChannel, maxCountChannelMpi
    real(8), allocatable :: perturbationVector(:)
    integer :: dfsCount, sizeSelect, nulDfs, nulSelec, nulHBHt
    integer, parameter :: stringLength=128
    character(len=stringLength) :: headerObs
    character(len=stringLength), allocatable :: headerObsForOutput(:)
    integer, allocatable :: stringInt(:), stringIntForOutput(:,:)
    !
    !- 1.  Initialization

    write(*,*)
    write(*,*) 'Computing HBHT from selected observations start'
    call utl_printTime()
    
    vco_anl => col_getVco(columnTrlOnAnlIncLev)
    !- 1.3 Create a gridstateVector to store the perturbations
    call gsv_allocate(stateVector, tim_nstepobsinc, hco_anl, vco_anl, &
                      dataKind_opt = pre_incrReal, mpi_local_opt = .true.)

    !- 1.4 Create column vectors to store the perturbation interpolated to obs horizontal locations
    call col_setVco(columnAnlInc, vco_anl)
    call col_allocate(columnAnlInc, col_getNumCol(columnTrlOnAnlIncLev))

    !- 1.6
    call oti_timeBinning(obsSpaceData, tim_nstepobsinc)
    
    numHeader = obs_numHeader(obsSpaceData)
    call rpn_comm_allReduce(numHeader, numHeaderMaxMpi, 1, 'mpi_integer', 'mpi_max', 'grid', ierr)

    allocate(headerIndexList(numHeaderMaxMpi))
    allocate(levelList(numHeaderMaxMpi,nLevelsDfs))
    allocate(bodyIndexList(numHeaderMaxMpi,nLevelsDfs))
    allocate(stdDevList(numHeaderMaxMpi,nLevelsDfs))

    headerIndexList(:) = MPC_missingValue_INT
    levelList(:,:) = MPC_missingValue_INT
    bodyIndexList(:,:) = MPC_missingValue_INT
    stdDevList(:,:) = MPC_missingValue_R8

    ! Thinning and modifying the flag associated
    if (doThinning) then

      call thn_thinHyper(obsSpaceData)
      
      call obs_set_current_header_list(obsSpaceData,trim(familyType))
      HEADER0: do
        headerIndex = obs_getHeaderIndex(obsSpaceData)
        if (headerIndex < 0) exit HEADER0
   
        bodyIndexBeg = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
        bodyIndexEnd = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) + bodyIndexBeg - 1
        BODY0:do bodyIndex1 = bodyIndexBeg, bodyIndexEnd

          if ( btest(obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex1),11) ) then            
            call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex1, obs_notAssimilated)
          end if
          
        end do BODY0
      end do HEADER0
    end if ! if (doThinning)

    ! First step count the number of selected observation for each MPI task
    countObs = 0
    countChannel = 0 ! necessary in the case where no obs in the file
    
    call obs_set_current_header_list(obsSpaceData,trim(familyType))
    HEADER1: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER1
      if (.not. isSelected(headerIndex)) cycle HEADER1
      bodyIndexBeg = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
      bodyIndexEnd = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) + bodyIndexBeg - 1
      countChannel = 0
      BODY1:do bodyIndex1 = bodyIndexBeg, bodyIndexEnd
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
        write(*,*) 'found one observation with all requested ', nlevelsDfs, 'levels/channels available ', headerIndex
        countObs = countObs + 1
        headerIndexList(countObs) = headerIndex
        levelList(countObs,:) = vCoordList(1:nLevelsDfs)
        countChannel = 0
        BODY2:do bodyIndex1 = bodyIndexBeg, bodyIndexEnd
          llok = (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex1) == obs_assimilated)
          if (llok) then
            countChannel = countChannel + 1
            bodyIndexList(countObs, countChannel) = bodyIndex1
            stdDevList(countObs, countChannel) = obs_bodyElem_r(obsSpaceData, OBS_OER, bodyIndex1)
          end if
        end do BODY2
      end if
    end do HEADER1
    
    call rpn_comm_allReduce(countObs, sumCountObsMpi, 1, 'mpi_integer', 'mpi_sum', 'grid', ierr)

    call rpn_comm_allReduce(countObs, maxCountObsMpi, 1, 'mpi_integer', 'mpi_max', 'grid', ierr)

    call rpn_comm_allReduce(countChannel, maxCountChannelMpi, 1, 'mpi_integer', 'mpi_max', 'grid', ierr)

    if (.not. computeInParallel) then
      allocate(headerIndexListMpi(maxCountObsMpi, mmpi_nprocs))
      allocate(levelListMpi(maxCountObsMpi, maxCountChannelMpi, mmpi_nprocs))
      allocate(bodyIndexListMpi(maxCountObsMpi, maxCountChannelMpi, mmpi_nprocs))
      allocate(stdDevListMpi(maxCountObsMpi, maxCountChannelMpi, mmpi_nprocs))
    
      headerIndexListMpi(:,:) = MPC_missingValue_INT
      levelListMpi(:,:,:) = MPC_missingValue_INT
      bodyIndexListMpi(:,:,:) = MPC_missingValue_INT
      stdDevListMpi(:,:,:) = MPC_missingValue_R8
    
      call rpn_comm_allgather(headerIndexList(1:maxCountObsMpi), maxCountObsMpi, 'mpi_integer', &
                              headerIndexListMpi, maxCountObsMpi, 'mpi_integer', 'grid', ierr)

      call rpn_comm_allgather(levelList(1:maxCountObsMpi,1:maxCountChannelMpi), &
                              maxCountObsMpi * maxCountChannelMpi, 'mpi_integer',  &
                              levelListMpi, maxCountObsMpi * maxCountChannelMpi, 'mpi_integer', 'grid', ierr)

      call rpn_comm_allgather(bodyIndexList(1:maxCountObsMpi,1:maxCountChannelMpi), &
                              maxCountObsMpi * maxCountChannelMpi, 'mpi_integer',  &
                              bodyIndexListMpi, maxCountObsMpi * maxCountChannelMpi, 'mpi_integer', 'grid', ierr)

      call rpn_comm_allgather(stdDevList(1:maxCountObsMpi,1:maxCountChannelMpi), &
                              maxCountObsMpi * maxCountChannelMpi, 'mpi_real8',  &
                              stdDevListMpi, maxCountObsMpi * maxCountChannelMpi, 'mpi_real8', 'grid', ierr)

      call rpn_comm_barrier('GRID', ierr)
    
      deallocate(bodyIndexList)
      deallocate(levelList)
      deallocate(headerIndexList)
      deallocate(stdDevList)
    end if
    
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
    
    allocate(perturbationVector(cvm_nvadim))
    allocate(HBHtMatrix(nLevelsDfs,nLevelsDfs))
    HBHtMatrix(:,:) = MPC_missingValue_R8

    if (doChannelSelection) then
      sizeSelect = nLevelsDfs
      if (maxSelect > 0) then
        sizeSelect = maxSelect
      end if
      allocate(dfsIncremental(sizeSelect))
      allocate(order(sizeSelect))
    end if
    
    dfsCount = 0
    
    if (computeInParallel) then
      observationLoop1:do obsIndex = 1, maxCountObsMpi
        channelLoop1:do channelIndex1 = 1, maxCountChannelMpi
          ! We need to initialize the full OBS_WORK column to zero 
          do bodyIndex2 = 1, obs_numBody(obsSpaceData)
            call obs_bodySet_r(obsSpaceData, OBS_WORK, bodyIndex2, 0.d0)
          end do
          headerIndex = headerIndexList(obsIndex)
          if (headerIndex /= MPC_missingValue_INT) then
            bodyIndex1 = bodyIndexList(obsIndex,channelIndex1)
            call obs_bodySet_r(obsSpaceData, OBS_WORK, bodyIndex1, 1.d0)
          end if          
          call applyHBHtOperator(columnAnlInc, columnTrlOnAnlIncLev, stateVector, perturbationVector, obsSpaceData) 
         
          do channelIndex2 = 1, maxCountChannelMpi
            bodyIndex2 = bodyIndexList(obsIndex,channelIndex2)
            if (bodyIndex2 /= MPC_missingValue_INT) then !useless ?
              HBHtMatrix(channelIndex1,channelIndex2) = obs_bodyElem_r(obsSpaceData, OBS_WORK, bodyIndex2)
            end if
          end do
          write(*,*) 'diagDfs: computed column ', channelIndex1, 'of HBHt'
          call utl_printTime()
        end do channelLoop1
      
        headerObs = ''
        dfs = 0.d0
        if (headerIndex /= MPC_missingValue_INT) then
          dfsCount = dfsCount + 1
          call createHeaderString(obsSpaceData, headerIndex, familyType, headerObs)
          
          ! Extraction of the R matrix
          allocate(Rsub(nLevelsDfs,nLevelsDfs))
          Rsub(:,:) = MPC_missingValue_R8
          if (familyType == 'TO') then
            sensorIndex = tvs_lsensor(headerIndex)
            call rmat_getRmatrix(sensorIndex, &
                levelList(obsIndex,:),        &
                stdDevList(obsIndex,:),       &
                Rsub)
          else
            Rsub(:,:) = 0.d0
            do channelIndex2 = 1, nLevelsDfs
              Rsub(channelIndex2,channelIndex2) = stdDevList(obsIndex,channelIndex2) ** 2
            end do
          end if

          ! Computation of total dfs
          dfs = computeDfs(HBHtMatrix, Rsub)
         
          ! Calculate the selection of channels 
          if (doChannelSelection) then
            call selectChannels(HBHtMatrix, Rsub, dfsIncremental, order, maxSelect)
          end if
          deallocate(Rsub)
        end if
        
        allocate(headerObsForOutput(mmpi_nprocs))
        allocate(stringInt(stringLength))
        do stringIndex = 1, stringLength
          stringInt(stringIndex) = iachar(headerObs(stringIndex:stringIndex))
        end do
        allocate(stringIntForOutput(stringLength,mmpi_nprocs))
        call rpn_comm_gather(stringInt, stringLength, 'mpi_integer', &
            stringIntForOutput, stringLength, 'mpi_integer', 0, 'grid', ierr)
        do stringIndex = 1, stringLength
          do outTaskIndex = 1, mmpi_nprocs
            headerObsForOutput(outTaskIndex)(stringIndex:stringIndex) = achar(stringIntForOutput(stringIndex,outTaskIndex))
          end do
        end do
        
        if (outputHBHt) then
          allocate(HBHtMatrixForOutput(nLevelsDfs,nLevelsDfs,mmpi_nprocs))
          call rpn_comm_gather(HBHtMatrix, nLevelsDfs*nLevelsDfs, 'mpi_real8', &
              HBHtMatrixForOutput, nLevelsDfs*nLevelsDfs, 'mpi_real8', 0, 'grid', ierr)
          
          if (mmpi_myId == 0) then
            do outTaskIndex = 1, mmpi_nprocs
              if (len_trim(headerObsForOutput(outTaskIndex)) == 0) cycle
              write(nulHBHt,'(A)') trim(headerObsForOutput(outTaskIndex))
              do channelIndex2 = 1, nLevelsDfs
                channelNumber2 = vcoordList(channelIndex2)
                do channelIndex1 = 1, nLevelsDfs
                  channelNumber1 = vcoordList(channelIndex1)
                  write(nulHBHt,'(A4,1x,2(i12,1x),e14.6)') 'HBHt', channelNumber1, &
                        channelNumber2, HBHtMatrixForOutput(channelIndex1,channelIndex2,outTaskIndex)
                end do
              end do
              write(nulHBHt,*)
            end do
          end if
        
        end if ! if (outputHBHt)
        
        allocate(dfsForOutput(mmpi_nprocs))
        call rpn_comm_gather(dfs, 1, 'mpi_real8', dfsForOutput, 1, 'mpi_real8', 0, 'grid', ierr)
        if (mmpi_myId == 0) then
          do outTaskIndex = 1, mmpi_nprocs
            if (len_trim(headerObsForOutput(outTaskIndex)) > 0) then
              write(nulDfs,'(A,1x,e14.6)') trim(headerObsForOutput(outTaskIndex)), dfsForOutput(outTaskIndex)
            end if
          end do
        end if
      
        if (doChannelSelection) then
          allocate(dfsIncrementalForOutput(sizeSelect,mmpi_nprocs))
          call rpn_comm_gather(dfsIncremental, sizeSelect, 'mpi_real8', &
              dfsIncrementalForOutput, sizeSelect, 'mpi_real8', 0, 'grid', ierr)
          allocate(orderForOutput(sizeSelect,mmpi_nprocs))
          call rpn_comm_gather(order, sizeSelect, 'mpi_integer', &
              orderForOutput, sizeSelect, 'mpi_integer', 0, 'grid', ierr)
          if (mmpi_myId == 0) then
            do outTaskIndex = 1, mmpi_nprocs
              if (len_trim(headerObsForOutput(outTaskIndex)) > 0) then
                write(nulSelec,'(A)') trim(headerObsForOutput(outTaskIndex))
                do channelIndex1 = 1, sizeSelect
                  write(nulSelec,'(3(i5,1x),e14.6)')  channelIndex1, orderForOutput(channelIndex1,outTaskIndex), &
                      vcoordList(orderForOutput(channelIndex1,outTaskIndex)), dfsIncrementalForOutput(channelIndex1,outTaskIndex)
                end do
                write(nulSelec,*)
              end if
            end do
          end if
        end if ! if (doChannelSelection)
        deallocate(headerObsForOutput)
        if (allocated(HBHtMatrixForOutput)) deallocate(HBHtMatrixForOutput)
        if (allocated(dfsIncrementalForOutput)) deallocate(dfsIncrementalForOutput)
        if (allocated(orderForOutput)) deallocate(orderForOutput)
        deallocate(dfsForOutput)
        deallocate(stringInt)
        deallocate(stringIntforOutput)
        if (dfsCount == nDfsMax) exit observationLoop1
      end do observationLoop1
     
     else ! if (computeInParallel)
      
      mpiTaskLoop:do procIndex = 1, mmpi_nprocs
        observationLoop2:do obsIndex = 1, maxCountObsMpi
          headerIndex = headerIndexListMpi(obsIndex,procIndex)
          if (headerIndex == MPC_missingValue_INT) cycle
          taskIndex = procIndex - 1
          channelLoop2:do channelIndex1 = 1, maxCountChannelMpi
            bodyIndex1 = bodyIndexListMpi(obsIndex,channelIndex1,procIndex)
            channelNumber1 = levelListMpi(obsIndex,channelIndex1,procIndex)
            ! We need to initialize the full OBS_WORK column to zero 
            do bodyIndex2 = 1, obs_numBody(obsSpaceData)
              call obs_bodySet_r(obsSpaceData, OBS_WORK, bodyIndex2, 0.d0)
            end do
            if (mmpi_myId == taskIndex) call obs_bodySet_r(obsSpaceData, OBS_WORK, bodyIndex1, 1.d0)
            call applyHBHtOperator(columnAnlInc, columnTrlOnAnlIncLev, stateVector, perturbationVector, obsSpaceData)
              
            do channelIndex2 = 1, maxCountChannelMpi
              bodyIndex2 = bodyIndexListMpi(obsIndex,channelIndex2,procIndex)
              if (mmpi_myId == taskIndex .and. bodyIndex2 /= MPC_missingValue_INT) then
                HBHtMatrix(channelIndex1,channelIndex2) = obs_bodyElem_r(obsSpaceData, OBS_WORK, bodyIndex2)
              end if
            end do
            write(*,*) 'diagDfs: computed column ', channelIndex1, 'of HBHt'
            call utl_printTime()
          end do channelLoop2
          dfsCount = dfsCount + 1

          if (mmpi_myId == taskIndex) then
            call createHeaderString(obsSpaceData, headerIndex, familyType, headerObs)
            
            ! Extraction of the R matrix
            allocate(Rsub(nLevelsDfs,nLevelsDfs))
            if (familyType == 'TO') then
              sensorIndex = tvs_lsensor(headerIndex)
              call rmat_getRmatrix(sensorIndex,        &
                  levelListMpi(obsIndex,:,procIndex),  &
                  stdDevListMpi(obsIndex,:,procIndex), &
                  Rsub)
            else
              Rsub(:,:) = 0.d0
              do channelIndex2 = 1, nLevelsDfs
                Rsub(channelIndex2,channelIndex2) = stdDevListMpi(obsIndex,channelIndex2,procIndex) ** 2
              end do
            end if

            ! Computation of total dfs
            dfs = computeDfs(HBHtMatrix, Rsub)
         
            ! Calculate the selection of channels 
            if (doChannelSelection) then
              call selectChannels(HBHtMatrix, Rsub, dfsIncremental, order, maxSelect)
            end if
            deallocate(Rsub)
          end if
        
          call rpn_comm_bcastc(headerObs, stringLength, 'MPI_CHARACTER', taskIndex, 'GRID', ierr)

          if (outputHBHt) then
            call rpn_comm_bcast(HBHtMatrix, nLevelsDfs*nLevelsDfs, 'MPI_REAL8', taskIndex, 'GRID', ierr)
            if (mmpi_myId == 0) then
              write(nulHBHt,'(A)') trim(headerObs)
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
          if (mmpi_myId == 0) write(nulDfs,'(A,1x,e14.6)') trim(headerObs), dfs
          if (doChannelSelection) then
            call rpn_comm_bcast(dfsIncremental, sizeSelect, 'MPI_REAL8', taskIndex, 'GRID', ierr)
            call rpn_comm_bcast(order, sizeSelect, 'MPI_INTEGER', taskIndex, 'GRID', ierr)
            if (mmpi_myId == 0) then
              write(nulSelec,'(A)') trim(headerObs)
              do channelIndex1 = 1, size(order)
                write(nulSelec,'(3(i5,1x),e14.6)')  channelIndex1, order(channelIndex1), &
                    levelListMpi(obsIndex,order(channelIndex1),procIndex), dfsIncremental(channelIndex1)
              end do
              write(nulSelec,*)
            end if
          end if
          if (dfsCount == nDfsMax) exit mpiTaskLoop
        end do observationLoop2
      end do mpiTaskLoop
    end if ! if (computeInParallel)

    if (mmpi_myId == 0) then
      if (outputHbHt) then
        ierr = fclos(nulhbht)
      end if
      if (doChannelSelection) then
        ierr = fclos(nulSelec)
      end if
      ierr = fclos(nuldfs)
    end if

    
    if (allocated(bodyIndexList))  deallocate(bodyIndexList)
    if (allocated(levelList))  deallocate(levelList)
    if (allocated(headerIndexList))  deallocate(headerIndexList)
    if (allocated(stdDevList))  deallocate(stdDevList)
    
    
    if (allocated(Rsub)) deallocate(Rsub)
    if (allocated(dfsIncremental)) deallocate(dfsIncremental)
    if (allocated(order)) deallocate(order)
    if (allocated(HBHtMatrix)) deallocate(HBHtMatrix)
    
    deallocate(perturbationVector)
    call col_deallocate(columnAnlInc)
    call gsv_deallocate(stateVector)
    
    call msg_memUsage('midas-dfs')
    write(*,*)
    write(*,*) 'Computing DFS from selected observations end'
    
  end subroutine diagDFS


  !--------------------------------------------------------------------------
  ! applyHBHtOperator
  !--------------------------------------------------------------------------
  subroutine applyHBHtOperator(columnAnlInc, columnTrlOnAnlIncLev, stateVector, perturbationVector, obsSpaceData)
    !
    !:Purpose: apply chain of operators to apply HBHt (input and output in obsSpaceData OBS_WORK)
    !
    implicit none
    
    ! Arguments
    type(struct_columnData), intent(inout) :: columnAnlInc          ! Analysis increment as a column
    type(struct_columnData), intent(inout) :: columnTrlOnAnlIncLev  ! Trial field interpolated on the analysis increment as a column
    type(struct_gsv),        intent(inout) :: stateVector           ! State vector
    real(8),                 intent(inout) :: perturbationVector(:) ! Vector of perturbations            
    type(struct_obs),        intent(inout) :: obsSpaceData          ! Observation-related data structure

    ! Locals
    integer       :: localDimension
    logical, save :: firstCall=.true.
    
    call msg_memUsage('midas-dfs')
    
    localDimension = size(perturbationvector)
    call col_zero(columnAnlInc)
    call oop_Had(columnAnlInc, & !output
        columnTrlOnAnlIncLev,  &
        obsSpaceData,          & !input
        initializeLinearization_opt=firstCall)
    firstCall = .false.
    call gsv_zero(stateVector)
    call s2c_ad(stateVector,  & ! output
        columnAnlInc,         & ! input
        columnTrlOnAnlIncLev, &
        obsSpaceData)
    perturbationVector(:) = 0.d0
    call bmat_sqrtBT(perturbationVector, & !output
        localDimension,                  &  
        stateVector)                       !input
    call gsv_zero(stateVector)
    call bmat_sqrtB(perturbationVector, & !input
        localDimension,                 &
        stateVector)                      !output
    call s2c_tl(stateVector,  & !input
        columnAnlInc,         & !output
        columnTrlOnAnlIncLev, & 
        obsSpaceData)
    call oop_Htl(columnAnlInc, & !input
        columnTrlOnAnlIncLev,  &
        obsSpaceData,          & !output
        min_nsim = 1, initializeLinearization_opt = .false.)
    
  end subroutine applyHBHtOperator
  
 !--------------------------------------------------------------------------
  ! createHeaderString
  !--------------------------------------------------------------------------
  subroutine createHeaderString(obsSpaceData, headerIndex, familyType, headerString)
    !
    !:Purpose: create header string to be writen into output ascii files
    !
    implicit none
    
    ! Arguments
    type(struct_obs), intent(inout) :: obsSpaceData ! Observation-related data structure
    integer,          intent(in)    :: headerIndex  ! Current observation index in header obsSpacedata structure
    character(len=*), intent(in)    :: familyType   ! Observation family
    character(len=*), intent(out)   :: headerString ! Output header string

    ! Locals
    character(len=16) :: headerObs
   
    write(headerString,"('# ',A12,1x,2e14.6,1x,i8.8,1x,i4.4)")                          &
        obs_elem_c(obsSpaceData, 'STID', headerIndex),                                  &
        obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex) * MPC_DEGREES_PER_RADIAN_R8, &
        obs_headElem_r(obsSpaceData, OBS_LON, headerIndex) * MPC_DEGREES_PER_RADIAN_R8, &
        obs_headElem_i(obsSpaceData, OBS_DAT, headerIndex),                             &
        obs_headElem_i(obsSpaceData, OBS_ETM, headerIndex)
    if (familyType == 'TO')  then
      write(headerObs,'(1x,e14.6)') obs_headElem_r(obsSpaceData, OBS_SZA, headerIndex)
      headerString = trim(headerString) // trim(headerObs) 
    end if
    
  end subroutine createHeaderString
  
  !--------------------------------------------------------------------------
  ! subsetMatrix
  !--------------------------------------------------------------------------
  subroutine subsetMatrix(matrixInput, order, matrixOutput)
    !
    !:Purpose: extract sub-Matrix for a subset of channels/levels
    !
    implicit none
    
    ! Arguments:
    real(8), intent(in)  :: matrixInput(:,:)                         ! Initial matrix
    integer, intent(in)  :: order(:)                                 ! Ordered list of channels wanted to be extracted
    real(8), intent(out) :: matrixOutput(size(order), size(order))   ! Subset of the initial matrix with columns extracted

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
    real(8), intent(in) :: HBHt(:,:)  ! HBHt matrix
    real(8), intent(in) :: R(:,:)     ! Observation covariance matrix
    ! Result
    real(8)             :: dfs        ! Degrees of freedom for signal
    
    ! Locals
    integer :: nbLevels, levelIndex
    real(8), allocatable :: dMatrix(:,:), inverse(:,:), hk(:,:)

    nbLevels = size(R, dim = 1)
    allocate(dMatrix(nbLevels,nbLevels))
    allocate(inverse(nbLevels,nbLevels))
    allocate(hk(nbLevels,nbLevels))
    
    dMatrix(:,:) =  HBHt(:,:) + R(:,:)
    call utl_fastInverse(dMatrix, inverse)
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
    real(8),            intent(in)    :: HBHt(:,:)                 ! Matrix HBHt
    real(8),            intent(in)    :: R(:,:)                    ! Observation covariance matrix R
    real(8),            intent(out)   :: orderedDfs(:)             ! list of incremental DFS for each new added channel
    integer,            intent(out)   :: orderedChannelIndexes(:)  ! list of the channels selected in order
    integer, optional,  intent(in)    :: nChannelsOut_opt          ! number of channels wanted to be selected
    
    ! Locals
    integer, allocatable :: tmpOrder(:), freeIndexList(:), tmpFree(:)
    real(8)              :: optimalDfs, dfsTest
    integer              :: channelIndex1, channelIndex2, optimalIndex, numberOfFreeChannels
    real(8), allocatable :: RSubset(:,:), HBHtSubset(:,:)
    integer              :: nChannelsIn, nChannelsOut

    write(*,*) 'selectChannels: start'
    call msg_memUsage('selectChannels')
    call utl_printTime()
    nChannelsIn = size(R, dim = 1)
    
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
    
    orderedChannelIndexes(:) = -1
    allocate(tmpOrder(nChannelsOut))
    do channelIndex1 = 1, nChannelsOut
      optimalIndex = 0
      optimalDfs = 0.0d0
      allocate(RSubset(channelIndex1,channelIndex1), HBHtSubset(channelIndex1,channelIndex1))
      numberOfFreeChannels = nChannelsIn - channelIndex1 + 1
      write(*,*) 'selectChannels: step ', channelIndex1
      call utl_printTime()
      do channelIndex2 = 1, numberOfFreeChannels
        tmpOrder(1:channelIndex1-1) = orderedChannelIndexes(1:channelIndex1-1)
        tmpOrder(channelIndex1) = freeIndexList(channelIndex2)
        call subsetMatrix(R, tmpOrder(1:channelIndex1), RSubset)
        call subsetMatrix(HBHt, tmpOrder(1:channelIndex1), HBHtSubset)
        dfsTest = computeDfs(HBHtSubset, RSubset)
        if (dfsTest > optimalDfs) then
          optimalIndex = freeIndexList(channelIndex2)
          optimalDfs = dfsTest
        end if
      end do
      deallocate(RSubset, HBHtSubset)
      orderedDfs(channelIndex1) = optimalDfs
      orderedChannelIndexes(channelIndex1) = optimalIndex
      allocate(tmpFree(numberOfFreeChannels))
      !https://stackoverflow.com/questions/42140832/automatic-array-allocation-upon-assignment-in-fortran could help to simplify a bit this 
      !but it is needed to get rid of the (:) which in contradiction with our coding style
      tmpFree(:) = freeIndexList(:)
      call utl_reAllocate(freeIndexList, numberOfFreeChannels - 1)
      freeIndexList(:) = pack(tmpFree, tmpFree /= optimalIndex)
      deallocate(tmpFree)
    end do
    
    deallocate(freeIndexList)
    deallocate(tmpOrder)

  end subroutine selectChannels

  !--------------------------------------------------------------------------
  ! isSelected
  !--------------------------------------------------------------------------
  function isSelected(headerIndex) result(selected)
    !
    !:Purpose: is the observation in the list of specifically selected observations ?
    !          if empty list it is always
    !
    implicit none

    ! Arguments
    integer, intent(in) :: headerIndex ! Position in the header of obsSpaceData
    !Result
    logical             :: selected
    
    ! Locals
    integer            :: obsIndex
    real(8)            :: latitude(1), longitude(1), satelliteZenithAngle
    integer            :: date, hour
    real(8), parameter :: epsilon=0.01d0
    integer, parameter :: boxSize=7
    integer            :: definedConditions, satisfiedConditions
    integer            :: latIndex(1), lonIndex(1), observationIndex(1)
    real(8)            :: lowResWaterFraction(1), highResWaterFraction

    latitude(1) = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex) * MPC_DEGREES_PER_RADIAN_R8
    longitude(1) = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex) * MPC_DEGREES_PER_RADIAN_R8
    observationIndex(1)=headerIndex
    call tvs_interp_sfc (latIndex, lonIndex, 1, latitude, longitude, observationIndex, skipAlbedo_opt= .true.)  

    highResWaterFraction = tvs_surfaceParameters(headerIndex) % pcnt_wat

    ! Find the regional water fraction (here in a 15x15 pixel box centered on profile)

    call tvs_pcnt_box (lowResWaterFraction, 1 ,latIndex,lonIndex, boxSize)

    waterFractionTest : if (highResWaterFraction < highResWaterFractionThreshold .or. lowResWaterFraction(1) < lowResWaterFractionThreshold) then
      selected= .false.      
    else
      
      if (selectSpecificObservationsFromList) then
        selected = .false.
        satelliteZenithAngle = obs_headElem_r(obsSpaceData, OBS_SZA, headerIndex)
        date = obs_headElem_i(obsSpaceData, OBS_DAT, headerIndex)
        hour = obs_headElem_i(obsSpaceData, OBS_ETM, headerIndex)
        
        do obsIndex = 1, nObsMax
          definedConditions = 0
          satisfiedConditions = 0
          
          if (latList(obsindex) /= MPC_missingValue_R8) then
            definedConditions =  definedConditions + 1
            
            if (abs(latList(obsIndex)-latitude(1)) < epsilon) &
                satisfiedConditions = satisfiedConditions + 1
          end if
          
          if (lonList(obsindex) /= MPC_missingValue_R8) then
            definedConditions =  definedConditions + 1
            
            if (abs(lonList(obsIndex)-longitude(1)) < epsilon) &
                satisfiedConditions = satisfiedConditions + 1
          end if

          if (satZenList(obsindex) /= MPC_missingValue_R8) then
            definedConditions =  definedConditions + 1
            
            if (abs(satZenList(obsIndex)-satelliteZenithAngle) < epsilon) &
                satisfiedConditions = satisfiedConditions + 1
          end if
          
          if (dayList(obsindex) /= MPC_missingValue_INT) then
            definedConditions =  definedConditions + 1
            
            if (dayList(obsIndex) == date) &
                satisfiedConditions = satisfiedConditions + 1
          end if
          
          if (timeList(obsindex) /= MPC_missingValue_INT) then
            definedConditions =  definedConditions + 1
            if (timeList(obsIndex) == hour) &
                satisfiedConditions = satisfiedConditions + 1
          end if
          
          if (satisfiedConditions > 0 .and. definedConditions == satisfiedConditions) then
            selected = .true.
            return
          end if
          
        end do
      
      else
        selected = .true.
      end if
    end if waterFractionTest

  end function isSelected
    
end program midas_dfs

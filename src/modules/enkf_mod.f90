
module enkf_mod
  ! MODULE enkf_mod (prefix='enkf' category='1. High-level functionality')
  !
  !:Purpose:  Various routines that are useful for implementing
  !           an EnKF in MIDAS, including the LETKF.
  !
  use mpi                 ! this is the mpi library module
  use rpn_comm
  use midasMpi_mod
  use utilities_mod
  use message_mod
  use mathPhysConstants_mod
  use timeCoord_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use ensembleStateVector_mod
  use gridStateVector_mod
  use obsSpaceData_mod
  use tovs_mod
  use ensembleObservations_mod
  use gridVariableTransforms_mod
  use localizationFunction_mod
  use varNameList_mod
  use codePrecision_mod
  use codTyp_mod
  use calcHeightAndPressure_mod

  implicit none
  save
  private

  ! public types
  public :: struct_enkfInterpInfo
  public :: struct_enkfNML
  public :: struct_enkfDFS

  ! Public procedures
  public :: enkf_readNML
  public :: enkf_setupInterpInfo, enkf_LETKFanalyses, enkf_modifyAMSUBobsError
  public :: enkf_rejectHighLatIR, enkf_getModulatedState, enkf_setupModulationFactor
  public :: enkf_deallocateDFS

  ! Constants (private)(
  integer, parameter :: maxNumLocalize = 10 ! hLocalize size

  ! for weight interpolation
  type struct_enkfInterpInfo
    integer              :: latLonStep
    integer, allocatable :: numIndexes(:,:)
    integer, allocatable :: lonIndexes(:,:,:)
    integer, allocatable :: latIndexes(:,:,:)
    real(8), allocatable :: interpWeights(:,:,:)
    integer              :: myLonBegHalo
    integer              :: myLonEndHalo
    integer              :: myLatBegHalo
    integer              :: myLatEndHalo
    integer              :: myLonBeg
    integer              :: myLonEnd
    integer              :: myLatBeg
    integer              :: myLatEnd
  end type struct_enkfInterpInfo

  ! all namelist variables
  type struct_enkfNML
    character(len=20)  :: algorithm
    logical            :: ensPostProcessing
    logical            :: recenterInputEns
    integer            :: numSubEns
    character(len=256) :: ensPathName
    integer  :: nEns
    logical  :: randomShuffleSubEns
    logical  :: writeLocalEnsObsToFile
    integer  :: maxNumLocalObs
    integer  :: maxNumLocalObsPerType
    integer  :: weightLatLonStep
    real(8)  :: alphaRandomPertPrior
    integer  :: numRetainedEigen
    integer  :: myNumLatLonSendFactor
    logical  :: modifyAmsubObsError
    logical  :: backgroundCheck
    logical  :: huberize
    logical  :: rejectHighLatIR
    logical  :: rejectRadNearSfc
    logical  :: ignoreEnsDate
    logical  :: outputOnlyEnsMean
    logical  :: outputEnsObs
    integer  :: localSelectionOutput
    logical  :: outputDFS
    logical  :: debug
    logical  :: readEnsObsFromFile
    real(8)  :: hLocalize(maxNumLocalize)
    real(8)  :: hLocalizePressure(maxNumLocalize)
    logical  :: hLinearLoc
    real(8)  :: vLocalize
    real(8)  :: minDistanceToLand
    character(len=20) :: obsTimeInterpType
    character(len=20) :: mpiDistribution
    character(len=12) :: etiket_anl
    character(len=20) :: localObsSorting
    integer  :: fileMemberIndex1       = 1
    logical  :: readEnsMeanFromFile    = .false.
    integer  :: numFullEns             = 0
  end type struct_enkfNML

  ! for dfs calculation
  type struct_enkfDFS
    logical              :: allocated = .false.
    real(8), allocatable :: locFun(:,:)
    integer, allocatable :: bodyIndex(:,:)
    real(8), allocatable :: lat(:)
    real(8), allocatable :: lon(:)
    real(8), allocatable :: lnp(:)
    real(8), allocatable :: dfs(:)
  end type struct_enkfDFS

  integer, external :: get_max_rss

contains

  !----------------------------------------------------------------------
  ! enkf_readNML
  !----------------------------------------------------------------------
  subroutine enkf_readNML(enkfNML)
    !
    !:Purpose: Read the namelist and put variables in a the derived type
    !          variable supplied as an argument. Also, a few checks and
    !          modifications are done.
    !
    implicit none

    ! Arguments:
    type(struct_enkfNML) :: enkfNML  ! Derived type variable with namelist variables

    ! Locals:
    integer :: ierr, locIndex

    ! namelist variables
    character(len=20)  :: algorithm  ! name of the chosen LETKF algorithm: 'LETKF', 'CVLETKF'
    logical            :: ensPostProcessing ! do all post-processing of analysis ensemble
    logical            :: recenterInputEns  ! read a deterministic state to recenter ensemble
    integer            :: numSubEns  ! number of sub-ensembles to split the full ensemble
    character(len=256) :: ensPathName ! absolute or relative path to ensemble directory
    integer  :: nEns                 ! ensemble size
    logical  :: randomShuffleSubEns  ! choose to randomly shuffle members into subensembles
    logical  :: writeLocalEnsObsToFile ! Controls writing the ensObs to file.
    integer  :: maxNumLocalObs       ! maximum number of obs in each local volume to assimilate
    integer  :: maxNumLocalObsPerType ! maximum number of obs of each type in each local volume to assimilate
    integer  :: weightLatLonStep     ! separation of lat-lon grid points for weight calculation
    real(8)  :: alphaRandomPertPrior ! Random perturbation additive inflation coeff applied to trials (0->1)
    integer  :: numRetainedEigen     ! number of retained eigenValues/Vectors of vertical localization matrix
    integer  :: myNumLatLonSendFactor ! factor to obtain max number of grid points computed on each mpi task
    logical  :: modifyAmsubObsError  ! reduce AMSU-B obs error stddev in tropics
    logical  :: backgroundCheck      ! apply additional background check using ensemble spread
    logical  :: huberize             ! apply huber norm quality control procedure
    logical  :: rejectHighLatIR      ! reject all IR observations at high latitudes
    logical  :: rejectRadNearSfc     ! reject radiance observations near the surface
    logical  :: ignoreEnsDate        ! when reading ensemble, ignore the date
    logical  :: outputOnlyEnsMean    ! when writing ensemble, can choose to only write member zero
    logical  :: outputEnsObs         ! to write trial and analysis ensemble members in observation space to sqlite
    integer  :: localSelectionOutput ! write output about the local selection of observations
    logical  :: outputDFS            ! write output about the DFS
    logical  :: debug                ! debug option to print values to the listings.
    logical  :: readEnsObsFromFile   ! instead of computing innovations, read ensObs%Yb from file.
    real(8)  :: hLocalize(maxNumLocalize)         ! horizontal localization radius (in km)
    real(8)  :: hLocalizePressure(maxNumLocalize) ! pressures where horizontal localization changes (in hPa)
    logical  :: hLinearLoc           ! apply piece-wise linear vertical interpolation for the localization radius
    real(8)  :: vLocalize            ! vertical localization radius (units: ln(Pressure in Pa) or meters)
    real(8)  :: minDistanceToLand    ! for ice/ocean DA: minimum distance to land for assimilating obs
    character(len=20) :: obsTimeInterpType ! type of time interpolation to obs time
    character(len=20) :: mpiDistribution   ! type of mpiDistribution for weight calculation ('ROUNDROBIN' or 'TILES')
    character(len=12) :: etiket_anl  ! etiket for output files
    character(len=20) :: localObsSorting   ! method to sort observations in eob_getLocalBodyIndices() ('HORIZONTAL' (default), 'LOCFUN' or 'MINTRACE')
    integer  :: fileMemberIndex1     ! first member index in ensemble set to be read
    logical  :: readEnsMeanFromFile  ! choose to read ens mean from file (when reading subset of members)
    integer  :: numFullEns           ! number of full ensemble set (needed only for modulated ensemble)

    NAMELIST /NAMLETKF/algorithm, ensPostProcessing, recenterInputEns, nEns, numSubEns, &
                       ensPathName, randomShuffleSubEns,  &
                       hLocalize, hLocalizePressure, hLinearLoc, vLocalize, minDistanceToLand,  &
                       maxNumLocalObs, maxNumLocalObsPerType, weightLatLonStep, alphaRandomPertPrior,  &
                       modifyAmsubObsError, backgroundCheck, huberize, rejectHighLatIR, rejectRadNearSfc,  &
                       ignoreEnsDate, outputOnlyEnsMean, outputEnsObs, localSelectionOutput, outputDFS, &
                       obsTimeInterpType, mpiDistribution, etiket_anl, localObsSorting, &
                       readEnsObsFromFile, writeLocalEnsObsToFile, &
                       numRetainedEigen, myNumLatLonSendFactor, debug, &
                       fileMemberIndex1, readEnsMeanFromFile, numFullEns

    !- 1.1 Setting default namelist variable values
    algorithm                = 'LETKF'
    ensPostProcessing        = .false.
    recenterInputEns         = .false.
    ensPathName              = 'ensemble'
    nEns                     = 10
    numSubEns                = 2
    randomShuffleSubEns      = .false.
    maxNumLocalObs           = 1000
    maxNumLocalObsPerType    = 1000000
    weightLatLonStep         = 1
    alphaRandomPertPrior     = 0.0D0
    modifyAmsubObsError      = .false.
    backgroundCheck          = .false.
    huberize                 = .false.
    rejectHighLatIR          = .false.
    rejectRadNearSfc         = .false.
    ignoreEnsDate            = .false.
    outputOnlyEnsMean        = .false.
    outputEnsObs             = .false.
    localSelectionOutput     = 0
    outputDFS                = .false.
    hLocalize(:)             = -1.0D0
    hLocalizePressure(:)     = -1.0D0
    hLinearLoc               = .false.
    vLocalize                = -1.0D0
    minDistanceToLand        = -1.0D0
    obsTimeInterpType        = 'LINEAR'
    mpiDistribution          = 'ROUNDROBIN'
    etiket_anl               = 'ENS_ANL'
    localObsSorting          = 'HORIZONTAL'
    readEnsObsFromFile       = .false.
    writeLocalEnsObsToFile   = .false.
    numRetainedEigen         = 0
    myNumLatLonSendFactor    = 10
    debug                    = .false.
    fileMemberIndex1         = 1
    readEnsMeanFromFile      = .false.
    numFullEns               = 0

    if ( .not. utl_isNamelistPresent('NAMLETKF','./flnml') ) then
      call utl_abort('enkf_readNML: namletkf is missing in the namelist. ')
    else
      call utl_tmg_start(181,'low-level--readNML')
      read(utl_flnml, nml=namletkf, iostat=ierr)
      if ( ierr /= 0) call utl_abort('enkf_readNML: Error reading namelist')
      if ( mmpi_myid == 0 ) write(*,nml=namletkf)
      call utl_tmg_stop(181)
    end if

    ! Some minor modifications of namelist values

    ! default hLocalizePressure values assume four radii, or just one
    if (hLocalizePressure(1) < 0.0d0) then
      write(*,*) 'enkf_readNML: hLocalizePressure not set in namelist. Setting default values.'
      if (hLinearLoc) then
        hLocalizePressure(1:4)   = (/6.0d0, 144.0d0, 237.0d0, 700.0d0/) ! midpoints
      else
        hLocalizePressure(1:3)   = (/14.0d0, 140.0d0, 400.0d0/) ! transition values
      endif
    endif

    ! if only 1 value given for hLocalize, use it for the entire column
    if (hLocalize(1) > 0.0d0 .and. hLocalize(2) < 0.0d0) then
      hLocalize(:) = hLocalize(1)
      if ( mmpi_myid == 0 ) write(*,*) 'enkf_readNML: hLocalize is modified after reading namelist. ' // &
           'hLocalize(:)=', hLocalize(1)
      ! if no value give for hLocalize, abort
    else if ( hLocalize(1) > 0.0d0 ) then
      ! Check hLocalizePressure and hLocalize lengths consistency
      ! For a linearly varying localization radius, the radius is set for the hLocalizePressure values
      ! Therefore, hLocalizePressure has the same length as hLocalize
      ! For a step varying localization radius, hLocalizePressure values are the transition values between hLocalize values
      ! Therefore, hLocalizePressurec has one less value than hLocalize
      if ( (count(hLocalize > 0.0d0) /= count(hLocalizePressure > 0.0d0)     .and.       hLinearLoc) .or. &
           (count(hLocalize > 0.0d0) /= count(hLocalizePressure > 0.0d0) + 1 .and. .not. hLinearLoc) ) then
        write(*,*) 'enkf_readNML: hLocalize and hLocalizePressure have inconsistent lengths.'
        write(*,*) 'enkf_readNML: hLocalize has',count(hLocalize > 0.0d0),'positive values'
        write(*,*) 'enkf_readNML: hLocalizePressure has',count(hLocalizePressure > 0.0d0),'positive values'
        write(*,*) 'enkf_readNML: hLocalize = ',hLocalize(:)
        write(*,*) 'enkf_readNML: hLocalizePressure = ',hLocalizePressure(:)
        call utl_abort('enkf_readNML: hLocalize and hLocalizePressure inconsistency')
      endif
    end if

    do locIndex = 1,maxNumLocalize-1
      ! check if hLocalizePressure positive values decrease
      if ((hLocalizePressure(locIndex) >= hLocalizePressure(locIndex+1)) .and. &
           hLocalizePressure(locIndex+1) > 0.0d0) then
        write(*,*) 'enkf_readNML: hLocalizePressure = ',hLocalizePressure(:)
        call utl_abort('enkf_readNML: hLocalizePressure does not decrease')
      end if
    enddo

    do locIndex=1,maxNumLocalize
      ! convert localization radius from km to m
      if (hLocalize(locIndex) > 0.0d0) then
        hLocalize(locIndex) = hLocalize(locIndex) * 1000.0d0
      end if
      ! convert pressure from hPa to log(Pa)
      if (hLocalizePressure(locIndex) > 0.0d0) then
        hLocalizePressure(locIndex) = log(hLocalizePressure(locIndex) * MPC_PA_PER_MBAR_R8)
      end if
    enddo

    if ( mmpi_myid == 0 ) then
      write(*,*) 'enkf_readNML: hLocalize (meters):',hlocalize
      write(*,*) 'enkf_readNML: hLocalizePressure (log(Pa)):',hlocalizePressure
    endif

    if (minDistanceToLand > 0.0D0) then
      minDistanceToLand = minDistanceToLand * 1000.0D0 ! convert from km to m
    end if

    if (trim(algorithm) /= 'LETKF'           .and. &
        trim(algorithm) /= 'CVLETKF'         .and. &
        trim(algorithm) /= 'CVLETKF-PERTOBS' .and. &
        trim(algorithm) /= 'LETKF-Gain'      .and. &
        trim(algorithm) /= 'LETKF-Gain-ME'   .and. &
        trim(algorithm) /= 'CVLETKF-ME') then
      call utl_abort('enkf_readNML: unknown LETKF algorithm: ' // trim(algorithm))
    end if

    if (numRetainedEigen < 0) then
      call utl_abort('enkf_readNML: numRetainedEigen should be ' // &
                     'equal or greater than zero')
    end if

    ! Copy all variables to the dervied type variable enkfNML
    enkfNML%algorithm              = algorithm
    enkfNML%ensPostProcessing      = ensPostProcessing
    enkfNML%recenterInputEns       = recenterInputEns
    enkfNML%ensPathName            = ensPathName
    enkfNML%nEns                   = nEns
    enkfNML%numSubEns              = numSubEns
    enkfNML%randomShuffleSubEns    = randomShuffleSubEns
    enkfNML%maxNumLocalObs         = maxNumLocalObs
    enkfNML%maxNumLocalObsPerType  = maxNumLocalObsPerType
    enkfNML%weightLatLonStep       = weightLatLonStep
    enkfNML%alphaRandomPertPrior   = alphaRandomPertPrior
    enkfNML%modifyAmsubObsError    = modifyAmsubObsError
    enkfNML%backgroundCheck        = backgroundCheck
    enkfNML%huberize               = huberize
    enkfNML%rejectHighLatIR        = rejectHighLatIR
    enkfNML%rejectRadNearSfc       = rejectRadNearSfc
    enkfNML%ignoreEnsDate          = ignoreEnsDate
    enkfNML%outputOnlyEnsMean      = outputOnlyEnsMean
    enkfNML%outputEnsObs           = outputEnsObs
    enkfNML%localSelectionOutput   = localSelectionOutput
    enkfNML%outputDFS              = outputDFS
    enkfNML%hLocalize(:)           = hLocalize(:)
    enkfNML%hLocalizePressure(:)   = hLocalizePressure(:)
    enkfNML%vLocalize              = vLocalize
    enkfNML%minDistanceToLand      = minDistanceToLand
    enkfNML%obsTimeInterpType      = obsTimeInterpType
    enkfNML%mpiDistribution        = mpiDistribution
    enkfNML%etiket_anl             = etiket_anl
    enkfNML%localObsSorting        = localObsSorting
    enkfNML%readEnsObsFromFile     = readEnsObsFromFile
    enkfNML%writeLocalEnsObsToFile = writeLocalEnsObsToFile
    enkfNML%numRetainedEigen       = numRetainedEigen
    enkfNML%myNumLatLonSendFactor  = myNumLatLonSendFactor
    enkfNML%debug                  = debug
    enkfNML%fileMemberIndex1       = fileMemberIndex1
    enkfNML%readEnsMeanFromFile    = readEnsMeanFromFile
    enkfNML%numFullEns             = numFullEns

  end subroutine enkf_readNML

  !----------------------------------------------------------------------
  ! enkf_LETKFanalyses
  !----------------------------------------------------------------------
  subroutine enkf_LETKFanalyses(enkfNML,  &
                                ensembleAnl, ensembleTrl, &
                                ensObs_mpiglobal, ensObsGain_mpiglobal, &
                                stateVectorMeanAnl, &
                                wInterpInfo, enkfDFS)
    !
    !:Purpose: Local subroutine containing the code for computing
    !          the LETKF analyses for all ensemble members, ensemble
    !          mean. Two MPI distributions are currently available for
    !          distributing the calculation of the LETKF weights at each
    !          grid point:
    !
    !           - 'MASTERWORKER': Dynamic load balancing with the first MPI task acting
    !             as the 'MASTER' with the sole job of assigning work to all of the other
    !             MPI tasks. When all of the LETKF weights are computed, then the MASTER
    !             send a signal to all workers that we are done and can move on to the next
    !             part of the algorithm.
    !           - 'ROUNDROBIN': Deal out all grid points where weights are computed to all
    !             mpi tasks without considering where they are needed for use. This should
    !             provide the best load balance in terms of the cost of covariance
    !             calculation which can vary depending on the local density of observations,
    !             but requires the most communication.
    !
    implicit none

    ! Arguments:
    type(struct_enkfNML),        intent(in)    :: enkfNML              ! Derived type variable with namelist variables
    type(struct_ens), pointer,   intent(inout) :: ensembleTrl          ! Trial ensemble
    type(struct_ens),            intent(inout) :: ensembleAnl          ! Analysis ensemble
    type(struct_eob), target,    intent(in)    :: ensObs_mpiglobal     ! Ensemble observations for original ensemble
    type(struct_eob),            intent(in)    :: ensObsGain_mpiglobal ! Ensemble observations for computing gain matrix
    type(struct_gsv),            intent(inout) :: stateVectorMeanAnl   ! Ensemble mean state vector
    type(struct_enkfInterpInfo), intent(in)    :: wInterpInfo          ! Derived type variable weight interpolation info
    type(struct_enkfDFS),        intent(inout) :: enkfDFS

    ! Locals:
    character :: readySignal
    integer :: workerProcID, finishedSignal, assignmentTag, readyTag, numFinished
    integer :: mpiStatus(MPI_STATUS_SIZE)
    integer, allocatable :: waitStatusesSend(:,:), waitStatusesRecv(:,:)
    integer :: latLonIndex, nEnsGain, nLev_weights, ierr
    integer :: procIndex, procIndexSend, latLonIndexMpiGlobal
    integer :: latIndex, lonIndex, levIndex
    integer :: countMaxExceeded, maxCountMaxExceeded, numGridPointWeights
    integer :: myNumLatLonRecv, numLatLonMpiGlobal, myNumLatLonCalcMax, myNumLatLonSendMax
    integer :: sendTag, recvTag, nsize, numRecv, numSend
    integer :: myLonBegHalo, myLonEndHalo, myLatBegHalo, myLatEndHalo
    integer :: requestIdRecvFinished(mmpi_nprocs-1), requestIdSendFinished(mmpi_nprocs-1)
    integer :: requestIdSignal
    integer, allocatable :: myLatIndexesRecv(:), myLonIndexesRecv(:)
    integer, allocatable :: latIndexesSendMpiGlobal(:), lonIndexesSendMpiGlobal(:)
    integer, allocatable :: numProcsSendMpiGlobal(:)
    integer, allocatable :: procIndexesSendMpiGlobal(:,:)
    integer, allocatable :: requestIdRecv(:), requestIdSend(:)
    integer, allocatable :: latLonTagMpiGlobal(:,:)
    real(8), allocatable :: weightsMembers(:,:,:,:), weightsMembersLatLon(:,:,:)
    real(8), allocatable :: weightsMean(:,:,:,:), weightsMeanLatLon(:,:,:)
    real(4), allocatable :: vertLocation_r4(:,:,:)
    real(8), allocatable :: weightsSendCombined(:,:,:), weightsRecvCombined(:,:,:)
    type(struct_hco), pointer :: hco_ens
    type(struct_vco), pointer :: vco_ens
    type(struct_gsv)          :: stateVectorMeanInc
    type(struct_gsv)          :: stateVectorMeanTrl
    logical :: procAlreadyFinished(mmpi_nprocs)
    logical :: useModulatedEns, masterIsMe

    call utl_tmg_start(131,'--LETKFanalysis')

    write(*,*)
    write(*,*) 'enkf_LETKFanalyses: starting'
    call msg_memUsage('enkf_LETKFanalyses')

    useModulatedEns = (enkfNML%numRetainedEigen > 0)
    if ( useModulatedEns ) then
      nEnsGain   = enkfNML%nEns * enkfNML%numRetainedEigen
    else
      nEnsGain   = enkfNML%nEns
    end if

    ! Set things up for the redistribution of work across mpi tasks
    call enkf_LETKFsetupMpiDistribution(numLatLonMpiGlobal, myNumLatLonRecv,  &
                                        myLatIndexesRecv, myLonIndexesRecv,   &
                                        latIndexesSendMpiGlobal, lonIndexesSendMpiGlobal,   &
                                        procIndexesSendMpiGlobal, &
                                        numProcsSendMpiGlobal, wInterpInfo)

    write(*,*) 'enkf_LETKFanalyses: numLatLonMpiGlobal, maxval(numProcsSendMpiGlobal) = ', &
                                    numLatLonMpiGlobal, maxval(numProcsSendMpiGlobal)

    ! Compute maximum expected number of grid points where weights computed on each mpi task
    if (trim(enkfNML%mpiDistribution) == 'ROUNDROBIN') then
      myNumLatLonCalcMax = ceiling(real(numLatLonMpiGlobal)/real(mmpi_nprocs))
    else if (trim(enkfNML%mpiDistribution) == 'MASTERWORKER') then
      myNumLatLonCalcMax = enkfNML%myNumLatLonSendFactor*ceiling(real(numLatLonMpiGlobal)/real(mmpi_nprocs))
    else
      write(*,*) 'mpiDistribution = ', trim(enkfNML%mpiDistribution)
      call utl_abort('enkf_LETKFanalyses: unknown mpiDistribution')
    end if
    myNumLatLonSendMax = ceiling(real(maxval(numProcsSendMpiGlobal))*real(myNumLatLonCalcMax))
    write(*,*) 'enkf_LETKFanalyses: myNumLatLonCalcMax, myNumLatLonSendMax = ',  &
                                    myNumLatLonCalcMax, myNumLatLonSendMax

    ! Allocate and initialize various arrays needed for MPI communication
    allocate(requestIdSend(myNumLatLonSendMax))
    allocate(requestIdRecv(myNumLatLonRecv))
    allocate(waitStatusesSend(MPI_STATUS_SIZE,myNumLatLonSendMax))
    allocate(waitStatusesRecv(MPI_STATUS_SIZE,myNumLatLonRecv))
    allocate(weightsRecvCombined(nEnsGain, enkfNML%nEns+1, myNumLatLonRecv))
    allocate(weightsSendCombined(nEnsGain, enkfNML%nEns+1, myNumLatLonCalcMax))
    requestIdRecv(:) = 0
    requestIdSend(:) = 0
    waitStatusesRecv(:,:) = 0
    waitStatusesSend(:,:) = 0
    weightsRecvCombined(:,:,:) = 0.0d0
    weightsSendCombined(:,:,:) = 0.0d0

    ! Initialize dimension related variables
    nLev_weights = max(ens_getNumLev(ensembleAnl,'MM'), ens_getNumLev(ensembleAnl,'DP'))
    if ( useModulatedEns ) nLev_weights = 1
    hco_ens => ens_getHco(ensembleAnl)
    vco_ens => ens_getVco(ensembleAnl)
    myLonBegHalo = wInterpInfo%myLonBegHalo
    myLonEndHalo = wInterpInfo%myLonEndHalo
    myLatBegHalo = wInterpInfo%myLatBegHalo
    myLatEndHalo = wInterpInfo%myLatEndHalo

    ! Allocate weights for mean analysis
    allocate(weightsMean(nEnsGain,1,myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
    weightsMean(:,:,:,:) = 0.0d0
    allocate(weightsMeanLatLon(nEnsGain,1,myNumLatLonCalcMax))
    weightsMeanLatLon(:,:,:) = 0.0d0
    ! Allocate weights for member analyses
    allocate(weightsMembers(nEnsGain,enkfNML%nEns,myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
    weightsMembers(:,:,:,:) = 0.0d0
    allocate(weightsMembersLatLon(nEnsGain,enkfNML%nEns,myNumLatLonCalcMax))
    weightsMembersLatLon(:,:,:) = 0.0d0

    if (enkfNML%outputDFS) call enkf_allocateDFS(enkfDFS, nLev_weights, myNumLatLonCalcMax, enkfNML%maxNumLocalObs)

    ! Allocate and initialize state vectors for mean trial and increment
    call gsv_allocate( stateVectorMeanTrl, tim_nstepobsinc, hco_ens, vco_ens, &
                       dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_zero(stateVectorMeanTrl)
    call gsv_allocate( stateVectorMeanInc, tim_nstepobsinc, hco_ens, vco_ens, &
                       dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_zero(stateVectorMeanInc)

    ! Compute trial ensemble mean and copy to state vector
    call ens_computeMean(ensembleTrl)
    call ens_copyEnsMean(ensembleTrl, stateVectorMeanTrl)

    ! Compute 3D field of vertical location needed for localization
    call lfn_Setup(LocFunctionWanted='FifthOrder')
    if (enkfNML%vLocalize > 0.0d0) then
      call enkf_computeVertLocation(vertLocation_r4,stateVectorMeanTrl)
    end if

    call utl_tmg_start(141,'----Barr')
    call rpn_comm_barrier('GRID',ierr)
    call utl_tmg_stop(141)

    ! Get mpi global list of tags used for mpi send/recv
    call utl_tmg_start(142, '----GetGlobalTags')
    allocate(latLonTagMpiGlobal(stateVectorMeanAnl%ni,stateVectorMeanAnl%nj))
    call enkf_LETKFgetMpiGlobalTags(latLonTagMpiGlobal,myLatIndexesRecv,myLonIndexesRecv)
    call utl_tmg_stop(142)

    ! Start of major loop for computing weights for ensemble mean and members
    countMaxExceeded = 0
    maxCountMaxExceeded = 0
    numGridPointWeights = 0
    LEV_LOOP: do levIndex = 1, nLev_weights
      write(*,*) 'Computing ensemble updates for vertical level = ', levIndex
      call utl_printTime(reset_opt = (levIndex==1))

      ! First post all recv instructions for communication of weights
      call utl_tmg_start(132,'----CommWeights')
      call utl_tmg_start(148,'----CommWeights-irecv')
      numSend = 0
      numRecv = 0
      do latLonIndex = 1, myNumLatLonRecv
        latIndex = myLatIndexesRecv(latLonIndex)
        lonIndex = myLonIndexesRecv(latLonIndex)
        recvTag = latLonTagMpiGlobal(lonIndex,latIndex) + (levIndex-1)*maxval(latLonTagMpiGlobal)
        recvTag = 1 + mod(recvTag-1, mmpi_maxTagValue - 10)

        nsize = nEnsGain * (enkfNML%nEns+1)
        numRecv = numRecv + 1
        weightsRecvCombined(:,:,latLonIndex) = -999.8d0
        call mpi_irecv( weightsRecvCombined(:,:,latLonIndex),  &
                        nsize, mmpi_datyp_real8, mpi_any_source, recvTag,  &
                        mmpi_comm_grid, requestIdRecv(numRecv), ierr )
      end do
      call utl_tmg_stop(148)
      call utl_tmg_stop(132)

      ! Set tag values for sending signals
      assignmentTag = mmpi_maxTagValue - 1
      readyTag      = mmpi_maxTagValue - 2

      ! Determine if I am the 'Master' (only for MASTERWORKER mpi distribution)
      masterIsMe = (mmpi_myid == 0) .and. (trim(enkfNML%mpiDistribution) == 'MASTERWORKER')

      MASTER_WORKER: if (masterIsMe) then  ! I am the master, I do no actual computations

        call utl_tmg_start(132,'----CommWeights')
        call utl_tmg_start(145,'----CommWeights-signals')

        procAlreadyFinished(:) = .false.

        write(*,*) 'Send all assignments'
        call utl_printTime()

        ! Loop over all gridpoints where calculations need to be performed
        write(*,*) 'Start of loop over all global grid points where weights computed'
        do latLonIndexMpiGlobal = 1, numLatLonMpiGlobal

          ReadyLoop: do ! Loop until we get a valid read signal

            ! Determine which MPI task is ready for a new work assignment
            call MPI_RECV(readySignal, 1, MPI_CHARACTER, mpi_any_source, readyTag, &
                          mmpi_comm_grid, mpiStatus, ierr)
            workerProcID = mpiStatus(MPI_SOURCE)

            if (readySignal == 'N') then
              procAlreadyFinished(workerProcID+1) = .true.
            end if

            ! Assign this MPI task the next gridpoint to be calculated
            ! NOTE: this assignment will be ignored when readySignal is 'N'
            call MPI_RSEND(latLonIndexMpiGlobal, 1, MPI_INTEGER, workerProcID, assignmentTag, &
                           mmpi_comm_grid, ierr)

            if (readySignal == 'Y') then
              exit ReadyLoop
            end if

          end do ReadyLoop

        end do

        ! Now that all work is done, we need to inform all workers
        write(*,*) 'Finished all grid points, send *finished* signal to all mpi tasks'
        call utl_printTime()

        numFinished = 0
        do procIndex = 2, mmpi_nprocs

          ! Skip this proc if he already told us that he is finished
          if (procAlreadyFinished(procIndex)) cycle

          numFinished = numFinished + 1

          ! Wait for signal for every task that last assignment is complete
          call MPI_IRECV(readySignal, 1, MPI_CHARACTER, procIndex-1, readyTag, &
                         mmpi_comm_grid, requestIdRecvFinished(numFinished), ierr)

          ! Tell this worker we are done
          finishedSignal = 0
          call MPI_ISEND(finishedSignal, 1, MPI_INTEGER, procIndex-1, assignmentTag, &
                         mmpi_comm_grid, requestIdSendFinished(numFinished), ierr)

        end do
        call mpi_waitAll(numFinished, requestIdSendFinished(1:numFinished), MPI_STATUSES_IGNORE, ierr)
        call mpi_waitAll(numFinished, requestIdRecvFinished(1:numFinished), MPI_STATUSES_IGNORE, ierr)

        write(*,*) 'All *finished* signals received and acknowledged'
        call utl_printTime()

        call utl_tmg_stop(145)
        call utl_tmg_stop(132)

      else  ! Not the master, I am a worker and therefore I do actual computations

        ! Main loop over grid points for computing analysis weights
        latLonIndex = 0
        latLonIndexMpiGlobal = 0
        LATLON_LOOP: do

          MPI_DISTRIBUTION: if (trim(enkfNML%mpiDistribution) == 'MASTERWORKER') then

            ! I am a worker and there exists a master, therefore need to get assignment from him

            call utl_tmg_start(132,'----CommWeights')
            call utl_tmg_start(145,'----CommWeights-signals')

            write(*,*) 'Requesting an assignment'
            call utl_printTime()

            ! Check if I have not yet reached max number of weight calculations allowed
            if ( (latLonIndex+1) <= myNumLatLonCalcMax) then

              ! Post recv to obtain assignment or signal that we are done
              call MPI_IRECV(latLonIndexMpiGlobal, 1, MPI_INTEGER, 0, assignmentTag, &
                             mmpi_comm_grid, requestIdSignal, ierr)

              ! Signal that I am ready for assignment
              readySignal = 'Y'
              call MPI_SEND(readySignal, 1, MPI_CHAR, 0, readyTag, &
                            mmpi_comm_grid, ierr)

              call mpi_wait(requestIdSignal, MPI_STATUS_IGNORE, ierr)

            else ! Reached the maximum calculations allowed, inform master

              ! Post recv to obtain assignment that WILL BE IGNORED
              call MPI_IRECV(latLonIndexMpiGlobal, 1, MPI_INTEGER, 0, assignmentTag, &
                             mmpi_comm_grid, requestIdSignal, ierr)

              ! Signal that I reached my maximum number of assignments
              write(*,*) 'Reached maximum number of assignments, wait until others are finished'
              readySignal = 'N'
              call MPI_SEND(readySignal, 1, MPI_CHAR, 0, readyTag, &
                            mmpi_comm_grid, ierr)

              call mpi_wait(requestIdSignal, MPI_STATUS_IGNORE, ierr)

              ! No assignment, so set latLonIndex to 0 to exit LATLON_LOOP
              latLonIndexMpiGlobal = 0

            end if ! check if max allowed calculations not yet reached

            write(*,*) 'Received assignment: ', latLonIndexMpiGlobal
            call utl_printTime()

            call utl_tmg_stop(145)
            call utl_tmg_stop(132)

            ! Check if we are done
            if (latLonIndexMpiGlobal == 0) then
              write(*,*) 'Received the *finished* signal, after computing weights for this many gridpoints: ', latLonIndex
              exit LATLON_LOOP
            end if

          else if (trim(enkfNML%mpiDistribution) == 'ROUNDROBIN') then
            ! There is no master, so just do the next grid point calculations in my list

            ! Find the next value of latLonIndexMpiGlobal that I am responsible for
            do
              latLonIndexMpiGlobal = latLonIndexMpiGlobal + 1
              if (latLonIndexMpiGlobal > numLatLonMpiGlobal) then
                write(*,*) 'Finished processing this many gridpoints: ', latLonIndex
                exit LATLON_LOOP
              end if
              if (mod(latLonIndexMpiGlobal-1,mmpi_nprocs)==mmpi_myid) then
                exit
              end if
            end do

          else ! mpiDistribution is not MASTERWORKER and not ROUNDROBIN

            write(*,*) 'mpiDistribution = ', trim(enkfNML%mpiDistribution)
            call utl_abort('enkf_LETKFanalyses: unknown value of mpiDistribution')

          end if MPI_DISTRIBUTION

          ! Increment main loop index
          latLonIndex = latLonIndex + 1

          ! Check if latLonIndex is larger than expected for allocations
          if (latLonIndex > myNumLatLonCalcMax) then
            write(*,*) 'enkf_LETKFanalyses: We encountered more latLonIndex values on this mpi task than expected!'
            write(*,*) '                    You should probably increase the value of namelist variable: myNumLatLonSendFactor'
            call utl_abort('enkf_LETKFanalyses: increase value of myNumLatLonSendFactor')
          end if

          ! The latitude/longitude index where I will now compute weights
          latIndex = latIndexesSendMpiGlobal(latLonIndexMpiGlobal)
          lonIndex = lonIndexesSendMpiGlobal(latLonIndexMpiGlobal)

          numGridPointWeights = numGridPointWeights + 1

          ! Call main subroutine for computing the LETKF analysis weights for this lev/lat/lon
          call enkf_LETKFcomputeWeights(enkfNML, weightsMeanLatLon(:,:,latLonIndex), &
                                        weightsMembersLatLon(:,:,latLonIndex), &
                                        ensembleAnl, levIndex, latIndex, lonIndex, numGridPointWeights, &
                                        vertLocation_r4, countMaxExceeded, maxCountMaxExceeded, &
                                        ensObs_mpiglobal, ensObsGain_mpiglobal, enkfDFS)

          ! Now post all send instructions (each lat-lon may be sent to multiple tasks)
          call utl_tmg_start(132,'----CommWeights')
          call utl_tmg_start(149,'----CommWeights-isend')
          latIndex = latIndexesSendMpiGlobal(latLonIndexMpiGlobal)
          lonIndex = lonIndexesSendMpiGlobal(latLonIndexMpiGlobal)
          weightsSendCombined(:,1:enkfNML%nEns,latLonIndex) = weightsMembersLatLon(:,:,latLonIndex)
          weightsSendCombined(:,(enkfNML%nEns+1),latLonIndex) = weightsMeanLatLon(:,1,latLonIndex)

          ! Loop over the tasks where I need to send the weights
          do procIndex = 1, numProcsSendMpiGlobal(latLonIndexMpiGlobal)
            sendTag = latLonTagMpiGlobal(lonIndex,latIndex) + (levIndex-1)*maxval(latLonTagMpiGlobal)
            sendTag = 1 + mod(sendTag-1, mmpi_maxTagValue - 10)
            procIndexSend = procIndexesSendMpiGlobal(latLonIndexMpiGlobal, procIndex)

            nsize = nEnsGain * (enkfNML%nEns+1)
            numSend = numSend + 1
            if (numSend > myNumLatLonSendMax) then
              call utl_abort('numSend larger than allowed limit')
            end if
            call mpi_isend( weightsSendCombined(:,:,latLonIndex),  &
                            nsize, mmpi_datyp_real8, procIndexSend-1, sendTag,  &
                            mmpi_comm_grid, requestIdSend(numSend), ierr )
          end do
          call utl_tmg_stop(149)
          call utl_tmg_stop(132)

        end do LATLON_LOOP

      end if MASTER_WORKER

      ! Wait for all previous RECV communications to finish before continuing
      call utl_tmg_start(132,'----CommWeights')
      call utl_tmg_start(147,'----CommWeights-waitRecv')
      if ( numRecv > 0 ) then
        call mpi_waitAll(numRecv, requestIdRecv(1:numRecv), waitStatusesRecv(:,1:numRecv), ierr)
        if (ierr == mpi_err_in_status) call utl_abort('Error code return by mpi_waitAll for IRECV')
      end if
      call utl_tmg_stop(147)

      ! Copy from combined to main separate arrays
      call utl_tmg_start(148,'----CommWeights-irecv')
      do latLonIndex = 1, myNumLatLonRecv
        latIndex = myLatIndexesRecv(latLonIndex)
        lonIndex = myLonIndexesRecv(latLonIndex)
        weightsMembers(:,:,lonIndex,latIndex) = weightsRecvCombined(:,1:enkfNML%nEns,latLonIndex)
        weightsMean(:,1,lonIndex,latIndex)    = weightsRecvCombined(:,(enkfNML%nEns+1),latLonIndex)
        if (any(weightsMembers(:,:,lonIndex,latIndex) < -999.0d0)) then
          write(*,*) 'latLonIndex, latIndex, lonIndex = ', latLonIndex, latIndex, lonIndex
          write(*,*) 'weightsMembers = ', weightsMembers(:,:,lonIndex,latIndex)
          call utl_abort('Invalid weight value received!!!')
        end if
        if (any(weightsMean(:,:,lonIndex,latIndex) < -999.0d0)) then
          write(*,*) 'latLonIndex, latIndex, lonIndex = ', latLonIndex, latIndex, lonIndex
          write(*,*) 'weightsMean = ', weightsMean(:,:,lonIndex,latIndex)
          call utl_abort('Invalid weight value received!!!')
        end if
      end do
      call utl_tmg_stop(148)
      call utl_tmg_stop(132)

      ! Interpolate weights from coarse to full resolution
      call utl_tmg_start(140,'----InterpolateWeights')
      if (wInterpInfo%latLonStep > 1) then
        call enkf_interpWeights(wInterpInfo, weightsMean)
        call enkf_interpWeights(wInterpInfo, weightsMembers)
      end if
      call utl_tmg_stop(140)

      ! Apply the weights to compute the ensemble mean and members
      call utl_tmg_start(143,'----ApplyWeights')
      call enkf_applyEnsWeights(enkfNML, stateVectorMeanInc, stateVectorMeanTrl, stateVectorMeanAnl, &
                                ensembleTrl, ensembleAnl, levIndex,  &
                                weightsMean, weightsMembers, useModulatedEns, &
                                myLonBegHalo, myLatBegHalo)
      call utl_tmg_stop(143)

      ! Wait for SEND communications to finish before continuing to the next level
      call utl_tmg_start(132,'----CommWeights')
      call utl_tmg_start(146,'----CommWeights-waitSend')
      if ( numSend > 0 ) then
        call mpi_waitAll(numSend, requestIdSend(1:numSend), waitStatusesSend(:,1:numSend), ierr)
        if (ierr == mpi_err_in_status) call utl_abort('Error code return by mpi_waitAll for ISEND')
      end if
      call utl_tmg_stop(146)
      call utl_tmg_stop(132)
      ! Should be safe to reuse buffer used for mpi_isend and the call to mpi_waitAll
      weightsSendCombined(:,:,:) = -999.9d0

    end do LEV_LOOP

    if (countMaxExceeded > 0) then
      write(*,*) 'enkf_LETKFanalyses: WARNING: Found more local obs than specified max number at ', &
                 real(100*countMaxExceeded)/real(numGridPointWeights), '% of grid points.'
      write(*,*) '                      Maximum number found was ', maxCountMaxExceeded,  &
                 ' which is greater than specified number ', enkfNML%maxNumLocalObs
      write(*,*) '                      Therefore will keep closest obs only.'
    end if

    call utl_tmg_start(141,'----Barr')
    call rpn_comm_barrier('GRID',ierr)
    call utl_tmg_stop(141)

    call gsv_deallocate(stateVectorMeanInc)
    call gsv_deallocate(stateVectorMeanTrl)

    write(*,*) 'enkf_LETKFanalyses: done'
    call msg_memUsage('enkf_LETKFanalyses')
    write(*,*)

    call utl_tmg_stop(131)

  end subroutine enkf_LETKFanalyses

  !----------------------------------------------------------------------
  ! enkf_LETKFcomputeWeights (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_LETKFcomputeWeights(enkfNML, weightsMeanLatLon, weightsMembersLatLon, &
                                      ensembleAnl, levIndex, latIndex, lonIndex, dfsIndex, &
                                      vertLocation_r4, countMaxExceeded, maxCountMaxExceeded, &
                                      ensObs_mpiglobal, ensObsGain_mpiglobal, enkfDFS)
    !
    !:Purpose: Main routine for computing the LETKF weights for 1 lat/lon/lev. The
    !          first part of the calculation is done in this routine with the rest
    !          done within separate subroutines for each specific algorithm.
    !
    implicit none

    ! Arguments:
    type(struct_enkfNML),     intent(in)    :: enkfNML                   ! Derived type variable with namelist variables
    real(8),                  intent(inout) :: weightsMeanLatLon(:,:)    ! Ens mean weights at one grid point
    real(8),                  intent(inout) :: weightsMembersLatLon(:,:) ! Ens member weights at one grid point
    type(struct_ens),         intent(in)    :: ensembleAnl               ! Analysis ensemble
    integer,                  intent(in)    :: levIndex                  ! Level index
    integer,                  intent(in)    :: latIndex                  ! Latitude index
    integer,                  intent(in)    :: lonIndex                  ! Longitude index
    integer,                  intent(in)    :: dfsIndex                  ! Grid index
    real(4),                  intent(in)    :: vertLocation_r4(:,:,:)    ! Vertical location for each grid point location
    integer,                  intent(inout) :: countMaxExceeded          ! Number of gridpts maxNumLocalObs exceeded
    integer,                  intent(inout) :: maxCountMaxExceeded       ! Maximum number of local obs found
    type(struct_eob), target, intent(in)    :: ensObs_mpiglobal          ! Ens observations for original ensemble
    type(struct_eob),         intent(in)    :: ensObsGain_mpiglobal      ! Ens observations for computing gain matrix
    type(struct_enkfDFS),     intent(inout) :: enkfDFS

    ! Locals:
    type(struct_hco), pointer :: hco_ens
    integer :: nEnsGain
    integer :: numLocalObs, numLocalObsFound, localObsIndex
    integer :: bodyIndex, memberIndex, subEnsIndex
    integer :: nEnsIndependentPerSubEns, nEnsPerSubEns, nEnsPerSubEns_mod
    integer :: eigenVectorColumnIndex, memberIndexInModEns
    logical :: useModulatedEns
    real(8) :: hLoc, anlLat, anlLon, anlVertLocation
    integer, allocatable,         save :: localBodyIndices(:)
    integer, allocatable,         save :: memberIndexSubEns(:,:), memberIndexSubEns_mod(:,:)
    integer, allocatable,         save :: memberIndexSubEnsComp(:,:)
    real(8), allocatable,         save :: locFun(:)
    real(8), pointer,             save :: YbTinvRYb_mean(:,:), YbTinvR_mean(:,:)
    real(8), allocatable, target, save :: YbTinvRYb_pert(:,:)
    real(8), allocatable, target, save :: YbTinvR_pert(:,:)
    real(8), allocatable,         save :: YbTinvRYb_mod(:,:)
    logical, save :: firstCall = .true.

    hco_ens => ens_getHco(ensembleAnl)

    useModulatedEns = (enkfNML%numRetainedEigen > 0)
    if ( useModulatedEns ) then
      nEnsGain   = enkfNML%nEns * enkfNML%numRetainedEigen
    else
      nEnsGain   = enkfNML%nEns
    end if

    if (firstCall) then
      allocate(YbTinvRYb_pert(nEnsGain,nEnsGain))
    end if

    ! Quantities needed for cross validation (CVLETKF, CVLETKF-PERTOBS, CVLETKF-ME)
    if ( trim(enkfNML%algorithm) == 'CVLETKF' .or. trim(enkfNML%algorithm) == 'CVLETKF-PERTOBS' .or. &
         trim(enkfNML%algorithm) == 'CVLETKF-ME' ) then
      nEnsPerSubEns = enkfNML%nEns / enkfNML%numSubEns
      if ( (nEnsPerSubEns * enkfNML%numSubEns) /= enkfNML%nEns ) then
        call utl_abort('enkf_LETKFanalyses: ensemble size not divisible by numSubEnsembles')
      end if
      if (enkfNML%numSubEns <= 1) then
        call utl_abort('enkf_LETKFanalyses: for CVLETKF(-PERTOBS)(-ME) algorithm, numSubEns must be greater than 1')
      end if
      if (useModulatedEns) then
        nEnsPerSubEns_mod = nEnsPerSubEns * enkfNML%numRetainedEigen
        nEnsIndependentPerSubEns = nEnsGain - nEnsPerSubEns_mod
      else
        nEnsIndependentPerSubEns = enkfNML%nEns - nEnsPerSubEns
      end if

      ! Define the subensembles: memberIndexSubEns, memberIndexSubEns_mod, memberIndexSubEnsComp
      if (firstCall) then
        allocate(memberIndexSubEns(nEnsPerSubEns,enkfNML%numSubEns))
        allocate(memberIndexSubEnsComp(nEnsIndependentPerSubEns,enkfNML%numSubEns))
        if ( useModulatedEns ) then
          allocate(memberIndexSubEns_mod(nEnsPerSubEns_mod,enkfNML%numSubEns))
        else
          allocate(memberIndexSubEns_mod(1,1))
        end if

        call enkf_defineSubEnsembles(enkfNML, memberIndexSubEns, memberIndexSubEns_mod, &
                                     memberIndexSubEnsComp, &
                                     nEnsPerSubEns, nEnsPerSubEns_mod, &
                                     useModulatedEns)

        if (mmpi_myid == 0) then
          write(*,*) 'nEns, numSubEns, nEnsPerSubEns, nEnsIndependentPerSubEns = ',  &
                      enkfNML%nEns, enkfNML%numSubEns, nEnsPerSubEns, nEnsIndependentPerSubEns
          do subEnsIndex = 1, enkfNML%numSubEns
            write(*,*) 'memberIndexSubEns = '
            write(*,*) memberIndexSubEns(:,subEnsIndex)
            if ( useModulatedEns ) then
              write(*,*) 'memberIndexSubEns_mod = '
              write(*,*) memberIndexSubEns_mod(:,subEnsIndex)
            end if
            write(*,*) 'memberIndexSubEnsComp = '
            write(*,*) memberIndexSubEnsComp(:,subEnsIndex)
          end do
        end if

      end if ! firstCall

    end if ! if CVLETKF(-PERTOBS)(-ME) algorithm

    ! Allocate arrays only on first call
    if (firstCall) then
      allocate(localBodyIndices(enkfNML%maxNumLocalObs))
      allocate(locFun(enkfNML%maxNumLocalObs))
      allocate(YbTinvR_pert(nEnsGain,enkfNML%maxNumLocalObs))
      if ( trim(enkfNML%algorithm) == 'CVLETKF-ME' .or. &
           trim(enkfNML%algorithm) == 'LETKF-Gain-ME' ) then
        allocate(YbTinvRYb_mod(nEnsGain,enkfNML%nEns))
      end if

      ! Only when observation is "simulated" separate quantities for ens mean needed
      if (eob_simObsAssim) then
        allocate(YbTinvR_mean(nEnsGain,enkfNML%maxNumLocalObs))
        allocate(YbTinvRYb_mean(nEnsGain,nEnsGain))
      else
        YbTinvR_mean => YbTinvR_pert
        YbTinvRYb_mean => YbTinvRYb_pert
      end if
    end if

    ! The lat-lon of the grid point for which we are computing the weights
    anlLat = hco_ens%lat2d_4(lonIndex,latIndex)
    anlLon = hco_ens%lon2d_4(lonIndex,latIndex)
    ! if there is vertical localization
    if (enkfNML%vLocalize > 0.0d0) then
      anlVertLocation = real(vertLocation_r4(lonIndex,latIndex,levIndex),8)
    end if

    ! Get list of nearby observations and localization functions to gridpoint.
    ! With modulated-ensembles, we get observations in entire column.
    call utl_tmg_start(133,'----GetLocalBodyIndices')
    if ( useModulatedEns ) anlVertLocation = MPC_missingValue_R8

    if (enkfDFS%allocated) then
      enkfDFS%lat(dfsIndex) = anlLat
      enkfDFS%lon(dfsIndex) = anlLon
      enkfDFS%lnp(dfsIndex) = anlVertLocation
    end if

    ! Find horizontal localization value for this vertical level
    call enkf_getLocalizationRadius(enkfNML%hLocalize, enkfNML%hLocalizePressure, anlVertLocation, enkfNML%hLinearLoc, hLoc)

    numLocalObs = eob_getLocalBodyIndices(ensObs_mpiglobal, localBodyIndices,     &
                                          locFun, anlLat, anlLon, anlVertLocation,  &
                                          hloc, enkfNML%vLocalize, &
                                          numLocalObsFound, enkfNML%maxNumLocalObsPerType, &
                                          enkfNML%localSelectionOutput, enkfNML%localObsSorting)
    if (numLocalObsFound > enkfNML%maxNumLocalObs) then
      countMaxExceeded = countMaxExceeded + 1
      maxCountMaxExceeded = max(maxCountMaxExceeded, numLocalObsFound)
    end if
    call utl_tmg_stop(133)

    call utl_tmg_start(134,'----CalculateWeights')

    ! Compute YbTinvR (assumes R is diagonal)
    do localObsIndex = 1, numLocalObs
      bodyIndex = localBodyIndices(localObsIndex)
      if (enkfDFS%allocated) then
        enkfDFS%locFun(dfsIndex,localObsIndex) = locFun(localObsIndex)
        enkfDFS%bodyIndex(dfsIndex,localObsIndex) = bodyIndex
      end if

      do memberIndex = 1, nEnsGain
        YbTinvR_pert(memberIndex,localObsIndex) =  &
             ensObsGain_mpiglobal%Yb_r4(memberIndex, bodyIndex) * &
             locFun(localObsIndex) * ensObsGain_mpiglobal%obsErrInv(bodyIndex)
      end do
      if (eob_simObsAssim) then
        do memberIndex = 1, nEnsGain
          ! Compute YbTinvR for the ensemble mean update for simulated observations
          YbTinvR_mean(memberIndex,localObsIndex) =  &
               ensObsGain_mpiglobal%Yb_r4(memberIndex, bodyIndex) * &
               locFun(localObsIndex) * ensObsGain_mpiglobal%obsErrInv_sim(bodyIndex)
        end do
      end if

    end do ! localObsIndex

    ! Compute covariance matrix in ensemble space: YbTinvRYb (this is expensive!)
    call utl_tmg_start(136,'------CalcYbTinvRYb')

    ! Here is the actual calculation
    call utl_tmg_start(138,'--------YbTinvRYb1')
    call enkf_calcYbTinvRYb(nEnsGain, nEnsGain, enkfNML%maxNumLocalObs, numLocalObs, &
                            YbTinvRYb_pert, YbTinvR_pert, &
                            ensObsGain_mpiglobal, localBodyIndices, &
                            YbTinvRYb_mean, YbTinvR_mean)
    call utl_tmg_stop(138)

    ! When using modulated ensemble, also compute YbTinvRYb for perturbation update
    if ( trim(enkfNML%algorithm) == 'CVLETKF-ME' .or. &
         trim(enkfNML%algorithm) == 'LETKF-Gain-ME' ) then

      call utl_tmg_start(139,'--------YbTinvRYb2')
      call enkf_calcYbTinvRYb(nEnsGain, enkfNML%nEns, enkfNML%maxNumLocalObs, numLocalObs, &
                              YbTinvRYb_mod, YbTinvR_pert, &
                              ensObs_mpiglobal, localBodyIndices)
      call utl_tmg_stop(139)

    end if !CVLETKF-ME or LETKF-GAIN-ME

    call utl_tmg_stop(136)

    ! Rest of computation of local weights at this grid point: separate routine for each algorithm
    localObsExist: if (numLocalObs > 0) then

      if (trim(enkfNML%algorithm) == 'LETKF') then

        call enkf_algorithmLETKF(enkfNML, weightsMeanLatLon, weightsMembersLatLon, &
                                 YbTinvRYb_mean, YbTinvRYb_pert, YbTinvR_mean, &
                                 numLocalObs, localBodyIndices, &
                                 ensObs_mpiglobal)

      else if (trim(enkfNML%algorithm) == 'LETKF-Gain') then

        call enkf_algorithmLETKFgain(enkfNML, weightsMeanLatLon, weightsMembersLatLon, &
                                     YbTinvRYb_mean, YbTinvRYb_pert, YbTinvR_mean, &
                                     numLocalObs, localBodyIndices, &
                                     ensObs_mpiglobal)

      else if (trim(enkfNML%algorithm) == 'LETKF-Gain-ME') then

        call enkf_algorithmLETKFgainME(enkfNML, weightsMeanLatLon, weightsMembersLatLon, &
                                       YbTinvRYb_mean, YbTinvRYb_pert, YbTinvRYb_mod, YbTinvR_mean, &
                                       nEnsGain, numLocalObs, localBodyIndices, &
                                       ensObs_mpiglobal)

      else if (trim(enkfNML%algorithm) == 'CVLETKF') then

        call enkf_algorithmCVLETKF(enkfNML, weightsMeanLatLon, weightsMembersLatLon, &
                                   YbTinvRYb_mean, YbTinvRYb_pert, YbTinvR_mean, &
                                   nEnsPerSubEns, nEnsIndependentPerSubEns, &
                                   memberIndexSubEns, memberIndexSubEnsComp, &
                                   numLocalObs, localBodyIndices, ensObs_mpiglobal)

      else if (trim(enkfNML%algorithm) == 'CVLETKF-ME') then

        call enkf_algorithmCVLETKFME(enkfNML, weightsMeanLatLon, weightsMembersLatLon, &
                                     YbTinvRYb_mean, YbTinvRYb_pert, YbTinvRYb_mod, YbTinvR_mean, &
                                     nEnsGain, nEnsPerSubEns, nEnsIndependentPerSubEns, &
                                     memberIndexSubEns, memberIndexSubEnsComp, &
                                     numLocalObs, localBodyIndices, ensObs_mpiglobal)

      else if (trim(enkfNML%algorithm) == 'CVLETKF-PERTOBS') then

        call enkf_algorithmCVLETKFPO(enkfNML, weightsMeanLatLon, weightsMembersLatLon, &
                                     YbTinvRYb_mean, YbTinvRYb_pert, YbTinvR_mean, YbTinvR_pert, &
                                     nEnsPerSubEns, nEnsIndependentPerSubEns, &
                                     memberIndexSubEns, memberIndexSubEnsComp, &
                                     numLocalObs, localBodyIndices, ensObs_mpiglobal)

      else

        call utl_abort('UNKNOWN LETKF ALGORITHM')

      end if

    else ! if numLocalObs > 0

      ! no obs near this grid point, mean weights zero, member weights identity
      weightsMeanLatLon(:,1) = 0.0d0
      weightsMembersLatLon(:,:) = 0.0d0
      do memberIndex = 1, enkfNML%nEns
        if ( useModulatedEns ) then
          do eigenVectorColumnIndex = 1, enkfNML%numRetainedEigen
            memberIndexInModEns = (eigenVectorColumnIndex - 1) * enkfNML%nEns + memberIndex
            weightsMembersLatLon(memberIndexInModEns,memberIndex) = 1.0d0
          end do
        else
          weightsMembersLatLon(memberIndex,memberIndex) = 1.0d0
        end if
      end do

    end if localObsExist

    call utl_tmg_stop(134)

    firstCall = .false.

  end subroutine enkf_LETKFcomputeWeights

  !----------------------------------------------------------------------
  ! enkf_algorithmLETKF (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_algorithmLETKF(enkfNML, weightsMeanLatLon, weightsMembersLatLon, &
                                 YbTinvRYb_mean, YbTinvRYb_pert, YbTinvR_mean, &
                                 numLocalObs, localBodyIndices, ensObs_mpiglobal)
    !
    !:Purpose: Weight calculation for standard LETKF algorithm.
    !
    implicit none

    ! Arguments:
    type(struct_enkfNML), intent(in)  :: enkfNML                   ! Derived type variable with namelist variables
    real(8),              intent(out) :: weightsMeanLatLon(:,:)    ! Ens mean weights at one grid point
    real(8),              intent(out) :: weightsMembersLatLon(:,:) ! Ens member weights at one grid point
    real(8),              intent(in)  :: YbTinvRYb_mean(:,:)       ! Cov matrix in ensemble space for ens mean calculation
    real(8),              intent(in)  :: YbTinvRYb_pert(:,:)       ! Cov matrix in ensemble space for ens member calculation
    real(8),              intent(in)  :: YbTinvR_mean(:,:)         ! Product of Yb and inv(R) needed for LETKF calculation
    integer,              intent(in)  :: numLocalObs               ! Number of local observations for computing weights
    integer,              intent(in)  :: localBodyIndices(:)       ! List of body indices of local observations
    type(struct_eob),     intent(in)  :: ensObs_mpiglobal          ! Ens observations for original ensemble

    ! Locals:
    integer :: memberIndex, memberIndex1, memberIndex2, bodyIndex, localObsIndex
    real(8) :: PaInv(enkfNML%nEns,enkfNML%nEns)
    real(8) :: Pa_pert(enkfNML%nEns,enkfNML%nEns), Pa_mean(enkfNML%nEns,enkfNML%nEns)
    real(8) :: PaSqrt_pert(enkfNML%nEns,enkfNML%nEns)
    real(8) :: weightsTemp(enkfNML%nEns)

    ! Add second term of PaInv
    PaInv(:,:) = YbTinvRYb_pert(:,:)
    do memberIndex = 1, enkfNML%nEns
      PaInv(memberIndex,memberIndex) = PaInv(memberIndex,memberIndex) + real(enkfNML%nEns - 1,8)
    end do
    ! Compute Pa and sqrt(Pa) matrices from PaInv
    Pa_pert(:,:) = PaInv(:,:)
    call utl_tmg_start(135,'------EigenDecomp')
    call utl_matInverse(Pa_pert, enkfNML%nEns, inverseSqrt_opt=PaSqrt_pert)
    call utl_tmg_stop(135)

    if (eob_simObsAssim) then
      PaInv(:,:) = YbTinvRYb_mean(:,:)
      do memberIndex = 1, enkfNML%nEns
        PaInv(memberIndex,memberIndex) = PaInv(memberIndex,memberIndex) + real(enkfNML%nEns - 1,8)
      end do
      Pa_mean(:,:) = PaInv(:,:)
      call utl_tmg_start(135,'------EigenDecomp')
      call utl_matInverse(Pa_mean, enkfNML%nEns)
      call utl_tmg_stop(135)
    else
      Pa_mean(:,:) = Pa_pert(:,:)
    end if

    ! Compute ensemble mean local weights as Pa * YbTinvR * (obs - meanYb)
    weightsTemp(:) = 0.0d0
    do localObsIndex = 1, numLocalObs
      bodyIndex = localBodyIndices(localObsIndex)
      do memberIndex = 1, enkfNML%nEns
        weightsTemp(memberIndex) = weightsTemp(memberIndex) +   &
                                   YbTinvR_mean(memberIndex,localObsIndex) *  &
                                   ( ensObs_mpiglobal%obsValue(bodyIndex) - &
                                     ensObs_mpiglobal%meanYb(bodyIndex) )
      end do
    end do

    weightsMeanLatLon(:,1) = 0.0d0
    do memberIndex2 = 1, enkfNML%nEns
      do memberIndex1 = 1, enkfNML%nEns
        weightsMeanLatLon(memberIndex1,1) =  &
             weightsMeanLatLon(memberIndex1,1) +  &
             Pa_mean(memberIndex1,memberIndex2)*weightsTemp(memberIndex2)
      end do
    end do

    ! Compute ensemble perturbation weights: [(Nens-1)^1/2*PaSqrt]
    weightsMembersLatLon(:,:) = sqrt(real(enkfNML%nEns - 1,8)) * PaSqrt_pert(:,:)

  end subroutine enkf_algorithmLETKF

  !----------------------------------------------------------------------
  ! enkf_algorithmLETKFgain (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_algorithmLETKFgain(enkfNML, weightsMeanLatLon, weightsMembersLatLon, &
                                     YbTinvRYb_mean, YbTinvRYb_pert, YbTinvR_mean, &
                                     numLocalObs, localBodyIndices, ensObs_mpiglobal)
    !
    !:Purpose: Weight calculation for gain-form LETKF algorithm.
    !
    implicit none

    ! Arguments:
    type(struct_enkfNML), intent(in)  :: enkfNML                   ! Derived type variable with namelist variables
    real(8),              intent(out) :: weightsMeanLatLon(:,:)    ! Ens mean weights at one grid point
    real(8),              intent(out) :: weightsMembersLatLon(:,:) ! Ens member weights at one grid point
    real(8),              intent(in)  :: YbTinvRYb_mean(:,:)       ! Cov matrix in ensemble space for ens mean calculation
    real(8),              intent(in)  :: YbTinvRYb_pert(:,:)       ! Cov matrix in ensemble space for ens member calculation
    real(8),              intent(in)  :: YbTinvR_mean(:,:)         ! Product of Yb and inv(R) needed for LETKF calculation
    integer,              intent(in)  :: numLocalObs               ! Number of local observations for computing weights
    integer,              intent(in)  :: localBodyIndices(:)       ! List of body indices of local observations
    type(struct_eob),     intent(in)  :: ensObs_mpiglobal          ! Ens observations for original ensemble

    ! Locals:
    integer :: memberIndex, memberIndex1, memberIndex2, bodyIndex, localObsIndex
    integer :: matrixRank
    real(8) :: tolerance
    real(8) :: eigenValues_mean(enkfNML%nEns), eigenVectors_mean(enkfNML%nEns,enkfNML%nEns)
    real(8) :: eigenValues_pert(enkfNML%nEns), eigenVectors_pert(enkfNML%nEns,enkfNML%nEns)
    real(8) :: weightsTemp(enkfNML%nEns), weightsTemp2(enkfNML%nEns)

    ! Compute eigenValues/Vectors of Yb^T R^-1 Yb = E * Lambda * E^T
    call utl_tmg_start(135,'------EigenDecomp')
    tolerance = 1.0D-50
    call utl_eigenDecomp(YbTinvRYb_pert, eigenValues_pert, eigenVectors_pert, tolerance, matrixRank)
    if (eob_simObsAssim) then
      call utl_eigenDecomp(YbTinvRYb_mean, eigenValues_mean, eigenVectors_mean, tolerance, matrixRank)
    else
      eigenValues_mean(:)    = eigenValues_pert(:)
      eigenVectors_mean(:,:) = eigenVectors_pert(:,:)
    end if
    call utl_tmg_stop(135)

    if (enkfNML%localSelectionOutput /= 0) call enkf_writeEdim(eigenValues_pert, enkfNML%nEns)

    ! Compute ensemble mean local weights as E * (Lambda + (Nens-1)*I)^-1 * E^T * YbTinvR * (obs - meanYb)
    weightsTemp(:) = 0.0d0
    do localObsIndex = 1, numLocalObs
      bodyIndex = localBodyIndices(localObsIndex)
      do memberIndex = 1, enkfNML%nEns
        weightsTemp(memberIndex) = weightsTemp(memberIndex) +   &
                                   YbTinvR_mean(memberIndex,localObsIndex) *  &
                                   ( ensObs_mpiglobal%obsValue(bodyIndex) - &
                                     ensObs_mpiglobal%meanYb(bodyIndex) )
      end do
    end do
    weightsTemp2(:) = 0.0d0
    do memberIndex2 = 1, matrixRank
      do memberIndex1 = 1, enkfNML%nEns
        weightsTemp2(memberIndex2) = weightsTemp2(memberIndex2) +   &
                                     eigenVectors_mean(memberIndex1,memberIndex2) *  &
                                     weightsTemp(memberIndex1)
      end do
    end do
    do memberIndex = 1, matrixRank
      weightsTemp2(memberIndex) = weightsTemp2(memberIndex) *  &
                                  1.0D0/(eigenValues_mean(memberIndex) + real(enkfNML%nEns - 1,8))
    end do
    weightsMeanLatLon(:,1) = 0.0d0
    do memberIndex2 = 1, matrixRank
      do memberIndex1 = 1, enkfNML%nEns
        weightsMeanLatLon(memberIndex1,1) =  &
             weightsMeanLatLon(memberIndex1,1) +   &
             eigenVectors_mean(memberIndex1,memberIndex2) *  &
             weightsTemp2(memberIndex2)
      end do
    end do

    ! Compute ensemble perturbation weights:
    ! Wa = [ - (Nens-1)^1/2 * E *
    !        {(Nens-1)^-1/2*I - (Lambda + (Nens-1)*I)^-1/2} * Lambda^-1 *
    !        E^T * YbTinvRYb ]
    ! Loop over members within the current sub-ensemble being updated
    do memberIndex = 1, enkfNML%nEns

      ! E^T * YbTinvRYb
      weightsTemp(:) = 0.0d0
      do memberIndex2 = 1, matrixRank
        do memberIndex1 = 1, enkfNML%nEns
          weightsTemp(memberIndex2) = weightsTemp(memberIndex2) +  &
                                      eigenVectors_pert(memberIndex1,memberIndex2) *  &
                                      YbTinvRYb_pert(memberIndex1,memberIndex)
        end do
      end do

      ! {(Nens-1)^-1/2*I - (Lambda + (Nens-1)*I)^-1/2} Lambda^-1 * previous_result

      do memberIndex1 = 1, matrixRank
        weightsTemp(memberIndex1) = weightsTemp(memberIndex1) *  &
                                    ( 1.0D0/sqrt(real(enkfNML%nEns - 1,8)) -   &
                                      1.0D0/sqrt(eigenValues_pert(memberIndex1) +  &
                                                 real(enkfNML%nEns - 1,8)) )
        weightsTemp(memberIndex1) = weightsTemp(memberIndex1) /  &
                                    eigenValues_pert(memberIndex1)
      end do

      ! E * previous_result
      weightsMembersLatLon(:,memberIndex) = 0.0d0
      do memberIndex2 = 1, matrixRank
        do memberIndex1 = 1, enkfNML%nEns
          weightsMembersLatLon(memberIndex1,memberIndex) =   &
               weightsMembersLatLon(memberIndex1,memberIndex) +   &
               eigenVectors_pert(memberIndex1,memberIndex2) *  &
               weightsTemp(memberIndex2)
        end do
      end do

      ! -1 * (Nens-1)^1/2 * previous_result
      weightsMembersLatLon(:,memberIndex) =  &
           -1.0D0 * sqrt(real(enkfNML%nEns - 1,8)) *  &
           weightsMembersLatLon(:,memberIndex)

      ! I + previous_result
      weightsMembersLatLon(memberIndex,memberIndex) =  &
           1.0D0 + weightsMembersLatLon(memberIndex,memberIndex)

    end do

    ! Remove the weights mean computed over the columns
    do memberIndex = 1, enkfNML%nEns
      weightsMembersLatLon(memberIndex,:) =  &
           weightsMembersLatLon(memberIndex,:) - &
           sum(weightsMembersLatLon(memberIndex,:))/real(enkfNML%nEns,8)
    end do

  end subroutine enkf_algorithmLETKFgain

  !----------------------------------------------------------------------
  ! enkf_algorithmLETKFgainME (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_algorithmLETKFgainME(enkfNML, weightsMeanLatLon, weightsMembersLatLon, &
                                       YbTinvRYb_mean, YbTinvRYb_pert, YbTinvRYb_mod, YbTinvR_mean, &
                                       nEnsGain, numLocalObs, localBodyIndices, ensObs_mpiglobal)
    !
    !:Purpose: Weight calculation for gain-form LETKF using modulated ensemble algorithm.
    !
    implicit none

    ! Arguments:
    type(struct_enkfNML), intent(in)  :: enkfNML                   ! Derived type variable with namelist variables
    real(8),              intent(out) :: weightsMeanLatLon(:,:)    ! Ens mean weights at one grid point
    real(8),              intent(out) :: weightsMembersLatLon(:,:) ! Ens member weights at one grid point
    real(8),              intent(in)  :: YbTinvRYb_mean(:,:)       ! Cov matrix in ensemble space for ens mean calculation
    real(8),              intent(in)  :: YbTinvRYb_pert(:,:)       ! Cov matrix in ensemble space for ens member calculation
    real(8),              intent(in)  :: YbTinvRYb_mod(:,:)        ! Cov matrix in ensemble space for modulated ensemble
    real(8),              intent(in)  :: YbTinvR_mean(:,:)         ! Product of Yb and inv(R) needed for LETKF calculation
    integer,              intent(in)  :: nEnsGain                  ! Number of members used to compute the gain matrix
    integer,              intent(in)  :: numLocalObs               ! Number of local observations for computing weights
    integer,              intent(in)  :: localBodyIndices(:)       ! List of body indices of local observations
    type(struct_eob),     intent(in)  :: ensObs_mpiglobal          ! Ens observations for original ensemble

    ! Locals:
    integer :: memberIndex, memberIndex1, memberIndex2, bodyIndex, localObsIndex
    integer :: matrixRank
    real(8) :: tolerance
    real(8) :: eigenValues_mean(nEnsGain), eigenVectors_mean(nEnsGain,nEnsGain)
    real(8) :: eigenValues_pert(nEnsGain), eigenVectors_pert(nEnsGain,nEnsGain)
    real(8) :: weightsTemp(nEnsGain), weightsTemp2(nEnsGain)

    ! Compute eigenValues/Vectors of Yb^T R^-1 Yb = E * Lambda * E^T
    call utl_tmg_start(135,'------EigenDecomp')
    tolerance = 1.0D-50
    call utl_eigenDecomp(YbTinvRYb_pert, eigenValues_pert, eigenVectors_pert, tolerance, matrixRank)
    if (eob_simObsAssim) then
      call utl_eigenDecomp(YbTinvRYb_mean, eigenValues_mean, eigenVectors_mean, tolerance, matrixRank)
    else
      eigenValues_mean(:)    = eigenValues_pert(:)
      eigenVectors_mean(:,:) = eigenVectors_pert(:,:)
    end if
    call utl_tmg_stop(135)

    if (enkfNML%localSelectionOutput /= 0) call enkf_writeEdim(eigenValues_pert, enkfNML%nEns)

    ! Compute ensemble mean local weights as E * (Lambda + (Nens-1)*I)^-1 * E^T * YbTinvR * (obs - meanYb)
    weightsTemp(:) = 0.0d0
    do localObsIndex = 1, numLocalObs
      bodyIndex = localBodyIndices(localObsIndex)
      do memberIndex = 1, nEnsGain
        weightsTemp(memberIndex) = weightsTemp(memberIndex) +   &
                                   YbTinvR_mean(memberIndex,localObsIndex) *  &
                                   (ensObs_mpiglobal%obsValue(bodyIndex) - &
                                    ensObs_mpiglobal%meanYb(bodyIndex))
      end do
    end do
    weightsTemp2(:) = 0.0d0
    do memberIndex2 = 1, matrixRank
      do memberIndex1 = 1, nEnsGain
        weightsTemp2(memberIndex2) = weightsTemp2(memberIndex2) +   &
                                     eigenVectors_mean(memberIndex1,memberIndex2) *  &
                                     weightsTemp(memberIndex1)
      end do
    end do
    do memberIndex = 1, matrixRank
      weightsTemp2(memberIndex) = weightsTemp2(memberIndex) *  &
                                  1.0D0/(eigenValues_mean(memberIndex) + real(nEnsGain - 1,8))
    end do
    weightsMeanLatLon(:,1) = 0.0d0
    do memberIndex2 = 1, matrixRank
      do memberIndex1 = 1, nEnsGain
        weightsMeanLatLon(memberIndex1,1) =  &
                 weightsMeanLatLon(memberIndex1,1) +   &
                 eigenVectors_mean(memberIndex1,memberIndex2) *  &
                 weightsTemp2(memberIndex2)
      end do
    end do

    ! Compute ensemble perturbation weights:
    ! Wa = [ - (Nens-1)^1/2 * E *
    !        {(Nens-1)^-1/2*I - (Lambda + (Nens-1)*I)^-1/2} * Lambda^-1 *
    !        E^T * YbTinvRYb_mod ]
    ! Loop over members within the current sub-ensemble being updated
    do memberIndex = 1, enkfNML%nEns

      ! E^T * YbTinvRYb_mod
      weightsTemp(:) = 0.0d0
      do memberIndex2 = 1, matrixRank
        do memberIndex1 = 1, nEnsGain
          weightsTemp(memberIndex2) = weightsTemp(memberIndex2) +  &
                                      eigenVectors_pert(memberIndex1,memberIndex2) *  &
                                      YbTinvRYb_mod(memberIndex1,memberIndex)
        end do
      end do

      ! {(Nens-1)^-1/2*I - (Lambda + (Nens-1)*I)^-1/2} Lambda^-1 * previous_result

      do memberIndex1 = 1, matrixRank
        weightsTemp(memberIndex1) = weightsTemp(memberIndex1) *  &
                                    ( 1.0D0/sqrt(real(nEnsGain - 1,8)) -   &
                                      1.0D0/sqrt(eigenValues_pert(memberIndex1) +  &
                                                 real(nEnsGain - 1,8)) )
        weightsTemp(memberIndex1) = weightsTemp(memberIndex1) /  &
                                    eigenValues_pert(memberIndex1)
      end do

      ! E * previous_result
      weightsMembersLatLon(:,memberIndex) = 0.0d0
      do memberIndex2 = 1, matrixRank
        do memberIndex1 = 1, nEnsGain
          weightsMembersLatLon(memberIndex1,memberIndex) =   &
               weightsMembersLatLon(memberIndex1,memberIndex) +   &
               eigenVectors_pert(memberIndex1,memberIndex2) *  &
               weightsTemp(memberIndex2)
        end do
      end do

      ! -1 * (Nens-1)^1/2 * previous_result
      weightsMembersLatLon(:,memberIndex) =  &
           -1.0D0 * sqrt(real(nEnsGain - 1,8)) *  &
           weightsMembersLatLon(:,memberIndex)

    end do

    ! Remove the weights mean computed over the columns
    do memberIndex = 1, nEnsGain
      weightsMembersLatLon(memberIndex,:) =  &
           weightsMembersLatLon(memberIndex,:) - &
           sum(weightsMembersLatLon(memberIndex,:))/real(enkfNML%nEns,8)
    end do

  end subroutine enkf_algorithmLETKFgainME

  !----------------------------------------------------------------------
  ! enkf_algorithmCVLETKF (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_algorithmCVLETKF(enkfNML, weightsMeanLatLon, weightsMembersLatLon, &
                                   YbTinvRYb_mean, YbTinvRYb_pert, YbTinvR_mean, &
                                   nEnsPerSubEns, nEnsIndependentPerSubEns, &
                                   memberIndexSubEns, memberIndexSubEnsComp, &
                                   numLocalObs, localBodyIndices, ensObs_mpiglobal)
    !
    !:Purpose: Weight calculation for LETKF with cross validation algorithm.
    !
    implicit none

    ! Arguments:
    type(struct_enkfNML), intent(in)  :: enkfNML                    ! Derived type variable with namelist variables
    real(8),              intent(out) :: weightsMeanLatLon(:,:)     ! Ens mean weights at one grid point
    real(8),              intent(out) :: weightsMembersLatLon(:,:)  ! Ens member weights at one grid point
    real(8),              intent(in)  :: YbTinvRYb_mean(:,:)        ! Cov matrix in ensemble space for ens mean calculation
    real(8),              intent(in)  :: YbTinvRYb_pert(:,:)        ! Cov matrix in ensemble space for ens member calculation
    real(8),              intent(in)  :: YbTinvR_mean(:,:)          ! Product of Yb and inv(R) needed for LETKF calculation
    integer,              intent(in)  :: nEnsPerSubEns              ! Number of members per sub-ensemble
    integer,              intent(in)  :: nEnsIndependentPerSubEns   ! Number of members in all other sub-ensembles
    integer,              intent(in)  :: memberIndexSubEns(:,:)     ! Member indexes in each sub-ensemble
    integer,              intent(in)  :: memberIndexSubEnsComp(:,:) ! Member indexes in all other sub-ensembles
    integer,              intent(in)  :: numLocalObs                ! Number of local observations for computing weights
    integer,              intent(in)  :: localBodyIndices(:)        ! List of body indices of local observations
    type(struct_eob),     intent(in)  :: ensObs_mpiglobal           ! Ens observations for original ensemble

    ! Locals:
    integer :: memberIndex, memberIndex1, memberIndex2, bodyIndex, localObsIndex
    integer :: matrixRank, subEnsIndex, memberIndexCV, memberIndexCV1, memberIndexCV2
    real(8) :: tolerance
    real(8) :: eigenValues_mean(enkfNML%nEns), eigenVectors_mean(enkfNML%nEns,enkfNML%nEns)
    real(8) :: weightsTemp(enkfNML%nEns), weightsTemp2(enkfNML%nEns)
    real(8) :: YbTinvRYb_CV(nEnsIndependentPerSubEns,nEnsIndependentPerSubEns)
    real(8) :: eigenValues_CV(nEnsIndependentPerSubEns)
    real(8) :: eigenVectors_CV(nEnsIndependentPerSubEns,nEnsIndependentPerSubEns)

    ! Compute eigenValues/Vectors of Yb^T R^-1 Yb = E * Lambda * E^T
    call utl_tmg_start(135,'------EigenDecomp')
    tolerance = 1.0D-50
    call utl_eigenDecomp(YbTinvRYb_mean, eigenValues_mean, eigenVectors_mean, tolerance, matrixRank)
    call utl_tmg_stop(135)

    if (enkfNML%localSelectionOutput /= 0) call enkf_writeEdim(eigenValues_mean, enkfNML%nEns)

    ! Compute ensemble mean local weights as E * (Lambda + (Nens-1)*I)^-1 * E^T * YbTinvR * (obs - meanYb)
    weightsTemp(:) = 0.0d0
    do localObsIndex = 1, numLocalObs
      bodyIndex = localBodyIndices(localObsIndex)
      do memberIndex = 1, enkfNML%nEns
        weightsTemp(memberIndex) = weightsTemp(memberIndex) +   &
                                   YbTinvR_mean(memberIndex,localObsIndex) *  &
                                   ( ensObs_mpiglobal%obsValue(bodyIndex) - &
                                     ensObs_mpiglobal%meanYb(bodyIndex) )
      end do
    end do
    weightsTemp2(:) = 0.0d0
    do memberIndex2 = 1, matrixRank
      do memberIndex1 = 1, enkfNML%nEns
        weightsTemp2(memberIndex2) = weightsTemp2(memberIndex2) +   &
                                     eigenVectors_mean(memberIndex1,memberIndex2) *  &
                                     weightsTemp(memberIndex1)
      end do
    end do
    do memberIndex = 1, matrixRank
      weightsTemp2(memberIndex) = weightsTemp2(memberIndex) *  &
                                  1.0D0/(eigenValues_mean(memberIndex) + real(enkfNML%nEns - 1,8))
    end do
    weightsMeanLatLon(:,1) = 0.0d0
    do memberIndex2 = 1, matrixRank
      do memberIndex1 = 1, enkfNML%nEns
        weightsMeanLatLon(memberIndex1,1) =  &
             weightsMeanLatLon(memberIndex1,1) +   &
             eigenVectors_mean(memberIndex1,memberIndex2) *  &
             weightsTemp2(memberIndex2)
      end do
    end do

    ! Compute ensemble perturbation weights:
    ! Wa = [ I - (Nens-1)^1/2 * E *
    !        {(Nens-1)^-1/2*I - (Lambda + (Nens-1)*I)^-1/2} * Lambda^-1 *
    !        E^T * YbTinvRYb ]
    ! Loop over sub-ensembles

    !$OMP PARALLEL DO PRIVATE(subEnsIndex, memberIndexCV, memberIndexCV1, memberIndexCV2, &
    !$OMP                     memberIndex, memberIndex1, memberIndex2, weightsTemp, tolerance, &
    !$OMP                     YbTinvRYb_CV, eigenValues_CV, eigenVectors_CV, matrixRank)
    do subEnsIndex = 1, enkfNML%numSubEns

      ! Use complement (independent) ens to get eigenValues/Vectors of Yb^T R^-1 Yb = E*Lambda*E^T
      call utl_tmg_start(135,'------EigenDecomp')
      do memberIndexCV2 = 1, nEnsIndependentPerSubEns
        memberIndex2 = memberIndexSubEnsComp(memberIndexCV2, subEnsIndex)
        do memberIndexCV1 = 1, nEnsIndependentPerSubEns
          memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
          YbTinvRYb_CV(memberIndexCV1,memberIndexCV2) = YbTinvRYb_pert(memberIndex1,memberIndex2)
        end do
      end do
      tolerance = 1.0D-50
      call utl_eigenDecomp(YbTinvRYb_CV, eigenValues_CV, eigenVectors_CV, tolerance, matrixRank)
      call utl_tmg_stop(135)

      ! Loop over members within the current sub-ensemble being updated
      do memberIndexCV = 1, nEnsPerSubEns

        ! This is index of member being updated
        memberIndex = memberIndexSubEns(memberIndexCV, subEnsIndex)

        ! E^T * YbTinvRYb
        weightsTemp(:) = 0.0d0
        do memberIndex2 = 1, matrixRank
          do memberIndexCV1 = 1, nEnsIndependentPerSubEns
            memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
            weightsTemp(memberIndex2) = weightsTemp(memberIndex2) +  &
                                        eigenVectors_CV(memberIndexCV1,memberIndex2) *  &
                                        YbTinvRYb_pert(memberIndex1,memberIndex)
          end do
        end do

        ! {(Nens-1)^-1/2*I - (Lambda + (Nens-1)*I)^-1/2} Lambda^-1 * previous_result

        do memberIndex1 = 1, matrixRank
          weightsTemp(memberIndex1) = weightsTemp(memberIndex1) *  &
                                      ( 1.0D0/sqrt(real(nEnsIndependentPerSubEns - 1,8)) -   &
                                        1.0D0/sqrt(eigenValues_CV(memberIndex1) +  &
                                                   real(nEnsIndependentPerSubEns - 1,8)) )
          weightsTemp(memberIndex1) = weightsTemp(memberIndex1) /  &
                                      eigenValues_CV(memberIndex1)
        end do

        ! E * previous_result
        weightsMembersLatLon(:,memberIndex) = 0.0d0
        do memberIndex2 = 1, matrixRank
          do memberIndexCV1 = 1, nEnsIndependentPerSubEns
            memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
            weightsMembersLatLon(memberIndex1,memberIndex) =   &
                 weightsMembersLatLon(memberIndex1,memberIndex) +   &
                 eigenVectors_CV(memberIndexCV1,memberIndex2) *  &
                 weightsTemp(memberIndex2)
          end do
        end do

        ! -1 * (Nens-1)^1/2 * previous_result
        weightsMembersLatLon(:,memberIndex) =  &
             -1.0D0 * sqrt(real(nEnsIndependentPerSubEns - 1,8)) *  &
             weightsMembersLatLon(:,memberIndex)

        ! I + previous_result
        weightsMembersLatLon(memberIndex,memberIndex) =  &
             1.0D0 + weightsMembersLatLon(memberIndex,memberIndex)

      end do ! memberIndexCV
    end do ! subEnsIndex
    !$OMP END PARALLEL DO

    ! Remove the weights mean computed over the columns
    do memberIndex = 1, enkfNML%nEns
      weightsMembersLatLon(memberIndex,:) =  &
           weightsMembersLatLon(memberIndex,:) - &
           sum(weightsMembersLatLon(memberIndex,:))/real(enkfNML%nEns,8)
    end do

  end subroutine enkf_algorithmCVLETKF

  !----------------------------------------------------------------------
  ! enkf_algorithmCVLETKFME (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_algorithmCVLETKFME(enkfNML, weightsMeanLatLon, weightsMembersLatLon, &
                                     YbTinvRYb_mean, YbTinvRYb_pert, YbTinvRYb_mod, YbTinvR_mean, &
                                     nEnsGain, nEnsPerSubEns, nEnsIndependentPerSubEns, &
                                     memberIndexSubEns, memberIndexSubEnsComp, &
                                     numLocalObs, localBodyIndices, ensObs_mpiglobal)
    !
    !:Purpose: Weight calculation for LETKF with cross validation and
    !          modulated ensemble algorithm.
    !
    implicit none

    ! Arguments:
    type(struct_enkfNML), intent(in)  :: enkfNML                    ! Derived type variable with namelist variables
    real(8),              intent(out) :: weightsMeanLatLon(:,:)     ! Ens mean weights at one grid point
    real(8),              intent(out) :: weightsMembersLatLon(:,:)  ! Ens member weights at one grid point
    real(8),              intent(in)  :: YbTinvRYb_mean(:,:)        ! Cov matrix in ensemble space for ens mean calculation
    real(8),              intent(in)  :: YbTinvRYb_pert(:,:)        ! Cov matrix in ensemble space for ens member calculation
    real(8),              intent(in)  :: YbTinvRYb_mod(:,:)         ! Cov matrix in ensemble space for modulated ensemble
    real(8),              intent(in)  :: YbTinvR_mean(:,:)          ! Product of Yb and inv(R) needed for LETKF calculation
    integer,              intent(in)  :: nEnsGain                   ! Number of members used to compute the gain matrix
    integer,              intent(in)  :: nEnsPerSubEns              ! Number of members per sub-ensemble
    integer,              intent(in)  :: nEnsIndependentPerSubEns   ! Number of members in all other sub-ensembles
    integer,              intent(in)  :: memberIndexSubEns(:,:)     ! Member indexes in each sub-ensemble
    integer,              intent(in)  :: memberIndexSubEnsComp(:,:) ! Member indexes in all other sub-ensembles
    integer,              intent(in)  :: numLocalObs                ! Number of local observations for computing weights
    integer,              intent(in)  :: localBodyIndices(:)        ! List of body indices of local observations
    type(struct_eob),     intent(in)  :: ensObs_mpiglobal           ! Ens observations for original ensemble

    ! Locals:
    integer :: memberIndex, memberIndex1, memberIndex2, bodyIndex, localObsIndex
    integer :: matrixRank, subEnsIndex, memberIndexCV, memberIndexCV1, memberIndexCV2
    real(8) :: tolerance
    real(8) :: eigenValues_mean(nEnsGain), eigenVectors_mean(nEnsGain,nEnsGain)
    real(8) :: weightsTemp(nEnsGain), weightsTemp2(nEnsGain)
    real(8) :: YbTinvRYb_CV(nEnsIndependentPerSubEns,nEnsIndependentPerSubEns)
    real(8) :: eigenValues_CV(nEnsIndependentPerSubEns)
    real(8) :: eigenVectors_CV(nEnsIndependentPerSubEns,nEnsIndependentPerSubEns)

    ! Compute eigenValues/Vectors of Yb^T R^-1 Yb = E * Lambda * E^T
    call utl_tmg_start(135,'------EigenDecomp')
    tolerance = 1.0D-50
    call utl_eigenDecomp(YbTinvRYb_mean, eigenValues_mean, eigenVectors_mean, tolerance, matrixRank)
    call utl_tmg_stop(135)
    !if (matrixRank < (nEns-1)) then
    !  write(*,*) 'YbTinvRYb is rank deficient =', matrixRank, nEns, numLocalObs
    !end if

    if (enkfNML%localSelectionOutput /= 0) call enkf_writeEdim(eigenValues_mean, enkfNML%nEns)

    ! Compute ensemble mean local weights as E * (Lambda + (Nens-1)*I)^-1 * E^T * YbTinvR * (obs - meanYb)
    weightsTemp(:) = 0.0d0
    do localObsIndex = 1, numLocalObs
      bodyIndex = localBodyIndices(localObsIndex)
      do memberIndex = 1, nEnsGain
        weightsTemp(memberIndex) = weightsTemp(memberIndex) +   &
                                   YbTinvR_mean(memberIndex,localObsIndex) *  &
                                   ( ensObs_mpiglobal%obsValue(bodyIndex) - &
                                     ensObs_mpiglobal%meanYb(bodyIndex) )
      end do
    end do
    weightsTemp2(:) = 0.0d0
    do memberIndex2 = 1, matrixRank
      do memberIndex1 = 1, nEnsGain
        weightsTemp2(memberIndex2) = weightsTemp2(memberIndex2) +   &
                                     eigenVectors_mean(memberIndex1,memberIndex2) *  &
                                     weightsTemp(memberIndex1)
      end do
    end do
    do memberIndex = 1, matrixRank
      weightsTemp2(memberIndex) = weightsTemp2(memberIndex) *  &
                                  1.0D0/(eigenValues_mean(memberIndex) + real(nEnsGain - 1,8))
    end do
    weightsMeanLatLon(:,1) = 0.0d0
    do memberIndex2 = 1, matrixRank
      do memberIndex1 = 1, nEnsGain
        weightsMeanLatLon(memberIndex1,1) =  &
             weightsMeanLatLon(memberIndex1,1) +   &
             eigenVectors_mean(memberIndex1,memberIndex2) *  &
             weightsTemp2(memberIndex2)
      end do
    end do

    ! Compute ensemble perturbation weights:
    ! Wa = [ - (Nens-1)^1/2 * E *
    !        {(Nens-1)^-1/2*I - (Lambda + (Nens-1)*I)^-1/2} * Lambda^-1 *
    !        E^T * YbTinvRYb_mod ]
    ! Loop over sub-ensembles
    !$OMP PARALLEL DO PRIVATE(subEnsIndex, memberIndexCV, memberIndexCV1, memberIndexCV2, &
    !$OMP                     memberIndex, memberIndex1, memberIndex2, weightsTemp, tolerance, &
    !$OMP                     YbTinvRYb_CV, eigenValues_CV, eigenVectors_CV, matrixRank)
    do subEnsIndex = 1, enkfNML%numSubEns

      ! Use complement (independent) ens to get eigenValues/Vectors of Yb^T R^-1 Yb = E*Lambda*E^T
      call utl_tmg_start(135,'------EigenDecomp')
      do memberIndexCV2 = 1, nEnsIndependentPerSubEns
        memberIndex2 = memberIndexSubEnsComp(memberIndexCV2, subEnsIndex)
        do memberIndexCV1 = 1, nEnsIndependentPerSubEns
          memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
          YbTinvRYb_CV(memberIndexCV1,memberIndexCV2) = YbTinvRYb_pert(memberIndex1,memberIndex2)
        end do
      end do
      tolerance = 1.0D-50
      call utl_eigenDecomp(YbTinvRYb_CV, eigenValues_CV, eigenVectors_CV, tolerance, matrixRank)
      call utl_tmg_stop(135)

      ! Loop over members within the current sub-ensemble being updated
      do memberIndexCV = 1, nEnsPerSubEns

        ! This is index of member being updated
        memberIndex = memberIndexSubEns(memberIndexCV, subEnsIndex)

        ! E^T * YbTinvRYb
        weightsTemp(:) = 0.0d0
        do memberIndex2 = 1, matrixRank
          do memberIndexCV1 = 1, nEnsIndependentPerSubEns
            memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
            weightsTemp(memberIndex2) = weightsTemp(memberIndex2) +  &
                                        eigenVectors_CV(memberIndexCV1,memberIndex2) *  &
                                        YbTinvRYb_mod(memberIndex1,memberIndex)
          end do
        end do

        ! {(Nens-1)^-1/2*I - (Lambda + (Nens-1)*I)^-1/2} Lambda^-1 * previous_result

        do memberIndex1 = 1, matrixRank
          weightsTemp(memberIndex1) = weightsTemp(memberIndex1) *  &
                                      ( 1.0D0/sqrt(real(nEnsIndependentPerSubEns - 1,8)) -   &
                                        1.0D0/sqrt(eigenValues_CV(memberIndex1) +  &
                                                   real(nEnsIndependentPerSubEns - 1,8)) )
          weightsTemp(memberIndex1) = weightsTemp(memberIndex1) /  &
                                      eigenValues_CV(memberIndex1)
        end do

        ! E * previous_result
        weightsMembersLatLon(:,memberIndex) = 0.0d0
        do memberIndex2 = 1, matrixRank
          do memberIndexCV1 = 1, nEnsIndependentPerSubEns
            memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
            weightsMembersLatLon(memberIndex1,memberIndex) =   &
                 weightsMembersLatLon(memberIndex1,memberIndex) +   &
                 eigenVectors_CV(memberIndexCV1,memberIndex2) *  &
                 weightsTemp(memberIndex2)
          end do
        end do

        ! -1 * (Nens-1)^1/2 * previous_result
        weightsMembersLatLon(:,memberIndex) =  &
             -1.0D0 * sqrt(real(nEnsIndependentPerSubEns - 1,8)) *  &
             weightsMembersLatLon(:,memberIndex)

      end do ! memberIndexCV
    end do ! subEnsIndex
    !$OMP END PARALLEL DO

    ! Remove the weights mean computed over the columns
    do memberIndex = 1, nEnsGain
      weightsMembersLatLon(memberIndex,:) =  &
           weightsMembersLatLon(memberIndex,:) - &
           sum(weightsMembersLatLon(memberIndex,:))/real(enkfNML%nEns,8)
    end do

  end subroutine enkf_algorithmCVLETKFME

  !----------------------------------------------------------------------
  ! enkf_algorithmCVLETKFPO (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_algorithmCVLETKFPO(enkfNML, weightsMeanLatLon, weightsMembersLatLon, &
                                     YbTinvRYb_mean, YbTinvRYb_pert, YbTinvR_mean, YbTinvR_pert, &
                                     nEnsPerSubEns, nEnsIndependentPerSubEns, &
                                     memberIndexSubEns, memberIndexSubEnsComp, &
                                     numLocalObs, localBodyIndices, ensObs_mpiglobal)
    !
    !:Purpose: Weight calculation for LETKF with cross validation with
    !          perturbed observations algorithm.
    !
    implicit none

    ! Arguments:
    type(struct_enkfNML), intent(in)  :: enkfNML                    ! Derived type variable with namelist variables
    real(8),              intent(out) :: weightsMeanLatLon(:,:)     ! Ens mean weights at one grid point
    real(8),              intent(out) :: weightsMembersLatLon(:,:)  ! Ens member weights at one grid point
    real(8),              intent(in)  :: YbTinvRYb_mean(:,:)        ! Cov matrix in ensemble space for ens mean calculation
    real(8),              intent(in)  :: YbTinvRYb_pert(:,:)        ! Cov matrix in ensemble space for ens member calculation
    real(8),              intent(in)  :: YbTinvR_mean(:,:)          ! Product of Yb and inv(R) needed for LETKF calculation
    real(8),              intent(in)  :: YbTinvR_pert(:,:)          ! Product of Yb and inv(R) needed for LETKF calculation
    integer,              intent(in)  :: nEnsPerSubEns              ! Number of members per sub-ensemble
    integer,              intent(in)  :: nEnsIndependentPerSubEns   ! Number of members in all other sub-ensembles
    integer,              intent(in)  :: memberIndexSubEns(:,:)     ! Member indexes in each sub-ensemble
    integer,              intent(in)  :: memberIndexSubEnsComp(:,:) ! Member indexes in all other sub-ensembles
    integer,              intent(in)  :: numLocalObs                ! Number of local observations for computing weights
    integer,              intent(in)  :: localBodyIndices(:)        ! List of body indices of local observations
    type(struct_eob),     intent(in)  :: ensObs_mpiglobal           ! Ens observations for original ensemble

    ! Locals:
    integer :: memberIndex, memberIndex1, memberIndex2, bodyIndex, localObsIndex
    integer :: matrixRank, subEnsIndex, memberIndexCV, memberIndexCV1, memberIndexCV2
    real(8) :: tolerance
    real(8) :: eigenValues_mean(enkfNML%nEns), eigenVectors_mean(enkfNML%nEns,enkfNML%nEns)
    real(8) :: weightsTemp(enkfNML%nEns), weightsTemp2(enkfNML%nEns)
    real(8) :: YbTinvRYb_CV(nEnsIndependentPerSubEns,nEnsIndependentPerSubEns)
    real(8) :: eigenValues_CV(nEnsIndependentPerSubEns)
    real(8) :: eigenVectors_CV(nEnsIndependentPerSubEns,nEnsIndependentPerSubEns)

    ! Compute eigenValues/Vectors of Yb^T R^-1 Yb = E * Lambda * E^T
    call utl_tmg_start(135,'------EigenDecomp')
    tolerance = 1.0D-50
    call utl_eigenDecomp(YbTinvRYb_mean, eigenValues_mean, eigenVectors_mean, tolerance, matrixRank)
    call utl_tmg_stop(135)

    if (enkfNML%localSelectionOutput /= 0) call enkf_writeEdim(eigenValues_mean, enkfNML%nEns)

    ! Compute ensemble mean local weights as E * (Lambda + (Nens-1)*I)^-1 * E^T * YbTinvR * (obs - meanYb)
    weightsTemp(:) = 0.0d0
    do localObsIndex = 1, numLocalObs
      bodyIndex = localBodyIndices(localObsIndex)
      do memberIndex = 1, enkfNML%nEns
        weightsTemp(memberIndex) = weightsTemp(memberIndex) +   &
                                   YbTinvR_mean(memberIndex,localObsIndex) *  &
                                   ( ensObs_mpiglobal%obsValue(bodyIndex) - &
                                     ensObs_mpiglobal%meanYb(bodyIndex) )
      end do
    end do
    weightsTemp2(:) = 0.0d0
    do memberIndex2 = 1, matrixRank
      do memberIndex1 = 1, enkfNML%nEns
        weightsTemp2(memberIndex2) = weightsTemp2(memberIndex2) +   &
                                     eigenVectors_mean(memberIndex1,memberIndex2) *  &
                                     weightsTemp(memberIndex1)
      end do
    end do
    do memberIndex = 1, matrixRank
      weightsTemp2(memberIndex) = weightsTemp2(memberIndex) *  &
                                  1.0D0/(eigenValues_mean(memberIndex) + real(enkfNML%nEns - 1,8))
    end do
    weightsMeanLatLon(:,1) = 0.0d0
    do memberIndex2 = 1, matrixRank
      do memberIndex1 = 1, enkfNML%nEns
        weightsMeanLatLon(memberIndex1,1) =  &
             weightsMeanLatLon(memberIndex1,1) +   &
             eigenVectors_mean(memberIndex1,memberIndex2) *  &
             weightsTemp2(memberIndex2)
      end do
    end do

    ! Compute ensemble perturbation weights using mean increment weights
    ! formula, but with subset of members:
    ! wa_i = I_i + E * (Lambda + (Nens-1)*I)^-1 * E^T * YbTinvR * (obs + randpert_i - Yb_i)
    ! Wa   = wa_i - mean_over_i(wa_i)
    !
    ! Loop over sub-ensembles
    !$OMP PARALLEL DO PRIVATE(subEnsIndex, memberIndexCV, memberIndexCV1, memberIndexCV2, &
    !$OMP                     memberIndex, memberIndex1, memberIndex2, weightsTemp, tolerance, &
    !$OMP                     YbTinvRYb_CV, eigenValues_CV, eigenVectors_CV, matrixRank, &
    !$OMP                     weightsTemp2, localObsIndex, bodyIndex)
    do subEnsIndex = 1, enkfNML%numSubEns

      ! Use complement (independent) ens to get eigenValues/Vectors of Yb^T R^-1 Yb = E*Lambda*E^T
      call utl_tmg_start(135,'------EigenDecomp')
      do memberIndexCV2 = 1, nEnsIndependentPerSubEns
        memberIndex2 = memberIndexSubEnsComp(memberIndexCV2, subEnsIndex)
        do memberIndexCV1 = 1, nEnsIndependentPerSubEns
          memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
          YbTinvRYb_CV(memberIndexCV1,memberIndexCV2) = YbTinvRYb_pert(memberIndex1,memberIndex2)
        end do
      end do
      tolerance = 1.0D-50
      call utl_eigenDecomp(YbTinvRYb_CV, eigenValues_CV, eigenVectors_CV, tolerance, matrixRank)
      call utl_tmg_stop(135)

      ! Loop over members within the current sub-ensemble being updated
      do memberIndexCV = 1, nEnsPerSubEns

        ! This is index of member being updated (i'th member)
        memberIndex = memberIndexSubEns(memberIndexCV, subEnsIndex)

        ! YbTinvRYb * (obsValue + randPert_i - Yb_i)
        weightsTemp(:) = 0.0d0
        do localObsIndex = 1, numLocalObs
          bodyIndex = localBodyIndices(localObsIndex)
          do memberIndexCV1 = 1, nEnsIndependentPerSubEns
            memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
            weightsTemp(memberIndexCV1) =  &
                 weightsTemp(memberIndexCV1) +   &
                 YbTinvR_pert(memberIndex1,localObsIndex) *  &
                 ( ensObs_mpiglobal%obsValue(bodyIndex) +  &
                   ensObs_mpiglobal%randPert_r4(memberIndex,bodyIndex) -  &
                   ( ensObs_mpiglobal%meanYb(bodyIndex) +  &
                     ensObs_mpiglobal%Yb_r4(memberIndex,bodyIndex) ) )
          end do
        end do

        ! E^T * previous_result
        weightsTemp2(:) = 0.0d0
        do memberIndex2 = 1, matrixRank
          do memberIndex1 = 1, nEnsIndependentPerSubEns
            weightsTemp2(memberIndex2) = weightsTemp2(memberIndex2) +   &
                                         eigenVectors_CV(memberIndex1,memberIndex2) *  &
                                         weightsTemp(memberIndex1)
          end do
        end do

        ! [lambda + (N_indep-1)*I]^-1 * previous_result
        do memberIndex1 = 1, matrixRank
          weightsTemp2(memberIndex1) =  &
               weightsTemp2(memberIndex1) *  &
               1.0D0/(eigenValues_CV(memberIndex1) + real(nEnsIndependentPerSubEns - 1,8))
        end do

        ! E * previous_result
        weightsMembersLatLon(:,memberIndex) = 0.0d0
        do memberIndex2 = 1, matrixRank
          do memberIndexCV1 = 1, nEnsIndependentPerSubEns
            memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
            weightsMembersLatLon(memberIndex1,memberIndex) =  &
                 weightsMembersLatLon(memberIndex1,memberIndex) +   &
                 eigenVectors_CV(memberIndexCV1,memberIndex2) *  &
                 weightsTemp2(memberIndex2)
          end do
        end do

        ! I + previous_result
        weightsMembersLatLon(memberIndex,memberIndex) =  &
             1.0D0 + weightsMembersLatLon(memberIndex,memberIndex)

      end do ! memberIndexCV
    end do ! subEnsIndex
    !$OMP END PARALLEL DO

    ! Remove the weights mean computed over the columns
    do memberIndex = 1, enkfNML%nEns
      weightsMembersLatLon(memberIndex,:) =  &
           weightsMembersLatLon(memberIndex,:) - &
           sum(weightsMembersLatLon(memberIndex,:))/real(enkfNML%nEns,8)
    end do

  end subroutine enkf_algorithmCVLETKFPO

  !----------------------------------------------------------------------
  ! enkf_defineSubEnsembles (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_defineSubEnsembles(enkfNML, memberIndexSubEns, memberIndexSubEns_mod, &
                                     memberIndexSubEnsComp, &
                                     nEnsPerSubEns, nEnsPerSubEns_mod, &
                                     useModulatedEns)
    !
    !:Purpose: Compute the list of member indexes for each sub-ensemble. Also computes
    !          a second list of member indexes for each sub-ensemble for the modulated
    !          ensemble (when this is being used).
    !
    implicit none

    ! Arguments:
    type(struct_enkfNML), intent(in)  :: enkfNML                    ! Derived type variable with namelist variables
    integer,              intent(out) :: memberIndexSubEns(:,:)     ! Array of member indexes for each subEns
    integer,              intent(out) :: memberIndexSubEns_mod(:,:) ! Array of modulated member indexes for each subEns
    integer,              intent(out) :: memberIndexSubEnsComp(:,:) ! Array of member indexes complementary to each subEns
    integer,              intent(in)  :: nEnsPerSubEns              ! Number of members per sub-ensemble
    integer,              intent(in)  :: nEnsPerSubEns_mod          ! Number of modulated members per sub-ensemble
    logical,              intent(in)  :: useModulatedEns            ! Switch to control if using modulated ensemble

    ! Locals:
    integer :: subEnsIndex, subEnsIndex2, memberIndex, memberIndex2
    integer :: imode, dateStamp, ierr, timePrint, datePrint, randomSeed
    integer :: eigenVectorColumnIndex, memberIndexInModEns, newDate
    integer, allocatable, save :: randomMemberIndexArray(:)
    logical, save :: firstCall = .true.

    if (.not.enkfNML%randomShuffleSubEns) then
      ! form subensembles with contiguous sequential groups of members
      do subEnsIndex = 1, enkfNML%numSubEns
        do memberIndex = 1, nEnsPerSubEns
          memberIndexSubEns(memberIndex,subEnsIndex) =  &
               (subEnsIndex-1)*nEnsPerSubEns + memberIndex
        end do
      end do
      if ( useModulatedEns ) then
        do subEnsIndex = 1, enkfNML%numSubEns
          memberIndex2 = 0
          do memberIndex = 1, nEnsPerSubEns
            do eigenVectorColumnIndex = 1, enkfNML%numRetainedEigen
              memberIndex2 = memberIndex2 + 1
              memberIndexInModEns = (eigenVectorColumnIndex - 1) * enkfNML%nEns + &
                   memberIndex
              memberIndexSubEns_mod(memberIndex2,subEnsIndex) =  &
                   (subEnsIndex-1)*nEnsPerSubEns + memberIndexInModEns
            end do
          end do
        end do
      end if
    else
      ! compute random seed from the date for randomly forming subensembles
      imode = -3 ! stamp to printable date and time: YYYYMMDD, HHMMSShh
      dateStamp = tim_getDateStamp()
      ierr = newdate(dateStamp, datePrint, timePrint, imode)
      timePrint = timePrint/1000000
      datePrint =  datePrint*100 + timePrint
      ! Remove the century, keeping 2 digits of the year
      randomSeed = datePrint - 100000000*(datePrint/100000000)
      if (firstCall) allocate(randomMemberIndexArray(enkfNML%nEns))
      do memberIndex = 1, enkfNML%nEns
        randomMemberIndexArray(memberIndex) = memberIndex
      end do
      call utl_randomOrderInt(randomMemberIndexArray,randomSeed)
      if (firstCall) then
        write(*,*) 'enkf_LETKFcomputeWeights: seed for random shuffle of sub ens = ', randomSeed
        write(*,*) 'enkf_LETKFcomputeWeights: randomOrder = ', randomMemberIndexArray(:)
      end if
      do subEnsIndex = 1, enkfNML%numSubEns
        do memberIndex = 1, nEnsPerSubEns
          memberIndexSubEns(memberIndex,subEnsIndex) =  &
               randomMemberIndexArray((subEnsIndex-1)*nEnsPerSubEns + memberIndex)
        end do
      end do
      if ( useModulatedEns ) then
        do subEnsIndex = 1, enkfNML%numSubEns
          memberIndex2 = 0
          do memberIndex = 1, nEnsPerSubEns
            do eigenVectorColumnIndex = 1, enkfNML%numRetainedEigen
              memberIndex2 = memberIndex2 + 1
              memberIndexSubEns_mod(memberIndex2,subEnsIndex) =  &
                   randomMemberIndexArray((subEnsIndex-1)*nEnsPerSubEns + memberIndex) + &
                   (eigenVectorColumnIndex - 1) * enkfNML%nEns
            end do
          end do
        end do
      end if
    end if

    do subEnsIndex = 1, enkfNML%numSubEns
      memberIndex = 1
      do subEnsIndex2 = 1, enkfNML%numSubEns
        if (subEnsIndex2 == subEnsIndex) cycle

        if ( .not. useModulatedEns ) then
          memberIndexSubEnsComp(memberIndex:memberIndex+nEnsPerSubEns-1,subEnsIndex) =  &
               memberIndexSubEns(:,subEnsIndex2)
          memberIndex = memberIndex + nEnsPerSubEns
        else
          memberIndexSubEnsComp(memberIndex:memberIndex+nEnsPerSubEns_mod-1,subEnsIndex) =  &
               memberIndexSubEns_mod(:,subEnsIndex2)
          memberIndex = memberIndex + nEnsPerSubEns_mod
        end if
      end do
    end do

    firstCall = .false.

  end subroutine enkf_defineSubEnsembles

  !----------------------------------------------------------------------
  ! enkf_applyEnsWeights (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_applyEnsWeights(enkfNML, stateVectorMeanInc, stateVectorMeanTrl, stateVectorMeanAnl, &
                                  ensembleTrl, ensembleAnl, levIndex,  &
                                  weightsMean, weightsMembers, useModulatedEns, &
                                  myLonBegHalo, myLatBegHalo)
    !
    !:Purpose: Use the computed weights and the background ensemble to compute analysis
    !          ensemble and ensemble mean. This calculation is done for just 1 level per call.
    !
    implicit none

    ! Arguments:
    type(struct_enkfNML), intent(in)    :: enkfNML            ! Derived type variable with namelist variables
    type(struct_gsv),     intent(inout) :: stateVectorMeanInc ! Ensemble mean increment
    type(struct_gsv),     intent(inout) :: stateVectorMeanTrl ! Ensemble mean trial
    type(struct_gsv),     intent(inout) :: stateVectorMeanAnl ! Ensemble mean analysis
    type(struct_ens),     intent(inout) :: ensembleTrl        ! Trial ensemble
    type(struct_ens),     intent(inout) :: ensembleAnl        ! Analysis ensemble
    integer,              intent(in)    :: levIndex           ! The `levIndex` being processed
    integer,              intent(in)    :: myLonBegHalo       ! First lon index of weights (for array indexing)
    integer,              intent(in)    :: myLatBegHalo       ! First lat index of weights (for array indexing)
    real(8),              intent(in)    :: weightsMean(1:,1:,myLonBegHalo:,myLatBegHalo:)    ! Weights for ens mean
    real(8),              intent(in)    :: weightsMembers(1:,1:,myLonBegHalo:,myLatBegHalo:) ! Weights for members
    logical,              intent(in)    :: useModulatedEns    ! Indicates if modulated ens is used

    ! Locals:
    character(len=4) :: varName
    integer :: stepIndex, memberIndex, memberIndex1, memberIndex2, memberIndexInModEns
    integer :: eigenVectorColumnIndex, latIndex, lonIndex, levIndex2, varLevIndex
    integer :: myLonBeg, myLonEnd, myLatBeg, myLatEnd, nLev_M, nLev_depth, numVarLev
    real(4) :: modulationFactor_r4, pert_r4
    real(4), pointer     :: meanTrl_ptr_r4(:,:,:,:), meanAnl_ptr_r4(:,:,:,:), meanInc_ptr_r4(:,:,:,:)
    real(4), pointer     :: memberTrl_ptr_r4(:,:,:,:), memberAnl_ptr_r4(:,:,:,:)
    real(8), allocatable :: memberAnlPert(:)
    integer,          allocatable, save :: levIndex2FromVarLevIndex(:)
    character(len=4), allocatable, save :: varLevelFromVarLev(:)
    character(len=2), allocatable, save :: varKindFromVarLev(:)
    integer,          allocatable, save :: levFromVarLev(:)
    logical, save :: firstCall = .true.

    myLonBeg = stateVectorMeanAnl%myLonBeg
    myLonEnd = stateVectorMeanAnl%myLonEnd
    myLatBeg = stateVectorMeanAnl%myLatBeg
    myLatEnd = stateVectorMeanAnl%myLatEnd
    nLev_M     = ens_getNumLev(ensembleAnl, 'MM')
    nLev_depth = ens_getNumLev(ensembleAnl, 'DP')
    numVarLev  = stateVectorMeanAnl%numVarLev

    allocate(memberAnlPert(enkfNML%nEns))
    call gsv_getField(stateVectorMeanInc,meanInc_ptr_r4)
    call gsv_getField(stateVectorMeanTrl,meanTrl_ptr_r4)
    call gsv_getField(stateVectorMeanAnl,meanAnl_ptr_r4)

    if (firstCall) then
      allocate(varLevelFromVarLev(numVarLev))
      allocate(levFromVarLev(numVarLev))
      allocate(varKindFromVarLev(numVarLev))
      do varLevIndex = 1, numVarLev
        varName = gsv_getVarNameFromVarLev(stateVectorMeanInc,varLevIndex)
        varLevelFromVarLev(varLevIndex) = vnl_varLevelFromVarname(varName)
        levFromVarLev(varLevIndex) = gsv_getLevFromVarLev(stateVectorMeanInc,varLevIndex)
        varKindFromVarLev(varLevIndex) = vnl_varKindFromVarname(varName)
      end do

      allocate(levIndex2FromVarLevIndex(numVarLev))
      do varLevIndex = 1, numVarLev
        ! Only treat varLevIndex values that correspond with current levIndex
        if (varLevelFromVarLev(varLevIndex) == 'SF'   .or. varLevelFromVarLev(varLevIndex) == 'SFMM' .or. &
            varLevelFromVarLev(varLevIndex) == 'SFTH' .or. varLevelFromVarLev(varLevIndex) == 'SS') then
          if (varKindFromVarLev(varLevIndex) == 'OC') then
            levIndex2 = 1
          else
            levIndex2 = max(nLev_M,nLev_depth)
          end if
        else if (varLevelFromVarLev(varLevIndex) == 'MM' .or. varLevelFromVarLev(varLevIndex) == 'TH' .or. &
                 varLevelFromVarLev(varLevIndex) == 'DP') then
          levIndex2 = levFromVarLev(varLevIndex)
        else if (varLevelFromVarLev(varLevIndex) == 'OT') then
          ! Most (all?) variables using the 'other' coordinate are surface
          levIndex2 = max(nLev_M,nLev_depth)
        else
          write(*,*) 'varLevel = ', varLevelFromVarLev(varLevIndex)
          call utl_abort('enkf_LETKFanalyses: unknown varLevel')
        end if
        levIndex2FromVarLevIndex(varLevIndex) = levIndex2
      end do
    end if ! firstCall

    !$OMP PARALLEL DO PRIVATE(latIndex, lonIndex, varLevIndex, levIndex2, memberTrl_ptr_r4, memberAnl_ptr_r4), &
    !$OMP PRIVATE(memberAnlPert, stepIndex, memberIndex, memberIndex2, memberIndex1, eigenVectorColumnIndex, pert_r4), &
    !$OMP PRIVATE(memberIndexInModEns, modulationFactor_r4)
    do latIndex = myLatBeg, myLatEnd
      LON_LOOP5: do lonIndex = myLonBeg, myLonEnd

        ! skip this grid point if all weights zero (no nearby obs)
        if (all(weightsMean(:,1,lonIndex,latIndex) == 0.0d0)) cycle LON_LOOP5

        ! Compute the ensemble mean increment and analysis
        do varLevIndex = 1, numVarLev
          levIndex2 = levIndex2FromVarLevIndex(varLevIndex)
          if (levIndex2 /= levIndex .and. .not. useModulatedEns) cycle
          memberTrl_ptr_r4 => ens_getOneLev_r4(ensembleTrl,varLevIndex)
          do stepIndex = 1, tim_nstepobsinc
            ! mean increment
            if ( useModulatedEns ) then
              do eigenVectorColumnIndex = 1, enkfNML%numRetainedEigen
                call getModulationFactor( enkfNML, stateVectorMeanInc%vco, levIndex2, &
                                          eigenVectorColumnIndex, &
                                          modulationFactor_r4 )

                do memberIndex = 1, enkfNML%nEns
                  pert_r4 = modulationFactor_r4 * ( memberTrl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) -  &
                                                    meanTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) )

                  ! Index of the modulated ensemble member corresponding to original
                  ! ensemble member index (memberIndex1) and eigenVectorColumnIndex.
                  memberIndexInModEns = (eigenVectorColumnIndex - 1) * enkfNML%nEns + memberIndex

                  meanInc_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) =  &
                       meanInc_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) +  &
                       weightsMean(memberIndexInModEns,1,lonIndex,latIndex) * pert_r4
                end do
              end do
            else
              do memberIndex = 1, enkfNML%nEns
                meanInc_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) =  &
                     meanInc_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) +  &
                     weightsMean(memberIndex,1,lonIndex,latIndex) *  &
                     (memberTrl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) -  &
                      meanTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex))
              end do
            end if

            ! mean analysis
            meanAnl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) =  &
                 meanTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) +  &
                 meanInc_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex)
          end do ! stepIndex
        end do ! varLevIndex

        ! Compute the ensemble member analyses
        call utl_tmg_start(144,'------ApplyWeightsMember')
        do varLevIndex = 1, numVarLev
          levIndex2 = levIndex2FromVarLevIndex(varLevIndex)
          if (levIndex2 /= levIndex .and. .not. useModulatedEns) cycle
          memberTrl_ptr_r4 => ens_getOneLev_r4(ensembleTrl,varLevIndex)
          memberAnl_ptr_r4 => ens_getOneLev_r4(ensembleAnl,varLevIndex)
          do stepIndex = 1, tim_nstepobsinc

            ! Compute analysis member perturbation
            memberAnlPert(:) = 0.0d0

            if ( useModulatedEns ) then
              do memberIndex2 = 1, enkfNML%nEns
                do eigenVectorColumnIndex = 1, enkfNML%numRetainedEigen
                  call getModulationFactor( enkfNML, stateVectorMeanInc%vco, levIndex2, &
                                            eigenVectorColumnIndex, &
                                            modulationFactor_r4 )

                  do memberIndex1 = 1, enkfNML%nEns
                    ! Compute background ensemble perturbations for the modulated ensemble (Xb_Mod)
                    pert_r4 = modulationFactor_r4 * ( memberTrl_ptr_r4(memberIndex1,stepIndex,lonIndex,latIndex) -  &
                                                      meanTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) )

                    ! Index of the modulated ensemble member corresponding to original
                    ! ensemble member index (memberIndex1) and eigenVectorColumnIndex.
                    memberIndexInModEns = (eigenVectorColumnIndex - 1) * enkfNML%nEns + memberIndex1

                    ! sum Xb_Mod * Wa over all modulated ensembles to get member perturbations for
                    !   original ensemble (memberIndex2)
                    memberAnlPert(memberIndex2) = memberAnlPert(memberIndex2) + &
                         weightsMembers(memberIndexInModEns,memberIndex2,lonIndex,latIndex) *  pert_r4
                  end do
                end do

                ! Compute final member perturbations by removing background original ensemble perturbations
                memberAnlPert(memberIndex2) = (memberTrl_ptr_r4(memberIndex2,stepIndex,lonIndex,latIndex) -  &
                                               meanTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex)) + &
                                               memberAnlPert(memberIndex2)

              end do ! memberIndex2
            else
              do memberIndex2 = 1, enkfNML%nEns
                do memberIndex1 = 1, enkfNML%nEns
                  memberAnlPert(memberIndex2) = memberAnlPert(memberIndex2) + &
                       weightsMembers(memberIndex1,memberIndex2,lonIndex,latIndex) *  &
                       ( memberTrl_ptr_r4(memberIndex1,stepIndex,lonIndex,latIndex) -  &
                         meanTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) )
                end do ! memberIndex1
              end do ! memberIndex2
            end if

            ! Add analysis member perturbation to mean analysis
            memberAnl_ptr_r4(:,stepIndex,lonIndex,latIndex) =  &
                 meanAnl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) + memberAnlPert(:)


          end do ! stepIndex
        end do ! varLevIndex
        call utl_tmg_stop(144)

      end do LON_LOOP5
    end do
    !$OMP END PARALLEL DO

    firstCall = .false.

  end subroutine enkf_applyEnsWeights

  !----------------------------------------------------------------------
  ! enkf_calcYbTinvRYb (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_calcYbTinvRYb(nEns1, nEns2, maxNumLocalObs, numLocalObs, &
                                YbTinvRYb_pert, YbTinvR_pert, &
                                ensObs_mpiglobal, localBodyIndices,  &
                                YbTinvRYb_mean, YbTinvR_mean)
    !
    !:Purpose: Compute the background covariance in ensemble space.
    !
    implicit none

    ! Arguments:
    integer,          intent(in)            :: nEns1                              ! First dimension of cov matrix in ens space
    integer,          intent(in)            :: nEns2                              ! Second dimension of cov matrix in ens space
    integer,          intent(in)            :: maxNumLocalObs                     ! Maximum number of local obs
    integer,          intent(in)            :: numLocalObs                        ! Actual number of local obs
    type(struct_eob), intent(in)            :: ensObs_mpiglobal                   ! Ensemble observations
    integer,          intent(in)            :: localBodyIndices(maxNumLocalObs)   ! List of body indexes for local obs
    real(8),          intent(out)           :: YbTinvRYb_pert(nEns1,nEns2)        ! Output cov matrix in ens space
    real(8),          intent(in)            :: YbTinvR_pert(nEns1,maxNumLocalObs) ! Input product of Yb and inv(R)
    real(8),          intent(out), optional :: YbTinvRYb_mean(nEns1,nEns2)        ! Output cov matrix in ens space for mean
    real(8),          intent(in),  optional :: YbTinvR_mean(nEns1,maxNumLocalObs) ! Input product of Yb and inv(R) for mean

    ! Locals:
    integer :: memberIndex, localObsIndex, bodyIndex
    logical :: isSymmetric
    real(8), allocatable :: YbCopy_r8(:,:)

    if ( numLocalObs == 0 ) then
      write(*,*) 'enkf_calcYbTinvRYb called with numLocalObs = 0'
      YbTinvRYb_pert(:,:) = 0.0d0
      if (present(YbTinvRYb_mean)) then
        YbTinvRYb_mean(:,:) = 0.0d0
      end if
      return
    end if

    allocate(YbCopy_r8(nEns2,numLocalObs))

    call utl_tmg_start(137,'--------YbArraysCopy')
    !$OMP PARALLEL DO PRIVATE (localObsIndex, bodyIndex, memberIndex)
    do localObsIndex = 1, numLocalObs
      bodyIndex = localBodyIndices(localObsIndex)
      do memberIndex = 1, nEns2
        YbCopy_r8(memberIndex,localObsIndex) = real(ensObs_mpiglobal%Yb_r4(memberIndex,bodyIndex), 8)
      end do
    end do
    !$OMP END PARALLEL DO
    call utl_tmg_stop(137)

    ! When nEns1 equals nEns2 we can assume the output matrix will be symmetric
    isSymmetric = nEns1 == nEns2

    call utl_fastMatMul(YbTinvR_pert, YbCopy_r8, YbTinvRYb_pert,                 &
                        isATransposed_opt = .false., isBTransposed_opt = .true., &
                        isCSymmetric_opt = isSymmetric, summationDim_opt = numLocalObs)

    if (eob_simObsAssim .and. present(YbTinvRYb_mean)) then
      call utl_fastMatMul(YbTinvR_mean, YbCopy_r8, YbTinvRYb_mean,                 &
                          isATransposed_opt = .false., isBTransposed_opt = .true., &
                          isCSymmetric_opt = isSymmetric, summationDim_opt = numLocalObs)
    end if

    deallocate(YbCopy_r8)

  end subroutine enkf_calcYbTinvRYb

  !----------------------------------------------------------------------
  ! enkf_computeVertLocation (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_computeVertLocation(vertLocation_r4,stateVectorMeanTrl)
    !
    !:Purpose:  Compute extract global 3D vertical location field from supplied
    !           stateVector. Can be either logPressure or depth levels.
    !
    implicit none

    ! Arguments:
    real(4), allocatable, intent(inout) :: vertLocation_r4(:,:,:) ! Vertical locations for all grid points
    type(struct_gsv),     intent(inout) :: stateVectorMeanTrl     ! Ensemble mean state vector

    ! Locals:
    integer          :: nLev_M, nLev_depth, nLev_vertLocation, levIndex, nsize, ierr
    real(4), pointer :: vertLocation_ptr_r4(:,:,:)
    type(struct_gsv) :: stateVectorMeanTrlPressure
    type(struct_gsv) :: stateVectorMeanTrlPressure_1step

    write(*,*) 'enkf_computeVertLocation: starting'

    nLev_M = gsv_getNumLev(stateVectorMeanTrl, 'MM')
    nLev_depth = gsv_getNumLev(stateVectorMeanTrl, 'DP')
    if ( nLev_M > 0 .and. nLev_depth > 0 ) then
      call utl_abort('enkf_computeVertLocation: both momentum and depth levels exist.')
    else if ( nLev_M == 0 .and. nLev_depth == 0 ) then
      call utl_abort('enkf_computeVertLocation: neither momentum nor depth levels exist.')
    end if
    nLev_vertLocation = max(nLev_M, nLev_depth)

    allocate(vertLocation_r4(stateVectorMeanTrl%hco%ni, &
                             stateVectorMeanTrl%hco%nj, &
                             nLev_vertLocation))

    if ( nLev_M > 0 ) then ! log pressure for NWP fields

      ! Compute background ens mean 3D log pressure and make mpiglobal for vertical localization
      call gsv_allocate( stateVectorMeanTrlPressure, tim_nstepobsinc,  &
                         stateVectorMeanTrl%hco, stateVectorMeanTrl%vco, dateStamp_opt=tim_getDateStamp(),  &
                         mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                         dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0 ','P_M','P_T','Z_M','Z_T','TT ','HU '/) )
      call gsv_zero(stateVectorMeanTrlPressure)
      call gsv_copy(stateVectorMeanTrl, stateVectorMeanTrlPressure, allowVarMismatch_opt=.true.)
      call gvt_transform(stateVectorMeanTrlPressure,'ZandP_nl')
      if (mmpi_myid == 0) then
        call gsv_allocate( stateVectorMeanTrlPressure_1step, 1,  &
                           stateVectorMeanTrl%hco, stateVectorMeanTrl%vco, dateStamp_opt=tim_getDateStamp(),  &
                           mpi_local_opt=.false., &
                           dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0 ','P_M','P_T','Z_M','Z_T','TT ','HU '/) )
      end if
      call gsv_transposeTilesToStep(stateVectorMeanTrlPressure_1step, stateVectorMeanTrlPressure, (tim_nstepobsinc+1)/2)
      call gsv_deallocate(stateVectorMeanTrlPressure)
      if (mmpi_myid == 0) then
        call gsv_getField(stateVectorMeanTrlPressure_1step,vertLocation_ptr_r4,'P_M')
        vertLocation_r4(:,:,:) = log(vertLocation_ptr_r4(:,:,:))
        write(*,*) 'enkf_computeVertLocation: vertLocation min/max = ', minval(vertLocation_r4), maxval(vertLocation_r4)
      end if
      nsize = stateVectorMeanTrlPressure%ni * stateVectorMeanTrlPressure%nj * nLev_M
      call rpn_comm_bcast(vertLocation_r4, nsize, 'mpi_real4', 0, 'GRID', ierr)

    else if ( nLev_depth > 0 ) then ! depth for ocean fields

      ! fill in all horizontal grid points with the same profile of depth values
      do levIndex = 1, nLev_depth
        write(*,*) 'setting vertLocation for levIndex =', levIndex, &
                   ', depth = ', stateVectorMeanTrl%vco%depths(levIndex)
        vertLocation_r4(:,:,levIndex) = stateVectorMeanTrl%vco%depths(levIndex)
      end do

    end if

    write(*,*) 'enkf_computeVertLocation: finished'

  end subroutine enkf_computeVertLocation

  !----------------------------------------------------------------------
  ! enkf_LETKFsetupMpiDistribution (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_LETKFsetupMpiDistribution(numLatLonMpiGlobal, myNumLatLonRecv,  &
                                            myLatIndexesRecv, myLonIndexesRecv, &
                                            latIndexesSendMpiGlobal, lonIndexesSendMpiGlobal, &
                                            procIndexesSendMpiGlobal, &
                                            numProcsSendMpiGlobal, wInterpInfo)
    !
    ! :Purpose: Setup for distribution of grid points over mpi tasks by building a global
    !           list of which MPI tasks need the computed weights for which grid points.
    !           Also build local list for grid points required to recv.
    !
    implicit none

    ! Arguments:
    integer,                     intent(out) :: numLatLonMpiGlobal            ! Total number of weights to be calculated
    integer,                     intent(out) :: myNumLatLonRecv               ! Number of weights needed locally (including halo)
    integer, allocatable,        intent(out) :: myLatIndexesRecv(:)           ! List of latIndex for locally needed weights
    integer, allocatable,        intent(out) :: myLonIndexesRecv(:)           ! List of lonIndex for locally needed weights
    integer, allocatable,        intent(out) :: latIndexesSendMpiGlobal(:)    ! Global list of latIndex
    integer, allocatable,        intent(out) :: lonIndexesSendMpiGlobal(:)    ! Global list of lonIndex
    integer, allocatable,        intent(out) :: procIndexesSendMpiGlobal(:,:) ! MPI task indexes where each weight is needed
    integer, allocatable,        intent(out) :: numProcsSendMpiGlobal(:)      ! Number of MPI tasks where each weight is needed
    type(struct_enkfInterpInfo), intent(in)  :: wInterpInfo                   ! LETKF weight interpolation info

    ! Locals:
    integer :: latIndex, lonIndex, procIndex, latLonIndex, myLatLonIndex, latLonIndexMpiGlobal
    integer :: myLonBeg, myLonEnd, myLatBeg, myLatEnd
    integer :: myLonBegHalo, myLonEndHalo, myLatBegHalo, myLatEndHalo
    integer :: numLatLonRecvMax, myNumLatLon, numLatLonMax, ierr
    integer, allocatable :: allLatIndexesRecv(:,:), allLonIndexesRecv(:,:)
    integer, allocatable :: allLatIndexesSend(:,:), allLonIndexesSend(:,:)
    integer, allocatable :: localLatIndexesSend(:), localLonIndexesSend(:)
    integer              :: allNumLatLonRecv(mmpi_nprocs), allNumLatLon(mmpi_nprocs)

    myLonBeg = wInterpInfo%myLonBeg
    myLonEnd = wInterpInfo%myLonEnd
    myLatBeg = wInterpInfo%myLatBeg
    myLatEnd = wInterpInfo%myLatEnd

    myLonBegHalo = wInterpInfo%myLonBegHalo
    myLonEndHalo = wInterpInfo%myLonEndHalo
    myLatBegHalo = wInterpInfo%myLatBegHalo
    myLatEndHalo = wInterpInfo%myLatEndHalo

    write(*,*) 'enkf_LETKFsetupMpiDistribution: starting'
    call msg_memUsage('enkf_LETKFsetupMpiDistribution')

    ! First, count the number of grid points where weights needed locally (for recv-ing)
    ! Note, this includes the "halo" needed for interpolation
    myNumLatLonRecv = 0
    do latIndex = myLatBegHalo, myLatEndHalo
      LON_LOOP1: do lonIndex = myLonBegHalo, myLonEndHalo
        ! If this lat-lon is to be interpolated, then skip calculation
        if (wInterpInfo%numIndexes(lonIndex,latIndex) > 0) cycle LON_LOOP1
        myNumLatLonRecv = myNumLatLonRecv + 1
      end do LON_LOOP1
    end do

    ! Communicate to all mpi tasks
    call rpn_comm_allgather(myNumLatLonRecv, 1, "mpi_integer",  &
                            allNumLatLonRecv, 1,"mpi_integer", "GRID", ierr)
    numLatLonRecvMax = maxval(allNumLatLonRecv)
    write(*,*) 'enkf_LETKFsetupMpiDistribution: allNumLatLonRecv =', allNumLatLonRecv(:)
    write(*,*) 'enkf_LETKFsetupMpiDistribution: numLatLonRecvSum =', sum(allNumLatLonRecv)
    write(*,*) 'enkf_LETKFsetupMpiDistribution: numLatLonRecvMax =', numLatLonRecvMax

    ! Now create a list of grid point indexes where weights needed locally (for recv-ing)
    allocate(myLatIndexesRecv(numLatLonRecvMax))
    allocate(myLonIndexesRecv(numLatLonRecvMax))
    myLatIndexesRecv(:) = -1
    myLonIndexesRecv(:) = -1
    myNumLatLonRecv = 0
    do latIndex = myLatBegHalo, myLatEndHalo
      LON_LOOP2: do lonIndex = myLonBegHalo, myLonEndHalo
        ! If this lat-lon is to be interpolated, then skip calculation
        if (wInterpInfo%numIndexes(lonIndex,latIndex) > 0) cycle LON_LOOP2
        myNumLatLonRecv = myNumLatLonRecv + 1

        myLatIndexesRecv(myNumLatLonRecv) = latIndex
        myLonIndexesRecv(myNumLatLonRecv) = lonIndex
      end do LON_LOOP2
    end do

    ! Communicate to all mpi tasks this list of grid point lat-lon indexes
    allocate(allLatIndexesRecv(numLatLonRecvMax, mmpi_nprocs))
    allocate(allLonIndexesRecv(numLatLonRecvMax, mmpi_nprocs))
    call rpn_comm_allgather(myLatIndexesRecv, numLatLonRecvMax, "mpi_integer",  &
                            allLatIndexesRecv, numLatLonRecvMax, "mpi_integer",  &
                            "GRID", ierr)
    call rpn_comm_allgather(myLonIndexesRecv, numLatLonRecvMax, "mpi_integer",  &
                            allLonIndexesRecv, numLatLonRecvMax, "mpi_integer",  &
                            "GRID", ierr)

    ! Now count number of local grid points without the halo
    myNumLatLon = 0
    do latIndex = myLatBeg, myLatEnd
      LON_LOOP3: do lonIndex = myLonBeg, myLonEnd
        ! If this lat-lon is to be interpolated, then skip calculation
        if (wInterpInfo%numIndexes(lonIndex,latIndex) > 0) cycle LON_LOOP3
        myNumLatLon = myNumLatLon + 1
      end do LON_LOOP3
    end do

    ! Communicate to all mpi tasks
    call rpn_comm_allreduce(myNumLatLon, numLatLonMpiGlobal, &
                            1,"mpi_integer","mpi_sum","GRID",ierr)
    call rpn_comm_allgather(myNumLatLon,  1, "mpi_integer",  &
                            allNumLatLon, 1, "mpi_integer",  &
                            "GRID", ierr)
    numLatLonMax = maxval(allNumLatLon)

    ! Build global lists of lat-lon indexes and list of mpi tasks where each needs to be sent
    allocate(allLatIndexesSend(numLatLonMax,mmpi_nprocs))
    allocate(allLonIndexesSend(numLatLonMax,mmpi_nprocs))
    allocate(latIndexesSendMpiGlobal(numLatLonMpiGlobal))
    allocate(lonIndexesSendMpiGlobal(numLatLonMpiGlobal))
    allocate(localLatIndexesSend(numLatLonMax))
    allocate(localLonIndexesSend(numLatLonMax))
    latIndexesSendMpiGlobal(:) = -1
    lonIndexesSendMpiGlobal(:) = -1
    localLatIndexesSend(:) = -1
    localLonIndexesSend(:) = -1

    myNumLatLon = 0
    do latIndex = myLatBeg, myLatEnd
      LON_LOOP4: do lonIndex = myLonBeg, myLonEnd
        ! If this lat-lon is to be interpolated, then skip calculation
        if (wInterpInfo%numIndexes(lonIndex,latIndex) > 0) cycle LON_LOOP4
        myNumLatLon = myNumLatLon + 1

        localLatIndexesSend(myNumLatLon) = latIndex
        localLonIndexesSend(myNumLatLon) = lonIndex
      end do LON_LOOP4
    end do

    call rpn_comm_allgather(localLatIndexesSend, numLatLonMax, "mpi_integer",  &
                            allLatIndexesSend,   numLatLonMax, "mpi_integer",  &
                            "GRID", ierr)
    call rpn_comm_allgather(localLonIndexesSend, numLatLonMax, "mpi_integer",  &
                            allLonIndexesSend,   numLatLonMax, "mpi_integer",  &
                            "GRID", ierr)

    ! Reorganize into single dimension list
    latLonIndexMpiGlobal = 0
    do procIndex = 1, mmpi_nprocs
      do myLatLonIndex = 1, allNumLatLon(procIndex)
        latLonIndexMpiGlobal = latLonIndexMpiGlobal + 1

        latIndexesSendMpiGlobal(latLonIndexMpiGlobal) = allLatIndexesSend(myLatLonIndex, procIndex)
        lonIndexesSendMpiGlobal(latLonIndexMpiGlobal) = allLonIndexesSend(myLatLonIndex, procIndex)
      end do
    end do

    ! Figure out which mpi tasks I will need to send my results to
    ! -- First just do counting for allocating procIndexesSendMpiGlobal
    allocate(numProcsSendMpiGlobal(numLatLonMpiGlobal))
    numProcsSendMpiGlobal(:) = 0
    do latLonIndex = 1, numLatLonMpiGlobal
      do procIndex = 1, mmpi_nprocs
        if ( any( (latIndexesSendMpiGlobal(latLonIndex) == allLatIndexesRecv(1:allNumLatLonRecv(procIndex), procIndex)) .and.  &
                  (lonIndexesSendMpiGlobal(latLonIndex) == allLonIndexesRecv(1:allNumLatLonRecv(procIndex), procIndex)) ) ) then
          numProcsSendMpiGlobal(latLonIndex) = numProcsSendMpiGlobal(latLonIndex) + 1
        end if
      end do
    end do
    ! -- Now allocate and store values in procIndexesSendMpiGlobal
    allocate(procIndexesSendMpiGlobal(numLatLonMpiGlobal,maxval(numProcsSendMpiGlobal(:))))
    procIndexesSendMpiGlobal(:,:) = -1
    numProcsSendMpiGlobal(:) = 0
    do latLonIndex = 1, numLatLonMpiGlobal
      do procIndex = 1, mmpi_nprocs
        if ( any( (latIndexesSendMpiGlobal(latLonIndex) == allLatIndexesRecv(1:allNumLatLonRecv(procIndex), procIndex)) .and.  &
                  (lonIndexesSendMpiGlobal(latLonIndex) == allLonIndexesRecv(1:allNumLatLonRecv(procIndex), procIndex)) ) ) then
          numProcsSendMpiGlobal(latLonIndex) = numProcsSendMpiGlobal(latLonIndex) + 1
          procIndexesSendMpiGlobal(latLonIndex,numProcsSendMpiGlobal(latLonIndex)) = procIndex
        end if
      end do
    end do

    write(*,*) 'enkf_LETKFsetupMpiDistribution: lat/lon/proc indexes I need to receive:'
    do latLonIndex = 1, myNumLatLonRecv
      write(*,*) myLatIndexesRecv(latLonIndex), myLonIndexesRecv(latLonIndex)
    end do

    if (mmpi_myid == 0) then
      write(*,*) 'enkf_LETKFsetupMpiDistribution: number of lat/lon indexes globally =', numLatLonMpiGlobal
      write(*,*) 'enkf_LETKFsetupMpiDistribution: global list of indexMpiGlobal/lat/lon/proc indexes:'
      do latLonIndexMpiGlobal = 1, numLatLonMpiGlobal
        write(*,*) latLonIndexMpiGlobal,  &
                   latIndexesSendMpiGlobal(latLonIndexMpiGlobal),  &
                   lonIndexesSendMpiGlobal(latLonIndexMpiGlobal),  &
                   procIndexesSendMpiGlobal(latLonIndexMpiGlobal,1:numProcsSendMpiGlobal(latLonIndexMpiGlobal))
      end do
    end if

    write(*,*) 'enkf_LETKFsetupMpiDistribution: done'
    call msg_memUsage('enkf_LETKFsetupMpiDistribution')

  end subroutine enkf_LETKFsetupMpiDistribution

  !----------------------------------------------------------------------
  ! enkf_LETKFgetMpiGlobalTags (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_LETKFgetMpiGlobalTags(latLonTagMpiGlobal,myLatIndexesRecv,myLonIndexesRecv)
    !
    !:Purpose: Define a set of unique MPI tags used for doing MPI communication
    !          of the LETKF weights.
    !
    implicit none

    ! Arguments:
    integer, intent(out) :: latLonTagMpiGlobal(:,:) ! Output list of unique MPI tags for weight communication
    integer, intent(in)  :: myLatIndexesRecv(:)     ! Input latIndex list for locally needed weights
    integer, intent(in)  :: myLonIndexesRecv(:)     ! Input lonIndex list for locally needed weights

    ! Locals:
    integer :: ierr, ni, nj, lonIndex, latIndex
    integer :: latPerPE, latPerPEmax, myLatBeg, myLatEnd
    integer :: lonPerPE, lonPerPEmax, myLonBeg, myLonEnd
    integer :: countTags, myNumLatLonRecv, numLatLonRecvMax
    integer, allocatable :: allNumLatLonRecv(:)
    integer, allocatable :: allLatIndexesRecv(:,:), allLonIndexesRecv(:,:)
    logical, allocatable :: tagNeededMpiLocal(:,:), tagNeededMpiGlobal(:,:)

    write(*,*) 'enkf_LETKFgetMpiGlobalTags: Starting'

    ni = size(latLonTagMpiGlobal,1)
    nj = size(latLonTagMpiGlobal,2)

    ! Set up local MPI tiles to speed up calculation
    call mmpi_setup_latbands(nj, latPerPE, latPerPEmax, myLatBeg, myLatEnd)
    call mmpi_setup_lonbands(ni, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd)

    myNumLatLonRecv = size(myLatIndexesRecv)
    allocate(allNumLatLonRecv(mmpi_nprocs))
    call rpn_comm_allgather(myNumLatLonRecv, 1, "mpi_integer",  &
                            allNumLatLonRecv, 1,"mpi_integer", "GRID", ierr)
    numLatLonRecvMax = maxval(allNumLatLonRecv)

    ! Communicate to all mpi tasks this list of grid point lat-lon indexes
    allocate(allLatIndexesRecv(numLatLonRecvMax, mmpi_nprocs))
    allocate(allLonIndexesRecv(numLatLonRecvMax, mmpi_nprocs))
    call rpn_comm_allgather(myLatIndexesRecv, numLatLonRecvMax, "mpi_integer",  &
                            allLatIndexesRecv, numLatLonRecvMax, "mpi_integer",  &
                            "GRID", ierr)
    call rpn_comm_allgather(myLonIndexesRecv, numLatLonRecvMax, "mpi_integer",  &
                            allLonIndexesRecv, numLatLonRecvMax, "mpi_integer",  &
                            "GRID", ierr)

    ! Determine grid points where weights are calculated - split work over MPI tasks

    allocate(tagNeededMpiLocal(ni,nj))
    allocate(tagNeededMpiGlobal(ni,nj))
    tagNeededMpiLocal(:,:) = .false.
    tagNeededMpiGlobal(:,:) = .false.

    !$OMP PARALLEL DO PRIVATE(latIndex, lonIndex)
    do lonIndex = myLonBeg, myLonEnd
      do latIndex = myLatBeg, myLatEnd
        if (any(lonIndex == allLonIndexesRecv(:,:) .and. latIndex == allLatIndexesRecv(:,:))) then
          tagNeededMpiLocal(lonIndex,latIndex) = .true.
        end if
      end do
    end do
    !$OMP END PARALLEL DO
    call rpn_comm_allreduce(tagNeededMpiLocal, tagNeededMpiGlobal, ni*nj, &
                            'mpi_logical','mpi_lor','GRID',ierr)

    ! Loop over global grid points with calculated weights to determine unique tag values
    countTags = 0
    latLonTagMpiGlobal(:,:) = 0
    do lonIndex = 1, ni
      do latIndex = 1, nj
        if (tagNeededMpiGlobal(lonIndex,latIndex)) then
          countTags = countTags + 1
          latLonTagMpiGlobal(lonIndex,latIndex) = countTags
        end if
      end do
    end do
    write(*,*) 'number of Recv grid points found = ', maxval(latLonTagMpiGlobal(:,:))

    if (maxval(latLonTagmpiGlobal(:,:)) > (mmpi_maxTagValue - 2)) then
      write(*,*) 'maximum allowable tag value = ', mmpi_maxTagValue - 2
      call utl_abort('enkf_LETKFgetMpiGlobalTags: mpi tag values exceeded max allowable value')
    end if

    deallocate(tagNeededMpiLocal)
    deallocate(tagNeededMpiGlobal)

    write(*,*) 'enkf_LETKFgetMpiGlobalTags: Finished'

  end subroutine enkf_LETKFgetMpiGlobalTags

  !--------------------------------------------------------------------------
  ! enkf_setupInterpInfo
  !--------------------------------------------------------------------------
  subroutine enkf_setupInterpInfo(wInterpInfo, hco, weightLatLonStep,  &
                                  myLonBeg,myLonEnd,myLatBeg,myLatEnd)
    !
    ! :Purpose: Setup the weights and lat/lon indices needed to bilinearly
    !           interpolate the LETKF weights from a coarse grid to the full
    !           resolution grid. The coarseness of the grid is specified by
    !           the weightLatLonStep argument.
    !
    implicit none

    ! Arguments:
    type(struct_enkfInterpInfo), intent(out) :: wInterpInfo      ! Output LETKF weight interpolation info
    type(struct_hco),            intent(in)  :: hco              ! Horizontal coordinate definition
    integer,                     intent(in)  :: weightLatLonStep ! Grid-point spacing of weight calculation
    integer,                     intent(in)  :: myLonBeg         ! Limits of local lat-lon tile
    integer,                     intent(in)  :: myLonEnd         ! Limits of local lat-lon tile
    integer,                     intent(in)  :: myLatBeg         ! Limits of local lat-lon tile
    integer,                     intent(in)  :: myLatEnd         ! Limits of local lat-lon tile

    ! Locals:
    integer :: lonIndex, latIndex, ni, nj
    integer :: myLonBegHalo, myLonEndHalo, myLatBegHalo, myLatEndHalo
    real(8) :: interpWeightLon, interpWeightLat
    logical :: includesYinYangBndry

    ni = hco%ni
    nj = hco%nj

    myLonBegHalo = 1 + weightLatLonStep * floor(real(myLonBeg - 1)/real(weightLatLonStep))
    myLonEndHalo = min(ni, 1 + weightLatLonStep * ceiling(real(myLonEnd - 1)/real(weightLatLonStep)))
    myLatBegHalo = 1 + weightLatLonStep * floor(real(myLatBeg - 1)/real(weightLatLonStep))
    myLatEndHalo = min(nj, 1 + weightLatLonStep * ceiling(real(myLatEnd - 1)/real(weightLatLonStep)))
    write(*,*) 'enkf_setupInterpInfo: myLonBeg/End, myLatBeg/End (original)  = ',  &
               myLonBeg, myLonEnd, myLatBeg, myLatEnd
    write(*,*) 'enkf_setupInterpInfo: myLonBeg/End, myLatBeg/End (with Halo) = ',  &
               myLonBegHalo, myLonEndHalo, myLatBegHalo, myLatEndHalo
    write(*,*) 'enkf_setupInterpInfo: myLonCount, myLatCount (with Halo) = ', &
               myLonEndHalo-myLonBegHalo+1, myLatEndHalo-myLatBegHalo+1
    write(*,*) 'enkf_setupInterpInfo: number of local gridpts where weights computed = ',  &
               ( 1 + ceiling(real(myLonEndHalo - myLonBegHalo) / real(weightLatLonStep)) ) *  &
               ( 1 + ceiling(real(myLatEndHalo - myLatBegHalo) / real(weightLatLonStep)) )
    call msg_memUsage('enkf_setupInterpInfo')

    wInterpInfo%latLonStep   = weightLatLonStep
    wInterpInfo%myLonBegHalo = myLonBegHalo
    wInterpInfo%myLonEndHalo = myLonEndHalo
    wInterpInfo%myLatBegHalo = myLatBegHalo
    wInterpInfo%myLatEndHalo = myLatEndHalo
    wInterpInfo%myLonBeg = myLonBeg
    wInterpInfo%myLonEnd = myLonEnd
    wInterpInfo%myLatBeg = myLatBeg
    wInterpInfo%myLatEnd = myLatEnd

    allocate(wInterpInfo%numIndexes(myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
    if (weightLatLonStep > 1) then
      ! Figure out if this tile straddles Yin-Yang boundary
      if (hco%grtyp == 'U' .and. myLatBegHalo <= nj/2 .and. myLatEndHalo >= ((nj/2)+1)) then
        includesYinYangBndry = .true.
      else
        includesYinYangBndry = .false.
      end if
      allocate(wInterpInfo%lonIndexes(4,myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
      allocate(wInterpInfo%latIndexes(4,myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
      allocate(wInterpInfo%interpWeights(4,myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
      wInterpInfo%lonIndexes(:,:,:) = 0
      wInterpInfo%latIndexes(:,:,:) = 0
      wInterpInfo%interpWeights(:,:,:) = 0.0D0
      ! Determine which lat-lon are interpolated (wInterpInfo%numIndexes>0)
      wInterpInfo%numIndexes(:,:) = 4
      do latIndex = myLatBegHalo, myLatEndHalo, weightLatLonStep
        do lonIndex = myLonBegHalo, myLonEndHalo, weightLatLonStep
          wInterpInfo%numIndexes(lonIndex,latIndex) = 0
        end do
      end do
      ! Ensure weights are computed along edge of domain
      if (myLonEndHalo == ni .and.  &
          myLatEndHalo == nj) then
        wInterpInfo%numIndexes(myLonEndHalo,myLatEndHalo) = 0
      end if
      if (myLonEndHalo == ni) then
        do latIndex = myLatBegHalo, myLatEndHalo, weightLatLonStep
          wInterpInfo%numIndexes(myLonEndHalo,latIndex) = 0
        end do
        ! Ensure weights are computed along both sides of Yin-Yang boundary
        if (includesYinYangBndry) then
          wInterpInfo%numIndexes(ni,nj/2) = 0
          wInterpInfo%numIndexes(ni,(nj/2)+1) = 0
          write(*,*) 'enkf_setupInterpInfo: Yin-Yang boundary (lon,lat1,lat2) =',  &
                     ni, nj/2, (nj/2)+1
        end if
      end if
      if (myLatEndHalo == nj) then
        do lonIndex = myLonBegHalo, myLonEndHalo, weightLatLonStep
          wInterpInfo%numIndexes(lonIndex,myLatEndHalo) = 0
        end do
      end if
      ! Ensure weights are computed along both sides of Yin-Yang boundary
      if (includesYinYangBndry) then
        do lonIndex = myLonBegHalo, myLonEndHalo, weightLatLonStep
          wInterpInfo%numIndexes(lonIndex,nj/2) = 0
          wInterpInfo%numIndexes(lonIndex,(nj/2)+1) = 0
          write(*,*) 'enkf_setupInterpInfo: Yin-Yang boundary (lon,lat1,lat2) =',  &
                     lonIndex, nj/2, (nj/2)+1
        end do
      end if
      ! For lon-only interpolation
      do latIndex = myLatBegHalo, myLatEndHalo
        if (wInterpInfo%numIndexes(myLonBegHalo,latIndex) > 0) cycle
        do lonIndex = myLonBegHalo, myLonEndHalo
          if (wInterpInfo%numIndexes(lonIndex,latIndex) == 0) cycle
          ! Find nearest grid point with a value towards left
          wInterpInfo%numIndexes(lonIndex,latIndex) = 2
          wInterpInfo%lonIndexes(1,lonIndex,latIndex) = myLonBegHalo +  &
               weightLatLonStep * floor(real(lonIndex - myLonBegHalo)/real(weightLatLonStep))
          wInterpInfo%lonIndexes(2,lonIndex,latIndex) = min(ni,  &
               wInterpInfo%lonIndexes(1,lonIndex,latIndex) + weightLatLonStep)
          wInterpInfo%latIndexes(1,lonIndex,latIndex) = latIndex
          wInterpInfo%latIndexes(2,lonIndex,latIndex) = latIndex
          wInterpInfo%interpWeights(1,lonIndex,latIndex) =   &
               real(wInterpInfo%lonIndexes(2,lonIndex,latIndex) - lonIndex, 8)/real(weightLatLonStep, 8)
          wInterpInfo%interpWeights(2,lonIndex,latIndex) = 1.0D0 -  &
               wInterpInfo%interpWeights(1,lonIndex,latIndex)
        end do
      end do
      ! For lat-only interpolation
      do latIndex = myLatBegHalo, myLatEndHalo
        do lonIndex = myLonBegHalo, myLonEndHalo, weightLatLonStep
          if (wInterpInfo%numIndexes(lonIndex,myLatBegHalo) > 0) cycle
          if (wInterpInfo%numIndexes(lonIndex,latIndex) == 0) cycle
          ! Find nearest grid point with a value towards bottom
          wInterpInfo%numIndexes(lonIndex,latIndex) = 2
          wInterpInfo%lonIndexes(1,lonIndex,latIndex) = lonIndex
          wInterpInfo%lonIndexes(2,lonIndex,latIndex) = lonIndex
          wInterpInfo%latIndexes(1,lonIndex,latIndex) = myLatBegHalo +  &
               weightLatLonStep * floor(real(latIndex - myLatBegHalo)/real(weightLatLonStep))
          wInterpInfo%latIndexes(2,lonIndex,latIndex) = min(nj,  &
               wInterpInfo%latIndexes(1,lonIndex,latIndex) + weightLatLonStep)
          ! Ensure we do not interpolate values across Yin-Yang boundary
          if (includesYinYangBndry) then
            if (latIndex <= nj/2) then
              wInterpInfo%latIndexes(2,lonIndex,latIndex) = min(nj/2, wInterpInfo%latIndexes(2,lonIndex,latIndex))
            else if(latIndex >= (nj/2)+1) then
              wInterpInfo%latIndexes(1,lonIndex,latIndex) = max((nj/2)+1, wInterpInfo%latIndexes(1,lonIndex,latIndex))
            end if
          end if
          wInterpInfo%interpWeights(1,lonIndex,latIndex) =  &
               real(wInterpInfo%latIndexes(2,lonIndex,latIndex) - latIndex, 8)/real(weightLatLonStep, 8)
          wInterpInfo%interpWeights(2,lonIndex,latIndex) = 1.0D0 -  &
               wInterpInfo%interpWeights(1,lonIndex,latIndex)
        end do
      end do
      ! For interior points needing 2D interpolation
      do latIndex = myLatBegHalo, myLatEndHalo
        do lonIndex = myLonBegHalo, myLonEndHalo
          if (wInterpInfo%numIndexes(lonIndex,latIndex) == 0) cycle ! no interpolation
          if (wInterpInfo%lonIndexes(1,lonIndex,latIndex) /= 0) cycle ! already set up
          wInterpInfo%numIndexes(lonIndex,latIndex) = 4
          ! 1. bottom-left indexes
          wInterpInfo%lonIndexes(1,lonIndex,latIndex) = myLonBegHalo +  &
               weightLatLonStep * floor(real(lonIndex - myLonBegHalo)/real(weightLatLonStep))
          wInterpInfo%latIndexes(1,lonIndex,latIndex) = myLatBegHalo +  &
               weightLatLonStep * floor(real(latIndex - myLatBegHalo)/real(weightLatLonStep))
          ! 2. bottom-right indexes
          wInterpInfo%lonIndexes(2,lonIndex,latIndex) = min(ni,  &
               wInterpInfo%lonIndexes(1,lonIndex,latIndex) + weightLatLonStep)
          wInterpInfo%latIndexes(2,lonIndex,latIndex) = wInterpInfo%latIndexes(1,lonIndex,latIndex)
          ! 3. top-left indexes
          wInterpInfo%lonIndexes(3,lonIndex,latIndex) = wInterpInfo%lonIndexes(1,lonIndex,latIndex)
          wInterpInfo%latIndexes(3,lonIndex,latIndex) = min(nj,  &
               wInterpInfo%latIndexes(1,lonIndex,latIndex) + weightLatLonStep)
          ! 4. top-right indexes
          wInterpInfo%lonIndexes(4,lonIndex,latIndex) = wInterpInfo%lonIndexes(2,lonIndex,latIndex)
          wInterpInfo%latIndexes(4,lonIndex,latIndex) = wInterpInfo%latIndexes(3,lonIndex,latIndex)
          ! Ensure we do not interpolate values across Yin-Yang boundary
          if (includesYinYangBndry) then
            if (latIndex <= nj/2) then
              wInterpInfo%latIndexes(3,lonIndex,latIndex) = min(nj/2, wInterpInfo%latIndexes(3,lonIndex,latIndex))
              wInterpInfo%latIndexes(4,lonIndex,latIndex) = min(nj/2, wInterpInfo%latIndexes(4,lonIndex,latIndex))
            else if(latIndex >= (nj/2)+1) then
              wInterpInfo%latIndexes(1,lonIndex,latIndex) = max((nj/2)+1, wInterpInfo%latIndexes(1,lonIndex,latIndex))
              wInterpInfo%latIndexes(2,lonIndex,latIndex) = max((nj/2)+1, wInterpInfo%latIndexes(2,lonIndex,latIndex))
            end if
          end if
          ! one-dimensional weights in lon and lat directions
          interpWeightLon = real(wInterpInfo%lonIndexes(4,lonIndex,latIndex) - lonIndex, 8) /  &
                            real(weightLatLonStep, 8)
          interpWeightLat = real(wInterpInfo%latIndexes(4,lonIndex,latIndex) - latIndex, 8) /  &
                            real(weightLatLonStep, 8)
          ! four interpolation weights
          wInterpInfo%interpWeights(1,lonIndex,latIndex) = interpWeightLon * interpWeightLat
          wInterpInfo%interpWeights(2,lonIndex,latIndex) = (1.0D0 - interpWeightLon) * interpWeightLat
          wInterpInfo%interpWeights(3,lonIndex,latIndex) = interpWeightLon * (1.0D0 - interpWeightLat)
          wInterpInfo%interpWeights(4,lonIndex,latIndex) = (1.0D0 - interpWeightLon) * (1.0D0 - interpWeightLat)
        end do
      end do
    else
      ! no interpolation, all weights are computed
      wInterpInfo%numIndexes(:,:) = 0
    end if
    call msg_memUsage('enkf_setupInterpInfo')

  end subroutine enkf_setupInterpInfo

  !--------------------------------------------------------------------------
  ! enkf_interpWeights
  !--------------------------------------------------------------------------
  subroutine enkf_interpWeights(wInterpInfo, weights)
    !
    ! :Purpose: Perform the bilinear interpolation of the weights
    !           using the precalculated interpolation info.
    !
    implicit none

    ! Arguments:
    type(struct_enkfInterpInfo), intent(in)  :: wInterpInfo  ! LETKF weight interpolation info
    real(8),                     intent(out) :: weights(1:,1:,wInterpInfo%myLonBegHalo:,wInterpInfo%myLatBegHalo:) ! Interpolated weights

    ! Locals:
    integer :: myLonBegHalo, myLonEndHalo, myLatBegHalo, myLatEndHalo
    integer :: lonIndex, latIndex, memberIndex1, memberIndex2, interpIndex
    integer :: interpLonIndex, interpLatIndex, numMembers1, numMembers2
    integer :: totalCount(mmpi_numthread)
    integer, external :: omp_get_thread_num
    logical, save :: firstCall = .true.

    myLonBegHalo = wInterpInfo%myLonBegHalo
    myLonEndHalo = wInterpInfo%myLonEndHalo
    myLatBegHalo = wInterpInfo%myLatBegHalo
    myLatEndHalo = wInterpInfo%myLatEndHalo
    numMembers1 = size(weights,1)
    numMembers2 = size(weights,2)
    totalCount(:) = 0

    !$OMP PARALLEL DO PRIVATE(latIndex, lonIndex, interpIndex, interpLatIndex, interpLonIndex, memberIndex1, memberIndex2)
    do latIndex = wInterpInfo%myLatBeg, wInterpInfo%myLatEnd
      do lonIndex = wInterpInfo%myLonBeg, wInterpInfo%myLonEnd
        if (wInterpInfo%numIndexes(lonIndex,latIndex) <= 0) cycle
        weights(:,:,lonIndex,latIndex) = 0.0D0
        if (wInterpInfo%lonIndexes(1,lonIndex,latIndex) == 0) cycle

        totalCount(omp_get_thread_num()+1) = totalCount(omp_get_thread_num()+1) + wInterpInfo%numIndexes(lonIndex,latIndex)

        do interpIndex = 1, wInterpInfo%numIndexes(lonIndex,latIndex)
          interpLonIndex = wInterpInfo%lonIndexes(interpIndex,lonIndex,latIndex)
          interpLatIndex = wInterpInfo%latIndexes(interpIndex,lonIndex,latIndex)

          do memberIndex2 = 1, numMembers2
            do memberIndex1 = 1, numMembers1
              weights(memberIndex1,memberIndex2,lonIndex,latIndex) =  &
                   weights(memberIndex1,memberIndex2,lonIndex,latIndex) + &
                   wInterpInfo%interpWeights(interpIndex,lonIndex,latIndex) *  &
                   weights(memberIndex1,memberIndex2,interpLonIndex,interpLatIndex)
            end do
          end do

        end do ! interpIndex

      end do ! lonIndex
    end do ! latIndex
    !$OMP END PARALLEL DO

    if (firstCall) write(*,*) 'enkf_interpWeights: totalCount = ', totalCount(:)
    firstCall = .false.

  end subroutine enkf_interpWeights

  !--------------------------------------------------------------------------
  ! enkf_modifyAMSUBobsError
  !--------------------------------------------------------------------------
  subroutine enkf_modifyAMSUBobsError(obsSpaceData)
    !
    !:Purpose: Ad-hoc modification of the observation error stddev of
    !          AMSUB (MHS, MWHS2) observations in the region equatorward
    !          of 40 degrees latitude. This was inherited from the
    !          original EnKF system.
    !
    implicit none

    ! Arguments:
    type(struct_obs), target, intent(inout) :: obsSpaceData ! Observation space information

    ! Locals:
    real(pre_obsReal), parameter :: AMSUB_trop_oer = 1.0 ! assumed value for AMSU-B obs error in tropics
    integer            :: headerIndex, bodyIndex, bodyIndexBeg, bodyIndexEnd, codeType
    real(pre_obsReal)  :: lat_obs

    ! for AMSUB observations set the observation error std dev equal to 1.0
    ! in the larger tropical area where the spread-skill correlation suggests
    ! that the data are accurate (.i.e |lat|<40. ). Otherwise don't reduce the
    ! observational error.
    do headerIndex = 1, obs_numheader(obsSpaceData)
      lat_obs = obs_headElem_r(obsSpaceData, obs_lat, headerIndex)
      codeType = obs_headElem_i(obsSpaceData, obs_ity, headerIndex)
      lat_obs = lat_obs * MPC_DEGREES_PER_RADIAN_R8
      if ( abs(lat_obs) < 40. .and. (codeType == codtyp_get_codtyp('amsub') .or.  &
                                     codeType == codtyp_get_codtyp('mhs') .or.  &
                                     codeType == codtyp_get_codtyp('mwhs2')) ) then
        bodyIndexBeg = obs_headElem_i(obsSpaceData, obs_rln, headerIndex)
        bodyIndexEnd = obs_headElem_i(obsSpaceData, obs_nlv, headerIndex) + bodyIndexBeg - 1
        do bodyIndex = bodyIndexBeg, bodyIndexEnd
          call obs_bodySet_r(obsSpaceData, obs_oer, bodyIndex, AMSUB_trop_oer)
        end do
      end if
    end do

  end subroutine enkf_modifyAMSUBobsError

  !--------------------------------------------------------------------------
  ! enkf_rejectHighLatIR
  !--------------------------------------------------------------------------
  subroutine enkf_rejectHighLatIR(obsSpaceData)
    !
    !:Purpose: Reject hyperspectral IR radiance observations at latitudes
    !          poleward of 60 degrees latitude.
    !
    implicit none

    ! Arguments:
    type(struct_obs), target, intent(inout) :: obsSpaceData ! Observation space information

    ! Locals:
    integer           :: headerIndex, bodyIndex, bodyIndexBeg, bodyIndexEnd, codeType
    real(pre_obsReal) :: lat_obs

    ! reject all HIR radiance observation in arctic and antarctic (.i.e |lat|>60. )
    do headerIndex = 1, obs_numheader(obsSpaceData)
      lat_obs = obs_headElem_r(obsSpaceData, obs_lat, headerIndex)
      codeType = obs_headElem_i(obsSpaceData, obs_ity, headerIndex)
      lat_obs = lat_obs * MPC_DEGREES_PER_RADIAN_R8
      if ( abs(lat_obs) > 60. .and. tvs_isIdBurpHyperSpectral(codeType) ) then
        write(*,*) 'enkf_rejectHighLatIR: !!!!!!!!--------WARNING--------!!!!!!!!'
        write(*,*) 'enkf_rejectHighLatIR: This HIR radiance profile was rejected because |lat|>60.'
        write(*,*) 'enkf_rejectHighLatIR: latidude= ', lat_obs, 'codtyp= ', codeType
        bodyIndexBeg = obs_headElem_i(obsSpaceData, obs_rln, headerIndex)
        bodyIndexEnd = obs_headElem_i(obsSpaceData, obs_nlv, headerIndex) + bodyIndexBeg - 1
        do bodyIndex = bodyIndexBeg, bodyIndexEnd
          call obs_bodySet_i(obsSpaceData, obs_ass, bodyIndex, obs_notAssimilated)
          ! also set the 'rejected by selection process' flag (bit 11)
          call obs_bodySet_i( obsSpaceData, obs_flg, bodyIndex,  &
                              ibset( obs_bodyElem_i( obsSpaceData, obs_flg, bodyIndex ), 11) )
        end do
      end if
    end do

  end subroutine enkf_rejectHighLatIR

  !--------------------------------------------------------------------------
  ! enkf_getModulatedState
  !--------------------------------------------------------------------------
  subroutine enkf_getModulatedState( enkfNML, stateVector_in, stateVectorMeanTrl, &
                                     eigenVectorColumnIndex, stateVector_out, &
                                     beSilent )
    !
    !:Purpose: Compute vertical localization matrix, and the corresponding
    !          eigenvectors/eigenvalues, to obtain modulated stateVector.
    !
    implicit none

    ! Arguments:
    type(struct_enkfNML), intent(in)    :: enkfNML                ! Derived type variable with namelist variables
    type(struct_gsv),     intent(in)    :: stateVector_in         ! State vector with ens member to be processed
    type(struct_gsv),     intent(in)    :: stateVectorMeanTrl     ! Ensemble mean state vector
    integer,              intent(in)    :: eigenVectorColumnIndex ! Index of eigenvector to use
    type(struct_gsv),     intent(inout) :: stateVector_out        ! Resulting modulated ensemble member
    logical,              intent(in)    :: beSilent               ! Control of verbose output

    ! Locals:
    real(4)          :: modulationFactor_r4
    real(4), pointer :: field_out_r4(:,:,:,:)
    integer :: nLev, nlev_out, levIndex, latIndex, lonIndex
    integer :: lon1, lon2, lat1, lat2
    integer :: varIndex, stepIndex, eigenVectorLevelIndex
    character(len=4) :: varName

    call utl_tmg_start(130,'--getModulatedState')

    if ( .not. beSilent ) write(*,*) 'enkf_getModulatedState: START'

    if ( stateVector_in%dataKind /= 4 ) then
      call utl_abort('enkf_getModulatedState: only dataKind=4 is implemented')
    end if

    nLev = stateVector_in%vco%nLev_M
    if ( enkfNML%vLocalize <= 0.0d0 .or. nLev <= 1 ) then
      call utl_abort('enkf_getModulatedState: no vertical localization')
    end if

    ! Compute perturbation by subtracting ensMean
    call gsv_copy(stateVector_in, stateVector_out, beSilent_opt=beSilent)
    call gsv_add(stateVectorMeanTrl, stateVector_out, scaleFactor_opt=-1.0d0)

    lon1 = stateVector_out%myLonBeg
    lon2 = stateVector_out%myLonEnd
    lat1 = stateVector_out%myLatBeg
    lat2 = stateVector_out%myLatEnd

    ! Compute modulated member perturbation from original member perturbation:
    !   v'_k = (Nens*nLamda/(Nens - 1))^1/2 * Lambda^1/2 * E * x'_k
    step_loop: do stepIndex = 1, stateVector_out%numStep
      var_loop: do varIndex = 1, vnl_numvarmax
        varName = vnl_varNameList(varIndex)
        if ( .not. gsv_varExist(stateVector_out,varName) ) cycle var_loop

        nlev_out  = stateVector_out%varNumLev(varIndex)

        call gsv_getField(statevector_out,field_out_r4,varName)

        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2
            do levIndex = 1, nlev_out
              if ( nlev_out == 1 ) then
                eigenVectorLevelIndex = nLev
              else
                eigenVectorLevelIndex = levIndex
              end if

              call getModulationFactor( enkfNML, stateVector_in%vco, eigenVectorLevelIndex, &
                                        eigenVectorColumnIndex, &
                                        modulationFactor_r4, beSilent_opt=beSilent )

              field_out_r4(lonIndex,latIndex,levIndex,stepIndex) = &
                                 field_out_r4(lonIndex,latIndex,levIndex,stepIndex) * &
                                 modulationFactor_r4
            end do
          end do
        end do

      end do var_loop
    end do step_loop

    ! Now add to ensMean to get modulated member
    ! v_k = v'_k + v_mean
    call gsv_add(stateVectorMeanTrl, stateVector_out)

    if ( .not. beSilent ) write(*,*) 'enkf_getModulatedState: END'

    call utl_tmg_stop(130)

  end subroutine enkf_getModulatedState

  !--------------------------------------------------------------------------
  ! enkf_setupModulationFactor
  !--------------------------------------------------------------------------
  subroutine enkf_setupModulationFactor(enkfNML, vco, beSilent)
    !
    !:Purpose: Setup modulationFactorArray by calling getModulationFactor for first time.
    !
    implicit none

    ! Arguments:
    type(struct_enkfNML),      intent(in) :: enkfNML  ! Derived type variable with namelist variables
    type(struct_vco), pointer, intent(in) :: vco      ! Vertical coordinate definition
    logical,                   intent(in) :: beSilent ! Control of verbose output

    ! Locals:
    integer :: eigenVectorColumnIndex
    integer :: eigenVectorLevelIndex
    real(4) :: modulationFactor_r4

    eigenVectorColumnIndex = 1
    eigenVectorLevelIndex = 1
    call getModulationFactor(enkfNML, vco, eigenVectorLevelIndex, &
                             eigenVectorColumnIndex, &
                             modulationFactor_r4, beSilent_opt=beSilent)

  end subroutine enkf_setupModulationFactor

  !--------------------------------------------------------------------------
  ! getModulationFactor
  !--------------------------------------------------------------------------
  subroutine getModulationFactor( enkfNML, vco, eigenVectorLevelIndex, &
                                  eigenVectorColumnIndex, &
                                  modulationFactor_r4, beSilent_opt )
    !
    !:Purpose: Compute modulation factor needed to multiply ensemble
    !          perturbation to get the modulated perturbation:
    !          (Nens*nLambda/(Nens - 1))^1/2 * Lambda^1/2
    !
    implicit none

    ! Arguments:
    type(struct_enkfNML),      intent(in)  :: enkfNML                ! Derived type variable with namelist variables
    type(struct_vco), pointer, intent(in)  :: vco                    ! Vertical coordinate definition
    integer,                   intent(in)  :: eigenVectorLevelIndex  ! Index of vertical level to use
    integer,                   intent(in)  :: eigenVectorColumnIndex ! Index of eigenvector to use
    real(4),                   intent(out) :: modulationFactor_r4    ! Resulting modulation factor
    logical, optional,         intent(in)  :: beSilent_opt           ! Control of verbose output

    ! Locals:
    integer             :: levIndex1, levIndex2, eigenIndex
    integer             :: nLev, nLev_M, nLev_depth, matrixRank
    real(8)             :: zr, zcorr, pSurfRef, hSurfRef
    real(8)             :: tolerance
    real(8), pointer    :: pressureProfile(:)
    real(8), allocatable, save :: eigenValues(:)
    real(8), allocatable, save :: eigenVectors(:,:)
    real(8), allocatable, save :: verticalLocalizationMat(:,:)
    real(8), allocatable, save :: verticalLocalizationMatLowRank(:,:)
    real(4), allocatable, save :: modulationFactorArray_r4(:,:)
    logical :: beSilent

    logical, save :: firstCall = .true.

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .false.
    end if

    ! Compute vertical localization matrix and its eigenValues/Vectors on first call
    if ( firstCall ) then
      firstCall = .false.
      if ( mmpi_myid == 0 .and. .not. beSilent ) then
        write(*,*) 'getModulationFactor: computing eigenValues/Vectors'
      end if

      nLev_M = vco%nLev_M
      nLev_depth = vco%nlev_depth
      nLev = max(nLev_M,nLev_depth)

      allocate(eigenValues(nLev))
      allocate(eigenVectors(nLev,nLev))
      allocate(verticalLocalizationMat(nLev,nLev))
      allocate(verticalLocalizationMatLowRank(nLev,nLev))
      allocate(modulationFactorArray_r4(enkfNML%numRetainedEigen,nLev))
      verticalLocalizationMatLowRank(:,:) = 0.0d0

      nullify(pressureProfile)
      if (vco%Vcode == 21001) then
        hSurfRef = 0.D0
        call czp_fetch1DLevels(vco, sfcValue=hSurfRef, sfcValueLS_opt=hSurfRef, &
                               profM_opt=pressureProfile)
        call czp_calcPressureProfileUsingStdAtm(pressureProfile,               & ! INOUT
                                                nLev_M)   ! IN
      else if(vco%Vcode == 5002 .or. vco%Vcode == 5005 .or. vco%Vcode == 5100) then
        pSurfRef = 101000.D0
        call czp_fetch1DLevels(vco, pSurfRef, sfcValueLS_opt=pSurfRef, &
                               profM_opt=pressureProfile)
      else
        call utl_abort('getModulationFactor: Unknown value of vcode')
      end if

      call lfn_Setup(LocFunctionWanted='FifthOrder')

      ! Calculate 5'th order function
      do levIndex1 = 1, nLev
        do levIndex2 = 1, nLev
          zr = abs(log(pressureProfile(levIndex2)) - log(pressureProfile(levIndex1)))
          zcorr = lfn_response(zr,enkfNML%vLocalize)
          verticalLocalizationMat(levIndex1,levIndex2) = zcorr
        end do
      end do

      ! Compute eigenValues/Vectors of vertical localization matrix
      tolerance = 1.0D-50
      call utl_eigenDecomp(verticalLocalizationMat, eigenValues, eigenVectors, &
                           tolerance, matrixRank)
      if ( matrixRank < enkfNML%numRetainedEigen ) then
        write(*,*) 'matrixRank=', matrixRank
        call utl_abort('getModulationFactor: verticalLocalizationMat is rank deficient=')
      end if

      ! Compute low-ranked vertical localization matrix
      do levIndex1 = 1, nLev
        do levIndex2 = 1, nLev
          do eigenIndex = 1, enkfNML%numRetainedEigen
            verticalLocalizationMatLowRank(levIndex1,levIndex2) = verticalLocalizationMatLowRank(levIndex1,levIndex2) + &
                                                                  eigenVectors(levIndex1,eigenIndex) * &
                                                                  eigenVectors(levIndex2,eigenIndex) * &
                                                                  eigenValues(eigenIndex)
          end do
        end do
      end do

      ! now compute the 2D modulationFactor array
      do levIndex1 = 1, nLev
        do eigenIndex = 1, enkfNML%numRetainedEigen
          modulationFactorArray_r4(eigenIndex,levIndex1) = real( &
                        1 / sqrt(verticalLocalizationMatLowRank(levIndex1,levIndex1)) * &
                        eigenVectors(levIndex1,eigenIndex) * &
                        eigenValues(eigenIndex) ** 0.5 * &
                        (enkfNML%nEns * enkfNML%numRetainedEigen / (enkfNML%nEns - 1)) ** 0.5,4)
        end do
      end do

      if ( mmpi_myid == 0 .and. .not. beSilent ) then
        do levIndex1 = 1, enkfNML%numRetainedEigen
          write(*,*) 'getModulationFactor: eigen mode=', levIndex1, ', eigenVectors=', eigenVectors(:,levIndex1)
        end do
        write(*,*) 'getModulationFactor: eigenValues=', eigenValues(1:enkfNML%numRetainedEigen)

        do levIndex1 = 1, nLev
          write(*,*) 'getModulationFactor: verticalLocalizationMat for lev ', levIndex1, '=', verticalLocalizationMat(levIndex1,:)
          write(*,*) 'getModulationFactor: verticalLocalizationMatLowRank for lev ', levIndex1, '=', verticalLocalizationMatLowRank(levIndex1,:)
        end do
      end if
    end if

    modulationFactor_r4 = modulationFactorArray_r4(eigenVectorColumnIndex,eigenVectorLevelIndex)

  end subroutine getModulationFactor

  !--------------------------------------------------------------------------
  ! enkf_getLocalizationRadius
  !--------------------------------------------------------------------------
  subroutine enkf_getLocalizationRadius(hLocalize, hLocalizePressure, &
                                        anlVertLocation, hLinearLoc, hLoc)
    !
    !:Purpose: get the localization radius, interpolated or not, at a given pressure
    !
    implicit none

    ! Arguments:
    real(8), intent(in)  :: hLocalize(:)         ! the list of localization radii (m)
    real(8), intent(in)  :: hLocalizePressure(:) ! the pressures where the radius changes (log(P))
    real(8), intent(in)  :: anlVertLocation      ! the gridpoint vertical coordinate in log(P)
    logical, intent(in)  :: hLinearLoc           ! apply linear vertical interpolation for the localization radius
    real(8), intent(out) :: hLoc                 ! the gridpoint localization radius

    ! Locals:
    integer              :: hLocIndex, numPresValues

    numPresValues = count(hLocalizePressure > 0.0d0)

    if ( hLocalize(1) < 0.0d0 ) then
      call utl_abort('enkf_getLocalizationRadius: hLocalize(1) < 0.0d0')
    end if

    ! radius is constant
    if ( all(hLocalize(:) == hLocalize(1)) ) then
      hLoc = hLocalize(1)

    ! radius varies vertically, and is linearly interpolated with log(P)
    else if (hLinearLoc) then
      hLocIndex = 1 + count(anlVertLocation >= hLocalizePressure(1:numPresValues))
      ! constant radius value near the top of the atmosphere
      if (hLocIndex == 1) then
        hLoc = hLocalize(1)
      ! constant radius value near the bottom of the atmosphere
      else if (hLocIndex == numPresValues+1) then
        hLoc = hLocalize(numPresValues)
      ! piece-wise linear interpolation
      else
        hLoc = hLocalize(hLocIndex-1) + &
               (anlVertLocation - hLocalizePressure(hLocIndex-1)) * &
               (hLocalize(hLocIndex) - hLocalize(hLocIndex-1))  / &
               (hLocalizePressure(hLocIndex) - hLocalizePressure(hLocIndex-1))
      end if

    ! radius varies vertically, but is not interpolated
    else
      hLocIndex = 1 + count(anlVertLocation > hLocalizePressure(1:numPresValues))
      hLoc = hLocalize(hLocIndex)

    end if

  end subroutine enkf_getLocalizationRadius

  !--------------------------------------------------------------------------
  ! enkf_writeEdim
  !--------------------------------------------------------------------------
  subroutine enkf_writeEdim(eigenValues,nEns)
    !
    !:Purpose: Compute and output ensemble dimension, degrees of freedom,
    !          and trace of Pa from the eigenvalues of Pa
    implicit none

    ! Arguments:
    integer, intent(in) :: nEns              ! number of members
    real(8), intent(in) :: eigenValues(nEns) ! eigenvalues of Pa in ensemble space

    ! Locals:
    real(8)             :: eDim, dof, trace
    character(len=50)   :: outfilename
    integer             :: memberIndex
    integer             :: fclos, funit, ierr

    eDim  = 0.0d0
    dof   = 0.0d0
    trace = 0.0d0
    do memberIndex = 1, nEns
      eDim  = eDim  + 1./sqrt((eigenValues(memberIndex)+nEns-1))
      dof   = dof   + 1./(1.+eigenValues(memberIndex)+nEns-1)
      trace = trace + 1./(eigenValues(memberIndex)+nEns-1)
    end do
    eDim = eDim**2/trace
    write(outfilename, '(I5.5)') mmpi_myid ! we assume there are less than 100 000 mpi tasks...
    outfilename = './eob_glbi_'//trim(adjustl(outfilename))
    call utl_open_asciifile(outfilename,funit)
    write(funit,'(A,2(1X,ES11.4),1X,ES13.6)') 'edim', eDim, dof, trace
    ierr = fclos(funit)

  end subroutine enkf_writeEdim

  !--------------------------------------------------------------------------
  ! enkf_allocateDFS
  !--------------------------------------------------------------------------
  subroutine enkf_allocateDFS(enkfDFS, numLev, numLatLon, maxNumLocalObs)
    !
    !:Purpose: allocate and initialize a DFS output structure
    !
    implicit none

    ! Arguments:
    type(struct_enkfDFS), intent(inout) :: enkfDFS        ! DFS structure
    integer,              intent(in)    :: numLev         ! number of vertical levels
    integer,              intent(in)    :: numLatLon      ! number of grid points
    integer,              intent(in)    :: maxNumLocalObs ! maximum number of obs in local volume

    ! Locals:
    integer :: numDFSIndex ! grid index

    numDFSIndex = numLev * numLatLon

    write(*,*) 'enkf_allocateDFS: enkfDFS dimensions:', numDFSIndex, maxNumLocalObs
    allocate(enkfDFS%locFun(numDFSIndex,maxNumLocalObs))
    allocate(enkfDFS%bodyIndex(numDFSIndex,maxNumLocalObs))
    allocate(enkfDFS%lat(numDFSIndex))
    allocate(enkfDFS%lon(numDFSIndex))
    allocate(enkfDFS%lnp(numDFSIndex))
    allocate(enkfDFS%dfs(numDFSIndex))
    write(*,*) 'enkf_allocateDFS: enkfDFS allocated'

    enkfDFS%allocated = .true.
    enkfDFS%locFun(:,:) = 0.0d0
    enkfDFS%bodyIndex(:,:) = 0
    enkfDFS%lat(:) = 0.0d0
    enkfDFS%lon(:) = 0.0d0
    enkfDFS%lnp(:) = 0.0d0
    enkfDFS%dfs(:) = 0.0d0

  end subroutine enkf_allocateDFS

  !--------------------------------------------------------------------------
  ! enkf_deallocateDFS
  !--------------------------------------------------------------------------
  subroutine enkf_deallocateDFS(enkfDFS)
    !
    !:Purpose: deallocate a DFS output structure
    !
    implicit none

    ! Arguments:
    type(struct_enkfDFS), intent(inout) :: enkfDFS ! DFS structure

    deallocate(enkfDFS%dfs)
    deallocate(enkfDFS%lnp)
    deallocate(enkfDFS%lon)
    deallocate(enkfDFS%lat)
    deallocate(enkfDFS%bodyIndex)
    deallocate(enkfDFS%locFun)
    enkfDFS%allocated = .false.

  end subroutine enkf_deallocateDFS

end module enkf_mod

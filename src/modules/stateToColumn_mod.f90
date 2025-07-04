
module stateToColumn_mod
  ! MODULE stateToColumn_mod (prefix='s2c' category='4. Data Object transformations')
  !
  !:Purpose:  Non-linear, tangent-linear and adjoint versions of
  !           horizontal-temporal interpolation between a gridStateVector object
  !           and a columnData object.
  !
  use mpi_f08 ! this is the Fortran 2008 MPI library module
  use midasMpi_mod
  use mathPhysConstants_mod
  use earthConstants_mod
  use message_mod
  use codePrecision_mod
  use gridStateVector_mod
  use obsSpaceData_mod
  use columnData_mod
  use horizontalCoord_mod
  use obsTimeInterp_mod
  use windRotation_mod
  use utilities_mod
  use gridVariableTransforms_mod
  use varNameList_mod
  use slantProfileLatLon_mod
  use tovs_mod
  use codtyp_mod
  use getGridPosition_mod
  use kdTree2_mod
  use calcHeightAndPressure_mod
  use humidityLimits_mod
  use obsFlags_mod

  implicit none
  save
  private

  ! Public routines
  public :: s2c_tl, s2c_ad, s2c_nl
  public :: s2c_bgcheck_bilin, s2c_getFootprintRadius, s2c_getWeightsAndGridPointIndexes
  public :: s2c_deallocInterpInfo

  ! private module variables and derived types

  ! Basic real48 types that contain allocatable arrays of both real(4) and real(8) types
  type real48_2d
    integer              :: dataKind = 0
    real(4), allocatable :: r4(:,:)
    real(8), allocatable :: r8(:,:)
    real(8), allocatable :: GZsfc(:)
  end type real48_2d
  type real48_3d
    integer              :: dataKind = 0
    real(4), allocatable :: r4(:,:,:)
    real(8), allocatable :: r8(:,:,:)
    real(8), allocatable :: GZsfc(:)
  end type real48_3d

  ! Secondary derived type for using the "2DFIELDS" MPI strategy
  type struct_stepProcData
    ! lat-lon location of observations to be interpolated
    real(8), pointer          :: allLat(:,:) => null()         ! (headerUsed, varLevIndex)
    real(8), pointer          :: allLon(:,:) => null()         ! (headerUsed, varLevIndex)
    ! lat-lon location on rotated grid of observations to be interpolated
    real(8), pointer          :: allLatRot(:,:,:) => null()    ! (subGrid, headerUsed, varLevIndex)
    real(8), pointer          :: allLonRot(:,:,:) => null()    ! (subGrid, headerUsed, varLevIndex)
    ! actual headerIndex, since the headerUsed is only for those obs with a non-zero interp weight
    integer, pointer          :: allHeaderIndex(:) => null()   ! (headerUsed)
    ! depotIndexBeg/End contain first/last indices into depots of interpolation weights and lat/lon indices
    integer, pointer          :: depotIndexBeg(:,:,:) => null() ! (subGrid, headerUsed, varLevIndex)
    integer, pointer          :: depotIndexEnd(:,:,:) => null() ! (subGrid, headerUsed, varLevIndex)
  end type struct_stepProcData

  ! Main derived type for using the "2DFIELDS" MPI strategy
  type struct_interpInfo
    logical                   :: initialized = .false.
    type(struct_hco), pointer :: hco => null() ! horizontal grid object
    type(struct_uvr), pointer :: uvr => null() ! windRotation object
    type(struct_oti), pointer :: oti => null() ! obsTimeInterp object

    ! number of obs headers on each proc having a non-zero interp weight for each stepIndex (headerUsed)
    integer, pointer          :: allNumHeaderUsed(:,:) => null()    ! (step, proc)

    ! structure containing all interpolation information that depends on (proc, step)
    type(struct_stepProcData), allocatable :: stepProcData(:,:) ! (proc, step)

    ! interpolation weights and lat/lon indices are accessed via the 'stepProcData%depotIndexBeg/End'
    real(8), allocatable      :: interpWeightDepot(:)                ! (depotIndex)
    integer, pointer          :: latIndexDepot(:)                    ! (depotIndex)
    integer, pointer          :: lonIndexDepot(:)                    ! (depotIndex)
    character(len=2)          :: inputStateVectorType
  end type struct_interpInfo

  ! This parameter is the value of latIndexDepot that indicates the value is
  ! located in the halo and the lonIndexDepot value indicates the value of
  ! haloIndex for this location.
  integer, parameter :: valueIsInHalo = -99

  ! Derived type for using the "TILES" MPI strategy
  type struct_interpInfoTiles
    logical                       :: initialized = .false.
    logical                       :: rotatedWinds = .false.
    integer                       :: numVarLevState
    integer, allocatable          :: varLevColFromVarLevState(:)
    integer                       :: sfcVarLevIndex
    character(len=2)              :: inputStateVectorType
    type(struct_hco), pointer     :: hco => null() ! horizontal grid object
    type(struct_uvr), pointer     :: uvr => null() ! windRotation object
    type(struct_oti), pointer     :: oti => null() ! obsTimeInterp object

    ! Number of obs headers on each proc
    integer, allocatable      :: allNumHeader(:)      ! (proc)

    ! MPI proc id for all observations wrt lat-lon tiles
    integer, allocatable      :: allObsTileMpiId(:,:) ! (headerIndex, proc)

    ! MPI proc id for all observations wrt original obsSpaceData
    integer              :: yourNumHeader
    integer, allocatable :: yourObsTileMpiId(:)      ! (yourHeaderIndex)
    integer, allocatable :: yourObsSubGridIndex(:,:) ! (varLevIndex,yourHeaderIndex)
    real(8), allocatable :: yourObsLat(:,:)          ! (varLevIndex,yourHeaderIndex)
    real(8), allocatable :: yourObsLon(:,:)          ! (varLevIndex,yourHeaderIndex)

    ! Information about obs on my Tile that I am responsible for interpolating
    integer              :: myInterpNumHeader            ! Number columns/headers in my tile
    real(8), allocatable :: myInterpObsLat(:,:)          ! (varLevIndex,myHeaderIndex)
    real(8), allocatable :: myInterpObsLon(:,:)          ! (varLevIndex,myHeaderIndex)
    real(8), allocatable :: myInterpObsLatRot(:,:)       ! (varLevIndex,myHeaderIndex)
    real(8), allocatable :: myInterpObsLonRot(:,:)       ! (varLevIndex,myHeaderIndex)
    real(4), allocatable :: myInterpObsXpos_r4(:,:)      ! (varLevIndex,myHeaderIndex)
    real(4), allocatable :: myInterpObsYpos_r4(:,:)      ! (varLevIndex,myHeaderIndex)
    integer, allocatable :: myInterpObsSubGridIndex(:,:) ! (varLevIndex,myHeaderIndex)
    real(4), allocatable :: myInterpObsFootprint_r4(:,:) ! (varLevIndex,myHeaderIndex)
    integer, allocatable :: myInterpObsMpiIdSrc(:)       ! (myHeaderIndex)
    integer, allocatable :: myInterpObsHeaderIndex(:)    ! (myHeaderIndex)

    ! Interpolation weights and grid point indexes, stored within a "depot"
    integer, allocatable :: depotIndexBeg(:,:)     ! (varLevIndex,myHeaderIndex)
    integer, allocatable :: depotIndexEnd(:,:)     ! (varLevIndex,myHeaderIndex)
    integer, allocatable :: latIndexDepot(:)       ! (depotIndex)
    integer, allocatable :: lonIndexDepot(:)       ! (depotIndex)
    real(8), allocatable :: interpWeightDepot(:)   ! (depotIndex)

    ! Halo information
    logical              :: periodic
    integer              :: myHaloSize
    integer              :: minLon                ! needed for dimensioning
    integer              :: minLat
    integer, allocatable :: myHaloLatIndex(:)     ! (haloIndex)
    integer, allocatable :: myHaloLonIndex(:)     ! (haloIndex)
    integer, allocatable :: myHaloMpiIdSrc(:)     ! (haloIndex)
    integer, allocatable :: myHaloMpiTag(:)       ! (haloIndex)
    integer              :: yourHaloSize
    integer, allocatable :: yourHaloLatIndex(:)   ! (haloIndex)
    integer, allocatable :: yourHaloLonIndex(:)   ! (haloIndex)
    integer, allocatable :: yourHaloMpiIdDst(:)   ! (haloIndex)
    integer, allocatable :: yourHaloMpiTag(:)     ! (haloIndex)

  end type struct_interpInfoTiles

  type(struct_interpInfo),      target :: interpInfo_tlad, interpInfo_nl
  type(struct_interpInfoTiles), target :: interpInfoTiles_tlad, interpInfoTiles_nl
  type(kdtree2), pointer  :: tree_nl => null()
  type(kdtree2), pointer  :: tree_tlad => null()

  character(len=20), parameter :: timeInterpType_tlad = 'LINEAR' ! hardcoded type of time interpolation for increment

  integer, parameter :: maxNumWrites = 50
  logical, parameter :: verbose = .false.

  ! "special" values of the footprint radius
  real(4), parameter :: nearestNeighbourFootprint = -3.0
  real(4), parameter ::             lakeFootprint = -2.0
  real(4), parameter ::         bilinearFootprint = -1.0
  integer, parameter :: maxNumLocalGridptsSearch = 1000
  integer, parameter :: minNumLocalGridptsSearch = 8

  ! namelist variables
  character(len=20) :: mpiMode             ! select mode for MPI distribution: 'TILES','2DFIELDS'
  logical :: slantPath_TO_nl               ! choose to use slant path for non-linear radiance operator
  logical :: slantPath_TO_tlad             ! choose to use slant path for linearized radiance operators
  logical :: slantPath_RO_nl               ! choose to use slant path for non-linear GPS-RO operator
  logical :: slantPath_RA_nl               ! choose to use slant path for non-linear radar operator
  logical :: calcHeightPressIncrOnColumn   ! choose to compute Z and P increment on column, instead of grid
  logical :: useFootprintForTovs           ! choose to use a horizontal footprint operator for radiance obs
  logical :: rejectObsNonMonotonicPressure ! choose to reject obs when interpolated column pressure is non-monotonic
  logical :: rejectObsOutsideGlobalGrid    ! choose to reject obs outside a global domain, currently employed for ORCA025 global grid
  logical :: NNInterpForCloudVars          ! to perform nearest neighbour horizontal interpolation for cloudy variables
  logical :: NNInterpForAllVars            ! to perform nearest neighbour horizontal interpolation for selected variablles

contains

  !---------------------------------------------------------
  ! readNml
  !---------------------------------------------------------
  subroutine readNml()
    !
    ! :Purpose: Read the namelist
    !
    implicit none

    ! Locals:
    integer :: ierr
    logical, save :: nmlAlreadyRead = .false.

    namelist /nams2c/ mpiMode, slantPath_TO_nl, slantPath_TO_tlad, slantPath_RO_nl, slantPath_RA_nl
    namelist /nams2c/ useFootprintForTovs, rejectObsNonMonotonicPressure, rejectObsOutsideGlobalGrid
    namelist /nams2c/ calcHeightPressIncrOnColumn
    namelist /nams2c/ NNInterpForCloudVars, NNInterpForAllVars

    if (nmlAlreadyRead) return

    write(*,*) 'readNml (s2c): STARTING'

    nmlAlreadyRead = .true.
    write(*,*) 'readNml (s2c): Reading namelist'

    ! default values
    mpiMode           = '2DFIELDS'
    slantPath_TO_nl   = .false.
    slantPath_TO_tlad = .false.
    slantPath_RO_nl   = .false.
    slantPath_RA_nl   = .false.
    calcHeightPressIncrOnColumn = .false.
    useFootprintForTovs = .false.
    rejectObsNonMonotonicPressure =.true.
    rejectObsOutsideGlobalGrid = .false.
    NNInterpForCloudVars = .false.
    NNInterpForAllVars = .false.

    if (.not. utl_isNamelistPresent('NAMS2C','./flnml') ) then

      if ( mmpi_myid == 0 ) then
        write(*,*) 'readNml (s2c): nams2c is missing in the namelist.'
        write(*,*) '               The default values will be taken.'
      end if

    else

      ! reading namelist variables
      call utl_tmg_start(181,'low-level--readNML')
      read(utl_flnml, nml = nams2c, iostat = ierr)
      if ( ierr /= 0 ) call utl_abort('readNml (s2c): Error reading namelist')
      call utl_tmg_stop(181)

    end if

    if (mmpi_myid == 0) write(*, nml = nams2c)

    ! Check for invalid combinations of namelist variables

    if (trim(mpiMode) /= '2DFIELDS' .and. trim(mpiMode) /= 'TILES') then
      call utl_abort('readNml (s2c): invalid value of mpiMode = '//trim(mpiMode))
    end if

    write(*,*) 'readNml (s2c): FINISHED'

  end subroutine readNml

  !-------------------------------------------------------------
  ! The following subroutines are for the newer "TILES" mpiMode
  !-------------------------------------------------------------

  !---------------------------------------------------------
  ! setupInterpInfoTiles
  !---------------------------------------------------------
  subroutine setupInterpInfoTiles(intInfo, obsSpaceData, stateVector,  &
                                  column, timeInterpType, rejectOutsideObs, &
                                  inputStateVectorType, beSilent)
    !
    ! :Purpose: Setup all of the information needed to quickly
    !           perform the horizontal interpolation to the observation
    !           locations using the mpi strategy of keeping most of
    !           the gridded data as lat-lon tiles to avoid lots of
    !           communication.
    !
    !           Note: this first version is simplified by not supporting
    !           multiple obs batches. Also, may still need to introduce
    !           some latlon checks.
    !
    implicit none

    ! Arguments
    type(struct_interpInfoTiles), intent(out)   :: intInfo              ! Interpolation info structure
    type(struct_obs)            , intent(inout) :: obsSpaceData         ! obs space object
    type(struct_gsv), target    , intent(in)    :: stateVector          ! stateVector object
    type(struct_columnData)     , intent(in)    :: column               ! column object
    logical                     , intent(inout) :: rejectOutsideObs     ! choose to reject obs outside domain
    character(len=*)            , intent(in)    :: timeInterpType       ! type of temporal interpolation
    character(len=*)            , intent(in)    :: inputStateVectorType ! type of state vector "nl", "tl" or "ad"
    logical,                      intent(in)    :: beSilent             ! choose to print nothing

    ! Locals:
    integer :: numStep, numHeader, varLevIndex
    character(len=4) :: varName

    if (.not. beSilent) then
      write(*,*) 'setupInterpInfoTiles (s2c): STARTING'
      call msg_memUsage('setupInterpInfoTiles')
    end if

    intInfo%initialized = .true.
    intInfo%numVarLevState = stateVector%numVarLev
    intInfo%yourNumHeader = obs_numHeader(obsSpaceData)
    intInfo%inputStateVectorType = trim(inputStateVectorType)
    allocate(intInfo%varLevColFromVarLevState(intInfo%numVarLevState))
    do varLevIndex = 1, intInfo%numVarLevState
      varName = gsv_getVarNameFromVarLev(statevector,varLevIndex)
      intInfo%varLevColFromVarLevState(varLevIndex) =  &
           col_getOffsetFromVarName(column,varName) + &
           gsv_getLevFromVarLev(statevector,varLevIndex)
    end do

    if (inputStateVectorType == 'nl' .and. rejectObsOutsideGlobalGrid) then
      rejectOutsideObs = .true.
    end if

    numStep = stateVector%numStep
    numHeader = obs_numHeader(obsSpaceData)

    call oti_setup(intInfo%oti, obsSpaceData, numStep,  &
                   1, numHeader, &
                   interpType_opt=timeInterpType, flagObsOutside_opt=.true.)

    ! copy the horizontal grid object
    intInfo%hco => stateVector%hco

    ! setup the information for wind rotation
    if (gsv_varExist(varName='UU') .and. &
        gsv_varExist(varName='VV') .and.  &
        stateVector%hco%rotated) then
      intInfo%rotatedWinds = .true.
      call uvr_Setup(intInfo%uvr, & ! INOUT
                     stateVector%hco)  ! IN
    end if

    if ( stateVector%hco%grtyp == 'G' .or. &
        (stateVector%hco%grtyp == 'Z' .and. stateVector%hco%global) ) then
      intInfo%periodic = .true.
    else
      intInfo%periodic = .false.
    end if

    if (.not. beSilent) then
      write(*,*) 'setupInterpInfoTiles (s2c): Set the varLevIndex for lat-lon of GZsfc'
      call msg_memUsage('setupInterpInfoTiles')
    end if
    call setSfcVarLevIndex(intInfo,stateVector)

    if (.not. beSilent) then
      write(*,*) 'setupInterpInfoTiles (s2c): Obtain profile of lat-lon for each obs'
      call msg_memUsage('setupInterpInfoTiles')
    end if
    call getObsLatLon(intInfo, obsSpaceData, stateVector)

    if (.not. beSilent) then
      write(*,*) 'setupInterpInfoTiles (s2c): Communicate all obs tile ID and locations globally'
      call msg_memUsage('setupInterpInfoTiles')
    end if
    call getObsTileMpiId(intInfo, obsSpaceData, stateVector)

    if (.not. beSilent) then
      write(*,*) 'setupInterpInfoTiles (s2c): Collect information on obs in my tile'
      call msg_memUsage('setupInterpInfoTiles')
    end if
    call getMyInterpObsLatLon(intInfo, stateVector)
    call getMyInterpObsXYposSubGridIndex(intInfo, obsSpaceData, stateVector)
    if (intInfo%rotatedWinds) then
      call getMyInterpObsRotLatLon(intInfo, stateVector)
    end if

    if (.not. beSilent) then
      write(*,*) 'setupInterpInfoTiles (s2c): Compute interpolation weights and indices for my obs'
      call msg_memUsage('setupInterpInfoTiles')
    end if
    call utl_tmg_start(33,'------s2c_SetupWeights')
    call getMyInterpWeights(intInfo, stateVector)
    call utl_tmg_stop(33)

    if (.not. beSilent) then
      write(*,*) 'setupInterpInfoTiles (s2c): Figure out which extra grid points I need (i.e. halo)'
    end if
    call getMyHalo(intInfo, stateVector)

    if (.not. beSilent) then
      write(*,*) 'setupInterpInfoTiles (s2c): Adjust lat/lonIndexDepot for halo locations'
      call msg_memUsage('setupInterpInfoTiles')
    end if
    call utl_tmg_start(33,'------s2c_SetupWeights')
    call getMyLatLonIndexHalo(intInfo, stateVector)
    call utl_tmg_stop(33)

    ! reject obs in obsSpaceData if any processor has zero weight
    ! called when a mask exists to catch land contaminated ocean obs
    if ( stateVector%oceanMask%maskPresent ) then
      call rejectZeroWeightObs(intInfo,obsSpaceData)
    end if

    if (.not. beSilent) call msg_memUsage('setupInterpInfoTiles')
    if (.not. beSilent) write(*,*) 'setupInterpInfoTiles (s2c): FINISHED'

  end subroutine setupInterpInfoTiles

  !--------------------------------------------------------------------------
  ! rejectZeroWeightObs (called by setupInterpInfoTiles)
  !--------------------------------------------------------------------------
  subroutine rejectZeroWeightObs(intInfo, obsSpaceData)
    !
    !:Purpose: To flag an observation in obsSpaceData as being rejected if
    !          it has zero interpolation weight (usually because an ocean
    !          obs is touching land) on any mpi task.
    !
    implicit none

    ! Arguments:
    type(struct_interpInfoTiles), intent(inout) :: intInfo
    type(struct_obs)       , intent(inout) :: obsSpaceData

    ! Locals:
    integer :: myHeaderIndex, headerIndex, varLevIndex, procIndex, depotIndex
    integer :: numVarLev, myNumHeader, myNumHeaderMax, bodyIndexBeg, bodyIndexEnd, bodyIndex
    integer :: headerCount(mmpi_nprocs), myNumHeaders(mmpi_nprocs)
    integer, save :: numWrites = 0
    logical, allocatable :: rejectObsSend(:,:), rejectObsRecv(:,:)
    integer, allocatable :: headIndexSend(:,:), headIndexRecv(:,:)

    write(*,*) 'rejectZeroWeightObs (s2c): Starting'

    myNumHeader = intInfo%myInterpNumHeader
    numVarLev = intInfo%numVarLevState

    ! Determine max number of headers on this mpi task from each mpi task
    myNumHeaders(:) = 0
    do myHeaderIndex = 1, myNumHeader
      procIndex = intInfo%myInterpObsMpiIdSrc(myHeaderIndex) + 1
      myNumHeaders(procIndex) = myNumHeaders(procIndex) + 1
    end do
    call mmpi_allReduce(maxval(myNumHeaders(:)), myNumHeaderMax, mmpi_max)

    allocate(rejectObsSend(myNumHeaderMax,mmpi_nprocs))
    allocate(rejectObsRecv(myNumHeaderMax,mmpi_nprocs))
    rejectObsSend(:,:) = .true.
    rejectObsRecv(:,:) = .true.
    allocate(headIndexSend(myNumHeaderMax,mmpi_nprocs))
    allocate(headIndexRecv(myNumHeaderMax,mmpi_nprocs))
    headIndexSend(:,:) = 0
    headIndexRecv(:,:) = 0

    headerCount(:) = 0
    do myHeaderIndex = 1, myNumHeader
      procIndex = intInfo%myInterpObsMpiIdSrc(myHeaderIndex) + 1
      headerCount(procIndex) = headerCount(procIndex) + 1

      headIndexSend(headerCount(procIndex),procIndex) = intInfo%myInterpObsHeaderIndex(myHeaderIndex)

      do varLevIndex = 1, numVarLev
        do depotIndex = intInfo%depotIndexBeg(varLevIndex,myHeaderIndex),  &
                        intInfo%depotIndexEnd(varLevIndex,myHeaderIndex)

          if (intInfo%interpWeightDepot(depotIndex) > 0.0d0) then
            rejectObsSend(headerCount(procIndex),procIndex) = .false.
          end if

        end do ! depotIndex
      end do ! varLevIndex

    end do ! myHeaderIndex

    ! do global communication of reject flags and original headerIndex
    call mmpi_alltoall(rejectObsSend, rejectObsRecv)
    call mmpi_alltoall(headIndexSend, headIndexRecv)

    ! modify obsSpaceData based on reject flags
    do procIndex = 1, mmpi_nprocs
      do myHeaderIndex = 1, myNumHeaderMax
        ! Get the original headerIndex for obsSpaceData
        headerIndex = headIndexRecv(myHeaderIndex,procIndex)

        ! Check if we already reached the last header from the mpi task procIndex
        if (headerIndex == 0) cycle

        if (rejectObsRecv(myHeaderIndex,procIndex)) then

          numWrites = numWrites + 1
          if (numWrites < maxNumWrites) then
            write(*,*) 'rejectZeroWeightObs (s2c): Rejecting OBS with zero weight, index ', headerIndex
          else if (numWrites == maxNumWrites) then
            write(*,*) 'rejectZeroWeightObs (s2c): More rejects, but reached maximum number of writes to the listing.'
          end if

          bodyIndexBeg = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
          bodyIndexEnd = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) + bodyIndexBeg -1
          do bodyIndex = bodyIndexBeg, bodyIndexEnd
            call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
          end do
          call obs_headSet_i(obsSpaceData, OBS_ST1, headerIndex,  &
               ibset( obs_headElem_i(obsSpaceData, OBS_ST1, headerIndex), 05))
        end if

      end do
    end do

    deallocate(rejectObsSend)
    deallocate(rejectObsRecv)
    deallocate(headIndexSend)
    deallocate(headIndexRecv)

    write(*,*) 'rejectZeroWeightObs (s2c): Finished'

  end subroutine rejectZeroWeightObs

  !---------------------------------------------------------
  ! getMyLatLonIndexHalo (called by setupInterpInfoTiles)
  !---------------------------------------------------------
  subroutine getMyLatLonIndexHalo(intInfo, stateVector)
    !
    ! :Purpose: Adjust the lat/lonIndexDepot for locations in the halo so
    !           that they point to the haloIndex instead of lat/lonIndex.
    !
    implicit none

    ! Arguments:
    type(struct_interpInfoTiles), intent(inout) :: intInfo     ! Interpolation info structure
    type(struct_gsv), target    , intent(in)    :: stateVector ! stateVector object

    ! Locals:
    integer :: myHeaderIndex, depotIndex, haloIndex
    integer :: lonIndex, latIndex
    integer :: varLevIndex, numVarLev

    numVarLev = intInfo%numVarLevState

    ! Go through all grid points in the depot
    do myHeaderIndex = 1, intInfo%myInterpNumHeader
      do varLevIndex = 1, numVarLev
        depotLoop: do depotIndex = intInfo%depotIndexBeg(varLevIndex,myHeaderIndex),  &
                                   intInfo%depotIndexEnd(varLevIndex,myHeaderIndex)

          ! Identify if this grid point in the depot is outside the tile interior
          latIndex = intInfo%latIndexDepot(depotIndex)
          lonIndex = intInfo%lonIndexDepot(depotIndex)
          if (latIndex < stateVector%myLatBeg .or. latIndex > stateVector%myLatEnd .or. &
              lonIndex < stateVector%myLonBeg .or. lonIndex > stateVector%myLonEnd) then

            ! Grid is outside, now find haloIndex value for this lat/lonIndex
            haloLoop: do haloIndex = 1, intInfo%myHaloSize
              if (latIndex == intInfo%myHaloLatIndex(haloIndex) .and. &
                  lonIndex == intInfo%myHaloLonIndex(haloIndex)) then

                ! Modify lat/lonIndexDepot values to record haloIndex value and exit loop
                intInfo%latIndexDepot(depotIndex) = valueIsInHalo
                intInfo%lonIndexDepot(depotIndex) = haloIndex
                exit haloLoop

              end if
            end do haloLoop

          end if

        end do depotLoop
      end do
    end do

  end subroutine getMyLatLonIndexHalo

  !---------------------------------------------------------
  ! getMyHalo (called by setupInterpInfoTiles)
  !---------------------------------------------------------
  subroutine getMyHalo(intInfo, stateVector)
    !
    ! :Purpose: Determine what grid points need to be communicated between
    !           lat-lon tiles, i.e. the halo.
    !
    implicit none

    ! Arguments:
    type(struct_interpInfoTiles), intent(inout) :: intInfo     ! Interpolation info structure
    type(struct_gsv), target    , intent(in)    :: stateVector ! stateVector object

    ! Locals:
    integer, parameter :: maxTempHaloSize = 10000000
    integer :: myHeaderIndex, depotIndex, haloIndex, haloIndex2, tileMpiIdX, tileMpiIdY
    integer :: procIndex, msgIndex, mpiMsgSize, lonIndex, latIndex
    integer :: varLevIndex, numVarLev
    integer, allocatable :: latVecSend(:,:), lonVecSend(:,:), latVecRecv(:,:), lonVecRecv(:,:)
    integer :: myHaloSizeMax, myHaloSizeMpi(mmpi_nprocs), yourHaloSizeMpi(mmpi_nprocs)

    numVarLev = intInfo%numVarLevState

    ! Temporarily allocate lat-lon arrays just for counting halo size
    allocate(intInfo%myHaloLatIndex(maxTempHaloSize))
    allocate(intInfo%myHaloLonIndex(maxTempHaloSize))

    ! Count number of grid points that are outside the lat-lon tile
    haloIndex = 0
    do myHeaderIndex = 1, intInfo%myInterpNumHeader
      do varLevIndex = 1, numVarLev
        depotLoop0: do depotIndex = intInfo%depotIndexBeg(varLevIndex,myHeaderIndex),  &
                                    intInfo%depotIndexEnd(varLevIndex,myHeaderIndex)

          ! Identify if this grid point is outside the tile interior
          latIndex = intInfo%latIndexDepot(depotIndex)
          lonIndex = intInfo%lonIndexDepot(depotIndex)
          if (latIndex < stateVector%myLatBeg .or. latIndex > stateVector%myLatEnd .or. &
              lonIndex < stateVector%myLonBeg .or. lonIndex > stateVector%myLonEnd) then

            ! Check if this location is already recorded as part of the halo
            if (haloIndex > 0) then
              do haloIndex2 = 1, haloIndex
                if (latIndex == intInfo%myHaloLatIndex(haloIndex2) .and. &
                    lonIndex == intInfo%myHaloLonIndex(haloIndex2)) then
                  cycle depotLoop0
                end if
              end do
            end if

            haloIndex = haloIndex + 1

            ! Lat-lon of halo grid point
            if (haloIndex > size(intInfo%myHaloLatIndex)) then
              call utl_abort('getMyHalo (s2c): temporary allocation size for halo not big enough')
            end if
            intInfo%myHaloLatIndex(haloIndex) = intInfo%latIndexDepot(depotIndex)
            intInfo%myHaloLonIndex(haloIndex) = intInfo%lonIndexDepot(depotIndex)

          end if
        end do depotLoop0
      end do
    end do
    intInfo%myHaloSize = haloIndex
    write(*,*) 'getMyHalo (s2c): myHaloSize    = ', intInfo%myHaloSize

    ! Deallocate temporarily allocated arrays
    deallocate(intInfo%myHaloLatIndex)
    deallocate(intInfo%myHaloLonIndex)

    ! Allocate arrays to store the halo information
    allocate(intInfo%myHaloLatIndex(intInfo%myHaloSize))
    allocate(intInfo%myHaloLonIndex(intInfo%myHaloSize))
    allocate(intInfo%myHaloMpiIdSrc(intInfo%myHaloSize))
    allocate(intInfo%myHaloMpiTag(intInfo%myHaloSize))

    ! Store my halo information
    haloIndex = 0
    do myHeaderIndex = 1, intInfo%myInterpNumHeader
      do varLevIndex = 1, numVarLev
        depotLoop1: do depotIndex = intInfo%depotIndexBeg(varLevIndex,myHeaderIndex),  &
                                    intInfo%depotIndexEnd(varLevIndex,myHeaderIndex)

          ! Identify if this grid point is outside the tile interior
          latIndex = intInfo%latIndexDepot(depotIndex)
          lonIndex = intInfo%lonIndexDepot(depotIndex)
          if (latIndex < stateVector%myLatBeg .or. latIndex > stateVector%myLatEnd .or. &
              lonIndex < stateVector%myLonBeg .or. lonIndex > stateVector%myLonEnd) then

            ! Check if this location is already recorded as part of the halo
            if (haloIndex > 0) then
              do haloIndex2 = 1, haloIndex
                if (latIndex == intInfo%myHaloLatIndex(haloIndex2) .and. &
                    lonIndex == intInfo%myHaloLonIndex(haloIndex2)) then
                  cycle depotLoop1
                end if
              end do
            end if

            ! Add this horizontal location to the halo
            haloIndex = haloIndex + 1
            if (haloIndex > intInfo%myHaloSize) then
              call utl_abort('getMyHalo (s2c): haloIndex is too big')
            end if

            ! Lat-lon of halo grid point
            intInfo%myHaloLatIndex(haloIndex) = intInfo%latIndexDepot(depotIndex)
            intInfo%myHaloLonIndex(haloIndex) = intInfo%lonIndexDepot(depotIndex)

            ! Adjust lonIndex used to determine source tile in special cases with periodic domain
            if (intInfo%periodic .and. intInfo%myHaloLonIndex(haloIndex) > stateVector%ni) then
              lonIndex = 1
            else if (intInfo%periodic .and. intInfo%myHaloLonIndex(haloIndex) < 1) then
              lonIndex = stateVector%ni
            else
              lonIndex = intInfo%myHaloLonIndex(haloIndex)
            end if
            latIndex = intInfo%myHaloLatIndex(haloIndex)

            ! MPI tile id in X, Y directions and global tile id
            tileMpiIdX = utl_findloc( &
                 lonIndex >= stateVector%allLonBeg(:) .and. &
                 lonIndex <= stateVector%allLonEnd(:)) - 1
            tileMpiIdY = utl_findloc( &
                 latIndex >= stateVector%allLatBeg(:) .and. &
                 latIndex <= stateVector%allLatEnd(:)) - 1
            intInfo%myHaloMpiIdSrc(haloIndex) = tileMpiIdX + tileMpiIdY*mmpi_npex

            ! Set mpi tag to distinguish between multiple messages received from the same mpi task
            intInfo%myHaloMpiTag(haloIndex) = count(intInfo%myHaloMpiIdSrc(1:haloIndex) == &
                                                    intInfo%myHaloMpiIdSrc(haloIndex))

          end if

        end do depotLoop1

      end do ! varLevIndex
    end do ! myHeaderIndex

    ! Send halo information to source MPI tiles so they will know what to send
    call mmpi_allReduce(intInfo%myHaloSize, myHaloSizeMax, mmpi_max)

    ! Number of grid points needed from each mpi task
    myHaloSizeMpi(:) = 0
    do procIndex = 1, mmpi_nprocs
      myHaloSizeMpi(procIndex) = count(intInfo%myHaloMpiIdSrc(:) == procIndex-1)
    end do
    call mmpi_alltoall(myHaloSizeMpi, yourHaloSizeMpi)

    ! Maximum size of MPI message for halo
    call mmpi_allReduce(maxval(myHaloSizeMpi), mpiMsgSize, mmpi_max)
    write(*,*) 'getMyHalo (s2c): mpiMsgSize      = ', mpiMsgSize

    allocate(latVecSend(mpiMsgSize,mmpi_nprocs))
    allocate(lonVecSend(mpiMsgSize,mmpi_nprocs))
    allocate(latVecRecv(mpiMsgSize,mmpi_nprocs))
    allocate(lonVecRecv(mpiMsgSize,mmpi_nprocs))
    latVecSend(:,:) = -999
    lonVecSend(:,:) = -999

    do haloIndex = 1, intInfo%myHaloSize
      procIndex = intInfo%myHaloMpiIdSrc(haloIndex) + 1
      msgIndex = count(latVecSend(:,procIndex) >= 0) + 1
      latVecSend(msgIndex,procIndex) = intInfo%myHaloLatIndex(haloIndex)
      lonVecSend(msgIndex,procIndex) = intInfo%myHaloLonIndex(haloIndex)
    end do

    call mmpi_alltoall(latVecSend, latVecRecv)
    call mmpi_alltoall(lonVecSend, lonVecRecv)

    ! Allocate arrays to store the halo information
    intInfo%yourHaloSize = count(latVecRecv(:,:) > 0)
    allocate(intInfo%yourHaloLatIndex(intInfo%yourHaloSize))
    allocate(intInfo%yourHaloLonIndex(intInfo%yourHaloSize))
    allocate(intInfo%yourHaloMpiIdDst(intInfo%yourHaloSize))
    allocate(intInfo%yourHaloMpiTag(intInfo%yourHaloSize))

    haloIndex = 0
    do procIndex = 1, mmpi_nprocs
      do msgIndex = 1, yourHaloSizeMpi(procIndex)
        haloIndex = haloIndex + 1

        ! Set the mpi id
        intInfo%yourHaloMpiIdDst(haloIndex) = procIndex - 1

        ! Lat and lon received from allToAll communication
        intInfo%yourHaloLatIndex(haloIndex) = latVecRecv(msgIndex,procIndex)
        intInfo%yourHaloLonIndex(haloIndex) = lonVecRecv(msgIndex,procIndex)

        ! Set mpi tag to distinguish between multiple messages received from the same mpi task
        intInfo%yourHaloMpiTag(haloIndex) = count(intInfo%yourHaloMpiIdDst(1:haloIndex) == &
                                                  intInfo%yourHaloMpiIdDst(haloIndex))

      end do
    end do

    deallocate(latVecSend)
    deallocate(lonVecSend)
    deallocate(latVecRecv)
    deallocate(lonVecRecv)

  end subroutine getMyHalo

  !---------------------------------------------------------
  ! getMyInterpWeights (called by setupInterpInfoTiles)
  !---------------------------------------------------------
  subroutine getMyInterpWeights(intInfo, stateVector)
    !
    ! :Purpose: Compute the interpolation weights and lat-lon indexes
    !           for each grid point involved in the interpolation to each
    !           observation location within my tile.
    !
    implicit none

    ! Arguments:
    type(struct_interpInfoTiles), intent(inout) :: intInfo     ! Interpolation info structure
    type(struct_gsv), target    , intent(in)    :: stateVector ! stateVector object

    ! Locals:
    integer :: myHeaderIndex, varLevIndex, depotIndex, depotSize
    integer :: gridIndex, latIndex, lonIndex, numGridpt, numVarLev
    real(8) :: lat, lon
    real(4) :: footprintRadius_r4
    real(kdkind), allocatable :: positionArray(:,:)
    type(kdtree2), pointer  :: tree

    numVarLev = intInfo%numVarLevState

    ! create kdtree to use in footprint operator, if any footprint radius > 0.
    if ( any(intInfo%myInterpObsFootprint_r4(:,:) > 0.0) ) then
      write(*,*) 'getMyInterpWeights (s2c): footPrint operator is used for inputStateVectorType=', &
                 intInfo%inputStateVectorType

      if ( (intInfo%inputStateVectorType == 'nl' .and. .not. associated(tree_nl))   .or. &
           (intInfo%inputStateVectorType == 'tl' .and. .not. associated(tree_tlad)) .or. &
           (intInfo%inputStateVectorType == 'ad' .and. .not. associated(tree_tlad)) ) then

        write(*,*) 'getMyInterpWeights (s2c): start creating kdtree for inputStateVectorType=', &
                   intInfo%inputStateVectorType
        call msg_memUsage('getMyInterpWeights')

        allocate(positionArray(3,statevector%hco%ni*statevector%hco%nj))

        gridIndex = 0
        do latIndex = 1, statevector%hco%nj
          do lonIndex = 1, statevector%hco%ni
            gridIndex = gridIndex + 1
            lat = real(stateVector % hco % lat2d_4(lonIndex,latIndex), 8)
            lon = real(stateVector % hco % lon2d_4(lonIndex,latIndex), 8)

            positionArray(:,gridIndex) = kdtree2_3dPosition(lon, lat)

          end do
        end do

        nullify(tree)
        tree => kdtree2_create(positionArray, sort=.false., rearrange=.true.)

        if ( intInfo%inputStateVectorType == 'nl' ) then
          tree_nl => tree
        else
          tree_tlad => tree
        end if

        deallocate(positionArray)

        write(*,*) 'getMyInterpWeights (s2c): done creating kdtree for inputStateVectorType = ', &
                   intInfo%inputStateVectorType
        call msg_memUsage('getMyInterpWeights')

      end if
    end if

    ! Allocate index pointer arrays into depot
    allocate(intInfo%depotIndexBeg(numVarLev,intInfo%myInterpNumHeader))
    allocate(intInfo%depotIndexEnd(numVarLev,intInfo%myInterpNumHeader))

    ! First determine the size of the depot and set the index pointer arrays
    depotSize = 0
    do myHeaderIndex = 1, intInfo%myInterpNumHeader

      do varLevIndex = 1, numVarLev
        ! Set the starting index for this myHeaderIndex
        intInfo%depotIndexBeg(varLevIndex,myHeaderIndex) = depotSize + 1

        footprintRadius_r4 = intInfo%myInterpObsFootprint_r4(varLevIndex,myHeaderIndex)
        if ( footprintRadius_r4 > 0.0 ) then

          call getMyInterpWeightsFootprint(footPrintRadius_r4, myHeaderIndex, &
                                           varLevIndex, numGridpt)

        else if ( utl_isEqual(footprintRadius_r4, bilinearFootprint) ) then

          call getMyInterpWeightsBilinear(myHeaderIndex, varLevIndex, numGridpt)

        else if ( utl_isEqual(footprintRadius_r4, lakeFootprint) ) then

          call getMyInterpWeightsLake(myHeaderIndex, varLevIndex, numGridpt)

        else if ( utl_isEqual(footprintRadius_r4, nearestNeighbourFootprint) ) then

          call getMyInterpWeightsNearestNeighbor(myHeaderIndex, varLevIndex, numGridpt)

        else

          write(*,*) 'footPrintRadius_r4 = ', footPrintRadius_r4
          call utl_abort('getMyInterpWeights (s2c): this type of interpolation not implemented')

        end if
        depotSize = depotSize + numGridpt

        ! Set the ending index for this myHeaderIndex
        intInfo%depotIndexEnd(varLevIndex,myHeaderIndex) = depotSize
      end do

    end do
    write(*,*) 'getMyInterpWeights (s2c): myInterpNumHeader, depotSize = ',  &
                                  intInfo%myInterpNumHeader, depotSize

    ! Allocate the depot variables
    allocate(intInfo%latIndexDepot(depotSize))
    allocate(intInfo%lonIndexDepot(depotSize))
    allocate(intInfo%interpWeightDepot(depotSize))

    ! Now assign the values in the depot
    do myHeaderIndex = 1, intInfo%myInterpNumHeader

      do varLevIndex = 1, numVarLev

        footprintRadius_r4 = intInfo%myInterpObsFootprint_r4(varLevIndex,myHeaderIndex)
        if ( footprintRadius_r4 > 0.0 ) then

          call getMyInterpWeightsFootprint(footPrintRadius_r4, myHeaderIndex, &
                                           varLevIndex, numGridpt)

        else if ( utl_isEqual(footprintRadius_r4, bilinearFootprint) ) then

          call getMyInterpWeightsBilinear(myHeaderIndex, varLevIndex, numGridpt)

        else if ( utl_isEqual(footprintRadius_r4, lakeFootprint) ) then

          call getMyInterpWeightsLake(myHeaderIndex, varLevIndex, numGridpt)

        else if ( utl_isEqual(footprintRadius_r4, nearestNeighbourFootprint) ) then

          call getMyInterpWeightsNearestNeighbor(myHeaderIndex, varLevIndex, numGridpt)

        else

          write(*,*) 'footPrintRadius_r4 = ', footPrintRadius_r4
          call utl_abort('getMyInterpWeights (s2c): this type of interpolation not implemented')

        end if

      end do ! varLevIndex
    end do ! myHeaderIndex

  contains

    !--------------------------------------------------------------------------
    ! getMyInterpWeightsBilinear (contained in getMyInterpWeights)
    !--------------------------------------------------------------------------
    subroutine getMyInterpWeightsBilinear(myHeaderIndex, varLevIndex, numGridpt)
      !
      ! :Purpose: Either just count or also assign interpolation weights and indexes.
      !
      implicit none

      ! Arguments:
      integer, intent(in)  :: myHeaderIndex ! headerIndex to be treated
      integer, intent(in)  :: varLevIndex   ! varLevIndex to be treated
      integer, intent(out) :: numGridpt     ! return total number of grid points

      ! Locals:
      integer, parameter :: leftIndex = 1, rightIndex = 2, bottomIndex = 1, topIndex = 2
      integer :: niP1, lonIndex, latIndex, lonIndexP1, lonIndexP1ForMask, iPoint
      integer :: latIndexVec(4), lonIndexVec(4)
      real(4) :: xpos_r4, ypos_r4
      real(8) :: WeightVec(4), dldx, dldy, weightsSum
      logical :: mask(2,2)

      xpos_r4 = intInfo%myInterpObsXpos_r4(varLevIndex,myHeaderIndex)
      ypos_r4 = intInfo%myInterpObsYpos_r4(varLevIndex,myHeaderIndex)

      ! Allow for periodicity in Longitude for global Gaussian grid
      if (intInfo%periodic) then
        niP1 = statevector%ni + 1
      else
        niP1 = statevector%ni
      end if

      ! Find the lower-left grid point next to the observation
      if ( xpos_r4 >= real(niP1) ) then
        xpos_r4 = real(niP1)
        lonIndex = niP1 - 1
      else if ( xpos_r4 < 1.0 ) then
        xpos_r4 = 1.0
        lonIndex = 1
      else
        lonIndex = floor(xpos_r4)
      end if

      if ( ypos_r4 >= real(statevector%nj) ) then
        ypos_r4 = real(statevector%nj)
        latIndex = statevector%nj - 1
      else if ( ypos_r4 < 1.0 ) then
        ypos_r4 = 1.0
        latIndex = 1
      else
        latIndex = floor(ypos_r4)
      end if

      if ( stateVector%hco%grtyp == 'U' ) then
        if ( utl_isEqual(ypos_r4, real(stateVector%nj/2)) ) then
          latIndex = floor(ypos_r4) - 1
        end if
      end if

      lonIndexP1 = lonIndex + 1

      ! Check if location is in between Yin and Yang (should not happen)
      if ( stateVector%hco%grtyp == 'U' ) then
        if ( ypos_r4 > real(stateVector%nj/2) .and.  &
             ypos_r4 < real((stateVector%nj/2)+1) ) then
          write(*,*) 'getMyInterpWeightsBilinear (s2c): WARNING, obs position between Yin and Yang!'
          write(*,*) '   xpos, ypos = ', xpos_r4, ypos_r4
        end if
      end if

      if ( stateVector%oceanMask%maskPresent ) then
        ! abort if 3D mask is present, since we may not handle this situation correctly
        if ( stateVector%oceanMask%nLev > 1 ) then
          call utl_abort('getMyInterpWeightsBilinear (s2c): 3D mask present - this case not properly handled')
        end if
        ! Handle periodicity in longitude for ocean mask value
        if (lonIndexP1 == stateVector%ni + 1) then
          lonIndexP1ForMask = 1
        else
          lonIndexP1ForMask = lonIndexP1
        end if
        mask(leftIndex ,bottomIndex) = stateVector%oceanMask%mask(lonIndex         ,latIndex    ,1)
        mask(rightIndex,bottomIndex) = stateVector%oceanMask%mask(lonIndexP1ForMask,latIndex    ,1)
        mask(leftIndex ,topIndex   ) = stateVector%oceanMask%mask(lonIndex         ,latIndex + 1,1)
        mask(rightIndex,topIndex   ) = stateVector%oceanMask%mask(lonIndexP1ForMask,latIndex + 1,1)
      else
        mask(:,:) = .true.
      end if

      WeightVec(:) = 0
      numGridpt = 0

      ! Compute the 4 weights of the bilinear interpolation
      dldx = real(xpos_r4,8) - real(lonIndex,8)
      dldy = real(ypos_r4,8) - real(latIndex,8)
      if (NNInterpForCloudVars) then
        call utl_abort('getMyInterpWeightsBilinear (s2c): NNInterpForCloudVars true is not supported')
      end if

      if ( mask(leftIndex ,bottomIndex) ) then
        numGridpt = numGridpt + 1
        latIndexVec(numGridpt) = latIndex
        lonIndexVec(numGridpt) = lonIndex
        WeightVec(numGridpt) = (1.d0-dldx) * (1.d0-dldy)
      end if

      if ( mask(rightIndex,bottomIndex) ) then
        numGridpt = numGridpt + 1
        latIndexVec(numGridpt) = latIndex
        lonIndexVec(numGridpt) = lonIndexP1
        WeightVec(numGridpt) =       dldx  * (1.d0-dldy)
      end if

      if ( mask(leftIndex ,topIndex   ) ) then
        numGridpt = numGridpt + 1
        latIndexVec(numGridpt) = latIndex + 1
        lonIndexVec(numGridpt) = lonIndex
        WeightVec(numGridpt) = (1.d0-dldx) *       dldy
      end if

      if ( mask(rightIndex,topIndex   ) ) then
        numGridpt = numGridpt + 1
        latIndexVec(numGridpt) = latIndex + 1
        lonIndexVec(numGridpt) = lonIndexP1
        WeightVec(numGridpt) =       dldx  *       dldy
      end if

      weightsSum = sum(WeightVec(1:numGridpt))
      if ( weightsSum > 0.d0 ) then
        WeightVec(1:numGridpt) = WeightVec(1:numGridpt) / weightsSum
      end if

      ! If the depot is allocated, then we fill in the values
      if ( allocated(intInfo%interpWeightDepot) ) then

        depotIndex = intInfo%depotIndexBeg(varLevIndex,myHeaderIndex)

        do ipoint = 1, numGridpt

          intInfo%interpWeightDepot(depotIndex) = WeightVec(ipoint)
          intInfo%latIndexDepot(depotIndex)     = latIndexVec(ipoint)
          intInfo%lonIndexDepot(depotIndex)     = lonIndexVec(ipoint)

          depotIndex = depotIndex + 1

        end do

      end if

    end subroutine getMyInterpWeightsBilinear

    !--------------------------------------------------------------------------
    ! getMyInterpWeightsFootprint (contained in getMyInterpWeights)
    !--------------------------------------------------------------------------
    subroutine getMyInterpWeightsFootprint(fpr, myHeaderIndex, varLevIndex, numGridpt)
      !
      !:Purpose: To determine the grid points and their associated weights
      !          for the footprint horizontal interpolation.
      !
      implicit none

      ! Arguments:
      real(4)                , intent(in)    :: fpr            ! footprint radius (metres)
      integer                , intent(in)    :: myHeaderIndex  ! headerIndex to be treated
      integer                , intent(in)    :: varLevIndex    ! varLevIndex to be treated
      integer                , intent(out)   :: numGridpt      ! return total number of grid points

      ! Locals:
      integer :: depotIndex, latIndexCentre, lonIndexCentre
      integer :: subGridIndex, numLocalGridptsFoundSearch
      integer :: ipoint, gridptCount
      integer :: lonIndex, latIndex, resultsIndex, gridIndex
      integer :: lonIndexVec(maxNumLocalGridptsSearch), latIndexVec(maxNumLocalGridptsSearch)
      real(8) :: lonObs, latObs
      real(4) :: xpos_r4, ypos_r4, lonObs_deg_r4, latObs_deg_r4
      type(kdtree2_result)   :: searchResults(maxNumLocalGridptsSearch)
      real(kdkind)           :: refPosition(3), maxRadiusSquared
      type(kdtree2), pointer :: tree

      numGridpt = 0

      ! Determine the grid point nearest the observation.

      latObs       = intInfo%myInterpObsLat(varLevIndex,myHeaderIndex)
      lonObs       = intInfo%myInterpObsLon(varLevIndex,myHeaderIndex)
      xpos_r4      = intInfo%myInterpObsXpos_r4(varLevIndex,myHeaderIndex)
      ypos_r4      = intInfo%myInterpObsYpos_r4(varLevIndex,myHeaderIndex)
      subGridIndex = intInfo%myInterpObsSubGridIndex(varLevIndex,myHeaderIndex)

      lonIndexCentre = nint(xpos_r4)
      latIndexCentre = nint(ypos_r4)

      if ( subGridIndex == 3 ) then
        call utl_abort('getMyInterpWeightsFootprint: two subGrids involved is not supported')
      end if

      ! Return if observation is not on the grid, or masked.
      if ( lonIndexCentre < 1 .or. lonIndexCentre > statevector%hco%ni .or.  &
           latIndexCentre < 1 .or. latIndexCentre > statevector%hco%nj ) return

      if ( stateVector%oceanMask%maskPresent ) then
        ! abort if 3D mask is present, since we may not handle this situation correctly
        if ( stateVector%oceanMask%nLev > 1 ) then
          call utl_abort('getMyInterpWeightsFootprint: 3D mask present - this case not properly handled')
        end if

        if ( .not. stateVector%oceanMask%mask(lonIndexCentre,latIndexCentre,1) ) return
      end if

      ! do the search
      maxRadiusSquared = real(fpr,8) ** 2
      refPosition(:) = kdtree2_3dPosition(lonObs, latObs)
      nullify(tree)
      if ( intInfo%inputStateVectorType == 'nl' ) then
        if ( associated(tree_nl) ) then
          tree => tree_nl
        else
          call utl_abort('getMyInterpWeightsFootprint: tree_nl is not allocated!')
        end if
      else if ( intInfo%inputStateVectorType == 'tl' .or. &
           intInfo%inputStateVectorType == 'ad' ) then
        if ( associated(tree_tlad) ) then
          tree => tree_tlad
        else
          call utl_abort('getMyInterpWeightsFootprint: tree_tlad is not allocated!')
        end if
      end if
      call kdtree2_r_nearest(tp=tree, qv=refPosition, r2=maxRadiusSquared, &
                             nfound=numLocalGridptsFoundSearch, &
                             nalloc=maxNumLocalGridptsSearch, &
                             results=searchResults)

      if (numLocalGridptsFoundSearch > maxNumLocalGridptsSearch ) then
        call utl_abort('getMyInterpWeightsFootprint: the parameter maxNumLocalGridptsSearch must be increased')
      else if ( numLocalGridptsFoundSearch < minNumLocalGridptsSearch .and. useFootprintForTovs ) then
        write(*,*) 'getMyInterpWeightsFootprint: Warning! For TOVS headerIndex=', myHeaderIndex, &
                   ' number of grid points found within footprint radius=', fpr, ' is less than ', &
                   minNumLocalGridptsSearch
      end if

      ! ensure at least the nearest neighbor is included in lonIndexVec/latIndexVec
      ! if footprint size is smaller than the grid spacing.
      gridptCount = 1
      lonIndexVec(gridptCount) = lonIndexCentre
      latIndexVec(gridptCount) = latIndexCentre

      ! fill the rest of lonIndexVec/latIndexVec
      gridLoop1: do resultsIndex = 1, numLocalGridptsFoundSearch

        gridIndex = searchResults(resultsIndex)%idx
        if ( gridIndex < 1 .or. gridIndex > statevector%hco%ni * statevector%hco%nj ) then
          write(*,*) 'getMyInterpWeightsFootprint: gridIndex=', gridIndex
          call utl_abort('getMyInterpWeightsFootprint: gridIndex out of bound.')
        end if

        latIndex = (gridIndex - 1) / statevector%hco%ni + 1
        lonIndex = gridIndex - (latIndex - 1) * statevector%hco%ni
        if ( lonIndex < 1 .or. lonIndex > statevector%hco%ni .or. &
             latIndex < 1 .or. latIndex > statevector%hco%nj ) then
          write(*,*) 'getMyInterpWeightsFootprint: lonIndex=', lonIndex, ',latIndex=', latIndex
          call utl_abort('getMyInterpWeightsFootprint: lonIndex/latIndex out of bound.')
        end if

        if ( stateVector%oceanMask%maskPresent ) then
          if ( .not. stateVector%oceanMask%mask(lonIndex,latIndex,1) ) cycle gridLoop1
        end if

        if ( lonIndex == lonIndexCentre .and. latIndex == latIndexCentre ) cycle gridLoop1

        gridptCount = gridptCount + 1
        lonIndexVec(gridptCount) = lonIndex
        latIndexVec(gridptCount) = latIndex

      end do gridLoop1

      ! If the depot is allocated, then we fill in the values
      if ( allocated(intInfo%interpWeightDepot) ) then

        depotIndex = intInfo%depotIndexBeg(varLevIndex,myHeaderIndex)

        do ipoint = 1, gridptCount

          intInfo%interpWeightDepot(depotIndex) = 1.0d0 / real(gridptCount,8)
          intInfo%latIndexDepot(depotIndex)     = latIndexVec(ipoint)
          intInfo%lonIndexDepot(depotIndex)     = lonIndexVec(ipoint)
          depotIndex = depotIndex + 1

        end do

      end if

      numGridpt = gridptCount

    end subroutine getMyInterpWeightsFootprint

    !--------------------------------------------------------------------------
    ! getMyInterpWeightsLake (contained in getMyInterpWeights)
    !--------------------------------------------------------------------------
    subroutine getMyInterpWeightsLake(myHeaderIndex, varLevIndex, numGridpt)
      !
      !:Purpose: To determine the grid points and their associated weights
      !          for the lake horizontal interpolation.
      !
      implicit none

      ! Arguments:
      integer                , intent(in)    :: myHeaderIndex
      integer                , intent(in)    :: varLevIndex
      integer                , intent(out)   :: numGridpt

      ! Locals:
      integer :: depotIndex, latIndexCentre, lonIndexCentre, subGridIndex
      integer :: ipoint, gridptCount
      integer :: lakeCount, latIndexCurrent, lonIndexCurrent
      integer :: lonIndex, latIndex, lakeIndex
      integer :: lonIndexVec(statevector%ni*statevector%nj), latIndexVec(statevector%ni*statevector%nj)
      real(8) :: lonObs, latObs
      real(4) :: xpos_r4, ypos_r4
      logical :: lake(statevector%ni,statevector%nj)

      if ( stateVector%hco%grtyp == 'U' ) then
        call utl_abort('getMyInterpWeightsLake (s2c): Yin-Yang grid not supported')
      end if

      if ( .not.stateVector%oceanMask%maskPresent ) then
        call utl_abort('getMyInterpWeightsLake (s2c): Only compatible when mask present')
      end if

      numGridpt = 0

      ! Determine the grid point nearest the observation.

      latObs       = intInfo%myInterpObsLat(varLevIndex,myHeaderIndex)
      lonObs       = intInfo%myInterpObsLon(varLevIndex,myHeaderIndex)
      xpos_r4      = intInfo%myInterpObsXpos_r4(varLevIndex,myHeaderIndex)
      ypos_r4      = intInfo%myInterpObsYpos_r4(varLevIndex,myHeaderIndex)
      subGridIndex = intInfo%myInterpObsSubGridIndex(varLevIndex,myHeaderIndex)

      lonIndexCentre = nint(xpos_r4)
      latIndexCentre = nint(ypos_r4)

      if ( subGridIndex == 3 ) then
        call utl_abort('getMyInterpWeightsLake (s2c): two subGrids involved is not supported')
      end if

      gridptCount = 0

      ! It can happen that the lake location is closest to a grid point
      ! where MASK(I,J) = .false. while there are other grid points for the
      ! same lake where MASK(I,J) = .true.. Code needs modifications
      ! for this case.

      ! If observation is not on the grid, don't use it.
      if ( lonIndexCentre < 1 .or. lonIndexCentre > statevector%ni .or.  &
           latIndexCentre < 1 .or. latIndexCentre > statevector%nj ) return

      if ( .not. stateVector%oceanMask%mask(lonIndexCentre,latIndexCentre,1) ) return

      lake(:,:) = .false.
      lake(lonIndexCentre,latIndexCentre) = .true.
      gridptCount = 1
      lonIndexVec(gridptCount) = lonIndexCentre
      latIndexVec(gridptCount) = latIndexCentre

      lakeCount = 0

      do while(lakeCount /= gridptCount)

        do lakeIndex = lakeCount+1, gridptCount

          if(lakeIndex == lakeCount+1) lakeCount = gridptCount

          lonIndexCurrent = lonIndexVec(lakeIndex)
          latIndexCurrent = latIndexVec(lakeIndex)

          do latIndex = max(1,latIndexCurrent-1), min(latIndexCurrent+1,statevector%nj)
            do lonIndex = max(1,lonIndexCurrent-1), min(lonIndexCurrent+1,statevector%ni)
              if(stateVector%oceanMask%mask(lonIndex,latIndex,1) .and. .not. lake(lonIndex,latIndex)) then
                lake(lonIndex,latIndex) = .true.
                gridptCount = gridptCount + 1
                lonIndexVec(gridptCount) = lonIndex
                latIndexVec(gridptCount) = latIndex
              end if
            end do
          end do

        end do

      end do

      ! If the depot is allocated, then we fill in the values
      if ( allocated(intInfo%interpWeightDepot) ) then

        depotIndex = intInfo%depotIndexBeg(varLevIndex,myHeaderIndex)

        do ipoint=1,gridptCount

          intInfo%interpWeightDepot(depotIndex) = 1.0d0 / real(gridptCount,8)
          intInfo%latIndexDepot(depotIndex)     = latIndexVec(ipoint)
          intInfo%lonIndexDepot(depotIndex)     = lonIndexVec(ipoint)
          depotIndex = depotIndex + 1

        end do

      end if

      numGridpt = gridptCount

    end subroutine getMyInterpWeightsLake

    !--------------------------------------------------------------------------
    ! getMyInterpWeightsNearestNeighbor (contained in getMyInterpWeights)
    !--------------------------------------------------------------------------
    subroutine getMyInterpWeightsNearestNeighbor(myHeaderIndex, varLevIndex, &
                                                 numGridpt)
      !
      !:Purpose: Determine the nearest grid points to the observations location
      !
      implicit none

      ! Arguments:
      integer                , intent(in)    :: myHeaderIndex
      integer                , intent(in)    :: varLevIndex
      integer                , intent(out)   :: numGridpt

      ! Locals:
      integer :: depotIndex
      integer :: latIndex, lonIndex
      integer :: subGridIndex
      real(4) :: xpos_r4, ypos_r4

      if ( stateVector%hco%grtyp == 'U' ) then
        call utl_abort('getMyInterpWeightsNearestNeighbor (s2c): Yin-Yang grid not supported')
      end if

      numGridpt = 0

      xpos_r4      = intInfo%myInterpObsXpos_r4(varLevIndex,myHeaderIndex)
      ypos_r4      = intInfo%myInterpObsYpos_r4(varLevIndex,myHeaderIndex)

      latIndex = nint(ypos_r4)
      lonIndex = nint(xpos_r4)

      ! Handle periodicity in longitude
      if ( lonIndex == statevector%ni+1 .and. intInfo%periodic ) lonIndex = 1

      ! Test bounds
      if ( lonIndex < 1 .or. lonIndex > statevector%ni .or. &
           latIndex < 1 .or. latIndex > statevector%nj  ) then

        write(*,*) 'getMyInterpWeightsNearestNeighbor (s2c): observation out of bounds'
        write(*,*) 'lonIndex. latIndex = ', lonIndex, latIndex

      else

        ! If the depot is allocated, then we fill in the value
        if ( allocated(intInfo%interpWeightDepot) ) then

          depotIndex = intInfo%depotIndexBeg(varLevIndex,myHeaderIndex)

          intInfo%interpWeightDepot(depotIndex) = 1.0d0
          intInfo%latIndexDepot(depotIndex)     = latIndex
          intInfo%lonIndexDepot(depotIndex)     = lonIndex

        end if

        numGridpt = 1

      end if

    end subroutine getMyInterpWeightsNearestNeighbor

  end subroutine getMyInterpWeights

  !---------------------------------------------------------
  ! getMyInterpObsLatLon (called by setupInterpInfoTiles)
  !---------------------------------------------------------
  subroutine getMyInterpObsLatLon(intInfo, stateVector)
    !
    ! :Purpose: Define lat-lon, mpiIdSrc and headerIndex of observations on my
    !           lat-lon tile where the interpolation will be performed.
    !
    implicit none

    ! Arguments:
    type(struct_interpInfoTiles), intent(inout) :: intInfo       ! Interpolation info structure
    type(struct_gsv), target    , intent(in)    :: stateVector   ! stateVector object

    ! Locals:
    integer :: headerIndex, myHeaderIndex, yourNumHeader
    integer :: procIndex
    integer :: varLevIndex, numVarLev
    integer :: offset
    integer :: sendCounts(mmpi_nprocs), recvCounts(mmpi_nprocs)
    integer :: sendDispls(mmpi_nprocs), recvDispls(mmpi_nprocs)
    real(8), allocatable :: sendObsLat(:), sendObsLon(:), recvObsLat(:), recvObsLon(:)

    yourNumHeader = intInfo%yourNumHeader
    numVarLev = intInfo%numVarLevState

    ! Allocate myInterpObs in object: lat-lon, mpiIdSrc, headerIndex
    allocate(intInfo%myInterpObsLat(numVarLev,intInfo%myInterpNumHeader))
    allocate(intInfo%myInterpObsLon(numVarLev,intInfo%myInterpNumHeader))
    allocate(intInfo%myInterpObsMpiIdSrc(intInfo%myInterpNumHeader))
    allocate(intInfo%myInterpObsHeaderIndex(intInfo%myInterpNumHeader))

    ! Determine the number of observations each process will send
    sendCounts(:) = 0
    do headerIndex = 1, yourNumHeader
      procIndex = intInfo%allObsTileMpiId(headerIndex, mmpi_myid+1) + 1
      sendCounts(procIndex) = sendCounts(procIndex) + numVarLev
    end do

    ! Compute displacements for send buffers
    sendDispls(1) = 0
    do procIndex = 2, mmpi_nprocs
      sendDispls(procIndex) = sendDispls(procIndex-1) + sendCounts(procIndex-1)
    end do

    ! Allocate arrays needed for MPI communication
    allocate(sendObsLat(sum(sendCounts)))
    allocate(sendObsLon(sum(sendCounts)))

    sendCounts(:) = 0
    do headerIndex = 1, yourNumHeader
      procIndex = intInfo%allObsTileMpiId(headerIndex, mmpi_myid+1) + 1
      offset = sendDispls(procIndex) + sendCounts(procIndex)
      do varLevIndex = 1, numVarLev
        sendObsLat(offset + varLevIndex) = intInfo%yourObsLat(varLevIndex, headerIndex)
        sendObsLon(offset + varLevIndex) = intInfo%yourObsLon(varLevIndex, headerIndex)
      end do
      sendCounts(procIndex) = sendCounts(procIndex) + numVarLev
    end do

    ! Exchange observation counts
    call mmpi_alltoall(sendCounts, recvCounts)

    ! Compute displacements for receive buffers
    recvDispls(1) = 0
    do procIndex = 2, mmpi_nprocs
      recvDispls(procIndex) = recvDispls(procIndex-1) + recvCounts(procIndex-1)
    end do

    ! Allocate arrays needed for MPI communication
    allocate(recvObsLat(sum(recvCounts)))
    allocate(recvObsLon(sum(recvCounts)))

    ! Perform 'mmpi_alltoallv' communication
    call mmpi_alltoallv(sendObsLat, sendCounts, sendDispls, &
                        recvObsLat, recvCounts, recvDispls)
    call mmpi_alltoallv(sendObsLon, sendCounts, sendDispls, &
                        recvObsLon, recvCounts, recvDispls)

    ! Shuffle received lat-lon values into correct arrays
    recvCounts(:) = 0
    myHeaderIndex = 0
    do procIndex = 1, mmpi_nprocs
      do headerIndex = 1, intInfo%allNumHeader(procIndex)

        if (intInfo%allObsTileMpiId(headerIndex, procIndex) /= mmpi_myid) cycle
        myHeaderIndex = myHeaderIndex + 1

        intInfo%myInterpObsMpiIdSrc(myHeaderIndex) = procIndex - 1
        intInfo%myInterpObsHeaderIndex(myHeaderIndex) = headerIndex

        offset = recvDispls(procIndex) + recvCounts(procIndex)
        do varLevIndex = 1, numVarLev
          intInfo%myInterpObsLat(varLevIndex, myHeaderIndex) = recvObsLat(offset + varLevIndex)
          intInfo%myInterpObsLon(varLevIndex, myHeaderIndex) = recvObsLon(offset + varLevIndex)
        end do
        recvCounts(procIndex) = recvCounts(procIndex) + numVarLev

      end do
    end do

    ! Deallocate local arrays
    deallocate(sendObsLat)
    deallocate(sendObsLon)
    deallocate(recvObsLat)
    deallocate(recvObsLon)

  end subroutine getMyInterpObsLatLon

  !---------------------------------------------------------
  ! getMyInterpObsXYposSubGridIndex (called by setupInterpInfoTiles)
  !---------------------------------------------------------
  subroutine getMyInterpObsXYposSubGridIndex(intInfo, obsSpaceData, stateVector)
    !
    ! :Purpose: Define x/y position, subGridIndex and Footprint of observations on my
    !           lat-lon tile where the interpolation will be performed.
    !           Also store yourObsSubGridIndex for use elsewhere.
    !
    implicit none

    ! Arguments:
    type(struct_interpInfoTiles), intent(inout) :: intInfo       ! Interpolation info structure
    type(struct_obs)            , intent(inout) :: obsSpaceData  ! obs space object
    type(struct_gsv), target    , intent(in)    :: stateVector   ! stateVector object

    ! Locals:
    integer :: headerIndex, myHeaderIndex, myNumHeader, yourNumHeader
    integer :: procIndex
    integer :: varLevIndex, numVarLev
    integer :: offset, ierr
    integer :: sendCounts(mmpi_nprocs), recvCounts(mmpi_nprocs)
    integer :: sendDispls(mmpi_nprocs), recvDispls(mmpi_nprocs)
    real(8) :: lat, lon
    integer, allocatable :: sendSubGridIndex(:), recvSubGridIndex(:)
    real(4), allocatable :: lat_deg_r4(:), lon_deg_r4(:)
    real(4), allocatable :: xpos_r4(:,:), ypos_r4(:,:), xpos2_r4(:), ypos2_r4(:)
    real(4), allocatable :: sendXpos_r4(:), sendYpos_r4(:), recvXpos_r4(:), recvYpos_r4(:)
    real(4), allocatable :: sendFoot_r4(:), recvFoot_r4(:)

    yourNumHeader = intInfo%yourNumHeader
    myNumHeader = intInfo%myInterpNumHeader
    numVarLev = intInfo%numVarLevState

    ! Allocate myInterpObs in object: lat-lon, x/ypos, SubGridIndex, Footprint
    allocate(intInfo%myInterpObsXpos_r4(numVarLev,myNumHeader))
    allocate(intInfo%myInterpObsYpos_r4(numVarLev,myNumHeader))
    allocate(intInfo%myInterpObsSubGridIndex(numVarLev,myNumHeader))
    allocate(intInfo%myInterpObsFootprint_r4(numVarLev,myNumHeader))

    ! Allocate yourObsSubGridIndex needed for later (original obs distribution)
    allocate(intInfo%yourObsSubGridIndex(numVarLev,yourNumHeader))

    ! Determine the number of observations each process will send
    sendCounts(:) = 0
    do headerIndex = 1, yourNumHeader
      procIndex = intInfo%allObsTileMpiId(headerIndex, mmpi_myid+1) + 1
      sendCounts(procIndex) = sendCounts(procIndex) + numVarLev
    end do

    ! Compute displacements for send buffers
    sendDispls(1) = 0
    do procIndex = 2, mmpi_nprocs
      sendDispls(procIndex) = sendDispls(procIndex-1) + sendCounts(procIndex-1)
    end do

    sendCounts(:) = 0
    do headerIndex = 1, yourNumHeader
      procIndex = intInfo%allObsTileMpiId(headerIndex, mmpi_myid+1) + 1
      sendCounts(procIndex) = sendCounts(procIndex) + numVarLev
    end do

    ! Exchange observation counts
    call mmpi_alltoall(sendCounts, recvCounts)

    ! Compute displacements for receive buffers
    recvDispls(1) = 0
    do procIndex = 2, mmpi_nprocs
      recvDispls(procIndex) = recvDispls(procIndex-1) + recvCounts(procIndex-1)
    end do

    ! Now compute x/y position, subGridIndex

    allocate(lat_deg_r4(numVarLev))
    allocate(lon_deg_r4(numVarLev))
    allocate(xpos2_r4(numVarLev)) ! not used
    allocate(ypos2_r4(numVarLev)) ! not used
    allocate(xpos_r4(numVarLev,yourNumHeader))
    allocate(ypos_r4(numVarLev,yourNumHeader))

    ! Compute and define x/ypos, subGridIndex and rotated lat-lon
    do headerIndex = 1, yourNumHeader

      do varLevIndex = 1, numVarLev

        ! Compute x-y position and subGridIndex
        lat = intInfo%yourObsLat(varLevIndex, headerIndex)
        lon = intInfo%yourObsLon(varLevIndex, headerIndex)
        lat_deg_r4(varLevIndex) = real(lat * MPC_DEGREES_PER_RADIAN_R8, 4) ! Radian To Degree
        lon_deg_r4(varLevIndex) = real(lon * MPC_DEGREES_PER_RADIAN_R8, 4)

      end do

      ! Determine the xpos/ypos for this observation on the stateVector grid
      ierr = gpos_getPositionXY(stateVector%hco%EZscintID,   &              ! IN
                                xpos_r4(:,headerIndex),      &              ! OUT
                                ypos_r4(:,headerIndex),      &              ! OUT
                                xpos2_r4, ypos2_r4,          &              ! OUT
                                lat_deg_r4, lon_deg_r4,      &              ! IN
                                intInfo%yourObsSubGridIndex(:,headerIndex)) ! OUT

      do varLevIndex = 1, numVarLev
        ! Check returned value of subGridIndex
        if (intInfo%yourObsSubGridIndex(varLevIndex,headerIndex) /=1 .and.  &
            intInfo%yourObsSubGridIndex(varLevIndex,headerIndex) /=2) then
          call utl_abort('getMyInterpObs (s2c): invalid value of subGridIndex')
        end if

      end do

    end do ! myHeaderIndex

    ! Allocate arrays needed for MPI communication
    allocate(sendXpos_r4(sum(sendCounts)))
    allocate(sendYpos_r4(sum(sendCounts)))
    allocate(sendSubGridIndex(sum(sendCounts)))
    allocate(sendFoot_r4(sum(sendCounts)))
    allocate(recvXpos_r4(sum(recvCounts)))
    allocate(recvYpos_r4(sum(recvCounts)))
    allocate(recvSubGridIndex(sum(recvCounts)))
    allocate(recvFoot_r4(sum(recvCounts)))

    ! Shuffle rotated x/y position and subGridIndex arrays in preparation for mpi communication
    sendCounts(:) = 0
    do headerIndex = 1, yourNumHeader
      procIndex = intInfo%allObsTileMpiId(headerIndex, mmpi_myid+1) + 1
      offset = sendDispls(procIndex) + sendCounts(procIndex)
      do varLevIndex = 1, numVarLev
        sendXpos_r4(offset + varLevIndex)      = xpos_r4(varLevIndex, headerIndex)
        sendYpos_r4(offset + varLevIndex)      = ypos_r4(varLevIndex, headerIndex)
        sendSubGridIndex(offset + varLevIndex) = intInfo%yourObsSubGridIndex(varLevIndex, headerIndex)
        sendFoot_r4(offset + varLevIndex)      = s2c_getFootprintRadius( &
                                                              obsSpaceData, &
                                                              stateVector, headerIndex)
      end do
      sendCounts(procIndex) = sendCounts(procIndex) + numVarLev
    end do

    ! Perform 'mmpi_alltoallv' communication
    call mmpi_alltoallv(sendXpos_r4, sendCounts, sendDispls, &
                        recvXpos_r4, recvCounts, recvDispls)
    call mmpi_alltoallv(sendYpos_r4, sendCounts, sendDispls, &
                        recvYpos_r4, recvCounts, recvDispls)
    call mmpi_alltoallv(sendSubGridIndex, sendCounts, sendDispls, &
                        recvSubGridIndex, recvCounts, recvDispls)
    call mmpi_alltoallv(sendFoot_r4, sendCounts, sendDispls, &
                        recvFoot_r4, recvCounts, recvDispls)

    ! Shuffle received lat-lon values into correct arrays
    recvCounts(:) = 0
    myHeaderIndex = 0
    do procIndex = 1, mmpi_nprocs
      do headerIndex = 1, intInfo%allNumHeader(procIndex)

        if (intInfo%allObsTileMpiId(headerIndex, procIndex) /= mmpi_myid) cycle
        myHeaderIndex = myHeaderIndex + 1

        offset = recvDispls(procIndex) + recvCounts(procIndex)
        do varLevIndex = 1, numVarLev
          intInfo%myInterpObsXpos_r4(varLevIndex, myHeaderIndex) =  &
               recvXpos_r4(offset + varLevIndex)
          intInfo%myInterpObsYpos_r4(varLevIndex, myHeaderIndex) =  &
               recvYpos_r4(offset + varLevIndex)
          intInfo%myInterpObsSubGridIndex(varLevIndex, myHeaderIndex) =  &
               recvSubGridIndex(offset + varLevIndex)
          intInfo%myInterpObsFootprint_r4(varLevIndex, myHeaderIndex) =  &
               recvFoot_r4(offset + varLevIndex)
        end do
        recvCounts(procIndex) = recvCounts(procIndex) + numVarLev

      end do
    end do

    ! Deallocate local arrays
    deallocate(lat_deg_r4)
    deallocate(lon_deg_r4)
    deallocate(xpos_r4)
    deallocate(ypos_r4)
    deallocate(xpos2_r4)
    deallocate(ypos2_r4)

    deallocate(sendXpos_r4)
    deallocate(sendYpos_r4)
    deallocate(sendSubGridIndex)
    deallocate(sendFoot_r4)
    deallocate(recvXpos_r4)
    deallocate(recvYpos_r4)
    deallocate(recvSubGridIndex)
    deallocate(recvFoot_r4)

  end subroutine getMyInterpObsXYposSubGridIndex

  !---------------------------------------------------------
  ! getMyInterpObsRotLatLon (called by setupInterpInfoTiles)
  !---------------------------------------------------------
  subroutine getMyInterpObsRotLatLon(intInfo, stateVector)
    !
    ! :Purpose: Define rotated lat-lon (if needed) of observations on my
    !           lat-lon tile where the interpolation will be performed.
    !
    implicit none

    ! Arguments:
    type(struct_interpInfoTiles), intent(inout) :: intInfo       ! Interpolation info structure
    type(struct_gsv), target    , intent(in)    :: stateVector   ! stateVector object

    ! Locals:
    integer :: headerIndex, myHeaderIndex, yourNumHeader
    integer :: procIndex, offset
    integer :: varLevIndex, numVarLev
    integer :: sendCounts(mmpi_nprocs), recvCounts(mmpi_nprocs)
    integer :: sendDispls(mmpi_nprocs), recvDispls(mmpi_nprocs)
    real(8), allocatable :: latRot(:,:), lonRot(:,:)
    real(8), allocatable :: sendObsLat(:), sendObsLon(:), recvObsLat(:), recvObsLon(:)

    yourNumHeader = intInfo%yourNumHeader
    numVarLev = intInfo%numVarLevState

    ! Allocate myInterpObs in object: rotated lat-lon
    allocate(intInfo%myInterpObsLatRot(numVarLev,intInfo%myInterpNumHeader))
    allocate(intInfo%myInterpObsLonRot(numVarLev,intInfo%myInterpNumHeader))

    ! Determine the number of observations each process will send
    sendCounts(:) = 0
    do headerIndex = 1, yourNumHeader
      procIndex = intInfo%allObsTileMpiId(headerIndex, mmpi_myid+1) + 1
      sendCounts(procIndex) = sendCounts(procIndex) + numVarLev
    end do

    ! Compute displacements for send buffers
    sendDispls(1) = 0
    do procIndex = 2, mmpi_nprocs
      sendDispls(procIndex) = sendDispls(procIndex-1) + sendCounts(procIndex-1)
    end do

    sendCounts(:) = 0
    do headerIndex = 1, yourNumHeader
      procIndex = intInfo%allObsTileMpiId(headerIndex, mmpi_myid+1) + 1
      sendCounts(procIndex) = sendCounts(procIndex) + numVarLev
    end do

    ! Exchange observation counts
    call mmpi_alltoall(sendCounts, recvCounts)

    ! Compute displacements for receive buffers
    recvDispls(1) = 0
    do procIndex = 2, mmpi_nprocs
      recvDispls(procIndex) = recvDispls(procIndex-1) + recvCounts(procIndex-1)
    end do

    ! Allocate arrays needed for MPI communication
    allocate(sendObsLat(sum(sendCounts)))
    allocate(sendObsLon(sum(sendCounts)))
    allocate(recvObsLat(sum(recvCounts)))
    allocate(recvObsLon(sum(recvCounts)))

    ! Compute rotated lat/lon
    allocate(latRot(numVarLev,yourNumHeader))
    allocate(lonRot(numVarLev,yourNumHeader))

    ! Compute and define x/ypos, subGridIndex and rotated lat-lon
    do headerIndex = 1, yourNumHeader

      call uvr_RotateLatLonVec(intInfo%uvr,                                 & ! INOUT
                               intInfo%yourObsSubGridIndex(:,headerIndex),  & ! IN
                               latRot(:,headerIndex),                       & ! OUT (radians)
                               lonRot(:,headerIndex),                       & ! OUT (radians)
                               intInfo%yourObsLat(:,headerIndex),           & ! IN  (radians)
                               intInfo%yourObsLon(:,headerIndex),           & ! IN  (radians)
                               'ToLatLonRot')                                 ! IN

    end do ! myHeaderIndex

    ! Shuffle rotated lat/lon arrays in preparation for mpi communication
    sendCounts(:) = 0
    do headerIndex = 1, yourNumHeader
      procIndex = intInfo%allObsTileMpiId(headerIndex, mmpi_myid+1) + 1
      offset = sendDispls(procIndex) + sendCounts(procIndex)
      do varLevIndex = 1, numVarLev
        sendObsLat(offset + varLevIndex) = latRot(varLevIndex, headerIndex)
        sendObsLon(offset + varLevIndex) = lonRot(varLevIndex, headerIndex)
      end do
      sendCounts(procIndex) = sendCounts(procIndex) + numVarLev
    end do

    ! Perform 'mmpi_alltoallv' communication
    call mmpi_alltoallv(sendObsLat, sendCounts, sendDispls, &
                        recvObsLat, recvCounts, recvDispls)
    call mmpi_alltoallv(sendObsLon, sendCounts, sendDispls, &
                        recvObsLon, recvCounts, recvDispls)

    ! Shuffle received lat-lon values into correct arrays
    recvCounts(:) = 0
    myHeaderIndex = 0
    do procIndex = 1, mmpi_nprocs
      do headerIndex = 1, intInfo%allNumHeader(procIndex)

        if (intInfo%allObsTileMpiId(headerIndex, procIndex) /= mmpi_myid) cycle
        myHeaderIndex = myHeaderIndex + 1

        offset = recvDispls(procIndex) + recvCounts(procIndex)
        do varLevIndex = 1, numVarLev
          intInfo%myInterpObsLatRot(varLevIndex, myHeaderIndex) = recvObsLat(offset + varLevIndex)
          intInfo%myInterpObsLonRot(varLevIndex, myHeaderIndex) = recvObsLon(offset + varLevIndex)
        end do
        recvCounts(procIndex) = recvCounts(procIndex) + numVarLev

      end do
    end do

    ! Deallocate local arrays
    deallocate(sendObsLat)
    deallocate(sendObsLon)
    deallocate(recvObsLat)
    deallocate(recvObsLon)
    deallocate(latRot)
    deallocate(lonRot)

  end subroutine getMyInterpObsRotLatLon

  !---------------------------------------------------------
  ! getObsTileMpiId (called by setupInterpInfoTiles)
  !---------------------------------------------------------
  subroutine getObsTileMpiId(intInfo, obsSpaceData, stateVector)
    !
    ! :Purpose: Compute the mpi ID for each observation with respect to
    !           the lat-lon tiles - communicate globally to all MPI tasks.
    !           Note: the original lat-lon in obsSpaceData is used to
    !           determine the lat-lon tile mpi task and not the vertically
    !           varying lat-lon (slant path).
    !
    implicit none

    ! Arguments:
    type(struct_interpInfoTiles), intent(inout) :: intInfo      ! Interpolation info structure
    type(struct_obs),             intent(in)    :: obsSpaceData ! obs space object
    type(struct_gsv),    target,  intent(in)    :: stateVector  ! stateVector object

    ! Locals:
    integer :: numHeader, numHeaderMaxMpi, headerIndex, subGridIndex, procIndex
    integer :: obsTileMpiIdX, obsTileMpiIdY, ierr
    integer, allocatable :: obsTileMpiId(:)
    real(8) :: obsLat, obsLon
    real(4) :: lat_deg_r4, lon_deg_r4, xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    logical :: outsideObs

    ! Allocation
    numHeader = obs_numHeader(obsSpaceData)
    numHeaderMaxMpi = maxval(intInfo%allNumHeader(:))

    ! Allocate local arrays
    allocate(obsTileMpiId(numHeaderMaxMpi))

    ! Allocate arrays in structure
    allocate(intInfo%allObsTileMpiId(numHeaderMaxMpi,mmpi_nprocs))
    allocate(intInfo%yourObsTileMpiId(numHeader))
    intInfo%allObsTileMpiId(:,:) = -99
    intInfo%yourObsTileMpiId(:) = -99

    ! Loop over all obs local headers to determine the mpi tile of each obs location
    do headerIndex = 1, numHeader

      ! Get lat/lon in degrees (from obsSpaceData)
      obsLat = obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex)
      obsLon = obs_headElem_r(obsSpaceData, OBS_LON, headerIndex)
      if (obsLon <  0.0         ) then
        obsLon = obsLon + 2.0*MPC_PI_R4
      end if
      if (obsLon >= 2.*MPC_PI_R4) then
        obsLon = obsLon - 2.0*MPC_PI_R4
      end if
      lat_deg_r4 = real(obsLat * MPC_DEGREES_PER_RADIAN_R8) ! Radian To Degree
      lon_deg_r4 = real(obsLon * MPC_DEGREES_PER_RADIAN_R8)

      ! Determine the xpos/ypos for this observation on the stateVector grid
      ierr = gpos_getPositionXY(stateVector%hco%EZscintID,   &  ! IN
                                xpos_r4, ypos_r4,            &  ! OUT
                                xpos2_r4, ypos2_r4,          &  ! OUT
                                lat_deg_r4, lon_deg_r4,      &  ! IN
                                subGridIndex)                   ! OUT

      ! Check returned value of subGridIndex
      if (subGridIndex /=1 .and. subGridIndex /=2) then
        call utl_abort('getObsTileMpiId (s2c): invalid value of subGridIndex')
      end if

      ! Use findloc to find the tile id in X, Y directions and global tile id
      obsTileMpiIdX = utl_findloc(nint(xpos_r4) >= stateVector%allLonBeg(:) .and. &
                                  nint(xpos_r4) <= stateVector%allLonEnd(:)) - 1
      obsTileMpiIdY = utl_findloc(nint(ypos_r4) >= stateVector%allLatBeg(:) .and. &
                                  nint(ypos_r4) <= stateVector%allLatEnd(:)) - 1
      obsTileMpiId(headerIndex) = obsTileMpiIdX + obsTileMpiIdY*mmpi_npex

      ! Check values of mpiIdX/mpiIdY
      outsideObs = .false.
      if (obsTileMpiIdX < 0 .or. obsTileMpiIdX >= mmpi_npex) then
        outsideObs = .true.
      end if
      if (obsTileMpiIdY < 0 .or. obsTileMpiIdY >= mmpi_npey) then
        outsideObs = .true.
      end if

      if (outsideObs) then
        obsTileMpiId(headerIndex) = 0
      end if

      intInfo%yourObsTileMpiId(headerIndex) = obsTileMpiId(headerIndex)

    end do ! headerIndex

    ! Communicate tile MPI ID, lat-lon and X-Y positions to all MPI tasks
    call mmpi_allGather(obsTileMpiId, intInfo%allObsTileMpiId)
    call msg_memUsage('After allgather of obsTileMpiId')

    deallocate(obsTileMpiId)

    ! Compute myInterpNumHeader based on allObsTileMpiId
    intInfo%myInterpNumHeader = 0
    do procIndex = 1, mmpi_nprocs
      intInfo%myInterpNumHeader = intInfo%myInterpNumHeader + &
           count(intInfo%allObsTileMpiId(1:intInfo%allNumHeader(procIndex),procIndex) == mmpi_myid)
    end do

  end subroutine getObsTileMpiId

  !---------------------------------------------------------
  ! getObsLatLon (called by setupInterpInfoTiles)
  !---------------------------------------------------------
  subroutine getObsLatLon(intInfo, obsSpaceData, stateVector)
    !
    ! :Purpose: Compute the lat-lon for each observations and
    !           each level/variable, including for slant path.
    !
    implicit none

    ! Arguments:
    type(struct_interpInfoTiles), intent(inout) :: intInfo      ! Interpolation info structure
    type(struct_obs),             intent(inout) :: obsSpaceData ! obs space object
    type(struct_gsv),             intent(in)    :: stateVector  ! Reference stateVector object

    ! Locals:
    type(struct_gsv) :: stateVector3dHeights, stateVector3dHeights_mpiGlb
    integer :: numHeader, numVarLev, headerIndex, varLevIndex
    integer :: codeType, nlev_T, nlev_M, levIndex
    real(4) :: lon_r4, lat_r4
    real(4), pointer :: height3D_r4_ptr2(:,:,:)
    real(8), pointer :: height3D_r8_ptr1(:,:,:)
    real(4), pointer :: height3D_M_r4(:,:,:), height3D_T_r4(:,:,:)
    logical :: rejectOutsideObs
    logical :: doSlantPath, SlantTO, SlantRO, SlantRA, doSetup3dHeights
    real(8), allocatable :: latLev_T(:), lonLev_T(:), latLev_M(:), lonLev_M(:)
    real(8) :: latLev_S, lonLev_S
    logical :: firstHeaderSlantPathTO, firstHeaderSlantPathRO, firstHeaderSlantPathRA
    character(len=4)          :: varLevel

    if (intInfo%inputStateVectorType == 'nl' .and. rejectObsOutsideGlobalGrid) then
      rejectOutsideObs = .true.
    end if

    doSlantPath = .false.
    SlantTO     = .false.
    SlantRO     = .false.
    SlantRA     = .false.
    if (slantPath_TO_nl   .and. intInfo%inputStateVectorType == 'nl') then
      doSlantPath = .true.
      SlantTO     = .true.
    endif
    if (slantPath_TO_tlad .and. intInfo%inputStateVectorType /= 'nl') then
      doSlantPath = .true.
      SlantTO     = .true.
    endif
    if (slantPath_RO_nl   .and. intInfo%inputStateVectorType == 'nl') then
      doSlantPath = .true.
      SlantRO     = .true.
    endif
    if (slantPath_RA_nl   .and. intInfo%inputStateVectorType == 'nl') then
      doSlantPath = .true.
      SlantRA     = .true.
    endif
    write(*,*) 'getObsLatLon (s2c): doSlantPath, SlantTO, SlantRO, SlantRA = ', &
                                    doSlantPath, SlantTO, SlantRO, SlantRA

    doSetup3dHeights = doSlantPath .and.  &
                       gsv_varExist(stateVector,'Z_T') .and. &
                       gsv_varExist(stateVector,'Z_M')

    ! Extract 3D heights and make mpi global
    if (doSetup3dHeights) then

      ! Create state vector with only 3D heights on tiles
      if ( intInfo%inputStateVectorType == 'nl' ) then

        call gsv_allocate(stateVector3dHeights, 1, &
                          stateVector%hco, stateVector%vco, &
                          mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                          dataKind_opt=4, varNames_opt=['Z_M','Z_T'])

        call gsv_getField(stateVector,         height3D_r8_ptr1,'Z_T')
        call gsv_getField(stateVector3dHeights,height3D_r4_ptr2,'Z_T')
        height3D_r4_ptr2(:,:,:) = height3D_r8_ptr1(:,:,:)

        call gsv_getField(stateVector,         height3D_r8_ptr1,'Z_M')
        call gsv_getField(stateVector3dHeights,height3D_r4_ptr2,'Z_M')
        height3D_r4_ptr2(:,:,:) = height3D_r8_ptr1(:,:,:)

      else

        call utl_abort('getObsLatLon (s2c): only "nl" supported so far')

      end if

      ! Communicate 3D height fields onto all mpi tasks
      call gsv_allocate(stateVector3dHeights_mpiGlb, 1, &
                        stateVector%hco, stateVector%vco, &
                        mpi_local_opt=.false., &
                        dataKind_opt=4, varNames_opt=['Z_M','Z_T'])
      call gsv_transposeTilesToMpiGlobal(stateVector3dHeights_mpiGlb, stateVector3dHeights)
      call gsv_getField(stateVector3dHeights_mpiGlb,height3D_T_r4,'Z_T')
      call gsv_getField(stateVector3dHeights_mpiGlb,height3D_M_r4,'Z_M')

      write(*,*) 'getObsLatLon (s2c): height3D_T_r4='
      write(*,*) height3D_T_r4(1,1,:)
      write(*,*) 'getObsLatLon (s2c): height3D_M_r4='
      write(*,*) height3D_M_r4(1,1,:)

      call gsv_deallocate(stateVector3dHeights)

    end if ! doSetup3dHeights
    call msg_memUsage('After setup 3D heights')

    ! Compute lat-lon of each local obs on each varLevIndex

    ! Allocation
    numVarLev = intInfo%numVarLevState
    numHeader = obs_numHeader(obsSpaceData)
    allocate(intInfo%allNumHeader(mmpi_nprocs))
    call mmpi_allGather(numHeader, intInfo%allNumHeader)
    allocate(intInfo%yourObsLat(numVarLev,numHeader))
    allocate(intInfo%yourObsLon(numVarLev,numHeader))
    intInfo%yourObsLat(:,:) =0.0d0
    intInfo%yourObsLon(:,:) =0.0d0

    if (doSlantPath .and. &
        gsv_varExist(stateVector,'Z_T') .and. &
        gsv_varExist(stateVector,'Z_M')) then

      nlev_T = gsv_getNumLev(stateVector,'TH')
      nlev_M = gsv_getNumLev(stateVector,'MM')

      allocate(latLev_T(nlev_T))
      allocate(lonLev_T(nlev_T))
      allocate(latLev_M(nlev_M))
      allocate(lonLev_M(nlev_M))
      latLev_T(:) = 0.0d0
      lonLev_T(:) = 0.0d0
      latLev_M(:) = 0.0d0
      lonLev_M(:) = 0.0d0

      firstHeaderSlantPathTO = .true.
      firstHeaderSlantPathRO = .true.
      firstHeaderSlantPathRA = .true.
      header_loop: do headerIndex = 1, numHeader

        !- Get LatLon of observation location
        lat_r4 = real(obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex), 4)
        lon_r4 = real(obs_headElem_r(obsSpaceData, OBS_LON, headerIndex), 4)
        if (lon_r4 <  0.0            ) lon_r4 = lon_r4 + 2.0 * MPC_PI_R4
        if (lon_r4 >= 2.0 * MPC_PI_R4) lon_r4 = lon_r4 - 2.0 * MPC_PI_R4

        codeType = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)

        if (tvs_isIdBurpTovs(codeType) .and. SlantTO) then
          if (firstHeaderSlantPathTO) then
            write(*,*) 'getObsLatLon (s2c): start slant-path for TOVS. ', &
                 'numHeader = ',numHeader
            firstHeaderSlantPathTO = .false.
          end if

          ! calculate lat/lon along the line of sight
          call utl_tmg_start(32,'------s2c_Slant')
          call slp_calcLatLonTovs(obsSpaceData, stateVector%hco, headerIndex, & ! IN
                                  height3D_T_r4, height3D_M_r4,               & ! IN
                                  latLev_T, lonLev_T,                         & ! OUT
                                  latLev_M, lonLev_M,                         & ! OUT
                                  latLev_S, lonLev_S             )              ! OUT
          call utl_tmg_stop(32)

        else if (codeType == codtyp_get_codtyp('ro') .and. SlantRO ) then
          if (firstHeaderSlantPathRO) then
            write(*,*) 'getObsLatLon (s2c): start slant-path for RO. ', &
                 'numHeader = ',numHeader
            firstHeaderSlantPathRO = .false.
          end if

          ! Calculate lat/lon along the GPSRO obs
          call utl_tmg_start(32,'------s2c_Slant')
          call slp_calcLatLonRO(obsSpaceData, stateVector%hco, headerIndex, & ! IN
                                height3D_T_r4, height3D_M_r4,               & ! IN
                                latLev_T, lonLev_T,                         & ! OUT
                                latLev_M, lonLev_M,                         & ! OUT
                                latLev_S, lonLev_S                          ) ! OUT
          call utl_tmg_stop(32)
        else if (codeType == codtyp_get_codtyp('radar') .and. SlantRA ) then
          if ( firstHeaderSlantPathRA ) then
            write(*,*) 'getObsLatLon (s2c): start slant-path for RADAR. ', &
                 'numHeader=',numHeader
            firstHeaderSlantPathRA = .false.
          end if

          ! calculate lat/lon along the radar beam obs
          call slp_calcLatLonRadar(obsSpaceData, stateVector%hco, headerIndex, & ! IN
                                   height3D_T_r4, height3D_M_r4,                 & ! IN
                                   latLev_T, lonLev_T,                           & ! OUT
                                   latLev_M, lonLev_M,                           & ! OUT
                                   latLev_S, lonLev_S                           ) ! OUT
        else

          latLev_T(:) = real(lat_r4,8)
          lonLev_T(:) = real(lon_r4,8)
          latLev_M(:) = real(lat_r4,8)
          lonLev_M(:) = real(lon_r4,8)
          latLev_S = real(lat_r4,8)
          lonLev_S = real(lon_r4,8)

        end if

        ! check if the slanted lat/lon is inside the domain
        call latlonChecks (obsSpaceData, stateVector%hco, & ! IN
                           headerIndex, rejectOutsideObs, & ! IN
                           latLev_T, lonLev_T,            & ! IN/OUT
                           latLev_M, lonLev_M,            & ! IN/OUT
                           latLev_S, lonLev_S )             ! IN/OUT

        ! put the lat/lon from TH/MM levels to varLevIndex
        do varLevIndex = 1, numVarLev
          levIndex = gsv_getLevFromVarLev(stateVector,varLevIndex)
          varLevel = vnl_varLevelFromVarname(gsv_getVarNameFromVarLev(stateVector,varLevIndex))

          if ( varLevel == 'TH' ) then
            intInfo%yourObsLat(varLevIndex,headerIndex) = latLev_T(levIndex)
            intInfo%yourObsLon(varLevIndex,headerIndex) = lonLev_T(levIndex)
          else if ( varLevel == 'MM' ) then
            intInfo%yourObsLat(varLevIndex,headerIndex) = latLev_M(levIndex)
            intInfo%yourObsLon(varLevIndex,headerIndex) = lonLev_M(levIndex)
          else if ( varLevel == 'SF' ) then
            intInfo%yourObsLat(varLevIndex,headerIndex) = latLev_S
            intInfo%yourObsLon(varLevIndex,headerIndex) = lonLev_S
          else
            call utl_abort('getObsLatLon (s2c): unknown value of varLevel')
          end if

        end do

      end do header_loop

      write(*,*) 'getObsLatLon (s2c): min/max(yourObsLat) = ',  &
           minval(intInfo%yourObsLat)*MPC_DEGREES_PER_RADIAN_R8,  &
           maxval(intInfo%yourObsLat)*MPC_DEGREES_PER_RADIAN_R8
      write(*,*) 'getObsLatLon (s2c): min/max(yourObsLon) = ',  &
           minval(intInfo%yourObsLon)*MPC_DEGREES_PER_RADIAN_R8,  &
           maxval(intInfo%yourObsLon)*MPC_DEGREES_PER_RADIAN_R8

      call gsv_deallocate(stateVector3dHeights_mpiGlb)

    else ! doSlantPath

      do headerIndex = 1, numHeader

        ! Get obs lat/lon and copy to all levels
        do varLevIndex = 1, numVarLev
          intInfo%yourObsLat(varLevIndex,headerIndex) = obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex)
          intInfo%yourObsLon(varLevIndex,headerIndex) = obs_headElem_r(obsSpaceData, OBS_LON, headerIndex)
        end do

      end do

    end if ! doSlantPath

    ! Loop over all obs headers to adjust/check the lat-lon of each obs location
    do headerIndex = 1, numHeader

      ! Get obs lat/lon and copy to all levels
      do varLevIndex = 1, numVarLev

        if (intInfo%yourObsLon(varLevIndex,headerIndex) <  0.0         ) then
          intInfo%yourObsLon(varLevIndex,headerIndex) =  &
               intInfo%yourObsLon(varLevIndex,headerIndex) + 2.0*MPC_PI_R4
        end if
        if (intInfo%yourObsLon(varLevIndex,headerIndex) >= 2.*MPC_PI_R4) then
          intInfo%yourObsLon(varLevIndex,headerIndex) =  &
               intInfo%yourObsLon(varLevIndex,headerIndex) - 2.0*MPC_PI_R4
        end if
      end do

!      ! check if the lat/lon is inside the domain
!      call latlonChecks(obsSpaceData, stateVector%hco, & ! IN
!                        headerIndex, rejectOutsideObs, & ! IN
!                        latLev_T, lonLev_T,            & ! IN/OUT
!                        latLev_M, lonLev_M )             ! IN/OUT

    end do

  end subroutine getObsLatLon

  !---------------------------------------------------------
  ! setSfcVarLevIndex (called by setupInterpInfoTiles)
  !---------------------------------------------------------
  subroutine setSfcVarLevIndex(intInfo,stateVector)
    !
    ! :Purpose: Set the varLevIndex for specifying the lat-lon
    !           at the surface for GZsfc.
    !
    implicit none

    ! Arguments:
    type(struct_interpInfoTiles), intent(inout) :: intInfo      ! Interpolation info structure
    type(struct_gsv),             intent(in)    :: stateVector  ! Reference stateVector object

    ! Locals:
    integer :: varLevIndex

    if (stateVector%HeightSfcPresent) then
      if (gsv_varExist(stateVector,'P0')) then
        varLevLoop: do varLevIndex = 1, intInfo%numVarLevState

          ! Use varLevIndex for P0
          if (trim(gsv_getVarNameFromVarLev(stateVector,varLevIndex)) == 'P0') then
            intInfo%sfcVarLevIndex = varLevIndex
            exit varLevLoop
          end if

        end do varLevLoop

      else

        ! No P0, assume the last value of varLevIndex is at the surface
        intInfo%sfcVarLevIndex = intInfo%numVarLevState

      end if

    else

      ! No GZsfc, so probably don't need sfcVarLevIndex, set to 1
      intInfo%sfcVarLevIndex = 1

    end if

  end subroutine setSfcVarLevIndex

  !---------------------------------------------------------
  ! nlTiles (called by s2c_nl)
  !---------------------------------------------------------
  subroutine nlTiles(stateVector, obsSpaceData, column, timeInterpType, &
                     rejectOutsideObs, beSilent, dealloc_opt)
    !
    ! :Purpose: Non-linear version of the horizontal interpolation,
    !           used for a full field (usually the background state when computing
    !           the innovation vector).
    !
    !           This version uses the newer mpi strategy that leaves the
    !           gridded data mostly in the original latitude-longitude mpi
    !           distribution to minimize communication.
    !
    implicit none

    ! Arguments:
    type(struct_gsv),        intent(inout) :: stateVector       ! Input stateVector object
    type(struct_obs),        intent(inout) :: obsSpaceData      ! obs space data object
    type(struct_columnData), intent(inout) :: column            ! Output columnData object
    character(len=*),        intent(in)    :: timeInterpType    ! Type of temporal interpolation
    logical,                 intent(inout) :: rejectOutsideObs  ! Choose to reject obs outside domain
    logical,                 intent(in)    :: beSilent          ! Choose to print nothing
    logical,       optional, intent(in)    :: dealloc_opt

    ! Locals:
    logical :: dealloc
    type(real48_3d)       :: myHaloValues_r48
    type(real48_2d)       :: myColumnData

    call utl_tmg_start(36,'------s2c_barrier_NL')
    call mmpi_barrier
    call utl_tmg_stop(36)

    if (.not. beSilent) then
      write(*,*) 'nlTiles (s2c): Starting'
    end if

    if (present(dealloc_opt)) then
      dealloc = dealloc_opt
    else
      dealloc = .true.
    end if

    ! check the column and statevector have same nk/varNameList
    call checkColumnStatevectorMatch(column,statevector)

    ! calculate delP_T/delP_M and del Z_T/Z_M on the grid
    call gvt_transform(statevector, 'ZandP_nl')

    if (interpInfoTiles_nl%initialized) then
      if (.not. hco_equal(interpInfoTiles_nl%hco,stateVector%hco) .or.  &
          interpInfoTiles_nl%oti%numStep /= stateVector%numStep) then
        write(*,*) 'nlTiles: WARNING! Current hco grid parameters differ from allocated interpInfo!'
        write(*,*) 'nlTiles: InterpInfo will be deallocated.'
        call s2c_deallocInterpInfo(inputStateVectorType = 'nl')
      end if
    end if

    if (.not. interpInfoTiles_nl%initialized) then
      call utl_tmg_stop(34)
      call utl_tmg_start(31,'----s2c_Setups')

      ! distribute obs over the lat-lon tiles and set interpolation weights
      call setupInterpInfoTiles(interpInfoTiles_nl, obsSpaceData, stateVector, &
                                column, timeInterpType, rejectOutsideObs, &
                                inputStateVectorType = 'nl', beSilent=beSilent)

      call utl_tmg_stop(31)
      call utl_tmg_start(34,'----s2c_NL')
    else
      if (.not. beSilent) then
        write(*,*) 'nlTiles (s2c): using existing setup of interpInfoTiles_nl'
      end if
    end if

    ! Communicate the list of extra grid-points needed on each MPI tile
    call sendRecvHalo(stateVector, interpInfoTiles_nl, beSilent,  &
                      myHaloValues_r48)
    if (.not. beSilent) then
      call msg_memUsage('nlTiles')
    end if

    ! Do interpolation to compute myColumnData from stateVector and myHaloValues
    call utl_tmg_start(35,'------s2c_NL_Hinterp')
    call tileToColumn(stateVector, myHaloValues_r48, interpInfoTiles_nl,  &
                      myColumnData, beSilent)
    call utl_tmg_stop(35)

    ! Send columns to the original mpi tasks and put in column object
    call sendRecvColumns(myColumnData, column, interpInfoTiles_nl, beSilent)

    ! Impose limits on some variables (e.g. humidity)
    call imposeLimits(column, beSilent)

    ! Ensure pressure is monotonically increasing with level index in the columns
    if (slantPath_TO_nl) call pressureProfileMonotonicityCheck(obsSpaceData, column)

    ! Deallocate interpInfo structure, if requested
    if (dealloc) call s2c_deallocInterpInfo(inputStateVectorType='nl')

    ! Deallocate myHaloValues_r48
    if (allocated(myHaloValues_r48%r4)) deallocate(myHaloValues_r48%r4)
    if (allocated(myHaloValues_r48%r8)) deallocate(myHaloValues_r48%r8)
    if (allocated(myHaloValues_r48%GZsfc)) deallocate(myHaloValues_r48%GZsfc)

    ! Deallocate myColumnData
    if (allocated(myColumnData%r4)) deallocate(myColumnData%r4)
    if (allocated(myColumnData%r8)) deallocate(myColumnData%r8)
    if (allocated(myColumnData%GZsfc)) deallocate(myColumnData%GZsfc)

    if (.not. beSilent) then
      write(*,*) 'nlTiles (s2c): Finished'
    end if

  end subroutine nlTiles

  !---------------------------------------------------------
  ! tlTiles (called by s2c_tl)
  !---------------------------------------------------------
  subroutine tlTiles(stateVector, obsSpaceData, columnAnlInc, beSilent_opt)
    !
    ! :Purpose: Tangent-linear version of the horizontal interpolation,
    !           used for increments or perturbations.
    !
    implicit none

    ! Arguments:
    type(struct_gsv),           intent(inout) :: stateVector
    type(struct_obs),           intent(inout) :: obsSpaceData
    type(struct_columnData),    intent(inout) :: columnAnlInc
    logical,          optional, intent(in)    :: beSilent_opt

    ! Locals:
    logical              :: rejectOutsideObs, beSilent
    type(real48_3d)      :: myHaloValues_r48
    type(real48_2d)      :: myColumnData

    if (present(beSilent_opt)) then
      beSilent = beSilent_opt
    else
      beSilent = .false.
    end if

    if ( .not. interpInfoTiles_tlad%initialized ) then
      rejectOutsideObs = .false.
      call utl_tmg_stop(38)
      call utl_tmg_start(31,'----s2c_Setups')
      call setupInterpInfoTiles(interpInfoTiles_tlad, obsSpaceData, stateVector,  &
                                columnAnlInc, timeInterpType_tlad,  rejectOutsideObs, &
                                inputStateVectorType='tl', beSilent=beSilent)
      call utl_tmg_stop(31)
      call utl_tmg_start(38,'----s2c_TL')
    end if

    ! set contents of column to zero
    call col_zero(columnAnlInc)

    ! Communicate the list of extra grid-points needed on each MPI tile
    call sendRecvHalo(stateVector, interpInfoTiles_tlad, beSilent,  &
                      myHaloValues_r48)

    ! Do interpolation to compute myColumnData from stateVector and myHaloValues
    call utl_tmg_start(39,'------s2c_TL_Hinterp')
    call tileToColumn(stateVector, myHaloValues_r48, interpInfoTiles_tlad,  &
                      myColumnData, beSilent)
    call utl_tmg_stop(39)

    ! Send columns to the original mpi tasks and put in column object
    call sendRecvColumns(myColumnData, columnAnlInc, interpInfoTiles_tlad, beSilent)

  end subroutine tlTiles

  !---------------------------------------------------------
  ! adTiles (called by s2c_ad)
  !---------------------------------------------------------
  subroutine adTiles(stateVector, obsSpaceData, columnAnlInc, beSilent_opt)
    !
    ! :Purpose: Adjoint version of the horizontal interpolation,
    !           used for the cost function gradient with respect to the increment.
    !
    implicit none

    ! Arguments:
    type(struct_gsv),        intent(inout) :: stateVector
    type(struct_obs),        intent(inout) :: obsSpaceData
    type(struct_columnData), intent(inout) :: columnAnlInc
    logical,       optional, intent(in)    :: beSilent_opt

    ! Locals:
    logical              :: rejectOutsideObs, beSilent
    type(real48_3d)      :: myHaloValuesAd_r48
    type(real48_2d)      :: myColumnDataAd

    if (present(beSilent_opt)) then
      beSilent = beSilent_opt
    else
      beSilent = .false.
    end if

    if ( .not. interpInfoTiles_tlad%initialized ) then
      rejectOutsideObs = .false.
      call utl_tmg_stop(41)
      call utl_tmg_start(31,'----s2c_Setups')
      call setupInterpInfoTiles(interpInfoTiles_tlad, obsSpaceData, stateVector,  &
                                columnAnlInc, timeInterpType_tlad,  rejectOutsideObs, &
                                inputStateVectorType='tl', beSilent=beSilent)
      call utl_tmg_stop(31)
      call utl_tmg_start(41,'----s2c_AD')
    end if

    ! Set stateVector to zero
    call gsv_zero(stateVector)

    ! Send columns from column object on the original mpi tasks back to mpi task with tile
    call sendRecvColumnsAd(myColumnDataAd, columnAnlInc, interpInfoTiles_tlad, beSilent)

    ! Do adjoint of interpolation to update stateVector and myHaloValuesAd from myColumnDataAd
    call utl_tmg_start(42,'------s2c_AD_Hinterp')
    call tileToColumnAd(stateVector, myHaloValuesAd_r48, interpInfoTiles_tlad,  &
                        myColumnDataAd, beSilent)
    call utl_tmg_stop(42)

    ! Communicate the list of extra grid-points needed on each MPI tile
    call sendRecvHaloAd(stateVector, interpInfoTiles_tlad, beSilent,  &
                        myHaloValuesAd_r48)

  end subroutine adTiles

  !---------------------------------------------------------
  ! imposeLimits
  !---------------------------------------------------------
  subroutine imposeLimits(column, beSilent)
    implicit none

    ! Arguments:
    type(struct_columnData), intent(inout) :: column    ! columnData object
    logical,                 intent(in)    :: beSilent  ! Choose to print nothing

    ! Locals:
    integer          :: columnIndex, varNameIndex
    real(8), pointer :: column_ptr(:)

    if (.not. beSilent) write(*,*) 'imposeLimits (s2c): Starting'

    ! Impose a lower limit on HU
    if(col_varExist(column,'HU')) then
      do columnIndex = 1, col_getNumCol(column)
        column_ptr => col_getColumn(column, columnIndex, 'HU')
        column_ptr(:) = max(column_ptr(:), col_rhumin)
      end do
    end if

    ! Impose a lower/upper limits on ALL cloud variables
    do varNameIndex = 1, vnl_numvarmaxCloud
      if(.not. col_varExist(column, vnl_varNameListCloud(varNameIndex))) cycle
      do columnIndex = 1, col_getNumCol(column)
        column_ptr => col_getColumn(column, columnIndex, vnl_varNameListCloud(varNameIndex))
        column_ptr(:) = max(column_ptr(:), qlim_getMinValueCloud(vnl_varNameListCloud(varNameIndex)))
        column_ptr(:) = min(column_ptr(:), qlim_getMaxValueCloud(vnl_varNameListCloud(varNameIndex)))
      end do
    end do

    if (.not. beSilent) write(*,*) 'imposeLimits (s2c): Finished'

  end subroutine imposeLimits

  !---------------------------------------------------------
  ! sendRecvColumns
  !---------------------------------------------------------
  subroutine sendRecvColumns(myColumnData, column, intInfo, beSilent)
    !
    ! :Purpose: Send columns from lat-lon tiles back to original mpi tasks
    !           where the observations are located.
    !
    implicit none

    ! Arguments:
    type(real48_2d),              intent(in)    :: myColumnData  ! Column values for the lat-lon tiles
    type(struct_columnData),      intent(inout) :: column        ! columnData object
    type(struct_interpInfoTiles), intent(in)    :: intInfo       ! Interpolation info structure
    logical,                      intent(in)    :: beSilent      ! Choose to print nothing

    ! Locals:
    integer :: myNumCol, yourNumCol, numVarLev, varLevIndex, varLevIndexCol, columnIndex
    integer :: indexBeg, indexEnd, nsize, procSrc, procDest
    integer, allocatable :: sendCounts(:), recvCounts(:), sendDispls(:), recvDispls(:)
    real(8)              :: gzSfc
    real(8), allocatable :: sendBuffer(:), recvBuffer(:)
    real(8), pointer :: colPtr(:,:)

    if (.not. beSilent) write(*,*) 'sendRecvColumns (s2c): Starting'

    myNumCol   = col_getNumCol(column)
    yourNumCol = intInfo%myInterpNumHeader

    if (.not. beSilent) then
      write(*,*) 'sendRecvColumns (s2c): myNumCol, yourNumCol = ', myNumCol, yourNumCol
    end if

    HeightSfcPresent: if (allocated(myColumnData%GZsfc)) then
      nsize = 1

      ! Allocate send/recv metadata arrays
      allocate(sendCounts(mmpi_nprocs))
      allocate(recvCounts(mmpi_nprocs))
      allocate(sendDispls(mmpi_nprocs))
      allocate(recvDispls(mmpi_nprocs))
      sendCounts(:) = 0
      recvCounts(:) = 0

      ! Count how much data each process will send/receive
      do columnIndex = 1, myNumCol
        procSrc = intInfo%yourObsTileMpiId(columnIndex)
        recvCounts(procSrc+1) = recvCounts(procSrc+1) + nsize
      end do

      do columnIndex = 1, yourNumCol
        procDest = intInfo%myInterpObsMpiIdSrc(columnIndex)
        sendCounts(procDest+1) = sendCounts(procDest+1) + nsize
      end do

      sendDispls(1) = 0
      recvDispls(1) = 0
      do columnIndex = 2, mmpi_nprocs
        sendDispls(columnIndex) = sendDispls(columnIndex-1) +  &
                                       sendCounts(columnIndex-1)
        recvDispls(columnIndex) = recvDispls(columnIndex-1) +  &
                                       recvCounts(columnIndex-1)
      end do

      allocate(sendBuffer(sum(sendCounts)))
      allocate(recvBuffer(sum(recvCounts)))

      ! Fill the send buffer (need to reset counts)
      sendCounts(:) = 0
      do columnIndex = 1, yourNumCol
        ! Range of index values for this column
        procDest = intInfo%myInterpObsMpiIdSrc(columnIndex)
        indexBeg = 1 + sendDispls(procDest+1) + sendCounts(procDest+1)
        ! Copy this column into correct slot of the sendBuffer and increment count
        sendBuffer(indexBeg) = myColumnData%GZsfc(columnIndex)
        sendCounts(procDest+1) = sendCounts(procDest+1) + nsize
      end do

      ! Perform 'mmpi_alltoallv' communication for GZsfc
      call mmpi_alltoallv(sendBuffer, sendCounts, sendDispls, &
                          recvBuffer, recvCounts, recvDispls)

      ! Copy received GZsfc data (need to reset counts)
      recvCounts(:) = 0
      do columnIndex = 1, myNumCol
        ! Range of index values for this column
        procSrc = intInfo%yourObsTileMpiId(columnIndex)
        indexBeg = 1 + recvDispls(procSrc+1) + recvCounts(procSrc+1)
        ! Copy from correct slot of recvBuffer into column and increment count
        gzSfc = recvBuffer(indexBeg)
        call col_setHeightSfc(column, columnIndex, gzSfc)
        recvCounts(procSrc+1) = recvCounts(procSrc+1) + nsize
      end do

      ! Deallocate memory
      deallocate(sendCounts)
      deallocate(recvCounts)
      deallocate(sendDispls)
      deallocate(recvDispls)
      deallocate(sendBuffer)
      deallocate(recvBuffer)

    end if HeightSfcPresent

    numVarLev  = intInfo%numVarLevState
    nsize      = numVarLev

    ! Allocate send/recv metadata arrays
    allocate(sendCounts(mmpi_nprocs))
    allocate(recvCounts(mmpi_nprocs))
    allocate(sendDispls(mmpi_nprocs))
    allocate(recvDispls(mmpi_nprocs))
    sendCounts(:) = 0
    recvCounts(:) = 0
    sendDispls(:) = 0
    recvDispls(:) = 0

    ! Set the counts for send and receive
    do columnIndex = 1, myNumCol
      procSrc = intInfo%yourObsTileMpiId(columnIndex)
      recvCounts(procSrc+1) = recvCounts(procSrc+1) + nsize
    end do
    do columnIndex = 1, yourNumCol
      procDest = intInfo%myInterpObsMpiIdSrc(columnIndex)
      sendCounts(procDest+1) = sendCounts(procDest+1) + nsize
    end do

    ! Calculate displacements
    do columnIndex = 2, mmpi_nprocs
      sendDispls(columnIndex) = sendDispls(columnIndex-1) + sendCounts(columnIndex-1)
      recvDispls(columnIndex) = recvDispls(columnIndex-1) + recvCounts(columnIndex-1)
    end do

    ! Allocate send/recv buffers
    allocate(sendBuffer(sum(sendCounts)))
    allocate(recvBuffer(sum(recvCounts)))

    ! Fill the send buffer (need to reset counts)
    sendCounts(:) = 0
    do columnIndex = 1, yourNumCol
      ! Range of index values for this column
      procDest = intInfo%myInterpObsMpiIdSrc(columnIndex)
      indexBeg = 1 + sendDispls(procDest+1) + sendCounts(procDest+1)
      indexEnd = indexBeg + nsize - 1
      ! Copy this column into correct slot of the sendBuffer and increment count
      sendBuffer(indexBeg:indexEnd) = myColumnData%r8(:,columnIndex)
      sendCounts(procDest+1) = sendCounts(procDest+1) + nsize
    end do

    ! Perform 'mmpi_alltoallv' communication for main data
    call mmpi_alltoallv(sendBuffer, sendCounts, sendDispls, &
                        recvBuffer, recvCounts, recvDispls)

    ! Copy received data to output structure (need to reset counts
    colPtr => col_getAllColumns(column)
    recvCounts(:) = 0
    do columnIndex = 1, myNumCol
      ! Range of index values for this column
      procSrc = intInfo%yourObsTileMpiId(columnIndex)
      indexBeg = 1 + recvDispls(procSrc+1) + recvCounts(procSrc+1)
      ! Copy from correct slot of recvBuffer into column and increment count
      do varLevIndex = 1, numVarLev
        varLevIndexCol = intInfo%varLevColFromVarLevState(varLevIndex)
        colPtr(varLevIndexCol,columnIndex) = recvBuffer(indexBeg+varLevIndex-1)
      end do
      recvCounts(procSrc+1) = recvCounts(procSrc+1) + nsize
    end do

    ! Deallocate memory
    deallocate(sendCounts)
    deallocate(recvCounts)
    deallocate(sendDispls)
    deallocate(recvDispls)
    deallocate(sendBuffer)
    deallocate(recvBuffer)

    if (.not. beSilent) write(*,*) 'sendRecvColumns (s2c): Finished'

  end subroutine sendRecvColumns

  !---------------------------------------------------------
  ! sendRecvColumnsAd
  !---------------------------------------------------------
  subroutine sendRecvColumnsAd(myColumnDataAd, column, intInfo, beSilent)
    !
    ! :Purpose: Send columns from original mpi tasks where the observations
    !           are located to lat-lon tiles.
    !
    implicit none

    ! Arguments:
    type(real48_2d),              intent(inout) :: myColumnDataAd  ! Column values for the lat-lon tiles
    type(struct_columnData),      intent(inout) :: column          ! columnData object
    type(struct_interpInfoTiles), intent(in)    :: intInfo         ! Interpolation info structure
    logical,                      intent(in)    :: beSilent        ! Choose to print nothing

    ! Locals:
    integer :: myNumCol, yourNumCol, numVarLev, varLevIndex, varLevIndexCol, columnIndex
    integer :: indexBeg, indexEnd, nsize, procSrc, procDest
    integer, allocatable :: sendCounts(:), recvCounts(:), sendDispls(:), recvDispls(:)
    real(8)              :: gzSfc
    real(8), allocatable :: sendBuffer(:), recvBuffer(:)
    real(8), pointer :: colPtr(:,:)

    if (.not. beSilent) write(*,*) 'sendRecvColumnsAd (s2c): Starting'

    myNumCol   = col_getNumCol(column)
    yourNumCol = intInfo%myInterpNumHeader

    if (.not. beSilent) then
      write(*,*) 'sendRecvColumnsAd (s2c): myNumCol, yourNumCol = ', myNumCol, yourNumCol
    end if

    HeightSfcPresent: if (allocated(myColumnDataAd%GZsfc)) then
      nsize = 1

      if (.not. allocated(myColumnDataAd%GZsfc)) then
        allocate(myColumnDataAd%GZsfc(yourNumCol))
      end if
      myColumnDataAd%GZsfc(:) = 0.0d0

      ! Allocate send/recv metadata arrays
      allocate(sendCounts(mmpi_nprocs))
      allocate(recvCounts(mmpi_nprocs))
      allocate(sendDispls(mmpi_nprocs))
      allocate(recvDispls(mmpi_nprocs))
      sendCounts(:) = 0
      recvCounts(:) = 0

      ! Count how much data each process will send/receive
      do columnIndex = 1, myNumCol
        procDest = intInfo%yourObsTileMpiId(columnIndex)
        sendCounts(procSrc+1) = sendCounts(procDest+1) + nsize
      end do

      do columnIndex = 1, yourNumCol
        procSrc = intInfo%myInterpObsMpiIdSrc(columnIndex)
        recvCounts(procDest+1) = recvCounts(procSrc+1) + nsize
      end do

      sendDispls(1) = 0
      recvDispls(1) = 0
      do columnIndex = 2, mmpi_nprocs
        sendDispls(columnIndex) = sendDispls(columnIndex-1) +  &
                                       sendCounts(columnIndex-1)
        recvDispls(columnIndex) = recvDispls(columnIndex-1) +  &
                                       recvCounts(columnIndex-1)
      end do

      allocate(sendBuffer(sum(sendCounts)))
      allocate(recvBuffer(sum(recvCounts)))
      sendBuffer(:) = 0.0d0
      recvBuffer(:) = 0.0d0

      ! Fill the send buffer (need to reset counts)
      sendCounts(:) = 0
      do columnIndex = 1, myNumCol
        ! Range of index values for this column
        procDest = intInfo%yourObsTileMpiId(columnIndex)
        indexBeg = 1 + sendDispls(procDest+1) + sendCounts(procDest+1)
        ! Copy from column into correct slot of sendBuffer and increment count
        gzSfc = col_getHeight(column, 1, columnIndex, 'SF')
        sendBuffer(indexBeg) = gzSfc
        sendCounts(procSrc+1) = sendCounts(procDest+1) + nsize
      end do

      ! Perform 'mmpi_alltoallv' communication for GZsfc
      call mmpi_alltoallv(sendBuffer, sendCounts, sendDispls, &
                          recvBuffer, recvCounts, recvDispls)

      ! Copy received GZsfc data (need to reset counts)
      recvCounts(:) = 0
      do columnIndex = 1, yourNumCol
        ! Range of index values for this column
        procSrc = intInfo%myInterpObsMpiIdSrc(columnIndex)
        indexBeg = 1 + recvDispls(procDest+1) + recvCounts(procSrc+1)
        ! Copy from correct slot of recvBuffer into column and increment count
        myColumnDataAd%GZsfc(columnIndex) = recvBuffer(indexBeg)
        recvCounts(procDest+1) = recvCounts(procSrc+1) + nsize
      end do

      ! Deallocate memory
      deallocate(sendCounts)
      deallocate(recvCounts)
      deallocate(sendDispls)
      deallocate(recvDispls)
      deallocate(sendBuffer)
      deallocate(recvBuffer)

    end if HeightSfcPresent

    numVarLev  = intInfo%numVarLevState
    nsize      = numVarLev

    if (.not. allocated(myColumnDataAd%r8)) then
      allocate(myColumnDataAd%r8(numVarLev,yourNumCol))
    end if
    myColumnDataAd%r8(:,:) = 0.0d0

    ! Allocate send/recv metadata arrays
    allocate(sendCounts(mmpi_nprocs))
    allocate(recvCounts(mmpi_nprocs))
    allocate(sendDispls(mmpi_nprocs))
    allocate(recvDispls(mmpi_nprocs))
    sendCounts(:) = 0
    recvCounts(:) = 0
    sendDispls(:) = 0
    recvDispls(:) = 0

    ! Set the counts for send and receive
    do columnIndex = 1, myNumCol
      procDest = intInfo%yourObsTileMpiId(columnIndex)
      sendCounts(procDest+1) = sendCounts(procDest+1) + nsize
    end do
    do columnIndex = 1, yourNumCol
      procSrc = intInfo%myInterpObsMpiIdSrc(columnIndex)
      recvCounts(procSrc+1) = recvCounts(procSrc+1) + nsize
    end do

    ! Calculate displacements
    do columnIndex = 2, mmpi_nprocs
      sendDispls(columnIndex) = sendDispls(columnIndex-1) + sendCounts(columnIndex-1)
      recvDispls(columnIndex) = recvDispls(columnIndex-1) + recvCounts(columnIndex-1)
    end do

    ! Allocate send/recv buffers
    allocate(sendBuffer(sum(sendCounts)))
    allocate(recvBuffer(sum(recvCounts)))
    sendBuffer(:) = 0.0d0
    recvBuffer(:) = 0.0d0

    ! Fill the send buffer (need to reset counts)
    colPtr => col_getAllColumns(column)
    sendCounts(:) = 0
    do columnIndex = 1, myNumCol
      ! Range of index values for this column
      procDest = intInfo%yourObsTileMpiId(columnIndex)
      indexBeg = 1 + sendDispls(procDest+1) + sendCounts(procDest+1)
      ! Copy from column into correct slot of sendBuffer and increment count
      do varLevIndex = 1, numVarLev
        varLevIndexCol = intInfo%varLevColFromVarLevState(varLevIndex)
        sendBuffer(indexBeg+varLevIndex-1) = colPtr(varLevIndexCol,columnIndex)
      end do
      sendCounts(procDest+1) = sendCounts(procDest+1) + nsize
    end do

    ! Perform 'mmpi_alltoallv' communication for main data
    call mmpi_alltoallv(sendBuffer, sendCounts, sendDispls, &
                        recvBuffer, recvCounts, recvDispls)

    ! Copy received data to output structure (need to reset counts
    recvCounts(:) = 0
    do columnIndex = 1, yourNumCol
      ! Range of index values for this column
      procSrc = intInfo%myInterpObsMpiIdSrc(columnIndex)
      indexBeg = 1 + recvDispls(procSrc+1) + recvCounts(procSrc+1)
      indexEnd = indexBeg + nsize - 1
      ! Copy from correct slot of recvBuffer into column and increment count
      myColumnDataAd%r8(:,columnIndex) = recvBuffer(indexBeg:indexEnd)
      recvCounts(procSrc+1) = recvCounts(procSrc+1) + nsize
    end do

    ! Deallocate memory
    deallocate(sendCounts)
    deallocate(recvCounts)
    deallocate(sendDispls)
    deallocate(recvDispls)
    deallocate(sendBuffer)
    deallocate(recvBuffer)

    if (.not. beSilent) write(*,*) 'sendRecvColumnsAd (s2c): Finished'

  end subroutine sendRecvColumnsAd

  !---------------------------------------------------------
  ! tileToColumn
  !---------------------------------------------------------
  subroutine tileToColumn(stateVector, myHaloValues_r48, intInfo,  &
                          myColumnData, beSilent)
    !
    ! :Purpose: Interpolate horizontally and in time from 4D tile to 1D columns.
    !
    implicit none

    ! Arguments:
    type(struct_gsv)       ,        intent(in)    :: stateVector     ! stateVector object
    type(real48_3d),                intent(in)    :: myHaloValues_r48 ! Halo values
    type(struct_interpInfoTiles),   intent(in)    :: intInfo         ! Interpolation info structure
    type(real48_2d),                intent(inout) :: myColumnData    ! Column values for this lat-lon tile
    logical,                        intent(in)    :: beSilent        ! Choose to print nothing

    ! Locals:
    integer :: numColumn, numStep, numVarLev

    if (.not. beSilent) write(*,*) 'tileToColumn (s2c): Starting'

    numColumn = intInfo%myInterpNumHeader
    numStep   = stateVector%numStep
    numVarLev = intInfo%numVarLevState

    heightSfcPresent: if (stateVector%HeightSfcPresent) then

      if (.not. allocated(myColumnData%GZsfc)) then
        allocate(myColumnData%GZsfc(numColumn))
      end if
      myColumnData%GZsfc(:) = 0.0d0

      ! Interpolate horizontally (input r8, output r8)
      call hInterpGZ_nl(intInfo, myColumnData, &
                        stateVector, myHaloValues_r48)

    end if heightSfcPresent

    if (.not. allocated(myColumnData%r8)) then
      allocate(myColumnData%r8(numVarLev,numColumn))
    end if
    myColumnData%r8(:,:) = 0.0d0

    ! Interpolate horizontally (input r4/r8, output r8)
    call interp_nl(intInfo, myColumnData, &
                   stateVector, myHaloValues_r48)

    if (.not. beSilent) write(*,*) 'tileToColumn (s2c): Finished'

  end subroutine tileToColumn

  !---------------------------------------------------------
  ! tileToColumnAd
  !---------------------------------------------------------
  subroutine tileToColumnAd(stateVector, myHaloValuesAd_r48, intInfo,  &
                            myColumnDataAd, beSilent)
    !
    ! :Purpose: Adjoint of interpolate horizontally and in time from 4D tile to 1D columns.
    !
    implicit none

    ! Arguments:
    type(struct_gsv),             intent(inout) :: stateVector       ! stateVector object
    type(real48_3d),              intent(inout) :: myHaloValuesAd_r48 ! Halo values
    type(struct_interpInfoTiles), intent(in)    :: intInfo           ! Interpolation info structure
    type(real48_2d),              intent(inout) :: myColumnDataAd    ! Column values for this lat-lon tile
    logical,                      intent(in)    :: beSilent          ! Choose to print nothing

    ! Locals:
    integer :: numColumn, numStep, numVarLev, numHalo

    if (.not. beSilent) write(*,*) 'tileToColumnAd (s2c): Starting'

    numColumn = intInfo%myInterpNumHeader
    numHalo   = intInfo%myHaloSize
    numStep   = stateVector%numStep
    numVarLev = intInfo%numVarLevState

    heightSfcPresent: if (stateVector%HeightSfcPresent) then

      if (.not. allocated(myHaloValuesAd_r48%GZsfc)) then
        allocate(myHaloValuesAd_r48%GZsfc(numHalo))
      end if
      myHaloValuesAd_r48%GZsfc(:) = 0.0d0

      ! Interpolate horizontally (input r8, output r8)
      call hInterpGZ_ad(intInfo, myColumnDataAd, &
                        stateVector, myHaloValuesAd_r48)

    end if heightSfcPresent

    if (stateVector%dataKind == 4) then
      if (.not. allocated(myHaloValuesAd_r48%r4)) then
        myHaloValuesAd_r48%dataKind = 4
        allocate(myHaloValuesAd_r48%r4(numVarLev,numStep,numHalo))
      end if
      myHaloValuesAd_r48%r4(:,:,:) = 0.0
    else
      if (.not. allocated(myHaloValuesAd_r48%r8)) then
        myHaloValuesAd_r48%dataKind = 8
        allocate(myHaloValuesAd_r48%r8(numVarLev,numStep,numHalo))
      end if
      myHaloValuesAd_r48%r8(:,:,:) = 0.0d0
    end if

    ! Interpolate horizontally (input r4/r8, output r8)
    call interp_ad(intInfo, myColumnDataAd, &
                   stateVector, myHaloValuesAd_r48)

    if (.not. beSilent) write(*,*) 'tileToColumnAd (s2c): Finished'

  end subroutine tileToColumnAd

  !---------------------------------------------------------
  ! interp_nl:
  !---------------------------------------------------------
  subroutine interp_nl(intInfo, myColumnData, stateVector, myHaloValues_r48)
    !
    ! :Purpose: Horizontal and temporal interpolation of both scalar
    !           and vector fields.
    !
    implicit none

    ! Arguments:
    type(struct_interpInfoTiles), intent(in)    :: intInfo          ! Interpolation info structure
    type(real48_2d),              intent(inout) :: myColumnData     ! Output column data
    type(struct_gsv),             intent(in)    :: stateVector      ! stateVector object
    type(real48_3d),              intent(in)    :: myHaloValues_r48 ! Halo values

    ! Locals:
    integer :: columnIndex, headerIndex, procIndex, numStep, stepIndex
    integer :: numColumn, numVarLev, lonIndex, latIndex, haloIndex
    integer :: varLevIndex, levIndex, depotIndex
    integer :: offsetUU, offsetVV, levIndexUU, levIndexVV, subGridIndex
    real(8) :: timeWeight, lat, lon, latRot, lonRot, gridWeight
    real(4), pointer :: ptr4d_r4(:,:,:,:)
    real(8), pointer :: ptr4d_r8(:,:,:,:)

    numColumn = intInfo%myInterpNumHeader
    numStep   = stateVector%numStep
    numVarLev = intInfo%numVarLevState

    if (stateVector%dataKind == 4) then
      call gsv_getField(stateVector, ptr4d_r4)
    else
      call gsv_getField(stateVector, ptr4d_r8)
    end if

    !$OMP PARALLEL DO PRIVATE (columnIndex, stepIndex, headerIndex, procIndex,  &
    !$OMP timeWeight, varLevIndex, depotIndex, lonIndex, latIndex, gridWeight, haloIndex)
    do columnIndex = 1, numColumn

      ! For this columnIndex, interpolate gridded state to obs location/time
      myColumnData%r8(:,columnIndex) = 0.0d0

      stepLoop: do stepIndex = 1, numStep

        ! Extract the temporal interpolation weight
        headerIndex = intInfo%myInterpObsHeaderIndex(columnIndex)
        procIndex   = intInfo%myInterpObsMpiIdSrc(columnIndex) + 1
        timeWeight = oti_getTimeInterpWeightMpiGlobal(intInfo%oti, headerIndex,  &
                                                      stepIndex,procIndex)

        ! skip this time step if weight is zero
        if ( utl_isEqual(timeWeight, 0.0d0) ) cycle

        do varLevIndex = 1, numVarLev
          do depotIndex = intInfo%depotIndexBeg(varLevIndex,columnIndex), &
                          intInfo%depotIndexEnd(varLevIndex,columnIndex)

            lonIndex   = intInfo%lonIndexDepot(depotIndex)
            latIndex   = intInfo%latIndexDepot(depotIndex)
            gridWeight = intInfo%interpWeightDepot(depotIndex)

            if (stateVector%dataKind == 4) then

              if (latIndex == valueIsInHalo) then
                ! This location is in the halo
                haloIndex = lonIndex
                myColumnData%r8(varLevIndex,columnIndex) =  &
                     myColumnData%r8(varLevIndex,columnIndex) +  &
                     timeWeight * gridWeight * &
                     real(myHaloValues_r48%r4(varLevIndex,stepIndex,haloIndex),8)
              else
                ! This location is in the tile interior
                myColumnData%r8(varLevIndex,columnIndex) =  &
                     myColumnData%r8(varLevIndex,columnIndex) +  &
                     timeWeight * gridWeight * &
                     real(ptr4d_r4(lonIndex, latIndex, varLevIndex, stepIndex),8)
              end if

            else

              if (latIndex == valueIsInHalo) then
                ! This location is in the halo
                haloIndex = lonIndex
                myColumnData%r8(varLevIndex,columnIndex) =  &
                     myColumnData%r8(varLevIndex,columnIndex) +  &
                     timeWeight * gridWeight * &
                     myHaloValues_r48%r8(varLevIndex,stepIndex,haloIndex)
              else
                ! This location is in the tile interior
                myColumnData%r8(varLevIndex,columnIndex) =  &
                     myColumnData%r8(varLevIndex,columnIndex) +  &
                     timeWeight * gridWeight * &
                     ptr4d_r8(lonIndex, latIndex, varLevIndex, stepIndex)
              end if

            end if

          end do ! depotIndex
        end do ! varLevIndex

      end do stepLoop

    end do ! columnIndex
    !$OMP END PARALLEL DO

    ! Rotate any vector variables
    if (intInfo%rotatedWinds) then
      offsetUU = gsv_getOffsetFromVarName(statevector,'UU')
      offsetVV = gsv_getOffsetFromVarName(statevector,'VV')

      !$OMP PARALLEL DO PRIVATE (columnIndex, levIndex, levIndexUU, levIndexVV, lat, lon, latRot, lonRot, subGridIndex)
      do columnIndex = 1, numColumn
        do levIndex = 1, gsv_getNumLev(statevector,'MM')
          levIndexUU = offsetUU + levIndex
          levIndexVV = offsetVV + levIndex

          lat          = intInfo%myInterpObsLat(levIndexUU,columnIndex)
          lon          = intInfo%myInterpObsLon(levIndexUU,columnIndex)
          latRot       = intInfo%myInterpObsLatRot(levIndexUU,columnIndex)
          lonRot       = intInfo%myInterpObsLonRot(levIndexUU,columnIndex)
          subGridIndex = intInfo%myInterpObsSubGridIndex(levIndexUU,columnIndex)

          call uvr_rotateWind_nl(intInfo%uvr,                              & ! IN
                                 subGridIndex,                             & ! IN
                                 myColumnData%r8(levIndexUU,columnIndex),  & ! INOUT
                                 myColumnData%r8(levIndexVV,columnIndex),  & ! INOUT
                                 lat, lon, latRot, lonRot,                 & ! IN
                                 'ToMetWind' )                               ! IN

        end do ! levIndex
      end do ! columnIndex
      !$OMP END PARALLEL DO

    end if

  end subroutine interp_nl

  !---------------------------------------------------------
  ! interp_ad:
  !---------------------------------------------------------
  subroutine interp_ad(intInfo, myColumnDataAd, stateVector, myHaloValuesAd_r48)
    !
    ! :Purpose: Adjoint of horizontal and temporal interpolation of both scalar
    !           and vector fields.
    !
    implicit none

    ! Arguments:
    type(struct_interpInfoTiles), intent(in)    :: intInfo            ! Interpolation info structure
    type(real48_2d),              intent(inout) :: myColumnDataAd     ! Output column data
    type(struct_gsv),             intent(inout) :: stateVector        ! stateVector object
    type(real48_3d),              intent(inout) :: myHaloValuesAd_r48 ! Halo values

    ! Locals:
    integer :: columnIndex, headerIndex, procIndex, numStep, stepIndex
    integer :: numColumn, numVarLev, lonIndex, latIndex, haloIndex
    integer :: varLevIndex, levIndex, depotIndex
    integer :: offsetUU, offsetVV, levIndexUU, levIndexVV, subGridIndex
    real(8) :: timeWeight, lat, lon, latRot, lonRot, gridWeight
    real(4), pointer :: ptr4d_r4(:,:,:,:)
    real(8), pointer :: ptr4d_r8(:,:,:,:)

    numColumn = intInfo%myInterpNumHeader
    numStep   = stateVector%numStep
    numVarLev = intInfo%numVarLevState

    ! Adjoint of rotate any vector variables
    if (intInfo%rotatedWinds) then
      offsetUU = gsv_getOffsetFromVarName(statevector,'UU')
      offsetVV = gsv_getOffsetFromVarName(statevector,'VV')

      !$OMP PARALLEL DO PRIVATE (columnIndex, levIndex, levIndexUU, levIndexVV, lat, lon, latRot, lonRot, subGridIndex)
      do columnIndex = 1, numColumn
        do levIndex = 1, gsv_getNumLev(statevector,'MM')
          levIndexUU = offsetUU + levIndex
          levIndexVV = offsetVV + levIndex

          lat          = intInfo%myInterpObsLat(levIndexUU,columnIndex)
          lon          = intInfo%myInterpObsLon(levIndexUU,columnIndex)
          latRot       = intInfo%myInterpObsLatRot(levIndexUU,columnIndex)
          lonRot       = intInfo%myInterpObsLonRot(levIndexUU,columnIndex)
          subGridIndex = intInfo%myInterpObsSubGridIndex(levIndexUU,columnIndex)

          call uvr_rotateWind_ad(intInfo%uvr,                                & ! IN
                                 subGridIndex,                               & ! IN
                                 myColumnDataAd%r8(levIndexUU,columnIndex),  & ! INOUT
                                 myColumnDataAd%r8(levIndexVV,columnIndex),  & ! INOUT
                                 lat, lon, latRot, lonRot,                   & ! IN
                                 'ToMetWind' )                                 ! IN

        end do ! levIndex
      end do ! columnIndex
      !$OMP END PARALLEL DO

    end if

    if (stateVector%dataKind == 4) then
      call gsv_getField(stateVector, ptr4d_r4)
    else
      call gsv_getField(stateVector, ptr4d_r8)
    end if

    !$OMP PARALLEL DO PRIVATE (stepIndex, columnIndex, headerIndex, procIndex,  &
    !$OMP timeWeight, varLevIndex, depotIndex, lonIndex, latIndex, gridWeight, haloIndex)
    stepLoop: do stepIndex = 1, numStep

      ! For this columnIndex, adjoint of interpolate gridded state to obs location/time
      do columnIndex = 1, numColumn

        ! Extract the temporal interpolation weight
        headerIndex = intInfo%myInterpObsHeaderIndex(columnIndex)
        procIndex   = intInfo%myInterpObsMpiIdSrc(columnIndex) + 1
        timeWeight = oti_getTimeInterpWeightMpiGlobal(intInfo%oti, headerIndex,  &
                                                      stepIndex,procIndex)

        ! skip this time step if weight is zero
        if ( utl_isEqual(timeWeight, 0.0d0) ) cycle

        do varLevIndex = 1, numVarLev
          do depotIndex = intInfo%depotIndexBeg(varLevIndex,columnIndex), &
                          intInfo%depotIndexEnd(varLevIndex,columnIndex)

            lonIndex   = intInfo%lonIndexDepot(depotIndex)
            latIndex   = intInfo%latIndexDepot(depotIndex)
            gridWeight = intInfo%interpWeightDepot(depotIndex)

            if (stateVector%dataKind == 4) then

              if (latIndex == valueIsInHalo) then
                ! This location is in the halo
                haloIndex = lonIndex
                myHaloValuesAd_r48%r4(varLevIndex,stepIndex,haloIndex) =  &
                     real(myHaloValuesAd_r48%r4(varLevIndex,stepIndex,haloIndex),8) +  &
                     timeWeight * gridWeight *  &
                     myColumnDataAd%r8(varLevIndex,columnIndex)
              else
                ! This location is in the tile interior
                ptr4d_r4(lonIndex, latIndex, varLevIndex, stepIndex) =  &
                     real(ptr4d_r4(lonIndex, latIndex, varLevIndex, stepIndex),8) +  &
                     timeWeight * gridWeight * &
                     myColumnDataAd%r8(varLevIndex,columnIndex)
              end if

            else

              if (latIndex == valueIsInHalo) then
                ! This location is in the halo
                haloIndex = lonIndex
                myHaloValuesAd_r48%r8(varLevIndex,stepIndex,haloIndex) =  &
                     myHaloValuesAd_r48%r8(varLevIndex,stepIndex,haloIndex) +  &
                     timeWeight * gridWeight *  &
                     myColumnDataAd%r8(varLevIndex,columnIndex)
              else
                ! This location is in the tile interior
                ptr4d_r8(lonIndex, latIndex, varLevIndex, stepIndex) =  &
                     ptr4d_r8(lonIndex, latIndex, varLevIndex, stepIndex) +  &
                     timeWeight * gridWeight * &
                     myColumnDataAd%r8(varLevIndex,columnIndex)
              end if

            end if

          end do ! depotIndex
        end do ! varLevIndex

      end do ! columnIndex

    end do stepLoop
    !$OMP END PARALLEL DO

  end subroutine interp_ad

  !---------------------------------------------------------
  ! hInterpGZ_nl:
  !---------------------------------------------------------
  subroutine hInterpGZ_nl(intInfo, myColumnData, stateVector, myHaloValues_r48)
    !
    ! :Purpose: Scalar horizontal interpolation of GZsfc field.
    !
    implicit none

    ! Arguments:
    type(struct_interpInfoTiles), intent(in)    :: intInfo          ! Interpolation info structure
    type(real48_2d),              intent(inout) :: myColumnData     ! Output column data
    type(struct_gsv),             intent(in)    :: stateVector      ! stateVector object
    type(real48_3d),              intent(in)    :: myHaloValues_r48 ! Halo values

    ! Locals:
    integer :: numColumn, columnIndex, lonIndex, latIndex, depotIndex, varLevIndex, haloIndex
    real(8) :: interpValue, weight
    real(8), pointer :: ptr2d_r8(:,:)

    numColumn = intInfo%myInterpNumHeader
    varLevIndex = intInfo%sfcVarLevIndex
    ptr2d_r8 => gsv_getHeightSfc(stateVector)

    do columnIndex = 1, numColumn

      ! Interpolate the model state to the obs point
      interpValue = 0.0d0

      do depotIndex = intInfo%depotIndexBeg(varLevIndex,columnIndex), &
                      intInfo%depotIndexEnd(varLevIndex,columnIndex)

        lonIndex = intInfo%lonIndexDepot(depotIndex)
        latIndex = intInfo%latIndexDepot(depotIndex)
        weight   = intInfo%interpWeightDepot(depotIndex)

        if (latIndex == valueIsInHalo) then
          ! This location is in the halo
          haloIndex = lonIndex
          interpValue = interpValue + weight *  &
                        myHaloValues_r48%GZsfc(haloIndex)
        else
          ! This location is in the tile interior
          interpValue = interpValue + weight *  &
                        ptr2d_r8(lonIndex, latIndex)
        end if

      end do ! depotIndex
      myColumnData%GZsfc(columnIndex) = interpValue

    end do ! columnIndex

  end subroutine hInterpGZ_nl

  !---------------------------------------------------------
  ! hInterpGZ_ad:
  !---------------------------------------------------------
  subroutine hInterpGZ_ad(intInfo, myColumnDataAd, stateVector, myHaloValuesAd_r48)
    !
    ! :Purpose: Scalar horizontal interpolation of GZsfc field.
    !
    implicit none

    ! Arguments:
    type(struct_interpInfoTiles), intent(in)    :: intInfo            ! Interpolation info structure
    type(real48_2d),              intent(inout) :: myColumnDataAd     ! Output column data
    type(struct_gsv),             intent(inout) :: stateVector        ! stateVector object
    type(real48_3d),              intent(inout) :: myHaloValuesAd_r48 ! Halo values

    ! Locals:
    integer :: numColumn, columnIndex, lonIndex, latIndex, depotIndex, varLevIndex, haloIndex
    real(8) :: weight
    real(8), pointer :: ptr2d_r8(:,:)

    numColumn = intInfo%myInterpNumHeader
    varLevIndex = intInfo%sfcVarLevIndex
    ptr2d_r8 => gsv_getHeightSfc(stateVector)

    do columnIndex = 1, numColumn

      ! Adjoint of interpolate the model state to the obs point

      do depotIndex = intInfo%depotIndexBeg(varLevIndex,columnIndex), &
                      intInfo%depotIndexEnd(varLevIndex,columnIndex)

        lonIndex = intInfo%lonIndexDepot(depotIndex)
        latIndex = intInfo%latIndexDepot(depotIndex)
        weight   = intInfo%interpWeightDepot(depotIndex)

        if (latIndex == valueIsInHalo) then
          ! This location is in the halo
          haloIndex = lonIndex
          myHaloValuesAd_r48%GZsfc(haloIndex) = myHaloValuesAd_r48%GZsfc(haloIndex) +  &
                                                weight * myColumnDataAd%GZsfc(columnIndex)
        else
          ! This location is in the tile interior
          ptr2d_r8(lonIndex, latIndex) = ptr2d_r8(lonIndex, latIndex) +  &
                                         weight * myColumnDataAd%GZsfc(columnIndex)
        end if

      end do ! depotIndex

    end do ! columnIndex

  end subroutine hInterpGZ_ad

  !---------------------------------------------------------
  ! sendRecvHalo
  !---------------------------------------------------------
  subroutine sendRecvHalo(stateVector, intInfo, beSilent, &
                          myHaloValues_r48)
    !
    ! :Purpose: Exchange data between mpi tasks to provide a list of the
    !           data values within the halo.
    !
    implicit none

    ! Arguments:
    type(struct_gsv)       ,        intent(in)  :: stateVector                ! stateVector object
    type(struct_interpInfoTiles),   intent(in)  :: intInfo                    ! Interpolation info structure
    logical,                        intent(in)  :: beSilent                   ! Choose to print nothing
    type(real48_3d),                intent(out) :: myHaloValues_r48           ! Halo values

    ! Locals:
    integer :: haloIndex, numRecv, numSend, recvTag, sendTag, nsize, ierr
    integer :: lonIndex, latIndex
    integer :: procSrc, procDest
    type(mpi_request), allocatable :: requestIdRecv(:), requestIdSend(:)
    real(4), allocatable :: haloSend_r4(:,:,:)
    real(8), allocatable :: haloSend_r8(:,:,:)
    real(8), allocatable :: haloSendGZsfc(:)
    real(4), pointer     :: ptr4d_r4(:,:,:,:)
    real(8), pointer     :: ptr4d_r8(:,:,:,:)
    real(8), pointer     :: ptr2d_r8(:,:)

    if (.not. beSilent) write(*,*) 'sendRecvHalo (s2c): Starting'

    allocate(requestIdSend(intInfo%yourHaloSize))
    allocate(requestIdRecv(intInfo%myHaloSize))

    if (.not. beSilent) then
      write(*,*) 'sendRecvHalo (s2c): myHaloSize, yourHaloSize = ', intInfo%myHaloSize, intInfo%yourHaloSize
    end if

    heightSfcPresent: if (stateVector%HeightSfcPresent) then

      ! Allocate halo and other arrays
      allocate(myHaloValues_r48%GZsfc(intInfo%myHaloSize))
      allocate(haloSendGZsfc(intInfo%yourHaloSize))
      myHaloValues_r48%GZsfc(:) = 0.0d0
      haloSendGZsfc(:) = 0.0d0

      nsize   = 1

      ! First post the recv commands
      numRecv = 0
      do haloIndex = 1, intInfo%myHaloSize
        recvTag = intInfo%myHaloMpiTag(haloIndex)
        procSrc = intInfo%myHaloMpiIdSrc(haloIndex)
        numRecv = numRecv + 1
        call mpi_irecv(myHaloValues_r48%GZsfc(haloIndex),    &
                       nsize, mmpi_real8, procSrc, recvTag,  &
                       mmpi_comm_grid, requestIdRecv(numRecv))
      end do

      ! Prepare the data for sending to other tiles
      numSend = 0
      ptr2d_r8 => gsv_getHeightSfc(stateVector)
      do haloIndex = 1, intInfo%yourHaloSize
        lonIndex = intInfo%yourHaloLonIndex(haloIndex)
        latIndex = intInfo%yourHaloLatIndex(haloIndex)
        ! Handle periodicity in longitude
        if (lonIndex > stateVector%ni) lonIndex = 1
        haloSendGZsfc(haloIndex) = ptr2d_r8(lonIndex,latIndex)

        ! Post the send commands
        sendTag  = intInfo%yourHaloMpiTag(haloIndex)
        procDest = intInfo%yourHaloMpiIdDst(haloIndex)
        numSend = numSend + 1
        call mpi_isend(haloSendGZsfc(haloIndex),             &
                       nsize, mmpi_real8, procDest, sendTag, &
                       mmpi_comm_grid, requestIdSend(numSend))
      end do

      ! Wait for all previous RECV communications to finish before continuing
      if (numRecv > 0) then
        call mpi_waitAll(numRecv, requestIdRecv(1:numRecv), MPI_STATUSES_IGNORE, ierr)
      end if
      if (numSend > 0) then
        call mpi_waitAll(numSend, requestIdSend(1:numSend), MPI_STATUSES_IGNORE, ierr)
      end if

      deallocate(haloSendGZsfc)

    end if HeightSfcPresent

    nsize   = intInfo%numVarLevState * stateVector%numStep

    if (stateVector%dataKind == 4) then

      ! Allocate halo and other arrays
      myHaloValues_r48%dataKind = 4
      allocate(myHaloValues_r48%r4(intInfo%numVarLevState,stateVector%numStep,intInfo%myHaloSize))
      allocate(haloSend_r4(intInfo%numVarLevState,stateVector%numStep,intInfo%yourHaloSize))
      myHaloValues_r48%r4(:,:,:) = 0.0
      haloSend_r4(:,:,:) = 0.0

      ! First post the recv commands
      numRecv = 0
      do haloIndex = 1, intInfo%myHaloSize
        recvTag = intInfo%myHaloMpiTag(haloIndex)
        procSrc = intInfo%myHaloMpiIdSrc(haloIndex)
        numRecv = numRecv + 1
        call mpi_irecv(myHaloValues_r48%r4(:,:,haloIndex),   &
                       nsize, mmpi_real4, procSrc, recvTag,  &
                       mmpi_comm_grid, requestIdRecv(numRecv))
      end do

      ! Send the data
      numSend = 0
      call gsv_getField(stateVector, ptr4d_r4)
      do haloIndex = 1, intInfo%yourHaloSize
        ! Prepare the data for sending to other tiles
        lonIndex = intInfo%yourHaloLonIndex(haloIndex)
        latIndex = intInfo%yourHaloLatIndex(haloIndex)
        ! Handle periodicity in longitude
        if (lonIndex > stateVector%ni) lonIndex = 1
        haloSend_r4(:,:,haloIndex) = ptr4d_r4(lonIndex,latIndex,:,:)

        ! Post the send commands
        sendTag  = intInfo%yourHaloMpiTag(haloIndex)
        procDest = intInfo%yourHaloMpiIdDst(haloIndex)
        numSend = numSend + 1
        call mpi_isend(haloSend_r4(:,:,haloIndex),           &
                       nsize, mmpi_real4, procDest, sendTag, &
                       mmpi_comm_grid, requestIdSend(numSend))
      end do

    else ! For dataKind=8

      ! Allocate halo and other arrays
      myHaloValues_r48%dataKind = 8
      allocate(myHaloValues_r48%r8(intInfo%numVarLevState,stateVector%numStep,intInfo%myHaloSize))
      allocate(haloSend_r8(intInfo%numVarLevState,stateVector%numStep,intInfo%yourHaloSize))
      myHaloValues_r48%r8(:,:,:) = 0.0d0
      haloSend_r8(:,:,:) = 0.0d0

      ! First post the recv commands
      numRecv = 0
      do haloIndex = 1, intInfo%myHaloSize
        recvTag = intInfo%myHaloMpiTag(haloIndex)
        procSrc = intInfo%myHaloMpiIdSrc(haloIndex)
        numRecv = numRecv + 1
        call mpi_irecv(myHaloValues_r48%r8(:,:,haloIndex),   &
                       nsize, mmpi_real8, procSrc, recvTag,  &
                       mmpi_comm_grid, requestIdRecv(numRecv))
      end do

      ! Send the data
      numSend = 0
      call gsv_getField(stateVector, ptr4d_r8)
      do haloIndex = 1, intInfo%yourHaloSize
        ! Prepare the data for sending to other tiles
        lonIndex = intInfo%yourHaloLonIndex(haloIndex)
        latIndex = intInfo%yourHaloLatIndex(haloIndex)
        ! Handle periodicity in longitude
        if (lonIndex > stateVector%ni) lonIndex = 1
        haloSend_r8(:,:,haloIndex) = ptr4d_r8(lonIndex,latIndex,:,:)

        ! Post the send commands
        sendTag  = intInfo%yourHaloMpiTag(haloIndex)
        procDest = intInfo%yourHaloMpiIdDst(haloIndex)
        numSend = numSend + 1
        call mpi_isend(haloSend_r8(:,:,haloIndex),           &
                       nsize, mmpi_real8, procDest, sendTag, &
                       mmpi_comm_grid, requestIdSend(numSend))
      end do

    end if

    ! Wait for all previous RECV communications to finish before continuing
    if (numRecv > 0) then
      call mpi_waitAll(numRecv, requestIdRecv(1:numRecv), MPI_STATUSES_IGNORE, ierr)
    end if
    if (numSend > 0) then
      call mpi_waitAll(numSend, requestIdSend(1:numSend), MPI_STATUSES_IGNORE, ierr)
    end if

    if (.not. beSilent) write(*,*) 'sendRecvHalo (s2c): Finished'

  end subroutine sendRecvHalo

  !---------------------------------------------------------
  ! sendRecvHaloAd
  !---------------------------------------------------------
  subroutine sendRecvHaloAd(stateVector, intInfo, beSilent, &
                            myHaloValuesAd_r48)
    !
    ! :Purpose: Adjoint of exchange data between mpi tasks to provide a list of the
    !           data values within the halo. That is, the values in the halo list
    !           are communicated back to the tiles where they originate and added to
    !           the interior tile values in stateVector.
    !
    implicit none

    ! Arguments:
    type(struct_gsv)       ,      intent(inout) :: stateVector         ! stateVector object
    type(struct_interpInfoTiles), intent(in)    :: intInfo             ! Interpolation info structure
    logical,                      intent(in)    :: beSilent            ! Choose to print nothing
    type(real48_3d),              intent(in)    :: myHaloValuesAd_r48  ! Halo values

    ! Locals:
    integer :: haloIndex, numRecv, numSend, recvTag, sendTag, nsize, ierr
    integer :: lonIndex, latIndex
    integer :: procSrc, procDest
    type(mpi_request), allocatable :: requestIdRecv(:), requestIdSend(:)
    real(4), allocatable :: haloRecv_r4(:,:,:)
    real(8), allocatable :: haloRecv_r8(:,:,:)
    real(8), allocatable :: haloRecvGZsfc(:)
    real(4), pointer     :: ptr4d_r4(:,:,:,:)
    real(8), pointer     :: ptr4d_r8(:,:,:,:)
    real(8), pointer     :: ptr2d_r8(:,:)

    if (.not. beSilent) write(*,*) 'sendRecvHaloAd (s2c): Starting'

    allocate(requestIdSend(intInfo%myHaloSize))
    allocate(requestIdRecv(intInfo%yourHaloSize))

    if (.not. beSilent) then
      write(*,*) 'sendRecvHaloAd (s2c): myHaloSize, yourHaloSize = ', intInfo%myHaloSize, intInfo%yourHaloSize
    end if

    heightSfcPresent: if (stateVector%HeightSfcPresent) then

      ! Allocate halo and other arrays
      allocate(haloRecvGZsfc(intInfo%yourHaloSize))

      nsize   = 1

      ! First post the recv commands
      numRecv = 0
      do haloIndex = 1, intInfo%yourHaloSize
        recvTag = intInfo%yourHaloMpiTag(haloIndex)
        procSrc = intInfo%yourHaloMpiIdDst(haloIndex)
        numRecv = numRecv + 1
        call mpi_irecv(haloRecvGZsfc(haloIndex),             &
                       nsize, mmpi_real8, procSrc, recvTag,  &
                       mmpi_comm_grid, requestIdRecv(numSend))
      end do

      ! Post the send commands
      numSend = 0
      do haloIndex = 1, intInfo%myHaloSize
        sendTag  = intInfo%myHaloMpiTag(haloIndex)
        procDest = intInfo%myHaloMpiIdSrc(haloIndex)
        numRecv = numRecv + 1
        call mpi_isend(myHaloValuesAd_r48%GZsfc(haloIndex),  &
                       nsize, mmpi_real8, procDest, sendTag, &
                       mmpi_comm_grid, requestIdSend(numRecv))
      end do

      ! Wait for all previous RECV communications to finish before continuing
      if (numRecv > 0) then
        call mpi_waitAll(numRecv, requestIdRecv(1:numRecv), MPI_STATUSES_IGNORE, ierr)
      end if
      if (numSend > 0) then
        call mpi_waitAll(numSend, requestIdSend(1:numSend), MPI_STATUSES_IGNORE, ierr)
      end if

      ! Copy recv'ed data into stateVector object
      ptr2d_r8 => gsv_getHeightSfc(stateVector)
      do haloIndex = 1, intInfo%yourHaloSize
        lonIndex = intInfo%yourHaloLonIndex(haloIndex)
        latIndex = intInfo%yourHaloLatIndex(haloIndex)
        ! Handle periodicity in longitude
        if (lonIndex > stateVector%ni) lonIndex = 1
        ptr2d_r8(lonIndex,latIndex) = ptr2d_r8(lonIndex,latIndex) +  &
                                      haloRecvGZsfc(haloIndex)
      end do

      deallocate(haloRecvGZsfc)

    end if HeightSfcPresent

    nsize   = intInfo%numVarLevState * stateVector%numStep

    if (stateVector%dataKind == 4) then

      ! Allocate halo and other arrays
      allocate(haloRecv_r4(intInfo%numVarLevState,stateVector%numStep,intInfo%yourHaloSize))
      haloRecv_r4(:,:,:) = 0.0

      ! First post the recv commands
      numRecv = 0
      do haloIndex = 1, intInfo%yourHaloSize
        recvTag = intInfo%yourHaloMpiTag(haloIndex)
        procSrc = intInfo%yourHaloMpiIdDst(haloIndex)
        numRecv = numRecv + 1
        call mpi_irecv(haloRecv_r4(:,:,haloIndex),           &
                       nsize, mmpi_real4, procSrc, recvTag,  &
                       mmpi_comm_grid, requestIdRecv(numRecv))
      end do

      ! Send the data
      numSend = 0
      do haloIndex = 1, intInfo%myHaloSize
        sendTag = intInfo%myHaloMpiTag(haloIndex)
        procDest = intInfo%myHaloMpiIdSrc(haloIndex)
        numSend = numSend + 1
        call mpi_isend(myHaloValuesAd_r48%r4(:,:,haloIndex), &
                       nsize, mmpi_real4, procDest, sendTag, &
                       mmpi_comm_grid, requestIdSend(numSend))
      end do

    else ! For dataKind=8

      ! Allocate halo and other arrays
      allocate(haloRecv_r8(intInfo%numVarLevState,stateVector%numStep,intInfo%yourHaloSize))
      haloRecv_r8(:,:,:) = 0.0d0

      ! First post the recv commands
      numRecv = 0
      do haloIndex = 1, intInfo%yourHaloSize
        recvTag = intInfo%yourHaloMpiTag(haloIndex)
        procSrc = intInfo%yourHaloMpiIdDst(haloIndex)
        numRecv = numRecv + 1
        call mpi_irecv(haloRecv_r8(:,:,haloIndex),           &
                       nsize, mmpi_real8, procSrc, recvTag,  &
                       mmpi_comm_grid, requestIdRecv(numRecv))
      end do

      ! Send the data
      numSend = 0
      do haloIndex = 1, intInfo%myHaloSize
        sendTag = intInfo%myHaloMpiTag(haloIndex)
        procDest = intInfo%myHaloMpiIdSrc(haloIndex)
        numSend = numSend + 1
        call mpi_isend(myHaloValuesAd_r48%r8(:,:,haloIndex), &
                       nsize, mmpi_real8, procDest, sendTag, &
                       mmpi_comm_grid, requestIdSend(numSend))
      end do

    end if

    ! Wait for all previous RECV communications to finish before continuing
    if (numRecv > 0) then
      call mpi_waitAll(numRecv, requestIdRecv(1:numRecv), MPI_STATUSES_IGNORE, ierr)
    end if
    if (numSend > 0) then
      call mpi_waitAll(numSend, requestIdSend(1:numSend), MPI_STATUSES_IGNORE, ierr)
    end if

    ! Copy recv'ed data into stateVector object
    if (stateVector%dataKind == 4) then

      call gsv_getField(stateVector, ptr4d_r4)
      do haloIndex = 1, intInfo%yourHaloSize
        lonIndex = intInfo%yourHaloLonIndex(haloIndex)
        latIndex = intInfo%yourHaloLatIndex(haloIndex)
        ! Handle periodicity in longitude
        if (lonIndex > stateVector%ni) lonIndex = 1
        ptr4d_r4(lonIndex,latIndex,:,:) = ptr4d_r4(lonIndex,latIndex,:,:) +  &
                                          haloRecv_r4(:,:,haloIndex)
      end do

    else

      call gsv_getField(stateVector, ptr4d_r8)
      do haloIndex = 1, intInfo%yourHaloSize
        lonIndex = intInfo%yourHaloLonIndex(haloIndex)
        latIndex = intInfo%yourHaloLatIndex(haloIndex)
        ! Handle periodicity in longitude
        if (lonIndex > stateVector%ni) lonIndex = 1
        ptr4d_r8(lonIndex,latIndex,:,:) = ptr4d_r8(lonIndex,latIndex,:,:) +  &
                                          haloRecv_r8(:,:,haloIndex)
      end do

    end if

    if (.not. beSilent) write(*,*) 'sendRecvHaloAd (s2c): Finished'

  end subroutine sendRecvHaloAd

  !-------------------------------------------------------------
  ! End of subroutines for the newer "TILES" mpiMode
  !-------------------------------------------------------------

  !---------------------------------------------------------
  ! pressureProfileMonotonicityCheck
  !---------------------------------------------------------
  subroutine pressureProfileMonotonicityCheck(obsSpaceData, column)
    !
    ! :Purpose: Check for non monotonic pressure profiles that can be computed in slantpathmode
    !
    implicit none

    ! Arguments:
    type(struct_obs),        intent(inout) :: obsSpaceData
    type(struct_columnData), intent(inout) :: column

    ! Locals:
    integer, parameter :: numWriteMax = 10
    integer :: headerIndex, bodyIndex, iterationCount, singularIndex, levelIndex
    integer :: pressureVarIndex
    integer :: nlv
    integer :: numWrites
    real(8), pointer :: pressureProfile(:)
    logical :: monotonicProfile
    integer, parameter :: nPressureVar =2
    character(len=4), parameter :: pressureVarList(nPressureVar)=['P_T ', 'P_M ']

    write(*,*) ' '
    write(*,*) 'pressureProfileMonotonicityCheck: START'
    write(*,*) ' '

    numWrites = 0

    call obs_set_current_header_list(obsSpaceData, 'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER
      do pressureVarIndex = 1, nPressureVar
        pressureProfile => col_getColumn(column, headerIndex, pressureVarList(pressureVarIndex))
        nlv = size(pressureProfile)
        monotonicProfile = .true.
        iterationCount = 0
        iterationLoop:do
          singularIndex = -1
          levelSearch:do levelIndex = 1, nlv - 1
            if ( pressureProfile(levelIndex) > pressureProfile(levelIndex+1)) then
              singularIndex = levelIndex
              exit levelSearch
            end if
          end do levelSearch
          if ( singularIndex == -1 ) exit iterationLoop !regular profile or correction OK
          iterationCount = iterationCount + 1
          if (iterationCount == 1) then
            monotonicProfile = .false.
            if (numWrites < numWriteMax) then
              numWrites = numWrites + 1
              write(*,*) 'pressureProfileMonotonicityCheck: found non monotonic pressure profile:', &
                   pressureVarList(pressureVarIndex), pressureProfile
            end if
          end if
          if (singularIndex == 1) then !should never happen
            write(*,*) 'pressureProfileMonotonicityCheck: ', pressureProfile(1:2)
            call utl_abort('pressureProfileMonotonicityCheck: profile in the wrong order ?' &
                 // pressureVarList(pressureVarIndex))
          end if
          pressureProfile(singularIndex) = 0.5d0 * ( pressureProfile(singularIndex - 1) + pressureProfile(singularIndex + 1) )
          write(*,*) 'pressureProfileMonotonicityCheck: profile iteration', &
               pressureVarList(pressureVarIndex), iterationCount
        end do iterationLoop

        ! if requested reject the corrected profile
        if (.not. monotonicProfile .and. rejectObsNonMonotonicPressure) then
          call obs_headSet_i(obsSpaceData,OBS_ST1,headerIndex, &
               ibset(obs_headElem_i(obsSpaceData,OBS_ST1,headerIndex), 05))
          call obs_set_current_body_list(obsSpaceData, headerIndex)
          BODY: do
            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if (bodyIndex < 0) exit BODY
            if (rejectObsNonMonotonicPressure) then
              call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
              call flg_setFlag(obsSpaceData, bodyIndex, flg_09rejBgck)
            end if
          end do BODY
        end if
      end do ! loop on pressure variables

    end do HEADER

    write(*,*) 'pressureProfileMonotonicityCheck: END'

  end subroutine pressureProfileMonotonicityCheck

  !---------------------------------------------------------
  ! latlonChecksAnlGrid
  !---------------------------------------------------------
  subroutine latlonChecksAnlGrid(obsSpaceData, hco_core, moveObsAtPole)
    !
    ! :Purpose: Check the lat/lon of observations and modify if necessary
    !
    implicit none

    ! Arguments:
    type(struct_obs)         , intent(inout) :: obsSpaceData
    type(struct_hco), pointer, intent(in)    :: hco_core
    logical                  , intent(in)    :: moveObsAtPole

    ! Locals:
    integer :: headerIndex, ierr
    integer :: idata, idatend, jdata, subGridIndex
    real(4) :: lat_r4, lon_r4, lat_deg_r4, lon_deg_r4
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    real(4) :: xposLowerBoundAnl_r4, xposUpperBoundAnl_r4
    real(8) :: lat_r8, lon_r8
    integer, save :: numWrites = 0
    integer :: gdllfxy

    write(*,*) ' '
    write(*,*) 'latlonChecksAnlGrid: STARTING'
    write(*,*) ' '
    call msg_memUsage('latlonChecksAnlGrid')

    !
    !-    Get the Analysis Grid structure
    !
    if ( hco_core % global ) then
       xposLowerBoundAnl_r4 = - huge(1.0) ! no limit since grid is global (periodic)
       xposUpperBoundAnl_r4 = + huge(1.0) ! no limit since grid is global (periodic)
    else
       xposLowerBoundAnl_r4 = 1.0
       xposUpperBoundAnl_r4 = real(hco_core % ni)
    end if

    header_loop: do headerIndex=1, obs_numheader(obsSpaceData)

      !- Get LatLon of observation location
      lat_r8 = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)
      lon_r8 = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex)
      lat_r4 = real(lat_r8,4)
      lon_r4 = real(lon_r8,4)
      if (lon_r4.lt.0.0         ) lon_r4 = lon_r4 + 2.0*MPC_PI_R4
      if (lon_r4.ge.2.*MPC_PI_R4) lon_r4 = lon_r4 - 2.0*MPC_PI_R4

      lat_deg_r4 = lat_r4 * MPC_DEGREES_PER_RADIAN_R4 ! Radian To Degree
      lon_deg_r4 = lon_r4 * MPC_DEGREES_PER_RADIAN_R4

      !
      !- Find the position in the analysis grid
      !
      ierr = gpos_getPositionXY(hco_core % EZscintID,  &
                                xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                                lat_deg_r4, lon_deg_r4, subGridIndex)

      !- Test if the obs is outside the analysis grid
      if ( xpos_r4 < xposLowerBoundAnl_r4  .or. &
           xpos_r4 > xposUpperBoundAnl_r4  .or. &
           ypos_r4 < 1.0                   .or. &
           ypos_r4 > real(hco_core % nj) ) then

        if ( hco_core % global ) then

          if ( moveObsAtPole ) then
            ! Modify latitude if we have an observation at or near the poles
            write(*,*) ''
            write(*,*) 'latlonChecksAnlGrid: Moving OBS inside the GLOBAL ANALYSIS grid, ', headerIndex
            write(*,*) '  true position : ', lat_deg_r4, lon_deg_r4, ypos_r4, xpos_r4

            !- Move the observation to the nearest grid point
            if ( ypos_r4 < 1.0 )                ypos_r4 = 1.0
            if ( ypos_r4 > real(hco_core % nj) ) ypos_r4 = real(hco_core % nj)

            ierr = gdllfxy( hco_core % EZscintID, &    ! IN
                            lat_deg_r4, lon_deg_r4, & ! OUT
                            xpos_r4, ypos_r4, 1)      ! IN

            write(*,*) '  new  position : ', lat_deg_r4, lon_deg_r4, ypos_r4, xpos_r4

            lat_r8 = real(lat_deg_r4,8) * MPC_RADIANS_PER_DEGREE_R8
            lon_r8 = real(lon_deg_r4,8) * MPC_RADIANS_PER_DEGREE_R8
            call obs_headSet_r(obsSpaceData,OBS_LAT,headerIndex, lat_r8) ! IN
            call obs_headSet_r(obsSpaceData,OBS_LON,headerIndex, lon_r8) ! IN
          else
            write(*,*)
            write(*,*) 'latlonChecksAnlGrid: OBS outside the GLOBAL ANALYSIS grid, but NOT moved, ', headerIndex
          end if

        else
          ! The observation is outside the domain
          ! In LAM Analysis mode we must discard this observation
          numWrites = numWrites + 1
          if (numWrites < maxNumWrites) then
            write(*,*) 'latlonChecksAnlGrid: Rejecting OBS outside the LAM ANALYSIS grid domain, ', headerIndex
            write(*,*) '  position : ', lat_deg_r4, lon_deg_r4, ypos_r4, xpos_r4
          else if (numWrites == maxNumWrites) then
            write(*,*) 'latlonChecksAnlGrid: More rejects, but reached maximum number of writes to the listing.'
          end if

          idata   = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
          idatend = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + idata -1
          do jdata = idata, idatend
            call obs_bodySet_i(obsSpaceData,OBS_ASS,JDATA, obs_notAssimilated)
          end do
          call obs_headSet_i(obsSpaceData,OBS_ST1,headerIndex,  &
                             ibset( obs_headElem_i(obsSpaceData,OBS_ST1,headerIndex), 05))
        end if

      end if

    end do header_loop

    write(*,*) 'latlonChecksAnlGrid: END'

  end subroutine latlonChecksAnlGrid

  !---------------------------------------------------------
  ! setupInterpInfo2dFields
  !---------------------------------------------------------
  subroutine setupInterpInfo2dFields(interpInfo, obsSpaceData, stateVector,  &
                                     headerIndexBeg, headerIndexEnd, &
                                     timeInterpType, rejectOutsideObs, &
                                     inputStateVectorType, lastCall_opt)
    !
    ! :Purpose: Setup all of the information needed to quickly
    !           perform the horizontal interpolation to the observation
    !           locations.
    !
    implicit none

    ! Arguments
    type(struct_interpInfo),   intent(out)   :: interpInfo
    type(struct_obs)        ,  intent(inout) :: obsSpaceData
    type(struct_gsv), target,  intent(in)    :: stateVector
    integer                 ,  intent(in)    :: headerIndexBeg
    integer                 ,  intent(in)    :: headerIndexEnd
    logical                 ,  intent(inout) :: rejectOutsideObs
    character(len=*)        ,  intent(in)    :: timeInterpType
    character(len=*)        ,  intent(in)    :: inputStateVectorType
    logical,       optional ,  intent(in)    :: lastCall_opt

    ! Locals:
    type(struct_gsv)          :: stateVector_VarsLevs_1Step, stateVector_Tiles_allVar_1Step
    type(struct_gsv)          :: stateVector_Tiles_1Step
    type(struct_gsv), save    :: stateVector_1Step
    type(struct_gsv), pointer :: stateVector_Tiles_ptr
    integer :: numHeader, numHeaderUsedMax, headerIndex, headerUsedIndex
    integer :: varLevIndex, varLevIndexCount, myVarLevBeg
    integer :: numStep, stepIndex, ierr
    integer :: procIndex, niP1, numGridptTotal, numGridptTotalMpi, numHeaderUsed
    integer :: subGridIndex, subGridForInterp, numSubGridsForInterp
    real(8) :: latRot, lonRot, lat, lon
    real(4) :: lon_r4, lat_r4, lon_deg_r4, lat_deg_r4
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    real(4) :: footprintRadius_r4 ! (metres)
    integer, allocatable :: numGridpt(:), allNumHeaderUsed(:,:)
    integer, allocatable :: allHeaderIndex(:,:,:), headerIndexVec(:,:)
    real(8), allocatable :: lat_send_r8(:,:), lat_recv_r8(:,:), lon_send_r8(:,:), lon_recv_r8(:,:)
    real(4), allocatable :: footprintRadiusVec_r4(:), allFootprintRadius_r4(:,:,:)
    real(8), allocatable :: allLatOneLev(:,:)
    real(8), allocatable :: allLonOneLev(:,:)
    character(len=4), pointer :: varNames(:)
    character(len=4)          :: varLevel
    real(8), allocatable :: latColumn(:,:), lonColumn(:,:)
    real(8), allocatable :: latLev_T(:), lonLev_T(:), latLev_M(:), lonLev_M(:)
    real(8) :: latLev_S, lonLev_S
    real(4), pointer :: height3D_r4_ptr1(:,:,:), height3D_r4_ptr2(:,:,:)
    real(4), save, pointer :: height3D_T_r4(:,:,:), height3D_M_r4(:,:,:)
    real(8), pointer :: height3D_r8_ptr1(:,:,:)
    real(kdkind), allocatable :: positionArray(:,:)
    integer :: sendsizes(mmpi_nprocs), recvsizes(mmpi_nprocs), senddispls(mmpi_nprocs)
    integer :: recvdispls(mmpi_nprocs), allVarLevBeg(mmpi_nprocs)
    integer :: codeType, nlev_T, nlev_M, levIndex
    integer :: lonIndex, latIndex, gridIndex
    integer :: maxkcount, numVarLevToSend, numTovsUsingFootprint, numAllTovs
    logical :: doSlantPath, SlantTO, SlantRO, SlantRA
    logical :: firstHeaderSlantPathTO, firstHeaderSlantPathRO, firstHeaderSlantPathRA
    logical :: doSetup3dHeights, lastCall
    type(kdtree2), pointer  :: tree

    write(*,*) 'setupInterpInfo2dFields: STARTING'
    call msg_memUsage('setupInterpInfo2dFields')

    if (present(lastCall_opt)) then
      lastCall = lastCall_opt
    else
      lastCall = .false.
    end if

    if (inputStateVectorType == 'nl' .and. rejectObsOutsideGlobalGrid) then
      rejectOutsideObs = .true.
    end if

    doSlantPath = .false.
    SlantTO     = .false.
    SlantRO     = .false.
    SlantRA     = .false.
    if (slantPath_TO_nl   .and. inputStateVectorType == 'nl') then
      doSlantPath = .true.
      SlantTO     = .true.
    endif
    if (slantPath_TO_tlad .and. inputStateVectorType /= 'nl') then
      doSlantPath = .true.
      SlantTO     = .true.
    endif
    if (slantPath_RO_nl   .and. inputStateVectorType == 'nl') then
      doSlantPath = .true.
      SlantRO     = .true.
    endif
    if (slantPath_RA_nl   .and. inputStateVectorType == 'nl') then
      doSlantPath = .true.
      SlantRA     = .true.
    endif
    write(*,*) 'setupInterpInfo2dFields: doSlantPath, SlantTO, SlantRO, SlantRA = ', &
               doSlantPath, SlantTO, SlantRO, SlantRA

    numStep = stateVector%numStep
    numHeader = headerIndexEnd - headerIndexBeg + 1

    call oti_setup(interpInfo%oti, obsSpaceData, numStep,  &
                   headerIndexBeg, headerIndexEnd, &
                   interpType_opt=timeInterpType, flagObsOutside_opt=.true.)

    if ((stateVector%heightSfcPresent) .and. (mmpi_myid == 0)) then
      myVarLevBeg = 0
    else
      myVarLevBeg = stateVector%myVarLevBeg
    end if

    call mmpi_allGather(myVarLevBeg, allVarLevBeg)

    ! Allow for periodicity in Longitude for global Gaussian grid
    if (stateVector%hco%grtyp == 'G' .or. &
        (stateVector%hco%grtyp == 'Z' .and. stateVector%hco%global)) then
      niP1 = stateVector%ni + 1
    else
      niP1 = stateVector%ni
    end if

    ! First count the number of headers for each stepIndex
    allocate(allNumHeaderUsed(numStep,mmpi_nprocs))
    allNumHeaderUsed(:,:) = 0
    do stepIndex = 1, numStep
      numHeaderUsed = 0

      header_loop1: do headerIndex = headerIndexBeg, headerIndexEnd

        ! if obs inside window, but zero weight for current stepIndex then skip it
        if ( utl_isEqual(oti_getTimeInterpWeight(interpInfo%oti,headerIndex,stepIndex), 0.0d0) ) then
          cycle header_loop1
        end if

        numHeaderUsed = numHeaderUsed + 1

      end do header_loop1

      ! gather the number of obs over all processors for each timestep
      call mmpi_allGather(numHeaderUsed, allNumHeaderUsed(stepIndex,:))
    end do

    numHeaderUsedMax = maxval(allNumHeaderUsed(:,:))
    write(*,*) 'setupInterpInfo2dFields: numHeaderUsedMax = ', numHeaderUsedMax

    ! temporary arrays
    allocate(headerIndexVec(numHeaderUsedMax,numStep))
    allocate(footprintRadiusVec_r4(numHeaderUsedMax))
    headerIndexVec(:,:) = 0

    ! copy the horizontal grid object
    interpInfo%hco => stateVector%hco

    ! setup the information for wind rotation
    if (gsv_varExist(varName='UU') .and. &
        gsv_varExist(varName='VV') .and. &
        interpInfo%hco%rotated) then
      call uvr_Setup(interpInfo%uvr, & ! INOUT
                     stateVector%hco)  ! IN
    end if

    allocate(interpInfo%stepProcData(mmpi_nprocs,numStep))
    do stepIndex = 1,numStep
      do procIndex = 1, mmpi_nprocs
        allocate(interpInfo%stepProcData(procIndex,stepIndex)%allLat(allNumHeaderUsed(stepIndex,procIndex),myVarLevBeg:stateVector%myVarLevEnd))
        allocate(interpInfo%stepProcData(procIndex,stepIndex)%allLon(allNumHeaderUsed(stepIndex,procIndex),myVarLevBeg:stateVector%myVarLevEnd))
        interpInfo%stepProcData(procIndex,stepIndex)%allLat(:,:) = 0.0d0
        interpInfo%stepProcData(procIndex,stepIndex)%allLon(:,:) = 0.0d0

        allocate(interpInfo%stepProcData(procIndex,stepIndex)%allHeaderIndex(allNumHeaderUsed(stepIndex,procIndex)))
        interpInfo%stepProcData(procIndex,stepIndex)%allHeaderIndex(:) = 0

        allocate(interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(interpInfo%hco%numSubGrid,numHeaderUsedMax,myVarLevBeg:stateVector%myVarLevEnd))
        allocate(interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(interpInfo%hco%numSubGrid,numHeaderUsedMax,myVarLevBeg:stateVector%myVarLevEnd))
        interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(:,:,:) = 0
        interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(:,:,:) = -1
      end do
    end do

    ! allocate arrays that will be returned
    allocate(interpInfo%allNumHeaderUsed(numStep,mmpi_nprocs))
    allocate(allLatOneLev(numHeaderUsedMax,mmpi_nprocs))
    allocate(allLonOneLev(numHeaderUsedMax,mmpi_nprocs))
    allocate(allFootprintRadius_r4(numHeaderUsedMax,numStep,mmpi_nprocs))
    allocate(numGridpt(interpInfo%hco%numSubGrid))
    allFootprintRadius_r4(:,:,:) = bilinearFootprint
    interpInfo%allNumHeaderUsed(:,:) = allNumHeaderUsed(:,:)

    if (interpInfo%hco%rotated) then
      do stepIndex = 1, numStep
        do procIndex = 1, mmpi_nprocs
          allocate(interpInfo%stepProcData(procIndex,stepIndex)%allLatRot(interpInfo%hco%numSubGrid,allNumHeaderUsed(stepIndex,procIndex),myVarLevBeg:stateVector%myVarLevEnd))
          allocate(interpInfo%stepProcData(procIndex,stepIndex)%allLonRot(interpInfo%hco%numSubGrid,allNumHeaderUsed(stepIndex,procIndex),myVarLevBeg:stateVector%myVarLevEnd))
          interpInfo%stepProcData(procIndex,stepIndex)%allLatRot(:,:,:) = 0.0d0
          interpInfo%stepProcData(procIndex,stepIndex)%allLonRot(:,:,:) = 0.0d0
        end do
      end do
    end if

    nlev_T = gsv_getNumLev(stateVector,'TH')
    nlev_M = gsv_getNumLev(stateVector,'MM')

    doSetup3dHeights = doSlantPath .and.  &
                       .not. gsv_isAllocated(stateVector_1Step) .and. &
                       gsv_varExist(stateVector,'Z_T') .and. &
                       gsv_varExist(stateVector,'Z_M')

    ! prepare for extracting the 3D height for slant-path calculation
    if ( doSetup3dHeights ) then

      write(*,*) 'setupInterpInfo2dFields: extracting 3D heights for slant-path for ', inputStateVectorType

      if ( inputStateVectorType == 'nl' ) then
        nullify(varNames)
        call gsv_varNamesList(varNames, stateVector)
        call gsv_allocate(stateVector_VarsLevs_1Step, 1, &
                          stateVector%hco, stateVector%vco, &
                          mpi_local_opt=.true., mpi_distribution_opt='VarsLevs', &
                          dataKind_opt=4, varNames_opt=varNames)

        call gsv_getField(stateVector,height3D_r4_ptr1)
        call gsv_getField(stateVector_VarsLevs_1Step,height3D_r4_ptr2)
        height3D_r4_ptr2(:,:,:) = height3D_r4_ptr1(:,:,:)

        call gsv_allocate(stateVector_Tiles_allVar_1Step, 1, &
                          stateVector%hco, stateVector%vco, &
                          mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                          dataKind_opt=4, varNames_opt=varNames)

        call gsv_transposeVarsLevsToTiles( stateVector_VarsLevs_1Step, stateVector_Tiles_allVar_1Step )
        call gsv_deallocate(statevector_VarsLevs_1Step)

        call gsv_allocate( stateVector_Tiles_1Step, 1, &
                           stateVector%hco, stateVector%vco, &
                           mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                           dataKind_opt=4, varNames_opt=(/'Z_M','Z_T'/) )

        call gsv_getField(stateVector_Tiles_allVar_1Step,height3D_r4_ptr1,'Z_T')
        call gsv_getField(stateVector_Tiles_1Step,height3D_r4_ptr2,'Z_T')
        height3D_r4_ptr2(:,:,:) = height3D_r4_ptr1(:,:,:)

        call gsv_getField(stateVector_Tiles_allVar_1Step,height3D_r4_ptr1,'Z_M')
        call gsv_getField(stateVector_Tiles_1Step,height3D_r4_ptr2,'Z_M')
        height3D_r4_ptr2(:,:,:) = height3D_r4_ptr1(:,:,:)

        call gsv_deallocate(stateVector_Tiles_allVar_1Step)

      else
        stateVector_Tiles_ptr => gvt_getStateVectorTrial('height')

        call gsv_allocate(stateVector_Tiles_1Step, 1, &
                          stateVector%hco, stateVector%vco, &
                          mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                          dataKind_opt=4, varNames_opt=(/'Z_M','Z_T'/))

        call gsv_getField(stateVector_Tiles_ptr,height3D_r8_ptr1,'Z_T')
        call gsv_getField(stateVector_Tiles_1Step,height3D_r4_ptr2,'Z_T')
        height3D_r4_ptr2(:,:,:) = real(height3D_r8_ptr1(:,:,:),4)

        call gsv_getField(stateVector_Tiles_ptr,height3D_r8_ptr1,'Z_M')
        call gsv_getField(stateVector_Tiles_1Step,height3D_r4_ptr2,'Z_M')
        height3D_r4_ptr2(:,:,:) = real(height3D_r8_ptr1(:,:,:),4)

      end if ! inputStateVectorType

      ! Communicate 3D height fields onto all mpi tasks
      call gsv_allocate(stateVector_1Step, 1, &
                        stateVector%hco, stateVector%vco, &
                        mpi_local_opt=.false., &
                        dataKind_opt=4, varNames_opt=(/'Z_M','Z_T'/))
      call utl_tmg_start(32,'------s2c_Slant')
      call gsv_transposeTilesToMpiGlobal(stateVector_1Step, stateVector_Tiles_1Step)
      call utl_tmg_stop(32)
      call gsv_getField(stateVector_1Step,height3D_T_r4,'Z_T')
      call gsv_getField(stateVector_1Step,height3D_M_r4,'Z_M')

      write(*,*) 'setupInterpInfo2dFields, height3D_T_r4='
      write(*,*) height3D_T_r4(1,1,:)
      write(*,*) 'setupInterpInfo2dFields, height3D_M_r4='
      write(*,*) height3D_M_r4(1,1,:)

      call msg_memUsage('setupInterpInfo2dFields')
    end if ! doSetup3dHeights

    ! get observation lat-lon and footprint radius onto all mpi tasks
    step_loop2: do stepIndex = 1, numStep
      numHeaderUsed = 0

      footprintRadiusVec_r4(:) = bilinearFootprint

      header_loop2: do headerIndex = headerIndexBeg, headerIndexEnd

        ! if obs inside window, but zero weight for current stepIndex then skip it
        if ( utl_isEqual(oti_getTimeInterpWeight(interpInfo%oti, headerIndex, stepIndex), 0.0d0) ) then
          cycle header_loop2
        end if

        numHeaderUsed = numHeaderUsed + 1
        headerIndexVec(numHeaderUsed,stepIndex) = headerIndex

        footprintRadiusVec_r4(numHeaderUsed) = s2c_getFootprintRadius(obsSpaceData, &
                                                                stateVector, headerIndex)

      end do header_loop2

      call mmpi_allGather(footprintRadiusVec_r4, allFootprintRadius_r4(:,stepIndex,:))

      allocate(latColumn(numHeaderUsedMax,allVarLevBeg(1):stateVector%numVarLev))
      allocate(lonColumn(numHeaderUsedMax,allVarLevBeg(1):stateVector%numVarLev))
      latColumn(:,:) = 0.0d0
      lonColumn(:,:) = 0.0d0

      if (doSlantPath .and. &
          gsv_varExist(stateVector,'Z_T') .and. &
          gsv_varExist(stateVector,'Z_M')) then

        allocate(latLev_T(nlev_T))
        allocate(lonLev_T(nlev_T))
        allocate(latLev_M(nlev_M))
        allocate(lonLev_M(nlev_M))
        latLev_T(:) = 0.0d0
        lonLev_T(:) = 0.0d0
        latLev_M(:) = 0.0d0
        lonLev_M(:) = 0.0d0

        firstHeaderSlantPathTO = .true.
        firstHeaderSlantPathRO = .true.
        firstHeaderSlantPathRA = .true.
        header_loop3: do headerUsedIndex = 1, numHeaderUsed
          headerIndex = headerIndexVec(headerUsedIndex,stepIndex)

          !- Get LatLon of observation location
          lat_r4 = real(obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex), 4)
          lon_r4 = real(obs_headElem_r(obsSpaceData, OBS_LON, headerIndex), 4)
          if (lon_r4 <  0.0            ) lon_r4 = lon_r4 + 2.0 * MPC_PI_R4
          if (lon_r4 >= 2.0 * MPC_PI_R4) lon_r4 = lon_r4 - 2.0 * MPC_PI_R4

          codeType = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)

          if (tvs_isIdBurpTovs(codeType) .and. SlantTO) then
            if (firstHeaderSlantPathTO) then
              write(*,'(a,i3,a,i8)') 'setupInterpInfo2dFields: start slant-path for TOVS. stepIndex = ', &
                   stepIndex,' and numHeaderUsed = ',numHeaderUsed
              firstHeaderSlantPathTO = .false.
            end if

            ! calculate lat/lon along the line of sight
            call utl_tmg_start(32,'------s2c_Slant')
            call slp_calcLatLonTovs(obsSpaceData, stateVector%hco, headerIndex, & ! IN
                                    height3D_T_r4, height3D_M_r4,               & ! IN
                                    latLev_T, lonLev_T,                         & ! OUT
                                    latLev_M, lonLev_M,                         & ! OUT
                                    latLev_S, lonLev_S             )              ! OUT
            call utl_tmg_stop(32)

          else if (codeType == codtyp_get_codtyp('ro') .and. SlantRO ) then
            if (firstHeaderSlantPathRO) then
              write(*,'(a,i3,a,i8)') 'setupInterpInfo2dFields: start slant-path for RO. stepIndex = ', &
                   stepIndex,' and numHeaderUsed = ',numHeaderUsed
              firstHeaderSlantPathRO = .false.
            end if

            ! Calculate lat/lon along the GPSRO obs
            call utl_tmg_start(32,'------s2c_Slant')
            call slp_calcLatLonRO(obsSpaceData, stateVector%hco, headerIndex, & ! IN
                                  height3D_T_r4, height3D_M_r4,               & ! IN
                                  latLev_T, lonLev_T,                         & ! OUT
                                  latLev_M, lonLev_M,                         & ! OUT
                                  latLev_S, lonLev_S                          ) ! OUT
            call utl_tmg_stop(32)
          else if (codeType == codtyp_get_codtyp('radar') .and. SlantRA ) then
            if ( firstHeaderSlantPathRA ) then
              write(*,'(a,i3,a,i8)') 'setupInterpInfo2dFields: start slant-path for RADAR. stepIndex=', &
                   stepIndex,' and numHeaderUsed=',numHeaderUsed
              firstHeaderSlantPathRA = .false.
            end if

             ! calculate lat/lon along the radar beam obs
             call slp_calcLatLonRadar(obsSpaceData, stateVector%hco, headerIndex, & ! IN
                                      height3D_T_r4, height3D_M_r4,                 & ! IN
                                      latLev_T, lonLev_T,                           & ! OUT
                                      latLev_M, lonLev_M,                           & ! OUT
                                      latLev_S, lonLev_S                           ) ! OUT
          else

            latLev_T(:) = real(lat_r4,8)
            lonLev_T(:) = real(lon_r4,8)
            latLev_M(:) = real(lat_r4,8)
            lonLev_M(:) = real(lon_r4,8)
            latLev_S = real(lat_r4,8)
            lonLev_S = real(lon_r4,8)
          end if

          ! check if the slanted lat/lon is inside the domain
          call latlonChecks (obsSpaceData, stateVector%hco, & ! IN
                             headerIndex, rejectOutsideObs, & ! IN
                             latLev_T, lonLev_T,            & ! IN/OUT
                             latLev_M, lonLev_M,            & ! IN/OUT
                             latLev_S, lonLev_S )             ! IN/OUT

          ! put the lat/lon from TH/MM levels to varLevIndex
          do varLevIndex = allVarLevBeg(1), stateVector%numVarLev
            if (varLevIndex == 0) then
              varLevel = 'SF'
            else
              levIndex = gsv_getLevFromVarLev(stateVector,varLevIndex)
              varLevel = vnl_varLevelFromVarname(gsv_getVarNameFromVarLev(stateVector,varLevIndex))
            end if

            if ( varLevel == 'TH' ) then
              latColumn(headerUsedIndex,varLevIndex) = latLev_T(levIndex)
              lonColumn(headerUsedIndex,varLevIndex) = lonLev_T(levIndex)
            else if ( varLevel == 'MM' ) then
              latColumn(headerUsedIndex,varLevIndex) = latLev_M(levIndex)
              lonColumn(headerUsedIndex,varLevIndex) = lonLev_M(levIndex)
            else if ( varLevel == 'SF' ) then
              latColumn(headerUsedIndex,varLevIndex) = latLev_S
              lonColumn(headerUsedIndex,varLevIndex) = lonLev_S
            else
              call utl_abort('setupInterpInfo2dFields: unknown value of varLevel')
            end if

          end do

        end do header_loop3

        ! MPI communication for the slant-path lat/lon

        maxkCount = maxval(stateVector%allVarLevCount(:) + stateVector%allVarLevBeg(:) - allVarLevBeg(:))
        numVarLevToSend = min(mmpi_nprocs,stateVector%numVarLev)

        allocate(lat_recv_r8(numHeaderUsedMax,mmpi_nprocs))
        lat_recv_r8(:,:) = 0.0d0
        allocate(lat_send_r8(numHeaderUsedMax,mmpi_nprocs))
        lat_send_r8(:,:) = 0.0d0
        allocate(lon_recv_r8(numHeaderUsedMax,mmpi_nprocs))
        lon_recv_r8(:,:) = 0.0d0
        allocate(lon_send_r8(numHeaderUsedMax,mmpi_nprocs))
        lon_send_r8(:,:) = 0.0d0

        ! only send the data from tasks with data, same amount to all
        sendsizes(:) = 0
        do procIndex = 1, numVarLevToSend
          sendsizes(procIndex) = numHeaderUsed
        end do
        senddispls(1) = 0
        do procIndex = 2, mmpi_nprocs
          senddispls(procIndex) = senddispls(procIndex-1) + numHeaderUsedMax
        end do

        recvdispls(1) = 0
        do procIndex = 2, mmpi_nprocs
          recvdispls(procIndex) = recvdispls(procIndex-1) + numHeaderUsedMax
        end do

        ! loop to send (at most) 1 level to (at most) all other mpi tasks
        do varLevIndexCount = 1, maxkCount

          sendsizes(:) = 0
          do procIndex = 1, mmpi_nprocs
            ! compute varLevIndex value being sent
            varLevIndex = varLevIndexCount + allVarLevBeg(procIndex) - 1
            if ( varLevIndex <= stateVector%allVarLevEnd(procIndex) ) then
              if( procIndex > numVarLevToSend ) then
                write(*,*) 'procIndex, numVarLevToSend = ', procIndex, numVarLevToSend
                call utl_abort('ERROR: with numVarLevToSend?')
              end if

              lat_send_r8(1:numHeaderUsed,procIndex) = latColumn(1:numHeaderUsed,varLevIndex)
              lon_send_r8(1:numHeaderUsed,procIndex) = lonColumn(1:numHeaderUsed,varLevIndex)
              sendsizes(procIndex) = numHeaderUsed
            else
              sendsizes(procIndex) = 0
            end if
          end do

          ! all tasks recv only from those with data
          varLevIndex = varLevIndexCount + myVarLevBeg - 1
          if ( varLevIndex <= stateVector%myVarLevEnd ) then
            do procIndex = 1, mmpi_nprocs
              recvsizes(procIndex) = allNumHeaderUsed(stepIndex,procIndex)
            end do
          else
            recvsizes(:) = 0
          end if

          call mmpi_alltoallv(lat_send_r8, sendsizes, senddispls, &
                              lat_recv_r8, recvsizes, recvdispls)
          call mmpi_alltoallv(lon_send_r8, sendsizes, senddispls, &
                              lon_recv_r8, recvsizes, recvdispls)

          do procIndex = 1, mmpi_nprocs
            ! all tasks copy the received step data into correct slot
            varLevIndex = varLevIndexCount + myVarLevBeg - 1
            if ( varLevIndex <= stateVector%myVarLevEnd ) then
              interpInfo%stepProcData(procIndex,stepIndex)%allLat(:,varLevIndex) = &
                   lat_recv_r8(1:allNumHeaderUsed(stepIndex,procIndex),procIndex)
              interpInfo%stepProcData(procIndex,stepIndex)%allLon(:,varLevIndex) = &
                   lon_recv_r8(1:allNumHeaderUsed(stepIndex,procIndex),procIndex)
            end if
          end do

        end do ! varLevIndexCount

        deallocate(lon_send_r8)
        deallocate(lon_recv_r8)
        deallocate(lat_send_r8)
        deallocate(lat_recv_r8)

        deallocate(latLev_T)
        deallocate(lonLev_T)
        deallocate(latLev_M)
        deallocate(lonLev_M)

      else ! not doSlantPath

        allocate(latLev_T(1))
        allocate(lonLev_T(1))
        allocate(latLev_M(1))
        allocate(lonLev_M(1))
        latLev_T(:) = 0.0d0
        lonLev_T(:) = 0.0d0
        latLev_M(:) = 0.0d0
        lonLev_M(:) = 0.0d0

        do headerUsedIndex = 1, numHeaderUsed
          headerIndex = headerIndexVec(headerUsedIndex,stepIndex)

          !- Get LatLon of observation location
          lat_r4 = real(obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex), 4)
          lon_r4 = real(obs_headElem_r(obsSpaceData, OBS_LON, headerIndex), 4)
          if (lon_r4 <  0.0          ) lon_r4 = lon_r4 + 2.0*MPC_PI_R4
          if (lon_r4 >= 2.0*MPC_PI_R4) lon_r4 = lon_r4 - 2.0*MPC_PI_R4

          latLev_T(:) = real(lat_r4,8)
          lonLev_T(:) = real(lon_r4,8)
          latLev_M(:) = real(lat_r4,8)
          lonLev_M(:) = real(lon_r4,8)

          ! check if the lat/lon is inside the domain
          call latlonChecks (obsSpaceData, stateVector%hco, & ! IN
                             headerIndex, rejectOutsideObs, & ! IN
                             latLev_T, lonLev_T,            & ! IN/OUT
                             latLev_M, lonLev_M )             ! IN/OUT

          latColumn(headerUsedIndex,allVarLevBeg(1)) = latLev_T(1)
          lonColumn(headerUsedIndex,allVarLevBeg(1)) = lonLev_T(1)
        end do

        ! gather geographical lat, lon positions of observations from all processors
        call mmpi_allGather(latColumn(:,allVarLevBeg(1)), allLatOneLev(:,:))
        call mmpi_allGather(lonColumn(:,allVarLevBeg(1)), allLonOneLev(:,:))

        k_loop: do varLevIndex = myVarLevBeg, statevector%myVarLevEnd
          do procIndex = 1, mmpi_nprocs
            interpInfo%stepProcData(procIndex,stepIndex)%allLat(:,varLevIndex) = &
                 allLatOneLev(1:allNumHeaderUsed(stepIndex,procIndex),procIndex)
            interpInfo%stepProcData(procIndex,stepIndex)%allLon(:,varLevIndex) = &
                 allLonOneLev(1:allNumHeaderUsed(stepIndex,procIndex),procIndex)
          end do
        end do k_loop

        deallocate(latLev_T)
        deallocate(lonLev_T)
        deallocate(latLev_M)
        deallocate(lonLev_M)

      end if ! doSlantPath

      deallocate(lonColumn)
      deallocate(latColumn)

    end do step_loop2

    if ( gsv_isAllocated(stateVector_1Step) .and. lastCall ) then
      write(*,*) 'setupInterpInfo2dFields: deallocate height3D fields'
      call gsv_deallocate(stateVector_1Step)
    end if
    deallocate(footprintRadiusVec_r4)

    write(*,*) 'setupInterpInfo2dFields: latlonChecks and lat/lon MPI comm finished.'

    allocate(allHeaderIndex(numHeaderUsedMax,numStep,mmpi_nprocs))
    ! gather the headerIndexVec arrays onto all processors
    call mmpi_allGather(headerIndexVec, allHeaderIndex)

    do procIndex = 1, mmpi_nprocs
      do stepIndex = 1, numStep
        do headerIndex = 1, allNumHeaderUsed(stepIndex,procIndex)
          interpInfo%stepProcData(procIndex, stepIndex)%allHeaderIndex(headerIndex) = allHeaderIndex(headerIndex,stepIndex,procIndex)
        end do
      end do
    end do

    ! create kdtree to use in footprint operator, if any footprint radius > 0.
    interpInfo%inputStateVectorType = inputStateVectorType
    if ( any(allFootprintRadius_r4(:,:,:) > 0.0) ) then
      if ( (inputStateVectorType == 'nl' .and. .not. associated(tree_nl))   .or. &
           (inputStateVectorType == 'tl' .and. .not. associated(tree_tlad)) .or. &
           (inputStateVectorType == 'ad' .and. .not. associated(tree_tlad)) ) then

        write(*,*) 'setupInterpInfo2dFields: start creating kdtree for inputStateVectorType=', &
                   inputStateVectorType
        call msg_memUsage('setupInterpInfo2dFields')

        allocate(positionArray(3,statevector%hco%ni*statevector%hco%nj))

        gridIndex = 0
        do latIndex = 1, statevector%hco%nj
          do lonIndex = 1, statevector%hco%ni
            gridIndex = gridIndex + 1
            lat = real(stateVector % hco % lat2d_4(lonIndex,latIndex), 8)
            lon = real(stateVector % hco % lon2d_4(lonIndex,latIndex), 8)

            positionArray(:,gridIndex) = kdtree2_3dPosition(lon, lat)

          end do
        end do

        nullify(tree)
        tree => kdtree2_create(positionArray, sort=.false., rearrange=.true.)

        if ( inputStateVectorType == 'nl' ) then
          tree_nl => tree
        else
          tree_tlad => tree
        end if

        deallocate(positionArray)

        write(*,*) 'setupInterpInfo2dFields: done creating kdtree for inputStateVectorType = ', &
                   inputStateVectorType
        call msg_memUsage('setupInterpInfo2dFields')

      end if
    end if

    do stepIndex = 1, numStep
      !$OMP PARALLEL DO PRIVATE (procIndex, varLevIndex, headerIndex, lat_deg_r4, lon_deg_r4, ierr, &
      !$OMP xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, subGridIndex, numSubGridsForInterp, subGridForInterp, &
      !$OMP lat, lon, latRot, lonRot, footprintRadius_r4, numGridpt)
      do procIndex = 1, mmpi_nprocs
        do varLevIndex = myVarLevBeg, statevector%myVarLevEnd
          do headerIndex = 1, allNumHeaderUsed(stepIndex,procIndex)

            ! Compute the rotated lat/lon, needed for the winds

            lat_deg_r4 = real(interpInfo%stepProcData(procIndex,stepIndex)%allLat(headerIndex,varLevIndex) *  &
                         MPC_DEGREES_PER_RADIAN_R8)
            lon_deg_r4 = real(interpInfo%stepProcData(procIndex,stepIndex)%allLon(headerIndex,varLevIndex) *  &
                         MPC_DEGREES_PER_RADIAN_R8)
            ierr = gpos_getPositionXY(stateVector%hco%EZscintID,   &
                                      xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                                      lat_deg_r4, lon_deg_r4, subGridIndex)

            if ( subGridIndex == 3 ) then
              ! both subGrids involved in interpolation, so first treat subGrid 1
              numSubGridsForInterp = 2
              subGridIndex = 1
            else
              ! only 1 subGrid involved in interpolation
              numSubGridsForInterp = 1
            end if

            do subGridForInterp = 1, numSubGridsForInterp

              if ( subGridForInterp == 1 ) then
                ! when only 1 subGrid involved, subGridIndex can be 1 or 2
              else
                ! when 2 subGrids, subGridIndex is set to 1 for 1st iteration, 2 for second
                subGridIndex = 2
              end if

              if (interpInfo%hco%rotated .and.  &
                  gsv_varExist(varName='UU') .and.  &
                  gsv_varExist(varName='VV')) then
                lat = interpInfo%stepProcData(procIndex,stepIndex)%allLat(headerIndex,varLevIndex)
                lon = interpInfo%stepProcData(procIndex,stepIndex)%allLon(headerIndex,varLevIndex)
                call uvr_RotateLatLon( interpInfo%uvr,   & ! INOUT
                                       subGridIndex,     & ! IN
                                       latRot, lonRot,   & ! OUT (radians)
                                       lat, lon,         & ! IN  (radians)
                                       'ToLatLonRot')      ! IN
                interpInfo%stepProcData(procIndex,stepIndex)%allLatRot(subGridIndex,headerIndex,varLevIndex) = latRot
                interpInfo%stepProcData(procIndex,stepIndex)%allLonRot(subGridIndex,headerIndex,varLevIndex) = lonRot
              end if

            end do ! subGridForInterp

            ! Count total number of grid points for allocating interp depot

            footprintRadius_r4 = allFootprintRadius_r4(headerIndex,stepIndex,procIndex)

            call s2c_setupHorizInterp(footprintRadius_r4, interpInfo, &
                                      stateVector, headerIndex, varLevIndex, stepIndex, &
                                      procIndex, numGridpt)

            ! for now, just store the number of gridpts for each obs in depotIndexEnd
            if ( (subGridIndex == 1) .or. (subGridIndex == 2) ) then
              ! indices for only 1 subgrid, other will have zeros
              interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex,headerIndex,varLevIndex) = &
                   numGridpt(subGridIndex)
            else
              ! locations on both subGrids will be averaged
              interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(1,headerIndex,varLevIndex) = numGridpt(1)
              interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(2,headerIndex,varLevIndex) = numGridpt(2)
            end if

          end do ! headerIndex
        end do ! varLevIndex
      end do ! procIndex
      !$OMP END PARALLEL DO
    end do ! stepIndex

    numGridptTotal = 0
    do stepIndex = 1, numStep
      do procIndex = 1, mmpi_nprocs
        do varLevIndex = myVarLevBeg, statevector%myVarLevEnd
          do headerIndex = 1, allNumHeaderUsed(stepIndex,procIndex)
            do subGridIndex = 1, interpInfo%hco%numSubGrid
              if ( interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex,headerIndex,varLevIndex) /= -1 ) then
                interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex,headerIndex,varLevIndex) = &
                     numGridptTotal + 1
                numGridptTotal = numGridptTotal + &
                     interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex,headerIndex,varLevIndex)
                interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex,headerIndex,varLevIndex) = &
                     numGridptTotal
              end if
            end do ! subGridIndex
          end do ! headerIndex
        end do ! varLevIndex
      end do ! procIndex
    end do ! stepIndex

    deallocate(allHeaderIndex)

    call msg_memUsage('setupInterpInfo2dFields')

    ! now that we know the size, allocate main arrays for storing interpolation information
    write(*,*) 'setupInterpInfo2dFields: numGridptTotal = ', numGridptTotal
    call mmpi_allReduce(numGridptTotal, numGridptTotalMpi, mmpi_sum)
    write(*,*) 'setupInterpInfo2dFields: numGridptTotal, numGridptTotalMpi = ',  &
                                         numGridptTotal, numGridptTotalMpi

    allocate(interpInfo%latIndexDepot(numGridptTotal))
    allocate(interpInfo%lonIndexDepot(numGridptTotal))
    allocate(interpInfo%interpWeightDepot(numGridptTotal))

    call utl_tmg_start(33,'------s2c_SetupWeights')
    !$OMP PARALLEL DO PRIVATE (procIndex, stepIndex, varLevIndex, headerIndex, footprintRadius_r4, numGridpt)
    do procIndex = 1, mmpi_nprocs
      do stepIndex = 1, numStep
        do varLevIndex = myVarLevBeg, statevector%myVarLevEnd

          do headerIndex = 1, allNumHeaderUsed(stepIndex,procIndex)

            footprintRadius_r4 = allFootprintRadius_r4(headerIndex, stepIndex, procIndex)

            call s2c_setupHorizInterp(footprintRadius_r4, interpInfo, stateVector, &
                                      headerIndex, varLevIndex, stepIndex, procIndex, numGridpt)

          end do ! headerIndex

        end do ! varLevIndex
      end do ! stepIndex
    end do ! procIndex
    !$OMP END PARALLEL DO
    call utl_tmg_stop(33)

    ! reject obs in obsSpaceData if any processor has zero weight
    ! called when a mask exists to catch land contaminated ocean obs
    if ( stateVector%oceanMask%maskPresent ) then
      call s2c_rejectZeroWeightObs(interpInfo,obsSpaceData,myVarLevBeg,stateVector%myVarLevEnd)
    end if

    ! on the last call, deallocate the tree_nl/tree_tlad
    if ( lastCall ) then
      if ( inputStateVectorType == 'nl' .and. associated(tree_nl) ) then
        call kdtree2_destroy(tree_nl)
      else if ( (inputStateVectorType == 'tl' .or. inputStateVectorType == 'ad') .and. &
                associated(tree_tlad) ) then
        call kdtree2_destroy(tree_tlad)
      end if
    end if

    ! Count the number of TOVS using footprint operator on one level
    if ( useFootprintForTovs ) then
      numTovsUsingFootprint = 0
      numAllTovs = 0
      procIndex = mmpi_myid + 1
      do stepIndex = 1, numStep
        do headerUsedIndex = 1, allNumHeaderUsed(stepIndex,procIndex)
          footprintRadius_r4 = allFootprintRadius_r4(headerUsedIndex, stepIndex, procIndex)
          headerIndex = headerIndexVec(headerUsedIndex,stepIndex)
          codeType = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)

          if ( tvs_isIdBurpTovs(codeType) ) then
            if ( footprintRadius_r4 > 0.0 ) numTovsUsingFootprint = numTovsUsingFootprint + 1
            numAllTovs = numAllTovs + 1
          end if
        end do
      end do

      if ( numAllTovs > 0 ) then
        write(*,'(A,2(I5,A2),F5.1,A)') 'setupInterpInfo2dFields: numTovsUsingFootprint/numAllTovs=', &
                       numTovsUsingFootprint, ' /', numAllTovs, ' (', &
                       real(numTovsUsingFootprint) / real(numAllTovs) * 100.0, '%)'
      end if
    end if

    deallocate(allFootprintRadius_r4)
    deallocate(allLonOneLev)
    deallocate(allLatOneLev)

    deallocate(headerIndexVec)
    deallocate(allNumHeaderUsed)

    interpInfo%initialized = .true.

    call msg_memUsage('setupInterpInfo2dFields')
    write(*,*) 'setupInterpInfo2dFields: FINISHED'

  end subroutine setupInterpInfo2dFields

  !---------------------------------------------------------
  ! s2c_tl
  !---------------------------------------------------------
  subroutine s2c_tl(statevector_in, columnAnlInc, columnTrlOnAnlIncLev,  &
                    obsSpaceData, beSilent_opt)
    !
    ! :Purpose: Tangent linear version of the horizontal
    !           interpolation, used for the increment (or perturbations).
    !
    implicit none

    ! Arguments:
    type(struct_gsv), target, intent(in)    :: stateVector_in
    type(struct_obs)        , intent(inout) :: obsSpaceData
    type(struct_columnData) , intent(inout) :: columnAnlInc
    type(struct_columnData) , intent(inout) :: columnTrlOnAnlIncLev
    logical,        optional, intent(in)    :: beSilent_opt

    ! Locals:
    type(struct_gsv), pointer  :: stateVector

    call utl_tmg_start(30,'--StateToColumn')
    call utl_tmg_start(38,'----s2c_TL')

    call utl_tmg_start(40,'------s2c_barrier_TL')
    call mmpi_barrier
    call utl_tmg_stop(40)

    if ( mmpi_myid == 0 ) write(*,*) 's2c_tl: Horizontal interpolation StateVector --> ColumnData'

    call readNml()

    if ( .not. gsv_isAllocated(stateVector_in) ) then
      call utl_abort('s2c_tl: stateVector must be allocated')
    end if

    if (trim(mpiMode) == '2DFIELDS' .and. interpInfo_tlad%initialized) then

      if (.not. hco_equal(interpInfo_tlad%hco,stateVector_in%hco)) then
        write(*,*) 's2c_tl: WARNING! Current hco grid parameters differ from allocated interpInfo_tlad!'
        write(*,*) 's2c_tl: InterpInfo_tlad will be deallocated.'
        call s2c_deallocInterpInfo(inputStateVectorType='tlad')
      end if

    else if(trim(mpiMode) == 'TILES' .and. interpInfoTiles_tlad%initialized) then

      if (.not. hco_equal(interpInfoTiles_tlad%hco,stateVector_in%hco)) then
        write(*,*) 's2c_tl: WARNING! Current hco grid parameters differ from allocated interpInfoTiles_tlad!'
        write(*,*) 's2c_tl: InterpInfoTiles_tlad will be deallocated.'
        call s2c_deallocInterpInfo(inputStateVectorType='tlad')
      end if

    end if

    ! check the column and statevector have same nk/varNameList
    call checkColumnStatevectorMatch(columnAnlInc,statevector_in)

    ! if we only compute Height and Pressure on column, make copy without them
    if ( calcHeightPressIncrOnColumn ) then
      allocate(stateVector)
      call gsv_allocate( stateVector, statevector_in%numstep, &
                         statevector_in%hco, statevector_in%vco, &
                         mpi_local_opt=.true., &
                         dataKind_opt=gsv_getDataKind(statevector_in), &
                         allocHeight_opt=.false., allocPressure_opt=.false. )
      call gsv_copy(stateVector_in, stateVector, allowVarMismatch_opt=.true.)
    else
      stateVector => stateVector_in

      ! calculate delP_T/delP_M and del Z_T/Z_M on the grid
      call gvt_transform( statevector, 'ZandP_tl' )
    end if

    if (trim(mpiMode) == '2DFIELDS') then
      call tl2dFields(stateVector, obsSpaceData, columnAnlInc)
    else if (trim(mpiMode) == 'TILES') then
      call tlTiles(stateVector, obsSpaceData, columnAnlInc, beSilent_opt)
    else
      call utl_abort('s2c_tl: invalid value of mpiMode = '//trim(mpiMode))
    end if

    if (calcHeightPressIncrOnColumn) then
      ! calculate delP_T/delP_M and  del Z_T/Z_M on the columns
      call czp_calcZandP_tl(columnAnlInc, columnTrlOnAnlIncLev)
    end if

    if (calcHeightPressIncrOnColumn) then
      call gsv_deallocate(stateVector)
      deallocate(stateVector)
    end if

    if (slantPath_TO_tlad) then
      call pressureProfileMonotonicityCheck(obsSpaceData, columnTrlOnAnlIncLev)
    end if

    call utl_tmg_stop(38)
    call utl_tmg_stop(30)

  end subroutine s2c_tl

  !---------------------------------------------------------
  ! tl2dFields
  !---------------------------------------------------------
  subroutine tl2dFields(stateVector, obsSpaceData, columnAnlInc)
    !
    ! :Purpose: Tangent-linear version of the horizontal interpolation,
    !           used for increments or perturbations.
    !
    implicit none

    ! Arguments:
    type(struct_gsv),           intent(inout) :: stateVector
    type(struct_obs),           intent(inout) :: obsSpaceData
    type(struct_columnData),    intent(inout) :: columnAnlInc

    ! Locals:
    type(struct_gsv) :: stateVector_VarsLevs
    integer :: varLevIndex, varLevIndex2, levIndex, kCount, stepIndex, numStep, myVarLevEndExtended
    integer :: headerIndex, numHeader, numHeaderMax, yourNumHeader
    integer :: procIndex, headerUsedIndex
    real(8) :: weight
    real(8), pointer     :: allCols_ptr(:,:)
    real(pre_incrReal), pointer :: ptr4d(:,:,:,:)
    real(pre_incrReal), pointer :: ptr3d_UV(:,:,:)
    real(8), allocatable :: cols_hint(:,:,:)
    real(8), allocatable :: cols_send(:,:)
    real(8), allocatable :: cols_recv(:,:)
    real(8), allocatable :: cols_send_1proc(:)
    logical              :: rejectOutsideObs
    character(len=4)     :: varName
    character(len=4), pointer :: varNames(:)

    nullify(varNames)
    call gsv_varNamesList(varNames, statevector)
    if (statevector%mpi_distribution == 'None') then
      call gsv_allocate( statevector_VarsLevs, statevector%numstep, &
                         statevector%hco, statevector%vco,          &
                         mpi_local_opt=.false., mpi_distribution_opt='None', &
                         dataKind_opt=gsv_getDataKind(statevector), &
                         varNames_opt=varNames )
    else
      call gsv_allocate( statevector_VarsLevs, statevector%numstep, &
                         statevector%hco, statevector%vco,          &
                         mpi_local_opt=.true., mpi_distribution_opt='VarsLevs', &
                         dataKind_opt=gsv_getDataKind(statevector), &
                         varNames_opt=varNames )
    end if
    deallocate(varNames)
    call gsv_transposeTilesToVarsLevs( statevector, statevector_VarsLevs )

    numStep = stateVector_VarsLevs%numStep
    numHeader = obs_numheader(obsSpaceData)
    call mmpi_allReduce(numHeader, numHeaderMax, mmpi_max)

    if ( .not. interpInfo_tlad%initialized ) then
      rejectOutsideObs = .false.
      call utl_tmg_stop(38)
      call utl_tmg_start(31,'----s2c_Setups')
      call setupInterpInfo2dFields( interpInfo_tlad, obsSpaceData, stateVector_VarsLevs,  &
                                    1, numHeader, timeInterpType_tlad,  rejectOutsideObs, &
                                    inputStateVectorType='tl' )
      call utl_tmg_stop(31)
      call utl_tmg_start(38,'----s2c_TL')
    end if

    ! arrays for interpolated column for 1 level/variable and each time step
    allocate(cols_hint(maxval(interpInfo_tlad%allNumHeaderUsed),numStep,mmpi_nprocs))
    cols_hint(:,:,:) = 0.0d0

    ! arrays for sending/receiving time interpolated column for 1 level/variable
    allocate(cols_send(numHeaderMax,mmpi_nprocs))
    cols_send(:,:) = 0.0d0

    allocate(cols_recv(numHeaderMax,mmpi_nprocs))
    cols_recv(:,:) = 0.0d0

    allocate(cols_send_1proc(numHeaderMax))

    ! set contents of column to zero
    call col_zero(columnAnlInc)

    call gsv_getField(stateVector_VarsLevs, ptr4d)

    myVarLevEndExtended = stateVector_VarsLevs%myVarLevBeg + maxval(stateVector_VarsLevs%allVarLevCount(:)) - 1

    kCount = 0
    k_loop: do varLevIndex = stateVector_VarsLevs%myVarLevBeg, myVarLevEndExtended

      kCount = kCount + 1

      if ( varLevIndex <= stateVector_VarsLevs%myVarLevEnd ) then
        varName = gsv_getVarNameFromVarLev(statevector,varLevIndex)

        if ( varName == 'UU' .or. varName == 'VV' ) then
          call gsv_getFieldUV(stateVector_VarsLevs,ptr3d_UV,varLevIndex)
        end if

        call utl_tmg_start(39,'------s2c_TL_Hinterp')
        !$OMP PARALLEL DO PRIVATE (stepIndex, procIndex, yourNumHeader, headerIndex)
        step_loop: do stepIndex = 1, numStep
          if ( maxval(interpInfo_tlad%allNumHeaderUsed(stepIndex,:)) == 0 ) cycle step_loop

          ! interpolate to the columns destined for all procs for all steps and one lev/var
          do procIndex = 1, mmpi_nprocs
            yourNumHeader = interpInfo_tlad%allNumHeaderUsed(stepIndex,procIndex)
            if ( yourNumHeader > 0 ) then
              if ( varName == 'UU' ) then
                call myezuvint_tl( cols_hint(1:yourNumHeader,stepIndex,procIndex), 'UU',  &
                                   ptr4d(:,:,varLevIndex,stepIndex), ptr3d_UV(:,:,stepIndex),  &
                                   interpInfo_tlad, varLevIndex, stepIndex, procIndex )
              else if ( varName == 'VV' ) then
                call myezuvint_tl( cols_hint(1:yourNumHeader,stepIndex,procIndex), 'VV',  &
                                   ptr3d_UV(:,:,stepIndex), ptr4d(:,:,varLevIndex,stepIndex),  &
                                   interpInfo_tlad, varLevIndex, stepIndex, procIndex )
              else
                call myezsint_tl( cols_hint(1:yourNumHeader,stepIndex,procIndex),  &
                                  ptr4d(:,:,varLevIndex,stepIndex), interpInfo_tlad, varLevIndex, &
                                  stepIndex, procIndex )
              end if
            end if
          end do

        end do step_loop
        !$OMP END PARALLEL DO
        call utl_tmg_stop(39)

        ! interpolate in time to the columns destined for all procs and one level/variable
        do procIndex = 1, mmpi_nprocs
          cols_send_1proc(:) = 0.0d0
          do stepIndex = 1, numStep
            !$OMP PARALLEL DO PRIVATE (headerUsedIndex, headerIndex, weight)
            do headerUsedIndex = 1, interpInfo_tlad%allNumHeaderUsed(stepIndex, procIndex)
              headerIndex = interpInfo_tlad%stepProcData(procIndex, stepIndex)%allHeaderIndex(headerUsedIndex)
              weight = oti_getTimeInterpWeightMpiGlobal(interpInfo_tlad%oti,  &
                                                        headerIndex,stepIndex,procIndex)
              cols_send_1proc(headerIndex) = cols_send_1proc(headerIndex) &
                            + weight * cols_hint(headerUsedIndex,stepIndex,procIndex)

            end do
            !$OMP END PARALLEL DO
          end do
          cols_send(:,procIndex) = cols_send_1proc(:)
        end do

      else

        ! this value of k does not exist on this mpi task
        cols_send(:,:) = 0.0

      end if ! if varLevIndex <= myVarLevEnd

      ! mpi communication: alltoall for one level/variable
      if(mmpi_nprocs > 1) then
        call mmpi_alltoall(cols_send, cols_recv)
      else
        cols_recv(:,1) = cols_send(:,1)
      end if

      ! reorganize ensemble of distributed columns
      !$OMP PARALLEL DO PRIVATE (procIndex, varLevIndex2, varName, levIndex, allCols_ptr, headerIndex)
      proc_loop: do procIndex = 1, mmpi_nprocs
        ! This is varLevIndex value of source (can be different for destination)
        varLevIndex2 = statevector_VarsLevs%allVarLevBeg(procIndex) + kCount - 1
        if ( varLevIndex2 > stateVector_VarsLevs%allVarLevEnd(procIndex) ) cycle proc_loop

        ! Figure out which variable/level of destination
        varName = gsv_getVarNameFromVarLev(statevector,varLevIndex2)
        levIndex = gsv_getLevFromVarLev(statevector,varLevIndex2)
        allCols_ptr => col_getAllColumns(columnAnlInc,varName)

        do headerIndex = 1, numHeader
          allCols_ptr(levIndex,headerIndex) = cols_recv(headerIndex,procIndex)
        end do
      end do proc_loop
      !$OMP END PARALLEL DO

    end do k_loop

    deallocate(cols_hint)
    deallocate(cols_send)
    deallocate(cols_recv)
    deallocate(cols_send_1proc)

    call gsv_deallocate( statevector_VarsLevs )

  end subroutine tl2dFields

  !---------------------------------------------------------
  ! s2c_ad
  !---------------------------------------------------------
  subroutine s2c_ad(statevector_out, columnAnlInc, columnTrlOnAnlIncLev,  &
                    obsSpaceData, beSilent_opt)
    !
    ! :Purpose: Adjoint version of the horizontal interpolation,
    !           used for the cost function gradient with respect to the increment.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), target, intent(inout) :: stateVector_out
    type(struct_obs)        , intent(inout) :: obsSpaceData
    type(struct_columnData) , intent(inout) :: columnAnlInc
    type(struct_columnData) , intent(inout) :: columnTrlOnAnlIncLev
    logical,        optional, intent(in)    :: beSilent_opt

    ! Locals:
    type(struct_gsv), pointer  :: stateVector

    call utl_tmg_start(30,'--StateToColumn')
    call utl_tmg_start(41,'----s2c_AD')

    call utl_tmg_start(43,'------s2c_barrier_AD')
    call mmpi_barrier
    call utl_tmg_stop(43)

    if(mmpi_myid == 0) write(*,*) 's2c_ad: Adjoint of horizontal interpolation StateVector --> ColumnData'

    if ( .not. gsv_isAllocated(stateVector_out) ) then
      call utl_abort('s2c_ad: stateVector must be allocated')
    end if

    if (trim(mpiMode) == '2DFIELDS' .and. interpInfo_tlad%initialized) then

      if (.not. hco_equal(interpInfo_tlad%hco,stateVector_out%hco)) then
        write(*,*) 's2c_ad: WARNING! Current hco grid parameters differ from allocated interpInfo_tlad!'
        write(*,*) 's2c_ad: InterpInfo_tlad will be deallocated.'
        call s2c_deallocInterpInfo(inputStateVectorType='tlad')
      end if

    else if(trim(mpiMode) == 'TILES' .and. interpInfoTiles_tlad%initialized) then

      if (.not. hco_equal(interpInfoTiles_tlad%hco,stateVector_out%hco)) then
        write(*,*) 's2c_ad: WARNING! Current hco grid parameters differ from allocated interpInfoTiles_tlad!'
        write(*,*) 's2c_ad: InterpInfoTiles_tlad will be deallocated.'
        call s2c_deallocInterpInfo(inputStateVectorType='tlad')
      end if

    end if

    ! if we only compute Height and Pressure on column, make copy without them
    if (calcHeightPressIncrOnColumn) then
      allocate(stateVector)
      call gsv_allocate(stateVector, statevector_out%numstep, &
                        statevector_out%hco, statevector_out%vco, &
                        mpi_local_opt=.true., &
                        dataKind_opt=gsv_getDataKind(statevector_out), &
                        allocHeight_opt=.false., allocPressure_opt=.false.)
      ! Adjoint of calculate del Z_T/Z_M and delP_T/delP_M on the columns
      call czp_calcZandP_ad(columnAnlInc, columnTrlOnAnlIncLev)
    else
      stateVector => stateVector_out
    end if

    if (trim(mpiMode) == '2DFIELDS') then
      call ad2dFields(stateVector, obsSpaceData, columnAnlInc)
    else if (trim(mpiMode) == 'TILES') then
      call adTiles(stateVector, obsSpaceData, columnAnlInc, beSilent_opt)
    else
      call utl_abort('s2c_ad: invalid value of mpiMode = '//trim(mpiMode))
    end if

    if (calcHeightPressIncrOnColumn) then
      call gsv_zero(statevector_out)
      call gsv_copy(stateVector, stateVector_out, allowVarMismatch_opt=.true.)
    else
      ! Adjoint of calculate del Z_T/Z_M and delP_T/delP_M on the grid
      call gvt_transform( statevector, 'ZandP_ad' )
    end if

    if (slantPath_TO_tlad) then
      call pressureProfileMonotonicityCheck(obsSpaceData, columnTrlOnAnlIncLev)
    end if

    if (calcHeightPressIncrOnColumn) then
      call gsv_deallocate( statevector )
      deallocate(stateVector)
    end if

    call utl_tmg_stop(41)
    call utl_tmg_stop(30)

  end subroutine s2c_ad

  !---------------------------------------------------------
  ! ad2dFields
  !---------------------------------------------------------
  subroutine ad2dFields(stateVector, obsSpaceData, columnAnlInc)
    !
    ! :Purpose: Adjoint version of the horizontal interpolation,
    !           used for the cost function gradient with respect to the increment.
    !
    implicit none

    ! Arguments:
    type(struct_gsv),        intent(inout) :: stateVector
    type(struct_obs),        intent(inout) :: obsSpaceData
    type(struct_columnData), intent(inout) :: columnAnlInc

    ! Locals:
    type(struct_gsv) :: stateVector_VarsLevs
    integer :: varLevIndex, varLevIndex2, kCount, levIndex, stepIndex, numStep, myVarLevEndExtended
    integer :: headerIndex, numHeader, numHeaderMax, yourNumHeader
    integer :: procIndex, headerUsedIndex
    character(len=4)     :: varName
    real(8) :: weight
    real(pre_incrReal), pointer :: ptr4d(:,:,:,:), ptr3d_UV(:,:,:)
    real(8), pointer     :: allCols_ptr(:,:)
    real(8), allocatable :: cols_hint(:,:,:)
    real(8), allocatable :: cols_send(:,:)
    real(8), allocatable :: cols_recv(:,:)
    logical              :: rejectOutsideObs
    character(len=4), pointer :: varNames(:)

    nullify(varNames)
    call gsv_varNamesList(varNames, statevector)
    call gsv_allocate( statevector_VarsLevs, statevector%numstep, &
                       statevector%hco, statevector%vco,          &
                       mpi_local_opt=.true., mpi_distribution_opt='VarsLevs', &
                       dataKind_opt=gsv_getDataKind(statevector), &
                       varNames_opt=varNames )
    deallocate(varNames)
    call gsv_zero( statevector_VarsLevs )

    numStep = stateVector_VarsLevs%numStep
    numHeader = obs_numheader(obsSpaceData)
    call mmpi_allReduce(numHeader, numHeaderMax, mmpi_max)

    if ( .not. interpInfo_tlad%initialized ) then
      rejectOutsideObs = .false.
      call utl_tmg_stop(41)
      call utl_tmg_start(31,'----s2c_Setups')
      call setupInterpInfo2dFields( interpInfo_tlad, obsSpaceData, stateVector_VarsLevs,  &
                                    1, numHeader, timeInterpType_tlad, rejectOutsideObs,  &
                                    inputStateVectorType='ad' )
      call utl_tmg_stop(31)
      call utl_tmg_start(41,'----s2c_AD')
    end if

    ! arrays for interpolated column for 1 level/variable and each time step
    allocate(cols_hint(maxval(interpInfo_tlad%allNumHeaderUsed),numStep,mmpi_nprocs))
    cols_hint(:,:,:) = 0.0d0

    ! arrays for sending/receiving time interpolated column for 1 level/variable
    allocate(cols_send(numHeaderMax,mmpi_nprocs))
    cols_send(:,:) = 0.0d0

    allocate(cols_recv(numHeaderMax,mmpi_nprocs))
    cols_recv(:,:) = 0.0d0

    call gsv_getField(stateVector_VarsLevs,ptr4d)
    myVarLevEndExtended = stateVector_VarsLevs%myVarLevBeg + maxval(stateVector_VarsLevs%allVarLevCount(:)) - 1

    kCount = 0
    k_loop: do varLevIndex = stateVector_VarsLevs%myVarLevBeg, myVarLevEndExtended

      kCount = kCount + 1

      ! reorganize ensemble of distributed columns
      !$OMP PARALLEL DO PRIVATE (procIndex, varLevIndex2, varName, levIndex, allCols_ptr, headerIndex)
      proc_loop: do procIndex = 1, mmpi_nprocs
        ! This is varLevIndex value of destination (can be different for source)
        varLevIndex2 = statevector_VarsLevs%allVarLevBeg(procIndex) + kCount - 1
        if ( varLevIndex2 > stateVector_VarsLevs%allVarLevEnd(procIndex) ) cycle proc_loop

        ! Figure out which variable/level of source
        varName = gsv_getVarNameFromVarLev(statevector,varLevIndex2)
        levIndex = gsv_getLevFromVarLev(statevector,varLevIndex2)
        allCols_ptr => col_getAllColumns(columnAnlInc,varName)

        do headerIndex = 1, numHeader
          cols_send(headerIndex,procIndex) = allCols_ptr(levIndex,headerIndex)
        end do
      end do proc_loop
      !$OMP END PARALLEL DO

      ! mpi communication: alltoall for one level/variable
      if(mmpi_nprocs > 1) then
        call mmpi_alltoall(cols_send, cols_recv)
      else
        cols_recv(:,1) = cols_send(:,1)
      end if

      if ( varLevIndex <= stateVector_VarsLevs%myVarLevEnd ) then
        varName = gsv_getVarNameFromVarLev(statevector,varLevIndex)

        if ( varName == 'UU' .or. varName == 'VV' ) then
          call gsv_getFieldUV(stateVector_VarsLevs, ptr3d_UV, varLevIndex)
        end if

        ! interpolate in time to the columns destined for all procs and one level/variable
        do procIndex = 1, mmpi_nprocs
          do stepIndex = 1, numStep
            !$OMP PARALLEL DO PRIVATE (headerIndex, headerUsedIndex, weight)
            do headerUsedIndex = 1, interpInfo_tlad%allNumHeaderUsed(stepIndex, procIndex)

              headerIndex = interpInfo_tlad%stepProcData(procIndex, stepIndex)%allHeaderIndex(headerUsedIndex)
              weight = oti_getTimeInterpWeightMpiGlobal(interpInfo_tlad%oti,  &
                                                        headerIndex,stepIndex,procIndex)

              cols_hint(headerUsedIndex,stepIndex,procIndex) =  &
                   weight * cols_recv(headerIndex,procIndex)

            end do
            !$OMP END PARALLEL DO
          end do
        end do

        call utl_tmg_start(42,'------s2c_AD_Hinterp')
        !$OMP PARALLEL DO PRIVATE (stepIndex, procIndex, yourNumHeader)
        step_loop: do stepIndex = 1, numStep
          if ( maxval(interpInfo_tlad%allNumHeaderUsed(stepIndex,:)) == 0 ) cycle step_loop

          ! interpolate to the columns destined for all procs for all steps and one lev/var
          do procIndex = 1, mmpi_nprocs
            yourNumHeader = interpInfo_tlad%allNumHeaderUsed(stepIndex,procIndex)
            if ( yourNumHeader > 0 ) then
              if ( varName == 'UU' ) then
                call myezuvint_ad( cols_hint(1:yourNumHeader,stepIndex,procIndex), 'UU',  &
                                   ptr4d(:,:,varLevIndex,stepIndex), ptr3d_UV(:,:,stepIndex),  &
                                   interpInfo_tlad, varLevIndex, stepIndex, procIndex )
              else if ( varName == 'VV' ) then
                call myezuvint_ad( cols_hint(1:yourNumHeader,stepIndex,procIndex), 'VV',  &
                                   ptr3d_UV(:,:,stepIndex), ptr4d(:,:,varLevIndex,stepIndex),  &
                                   interpInfo_tlad, varLevIndex, stepIndex, procIndex )
              else
                call myezsint_ad( cols_hint(1:yourNumHeader,stepIndex,procIndex), &
                                  ptr4d(:,:,varLevIndex,stepIndex), interpInfo_tlad, varLevIndex, &
                                  stepIndex, procIndex )
              end if
            end if
          end do

        end do step_loop
        !$OMP END PARALLEL DO
        call utl_tmg_stop(42)

      end if ! if varLevIndex <= myVarLevEnd

    end do k_loop

    deallocate(cols_hint)
    deallocate(cols_send)
    deallocate(cols_recv)

    call gsv_transposeTilesToVarsLevsAd( statevector_VarsLevs, statevector )

    call gsv_deallocate( statevector_VarsLevs )

  end subroutine ad2dFields

  !---------------------------------------------------------
  ! s2c_nl
  !---------------------------------------------------------
  subroutine s2c_nl(stateVector, obsSpaceData, column, hco_core, timeInterpType, &
                    numObsBatches_opt, dealloc_opt, moveObsAtPole_opt, &
                    beSilent_opt)
    !
    ! :Purpose: Non-linear version of the horizontal interpolation,
    !           used for a full field (usually the background state when computing
    !           the innovation vector).
    !
    implicit none

    ! Arguments:
    type(struct_gsv)       ,    intent(inout) :: stateVector
    type(struct_obs)       ,    intent(inout) :: obsSpaceData
    type(struct_columnData),    intent(inout) :: column
    type(struct_hco), pointer,  intent(in)    :: hco_core
    character(len=*)          , intent(in)    :: timeInterpType
    integer,          optional, intent(in)    :: numObsBatches_opt
    logical,          optional, intent(in)    :: dealloc_opt
    logical,          optional, intent(in)    :: moveObsAtPole_opt
    logical,          optional, intent(in)    :: beSilent_opt

    ! Locals:
    logical :: moveObsAtPole, rejectOutsideObs, beSilent

    call utl_tmg_start(30,'--StateToColumn')
    call utl_tmg_start(34,'----s2c_NL')

    ! Read the namelist
    call readNml()

    if (present(moveObsAtPole_opt)) then
      moveObsAtPole = moveObsAtPole_opt
    else
      moveObsAtPole = .false.
    end if

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .false.
    end if

    if (.not. interpInfo_nl%initialized .and. &
        .not. interpInfoTiles_nl%initialized) then

      call utl_tmg_start(31,'----s2c_Setups')
      ! Reject obs outside (LAM) domain and optionally move obs near
      ! numerical pole to first/last analysis grid latitude
      call latlonChecksAnlGrid(obsSpaceData, hco_core, moveObsAtPole)

      ! Do not reject obs for global domain
      rejectOutsideObs = .not. stateVector%hco%global
      write(*,*) 's2c_nl: rejectOutsideObs = ', rejectOutsideObs
      call utl_tmg_stop(31)

    end if

    if (.not. beSilent) then
      write(*,*) 's2c_nl: oceanMaskPresent = ', stateVector%oceanMask%maskPresent
    end if

    if (trim(mpiMode) == '2DFIELDS') then
      call nl2dFields(stateVector, obsSpaceData, column, timeInterpType, &
                      rejectOutsideObs, beSilent, numObsBatches_opt, dealloc_opt)
    else if (trim(mpiMode) == 'TILES') then
      call nlTiles(stateVector, obsSpaceData, column, timeInterpType, &
                   rejectOutsideObs, beSilent, dealloc_opt)
    else
      call utl_abort('s2c_nl: invalid value of mpiMode = '//trim(mpiMode))
    end if

    call utl_tmg_stop(34)
    call utl_tmg_stop(30)

  end subroutine s2c_nl

  !---------------------------------------------------------
  ! nl2dFields
  !---------------------------------------------------------
  subroutine nl2dFields(stateVector, obsSpaceData, column, timeInterpType, &
                        rejectOutsideObs, beSilent, numObsBatches_opt, dealloc_opt)
    !
    ! :Purpose: Non-linear version of the horizontal interpolation,
    !           used for a full field (usually the background state when computing
    !           the innovation vector).
    !
    implicit none

    ! Arguments:
    type(struct_gsv)       ,    intent(inout) :: stateVector
    type(struct_obs)       ,    intent(inout) :: obsSpaceData
    type(struct_columnData),    intent(inout) :: column
    character(len=*)          , intent(in)    :: timeInterpType
    logical,                    intent(inout) :: rejectOutsideObs
    logical,                    intent(in)    :: beSilent
    integer,          optional, intent(in)    :: numObsBatches_opt
    logical,          optional, intent(in)    :: dealloc_opt

    ! Locals:
    type(struct_gsv), save :: stateVector_VarsLevs
    integer :: varLevIndex, varLevIndex2, kCount, stepIndex, numStep, myVarLevEndExtended, levIndex
    integer :: headerIndex, headerIndex2, numHeader, numHeaderMax, yourNumHeader
    integer :: headerIndexBeg, headerIndexEnd, obsBatchIndex, numObsBatches
    integer :: procIndex, headerUsedIndex, allHeaderIndexBeg(mmpi_nprocs)
    integer :: varLevIndexHeightSfc, varNameIndex, allNumHeader(mmpi_nprocs)
    real(8) :: weight
    character(len=4)     :: varName
    real(8), pointer     :: column_ptr(:), ptr2d_r8(:,:), allCols_ptr(:,:)
    real(4), pointer     :: ptr4d_r4(:,:,:,:), ptr3d_UV_r4(:,:,:)
    real(8), allocatable :: cols_hint(:,:,:)
    real(8), allocatable :: cols_send(:,:)
    real(8), allocatable :: cols_recv(:,:)
    real(8), allocatable :: cols_send_1proc(:)
    integer              :: displs(mmpi_nprocs), nsizes(mmpi_nprocs)
    integer              :: senddispls(mmpi_nprocs), sendsizes(mmpi_nprocs)
    integer              :: recvdispls(mmpi_nprocs), recvsizes(mmpi_nprocs)
    logical              :: dealloc
    logical, save        :: firstCall = .true.
    character(len=4), pointer :: varNames(:)

    if (.not. beSilent) then
      write(*,*) 'nl2dFields: STARTING'
      call msg_memUsage('s2c_nl')
    end if

    if (.not. gsv_isAllocated(stateVector)) then
      call utl_abort('nl2dFields: stateVector must be allocated')
    end if

    if (present(dealloc_opt)) then
      dealloc = dealloc_opt
    else
      dealloc = .true.
    end if

    ! determine number of obs batches (to reduce memory usage)
    if (present(numObsBatches_opt)) then
      numObsBatches = numObsBatches_opt
    else
      numObsBatches = 1
    end if

    if (.not. dealloc) then
      if (numObsBatches /= 1) then
        write(*,*) 's2c_nl: WARNING! numObsBatches=', numObsBatches, ' will be set to 1.'
      end if

      numObsBatches = 1 ! multiple batches only possible if dealloc=.true.
    end if

    if (interpInfo_nl%initialized) then
      if (.not. hco_equal(interpInfo_nl%hco,stateVector%hco) .or. numObsBatches > 1) then
        write(*,*) 'nl2dFields: WARNING! Current hco grid parameters differ from allocated interpInfo!'
        write(*,*) 'nl2dFields: InterpInfo will be deallocated.'
        call s2c_deallocInterpInfo(inputStateVectorType = 'nl')
      end if
    end if

    if (stateVector%mpi_distribution /= 'Tiles') then
      call utl_abort('nl2dFields: stateVector must by Tiles distributed')
    end if

    ! check the column and statevector have same nk/varNameList
    call checkColumnStatevectorMatch(column,statevector)

    ! calculate delP_T/delP_M and del Z_T/Z_M on the grid
    call gvt_transform(statevector, 'ZandP_nl')

    if (dealloc .or. firstCall) then
      nullify(varNames)
      call gsv_varNamesList(varNames, statevector)
      call gsv_allocate(statevector_VarsLevs, stateVector%numstep, &
                        stateVector%hco, stateVector%vco, mpi_local_opt = .true., &
                        mpi_distribution_opt = 'VarsLevs', dataKind_opt = 4, &
                        allocHeightSfc_opt = .true., varNames_opt = varNames )
      deallocate(varNames)
    else
      if (mmpi_myid == 0 .and. .not. beSilent) write(*,*) 'nl2dFields: avoid re-allocating statevector_VarsLevs'
      call gsv_zero(statevector_VarsLevs)
    end if

    call gsv_transposeTilesToVarsLevs(stateVector, stateVector_VarsLevs, &
                                      beSilent_opt=beSilent)

    numStep = stateVector_VarsLevs%numStep

    ! set contents of column to zero
    call col_zero(column)

    OBSBATCH: do obsBatchIndex = 1, numObsBatches
      headerIndexBeg = 1 + (obsBatchIndex - 1) * (obs_numheader(obsSpaceData) / numObsBatches)
      if (obsBatchIndex == numObsBatches) then
        headerIndexEnd = obs_numheader(obsSpaceData)
      else
        headerIndexEnd = headerIndexBeg + (obs_numheader(obsSpaceData) / numObsBatches) - 1
      end if
      numHeader = headerIndexEnd - headerIndexBeg + 1
      call mmpi_allReduce(numHeader, numHeaderMax, mmpi_max)

      call mmpi_allGather(numHeader,      allNumHeader)
      call mmpi_allGather(headerIndexBeg, allHeaderIndexBeg)
      if ( .not. beSilent ) then
        write(*,*) 'nl2dFields: headerIndexBeg/End, numHeader, numHeaderMax = ',  &
             headerIndexBeg, headerIndexEnd, numHeader, numHeaderMax
        if (mmpi_myid == 0) then
           write(*,*) 'nl2dFields: min/max of allNumHeader = ', minval(allNumHeader), maxval(allNumHeader)
        end if
      end if

      if (.not. interpInfo_nl%initialized) then
        call utl_tmg_stop(34)
        call utl_tmg_start(31,'----s2c_Setups')

        ! compute and collect all obs grids onto all mpi tasks
        call setupInterpInfo2dFields(interpInfo_nl, obsSpaceData, stateVector_VarsLevs, &
                                     headerIndexBeg, headerIndexEnd, &
                                     timeInterpType, rejectOutsideObs, &
                                     inputStateVectorType = 'nl', &
                                     lastCall_opt = (obsBatchIndex == numObsBatches))
        if (mmpi_myid == 0 .and. verbose) then
          do stepIndex = 1, numStep
            write(*,*) 'nl2dFields: stepIndex, allNumHeaderUsed = ',  &
                       stepIndex, interpInfo_nl%allNumHeaderUsed(stepIndex,:)
          end do
        end if
        call utl_tmg_stop(31)
        call utl_tmg_start(34,'----s2c_NL')
      end if

      ! arrays for interpolated column for 1 level/variable and each time step
      allocate(cols_hint(maxval(interpInfo_nl%allNumHeaderUsed), numStep, mmpi_nprocs))
      cols_hint(:,:,:) = 0.0d0

      ! arrays for sending/receiving time interpolated column for 1 level/variable
      allocate(cols_send(numHeaderMax, mmpi_nprocs))
      cols_send(:,:) = 0.0d0

      allocate(cols_recv(numHeaderMax, mmpi_nprocs))
      cols_recv(:,:) = 0.0d0

      allocate(cols_send_1proc(numHeaderMax))

      call gsv_getField(stateVector_VarsLevs, ptr4d_r4)

      myVarLevEndExtended = stateVector_VarsLevs%myVarLevBeg + maxval(stateVector_VarsLevs%allVarLevCount(:)) - 1

      kCount = 0
      k_loop: do varLevIndex = stateVector_VarsLevs%myVarLevBeg, myVarLevEndExtended
        kCount = kCount + 1

        if ( varLevIndex <= stateVector_VarsLevs%myVarLevEnd ) then
          varName = gsv_getVarNameFromVarLev(stateVector_VarsLevs, varLevIndex)

          call utl_tmg_start(35,'------s2c_NL_Hinterp')
          if ( varName == 'UU' .or. varName == 'VV' ) then
            call gsv_getFieldUV(stateVector_VarsLevs,ptr3d_UV_r4,varLevIndex)
          end if

          step_loop: do stepIndex = 1, numStep
            if (maxval(interpInfo_nl%allNumHeaderUsed(stepIndex,:)) == 0) cycle step_loop

            ! interpolate to the columns destined for all procs for all steps and one lev/var
            !$OMP PARALLEL DO PRIVATE (procIndex, yourNumHeader)
            do procIndex = 1, mmpi_nprocs
              yourNumHeader = interpInfo_nl%allNumHeaderUsed(stepIndex,procIndex)
              if (yourNumHeader > 0) then
                if (varName == 'UU') then
                  call myezuvint_nl(cols_hint(1:yourNumHeader, stepIndex, procIndex), 'UU',  &
                                    ptr4d_r4(:,:, varLevIndex, stepIndex), ptr3d_UV_r4(:,:,stepIndex), &
                                    interpInfo_nl, varLevIndex, stepIndex, procIndex)
                else if ( varName == 'VV' ) then
                  call myezuvint_nl(cols_hint(1:yourNumHeader, stepIndex, procIndex), 'VV',  &
                                    ptr3d_UV_r4(:,:,stepIndex), ptr4d_r4(:,:, varLevIndex, stepIndex), &
                                    interpInfo_nl, varLevIndex, stepIndex, procIndex)
                else
                  call myezsint_nl(cols_hint(1:yourNumHeader, stepIndex, procIndex), &
                                   ptr4d_r4(:,:, varLevIndex, stepIndex),  &
                                   interpInfo_nl, varLevIndex, stepIndex, procIndex)
                end if
              end if
            end do
            !$OMP END PARALLEL DO

          end do step_loop
          call utl_tmg_stop(35)

          ! interpolate in time to the columns destined for all procs and one level/variable
          do procIndex = 1, mmpi_nprocs
            cols_send_1proc(:) = 0.0d0
            do stepIndex = 1, numStep
              !$OMP PARALLEL DO PRIVATE (headerIndex, headerIndex2, headerUsedIndex, weight)
              do headerUsedIndex = 1, interpInfo_nl%allNumHeaderUsed(stepIndex, procIndex)
                headerIndex = interpInfo_nl%stepProcData(procIndex, stepIndex)%allHeaderIndex(headerUsedIndex)
                headerIndex2 = headerIndex - allHeaderIndexBeg(procIndex) + 1
                weight = oti_getTimeInterpWeightMpiGlobal(interpInfo_nl%oti,  &
                                                          headerIndex2, stepIndex, procIndex)
                cols_send_1proc(headerIndex2) = cols_send_1proc(headerIndex2) + &
                                                weight * cols_hint(headerUsedIndex, stepIndex, procIndex)
              end do
              !$OMP END PARALLEL DO
            end do
            cols_send(:,procIndex) = cols_send_1proc(:)
          end do

        else

          ! this value of k does not exist on this mpi task
          cols_send(:,:) = 0.0d0

        end if ! if varLevIndex <= myVarLevEnd

        ! mpi communication: alltoallv for one level/variable

        ! only receive the data from tasks with data, same amount from all of those
        recvsizes(:) = 0
        do procIndex = 1, mmpi_nprocs
          varLevIndex2 = stateVector_VarsLevs%allVarLevBeg(procIndex) + kCount - 1
          if (varLevIndex2 > stateVector_VarsLevs%allVarLevEnd(procIndex)) cycle
          recvsizes(procIndex) = numHeader
        end do
        recvdispls(1) = 0
        do procIndex = 2, mmpi_nprocs
          recvdispls(procIndex) = recvdispls(procIndex-1) + numHeaderMax
        end do

        ! tasks send only from those with data
        if (varLevIndex <= stateVector_VarsLevs%myVarLevEnd) then
          do procIndex = 1, mmpi_nprocs
            sendsizes(procIndex) = allNumHeader(procIndex)
          end do
        else
          sendsizes(:) = 0
        end if
        senddispls(1) = 0
        do procIndex = 2, mmpi_nprocs
          senddispls(procIndex) = senddispls(procIndex-1) + numHeaderMax
        end do

        if(mmpi_nprocs > 1) then
          call mmpi_alltoallv(cols_send, sendsizes, senddispls, &
                              cols_recv, recvsizes, recvdispls)
        else
          cols_recv(:,1) = cols_send(:,1)
        end if

        ! reorganize ensemble of distributed columns
        !$OMP PARALLEL DO PRIVATE (procIndex, varLevIndex2, headerIndex, headerIndex2, varName, &
        !$OMP levIndex, allCols_ptr)
        proc_loop: do procIndex = 1, mmpi_nprocs

          varLevIndex2 = stateVector_VarsLevs%allVarLevBeg(procIndex) + kCount - 1
          if ( varLevIndex2 > stateVector_VarsLevs%allVarLevEnd(procIndex) ) cycle proc_loop

          varName = gsv_getVarNameFromVarLev(stateVector_VarsLevs, varLevIndex2)
          levIndex = gsv_getLevFromVarLev(statevector, varLevIndex2)
          allCols_ptr => col_getAllColumns(column, varName)

          do headerIndex = 1, numHeader
            headerIndex2 = headerIndex + headerIndexBeg - 1
            allCols_ptr(levIndex, headerIndex2) = cols_recv(headerIndex, procIndex)
          end do

        end do proc_loop
        !$OMP END PARALLEL DO

      end do k_loop

      ! impose a lower limit on HU
      if(col_varExist(column,'HU')) then
        do headerIndex = headerIndexBeg, headerIndexEnd
          column_ptr => col_getColumn(column, headerIndex,'HU')
          column_ptr(:) = max(column_ptr(:), col_rhumin)
        end do
      end if

      ! impose a lower/upper limits on ALL cloud variables
      do varNameIndex = 1, vnl_numvarmaxCloud
        if(col_varExist(column, vnl_varNameListCloud(varNameIndex))) then
          do headerIndex = headerIndexBeg, headerIndexEnd
            column_ptr => col_getColumn(column,headerIndex, vnl_varNameListCloud(varNameIndex))
            column_ptr(:) = max(column_ptr(:), qlim_getMinValueCloud(vnl_varNameListCloud(varNameIndex)))
            column_ptr(:) = min(column_ptr(:), qlim_getMaxValueCloud(vnl_varNameListCloud(varNameIndex)))
          end do
        end if
      end do

      ! Interpolate surface height separately, only exists on mpi task 0
      HeightSfcPresent: if ( stateVector_VarsLevs%HeightSfcPresent ) then

        if (mmpi_myid == 0) then
          !varName = 'GZ'
          varLevIndexHeightSfc = 0

          step_loop_height: do stepIndex = 1, numStep

            if (maxval(interpInfo_nl%allNumHeaderUsed(stepIndex,:)) == 0) cycle step_loop_height

            ! interpolate to the columns destined for all procs for all steps and one lev/var
            !$OMP PARALLEL DO PRIVATE (procIndex, yourNumHeader, ptr2d_r8)
            do procIndex = 1, mmpi_nprocs
              yourNumHeader = interpInfo_nl%allNumHeaderUsed(stepIndex, procIndex)
              if ( yourNumHeader > 0 ) then
                ptr2d_r8 => gsv_getHeightSfc(stateVector_VarsLevs)
                call myezsint_r8_nl( cols_hint(1:yourNumHeader, stepIndex, procIndex), &
                                     ptr2d_r8(:,:), interpInfo_nl, varLevIndexHeightSfc, stepIndex, procIndex )
              end if
            end do
            !$OMP END PARALLEL DO

          end do step_loop_height

          ! interpolate in time to the columns destined for all procs and one level/variable
          do procIndex = 1, mmpi_nprocs
            cols_send(:,procIndex) = 0.0d0
            do stepIndex = 1, numStep
              !$OMP PARALLEL DO PRIVATE (headerIndex, headerIndex2, headerUsedIndex)
              do headerUsedIndex = 1, interpInfo_nl%allNumHeaderUsed(stepIndex, procIndex)
                headerIndex = interpInfo_nl%stepProcData(procIndex, stepIndex)%allHeaderIndex(headerUsedIndex)
                ! just copy, since surface height same for all time steps
                headerIndex2 = headerIndex - allHeaderIndexBeg(procIndex) + 1
                cols_send(headerIndex2,procIndex) = cols_hint(headerUsedIndex, stepIndex, procIndex)
              end do
              !$OMP END PARALLEL DO
            end do
          end do

        end if

        ! mpi communication: scatter data from task 0
        if(mmpi_nprocs > 1) then
          do procIndex = 1, mmpi_nprocs
            displs(procIndex) = (procIndex - 1) * numHeaderMax
            nsizes(procIndex) = allNumHeader(procIndex)
          end do
          call mmpi_scatterv(cols_send, cols_recv, nsizes, displs, numHeader)
        else
          cols_recv(:,1) = cols_send(:,1)
        end if

        do headerIndex = headerIndexBeg, headerIndexEnd
          headerIndex2 = headerIndex - headerIndexBeg + 1
          call col_setHeightSfc(column, headerIndex, cols_recv(headerIndex2,1))
        end do

      end if HeightSfcPresent

      ! Interpolate surface height LS separately, only exists on mpi task 0
      HeightSfcLsPresent: if ( stateVector_VarsLevs%HeightSfcLsPresent ) then

        if (mmpi_myid == 0) then
          !varName = 'MELS'
          varLevIndexHeightSfc = 0
          step_loop_heightLs: do stepIndex = 1, numStep

            if (maxval(interpInfo_nl%allNumHeaderUsed(stepIndex,:)) == 0) cycle step_loop_heightLs

            ! interpolate to the columns destined for all procs for all steps and one lev/var
            !$OMP PARALLEL DO PRIVATE (procIndex, yourNumHeader, ptr2d_r8)
            do procIndex = 1, mmpi_nprocs
              yourNumHeader = interpInfo_nl%allNumHeaderUsed(stepIndex, procIndex)
              if ( yourNumHeader > 0 ) then
                ptr2d_r8 => gsv_getHeightSfcLS(stateVector_VarsLevs)
                call myezsint_r8_nl( cols_hint(1:yourNumHeader, stepIndex, procIndex), &
                                     ptr2d_r8(:,:), interpInfo_nl, varLevIndexHeightSfc, stepIndex, procIndex )
              end if
            end do
            !$OMP END PARALLEL DO

          end do step_loop_heightLs

          ! interpolate in time to the columns destined for all procs and one level/variable
          do procIndex = 1, mmpi_nprocs
            cols_send(:,procIndex) = 0.0d0
            do stepIndex = 1, numStep
              !$OMP PARALLEL DO PRIVATE (headerIndex, headerIndex2, headerUsedIndex)
              do headerUsedIndex = 1, interpInfo_nl%allNumHeaderUsed(stepIndex, procIndex)
                headerIndex = interpInfo_nl%stepProcData(procIndex, stepIndex)%allHeaderIndex(headerUsedIndex)
                ! just copy, since surface height same for all time steps
                headerIndex2 = headerIndex - allHeaderIndexBeg(procIndex) + 1
                cols_send(headerIndex2,procIndex) = cols_hint(headerUsedIndex, stepIndex, procIndex)
              end do
              !$OMP END PARALLEL DO
            end do
          end do

        end if

        ! mpi communication: scatter data from task 0
        if(mmpi_nprocs > 1) then
          do procIndex = 1, mmpi_nprocs
            displs(procIndex) = (procIndex - 1) * numHeaderMax
            nsizes(procIndex) = allNumHeader(procIndex)
          end do
          call mmpi_scatterv(cols_send, cols_recv, nsizes, displs, numHeader)
        else
          cols_recv(:,1) = cols_send(:,1)
        end if

        do headerIndex = headerIndexBeg, headerIndexEnd
          headerIndex2 = headerIndex - headerIndexBeg + 1
          call col_setHeightSfcLs(column, headerIndex, cols_recv(headerIndex2,1))
        end do

      end if HeightSfcLsPresent

      deallocate(cols_hint)
      deallocate(cols_send)
      deallocate(cols_recv)
      deallocate(cols_send_1proc)

      if ( dealloc ) call s2c_deallocInterpInfo( inputStateVectorType='nl' )

    end do OBSBATCH

    if (dealloc) call gsv_deallocate( statevector_VarsLevs )

    if (slantPath_TO_nl) call pressureProfileMonotonicityCheck(obsSpaceData, column)

    firstCall = .false.

    if ( .not. beSilent ) then
      write(*,*) 'nl2dFields: FINISHED'
      call msg_memUsage('s2c_nl')
    end if

  end subroutine nl2dFields

  ! -------------------------------------------------
  ! myezsint_nl: Scalar field horizontal interpolation
  ! -------------------------------------------------
  subroutine myezsint_nl( column_out, field_in, interpInfo, varLevIndex, stepIndex, procIndex )
    !
    ! :Purpose: Scalar horizontal interpolation, replaces the
    !           ezsint routine from rmnlib.
    !
    implicit none

    ! Arguments:
    real(8)                , intent(out) :: column_out(:)
    real(4)                , intent(in)  :: field_in(:,:)
    type(struct_interpInfo), intent(in)  :: interpInfo
    integer                , intent(in)  :: stepIndex
    integer                , intent(in)  :: procIndex
    integer                , intent(in)  :: varLevIndex

    ! Locals:
    integer :: lonIndex, latIndex, gridptIndex, headerIndex, subGridIndex, numColumn
    real(8) :: interpValue, weight

    numColumn = size( column_out )

    do headerIndex = 1, numColumn

      ! Interpolate the model state to the obs point
      interpValue = 0.0d0

      do subGridIndex = 1, interpInfo%hco%numSubGrid

        do gridptIndex =  &
             interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, varLevIndex), &
             interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex, headerIndex, varLevIndex)

          lonIndex = interpInfo%lonIndexDepot(gridptIndex)
          latIndex = interpInfo%latIndexDepot(gridptIndex)
          weight = interpInfo%interpWeightDepot(gridptIndex)

          interpValue = interpValue + weight * real(field_in(lonIndex, latIndex),8)

        end do

      end do
      column_out(headerIndex) = interpValue

    end do

  end subroutine myezsint_nl

  ! -------------------------------------------------
  ! myezsint_r8_nl: Scalar field horizontal interpolation
  ! -------------------------------------------------
  subroutine myezsint_r8_nl( column_out, field_in, interpInfo, varLevIndex, stepIndex, procIndex )
    !
    ! :Purpose: Scalar horizontal interpolation, replaces the
    !           ezsint routine from rmnlib.
    !
    implicit none

    ! Arguments:
    real(8)                , intent(out) :: column_out(:)
    real(8)                , intent(in)  :: field_in(:,:)
    type(struct_interpInfo), intent(in)  :: interpInfo
    integer                , intent(in)  :: stepIndex
    integer                , intent(in)  :: procIndex
    integer                , intent(in)  :: varLevIndex

    ! Locals:
    integer :: lonIndex, latIndex, gridptIndex, headerIndex, subGridIndex, numColumn
    real(8) :: interpValue, weight

    numColumn = size( column_out )

    do headerIndex = 1, numColumn

      ! Interpolate the model state to the obs point
      interpValue = 0.0d0

      do subGridIndex = 1, interpInfo%hco%numSubGrid

        do gridptIndex =  &
             interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, varLevIndex), &
             interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex, headerIndex, varLevIndex)

          lonIndex = interpInfo%lonIndexDepot(gridptIndex)
          latIndex = interpInfo%latIndexDepot(gridptIndex)
          weight = interpInfo%interpWeightDepot(gridptIndex)

          interpValue = interpValue + weight * field_in(lonIndex, latIndex)

        end do

      end do
      column_out(headerIndex) = interpValue

    end do

  end subroutine myezsint_r8_nl

  ! -------------------------------------------------
  ! myezsint_tl: Scalar field horizontal interpolation
  ! -------------------------------------------------
  subroutine myezsint_tl( column_out, field_in, interpInfo, varLevIndex, stepIndex, procIndex )
    !
    ! :Purpose: Scalar horizontal interpolation, replaces the
    !           ezsint routine from rmnlib.
    !
    implicit none

    ! Arguments:
    real(8)                , intent(out) :: column_out(:)
    real(pre_incrReal)     , intent(in)  :: field_in(:,:)
    type(struct_interpInfo), intent(in)  :: interpInfo
    integer                , intent(in)  :: stepIndex
    integer                , intent(in)  :: procIndex
    integer                , intent(in)  :: varLevIndex

    ! Locals:
    integer :: lonIndex, latIndex, gridptIndex, headerIndex, subGridIndex, numColumn
    real(8) :: interpValue, weight

    numColumn = size( column_out )

    do headerIndex = 1, numColumn

      ! Interpolate the model state to the obs point
      interpValue = 0.0d0

      do subGridIndex = 1, interpInfo%hco%numSubGrid

        do gridptIndex =  &
             interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, varLevIndex), &
             interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex, headerIndex, varLevIndex)

          lonIndex = interpInfo%lonIndexDepot(gridptIndex)
          latIndex = interpInfo%latIndexDepot(gridptIndex)
          weight = interpInfo%interpWeightDepot(gridptIndex)

          interpValue = interpValue + weight * field_in(lonIndex, latIndex)

        end do

      end do
      column_out(headerIndex) = interpValue

    end do

  end subroutine myezsint_tl

  ! -------------------------------------------------------------
  ! myezsint_ad: Adjoint of scalar field horizontal interpolation
  ! -------------------------------------------------------------
  subroutine myezsint_ad( column_in, field_out, interpInfo, varLevIndex, stepIndex, procIndex )
    !
    ! :Purpose: Adjoint of the scalar horizontal interpolation.
    !
    implicit none

    ! Arguments:
    real(8)                , intent(in)    :: column_in(:)
    real(pre_incrReal)     , intent(inout) :: field_out(:,:)
    type(struct_interpInfo), intent(in)    :: interpInfo
    integer                , intent(in)    :: stepIndex
    integer                , intent(in)    :: procIndex
    integer                , intent(in)    :: varLevIndex

    ! Locals:
    integer :: lonIndex, latIndex, gridptIndex, headerIndex, subGridIndex, numColumn
    real(8) :: weight

    numColumn = size( column_in )

    do headerIndex = 1, numColumn

      ! Interpolate the model state to the obs point

      do subGridIndex = 1, interpInfo%hco%numSubGrid

        do gridptIndex =  &
             interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, varLevIndex), &
             interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex, headerIndex, varLevIndex)

          lonIndex = interpInfo%lonIndexDepot(gridptIndex)
          latIndex = interpInfo%latIndexDepot(gridptIndex)
          weight = interpInfo%interpWeightDepot(gridptIndex)

          field_out(lonIndex, latIndex) = field_out(lonIndex, latIndex) +  &
                                          weight * column_in(headerIndex)

        end do

      end do

    end do

  end subroutine myezsint_ad

  ! -------------------------------------------------------------
  ! myezuvint_nl: Vector field horizontal interpolation
  ! -------------------------------------------------------------
  subroutine myezuvint_nl( column_out, varName, fieldUU_in, fieldVV_in,  &
                           interpInfo, varLevIndex, stepIndex, procIndex )
    !
    ! :Purpose: Vector horizontal interpolation, replaces the
    !           ezuvint routine from rmnlib.
    !
    implicit none

    ! Arguments:
    real(8)                , intent(out) :: column_out(:)
    character(len=*)       , intent(in)  :: varName
    real(4)                , intent(in)  :: fieldUU_in(:,:)
    real(4)                , intent(in)  :: fieldVV_in(:,:)
    type(struct_interpInfo), intent(in)  :: interpInfo
    integer                , intent(in)  :: stepIndex
    integer                , intent(in)  :: procIndex
    integer                , intent(in)  :: varLevIndex

    ! Locals:
    integer :: lonIndex, latIndex, indexBeg, indexEnd, gridptIndex, headerIndex
    integer :: numColumn, subGridIndex
    real(8) :: interpUU(interpInfo%hco%numSubGrid), interpVV(interpInfo%hco%numSubGrid)
    real(8) :: lat, lon, latRot, lonRot, weight
    logical :: doUU, doVV

    numColumn = size( column_out )

    doUU = (trim(varName) == 'UU' .or. interpInfo%hco%rotated)
    doVV = (trim(varName) == 'VV' .or. interpInfo%hco%rotated)

    header_loop: do headerIndex = 1, numColumn

      interpUU(:) = 0.0d0
      interpVV(:) = 0.0d0

      subGrid_loop: do subGridIndex = 1, interpInfo%hco%numSubGrid

        indexBeg = interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, varLevIndex)
        indexEnd = interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex, headerIndex, varLevIndex)

        if ( indexEnd < IndexBeg ) cycle subGrid_loop

        ! Interpolate the model UU to the obs point
        do gridptIndex = indexBeg, indexEnd

          lonIndex = interpInfo%lonIndexDepot(gridptIndex)
          latIndex = interpInfo%latIndexDepot(gridptIndex)
          weight = interpInfo%interpWeightDepot(gridptIndex)

          if ( doUU ) interpUU(subGridIndex) = interpUU(subGridIndex) +  &
                      weight * real(fieldUU_in(lonIndex, latIndex),8)
          if ( doVV ) interpVV(subGridIndex) = interpVV(subGridIndex) +  &
                      weight * real(fieldVV_in(lonIndex, latIndex),8)

        end do
        ! now rotate the wind vector
        if ( interpInfo%hco%rotated ) then
          lat = interpInfo%stepProcData(procIndex, stepIndex)%allLat(headerIndex, varLevIndex)
          lon = interpInfo%stepProcData(procIndex, stepIndex)%allLon(headerIndex, varLevIndex)
          latRot = interpInfo%stepProcData(procIndex,stepIndex)%allLatRot(subGridIndex, headerIndex, varLevIndex)
          lonRot = interpInfo%stepProcData(procIndex,stepIndex)%allLonRot(subGridIndex, headerIndex, varLevIndex)

          call uvr_rotateWind_nl( interpInfo%uvr,            & ! IN
                                  subGridIndex,              & ! IN
                                  interpUU(subGridIndex),    & ! INOUT
                                  interpVV(subGridIndex),    & ! INOUT
                                  lat, lon, latRot, lonRot,  & ! IN
                                  'ToMetWind' )                ! IN
        end if

      end do subGrid_loop

      ! return only the desired component
      if ( trim(varName) == 'UU' ) then
        column_out(headerIndex) = sum(interpUU(:))
      else
        column_out(headerIndex) = sum(interpVV(:))
      end if

    end do header_loop

  end subroutine myezuvint_nl

  ! -------------------------------------------------------------
  ! myezuvint_tl: Vector field horizontal interpolation
  ! -------------------------------------------------------------
  subroutine myezuvint_tl( column_out, varName, fieldUU_in, fieldVV_in,  &
                           interpInfo, varLevIndex, stepIndex, procIndex )
    ! :Purpose: Vector horizontal interpolation, replaces the
    !           ezuvint routine from rmnlib.
    !
    implicit none

    ! Arguments:
    real(8)                , intent(out) :: column_out(:)
    character(len=*)       , intent(in)  :: varName
    real(pre_incrReal)     , intent(in)  :: fieldUU_in(:,:)
    real(pre_incrReal)     , intent(in)  :: fieldVV_in(:,:)
    type(struct_interpInfo), intent(in)  :: interpInfo
    integer                , intent(in)  :: stepIndex
    integer                , intent(in)  :: procIndex
    integer                , intent(in)  :: varLevIndex

    ! Locals:
    integer :: lonIndex, latIndex, indexBeg, indexEnd, gridptIndex, headerIndex
    integer :: numColumn, subGridIndex
    real(8) :: interpUU(interpInfo%hco%numSubGrid), interpVV(interpInfo%hco%numSubGrid)
    real(8) :: lat, lon, latRot, lonRot, weight
    logical :: doUU, doVV

    numColumn = size( column_out )

    doUU = (trim(varName) == 'UU' .or. interpInfo%hco%rotated)
    doVV = (trim(varName) == 'VV' .or. interpInfo%hco%rotated)

    header_loop: do headerIndex = 1, numColumn

      interpUU(:) = 0.0d0
      interpVV(:) = 0.0d0

      subGrid_loop: do subGridIndex = 1, interpInfo%hco%numSubGrid

        indexBeg = interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, varLevIndex)
        indexEnd = interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex, headerIndex, varLevIndex)

        if ( indexEnd < IndexBeg ) cycle subGrid_loop

        ! Interpolate the model UU to the obs point
        do gridptIndex = indexBeg, indexEnd

          lonIndex = interpInfo%lonIndexDepot(gridptIndex)
          latIndex = interpInfo%latIndexDepot(gridptIndex)
          weight = interpInfo%interpWeightDepot(gridptIndex)

          if ( doUU ) interpUU(subGridIndex) = interpUU(subGridIndex) +  &
                      weight * fieldUU_in(lonIndex, latIndex)
          if ( doVV ) interpVV(subGridIndex) = interpVV(subGridIndex) +  &
                      weight * fieldVV_in(lonIndex, latIndex)

        end do
        ! now rotate the wind vector
        if ( interpInfo%hco%rotated ) then
          lat = interpInfo%stepProcData(procIndex, stepIndex)%allLat(headerIndex, varLevIndex)
          lon = interpInfo%stepProcData(procIndex, stepIndex)%allLon(headerIndex, varLevIndex)
          latRot = interpInfo%stepProcData(procIndex,stepIndex)%allLatRot(subGridIndex, headerIndex, varLevIndex)
          lonRot = interpInfo%stepProcData(procIndex,stepIndex)%allLonRot(subGridIndex, headerIndex, varLevIndex)

          call uvr_rotateWind_tl( interpInfo%uvr,            & ! IN
                                  subGridIndex,              & ! IN
                                  interpUU(subGridIndex),    & ! INOUT
                                  interpVV(subGridIndex),    & ! INOUT
                                  lat, lon, latRot, lonRot,  & ! IN
                                  'ToMetWind' )                ! IN
        end if

      end do subGrid_loop

      ! return only the desired component
      if ( trim(varName) == 'UU' ) then
        column_out(headerIndex) = sum(interpUU(:))
      else
        column_out(headerIndex) = sum(interpVV(:))
      end if

    end do header_loop

  end subroutine myezuvint_tl

  ! -------------------------------------------------------------
  ! myezuvint_ad: Adjoint of vector field horizontal interpolation
  ! -------------------------------------------------------------
  subroutine myezuvint_ad( column_in, varName, fieldUU_out, fieldVV_out, &
                           interpInfo, varLevIndex, stepIndex, procIndex )
    !
    ! :Purpose: Adjoint of the vector horizontal interpolation.
    !
    implicit none

    ! Arguments:
    real(8)                , intent(in)    :: column_in(:)
    character(len=*)       , intent(in)    :: varName
    real(pre_incrReal)     , intent(inout) :: fieldUU_out(:,:)
    real(pre_incrReal)     , intent(inout) :: fieldVV_out(:,:)
    type(struct_interpInfo), intent(in)    :: interpInfo
    integer                , intent(in)    :: stepIndex
    integer                , intent(in)    :: procIndex
    integer                , intent(in)    :: varLevIndex

    ! Locals:
    integer :: lonIndex, latIndex, indexBeg, indexEnd, gridptIndex, headerIndex
    integer :: numColumn, subGridIndex
    real(8) :: interpUU(interpInfo%hco%numSubGrid), interpVV(interpInfo%hco%numSubGrid)
    real(8) :: lat, lon, latRot, lonRot, weight
    logical :: doUU, doVV

    numColumn = size( column_in )

    doUU = (trim(varName) == 'UU' .or. interpInfo%hco%rotated)
    doVV = (trim(varName) == 'VV' .or. interpInfo%hco%rotated)

    header_loop: do headerIndex = 1, numColumn

      if ( trim(varName) == 'UU' ) then
        interpUU(:) = column_in(headerIndex)
        interpVV(:) = 0.0d0
      else
        interpUU(:) = 0.0d0
        interpVV(:) = column_in(headerIndex)
      end if

      subGrid_loop: do subGridIndex = 1, interpInfo%hco%numSubGrid

        indexBeg = interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, varLevIndex)
        indexEnd = interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex, headerIndex, varLevIndex)

        if ( indexEnd < IndexBeg ) cycle subGrid_loop

        ! now rotate the wind vector and return the desired component
        if ( interpInfo%hco%rotated ) then
          lat = interpInfo%stepProcData(procIndex, stepIndex)%allLat(headerIndex, varLevIndex)
          lon = interpInfo%stepProcData(procIndex, stepIndex)%allLon(headerIndex, varLevIndex)
          latRot = interpInfo%stepProcData(procIndex,stepIndex)%allLatRot(subGridIndex, headerIndex, varLevIndex)
          lonRot = interpInfo%stepProcData(procIndex,stepIndex)%allLonRot(subGridIndex, headerIndex, varLevIndex)

          call uvr_rotateWind_ad( interpInfo%uvr,           & ! IN
                                  subGridIndex,             & ! IN
                                  interpUU(subGridIndex),   & ! INOUT
                                  interpVV(subGridIndex),   & ! INOUT
                                  lat, lon, latRot, lonRot, & ! IN
                                  'ToMetWind' )               ! IN
        end if

        ! Interpolate the model VV to the obs point
        do gridptIndex = indexBeg, indexEnd

          lonIndex = interpInfo%lonIndexDepot(gridptIndex)
          latIndex = interpInfo%latIndexDepot(gridptIndex)
          weight = interpInfo%interpWeightDepot(gridptIndex)

          if ( doUU ) fieldUU_out(lonIndex, latIndex) =  &
                      fieldUU_out(lonIndex, latIndex) + weight * interpUU(subGridIndex)
          if ( doVV ) fieldVV_out(lonIndex, latIndex) =  &
                      fieldVV_out(lonIndex, latIndex) + weight * interpVV(subGridIndex)

        end do

      end do subGrid_loop

    end do header_loop

  end subroutine myezuvint_ad

  !---------------------------------------------------------
  ! s2c_bgcheck_bilin
  !---------------------------------------------------------
  subroutine s2c_bgcheck_bilin(column,statevector,obsSpaceData)
    !
    ! :Purpose: Special version of s2c_tl used for background check. This should
    !           be replaced by direct call to s2c_tl. It is not general enough to
    !           be used for new analysis variables.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(inout) :: column
    type(struct_gsv),        intent(in)    :: statevector
    type(struct_obs),        intent(in)    :: obsSpaceData

    ! Locals:
    integer :: jk, jk2, jgl, headerIndex
    integer :: lonIndex, ila, ierr, subGridIndex
    integer :: extraLongitude
    real(8) :: lat, lon
    real(4) :: lat_r4, lon_r4, lat_deg_r4, lon_deg_r4, xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    real(8) :: dldy, dlw1, dlw2, dlw3, dlw4, dldx, ypos, xpos
    real(8), allocatable ::zgd(:,:,:)
    real(8), pointer :: field_ptr(:,:,:)
    real(8), pointer :: varColumn(:)
    integer :: varIndex
    character(len=4) :: varName

    call utl_tmg_start(30,'--StateToColumn')

    ! Note: We assume here the all the obs between the poles and the last grid points
    !       (i.e. outside the grid) have been moved within the grid by suprep

    if (statevector%hco%global) then
      extraLongitude = 1
    else
      extraLongitude = 0
    end if

    allocate(zgd(statevector%ni+extraLongitude,statevector%nj,statevector%numVarLev))

    zgd(:,:,:)=0.0d0
    call gsv_getField(statevector,field_ptr)
    zgd(1:statevector%ni,1:statevector%nj,1:statevector%numVarLev)= &
         field_ptr(1:statevector%ni,1:statevector%nj,1:statevector%numVarLev)

    !
    !- 1.  Expand field by repeating meridian 1 into into meridian ni+1
    !
    if (extraLongitude == 1) then
      do jk = 1, statevector%numVarLev
        do jgl = 1, statevector%nj
          zgd(statevector%ni+1,jgl,jk) = zgd( 1,jgl,jk)
        end do
      end do
    end if

    !
    !- 2.  Loop over all the headers
    !
    do headerIndex = 1, col_getNumCol(column)

      !- 2.1 Find the obs positin within the analysis grid
      lat    = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)
      lon    = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex)
      lat_r4 = real(lat,4)
      lon_r4 = real(lon,4)
      if (lon_r4.lt.0.0         ) lon_r4 = lon_r4 + 2.0*MPC_PI_R4
      if (lon_r4.ge.2.*MPC_PI_R4) lon_r4 = lon_r4 - 2.0*MPC_PI_R4
      lat_deg_r4 = lat_r4 * MPC_DEGREES_PER_RADIAN_R4 ! Radian To Degree
      lon_deg_r4 = lon_r4 * MPC_DEGREES_PER_RADIAN_R4
      ierr = gpos_getPositionXY(stateVector % hco % EZscintID,   &
                                xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                                lat_deg_r4, lon_deg_r4, subGridIndex)
      xpos = real(xpos_r4,8)
      ypos = real(ypos_r4,8)

      !- Make sure we are within bounds
      if ( ypos < 1.d0                        .or. &
           ypos > real(statevector%nj    , 8) .or. &
           xpos < 1.d0                        .or. &
           xpos > real(statevector%ni + extraLongitude, 8) ) then
        write(*,*) 's2c_bgcheck_bilin: Obs outside local domain for headerIndex = ', &
                   headerIndex
        write(*,*) '  obs    lat, lon position            = ', &
                   Lat*MPC_DEGREES_PER_RADIAN_R8, Lon*MPC_DEGREES_PER_RADIAN_R8
        write(*,*) '  obs    x, y     position            = ', &
                   xpos, ypos
        write(*,*) '  domain x_end, y_end bounds          = ', &
                   statevector%ni + extraLongitude, statevector%nj
        call utl_abort('s2c_bgcheck_bilin')
      end if

      !- 2.2 Find the lower-left grid point next to the observation
      if ( utl_isEqual(xpos, real(statevector%ni + extraLongitude,8)) ) then
        lonIndex = floor(xpos) - 1
      else
        lonIndex = floor(xpos)
      end if

      if ( utl_isEqual(ypos, real(statevector%nj,8)) ) then
        ILA = floor(ypos) - 1
      else
        ILA = floor(ypos)
      end if

      !- 2.3 Compute the 4 weights of the bilinear interpolation
      dldx = xpos - real(lonIndex,8)
      dldy = ypos - real(ILA,8)

      dlw1 = (1.d0-dldx) * (1.d0-dldy)
      dlw2 =       dldx  * (1.d0-dldy)
      dlw3 = (1.d0-dldx) *       dldy
      dlw4 =       dldx  *       dldy

      !- 2.4 Interpolate the model state to the obs point

      do varIndex = 1, vnl_numvarmax
        if (.not. col_varExist(column,trim(vnl_varNameList(varIndex)))) cycle
        varName=trim(vnl_varNameList(varIndex))
        varColumn => col_getColumn(column,headerIndex,varName)

        if(gsv_varExist(statevector,varName)) then
          do jk = 1, gsv_getNumLevFromVarName(statevector,varName)
              jk2=jk+gsv_getOffsetFromVarName(statevector,varName)
              varColumn(jk) =   dlw1*zgd(lonIndex  ,ila,jk2)  &
                                + dlw2*zgd(lonIndex+1,ila,jk2)  &
                                + dlw3*zgd(lonIndex  ,ila+1,jk2)  &
                                + dlw4*zgd(lonIndex+1,ila+1,jk2)
          end do
        end if

        nullify(varColumn)
      end do

    end do

    deallocate(zgd)

    call utl_tmg_stop(30)

  end subroutine s2c_bgcheck_bilin

  !--------------------------------------------------------------------------
  ! s2c_setupHorizInterp
  !--------------------------------------------------------------------------
  subroutine s2c_setupHorizInterp(footprintRadius_r4, interpInfo, &
                                  stateVector, headerIndex, varLevIndex, stepIndex, &
                                  procIndex, numGridpt)
    !
    !:Purpose: To identify the appropriate horizontal interpolation scheme based
    !          on footprint radius value. Then to call the corresponding
    !          subroutine to determine the grid points and their associated
    !          weights.
    !
    implicit none

    ! Arguments:
    real(4)                , intent(in)    :: footprintRadius_r4 ! (metres)
    type(struct_interpInfo), intent(inout) :: interpInfo
    type(struct_gsv)       , intent(in)    :: stateVector
    integer                , intent(in)    :: headerIndex, varLevIndex, stepIndex
    integer                , intent(in)    :: procIndex
    integer                , intent(out)   :: numGridpt(interpInfo%hco%numSubGrid)

    if ( footprintRadius_r4 > 0.0 ) then

      call s2c_setupFootprintInterp(footprintRadius_r4, interpInfo, stateVector, &
                                    headerIndex, varLevIndex, stepIndex, procIndex, numGridpt)

    else if ( utl_isEqual(footprintRadius_r4, bilinearFootprint) ) then

      call s2c_setupBilinearInterp(interpInfo, stateVector, headerIndex, varLevIndex, stepIndex, &
                                   procIndex, numGridpt)

    else if ( utl_isEqual(footprintRadius_r4, lakeFootprint) ) then

      call s2c_setupLakeInterp(interpInfo, stateVector, headerIndex, varLevIndex, stepIndex, &
                               procIndex, numGridpt)

    else if ( utl_isEqual(footprintRadius_r4, nearestNeighbourFootprint) ) then

      call s2c_setupNearestNeighbor(interpInfo, stateVector, headerIndex, varLevIndex, stepIndex, &
                                    procIndex, numGridpt)

    else

      write(*,*) 'footprint radius = ',footprintRadius_r4
      call utl_abort('s2c_setupHorizInterp: footprint radius not permitted')

    end if

  end subroutine s2c_setupHorizInterp

  !------------------------------------------------------------------
  ! s2c_getFootprintRadius
  !------------------------------------------------------------------
  function s2c_getFootprintRadius( obsSpaceData, stateVector, headerIndex ) result(fpr)
    !
    !:Purpose: To determine the footprint radius (metres) of the observation.
    !          In the case of bilinear horizontal interpolation,
    !          the returned footprint is zero (default).
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(in)  :: obsSpaceData
    type(struct_gsv), intent(in)  :: stateVector
    integer         , intent(in)  :: headerIndex
    ! Result:
    real(4)                       :: fpr

    ! Locals:
    character(len=2)  :: obsFamily
    character(len=obs_stnidLength) :: cstnid
    integer           :: codeType

    fpr = bilinearFootprint

    obsFamily = obs_getFamily ( obsSpaceData, headerIndex_opt=headerIndex )
    if ( obsFamily == 'GL' ) then

      cstnid = obs_elem_c ( obsSpaceData, 'STID' , headerIndex )
      codeType = obs_headElem_i( obsSpaceData, OBS_ITY, headerIndex )

      if (index(cstnid,'DMSP') == 1) then

        select case(cstnid)
        case('DMSP15')
          fpr = 27.5e3
        case('DMSP16','DMSP17','DMSP18')
          fpr = 29.0e3
        case DEFAULT
          call utl_abort('s2c_getFootprintRadius: UNKNOWN station id: '//cstnid)
        end select

      else if (cstnid == 'GCOM-W1') then

        fpr = 11.0e3

      else if (cstnid(1:6) == 'METOP-') then

        fpr = 25.0e3

      else if (cstnid == 'noaa-19') then

        fpr = 2.75e3

      else if (cstnid == 'CIS_DAILY') then

        fpr = bilinearFootprint

      else if (cstnid == 'RS1_IMG') then

        fpr = bilinearFootprint

      else if (codtyp_get_name(codeType) == 'iceclake') then

        fpr = lakeFootprint

      else if (cstnid == 'CIS_REGIONAL') then

        fpr = bilinearFootprint

      else if (cstnid(1:3) == 'RCM') then

        fpr = 0.8e3

      else

        call utl_abort('s2c_getFootprintRadius: UNKNOWN station id: '//cstnid)

      end if

    else if (obsFamily == 'HY') then

      fpr = nearestNeighbourFootprint

    else if (obsFamily == 'TO' .and. useFootprintForTovs ) then

      fpr = getTovsFootprintRadius(obsSpaceData, headerIndex, beSilent_opt=.true.)

      ! As safety margin, add 10% to maxGridSpacing before comparing to the footprint radius.
      if ( fpr < 1.1 * real(stateVector%hco%maxGridSpacing,4) ) fpr = bilinearFootprint

    else

      fpr = bilinearFootprint

    end if

  end function s2c_getFootprintRadius

  !--------------------------------------------------------------------------
  ! s2c_rejectZeroWeightObs
  !--------------------------------------------------------------------------
  subroutine s2c_rejectZeroWeightObs(interpInfo, obsSpaceData, myVarLevBeg, myVarLevEnd)
    !
    !:Purpose: To flag an observation in obsSpaceData as being rejected if
    !          it has zero interpolation weight (usually because an ocean
    !          obs is touching land) on any mpi task.
    !
    implicit none

    ! Arguments:
    type(struct_interpInfo), intent(inout) :: interpInfo
    type(struct_obs)       , intent(inout) :: obsSpaceData
    integer                , intent(in)    :: myVarLevBeg
    integer                , intent(in)    :: myVarLevEnd

    ! Locals:
    integer :: numStep, procIndex, stepIndex, headerUsedIndex, headerIndex, varLevIndex
    integer :: numHeader, numHeaderMax, bodyIndexBeg, bodyIndexEnd, bodyIndex
    integer :: subGridIndex, gridptIndex
    integer, save :: numWrites = 0
    logical, allocatable :: allRejectObs(:,:), allRejectObsMpiGlobal(:,:)

    write(*,*) 's2c_rejectZeroWeightObs: Starting'

    numHeader = obs_numheader(obsSpaceData)
    call mmpi_allReduce(numHeader, numHeaderMax, mmpi_max)

    allocate(allRejectObs(numHeaderMax,mmpi_nprocs))
    allocate(allRejectObsMpiGlobal(numHeaderMax,mmpi_nprocs))
    allRejectObs(:,:) = .false.
    allRejectObsMpiGlobal(:,:) = .false.

    numStep = size(interpInfo%stepProcData(1,:))
    do procIndex = 1, mmpi_nprocs
      do stepIndex = 1, numStep
        do headerUsedIndex = 1, interpInfo%allNumHeaderUsed(stepIndex,procIndex)
          headerIndex = interpInfo%stepProcData(procIndex,stepIndex)%allHeaderIndex(headerUsedIndex)
          do varLevIndex = myVarLevBeg, myVarLevEnd
            if (varLevIndex == myVarLevBeg) allRejectObs(headerIndex,procIndex) = .true.
            do subGridIndex = 1, interpInfo%hco%numSubGrid
              do gridptIndex =  &
                   interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerUsedIndex, varLevIndex), &
                   interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex, headerUsedIndex, varLevIndex)
                if (interpInfo%interpWeightDepot(gridptIndex) > 0.0d0) then
                  allRejectObs(headerIndex,procIndex) = .false.
                end if
              end do
            end do
          end do ! varLevIndex
        end do ! headerUsedIndex
      end do ! stepIndex
    end do ! procIndex

    ! do global communication of reject flags
    call mmpi_allReduce(allRejectObs, allRejectObsMpiGlobal, mmpi_lor)

    ! modify obsSpaceData based on reject flags
    do headerIndex = 1, obs_numHeader(obsSpaceData)
      if (allRejectObsMpiGlobal(headerIndex,mmpi_myid+1)) then

        numWrites = numWrites + 1
        if (numWrites < maxNumWrites) then
          write(*,*) 's2c_rejectZeroWeightObs: Rejecting OBS with zero weight, index ', headerIndex
        else if (numWrites == maxNumWrites) then
          write(*,*) 's2c_rejectZeroWeightObs: More rejects, but reached maximum number of writes to the listing.'
        end if

        bodyIndexBeg = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
        bodyIndexEnd = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) + bodyIndexBeg -1
        do bodyIndex = bodyIndexBeg, bodyIndexEnd
          call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
        end do
        call obs_headSet_i(obsSpaceData, OBS_ST1, headerIndex,  &
                    ibset( obs_headElem_i(obsSpaceData, OBS_ST1, headerIndex), 05))
      end if
    end do

    deallocate(allRejectObs)
    deallocate(allRejectObsMpiGlobal)

    write(*,*) 's2c_rejectZeroWeightObs: Finished'

  end subroutine s2c_rejectZeroWeightObs

  !--------------------------------------------------------------------------
  ! s2c_setupBilinearInterp
  !--------------------------------------------------------------------------
  subroutine s2c_setupBilinearInterp(interpInfo, stateVector, headerIndex, varLevIndex, &
                                     stepIndex, procIndex, numGridpt)
    !
    !:Purpose: To determine the grid points and their associated weights
    !          for the bilinear horizontal interpolation. If mask is present
    !          we currently can only handle a single 2D mask (like for sea
    !          ice or SST analysis). Will abort if multiple ocean levels present.
    !
    implicit none

    ! Arguments:
    type(struct_interpInfo), intent(inout) :: interpInfo
    type(struct_gsv)       , intent(in)    :: stateVector
    integer                , intent(in)    :: headerIndex, varLevIndex, stepIndex
    integer                , intent(in)    :: procIndex
    integer                , intent(out)   :: numGridpt(interpInfo%hco%numSubGrid)

    ! Locals:
    integer :: depotIndex
    integer :: ierr, niP1
    integer :: latIndex, lonIndex, latIndex2, lonIndex2, lonIndexP1
    integer :: subGridIndex, subGridForInterp, numSubGridsForInterp
    integer :: ipoint, gridptCount
    integer :: latIndexVec(4), lonIndexVec(4)
    logical :: mask(2,2)
    real(8) :: WeightVec(4)
    real(8) :: dldx, dldy
    real(8) :: weightsSum
    real(4) :: lon_deg_r4, lat_deg_r4
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    integer, parameter :: leftIndex = 1, rightIndex = 2, bottomIndex = 1, topIndex = 2
    logical :: isCloudVariable

    numGridpt(:) = 0

    lat_deg_r4 = real(interpInfo%stepProcData(procIndex, stepIndex)%allLat(headerIndex, varLevIndex) *  &
                 MPC_DEGREES_PER_RADIAN_R8)
    lon_deg_r4 = real(interpInfo%stepProcData(procIndex, stepIndex)%allLon(headerIndex, varLevIndex) *  &
                 MPC_DEGREES_PER_RADIAN_R8)
    ierr = gpos_getPositionXY(stateVector%hco%EZscintID,   &
                              xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                              lat_deg_r4, lon_deg_r4, subGridIndex)

    ! Allow for periodicity in Longitude for global Gaussian grid
    if ( stateVector%hco%grtyp == 'G' .or. &
         (stateVector%hco%grtyp == 'Z' .and. stateVector%hco%global) ) then
      niP1 = statevector%ni + 1
    else
      niP1 = statevector%ni
    end if

    ! Find the lower-left grid point next to the observation
    if ( xpos_r4 >= real(niP1) ) then
      xpos_r4 = real(niP1)
      lonIndex = niP1 - 1
    else if ( xpos_r4 < 1.0 ) then
      xpos_r4 = 1.0
      lonIndex = 1
    else
      lonIndex = floor(xpos_r4)
    end if
    if ( xpos2_r4 >= real(niP1) ) then
      xpos2_r4 = real(niP1)
      lonIndex2 = niP1 - 1
    else if ( xpos2_r4 < 1.0 ) then
      xpos2_r4 = 1.0
      lonIndex2 = 1
    else
      lonIndex2 = floor(xpos2_r4)
    end if

    if ( ypos_r4 >= real(statevector%nj) ) then
      ypos_r4 = real(statevector%nj)
      latIndex = statevector%nj - 1
    else if ( ypos_r4 < 1.0 ) then
      ypos_r4 = 1.0
      latIndex = 1
    else
      latIndex = floor(ypos_r4)
    end if
    if ( ypos2_r4 >= real(statevector%nj) ) then
      ypos2_r4 = real(statevector%nj)
      latIndex2 = statevector%nj - 1
    else if ( ypos2_r4 < 1.0 ) then
      ypos2_r4 = 1.0
      latIndex2 = 1
    else
      latIndex2 = floor(ypos2_r4)
    end if

    if ( stateVector%hco%grtyp == 'U' ) then
      if ( utl_isEqual(ypos_r4, real(stateVector%nj/2,4)) ) then
        latIndex = floor(ypos_r4) - 1
      end if
      if ( utl_isEqual(ypos2_r4, real(stateVector%nj/2,4)) ) then
        latIndex2 = floor(ypos2_r4) - 1
      end if
    end if

    ! Handle periodicity in longitude
    lonIndexP1 = lonIndex + 1
    if ( lonIndexP1 == statevector%ni + 1 ) lonIndexP1 = 1

    ! Check if location is in between Yin and Yang (should not happen)
    if ( stateVector%hco%grtyp == 'U' ) then
      if ( ypos_r4 > real(stateVector%nj/2) .and.  &
           ypos_r4 < real((stateVector%nj/2)+1) ) then
        write(*,*) 's2c_setupBilinearInterp: WARNING, obs position in between Yin and Yang!!!'
        write(*,*) '   xpos, ypos = ', xpos_r4, ypos_r4
      end if
      if ( ypos2_r4 > real(stateVector%nj/2) .and.  &
           ypos2_r4 < real((stateVector%nj/2)+1) ) then
        write(*,*) 's2c_setupBilinearInterp: WARNING, obs position in between Yin and Yang!!!'
        write(*,*) '   xpos2, ypos2 = ', xpos2_r4, ypos2_r4
      end if
    end if

    if ( subGridIndex == 3 ) then
      ! both subGrids involved in interpolation, so first treat subGrid 1
      numSubGridsForInterp = 2
      subGridIndex = 1
    else
      ! only 1 subGrid involved in interpolation
      numSubGridsForInterp = 1
    end if

    if ( stateVector%oceanMask%maskPresent ) then
      ! abort if 3D mask is present, since we may not handle this situation correctly
      if ( stateVector%oceanMask%nLev > 1 ) then
        call utl_abort('s2c_setupBilinearInterp: 3D mask present - this case not properly handled')
      end if
      ! extract the ocean mask
      mask(leftIndex ,bottomIndex) = stateVector%oceanMask%mask(lonIndex  ,latIndex    ,1)
      mask(rightIndex,bottomIndex) = stateVector%oceanMask%mask(lonIndexP1,latIndex    ,1)
      mask(leftIndex ,topIndex   ) = stateVector%oceanMask%mask(lonIndex  ,latIndex + 1,1)
      mask(rightIndex,topIndex   ) = stateVector%oceanMask%mask(lonIndexP1,latIndex + 1,1)
    else
      mask(:,:) = .true.
    end if

    do subGridForInterp = 1, numSubGridsForInterp

      WeightVec(:) = 0
      gridptCount = 0

      ! Compute the 4 weights of the bilinear interpolation
      if ( subGridForInterp == 1 ) then
        ! when only 1 subGrid involved, subGridIndex can be 1 or 2
        dldx = real(xpos_r4,8) - real(lonIndex,8)
        dldy = real(ypos_r4,8) - real(latIndex,8)
      else
        ! when 2 subGrids, subGridIndex is set to 1 for 1st iteration, 2 for second
        subGridIndex = 2
        lonIndex = lonIndex2
        latIndex = latIndex2
        lonIndexP1 = lonIndex2 + 1
        dldx = real(xpos2_r4,8) - real(lonIndex,8)
        dldy = real(ypos2_r4,8) - real(latIndex,8)
      end if

      if (NNInterpForCloudVars) then
        if (varLevIndex > 0) then
          isCloudVariable = vnl_isCloudVar( gsv_getVarNameFromVarLev(stateVector,varLevIndex) )
        else
          isCloudVariable = .false.
        end if
        if (isCloudVariable) then
          dldx = real(nint(dldx), 8)
          dldy = real(nint(dldy), 8)
        end if
      end if

      if (NNInterpForAllVars) then
        dldx = real(nint(dldx), 8)
        dldy = real(nint(dldy), 8)
      end if

      if ( mask(leftIndex ,bottomIndex) ) then
        gridptCount = gridptCount + 1
        latIndexVec(gridptCount) = latIndex
        lonIndexVec(gridptCount) = lonIndex
        WeightVec(gridptCount) = (1.d0-dldx) * (1.d0-dldy)
      end if

      if ( mask(rightIndex,bottomIndex) ) then
        gridptCount = gridptCount + 1
        latIndexVec(gridptCount) = latIndex
        lonIndexVec(gridptCount) = lonIndexP1
        WeightVec(gridptCount) =       dldx  * (1.d0-dldy)
      end if

      if ( mask(leftIndex ,topIndex   ) ) then
        gridptCount = gridptCount + 1
        latIndexVec(gridptCount) = latIndex + 1
        lonIndexVec(gridptCount) = lonIndex
        WeightVec(gridptCount) = (1.d0-dldx) *       dldy
      end if

      if ( mask(rightIndex,topIndex   ) ) then
        gridptCount = gridptCount + 1
        latIndexVec(gridptCount) = latIndex + 1
        lonIndexVec(gridptCount) = lonIndexP1
        WeightVec(gridptCount) =       dldx  *       dldy
      end if

      weightsSum = sum(WeightVec(1:gridptCount))
      if ( weightsSum > 0.d0 ) then
        WeightVec(1:gridptCount) = WeightVec(1:gridptCount) / weightsSum
      end if

      ! divide weight by number of subGrids
      WeightVec(1:gridptCount) = WeightVec(1:gridptCount) / real(numSubGridsForInterp,8)

      if ( allocated(interpInfo%interpWeightDepot) ) then

        depotIndex = interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, varLevIndex)

        do ipoint=1,gridptCount

          interpInfo%interpWeightDepot(depotIndex) = WeightVec(ipoint)
          interpInfo%latIndexDepot(depotIndex)     = latIndexVec(ipoint)
          interpInfo%lonIndexDepot(depotIndex)     = lonIndexVec(ipoint)
          depotIndex = depotIndex + 1

        end do

      end if

      numGridpt(subGridIndex) = gridptCount

    end do ! subGrid

  end subroutine s2c_setupBilinearInterp

  !--------------------------------------------------------------------------
  ! s2c_setupFootprintInterp
  !--------------------------------------------------------------------------
  subroutine s2c_setupFootprintInterp(fpr, interpInfo, stateVector, headerIndex, &
                                      varLevIndex, stepIndex, procIndex, numGridpt)
    !
    !:Purpose: To determine the grid points and their associated weights
    !          for the footprint horizontal interpolation.
    !
    implicit none

    ! Arguments:
    real(4)                , intent(in)    :: fpr ! footprint radius (metres)
    type(struct_interpInfo), intent(inout) :: interpInfo
    type(struct_gsv)       , intent(in)    :: stateVector
    integer                , intent(in)    :: headerIndex, varLevIndex, stepIndex
    integer                , intent(in)    :: procIndex
    integer                , intent(out)   :: numGridpt(interpInfo%hco%numSubGrid)

    ! Locals:
    integer :: depotIndex
    integer :: ierr
    integer :: latIndexCentre, lonIndexCentre
    integer :: subGridIndex, numLocalGridptsFoundSearch
    real(4) :: lonObs_deg_r4, latObs_deg_r4
    real(8) :: lonObs, latObs
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    integer :: ipoint, gridptCount
    integer :: lonIndex, latIndex, resultsIndex, gridIndex
    integer :: lonIndexVec(maxNumLocalGridptsSearch), latIndexVec(maxNumLocalGridptsSearch)
    type(kdtree2_result)      :: searchResults(maxNumLocalGridptsSearch)
    real(kdkind)              :: refPosition(3), maxRadiusSquared
    type(kdtree2), pointer    :: tree

    numGridpt(:) = 0

    ! Determine the grid point nearest the observation.

    latObs = interpInfo % stepProcData(procIndex, stepIndex) % allLat(headerIndex, varLevIndex)
    lonObs = interpInfo % stepProcData(procIndex, stepIndex) % allLon(headerIndex, varLevIndex)

    latObs_deg_r4 = real(latObs * MPC_DEGREES_PER_RADIAN_R8)
    lonObs_deg_r4 = real(lonObs * MPC_DEGREES_PER_RADIAN_R8)
    ierr = gpos_getPositionXY(stateVector%hco%EZscintID,   &
                              xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                              latObs_deg_r4, lonObs_deg_r4, subGridIndex)

    lonIndexCentre = nint(xpos_r4)
    latIndexCentre = nint(ypos_r4)

    if ( subGridIndex == 3 ) then
      write(*,*) 's2c_setupFootprintInterp: revise code'
      call utl_abort('s2c_setupFootprintInterp: both subGrids involved in interpolation.')
    end if

    ! Return if observation is not on the grid, or masked.
    if ( lonIndexCentre < 1 .or. lonIndexCentre > statevector%hco%ni .or.  &
         latIndexCentre < 1 .or. latIndexCentre > statevector%hco%nj ) return

    if ( stateVector%oceanMask%maskPresent ) then
      ! abort if 3D mask is present, since we may not handle this situation correctly
      if ( stateVector%oceanMask%nLev > 1 ) then
        call utl_abort('s2c_setupFootprintInterp: 3D mask present - this case not properly handled')
      end if

      if ( .not. stateVector%oceanMask%mask(lonIndexCentre,latIndexCentre,1) ) return
    end if

    ! do the search
    maxRadiusSquared = real(fpr,8) ** 2
    refPosition(:) = kdtree2_3dPosition(lonObs, latObs)
    nullify(tree)
    if ( interpInfo%inputStateVectorType == 'nl' ) then
      if ( associated(tree_nl) ) then
        tree => tree_nl
      else
        call utl_abort('s2c_setupFootprintInterp: tree_nl is not allocated!')
      end if
    else if ( interpInfo%inputStateVectorType == 'tl' .or. &
              interpInfo%inputStateVectorType == 'ad' ) then
      if ( associated(tree_tlad) ) then
        tree => tree_tlad
      else
        call utl_abort('s2c_setupFootprintInterp: tree_tlad is not allocated!')
      end if
    end if
    call kdtree2_r_nearest(tp=tree, qv=refPosition, r2=maxRadiusSquared, &
                           nfound=numLocalGridptsFoundSearch, &
                           nalloc=maxNumLocalGridptsSearch, &
                           results=searchResults)

    if (numLocalGridptsFoundSearch > maxNumLocalGridptsSearch ) then
      call utl_abort('s2c_setupFootprintInterp: the parameter maxNumLocalGridptsSearch must be increased')
    else if ( numLocalGridptsFoundSearch < minNumLocalGridptsSearch .and. useFootprintForTovs ) then
      write(*,*) 's2c_setupFootprintInterp: Warning! For TOVS headerIndex=', headerIndex, &
                 ' number of grid points found within footprint radius=', fpr, ' is less than ', &
                 minNumLocalGridptsSearch
    end if

    ! ensure at least the nearest neighbor is included in lonIndexVec/latIndexVec
    ! if footprint size is smaller than the grid spacing.
    gridptCount = 1
    lonIndexVec(gridptCount) = lonIndexCentre
    latIndexVec(gridptCount) = latIndexCentre

    ! fill the rest of lonIndexVec/latIndexVec
    gridLoop1: do resultsIndex = 1, numLocalGridptsFoundSearch
      gridIndex = searchResults(resultsIndex)%idx
      if ( gridIndex < 1 .or. gridIndex > statevector%hco%ni * statevector%hco%nj ) then
        write(*,*) 's2c_setupFootprintInterp: gridIndex=', gridIndex
        call utl_abort('s2c_setupFootprintInterp: gridIndex out of bound.')
      end if

      latIndex = (gridIndex - 1) / statevector%hco%ni + 1
      lonIndex = gridIndex - (latIndex - 1) * statevector%hco%ni
      if ( lonIndex < 1 .or. lonIndex > statevector%hco%ni .or. &
           latIndex < 1 .or. latIndex > statevector%hco%nj ) then
        write(*,*) 's2c_setupFootprintInterp: lonIndex=', lonIndex, ',latIndex=', latIndex
        call utl_abort('s2c_setupFootprintInterp: lonIndex/latIndex out of bound.')
      end if

      if ( stateVector%oceanMask%maskPresent ) then
        if ( .not. stateVector%oceanMask%mask(lonIndex,latIndex,1) ) cycle gridLoop1
      end if

      if ( lonIndex == lonIndexCentre .and. latIndex == latIndexCentre ) cycle gridLoop1

      gridptCount = gridptCount + 1
      lonIndexVec(gridptCount) = lonIndex
      latIndexVec(gridptCount) = latIndex
    end do gridLoop1

    if ( allocated(interpInfo%interpWeightDepot) ) then

      depotIndex = interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, varLevIndex)

      do ipoint = 1, gridptCount

        interpInfo%interpWeightDepot(depotIndex) = 1.0d0 / real(gridptCount,8)
        interpInfo%latIndexDepot(depotIndex)     = latIndexVec(ipoint)
        interpInfo%lonIndexDepot(depotIndex)     = lonIndexVec(ipoint)
        depotIndex = depotIndex + 1

      end do

    end if

    numGridpt(subGridIndex) = gridptCount

  end subroutine s2c_setupFootprintInterp

  !--------------------------------------------------------------------------
  ! s2c_setupLakeInterp
  !--------------------------------------------------------------------------
  subroutine s2c_setupLakeInterp(interpInfo, stateVector, headerIndex, varLevIndex, &
                                 stepIndex, procIndex, numGridpt)
    !
    !:Purpose: To determine the grid points and their associated weights
    !          for the lake horizontal interpolation.
    !
    implicit none

    ! Arguments:
    type(struct_interpInfo), intent(inout) :: interpInfo
    type(struct_gsv)       , intent(in)    :: stateVector
    integer                , intent(in)    :: headerIndex, varLevIndex, stepIndex
    integer                , intent(in)    :: procIndex
    integer                , intent(out)   :: numGridpt(interpInfo%hco%numSubGrid)

    ! Locals:
    integer :: depotIndex
    integer :: ierr
    integer :: latIndexCentre, lonIndexCentre
    integer :: subGridIndex, subGridForInterp, numSubGridsForInterp
    real(4) :: lon_deg_r4, lat_deg_r4
    real(8) :: lon_rad, lat_rad
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    integer :: ipoint, gridptCount
    integer :: lakeCount
    integer :: lonIndex, latIndex, lakeIndex
    integer :: lonIndexVec(statevector%ni*statevector%nj), latIndexVec(statevector%ni*statevector%nj)
    logical :: reject, lake(statevector%ni,statevector%nj)
    integer :: k, l

    if ( stateVector%hco%grtyp == 'U' ) then
      call utl_abort('s2c_setupLakeInterp: Yin-Yang grid not supported')
    end if

    if ( .not.stateVector%oceanMask%maskPresent ) then
      call utl_abort('s2c_setupLakeInterp: Only compatible when mask present')
    end if

    numGridpt(:) = 0

    reject = .false.

    numGridpt(:) = 0

    ! Determine the grid point nearest the observation.

    lat_rad = interpInfo%stepProcData(procIndex, stepIndex)%allLat(headerIndex, varLevIndex)
    lon_rad = interpInfo%stepProcData(procIndex, stepIndex)%allLon(headerIndex, varLevIndex)
    lat_deg_r4 = real(lat_rad * MPC_DEGREES_PER_RADIAN_R8)
    lon_deg_r4 = real(lon_rad * MPC_DEGREES_PER_RADIAN_R8)
    ierr = gpos_getPositionXY(stateVector%hco%EZscintID,   &
                              xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                              lat_deg_r4, lon_deg_r4, subGridIndex)

    lonIndexCentre = nint(xpos_r4)
    latIndexCentre = nint(ypos_r4)

    if ( subGridIndex == 3 ) then
      write(*,*) 's2c_setupLakeInterp: revise code'
      call utl_abort('s2c_setupLakeInterp: both subGrids involved in interpolation.')
      numSubGridsForInterp = 2
      subGridIndex = 1
    else
      ! only 1 subGrid involved in interpolation
      numSubGridsForInterp = 1
    end if

    do subGridForInterp = 1, numSubGridsForInterp

      gridptCount = 0

      ! It can happen that the lake location is closest to a grid point
      ! where MASK(I,J) = .false. while there are other grid points for the
      ! same lake where MASK(I,J) = .true.. Code needs modifications
      ! for this case.

      ! If observation is not on the grid, don't use it.
      if ( lonIndexCentre < 1 .or. lonIndexCentre > statevector%ni .or.  &
           latIndexCentre < 1 .or. latIndexCentre > statevector%nj ) reject = .true.

      if ( .not. stateVector%oceanMask%mask(lonIndexCentre,latIndexCentre,1) ) reject = .true.

      if ( .not. reject ) then

        lake(:,:) = .false.
        lake(lonIndexCentre,latIndexCentre) = .true.
        gridptCount = 1
        lonIndexVec(gridptCount) = lonIndexCentre
        latIndexVec(gridptCount) = latIndexCentre

        lakeCount = 0

        do while(lakeCount /= gridptCount)

          do lakeIndex = lakeCount+1, gridptCount

            if(lakeIndex == lakeCount+1) lakeCount = gridptCount

            k = lonIndexVec(lakeIndex)
            l = latIndexVec(lakeIndex)

            do latIndex = max(1,l-1),min(l+1,statevector%nj)
              do lonIndex = max(1,k-1),min(k+1,statevector%ni)
                if(stateVector%oceanMask%mask(lonIndex,latIndex,1) .and. .not. lake(lonIndex,latIndex)) then
                  lake(lonIndex,latIndex) = .true.
                  gridptCount = gridptCount + 1
                  lonIndexVec(gridptCount) = lonIndex
                  latIndexVec(gridptCount) = latIndex
                end if
              end do
            end do

          end do

        end do

        if ( allocated(interpInfo%interpWeightDepot) ) then

          depotIndex = interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, varLevIndex)

          do ipoint=1,gridptCount

            interpInfo%interpWeightDepot(depotIndex) = 1.0d0 / real(gridptCount,8)
            interpInfo%latIndexDepot(depotIndex)     = latIndexVec(ipoint)
            interpInfo%lonIndexDepot(depotIndex)     = lonIndexVec(ipoint)
            depotIndex = depotIndex + 1

          end do

        end if

        numGridpt(subGridIndex) = gridptCount

      end if ! not reject

    end do ! subGrid

  end subroutine s2c_setupLakeInterp

  !--------------------------------------------------------------------------
  ! s2c_setupNearestNeighbor
  !--------------------------------------------------------------------------
  subroutine s2c_setupNearestNeighbor(interpInfo, stateVector, headerIndex, varLevIndex, &
                                      stepIndex, procIndex, numGridpt)
    !
    !:Purpose: Determine the nearest grid points to the observations location
    !
    implicit none

    ! Arguments:
    type(struct_interpInfo), intent(inout) :: interpInfo
    type(struct_gsv)       , intent(in)    :: stateVector
    integer                , intent(in)    :: headerIndex, varLevIndex, stepIndex, procIndex
    integer                , intent(out)   :: numGridpt(interpInfo%hco%numSubGrid)

    ! Locals:
    integer :: depotIndex
    integer :: ierr
    integer :: latIndex, lonIndex
    integer :: subGridIndex
    real(4) :: lon_deg_r4, lat_deg_r4
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4

    if ( stateVector%hco%grtyp == 'U' ) then
      call utl_abort('s2c_setupNearestNeighbor: Yin-Yang grid not supported')
    end if

    numGridpt(:) = 0

    lat_deg_r4 = real(interpInfo%stepProcData(procIndex, stepIndex)%allLat(headerIndex, varLevIndex) *  &
                 MPC_DEGREES_PER_RADIAN_R8)
    lon_deg_r4 = real(interpInfo%stepProcData(procIndex, stepIndex)%allLon(headerIndex, varLevIndex) *  &
                 MPC_DEGREES_PER_RADIAN_R8)

    ierr = gpos_getPositionXY(stateVector%hco%EZscintID,   &
                              xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                              lat_deg_r4, lon_deg_r4, subGridIndex)

    latIndex = nint(ypos_r4)
    lonIndex = nint(xpos_r4)

    ! Handle periodicity in longitude
    if ( lonIndex == statevector%ni+1 .and. stateVector%hco%grtyp == 'G' ) lonIndex = 1

    ! Test bounds
    if ( lonIndex < 1 .or. lonIndex > statevector%ni .or. &
         latIndex < 1 .or. latIndex > statevector%nj  ) then

      write(*,*) 's2c_setupNearestNeighbor: observation out of bounds'

    else

      if ( allocated(interpInfo%interpWeightDepot) ) then

        depotIndex = interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, varLevIndex)

        interpInfo%interpWeightDepot(depotIndex) = 1.d0
        interpInfo%latIndexDepot    (depotIndex) = latIndex
        interpInfo%lonIndexDepot    (depotIndex) = lonIndex

      end if

      numGridpt(subGridIndex) = 1

    end if

  end subroutine s2c_setupNearestNeighbor

  !--------------------------------------------------------------------------
  ! checkColumnStatevectorMatch
  !--------------------------------------------------------------------------
  subroutine checkColumnStatevectorMatch(column,statevector)
    !
    !:Purpose: To check column and statevector have identical nk and variables.
    !
    implicit none

    ! Arguments:
    type(struct_gsv)       , intent(in) :: statevector
    type(struct_columnData), intent(in) :: column

    ! Locals:
    integer :: varLevIndex

    ! check column/statevector have same nk
    if ( col_getNumVarLev(column) /= gsv_getNumVarLev(statevector) ) then
      write(*,*) 'checkColumnStatevectorMatch: col_getNumVarLev(column), gsv_getNumVarLev(statevector)', &
                 col_getNumVarLev(column), gsv_getNumVarLev(statevector)
      call utl_abort('checkColumnStatevectorMatch: col_getNumVarLev(column) /= gsv_getNumVarLev(statevector)')
    end if

    ! loop through k and check varNames are same between column/statevector
    do varLevIndex = 1, col_getNumVarLev(column)
      if (gsv_getVarNameFromVarLev(statevector,varLevIndex) /= col_getVarNameFromVarLev(column,varLevIndex)) then
        write(*,*) 'checkColumnStatevectorMatch: varLevIndex, varname in statevector and column: ', varLevIndex, &
                   gsv_getVarNameFromVarLev(statevector,varLevIndex), col_getVarNameFromVarLev(column,varLevIndex)
        call utl_abort('checkColumnStatevectorMatch: varname in column and statevector do not match')
      end if
    end do

  end subroutine checkColumnStatevectorMatch

  !--------------------------------------------------------------------------
  ! latlonChecks
  !--------------------------------------------------------------------------
  subroutine latlonChecks( obsSpaceData, hco, headerIndex, rejectOutsideObs, &
    latLev_T, lonLev_T, latLev_M, lonLev_M, latLev_S_opt, lonLev_S_opt )
    !
    !:Purpose: To check if the obs are inside the domain.
    !
    implicit none

    ! Arguments:
    type(struct_obs),  intent(inout) :: obsSpaceData
    type(struct_hco),  intent(in)    :: hco
    integer,           intent(in)    :: headerIndex
    logical,           intent(in)    :: rejectOutsideObs
    real(8),           intent(inout) :: latLev_T(:)
    real(8),           intent(inout) :: lonLev_T(:)
    real(8),           intent(inout) :: latLev_M(:)
    real(8),           intent(inout) :: lonLev_M(:)
    real(8), optional, intent(inout) :: latLev_S_opt
    real(8), optional, intent(inout) :: lonLev_S_opt

    ! Locals:
    integer :: ierr
    integer :: bodyIndex, bodyIndexBeg, bodyIndexEnd, niP1, subGridIndex
    integer :: nlev_T, nlev_M
    real(4) :: lon_r4, lat_r4, lon_deg_r4, lat_deg_r4
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    logical :: latlonOutsideGrid, rejectHeader
    integer :: gdllfxy

    ! Allow for periodicity in Longitude for global Gaussian grid
    if ( hco%grtyp == 'G' .or. (hco%grtyp == 'Z' .and. hco%global) ) then
      niP1 = hco%ni + 1
    else
      niP1 = hco%ni
    end if

    nlev_T = size(latLev_T)
    nlev_M = size(latLev_M)

    ! check if lat/lon of last thermo level is outside domain.
    rejectHeader = .false.
    lat_r4 = real(latLev_T(nlev_T),4)
    lon_r4 = real(lonLev_T(nlev_T),4)

    ! TODO: simplify the floating point precision conversions
    !     lat_deg_r4 = lat_r4 * MPC_DEGREES_PER_RADIAN_R4
    !     lon_deg_r4 = lon_r4 * MPC_DEGREES_PER_RADIAN_R4
    lat_deg_r4 = real( real(lat_r4,8) * MPC_DEGREES_PER_RADIAN_R8, 4)
    lon_deg_r4 = real( real(lon_r4,8) * MPC_DEGREES_PER_RADIAN_R8, 4)
    ierr = gpos_getPositionXY(hco%EZscintID,   &
                              xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                              lat_deg_r4, lon_deg_r4, subGridIndex)

    latlonOutsideGrid = ( xpos_r4 < 1.0        .or. &
                          xpos_r4 > real(niP1) .or. &
                          ypos_r4 < 1.0        .or. &
                          ypos_r4 > real(hco%nj) )

    if ( latlonOutsideGrid .and. rejectOutsideObs ) then
      rejectHeader = .true.
    end if

    !  check if lat/lon of last momentum level is outside domain.
    if ( .not. rejectHeader ) then
      lat_r4 = real(latLev_M(nlev_M),4)
      lon_r4 = real(lonLev_M(nlev_M),4)

      ! TODO: simplify the floating point precision conversions
      !     lat_deg_r4 = lat_r4 * MPC_DEGREES_PER_RADIAN_R4
      !     lon_deg_r4 = lon_r4 * MPC_DEGREES_PER_RADIAN_R4
      lat_deg_r4 = real( real(lat_r4,8) * MPC_DEGREES_PER_RADIAN_R8, 4)
      lon_deg_r4 = real( real(lon_r4,8) * MPC_DEGREES_PER_RADIAN_R8, 4)
      ierr = gpos_getPositionXY(hco%EZscintID,   &
                                xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                                lat_deg_r4, lon_deg_r4, subGridIndex)

      latlonOutsideGrid = ( xpos_r4 < 1.0        .or. &
                            xpos_r4 > real(niP1) .or. &
                            ypos_r4 < 1.0        .or. &
                            ypos_r4 > real(hco%nj) )

      if ( latlonOutsideGrid .and. rejectOutsideObs ) then
        rejectHeader = .true.
      end if
    end if

    !  check if lat/lon of surface level is outside domain.
    if ( present(latLev_S_opt) .and. present(lonLev_S_opt) .and. .not. rejectHeader ) then
      lat_r4 = real(latLev_S_opt,4)
      lon_r4 = real(lonLev_S_opt,4)

      ! TODO: simplify the floating point precision conversions
      !     lat_deg_r4 = lat_r4 * MPC_DEGREES_PER_RADIAN_R4
      !     lon_deg_r4 = lon_r4 * MPC_DEGREES_PER_RADIAN_R4
      lat_deg_r4 = real( real(lat_r4,8) * MPC_DEGREES_PER_RADIAN_R8, 4)
      lon_deg_r4 = real( real(lon_r4,8) * MPC_DEGREES_PER_RADIAN_R8, 4)
      ierr = gpos_getPositionXY(hco%EZscintID,   &
                                xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                                lat_deg_r4, lon_deg_r4, subGridIndex)

      latlonOutsideGrid = ( xpos_r4 < 1.0        .or. &
                            xpos_r4 > real(niP1) .or. &
                            ypos_r4 < 1.0        .or. &
                            ypos_r4 > real(hco%nj) )

      if ( latlonOutsideGrid .and. rejectOutsideObs ) then
        rejectHeader = .true.
      end if
    end if

    if ( rejectHeader ) then
      ! The observation is outside the domain.
      ! With a LAM trial field we must discard this observation
      write(*,*) 'latlonChecks: Rejecting OBS outside the hco domain, ', headerIndex
      write(*,*) '  position : ', lat_deg_r4, lon_deg_r4, ypos_r4, xpos_r4

      bodyIndexBeg = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
      bodyIndexEnd = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) + bodyIndexBeg -1
      do bodyIndex = bodyIndexBeg, bodyIndexEnd
        call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
      end do
      call obs_headSet_i(obsSpaceData, OBS_ST1, headerIndex,  &
           ibset( obs_headElem_i(obsSpaceData, OBS_ST1, headerIndex), 05))

      ! Assign domain mid-point lat-lon to this header
      if ( hco%grtyp == 'Y' ) then
        lat_deg_r4 = hco%lat2d_4(hco%ni/2,hco%nj/2)
        lon_deg_r4 = hco%lon2d_4(hco%ni/2,hco%nj/2)
      else
        xpos_r4 = real(hco%ni)/2.0
        ypos_r4 = real(hco%nj)/2.0
        ierr = gdllfxy(hco%EZscintID, lat_deg_r4, lon_deg_r4, &
                       xpos_r4, ypos_r4, 1)
      end if

      lonLev_T(:) = real(lon_deg_r4 * MPC_RADIANS_PER_DEGREE_R4,8)
      latLev_T(:) = real(lat_deg_r4 * MPC_RADIANS_PER_DEGREE_R4,8)
      lonLev_M(:) = real(lon_deg_r4 * MPC_RADIANS_PER_DEGREE_R4,8)
      latLev_M(:) = real(lat_deg_r4 * MPC_RADIANS_PER_DEGREE_R4,8)
      if (present(lonLev_S_opt) .and. present(latLev_S_opt)) then
        lonLev_S_opt = real(lon_deg_r4 * MPC_RADIANS_PER_DEGREE_R4,8)
        latLev_S_opt = real(lat_deg_r4 * MPC_RADIANS_PER_DEGREE_R4,8)
      end if

    end if

  end subroutine latlonChecks

  !--------------------------------------------------------------------------
  ! getTovsFootprintRadius
  !--------------------------------------------------------------------------
  function getTovsFootprintRadius(obsSpaceData, headerIndex, beSilent_opt) result(footPrintRadius_r4)
    !
    !:Purpose: calculate foot-print radius for TOVS observations
    !
    implicit none

    ! Arguments:
    type(struct_obs),  intent(in) :: obsSpaceData
    integer         ,  intent(in) :: headerIndex
    logical, optional, intent(in) :: beSilent_opt
    ! Result:
    real(4)                       :: footPrintRadius_r4

    ! Locals:
    integer :: codtyp, sensorIndex
    real(8) :: fovAngularDiameter, satHeight, footPrintRadius
    character(len=codtyp_name_length) :: instrumName
    logical :: beSilent

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .true.
    end if

    ! get nominal satellite height
    sensorIndex = tvs_lsensor(headerIndex)
    satHeight = tvs_coefs(sensorIndex)%coef%fc_sat_height

    ! FOV angular diameter
    codtyp = obs_headElem_i( obsSpaceData, OBS_ITY, headerIndex )
    instrumName = codtyp_get_name(codtyp)
    select case(trim(instrumName))
    case('amsua')
      fovAngularDiameter = 3.3d0
    case('amsub')
      fovAngularDiameter = 1.1d0
    case('mhs')
      fovAngularDiameter = 10.0d0 / 9.0d0
    case('airs')
      fovAngularDiameter = 1.1d0
    case('iasi')
      fovAngularDiameter = 14.65d0 / 1000.0d0 * MPC_DEGREES_PER_RADIAN_R8
    case('radianceclear')
      fovAngularDiameter = 0.125d0
    case('ssmis')
      fovAngularDiameter = 1.2d0
    case('atms')
      fovAngularDiameter = 2.2d0
    case('cris')
      fovAngularDiameter = 14.0d0 / 824.0d0 * MPC_DEGREES_PER_RADIAN_R8
    case default
      fovAngularDiameter = -1.0d0
    end select

    if ( fovAngularDiameter < 0.0d0 ) then
      footPrintRadius_r4 = bilinearFootprint
    else
      ! get foot print radius (meter) from angular diameter
      footPrintRadius = 0.5d0 * fovAngularDiameter * MPC_RADIANS_PER_DEGREE_R8 * satHeight * 1000
      footPrintRadius_r4 = real(footPrintRadius,4)
    end if

    if ( .not. beSilent ) then
      write(*,*) 'getTovsFootprintRadius: sensorIndex=', sensorIndex, &
                ',satHeight=', satHeight, ',fovAngularDiameter=', fovAngularDiameter, ',codtyp=', codtyp, &
                ',footPrintRadius=', footPrintRadius_r4
    end if

  end function getTovsFootprintRadius

  ! -------------------------------------------------------------
  ! s2c_getWeightsAndGridPointIndexes
  ! -------------------------------------------------------------
  subroutine s2c_getWeightsAndGridPointIndexes(headerIndex, varLevIndex, stepIndex, procIndex, &
                                               interpWeight, latIndex, lonIndex, gridptCount)
    ! :Purpose: Returns the weights and grid point indexes for a single observation.
    !
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: headerIndex
    integer, intent(in)  :: varLevIndex
    integer, intent(in)  :: stepIndex
    integer, intent(in)  :: procIndex
    real(8), intent(out) :: interpWeight(:)
    integer, intent(out) :: latIndex(:)
    integer, intent(out) :: lonIndex(:)
    integer, intent(out) :: gridptCount

    ! Locals:
    integer :: indexBeg, indexEnd, gridptIndex
    integer :: subGridIndex, maxGridpt

    call utl_tmg_start(30,'--StateToColumn')

    maxGridpt = size( interpWeight )

    gridptCount = 0

    if ( interpInfo_tlad%stepProcData(procIndex, stepIndex)%allHeaderIndex(headerIndex) /= headerIndex ) then
      call utl_abort('s2c_getWeightsAndGridPointIndexes: headerUsedIndex and headerIndex differ.'//    &
                     ' If using multiple time steps in the assimilation window,'//                     &
                     ' the code needs to be modified to convert values of headerIndex into headerUsedIndex.')
    end if

    subGrid_loop: do subGridIndex = 1, interpInfo_tlad%hco%numSubGrid

      indexBeg = interpInfo_tlad%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, varLevIndex)
      indexEnd = interpInfo_tlad%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex, headerIndex, varLevIndex)

      if ( indexEnd < IndexBeg ) cycle subGrid_loop

      do gridptIndex = indexBeg, indexEnd

        gridptCount = gridptCount + 1

        if ( gridptCount > maxGridpt ) then
          call utl_abort('s2c_getWeightsAndGridPointIndexes: maxGridpt must be increased')
        end if

        lonIndex(gridptCount) = interpInfo_tlad%lonIndexDepot(gridptIndex)
        latIndex(gridptCount) = interpInfo_tlad%latIndexDepot(gridptIndex)
        interpWeight(gridptCount) = interpInfo_tlad%interpWeightDepot(gridptIndex)

      end do

    end do subGrid_loop

    call utl_tmg_stop(30)

  end subroutine s2c_getWeightsAndGridPointIndexes

  ! -------------------------------------------------------------
  ! s2c_deallocInterpInfo
  ! -------------------------------------------------------------
  subroutine s2c_deallocInterpInfo( inputStateVectorType )
    ! :Purpose: Deallocate interpInfo_nl/tlad object.
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: inputStateVectorType

    ! Locals:
    type(struct_interpInfo),      pointer :: interpInfo
    type(struct_interpInfoTiles), pointer :: intInfoTiles
    integer :: stepIndex, procIndex, numStep

    if (trim(mpiMode) == '2DFIELDS') then

      select case( trim(inputStateVectorType) )
        case('nl')
          interpInfo => interpInfo_nl
        case('tlad')
          interpInfo => interpInfo_tlad
        case default
          call utl_abort('s2c_deallocInterpInfo: invalid input argument' // inputStateVectorType)
      end select

      if ( .not. interpInfo%initialized ) return

      write(*,*) 's2c_deallocInterpInfo: deallocating interpInfo for inputStateVectorType=', &
           inputStateVectorType

      numStep = size(interpInfo%stepProcData,2)

      deallocate(interpInfo%interpWeightDepot)
      deallocate(interpInfo%latIndexDepot)
      deallocate(interpInfo%lonIndexDepot)
      do stepIndex = 1, numStep
        do procIndex = 1, mmpi_nprocs
          deallocate(interpInfo%stepProcData(procIndex,stepIndex)%allLat)
          deallocate(interpInfo%stepProcData(procIndex,stepIndex)%allLon)
          deallocate(interpInfo%stepProcData(procIndex,stepIndex)%allHeaderIndex)
          deallocate(interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg)
          deallocate(interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd)
          if ( interpInfo%hco%rotated ) then
            deallocate(interpInfo%stepProcData(procIndex,stepIndex)%allLonRot)
            deallocate(interpInfo%stepProcData(procIndex,stepIndex)%allLatRot)
          end if
        end do
      end do
      deallocate(interpInfo%stepProcData)
      deallocate(interpInfo%allNumHeaderUsed)
      call oti_deallocate(interpInfo%oti)

      interpInfo%initialized = .false.

    else if(trim(mpiMode) == 'TILES') then

      select case( trim(inputStateVectorType) )
        case('nl')
          intInfoTiles => interpInfoTiles_nl
        case('tlad')
          intInfoTiles => interpInfoTiles_tlad
        case default
          call utl_abort('s2c_deallocInterpInfo: invalid input argument' // inputStateVectorType)
      end select

      if ( .not. intInfoTiles%initialized ) return

      write(*,*) 's2c_deallocInterpInfo: deallocating interpInfo for inputStateVectorType=', &
           inputStateVectorType

      deallocate(intInfoTiles%varLevColFromVarLevState)

      deallocate(intInfoTiles%allObsTileMpiId)
      deallocate(intInfoTiles%allNumHeader)
      deallocate(intInfoTiles%yourObsTileMpiId)
      deallocate(intInfoTiles%yourObsSubGridIndex)
      deallocate(intInfoTiles%yourObsLat)
      deallocate(intInfoTiles%yourObsLon)

      deallocate(intInfoTiles%myInterpObsLat)
      deallocate(intInfoTiles%myInterpObsLon)
      deallocate(intInfoTiles%myInterpObsXpos_r4)
      deallocate(intInfoTiles%myInterpObsYpos_r4)
      deallocate(intInfoTiles%myInterpObsSubGridIndex)
      deallocate(intInfoTiles%myInterpObsFootprint_r4)
      deallocate(intInfoTiles%myInterpObsMpiIdSrc)
      deallocate(intInfoTiles%myInterpObsHeaderIndex)
      if (intInfoTiles%rotatedWinds) then
        deallocate(intInfoTiles%myInterpObsLatRot)
        deallocate(intInfoTiles%myInterpObsLonRot)
      end if

      deallocate(intInfoTiles%depotIndexBeg)
      deallocate(intInfoTiles%depotIndexEnd)
      deallocate(intInfoTiles%latIndexDepot)
      deallocate(intInfoTiles%lonIndexDepot)
      deallocate(intInfoTiles%interpWeightDepot)

      deallocate(intInfoTiles%myHaloLatIndex)
      deallocate(intInfoTiles%myHaloLonIndex)
      deallocate(intInfoTiles%myHaloMpiIdSrc)
      deallocate(intInfoTiles%myHaloMpiTag)
      deallocate(intInfoTiles%yourHaloLatIndex)
      deallocate(intInfoTiles%yourHaloLonIndex)
      deallocate(intInfoTiles%yourHaloMpiIdDst)
      deallocate(intInfoTiles%yourHaloMpiTag)

      call oti_deallocate(intInfoTiles%oti)

      intInfoTiles%initialized = .false.

    end if

  end subroutine s2c_deallocInterpInfo

end module stateToColumn_mod

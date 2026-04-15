
module oceanObservations_mod
  ! MODULE oceanObservations_mod (prefix='oobs' category='1. High-level functionality')
  !
  !:Purpose: storage for ocean observations related subroutines.
  !
  use midasMpi_mod
  use utilities_mod
  use runtimeInfo_mod
  use obsSpaceData_mod
  use codtyp_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use oceanMask_mod
  use timeCoord_mod
  use codePrecision_mod
  use sqliteRead_mod
  use bufr_mod
  use physicsFunctions_mod

  implicit none
  save
  private

  ! Public functions/subroutines
  public :: oobs_pseudoSST, oobs_pseudoSIC

  ! mpi topology
  integer :: myLatBeg, myLatEnd
  integer :: myLonBeg, myLonEnd
  integer :: latPerPE, latPerPEmax, lonPerPE, lonPerPEmax

  contains

  !----------------------------------------------------------------------------------------
  ! oobs_pseudoSST
  !----------------------------------------------------------------------------------------
  subroutine oobs_pseudoSST(hco, vco, iceFractionThreshold, outputSST, outputFreshWaterST, &
                            iceThinning, outputFileName, seaWaterThreshold, useSalinity)
    !
    !:Purpose: to generate pseudo SST data
    !
    implicit none

    ! Arguments:
    type(struct_hco), pointer , intent(inout) :: hco                  ! horizontal grid structure
    type(struct_vco), pointer , intent(in)    :: vco                  ! vertical grid structure
    real(8)                   , intent(in)    :: iceFractionThreshold ! consider no ice condition below this threshold
    real(8)                   , intent(in)    :: outputSST            ! output SST value for pseudo observations
    real(8)                   , intent(in)    :: outputFreshWaterST   ! output fresh water surface temperature for pseudo obs
    integer                   , intent(in)    :: iceThinning          ! generate pseudo obs in every 'iceThinning' points
    character(len=*)          , intent(in)    :: outputFileName
    real(8)                   , intent(in)    :: seaWaterThreshold    ! to distinguish inland water from sea water
    logical                   , intent(in)    :: useSalinity          ! to use or not NEMO salinity field to compute freezing point temperature

    ! Locals:
    type(struct_gsv)     :: stateVector_ice, stateVector_seaWater, stateVector_salinity
    real(8), pointer     :: seaIce_ptr(:, :, :), seaWater_ptr(:, :, :)
    real(4), pointer     :: salinity_ptr(:, :, :)
    type(struct_ocm)     :: oceanMask
    integer              :: numberIceCoveredPoints, lonIndex, latIndex, dateStamp, inlandWaterPoints
    integer              :: datePrint, timePrint, seaWaterPoints
    integer              :: randomSeed
    integer, allocatable :: iceDomainIndexesAux(:), iceDomainIndexes(:)
    real(8), allocatable :: seaWaterFractionAux(:), iceLonsAux(:), iceLatsAux(:) , salinityAux(:)
    real(8), allocatable :: seaWaterFraction(:), iceLons(:), iceLats(:), salinity(:)
    type(struct_obs)     :: obsData

    call rti_tmg_start(185,'--oobs_pseudoSST')

    ! get mpi topology
    call mmpi_setup_lonbands(hco%ni, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd)
    call mmpi_setup_latbands(hco%nj, latPerPE, latPerPEmax, myLatBeg, myLatEnd)

    ! get latest sea-ice analysis
    call gsv_allocate(stateVector_ice, 1, hco, vco, dataKind_opt = 8, &
                      datestamp_opt = -1, mpi_local_opt = .false., varNames_opt = (/'LG'/), &
                      hInterpolateDegree_opt = 'LINEAR')
    call gio_readFromFile(stateVector_ice, './seaice_analysis', ' ','A', &
                          unitConversion_opt=.false., containsFullField_opt=.true.)
    call gsv_getField(stateVector_ice, seaIce_ptr)

    ! read sea water fraction
    call gsv_allocate(stateVector_seaWater, 1, hco, vco, dataKind_opt = 8, &
                      datestamp_opt = -1, mpi_local_opt = .false., varNames_opt = (/'VF'/), &
                      hInterpolateDegree_opt = 'LINEAR')
    call gio_readFromFile(stateVector_seaWater, './seaice_analysis', ' ','A', &
                          unitConversion_opt = .false., containsFullField_opt = .true.)
    call gsv_getField(stateVector_seaWater, seaWater_ptr)

    ! Get land mask from analysisgrid file (1=water, 0=land)
    call ocm_readMaskFromFile(oceanMask, hco, vco, './analysisgrid')

    if (useSalinity) then
      ! get latest sea-ice analysis
      call gsv_allocate(stateVector_salinity, 1, hco, vco, dataKind_opt = 4, &
                        datestamp_opt = -1, mpi_local_opt = .false., varNames_opt = (/'SSS'/))
      call gio_readFromFile(stateVector_salinity, './restart', ' ', 'A', &
                            unitConversion_opt = .false., containsFullField_opt = .true.)
      call gsv_getField(stateVector_salinity, salinity_ptr)
      allocate(salinityAux((myLonEnd - myLonBeg + 1) * (myLatEnd - myLatBeg + 1)))
    end if

    allocate(iceDomainIndexesAux((myLonEnd - myLonBeg + 1) * (myLatEnd - myLatBeg + 1)))
    allocate(iceLonsAux((myLonEnd - myLonBeg + 1) * (myLatEnd - myLatBeg + 1)))
    allocate(iceLatsAux((myLonEnd - myLonBeg + 1) * (myLatEnd - myLatBeg + 1)))
    allocate(seaWaterFractionAux((myLonEnd - myLonBeg + 1) * (myLatEnd - myLatBeg + 1)))

    numberIceCoveredPoints = 0
    inlandWaterPoints = 0
    seaWaterPoints = 0

    do lonIndex = myLonBeg, myLonEnd
      do latIndex = myLatBeg, myLatEnd
        if (oceanMask%mask(lonIndex, latIndex, 1)) then
          if (seaice_ptr(lonIndex, latIndex, 1) > iceFractionThreshold) then

            numberIceCoveredPoints = numberIceCoveredPoints + 1

            iceDomainIndexesAux(numberIceCoveredPoints) = numberIceCoveredPoints
            seaWaterFractionAux(numberIceCoveredPoints) = seaWater_ptr(lonIndex, latIndex, 1)

            if(seaWater_ptr(lonIndex, latIndex, 1) <= seaWaterThreshold) then
              inlandWaterPoints = inlandWaterPoints + 1
            else
              seaWaterPoints = seaWaterPoints + 1
            end if

            iceLonsAux(numberIceCoveredPoints) = hco%lon2d_4 (lonIndex, latIndex)
            iceLatsAux(numberIceCoveredPoints) = hco%lat2d_4 (lonIndex, latIndex)
            iceDomainIndexesAux(numberIceCoveredPoints) = numberIceCoveredPoints
            seaWaterFractionAux(numberIceCoveredPoints) = seaWater_ptr(lonIndex, latIndex, 1)
            if (useSalinity) salinityAux(numberIceCoveredPoints) = salinity_ptr(lonIndex, latIndex, 1)

          end if
        end if
      end do
    end do
    call ocm_deallocate(oceanMask)
    call gsv_deallocate(stateVector_ice)
    call gsv_deallocate(stateVector_seaWater)
    if (useSalinity) call gsv_deallocate(stateVector_salinity)

    write(*,*) 'oobs_pseudoSST: ', numberIceCoveredPoints, ' ice-covered points found'
    write(*,*) 'oobs_pseudoSST: where ', inlandWaterPoints, ' are inland water points'
    write(*,*) 'oobs_pseudoSST: ', seaWaterPoints, ' sea water points found'

    if (numberIceCoveredPoints > 0) then

      allocate(iceDomainIndexes(1:numberIceCoveredPoints))
      iceDomainIndexes(:) = iceDomainIndexesAux(1:numberIceCoveredPoints)
      allocate(seaWaterFraction(1:numberIceCoveredPoints))
      seaWaterFraction(:) = seaWaterFractionAux(1:numberIceCoveredPoints)
      allocate(iceLons(1:numberIceCoveredPoints))
      allocate(iceLats(1:numberIceCoveredPoints))
      iceLons(:) = iceLonsAux(1:numberIceCoveredPoints)
      iceLats(:) = iceLatsAux(1:numberIceCoveredPoints)

      if (useSalinity) then
        allocate(salinity(1:numberIceCoveredPoints))
        salinity(:) = salinityAux(1:numberIceCoveredPoints)
      end if

    end if

    deallocate(iceLonsAux)
    deallocate(iceLatsAux)
    deallocate(iceDomainIndexesAux)
    deallocate(seaWaterFractionAux)
    if (useSalinity) deallocate(salinityAux)

    dateStamp = tim_getDatestampFromFile('./seaice_analysis', varNameForDate_opt = 'LG')
    write(*,*) 'oobs_pseudoSST: datestamp: ', dateStamp
    ! compute random seed from the date for randomly forming sea-ice subdomain
    call tim_dateStampToYYYYMMDDHH(dateStamp, datePrint, timePrint)
    timePrint = timePrint / 1000000
    datePrint =  datePrint * 100 + timePrint

    ! Remove the century, keeping 2 digits of the year
    randomSeed = datePrint - 100000000 * (datePrint / 100000000)
    write(*,*) 'oobs_pseudoSST: datePrint, timePrint: ', datePrint, timePrint

    if (numberIceCoveredPoints > 0) then

      call utl_randomOrderInt(iceDomainIndexes, randomSeed)
      write(*,*) 'oobs_pseudoSST: seed for random shuffle of sea-ice points: ', randomSeed

      if (useSalinity) then
        call oobs_computeObsData(obsData, iceDomainIndexes, iceLons, iceLats, &
                                 iceThinning, outputSST, outputFreshWaterST, &
                                 outputFileName, datePrint, timePrint, &
                                 seaWaterFraction, seaWaterThreshold, &
                                 inlandWaterPoints, salinity_opt = salinity)
      else
        call oobs_computeObsData(obsData, iceDomainIndexes, iceLons, iceLats, &
                                 iceThinning, outputSST, outputFreshWaterST, &
                                 outputFileName, datePrint, timePrint, &
                                 seaWaterFraction, seaWaterThreshold, &
                                 inlandWaterPoints)
      end if

    else

      call obs_initialize(obsData, numHeader_max_opt = 0, numBody_max_opt = 0, mpi_local_opt = .true.)
      call sqlr_writePseudoOceanIceObs(obsData, 'OS', outputFileName)

    end if

    if (numberIceCoveredPoints > 0) then
      deallocate(iceLons)
      deallocate(iceLats)
      deallocate(iceDomainIndexes)
      deallocate(seaWaterFraction)
      if (useSalinity) deallocate(salinity)
    end if

    write(*,*) 'oobs_pseudoSST: done'
    call rti_tmg_stop(185)

  end subroutine oobs_pseudoSST

  !----------------------------------------------------------------------------------------
  ! oobs_pseudoSIC
  !----------------------------------------------------------------------------------------
  subroutine oobs_pseudoSIC(hco, vco, iceFractionThreshold, outputFileName, seaIceBand)
    !
    !:Purpose: to generate pseudo Sea Ice Concentration (SIC) data to preserve sharp
    !          SIC gradients at the boundary between ice-covered regions and open water,
    !          which is essential for maintaining a realistic ice edge during
    !          strongly coupled sea-ice-ocean data assimilation.
    !
    !          The computation of the pseudo SIC observations is performed globally
    !          on the MPI task 0 (the global sea-ice field is read by the MPI task 0).
    !          Other MPI tasks produce empty pseudo SIC observation files.
    !
    implicit none

    ! Arguments
    type(struct_hco), pointer , intent(inout) :: hco
    type(struct_vco), pointer , intent(in)    :: vco
    real(8)                   , intent(in)    :: iceFractionThreshold
    character(len=*)          , intent(in)    :: outputFileName
    real(8)                   , intent(in)    :: seaIceBand

    ! Locals
    type(struct_gsv)     :: stateVector_ice
    real(8), pointer     :: seaIce_ptr(:, :, :)
    type(struct_ocm)     :: oceanMask
    logical, allocatable :: isIce(:,:), isEdge(:,:), isBand(:,:), isGap(:,:), isOcean(:,:)
    integer              :: newDate
    integer              :: ni, nj, bandCells
    integer              :: lonIndex, latIndex, nIceNeighbors
    type(struct_obs)     :: obsData
    integer              :: numberPseudoSICPoints, headerIndex
    integer              :: dateStamp, datePrint, timePrint, imode, ierr
    integer              :: codeType
    real(pre_obsReal)    :: obsLon, obsLat, obsValue
    real(8), allocatable :: iceLons(:), iceLats(:)
    real(8)              :: gridResolution   ! grid resolution in km
    ! BFS queue (Breadth-First Search)
    ! It's an algorithm used to explore faster a grid layer by layer
    integer, allocatable :: queueLonIndex(:), queueLatIndex(:), queueDistance(:)
    integer :: queueHead, queueTail

    call rti_tmg_start(186,'--oobs_pseudoSIC')

    if (mmpi_myid /= 0) then
      call obs_initialize(obsData, numHeader_max_opt = 0, numBody_max_opt = 0, mpi_local_opt = .true.)
      call sqlr_writePseudoOceanIceObs(obsData, 'GL', outputFileName)
      call rti_tmg_stop(186)
      return
    end if

    ni = hco%ni
    nj = hco%nj

    ! get latest sea-ice analysis
    call gsv_allocate(stateVector_ice, 1, hco, vco, dataKind_opt = 8, &
                      datestamp_opt = -1, mpi_local_opt = .false., varNames_opt = (/'LG'/), &
                      hInterpolateDegree_opt = 'LINEAR')
    call gio_readFromFile(stateVector_ice, './seaice_analysis', ' ','A', &
                          unitConversion_opt = .false., containsFullField_opt = .true.)
    call gsv_getField(stateVector_ice, seaIce_ptr)

    ! Get land mask from analysisgrid file (1=water, 0=land)
    call ocm_readMaskFromFile(oceanMask, hco, vco, './analysisgrid')

    !---------------------------------------------------------
    ! Allocate
    !---------------------------------------------------------
    allocate(isIce(ni,nj), isEdge(ni,nj), isBand(ni,nj), isGap(ni,nj), isOcean(ni,nj))

    isIce(:,:)  = .false.
    isEdge(:,:) = .false.
    isBand(:,:) = .false.
    isGap(:,:)  = .false.

    isOcean(:,:) = oceanMask%mask(:,:,1)

    !---------------------------------------------------------
    ! 1. Ice mask
    !---------------------------------------------------------
    do latIndex = 1, nj
      do lonIndex = 1, ni
        if (isOcean(lonIndex, latIndex)) then
          if (seaIce_ptr(lonIndex, latIndex, 1) > iceFractionThreshold) then
            isIce(lonIndex, latIndex) = .true.
          end if
        end if
      end do
    end do

    !---------------------------------------------------------
    ! 2. Edge detection
    !---------------------------------------------------------
    do latIndex = 2, nj - 1
      do lonIndex = 2, ni - 1

        if (.not. isOcean(lonIndex, latIndex)) cycle

        if (isIce(lonIndex, latIndex) .neqv. isIce(lonIndex - 1, latIndex) .or. &
            isIce(lonIndex, latIndex) .neqv. isIce(lonIndex + 1, latIndex) .or. &
            isIce(lonIndex, latIndex) .neqv. isIce(lonIndex, latIndex - 1) .or. &
            isIce(lonIndex, latIndex) .neqv. isIce(lonIndex, latIndex + 1)) then

          isEdge(lonIndex, latIndex) = .true.

        end if

      end do
    end do

    gridResolution = hco%maxGridSpacing / 1000.d0
    bandCells = max(1, nint(seaIceBand / gridResolution))
    write(*,*) 'oobs_pseudoSIC: grid resolution ', gridResolution
    write(*,*) 'oobs_pseudoSIC: number of band cells: ', bandCells

    !---------------------------------------------------------
    ! 3. BFS band expansion
    !---------------------------------------------------------

    allocate(queueLonIndex(ni*nj), queueLatIndex(ni*nj), queueDistance(ni*nj))
    queueHead = 1
    queueTail = 0

    ! Initialize queue with edges
    do latIndex = 1, nj
      do lonIndex = 1, ni
        if (isEdge(lonIndex, latIndex)) then
          queueTail = queueTail + 1
          queueLonIndex(queueTail) = lonIndex
          queueLatIndex(queueTail) = latIndex
          queueDistance(queueTail) = 0
          isBand(lonIndex, latIndex) = .true.
        end if
      end do
    end do

    ! BFS propagation
    do while (queueHead <= queueTail)

      lonIndex = queueLonIndex(queueHead)
      latIndex = queueLatIndex(queueHead)

      if (queueDistance(queueHead) < bandCells) then

        if (lonIndex > 1) then
          if (.not. isBand(lonIndex - 1, latIndex)) then
            queueTail = queueTail + 1
            queueLonIndex(queueTail) = lonIndex - 1
            queueLatIndex(queueTail) = latIndex
            queueDistance(queueTail) = queueDistance(queueHead) + 1
            isBand(lonIndex - 1, latIndex) = .true.
          end if
        end if

        if (lonIndex < ni) then
          if (.not. isBand(lonIndex + 1, latIndex)) then
            queueTail = queueTail + 1
            queueLonIndex(queueTail) = lonIndex + 1
            queueLatIndex(queueTail) = latIndex
            queueDistance(queueTail) = queueDistance(queueHead) + 1
            isBand(lonIndex + 1, latIndex) = .true.
          end if
        end if

        if (latIndex > 1) then
          if (.not. isBand(lonIndex, latIndex - 1)) then
            queueTail = queueTail + 1
            queueLonIndex(queueTail) = lonIndex
            queueLatIndex(queueTail) = latIndex - 1
            queueDistance(queueTail) = queueDistance(queueHead) + 1
            isBand(lonIndex, latIndex - 1) = .true.
          end if
        end if

        if (latIndex < nj) then
          if (.not. isBand(lonIndex, latIndex + 1)) then
            queueTail = queueTail + 1
            queueLonIndex(queueTail) = lonIndex
            queueLatIndex(queueTail) = latIndex + 1
            queueDistance(queueTail) = queueDistance(queueHead) + 1
            isBand(lonIndex, latIndex + 1) = .true.
          end if
        end if

      end if ! queueDistance(queueHead) < bandCells

      queueHead = queueHead + 1

    end do ! BFS propagation

    deallocate(queueLonIndex, queueLatIndex, queueDistance)

    !---------------------------------------------------------
    ! 4. Gap detection
    !---------------------------------------------------------
    do latIndex = 2, nj-1
      do lonIndex = 2, ni-1

        if (.not. isOcean(lonIndex, latIndex)) cycle
        if (isIce(lonIndex, latIndex)) cycle

        nIceNeighbors = 0

        if (isIce(lonIndex - 1, latIndex)) nIceNeighbors = nIceNeighbors + 1
        if (isIce(lonIndex + 1, latIndex)) nIceNeighbors = nIceNeighbors + 1
        if (isIce(lonIndex, latIndex - 1)) nIceNeighbors = nIceNeighbors + 1
        if (isIce(lonIndex, latIndex + 1)) nIceNeighbors = nIceNeighbors + 1

        if (nIceNeighbors >= 3) isGap(lonIndex, latIndex) = .true.

      end do
    end do

    !---------------------------------------------------------
    ! 5. Collect obs points
    !---------------------------------------------------------
    numberPseudoSICPoints = 0
    allocate(iceLons(ni*nj), iceLats(ni*nj))

    do latIndex = 1, nj
      do lonIndex = 1, ni

        if (.not. isOcean(lonIndex, latIndex)) cycle

        if (.not. isIce(lonIndex, latIndex)) then
          if (isBand(lonIndex, latIndex) .or. isGap(lonIndex, latIndex)) then
            numberPseudoSICPoints = numberPseudoSICPoints + 1
            iceLons(numberPseudoSICPoints) = hco%lon2d_4(lonIndex, latIndex)
            iceLats(numberPseudoSICPoints) = hco%lat2d_4(lonIndex, latIndex)
          end if
        end if

      end do
    end do

    !---------------------------------------------------------
    ! Date
    !---------------------------------------------------------
    dateStamp = tim_getDatestampFromFile('./seaice_analysis', varNameForDate_opt = 'LG')
    imode = -3
    ierr = newdate(dateStamp, datePrint, timePrint, imode)
    timePrint = timePrint / 1000000
    datePrint = datePrint * 100 + timePrint

    !---------------------------------------------------------
    ! Write obs
    !---------------------------------------------------------

    write(*,*) 'oobs_pseudoSIC: pseudo SIC points: ', numberPseudoSICPoints

    if (numberPseudoSICPoints > 0) then

      call obs_initialize(obsData, numHeader_max_opt = numberPseudoSICPoints, &
                          numBody_max_opt = numberPseudoSICPoints, mpi_local_opt = .true.)
      codeType = codtyp_get_codtyp('pseudosfc')
      obsValue = 0.0_pre_obsReal

      do headerIndex = 1, numberPseudoSICPoints

        obsLon = iceLons(headerIndex)
        obsLat = iceLats(headerIndex)

        call obs_setFamily(obsData, 'GL', headerIndex)
        call obs_headSet_i(obsData, OBS_ONM, headerIndex, headerIndex)
        call obs_headSet_i(obsData, OBS_ITY, headerIndex, codeType)
        call obs_headSet_r(obsData, OBS_LAT, headerIndex, obsLat)
        call obs_headSet_r(obsData, OBS_LON, headerIndex, obsLon)
        call obs_bodySet_r(obsData, OBS_VAR, headerIndex, obsValue)
        call obs_bodySet_i(obsData, OBS_VNM, headerIndex, bufr_iceBogus)
        call obs_set_c(obsData, 'STID', headerIndex, 'BOGUS')
        call obs_headSet_i(obsData, OBS_NLV, headerIndex, 1)
        call obs_headSet_i(obsData, OBS_RLN, headerIndex, headerIndex)
        call obs_headSet_i(obsData, OBS_DAT, headerIndex, datePrint / 100)
        call obs_headSet_i(obsData, OBS_ETM, headerIndex, timePrint)

      end do

      call sqlr_writePseudoOceanIceObs(obsData, 'GL', outputFileName)

    else

      write(*,*) 'oobs_pseudoSIC: WARNING: No pseudo SIC observations were generated.'
      write(*,*) '                Check namelist parameters: iceFractionThreshold and seaIceBand.'
      call obs_initialize(obsData, numHeader_max_opt = 0, numBody_max_opt = 0, mpi_local_opt = .true.)
      call sqlr_writePseudoOceanIceObs(obsData, 'GL', outputFileName)

    end if

    call obs_finalize(obsData)

    !---------------------------------------------------------
    ! Cleanup
    !---------------------------------------------------------
    deallocate(isIce, isEdge, isBand, isGap, isOcean)
    deallocate(iceLons, iceLats)

    call ocm_deallocate(oceanMask)
    call gsv_deallocate(stateVector_ice)
    call rti_tmg_stop(186)

  end subroutine oobs_pseudoSIC

  !--------------------------------------------------------------------------
  ! oobs_computeObsData
  !--------------------------------------------------------------------------
  subroutine oobs_computeObsData(obsData, iceDomainIndexes, iceLons, iceLats, iceThinning, &
                                 outputSST, outputFreshWaterST, outputFileName, &
                                 datePrint, timePrint, seaWaterFraction, &
                                 seaWaterThreshold, inlandWaterPoints, salinity_opt)
    !
    !:Purpose: pseudo SST data are put into obsSpaceData
    !          and written into an SQLite file

    implicit none

    ! Arguments:
    type(struct_obs) , intent(inout) :: obsData            ! obsSpaceData
    integer          , intent(in)    :: iceDomainIndexes(:)! array of the ice-covered point indexes
    real(8)          , intent(in)    :: iceLons(:)         ! longitudes of sea ice
    real(8)          , intent(in)    :: iceLats(:)         ! latitudes of sea ice
    integer          , intent(in)    :: iceThinning        ! generate pseudo obs in every 'iceThinning' points
    real(8)          , intent(in)    :: outputSST          ! output SST value for pseudo observations
    real(8)          , intent(in)    :: outputFreshWaterST ! output fresh water surface temperature for pseudo obs
    character(len=*) , intent(in)    :: outputFileName
    integer          , intent(in)    :: datePrint
    integer          , intent(in)    :: timePrint
    real(8)          , intent(in)    :: seaWaterFraction(:)! sea water fraction data: 0: fresh water; 1: sea water
    real(8)          , intent(in)    :: seaWaterThreshold  ! to distinguish inland water from sea water
    integer          , intent(in)    :: inlandWaterPoints  ! number of inland water points
    real(8), optional, intent(in)    :: salinity_opt(:)    ! to use or not NEMO salinity field to compute freezing point temperature

    ! Locals:
    real(pre_obsReal) :: obsLon, obsLat, obsValue
    integer           :: iceIndex, iceDomainDimension, pseudoObsDimension
    integer           :: codeType, headerIndex
    integer           :: coordinatesIndex, counterThinning, checkInlandWatersCount, checkSeaWatersCount

    iceDomainDimension = size(iceDomainIndexes)
    pseudoObsDimension = floor(real((iceDomainDimension - inlandWaterPoints) / iceThinning)) + inlandWaterPoints

    write(*,*) 'oobs_computeObsData: sea-ice domain dimension: ', iceDomainDimension
    write(*,*) 'oobs_computeObsData: pseudo obs vector dimension: ', pseudoObsDimension
    write(*,*) 'oobs_computeObsData: pseudo SST obs will be generated in every ', iceThinning, &
    ' points of the sea-ice field for sea water, '
    write(*,*) 'oobs_computeObsData: and in every point for inland waters, where sea water fraction <= ', seaWaterThreshold
    write(*,*) 'oobs_computeObsData: number of inland waters points: ', inlandWaterPoints

    call obs_initialize(obsData, numHeader_max_opt = pseudoObsDimension, &
                        numBody_max_opt = pseudoObsDimension, mpi_local_opt = .true.)
    codeType = codtyp_get_codtyp('pseudosfc')

    headerIndex = 1
    counterThinning = iceThinning
    checkInlandWatersCount = 0
    checkSeaWatersCount = 0

    do iceIndex = 1, iceDomainDimension

      if (headerIndex > pseudoObsDimension ) cycle

      coordinatesIndex = iceDomainIndexes(iceIndex)
      obsLon   = iceLons(coordinatesIndex)
      obsLat   = iceLats(coordinatesIndex)

      if (seaWaterFraction(coordinatesIndex) <= seaWaterThreshold) then
        obsValue = real((1.0d0 - seaWaterFraction(coordinatesIndex)) * outputFreshWaterST + &
                        seaWaterFraction(coordinatesIndex)* outputSST, pre_obsReal)
        checkInlandWatersCount =  checkInlandWatersCount + 1
      else
        if (counterThinning == iceThinning) then
          obsValue = real(outputSST, pre_obsReal)
          checkSeaWatersCount = checkSeaWatersCount + 1
          counterThinning = 1
        else
          counterThinning = counterThinning + 1
          cycle
        end if
      end if

      if (present(salinity_opt)) then
        ! observation value is set to freezing point temperature (in Kelvin),
        ! pressure is set to zero as in the current operational system
        obsValue = real(phf_getFreezingPoint(salinity_opt(coordinatesIndex), 0.0d0), pre_obsReal)
      end if

      call obs_setFamily(obsData, 'OS'   , headerIndex)
      call obs_headSet_i(obsData, OBS_ONM, headerIndex, headerIndex)
      call obs_headSet_i(obsData, OBS_ITY, headerIndex, codeType)
      call obs_headSet_r(obsData, OBS_LAT, headerIndex, obsLat)
      call obs_headSet_r(obsData, OBS_LON, headerIndex, obsLon)
      call obs_bodySet_r(obsData, OBS_VAR, headerIndex, obsValue)
      call obs_bodySet_i(obsData, OBS_VNM, headerIndex, bufr_sst)
      call     obs_set_c(obsData, 'STID' , headerIndex, 'ABOG')
      call obs_headSet_i(obsData, OBS_NLV, headerIndex, 1)
      call obs_headSet_i(obsData, OBS_RLN, headerIndex, headerIndex)
      call obs_headSet_i(obsData, OBS_DAT, headerIndex, datePrint / 100)
      call obs_headSet_i(obsData, OBS_ETM, headerIndex, timePrint)
      headerIndex = headerIndex + 1
    end do

    call sqlr_writePseudoOceanIceObs(obsData, 'OS', outputFileName)

    ! Deallocate obsSpaceData
    call obs_finalize(obsData)

  end subroutine oobs_computeObsData

end module oceanObservations_mod

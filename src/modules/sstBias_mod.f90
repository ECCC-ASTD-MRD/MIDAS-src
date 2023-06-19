
module sstBias_mod

  ! MODULE sstBias (prefix='sstb' category='1. High-level functionality')
  !
  !:Purpose: Compute SST satellite data bias estimation and correction
  !
  use obsSpaceData_mod  
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use kdTree2_mod
  use codePrecision_mod
  use mathPhysConstants_mod
  use utilities_mod
  use midasMpi_mod
  use codtyp_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use oceanMask_mod
  use localizationFunction_mod
  use columnData_mod
  use statetocolumn_mod 
   
  implicit none
  save
  private

  ! public subroutines
  public :: sstb_computeBias, sstb_applySatelliteSSTBiasCorrection
  
  ! external 
  integer, external  :: fnom, fclos
  
  ! mpi topology
  integer            :: myLatBeg, myLatEnd
  integer            :: myLonBeg, myLonEnd
  integer            :: latPerPE, latPerPEmax, lonPerPE, lonPerPEmax

  integer, parameter :: maxNumberSensors = 10

  ! namelist variables
  real(8)            :: searchRadius                 ! horizontal search radius, in km, for obs gridding
  real(4)            :: maxBias                      ! max acceptable difference (insitu - satellite)
  real(4)            :: iceFractionThreshold         ! consider no ice condition below this threshold
  integer            :: numberPointsBG               ! parameter, number of matchups of the background bias estimation
  character(len=20)  :: timeInterpType_nl            ! 'NEAREST' or 'LINEAR'
  integer            :: numObsBatches                ! number of batches for calling interp setup
  integer            :: numberSensors                ! MUST NOT BE INCLUDED IN NAMELIST!
  character(len=10)  :: sensorList(maxNumberSensors) ! list of sensors
  logical            :: saveAuxFields                ! to store or not auxiliary fields: nobs and weight        
  real(4)            :: weightMin                    ! minimum value of weight for the current day bias
  real(4)            :: weightMax                    ! maximum value of weight for the current day bias
  real(4)            :: bgTermZeroBias               ! background term of zero bias estimate

  contains

  !--------------------------------------------------------------------------
  ! sstb_computeBias
  !--------------------------------------------------------------------------
  subroutine sstb_computeBias(obsData, hco, vco)
    !
    !:Purpose: compute bias for SST satellite data with respect to insitu data 
    !  
    implicit none
    
    ! Arguments: 
    type(struct_obs), intent(inout)          :: obsData              ! obsSpaceData
    type(struct_hco), intent(inout), pointer :: hco                  ! horizontal grid structure
    type(struct_vco), intent(in)   , pointer :: vco                  ! vertical grid structure
      
    ! Locals:
    integer                     :: sensorIndex, productIndex
    real(8)                     :: insituGrid   (hco%ni, hco%nj)
    real(8)                     :: satelliteGrid(hco%ni, hco%nj)
    logical                     :: mask(hco%ni, hco%nj), openWater(hco%ni, hco%nj) 
    type(struct_ocm)            :: oceanMask
    integer                     :: numberOpenWaterPoints, lonIndex, latIndex
    integer                     :: nobsFoundInsituGlob, nobsFoundInsituLoc
    integer                     :: nobsFoundSatGlob, nobsFoundSatLoc
    type(struct_gsv)            :: stateVector_ice
    real(4), pointer            :: seaice_ptr(:,:,:)
    integer         , parameter :: numberProducts = 2  ! day and night
    character(len=*), parameter :: listProducts(numberProducts)= (/'day', 'night'/)

    write(*,*) 'sstb_computeBias: Starting...'
    write(*,*) 'sstb_computeBias: Sea-ice Fraction threshold: ', iceFractionThreshold

    ! read the namelist
    call readNml()

    ! get mpi topology
    call mmpi_setup_lonbands(hco%ni, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd)
    call mmpi_setup_latbands(hco%nj, latPerPE, latPerPEmax, myLatBeg, myLatEnd)

    ! get latest sea-ice analysis
    call gsv_allocate(stateVector_ice, 1, hco, vco, dataKind_opt = 4, &
                      datestamp_opt = -1, mpi_local_opt = .false.,    &
                      varNames_opt = (/'LG'/), hInterpolateDegree_opt ='LINEAR')
    call gio_readFromFile(stateVector_ice, './seaice_analysis', ' ','A', &
                           unitConversion_opt=.false., containsFullField_opt=.true.)
    call gsv_getField(stateVector_ice, seaice_ptr)
    
    ! Get land mask from analysisgrid file (1=water, 0=land) 
    ! and the number of open water points
    call ocm_readMaskFromFile(oceanMask, hco, vco, './analysisgrid')

    numberOpenWaterPoints = 0
    mask(:, :) = .false.
    openWater(:, :) = .false.
    do latIndex = 1, hco%nj
      do lonIndex = 1, hco%ni
        if (oceanMask%mask(lonIndex, latIndex, 1)) then
          mask(lonIndex, latIndex) = .true.
          if (seaice_ptr(lonIndex, latIndex, 1) <= iceFractionThreshold) then
            openWater(lonIndex, latIndex) = .true.
            numberOpenWaterPoints = numberOpenWaterPoints + 1
          end if
        end if
      end do
    end do
    call ocm_deallocate(oceanMask)
    call gsv_deallocate(stateVector_ice)
    write(*,*) 'sstb_computeBias: computing bias for ', numberOpenWaterPoints, ' open water points'
    
    insituGrid(:, :) = MPC_missingValue_R8

    call sstb_getGriddedObs(obsData, insituGrid, nobsFoundInsituGlob, &
                            nobsFoundInsituLoc, hco, openWater, 'insitu')

    do sensorIndex = 1, numberSensors 
      do productIndex = 1, numberProducts
        if (nobsFoundInsituGlob > 0) then
          satelliteGrid(:, :) = MPC_missingValue_R8
          call sstb_getGriddedObs(obsData, satelliteGrid(:, :), nobsFoundSatGlob, &
                                  nobsFoundSatLoc, hco, openWater, sensorList(sensorIndex), &
                                  dayOrNight_opt = listProducts(productIndex))
          if (nobsFoundSatGlob > 0) then
            call sstb_getGriddedBias(satelliteGrid(:, :), insituGrid, nobsFoundSatGlob, &
                                     nobsFoundSatLoc, hco, vco, mask, openWater, &
                                     sensorList(sensorIndex), numberOpenWaterPoints, &
                                     listProducts(productIndex))
          else
            write(*,*) 'sstb_computeBias: WARNING: missing ', trim(sensorList(sensorIndex)), ' ', &
                       trim(listProducts(productIndex)),' data.' 
            write(*,*) 'Bias estimate will be read from the previous state...'
            call sstb_getBiasFromPreviousState(hco, vco, &
                                               sensorList(sensorIndex), &
                                               listProducts(productIndex)) 
          end if
        else
          write(*,*) 'sstb_computeBias: WARNING: missing insitu data.'
          write(*,*) 'sstb_computeBias: bias estimates for all sensors will be read from the previous state...'
          call sstb_getBiasFromPreviousState(hco, vco, &
                                             sensorList(sensorIndex), &
                                             listProducts(productIndex)) 
        end if
      end do
    end do
      
    write(*,*) 'sstb_computeBias: done.'
  
  end subroutine sstb_computeBias
  
  !--------------------------------------------------------------------------
  ! sstb_getGriddedObs
  !--------------------------------------------------------------------------
  subroutine sstb_getGriddedObs(obsData, obsGrid, countObsGlob, countObsLoc, &
                                hco, openWater, instrument, dayOrNight_opt)
    !
    !:Purpose: put observations of a given family on the grid
    !           
    implicit none
    
    ! Arguments: 
    type(struct_obs), intent(inout)          :: obsData        ! obsSpaceData
    real(8)         , intent(inout)          :: obsGrid(:,:)   ! observations on the grid
    integer         , intent(out)            :: countObsGlob   ! global number of data found (all procs)
    integer         , intent(out)            :: countObsLoc    ! number of data found (current MPI proc)
    type(struct_hco), intent(in)   , pointer :: hco            ! horizontal grid structure
    logical         , intent(in)             :: openWater(:,:) ! open water points (.true.)
    character(len=*), intent(in)             :: instrument     ! name of instrument
    character(len=*), intent(in), optional   :: dayOrNight_opt ! look for daytime or nighttime obs
    
    ! Locals:
    integer, parameter        :: maxPointsSearch = 200000
    real(4), parameter        :: solarZenithThreshold = 90.0 ! to distinguish day and night
    type(kdtree2), pointer    :: tree => null() 
    real(kdkind), allocatable :: positionArray(:,:)
    type(kdtree2_result)      :: searchResults(maxPointsSearch)
    real(kdkind)              :: refPosition(3)
    real(kdkind)              :: lon_grd, lat_grd
    real(pre_obsReal)         :: lon_obs, lat_obs
    integer                   :: lonIndex, latIndex
    integer                   :: bodyIndex, headerIndex, ierr, headerCounter, codtyp
    integer                   :: localObsIndex
    integer                   :: ndataFoundGridLoc(hco%ni, hco%nj)  ! kd-tree output: number of data found within the search radius for every grid point
    integer                   :: ndataFoundGridGlob(hco%ni, hco%nj) ! to compute mpi_allreduce of ndataFoundGridLoc 
    real(kdkind)              :: searchRadiusSquared
    integer, allocatable      :: headerIndexes(:)
    real(8)                   :: currentObs
    character(len=50)         :: instrumentString 

    countObsLoc = 0
    
    ! count local observations for the given instrument
    if (trim(instrument) == 'insitu') then
      instrumentString = trim(instrument)
      do headerIndex = 1, obs_numheader(obsData)
        codtyp = obs_headElem_i(obsData, obs_ity, headerIndex)
        if (codtyp == codtyp_get_codtyp('shipnonauto') .or. &
             codtyp == codtyp_get_codtyp('drifter')     .or. &
             codtyp == codtyp_get_codtyp('ashipauto')) then
          countObsLoc = countObsLoc + 1
        end if
      end do
    else
      instrumentString = trim(instrument)//' '//dayOrNight_opt
      do headerIndex = 1, obs_numheader(obsData)
        if (obs_elem_c(obsData, 'STID' , headerIndex) == trim(instrument)) then
          if (trim(dayOrNight_opt) == 'day') then
            if (obs_headElem_r(obsData, obs_sun, headerIndex) <  solarZenithThreshold) &
              countObsLoc = countObsLoc + 1
          else if (trim(dayOrNight_opt) == 'night') then
            if (obs_headElem_r(obsData, obs_sun, headerIndex) >= solarZenithThreshold) &
              countObsLoc = countObsLoc + 1
          end if
        end if
      end do 
    end if  
    write(*,*) ''
    write(*,"(a, i10, a)") 'sstb_getGriddedObs: found ', countObsLoc, ' '//trim(instrumentString)//' data'

    call rpn_comm_allreduce(countObsLoc, countObsGlob, 1, "mpi_integer", "mpi_sum", "grid", ierr)

    obsGrid(:, :) = 0.0d0
    ndataFoundGridLoc(:,:) = 0
    ndataFoundGridGlob(:,:) = 0
    
    POSITIVECOUNTOBSLOC: if (countObsLoc > 0) then

      write(*,*)'sstb_getGriddedObs: define kd-tree using data positions...'    
      allocate(positionArray(3, countObsLoc))
      allocate(headerIndexes(countObsLoc))

      headerCounter = 0
      do headerIndex = 1, obs_numheader(obsData)
        if (trim(instrument) == 'insitu') then
          codtyp = obs_headElem_i(obsData, obs_ity, headerIndex)
	  if (codtyp == codtyp_get_codtyp('shipnonauto') .or. &
              codtyp == codtyp_get_codtyp('drifter')     .or. &
              codtyp == codtyp_get_codtyp('ashipauto')) then
            lon_obs = obs_headElem_r(obsData, obs_lon, headerIndex)
            lat_obs = obs_headElem_r(obsData, obs_lat, headerIndex)
            headerCounter = headerCounter + 1
            positionArray(:, headerCounter) = kdtree2_3dPosition(lon_obs, lat_obs)
            headerIndexes(headerCounter) = headerIndex
          end if
        else
          if (obs_elem_c(obsData, 'STID' , headerIndex) == trim(instrument)) then
            if (trim(dayOrNight_opt) == 'day') then
              if (obs_headElem_r(obsData, obs_sun, headerIndex) <  solarZenithThreshold) then
                lon_obs = obs_headElem_r(obsData, obs_lon, headerIndex)
                lat_obs = obs_headElem_r(obsData, obs_lat, headerIndex)
                headerCounter = headerCounter + 1
                positionArray(:, headerCounter) = kdtree2_3dPosition(lon_obs, lat_obs)
                headerIndexes(headerCounter) = headerIndex
              end if
            else if (trim(dayOrNight_opt) == 'night') then 
              if (obs_headElem_r(obsData, obs_sun, headerIndex) >= solarZenithThreshold) then
                lon_obs = obs_headElem_r(obsData, obs_lon, headerIndex)
                lat_obs = obs_headElem_r(obsData, obs_lat, headerIndex)
                headerCounter = headerCounter + 1
                positionArray(:, headerCounter) = kdtree2_3dPosition(lon_obs, lat_obs)
                headerIndexes(headerCounter) = headerIndex
              end if
            end if  
          end if
        end if 
      end do
    
      nullify(tree)
      tree => kdtree2_create(positionArray, sort=.true., rearrange=.true.)
      
      ! do the search
      write(*, '(a,f5.1,a)') 'sstb_getGriddedObs: Collocation radius: ', searchRadius, ' km'
      searchRadiusSquared = (1.1d0 * searchRadius * 1000.d0)**2 ! convert from km to m2

      write(*,*) 'sstb_getGriddedObs: computing the sum of data values and '//&
                 'their number within the collocation radius for every grid point...'
      do latIndex = 1, hco%nj
        do lonIndex = 1, hco%ni 

          ! compute gridded obs for every open water point
          OPENWATERPTS: if (openWater(lonIndex, latIndex)) then 
        
            lon_grd = real(hco % lon2d_4 (lonIndex, latIndex), 8)
            lat_grd = real(hco % lat2d_4 (lonIndex, latIndex), 8)
            refPosition(:) = kdtree2_3dPosition(lon_grd, lat_grd)
            call kdtree2_r_nearest(tp = tree, qv = refPosition, r2 = searchRadiusSquared, &
                                   nfound = ndataFoundGridLoc(lonIndex, latIndex), & 
                                   nalloc = maxPointsSearch, results = searchResults)
            if (ndataFoundGridLoc(lonIndex, latIndex) > maxPointsSearch) &
            call utl_abort('sstb_getGriddedObs: the parameter maxPointsSearch must be increased')

            if (ndataFoundGridLoc(lonIndex, latIndex) > 0) then
              do localObsIndex = 1, ndataFoundGridLoc(lonIndex, latIndex)
                bodyIndex  = obs_headElem_i(obsData, obs_rln, headerIndexes(searchResults(localObsIndex)%idx))
                currentObs = obs_bodyElem_r(obsData, obs_var, bodyIndex)
                obsGrid(lonIndex, latIndex) = obsGrid(lonIndex, latIndex) + currentObs 
              end do
	    end if
	    
	  end if OPENWATERPTS
        end do
      end do

      deallocate(headerIndexes)
      deallocate(positionArray)
      
    end if POSITIVECOUNTOBSLOC

    POSITIVECOUNTOBSGLOB: if (countObsGlob > 0) then
      write(*,*) 'sstb_getGriddedObs: computing average values for every grid point...'
      do latIndex = 1, hco%nj
        do lonIndex = 1, hco%ni 
          ! summing the values over all mpi tasks and sending them back to all tasks preserving the order of summation
          call mmpi_allreduce_sumreal8scalar(obsGrid(lonIndex, latIndex), "grid")
          call rpn_comm_allreduce(ndataFoundGridLoc(lonIndex, latIndex), &
                                  ndataFoundGridGlob(lonIndex, latIndex), 1, &
                                  'mpi_integer', 'mpi_sum', 'grid', ierr)
	  if (ndataFoundGridGlob(lonIndex, latIndex) > 0) then
	    obsGrid(lonIndex, latIndex) = obsGrid(lonIndex, latIndex) / &
                                          real(ndataFoundGridGlob(lonIndex, latIndex))
	  else
	    obsGrid(lonIndex, latIndex) = MPC_missingValue_R8
	  end if
          if (.not.openWater(lonIndex, latIndex)) obsGrid(lonIndex, latIndex) = MPC_missingValue_R8
        end do
      end do
      call rpn_comm_barrier('GRID', ierr)
      write(*,*) 'sstb_getGriddedObs: gridding for '//trim(instrumentString)//' data completed'
    end if POSITIVECOUNTOBSGLOB

  end subroutine sstb_getGriddedObs

  !--------------------------------------------------------------------------
  ! sstb_getGriddedBias
  !--------------------------------------------------------------------------
  subroutine sstb_getGriddedBias(satelliteGrid, insituGrid, nobsGlob, nobsLoc, &
                                 hco, vco,  mask, openWater, sensor, &
                                 numberOpenWaterPoints, dayOrNight)
    !
    !:Purpose: compute the satellite SST data bias estimation field on a grid
    !           
    implicit none
    
    ! Arguments: 
    real(8)         , intent(inout)          :: satelliteGrid(:,:)   ! gridded satellite data
    real(8)         , intent(inout)          :: insituGrid(:,:)      ! gridded insitu data
    integer         , intent(in)             :: nobsGlob             ! number of data on all procs 
    integer         , intent(in)             :: nobsLoc              ! number of data on the current MPI proc 
    type(struct_hco), intent(in)   , pointer :: hco                  ! horizontal grid structure
    type(struct_vco), intent(in)   , pointer :: vco                  ! vertical grid structure
    logical         , intent(in)             :: mask(:,:)            ! land-ocean mask
    logical         , intent(in)             :: openWater(:,:)       ! open water points (.true.)
    character(len=*), intent(in)             :: sensor               ! current sensor name
    integer         , intent(in)             :: numberOpenWaterPoints! number of open water points to allocate kd-tree work arrays 
    character(len=*), intent(in)             :: dayOrNight           ! look for daytime or nighttime obs

    ! Locals:
    real(4), parameter          :: solarZenithThreshold = 90.0       ! to distinguish day and night
    type(kdtree2), pointer      :: tree => null() 
    real(kdkind), allocatable   :: positionArray(:,:)
    integer, parameter          :: maxPointsSearch = 200000    
    type(kdtree2_result)        :: searchResults(maxPointsSearch)
    real(kdkind)                :: refPosition(3)
    real(kdkind)                :: lon_grd, lat_grd
    integer                     :: lonIndex, latIndex
    integer                     :: ierr
    integer                     :: localIndex, indexCounter, localLonIndex, localLatIndex
    integer                     :: numPointsFound
    real(kdkind)                :: searchRadiusSquared
    type(struct_gsv)            :: stateVector                      ! state vector containing bias estimation field
    type(struct_gsv)            :: stateVector_searchRadius, stateVector_previous
    type(struct_gsv)            :: stateVectorNobs, stateVectorWeight    
    real(4), pointer            :: griddedBias_r4_ptr(:,:,:), searchRadius_ptr(:,:,:)
    real(4), pointer            :: nobsField_r4_ptr(:,:,:), weightField_r4_ptr(:,:,:)
    real(4), pointer            :: griddedBias_r4_previous_ptr(:,:,:)
    integer, allocatable        :: gridPointIndexes(:,:)
    real(8)                     :: weight, distance, correlation, lengthscale, difference, numberPoints
    character(len=1)            :: extension
    character(len=*), parameter :: outputFileName = './satellite_bias.fst'
    character(len=*), parameter :: outputAuxFileName = './auxOutput.fst' ! auxiliary file name to store nobs^a and weight
    
    write(*,*) ''
    write(*,*) 'sstb_getGriddedBias: computing bias for: '//trim(sensor)//' '//trim(dayOrNight)
    write(*,*) 'sstb_getGriddedBias: the current processor contains ', nobsLoc, ' data out of ', nobsGlob  
    
    if (dayOrNight == 'day') then
      extension = 'D'
    else if (dayOrNight == 'night') then
      extension = 'N'	  
    else  
      call utl_abort('sstb_getGriddedBias: wrong extension: '//trim(extension)) 
    end if  
    
    POSITIVEOBSNUMBER: if (nobsLoc > 0) then

      write(*,*) 'sstb_getGriddedBias: compute kd-tree with all the grid points...'
      allocate(positionArray(3, numberOpenWaterPoints))
      allocate(gridPointIndexes(2, numberOpenWaterPoints))
    
      call lfn_setup('FifthOrder')

      indexCounter = 0
      do latIndex = 1, hco%nj
        do lonIndex = 1, hco%ni 
          if (openWater(lonIndex, latIndex)) then
            indexCounter = indexCounter + 1
            lon_grd = real(hco%lon2d_4(lonIndex, latIndex), 8)
            lat_grd = real(hco%lat2d_4(lonIndex, latIndex), 8)
            positionArray(:, indexCounter) = kdtree2_3dPosition(lon_grd, lat_grd)
            gridPointIndexes(1, indexCounter) = lonIndex
            gridPointIndexes(2, indexCounter) = latIndex
          end if  
        end do
      end do
    
      nullify(tree)
      tree => kdtree2_create(positionArray, sort=.true., rearrange=.true.)
 
    end if POSITIVEOBSNUMBER 

    ! get search radius field
    call gsv_allocate(stateVector_searchRadius, 1, hco, vco, dataKind_opt = 4, &
                      datestamp_opt = -1, mpi_local_opt = .true., varNames_opt = (/'TM'/))
    call gio_readFromFile(stateVector_searchRadius, './searchRadius', 'RADIUS','A', &
                          unitConversion_opt=.false., containsFullField_opt=.true.)
    call gsv_getField(stateVector_searchRadius, searchRadius_ptr)

    if (saveAuxFields) then
      ! output nobs state vector
      call gsv_allocate(stateVectorNobs, 1, hco, vco, dataKind_opt = 4, &
                        datestamp_opt = tim_getDateStamp(), mpi_local_opt = .true., &
                        varNames_opt = (/'TM'/))
      ! pointer for nobs stateVector
      call gsv_getField(stateVectorNobs, nobsField_r4_ptr)
      ! output weight state vector
      call gsv_allocate(stateVectorWeight, 1, hco, vco, dataKind_opt = 4, &
                        datestamp_opt = tim_getDateStamp(), mpi_local_opt = .true., &
                        varNames_opt = (/'TM'/))
      ! pointer for weight stateVector
      call gsv_getField(stateVectorWeight, weightField_r4_ptr)
      nobsField_r4_ptr(:,:,1) = 0.0d0
      weightField_r4_ptr(:,:,1) = 0.0d0
    end if
    
    ! previous bias estimation
    call gsv_allocate(stateVector_previous, 1, hco, vco, dataKind_opt = 4, &
                      datestamp_opt = -1, mpi_local_opt = .true., varNames_opt = (/'TM'/))
    call gio_readFromFile(stateVector_previous, './trlm_01', 'B_'//trim(sensor)//'_'//trim(extension), &
                          ' ', unitConversion_opt=.false., containsFullField_opt=.true.)
    call gsv_getField(stateVector_previous, griddedBias_r4_previous_ptr) 
       
    ! resulting bias estimation state vector
    call gsv_allocate(stateVector, 1, hco, vco, dataKind_opt = 4, &
                      datestamp_opt = tim_getDateStamp(), mpi_local_opt = .true., varNames_opt = (/'TM'/))
    ! pointer for bias estimation stateVector
    call gsv_getField(stateVector, griddedBias_r4_ptr)
    
    griddedBias_r4_ptr(:,:,1) = 0.0d0

    if (nobsLoc > 0) then
      write(*,*) 'sstb_getGriddedBias: do the search for '//trim(sensor)//' '//trim(dayOrNight)//'...'
    else
      write(*,*) 'sstb_getGriddedBias: no '//trim(sensor)//' '&
                 //trim(dayOrNight)//' data on the current processor' 
      write(*,*) 'sstb_getGriddedBias: previous estimation state will be used.'
    end if

    do latIndex = myLatBeg, myLatEnd
      do lonIndex = myLonBeg, myLonEnd 

        POSITIVEOBSNUMBER2: if (nobsLoc > 0) then
    
	  numberPoints = 0.0d0
          LANDSEAMASK: if (mask(lonIndex, latIndex)) then
            lon_grd = real(hco%lon2d_4(lonIndex, latIndex), 8)
            lat_grd = real(hco%lat2d_4(lonIndex, latIndex), 8)
            refPosition(:) = kdtree2_3dPosition(lon_grd, lat_grd)
            ! convert from km to m2
            searchRadiusSquared = (1.1d0 * searchRadius_ptr(lonIndex, latIndex, 1) * 1000.d0)**2
            call kdtree2_r_nearest(tp = tree, qv = refPosition, &
                                   r2 = searchRadiusSquared, &
                                   nfound = numPointsFound, &
                                   nalloc = maxPointsSearch, &
                                   results = searchResults)
            POSITIVENPTSFOUND: if (numPointsFound > 0) then
              do localIndex = 1, numPointsFound
                localLonIndex = gridPointIndexes(1, searchResults(localIndex)%idx)
                localLatIndex = gridPointIndexes(2, searchResults(localIndex)%idx)
                difference = satelliteGrid(localLonIndex, localLatIndex) - &
                             insituGrid(localLonIndex, localLatIndex)
                if (insituGrid   (localLonIndex, localLatIndex) /= MPC_missingValue_R8 .and. &
                    satelliteGrid(localLonIndex, localLatIndex) /= MPC_missingValue_R8 .and. &
                    abs(difference) < maxBias) then
                  distance = sqrt(searchResults(localIndex)%dis)
                  lengthscale = 1000.d0 * searchRadius_ptr(lonIndex, latIndex, 1)
                  correlation = lfn_response(distance, lengthscale)
                  griddedBias_r4_ptr(lonIndex, latIndex, 1) = griddedBias_r4_ptr(lonIndex, latIndex, 1) + &
                                                              correlation * difference 
                  numberPoints = numberPoints + correlation
                end if
              end do
    
              if (numberPoints > 0.0d0) &
              griddedBias_r4_ptr(lonIndex, latIndex, 1) = griddedBias_r4_ptr(lonIndex, latIndex, 1) / numberPoints
            end if POSITIVENPTSFOUND   
          end if LANDSEAMASK

          weight = numberPoints / (numberPoints + numberPointsBG)
          if (weight < weightMin) weight = weightMin
          if (weight > weightMax) weight = weightMax
       
	  if (saveAuxFields) then
            weightField_r4_ptr(lonIndex, latIndex, 1) = weight
	    nobsField_r4_ptr(lonIndex, latIndex, 1) = numberPoints
	  end if
	  
	  ! computation of the bias:
          griddedBias_r4_ptr(lonIndex, latIndex, 1) = (1.0d0 - weight) * bgTermZeroBias * &
                                                      griddedBias_r4_previous_ptr(lonIndex, latIndex, 1) + &
                                                      weight * griddedBias_r4_ptr(lonIndex, latIndex, 1)
        else ! no data on the current processor   
       
	  ! the bias estimation on the current processor is the estimation from previous state:
          griddedBias_r4_ptr(lonIndex, latIndex, 1) = griddedBias_r4_previous_ptr(lonIndex, latIndex, 1)

        end if POSITIVEOBSNUMBER2
      end do
    end do
    
    call rpn_comm_barrier('GRID', ierr)
    call gio_writeToFile(stateVector, outputFileName, 'B_'//trim(sensor)//'_'//trim(extension))

    if (nobsLoc > 0) then
      deallocate(gridPointIndexes)
      deallocate(positionArray)
    end if

    call gsv_deallocate(stateVector_searchRadius)
    call gsv_deallocate(stateVector_previous)
    call gsv_deallocate(stateVector)
    if (saveAuxFields) then
      call gio_writeToFile(stateVectorNobs, outputAuxFileName, 'N_'//trim(sensor)//'_'//trim(extension))
      call gio_writeToFile(stateVectorWeight, outputAuxFileName, 'W_'//trim(sensor)//'_'//trim(extension))
      call gsv_deallocate(stateVectorNobs)
      call gsv_deallocate(stateVectorWeight)
    end if
    
    write(*,*) 'sstb_getGriddedBias: completed for: '//trim(sensor)//' '//trim(dayOrNight)

  end subroutine sstb_getGriddedBias

  !--------------------------------------------------------------------------
  ! sstb_getBiasCorrection
  !--------------------------------------------------------------------------
  subroutine sstb_getBiasCorrection(stateVector, column, obsData, hco, sensor, dayOrNight)
    !
    !:Purpose: To compute bias correction and put it into obsSpace data. 
    !          Columns from input field are interpolated to obs location
    !     
    implicit none
    
    ! Arguments:
    type(struct_gsv)       , intent(inout)       :: stateVector       ! state vector containing bias estimation field    
    type(struct_columnData), intent(inout)       :: column            ! column data
    type(struct_obs)       , intent(inout)       :: obsData           ! obsSpaceData
    type(struct_hco)       , intent(in), pointer :: hco               ! horizontal grid
    character(len=*)       , intent(in)          :: sensor            ! current sensor name
    character(len=*)       , intent(in)          :: dayOrNight        ! look for daytime or nighttime obs

    ! Locals:
    real(4), parameter :: solarZenithThreshold = 90.0 ! to distinguish day and night
    integer            :: bodyIndex, headerIndex
    real(8)            :: currentObs

    write(*,*)
    write(*,*) 'sstb_getBiasCorrection: computing bias correction for ', trim(sensor), ' ', trim(dayOrNight), 'time ************'

    call s2c_nl(stateVector, obsData, column, hco, timeInterpType = timeInterpType_nl, &
                moveObsAtPole_opt = .true., numObsBatches_opt = numObsBatches, dealloc_opt = .false.)

    do headerIndex = 1, obs_numheader(obsData)
      
      if (obs_elem_c(obsData, 'STID' , headerIndex) == trim(sensor)) then
        
        bodyIndex  = obs_headElem_i(obsData, obs_rln, headerIndex)
	currentObs = obs_bodyElem_r(obsData, obs_var, bodyIndex)
	
	if (trim(dayOrNight) == 'day') then
          if (obs_headElem_r(obsData, obs_sun, headerIndex) <  solarZenithThreshold) then
            call obs_bodySet_r(obsData, obs_bcor, bodyIndex, col_getElem(column, 1, headerIndex, 'TM'))
            call obs_bodySet_r(obsData, obs_var , bodyIndex, currentObs - col_getElem(column, 1, headerIndex, 'TM'))
	  end if  
	else if (trim(dayOrNight) == 'night') then
          if (obs_headElem_r(obsData, obs_sun, headerIndex) >= solarZenithThreshold) then
            call obs_bodySet_r(obsData, obs_bcor, bodyIndex, col_getElem(column, 1, headerIndex, 'TM'))
            call obs_bodySet_r(obsData, obs_var , bodyIndex, currentObs - col_getElem(column, 1, headerIndex, 'TM'))
          end if
	end if
      end if
      
    end do 
 
    write(*,*) 'sstb_getBiasCorrection: END'

  end subroutine sstb_getBiasCorrection

  !--------------------------------------------------------------------------
  ! readNml (private subroutine)
  !--------------------------------------------------------------------------
  subroutine readNml()
    !
    !:Purpose: Read the namelist `namSSTbiasEstimate`
    !
    implicit none

    ! Locals:
    integer :: ierr, nulnam, sensorIndex
    
    namelist /namSSTbiasEstimate/ searchRadius, maxBias, iceFractionThreshold, numberPointsBG, &
                                  timeInterpType_nl, numObsBatches, numberSensors, sensorList, &
                                  weightMin, weightMax, saveAuxFields, bgTermZeroBias

    ! Setting default namelist variable values
    searchRadius = 10.            
    maxBias = 1.                  
    iceFractionThreshold   = 0.6 
    numberSensors = MPC_missingValue_INT
    numberPointsBG = 0            
    timeInterpType_nl = 'NEAREST'
    numObsBatches = 20
    sensorList(:) = ''
    weightMin = 0.0
    weightMax = 1.0
    saveAuxFields = .False.
    bgTermZeroBias = 1.0
    
    ! Read the namelist
    nulnam = 0
    ierr = fnom( nulnam, './flnml', 'FTN+SEQ+R/O', 0 )
    read(nulnam, nml = namSSTbiasEstimate, iostat = ierr )
    if (ierr /= 0) call utl_abort('readNml (sstb): Error reading namelist')
    if (mmpi_myid == 0) write(*, nml = namSSTbiasEstimate )
    ierr = fclos( nulnam )

    if (numberSensors /= MPC_missingValue_INT) then
      call utl_abort('readNml (sstb): check namSSTbiasEstimate namelist section: numberSensors should be removed')
    end if

    numberSensors = 0
    do sensorIndex = 1, maxNumberSensors
      if (trim(sensorList(sensorIndex))=='') exit
      numberSensors = numberSensors + 1
    end do
    if (numberSensors == 0) call utl_abort('readNml (sstb): check namSSTbiasEstimate namelist section: empty sensorList')

    if (mmpi_myid == 0) then
      write(*,*)
      write(*,*) 'readNml (sstb): sensors to treat: '
      do sensorIndex = 1, numberSensors
        write(*,*) 'readNml (sstb): sensor index: ', sensorIndex, ', sensor: ', sensorList( sensorIndex )
      end do
      write(*,*) 'readNml (sstb): interpolation type: ', timeInterpType_nl
      write(*,*) 'readNml (sstb): number obs batches: ', numObsBatches
      write(*,*) 'readNml (sstb): weight limits for current bias estimate: ', weightMin, weightMax
    end if

  end subroutine readNml

  !--------------------------------------------------------------------------
  ! sstb_applySatelliteSSTBiasCorrection
  !--------------------------------------------------------------------------
  subroutine sstb_applySatelliteSSTBiasCorrection(obsData, hco, vco, column)
    !
    !:Purpose: To apply bias satellite SST data bias correction and put it into obsSpace data. 
    !          Columns from input field are interpolated to obs location
    !     
    implicit none
    
    ! Arguments:
    type(struct_obs)       , intent(inout)          :: obsData ! obsSpaceData
    type(struct_hco)       , intent(inout), pointer :: hco     ! horizontal grid structure
    type(struct_vco)       , intent(in)   , pointer :: vco     ! vertical grid structure
    type(struct_columnData), intent(inout)          :: column  ! column data 

    ! Locals:
    type(struct_gsv)            :: stateVector
    integer         , parameter :: numberProducts = 2            ! day and night
    character(len=*), parameter :: listProducts(numberProducts)= (/'day', 'night'/) ! day and night biases
    integer                     :: sensorIndex, productIndex
    character(len=1)            :: extension
    character(len=*), parameter :: biasFileName = './satellite_bias.fst'

    ! read the namelist
    call readNml()

    ! allocate state vector for bias estimation field
    call gsv_allocate(stateVector, 1, hco, vco, dataKind_opt = 4, hInterpolateDegree_opt = 'LINEAR', &
                      datestamp_opt = -1, mpi_local_opt = .true., varNames_opt = (/'TM'/))

    do sensorIndex = 1, numberSensors 
      do productIndex = 1, numberProducts
      
        if (trim(listProducts(productIndex)) == 'day') then
          extension = 'D'
        else if (trim(listProducts(productIndex)) == 'night') then
          extension = 'N'
        else
          call utl_abort('sstb_applySatelliteSSTBiasCorrection: wrong extension: '//trim(extension)) 
        end if
	
        call gio_readFromFile(stateVector, biasFileName, 'B_'//trim(sensorList(sensorIndex))//'_'//trim(extension), &
                              ' ', unitConversion_opt = .false., containsFullField_opt = .true.)
        call sstb_getBiasCorrection(stateVector, column, obsData, hco, trim(sensorList(sensorIndex)), &
                                    trim(listProducts(productIndex)))
      end do
    end do
    				    
    call gsv_deallocate(stateVector)
			    
  end subroutine sstb_applySatelliteSSTBiasCorrection

  !--------------------------------------------------------------------------
  ! sstb_getBiasFromPreviousState
  !--------------------------------------------------------------------------
  subroutine sstb_getBiasFromPreviousState(hco, vco, sensor, dayOrNight)
    !
    !:Purpose: to get a satellite SST data bias estimate from the previous state if data is missing.
    !          or there are no insitu data for the current dateStamp,
    !          hence, unable to compute bias estimates for the satellite data.
    !
    implicit none
    
    ! Arguments:
    type(struct_hco), intent(inout), pointer :: hco        ! horizontal grid structure
    type(struct_vco), intent(in)   , pointer :: vco        ! vertical grid structure
    character(len=*), intent(in)             :: sensor     ! sensor name
    character(len=*), intent(in)             :: dayOrNight ! look for daytime or nighttime bias estimation

    ! Locals:
    type(struct_gsv)            :: stateVector                  ! state vector containing current  bias estimation field
    type(struct_gsv)            :: stateVector_previous         ! state vector containing previous bias estimation field
    real(4), pointer            :: griddedBias_r4_ptr(:, :, :)
    real(4), pointer            :: griddedBias_r4_previous_ptr(:, :, :)
    character(len=1)            :: extension
    character(len=*), parameter :: outputFileName = './satellite_bias.fst'

    write(*,*) ''
    write(*,*) 'sstb_getBiasFromPreviousState: for ', trim(sensor), ' ', trim(dayOrNight), ' data...'

    if (trim(dayOrNight) == 'day') then
      extension = 'D'
    else if (trim(dayOrNight) == 'night') then
      extension = 'N'	  
    else  
      call utl_abort('sstb_getBiasFromPreviousState: wrong extension: '//trim(extension))
    end if

    ! read previous bias estimation
    call gsv_allocate(stateVector_previous, 1, hco, vco, dataKind_opt = 4, &
                       datestamp_opt = -1, mpi_local_opt = .true., varNames_opt = (/'TM'/))
    call gio_readFromFile(stateVector_previous, './trlm_01', 'B_'//trim(sensor)//'_'//trim(extension), &
                          ' ', unitConversion_opt=.false., containsFullField_opt=.true.)
    call gsv_getField(stateVector_previous, griddedBias_r4_previous_ptr) 

    ! resulting bias estimation state vector
    call gsv_allocate(stateVector, 1, hco, vco, dataKind_opt = 4, &
                      datestamp_opt = tim_getDateStamp(), mpi_local_opt = .true., varNames_opt = (/'TM'/))
    ! pointer for bias estimation stateVector
    call gsv_getField(stateVector, griddedBias_r4_ptr)
        
    griddedBias_r4_ptr(:, :, :) = griddedBias_r4_previous_ptr(:, :, :)
    call gio_writeToFile(stateVector, outputFileName, 'B_'//trim(sensor)//'_'//trim(extension))

    call gsv_deallocate(stateVector)
    call gsv_deallocate(stateVector_previous)

    write(*,*) 'sstb_getBiasFromPreviousState: done ', trim(sensor), ' ', trim(dayOrNight), ' data.'

  end subroutine sstb_getBiasFromPreviousState
  
end module SSTbias_mod

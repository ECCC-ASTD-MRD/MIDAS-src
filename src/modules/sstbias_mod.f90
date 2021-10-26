!--------------------------------------- LICENCE BEGIN -----------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------

module SSTbias_mod

  ! MODULE SSTbias (prefix='sstb' category='1. High-level functionality')
  !
  ! :Purpose: Compute SST bias estimation and correction
  !
  use obsSpaceData_mod  
  use horizontalCoord_mod
  use verticalCoord_mod
  use kdtree2_mod
  use earthconstants_mod
  use codePrecision_mod
  use mathPhysConstants_mod
  use utilities_mod
  use mpi_mod
  use codtyp_mod
  use mpivar_mod
  use gridStateVector_mod
  use oceanMask_mod
  use timeCoord_mod
  use localizationFunction_mod
  
  implicit none
  save
  private

  ! public subroutines
  public :: sstb_computeBias
  
  ! mpi topology
  integer           :: myLatBeg, myLatEnd
  integer           :: myLonBeg, myLonEnd
  integer           :: latPerPE, latPerPEmax, lonPerPE, lonPerPEmax

  contains

  !--------------------------------------------------------------------------
  ! sstb_computeBias
  !--------------------------------------------------------------------------
  subroutine sstb_computeBias( obsData, hco, vco, iceFractionThreshold, searchRadius, &
                               numberSensors, sensorList, dateStamp, maxBias, numberPointsBG )
    !
    ! :Purpose: compute bias for SST satellite data with respect to insitu data 
    !  
    implicit none
    
    ! Arguments: 
    type(struct_obs), intent(inout)          :: obsData              ! obsSpaceData
    type(struct_hco), intent(inout), pointer :: hco                  ! horizontal grid structure
    type(struct_vco), intent(in)   , pointer :: vco                  ! vertical grid structure
    real(4)         , intent(in)             :: iceFractionThreshold ! for ice fraction below it, consider open water      
    real(8)         , intent(in)             :: searchRadius         ! horizontal search radius for obs gridding
    integer         , intent(in)             :: numberSensors        ! Current satellites: AMSR2, METO-B, METO-A, NOAA19, NPP
    integer         , intent(in)             :: dateStamp            ! date to put into output fst file
    character(len=*), intent(in)             :: sensorList(:)        ! list of satellite names
    real(4)         , intent(in)             :: maxBias              ! max insitu - satellite difference in degrees  
    integer         , intent(in)             :: numberPointsBG       ! namelist parameter: number of points used to compute
                                                                     ! the background state
    ! locals
    character(len=*), parameter :: myName = 'sstb_computeBias'
    integer                     :: headerIndex, sensorIndex, productIndex
    real(8)                     :: insituGrid   ( hco % ni, hco % nj )
    real(8)                     :: satelliteGrid( hco % ni, hco % nj )
    logical                     :: mask( hco % ni, hco % nj ), openWater( hco % ni, hco % nj ) 
    type(struct_ocm)            :: oceanMask
    integer                     :: numberOpenWaterPoints, lonIndex, latIndex, ierr
    type(struct_gsv)            :: stateVector_ice
    real(4), pointer            :: seaice_ptr( :, :, : )
    integer         , parameter :: numberProducts = 2  ! day and night
    character(len=*), parameter :: listProducts( numberProducts )= (/ 'day', 'night' /)

    write(*,*) 'Starting '//myName//'...'
    write(*,*) myName//': Current analysis date: ', dateStamp
    write(*,*) myName//': Sea-ice Fraction threshold: ', iceFractionThreshold

    ! get mpi topology
    call mpivar_setup_lonbands( hco % ni, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd )
    call mpivar_setup_latbands( hco % nj, latPerPE, latPerPEmax, myLatBeg, myLatEnd )

    ! get latest sea-ice analysis
    call gsv_allocate( stateVector_ice, 1, hco, vco, dataKind_opt = 4, &
                       datestamp_opt = -1, mpi_local_opt = .false.,    &
                       varNames_opt = (/'LG'/) )
    call gsv_readFromFile( stateVector_ice, './seaice_analysis', ' ','A', &
                           unitConversion_opt=.false., containsFullField_opt=.true. )
    call gsv_getField( stateVector_ice, seaice_ptr )
    
    ! Get land mask from analysisgrid file ( 1=water, 0=land ) 
    ! and the number of open water points
    call ocm_readMaskFromFile( oceanMask, hco, vco, './analysisgrid' )

    numberOpenWaterPoints = 0
    mask( :, : ) = .false.
    openWater( :, : ) = .false.
    do lonIndex = 1, hco % ni
      do latIndex = 1, hco % nj
        if ( oceanMask%mask( lonIndex, latIndex, 1 ) ) then
          mask( lonIndex, latIndex ) = .true.
          if ( seaice_ptr( lonIndex, latIndex, 1 ) <= iceFractionThreshold ) then
            openWater( lonIndex, latIndex ) = .true.
            numberOpenWaterPoints = numberOpenWaterPoints + 1
          end if
        end if
      end do
    end do
    call ocm_deallocate( oceanMask )
    call gsv_deallocate( stateVector_ice )
    write(*,*) myName//': computing bias for ', numberOpenWaterPoints, ' open water points'
    
    insituGrid( :, : ) = MPC_missingValue_R8

    call sstb_getGriddedObs( obsData, insituGrid, hco, vco, searchRadius, openWater, 'insitu' )

    do sensorIndex = 1, numberSensors 
      do productIndex = 1, numberProducts
        satelliteGrid( :, : ) = MPC_missingValue_R8
        call sstb_getGriddedObs( obsData, satelliteGrid ( :, : ), hco, vco, searchRadius, &
                                 openWater, trim( sensorList( sensorIndex )), &
                                 dayOrNight_opt = trim( listProducts(productIndex)) )
        call sstb_getGriddedBias( satelliteGrid ( :, : ), insituGrid, hco, vco, mask, openWater, &
                                  maxBias, trim( sensorList( sensorIndex )), numberOpenWaterPoints, &
                                  numberPointsBG, dateStamp, trim( listProducts(productIndex)) )
      end do
    end do
    
  end subroutine sstb_computeBias
  
  !--------------------------------------------------------------------------
  ! sstb_getGriddedObs
  !--------------------------------------------------------------------------
  subroutine sstb_getGriddedObs( obsData, obsGrid, hco, vco, searchRadius, &
                                 openWater, instrument, dayOrNight_opt )
    !
    ! :Purpose: put observations of a given family on the grid
    !           
    implicit none
    
    ! Arguments: 
    type(struct_obs), intent(inout)          :: obsData        ! obsSpaceData
    real(8)         , intent(inout)          :: obsGrid(:,:)   ! observations on the grid
    type(struct_hco), intent(in)   , pointer :: hco            ! horizontal grid structure
    type(struct_vco), intent(in)   , pointer :: vco            ! vertical grid structure
    real(8)         , intent(in)             :: searchRadius   ! horizontal search radius where to search obs
    logical         , intent(in)             :: openWater(:,:) ! open water points (.true.)
    character(len=*), intent(in)             :: instrument     ! name of instrument
    character(len=*), intent(in), optional   :: dayOrNight_opt ! look for daytime or nighttime obs
    
    ! locals
    integer, parameter          :: maxPointsSearch = 200000
    real, parameter             :: solarZenithThreshold = 90.0      ! to distinguish day and night
    type(kdtree2), pointer      :: tree => null() 
    real(kdkind), allocatable   :: positionArray(:,:)
    type(kdtree2_result)        :: searchResults( maxPointsSearch )
    real(kdkind)                :: refPosition(3)
    real(kdkind)                :: lon_grd, lat_grd
    real(pre_obsReal)           :: lon_obs, lat_obs
    integer                     :: lonIndex, latIndex
    integer                     :: bodyIndex, headerIndex, ierr, headerCounter, codtyp
    integer                     :: localObsIndex, countObs
    integer                     :: numObsFoundLoc, numObsFoundGlob
    real(kdkind)                :: searchRadiusSquared
    integer, allocatable        :: headerIndexes(:)
    real(8)                     :: currentObs
    type(struct_gsv)            :: stateVector
    real(4), pointer            :: obsGrid_r4_ptr( :, :, : )
    character(len=50)           :: instrumentString 
    character(len=*), parameter :: myName = 'sstb_getGriddedObs'

    countObs = 0
    
    if (instrument == 'insitu') then
      
      instrumentString = instrument
      do headerIndex = 1, obs_numheader( obsData )
        codtyp = obs_headElem_i( obsData, obs_ity, headerIndex )
        if ( codtyp == codtyp_get_codtyp( 'shipnonauto' ) .or. &
             codtyp == codtyp_get_codtyp( 'drifter' )     .or. &
             codtyp == codtyp_get_codtyp( 'ashipauto' )        ) then
          countObs = countObs + 1
        end if
      end do

    else
      
      instrumentString = instrument//' '//dayOrNight_opt
      do headerIndex = 1, obs_numheader( obsData )
        if ( obs_elem_c( obsData, 'STID' , headerIndex ) == trim(instrument) ) then
          if ( trim( dayOrNight_opt ) == 'day'   ) then
            if ( obs_headElem_r( obsData, obs_sun, headerIndex ) <  solarZenithThreshold ) &
              countObs = countObs + 1
          else if ( trim( dayOrNight_opt ) == 'night' ) then
            if ( obs_headElem_r( obsData, obs_sun, headerIndex ) >= solarZenithThreshold ) &
              countObs = countObs + 1
          end if
        end if
      end do 
   
    end if  
    
    write(*,*) ''
    write(*,"(a, i10, a)") myName//': found ', countObs, ' '//trim(instrumentString)//' data'
    
    if ( countObs > 0 ) then
    
      obsGrid( :, : ) = 0.0d0
      allocate( positionArray( 3, countObs ))
      allocate( headerIndexes( countObs ))
    
      headerCounter = 0
      do headerIndex = 1, obs_numheader( obsData )
        if (instrument == 'insitu') then
          codtyp = obs_headElem_i( obsData, obs_ity, headerIndex )
	  if ( codtyp == codtyp_get_codtyp( 'shipnonauto' ) .or. &
               codtyp == codtyp_get_codtyp( 'drifter' )     .or. &
               codtyp == codtyp_get_codtyp( 'ashipauto' )        ) then

            lon_obs = obs_headElem_r( obsData, obs_lon, headerIndex )
            lat_obs = obs_headElem_r( obsData, obs_lat, headerIndex )
            headerCounter = headerCounter + 1
            positionArray( :, headerCounter ) = kdtree2_3dPosition( lon_obs, lat_obs )
            headerIndexes( headerCounter ) = headerIndex
          end if
        else
          if ( obs_elem_c( obsData, 'STID' , headerIndex ) == trim(instrument) ) then
            if ( trim( dayOrNight_opt ) == 'day'   ) then
              if ( obs_headElem_r( obsData, obs_sun, headerIndex ) <  solarZenithThreshold ) then
                lon_obs = obs_headElem_r( obsData, obs_lon, headerIndex )
                lat_obs = obs_headElem_r( obsData, obs_lat, headerIndex )
                headerCounter = headerCounter + 1
                positionArray( :, headerCounter ) = kdtree2_3dPosition( lon_obs, lat_obs )
                headerIndexes( headerCounter ) = headerIndex
              end if
            else if ( trim( dayOrNight_opt ) == 'night'   ) then 
              if ( obs_headElem_r( obsData, obs_sun, headerIndex ) >= solarZenithThreshold ) then
                lon_obs = obs_headElem_r( obsData, obs_lon, headerIndex )
                lat_obs = obs_headElem_r( obsData, obs_lat, headerIndex )
                headerCounter = headerCounter + 1
                positionArray( :, headerCounter ) = kdtree2_3dPosition( lon_obs, lat_obs )
                headerIndexes( headerCounter ) = headerIndex
              end if
            end if  
          end if
        end if 
      end do
    
      nullify(tree)
      tree => kdtree2_create( positionArray, sort=.true., rearrange=.true. )
      
      ! do the search
      write(*, "(a, f5.1, a)" ) myName//': Collocation radius: ', searchRadius, ' km'
      searchRadiusSquared = ( 1.1d0 * searchRadius * 1000.d0 )**2 ! convert from km to m2

      do lonIndex = 1, hco % ni 
        do latIndex = 1, hco % nj

          ! compute gridded obs for every open water point
          if ( openWater( lonIndex, latIndex ) == .true. ) then 
        
            lon_grd = real( hco % lon2d_4 ( lonIndex, latIndex ), 8 )
            lat_grd = real( hco % lat2d_4 ( lonIndex, latIndex ), 8 )
            refPosition(:) = kdtree2_3dPosition( lon_grd, lat_grd )
            call kdtree2_r_nearest( tp = tree, qv = refPosition, r2 = searchRadiusSquared, &
                                    nfound = numObsFoundLoc, & 
                                    nalloc = maxPointsSearch, results = searchResults )
            if ( numObsFoundLoc > maxPointsSearch ) &
            call utl_abort( myName//': the parameter maxPointsSearch must be increased' )

            if ( numObsFoundLoc > 0 ) then
	    
              do localObsIndex = 1, numObsFoundLoc
                bodyIndex  = obs_headElem_i( obsData, obs_rln, headerIndexes( searchResults( localObsIndex ) % idx ))
                currentObs = obs_bodyElem_r( obsData, obs_var, bodyIndex )
                obsGrid( lonIndex, latIndex ) = obsGrid( lonIndex, latIndex ) + currentObs 
              end do
            
	    else  
            
	      obsGrid( lonIndex, latIndex ) = 0.0d0  
            
	    end if
  
            ! summing the values over all mpi tasks and sending them back to all tasks preserving the order of summation
            call mpi_allreduce_sumreal8scalar( obsGrid( lonIndex, latIndex ), "grid" )
            call rpn_comm_allreduce( numObsFoundLoc, numObsFoundGlob, 1, "mpi_integer", "mpi_sum", "grid", ierr )
  
            if ( numObsFoundGlob > 0 ) then
              obsGrid( lonIndex, latIndex ) = obsGrid( lonIndex, latIndex ) / real( numObsFoundGlob )
            else
              obsGrid( lonIndex, latIndex ) = MPC_missingValue_R8
            end if

          else

            obsGrid( lonIndex, latIndex ) = MPC_missingValue_R8
  
          end if   

        end do
      end do

      deallocate( headerIndexes )
      deallocate( positionArray )
    
    end if  
    
    call rpn_comm_barrier( 'GRID', ierr )
    write(*,*) myName//': gridding for '//trim(instrumentString)//' data completed'

  end subroutine sstb_getGriddedObs

  !--------------------------------------------------------------------------
  ! sstb_getGriddedBias
  !--------------------------------------------------------------------------
  subroutine sstb_getGriddedBias( satelliteGrid, insituGrid, hco, vco,  mask, openWater, &
                                  maxBias, sensor, numberOpenWaterPoints, &
                                  numberPointsBG, dateStamp, dayOrNight )
    !
    ! :Purpose: compute the satellite SST data bias estimation field on a grid
    !           
    implicit none
    
    ! Arguments: 
    real(8)         , intent(inout)          :: satelliteGrid(:,:)   ! gridded satellite data
    real(8)         , intent(inout)          :: insituGrid(:,:)      ! gridded insitu data
    type(struct_hco), intent(in)   , pointer :: hco                  ! horizontal grid structure
    type(struct_vco), intent(in)   , pointer :: vco                  ! vertical grid structure
    logical         , intent(in)             :: mask( :, : )         ! land-ocean mask
    logical         , intent(in)             :: openWater( :, : )    ! open water points (.true.)
    real(4)         , intent(in)             :: maxBias              ! maximum difference in degrees between satellite and insitu SST
    character(len=*), intent(in)             :: sensor               ! current sensor
    integer         , intent(in)             :: numberOpenWaterPoints! number of open water points to allocate kd-tree work arrays 
    integer         , intent(in)             :: numberPointsBG       ! namelist parameter: number of points 
                                                                     ! used to compute the previous (background) bias estimation
    integer         , intent(in)             :: dateStamp            ! date to put into output fst files
    character(len=*), intent(in)             :: dayOrNight           ! look for daytime or nighttime obs
    
    ! locals
    real, parameter             :: solarZenithThreshold = 90.0      ! to distinguish day and night
    type(kdtree2), pointer      :: tree => null() 
    real(kdkind), allocatable   :: positionArray(:,:)
    integer, parameter          :: maxPointsSearch = 200000    
    type(kdtree2_result)        :: searchResults( maxPointsSearch )
    real(kdkind)                :: refPosition(3)
    real(kdkind)                :: lon_grd, lat_grd
    real(pre_obsReal)           :: lon_obs, lat_obs
    integer                     :: lonIndex, latIndex
    integer                     :: bodyIndex, headerIndex, ierr, headerCounter, codtyp
    integer                     :: localIndex, indexCounter, localLonIndex, localLatIndex
    integer                     :: numPointsFound
    real(kdkind)                :: searchRadiusSquared
    integer, allocatable        :: headerIndexes(:)
    type(struct_gsv)            :: stateVector, stateVector_searchRadius, stateVector_previous
    real(4), pointer            :: griddedBias_r4_ptr( :, :, : ), searchRadius_ptr( :, :, : )
    real(4), pointer            :: griddedBias_r4_previous_ptr( :, :, : )
    integer, allocatable        :: gridPointIndexes(:,:)
    real(8)                     :: weight, distance, correlation, lengthscale, difference, numberPoints
    character(len=1)            :: extension
    character(len=*), parameter :: outputFileName = './satellite_bias.fst'
    character(len=*), parameter :: myName = 'sstb_getGriddedBias'
    
    write(*,*) ''
    write(*,*) myName//' computing bias for: '//trim(sensor)//' '//trim(dayOrNight)
    
    if ( dayOrNight == 'day' ) then
      extension = 'D'
    else if ( dayOrNight == 'night' ) then
      extension = 'N'
    end if  
    
    allocate( positionArray( 3, numberOpenWaterPoints ))
    allocate( gridPointIndexes( 2, numberOpenWaterPoints ))

    call lfn_setup('FifthOrder')
    
    indexCounter = 0
    do lonIndex = 1, hco % ni 
      do latIndex = 1, hco % nj
    
        if ( openWater( lonIndex, latIndex ) == .true. ) then

          indexCounter = indexCounter + 1
          lon_grd = real( hco % lon2d_4 ( lonIndex, latIndex ), 8 )
          lat_grd = real( hco % lat2d_4 ( lonIndex, latIndex ), 8 )
          positionArray( :, indexCounter ) = kdtree2_3dPosition( lon_grd, lat_grd )
          gridPointIndexes( 1, indexCounter ) = lonIndex
          gridPointIndexes( 2, indexCounter ) = latIndex
    
        end if  

      end do
    end do
    
    nullify(tree)
    tree => kdtree2_create( positionArray, sort=.true., rearrange=.true. ) 
    
    ! get search radius field
    call gsv_allocate( stateVector_searchRadius, 1, hco, vco, dataKind_opt = 4, &
                       datestamp_opt = -1, mpi_local_opt = .true., varNames_opt = (/'TM'/) )
    call gsv_readFromFile( stateVector_searchRadius, './searchRadius', 'RADIUS','A', &
                           unitConversion_opt=.false., containsFullField_opt=.true. )
    call gsv_getField( stateVector_searchRadius, searchRadius_ptr )    

    ! previous bias estimation
    call gsv_allocate( stateVector_previous, 1, hco, vco, dataKind_opt = 4, &
                       datestamp_opt = -1, mpi_local_opt = .true., varNames_opt = (/'TM'/) )
    call gsv_readFromFile( stateVector_previous, './trlm_01', 'B_'//sensor//'_'//extension, &
                           'R', unitConversion_opt=.false., containsFullField_opt=.true. )
    call gsv_getField( stateVector_previous, griddedBias_r4_previous_ptr ) 
       
    ! allocate state vector for bias estimation field
    call gsv_allocate( stateVector, 1, hco, vco, dataKind_opt = 4, &
                       datestamp_opt = dateStamp, mpi_local_opt = .true., varNames_opt = (/'TM'/) )
    call gsv_getField( stateVector, griddedBias_r4_ptr )

    ! do the search
    write(*,*) myName//': do the search for '//trim(sensor)//' '//trim(dayOrNight)//'...' 

    do lonIndex = myLonBeg, myLonEnd 
      do latIndex = myLatBeg, myLatEnd
    
        griddedBias_r4_ptr( lonIndex, latIndex, 1 ) = 0.0d0
        numberPoints = 0.0d0
    
        if ( mask( lonIndex, latIndex ) == .true. ) then
    
          lon_grd = real( hco % lon2d_4 ( lonIndex, latIndex ), 8 )
          lat_grd = real( hco % lat2d_4 ( lonIndex, latIndex ), 8 )
          refPosition(:) = kdtree2_3dPosition( lon_grd, lat_grd )
          searchRadiusSquared = ( 1.1d0 * searchRadius_ptr( lonIndex, latIndex, 1 ) * 1000.d0 )**2 ! convert from km to m2
          call kdtree2_r_nearest( tp = tree, qv = refPosition, r2 = searchRadiusSquared, &
                                  nfound = numPointsFound, nalloc = maxPointsSearch, results = searchResults )
          if ( numPointsFound > 0 ) then
            do localIndex = 1, numPointsFound
              localLonIndex = gridPointIndexes( 1, searchResults( localIndex ) % idx )
              localLatIndex = gridPointIndexes( 2, searchResults( localIndex ) % idx )
              difference = satelliteGrid( localLonIndex, localLatIndex ) - &
                           insituGrid( localLonIndex, localLatIndex )
              if ( insituGrid   ( localLonIndex, localLatIndex ) /= MPC_missingValue_R8 .and. &
                   satelliteGrid( localLonIndex, localLatIndex ) /= MPC_missingValue_R8 .and. &
                   abs( difference ) < maxBias ) then
                distance = sqrt( searchResults( localIndex ) % dis )
                lengthscale = 1000.d0 * searchRadius_ptr( lonIndex, latIndex, 1 )
                correlation = lfn_response( distance, lengthscale )
                griddedBias_r4_ptr( lonIndex, latIndex, 1 ) = griddedBias_r4_ptr( lonIndex, latIndex, 1 ) + &
                                                              correlation * difference 
                numberPoints = numberPoints + correlation
              end if
            end do
    
            if ( numberPoints > 0.0d0 ) &
            griddedBias_r4_ptr( lonIndex, latIndex, 1 ) = griddedBias_r4_ptr( lonIndex, latIndex, 1 ) / numberPoints
          end if    
  
        end if

        weight = numberPoints / ( numberPoints + numberPointsBG )
        griddedBias_r4_ptr( lonIndex, latIndex, 1 ) = ( 1.0d0 - weight )  * griddedBias_r4_previous_ptr( lonIndex, latIndex, 1 ) + &
                                                      weight              * griddedBias_r4_ptr( lonIndex, latIndex, 1 )

      end do
    end do
    
    call rpn_comm_barrier( 'GRID', ierr )
    call gsv_writeToFile( stateVector, outputFileName, 'B_'//sensor//'_'//extension )
    
    deallocate( gridPointIndexes )
    deallocate( positionArray )
    call gsv_deallocate( stateVector )
    call gsv_deallocate( stateVector_searchRadius )
    call gsv_deallocate( stateVector_previous )
    
    write(*,*) myName//' completed for: '//trim(sensor)//' '//trim(dayOrNight)

  end subroutine sstb_getGriddedBias

end module SSTbias_mod

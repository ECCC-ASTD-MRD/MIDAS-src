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

  ! MODULE SSTbiasEstimation_mod (prefix='sstb')
  !
  ! :Purpose: Estimate SST bias for satellite data
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
  use mpivar_mod
  use gridStateVector_mod
  use oceanMask_mod
  use timeCoord_mod
  use localizationFunction_mod
  
  implicit none
  save
  private

  ! Public Subroutines
  public :: sstb_computeBias
  integer, external :: get_max_rss
  integer           :: myLatBeg, myLatEnd
  integer           :: myLonBeg, myLonEnd
  integer           :: latPerPE, latPerPEmax, lonPerPE, lonPerPEmax

  contains
  
  subroutine sstb_computeBias( obsData, hco, vco, iceFractionThreshold, searchRadius, &
                               numberSensors, sensorList, dateStamp, &
                               maxBias, numberPointsBG )
    !
    ! :Purpose: compute bias for SST satellite data with respect to insitu data 
    !  
    implicit none
    
    ! Arguments: 
    type(struct_obs)                         :: obsData
    type(struct_hco), intent(inout), pointer :: hco                  ! horizontal grid structure
    type(struct_vco), intent(in)   , pointer :: vco                  ! vertical grid structure
    real(4)         , intent(in)             :: iceFractionThreshold ! for ice fraction below it, consider open water      
    real(8)         , intent(in)             :: searchRadius         ! horizontal search radius for obs gridding
    integer         , intent(in)             :: numberSensors        ! Current satellites: AMSR2, METO-B, METO-A, NOAA19, NPP
    integer         , intent(in)             :: dateStamp            ! date to put into output fst files
    character(len=*), intent(in)             :: sensorList(:)        ! list of satellite names
    real(4)         , intent(in)             :: maxBias              ! max insitu - satellite difference in degrees  
    integer         , intent(in)             :: numberPointsBG       ! namelist parameter: number of points used to compute
                                                                     ! the background state
    ! locals
    character(len=*), parameter :: myName = 'sstb_computeBias'
    integer                     :: headerIndex, sensorIndex
    real(8)                     :: insituGrid    ( hco % ni, hco % nj )
    real(8)                     :: satelliteGrid ( hco % ni, hco % nj, numberSensors, 2 )
    integer, allocatable        :: mask( :, : ) 
    type(struct_ocm)            :: oceanMask
    integer                     :: numberOpenWaterPoints, lonIndex, latIndex, ierr
    type(struct_gsv)            :: stateVector_ice
    real(4), pointer            :: seaice_ptr( :, :, : )
    
    write(*,*) 'Starting '//myName//'...'
    write(*,*) myName//': Current analysis date: ', dateStamp
    write(*,*) myName//': satellites to treat: '
    do sensorIndex = 1, numberSensors
      write(*,*) myName//': satellite index: ', sensorIndex, ', satellite: ', sensorList( sensorIndex )
    end do
    write(*,*) myName//': Sea-ice Fraction threshold: ', iceFractionThreshold

    ! get mpi topology
    call mpivar_setup_lonbands( hco % ni, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd )
    call mpivar_setup_latbands( hco % nj, latPerPE, latPerPEmax, myLatBeg, myLatEnd )

    ! get latest sea-ice analysis
    call gsv_allocate( stateVector_ice, 1, hco, vco, dataKind_opt = 4, &
                       datestamp_opt = -1, mpi_local_opt = .false.,    &
                       varNames_opt = (/'LG'/) )
    call gsv_readFromFile( stateVector_ice, './seaice_analysis', 'G6_1_2_1N','A', &
                           unitConversion_opt=.false., containsFullField_opt=.true. )
    call gsv_getField( stateVector_ice, seaice_ptr )
    
    ! Get land mask from analysisgrid file ( 1=water, 0=land ) 
    ! and the number of open water points
    allocate( mask( hco % ni, hco % nj) )
    write(*,*) myName//': reading ocean-land mask...'
    call ocm_readMaskFromFile( oceanMask, hco, vco, './analysisgrid' )
    write(*,*) myName//': ocean-land mask ok'

    numberOpenWaterPoints = 0
    do lonIndex = 1, hco % ni
      do latIndex = 1, hco % nj
        if ( oceanMask%mask ( lonIndex, latIndex, 1 ) ) then
          mask ( lonIndex, latIndex ) = 1
          if ( seaice_ptr ( lonIndex, latIndex, 1 ) <= iceFractionThreshold ) then
            numberOpenWaterPoints = numberOpenWaterPoints + 1
          end if 
        else
          mask ( lonIndex, latIndex ) = 0
        end if
      end do
    end do
    call ocm_deallocate( oceanMask )
    write(*,*) myName//': computing bias for ', numberOpenWaterPoints, ' open water points'
    
    insituGrid( :, : ) = MPC_missingValue_R8
    satelliteGrid( :, : , : , : ) = MPC_missingValue_R8

    call sstb_getGriddedObs( obsData, insituGrid, hco, vco, searchRadius, &
                             mask, seaice_ptr, iceFractionThreshold, 'insitu', dateStamp )

    do sensorIndex = 1, numberSensors 
    

      call sstb_getGriddedObs( obsData, satelliteGrid ( :, :, sensorIndex, 1 ), hco, vco, &
                               searchRadius, mask, seaice_ptr, iceFractionThreshold, &
                               trim( sensorList( sensorIndex )), dateStamp, dayOrNight = 'day' )

      call sstb_getGriddedObs( obsData, satelliteGrid ( :, :, sensorIndex, 2 ), hco, vco, &
                               searchRadius, mask, seaice_ptr, iceFractionThreshold, &
                               trim( sensorList( sensorIndex )), dateStamp, dayOrNight = 'night' )

      call sstb_getGriddedBias( satelliteGrid ( :, :, sensorIndex, 1 ), insituGrid, hco, vco, &
                                seaice_ptr, iceFractionThreshold, mask, maxBias, trim( sensorList( sensorIndex )), &
                                numberOpenWaterPoints, numberPointsBG, dateStamp, dayOrNight = 'day' )

      call sstb_getGriddedBias( satelliteGrid ( :, :, sensorIndex, 2 ), insituGrid, hco, vco, &
                                seaice_ptr, iceFractionThreshold, mask, maxBias, trim( sensorList( sensorIndex )), &
                                numberOpenWaterPoints, numberPointsBG, dateStamp, dayOrNight = 'night' )
    
    end do
    
    deallocate( mask )
    call gsv_deallocate( stateVector_ice )
    
  end subroutine sstb_computeBias
  

  subroutine sstb_getGriddedObs( obsData, obsGrid, hco, vco, searchRadius, mask, seaice_ptr, &
                                 iceFractionThreshold, instrument, dateStamp, dayOrNight )
    !
    ! :Purpose: put observations of a given family on the grid
    !           
    implicit none
    
    ! Arguments: 
    type(struct_obs), intent(inout)          :: obsData              ! observation structure
    real(8)         , intent(inout)          :: obsGrid(:,:)         ! observations on the grid
    type(struct_hco), intent(inout), pointer :: hco                  ! horizontal grid structure
    type(struct_vco), intent(in)   , pointer :: vco                  ! vertical grid structure
    real(8)         , intent(in)             :: searchRadius         ! horizontal search radius where to search obs
    integer         , intent(in)             :: mask( :, : )         ! ocean mask
    real(4)         , intent(in)   , pointer :: seaice_ptr( :, :, : )! sea-ice fraction
    real(4)         , intent(in)             :: iceFractionThreshold ! for ice fraction below it, open water      
    character(len=*), intent(in)             :: instrument           ! name of instrument
    integer         , intent(in)             :: dateStamp            ! date to put into output fst files
    character(len=*), intent(in), optional   :: dayOrNight           ! look for day or night obs
    
    ! locals
    integer, parameter          :: maxObsPointsSearch = 200000
    real, parameter             :: solarZenithThreshold = 90.0      ! to distinguish day and night
    type(kdtree2), pointer      :: tree => null() 
    real(kdkind), allocatable   :: positionArray(:,:)
    type(kdtree2_result)        :: searchResults( maxObsPointsSearch )
    real(kdkind)                :: refPosition(3)
    real(kdkind)                :: lon_grd, lat_grd
    real(pre_obsReal)           :: lon_obs, lat_obs
    integer                     :: lonIndex, latIndex
    integer                     :: bodyIndex, headerIndex, ierr, headerCounter, codtyp
    integer                     :: localObsIndex, countObs
    integer                     :: numObsFound, numObsFoundMPIGlobal
    real(kdkind)                :: searchRadiusSquared
    integer, allocatable        :: headerIndexes(:)
    real(8)                     :: currentObs
    type(struct_gsv)            :: stateVector
    real(4), pointer            :: obsGrid_r4_ptr( :, :, : )
    character(len=50)           :: instrumentString 
    character(len=1)            :: extention
    character(len=*), parameter :: myName = 'sstb_getGriddedObs'

    countObs = 0
    
    if (instrument == 'insitu') then
      
      instrumentString = instrument

      do headerIndex = 1, obs_numheader( obsData )
    
        codtyp = obs_headElem_i( obsData, obs_ity, headerIndex )
        if ( codtyp == 13 .or. codtyp == 18 .or. codtyp == 147 ) then
          countObs = countObs + 1
        end if
      
      end do

    else
      
      instrumentString = instrument//' '//dayOrNight
      
      do headerIndex = 1, obs_numheader( obsData )
      
        if ( obs_elem_c( obsData, 'STID' , headerIndex ) == trim(instrument) ) then

          if ( trim(dayOrNight) == 'day'   ) then
            if ( obs_headElem_r( obsData, obs_sun, headerIndex ) <  solarZenithThreshold ) &
              countObs = countObs + 1
            extention = 'D'
          else if ( trim(dayOrNight) == 'night' ) then
            if ( obs_headElem_r( obsData, obs_sun, headerIndex ) >= solarZenithThreshold ) &
              countObs = countObs + 1
            extention = 'N'
          end if
        end if

      end do 
   
    end if  
    
    write(*,*) ''
    write(*,*) myName//': found ', countObs, ' ', instrumentString, ' data'
    
    if ( countObs > 0 ) then
    
      obsGrid( :, : ) = 0.d0
      allocate( positionArray( 3, countObs ))
      allocate( headerIndexes( countObs ))
    
      headerCounter = 0
      do headerIndex = 1, obs_numheader( obsData )

        if (instrument == 'insitu') then

          codtyp = obs_headElem_i( obsData, obs_ity, headerIndex )

          if ( codtyp == 13 .or. codtyp == 18 .or. codtyp == 147 ) then
          
            lon_obs = obs_headElem_r( obsData, obs_lon, headerIndex )
            lat_obs = obs_headElem_r( obsData, obs_lat, headerIndex )
            headerCounter = headerCounter + 1
            positionArray( :, headerCounter ) = kdtree2_3dPosition( lon_obs, lat_obs )
            headerIndexes( headerCounter ) = headerIndex
  
          end if

        else
      
          if ( obs_elem_c( obsData, 'STID' , headerIndex ) == trim(instrument) ) then

            if ( trim(dayOrNight) == 'day'   ) then
            
              if ( obs_headElem_r( obsData, obs_sun, headerIndex ) <  solarZenithThreshold ) then
   
                lon_obs = obs_headElem_r( obsData, obs_lon, headerIndex )
                lat_obs = obs_headElem_r( obsData, obs_lat, headerIndex )
                headerCounter = headerCounter + 1
                positionArray( :, headerCounter ) = kdtree2_3dPosition( lon_obs, lat_obs )
                headerIndexes( headerCounter ) = headerIndex
           
              end if
    
            else if ( trim(dayOrNight) == 'night'   ) then 
   
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
      write(*,*) myName//': Collocation radius: ', searchRadius
      searchRadiusSquared = ( 1.1d0 * searchRadius * 1000.d0 )**2 ! convert from km to m2

      do lonIndex = 1, hco % ni 
        do latIndex = 1, hco % nj

          ! compute gridded obs for every open water point
          if ( mask( lonIndex, latIndex ) == 1 .and. &
               seaice_ptr( lonIndex, latIndex, 1 ) <= iceFractionThreshold ) then 
        
            lon_grd = real( hco % lon2d_4 ( lonIndex, latIndex ), 8 )
            lat_grd = real( hco % lat2d_4 ( lonIndex, latIndex ), 8 )
            refPosition(:) = kdtree2_3dPosition( lon_grd, lat_grd )

            call kdtree2_r_nearest( tp = tree, qv = refPosition, r2 = searchRadiusSquared, &
                                    nfound = numObsFound, & 
                                    nalloc = maxObsPointsSearch, results = searchResults )

            if ( numObsFound > maxObsPointsSearch ) &
              call utl_abort( myName//': the parameter maxObsPointsSearch must be increased' )

            if ( numObsFound > 0 ) then

              do localObsIndex = 1, numObsFound
                bodyIndex  = obs_headElem_i( obsData, obs_rln, headerIndexes( searchResults( localObsIndex ) % idx ))
                currentObs = obs_bodyElem_r( obsData, obs_var, bodyIndex )
                obsGrid( lonIndex, latIndex ) = obsGrid( lonIndex, latIndex ) + currentObs 
              end do
    
            else  

              obsGrid( lonIndex, latIndex ) = 0.0d0  

            end if
  
            ! summing the values over all mpi tasks and sending them back to all tasks preserving the order of summation
            call mpi_allreduce_sumreal8scalar( obsGrid( lonIndex, latIndex ), "grid" )
            call rpn_comm_allreduce( numObsFound, numObsFoundMPIGlobal, 1, "mpi_integer", "mpi_sum", "grid", ierr )
  
            if ( numObsFoundMPIGlobal > 0 ) then
              obsGrid( lonIndex, latIndex ) = obsGrid( lonIndex, latIndex ) / real( numObsFoundMPIGlobal )
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

  end subroutine sstb_getGriddedObs
  
  
  subroutine sstb_getGriddedBias( satelliteGrid, insituGrid, hco, vco, seaice_ptr, iceFractionThreshold, &
                                  mask, maxBias, instrument, numberOpenWaterPoints, &
                                  numberPointsBG, dateStamp, dayOrNight )
    !
    ! :Purpose: compute the satellite SST data bias estimation field on a grid
    !           
    implicit none
    
    ! Arguments: 
    real(8)         , intent(inout)          :: satelliteGrid(:,:)   ! satellite SST put on the grid
    real(8)         , intent(inout)          :: insituGrid(:,:)      ! insitu SST put on the grid
    type(struct_hco), intent(inout), pointer :: hco                  ! horizontal grid structure
    type(struct_vco), intent(in)   , pointer :: vco                  ! vertical grid structure
    real(4)         , intent(in)   , pointer :: seaice_ptr( :, :, : )! sea-ice fraction
    real(4)         , intent(in)             :: iceFractionThreshold ! for ice fraction below it, open water      
    integer         , intent(in)             :: mask( :, : )         ! ocean mask
    real(4)         , intent(in)             :: maxBias              ! maximum difference in degrees between satellite and insitu SST
    character(len=*), intent(in)             :: instrument           ! name of instrument
    integer         , intent(in)             :: numberOpenWaterPoints! number of open water points to allocate kd-tree work arrays 
    integer         , intent(in)             :: numberPointsBG       ! namelist parameter: number of points 
                                                                     ! used to compute the previous (backgrouind) bias estimation
    integer         , intent(in)             :: dateStamp            ! date to put into output fst files
    character(len=*), intent(in)             :: dayOrNight           ! look for day or night obs
    
    ! locals
    integer, parameter          :: maxObsPointsSearch = 200000
    real, parameter             :: solarZenithThreshold = 90.0      ! to distinguish day and night
    type(kdtree2), pointer      :: tree => null() 
    real(kdkind), allocatable   :: positionArray(:,:)
    integer, parameter          :: maxPointsSearch = 200000    
    type(kdtree2_result)        :: searchResults( maxObsPointsSearch )
    real(kdkind)                :: refPosition(3)
    real(kdkind)                :: lon_grd, lat_grd
    real(pre_obsReal)           :: lon_obs, lat_obs
    integer                     :: lonIndex, latIndex
    integer                     :: bodyIndex, headerIndex, ierr, headerCounter, codtyp, countFoundPointsInsitu
    integer                     :: localIndex, indexCounter, localLonIndex, localLatIndex
    integer                     :: numPointsFound
    real(kdkind)                :: searchRadiusSquared
    integer, allocatable        :: headerIndexes(:)
    character(len=*), parameter :: myName = 'sstb_getGriddedBias'
    type(struct_gsv)            :: stateVector, stateVector_npoints, stateVector_searchRadius, stateVector_previous
    real(4), pointer            :: griddedBias_r4_ptr( :, :, : ), npoints_ptr(:,:,:), searchRadius_ptr( :, :, : )
    real(4), pointer            :: griddedBias_r4_previous_ptr( :, :, : )
    integer, allocatable        :: gridPointIndexes(:,:)
    real(8)                     :: weight, distance, correlation, lengthscale, correlationSum, difference
    character(len=1)            :: extention
    character(len=*), parameter :: outputFileName = './satellite_bias.fst'
    
    write(*,*) ''
    write(*,*) myName//' computing bias for: ', instrument, ' ', dayOrNight
    
    if ( dayOrNight == 'day' ) then
      extention = 'D'
    else if ( dayOrNight == 'night' ) then
      extention = 'N'
    end if  
    
    allocate( positionArray( 3, numberOpenWaterPoints ))
    allocate( gridPointIndexes( 2, numberOpenWaterPoints ))

    call lfn_setup('FifthOrder')
    
    indexCounter = 0
    do lonIndex = 1, hco % ni 
      do latIndex = 1, hco % nj
    
        if (   mask( lonIndex, latIndex ) == 1 .and. & 
             seaice_ptr( lonIndex, latIndex, 1 ) <= iceFractionThreshold ) then

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
    call gsv_readFromFile( stateVector_previous, './trlm_01', 'B_'//instrument//'_'//extention, &
                           'R', unitConversion_opt=.false., containsFullField_opt=.true. )
    call gsv_getField( stateVector_previous, griddedBias_r4_previous_ptr ) 
       
    ! allocate state vector for bias estimation and number of points field
    call gsv_allocate( stateVector, 1, hco, vco, dataKind_opt = 4, &
                       datestamp_opt = dateStamp, mpi_local_opt = .true., varNames_opt = (/'TM'/) )
    call gsv_getField( stateVector, griddedBias_r4_ptr )

    call gsv_allocate( stateVector_npoints, 1, hco, vco, dataKind_opt = 4, &
                       datestamp_opt = dateStamp, mpi_local_opt = .true., varNames_opt = (/'TM'/) )
    call gsv_getField( stateVector_npoints, npoints_ptr )

    ! do the search
    write(*,*) myName//': do the search for ', instrument, ' ', dayOrNight,'...' 

    do lonIndex = myLonBeg, myLonEnd 
      do latIndex = myLatBeg, myLatEnd
    
        griddedBias_r4_ptr( lonIndex, latIndex, 1 ) = 0.0d0
        npoints_ptr( lonIndex, latIndex, 1 ) = 0.0d0
    
        if ( mask( lonIndex, latIndex ) == 1 ) then
    
          lon_grd = real( hco % lon2d_4 ( lonIndex, latIndex ), 8 )
          lat_grd = real( hco % lat2d_4 ( lonIndex, latIndex ), 8 )
          refPosition(:) = kdtree2_3dPosition( lon_grd, lat_grd )
  
          searchRadiusSquared = ( 1.1d0 * searchRadius_ptr( lonIndex, latIndex, 1 ) * 1000.d0 )**2 ! convert from km to m2
          
          call kdtree2_r_nearest( tp = tree, qv = refPosition, r2 = searchRadiusSquared, &
                                  nfound = numPointsFound, nalloc = maxPointsSearch, results = searchResults )

          if ( numPointsFound > 0 ) then
    
            correlationSum = 0.0d0
    
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
                npoints_ptr( lonIndex, latIndex, 1 ) = npoints_ptr( lonIndex, latIndex, 1 ) + correlation
                correlationSum = correlationSum + correlation

              end if
    
            end do
    
            if ( correlationSum > 0.0d0 ) then
              griddedBias_r4_ptr( lonIndex, latIndex, 1 ) = griddedBias_r4_ptr( lonIndex, latIndex, 1 ) / correlationSum
            end if    
    
          end if    
  
        end if

        weight = npoints_ptr( lonIndex, latIndex, 1 ) / ( npoints_ptr( lonIndex, latIndex, 1 ) + numberPointsBG )

        griddedBias_r4_ptr( lonIndex, latIndex, 1 ) = ( 1.0d0 - weight )  * griddedBias_r4_previous_ptr( lonIndex, latIndex, 1 ) + &
                                                      weight              * griddedBias_r4_ptr( lonIndex, latIndex, 1 )

      end do
    end do
    
    call rpn_comm_barrier( 'GRID', ierr )
    write(*,*) myName//': saving results...'
    call gsv_writeToFile( stateVector, outputFileName, 'B_'//instrument//'_'//extention )
    
    deallocate( gridPointIndexes )
    deallocate( positionArray )
    call gsv_deallocate( stateVector )
    call gsv_deallocate( stateVector_npoints )
    call gsv_deallocate( stateVector_searchRadius )
    call gsv_deallocate(stateVector_previous)
    
    write(*,*) myName//' completed for: ', instrument, ' ', dayOrNight

  end subroutine sstb_getGriddedBias

end module SSTbias_mod

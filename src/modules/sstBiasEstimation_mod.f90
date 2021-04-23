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

module SSTbiasEstimation_mod

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
  use gridStateVector_mod
  
  implicit none
  save
  private

  ! Public Subroutines
  public :: sstb_computeGriddedObservations
  integer, external :: get_max_rss

  contains
  
  subroutine sstb_computeGriddedObservations( obsData, hco, vco, horizontalSearchRadius, numberSatellites, satelliteList, dateStamp )

    implicit none
    
    ! Arguments: 
    type(struct_obs), intent(in)          :: obsData                ! satellite observations
    type(struct_hco), intent(in), pointer :: hco                    ! horizontal grid structure
    type(struct_vco), intent(in), pointer :: vco                    ! vertical grid structure
    real(8)         , intent(in)          :: horizontalSearchRadius ! horizontal search radius where to search observations
    integer         , intent(in)          :: numberSatellites       ! Current satellites: AMSR2, METO-B, METO-A, NOAA19, NPP
    integer         , intent(in)          :: dateStamp              ! date to put into output fst files
    character(len=*), intent(in)          :: satelliteList(:)       ! list of satellite names

    ! locals
    character(len=*), parameter :: myName = 'sstb_computeGriddedObservations'
    integer                     :: headerIndex, satelliteIndex
    
    write(*,*) 'Starting '//myName//'...'
    
    
    do satelliteIndex = 1, numberSatellites 
        
      write(*,*) myName//': treating satellite: ', satelliteIndex, trim( satelliteList( satelliteIndex ))
      call sstb_getMeanObs( obsData, hco, vco, horizontalSearchRadius, trim( satelliteList( satelliteIndex )), 'day'  , dateStamp )
      call sstb_getMeanObs( obsData, hco, vco, horizontalSearchRadius, trim( satelliteList( satelliteIndex )), 'night', dateStamp )
      
    end do	   
    
  end subroutine sstb_computeGriddedObservations
  

  subroutine sstb_getMeanObs( obsData, hco, vco, horizontalSearchRadius, instrument, dayOrNight, dateStamp )
    !
    ! :Purpose: for every horizontal point of the input horizontal grid,
    !           compute the mean of input SST observations inside a given horizontal search radius
    !           and to save the output in a standard file.
    !           The mean is computed for a given instrument and day / night time. 
    !
    implicit none
    
    ! Arguments: 
    type(struct_obs), intent(in)          :: obsData                ! satellite observations
    type(struct_hco), intent(in), pointer :: hco                    ! horizontal grid structure
    type(struct_vco), intent(in), pointer :: vco                    ! vertical grid structure
    real(8)         , intent(in)          :: horizontalSearchRadius ! horizontal search radius where to search obs
    character(len=*), intent(in)          :: instrument             ! name of instrument
    character(len=*), intent(in)          :: dayOrNight             ! look for day or night obs
    integer         , intent(in)          :: dateStamp              ! date to put into output fst files

    ! locals
    integer, parameter          :: maxObsPointsSearch = 200000
    real, parameter             :: solarZenithThreshold = 90.0      ! to distinguish day and night
    type(kdtree2), pointer      :: tree => null() 
    real(kdkind), allocatable   :: positionArray(:,:)
    type(kdtree2_result)        :: searchResults( maxObsPointsSearch )
    real(kdkind)                :: maxRadius
    real(kdkind)                :: refPosition(3)
    real(pre_obsReal)           :: lat_obs, lon_obs
    integer                     :: bodyIndex, headerIndex, ierr, countObs, headerCounter
    integer                     :: lonIndex, latIndex, localObsIndex
    real(kdkind)                :: lon_grd, lat_grd
    integer                     :: numObsFound, numObsFoundMPIGlobal
    real(kdkind)                :: searchRadiusSquared
    integer, allocatable        :: headerIndexes(:)
    type(struct_gsv)            :: stateVector
    real(8), pointer            :: meanObs_ptr( :, :, : )
    character(len=*), parameter :: myName = 'sstb_getSelectedObs'
      
    write(*,*) myName//': ##### computing mean '//trim(instrument)//' observations for ', trim(dayOrNight), ' time... #####'

    countObs = 0 
    do headerIndex = 1, obs_numheader( obsData )
      
      if ( obs_elem_c( obsData, 'STID' , headerIndex ) == trim(instrument) ) then
        
	if ( trim(dayOrNight) == 'day'   ) then
          if ( obs_headElem_r( obsData, obs_sun, headerIndex ) <  solarZenithThreshold ) countObs = countObs + 1
        else if ( trim(dayOrNight) == 'night' ) then
          if ( obs_headElem_r( obsData, obs_sun, headerIndex ) >= solarZenithThreshold ) countObs = countObs + 1
	end if  
	
      end if	
		
    end do
    
    call gsv_allocate( stateVector, 1, hco, vco, dataKind_opt=8, &
                       datestamp_opt = dateStamp, mpi_local_opt=.false., &
		       hInterpolateDegree_opt='LINEAR', varNames_opt=(/'TM'/) )
    call gsv_getField( stateVector, meanObs_ptr )
    meanObs_ptr = 0.0d0

    if ( trim(dayOrNight) == 'day'   ) then
      write(*,*) myName//': found ', countObs, ' day observations'
    else if ( trim(dayOrNight) == 'night' ) then
      write(*,*) myName//': found ', countObs, ' night observations'
    end if	
      	
    allocate( positionArray( 3, countObs ))
    allocate( headerIndexes( countObs ))

    headerCounter = 0
    do headerIndex = 1, obs_numheader( obsData )
      
      if ( obs_elem_c( obsData, 'STID' , headerIndex ) == trim(instrument) ) then
      
        lon_obs = obs_headElem_r( obsData, obs_lon, headerIndex )
        lat_obs = obs_headElem_r( obsData, obs_lat, headerIndex )
	
        if ( trim(dayOrNight) == 'day' ) then
	
          if ( obs_headElem_r( obsData, obs_sun, headerIndex ) < solarZenithThreshold ) then
	  
	    headerCounter = headerCounter + 1
            positionArray( :, headerCounter ) = kdtree2_3dPosition( lon_obs, lat_obs )
	    headerIndexes( headerCounter ) = headerIndex
	  
	  end if  
	    
        else if( trim(dayOrNight) == 'night' ) then
	
          if ( obs_headElem_r( obsData, obs_sun, headerIndex ) >= solarZenithThreshold ) then
	  
	    headerCounter = headerCounter + 1
            positionArray(:, headerCounter ) = kdtree2_3dPosition( lon_obs, lat_obs )
	    headerIndexes( headerCounter ) = headerIndex
	    
	  end if
	  
        end if
	
      end if
      	    
    end do
    
    nullify(tree)
    tree => kdtree2_create( positionArray, sort=.true., rearrange=.true. ) 
    write(*,*) 'Memory Used: ', get_max_rss() / 1024, 'Mb'
    
    ! do the search
    write(*,*) myName//': horizontal search radius, in km : ', horizontalSearchRadius
    searchRadiusSquared = ( 1.1d0 * horizontalSearchRadius * 1000.d0 )**2 ! convert from km to m2
    
    do lonIndex = 1, hco % ni
      do latIndex = 1, hco % nj
        
	lon_grd = real( hco % lon2d_4 ( lonIndex, latIndex ), 8 )
	lat_grd = real( hco % lat2d_4 ( lonIndex, latIndex ), 8 )
        refPosition(:) = kdtree2_3dPosition( lon_grd, lat_grd )
	
        call kdtree2_r_nearest( tp = tree, qv = refPosition, r2 = searchRadiusSquared, nfound = numObsFound, & 
				nalloc = maxObsPointsSearch, results = searchResults )
	
	if ( numObsFound > maxObsPointsSearch ) &
          call utl_abort( myName//': the parameter maxObsPointsSearch must be increased' )
	
        if ( numObsFound > 0 ) then
	
	  do localObsIndex = 1, numObsFound
	  
            bodyIndex  = obs_headElem_i( obsData, obs_rln, headerIndexes( searchResults( localObsIndex ) % idx ))
            meanObs_ptr( lonIndex, latIndex, 1 ) = meanObs_ptr( lonIndex, latIndex, 1 ) + &
	                                           obs_bodyElem_r( obsData, obs_var, bodyIndex ) 
	    
	  end do
	  
	end if 
	
	! summing the values over all mpi tasks and sending them back to all tasks preserving the order of summation
        call mpi_allreduce_sumreal8scalar( meanObs_ptr( lonIndex, latIndex, 1 ), "grid" )
	! doing the same for numObsFound, no need to preserve the order of summation 
	call rpn_comm_allreduce( numObsFound, numObsFoundMPIGlobal, 1, "mpi_integer", "mpi_sum", "grid", ierr )
	
        if ( numObsFoundMPIGlobal > 0 ) then
	  meanObs_ptr( lonIndex, latIndex, 1 ) = meanObs_ptr( lonIndex, latIndex, 1 ) / real( numObsFoundMPIGlobal )
	else  
	  meanObs_ptr( lonIndex, latIndex, 1 ) = MPC_missingValue_R8
	end if
	  
      end do
    end do
    
    if( mpi_myid == 0 ) then
    
      ! Save results
      call gsv_writeToFile( stateVector, './mean_observations.fst', &
                            trim(instrument)//'_'//trim(dayOrNight), &
	                    containsFullField_opt=.true. )
    end if
      
    call gsv_deallocate( stateVector )
    deallocate( headerIndexes )
    deallocate( positionArray )
    
  end subroutine sstb_getMeanObs  

end module SSTbiasEstimation_mod

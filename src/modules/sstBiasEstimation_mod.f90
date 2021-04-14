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
  use kdtree2_mod
  use earthconstants_mod
  use codePrecision_mod
  use mathPhysConstants_mod
  use utilities_mod
  
  implicit none
  save
  private

  ! Public Subroutines
  public :: sstb_computeGriddedObservations
  integer, external :: get_max_rss

  contains
  
  subroutine sstb_computeGriddedObservations( obsData, hco, horizontalSearchRadius )

    implicit none
    
    ! Arguments: 
    type(struct_obs), intent(in)          :: obsData                            ! satellite observations
    type(struct_hco), intent(in), pointer :: hco                                ! horizontal grid structure
    real(8), intent(in)                   :: horizontalSearchRadius             ! horizontal search radius where to search observations

    ! locals
    integer, parameter          :: maxObsPointsSearch = 200000
    character(len=*), parameter :: myName = 'sstb_computeGriddedObservations'
    integer                     :: bodyIndex, headerIndex, bodyCounter
    type(kdtree2), pointer      :: tree => null()
    real(kdkind), allocatable   :: positionArray(:,:)
    type(kdtree2_result)        :: searchResults( maxObsPointsSearch )
    real(kdkind)                :: maxRadius
    real(kdkind)                :: refPosition(3)
    real(pre_obsReal)           :: lat_obs, lon_obs
    integer                     :: lonIndex, latIndex, ni, nj, localObsIndex
    real(kdkind)                :: gridLon, gridLat
    integer                     :: numObsFound
    real(kdkind)                :: searchRadiusSquared
    real(8)                     :: obsMean, currentObs
    real(8), allocatable        :: biasEstimate(:,:)
    
    write(*,*) 'Starting '//myName//'...'
    
    ! create the kdtree on the first call
    if (.not. associated(tree)) then
    
      write(*,*) myName//': start creating kd-tree'
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
      
      write(*,*) myName//': obs_numheader(obsData) = ', obs_numheader(obsData)
      
      allocate(positionArray( 3, obs_numheader(obsData) ))
      
      do headerIndex = 1, obs_numheader( obsData)
      
        lat_obs = obs_headElem_r( obsData, obs_lat, headerIndex )
        lon_obs = obs_headElem_r( obsData, obs_lon, headerIndex )
        lon_obs = lon_obs * MPC_DEGREES_PER_RADIAN_R8
        lat_obs = lat_obs * MPC_DEGREES_PER_RADIAN_R8

        positionArray(:, headerIndex ) = kdtree2_3dPosition( lon_obs, lat_obs )
	
      end do
      
      tree => kdtree2_create( positionArray, sort=.true., rearrange=.true. ) 
      
      write(*,*) myName//': kd-tree done'
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
      
    end if
    
    ! do the search
    
    ni = hco % ni
    nj = hco % nj
    write(*,*) myName//': input grid dimensions, [ni, nj]: ', ni, nj
    
    allocate( biasEstimate( ni, nj ) )
    
    write(*,*) myName//': horizontal search radius, in km : ', horizontalSearchRadius
    searchRadiusSquared = (1.1d0 * horizontalSearchRadius * 1000.d0 )**2 ! also convert from km to m
    
    do lonIndex = 1, ni
      do latIndex = 1, nj
        
	gridLon = real( hco % lon2d_4 ( lonIndex, latIndex ), 8 )
	gridLat = real( hco % lat2d_4 ( lonIndex, latIndex ), 8 )
        refPosition(:) = kdtree2_3dPosition( gridLon, gridLat )
        call kdtree2_r_nearest( tp = tree, qv = refPosition, r2 = searchRadiusSquared, &
                                nfound = numObsFound, &
                                nalloc = maxObsPointsSearch, results = searchResults )
	if ( numObsFound > maxObsPointsSearch ) &
          call utl_abort( myName//': the parameter maxObsPointsSearch must be increased' )
	
	obsMean = 0.0d0
	LOCALOBS: do localObsIndex = 1, numObsFound
	
	  headerIndex = searchResults( localObsIndex ) % idx
          bodyIndex  = obs_headElem_i( obsData, OBS_RLN, headerIndex )
          currentObs = obs_bodyElem_r( obsData, OBS_VAR, bodyIndex   )
          obsMean = obsMean + currentObs 
	  
	end do LOCALOBS
	
	biasEstimate( lonIndex, latIndex ) = obsMean / numObsFound 
	
      end do
    end do  	
    
    deallocate( biasEstimate )
    
  end subroutine sstb_computeGriddedObservations

end module SSTbiasEstimation_mod

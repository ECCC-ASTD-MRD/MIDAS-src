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
  use gridStateVector_mod
  use oceanMask_mod
  use timeCoord_mod

  implicit none
  save
  private

  ! Public Subroutines
  public :: sstb_computeBias
  integer, external :: get_max_rss

  contains
  
  subroutine sstb_computeBias( obsData, hco, vco, numberDays, iceFractionThreshold, &
                               horizontalSearchRadius, numberSatellites , &
                               satelliteList, dateStamp )
    !
    ! :Purpose: compute OmA bias for SST satellite data 
    !  
    implicit none
    
    ! Arguments: 
    type(struct_obs), intent(in)             :: obsData                ! satellite observations
    type(struct_hco), intent(inout), pointer :: hco                    ! horizontal grid structure
    type(struct_vco), intent(in)   , pointer :: vco                    ! vertical grid structure
    integer         , intent(in)             :: numberDays             ! number of days to average temporally A and O 
    real(4)         , intent(in)             :: iceFractionThreshold   ! for ice fraction below it, consider open water      
    real(8)         , intent(in)             :: horizontalSearchRadius ! horizontal search radius where to search observations
    integer         , intent(in)             :: numberSatellites       ! Current satellites: AMSR2, METO-B, METO-A, NOAA19, NPP
    integer         , intent(in)             :: dateStamp              ! date to put into output fst files
    character(len=*), intent(in)             :: satelliteList(:)       ! list of satellite names

    ! locals
    character(len=*), parameter :: myName = 'sstb_computeBias'
    integer                     :: headerIndex, satelliteIndex
    real(8)                     :: averageAna( hco % ni, hco % nj, numberDays ), incrementInHours
    real(8)                     :: averageObs( hco % ni, hco % nj, numberDays, numberSatellites, 2 )
    integer                     :: mask( hco % ni, hco % nj ) 
    real(4)                     :: seaice( hco % ni, hco % nj )
    integer                     :: ierr, newdate
    character(len=*), parameter :: outputFileName = './satellite_bias.fst'
    integer                     :: currentDateStamp, currentDate, currentTime, recordIndex
    type(struct_ocm)            :: oceanMask
    integer                     :: numberOpenWaterPoints, lonIndex, latIndex
    type(struct_gsv)            :: stateVector_ice
    real(4), pointer            :: seaice_ptr( :, :, : )
    
    write(*,*) 'Starting '//myName//'...'
    write(*,*) myName//': Analyse and Observations will be averaged within: ', numberDays, ' days'
    write(*,*) myName//': Current analysis date: ', dateStamp
    write(*,*) myName//': satellites to treat: '
    do satelliteIndex = 1, numberSatellites
      write(*,*) myName//': satellite index: ', satelliteIndex, ', satellite: ', satelliteList( satelliteIndex )
    end do
    write(*,*) myName//': Sea-ice Fraction threshold: ', iceFractionThreshold
    call tim_checkAnchorAnalysesFile( './analysisgrid', numberDays, dateStamp )
 
    ! get latest sea-ice analysis
    call gsv_allocate( stateVector_ice, 1, hco, vco, dataKind_opt = 4, &
                       datestamp_opt = -1, mpi_local_opt = .false.,    &
     		       varNames_opt = (/'LG'/) )
    call gsv_zero( stateVector_ice )
    call gsv_readFromFile( stateVector_ice, './seaice_analysis', 'G6_1_2_1N','A', &
                           unitConversion_opt=.false., containsFullField_opt=.true. )
    call gsv_getField( stateVector_ice, seaice_ptr )
    seaice = seaice_ptr( :, :, 1 )
    call gsv_deallocate( stateVector_ice )

    ! Get land mask from analysisgrid file ( 1=water, 0=land ) 
    ! and the number of open water points
    call ocm_readMaskFromFile( oceanMask, hco, vco, './analysisgrid' )
    numberOpenWaterPoints = 0
    do latIndex = 1, hco % nj
      do lonIndex = 1, hco % ni
        if ( oceanMask%mask ( lonIndex, latIndex, 1 ) ) then
	  mask ( lonIndex, latIndex ) = 1
	  if ( seaice ( lonIndex, latIndex ) <= iceFractionThreshold ) then
	    numberOpenWaterPoints = numberOpenWaterPoints + 1
	  end if  	  
        else
	  mask ( lonIndex, latIndex ) = 0
	end if
      end do
    end do
    call ocm_deallocate( oceanMask )
    
    currentDateStamp = dateStamp
    TIME: do recordIndex = 1, numberDays
    
      ierr = newdate( currentDateStamp, currentDate, currentTime, -3)

      call sstb_averageAna ( averageAna( :, :, recordIndex ), horizontalSearchRadius, &
                             hco, vco, mask, numberOpenWaterPoints, seaice, &
			     iceFractionThreshold, currentDateStamp )

      do satelliteIndex = 1, numberSatellites 
      
        call sstb_averageObs( obsData, hco, vco, horizontalSearchRadius, mask, seaice, &
	                      iceFractionThreshold, trim( satelliteList( satelliteIndex )), 'day'  , &
			      currentDateStamp, currentDate, averageObs( :, : , recordIndex, satelliteIndex, 1 ))
        call sstb_averageObs( obsData, hco, vco, horizontalSearchRadius, mask, seaice, &
	                      iceFractionThreshold, trim( satelliteList( satelliteIndex )), 'night', &
			      currentDateStamp, currentDate, averageObs( :, : , recordIndex, satelliteIndex, 2 ))
      end do
      
      incrementInHours = -24.d0
      call incdatr( currentDateStamp, currentDateStamp, incrementInHours )
      
    end do TIME
    
    do satelliteIndex = 1, numberSatellites 
        
      call sstb_averageOmA( averageAna, averageObs( :, : , :, satelliteIndex, 1),  &
                            trim( satelliteList( satelliteIndex )), 'day', hco, vco, mask, &
                            seaice, iceFractionThreshold, dateStamp, outputFileName )
      call sstb_averageOmA( averageAna, averageObs( :, : , :, satelliteIndex, 2),  &
                            trim( satelliteList( satelliteIndex )), 'night', hco, vco, mask, &
                            seaice, iceFractionThreshold, dateStamp, outputFileName )
    end do			   
    
  end subroutine sstb_computeBias
  
!########################################################################

  subroutine sstb_averageAna ( averageAna, horizontalSearchRadius, &
                               hco, vco, mask, numberWaterPoints, seaice, &
			       iceFractionThreshold, dateStamp )
    !
    ! :Purpose: for every horizontal point of the input horizontal grid,
    !           compute the mean inside a given horizontal search radius
    ! 
    implicit none

    ! Arguments:    
    real(8)         , intent(inout)          :: averageAna( :, : )     ! average anchor analysis
    real(8)         , intent(in)             :: horizontalSearchRadius ! horizontal search radius where to search obs
    type(struct_hco), intent(inout), pointer :: hco                    ! horizontal grid structure
    type(struct_vco), intent(in)   , pointer :: vco                    ! vertical grid structure
    integer         , intent(in)             :: mask( :, : )           ! ocean mask
    integer         , intent(in)             :: numberWaterPoints      ! number of open water points
    real(4)         , intent(in)             :: seaice( :, : )         ! sea-ice fraction
    real(4)         , intent(in)             :: iceFractionThreshold   ! for ice fraction below it, consider open water   
    integer         , intent(in)             :: dateStamp              ! current date stamp
    
    ! locals
    type(struct_gsv)            :: stateVector_ana
    real(4), pointer            :: anchorAnalysis_r4_ptr( :, :, : )
    integer, parameter          :: maxPointsSearch = 200000
    type(kdtree2), pointer      :: tree => null() 
    real(kdkind), allocatable   :: positionArray(:,:)
    type(kdtree2_result)        :: searchResults( maxPointsSearch )
    real(kdkind)                :: maxRadius
    real(kdkind)                :: refPosition(3)
    integer                     :: lonIndex, latIndex, indexCounter
    integer                     :: localLonIndex, localLatIndex, localIndex, ierr
    real(kdkind)                :: lon_grd, lat_grd
    integer                     :: numPointsFound, numPointsFoundMPIGlobal
    real(kdkind)                :: searchRadiusSquared
    integer, allocatable        :: gridPointIndexes(:,:)
    character(len=*), parameter :: myName = 'sstb_averageAna'
     
    call gsv_allocate( stateVector_ana, 1, hco, vco, dataKind_opt = 4, &
                       datestamp_opt = dateStamp, mpi_local_opt = .false., &
		       varNames_opt = (/'TM'/) )
    call gsv_zero( stateVector_ana )
    call gsv_readFromFile( stateVector_ana, './analysisgrid', 'ANCHOR', 'P@', &
                           unitConversion_opt=.false., containsFullField_opt=.true. )
    call gsv_getField( stateVector_ana, anchorAnalysis_r4_ptr )

    averageAna = 0.0d0
    
    allocate( positionArray( 3,  numberWaterPoints ))
    allocate( gridPointIndexes( 2, numberWaterPoints ))
      
    indexCounter = 0
    do lonIndex = 1, hco % ni
      do latIndex = 1, hco % nj
      
        if ( mask( lonIndex, latIndex ) == 1 .and. &
	   seaice( lonIndex, latIndex ) <= iceFractionThreshold ) then 
        
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

    ! do the search
    searchRadiusSquared = ( 1.1d0 * horizontalSearchRadius * 1000.d0 )**2 ! convert from km to m2

    do lonIndex = 1, hco % ni
      do latIndex = 1, hco % nj
      
        ! do averaging for every water point
        if (   mask( lonIndex, latIndex ) == 1 .and. & 
             seaice( lonIndex, latIndex ) <= iceFractionThreshold ) then
	       
	  lon_grd = real( hco % lon2d_4 ( lonIndex, latIndex ), 8 )
	  lat_grd = real( hco % lat2d_4 ( lonIndex, latIndex ), 8 )
          refPosition(:) = kdtree2_3dPosition( lon_grd, lat_grd )
	
          call kdtree2_r_nearest( tp = tree, qv = refPosition, r2 = searchRadiusSquared, &
	                          nfound = numPointsFound, nalloc = maxPointsSearch, results = searchResults )
	  
          do localIndex = 1, numPointsFound
	  
	    localLonIndex = gridPointIndexes( 1, searchResults( localIndex ) % idx )
	    localLatIndex = gridPointIndexes( 2, searchResults( localIndex ) % idx )
	    averageAna( lonIndex, latIndex ) = averageAna( lonIndex, latIndex ) + &
	                                       anchorAnalysis_r4_ptr( localLonIndex, localLatIndex, 1 )
	  end do 

	  ! summing the values over all mpi tasks and sending them back to all tasks preserving the order of summation
          call mpi_allreduce_sumreal8scalar( averageAna( lonIndex, latIndex ), "grid" )
	  call rpn_comm_allreduce( numPointsFound, numPointsFoundMPIGlobal, 1, "mpi_integer", "mpi_sum", "grid", ierr )
	  
	  averageAna( lonIndex, latIndex ) = averageAna( lonIndex, latIndex ) / real( numPointsFoundMPIGlobal )
	  
        else if ( mask( lonIndex, latIndex ) == 0 .or. & 
                seaice( lonIndex, latIndex ) >= iceFractionThreshold ) then
	
	  averageAna( lonIndex, latIndex ) = MPC_missingValue_R8
	  
	end if 

      end do
    end do
       
    deallocate( positionArray )
    deallocate( gridPointIndexes )
     
    call gsv_deallocate( stateVector_ana )

  end subroutine sstb_averageAna
   
!########################################################################

  subroutine sstb_averageObs( obsData, hco, vco, horizontalSearchRadius, mask, &
                              seaice, iceFractionThreshold, instrument, dayOrNight, &
                              dateStamp, currentDate, averageObs )
    !
    ! :Purpose: for every horizontal point of the input horizontal grid,
    !           compute mean Obs within a given horizontal search radius,
    !           The mean is computed for a given instrument and day / night time. 
    !
    implicit none
    
    ! Arguments: 
    type(struct_obs), intent(in)             :: obsData                ! satellite observations
    type(struct_hco), intent(inout), pointer :: hco                    ! horizontal grid structure
    type(struct_vco), intent(in)   , pointer :: vco                    ! vertical grid structure
    real(8)         , intent(in)             :: horizontalSearchRadius ! horizontal search radius where to search obs
    integer         , intent(in)             :: mask( :, : )           ! ocean mask
    real(4)         , intent(in)             :: seaice( :, : )         ! sea-ice fraction
    real(4)         , intent(in)             :: iceFractionThreshold   ! for ice fraction below it, open water      
    character(len=*), intent(in)             :: instrument             ! name of instrument
    character(len=*), intent(in)             :: dayOrNight             ! look for day or night obs
    integer         , intent(in)             :: dateStamp              ! date to put into output fst files
    integer         , intent(in)             :: currentDate            ! current date
    real(8)         , intent(out)            :: averageObs( :, : )     ! average obs for a given day
    
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
    real(8)                     :: currentObs
    character(len=*), parameter :: myName = 'sstb_averageObs'
     
    countObs = 0 
    do headerIndex = 1, obs_numheader( obsData )
    
      if ( obs_headElem_i( obsData, obs_dat, headerIndex ) == currentDate .and. &
           obs_elem_c( obsData, 'STID' , headerIndex ) == trim(instrument) ) then
        
	if ( trim(dayOrNight) == 'day'   ) then
          if ( obs_headElem_r( obsData, obs_sun, headerIndex ) <  solarZenithThreshold ) &
	    countObs = countObs + 1
        else if ( trim(dayOrNight) == 'night' ) then
          if ( obs_headElem_r( obsData, obs_sun, headerIndex ) >= solarZenithThreshold ) &
	    countObs = countObs + 1
	end if  
	
      end if	
		
    end do
    
    averageObs = 0.0d0
    
    if ( trim(dayOrNight) == 'day'   ) then
      write(*,*) myName//': found ', countObs, ' ', trim(instrument), ' day observations'
    else if ( trim(dayOrNight) == 'night' ) then
      write(*,*) myName//': found ', countObs, ' ', trim(instrument), ' night observations'
    end if	
      	
    allocate( positionArray( 3, countObs ))
    allocate( headerIndexes( countObs ))

    headerCounter = 0
    HEADER: do headerIndex = 1, obs_numheader( obsData )
      
      if ( obs_headElem_i( obsData, obs_dat, headerIndex ) == currentDate .and. &
           obs_elem_c( obsData, 'STID' , headerIndex ) == trim(instrument) ) then
      
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
      	    
    end do HEADER
    
    nullify(tree)
    tree => kdtree2_create( positionArray, sort=.true., rearrange=.true. ) 
    
    ! do the search
    searchRadiusSquared = ( 1.1d0 * horizontalSearchRadius * 1000.d0 )**2 ! convert from km to m2
    
    do lonIndex = 1, hco % ni
      do latIndex = 1, hco % nj
      
        ! compute bias for every open water point
        if (   mask( lonIndex, latIndex ) == 1 .and. &
	     seaice( lonIndex, latIndex ) <= iceFractionThreshold ) then 
        
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
 	      currentObs = obs_bodyElem_r( obsData, obs_var, bodyIndex ) - MPC_K_C_DEGREE_OFFSET_R8
              averageObs( lonIndex, latIndex ) = averageObs( lonIndex, latIndex ) + currentObs 
	      
	    end do
	    
	  else
	    
	    averageObs( lonIndex, latIndex ) = 0.d0
	  
	  end if
	  
	  ! summing the values over all mpi tasks and sending them back to all tasks preserving the order of summation
          call mpi_allreduce_sumreal8scalar( averageObs( lonIndex, latIndex ), "grid" )
	  call rpn_comm_allreduce( numObsFound, numObsFoundMPIGlobal, 1, "mpi_integer", "mpi_sum", "grid", ierr )

          if ( numObsFoundMPIGlobal > 0 ) then

	    averageObs( lonIndex, latIndex ) = averageObs( lonIndex, latIndex ) / real( numObsFoundMPIGlobal )
	    
	  else
	  
	    averageObs( lonIndex, latIndex ) = MPC_missingValue_R8
	  
	  end if  
	 
	else if (   mask( lonIndex, latIndex ) == 0 .or. &
	          seaice( lonIndex, latIndex ) >= iceFractionThreshold ) then
		   
	  averageObs( lonIndex, latIndex ) = MPC_missingValue_R8
	
	end if 
	  
      end do
    end do
    
    deallocate( headerIndexes )
    deallocate( positionArray )
    
  end subroutine sstb_averageObs  


  subroutine sstb_averageOmA( averageAna, averageObs, instrument, dayOrNight, &
                              hco, vco, mask, seaice, &
			      iceFractionThreshold, dateStamp, outputFileName )
    !
    ! :Purpose: Compute satellite bias as time averaged OmA and save it into a std.file 
    ! 
    implicit none
  
    real(8)         , intent(inout)          :: averageAna( :, :, : ) ! average anchor analysis
    real(8)         , intent(inout)          :: averageObs( :, :, : ) ! average obs 
    character(len=*), intent(in)             :: instrument            ! name of instrument
    character(len=*), intent(in)             :: dayOrNight            ! look for day or night obs
    type(struct_hco), intent(inout), pointer :: hco                   ! horizontal grid structure
    type(struct_vco), intent(in)   , pointer :: vco                   ! vertical grid structure
    integer         , intent(in)             :: mask( :, : )          ! ocean mask
    real(4)         , intent(in)             :: seaice( :, : )        ! sea-ice fraction  
    real(4)         , intent(in)             :: iceFractionThreshold  ! for ice fraction below it, consider open water      
    integer         , intent(in)             :: dateStamp             ! current date stamp
    character(len=*), intent(in)             :: outputFileName        ! output file name
 
    ! locals
    character(len=*), parameter :: myName = 'sstb_averageOmA'
    type(struct_gsv) :: stateVectorBias
    real(8), pointer :: bias_ptr( :, :, : )
    integer          :: numberDays, recordIndex, countStates
    integer          :: lonIndex, latIndex, countStatesMPIGlobal, ierr

    write(*,*) myName//': computing OmA: ', instrument, ' ', dayOrNight

    numberDays = size( averageAna, dim = 3 )

    call gsv_allocate( stateVectorBias, 1, hco, vco, dataKind_opt = 8, &
                       datestamp_opt = dateStamp, mpi_local_opt=.false., &
		       varNames_opt=(/'TM'/) )
    call gsv_getField( stateVectorBias, bias_ptr )
    bias_ptr = 0.0d0
    
    do lonIndex = 1, hco % ni
      do latIndex = 1, hco % nj
     
        if ( mask( lonIndex, latIndex ) == 1 .and. &
	     seaice( lonIndex, latIndex ) <= iceFractionThreshold ) then
        
	  countStates = 0 
          do recordIndex = 1, numberDays
	  
    	    if ( averageObs( lonIndex, latIndex, recordIndex ) /= MPC_missingValue_R8 ) then
	  
              countStates = countStates + 1
              bias_ptr( lonIndex, latIndex, 1 ) = bias_ptr( lonIndex, latIndex, 1 ) + &
	                                          averageObs( lonIndex, latIndex, recordIndex ) - &
	      	      			          averageAna( lonIndex, latIndex, recordIndex )
	    end if
	    					  
	  end do
	     
	  ! summing the values over all mpi tasks and sending them back 
	  ! to all tasks preserving the order of summation
          call mpi_allreduce_sumreal8scalar( bias_ptr( lonIndex, latIndex, 1 ), "grid" )
	  call rpn_comm_allreduce( countStates, countStatesMPIGlobal, 1, "mpi_integer", "mpi_sum", "grid", ierr )
	
	  if ( countStatesMPIGlobal > 0 ) then
	    bias_ptr( lonIndex, latIndex, 1 ) = bias_ptr( lonIndex, latIndex, 1 ) / real(countStatesMPIGlobal)
	  else  
            bias_ptr( lonIndex, latIndex, 1 ) = MPC_missingValue_R8
	  end if
		    					  
	else
	  
          bias_ptr( lonIndex, latIndex, 1 ) = MPC_missingValue_R8
	 
	end if  

      end do 	  
    end do
    
    ! Save results
    if ( mpi_myid == 0 ) &
    call gsv_writeToFile( stateVectorBias, outputFileName, 'B_'//trim(instrument)//'_'//trim(dayOrNight))
   

    call gsv_deallocate( stateVectorBias )

  end subroutine sstb_averageOmA

end module SSTbias_mod

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
!-------------------------------------- LICENCE end --------------------------------------

module oceanObservations_mod
  ! MODULE oceobs_mod (prefix='oceobs' category='1. High-level functionality')
  !
  ! :Purpose: storage for ocean observations related subroutines
  !
  use mpi_mod
  use mpivar_mod
  use utilities_mod
  use obsSpaceData_mod
  use codtyp_mod
  use gridStateVector_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use oceanMask_mod
  use timeCoord_mod
  use codePrecision_mod  
  use sqliteRead_mod
  use bufr_mod
  use MathPhysConstants_mod
  
  implicit none

  save
  private

  ! Public functions/subroutines
  public :: oceobs_pseudoSST
 
  ! External functions
  integer, external :: fnom, fclos  

  ! mpi topology
  integer           :: myLatBeg, myLatEnd
  integer           :: myLonBeg, myLonEnd
  integer           :: latPerPE, latPerPEmax, lonPerPE, lonPerPEmax

  type, public :: exceptRegion_type
    character (len=12) :: name   = ''    ! name of the region
    real(4)            :: SST    = 271.4 ! output SST value for this region
    real(4)            :: lonMin = 0.0   ! 
    real(4)            :: lonMax = 0.0   !
    real(4)            :: latMin = 0.0   !
    real(4)            :: latMax = 0.0   !
  end type exceptRegion_type

  contains

  !----------------------------------------------------------------------------------------
  ! oceobs_pseudoSST
  !----------------------------------------------------------------------------------------
  subroutine oceobs_pseudoSST( hco, vco, iceFractionThreshold, outputSST, outputFreshWaterST, seaiceThinning, &
                               outputFileName, etiket, useSeaWaterFraction, numberExceptRegions, exceptRegion )
    !
    !: Purpose: to generate pseudo SST data  
    !           
    
    implicit none
    
    ! Arguments
    type(struct_hco)       , intent(inout), pointer :: hco                  ! horizontal grid structure
    type(struct_vco)       , intent(in)   , pointer :: vco                  ! vertical grid structure
    real(4)                , intent(in)             :: iceFractionThreshold ! consider no ice condition below this threshold
    real(4)                , intent(in)             :: outputSST            ! output SST value for pseudo observations
    real(4)                , intent(in)             :: outputFreshWaterST   ! output fresh water surface temperature for pseudo observations
    integer                , intent(in)             :: seaiceThinning       ! generate pseudo obs in every 'seaiceThinning' points   
    character(len=20)      , intent(in)             :: outputFileName    
    character(len=20)      , intent(in)             :: etiket    
    logical                , intent(in)             :: useSeaWaterFraction    
    integer                , intent(in)             :: numberExceptRegions
    type(exceptRegion_type), intent(in)             :: exceptRegion(1:numberExceptRegions)
    
    ! Locals:
    type(struct_gsv)            :: stateVector_ice, stateVector_salinity
    real(4), pointer            :: seaice_ptr( :, :, : ), salinity_ptr( :, :, : )
    type(struct_ocm)            :: oceanMask
    integer                     :: numberIceCoveredPoints, lonIndex, latIndex, dateStamp
    integer                     :: datePrint, timePrint, imode
    integer                     :: randomSeed, newDate, ierr 
    integer, allocatable        :: seaiceDomainIndexes(:), seaiceDomainIndexesTemp(:)
    real(4), allocatable        :: seaWaterFraction(:), seaWaterFractionTemp(:)
    real(4), allocatable        :: seaiceLons(:), seaiceLats(:), seaiceLonsTemp(:), seaiceLatsTemp(:)  
    type(struct_obs)            :: obsData              ! obsSpaceData   
    character(len=*), parameter :: myName = 'oceobs_pseudoSST'
    
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

    if ( useSeaWaterFraction ) then
      ! read salinity fraction
      call gsv_allocate( stateVector_salinity, 1, hco, vco, dataKind_opt = 4, &
                         datestamp_opt = -1, mpi_local_opt = .false.,    &
                         varNames_opt = (/'VF'/) )
      call gsv_readFromFile( stateVector_salinity, './seaice_analysis', ' ','A', &
                           unitConversion_opt=.false., containsFullField_opt=.true. )
      call gsv_getField( stateVector_salinity, salinity_ptr )
      
    end if  

    ! Get land mask from analysisgrid file ( 1=water, 0=land ) 
    ! and the number of open water points
    call ocm_readMaskFromFile( oceanMask, hco, vco, './analysisgrid' )

    allocate( seaiceDomainIndexesTemp( ( myLonEnd - myLonBeg + 1 ) * ( myLatEnd - myLatBeg + 1 ) ))
    if ( useSeaWaterFraction ) &
      allocate( seaWaterFractionTemp( ( myLonEnd - myLonBeg + 1 ) * ( myLatEnd - myLatBeg + 1 ) ))
    allocate( seaiceLonsTemp( ( myLonEnd - myLonBeg + 1 ) * ( myLatEnd - myLatBeg + 1 ) ))
    allocate( seaiceLatsTemp( ( myLonEnd - myLonBeg + 1 ) * ( myLatEnd - myLatBeg + 1 ) ))
    
    numberIceCoveredPoints = 0
    do lonIndex = myLonBeg, myLonEnd 
      do latIndex = myLatBeg, myLatEnd
        if ( oceanMask%mask( lonIndex, latIndex, 1 ) ) then
          if ( seaice_ptr( lonIndex, latIndex, 1 ) > iceFractionThreshold ) then
            numberIceCoveredPoints = numberIceCoveredPoints + 1
	    seaiceDomainIndexesTemp( numberIceCoveredPoints ) = numberIceCoveredPoints
	    if ( useSeaWaterFraction ) &
	      seaWaterFractionTemp( numberIceCoveredPoints ) = salinity_ptr( lonIndex, latIndex, 1 )
	    seaiceLonsTemp( numberIceCoveredPoints ) = hco % lon2d_4 ( lonIndex, latIndex )
	    seaiceLatsTemp( numberIceCoveredPoints ) = hco % lat2d_4 ( lonIndex, latIndex )	
          end if
        end if
      end do
    end do
    call ocm_deallocate( oceanMask )
    call gsv_deallocate( stateVector_ice )
    if ( useSeaWaterFraction ) call gsv_deallocate( stateVector_salinity )
    write(*,*) myName//': ', numberIceCoveredPoints, ' ice-covered points found'
    
    allocate( seaiceDomainIndexes( numberIceCoveredPoints ))
    if ( useSeaWaterFraction ) allocate( seaWaterFraction( numberIceCoveredPoints ))
    allocate( seaiceLons( numberIceCoveredPoints ))
    allocate( seaiceLats( numberIceCoveredPoints ))
    
    seaiceDomainIndexes(:) = seaiceDomainIndexesTemp( 1 : numberIceCoveredPoints )
    if ( useSeaWaterFraction ) seaWaterFraction(:) = seaWaterFractionTemp( 1 : numberIceCoveredPoints )
    seaiceLons(:) = seaiceLonsTemp( 1 : numberIceCoveredPoints )
    seaiceLats(:) = seaiceLatsTemp( 1 : numberIceCoveredPoints )
    
    deallocate( seaiceDomainIndexesTemp )
    if ( useSeaWaterFraction ) deallocate( seaWaterFractionTemp )
    deallocate( seaiceLonsTemp )
    deallocate( seaiceLatsTemp )
    
        
    dateStamp = tim_getDatestampFromFile( './seaice_analysis' )
    write(*,*) myName//': datestamp: ', dateStamp 
    ! compute random seed from the date for randomly forming sea-ice subdomain
    imode = -3 ! stamp to printable date and time: YYYYMMDD, HHMMSShh
    ierr = newdate( dateStamp, datePrint, timePrint, imode )
    timePrint = timePrint / 1000000
    datePrint =  datePrint * 100 + timePrint
    
    ! Remove the century, keeping 2 digits of the year
    randomSeed = datePrint - 100000000 * ( datePrint / 100000000 )
    write(*,*) myName//': datePrint, timePrint: ', datePrint, timePrint 
    call utl_randomOrderInt( seaiceDomainIndexes, randomSeed )
    write(*,*) myName//': seed for random shuffle of sea-ice points: ', randomSeed

    if ( useSeaWaterFraction ) then
      call oceobs_computeObsData( obsData, seaiceDomainIndexes, seaiceLons, seaiceLats, &
                                  seaiceThinning, outputSST, outputFreshWaterST, exceptRegion, &
                                  outputFileName, etiket, datePrint, timePrint, &
                                  seaWaterFraction_opt = seaWaterFraction )
    else
      call oceobs_computeObsData( obsData, seaiceDomainIndexes, seaiceLons, seaiceLats, &
                                  seaiceThinning, outputSST, outputFreshWaterST, exceptRegion, &
                                  outputFileName, etiket, datePrint, timePrint )
    end if				  
    
    deallocate( seaiceDomainIndexes )
    if ( useSeaWaterFraction ) deallocate( seaWaterFraction )
    deallocate( seaiceLons )
    deallocate( seaiceLats )

    write(*,*) myName//': done'
    
  end subroutine oceobs_pseudoSST

  !--------------------------------------------------------------------------
  ! oceobs_computeObsData
  !--------------------------------------------------------------------------

  subroutine oceobs_computeObsData( obsData, seaiceDomainIndexes, seaiceLons, seaiceLats, seaiceThinning, &
                                    outputSST, outputFreshWaterST, exceptRegion, outputFileName, &
                                    etiket, datePrint, timePrint, seaWaterFraction_opt )
    !
    !: Purpose: pseudo SST data are put into obsSpaceData  
    !           and written into an SQLite file
    
    implicit none
    
    ! Arguments
    type(struct_obs)       , intent(inout)        :: obsData                ! obsSpaceData   
    integer                , intent(in)           :: seaiceDomainIndexes(:) ! array of the ice-covered point indexes
    real(4)                , intent(in)           :: seaiceLons(:)          ! longitudes of sea ice 
    real(4)                , intent(in)           :: seaiceLats(:)          ! latitudes of sea ice 
    integer                , intent(in)           :: seaiceThinning         ! generate pseudo obs in every 'seaiceThinning' points   
    real(4)                , intent(in)           :: outputSST              ! output SST value for pseudo observations
    real(4)                , intent(in)           :: outputFreshWaterST     ! output fresh water surface temperature for pseudo observations
    type(exceptRegion_type), intent(in)           :: exceptRegion(:)
    character(len=20)      , intent(in)           :: outputFileName    
    character(len=20)      , intent(in)           :: etiket    
    integer                , intent(in)           :: datePrint
    integer                , intent(in)           :: timePrint
    real(4)                , intent(in), optional :: seaWaterFraction_opt(:)! sea water fraction data:
                                                                            ! 0: fresh water
                                                                            ! 1: sea water
    
    ! Locals
    real(pre_obsReal)           :: obsLon, obsLat, obsLon_degrees, obsLat_degrees
    real(4)                     :: obsValue, seaWaterFraction
    integer                     :: seaiceIndex, seaiceDomainDimension, pseudoObsDimension
    integer                     :: codeType, headerIndex, regionIndex, numberExceptRegions
    integer                     :: coordinatesIndex
    character(len=*), parameter :: myName = 'oceobs_computeObsData'
    
    seaiceDomainDimension = size( seaiceDomainIndexes )
    pseudoObsDimension = floor( real( seaiceDomainDimension ) / real ( seaiceThinning ) ) + 1
    numberExceptRegions = size( exceptRegion )
    write(*,*) myName//': sea-ice domain dimension: ', seaiceDomainDimension
    write(*,*) myName//': pseudo obs vector dimension: ', pseudoObsDimension
    write(*,*) myName//': exception regions: '
    do regionIndex = 1, numberExceptRegions
      write(*,*) myName//': ', regionIndex, ', ', exceptRegion( regionIndex ) % name, &
                 ', SST: ', exceptRegion( regionIndex ) % SST
    end do   

    call obs_class_initialize( 'VAR' )
    call obs_initialize( obsData, numHeader_max = pseudoObsDimension, numBody_max = pseudoObsDimension, mpi_local= .true. )
    codeType = codtyp_get_codtyp( 'pseudosfc' )
    
    headerIndex = 1
    
    do seaiceIndex = 1, seaiceDomainDimension, seaiceThinning 
	  
      coordinatesIndex = seaiceDomainIndexes( seaiceIndex )
      
      obsLon   = seaiceLons( coordinatesIndex )
      obsLat   = seaiceLats( coordinatesIndex )
      if ( present( seaWaterFraction_opt ) ) then
        seaWaterFraction = seaWaterFraction_opt( coordinatesIndex )
        obsValue = ( 1.0d0 - seaWaterFraction ) * outputFreshWaterST + seaWaterFraction * outputSST
      else
        obsValue = outputSST
      end if	
      
      ! compute exception regions as in operational OI approach      
      do regionIndex = 1, numberExceptRegions
        obsLon_degrees = obsLon * MPC_DEGREES_PER_RADIAN_R8
        if ( obsLon_degrees > 180. ) obsLon_degrees = obsLon_degrees - 360.
        obsLat_degrees = obsLat * MPC_DEGREES_PER_RADIAN_R8
        if ( obsLon_degrees > exceptRegion( regionIndex ) % lonMin .and. &
             obsLon_degrees < exceptRegion( regionIndex ) % lonMax .and. &
             obsLat_degrees > exceptRegion( regionIndex ) % latMin .and. &
             obsLat_degrees < exceptRegion( regionIndex ) % latMax       ) then
          obsValue = exceptRegion( regionIndex ) % SST
        end if
      end do	  	     
    
      call obs_setFamily( obsData, 'SF'   , headerIndex )
      call obs_headSet_i( obsData, OBS_ONM, headerIndex, headerIndex )
      call obs_headSet_i( obsData, OBS_ITY, headerIndex, codeType    )
      call obs_headSet_r( obsData, OBS_LAT, headerIndex, obsLat      )
      call obs_headSet_r( obsData, OBS_LON, headerIndex, obsLon      )
      call obs_bodySet_r( obsData, OBS_VAR, headerIndex, obsValue    )
      call obs_bodySet_i( obsData, OBS_VNM, headerIndex, bufr_sst    )  
      call     obs_set_c( obsData, 'STID' , headerIndex, 'ABOG'      )
      call obs_headSet_i( obsData, OBS_NLV, headerIndex, 1           )
      call obs_headSet_i( obsData, OBS_RLN, headerIndex, headerIndex )
      
      headerIndex = headerIndex + 1
      
    end do  

    call sqlr_writeSqlDiagFile( obsData, 'SF', .false., .false., outputFileName, &
                                pseudoObs_opt = .true., etiket_opt = etiket, &
				datePrint_opt = datePrint, timePrint_opt = timePrint ) 
 
    ! Deallocate obsSpaceData
    call obs_finalize( obsData )

  end subroutine oceobs_computeObsData
  
end module oceanObservations_mod  

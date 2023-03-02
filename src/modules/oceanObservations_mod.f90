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
  ! MODULE oceanObservations_mod (prefix='oobs' category='1. High-level functionality')
  !
  ! :Purpose: storage for ocean observations related subroutines
  !
  use midasMpi_mod
  use utilities_mod
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
  
  implicit none

  save
  private

  ! Public functions/subroutines
  public :: oobs_pseudoSST
 
  ! External functions
  integer, external :: fnom, fclos  

  ! mpi topology
  integer :: myLatBeg, myLatEnd
  integer :: myLonBeg, myLonEnd
  integer :: latPerPE, latPerPEmax, lonPerPE, lonPerPEmax

  contains

  !----------------------------------------------------------------------------------------
  ! oobs_pseudoSST
  !----------------------------------------------------------------------------------------
  subroutine oobs_pseudoSST(hco, vco, iceFractionThreshold, outputSST, outputFreshWaterST, &
                            iceThinning, outputFileName, etiket, seaWaterThreshold)
    !
    !: Purpose: to generate pseudo SST data  
    !           
    
    implicit none
    
    ! Arguments
    type(struct_hco) , intent(inout), pointer :: hco                  ! horizontal grid structure
    type(struct_vco) , intent(in)   , pointer :: vco                  ! vertical grid structure
    real(4)          , intent(in)             :: iceFractionThreshold ! consider no ice condition below this threshold
    real(4)          , intent(in)             :: outputSST            ! output SST value for pseudo observations
    real(4)          , intent(in)             :: outputFreshWaterST   ! output fresh water surface temperature for pseudo observations
    integer          , intent(in)             :: iceThinning          ! generate pseudo obs in every 'iceThinning' points   
    character(len=*) , intent(in)             :: outputFileName    
    character(len=*) , intent(in)             :: etiket    
    real(4)          , intent(in)             :: seaWaterThreshold    ! to distinguish inland water from sea water  
    
    ! Locals:
    type(struct_gsv)            :: stateVector_ice, stateVector_seaWater
    real(4), pointer            :: seaIce_ptr(:, :, :), seaWater_ptr(:, :, :)
    type(struct_ocm)            :: oceanMask
    integer                     :: numberIceCoveredPoints, lonIndex, latIndex, dateStamp, inlandWaterPoints
    integer                     :: datePrint, timePrint, imode, seaWaterPoints
    integer                     :: randomSeed, newDate, ierr 
    integer, allocatable        :: iceDomainIndexesAux(:), iceDomainIndexes(:)
    real(4), allocatable        :: seaWaterFractionAux(:), iceLonsAux(:), iceLatsAux(:)  
    real(4), allocatable        :: seaWaterFraction(:), iceLons(:), iceLats(:)
    type(struct_obs)            :: obsData   
    
    ! get mpi topology
    call mmpi_setup_lonbands(hco%ni, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd)
    call mmpi_setup_latbands(hco%nj, latPerPE, latPerPEmax, myLatBeg, myLatEnd)

    ! get latest sea-ice analysis
    call gsv_allocate(stateVector_ice, 1, hco, vco, dataKind_opt = 4, &
                      datestamp_opt = -1, mpi_local_opt = .false., varNames_opt = (/'LG'/))
    call gio_readFromFile(stateVector_ice, './seaice_analysis', ' ','A', &
                          unitConversion_opt=.false., containsFullField_opt=.true.)
    call gsv_getField(stateVector_ice, seaIce_ptr)

    ! read sea water fraction
    call gsv_allocate(stateVector_seaWater, 1, hco, vco, dataKind_opt = 4, &
                      datestamp_opt = -1, mpi_local_opt = .false., varNames_opt = (/'VF'/))
    call gio_readFromFile(stateVector_seaWater, './seaice_analysis', ' ','A', &
                          unitConversion_opt=.false., containsFullField_opt=.true.)
    call gsv_getField(stateVector_seaWater, seaWater_ptr)

    ! Get land mask from analysisgrid file (1=water, 0=land)
    call ocm_readMaskFromFile(oceanMask, hco, vco, './analysisgrid')

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
          
	  end if
        end if
      end do
    end do
    call ocm_deallocate(oceanMask)
    call gsv_deallocate(stateVector_ice)
    call gsv_deallocate(stateVector_seaWater)
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
    end if  
    deallocate(iceLonsAux)
    deallocate(iceLatsAux)
    deallocate(iceDomainIndexesAux)
    deallocate(seaWaterFractionAux)
        
    dateStamp = tim_getDatestampFromFile('./seaice_analysis', varNameForDate_opt = 'LG')
    write(*,*) 'oobs_pseudoSST: datestamp: ', dateStamp 
    ! compute random seed from the date for randomly forming sea-ice subdomain
    imode = -3 ! stamp to printable date and time: YYYYMMDD, HHMMSShh
    ierr = newdate(dateStamp, datePrint, timePrint, imode)
    timePrint = timePrint / 1000000
    datePrint =  datePrint * 100 + timePrint
    
    ! Remove the century, keeping 2 digits of the year
    randomSeed = datePrint - 100000000 * (datePrint / 100000000)
    write(*,*) 'oobs_pseudoSST: datePrint, timePrint: ', datePrint, timePrint 

    if (numberIceCoveredPoints > 0) then

      call utl_randomOrderInt(iceDomainIndexes, randomSeed)
      write(*,*) 'oobs_pseudoSST: seed for random shuffle of sea-ice points: ', randomSeed
    
      call oobs_computeObsData(obsData, iceDomainIndexes, iceLons, iceLats, &
                               iceThinning, outputSST, outputFreshWaterST, &
                               outputFileName, etiket, datePrint, timePrint, &
                               seaWaterFraction, seaWaterThreshold, inlandWaterPoints)
    else 
    
      call obs_initialize(obsData, numHeader_max = 0, numBody_max = 0, mpi_local= .true.)
      call sqlr_writeEmptyPseudoSSTobsFile(obsData, 'SF', outputFileName, etiket, datePrint, timePrint)   

    end if

    if (numberIceCoveredPoints > 0) then
      deallocate(iceLons)
      deallocate(iceLats)
      deallocate(iceDomainIndexes)
      deallocate(seaWaterFraction)
    end if  
    
    write(*,*) 'oobs_pseudoSST: done'
    
  end subroutine oobs_pseudoSST

  !--------------------------------------------------------------------------
  ! oobs_computeObsData
  !--------------------------------------------------------------------------

  subroutine oobs_computeObsData(obsData, iceDomainIndexes, iceLons, iceLats, iceThinning, &
                                 outputSST, outputFreshWaterST, outputFileName, &
                                 etiket, datePrint, timePrint, &
                                 seaWaterFraction, seaWaterThreshold, inlandWaterPoints)
    !
    !: Purpose: pseudo SST data are put into obsSpaceData  
    !           and written into an SQLite file
    
    implicit none
    
    ! Arguments
    type(struct_obs) , intent(inout) :: obsData            ! obsSpaceData   
    integer          , intent(in)    :: iceDomainIndexes(:)! array of the ice-covered point indexes
    real(4)          , intent(in)    :: iceLons(:)         ! longitudes of sea ice 
    real(4)          , intent(in)    :: iceLats(:)         ! latitudes of sea ice 
    integer          , intent(in)    :: iceThinning        ! generate pseudo obs in every 'iceThinning' points   
    real(4)          , intent(in)    :: outputSST          ! output SST value for pseudo observations
    real(4)          , intent(in)    :: outputFreshWaterST ! output fresh water surface temperature for pseudo obs
    character(len=*) , intent(in)    :: outputFileName    
    character(len=*) , intent(in)    :: etiket    
    integer          , intent(in)    :: datePrint
    integer          , intent(in)    :: timePrint
    real(4)          , intent(in)    :: seaWaterFraction(:)! sea water fraction data: 0: fresh water; 1: sea water
    real(4)          , intent(in)    :: seaWaterThreshold  ! to distinguish inland water from sea water 
    integer          , intent(in)    :: inlandWaterPoints  ! number of inland water points 
     
    ! Locals
    real(pre_obsReal)           :: obsLon, obsLat
    real(4)                     :: obsValue
    integer                     :: iceIndex, iceDomainDimension, pseudoObsDimension
    integer                     :: codeType, headerIndex
    integer                     :: coordinatesIndex, counterThinning, checkInlandWatersCount, checkSeaWatersCount
    character(len=*), parameter :: myName = 'oobs_computeObsData'
    
    iceDomainDimension = size(iceDomainIndexes)
    pseudoObsDimension = floor(real((iceDomainDimension - inlandWaterPoints) / iceThinning)) + inlandWaterPoints
      
    write(*,*) myName//': sea-ice domain dimension: ', iceDomainDimension
    write(*,*) myName//': pseudo obs vector dimension: ', pseudoObsDimension
    write(*,*) myName//': pseudo SST obs will be generated in every ', iceThinning, &
    ' points of the sea-ice field for sea water, '
    write(*,*) myName//': and in every point for inland waters, where sea water fraction <= ', seaWaterThreshold
    write(*,*) myName//': number of inland waters points: ', inlandWaterPoints  
    
    call obs_initialize(obsData, numHeader_max = pseudoObsDimension, numBody_max = pseudoObsDimension, mpi_local= .true.)
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
        obsValue = (1.0d0 - seaWaterFraction(coordinatesIndex)) * outputFreshWaterST + &
                   seaWaterFraction(coordinatesIndex)* outputSST
        checkInlandWatersCount =  checkInlandWatersCount + 1
      else
        if (counterThinning == iceThinning) then  
          obsValue = outputSST
          checkSeaWatersCount = checkSeaWatersCount + 1	
          counterThinning = 1
        else
          counterThinning = counterThinning + 1
          cycle
        end if
      end if

      call obs_setFamily(obsData, 'SF'   , headerIndex)
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
    
    call sqlr_writePseudoSSTobs(obsData, 'SF', outputFileName, etiket, datePrint, timePrint) 
 
    ! Deallocate obsSpaceData
    call obs_finalize(obsData)

  end subroutine oobs_computeObsData

end module oceanObservations_mod  

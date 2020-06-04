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

module thinning_mod
  ! MODULE thinning (prefix='thn' category='1. High-level functionality')
  !
  ! :Purpose: Using observation-type-specific algorithms, set bit 11 of 'flag'
  !           on data that are not to be assimilated.
  !
  ! :Note:    This module is intended to group all of the thinning methods in a
  !           single fortran module.
  !
  use mpi_mod
  use bufr_mod
  use mathPhysConstants_mod
  use earthConstants_mod
  use obsSpaceData_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use gridStateVector_mod
  use codtyp_mod
  use physicsFunctions_mod
  use utilities_mod
  use kdtree2_mod
  implicit none
  private

  public :: thn_thinHyper, thn_thinTovs, thn_thinAircraft, thn_thinAladin

  integer, external :: get_max_rss

contains

  !--------------------------------------------------------------------------
  ! thn_thinAircraft
  !--------------------------------------------------------------------------
  subroutine thn_thinAircraft(obsdat)
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat

    ! Locals:
    integer :: nulnam
    integer :: fnom, fclos, ierr

    ! Namelist variables
    integer :: deltmax ! maximum time difference (in minutes)

    namelist /thin_aircraft/deltmax

    ! Default values for namelist variables
    deltmax = 90

    ! Read the namelist for Aircraft observations (if it exists)
    if (utl_isNamelistPresent('thin_aircraft','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_thinAircraft: Error opening file flnml')
      read(nulnam,nml=thin_aircraft,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinAircraft: Error reading namelist')
      if (mpi_myid == 0) write(*,nml=thin_aircraft)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'thn_thinAircraft: Namelist block thin_aircraft is missing in the namelist.'
      write(*,*) '                  The default value will be taken.'
      if (mpi_myid == 0) write(*,nml=thin_aircraft)
    end if

    call thn_aircraftByBoxes(obsdat, 'AI', deltmax)

  end subroutine thn_thinAircraft

  !--------------------------------------------------------------------------
  ! thn_thinAladin
  !--------------------------------------------------------------------------
  subroutine thn_thinAladin(obsdat)
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat

    ! Locals:
    integer :: nulnam
    integer :: fnom, fclos, ierr

    ! Namelist variables
    integer :: keepNthVertical ! keep every nth vertical datum

    namelist /thin_aladin/keepNthVertical

    ! Default values for namelist variables
    keepNthVertical=-1

    ! Read the namelist for Aladin observations (if it exists)
    if (utl_isNamelistPresent('thin_aladin','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_thinAladin: Error opening file flnml')
      read(nulnam,nml=thin_aladin,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinAladin: Error reading namelist')
      if (mpi_myid == 0) write(*,nml=thin_aladin)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'thn_thinAladin: Namelist block thin_aladin is missing in the namelist.'
      write(*,*) '                The default value will be taken.'
      if (mpi_myid == 0) write(*,nml=thin_aladin)
    end if

    if (keepNthVertical > 0) then
      call thn_keepNthObs(obsdat, 'AL', keepNthVertical)
    end if

  end subroutine thn_thinAladin

  !--------------------------------------------------------------------------
  ! thn_thinTovs
  !--------------------------------------------------------------------------
  subroutine thn_thinTovs(obsdat)
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat

    ! Locals:
    integer :: nulnam
    integer :: fnom, fclos, ierr

    ! Namelist variables
    integer :: delta    ! 

    namelist /thin_tovs/delta

    ! Default namelist values
    delta = 100

    ! Read the namelist for TOVS observations (if it exists)
    if (utl_isNamelistPresent('thin_tovs','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_thinTovs: Error opening file flnml')
      read(nulnam,nml=thin_tovs,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinTovs: Error reading namelist')
      if (mpi_myid == 0) write(*,nml=thin_tovs)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'thn_thinTovs: Namelist block thin_tovs is missing in the namelist.'
      write(*,*) '              The default value will be taken.'
      if (mpi_myid == 0) write(*,nml=thin_tovs)
    end if

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_tovsFilt(obsdat, delta, codtyp_get_codtyp('amsua'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_tovsFilt(obsdat, delta, codtyp_get_codtyp('amsub'), &
                      codtyp2_opt=codtyp_get_codtyp('mhs'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_tovsFilt(obsdat, delta, codtyp_get_codtyp('atms'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_tovsFilt(obsdat, delta, codtyp_get_codtyp('mwhs2'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine thn_thinTovs

  !--------------------------------------------------------------------------
  ! thn_thinHyper
  !--------------------------------------------------------------------------
  subroutine thn_thinHyper(obsdat)
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat

    ! Locals:
    integer :: nulnam
    integer :: fnom, fclos, ierr

    ! Namelist variables
    logical :: removeUnCorrected ! indicate if obs without bias correction should be removed
    integer :: deltmax           ! time window by bin (from bin center to bin edge) (in minutes)
    integer :: deltax            ! thinning (dimension of box sides) (in km)
    integer :: deltrad           ! radius around box center for chosen obs (in km)
    namelist /thin_hyper/removeUnCorrected, deltmax, deltax, deltrad

    ! Default namelist values
    removeUnCorrected = .true.
    deltmax = 22
    deltax  = 150
    deltrad = 45

    ! Read the namelist for Aladin observations (if it exists)
    if (utl_isNamelistPresent('thin_hyper','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_thinHyper: Error opening file flnml')
      read(nulnam,nml=thin_hyper,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinHyper: Error reading namelist')
      if (mpi_myid == 0) write(*,nml=thin_hyper)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'thn_thinHyper: Namelist block thin_hyper is missing in the namelist.'
      write(*,*) '               The default value will be taken.'
      if (mpi_myid == 0) write(*,nml=thin_hyper)
    end if

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_thinByLatLonBoxes(obsdat, removeUnCorrected, deltmax, deltax, deltrad, &
                               'TO', codtyp_get_codtyp('airs'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_thinByLatLonBoxes(obsdat, removeUnCorrected, deltmax, deltax, deltrad, &
                               'TO', codtyp_get_codtyp('iasi'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_thinByLatLonBoxes(obsdat, removeUnCorrected, deltmax, deltax, deltrad, &
                               'TO', codtyp_get_codtyp('cris'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_thinByLatLonBoxes(obsdat, removeUnCorrected, deltmax, deltax, deltrad, &
                               'TO', codtyp_get_codtyp('crisfsr'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine thn_thinHyper

!_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
!_/
!_/ The following methods are intended to be general algorithms that may be
!_/ called by any of the observation-type-specific thinning methods.
!_/
!_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/


  !--------------------------------------------------------------------------
  ! thn_aircraftByBoxes
  !--------------------------------------------------------------------------
  subroutine thn_aircraftByBoxes(obsdat, familyType, deltmax)
    !
    ! :Purpose: Original method for thinning aircraft data by lat-lon boxes.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in)    :: familyType
    integer,          intent(in)    :: deltmax

    ! Locals:
    character(len=20) :: trlmFileName
    character(len=2)  :: fileNumber
    type(struct_hco), pointer :: hco_thinning
    type(struct_vco), pointer :: vco_sfc
    type(struct_gsv)          :: stateVectorPsfc
    integer :: numLon, numLat, nsize, headerIndexBeg, headerIndexEnd
    integer :: nulnam, ierr, lonIndex, latIndex, levIndex, stepIndex, codtyp
    integer :: obsLonIndex, obsLatIndex, obsLevIndex, obsStepIndex
    integer :: numHeader, numHeaderMaxMpi, headerIndex, bodyIndex
    integer :: aiTypeCount(4), aiTypeCountMpi(4)
    integer :: obsFlag, obsVarno, obsDate, obsTime
    logical :: ttMissing, huMissing, uuMissing, vvMissing, ddMissing, ffMissing
    real(8) :: zpresa, zpresb, obsPressure, delMinutes, obsStepIndex_r8
    real(8) :: obsLonInRad, obsLatInRad, obsLonInDegrees, obsLatInDegrees
    real(8) :: deltaLon, deltaLat, deltaPress, midDistance, score
    real(4), allocatable :: gridLat(:), gridLon(:)
    integer, allocatable :: rejectCount(:), rejectCountMpi(:)
    real(8), allocatable :: gridPressure(:,:,:,:)
    real(8), pointer     :: surfPressure(:,:,:,:)
    logical, allocatable :: keepObs(:), keepObsMpi(:), isAircraft(:)
    integer, allocatable :: obsLatIndexVec(:), obsLonIndexVec(:)
    integer, allocatable :: obsTimeIndexVec(:), obsLevIndexVec(:)
    integer, allocatable :: obsLatIndexMpi(:), obsLonIndexMpi(:)
    integer, allocatable :: obsTimeIndexMpi(:), obsLevIndexMpi(:)
    logical, allocatable :: obsUVPresent(:), obsTTPresent(:)
    logical, allocatable :: obsUVPresentMpi(:), obsTTPresentMpi(:)
    real(4), allocatable :: obsDistance(:), obsUU(:), obsVV(:), obsTT(:)
    real(4), allocatable :: obsDistanceMpi(:), obsUUMpi(:), obsVVMpi(:), obsTTMpi(:)
    integer, allocatable :: handlesGrid(:,:,:), numObsGrid(:,:,:)
    real(4), allocatable :: minScoreGrid(:,:,:), minDistGrid(:,:,:), maxDistGrid(:,:,:)
    real(4), allocatable :: uuSumGrid(:,:,:), vvSumGrid(:,:,:), ttSumGrid(:,:,:)

    integer, external :: fnom, fclos
    integer, parameter :: maxLev = 500

    ! Namelist variables:
    real(8) :: rprefinc
    real(8) :: rptopinc
    real(8) :: rcoefinc
    real(4) :: vlev(maxLev)
    integer :: numlev
    namelist /namgem/rprefinc, rptopinc, rcoefinc, numlev, vlev

    write(*,*) 'thn_aircraftByBoxes: Starting'

    numHeader = obs_numHeader(obsdat)
    call rpn_comm_allReduce(numHeader, numHeaderMaxMpi, 1, 'mpi_integer', &
                            'mpi_max','grid',ierr)
    write(*,*) 'thn_aircraftByBoxes: numHeader, numHeaderMaxMpi = ', &
         numHeader, numHeaderMaxMpi

    allocate(keepObs(numHeaderMaxMpi))
    allocate(isAircraft(numHeader))
    keepObs(:) = .false.
    aiTypeCount(:) = 0

    call obs_set_current_header_list(obsdat,trim(familyType))
    HEADER0: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER0

      ! check observation type
      codtyp = obs_headElem_i(obsdat, OBS_ITY, headerIndex)
      if ( codtyp == codtyp_get_codtyp('airep') ) then
        aiTypeCount(1) = aiTypeCount(1) + 1
        keepObs(headerIndex) = .true.
      else if ( codtyp == codtyp_get_codtyp('amdar')  ) then
        aiTypeCount(2) = aiTypeCount(2) + 1
        keepObs(headerIndex) = .true.
      else if ( codtyp == codtyp_get_codtyp('acars') ) then
        aiTypeCount(3) = aiTypeCount(3) + 1
        keepObs(headerIndex) = .true.
      else if ( codtyp == codtyp_get_codtyp('ads') ) then
        aiTypeCount(4) = aiTypeCount(4) + 1
        keepObs(headerIndex) = .true.
      end if
    end do HEADER0
    isAircraft(:) = keepObs(:)

    ! Return if no aircraft obs to thin
    allocate(keepObsMpi(numHeaderMaxMpi*mpi_nprocs))
    nsize = numHeaderMaxMpi
    call rpn_comm_allgather(keepObs,    nsize, 'mpi_logical',  &
                            keepObsMpi, nsize, 'mpi_logical', 'grid', ierr)
    if (count(keepObsMpi(:)) == 0) then
      write(*,*) 'thn_aircraftByBoxes: no aircraft observations present'
      return
    end if

    write(*,*) 'thn_aircraftByBoxes: number of *local* obs initial = ', count(keepObs(:))

    ! Setup horizontal thinning grid
    nullify(hco_thinning)
    call hco_SetupFromFile(hco_thinning, './analysisgrid_thinning_ai', &
                           'ANALYSIS', 'Analysis')

    ! Default namelist values
    numlev = 80
    vlev(:) = -1
    rprefinc = 0.0d0
    rptopinc = 0.0d0
    rcoefinc = 0.0d0
    ! Read the namelist defining the vertical levels
    if (utl_isNamelistPresent('namgem','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_aircraftByBoxes: Error opening file flnml')
      read(nulnam,nml=namgem,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_aircraftByBoxes: Error reading namelist')
      write(*,nml=namgem)
      ierr = fclos(nulnam)
    else
      call utl_abort('thn_aircraftByBoxes: Namelist block namgem is missing in the namelist.')
    end if

    ! Setup thinning grid parameters
    numLon = hco_thinning%ni
    numLat = hco_thinning%nj
    allocate(gridLat(numLat))
    allocate(gridLon(numLon))
    gridLon(:) = hco_thinning%lon(:) * MPC_DEGREES_PER_RADIAN_R8
    gridLat(:) = hco_thinning%lat(:) * MPC_DEGREES_PER_RADIAN_R8
    write(*,*) 'thinning grid lats = '
    write(*,*) gridLat(:)
    write(*,*) 'thinning grid lons = '
    write(*,*) gridLon(:)
    write(*,*) 'thinning grid vlev = '
    write(*,*) vlev(1:numLev)

    ! Allocate vectors
    allocate(rejectCount(tim_nstepObs))
    allocate(obsLatIndexVec(numHeaderMaxMpi))
    allocate(obsLonIndexVec(numHeaderMaxMpi))
    allocate(obsLevIndexVec(numHeaderMaxMpi))
    allocate(obsTimeIndexVec(numHeaderMaxMpi))
    allocate(obsDistance(numHeaderMaxMpi))
    allocate(obsUU(numHeaderMaxMpi))
    allocate(obsVV(numHeaderMaxMpi))
    allocate(obsTT(numHeaderMaxMpi))
    allocate(obsUVPresent(numHeaderMaxMpi))
    allocate(obsTTPresent(numHeaderMaxMpi))

    ! Allocate mpi global vectors
    allocate(rejectCountMpi(tim_nstepObs))
    allocate(obsLatIndexMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsLonIndexMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsLevIndexMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsTimeIndexMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsDistanceMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsUUMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsVVMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsTTMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsUVPresentMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsTTPresentMpi(numHeaderMaxMpi*mpi_nprocs))

    ! Initialize vectors
    rejectCount(:) = 0
    obsLatIndexVec(:) = 0
    obsLonIndexVec(:) = 0
    obsLevIndexVec(:) = 0
    obsTimeIndexVec(:) = 0
    obsDistance(:) = 0.
    obsUU(:) = 0.
    obsVV(:) = 0.
    obsTT(:) = 0.
    obsUVPresent(:) = .true.
    obsTTPresent(:) = .true.

    ! Read and interpolate the trial surface pressure
    nullify(vco_sfc)
    trlmFileName = './trlm_01'
    call vco_setupFromFile(vco_sfc, trlmFileName)
    call gsv_allocate( stateVectorPsfc, tim_nstepobs, hco_thinning, vco_sfc, &
                       datestamp_opt=tim_getDatestamp(), mpi_local_opt=.false., &
                       varNames_opt=(/'P0'/), hInterpolateDegree_opt='LINEAR' )
    do stepIndex = 1, tim_nstepobs
      write(fileNumber,'(I2.2)') stepIndex
      trlmFileName = './trlm_' // trim(fileNumber)
      call gsv_readFromFile( stateVectorPsfc, trlmFileName, ' ', ' ',  &
                             stepIndex_opt=stepIndex, containsFullField_opt=.true. )
    end do

    ! compute pressure of each model level using p0
    allocate(gridPressure(numLon,numLat,numLev,tim_nstepobs))
    gridPressure(:,:,:,:) = -1.0
    call gsv_getField(stateVectorPsfc,surfPressure)
    do stepIndex = 1, tim_nstepobs
      do levIndex  = 1, numLev
        zpresb = ( (vlev(levIndex) - rptopinc/rprefinc) /  &
                   (1.0D0-rptopinc/rprefinc) )**rcoefinc
        zpresa = rprefinc * (vlev(levIndex)-zpresb)
        do latIndex = 1, numLat
          do lonIndex = 1, numLon
            gridPressure(lonIndex,latIndex,levIndex,stepIndex) =  &
                 zpresa + zpresb*surfPressure(lonIndex,latIndex,1,stepIndex)
          end do
        end do
      end do
    end do
    call gsv_deallocate(stateVectorPsfc)

    call obs_set_current_header_list(obsdat,trim(familyType))
    HEADER1: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER1

      ! find time difference
      obsDate = obs_headElem_i(obsdat, OBS_DAT, headerIndex)
      obsTime = obs_headElem_i(obsdat, OBS_ETM, headerIndex)
      call tim_getStepObsIndex(obsStepIndex_r8, tim_getDatestamp(), &
                               obsDate, obsTime, tim_nstepobs)
      obsStepIndex = nint(obsStepIndex_r8)
      delMinutes = abs(nint(60.0 * tim_dstepobs * abs(real(obsStepIndex) - obsStepIndex_r8)))

      ! check time window
      if ( delMinutes > deltmax ) then
        keepObs(headerIndex) = .false.
        rejectCount(obsStepIndex) = rejectCount(obsStepIndex) + 1
      end if

      ! Only accept obs below 175hPa
      obsPressure = -1.0d0
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY
        
        if (obsPressure < 0.0d0) then
          obsPressure = obs_bodyElem_r(obsdat, OBS_PPP, bodyIndex)
          if ( obsPressure < 17500.0 .or. obsPressure > 110000.0 ) keepObs(headerIndex) = .false.
        end if
      end do BODY

      ttMissing = .true.
      huMissing = .true.
      uuMissing = .true.
      vvMissing = .true.
      ddMissing = .true.
      ffMissing = .true.      
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY2: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY2
        
        obsFlag  = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
        obsVarno = obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex)

        ! find number of elements availables
        if (obsVarno == BUFR_NETT) then
          if ( .not. (btest(obsFlag,18) .or. btest(obsFlag,16) .or. &
                      btest(obsFlag,9)  .or. btest(obsFlag,8)  .or. &
                      btest(obsFlag,2)) ) then
            ttMissing = .false.
            obsTT(headerIndex) = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex)
          end if
        else if (obsVarno == BUFR_NEES) then
          if ( .not. (btest(obsFlag,18) .or. btest(obsFlag,16) .or. &
                      btest(obsFlag,9)  .or. btest(obsFlag,8)  .or. &
                      btest(obsFlag,2)) ) huMissing = .false.
        else if (obsVarno == BUFR_NEUU) then
          if ( .not. (btest(obsFlag,18) .or. btest(obsFlag,16) .or. &
                      btest(obsFlag,9)  .or. btest(obsFlag,8)  .or. &
                      btest(obsFlag,2)) ) then
            uuMissing = .false.
            obsUU(headerIndex) = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex)
          end if
        else if (obsVarno == BUFR_NEVV) then
          if ( .not. (btest(obsFlag,18) .or. btest(obsFlag,16) .or. &
                      btest(obsFlag,9)  .or. btest(obsFlag,8)  .or. &
                      btest(obsFlag,2)) ) then
            vvMissing = .false.
            obsVV(headerIndex) = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex)
          end if
        else if (obsVarno == BUFR_NEDD) then
          if ( .not. (btest(obsFlag,18) .or. btest(obsFlag,16) .or. &
                      btest(obsFlag,9)  .or. btest(obsFlag,8)  .or. &
                      btest(obsFlag,2)) ) ddMissing = .false.
        else if (obsVarno == BUFR_NEFF) then
          if ( .not. (btest(obsFlag,18) .or. btest(obsFlag,16) .or. &
                      btest(obsFlag,9)  .or. btest(obsFlag,8)  .or. &
                      btest(obsFlag,2)) ) ffMissing = .false.
        end if

      end do BODY2

      ! wind components are rejected if speed or direction 
      if (ddMissing .or. ffMissing) then
        uuMissing = .true.
        vvMissing = .true.
      end if

      ! eliminate records with nothing to assimilate
      if (ttMissing .and. huMissing .and. uuMissing .and. vvMissing) then
        rejectCount(obsStepIndex) = rejectCount(obsStepIndex) + 1
        keepObs(headerIndex) = .false.
      end if

      if( uuMissing .or. vvMissing ) obsUVPresent(headerIndex) = .false.
      if( ttMissing )                obsTTPresent(headerIndex) = .false.

      ! Lat and Lon for each observation
      obsLonInRad = obs_headElem_r(obsdat, OBS_LON, headerIndex)
      obsLatInRad = obs_headElem_r(obsdat, OBS_LAT, headerIndex)
      obsLonInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obsLonInRad
      obsLatInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obsLatInRad

      ! latitude index
      deltaLat = abs(gridLat(1) - obsLatInDegrees)
      obsLatIndex = 1
      do latIndex = 2, numLat
        if (abs(gridLat(latIndex) - obsLatInDegrees) < deltaLat) then
          deltaLat = abs(gridLat(latIndex) - obsLatInDegrees)
          obsLatIndex = latIndex
        end if
      end do

      ! longitude index
      deltaLon = abs(gridLon(1) - obsLonInDegrees)
      obsLonIndex = 1
      do lonIndex = 2, numLon
        if ( abs(gridLon(lonIndex) - obsLonInDegrees) < deltaLon ) then
          deltaLon = abs(gridLon(lonIndex) - obsLonInDegrees)
          obsLonIndex = lonIndex
        end if
      end do

      ! layer index
      deltaPress = abs(gridPressure(obsLonIndex,obsLatIndex,1,obsStepIndex) - obsPressure)
      obsLevIndex = 1
      do levIndex = 2, numLev
        if ( abs(gridPressure(obsLonIndex,obsLatIndex,levIndex,obsStepIndex) - obsPressure) < deltaPress ) then
          deltaPress = abs(gridPressure(obsLonIndex,obsLatIndex,levIndex,obsStepIndex) - obsPressure)
          obsLevIndex = levIndex
        end if
      end do

      obsLatIndexVec(headerIndex) = obsLatIndex
      obsLonIndexVec(headerIndex) = obsLonIndex
      obsLevIndexVec(headerIndex) = obsLevIndex
      obsTimeIndexVec(headerIndex) = obsStepIndex
      obsDistance(headerIndex) = sqrt((deltaLon*3)**2 + (deltaLat*3)**2 + (deltaPress/100.0)**2)

    end do HEADER1

    call rpn_comm_allReduce(aiTypeCount, aiTypeCountMpi, 4, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(rejectCount, rejectCountMpi, tim_nstepObs, 'mpi_integer', &
                            'mpi_sum','grid',ierr)

    write(*,*)
    write(*,'(a50,i10)') ' Total number of obs = ', sum(aiTypeCountMpi(:))
    write(*,*)
    do stepIndex = 1, tim_nstepobs
      write(*,'(a50,2i10)')' Number of rejects for bin = ', stepIndex, rejectCountMpi(stepIndex)
    enddo
    write(*,'(a50,i10)')' Total number of rejects = ', sum(rejectCountMpi)
    write(*,*)
    write(*,'(a50,i10)') '====nb AIREP = ', aiTypeCountMpi(1)
    write(*,'(a50,i10)') '====nb AMDAR = ', aiTypeCountMpi(2)
    write(*,'(a50,i10)') '====nb ACARS = ', aiTypeCountMpi(3)
    write(*,'(a50,i10)') '====nb ADS = ',   aiTypeCountMpi(4)
    write(*,*)

    allocate(handlesGrid(numLat,numLon,numLev))
    allocate(minScoreGrid(numLat,numLon,numLev))
    allocate(minDistGrid(numLat,numLon,numLev))
    allocate(maxDistGrid(numLat,numLon,numLev))
    allocate(numObsGrid(numLat,numLon,numLev))
    allocate(uuSumGrid(numLat,numLon,numLev))
    allocate(vvSumGrid(numLat,numLon,numLev))
    allocate(ttSumGrid(numLat,numLon,numLev))

    ! Make all inputs to the following tests mpiglobal
    nsize = numHeaderMaxMpi
    call rpn_comm_allgather(keepObs,    nsize, 'mpi_logical',  &
                            keepObsMpi, nsize, 'mpi_logical', 'grid', ierr)
    call rpn_comm_allgather(obsLatIndexVec, nsize, 'mpi_integer',  &
                            obsLatIndexMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsLonIndexVec, nsize, 'mpi_integer',  &
                            obsLonIndexMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsLevIndexVec, nsize, 'mpi_integer',  &
                            obsLevIndexMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsTimeIndexVec,nsize, 'mpi_integer',  &
                            obsTimeIndexMpi,nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsDistance,    nsize, 'mpi_real4',  &
                            obsDistanceMpi, nsize, 'mpi_real4', 'grid', ierr)
    call rpn_comm_allgather(obsUU,    nsize, 'mpi_real4',  &
                            obsUUMpi, nsize, 'mpi_real4', 'grid', ierr)
    call rpn_comm_allgather(obsVV,    nsize, 'mpi_real4',  &
                            obsVVMpi, nsize, 'mpi_real4', 'grid', ierr)
    call rpn_comm_allgather(obsTT,    nsize, 'mpi_real4',  &
                            obsTTMpi, nsize, 'mpi_real4', 'grid', ierr)
    call rpn_comm_allgather(obsUVPresent,    nsize, 'mpi_logical',  &
                            obsUVPresentMpi, nsize, 'mpi_logical', 'grid', ierr)
    call rpn_comm_allgather(obsTTPresent,    nsize, 'mpi_logical',  &
                            obsTTPresentMpi, nsize, 'mpi_logical', 'grid', ierr)

    STEP: do stepIndex = 1, tim_nstepobs
      write (*,'(a50,i10)' ) ' Process bin number = ', stepIndex

      handlesGrid(:,:,:) = -1
      minScoreGrid(:,:,:) = 1000000.
      minDistGrid(:,:,:) = 1000000.
      maxDistGrid(:,:,:) = 0.
      numObsGrid(:,:,:) = 0 
      uuSumGrid(:,:,:) = 0. 
      vvSumGrid(:,:,:) = 0. 
      ttSumGrid(:,:,:) = 0. 

      ! Calcul des distances min et max du centre la boite des rapports 
      ! contenus dans les boites
      do headerIndex = 1, numHeaderMaxMpi*mpi_nprocs
        if( .not. keepObsMpi(headerIndex) ) cycle
        if( obsTimeIndexMpi(headerIndex) /= stepIndex ) cycle
        latIndex = obsLatIndexMpi(headerIndex)
        lonIndex = obsLonIndexMpi(headerIndex)
        levIndex = obsLevIndexMpi(headerIndex)
        if ( obsDistanceMpi(headerIndex) < minDistGrid(latIndex,lonIndex,levIndex) ) then
          minDistGrid(latIndex,lonIndex,levIndex) = obsDistanceMpi(headerIndex)
        end if
        if ( obsDistanceMpi(headerIndex) > maxDistGrid(latIndex,lonIndex,levIndex) ) then
          maxDistGrid(latIndex,lonIndex,levIndex) = obsDistanceMpi(headerIndex)
        end if
      end do

      ! Calcul des sommes de u, v et t des observations situees a une distance midDistance
      ! du centre de la boite
      do headerIndex = 1, numHeaderMaxMpi*mpi_nprocs
        if( .not. keepObsMpi(headerIndex) ) cycle
        if( obsTimeIndexMpi(headerIndex) /= stepIndex ) cycle
        latIndex = obsLatIndexMpi(headerIndex)
        lonIndex = obsLonIndexMpi(headerIndex)
        levIndex = obsLevIndexMpi(headerIndex)
        midDistance = ( minDistGrid(latIndex,lonIndex,levIndex) + &
                        maxDistGrid(latIndex,lonIndex,levIndex) )/2.
        if ( (obsDistanceMpi(headerIndex) < midDistance) .and. &
             obsTTPresentMpi(headerIndex) .and. &
             obsUVPresentMpi(headerIndex) ) then
          numObsGrid(latIndex,lonIndex,levIndex) =  &
               numObsGrid(latIndex,lonIndex,levIndex) + 1
          uuSumGrid(latIndex,lonIndex,levIndex) =  &
               uuSumGrid(latIndex,lonIndex,levIndex) + obsUUMpi(headerIndex)
          vvSumGrid(latIndex,lonIndex,levIndex) =  &
               vvSumGrid(latIndex,lonIndex,levIndex) + obsVVMpi(headerIndex)
          ttSumGrid(latIndex,lonIndex,levIndex) =  &
               ttSumGrid(latIndex,lonIndex,levIndex) + obsTTMpi(headerIndex)
        end if
      end do

      ! Calcul la moyenne de u, v et t s'il y a plus de 3 rapports dans la boite
      do latIndex = 1, numLat
        do lonIndex = 1, numLon
          do levIndex = 1, numLev
            if(numObsGrid(latIndex,lonIndex,levIndex) >= 3) then
              uuSumGrid(latIndex,lonIndex,levIndex) =  &
                   uuSumGrid(latIndex,lonIndex,levIndex)/numObsGrid(latIndex,lonIndex,levIndex)
              vvSumGrid(latIndex,lonIndex,levIndex) =  &
                   vvSumGrid(latIndex,lonIndex,levIndex)/numObsGrid(latIndex,lonIndex,levIndex)
              ttSumGrid(latIndex,lonIndex,levIndex) =  &
                   ttSumGrid(latIndex,lonIndex,levIndex)/numObsGrid(latIndex,lonIndex,levIndex)
            end if
          end do
        end do
      end do

      ! S'il y a plus de 3 rapports dans la boite, le rapport dont le score est le plus
      ! petit est retenu. S'il y a 2 rapports ou moins, le rapport le plus pres du centre
      ! de la boite est retenu.
      do headerIndex = 1, numHeaderMaxMpi*mpi_nprocs
        if( .not. keepObsMpi(headerIndex) ) cycle
        if( obsTimeIndexMpi(headerIndex) /= stepIndex ) cycle
        latIndex = obsLatIndexMpi(headerIndex)
        lonIndex = obsLonIndexMpi(headerIndex)
        levIndex = obsLevIndexMpi(headerIndex)

        if(numObsGrid(latIndex,lonIndex,levIndex) >= 3) then

          midDistance = ( minDistGrid(latIndex,lonIndex,levIndex) + &
                          maxDistGrid(latIndex,lonIndex,levIndex) )/2.
          if ((obsDistanceMpi(headerIndex) < midDistance) .and. &
               obsTTPresentMpi(headerIndex) .and. &
               obsUVPresentMpi(headerIndex) ) then
            score = sqrt( (uuSumGrid(latIndex,lonIndex,levIndex) - &
                           obsUUMpi(headerIndex))**2/(1.4**2) +   &
                          (vvSumGrid(latIndex,lonIndex,levIndex) - &
                           obsVVMpi(headerIndex))**2/(1.4**2) ) + &
                    (ttSumGrid(latIndex,lonIndex,levIndex) - &
                     obsTTMpi(headerIndex))**2/(0.9**2)

            if ( handlesGrid(latIndex,lonIndex,levIndex) /= -1 ) then
              if ( score >= minScoreGrid(latIndex,lonIndex,levIndex) ) then
                keepObsMpi(headerIndex) = .false.
              end if
            end if
          
            if ( keepObsMpi(headerIndex) ) then
              if ( handlesGrid(latIndex,lonIndex,levIndex) /= -1 ) then
                keepObsMpi(handlesGrid(latIndex,lonIndex,levIndex)) = .false.
              end if
              minScoreGrid(latIndex,lonIndex,levIndex) = score
              keepObsMpi(headerIndex) = .true.
              handlesGrid(latIndex,lonIndex,levIndex) = headerIndex
            end if

          else
            keepObsMpi(headerIndex) = .false.
          end if

        else ! if(numObsGrid(latIndex,lonIndex,levIndex) < 3)

          if ( handlesGrid(latIndex,lonIndex,levIndex) /= -1 ) then
            if ( obsDistanceMpi(headerIndex) > minDistGrid(latIndex,lonIndex,levIndex) ) then
              keepObsMpi(headerIndex) = .false.
            end if
          end if

          if ( keepObsMpi(headerIndex) ) then
            if ( handlesGrid(latIndex,lonIndex,levIndex) /= -1 ) then
              keepObsMpi(handlesGrid(latIndex,lonIndex,levIndex)) = .false.
            end if
            keepObsMpi(headerIndex) = .true.
            handlesGrid(latIndex,lonIndex,levIndex) = headerIndex
          end if

        end if

      end do ! headerIndex

    end do  STEP

    ! Update local copy of keepObs from global mpi version
    headerIndexBeg = 1 + mpi_myid * numHeaderMaxMpi
    headerIndexEnd = headerIndexBeg + numHeaderMaxMpi - 1
    keepObs(:) = keepObsMpi(headerIndexBeg:headerIndexEnd)

    deallocate(handlesGrid)
    deallocate(minScoreGrid)
    deallocate(minDistGrid)
    deallocate(maxDistGrid)
    deallocate(numObsGrid)
    deallocate(uuSumGrid)
    deallocate(vvSumGrid)
    deallocate(ttSumGrid)

    write (*,*)
    write (*,'(a50,i10)') " Number of obs in  = ", sum(aiTypeCountMpi(:))
    write (*,'(a50,i10)') " Number of obs out = ", count(keepObsMpi(:))
    write (*,'(a50,i10)') " Number of obs not out = ", &
         sum(aiTypeCountMpi(:)) - count(keepObsMpi(:))
    write (*,*)

    ! Modify the flags for rejected observations
    do headerIndex = 1, numHeader
      ! skip observation if we're not supposed to consider it
      if (.not. isAirCraft(headerIndex)) cycle
     
      if (.not. keepObs(headerIndex)) then
        call obs_set_current_body_list(obsdat, headerIndex)
        BODY3: do 
          bodyIndex = obs_getBodyIndex(obsdat)
          if (bodyIndex < 0) exit BODY3
        
          obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
          call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,11))

        end do BODY3
      end if
    end do

    ! Deallocation
    deallocate(keepObs)
    deallocate(isAircraft)
    deallocate(keepObsMpi)
    deallocate(gridLat)
    deallocate(gridLon)
    deallocate(rejectCount)
    deallocate(obsLatIndexVec)
    deallocate(obsLonIndexVec)
    deallocate(obsLevIndexVec)
    deallocate(obsTimeIndexVec)
    deallocate(obsDistance)
    deallocate(obsUU)
    deallocate(obsVV)
    deallocate(obsTT)
    deallocate(obsUVPresent)
    deallocate(obsTTPresent)
    deallocate(rejectCountMpi)
    deallocate(obsLatIndexMpi)
    deallocate(obsLonIndexMpi)
    deallocate(obsLevIndexMpi)
    deallocate(obsTimeIndexMpi)
    deallocate(obsDistanceMpi)
    deallocate(obsUUMpi)
    deallocate(obsVVMpi)
    deallocate(obsTTMpi)
    deallocate(obsUVPresentMpi)
    deallocate(obsTTPresentMpi)
    deallocate(gridPressure)

    write(*,*) 'thn_aircraftByBoxes: Finished'

  end subroutine thn_aircraftByBoxes

  !--------------------------------------------------------------------------
  ! thn_keepNthObs
  !--------------------------------------------------------------------------
  subroutine thn_keepNthObs(obsdat, familyType, keepNthVertical)
    !
    ! :Purpose: Of the observations in a column that have not already been
    !           rejected, keep every nth observation and throw out the rest.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in) :: familyType
    integer,          intent(in) :: keepNthVertical

    ! Locals:
    integer, parameter :: PROFILE_NOT_FOUND=-1
    integer :: headerIndex, bodyIndex, bodyIndex2, bodyIndexStart, bodyIndexEnd
    integer :: flag
    integer :: countKeepN ! count to keep every Nth observation in the column
    integer :: newProfileId

    write(*,*) 'thn_keepNthObs: Starting'
    countKeepN=0

    ! Loop over all body indices (columns) of the family of interest and
    ! thin each column independently of the others
    call obs_set_current_body_list(obsdat, familyType)
    BODY: do 
      bodyIndex = obs_getBodyIndex(obsdat)
      if (bodyIndex < 0) exit BODY

      ! If the datum is not being assimilated, ignore it
      if (obs_bodyElem_i(obsdat, OBS_ASS, bodyIndex) == obs_notAssimilated) cycle BODY

      ! If datum already rejected, ignore it
      flag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
      if (btest(flag,9)) cycle BODY

      headerIndex  = obs_bodyElem_i(obsdat, OBS_HIND, bodyIndex  )
      newProfileId = obs_headElem_i(obsdat, OBS_PRFL, headerIndex)

      countKeepN=countKeepN + 1
      if ( countKeepN == keepNthVertical .or. &
           new_column() ) then
        ! Reset the counter and keep this observation
        countKeepN=0
      else
        ! Reject this observation
        call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(flag,11))
      end if

    end do BODY

    write(*,*) 'thn_keepNthObs: Finished'

  contains
    function new_column()
      ! Determine whether the current observation begins a new vertical column
      ! (Assume that observations are in chronological order)
      !
      ! Note:  This method has been written with aladin observations in mind.
      !        It might be necessary to generalize the method.
      implicit none
      logical :: new_column

      integer, save :: previousProfileId=huge(previousProfileId)

      if (newProfileId == PROFILE_NOT_FOUND) then
        ! The profile ID for this element is missing.
        ! Assume that it is the same as the previous element
        newProfileId = previousProfileId
      end if

      if (newProfileId /= previousProfileId) then
        previousProfileId = newProfileId
        new_column=.true.
      else
        new_column=.false.
      end if
    end function new_column

  end subroutine thn_keepNthObs

  !--------------------------------------------------------------------------
  ! thn_tovsFilt
  !--------------------------------------------------------------------------
  subroutine thn_tovsFilt(obsdat, delta, codtyp, codtyp2_opt)
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    integer, intent(in)             :: delta
    integer, intent(in)             :: codtyp
    integer, optional, intent(in)   :: codtyp2_opt

    ! Locals:
    integer :: numLat, numLon, headerIndex, headerIndexKeep, latIndex, lonIndex, latIndex2
    integer :: gridIndex, ngrdtot, obsTime, obsDate, numHeader, numHeaderMaxMpi, ierr
    integer :: bodyIndex, stepIndex, obsIndex, obsFov
    integer :: loscan, hiscan, obsFlag, nbstn2, nsize
    integer :: procIndex, procIndexKeep, minLonBurpFile, obsCount, obsCountMpi
    integer :: countQc, countKept, countOther, countKeptMpi, countQcMpi, countGridPoints
    integer :: allMinLonBurpFile(mpi_nprocs)
    real(4) :: zlength, zlatrad, zdlon, zlatchk, zlonchk, zlontmp
    real(4) :: obsLatInRad, obsLonInRad, obsLat, obsLon
    real(4) :: obsLatInDegrees, obsLonInDegrees
    real(4) :: minDistance
    real(4) :: allMinDistance(mpi_nprocs)
    real(4) :: zdatflgs, zlatmid, zlonmid, zlatobs, zlonobs
    real(4) :: percentTotal, percentQc, percentOther, percentKept
    real(8), allocatable :: stepObsIndex(:)
    real(4), allocatable :: zlatdeg(:), zlatg(:), zlong(:), zdobs(:)
    logical, allocatable :: valid(:), indexs(:)
    integer, allocatable :: ngrd(:), nbstn(:), igrds(:)
    integer, allocatable :: obsLonBurpFile(:), obsLatBurpFile(:), nassim(:)
    integer, allocatable :: bufref(:), link(:)
    integer, allocatable :: headerIndexList(:), headerIndexList2(:)

    ! Local parameters:
    integer, parameter :: xlat=10000, xlon=40000, xdelt=50
    integer, parameter :: mxscanamsua=30
    integer, parameter :: mxscanamsub=90
    integer, parameter :: mxscanatms =96
    integer, parameter :: mxscanmwhs2=98

    write(*,*) 'thn_tovsFilt: Starting, ', trim(codtyp_get_name(codtyp))

    numHeader = obs_numHeader(obsdat)
    call rpn_comm_allReduce(numHeader, numHeaderMaxMpi, 1, 'mpi_integer', &
                            'mpi_max','grid',ierr)
    write(*,*) 'thn_tovsFilt: numHeader, numHeaderMaxMpi = ', numHeader, numHeaderMaxMpi

    ! Check if we have any observations to process
    allocate(valid(numHeaderMaxMpi))
    valid(:) = .false.
    do headerIndex = 1, numHeader
      if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) == codtyp) then
        valid(headerIndex) = .true.
      else if (present(codtyp2_opt)) then
        if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) == codtyp2_opt) then
          valid(headerIndex) = .true.
        end if
      end if
    end do
    if (count(valid(:)) == 0) then
      write(*,*) 'thn_tovsFilt: no observations for this instrument'
      return
    end if

    write(*,*) 'thn_tovsFilt: obsCount initial = ', count(valid(:))

    ! Remove RARS obs that are also present from a global originating centre
    call thn_removeRarsDuplicates(obsdat, valid)

    write(*,*) 'thn_tovsFilt: obsCount after thn_removeRarsDuplicates = ', count(valid(:))

    numLat = 2*xlat/delta
    numLon = xlon/delta

    ! Allocations
    allocate(zlatdeg(numLat))
    allocate(zlatg(numLat*numLon))
    allocate(zlong(numLat*numLon))
    allocate(nbstn(numLat*numLon))
    allocate(bufref(numLat*numLon))
    allocate(link(numHeader))
    allocate(headerIndexList(numHeader))
    allocate(headerIndexList2(numHeader))
    allocate(ngrd(numLat))
    allocate(igrds(numHeader))
    allocate(indexs(numHeader))
    allocate(obsLonBurpFile(numHeader))
    allocate(obsLatBurpFile(numHeader))
    allocate(nassim(numHeader))
    allocate(zdobs(numHeader))
    allocate(stepObsIndex(numHeader))

    ! Initialize arrays
    zlatdeg(:) = 0.0
    zlatg(:) = 0.0
    zlong(:) = 0.0
    zdobs(:) = 0.0
    nbstn(:) = 0
    bufref(:) = 0
    link(:) = 0
    headerIndexList(:) = 0
    headerIndexList2(:) = 0
    ngrd(:) = 0
    igrds(:) = 0
    indexs(:) = 0
    obsLonBurpFile(:) = 0
    obsLatBurpFile(:) = 0
    nassim(:) = 0
    stepObsIndex(:) = 0.0d0

    ! Set up the grid used for thinning
    ngrdtot = 0
    do latIndex = 1, numLat
      zlatdeg(latIndex) = (latIndex*180./numLat) - 90.
      zlatrad           = zlatdeg(latIndex) * mpc_pi_r4 / 180.
      zlength           = xlon * cos(zlatrad)
      ngrd(latIndex)    = nint(zlength/delta)
      ngrd(latIndex)    = max(ngrd(latIndex),1)
      ngrdtot           = ngrdtot + ngrd(latIndex)
    end do

    gridIndex = 0
    do latIndex = 1, numLat
      zdlon = 360./ngrd(latIndex)
      do lonIndex = 1, ngrd(latIndex)
        zlatg(gridIndex+lonIndex) = zlatdeg(latIndex)
        zlong(gridIndex+lonIndex) = (lonIndex-1)*zdlon
      end do
      gridIndex = gridIndex + ngrd(latIndex)
    end do

    ! Loop over all observation locations
    do headerIndex = 1, numHeader
      if ( .not. valid(headerIndex) ) cycle

      ! Lat and Lon for each observation
      obsLonInRad = obs_headElem_r(obsdat, OBS_LON, headerIndex)
      obsLatInRad = obs_headElem_r(obsdat, OBS_LAT, headerIndex)

      obsLonInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obsLonInRad
      obsLatInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obsLatInRad
      obsLonBurpFile(headerIndex) = nint(100.0*(obsLonInDegrees - 180.0))
      if(obsLonBurpFile(headerIndex) < 0) then
        obsLonBurpFile(headerIndex) = obsLonBurpFile(headerIndex) + 36000
      end if
      obsLatBurpFile(headerIndex) = 9000+nint(100.0*obsLatInDegrees)

      obsDate = obs_headElem_i(obsdat, OBS_DAT, headerIndex)
      obsTime = obs_headElem_i(obsdat, OBS_ETM, headerIndex)
      call tim_getStepObsIndex(stepObsIndex(headerIndex), tim_getDatestamp(), &
                               obsDate, obsTime, tim_nstepobs)

      ! Associate each observation to a grid point
      zlatchk = (obsLatBurpFile(headerIndex) - 9000.) / 100.
      zlonchk = obsLonBurpFile(headerIndex) / 100.
      zlontmp = zlonchk
      if (zlonchk > 180.) zlonchk = zlonchk - 360.
      do latIndex = 1,numLat-1
        if (zlatchk <  zlatdeg(1))     zlatchk = zlatdeg(1)
        if (zlatchk >= zlatdeg(numLat)) zlatchk = zlatdeg(numLat) - 0.5
        if (zlatchk >= zlatdeg(latIndex) .and. zlatchk < zlatdeg(latIndex+1)) then
          gridIndex = 1
          do latIndex2 = 1, latIndex-1
            gridIndex = gridIndex + ngrd(latIndex2)
          end do
          zdlon = 360./ngrd(latIndex)
          gridIndex = gridIndex + ifix(zlontmp/zdlon)
          exit
        end if
      end do
      igrds(headerIndex) = gridIndex
      nbstn(gridIndex) = nbstn(gridIndex) + 1

    end do ! headerIndex

    if      ( codtyp == codtyp_get_codtyp('amsua') ) then
      loscan   = 1
      hiscan   = mxscanamsua
    else if ( codtyp == codtyp_get_codtyp('amsub') ) then
      loscan   = 1
      hiscan   = mxscanamsub
    else if ( codtyp == codtyp_get_codtyp('atms') ) then
      loscan   = 2
      hiscan   = mxscanatms - 1
    else if ( codtyp == codtyp_get_codtyp('mwhs2') ) then
      loscan   = 1
      hiscan   = mxscanmwhs2
    end if

    countQc = 0
    nassim(:) = 0
    do headerIndex = 1, numHeader

      if ( .not. valid(headerIndex) ) cycle

      ! Look at the obs flags
      zdatflgs = 0.

      call obs_set_current_body_list(obsdat, headerIndex)
      BODY: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY
        
        ! If not a blacklisted channel (note that bit 11 is set in 
        ! satqc_amsu*.f for blacklisted channels)
        obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
        if ( .not. btest(obsFlag,11) ) then
          nassim(headerIndex) = nassim(headerIndex) + 1
          if ( btest(obsFlag,9) ) then
            zdatflgs = zdatflgs + 1.0
          end if
        end if
      end do BODY

      ! fixer le % de rejets a 100% si aucun canal n'est assimilable         
      if ( zdatflgs == 0. .and. nassim(headerIndex) == 0 ) then
        zdatflgs = 1.
      else
        zdatflgs = zdatflgs / max(nassim(headerIndex),1)  
      end if

      obsFov = obs_headElem_i(obsdat, OBS_FOV, headerIndex)
      if ( zdatflgs >= 0.80 ) then
        countQc = countQc + 1
        valid(headerIndex) = .false.
      else if (obsFov < loscan .or.  &
               obsFov > hiscan ) then
        countQc = countQc + 1
        valid(headerIndex) = .false.
      end if

    end do

    write(*,*) 'thn_tovsFilt: obsCount after QC = ', count(valid(:))

    ! Calculate distance of obs. from center of its grid box center
    do headerIndex = 1, numHeader
      if ( .not. valid(headerIndex) ) cycle

      gridIndex = igrds(headerIndex)
      if (nbstn(gridIndex) /= 0) then
        latIndex = (zlatg(gridIndex)+90.)/(180./numLat)
        zdlon = 360./ngrd(latIndex)
        zlatmid = zlatg(gridIndex) + 0.5*(180./numLat)
        zlonmid = zlong(gridIndex) + 0.5*zdlon
        zlatobs = (obsLatBurpFile(headerIndex) - 9000.) / 100.
        zlonobs = obsLonBurpFile(headerIndex) / 100.
        zdobs(headerIndex) = separa(zlonobs,zlatobs,zlonmid,zlatmid) * &
             float(xlat) / 90.
      end if
    end do

    ! Create a linked list of observations (link to grid point)
    bufref(:) = 0
    link(:) = 0
    obsCount = 0
    do headerIndex = 1, numHeader
      if ( .not. valid(headerIndex) ) cycle

      gridIndex = igrds(headerIndex)
      if (nbstn(gridIndex) /= 0) then
        obsCount = obsCount + 1
        headerIndexList(obsCount) = headerIndex
        link(obsCount) = bufref(gridIndex)
        bufref(gridIndex) = obsCount
      end if
    end do

    ! Loop over stepObs
    do stepIndex = 1, tim_nstepobs

      ! Loop over all grid points
      countGridPoints = 0
      do gridIndex = 1, ngrdtot
        if (nbstn(gridIndex) /= 0) then
          countGridPoints = countGridPoints + 1
        end if
        nbstn2 = 0
        obsIndex = bufref(gridIndex)
        do
          if (obsIndex == 0) exit
          headerIndex = headerIndexList(obsIndex)
          if ( igrds(headerIndex) == gridIndex  .and. &
               valid(headerIndex)               .and. &
               nint(stepObsIndex(headerIndex)) == stepIndex ) then
            nbstn2 = nbstn2 + 1
            headerIndexList2(nbstn2) = headerIndex
          end if

          obsIndex = link(obsIndex)
        end do

        minDistance = 1000000.             

        ! Choose the obs closest to the grid point
        do obsIndex = 1, nbstn2
          if (zdobs(headerIndexList2(obsIndex)) < minDistance) then
            minDistance = zdobs(headerIndexList2(obsIndex))
            minLonBurpFile = obsLonBurpFile(headerIndexList2(obsIndex))
            headerIndexKeep = headerIndexList2(obsIndex)
          end if
        end do

        ! Check for multiple obs with same distance to grid point
        if (nbstn2 > 0) then
          if ( count(zdobs(headerIndexList2(1:nbstn2)) == minDistance) > 1 ) then
            ! resolve ambiguity by choosing obs with min value of lon
            minLonBurpFile = 10000000
            do obsIndex = 1, nbstn2
              if (zdobs(headerIndexList2(obsIndex)) == minDistance) then
                if (obsLonBurpFile(headerIndexList2(obsIndex)) < minLonBurpFile) then
                  minLonBurpFile = obsLonBurpFile(headerIndexList2(obsIndex))
                  headerIndexKeep = headerIndexList2(obsIndex)
                end if
              end if
            end do
          end if
        end if

        do obsIndex = 1, nbstn2
          valid(headerIndexList2(obsIndex)) = .false.
        end do
        if (nbstn2 > 0 .and. minDistance <= 75. ) then
          valid(headerIndexKeep) = .true.
        end if

        ! Communicate the distance of chosen observation among all mpi tasks
        call rpn_comm_allgather(minDistance,    1, 'mpi_real4',  &
                                allMinDistance, 1, 'mpi_real4', 'grid', ierr)

        ! Choose the closest to the center of the box among all mpi tasks
        minDistance = 1000000.
        do procIndex = 1, mpi_nprocs
          if (allMinDistance(procIndex) < minDistance) then
            minDistance = allMinDistance(procIndex)
            procIndexKeep = procIndex
          end if
        end do

        ! Adjust flags to only keep 1 observation among all mpi tasks
        if (minDistance < 1000000.) then
          if ( count(allMinDistance(:) == minDistance) > 1 ) then
            ! resolve ambiguity by choosing obs with min value of lon
            call rpn_comm_allgather(minLonBurpFile,    1, 'mpi_integer',  &
                                    allMinLonBurpFile, 1, 'mpi_integer', 'grid', ierr)
            minLonBurpFile = 10000000
            do procIndex = 1, mpi_nprocs
              if (allMinDistance(procIndex) == minDistance) then
                if (allMinLonBurpFile(procIndex) < minLonBurpFile) then
                  minLonBurpFile = allMinLonBurpFile(procIndex)
                  procIndexKeep = procIndex
                end if
              end if
            end do            
          end if
          if (nbstn2 > 0) then
            if (mpi_myid /= (procIndexKeep-1)) then
              valid(headerIndexKeep) = .false.
            end if
          end if
        end if

      end do ! gridIndex
    end do ! stepIndex

    write(*,*) 'thn_tovsFilt: obsCount after thinning = ', count(valid(:))

    ! modify the observation flags in obsSpaceData
    obsCount = 0
    do headerIndex = 1, numHeader
      ! skip observation if we're not supposed to consider it
      if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) /= codtyp) then
        if (present(codtyp2_opt)) then
          if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) /= codtyp2_opt) then
            cycle
          end if
        else
          cycle
        end if
      end if
     
      obsCount = obsCount + 1

      if (.not. valid(headerIndex)) then
        call obs_set_current_body_list(obsdat, headerIndex)
        BODY2: do 
          bodyIndex = obs_getBodyIndex(obsdat)
          if (bodyIndex < 0) exit BODY2
        
          obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
          call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,11))

        end do BODY2
      end if
    end do

    ! print a summary to the listing
    countKept = count(valid)
    call rpn_comm_allReduce(countKept, countKeptMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(countQc,   countQcMpi,   1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(obsCount, obsCountMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)

    countOther = obsCountMpi - countKeptMpi - countQcMpi

    percentTotal = 100.0
    percentQc    = (float(countQcMpi)   / float(obsCountMpi)) * 100.0
    percentOther = (float(countOther)   / float(obsCountMpi)) * 100.0
    percentKept  = (float(countKeptMpi) / float(obsCountMpi)) * 100.0
         
    write(*,100)
100 format(/,' SOMMAIRE DES RESULTATS',/)
    write(*,200) obsCountMpi, percentTotal, countQcMpi, percentQc, &
                 countOther, percentOther, countKeptMpi, percentKept, &
                 delta, ngrdtot, countGridPoints
200 format(' NB.STNS TOTAL AU DEBUT EN ENTREE        =    ',I7, &
           '  =  ',F6.1,'%',/, &
           ' NB.STNS REJETEES AU CONTROLE QUALITATIF =  - ',I7, &
           '  =  ',F6.1,'%',/, &
           ' NB.STNS REJETEES POUR LES AUTRES RAISONS=  - ',I7, &
           '  =  ',F6.1,'%',/, &
           '                                              ', &
           '-------- = -------',/, &
           ' NB.STNS GARDEES A LA FIN                =    ',I7, &
           '  =  ',F6.1,'%',/////, &
           ' NB.PTS GRILLE AVEC LA RESOLUTION ',I5,' KM   =    ', &
           I7,/, &
           ' NB.PTS GRILLE TRAITES                       =    ',I7)

    ! end of sommair

    deallocate(valid)
    deallocate(zlatdeg)
    deallocate(zlatg)
    deallocate(zlong)
    deallocate(nbstn)
    deallocate(bufref)
    deallocate(link)
    deallocate(headerIndexList)
    deallocate(headerIndexList2)
    deallocate(ngrd)
    deallocate(igrds)
    deallocate(indexs)
    deallocate(obsLonBurpFile)
    deallocate(obsLatBurpFile)
    deallocate(nassim)
    deallocate(zdobs)
    deallocate(stepObsIndex)

    write(*,*) 'thn_tovsFilt: finished'

  end subroutine thn_tovsFilt

  !--------------------------------------------------------------------------
  ! thn_removeRarsDuplicates
  !--------------------------------------------------------------------------
  subroutine thn_removeRarsDuplicates(obsdat, valid)
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    logical,          intent(inout) :: valid(:)

    ! Locals:
    integer :: nsize, ierr, lenStnId, headerIndex, headerIndex1, headerIndex2
    integer :: numHeader, numHeaderMaxMpi, charIndex, headerIndexBeg, headerIndexEnd
    integer :: obsDate, obsTime
    real(4) :: obsLatInRad, obsLonInRad
    real(8) :: dlhours
    logical :: global1, global2
    integer, allocatable :: centreOrig(:), allCentreOrig(:)
    integer, allocatable :: obsFov(:), allObsFov(:)
    integer, allocatable :: obsDateStamp(:), allObsDateStamp(:)
    integer, allocatable :: stnIdInt(:,:), allStnIdInt(:,:)
    logical, allocatable :: allValid(:)
    character(len=12)    :: stnId
    
    ! Locals related to kdtree2:
    type(kdtree2), pointer            :: tree
    integer, parameter                :: maxNumSearch = 100
    integer                           :: numFoundSearch, resultIndex
    type(kdtree2_result)              :: searchResults(maxNumSearch)
    real(kdkind)                      :: maxRadius = 100.d6
    real(kdkind)                      :: refPosition(3)
    real(kdkind), allocatable         :: obsPosition3d(:,:)
    real(kdkind), allocatable         :: allObsPosition3d(:,:)

    ! Local parameters:
    integer, parameter :: centreOrigGlobal(3)=(/53, 74, 160/)

    ! Externals:
    integer, external    :: newdate

    numHeader = obs_numHeader(obsdat)
    call rpn_comm_allReduce(numHeader, numHeaderMaxMpi, 1, 'mpi_integer', &
                            'mpi_max','grid',ierr)

    ! Allocations
    allocate(obsPosition3d(3,numHeaderMaxMpi))
    allocate(allObsPosition3d(3,numHeaderMaxMpi*mpi_nprocs))
    allocate(centreOrig(numHeaderMaxMpi))
    allocate(allCentreOrig(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsFov(numHeaderMaxMpi))
    allocate(allObsFov(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsDateStamp(numHeaderMaxMpi))
    allocate(allObsDateStamp(numHeaderMaxMpi*mpi_nprocs))
    allocate(allValid(numHeaderMaxMpi*mpi_nprocs))
    lenStnId = len(stnId)
    allocate(stnIdInt(lenStnId,numHeaderMaxMpi))
    allocate(allStnIdInt(lenStnId,numHeaderMaxMpi*mpi_nprocs))

    ! Some initializations
    centreOrig(:) = 0
    obsPosition3d(:,:) = 0.0

    ! Loop over all observation locations
    do headerIndex = 1, numHeader
      if ( .not. valid(headerIndex) ) cycle

      ! Originating centre of data
      centreOrig(headerIndex) = obs_headElem_i(obsdat, OBS_ORI, headerIndex)

      ! Station ID converted to integer array
      stnId = obs_elem_c(obsdat,'STID',headerIndex)
      do charIndex = 1, lenStnId
        stnIdInt(charIndex,headerIndex) = iachar(stnId(charIndex:charIndex))
      end do

      ! Date stamp for each observation
      obsDate = obs_headElem_i(obsdat, OBS_DAT, headerIndex)
      obsTime = obs_headElem_i(obsdat, OBS_ETM, headerIndex)
      ierr = newdate(obsDateStamp(headerIndex), obsDate, obsTime*10000+2900, 3)

      ! Field of View for each observation
      obsFov(headerIndex) = obs_headElem_i(obsdat, OBS_FOV, headerIndex)

      ! Lat and Lon for each observation
      obsLonInRad = obs_headElem_r(obsdat, OBS_LON, headerIndex)
      obsLatInRad = obs_headElem_r(obsdat, OBS_LAT, headerIndex)

      ! 3D location array for kdtree
      obsPosition3d(1,headerIndex) = RA * sin(obsLonInRad) * cos(obsLatInRad)
      obsPosition3d(2,headerIndex) = RA * cos(obsLonInRad) * cos(obsLatInRad)
      obsPosition3d(3,headerIndex) = RA *                    sin(obsLatInRad)
    end do

    nsize = 3 * numHeaderMaxMpi
    call rpn_comm_allgather(obsPosition3d,    nsize, 'mpi_real8',  &
                            allObsPosition3d, nsize, 'mpi_real8', 'grid', ierr)
    nsize = numHeaderMaxMpi
    call rpn_comm_allgather(valid,    nsize, 'mpi_logical',  &
                            allValid, nsize, 'mpi_logical', 'grid', ierr)
    call rpn_comm_allgather(centreOrig,    nsize, 'mpi_integer',  &
                            allCentreOrig, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsFov,    nsize, 'mpi_integer',  &
                            allObsFov, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsDateStamp,    nsize, 'mpi_integer',  &
                            allObsDateStamp, nsize, 'mpi_integer', 'grid', ierr)
    nsize = lenStnId * numHeaderMaxMpi
    call rpn_comm_allgather(stnIdInt,    nsize, 'mpi_integer',  &
                            allStnIdInt, nsize, 'mpi_integer', 'grid', ierr)
    nullify(tree)

    tree => kdtree2_create(allObsPosition3d, sort=.true., rearrange=.true.)
    HEADER1: do headerIndex1 = 1, mpi_nprocs*numHeaderMaxMpi
        
      if ( .not. allValid(headerIndex1) ) cycle HEADER1

      ! Find all obs within 10km
      refPosition(:) = allObsPosition3d(:,headerIndex1)
      call kdtree2_r_nearest(tp=tree, qv=refPosition, r2=maxRadius, nfound=numFoundSearch, &
                             nalloc=maxNumSearch, results=searchResults)
      if (numFoundSearch >= maxNumSearch) then
        call utl_abort('thn_tovsFilt: the parameter maxNumSearch must be increased')
      end if
      if (numFoundSearch == 0) then
        call utl_abort('thn_tovsFilt: no match found. This should not happen!!!')
      end if

      ! Loop over all of these nearby locations
      HEADER2: do resultIndex = 1, numFoundSearch
        headerIndex2 = searchResults(resultIndex)%idx

        if ( .not. allValid(headerIndex2) ) cycle HEADER2

        ! Certaines stations locales nous envoient 
        ! le mauvais numero d'orbite. On ne peut donc
        ! pas s'y fier.
        ! Il faut comparer le temps de la reception des
        ! donnees
        if ( allCentreOrig(headerIndex1) /= allCentreOrig(headerIndex2) ) then
          if ( allObsFov(headerIndex1) == allObsFov(headerIndex2) ) then
            if ( all(allStnIdInt(:,headerIndex1) ==  allStnIdInt(:,headerIndex2)) ) then
            
              ! Difference (in hours) between obs time
              call difdatr(allObsDateStamp(headerIndex1),allObsDateStamp(headerIndex2),dlhours)

              ! Si la difference est moins de 6 minutes,
              ! on peut avoir affaire a un rars

              if ( abs(dlhours) <= 0.1 ) then

                ! si l'element_i est global, on doit le garder et rejeter l'element_j
                global1 = any(centreOrigGlobal(:) == allCentreOrig(headerIndex1))
                if (global1) then 
                  allValid(headerIndex2) = .false.
                else
                  ! toutefois, ca ne signifie pas que l'element_j est un rars
                  ! VERIFIER SI LA STATION 2 EST RARS
                  global2 = any(centreOrigGlobal(:) == allCentreOrig(headerIndex2))

                  ! Si l'element_j est global, rejeter l'element_i
                  ! Si les 2 elements sont rars, garder le 1er
                  if (global2) then 
                    allValid(headerIndex1) = .false.
                    cycle HEADER1
                  else
                    allValid(headerIndex2) = .false.
                  end if
                end if

              end if ! abs(dlhours) <= 0.1
              
            end if ! STID1 == STID2
          end if ! FOV1 == FOV2
        end if ! centreOrig1 /= centreOrig2

      end do HEADER2
    end do HEADER1
    call kdtree2_destroy(tree)

    ! update local copy of 'valid' array
    headerIndexBeg = 1 + mpi_myid * numHeaderMaxMpi
    headerIndexEnd = headerIndexBeg + numHeaderMaxMpi - 1
    valid(:) = allValid(headerIndexBeg:headerIndexEnd)

    deallocate(obsPosition3d)
    deallocate(allObsPosition3d)
    deallocate(centreOrig)
    deallocate(allCentreOrig)
    deallocate(obsFov)
    deallocate(allObsFov)
    deallocate(obsDateStamp)
    deallocate(allObsDateStamp)
    deallocate(allValid)
    deallocate(stnIdInt)
    deallocate(allStnIdInt)

  end subroutine thn_removeRarsDuplicates

  function separa(xlon1,xlat1,xlon2,xlat2)
    implicit none

    ! Arguments:
    real(4) :: separa
    real(4) :: xlat1, xlat2, xlon1, xlon2

    ! Locals:
    real(4) :: cosval, degrad, raddeg

    raddeg = 180.0/3.14159265358979
    degrad = 1.0/raddeg
    cosval = sin(xlat1*degrad) * sin(xlat2*degrad) + &
             cos(xlat1*degrad) * cos(xlat2*degrad) * &
             cos((xlon1-xlon2) * degrad)

    if (cosval < -1.0d0) then
      cosval = -1.0d0
    else if (cosval > 1.0d0) then
      cosval = 1.0d0
    end if
    separa = acos(cosval) * raddeg

  end function separa


  !--------------------------------------------------------------------------
  ! thn_thinByLatLonBoxes
  !--------------------------------------------------------------------------
  subroutine thn_thinByLatLonBoxes(obsdat, removeUnCorrected, &
                                   deltmax, deltax, deltrad,  &
                                   familyType, codtyp)
    !
    ! :Purpose: Only keep the observation closest to the center of each
    !           lat-lon (and time) box.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    logical, intent(in)             :: removeUnCorrected
    integer, intent(in)             :: deltmax
    integer, intent(in)             :: deltax
    integer, intent(in)             :: deltrad
    character(len=*), intent(in)    :: familyType
    integer, intent(in)             :: codtyp

    ! Locals:
    integer :: headerIndex, bodyIndex, obsDate, obsTime, obsFlag
    integer :: numLat, numLon, latIndex, numChannels, delMinutes
    integer :: lonBinIndex, latBinIndex, timeBinIndex
    integer :: ierr, nsize, procIndex, countHeader
    real(4) :: latr, length, distance
    real(8) :: lonBoxCenterInDegrees, latBoxCenterInDegrees
    real(8) :: obsLatInRad, obsLonInRad, obsLat, obsLon
    real(8) :: obsLatInDegrees, obsLonInDegrees, stepObsIndex
    integer, parameter :: LAT_LENGTH = 10000
    integer, parameter :: LON_LENGTH = 40000
    real(4), allocatable :: latdeg(:)
    integer, allocatable :: ngrd(:)
    integer, allocatable :: headerIndexKeep(:,:,:), numChannelsKeep(:,:,:)
    integer, allocatable :: allHeaderIndex(:,:,:,:), allNumChannels(:,:,:,:)
    integer, allocatable :: delMinutesKeep(:,:,:)
    integer, allocatable :: allDelMinutes(:,:,:,:)
    integer, allocatable :: procIndexKeep(:,:,:)
    real(4), allocatable :: distanceKeep(:,:,:)
    real(4), allocatable :: allDistance(:,:,:,:)
    logical, allocatable :: rejectThisHeader(:)
    logical :: keepThisObs
    integer :: obsLonBurpFile, obsLatBurpFile
    character(len=9) :: stnid

    write(*,*) 'thn_thinByLatLonBoxes: Starting, ', trim(codtyp_get_name(codtyp))

    ! Initial setup
    numLat = nint( 2. * real(LAT_LENGTH) / real(deltax) )
    numLon = nint(      real(LON_LENGTH) / real(deltax) )
    allocate(headerIndexKeep(numLat,numLon,tim_nstepobs))
    allocate(numChannelsKeep(numLat,numLon,tim_nstepobs))
    allocate(distanceKeep(numLat,numLon,tim_nstepobs))
    allocate(delMinutesKeep(numLat,numLon,tim_nstepobs))
    allocate(latdeg(numLat))
    allocate(ngrd(numLat))
    headerIndexKeep(:,:,:) = -1
    numChannelsKeep(:,:,:) = 0
    distanceKeep(:,:,:)    = 0.0
    delMinutesKeep(:,:,:)  = deltmax
    latdeg(:)              = 0.0
    ngrd(:)                = 0

    allocate(allHeaderIndex(numLat,numLon,tim_nstepobs,mpi_nprocs))
    allocate(allNumChannels(numLat,numLon,tim_nstepobs,mpi_nprocs))
    allocate(allDistance(numLat,numLon,tim_nstepobs,mpi_nprocs))
    allocate(allDelMinutes(numLat,numLon,tim_nstepobs,mpi_nprocs))
    allocate(procIndexKeep(numLat,numLon,tim_nstepobs))
    procIndexKeep(:,:,:) = -1

    ! set spatial boxes properties
    ! latdeg(:) : latitude (deg) of northern side of the box
    ! ngrd(:)   : number of longitudinal boxes at this latitude
    do latIndex = 1, numLat
      latdeg(latIndex) = (latIndex*180./numLat) - 90.
      if ( latdeg(latIndex) <= 0.0 ) then
        latr = latdeg(latIndex) * MPC_PI_R8 / 180.
      else
        latr = latdeg(latIndex-1) * MPC_PI_R8 / 180.
      end if
      length = LON_LENGTH * cos(latr)
      ngrd(latIndex)   = nint(length/deltax)
    end do

    countHeader = 0

    ! loop over all header indices of the specified family
    call obs_set_current_header_list(obsdat,trim(familyType))
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER

      if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) /= codtyp) cycle HEADER

      countHeader = countHeader + 1

      obsLonInRad = obs_headElem_r(obsdat, OBS_LON, headerIndex)
      obsLatInRad = obs_headElem_r(obsdat, OBS_LAT, headerIndex)
      obsDate = obs_headElem_i(obsdat, OBS_DAT, headerIndex)
      obsTime = obs_headElem_i(obsdat, OBS_ETM, headerIndex)

      obsLonInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obsLonInRad
      obsLatInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obsLatInRad
      obsLonBurpFile = nint(100.0*(obsLonInDegrees - 180.0))
      if(obsLonBurpFile < 0) obsLonBurpFile = obsLonBurpFile + 36000
      obsLatBurpFile = 9000+nint(100.0*obsLatInDegrees)

      numChannels = 0

      ! loop over all body indices for this headerIndex
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY

        ! mark for rejection if not bias corrected (bit 6 not set)
        if (removeUnCorrected) then
          obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
          if (.not. btest(obsFlag,6)) then
            call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,11))
            call obs_bodySet_i(obsdat, OBS_ASS, bodyIndex, obs_notAssimilated)
          end if
        end if

        ! count the number of accepted channels
        if (obs_bodyElem_i(obsdat, OBS_ASS, bodyIndex) == obs_assimilated) then
          numChannels = numChannels + 1
        end if

      end do BODY

      ! Determine the lat and lon bin indexes
      do latIndex = 1, numLat
        if ( obsLatInDegrees <= (latdeg(latIndex)+0.000001) ) then
          latBinIndex = latIndex
          exit
        end if
      end do
      lonBinIndex = int( obsLonBurpFile/(36000.0/ngrd(latBinIndex)) ) + 1
      if ( lonBinIndex > ngrd(latBinIndex) ) lonBinIndex = ngrd(latBinIndex)

      ! Determine the time bin index
      call tim_getStepObsIndex(stepObsIndex, tim_getDatestamp(), &
                               obsDate, obsTime, tim_nstepobs)
      timeBinIndex = nint(stepObsIndex)
      delMinutes = nint(60.0 * tim_dstepobs * abs(real(timeBinIndex) - stepObsIndex))

      ! Determine distance from box center
      latBoxCenterInDegrees = latdeg(latBinIndex) - 0.5 * (180./numLat)
      lonBoxCenterInDegrees = (360. / ngrd(latBinIndex)) * (lonBinIndex - 0.5)
      obsLat = (obsLatBurpFile - 9000.) / 100.
      obsLon = obsLonBurpFile / 100.
      distance = 1.0d-3 * phf_calcDistance(MPC_RADIANS_PER_DEGREE_R8 * latBoxCenterInDegrees, &
                                           MPC_RADIANS_PER_DEGREE_R8 * lonBoxCenterInDegrees, &
                                           MPC_RADIANS_PER_DEGREE_R8 * obsLat, &
                                           MPC_RADIANS_PER_DEGREE_R8 * obsLon )

      stnid = obs_elem_c(obsdat,'STID',headerIndex)

      ! Apply thinning criteria
      keepThisObs = .false.

      ! keep if distance to box center smaller than limit and 
      ! time under limit fixed at input and maximise number of channels 
      if ( (distance < deltrad)                                              .and. &
           (numChannels >= numChannelsKeep(latBinIndex,lonBinIndex,timeBinIndex)) .and. &
           (delMinutes <= deltmax) ) keepThisObs = .true.
     
      ! keep the closest to bin central time
      if ( numChannels == numChannelsKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
        if ( headerIndexKeep(latBinIndex,lonBinIndex,timeBinIndex) == -1 ) then
          keepThisObs = .false.
        else
          if ( delMinutes > delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
            keepThisObs = .false.
          end if
         
          if ( delMinutes == delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
            if ( distance > distanceKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
              keepThisObs = .false.
            end if
          end if

        end if
      end if

      ! save the observation if thinning criteria is satisfied
      if ( keepThisObs ) then
        headerIndexKeep(latBinIndex,lonBinIndex,timeBinIndex) = headerIndex
        distanceKeep(latBinIndex,lonBinIndex,timeBinIndex)    = distance
        delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex)  = delMinutes
        numChannelsKeep(latBinIndex,lonBinIndex,timeBinIndex) = numChannels
      end if

    end do HEADER

    ! return if no observations for this instrument
    if (countHeader == 0) then
      call deallocLocals()
      write(*,*) 'thn_thinByLatLonBoxes: no observations for this instrument'
      return
    end if

    ! communicate results to all other mpi tasks
    nsize = numLat * numLon * tim_nstepobs
    call rpn_comm_allgather(distanceKeep, nsize, 'mpi_real4',  &
                            allDistance,  nsize, 'mpi_real4', 'grid', ierr)
    call rpn_comm_allgather(delMinutesKeep, nsize, 'mpi_integer',  &
                            allDelMinutes,  nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(numChannelsKeep, nsize, 'mpi_integer',  &
                            allNumChannels,  nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(headerIndexKeep, nsize, 'mpi_integer',  &
                            allheaderIndex,  nsize, 'mpi_integer', 'grid', ierr)

    ! reset arrays that store info about kept obs
    headerIndexKeep(:,:,:) = -1
    numChannelsKeep(:,:,:) = 0
    distanceKeep(:,:,:)    = 0.0
    delMinutesKeep(:,:,:)  = deltmax

    do timeBinIndex = 1, tim_nstepobs
      do lonBinIndex = 1, numLon
        do latBinIndex = 1, numLat

          ! Apply thinning criteria to results from all mpi tasks
          do procIndex = 1, mpi_nprocs

            headerIndex = allHeaderIndex(latBinIndex,lonBinIndex,timeBinIndex,procIndex)
            distance    = allDistance(latBinIndex,lonBinIndex,timeBinIndex,procIndex)
            delMinutes  = allDelMinutes(latBinIndex,lonBinIndex,timeBinIndex,procIndex)
            numChannels = allNumChannels(latBinIndex,lonBinIndex,timeBinIndex,procIndex)
            
            keepThisObs = .false.

            ! keep if distance to box center smaller than limit and 
            ! time under limit fixed at input and maximise number of channels 
            if ( (distance < deltrad)                                              .and. &
                 (numChannels >= numChannelsKeep(latBinIndex,lonBinIndex,timeBinIndex)) .and. &
                 (delMinutes <= deltmax) ) keepThisObs = .true.
     
            ! keep the closest to bin central time
            if ( numChannels == numChannelsKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
              if ( headerIndexKeep(latBinIndex,lonBinIndex,timeBinIndex) == -1 ) then
                keepThisObs = .false.
              else
                if ( delMinutes > delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
                  keepThisObs = .false.
                end if
         
                if ( delMinutes == delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
                  if ( distance > distanceKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
                    keepThisObs = .false.
                  end if
                end if

              end if
            end if

            ! save the observation if thinning criteria is satisfied
            if ( keepThisObs ) then
              procIndexKeep(latBinIndex,lonBinIndex,timeBinIndex)   = procIndex
              headerIndexKeep(latBinIndex,lonBinIndex,timeBinIndex) = headerIndex
              distanceKeep(latBinIndex,lonBinIndex,timeBinIndex)    = distance
              delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex)  = delMinutes
              numChannelsKeep(latBinIndex,lonBinIndex,timeBinIndex) = numChannels
            end if

          end do ! procIndex

        end do
      end do
    end do

    ! determine which headerIndex values are rejected
    allocate(rejectThisHeader(obs_numheader(obsdat)))
    rejectThisHeader(:) = .true.
    do timeBinIndex = 1, tim_nstepobs
      do lonBinIndex = 1, numLon
        do latBinIndex = 1, numLat
          if (procIndexKeep(latBinIndex,lonBinIndex,timeBinIndex) == mpi_myid+1) then
            headerIndex = headerIndexKeep(latBinIndex,lonBinIndex,timeBinIndex)
            rejectThisHeader(headerIndex) = .false.
          end if
        end do
      end do
    end do

    ! modify flags, by setting bit 11 for those rejected
    call obs_set_current_header_list(obsdat,trim(familyType))
    HEADER2: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER2

      if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) /= codtyp) cycle HEADER2

      if (.not. rejectThisHeader(headerIndex)) cycle HEADER2

      ! loop over all body indices for this headerIndex
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY2: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY2

        obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
        call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,11))
      end do BODY2

    end do HEADER2
    deallocate(rejectThisHeader)

    call deallocLocals()

    write(*,*) 'thn_thinByLatLonBoxes: Finished'

  contains

    subroutine deallocLocals()
      implicit none

      deallocate(headerIndexKeep)
      deallocate(numChannelsKeep)
      deallocate(distanceKeep)
      deallocate(delMinutesKeep)
      deallocate(latdeg)
      deallocate(ngrd)
      deallocate(allHeaderIndex)
      deallocate(allNumChannels)
      deallocate(allDistance)
      deallocate(allDelMinutes)
      deallocate(procIndexKeep)
      
    end subroutine deallocLocals

  end subroutine thn_thinByLatLonBoxes

end module thinning_mod

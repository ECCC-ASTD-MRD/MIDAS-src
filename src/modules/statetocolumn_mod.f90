!-------------------------------------- LICENCE BEGIN ------------------------------------
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

module stateToColumn_mod
  ! MODULE stateToColumn (prefix='s2c' category='4. Data Object transformations')
  !
  ! :Purpose: Non-linear, tangent-linear and adjoint versions of
  !           horizontal-temporal interpolation between a gridStateVector object
  !           and a columnData object.
  !
  use mathPhysConstants_mod
  use earthConstants_mod
  use mpi, only : mpi_status_size ! this is the mpi library module
  use midasMpi_mod
  use codePrecision_mod
  use gridstatevector_mod
  use obsSpaceData_mod
  use columnData_mod
  use horizontalCoord_mod
  use obsTimeInterp_mod
  use windRotation_mod
  use utilities_mod
  use gridVariableTransforms_mod
  use varNameList_mod
  use slantprofilelatlon_mod
  use tovs_nl_mod
  use codtyp_mod
  use getGridPosition_mod
  use kdtree2_mod
  use calcHeightAndPressure_mod
  use humidityLimits_mod

  implicit none
  save
  private
  
  ! public routines
  public :: s2c_tl, s2c_ad, s2c_nl
  public :: s2c_bgcheck_bilin, s2c_getFootprintRadius, s2c_getWeightsAndGridPointIndexes
  public :: s2c_deallocInterpInfo

  ! private module variables and derived types

  type struct_stepProcData
    ! lat-lon location of observations to be interpolated
    real(8), pointer          :: allLat(:,:) => null()         ! (headerUsed, kIndex)
    real(8), pointer          :: allLon(:,:) => null()         ! (headerUsed, kIndex)
    ! lat-lon location on rotated grid of observations to be interpolated
    real(8), pointer          :: allLatRot(:,:,:) => null()    ! (subGrid, headerUsed, kIndex)
    real(8), pointer          :: allLonRot(:,:,:) => null()    ! (subGrid, headerUsed, kIndex)
    ! actual headerIndex, since the headerUsed is only for those obs with a non-zero interp weight
    integer, pointer          :: allHeaderIndex(:) => null()   ! (headerUsed)
    ! depotIndexBeg/End contain first/last indices into depots of interpolation weights and lat/lon indices
    integer, pointer          :: depotIndexBeg(:,:,:) => null() ! (subGrid, headerUsed, kIndex)
    integer, pointer          :: depotIndexEnd(:,:,:) => null() ! (subGrid, headerUsed, kIndex)
  end type struct_stepProcData

  type struct_interpInfo
    logical                   :: initialized = .false.
    type(struct_hco), pointer :: hco => null() ! horizontal grid object
    type(struct_uvr), pointer :: uvr => null() ! windRotation object
    type(struct_oti), pointer :: oti => null() ! obsTimeInterp object

    ! number of obs headers on each proc having a non-zero interp weight for each stepIndex (headerUsed)
    integer, pointer          :: allNumHeaderUsed(:,:) => null()    ! (step, proc)

    ! structure containing all interpolation information that depends on (proc, step)
    type(struct_stepProcData), allocatable :: stepProcData(:,:) ! (proc, step)

    ! interpolation weights and lat/lon indices are accessed via the 'stepProcData%depotIndexBeg/End'
    real(8), allocatable      :: interpWeightDepot(:)                ! (depotIndex)
    integer, pointer          :: latIndexDepot(:)                    ! (depotIndex)
    integer, pointer          :: lonIndexDepot(:)                    ! (depotIndex)
    character(len=2)          :: inputStateVectorType
  end type struct_interpInfo

  type(struct_interpInfo), target :: interpInfo_tlad, interpInfo_nl
  type(kdtree2), pointer  :: tree_nl => null()
  type(kdtree2), pointer  :: tree_tlad => null()

  character(len=20), parameter :: timeInterpType_tlad = 'LINEAR' ! hardcoded type of time interpolation for increment

  integer, parameter :: maxNumWrites = 50
  logical, parameter :: verbose = .false.

  ! "special" values of the footprint radius
  real(4), parameter :: nearestNeighbourFootprint = -3.0
  real(4), parameter ::             lakeFootprint = -2.0
  real(4), parameter ::         bilinearFootprint = -1.0
  integer, parameter :: maxNumLocalGridptsSearch = 1000
  integer, parameter :: minNumLocalGridptsSearch = 8

  ! namelist variables:
  logical :: slantPath_TO_nl
  logical :: slantPath_TO_tlad
  logical :: slantPath_RO_nl
  logical :: slantPath_RA_nl
  logical :: calcHeightPressIncrOnColumn
  logical :: useFootprintForTovs
  logical :: rejectObsNonMonotonicPressure

  integer, external :: get_max_rss

contains 


  !---------------------------------------------------------
  ! pressureProfileMonotonicityCheck
  !---------------------------------------------------------
  subroutine pressureProfileMonotonicityCheck(obsSpaceData, column)
    !
    ! :Purpose: Check for non monotonic pressure profiles that can be computed in slantpathmode
    !
    implicit none

    ! arguments
    type(struct_obs), intent(inout)       :: obsSpaceData
    type(struct_columnData), intent(inout):: column

    ! locals
    integer, parameter :: numWriteMax = 10
    integer :: headerIndex, bodyIndex, iterationCount, singularIndex, levelIndex
    integer :: pressureVarIndex
    integer :: nlv
    integer :: numWrites
    real(8), pointer :: pressureProfile(:)
    logical :: monotonicProfile
    integer, parameter :: nPressureVar =2
    character(len=4), parameter :: pressureVarList(nPressureVar)=['P_T ', 'P_M ']

    write(*,*) ' '
    write(*,*) 'pressureProfileMonotonicityCheck: START'
    write(*,*) ' '
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    numWrites = 0

    call obs_set_current_header_list(obsSpaceData, 'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER
      do pressureVarIndex = 1, nPressureVar
        pressureProfile => col_getColumn(column, headerIndex, pressureVarList(pressureVarIndex))
        nlv = size(pressureProfile)
        monotonicProfile = .true.
        iterationCount = 0
        iterationLoop:do
          singularIndex = -1
          levelSearch:do levelIndex = 1, nlv - 1
            if ( pressureProfile(levelIndex) > pressureProfile(levelIndex+1)) then
              singularIndex = levelIndex
              exit levelSearch
            end if
          end do levelSearch
          if ( singularIndex == -1 ) exit iterationLoop !regular profile or correction OK
          iterationCount = iterationCount + 1
          if (iterationCount == 1) then
            monotonicProfile = .false.
            if (numWrites < numWriteMax) then
              numWrites = numWrites + 1
              write(*,*) 'pressureProfileMonotonicityCheck: found non monotonic pressure profile:', &
                   pressureVarList(pressureVarIndex), pressureProfile
            end if
          end if
          if (singularIndex == 1) then !should never happen
            write(*,*) 'pressureProfileMonotonicityCheck: ', pressureProfile(1:2)
            call utl_abort('pressureProfileMonotonicityCheck: profile in the wrong order ?' &
                 // pressureVarList(pressureVarIndex))
          end if
          pressureProfile(singularIndex) = 0.5d0 * ( pressureProfile(singularIndex - 1) + pressureProfile(singularIndex + 1) )
          write(*,*) 'pressureProfileMonotonicityCheck: profile iteration', &
               pressureVarList(pressureVarIndex), iterationCount
        end do iterationLoop

        ! if requested reject the corrected profile
        if (.not. monotonicProfile .and. rejectObsNonMonotonicPressure) then
          call obs_headSet_i(obsSpaceData,OBS_ST1,headerIndex, &
               ibset( obs_headElem_i(obsSpaceData,OBS_ST1,headerIndex), 05))
          call obs_set_current_body_list(obsSpaceData, headerIndex)
          BODY: do
            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if (bodyIndex < 0) exit BODY
            if (rejectObsNonMonotonicPressure) then
              call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
              call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex, &
                   ibset(obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex),9))
            end if
          end do BODY
        end if
      end do ! loop on pressure variables
     
    end do HEADER

    write(*,*) 'pressureProfileMonotonicityCheck: END'

  end subroutine pressureProfileMonotonicityCheck

  !---------------------------------------------------------
  ! latlonChecksAnlGrid
  !---------------------------------------------------------
  subroutine latlonChecksAnlGrid(obsSpaceData, hco_core, moveObsAtPole)
    !
    ! :Purpose: Check the lat/lon of observations and modify if necessary
    !
    implicit none

    ! arguments
    type(struct_obs)          :: obsSpaceData
    type(struct_hco), pointer :: hco_core
    logical                   :: moveObsAtPole

    ! locals
    integer :: headerIndex, ierr
    integer :: idata, idatend, jdata, subGridIndex
    real(4) :: lat_r4, lon_r4, lat_deg_r4, lon_deg_r4
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    real(4) :: xposLowerBoundAnl_r4, xposUpperBoundAnl_r4
    real(8) :: lat_r8, lon_r8
    integer, save :: numWrites = 0

    ! external functions
    integer :: gdllfxy

    write(*,*) ' '
    write(*,*) 'latlonChecksAnlGrid: STARTING'
    write(*,*) ' '
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !-    Get the Analysis Grid structure
    !
    if ( hco_core % global ) then
       xposLowerBoundAnl_r4 = - huge(1.0) ! no limit since grid is global (periodic)
       xposUpperBoundAnl_r4 = + huge(1.0) ! no limit since grid is global (periodic)
    else
       xposLowerBoundAnl_r4 = 1.0
       xposUpperBoundAnl_r4 = real(hco_core % ni)
    end if

    header_loop: do headerIndex=1, obs_numheader(obsSpaceData)

      !- Get LatLon of observation location
      lat_r8 = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)
      lon_r8 = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex)
      lat_r4 = real(lat_r8,4)
      lon_r4 = real(lon_r8,4)
      if (lon_r4.lt.0.0         ) lon_r4 = lon_r4 + 2.0*MPC_PI_R4
      if (lon_r4.ge.2.*MPC_PI_R4) lon_r4 = lon_r4 - 2.0*MPC_PI_R4

      lat_deg_r4 = lat_r4 * MPC_DEGREES_PER_RADIAN_R4 ! Radian To Degree
      lon_deg_r4 = lon_r4 * MPC_DEGREES_PER_RADIAN_R4

      !
      !- Find the position in the analysis grid
      !
      ierr = gpos_getPositionXY( hco_core % EZscintID,  &
                                xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                                lat_deg_r4, lon_deg_r4, subGridIndex )

      !- Test if the obs is outside the analysis grid
      if ( xpos_r4 < xposLowerBoundAnl_r4  .or. &
           xpos_r4 > xposUpperBoundAnl_r4  .or. &
           ypos_r4 < 1.0                   .or. &
           ypos_r4 > real(hco_core % nj) ) then

        if ( hco_core % global ) then

          if ( moveObsAtPole ) then
            ! Modify latitude if we have an observation at or near the poles
            write(*,*) ''
            write(*,*) 'latlonChecksAnlGrid: Moving OBS inside the GLOBAL ANALYSIS grid, ', headerIndex
            write(*,*) '  true position : ', lat_deg_r4, lon_deg_r4, ypos_r4, xpos_r4

            !- Move the observation to the nearest grid point
            if ( ypos_r4 < 1.0 )                ypos_r4 = 1.0
            if ( ypos_r4 > real(hco_core % nj) ) ypos_r4 = real(hco_core % nj)

            ierr = gdllfxy( hco_core % EZscintID, &    ! IN
                            lat_deg_r4, lon_deg_r4, & ! OUT
                            xpos_r4, ypos_r4, 1)      ! IN

            write(*,*) '  new  position : ', lat_deg_r4, lon_deg_r4, ypos_r4, xpos_r4

            lat_r8 = real(lat_deg_r4,8) * MPC_RADIANS_PER_DEGREE_R8
            lon_r8 = real(lon_deg_r4,8) * MPC_RADIANS_PER_DEGREE_R8
            call obs_headSet_r(obsSpaceData,OBS_LAT,headerIndex, lat_r8) ! IN
            call obs_headSet_r(obsSpaceData,OBS_LON,headerIndex, lon_r8) ! IN
          else
            write(*,*)
            write(*,*) 'latlonChecksAnlGrid: OBS outside the GLOBAL ANALYSIS grid, but NOT moved, ', headerIndex
          end if

        else
          ! The observation is outside the domain
          ! In LAM Analysis mode we must discard this observation
          numWrites = numWrites + 1
          if (numWrites < maxNumWrites) then
            write(*,*) 'latlonChecksAnlGrid: Rejecting OBS outside the LAM ANALYSIS grid domain, ', headerIndex
            write(*,*) '  position : ', lat_deg_r4, lon_deg_r4, ypos_r4, xpos_r4
          else if (numWrites == maxNumWrites) then
            write(*,*) 'latlonChecksAnlGrid: More rejects, but reached maximum number of writes to the listing.'
          end if

          idata   = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
          idatend = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + idata -1
          do jdata = idata, idatend
            call obs_bodySet_i(obsSpaceData,OBS_ASS,JDATA, obs_notAssimilated)
          end do
          call obs_headSet_i(obsSpaceData,OBS_ST1,headerIndex,  &
                             ibset( obs_headElem_i(obsSpaceData,OBS_ST1,headerIndex), 05))
        end if

      end if

    end do header_loop

    write(*,*) 'latlonChecksAnlGrid: END'

  end subroutine latlonChecksAnlGrid

  !---------------------------------------------------------
  ! s2c_setupInterpInfo
  !---------------------------------------------------------
  subroutine s2c_setupInterpInfo( interpInfo, obsSpaceData, stateVector,  &
                                  headerIndexBeg, headerIndexEnd, &
                                  timeInterpType, rejectOutsideObs, &
                                  inputStateVectorType, lastCall_opt )
    ! :Purpose: Setup all of the information needed to quickly
    !           perform the horizontal interpolation to the observation
    !           locations.
    !
    implicit none

    ! arguments
    type(struct_interpInfo)    :: interpInfo
    type(struct_obs)           :: obsSpaceData
    type(struct_gsv), target   :: stateVector
    integer                    :: headerIndexBeg
    integer                    :: headerIndexEnd
    logical                    :: rejectOutsideObs
    character(len=*)           :: timeInterpType
    character(len=*)           :: inputStateVectorType
    logical, optional          :: lastCall_opt

    ! locals
    type(struct_gsv)          :: stateVector_VarsLevs_1Step, stateVector_Tiles_allVar_1Step
    type(struct_gsv)          :: stateVector_Tiles_1Step
    type(struct_gsv), save    :: stateVector_1Step
    type(struct_gsv), pointer :: stateVector_Tiles_ptr
    integer :: numHeader, numHeaderUsedMax, headerIndex, headerUsedIndex
    integer :: kIndex, kIndexCount, myKBeg
    integer :: numStep, stepIndex, fnom, fclos, nulnam, ierr
    integer :: procIndex, niP1, numGridptTotal, numHeaderUsed
    integer :: subGridIndex, subGridForInterp, numSubGridsForInterp
    real(8) :: latRot, lonRot, lat, lon
    real(4) :: lon_r4, lat_r4, lon_deg_r4, lat_deg_r4
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    real(4) :: footprintRadius_r4 ! (metres)
    integer, allocatable :: numGridpt(:), allNumHeaderUsed(:,:)
    integer, allocatable :: allHeaderIndex(:,:,:), headerIndexVec(:,:)
    real(8), allocatable :: lat_send_r8(:,:), lat_recv_r8(:,:), lon_send_r8(:,:), lon_recv_r8(:,:)
    real(4), allocatable :: footprintRadiusVec_r4(:), allFootprintRadius_r4(:,:,:)
    real(8), allocatable :: allLatOneLev(:,:)
    real(8), allocatable :: allLonOneLev(:,:)
    character(len=4), pointer :: varNames(:)
    character(len=4)          :: varLevel
    real(8), allocatable :: latColumn(:,:), lonColumn(:,:)
    real(8), allocatable :: latLev_T(:), lonLev_T(:), latLev_M(:), lonLev_M(:)
    real(8) :: latLev_S, lonLev_S
    real(4), pointer :: height3D_r4_ptr1(:,:,:), height3D_r4_ptr2(:,:,:)
    real(4), save, pointer :: height3D_T_r4(:,:,:), height3D_M_r4(:,:,:)
    real(8), pointer :: height3D_r8_ptr1(:,:,:)
    real(kdkind), allocatable :: positionArray(:,:)
    integer :: sendsizes(mmpi_nprocs), recvsizes(mmpi_nprocs), senddispls(mmpi_nprocs)
    integer :: recvdispls(mmpi_nprocs), allkBeg(mmpi_nprocs)
    integer :: codeType, nlev_T, nlev_M, levIndex 
    integer :: lonIndex, latIndex, gridIndex
    integer :: maxkcount, numkToSend, numTovsUsingFootprint, numAllTovs
    logical :: doSlantPath, SlantTO, SlantRO, SlantRA, firstHeaderSlantPathTO, firstHeaderSlantPathRO, firstHeaderSlantPathRA
    logical :: doSetup3dHeights, lastCall
    logical, save :: nmlAlreadyRead = .false.

    type(kdtree2), pointer  :: tree

    namelist /nams2c/ slantPath_TO_nl, slantPath_TO_tlad, slantPath_RO_nl, slantPath_RA_nl, calcHeightPressIncrOnColumn
    namelist /nams2c/ useFootprintForTovs, rejectObsNonMonotonicPressure 

    write(*,*) 's2c_setupInterpInfo: STARTING'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    write(*,*) 's2c_setupInterpInfo: inputStateVectorType=', inputStateVectorType

    if ( present(lastCall_opt) ) then
      lastCall = lastCall_opt
    else
      lastCall = .false.
    end if

    if ( .not. nmlAlreadyRead ) then
      nmlAlreadyRead = .true.

      ! default values
      slantPath_TO_nl   = .false.
      slantPath_TO_tlad = .false.
      slantPath_RO_nl   = .false.
      slantPath_RA_nl   = .false.
      calcHeightPressIncrOnColumn = .false.
      useFootprintForTovs = .false.
      rejectObsNonMonotonicPressure =.true.

      if ( .not. utl_isNamelistPresent('NAMS2C','./flnml') ) then
        if ( mmpi_myid == 0 ) then
          write(*,*) 's2c_setupInterpInfo: nams2c is missing in the namelist.'
          write(*,*) '                     The default values will be taken.'
        end if

      else
        ! reading namelist variables
        nulnam = 0
        ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
        read(nulnam, nml = nams2c, iostat = ierr)
        if ( ierr /= 0 ) call utl_abort('s2c_setupInterpInfo: Error reading namelist')
        ierr = fclos(nulnam)
      end if
      if ( mmpi_myid == 0 ) write(*, nml = nams2c)
    end if

    doSlantPath = .false.
    SlantTO     = .false.
    SlantRO     = .false.
    if ( slantPath_TO_nl   .and. inputStateVectorType == 'nl' ) then
      doSlantPath = .true.
      SlantTO     = .true.
    endif
    if ( slantPath_TO_tlad .and. inputStateVectorType /= 'nl' ) then
      doSlantPath = .true.
      SlantTO     = .true.
    endif
    if ( slantPath_RO_nl   .and. inputStateVectorType == 'nl' ) then
      doSlantPath = .true.
      SlantRO     = .true.
    endif
    if ( slantPath_RA_nl   .and. inputStateVectorType == 'nl' ) then
      doSlantPath = .true.
      SlantRA     = .true.
    endif
    write(*,*) 's2c_setupInterpInfo: doSlantPath, SlantTO, SlantRO, SlantRA = ', &
               doSlantPath, SlantTO, SlantRO, SlantRA

    numStep = stateVector%numStep
    numHeader = headerIndexEnd - headerIndexBeg + 1

    call oti_setup(interpInfo%oti, obsSpaceData, numStep,  &
                   headerIndexBeg, headerIndexEnd, &
                   interpType_opt=timeInterpType, flagObsOutside_opt=.true.)

    if ((stateVector%heightSfcPresent) .and. ( mmpi_myid == 0)) then
      mykBeg = 0 
    else
      mykBeg = stateVector%mykBeg
    end if   

    call rpn_comm_allgather(mykBeg, 1,'mpi_integer',       &
                            allkBeg,1,'mpi_integer','grid',ierr)

    ! Allow for periodicity in Longitude for global Gaussian grid
    if ( stateVector%hco%grtyp == 'G' .or. &
         (stateVector%hco%grtyp == 'Z' .and. stateVector%hco%global) ) then
      niP1 = stateVector%ni + 1
    else
      niP1 = stateVector%ni
    end if

    ! First count the number of headers for each stepIndex
    allocate(allNumHeaderUsed(numStep,mmpi_nprocs))
    allNumHeaderUsed(:,:) = 0
    do stepIndex = 1, numStep
      numHeaderUsed = 0

      header_loop1: do headerIndex = headerIndexBeg, headerIndexEnd

        ! if obs inside window, but zero weight for current stepIndex then skip it
        if ( oti_getTimeInterpWeight(interpInfo%oti,headerIndex,stepIndex) == 0.0d0 ) then
          cycle header_loop1
        end if

        numHeaderUsed = numHeaderUsed + 1

      end do header_loop1
      ! gather the number of obs over all processors for each timestep
      call rpn_comm_allgather(numHeaderUsed,                 1, 'MPI_INTEGER', &
                              allNumHeaderUsed(stepIndex,:), 1, 'MPI_INTEGER', &
                              'GRID',ierr)

    end do

    numHeaderUsedMax = maxval(allNumHeaderUsed(:,:))
    write(*,*) 's2c_setupInterpInfo: numHeaderUsedMax = ', numHeaderUsedMax

    ! temporary arrays
    allocate(headerIndexVec(numHeaderUsedMax,numStep))
    allocate(footprintRadiusVec_r4(numHeaderUsedMax))
    headerIndexVec(:,:) = 0

    ! copy the horizontal grid object
    interpInfo%hco => stateVector%hco

    ! setup the information for wind rotation
    if ( (gsv_varExist(varName='UU') .or. gsv_varExist(varName='VV')) .and.  &
         stateVector%hco%rotated ) then
      call uvr_Setup( interpInfo%uvr, & ! INOUT
                      stateVector%hco ) ! IN 
    end if

    allocate(interpInfo%stepProcData(mmpi_nprocs,numStep))
    do stepIndex = 1,numStep
      do procIndex = 1, mmpi_nprocs
        allocate(interpInfo%stepProcData(procIndex,stepIndex)%allLat(allNumHeaderUsed(stepIndex,procIndex),mykBeg:stateVector%mykEnd))
        allocate(interpInfo%stepProcData(procIndex,stepIndex)%allLon(allNumHeaderUsed(stepIndex,procIndex),mykBeg:stateVector%mykEnd))
        interpInfo%stepProcData(procIndex,stepIndex)%allLat(:,:) = 0.0d0
        interpInfo%stepProcData(procIndex,stepIndex)%allLon(:,:) = 0.0d0

        allocate(interpInfo%stepProcData(procIndex,stepIndex)%allHeaderIndex(allNumHeaderUsed(stepIndex,procIndex)))
        interpInfo%stepProcData(procIndex,stepIndex)%allHeaderIndex(:) = 0

        allocate(interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(interpInfo%hco%numSubGrid,numHeaderUsedMax,mykBeg:stateVector%mykEnd))
        allocate(interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(interpInfo%hco%numSubGrid,numHeaderUsedMax,mykBeg:stateVector%mykEnd))
        interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(:,:,:) = 0
        interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(:,:,:) = -1
      end do
    end do

    ! allocate arrays that will be returned
    allocate(interpInfo%allNumHeaderUsed(numStep,mmpi_nprocs))
    allocate(allLatOneLev(numHeaderUsedMax,mmpi_nprocs))
    allocate(allLonOneLev(numHeaderUsedMax,mmpi_nprocs))
    allocate(allFootprintRadius_r4(numHeaderUsedMax,numStep,mmpi_nprocs))
    allocate(numGridpt(interpInfo%hco%numSubGrid))
    allFootprintRadius_r4(:,:,:) = bilinearFootprint
    interpInfo%allNumHeaderUsed(:,:) = allNumHeaderUsed(:,:)

    if ( interpInfo%hco%rotated ) then
      do stepIndex = 1, numStep
        do procIndex = 1, mmpi_nprocs
          allocate(interpInfo%stepProcData(procIndex,stepIndex)%allLatRot(interpInfo%hco%numSubGrid,allNumHeaderUsed(stepIndex,procIndex),mykBeg:stateVector%mykEnd))
          allocate(interpInfo%stepProcData(procIndex,stepIndex)%allLonRot(interpInfo%hco%numSubGrid,allNumHeaderUsed(stepIndex,procIndex),mykBeg:stateVector%mykEnd))
          interpInfo%stepProcData(procIndex,stepIndex)%allLatRot(:,:,:) = 0.0d0
          interpInfo%stepProcData(procIndex,stepIndex)%allLonRot(:,:,:) = 0.0d0
        end do
      end do
    end if

    nlev_T = gsv_getNumLev(stateVector,'TH')
    nlev_M = gsv_getNumLev(stateVector,'MM')

    doSetup3dHeights = doSlantPath .and.  &
                       .not. gsv_isAllocated(stateVector_1Step) .and. &
                       gsv_varExist(stateVector,'Z_T') .and. &
                       gsv_varExist(stateVector,'Z_M')

    ! prepare for extracting the 3D height for slant-path calculation
    if ( doSetup3dHeights ) then

      write(*,*) 's2c_setupInterpInfo: extracting 3D heights for slant-path for ', inputStateVectorType 

      if ( inputStateVectorType == 'nl' ) then
        nullify(varNames)
        call gsv_varNamesList(varNames, stateVector)
        call gsv_allocate( stateVector_VarsLevs_1Step, 1, &
                           stateVector%hco, stateVector%vco, &
                           mpi_local_opt=.true., mpi_distribution_opt='VarsLevs', &
                           dataKind_opt=4, varNames_opt=varNames )

        call gsv_getField(stateVector,height3D_r4_ptr1)
        call gsv_getField(stateVector_VarsLevs_1Step,height3D_r4_ptr2)
        height3D_r4_ptr2(:,:,:) = height3D_r4_ptr1(:,:,:)

        call gsv_allocate( stateVector_Tiles_allVar_1Step, 1, &
                           stateVector%hco, stateVector%vco, &
                           mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                           dataKind_opt=4, varNames_opt=varNames )

        call gsv_transposeVarsLevsToTiles( stateVector_VarsLevs_1Step, stateVector_Tiles_allVar_1Step )
        call gsv_deallocate(statevector_VarsLevs_1Step)

        call gsv_allocate( stateVector_Tiles_1Step, 1, &
                           stateVector%hco, stateVector%vco, &
                           mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                           dataKind_opt=4, varNames_opt=(/'Z_M','Z_T'/) )

        call gsv_getField(stateVector_Tiles_allVar_1Step,height3D_r4_ptr1,'Z_T')
        call gsv_getField(stateVector_Tiles_1Step,height3D_r4_ptr2,'Z_T')
        height3D_r4_ptr2(:,:,:) = height3D_r4_ptr1(:,:,:)

        call gsv_getField(stateVector_Tiles_allVar_1Step,height3D_r4_ptr1,'Z_M')
        call gsv_getField(stateVector_Tiles_1Step,height3D_r4_ptr2,'Z_M')
        height3D_r4_ptr2(:,:,:) = height3D_r4_ptr1(:,:,:)

        call gsv_deallocate(stateVector_Tiles_allVar_1Step)

      else
        stateVector_Tiles_ptr => gvt_getStateVectorTrial('height')

        call gsv_allocate( stateVector_Tiles_1Step, 1, &
                           stateVector%hco, stateVector%vco, &
                           mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                           dataKind_opt=4, varNames_opt=(/'Z_M','Z_T'/) )

        call gsv_getField(stateVector_Tiles_ptr,height3D_r8_ptr1,'Z_T')
        call gsv_getField(stateVector_Tiles_1Step,height3D_r4_ptr2,'Z_T')
        height3D_r4_ptr2(:,:,:) = height3D_r8_ptr1(:,:,:)

        call gsv_getField(stateVector_Tiles_ptr,height3D_r8_ptr1,'Z_M')
        call gsv_getField(stateVector_Tiles_1Step,height3D_r4_ptr2,'Z_M')
        height3D_r4_ptr2(:,:,:) = height3D_r8_ptr1(:,:,:)

      end if ! inputStateVectorType 

      ! Communicate 3D height fields onto all mpi tasks
      call gsv_allocate( stateVector_1Step, 1, &
                         stateVector%hco, stateVector%vco, &
                         mpi_local_opt=.false., &
                         dataKind_opt=4, varNames_opt=(/'Z_M','Z_T'/) )
      call utl_tmg_start(32,'------s2c_Slant')
      call gsv_transposeTilesToMpiGlobal(stateVector_1Step, stateVector_Tiles_1Step)
      call utl_tmg_stop(32)
      call gsv_getField(stateVector_1Step,height3D_T_r4,'Z_T')
      call gsv_getField(stateVector_1Step,height3D_M_r4,'Z_M')
    
      write(*,*) 's2c_setupInterpInfo, height3D_T_r4='
      write(*,*) height3D_T_r4(1,1,:)
      write(*,*) 's2c_setupInterpInfo, height3D_M_r4='
      write(*,*) height3D_M_r4(1,1,:)

      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    end if ! doSlantPath 

    ! get observation lat-lon and footprint radius onto all mpi tasks
    step_loop2: do stepIndex = 1, numStep
      numHeaderUsed = 0

      footprintRadiusVec_r4(:) = bilinearFootprint

      header_loop2: do headerIndex = headerIndexBeg, headerIndexEnd

        ! if obs inside window, but zero weight for current stepIndex then skip it
        if ( oti_getTimeInterpWeight(interpInfo%oti, headerIndex, stepIndex) == 0.0d0 ) then
          cycle header_loop2
        end if

        numHeaderUsed = numHeaderUsed + 1
        headerIndexVec(numHeaderUsed,stepIndex) = headerIndex

        footprintRadiusVec_r4(numHeaderUsed) = s2c_getFootprintRadius(obsSpaceData, &
                                                                stateVector, headerIndex)

      end do header_loop2

      call rpn_comm_allgather(footprintRadiusVec_r4,                numHeaderUsedMax, 'MPI_REAL4', &
                              allFootprintRadius_r4(:,stepIndex,:), numHeaderUsedMax, 'MPI_REAL4', &
                              'GRID', ierr)

      allocate(latColumn(numHeaderUsedMax,allkBeg(1):stateVector%nk))
      allocate(lonColumn(numHeaderUsedMax,allkBeg(1):stateVector%nk))
      latColumn(:,:) = 0.0d0
      lonColumn(:,:) = 0.0d0

      if ( doSlantPath .and. &
           gsv_varExist(stateVector,'Z_T') .and. &
           gsv_varExist(stateVector,'Z_M') ) then

        allocate(latLev_T(nlev_T))
        allocate(lonLev_T(nlev_T))
        allocate(latLev_M(nlev_M))
        allocate(lonLev_M(nlev_M))
        latLev_T(:) = 0.0d0
        lonLev_T(:) = 0.0d0
        latLev_M(:) = 0.0d0
        lonLev_M(:) = 0.0d0

        firstHeaderSlantPathTO = .true.
        firstHeaderSlantPathRO = .true.
        firstHeaderSlantPathRA = .true.
        header_loop3: do headerUsedIndex = 1, numHeaderUsed
          headerIndex = headerIndexVec(headerUsedIndex,stepIndex)

          !- Get LatLon of observation location
          lat_r4 = real(obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex), 4)
          lon_r4 = real(obs_headElem_r(obsSpaceData, OBS_LON, headerIndex), 4)
          if (lon_r4 <  0.0          ) lon_r4 = lon_r4 + 2.0*MPC_PI_R4
          if (lon_r4 >= 2.0*MPC_PI_R4) lon_r4 = lon_r4 - 2.0*MPC_PI_R4

          codeType = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)

          if ( tvs_isIdBurpTovs(codeType) .and. SlantTO ) then
            if ( firstHeaderSlantPathTO ) then
              write(*,'(a,i3,a,i8)') 's2c_setupInterpInfo: start slant-path for TOVS. stepIndex = ', &
                   stepIndex,' and numHeaderUsed = ',numHeaderUsed
              firstHeaderSlantPathTO = .false.
            end if

            ! calculate lat/lon along the line of sight
            call utl_tmg_start(32,'------s2c_Slant')
            call slp_calcLatLonTovs( obsSpaceData, stateVector%hco, headerIndex, & ! IN
                                     height3D_T_r4, height3D_M_r4,               & ! IN
                                     latLev_T, lonLev_T,                         & ! OUT
                                     latLev_M, lonLev_M,                         & ! OUT
                                     latLev_S, lonLev_S             )              ! OUT
            call utl_tmg_stop(32)

          else if (codeType == codtyp_get_codtyp('ro') .and. SlantRO ) then
            if ( firstHeaderSlantPathRO ) then
              write(*,'(a,i3,a,i8)') 's2c_setupInterpInfo: start slant-path for RO. stepIndex = ', &
                   stepIndex,' and numHeaderUsed = ',numHeaderUsed
              firstHeaderSlantPathRO = .false.
            end if

            ! Calculate lat/lon along the GPSRO obs
            call utl_tmg_start(32,'------s2c_Slant')
            call slp_calcLatLonRO( obsSpaceData, stateVector%hco, headerIndex, & ! IN
                                   height3D_T_r4, height3D_M_r4,               & ! IN
                                   latLev_T, lonLev_T,                         & ! OUT
                                   latLev_M, lonLev_M,                         & ! OUT
                                   latLev_S, lonLev_S                          ) ! OUT
            call utl_tmg_stop(32)
          else if (codeType == codtyp_get_codtyp('radar') .and. SlantRA ) then
            if ( firstHeaderSlantPathRA ) then
              write(*,'(a,i3,a,i8)') 's2c_setupInterpInfo: start slant-path for RADAR. stepIndex=', &
                   stepIndex,' and numHeaderUsed=',numHeaderUsed
              firstHeaderSlantPathRA = .false.
            end if
            
             ! calculate lat/lon along the radar beam obs
             call slp_calcLatLonRadar( obsSpaceData, stateVector%hco, headerIndex, & ! IN
                                     height3D_T_r4, height3D_M_r4,                 & ! IN
                                     latLev_T, lonLev_T,                           & ! OUT
                                     latLev_M, lonLev_M,                           & ! OUT
                                     latLev_S, lonLev_S                           ) ! OUT
          else

            latLev_T(:) = real(lat_r4,8)
            lonLev_T(:) = real(lon_r4,8)
            latLev_M(:) = real(lat_r4,8)
            lonLev_M(:) = real(lon_r4,8)
            latLev_S = real(lat_r4,8)
            lonLev_S = real(lon_r4,8)
          end if

          ! check if the slanted lat/lon is inside the domain
          call latlonChecks ( obsSpaceData, stateVector%hco, & ! IN
                              headerIndex, rejectOutsideObs, & ! IN
                              latLev_T, lonLev_T,            & ! IN/OUT
                              latLev_M, lonLev_M,            & ! IN/OUT
                              latLev_S, lonLev_S )             ! IN/OUT

          ! put the lat/lon from TH/MM levels to kIndex
          do kIndex = allkBeg(1), stateVector%nk
            if ( kIndex == 0 ) then
              varLevel = 'SF'
            else
              levIndex = gsv_getLevFromK(stateVector,kIndex)
              varLevel = vnl_varLevelFromVarname(gsv_getVarNameFromK(stateVector,kIndex))
            end if

            if ( varLevel == 'TH' ) then
              latColumn(headerUsedIndex,kIndex) = latLev_T(levIndex)
              lonColumn(headerUsedIndex,kIndex) = lonLev_T(levIndex)
            else if ( varLevel == 'MM' ) then
              latColumn(headerUsedIndex,kIndex) = latLev_M(levIndex)
              lonColumn(headerUsedIndex,kIndex) = lonLev_M(levIndex)
            else if ( varLevel == 'SF' ) then
              latColumn(headerUsedIndex,kIndex) = latLev_S
              lonColumn(headerUsedIndex,kIndex) = lonLev_S
            else
              call utl_abort('s2c_setupInterpInfo: unknown value of varLevel')
            end if

          end do

        end do header_loop3

        ! MPI communication for the slant-path lat/lon

        maxkCount = maxval(stateVector%allkCount(:) + stateVector%allkBeg(:) - allkBeg(:))
        numkToSend = min(mmpi_nprocs,stateVector%nk)

        allocate(lat_recv_r8(numHeaderUsedMax,mmpi_nprocs))
        lat_recv_r8(:,:) = 0.0d0
        allocate(lat_send_r8(numHeaderUsedMax,mmpi_nprocs))
        lat_send_r8(:,:) = 0.0d0
        allocate(lon_recv_r8(numHeaderUsedMax,mmpi_nprocs))
        lon_recv_r8(:,:) = 0.0d0
        allocate(lon_send_r8(numHeaderUsedMax,mmpi_nprocs))
        lon_send_r8(:,:) = 0.0d0

        ! only send the data from tasks with data, same amount to all
        sendsizes(:) = 0
        do procIndex = 1, numkToSend
          sendsizes(procIndex) = numHeaderUsed
        end do
        senddispls(1) = 0
        do procIndex = 2, mmpi_nprocs
          senddispls(procIndex) = senddispls(procIndex-1) + numHeaderUsedMax
        end do

        recvdispls(1) = 0
        do procIndex = 2, mmpi_nprocs
          recvdispls(procIndex) = recvdispls(procIndex-1) + numHeaderUsedMax
        end do

        ! loop to send (at most) 1 level to (at most) all other mpi tasks
        do kIndexCount = 1, maxkCount

          sendsizes(:) = 0
          do procIndex = 1, mmpi_nprocs
            ! compute kIndex value being sent
            kIndex = kIndexCount + allkBeg(procIndex) - 1
            if ( kIndex <= stateVector%allkEnd(procIndex) ) then
              if( procIndex > numkToSend ) then
                write(*,*) 'procIndex, numkToSend = ', procIndex, numkToSend
                call utl_abort('ERROR: with numkToSend?')
              end if

              lat_send_r8(1:numHeaderUsed,procIndex) = latColumn(1:numHeaderUsed,kIndex)
              lon_send_r8(1:numHeaderUsed,procIndex) = lonColumn(1:numHeaderUsed,kIndex)
              sendsizes(procIndex) = numHeaderUsed
            else
              sendsizes(procIndex) = 0
            end if
          end do

          ! all tasks recv only from those with data
          kIndex = kIndexCount + mykBeg - 1
          if ( kIndex <= stateVector%mykEnd ) then
            do procIndex = 1, mmpi_nprocs
              recvsizes(procIndex) = allNumHeaderUsed(stepIndex,procIndex)
            end do
          else
            recvsizes(:) = 0
          end if

          call mpi_alltoallv(lat_send_r8, sendsizes, senddispls, mmpi_datyp_real8,  &
                             lat_recv_r8, recvsizes, recvdispls, mmpi_datyp_real8,  &
                             mmpi_comm_grid, ierr)
          call mpi_alltoallv(lon_send_r8, sendsizes, senddispls, mmpi_datyp_real8,  &
                             lon_recv_r8, recvsizes, recvdispls, mmpi_datyp_real8,  &
                             mmpi_comm_grid, ierr)

          do procIndex = 1, mmpi_nprocs
            ! all tasks copy the received step data into correct slot
            kIndex = kIndexCount + mykBeg - 1
            if ( kIndex <= stateVector%mykEnd ) then
              interpInfo%stepProcData(procIndex,stepIndex)%allLat(:,kIndex) = &
                   lat_recv_r8(1:allNumHeaderUsed(stepIndex,procIndex),procIndex)
              interpInfo%stepProcData(procIndex,stepIndex)%allLon(:,kIndex) = &
                   lon_recv_r8(1:allNumHeaderUsed(stepIndex,procIndex),procIndex)
            end if
          end do

        end do ! kIndexCount

        deallocate(lon_send_r8)
        deallocate(lon_recv_r8)
        deallocate(lat_send_r8)
        deallocate(lat_recv_r8)

        deallocate(latLev_T)
        deallocate(lonLev_T)
        deallocate(latLev_M)
        deallocate(lonLev_M)

      else ! not doSlantPath

        allocate(latLev_T(1))
        allocate(lonLev_T(1))
        allocate(latLev_M(1))
        allocate(lonLev_M(1))
        latLev_T(:) = 0.0d0
        lonLev_T(:) = 0.0d0
        latLev_M(:) = 0.0d0
        lonLev_M(:) = 0.0d0

        do headerUsedIndex = 1, numHeaderUsed
          headerIndex = headerIndexVec(headerUsedIndex,stepIndex)

          !- Get LatLon of observation location
          lat_r4 = real(obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex), 4)
          lon_r4 = real(obs_headElem_r(obsSpaceData, OBS_LON, headerIndex), 4)
          if (lon_r4 <  0.0          ) lon_r4 = lon_r4 + 2.0*MPC_PI_R4
          if (lon_r4 >= 2.0*MPC_PI_R4) lon_r4 = lon_r4 - 2.0*MPC_PI_R4

          latLev_T(:) = real(lat_r4,8)
          lonLev_T(:) = real(lon_r4,8)
          latLev_M(:) = real(lat_r4,8)
          lonLev_M(:) = real(lon_r4,8)

          ! check if the lat/lon is inside the domain
          call latlonChecks ( obsSpaceData, stateVector%hco, & ! IN
                              headerIndex, rejectOutsideObs, & ! IN
                              latLev_T, lonLev_T,            & ! IN/OUT
                              latLev_M, lonLev_M )             ! IN/OUT 

          latColumn(headerUsedIndex,allkBeg(1)) = latLev_T(1)
          lonColumn(headerUsedIndex,allkBeg(1)) = lonLev_T(1)
        end do

        ! gather geographical lat, lon positions of observations from all processors
        call rpn_comm_allgather(latColumn(:,allkBeg(1)), numHeaderUsedMax, 'MPI_REAL8', &
                                allLatOneLev(:,:), numHeaderUsedMax, 'MPI_REAL8', 'GRID', ierr)
        call rpn_comm_allgather(lonColumn(:,allkBeg(1)), numHeaderUsedMax, 'MPI_REAL8', &
                                allLonOneLev(:,:), numHeaderUsedMax, 'MPI_REAL8', 'GRID', ierr)

        k_loop: do kIndex = mykBeg, statevector%mykEnd
          do procIndex = 1, mmpi_nprocs
            interpInfo%stepProcData(procIndex,stepIndex)%allLat(:,kIndex) = allLatOneLev(1:allNumHeaderUsed(stepIndex,procIndex),procIndex)
            interpInfo%stepProcData(procIndex,stepIndex)%allLon(:,kIndex) = allLonOneLev(1:allNumHeaderUsed(stepIndex,procIndex),procIndex)
          end do
        end do k_loop

        deallocate(latLev_T)
        deallocate(lonLev_T)
        deallocate(latLev_M)
        deallocate(lonLev_M)

      end if ! doSlantPath 

      deallocate(lonColumn)
      deallocate(latColumn)

    end do step_loop2

    if ( gsv_isAllocated(stateVector_1Step) .and. lastCall ) then
      write(*,*) 's2c_setupInterpInfo: deallocate height3D fields'
      call gsv_deallocate(stateVector_1Step)
    end if
    deallocate(footprintRadiusVec_r4)

    write(*,*) 's2c_setupInterpInfo: latlonChecks and lat/lon MPI comm finished.'

    allocate(allHeaderIndex(numHeaderUsedMax,numStep,mmpi_nprocs))
    ! gather the headerIndexVec arrays onto all processors
    call rpn_comm_allgather(headerIndexVec, numHeaderUsedMax*numStep, 'MPI_INTEGER', &
                            allHeaderIndex, numHeaderUsedMax*numStep, 'MPI_INTEGER', &
                            'GRID',ierr)

    do procIndex = 1, mmpi_nprocs
      do stepIndex = 1, numStep
        do headerIndex = 1, allNumHeaderUsed(stepIndex,procIndex)
          interpInfo%stepProcData(procIndex, stepIndex)%allHeaderIndex(headerIndex) = allHeaderIndex(headerIndex,stepIndex,procIndex)
        end do
      end do
    end do

    ! create kdtree to use in footprint operator, if any footprint radius > 0.
    interpInfo%inputStateVectorType = inputStateVectorType
    if ( any(allFootprintRadius_r4(:,:,:) > 0.0) ) then
      if ( (inputStateVectorType == 'nl' .and. .not. associated(tree_nl))   .or. &
           (inputStateVectorType == 'tl' .and. .not. associated(tree_tlad)) .or. &
           (inputStateVectorType == 'ad' .and. .not. associated(tree_tlad)) ) then

        write(*,*) 's2c_setupInterpInfo: start creating kdtree for inputStateVectorType=', &
                   inputStateVectorType
        write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

        allocate(positionArray(3,statevector%hco%ni*statevector%hco%nj))

        gridIndex = 0
        do latIndex = 1, statevector%hco%nj
          do lonIndex = 1, statevector%hco%ni
            gridIndex = gridIndex + 1
            lat = real(stateVector % hco % lat2d_4(lonIndex,latIndex), 8)
            lon = real(stateVector % hco % lon2d_4(lonIndex,latIndex), 8)

            positionArray(:,gridIndex) = kdtree2_3dPosition(lon, lat)

          end do
        end do

        nullify(tree)
        tree => kdtree2_create(positionArray, sort=.false., rearrange=.true.) 

        if ( inputStateVectorType == 'nl' ) then
          tree_nl => tree
        else 
          tree_tlad => tree
        end if

        deallocate(positionArray)

        write(*,*) 's2c_setupInterpInfo: done creating kdtree for inputStateVectorType=', &
                   inputStateVectorType
        write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      end if
    end if

    do stepIndex = 1, numStep
      !$OMP PARALLEL DO PRIVATE (procIndex, kIndex, headerIndex, lat_deg_r4, lon_deg_r4, ierr, xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
      !$OMP subGridIndex, numSubGridsForInterp, subGridForInterp, lat, lon, latRot, lonRot, footprintRadius_r4, numGridpt)
      do procIndex = 1, mmpi_nprocs
        do kIndex = mykBeg, statevector%mykEnd
          do headerIndex = 1, allNumHeaderUsed(stepIndex,procIndex)

            ! Compute the rotated lat/lon, needed for the winds

            lat_deg_r4 = real(interpInfo%stepProcData(procIndex,stepIndex)%allLat(headerIndex,kIndex) *  &
                         MPC_DEGREES_PER_RADIAN_R8)
            lon_deg_r4 = real(interpInfo%stepProcData(procIndex,stepIndex)%allLon(headerIndex,kIndex) *  &
                         MPC_DEGREES_PER_RADIAN_R8)
            ierr = gpos_getPositionXY( stateVector%hco%EZscintID,   &
                                      xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                                      lat_deg_r4, lon_deg_r4, subGridIndex )

            if ( subGridIndex == 3 ) then
              ! both subGrids involved in interpolation, so first treat subGrid 1
              numSubGridsForInterp = 2
              subGridIndex = 1
            else
              ! only 1 subGrid involved in interpolation
              numSubGridsForInterp = 1
            end if

            do subGridForInterp = 1, numSubGridsForInterp

              if ( subGridForInterp == 1 ) then
                ! when only 1 subGrid involved, subGridIndex can be 1 or 2
              else
                ! when 2 subGrids, subGridIndex is set to 1 for 1st iteration, 2 for second
                subGridIndex = 2
              end if

              if ( interpInfo%hco%rotated .and.  &
                   (gsv_varExist(varName='UU') .or.  &
                    gsv_varExist(varName='VV')) ) then
                lat = interpInfo%stepProcData(procIndex,stepIndex)%allLat(headerIndex,kIndex)
                lon = interpInfo%stepProcData(procIndex,stepIndex)%allLon(headerIndex,kIndex)
                call uvr_RotateLatLon( interpInfo%uvr,   & ! INOUT
                                       subGridIndex,     & ! IN
                                       latRot, lonRot,   & ! OUT (radians)
                                       lat, lon,         & ! IN  (radians)
                                       'ToLatLonRot')      ! IN
                interpInfo%stepProcData(procIndex,stepIndex)%allLatRot(subGridIndex,headerIndex,kIndex) = latRot
                interpInfo%stepProcData(procIndex,stepIndex)%allLonRot(subGridIndex,headerIndex,kIndex) = lonRot
              end if

            end do ! subGridForInterp

            ! Count total number of grid points for allocating interp depot

            footprintRadius_r4 = allFootprintRadius_r4(headerIndex,stepIndex,procIndex)

            call s2c_setupHorizInterp(footprintRadius_r4, interpInfo, &
                                      stateVector, headerIndex, kIndex, stepIndex, &
                                      procIndex, numGridpt)

            ! for now, just store the number of gridpts for each obs in depotIndexEnd
            if ( (subGridIndex == 1) .or. (subGridIndex == 2) ) then
              ! indices for only 1 subgrid, other will have zeros
              interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex,headerIndex,kIndex) = numGridpt(subGridIndex)
            else
              ! locations on both subGrids will be averaged
              interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(1,headerIndex,kIndex) = numGridpt(1)
              interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(2,headerIndex,kIndex) = numGridpt(2)
            end if

          end do ! headerIndex
        end do ! kIndex
      end do ! procIndex
      !$OMP END PARALLEL DO
    end do ! stepIndex

    numGridptTotal = 0
    do stepIndex = 1, numStep
      do procIndex = 1, mmpi_nprocs
        do kIndex = mykBeg, statevector%mykEnd
          do headerIndex = 1, allNumHeaderUsed(stepIndex,procIndex)
            do subGridIndex = 1, interpInfo%hco%numSubGrid
              if ( interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex,headerIndex,kIndex) /= -1 ) then
                interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex,headerIndex,kIndex) = numGridptTotal + 1
                numGridptTotal = numGridptTotal + interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex,headerIndex,kIndex)
                interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex,headerIndex,kIndex) = numGridptTotal
              end if
            end do ! subGridIndex
          end do ! headerIndex
        end do ! kIndex
      end do ! procIndex
    end do ! stepIndex

    deallocate(allHeaderIndex)

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! now that we know the size, allocate main arrays for storing interpolation information
    write(*,*) 's2c_setupInterpInfo: numGridptTotal = ', numGridptTotal
    allocate( interpInfo%latIndexDepot(numGridptTotal) )
    allocate( interpInfo%lonIndexDepot(numGridptTotal) )
    allocate( interpInfo%interpWeightDepot(numGridptTotal) )

    call utl_tmg_start(33,'------s2c_SetupWeights')
    !$OMP PARALLEL DO PRIVATE (procIndex, stepIndex, kIndex, headerIndex, footprintRadius_r4, numGridpt)
    do procIndex = 1, mmpi_nprocs
      do stepIndex = 1, numStep
        do kIndex = mykBeg, statevector%mykEnd

          do headerIndex = 1, allNumHeaderUsed(stepIndex,procIndex)

            footprintRadius_r4 = allFootprintRadius_r4(headerIndex, stepIndex, procIndex)

            call s2c_setupHorizInterp(footprintRadius_r4, interpInfo, stateVector, &
                                      headerIndex, kIndex, stepIndex, procIndex, numGridpt)

          end do ! headerIndex

        end do ! kIndex
      end do ! stepIndex
    end do ! procIndex
    !$OMP END PARALLEL DO
    call utl_tmg_stop(33)

    ! reject obs in obsSpaceData if any processor has zero weight
    ! called when a mask exists to catch land contaminated ocean obs
    if ( stateVector%oceanMask%maskPresent ) then
      call s2c_rejectZeroWeightObs(interpInfo,obsSpaceData,mykBeg,stateVector%mykEnd)
    end if

    ! on the last call, deallocate the tree_nl/tree_tlad
    if ( lastCall ) then
      if ( inputStateVectorType == 'nl' .and. associated(tree_nl) ) then
        call kdtree2_destroy(tree_nl)
      else if ( (inputStateVectorType == 'tl' .or. inputStateVectorType == 'ad') .and. &
                associated(tree_tlad) ) then
        call kdtree2_destroy(tree_tlad)
      end if
    end if

    ! Count the number of TOVS using footprint operator on one level
    if ( useFootprintForTovs ) then 
      numTovsUsingFootprint = 0
      numAllTovs = 0
      procIndex = mmpi_myid + 1
      do stepIndex = 1, numStep
        do headerUsedIndex = 1, allNumHeaderUsed(stepIndex,procIndex)
          footprintRadius_r4 = allFootprintRadius_r4(headerUsedIndex, stepIndex, procIndex)
          headerIndex = headerIndexVec(headerUsedIndex,stepIndex)
          codeType = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)

          if ( tvs_isIdBurpTovs(codeType) ) then
            if ( footprintRadius_r4 > 0.0 ) numTovsUsingFootprint = numTovsUsingFootprint + 1
            numAllTovs = numAllTovs + 1 
          end if
        end do
      end do

      if ( numAllTovs > 0 ) then 
        write(*,'(A,2(I5,A2),F5.1,A)') 's2c_setupInterpInfo: numTovsUsingFootprint/numAllTovs=', &
                       numTovsUsingFootprint, ' /', numAllTovs, ' (', &
                       real(numTovsUsingFootprint) / real(numAllTovs) * 100.0, '%)'
      end if
    end if
    
    deallocate(allFootprintRadius_r4)
    deallocate(allLonOneLev)
    deallocate(allLatOneLev)

    deallocate(headerIndexVec)
    deallocate(allNumHeaderUsed)

    interpInfo%initialized = .true.

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 's2c_setupInterpInfo: FINISHED'

  end subroutine s2c_setupInterpInfo

  !---------------------------------------------------------
  ! s2c_tl
  !---------------------------------------------------------
  subroutine s2c_tl(statevector_in, columnAnlInc, columnTrlOnAnlIncLev, obsSpaceData)
    !
    ! :Purpose: Tangent linear version of the horizontal
    !           interpolation, used for the increment (or perturbations).
    !
    implicit none

    ! arguments
    type(struct_gsv), target   :: stateVector_in
    type(struct_obs)           :: obsSpaceData
    type(struct_columnData)    :: columnAnlInc
    type(struct_columnData)    :: columnTrlOnAnlIncLev

    ! locals
    type(struct_gsv)           :: stateVector_VarsLevs
    type(struct_gsv), pointer  :: stateVector
    integer :: kIndex, kIndex2, levIndex, kCount, stepIndex, numStep, mykEndExtended
    integer :: headerIndex, numHeader, numHeaderMax, yourNumHeader
    integer :: procIndex, nsize, ierr, headerUsedIndex
    real(8) :: weight
    real(8), pointer     :: allCols_ptr(:,:)
    real(pre_incrReal), pointer :: ptr4d(:,:,:,:)
    real(pre_incrReal), pointer :: ptr3d_UV(:,:,:)
    real(8), allocatable :: cols_hint(:,:,:)
    real(8), allocatable :: cols_send(:,:)
    real(8), allocatable :: cols_recv(:,:)
    real(8), allocatable :: cols_send_1proc(:)
    logical              :: rejectOutsideObs
    character(len=4)     :: varName
    character(len=4), pointer :: varNames(:)

    call utl_tmg_start(30,'--StateToColumn')

    if ( mmpi_myid == 0 ) write(*,*) 's2c_tl: Horizontal interpolation StateVector --> ColumnData'
    call utl_tmg_start(37,'----s2c_TL')

    call rpn_comm_barrier('GRID',ierr)

    if ( .not. gsv_isAllocated(stateVector_in) ) then 
      call utl_abort('s2c_tl: stateVector must be allocated')
    end if

    if (interpInfo_tlad%initialized) then
      if (.not. hco_equal(interpInfo_tlad%hco,stateVector_in%hco)) then
        write(*,*) 's2c_tl: WARNING! Current hco grid parameters differ from allocated interpInfo_tlad!'
        write(*,*) 's2c_tl: InterpInfo_tlad will be deallocated.'
	call s2c_deallocInterpInfo(inputStateVectorType='tlad')
      end if
    end if

    ! check the column and statevector have same nk/varNameList
    call checkColumnStatevectorMatch(columnAnlInc,statevector_in)

    ! if we only compute Height and Pressure on column, make copy without them
    if ( calcHeightPressIncrOnColumn ) then
      allocate(stateVector)
      call gsv_allocate( stateVector, statevector_in%numstep, &
                         statevector_in%hco, statevector_in%vco, &
                         mpi_local_opt=.true., &
                         dataKind_opt=gsv_getDataKind(statevector_in), &
                         allocHeight_opt=.false., allocPressure_opt=.false. )
      call gsv_copy(stateVector_in, stateVector, allowVarMismatch_opt=.true.)
    else
      stateVector => stateVector_in

      ! calculate delP_T/delP_M and del Z_T/Z_M on the grid
      call gvt_transform( statevector, 'ZandP_tl' )
    end if

    nullify(varNames)
    call gsv_varNamesList(varNames, statevector)
    call gsv_allocate( statevector_VarsLevs, statevector%numstep, &
                       statevector%hco, statevector%vco,          &
                       mpi_local_opt=.true., mpi_distribution_opt='VarsLevs', &
                       dataKind_opt=gsv_getDataKind(statevector), &
                       varNames_opt=varNames )
    deallocate(varNames)
    call gsv_transposeTilesToVarsLevs( statevector, statevector_VarsLevs )

    numStep = stateVector_VarsLevs%numStep
    numHeader = obs_numheader(obsSpaceData)
    call rpn_comm_allreduce(numHeader, numHeaderMax, 1,  &
                            'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)

    if ( .not. interpInfo_tlad%initialized ) then
      rejectOutsideObs = .false.
      call utl_tmg_stop(37)
      call utl_tmg_start(31,'----s2c_Setups')
      call s2c_setupInterpInfo( interpInfo_tlad, obsSpaceData, stateVector_VarsLevs,  &
                                1, numHeader, timeInterpType_tlad,  rejectOutsideObs, &
                                inputStateVectorType='tl' )
      call utl_tmg_stop(31)
      call utl_tmg_start(37,'----s2c_TL')
    end if

    ! arrays for interpolated column for 1 level/variable and each time step
    allocate(cols_hint(maxval(interpInfo_tlad%allNumHeaderUsed),numStep,mmpi_nprocs))
    cols_hint(:,:,:) = 0.0d0

    ! arrays for sending/receiving time interpolated column for 1 level/variable
    allocate(cols_send(numHeaderMax,mmpi_nprocs))
    cols_send(:,:) = 0.0d0

    allocate(cols_recv(numHeaderMax,mmpi_nprocs))
    cols_recv(:,:) = 0.0d0

    allocate(cols_send_1proc(numHeaderMax))

    ! set contents of column to zero
    allCols_ptr => col_getAllColumns(columnAnlInc)
    if ( numHeader > 0 ) allCols_ptr(:,:) = 0.0d0

    call gsv_getField(stateVector_VarsLevs, ptr4d)

    mykEndExtended = stateVector_VarsLevs%mykBeg + maxval(stateVector_VarsLevs%allkCount(:)) - 1

    kCount = 0
    k_loop: do kIndex = stateVector_VarsLevs%mykBeg, mykEndExtended

      kCount = kCount + 1

      if ( kIndex <= stateVector_VarsLevs%mykEnd ) then
        varName = gsv_getVarNameFromK(statevector,kIndex)

        if ( varName == 'UU' .or. varName == 'VV' ) then
          call gsv_getFieldUV(stateVector_VarsLevs,ptr3d_UV,kIndex)
        end if

        call utl_tmg_start(38,'------s2c_TL_Hinterp')
        !$OMP PARALLEL DO PRIVATE (stepIndex, procIndex, yourNumHeader, headerIndex)
        step_loop: do stepIndex = 1, numStep
          if ( maxval(interpInfo_tlad%allNumHeaderUsed(stepIndex,:)) == 0 ) cycle step_loop

          ! interpolate to the columns destined for all procs for all steps and one lev/var
          do procIndex = 1, mmpi_nprocs
            yourNumHeader = interpInfo_tlad%allNumHeaderUsed(stepIndex,procIndex)
            if ( yourNumHeader > 0 ) then
              if ( varName == 'UU' ) then
                call myezuvint_tl( cols_hint(1:yourNumHeader,stepIndex,procIndex), 'UU',  &
                                   ptr4d(:,:,kIndex,stepIndex), ptr3d_UV(:,:,stepIndex),  &
                                   interpInfo_tlad, kIndex, stepIndex, procIndex )
              else if ( varName == 'VV' ) then
                call myezuvint_tl( cols_hint(1:yourNumHeader,stepIndex,procIndex), 'VV',  &
                                   ptr3d_UV(:,:,stepIndex), ptr4d(:,:,kIndex,stepIndex),  &
                                   interpInfo_tlad, kIndex, stepIndex, procIndex )
              else
                call myezsint_tl( cols_hint(1:yourNumHeader,stepIndex,procIndex),  &
                                  ptr4d(:,:,kIndex,stepIndex), interpInfo_tlad, kIndex, &
                                  stepIndex, procIndex )
              end if
            end if
          end do

        end do step_loop
        !$OMP END PARALLEL DO
        call utl_tmg_stop(38)

        ! interpolate in time to the columns destined for all procs and one level/variable
        do procIndex = 1, mmpi_nprocs
          cols_send_1proc(:) = 0.0d0
          do stepIndex = 1, numStep
            !$OMP PARALLEL DO PRIVATE (headerUsedIndex, headerIndex, weight)
            do headerUsedIndex = 1, interpInfo_tlad%allNumHeaderUsed(stepIndex, procIndex)
              headerIndex = interpInfo_tlad%stepProcData(procIndex, stepIndex)%allHeaderIndex(headerUsedIndex)
              weight = oti_getTimeInterpWeightMpiGlobal(interpInfo_tlad%oti,  &
                                                        headerIndex,stepIndex,procIndex)
              cols_send_1proc(headerIndex) = cols_send_1proc(headerIndex) &
                            + weight * cols_hint(headerUsedIndex,stepIndex,procIndex)

            end do
            !$OMP END PARALLEL DO
          end do
          cols_send(:,procIndex) = cols_send_1proc(:)
        end do

      else

        ! this value of k does not exist on this mpi task
        cols_send(:,:) = 0.0

      end if ! if kIndex <= mykEnd

      call rpn_comm_barrier('GRID',ierr)

      ! mpi communication: alltoall for one level/variable
      nsize = numHeaderMax
      if(mmpi_nprocs > 1) then
        call rpn_comm_alltoall(cols_send, nsize, 'MPI_REAL8',  &
                               cols_recv, nsize, 'MPI_REAL8', 'GRID', ierr)
      else
        cols_recv(:,1) = cols_send(:,1)
      end if

      ! reorganize ensemble of distributed columns
      !$OMP PARALLEL DO PRIVATE (procIndex, kIndex2, varName, levIndex, allCols_ptr, headerIndex)
      proc_loop: do procIndex = 1, mmpi_nprocs
        ! This is kIndex value of source (can be different for destination)
        kIndex2 = statevector_VarsLevs%allkBeg(procIndex) + kCount - 1
        if ( kIndex2 > stateVector_VarsLevs%allkEnd(procIndex) ) cycle proc_loop

        ! Figure out which variable/level of destination
        varName = gsv_getVarNameFromK(statevector,kIndex2)
        levIndex = gsv_getLevFromK(statevector,kIndex2)
        allCols_ptr => col_getAllColumns(columnAnlInc,varName)

        do headerIndex = 1, numHeader
          allCols_ptr(levIndex,headerIndex) = cols_recv(headerIndex,procIndex)
        end do
      end do proc_loop
      !$OMP END PARALLEL DO

    end do k_loop

    if (calcHeightPressIncrOnColumn) then
      ! calculate delP_T/delP_M and  del Z_T/Z_M on the columns
      call czp_calcZandP_tl(columnAnlInc, columnTrlOnAnlIncLev)
    end if

    deallocate(cols_hint)
    deallocate(cols_send)
    deallocate(cols_recv)
    deallocate(cols_send_1proc)

    call gsv_deallocate( statevector_VarsLevs )
    if (calcHeightPressIncrOnColumn) then
      call gsv_deallocate( stateVector )
    end if
    
    if (slantPath_TO_tlad) call pressureProfileMonotonicityCheck(obsSpaceData, columnTrlOnAnlIncLev)

    call utl_tmg_stop(37)

    call utl_tmg_stop(30)

  end subroutine s2c_tl

  !---------------------------------------------------------
  ! s2c_ad
  !---------------------------------------------------------
  subroutine s2c_ad(statevector_out, columnAnlInc, columnTrlOnAnlIncLev, obsSpaceData)
    !
    ! :Purpose: Adjoint version of the horizontal interpolation,
    !           used for the cost function gradient with respect to the increment.
    !
    implicit none

    ! arguments
    type(struct_gsv), target   :: stateVector_out
    type(struct_obs)           :: obsSpaceData
    type(struct_columnData)    :: columnAnlInc
    type(struct_columnData)    :: columnTrlOnAnlIncLev

    ! locals
    type(struct_gsv)           :: stateVector_VarsLevs
    type(struct_gsv), pointer  :: stateVector
    integer :: kIndex, kIndex2, kCount, levIndex, stepIndex, numStep, mykEndExtended
    integer :: headerIndex, numHeader, numHeaderMax, yourNumHeader
    integer :: procIndex, nsize, ierr, headerUsedIndex
    character(len=4)     :: varName
    real(8) :: weight
    real(8), pointer     :: allCols_ptr(:,:)
    real(pre_incrReal), pointer :: ptr4d(:,:,:,:), ptr3d_UV(:,:,:)
    real(8), allocatable :: cols_hint(:,:,:)
    real(8), allocatable :: cols_send(:,:)
    real(8), allocatable :: cols_recv(:,:)
    logical              :: rejectOutsideObs
    character(len=4), pointer :: varNames(:)

    call utl_tmg_start(30,'--StateToColumn')

    if(mmpi_myid == 0) write(*,*) 's2c_ad: Adjoint of horizontal interpolation StateVector --> ColumnData'
    call utl_tmg_start(39,'----s2c_AD')

    call rpn_comm_barrier('GRID',ierr)

    if ( .not. gsv_isAllocated(stateVector_out) ) then 
      call utl_abort('s2c_ad: stateVector must be allocated')
    end if

    if (interpInfo_tlad%initialized) then
      if (.not. hco_equal(interpInfo_tlad%hco,stateVector_out%hco)) then
        write(*,*) 's2c_ad: WARNING! Current hco grid parameters differ from allocated interpInfo_tlad!'
        write(*,*) 's2c_ad: InterpInfo_tlad will be deallocated.'
        call s2c_deallocInterpInfo(inputStateVectorType='tlad')
      end if
    end if


    ! if we only compute Height and Pressure on column, make copy without them
    if (calcHeightPressIncrOnColumn) then
      allocate(stateVector)
      call gsv_allocate( stateVector, statevector_out%numstep, &
                         statevector_out%hco, statevector_out%vco, &
                         mpi_local_opt=.true., &
                         dataKind_opt=gsv_getDataKind(statevector_out), &
                         allocHeight_opt=.false., allocPressure_opt=.false. )
      ! Adjoint of calculate del Z_T/Z_M and delP_T/delP_M on the columns
      call czp_calcZandP_ad(columnAnlInc, columnTrlOnAnlIncLev)
    else
      stateVector => stateVector_out
    end if

    nullify(varNames)
    call gsv_varNamesList(varNames, statevector)
    call gsv_allocate( statevector_VarsLevs, statevector%numstep, &
                       statevector%hco, statevector%vco,          &
                       mpi_local_opt=.true., mpi_distribution_opt='VarsLevs', &
                       dataKind_opt=gsv_getDataKind(statevector), &
                       varNames_opt=varNames )
    deallocate(varNames)
    call gsv_zero( statevector_VarsLevs )

    numStep = stateVector_VarsLevs%numStep
    numHeader = obs_numheader(obsSpaceData)
    call rpn_comm_allreduce(numHeader, numHeaderMax, 1,  &
                            'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)

    if ( .not. interpInfo_tlad%initialized ) then
      rejectOutsideObs = .false.
      call utl_tmg_stop(39)
      call utl_tmg_start(31,'----s2c_Setups')
      call s2c_setupInterpInfo( interpInfo_tlad, obsSpaceData, stateVector_VarsLevs,  &
                                1, numHeader, timeInterpType_tlad, rejectOutsideObs,  &
                                inputStateVectorType='ad' )
      call utl_tmg_stop(31)
      call utl_tmg_start(39,'----s2c_AD')
    end if

    ! arrays for interpolated column for 1 level/variable and each time step
    allocate(cols_hint(maxval(interpInfo_tlad%allNumHeaderUsed),numStep,mmpi_nprocs))
    cols_hint(:,:,:) = 0.0d0

    ! arrays for sending/receiving time interpolated column for 1 level/variable
    allocate(cols_send(numHeaderMax,mmpi_nprocs))
    cols_send(:,:) = 0.0d0

    allocate(cols_recv(numHeaderMax,mmpi_nprocs))
    cols_recv(:,:) = 0.0d0

    ! set contents of column to zero
    allCols_ptr => col_getAllColumns(columnAnlInc)

    call gsv_getField(stateVector_VarsLevs,ptr4d)
    mykEndExtended = stateVector_VarsLevs%mykBeg + maxval(stateVector_VarsLevs%allkCount(:)) - 1

    kCount = 0
    k_loop: do kIndex = stateVector_VarsLevs%mykBeg, mykEndExtended

      kCount = kCount + 1

      ! reorganize ensemble of distributed columns
      !$OMP PARALLEL DO PRIVATE (procIndex, kIndex2, varName, levIndex, allCols_ptr, headerIndex)
      proc_loop: do procIndex = 1, mmpi_nprocs
        ! This is kIndex value of destination (can be different for source)
        kIndex2 = statevector_VarsLevs%allkBeg(procIndex) + kCount - 1
        if ( kIndex2 > stateVector_VarsLevs%allkEnd(procIndex) ) cycle proc_loop
        
        ! Figure out which variable/level of source
        varName = gsv_getVarNameFromK(statevector,kIndex2)
        levIndex = gsv_getLevFromK(statevector,kIndex2)
        allCols_ptr => col_getAllColumns(columnAnlInc,varName)

        do headerIndex = 1, numHeader
          cols_send(headerIndex,procIndex) = allCols_ptr(levIndex,headerIndex)
        end do
      end do proc_loop
      !$OMP END PARALLEL DO

      call rpn_comm_barrier('GRID',ierr)

      ! mpi communication: alltoall for one level/variable
      nsize = numHeaderMax
      if(mmpi_nprocs > 1) then
        call rpn_comm_alltoall(cols_send, nsize, 'MPI_REAL8',  &
                               cols_recv, nsize, 'MPI_REAL8', 'GRID', ierr)
      else
        cols_recv(:,1) = cols_send(:,1)
      end if

      if ( kIndex <= stateVector_VarsLevs%mykEnd ) then
        varName = gsv_getVarNameFromK(statevector,kIndex)

        if ( varName == 'UU' .or. varName == 'VV' ) then
          call gsv_getFieldUV(stateVector_VarsLevs, ptr3d_UV, kIndex)
        end if

        ! interpolate in time to the columns destined for all procs and one level/variable
        do procIndex = 1, mmpi_nprocs
          do stepIndex = 1, numStep
            !$OMP PARALLEL DO PRIVATE (headerIndex, headerUsedIndex, weight)
            do headerUsedIndex = 1, interpInfo_tlad%allNumHeaderUsed(stepIndex, procIndex)

              headerIndex = interpInfo_tlad%stepProcData(procIndex, stepIndex)%allHeaderIndex(headerUsedIndex)
              weight = oti_getTimeInterpWeightMpiGlobal(interpInfo_tlad%oti,  &
                                                        headerIndex,stepIndex,procIndex)

              cols_hint(headerUsedIndex,stepIndex,procIndex) =  &
                   weight * cols_recv(headerIndex,procIndex)

            end do
            !$OMP END PARALLEL DO
          end do
        end do

        call utl_tmg_start(40,'------s2c_AD_Hinterp')
        !$OMP PARALLEL DO PRIVATE (stepIndex, procIndex, yourNumHeader)
        step_loop: do stepIndex = 1, numStep
          if ( maxval(interpInfo_tlad%allNumHeaderUsed(stepIndex,:)) == 0 ) cycle step_loop

          ! interpolate to the columns destined for all procs for all steps and one lev/var
          do procIndex = 1, mmpi_nprocs
            yourNumHeader = interpInfo_tlad%allNumHeaderUsed(stepIndex,procIndex)
            if ( yourNumHeader > 0 ) then
              if ( varName == 'UU' ) then
                call myezuvint_ad( cols_hint(1:yourNumHeader,stepIndex,procIndex), 'UU',  &
                                   ptr4d(:,:,kIndex,stepIndex), ptr3d_UV(:,:,stepIndex),  &
                                   interpInfo_tlad, kIndex, stepIndex, procIndex )
              else if ( varName == 'VV' ) then
                call myezuvint_ad( cols_hint(1:yourNumHeader,stepIndex,procIndex), 'VV',  &
                                   ptr3d_UV(:,:,stepIndex), ptr4d(:,:,kIndex,stepIndex),  &
                                   interpInfo_tlad, kIndex, stepIndex, procIndex )
              else
                call myezsint_ad( cols_hint(1:yourNumHeader,stepIndex,procIndex), &
                                  ptr4d(:,:,kIndex,stepIndex), interpInfo_tlad, kIndex, &
                                  stepIndex, procIndex )
              end if
            end if
          end do

        end do step_loop
        !$OMP END PARALLEL DO
        call utl_tmg_stop(40)

      end if ! if kIndex <= mykEnd

    end do k_loop

    deallocate(cols_hint)
    deallocate(cols_send)
    deallocate(cols_recv)

    call rpn_comm_barrier('GRID',ierr)

    call gsv_transposeTilesToVarsLevsAd( statevector_VarsLevs, statevector )

    if (calcHeightPressIncrOnColumn) then
      call gsv_zero(statevector_out)
      call gsv_copy(stateVector, stateVector_out, allowVarMismatch_opt=.true.)
    else
      ! Adjoint of calculate del Z_T/Z_M and delP_T/delP_M on the grid
      call gvt_transform( statevector, 'ZandP_ad' )
    end if

    call gsv_deallocate( statevector_VarsLevs )

    if (slantPath_TO_tlad) call pressureProfileMonotonicityCheck(obsSpaceData, columnTrlOnAnlIncLev)

    call utl_tmg_stop(39)

    call utl_tmg_stop(30)

  end subroutine s2c_ad

  !---------------------------------------------------------
  ! s2c_nl
  !---------------------------------------------------------
  subroutine s2c_nl(stateVector, obsSpaceData, column, hco_core, timeInterpType, &
                   varName_opt, numObsBatches_opt, dealloc_opt, moveObsAtPole_opt, &
                   beSilent_opt)
    !
    ! :Purpose: Non-linear version of the horizontal interpolation,
    !           used for a full field (usually the background state when computing
    !           the innovation vector).
    !
    implicit none

    ! arguments
    type(struct_gsv)           :: stateVector
    type(struct_obs)           :: obsSpaceData
    type(struct_columnData)    :: column
    type(struct_hco), pointer  :: hco_core
    character(len=*)           :: timeInterpType
    character(len=*), optional :: varName_opt
    integer, optional          :: numObsBatches_opt
    logical, optional          :: dealloc_opt
    logical, optional          :: moveObsAtPole_opt
    logical, optional          :: beSilent_opt

    ! locals
    type(struct_gsv), save :: stateVector_VarsLevs 
    integer :: kIndex, kIndex2, kCount, stepIndex, numStep, mykEndExtended
    integer :: headerIndex, headerIndex2, numHeader, numHeaderMax, yourNumHeader
    integer :: headerIndexBeg, headerIndexEnd, obsBatchIndex, numObsBatches
    integer :: procIndex, nsize, ierr, headerUsedIndex, allHeaderIndexBeg(mmpi_nprocs)
    integer :: kIndexHeightSfc
    real(8) :: weight
    character(len=4)     :: varName
    real(8), pointer     :: column_ptr(:), ptr2d_r8(:,:), allCols_ptr(:,:)
    real(4), pointer     :: ptr4d_r4(:,:,:,:), ptr3d_UV_r4(:,:,:)
    real(8), allocatable :: cols_hint(:,:,:)
    real(8), allocatable :: cols_send(:,:)
    real(8), allocatable :: cols_recv(:,:)
    real(8), allocatable :: cols_send_1proc(:)
    integer, allocatable :: displs(:), nsizes(:)
    logical              :: dealloc, moveObsAtPole, rejectOutsideObs, beSilent
    logical, save        :: firstCall = .true.
    character(len=4), pointer :: varNames(:)

    call utl_tmg_start(30,'--StateToColumn')
    call utl_tmg_start(34,'----s2c_NL')

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .false.
    end if

    if ( .not. beSilent ) then
      write(*,*) 's2c_nl: STARTING'
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    end if

    if ( .not. gsv_isAllocated(stateVector) ) then 
      call utl_abort('s2c_nl: stateVector must be allocated')
    end if

    if (present(dealloc_opt)) then
      dealloc = dealloc_opt
    else
      dealloc = .true.
    end if

    ! determine number of obs batches (to reduce memory usage)
    if (present(numObsBatches_opt)) then
      numObsBatches = numObsBatches_opt
    else
      numObsBatches = 1
    end if
    if (.not. dealloc) then
      numObsBatches = 1 ! multiple batches only possible if dealloc=.true.
    end if

    if (interpInfo_nl%initialized) then
      if (.not. hco_equal(interpInfo_nl%hco,stateVector%hco) .or. numObsBatches > 1) then
        write(*,*) 's2c_nl: WARNING! Current hco grid parameters differ from allocated interpInfo!'
        write(*,*) 's2c_nl: InterpInfo will be deallocated.'
        call s2c_deallocInterpInfo(inputStateVectorType='nl')
      end if	
    end if

    if ( stateVector%mpi_distribution /= 'Tiles' ) then 
      call utl_abort('s2c_nl: stateVector must by Tiles distributed')
    end if

    if ( present(moveObsAtPole_opt) ) then
      moveObsAtPole = moveObsAtPole_opt
    else
      moveObsAtPole = .false.
    end if

    ! check the column and statevector have same nk/varNameList
    call checkColumnStatevectorMatch(column,statevector)

    ! calculate delP_T/delP_M and del Z_T/Z_M on the grid
    call gvt_transform( statevector, 'ZandP_nl' )

    if ( dealloc .or. firstCall ) then
      nullify(varNames)
      call gsv_varNamesList(varNames, statevector)
      call gsv_allocate( statevector_VarsLevs, stateVector%numstep, &
                         stateVector%hco, stateVector%vco, mpi_local_opt=.true., &
                         mpi_distribution_opt='VarsLevs', dataKind_opt=4, &
                         allocHeightSfc_opt=.true., varNames_opt=varNames )
      deallocate(varNames)
    else
      if (mmpi_myid == 0 .and. .not. beSilent) write(*,*) 's2c_nl: avoid re-allocating statevector_VarsLevs'
      call gsv_zero(statevector_VarsLevs)
    end if

    call gsv_transposeTilesToVarsLevs( stateVector, stateVector_VarsLevs, &
                                        beSilent_opt=beSilent )

    numStep = stateVector_VarsLevs%numStep

    if ( .not. interpInfo_nl%initialized ) then
      call utl_tmg_stop(34)
      call utl_tmg_start(31,'----s2c_Setups')
      ! also reject obs outside (LAM) domain and optionally move obs near 
      ! numerical pole to first/last analysis grid latitude
      call latlonChecksAnlGrid( obsSpaceData, hco_core, moveObsAtPole )

      ! Do not reject obs for global domain
      rejectOutsideObs = .not. stateVector_VarsLevs%hco%global
      write(*,*) 's2c_nl: rejectOutsideObs = ', rejectOutsideObs
      call utl_tmg_stop(31)
      call utl_tmg_start(34,'----s2c_NL')

    end if

    ! set contents of column to zero (1 variable or all)
    allCols_ptr => col_getAllColumns(column,varName_opt)
    if ( obs_numHeader(obsSpaceData) > 0 ) allCols_ptr(:,:) = 0.0d0

    OBSBATCH: do obsBatchIndex = 1, numObsBatches
      headerIndexBeg = 1 + (obsBatchIndex - 1) * (obs_numheader(obsSpaceData) / numObsBatches)
      if (obsBatchIndex == numObsBatches) then
        headerIndexEnd = obs_numheader(obsSpaceData)
      else
        headerIndexEnd = headerIndexBeg + (obs_numheader(obsSpaceData) / numObsBatches) - 1
      end if
      numHeader = headerIndexEnd - headerIndexBeg + 1
      call rpn_comm_allreduce(numHeader, numHeaderMax, 1,  &
                              'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)
      if ( .not. beSilent ) write(*,*) 's2c_nl: headerIndexBeg/End, numHeader = ',  &
                                        headerIndexBeg, headerIndexEnd, numHeader
      call rpn_comm_allgather(headerIndexBeg,   1,'mpi_integer', &
                              allHeaderIndexBeg,1,'mpi_integer','grid',ierr)

      if ( .not. interpInfo_nl%initialized ) then
        call utl_tmg_stop(34)
        call utl_tmg_start(31,'----s2c_Setups')

        ! compute and collect all obs grids onto all mpi tasks
        call s2c_setupInterpInfo( interpInfo_nl, obsSpaceData, stateVector_VarsLevs,  &
                                  headerIndexBeg, headerIndexEnd, &
                                  timeInterpType, rejectOutsideObs, &
                                  inputStateVectorType='nl', &
                                  lastCall_opt=(obsBatchIndex==numObsBatches))
        if ( mmpi_myid == 0 .and. verbose ) then
          do stepIndex = 1, numStep
            write(*,*) 's2c_nl: stepIndex, allNumHeaderUsed = ',  &
                       stepIndex, interpInfo_nl%allNumHeaderUsed(stepIndex,:)
          end do
        end if
        call utl_tmg_stop(31)
        call utl_tmg_start(34,'----s2c_NL')
      end if

      ! arrays for interpolated column for 1 level/variable and each time step
      allocate(cols_hint(maxval(interpInfo_nl%allNumHeaderUsed),numStep,mmpi_nprocs))
      cols_hint(:,:,:) = 0.0d0

      ! arrays for sending/receiving time interpolated column for 1 level/variable
      allocate(cols_send(numHeaderMax,mmpi_nprocs))
      cols_send(:,:) = 0.0d0

      allocate(cols_recv(numHeaderMax,mmpi_nprocs))
      cols_recv(:,:) = 0.0d0

      allocate(cols_send_1proc(numHeaderMax))

      call gsv_getField(stateVector_VarsLevs,ptr4d_r4)

      mykEndExtended = stateVector_VarsLevs%mykBeg + maxval(stateVector_VarsLevs%allkCount(:)) - 1

      kCount = 0
      k_loop: do kIndex = stateVector_VarsLevs%mykBeg, mykEndExtended
        kCount = kCount + 1

        if ( kIndex <= stateVector_VarsLevs%mykEnd ) then
          varName = gsv_getVarNameFromK(stateVector_VarsLevs,kIndex)

          call utl_tmg_start(35,'------s2c_NL_Hinterp')
          if ( varName == 'UU' .or. varName == 'VV' ) then
            call gsv_getFieldUV(stateVector_VarsLevs,ptr3d_UV_r4,kIndex)
          end if

          step_loop: do stepIndex = 1, numStep
            if ( maxval(interpInfo_nl%allNumHeaderUsed(stepIndex,:)) == 0 ) cycle step_loop

            ! interpolate to the columns destined for all procs for all steps and one lev/var
            !$OMP PARALLEL DO PRIVATE (procIndex, yourNumHeader)
            do procIndex = 1, mmpi_nprocs
              yourNumHeader = interpInfo_nl%allNumHeaderUsed(stepIndex,procIndex)
              if ( yourNumHeader > 0 ) then
                if ( varName == 'UU' ) then
                  call myezuvint_nl( cols_hint(1:yourNumHeader,stepIndex,procIndex), 'UU',  &
                                     ptr4d_r4(:,:,kIndex,stepIndex), ptr3d_UV_r4(:,:,stepIndex), &
                                     interpInfo_nl, kindex, stepIndex, procIndex )
                else if ( varName == 'VV' ) then
                  call myezuvint_nl( cols_hint(1:yourNumHeader,stepIndex,procIndex), 'VV',  &
                                     ptr3d_UV_r4(:,:,stepIndex), ptr4d_r4(:,:,kIndex,stepIndex), &
                                     interpInfo_nl, kindex, stepIndex, procIndex )
                else
                  call myezsint_nl( cols_hint(1:yourNumHeader,stepIndex,procIndex), &
                                    ptr4d_r4(:,:,kIndex,stepIndex),  &
                                    interpInfo_nl, kindex, stepIndex, procIndex )
                end if
              end if
            end do
            !$OMP END PARALLEL DO

          end do step_loop
          call utl_tmg_stop(35)

          ! interpolate in time to the columns destined for all procs and one level/variable
          do procIndex = 1, mmpi_nprocs
            cols_send_1proc(:) = 0.0d0
            do stepIndex = 1, numStep
              !$OMP PARALLEL DO PRIVATE (headerIndex, headerIndex2, headerUsedIndex, weight)
              do headerUsedIndex = 1, interpInfo_nl%allNumHeaderUsed(stepIndex, procIndex)
                headerIndex = interpInfo_nl%stepProcData(procIndex,stepIndex)%allHeaderIndex(headerUsedIndex)
                headerIndex2 = headerIndex - allHeaderIndexBeg(procIndex) + 1
                weight = oti_getTimeInterpWeightMpiGlobal(interpInfo_nl%oti,  &
                                                          headerIndex2,stepIndex,procIndex)
                cols_send_1proc(headerIndex2) = cols_send_1proc(headerIndex2) + &
                                                weight * cols_hint(headerUsedIndex,stepIndex,procIndex)
              end do
              !$OMP END PARALLEL DO
            end do
            cols_send(:,procIndex) = cols_send_1proc(:)
          end do

        else

          ! this value of k does not exist on this mpi task
          cols_send(:,:) = 0.0d0

        end if ! if kIndex <= mykEnd

        call rpn_comm_barrier('GRID',ierr)

        call utl_tmg_start(36,'------s2c_NL_allToAll')
        ! mpi communication: alltoall for one level/variable
        nsize = numHeaderMax
        if(mmpi_nprocs > 1) then
          call rpn_comm_alltoall(cols_send, nsize, 'MPI_REAL8',  &
                                 cols_recv, nsize, 'MPI_REAL8', 'GRID', ierr)
        else
          cols_recv(:,1) = cols_send(:,1)
        end if
        call utl_tmg_stop(36)
	
        ! reorganize ensemble of distributed columns
        !$OMP PARALLEL DO PRIVATE (procIndex, kIndex2, headerIndex, headerIndex2)
        proc_loop: do procIndex = 1, mmpi_nprocs
          kIndex2 = stateVector_VarsLevs%allkBeg(procIndex) + kCount - 1
          if ( kIndex2 > stateVector_VarsLevs%allkEnd(procIndex) ) cycle proc_loop
          do headerIndex = 1, numHeader
            headerIndex2 = headerIndex + headerIndexBeg - 1
            allCols_ptr(kIndex2,headerIndex2) = cols_recv(headerIndex,procIndex)
          end do
        end do proc_loop
        !$OMP END PARALLEL DO

      end do k_loop

      ! impose a lower limit on HU
      if(col_varExist(column,'HU')) then
        do headerIndex = headerIndexBeg, headerIndexEnd
          column_ptr => col_getColumn(column,headerIndex,'HU')
          column_ptr(:) = max(column_ptr(:),col_rhumin)
        end do
      end if

      ! impose a lower/upper limits on LWCR
      if( col_varExist(column,'LWCR') ) then
        do headerIndex = headerIndexBeg, headerIndexEnd
          column_ptr => col_getColumn(column,headerIndex,'LWCR')
          column_ptr(:) = max(column_ptr(:),qlim_readMinClwValue())
          column_ptr(:) = min(column_ptr(:),qlim_readMaxClwValue())
        end do
      end if

      ! Interpolate surface height separately, only exists on mpi task 0
      HeightSfcPresent: if ( stateVector_VarsLevs%HeightSfcPresent ) then

        if ( mmpi_myid == 0 ) then
          varName = 'GZ'
          kIndexHeightSfc = 0    
          step_loop_height: do stepIndex = 1, numStep

            if ( maxval(interpInfo_nl%allNumHeaderUsed(stepIndex,:)) == 0 ) cycle step_loop_height

            ! interpolate to the columns destined for all procs for all steps and one lev/var
            !$OMP PARALLEL DO PRIVATE (procIndex, yourNumHeader, ptr2d_r8)
            do procIndex = 1, mmpi_nprocs
              yourNumHeader = interpInfo_nl%allNumHeaderUsed(stepIndex,procIndex)
              if ( yourNumHeader > 0 ) then
                ptr2d_r8 => gsv_getHeightSfc(stateVector_VarsLevs)
                call myezsint_r8_nl( cols_hint(1:yourNumHeader,stepIndex,procIndex), &
                                     ptr2d_r8(:,:), interpInfo_nl, kIndexHeightSfc, stepIndex, procIndex )
              end if
            end do
            !$OMP END PARALLEL DO

          end do step_loop_height

          ! interpolate in time to the columns destined for all procs and one level/variable
          do procIndex = 1, mmpi_nprocs
            cols_send(:,procIndex) = 0.0d0
            do stepIndex = 1, numStep
              !$OMP PARALLEL DO PRIVATE (headerIndex, headerIndex2, headerUsedIndex)
              do headerUsedIndex = 1, interpInfo_nl%allNumHeaderUsed(stepIndex, procIndex)
                headerIndex = interpInfo_nl%stepProcData(procIndex, stepIndex)%allHeaderIndex(headerUsedIndex)
                ! just copy, since surface height same for all time steps
                headerIndex2 = headerIndex - allHeaderIndexBeg(procIndex) + 1
                cols_send(headerIndex2,procIndex) = cols_hint(headerUsedIndex,stepIndex,procIndex)
              end do
              !$OMP END PARALLEL DO
            end do
          end do

        end if

        ! mpi communication: scatter data from task 0
        nsize = numHeaderMax
        if(mmpi_nprocs > 1) then
          allocate(displs(mmpi_nprocs))
          allocate(nsizes(mmpi_nprocs))
          do procIndex = 1, mmpi_nprocs
            displs(procIndex) = (procIndex - 1) * nsize
            nsizes(procIndex) = nsize
          end do
          call rpn_comm_scatterv(cols_send, nsizes, displs, 'MPI_REAL8', &
                                 cols_recv, nsize, 'MPI_REAL8', &
                                 0, 'GRID', ierr)
          deallocate(displs)
          deallocate(nsizes)

        else
          cols_recv(:,1) = cols_send(:,1)
        end if

        do headerIndex = headerIndexBeg, headerIndexEnd
          headerIndex2 = headerIndex - headerIndexBeg + 1
          call col_setHeightSfc(column, headerIndex, cols_recv(headerIndex2,1))
        end do

      end if HeightSfcPresent

      deallocate(cols_hint)
      deallocate(cols_send)
      deallocate(cols_recv)
      deallocate(cols_send_1proc)

      if ( dealloc ) call s2c_deallocInterpInfo( inputStateVectorType='nl' )

    end do OBSBATCH

    if ( dealloc) call gsv_deallocate( statevector_VarsLevs )

    if (slantPath_TO_nl) call pressureProfileMonotonicityCheck(obsSpaceData, column)

    firstCall = .false.

    if ( .not. beSilent ) then
      write(*,*) 's2c_nl: FINISHED'
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    end if

    call utl_tmg_stop(34)

    call utl_tmg_stop(30)

  end subroutine s2c_nl

  ! -------------------------------------------------
  ! myezsint_nl: Scalar field horizontal interpolation
  ! -------------------------------------------------
  subroutine myezsint_nl( column_out, field_in, interpInfo, kIndex, stepIndex, procIndex )
    !
    ! :Purpose: Scalar horizontal interpolation, replaces the
    !           ezsint routine from rmnlib.
    !
    implicit none

    ! arguments
    real(8)                 :: column_out(:)
    real(4)                 :: field_in(:,:)
    type(struct_interpInfo) :: interpInfo
    integer                 :: stepIndex
    integer                 :: procIndex
    integer                 :: kIndex

    ! locals
    integer :: lonIndex, latIndex, gridptIndex, headerIndex, subGridIndex, numColumn
    real(8) :: interpValue, weight

    numColumn = size( column_out )

    do headerIndex = 1, numColumn

      ! Interpolate the model state to the obs point
      interpValue = 0.0d0

      do subGridIndex = 1, interpInfo%hco%numSubGrid

        do gridptIndex =  &
             interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, kIndex), &
             interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex, headerIndex, kIndex)

          lonIndex = interpInfo%lonIndexDepot(gridptIndex)
          latIndex = interpInfo%latIndexDepot(gridptIndex)
          weight = interpInfo%interpWeightDepot(gridptIndex)

          interpValue = interpValue + weight * real(field_in(lonIndex, latIndex),8)

        end do

      end do
      column_out(headerIndex) = interpValue

    end do

  end subroutine myezsint_nl

  ! -------------------------------------------------
  ! myezsint_r8_nl: Scalar field horizontal interpolation
  ! -------------------------------------------------
  subroutine myezsint_r8_nl( column_out, field_in, interpInfo, kIndex, stepIndex, procIndex )
    !
    ! :Purpose: Scalar horizontal interpolation, replaces the
    !           ezsint routine from rmnlib.
    !
    implicit none

    ! arguments
    real(8)                 :: column_out(:)
    real(8)                 :: field_in(:,:)
    type(struct_interpInfo) :: interpInfo
    integer                 :: stepIndex
    integer                 :: procIndex
    integer                 :: kIndex

    ! locals
    integer :: lonIndex, latIndex, gridptIndex, headerIndex, subGridIndex, numColumn
    real(8) :: interpValue, weight

    numColumn = size( column_out )

    do headerIndex = 1, numColumn

      ! Interpolate the model state to the obs point
      interpValue = 0.0d0

      do subGridIndex = 1, interpInfo%hco%numSubGrid

        do gridptIndex =  &
             interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, kIndex), &
             interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex, headerIndex, kIndex)

          lonIndex = interpInfo%lonIndexDepot(gridptIndex)
          latIndex = interpInfo%latIndexDepot(gridptIndex)
          weight = interpInfo%interpWeightDepot(gridptIndex)

          interpValue = interpValue + weight * field_in(lonIndex, latIndex)

        end do

      end do
      column_out(headerIndex) = interpValue

    end do

  end subroutine myezsint_r8_nl

  ! -------------------------------------------------
  ! myezsint_tl: Scalar field horizontal interpolation
  ! -------------------------------------------------
  subroutine myezsint_tl( column_out, field_in, interpInfo, kIndex, stepIndex, procIndex )
    !
    ! :Purpose: Scalar horizontal interpolation, replaces the
    !           ezsint routine from rmnlib.
    !
    implicit none

    ! arguments
    real(8)                 :: column_out(:)
    real(pre_incrReal)      :: field_in(:,:)
    type(struct_interpInfo) :: interpInfo
    integer                 :: stepIndex
    integer                 :: procIndex
    integer                 :: kIndex

    ! locals
    integer :: lonIndex, latIndex, gridptIndex, headerIndex, subGridIndex, numColumn
    real(8) :: interpValue, weight

    numColumn = size( column_out )

    do headerIndex = 1, numColumn

      ! Interpolate the model state to the obs point
      interpValue = 0.0d0

      do subGridIndex = 1, interpInfo%hco%numSubGrid

        do gridptIndex =  &
             interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, kIndex), &
             interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex, headerIndex, kIndex)

          lonIndex = interpInfo%lonIndexDepot(gridptIndex)
          latIndex = interpInfo%latIndexDepot(gridptIndex)
          weight = interpInfo%interpWeightDepot(gridptIndex)

          interpValue = interpValue + weight * field_in(lonIndex, latIndex)

        end do

      end do
      column_out(headerIndex) = interpValue

    end do

  end subroutine myezsint_tl

  ! -------------------------------------------------------------
  ! myezsint_ad: Adjoint of scalar field horizontal interpolation
  ! -------------------------------------------------------------
  subroutine myezsint_ad( column_in, field_out, interpInfo, kIndex, stepIndex, procIndex )
    !
    ! :Purpose: Adjoint of the scalar horizontal interpolation.
    !
    implicit none

    ! Arguments:
    real(8)                 :: column_in(:)
    real(pre_incrReal)      :: field_out(:,:)
    type(struct_interpInfo) :: interpInfo
    integer                 :: stepIndex
    integer                 :: procIndex
    integer                 :: kIndex

    ! Locals:
    integer :: lonIndex, latIndex, gridptIndex, headerIndex, subGridIndex, numColumn
    real(8) :: weight

    numColumn = size( column_in )

    do headerIndex = 1, numColumn

      ! Interpolate the model state to the obs point

      do subGridIndex = 1, interpInfo%hco%numSubGrid

        do gridptIndex =  &
             interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, kIndex), &
             interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex, headerIndex, kIndex)

          lonIndex = interpInfo%lonIndexDepot(gridptIndex)
          latIndex = interpInfo%latIndexDepot(gridptIndex)
          weight = interpInfo%interpWeightDepot(gridptIndex)

          field_out(lonIndex, latIndex) = field_out(lonIndex, latIndex) +  &
                                          weight * column_in(headerIndex)

        end do

      end do

    end do

  end subroutine myezsint_ad

  ! -------------------------------------------------------------
  ! myezuvint_nl: Vector field horizontal interpolation
  ! -------------------------------------------------------------
  subroutine myezuvint_nl( column_out, varName, fieldUU_in, fieldVV_in,  &
                           interpInfo, kIndex, stepIndex, procIndex )
    !
    ! :Purpose: Vector horizontal interpolation, replaces the
    !           ezuvint routine from rmnlib.
    !
    implicit none

    ! arguments
    real(8)                 :: column_out(:)
    character(len=*)        :: varName
    real(4)                 :: fieldUU_in(:,:)
    real(4)                 :: fieldVV_in(:,:)
    type(struct_interpInfo) :: interpInfo
    integer                 :: stepIndex
    integer                 :: procIndex
    integer                 :: kIndex

    ! locals
    integer :: lonIndex, latIndex, indexBeg, indexEnd, gridptIndex, headerIndex
    integer :: numColumn, subGridIndex
    real(8) :: interpUU(interpInfo%hco%numSubGrid), interpVV(interpInfo%hco%numSubGrid)
    real(8) :: lat, lon, latRot, lonRot, weight
    logical :: doUU, doVV

    numColumn = size( column_out )

    doUU = (trim(varName) == 'UU' .or. interpInfo%hco%rotated)
    doVV = (trim(varName) == 'VV' .or. interpInfo%hco%rotated)

    header_loop: do headerIndex = 1, numColumn

      interpUU(:) = 0.0d0
      interpVV(:) = 0.0d0

      subGrid_loop: do subGridIndex = 1, interpInfo%hco%numSubGrid

        indexBeg = interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, kIndex)
        indexEnd = interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex, headerIndex, kIndex)

        if ( indexEnd < IndexBeg ) cycle subGrid_loop

        ! Interpolate the model UU to the obs point
        do gridptIndex = indexBeg, indexEnd

          lonIndex = interpInfo%lonIndexDepot(gridptIndex)
          latIndex = interpInfo%latIndexDepot(gridptIndex)
          weight = interpInfo%interpWeightDepot(gridptIndex)

          if ( doUU ) interpUU(subGridIndex) = interpUU(subGridIndex) +  &
                      weight * real(fieldUU_in(lonIndex, latIndex),8)
          if ( doVV ) interpVV(subGridIndex) = interpVV(subGridIndex) +  &
                      weight * real(fieldVV_in(lonIndex, latIndex),8)

        end do
        ! now rotate the wind vector
        if ( interpInfo%hco%rotated ) then
          lat = interpInfo%stepProcData(procIndex, stepIndex)%allLat(headerIndex, kIndex)
          lon = interpInfo%stepProcData(procIndex, stepIndex)%allLon(headerIndex, kIndex)
          latRot = interpInfo%stepProcData(procIndex,stepIndex)%allLatRot(subGridIndex, headerIndex, kIndex)
          lonRot = interpInfo%stepProcData(procIndex,stepIndex)%allLonRot(subGridIndex, headerIndex, kIndex)

          call uvr_rotateWind_nl( interpInfo%uvr,            & ! IN
                                  subGridIndex,              & ! IN
                                  interpUU(subGridIndex),    & ! INOUT
                                  interpVV(subGridIndex),    & ! INOUT
                                  lat, lon, latRot, lonRot,  & ! IN
                                  'ToMetWind' )                ! IN
        end if

      end do subGrid_loop

      ! return only the desired component
      if ( trim(varName) == 'UU' ) then
        column_out(headerIndex) = sum(interpUU(:))
      else
        column_out(headerIndex) = sum(interpVV(:))
      end if

    end do header_loop

  end subroutine myezuvint_nl

  ! -------------------------------------------------------------
  ! myezuvint_tl: Vector field horizontal interpolation
  ! -------------------------------------------------------------
  subroutine myezuvint_tl( column_out, varName, fieldUU_in, fieldVV_in,  &
                           interpInfo, kIndex, stepIndex, procIndex )
    ! :Purpose: Vector horizontal interpolation, replaces the
    !           ezuvint routine from rmnlib.
    !
    implicit none

    ! arguments
    real(8)                 :: column_out(:)
    character(len=*)        :: varName
    real(pre_incrReal)      :: fieldUU_in(:,:)
    real(pre_incrReal)      :: fieldVV_in(:,:)
    type(struct_interpInfo) :: interpInfo
    integer                 :: stepIndex
    integer                 :: procIndex
    integer                 :: kIndex

    ! locals
    integer :: lonIndex, latIndex, indexBeg, indexEnd, gridptIndex, headerIndex
    integer :: numColumn, subGridIndex
    real(8) :: interpUU(interpInfo%hco%numSubGrid), interpVV(interpInfo%hco%numSubGrid)
    real(8) :: lat, lon, latRot, lonRot, weight
    logical :: doUU, doVV

    numColumn = size( column_out )

    doUU = (trim(varName) == 'UU' .or. interpInfo%hco%rotated)
    doVV = (trim(varName) == 'VV' .or. interpInfo%hco%rotated)

    header_loop: do headerIndex = 1, numColumn

      interpUU(:) = 0.0d0
      interpVV(:) = 0.0d0

      subGrid_loop: do subGridIndex = 1, interpInfo%hco%numSubGrid

        indexBeg = interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, kIndex)
        indexEnd = interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex, headerIndex, kIndex)

        if ( indexEnd < IndexBeg ) cycle subGrid_loop

        ! Interpolate the model UU to the obs point
        do gridptIndex = indexBeg, indexEnd

          lonIndex = interpInfo%lonIndexDepot(gridptIndex)
          latIndex = interpInfo%latIndexDepot(gridptIndex)
          weight = interpInfo%interpWeightDepot(gridptIndex)

          if ( doUU ) interpUU(subGridIndex) = interpUU(subGridIndex) +  &
                      weight * fieldUU_in(lonIndex, latIndex)
          if ( doVV ) interpVV(subGridIndex) = interpVV(subGridIndex) +  &
                      weight * fieldVV_in(lonIndex, latIndex)

        end do
        ! now rotate the wind vector
        if ( interpInfo%hco%rotated ) then
          lat = interpInfo%stepProcData(procIndex, stepIndex)%allLat(headerIndex, kIndex)
          lon = interpInfo%stepProcData(procIndex, stepIndex)%allLon(headerIndex, kIndex)
          latRot = interpInfo%stepProcData(procIndex,stepIndex)%allLatRot(subGridIndex, headerIndex, kIndex)
          lonRot = interpInfo%stepProcData(procIndex,stepIndex)%allLonRot(subGridIndex, headerIndex, kIndex)

          call uvr_rotateWind_tl( interpInfo%uvr,            & ! IN
                                  subGridIndex,              & ! IN
                                  interpUU(subGridIndex),    & ! INOUT
                                  interpVV(subGridIndex),    & ! INOUT
                                  lat, lon, latRot, lonRot,  & ! IN
                                  'ToMetWind' )                ! IN
        end if

      end do subGrid_loop

      ! return only the desired component
      if ( trim(varName) == 'UU' ) then
        column_out(headerIndex) = sum(interpUU(:))
      else
        column_out(headerIndex) = sum(interpVV(:))
      end if

    end do header_loop

  end subroutine myezuvint_tl

  ! -------------------------------------------------------------
  ! myezuvint_ad: Adjoint of vector field horizontal interpolation
  ! -------------------------------------------------------------
  subroutine myezuvint_ad( column_in, varName, fieldUU_out, fieldVV_out, &
                           interpInfo, kIndex, stepIndex, procIndex )
    !
    ! :Purpose: Adjoint of the vector horizontal interpolation.
    !
    implicit none

    ! arguments
    real(8)                 :: column_in(:)
    character(len=*)        :: varName
    real(pre_incrReal)      :: fieldUU_out(:,:)
    real(pre_incrReal)      :: fieldVV_out(:,:)
    type(struct_interpInfo) :: interpInfo
    integer                 :: stepIndex
    integer                 :: procIndex
    integer                 :: kIndex

    ! locals
    integer :: lonIndex, latIndex, indexBeg, indexEnd, gridptIndex, headerIndex
    integer :: numColumn, subGridIndex
    real(8) :: interpUU(interpInfo%hco%numSubGrid), interpVV(interpInfo%hco%numSubGrid)
    real(8) :: lat, lon, latRot, lonRot, weight
    logical :: doUU, doVV

    numColumn = size( column_in )

    doUU = (trim(varName) == 'UU' .or. interpInfo%hco%rotated)
    doVV = (trim(varName) == 'VV' .or. interpInfo%hco%rotated)

    header_loop: do headerIndex = 1, numColumn

      if ( trim(varName) == 'UU' ) then
        interpUU(:) = column_in(headerIndex)
        interpVV(:) = 0.0d0
      else
        interpUU(:) = 0.0d0
        interpVV(:) = column_in(headerIndex)
      end if

      subGrid_loop: do subGridIndex = 1, interpInfo%hco%numSubGrid

        indexBeg = interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, kIndex)
        indexEnd = interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex, headerIndex, kIndex)

        if ( indexEnd < IndexBeg ) cycle subGrid_loop

        ! now rotate the wind vector and return the desired component
        if ( interpInfo%hco%rotated ) then
          lat = interpInfo%stepProcData(procIndex, stepIndex)%allLat(headerIndex, kIndex)
          lon = interpInfo%stepProcData(procIndex, stepIndex)%allLon(headerIndex, kIndex)
          latRot = interpInfo%stepProcData(procIndex,stepIndex)%allLatRot(subGridIndex, headerIndex, kIndex)
          lonRot = interpInfo%stepProcData(procIndex,stepIndex)%allLonRot(subGridIndex, headerIndex, kIndex)

          call uvr_rotateWind_ad( interpInfo%uvr,           & ! IN 
                                  subGridIndex,             & ! IN
                                  interpUU(subGridIndex),   & ! INOUT
                                  interpVV(subGridIndex),   & ! INOUT
                                  lat, lon, latRot, lonRot, & ! IN
                                  'ToMetWind' )               ! IN
        end if

        ! Interpolate the model VV to the obs point
        do gridptIndex = indexBeg, indexEnd

          lonIndex = interpInfo%lonIndexDepot(gridptIndex)
          latIndex = interpInfo%latIndexDepot(gridptIndex)
          weight = interpInfo%interpWeightDepot(gridptIndex)

          if ( doUU ) fieldUU_out(lonIndex, latIndex) =  &
                      fieldUU_out(lonIndex, latIndex) + weight * interpUU(subGridIndex)
          if ( doVV ) fieldVV_out(lonIndex, latIndex) =  &
                      fieldVV_out(lonIndex, latIndex) + weight * interpVV(subGridIndex)

        end do

      end do subGrid_loop

    end do header_loop

  end subroutine myezuvint_ad

  !---------------------------------------------------------
  ! s2c_bgcheck_bilin
  !---------------------------------------------------------
  subroutine s2c_bgcheck_bilin(column,statevector,obsSpaceData)
    !
    ! :Purpose: Special version of s2c_tl used for background check. This should
    !           be replaced by direct call to s2c_tl. It is not general enough to
    !           be used for new analysis variables.
    !
    implicit none

    ! arguments
    type(struct_columnData) :: column
    type(struct_gsv) :: statevector
    type(struct_obs) :: obsSpaceData

    ! locals
    integer :: jk, jk2, jgl, headerIndex
    integer :: lonIndex, ila, ierr, subGridIndex
    integer :: extraLongitude
    real(8) :: lat, lon
    real(4) :: lat_r4, lon_r4, lat_deg_r4, lon_deg_r4, xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    real(8) :: dldy, dlw1, dlw2, dlw3, dlw4, dldx, ypos, xpos
    real(8), allocatable ::zgd(:,:,:)
    real(8), pointer :: field_ptr(:,:,:)

    real(8), pointer :: varColumn(:)
    integer :: varIndex
    character(len=4) :: varName

    call utl_tmg_start(30,'--StateToColumn')

    ! Note: We assume here the all the obs between the poles and the last grid points
    !       (i.e. outside the grid) have been moved within the grid by suprep

    if (statevector%hco%global) then
      extraLongitude = 1
    else
      extraLongitude = 0
    end if

    allocate(zgd(statevector%ni+extraLongitude,statevector%nj,statevector%nk))
  
    zgd(:,:,:)=0.0d0
    call gsv_getField(statevector,field_ptr)
    zgd(1:statevector%ni,1:statevector%nj,1:statevector%nk)= &
         field_ptr(1:statevector%ni,1:statevector%nj,1:statevector%nk)

    !
    !- 1.  Expand field by repeating meridian 1 into into meridian ni+1
    !
    if (extraLongitude == 1) then
      do jk = 1, statevector%nk
        do jgl = 1, statevector%nj
          zgd(statevector%ni+1,jgl,jk) = zgd( 1,jgl,jk)
        end do
      end do
    end if

    !
    !- 2.  Loop over all the headers
    !
    do headerIndex = 1, col_getNumCol(column)

      !- 2.1 Find the obs positin within the analysis grid
      lat    = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)
      lon    = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex)
      lat_r4 = real(lat,4)
      lon_r4 = real(lon,4)
      if (lon_r4.lt.0.0         ) lon_r4 = lon_r4 + 2.0*MPC_PI_R4
      if (lon_r4.ge.2.*MPC_PI_R4) lon_r4 = lon_r4 - 2.0*MPC_PI_R4
      lat_deg_r4 = lat_r4 * MPC_DEGREES_PER_RADIAN_R4 ! Radian To Degree
      lon_deg_r4 = lon_r4 * MPC_DEGREES_PER_RADIAN_R4
      ierr = gpos_getPositionXY( stateVector % hco % EZscintID,   &
                                xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                                lat_deg_r4, lon_deg_r4, subGridIndex )
      xpos = real(xpos_r4,8)
      ypos = real(ypos_r4,8)

      !- Make sure we are within bounds
      if ( ypos < 1.d0                        .or. &
           ypos > real(statevector%nj    , 8) .or. &
           xpos < 1.d0                        .or. &
           xpos > real(statevector%ni + extraLongitude, 8) ) then
        write(*,*) 's2c_bgcheck_bilin: Obs outside local domain for headerIndex = ', &
                   headerIndex
        write(*,*) '  obs    lat, lon position            = ', &
                   Lat*MPC_DEGREES_PER_RADIAN_R8, Lon*MPC_DEGREES_PER_RADIAN_R8
        write(*,*) '  obs    x, y     position            = ', &
                   xpos, ypos
        write(*,*) '  domain x_end, y_end bounds          = ', &
                   statevector%ni + extraLongitude, statevector%nj
        call utl_abort('s2c_bgcheck_bilin')
      end if

      !- 2.2 Find the lower-left grid point next to the observation
      if ( xpos == real(statevector%ni + extraLongitude,8) ) then
        lonIndex = floor(xpos) - 1
      else
        lonIndex = floor(xpos)
      end if

      if ( ypos == real(statevector%nj,8) ) then
        ILA = floor(ypos) - 1
      else
        ILA = floor(ypos)
      end if

      !- 2.3 Compute the 4 weights of the bilinear interpolation
      dldx = xpos - real(lonIndex,8)
      dldy = ypos - real(ILA,8)

      dlw1 = (1.d0-dldx) * (1.d0-dldy)
      dlw2 =       dldx  * (1.d0-dldy)
      dlw3 = (1.d0-dldx) *       dldy
      dlw4 =       dldx  *       dldy

      !- 2.4 Interpolate the model state to the obs point
           
      do varIndex = 1, vnl_numvarmax
        if (.not. col_varExist(column,trim(vnl_varNameList(varIndex)))) cycle
        varName=trim(vnl_varNameList(varIndex))
        varColumn => col_getColumn(column,headerIndex,varName)
        
        if(gsv_varExist(statevector,varName)) then
          do jk = 1, gsv_getNumLevFromVarName(statevector,varName)      
              jk2=jk+gsv_getOffsetFromVarName(statevector,varName)
              varColumn(jk) =   dlw1*zgd(lonIndex  ,ila,jk2)  &
                                + dlw2*zgd(lonIndex+1,ila,jk2)  &
                                + dlw3*zgd(lonIndex  ,ila+1,jk2)  &
                                + dlw4*zgd(lonIndex+1,ila+1,jk2)
          end do
        end if
        
        nullify(varColumn)
      end do
      
    end do

    deallocate(zgd)

    call utl_tmg_stop(30)

  end subroutine s2c_bgcheck_bilin

  !--------------------------------------------------------------------------
  ! s2c_setupHorizInterp
  !--------------------------------------------------------------------------
  subroutine s2c_setupHorizInterp(footprintRadius_r4, interpInfo, &
                                  stateVector, headerIndex, kIndex, stepIndex, &
                                  procIndex, numGridpt)
    !
    !:Purpose: To identify the appropriate horizontal interpolation scheme based
    !          on footprint radius value. Then to call the corresponding
    !          subroutine to determine the grid points and their associated
    !          weights.
    !
    implicit none

    ! Arguments:
    real(4)                , intent(in)    :: footprintRadius_r4 ! (metres)
    type(struct_interpInfo), intent(inout) :: interpInfo
    type(struct_gsv)       , intent(in)    :: stateVector
    integer                , intent(in)    :: headerIndex, kIndex, stepIndex
    integer                , intent(in)    :: procIndex
    integer                , intent(out)   :: numGridpt(interpInfo%hco%numSubGrid)

    if ( footprintRadius_r4 > 0.0 ) then

      call s2c_setupFootprintInterp(footprintRadius_r4, interpInfo, stateVector, headerIndex, kIndex, stepIndex, procIndex, numGridpt)

    else if ( footprintRadius_r4 == bilinearFootprint ) then

      call s2c_setupBilinearInterp(interpInfo, stateVector, headerIndex, kIndex, stepIndex, procIndex, numGridpt)

    else if ( footprintRadius_r4 == lakeFootprint ) then

      call s2c_setupLakeInterp(interpInfo, stateVector, headerIndex, kIndex, stepIndex, procIndex, numGridpt)

    else if ( footprintRadius_r4 == nearestNeighbourFootprint ) then

      call s2c_setupNearestNeighbor(interpInfo, stateVector, headerIndex, kIndex, stepIndex, procIndex, numGridpt)

    else

      write(*,*) 'footprint radius = ',footprintRadius_r4
      call utl_abort('s2c_setupHorizInterp: footprint radius not permitted')

    end if

  end subroutine s2c_setupHorizInterp

  !------------------------------------------------------------------
  ! s2c_getFootprintRadius
  !------------------------------------------------------------------
  function s2c_getFootprintRadius( obsSpaceData, stateVector, headerIndex ) result(fpr)
    !
    !:Purpose: To determine the footprint radius (metres) of the observation.
    !          In the case of bilinear horizontal interpolation,
    !          the returned footprint is zero (default).
    !
    implicit none
    real(4)                       :: fpr

    ! Arguments:
    type(struct_obs), intent(in)  :: obsSpaceData
    type(struct_gsv), intent(in)  :: stateVector
    integer         , intent(in)  :: headerIndex

    ! locals
    character(len=2)  :: obsFamily
    character(len=12) :: cstnid
    integer           :: codeType

    fpr = bilinearFootprint

    obsFamily = obs_getFamily ( obsSpaceData, headerIndex )
    if ( obsFamily == 'GL' ) then

      cstnid = obs_elem_c ( obsSpaceData, 'STID' , headerIndex )
      codeType = obs_headElem_i( obsSpaceData, OBS_ITY, headerIndex )

      if (index(cstnid,'DMSP') == 1) then

        select case(cstnid)
        case('DMSP15')
          fpr = 27.5e3
        case('DMSP16','DMSP17','DMSP18')
          fpr = 29.0e3
        case DEFAULT
          call utl_abort('s2c_getFootprintRadius: UNKNOWN station id: '//cstnid)
        end select

      else if (cstnid == 'GCOM-W1') then

        fpr = 11.0e3

      else if (cstnid(1:6) == 'METOP-') then

        fpr = 25.0e3

      else if (cstnid == 'noaa-19') then

        fpr = 2.75e3

      else if (cstnid == 'CIS_DAILY') then

        fpr = bilinearFootprint

      else if (cstnid == 'RS1_IMG') then

        fpr = bilinearFootprint

      else if (codtyp_get_name(codeType) == 'iceclake') then

        fpr = lakeFootprint

      else if (cstnid == 'CIS_REGIONAL') then

        fpr = bilinearFootprint

      else

        call utl_abort('s2c_getFootprintRadius: UNKNOWN station id: '//cstnid)

      end if

    else if (obsFamily == 'HY') then

      fpr = nearestNeighbourFootprint

    else if (obsFamily == 'TO' .and. useFootprintForTovs ) then

      fpr = getTovsFootprintRadius(obsSpaceData, headerIndex, beSilent_opt=.true.)

      ! As safety margin, add 10% to maxGridSpacing before comparing to the footprint radius.
      if ( fpr < 1.1 * real(stateVector%hco%maxGridSpacing,4) ) fpr = bilinearFootprint

    else

      fpr = bilinearFootprint

    end if

  end function s2c_getFootprintRadius

  !--------------------------------------------------------------------------
  ! s2c_rejectZeroWeightObs
  !--------------------------------------------------------------------------
  subroutine s2c_rejectZeroWeightObs(interpInfo, obsSpaceData, mykBeg, mykEnd)
    !
    !:Purpose: To flag an observation in obsSpaceData as being rejected if
    !          it has zero interpolation weight (usually because an ocean
    !          obs is touching land) on any mpi task.
    !
    implicit none    

    ! Arguments:
    type(struct_interpInfo), intent(inout) :: interpInfo
    type(struct_obs)       , intent(inout) :: obsSpaceData
    integer                , intent(in)    :: mykBeg
    integer                , intent(in)    :: mykEnd

    ! Locals:
    integer :: numStep, procIndex, stepIndex, headerUsedIndex, headerIndex, kIndex
    integer :: numHeader, numHeaderMax, bodyIndexBeg, bodyIndexEnd, bodyIndex
    integer :: subGridIndex, gridptIndex, ierr, nsize
    integer, save :: numWrites = 0
    logical, allocatable :: allRejectObs(:,:), allRejectObsMpiGlobal(:,:)

    write(*,*) 's2c_rejectZeroWeightObs: Starting'

    numHeader = obs_numheader(obsSpaceData)
    call rpn_comm_allreduce(numHeader, numHeaderMax, 1,  &
                            'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)

    allocate(allRejectObs(numHeaderMax,mmpi_nprocs))
    allocate(allRejectObsMpiGlobal(numHeaderMax,mmpi_nprocs))
    allRejectObs(:,:) = .false.
    allRejectObsMpiGlobal(:,:) = .false.

    numStep = size(interpInfo%stepProcData(1,:))
    do procIndex = 1, mmpi_nprocs
      do stepIndex = 1, numStep
        do headerUsedIndex = 1, interpInfo%allNumHeaderUsed(stepIndex,procIndex)
          headerIndex = interpInfo%stepProcData(procIndex,stepIndex)%allHeaderIndex(headerUsedIndex)
          do kIndex = mykBeg, mykEnd
            if (kIndex == mykBeg) allRejectObs(headerIndex,procIndex) = .true.
            do subGridIndex = 1, interpInfo%hco%numSubGrid
              do gridptIndex =  &
                   interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerUsedIndex, kIndex), &
                   interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex, headerUsedIndex, kIndex)
                if (interpInfo%interpWeightDepot(gridptIndex) > 0.0d0) then
                  allRejectObs(headerIndex,procIndex) = .false.
                end if
              end do
            end do
          end do ! kIndex
        end do ! headerUsedIndex
      end do ! stepIndex
    end do ! procIndex

    ! do global communication of reject flags
    nsize = numHeaderMax*mmpi_nprocs
    call rpn_comm_allreduce(allRejectObs,allRejectObsMpiGlobal,nsize,'MPI_LOGICAL','MPI_LOR','GRID',ierr)

    ! modify obsSpaceData based on reject flags
    do headerIndex = 1, obs_numHeader(obsSpaceData)
      if (allRejectObsMpiGlobal(headerIndex,mmpi_myid+1)) then

        numWrites = numWrites + 1
        if (numWrites < maxNumWrites) then
          write(*,*) 's2c_rejectZeroWeightObs: Rejecting OBS with zero weight, index ', headerIndex
        else if (numWrites == maxNumWrites) then
          write(*,*) 's2c_rejectZeroWeightObs: More rejects, but reached maximum number of writes to the listing.'
        end if

        bodyIndexBeg = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
        bodyIndexEnd = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) + bodyIndexBeg -1
        do bodyIndex = bodyIndexBeg, bodyIndexEnd
          call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
        end do
        call obs_headSet_i(obsSpaceData, OBS_ST1, headerIndex,  &
                    ibset( obs_headElem_i(obsSpaceData, OBS_ST1, headerIndex), 05))
      end if
    end do

    deallocate(allRejectObs)
    deallocate(allRejectObsMpiGlobal)

    write(*,*) 's2c_rejectZeroWeightObs: Finished'

  end subroutine s2c_rejectZeroWeightObs

  !--------------------------------------------------------------------------
  ! s2c_setupBilinearInterp
  !--------------------------------------------------------------------------
  subroutine s2c_setupBilinearInterp(interpInfo, stateVector, headerIndex, kIndex, &
                                     stepIndex, procIndex, numGridpt)
    !
    !:Purpose: To determine the grid points and their associated weights
    !          for the bilinear horizontal interpolation. If mask is present
    !          we currently can only handle a single 2D mask (like for sea
    !          ice or SST analysis). Will abort if multiple ocean levels present.
    !
    implicit none

    ! Arguments:
    type(struct_interpInfo), intent(inout) :: interpInfo
    type(struct_gsv)       , intent(in)    :: stateVector
    integer                , intent(in)    :: headerIndex, kIndex, stepIndex
    integer                , intent(in)    :: procIndex
    integer                , intent(out)   :: numGridpt(interpInfo%hco%numSubGrid)

    ! Locals:
    integer :: depotIndex
    integer :: ierr, niP1
    integer :: latIndex, lonIndex, latIndex2, lonIndex2, lonIndexP1
    integer :: subGridIndex, subGridForInterp, numSubGridsForInterp
    integer :: ipoint, gridptCount
    integer :: latIndexVec(4), lonIndexVec(4)
    logical :: mask(2,2)
    real(8) :: WeightVec(4)
    real(8) :: dldx, dldy
    real(8) :: weightsSum
    real(4) :: lon_deg_r4, lat_deg_r4
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    integer, parameter :: leftIndex = 1, rightIndex = 2, bottomIndex = 1, topIndex = 2

    numGridpt(:) = 0

    lat_deg_r4 = real(interpInfo%stepProcData(procIndex, stepIndex)%allLat(headerIndex, kIndex) *  &
                 MPC_DEGREES_PER_RADIAN_R8)
    lon_deg_r4 = real(interpInfo%stepProcData(procIndex, stepIndex)%allLon(headerIndex, kIndex) *  &
                 MPC_DEGREES_PER_RADIAN_R8)
    ierr = gpos_getPositionXY( stateVector%hco%EZscintID,   &
                              xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                              lat_deg_r4, lon_deg_r4, subGridIndex )

    ! Allow for periodicity in Longitude for global Gaussian grid
    if ( stateVector%hco%grtyp == 'G' .or. &
         (stateVector%hco%grtyp == 'Z' .and. stateVector%hco%global) ) then
      niP1 = statevector%ni + 1
    else
      niP1 = statevector%ni
    end if

    ! Find the lower-left grid point next to the observation
    if ( xpos_r4 >= real(niP1) ) then
      xpos_r4 = real(niP1)
      lonIndex = niP1 - 1
    else if ( xpos_r4 < 1.0 ) then
      xpos_r4 = 1.0
      lonIndex = 1
    else
      lonIndex = floor(xpos_r4)
    end if
    if ( xpos2_r4 >= real(niP1) ) then
      xpos2_r4 = real(niP1)
      lonIndex2 = niP1 - 1
    else if ( xpos2_r4 < 1.0 ) then
      xpos2_r4 = 1.0
      lonIndex2 = 1
    else
      lonIndex2 = floor(xpos2_r4)
    end if

    if ( ypos_r4 >= real(statevector%nj) ) then
      ypos_r4 = real(statevector%nj)
      latIndex = statevector%nj - 1
    else if ( ypos_r4 < 1.0 ) then
      ypos_r4 = 1.0
      latIndex = 1
    else
      latIndex = floor(ypos_r4)
    end if
    if ( ypos2_r4 >= real(statevector%nj) ) then
      ypos2_r4 = real(statevector%nj)
      latIndex2 = statevector%nj - 1
    else if ( ypos2_r4 < 1.0 ) then
      ypos2_r4 = 1.0
      latIndex2 = 1
    else
      latIndex2 = floor(ypos2_r4)
    end if

    if ( stateVector%hco%grtyp == 'U' ) then
      if ( ypos_r4 == real(stateVector%nj/2) ) then
        latIndex = floor(ypos_r4) - 1
      end if
      if ( ypos2_r4 == real(stateVector%nj/2) ) then
        latIndex2 = floor(ypos2_r4) - 1
      end if
    end if

    ! Handle periodicity in longitude
    lonIndexP1 = lonIndex + 1
    if ( lonIndexP1 == statevector%ni + 1 ) lonIndexP1 = 1

    ! Check if location is in between Yin and Yang (should not happen)
    if ( stateVector%hco%grtyp == 'U' ) then
      if ( ypos_r4 > real(stateVector%nj/2) .and.  &
           ypos_r4 < real((stateVector%nj/2)+1) ) then
        write(*,*) 's2c_setupBilinearInterp: WARNING, obs position in between Yin and Yang!!!'
        write(*,*) '   xpos, ypos = ', xpos_r4, ypos_r4
      end if
      if ( ypos2_r4 > real(stateVector%nj/2) .and.  &
           ypos2_r4 < real((stateVector%nj/2)+1) ) then
        write(*,*) 's2c_setupBilinearInterp: WARNING, obs position in between Yin and Yang!!!'
        write(*,*) '   xpos2, ypos2 = ', xpos2_r4, ypos2_r4
      end if
    end if

    if ( subGridIndex == 3 ) then
      ! both subGrids involved in interpolation, so first treat subGrid 1
      numSubGridsForInterp = 2
      subGridIndex = 1
    else
      ! only 1 subGrid involved in interpolation
      numSubGridsForInterp = 1
    end if

    if ( stateVector%oceanMask%maskPresent ) then
      ! abort if 3D mask is present, since we may not handle this situation correctly
      if ( stateVector%oceanMask%nLev > 1 ) then
        call utl_abort('s2c_setupBilinearInterp: 3D mask present - this case not properly handled')
      end if
      ! extract the ocean mask
      mask(leftIndex ,bottomIndex) = stateVector%oceanMask%mask(lonIndex  ,latIndex    ,1)
      mask(rightIndex,bottomIndex) = stateVector%oceanMask%mask(lonIndexP1,latIndex    ,1)
      mask(leftIndex ,topIndex   ) = stateVector%oceanMask%mask(lonIndex  ,latIndex + 1,1)
      mask(rightIndex,topIndex   ) = stateVector%oceanMask%mask(lonIndexP1,latIndex + 1,1)
    else
      mask(:,:) = .true.
    end if

    do subGridForInterp = 1, numSubGridsForInterp

      WeightVec(:) = 0
      gridptCount = 0

      ! Compute the 4 weights of the bilinear interpolation
      if ( subGridForInterp == 1 ) then
        ! when only 1 subGrid involved, subGridIndex can be 1 or 2
        dldx = real(xpos_r4,8) - real(lonIndex,8)
        dldy = real(ypos_r4,8) - real(latIndex,8)
      else
        ! when 2 subGrids, subGridIndex is set to 1 for 1st iteration, 2 for second
        subGridIndex = 2
        lonIndex = lonIndex2
        latIndex = latIndex2
        lonIndexP1 = lonIndex2 + 1
        dldx = real(xpos2_r4,8) - real(lonIndex,8)
        dldy = real(ypos2_r4,8) - real(latIndex,8)
      end if

      if ( mask(leftIndex ,bottomIndex) ) then
        gridptCount = gridptCount + 1
        latIndexVec(gridptCount) = latIndex
        lonIndexVec(gridptCount) = lonIndex
        WeightVec(gridptCount) = (1.d0-dldx) * (1.d0-dldy)
      end if

      if ( mask(rightIndex,bottomIndex) ) then
        gridptCount = gridptCount + 1
        latIndexVec(gridptCount) = latIndex
        lonIndexVec(gridptCount) = lonIndexP1
        WeightVec(gridptCount) =       dldx  * (1.d0-dldy)
      end if

      if ( mask(leftIndex ,topIndex   ) ) then
        gridptCount = gridptCount + 1
        latIndexVec(gridptCount) = latIndex + 1
        lonIndexVec(gridptCount) = lonIndex
        WeightVec(gridptCount) = (1.d0-dldx) *       dldy
      end if

      if ( mask(rightIndex,topIndex   ) ) then
        gridptCount = gridptCount + 1
        latIndexVec(gridptCount) = latIndex + 1
        lonIndexVec(gridptCount) = lonIndexP1
        WeightVec(gridptCount) =       dldx  *       dldy
      end if

      weightsSum = sum(WeightVec(1:gridptCount))
      if ( weightsSum > 0.d0 ) then
        WeightVec(1:gridptCount) = WeightVec(1:gridptCount) / weightsSum
      end if

      ! divide weight by number of subGrids
      WeightVec(1:gridptCount) = WeightVec(1:gridptCount) / real(numSubGridsForInterp,8)

      if ( allocated(interpInfo%interpWeightDepot) ) then

        depotIndex = interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, kIndex)

        do ipoint=1,gridptCount

          interpInfo%interpWeightDepot(depotIndex) = WeightVec(ipoint)
          interpInfo%latIndexDepot(depotIndex)     = latIndexVec(ipoint)
          interpInfo%lonIndexDepot(depotIndex)     = lonIndexVec(ipoint)
          depotIndex = depotIndex + 1

        end do

      end if

      numGridpt(subGridIndex) = gridptCount

    end do ! subGrid

  end subroutine s2c_setupBilinearInterp

  !--------------------------------------------------------------------------
  ! s2c_setupFootprintInterp
  !--------------------------------------------------------------------------
  subroutine s2c_setupFootprintInterp(fpr, interpInfo, stateVector, headerIndex, &
                                      kIndex, stepIndex, procIndex, numGridpt)
    !
    !:Purpose: To determine the grid points and their associated weights
    !          for the footprint horizontal interpolation.
    !
    implicit none

    ! Arguments:
    real(4)                , intent(in)    :: fpr ! footprint radius (metres)
    type(struct_interpInfo), intent(inout) :: interpInfo
    type(struct_gsv)       , intent(in)    :: stateVector
    integer                , intent(in)    :: headerIndex, kIndex, stepIndex
    integer                , intent(in)    :: procIndex
    integer                , intent(out)   :: numGridpt(interpInfo%hco%numSubGrid)

    ! Locals:
    integer :: depotIndex
    integer :: ierr
    integer :: latIndexCentre, lonIndexCentre, latIndexCentre2, lonIndexCentre2
    integer :: subGridIndex, numLocalGridptsFoundSearch
    real(4) :: lonObs_deg_r4, latObs_deg_r4
    real(8) :: lonObs, latObs
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    integer :: ipoint, gridptCount
    integer :: lonIndex, latIndex, resultsIndex, gridIndex
    integer :: lonIndexVec(maxNumLocalGridptsSearch), latIndexVec(maxNumLocalGridptsSearch)

    type(kdtree2_result)      :: searchResults(maxNumLocalGridptsSearch)
    real(kdkind)              :: refPosition(3), maxRadiusSquared
    type(kdtree2), pointer    :: tree

    numGridpt(:) = 0

    ! Determine the grid point nearest the observation.

    latObs = interpInfo % stepProcData(procIndex, stepIndex) % allLat(headerIndex, kIndex)
    lonObs = interpInfo % stepProcData(procIndex, stepIndex) % allLon(headerIndex, kIndex)

    latObs_deg_r4 = real(latObs * MPC_DEGREES_PER_RADIAN_R8)
    lonObs_deg_r4 = real(lonObs * MPC_DEGREES_PER_RADIAN_R8)
    ierr = gpos_getPositionXY( stateVector%hco%EZscintID,   &
                              xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                              latObs_deg_r4, lonObs_deg_r4, subGridIndex )

    lonIndexCentre = nint(xpos_r4)
    latIndexCentre = nint(ypos_r4)
    lonIndexCentre2 = nint(xpos2_r4)
    latIndexCentre2 = nint(ypos2_r4)

    if ( subGridIndex == 3 ) then
      write(*,*) 's2c_setupFootprintInterp: revise code'
      call utl_abort('s2c_setupFootprintInterp: both subGrids involved in interpolation.')
    end if

    ! Return if observation is not on the grid, or masked.
    if ( lonIndexCentre < 1 .or. lonIndexCentre > statevector%hco%ni .or.  &
         latIndexCentre < 1 .or. latIndexCentre > statevector%hco%nj ) return

    if ( stateVector%oceanMask%maskPresent ) then
      ! abort if 3D mask is present, since we may not handle this situation correctly
      if ( stateVector%oceanMask%nLev > 1 ) then
        call utl_abort('s2c_setupFootprintInterp: 3D mask present - this case not properly handled')
      end if

      if ( .not. stateVector%oceanMask%mask(lonIndexCentre,latIndexCentre,1) ) return
    end if

    ! do the search
    maxRadiusSquared = real(fpr,8) ** 2
    refPosition(:) = kdtree2_3dPosition(lonObs, latObs)
    nullify(tree)
    if ( interpInfo%inputStateVectorType == 'nl' ) then
      if ( associated(tree_nl) ) then
        tree => tree_nl
      else
        call utl_abort('s2c_setupFootprintInterp: tree_nl is not allocated!')
      end if
    else if ( interpInfo%inputStateVectorType == 'tl' .or. &
              interpInfo%inputStateVectorType == 'ad' ) then
      if ( associated(tree_tlad) ) then
        tree => tree_tlad
      else
        call utl_abort('s2c_setupFootprintInterp: tree_tlad is not allocated!')
      end if
    end if
    call kdtree2_r_nearest(tp=tree, qv=refPosition, r2=maxRadiusSquared, &
                           nfound=numLocalGridptsFoundSearch, &
                           nalloc=maxNumLocalGridptsSearch, &
                           results=searchResults)

    if (numLocalGridptsFoundSearch > maxNumLocalGridptsSearch ) then
      call utl_abort('s2c_setupFootprintInterp: the parameter maxNumLocalGridptsSearch must be increased')
    else if ( numLocalGridptsFoundSearch < minNumLocalGridptsSearch .and. useFootprintForTovs ) then
      write(*,*) 's2c_setupFootprintInterp: Warning! For TOVS headerIndex=', headerIndex, &
                 ' number of grid points found within footprint radius=', fpr, ' is less than ', &
                 minNumLocalGridptsSearch 
    end if

    ! ensure at least the nearest neighbor is included in lonIndexVec/latIndexVec
    ! if footprint size is smaller than the grid spacing.
    gridptCount = 1
    lonIndexVec(gridptCount) = lonIndexCentre
    latIndexVec(gridptCount) = latIndexCentre

    ! fill the rest of lonIndexVec/latIndexVec
    gridLoop1: do resultsIndex = 1, numLocalGridptsFoundSearch
      gridIndex = searchResults(resultsIndex)%idx
      if ( gridIndex < 1 .or. gridIndex > statevector%hco%ni * statevector%hco%nj ) then
        write(*,*) 's2c_setupFootprintInterp: gridIndex=', gridIndex
        call utl_abort('s2c_setupFootprintInterp: gridIndex out of bound.')
      end if

      latIndex = (gridIndex - 1) / statevector%hco%ni + 1
      lonIndex = gridIndex - (latIndex - 1) * statevector%hco%ni
      if ( lonIndex < 1 .or. lonIndex > statevector%hco%ni .or. &
           latIndex < 1 .or. latIndex > statevector%hco%nj ) then
        write(*,*) 's2c_setupFootprintInterp: lonIndex=', lonIndex, ',latIndex=', latIndex
        call utl_abort('s2c_setupFootprintInterp: lonIndex/latIndex out of bound.')
      end if

      if ( stateVector%oceanMask%maskPresent ) then
        if ( .not. stateVector%oceanMask%mask(lonIndex,latIndex,1) ) cycle gridLoop1
      end if

      if ( lonIndex == lonIndexCentre .and. latIndex == latIndexCentre ) cycle gridLoop1

      gridptCount = gridptCount + 1
      lonIndexVec(gridptCount) = lonIndex
      latIndexVec(gridptCount) = latIndex
    end do gridLoop1

    if ( allocated(interpInfo%interpWeightDepot) ) then

      depotIndex = interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, kIndex)

      do ipoint = 1, gridptCount

        interpInfo%interpWeightDepot(depotIndex) = 1.0d0 / real(gridptCount,8)
        interpInfo%latIndexDepot(depotIndex)     = latIndexVec(ipoint)
        interpInfo%lonIndexDepot(depotIndex)     = lonIndexVec(ipoint)
        depotIndex = depotIndex + 1

      end do

    end if

    numGridpt(subGridIndex) = gridptCount

  end subroutine s2c_setupFootprintInterp

  !--------------------------------------------------------------------------
  ! s2c_setupLakeInterp
  !--------------------------------------------------------------------------
  subroutine s2c_setupLakeInterp(interpInfo, stateVector, headerIndex, kIndex, &
                                 stepIndex, procIndex, numGridpt)
    !
    !:Purpose: To determine the grid points and their associated weights
    !          for the lake horizontal interpolation.
    !
    implicit none

    ! Arguments:
    type(struct_interpInfo), intent(inout) :: interpInfo
    type(struct_gsv)       , intent(in)    :: stateVector
    integer                , intent(in)    :: headerIndex, kIndex, stepIndex
    integer                , intent(in)    :: procIndex
    integer                , intent(out)   :: numGridpt(interpInfo%hco%numSubGrid)

    ! Locals:
    integer :: depotIndex
    integer :: ierr
    integer :: latIndexCentre, lonIndexCentre, latIndexCentre2, lonIndexCentre2
    integer :: subGridIndex, subGridForInterp, numSubGridsForInterp
    real(4) :: lon_deg_r4, lat_deg_r4
    real(8) :: lon_rad, lat_rad
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    integer :: ipoint, gridptCount
    integer :: lakeCount
    integer :: lonIndex, latIndex, lakeIndex
    integer :: lonIndexVec(statevector%ni*statevector%nj), latIndexVec(statevector%ni*statevector%nj)
    logical :: reject, lake(statevector%ni,statevector%nj)
    integer :: k, l

    if ( stateVector%hco%grtyp == 'U' ) then
      call utl_abort('s2c_setupLakeInterp: Yin-Yang grid not supported')
    end if

    if ( .not.stateVector%oceanMask%maskPresent ) then
      call utl_abort('s2c_setupLakeInterp: Only compatible when mask present')
    end if

    numGridpt(:) = 0

    reject = .false.

    numGridpt(:) = 0

    ! Determine the grid point nearest the observation.

    lat_rad = interpInfo%stepProcData(procIndex, stepIndex)%allLat(headerIndex, kIndex)
    lon_rad = interpInfo%stepProcData(procIndex, stepIndex)%allLon(headerIndex, kIndex)
    lat_deg_r4 = real(lat_rad * MPC_DEGREES_PER_RADIAN_R8)
    lon_deg_r4 = real(lon_rad * MPC_DEGREES_PER_RADIAN_R8)
    ierr = gpos_getPositionXY( stateVector%hco%EZscintID,   &
                               xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                               lat_deg_r4, lon_deg_r4, subGridIndex )

    lonIndexCentre = nint(xpos_r4)
    latIndexCentre = nint(ypos_r4)
    lonIndexCentre2 = nint(xpos2_r4)
    latIndexCentre2 = nint(ypos2_r4)

    if ( subGridIndex == 3 ) then
      write(*,*) 's2c_setupLakeInterp: revise code'
      call utl_abort('s2c_setupLakeInterp: both subGrids involved in interpolation.')
      numSubGridsForInterp = 2
      subGridIndex = 1
    else
      ! only 1 subGrid involved in interpolation
      numSubGridsForInterp = 1
    end if

    do subGridForInterp = 1, numSubGridsForInterp

      gridptCount = 0

      ! It can happen that the lake location is closest to a grid point
      ! where MASK(I,J) = .false. while there are other grid points for the
      ! same lake where MASK(I,J) = .true.. Code needs modifications
      ! for this case.

      ! If observation is not on the grid, don't use it.
      if ( lonIndexCentre < 1 .or. lonIndexCentre > statevector%ni .or.  &
           latIndexCentre < 1 .or. latIndexCentre > statevector%nj ) reject = .true.

      if ( .not. stateVector%oceanMask%mask(lonIndexCentre,latIndexCentre,1) ) reject = .true.

      if ( .not. reject ) then

        lake(:,:) = .false.
        lake(lonIndexCentre,latIndexCentre) = .true.
        gridptCount = 1
        lonIndexVec(gridptCount) = lonIndexCentre
        latIndexVec(gridptCount) = latIndexCentre

        lakeCount = 0

        do while(lakeCount /= gridptCount)

          do lakeIndex = lakeCount+1, gridptCount

            if(lakeIndex == lakeCount+1) lakeCount = gridptCount

            k = lonIndexVec(lakeIndex)
            l = latIndexVec(lakeIndex)

            do latIndex = max(1,l-1),min(l+1,statevector%nj)
              do lonIndex = max(1,k-1),min(k+1,statevector%ni)
                if(stateVector%oceanMask%mask(lonIndex,latIndex,1) .and. .not. lake(lonIndex,latIndex)) then
                  lake(lonIndex,latIndex) = .true.
                  gridptCount = gridptCount + 1
                  lonIndexVec(gridptCount) = lonIndex
                  latIndexVec(gridptCount) = latIndex
                end if
              end do
            end do

          end do

        end do

        if ( allocated(interpInfo%interpWeightDepot) ) then

          depotIndex = interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, kIndex)

          do ipoint=1,gridptCount

            interpInfo%interpWeightDepot(depotIndex) = 1.0d0 / real(gridptCount,8)
            interpInfo%latIndexDepot(depotIndex)     = latIndexVec(ipoint)
            interpInfo%lonIndexDepot(depotIndex)     = lonIndexVec(ipoint)
            depotIndex = depotIndex + 1

          end do

        end if

        numGridpt(subGridIndex) = gridptCount

      end if ! not reject

    end do ! subGrid

  end subroutine s2c_setupLakeInterp

  !--------------------------------------------------------------------------
  ! s2c_setupNearestNeighbor
  !--------------------------------------------------------------------------
  subroutine s2c_setupNearestNeighbor(interpInfo, stateVector, headerIndex, kIndex, &
                                      stepIndex, procIndex, numGridpt)
    !
    !:Purpose: Determine the nearest grid points to the observations location
    !
    implicit none

    ! arguments
    type(struct_interpInfo), intent(inout) :: interpInfo
    type(struct_gsv)       , intent(in)    :: stateVector
    integer                , intent(in)    :: headerIndex, kIndex, stepIndex, procIndex
    integer                , intent(out)   :: numGridpt(interpInfo%hco%numSubGrid)

    ! locals
    integer :: depotIndex
    integer :: ierr
    integer :: latIndex, lonIndex
    integer :: subGridIndex
    real(4) :: lon_deg_r4, lat_deg_r4
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4

    if ( stateVector%hco%grtyp == 'U' ) then
      call utl_abort('s2c_setupNearestNeighbor: Yin-Yang grid not supported')
    end if

    numGridpt(:) = 0

    lat_deg_r4 = real(interpInfo%stepProcData(procIndex, stepIndex)%allLat(headerIndex, kIndex) *  &
                 MPC_DEGREES_PER_RADIAN_R8)
    lon_deg_r4 = real(interpInfo%stepProcData(procIndex, stepIndex)%allLon(headerIndex, kIndex) *  &
                 MPC_DEGREES_PER_RADIAN_R8)

    ierr = gpos_getPositionXY( stateVector%hco%EZscintID,   &
                              xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                              lat_deg_r4, lon_deg_r4, subGridIndex )

    latIndex = nint(ypos_r4)
    lonIndex = nint(xpos_r4)

    ! Handle periodicity in longitude
    if ( lonIndex == statevector%ni+1 .and. stateVector%hco%grtyp == 'G' ) lonIndex = 1

    ! Test bounds
    if ( lonIndex < 1 .or. lonIndex > statevector%ni .or. &
         latIndex < 1 .or. latIndex > statevector%nj  ) then

      write(*,*) 's2c_setupNearestNeighbor: observation out of bounds'

    else

      if ( allocated(interpInfo%interpWeightDepot) ) then
      
        depotIndex = interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, kIndex)

        interpInfo%interpWeightDepot(depotIndex) = 1.d0
        interpInfo%latIndexDepot    (depotIndex) = latIndex
        interpInfo%lonIndexDepot    (depotIndex) = lonIndex

      end if

      numGridpt(subGridIndex) = 1

    end if

  end subroutine s2c_setupNearestNeighbor

  !--------------------------------------------------------------------------
  ! checkColumnStatevectorMatch
  !--------------------------------------------------------------------------
  subroutine checkColumnStatevectorMatch(column,statevector)
    !
    !:Purpose: To check column and statevector have identical nk and variables.
    !
    implicit none

    ! Arguments:
    type(struct_gsv)       , intent(in) :: statevector
    type(struct_columnData), intent(in) :: column

    ! Locals:
    integer :: kIndex

    ! check column/statevector have same nk
    if ( column%nk /= gsv_getNumK(statevector) ) then
      write(*,*) 'checkColumnStatevectorMatch: column%nk, gsv_getNumK(statevector)', column%nk, gsv_getNumK(statevector)
      call utl_abort('checkColumnStatevectorMatch: column%nk /= gsv_getNumK(statevector)')
    end if
    
    ! loop through k and check varNames are same between column/statevector
    do kIndex = 1, column%nk
      if (gsv_getVarNameFromK(statevector,kIndex) /= col_getVarNameFromK(column,kIndex)) then
        write(*,*) 'checkColumnStatevectorMatch: kIndex, varname in statevector and column: ', kIndex, &
                   gsv_getVarNameFromK(statevector,kIndex), col_getVarNameFromK(column,kIndex) 
        call utl_abort('checkColumnStatevectorMatch: varname in column and statevector do not match')
      end if	
    end do

  end subroutine checkColumnStatevectorMatch

  !--------------------------------------------------------------------------
  ! latlonChecks
  !--------------------------------------------------------------------------
  subroutine latlonChecks( obsSpaceData, hco, headerIndex, rejectOutsideObs, &
    latLev_T, lonLev_T, latLev_M, lonLev_M, latLev_S, lonLev_S )
    !
    !:Purpose: To check if the obs are inside the domain.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout)  :: obsSpaceData
    type(struct_hco), intent(in)  :: hco
    integer,          intent(in)  :: headerIndex
    logical,          intent(in)  :: rejectOutsideObs
    real(8),          intent(inout)  :: latLev_T(:)
    real(8),          intent(inout)  :: lonLev_T(:)
    real(8),          intent(inout)  :: latLev_M(:)
    real(8),          intent(inout)  :: lonLev_M(:)
    real(8), intent(inout),optional  :: latLev_S
    real(8), intent(inout),optional  :: lonLev_S

    ! Locals:
    integer :: ierr
    integer :: bodyIndex, bodyIndexBeg, bodyIndexEnd, niP1, subGridIndex
    integer :: nlev_T, nlev_M
    real(4) :: lon_r4, lat_r4, lon_deg_r4, lat_deg_r4
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    logical :: latlonOutsideGrid, rejectHeader

    ! external functions
    integer :: gdllfxy

    ! Allow for periodicity in Longitude for global Gaussian grid
    if ( hco%grtyp == 'G' .or. (hco%grtyp == 'Z' .and. hco%global) ) then
      niP1 = hco%ni + 1
    else
      niP1 = hco%ni
    end if

    nlev_T = size(latLev_T)
    nlev_M = size(latLev_M)

    ! check if lat/lon of last thermo level is outside domain.
    rejectHeader = .false.
    lat_r4 = real(latLev_T(nlev_T),4)
    lon_r4 = real(lonLev_T(nlev_T),4)

    lat_deg_r4 = lat_r4 * MPC_DEGREES_PER_RADIAN_R8
    lon_deg_r4 = lon_r4 * MPC_DEGREES_PER_RADIAN_R8
    ierr = gpos_getPositionXY( hco%EZscintID,   &
                              xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                              lat_deg_r4, lon_deg_r4, subGridIndex )

    latlonOutsideGrid = ( xpos_r4 < 1.0        .or. &
                          xpos_r4 > real(niP1) .or. &
                          ypos_r4 < 1.0        .or. &
                          ypos_r4 > real(hco%nj) )

    if ( latlonOutsideGrid .and. rejectOutsideObs ) then
      rejectHeader = .true.
    end if

    !  check if lat/lon of last momentum level is outside domain.
    if ( .not. rejectHeader ) then
      lat_r4 = real(latLev_M(nlev_M),4)
      lon_r4 = real(lonLev_M(nlev_M),4)

      lat_deg_r4 = lat_r4 * MPC_DEGREES_PER_RADIAN_R8
      lon_deg_r4 = lon_r4 * MPC_DEGREES_PER_RADIAN_R8
      ierr = gpos_getPositionXY( hco%EZscintID,   &
                                xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                                lat_deg_r4, lon_deg_r4, subGridIndex )

      latlonOutsideGrid = ( xpos_r4 < 1.0        .or. &
                            xpos_r4 > real(niP1) .or. &
                            ypos_r4 < 1.0        .or. &
                            ypos_r4 > real(hco%nj) )

      if ( latlonOutsideGrid .and. rejectOutsideObs ) then
        rejectHeader = .true.
      end if
    end if

    !  check if lat/lon of surface level is outside domain.
    if ( present(latLev_S) .and. present(lonLev_S) .and. .not. rejectHeader ) then
      lat_r4 = real(latLev_S,4)
      lon_r4 = real(lonLev_S,4)

      lat_deg_r4 = lat_r4 * MPC_DEGREES_PER_RADIAN_R8
      lon_deg_r4 = lon_r4 * MPC_DEGREES_PER_RADIAN_R8
      ierr = gpos_getPositionXY( hco%EZscintID,   &
                                xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                                lat_deg_r4, lon_deg_r4, subGridIndex )

      latlonOutsideGrid = ( xpos_r4 < 1.0        .or. &
                            xpos_r4 > real(niP1) .or. &
                            ypos_r4 < 1.0        .or. &
                            ypos_r4 > real(hco%nj) )

      if ( latlonOutsideGrid .and. rejectOutsideObs ) then
        rejectHeader = .true.
      end if
    end if

    if ( rejectHeader ) then
      ! The observation is outside the domain.
      ! With a LAM trial field we must discard this observation
      write(*,*) 'latlonChecks: Rejecting OBS outside the hco domain, ', headerIndex
      write(*,*) '  position : ', lat_deg_r4, lon_deg_r4, ypos_r4, xpos_r4

      bodyIndexBeg = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
      bodyIndexEnd = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) + bodyIndexBeg -1
      do bodyIndex = bodyIndexBeg, bodyIndexEnd
        call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
      end do
      call obs_headSet_i(obsSpaceData, OBS_ST1, headerIndex,  &
           ibset( obs_headElem_i(obsSpaceData, OBS_ST1, headerIndex), 05))

      ! Assign domain mid-point lat-lon to this header
      if ( hco%grtyp == 'Y' ) then
        lat_deg_r4 = hco%lat2d_4(hco%ni/2,hco%nj/2)
        lon_deg_r4 = hco%lon2d_4(hco%ni/2,hco%nj/2)
      else
        xpos_r4 = real(hco%ni)/2.0
        ypos_r4 = real(hco%nj)/2.0
        ierr = gdllfxy(hco%EZscintID, lat_deg_r4, lon_deg_r4, &
                       xpos_r4, ypos_r4, 1)
      end if

      lonLev_T(:) = real(lon_deg_r4 * MPC_RADIANS_PER_DEGREE_R4,8)
      latLev_T(:) = real(lat_deg_r4 * MPC_RADIANS_PER_DEGREE_R4,8)
      lonLev_M(:) = real(lon_deg_r4 * MPC_RADIANS_PER_DEGREE_R4,8)
      latLev_M(:) = real(lat_deg_r4 * MPC_RADIANS_PER_DEGREE_R4,8)
      if (present(lonLev_S) .and. present(latLev_S)) then
        lonLev_S    = real(lon_deg_r4 * MPC_RADIANS_PER_DEGREE_R4,8)
        latLev_S    = real(lat_deg_r4 * MPC_RADIANS_PER_DEGREE_R4,8)
      end if

    end if

  end subroutine latlonChecks

  !--------------------------------------------------------------------------
  ! getTovsFootprintRadius
  !--------------------------------------------------------------------------
  function getTovsFootprintRadius(obsSpaceData, headerIndex, beSilent_opt) result(footPrintRadius_r4)
    !
    !:Purpose: calculate foot-print radius for TOVS observations
    !
    implicit none

    ! Arguments
    type(struct_obs), intent(in)  :: obsSpaceData
    real(4)                       :: footPrintRadius_r4
    integer         , intent(in)  :: headerIndex
    logical         , intent(in), optional :: beSilent_opt
    
    ! local
    integer :: codtyp, sensorIndex 
    real(8) :: fovAngularDiameter, satHeight, footPrintRadius
    character(len=codtyp_name_length) :: instrumName
    logical :: beSilent

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .true.
    end if

    ! get nominal satellite height
    sensorIndex = tvs_lsensor(tvs_tovsIndex(headerIndex))
    satHeight = tvs_coefs(sensorIndex)%coef%fc_sat_height

    ! FOV angular diameter  
    codtyp = obs_headElem_i( obsSpaceData, OBS_ITY, headerIndex )
    instrumName = codtyp_get_name(codtyp)
    select case(trim(instrumName))
    case('amsua')
      fovAngularDiameter = 3.3d0
    case('amsub')
      fovAngularDiameter = 1.1d0
    case('mhs')
      fovAngularDiameter = 10.0d0 / 9.0d0
    case('airs')
      fovAngularDiameter = 1.1d0
    case('iasi')
      fovAngularDiameter = 14.65d0 / 1000.0d0 * MPC_DEGREES_PER_RADIAN_R8
    case('radianceclear')
      fovAngularDiameter = 0.125d0
    case('ssmis')
      fovAngularDiameter = 1.2d0
    case('atms')
      fovAngularDiameter = 2.2d0
    case('cris')
      fovAngularDiameter = 14.0d0 / 824.0d0 * MPC_DEGREES_PER_RADIAN_R8
    case default
      fovAngularDiameter = -1.0d0
    end select

    if ( fovAngularDiameter < 0.0d0 ) then 
      footPrintRadius_r4 = bilinearFootprint
    else
      ! get foot print radius (meter) from angular diameter
      footPrintRadius = 0.5d0 * fovAngularDiameter * MPC_RADIANS_PER_DEGREE_R8 * satHeight * 1000
      footPrintRadius_r4 = real(footPrintRadius,4)
    end if

    if ( .not. beSilent ) then
      write(*,*) 'getTovsFootprintRadius: sensorIndex=', sensorIndex, &
                ',satHeight=', satHeight, ',fovAngularDiameter=', fovAngularDiameter, ',codtyp=', codtyp, &
                ',footPrintRadius=', footPrintRadius_r4
    end if

  end function getTovsFootprintRadius

  ! -------------------------------------------------------------
  ! s2c_getWeightsAndGridPointIndexes
  ! -------------------------------------------------------------
  subroutine s2c_getWeightsAndGridPointIndexes(headerIndex, kIndex, stepIndex, procIndex, &
       interpWeight, latIndex, lonIndex, gridptCount)
    ! :Purpose: Returns the weights and grid point indexes for a single observation.
    !           
    !
    implicit none

    ! arguments
    integer, intent(in)  :: headerIndex
    integer, intent(in)  :: kIndex
    integer, intent(in)  :: stepIndex
    integer, intent(in)  :: procIndex
    real(8), intent(out) :: interpWeight(:)
    integer, intent(out) :: latIndex(:)
    integer, intent(out) :: lonIndex(:)
    integer, intent(out) :: gridptCount

    ! locals
    integer :: indexBeg, indexEnd, gridptIndex
    integer :: subGridIndex, maxGridpt

    call utl_tmg_start(30,'--StateToColumn')

    maxGridpt = size( interpWeight )

    gridptCount = 0

    if ( interpInfo_tlad%stepProcData(procIndex, stepIndex)%allHeaderIndex(headerIndex) /= headerIndex ) then
      call utl_abort('s2c_getWeightsAndGridPointIndexes: headerUsedIndex and headerIndex differ.'//    &
                     ' If using multiple time steps in the assimilation window,'//                     &
                     ' the code needs to be modified to convert values of headerIndex into headerUsedIndex.')
    end if

    subGrid_loop: do subGridIndex = 1, interpInfo_tlad%hco%numSubGrid

      indexBeg = interpInfo_tlad%stepProcData(procIndex,stepIndex)%depotIndexBeg(subGridIndex, headerIndex, kIndex)
      indexEnd = interpInfo_tlad%stepProcData(procIndex,stepIndex)%depotIndexEnd(subGridIndex, headerIndex, kIndex)

      if ( indexEnd < IndexBeg ) cycle subGrid_loop

      do gridptIndex = indexBeg, indexEnd

        gridptCount = gridptCount + 1

        if ( gridptCount > maxGridpt ) then
          call utl_abort('s2c_getWeightsAndGridPointIndexes: maxGridpt must be increased')
        end if

        lonIndex(gridptCount) = interpInfo_tlad%lonIndexDepot(gridptIndex)
        latIndex(gridptCount) = interpInfo_tlad%latIndexDepot(gridptIndex)
        interpWeight(gridptCount) = interpInfo_tlad%interpWeightDepot(gridptIndex)

      end do

    end do subGrid_loop

    call utl_tmg_stop(30)

  end subroutine s2c_getWeightsAndGridPointIndexes

  ! -------------------------------------------------------------
  ! s2c_deallocInterpInfo
  ! -------------------------------------------------------------
  subroutine s2c_deallocInterpInfo( inputStateVectorType )
    ! :Purpose: Deallocate interpInfo_nl/tlad object.
    !
    implicit none

    ! arguments
    character(len=*), intent(in) :: inputStateVectorType

    ! locals
    type(struct_interpInfo), pointer :: interpInfo
    integer :: stepIndex, procIndex, numStep

    select case( trim(inputStateVectorType) )
      case('nl')
        interpInfo => interpInfo_nl
      case('tlad')
        interpInfo => interpInfo_tlad
      case default
        call utl_abort('s2c_deallocInterpInfo: invalid input argument' // inputStateVectorType)
    end select

    if ( .not. interpInfo%initialized ) return

    write(*,*) 's2c_deallocInterpInfo: deallocating interpInfo for inputStateVectorType=', &
                inputStateVectorType

    numStep = size(interpInfo%stepProcData,2)

    deallocate(interpInfo%interpWeightDepot)
    deallocate(interpInfo%latIndexDepot)
    deallocate(interpInfo%lonIndexDepot)
    do stepIndex = 1, numStep
      do procIndex = 1, mmpi_nprocs
        deallocate(interpInfo%stepProcData(procIndex,stepIndex)%allLat)
        deallocate(interpInfo%stepProcData(procIndex,stepIndex)%allLon)
        deallocate(interpInfo%stepProcData(procIndex,stepIndex)%allHeaderIndex)
        deallocate(interpInfo%stepProcData(procIndex,stepIndex)%depotIndexBeg)
        deallocate(interpInfo%stepProcData(procIndex,stepIndex)%depotIndexEnd)
        if ( interpInfo%hco%rotated ) then
          deallocate(interpInfo%stepProcData(procIndex,stepIndex)%allLonRot)
          deallocate(interpInfo%stepProcData(procIndex,stepIndex)%allLatRot)
        end if
      end do
    end do
    deallocate(interpInfo%stepProcData)
    deallocate(interpInfo%allNumHeaderUsed)
    call oti_deallocate(interpInfo%oti)

    interpInfo%initialized = .false.

  end subroutine s2c_deallocInterpInfo

end module stateToColumn_mod

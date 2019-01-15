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
  !
  ! MODULE stateToColumn (prefix='s2c' category='3. High-level transformations')
  !
  ! **Purpose:** 
  ! Non-linear, tangent-linear and adjoint versions of
  ! horizontal-temporal interpolation between a gridStateVector object
  ! and a columnData object.
  !
  use mathPhysConstants_mod
  use mpi_mod
  use mpivar_mod 
  use gridstatevector_mod
  use obsSpaceData_mod
  use columnData_mod
  use analysisgrid_mod
  use horizontalCoord_mod
  use obsTimeInterp_mod
  use tt2phi_mod
  use windRotation_mod
  use utilities_mod
  use variabletransforms_mod
  use calcPressure_mod
  
  implicit none
  save
  private
  
  ! public routines
  public :: s2c_tl, s2c_ad, s2c_nl
  public :: s2c_column_hbilin, s2c_bgcheck_bilin

  ! private module variables and derived types

  type struct_interpInfo
    logical                   :: initialized = .false.
    type(struct_hco), pointer :: hco => null() ! horizontal grid object
    type(struct_uvr), pointer :: uvr => null() ! windRotation object
    type(struct_oti), pointer :: oti => null() ! obsTimeInterp object

    ! number of obs headers on each proc having a non-zero interp weight for each stepIndex (headerUsed)
    integer, pointer          :: allNumHeaderUsed(:,:) => null()    ! (step, proc)

    ! actual headerIndex, since the headerUsed is only for those obs with a non-zero interp weight
    integer, pointer          :: allHeaderIndex(:,:,:) => null() ! (headerUsed, step, proc)

    ! lat-lon location of observations to be interpolated (only needed to rotate winds)
    real(8), pointer          :: allLat(:,:,:) => null()         ! (headerUsed, step, proc)
    real(8), pointer          :: allLon(:,:,:) => null()         ! (headerUsed, step, proc)
    real(8), pointer          :: allLatRot(:,:,:,:) => null()    ! (subGrid, headerUsed, step, proc)
    real(8), pointer          :: allLonRot(:,:,:,:) => null()    ! (subGrid, headerUsed, step, proc)

    ! interpolation weights and lat/lon indices are accessed via the 'depotIndexBeg/End'
    integer, pointer          :: depotIndexBeg(:,:,:,:) => null()    ! (subGrid, headerUsed, step, proc)
    integer, pointer          :: depotIndexEnd(:,:,:,:) => null()    ! (subGrid, headerUsed, step, proc)
    real(8), pointer          :: interpWeightDepot(:)                ! (depotIndex)
    integer, pointer          :: latIndexDepot(:)                    ! (depotIndex)
    integer, pointer          :: lonIndexDepot(:)                    ! (depotIndex)
  end type struct_interpInfo

  type(struct_interpInfo) :: interpInfo_tlad, interpInfo_nl

  character(len=20), parameter :: timeInterpType_tlad = 'LINEAR' ! hardcoded type of time interpolation for increment

  integer, external    :: get_max_rss
  character(len=4) :: varsToInterpolate(5) = (/ 'ALL ','GZ_T','GZ_M','P_T ','P_M ' /)


contains 

  !---------------------------------------------------------
  ! s2c_latLonChecks
  !---------------------------------------------------------
  subroutine s2c_latLonChecks(obsSpaceData, moveObsAtPole)
    ! **Purpose:** 
    ! Check the lat/lon of observations and modify if necessary
    !
    implicit none

    ! arguments
    type(struct_obs) :: obsSpaceData
    logical          :: moveObsAtPole

    ! locals
    type(struct_hco), pointer :: hco_anl
    integer :: headerIndex, ierr
    integer :: idata, idatend, jdata, subGridIndex
    real(4) :: lat_r4, lon_r4, lat_deg_r4, lon_deg_r4
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    real(4) :: xposLowerBoundAnl_r4, xposUpperBoundAnl_r4
    real(8) :: lat_r8, lon_r8

    ! external functions
    integer :: gdllfxy

    write(*,*) ' '
    write(*,*) 's2c_latLonChecks: STARTING'
    write(*,*) ' '
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !-    Get the Analysis Grid structure
    !
    hco_anl => agd_getHco('CoreGrid')

    if ( hco_anl % global ) then
       xposLowerBoundAnl_r4 = - huge(1.0) ! no limit since grid is global (periodic)
       xposUpperBoundAnl_r4 = + huge(1.0) ! no limit since grid is global (periodic)
    else
       xposLowerBoundAnl_r4 = 1.0
       xposUpperBoundAnl_r4 = real(hco_anl % ni)
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
      ierr = getPositionXY( hco_anl % EZscintID,  &
                            xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                            lat_deg_r4, lon_deg_r4, subGridIndex )

      !- Test if the obs is outside the analysis grid
      if ( xpos_r4 < xposLowerBoundAnl_r4  .or. &
           xpos_r4 > xposUpperBoundAnl_r4  .or. &
           ypos_r4 < 1.0                   .or. &
           ypos_r4 > real(hco_anl % nj) ) then

        if ( hco_anl % global ) then

          if ( moveObsAtPole ) then
            ! Modify latitude if we have an observation at or near the poles
            write(*,*) ''
            write(*,*) 's2c_latLonChecks: Moving OBS inside the GLOBAL ANALYSIS grid, ', headerIndex
            write(*,*) '  true position : ', lat_deg_r4, lon_deg_r4, ypos_r4, xpos_r4

            !- Move the observation to the nearest grid point
            if ( ypos_r4 < 1.0 )                ypos_r4 = 1.0
            if ( ypos_r4 > real(hco_anl % nj) ) ypos_r4 = real(hco_anl % nj)

            ierr = gdllfxy( hco_anl % EZscintID, &    ! IN
                            lat_deg_r4, lon_deg_r4, & ! OUT
                            xpos_r4, ypos_r4, 1)      ! IN

            write(*,*) '  new  position : ', lat_deg_r4, lon_deg_r4, ypos_r4, xpos_r4

            lat_r8 = real(lat_deg_r4,8) * MPC_RADIANS_PER_DEGREE_R8
            lon_r8 = real(lon_deg_r4,8) * MPC_RADIANS_PER_DEGREE_R8
            call obs_headSet_r(obsSpaceData,OBS_LAT,headerIndex, lat_r8) ! IN
            call obs_headSet_r(obsSpaceData,OBS_LON,headerIndex, lon_r8) ! IN
          else
            write(*,*)
            write(*,*) 's2c_latLonChecks: OBS outside the GLOBAL ANALYSIS grid, but NOT moved, ', headerIndex
          end if

        else
          ! The observation is outside the domain
          ! In LAM Analysis mode we must discard this observation
          write(*,*) 's2c_latLonChecks: Rejecting OBS outside the LAM ANALYSIS grid domain, ', headerIndex
          write(*,*) '  position : ', lat_deg_r4, lon_deg_r4, ypos_r4, xpos_r4

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

    write(*,*) 's2c_latLonChecks: END'

  end subroutine s2c_latLonChecks

  !---------------------------------------------------------
  ! s2c_setupInterpInfo
  !---------------------------------------------------------
  subroutine s2c_setupInterpInfo( interpInfo, obsSpaceData, stateVector,  &
                                  timeInterpType, rejectOutsideObs )
    ! **Purpose:**
    ! Setup all of the information needed to quickly
    ! perform the horizontal interpolation to the observation
    ! locations.
    !
    implicit none

    ! arguments
    type(struct_interpInfo)    :: interpInfo
    type(struct_obs)           :: obsSpaceData
    type(struct_gsv)           :: stateVector
    logical                    :: rejectOutsideObs
    character(len=*)           :: timeInterpType

    ! locals
    integer :: numHeader, numHeaderUsedMax, headerIndex, bodyIndex
    integer :: numStep, stepIndex, ierr, indexBeg, indexEnd
    integer :: bodyIndexBeg, bodyIndexEnd, procIndex, niP1, numGridpt, numGridptTotal, numHeaderUsed
    integer :: latIndex, lonIndex, latIndex2, lonIndex2, lonIndexP1
    integer :: subGridIndex, subGridForInterp, numSubGridsForInterp
    real(8) :: dldx, dldy, xpos, ypos
    real(8) :: latRot, lonRot, lat, lon, weightsSum
    real(4) :: lon_r4, lat_r4, lon_deg_r4, lat_deg_r4
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    integer, allocatable :: allNumHeaderUsed(:,:), headerIndexVec(:,:)
    real(4), allocatable :: lonVec_r4(:), latVec_r4(:)
    integer :: ezgdef, gdllfxy
    logical :: obsOutsideGrid

    write(*,*) 's2c_setupInterpInfo: STARTING'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    numStep = stateVector%numStep
    numHeader = obs_numheader(obsSpaceData)

    call oti_setup(interpInfo%oti, obsSpaceData, numStep, timeInterpType, flagObsOutside_opt=.true.)

    ! Allow for periodicity in Longitude for global Gaussian grid
    if ( stateVector%hco%grtyp == 'G' ) then
      niP1 = statevector%ni + 1
    else
      niP1 = statevector%ni
    end if

    ! First count the number of headers for each stepIndex
    allocate(allNumHeaderUsed(numStep,mpi_nprocs))
    allNumHeaderUsed(:,:) = 0
    do stepIndex = 1, numStep
      numHeaderUsed = 0

      header_loop1: do headerIndex = 1, numHeader

        ! if obs inside window, but zero weight for current stepIndex then skip it
        if ( oti_getTimeInterpWeight(interpInfo%oti,headerIndex,stepIndex) == 0.0d0 ) cycle header_loop1

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
    allocate(lonVec_r4(numHeaderUsedMax))
    allocate(latVec_r4(numHeaderUsedMax))
    allocate(headerIndexVec(numHeaderUsedMax,numStep))
    headerIndexVec(:,:) = 0

    ! copy the horizontal grid object
    interpInfo%hco => stateVector%hco

    ! setup the information for wind rotation
    if ( (gsv_varExist(varName='UU') .or. gsv_varExist(varName='VV')) .and.  &
         stateVector%hco%rotated ) then
      call uvr_Setup( interpInfo%uvr, & ! INOUT
                      stateVector%hco ) ! IN 
    end if

    ! allocate arrays that will be returned
    allocate(interpInfo%allNumHeaderUsed(numStep,mpi_nprocs))
    allocate(interpInfo%depotIndexBeg(interpInfo%hco%numSubGrid,numHeaderUsedMax,numStep,mpi_nprocs))
    allocate(interpInfo%depotIndexEnd(interpInfo%hco%numSubGrid,numHeaderUsedMax,numStep,mpi_nprocs))
    allocate(interpInfo%allHeaderIndex(numHeaderUsedMax,numStep,mpi_nprocs))
    allocate(interpInfo%allLat(numHeaderUsedMax,numStep,mpi_nprocs))
    allocate(interpInfo%allLon(numHeaderUsedMax,numStep,mpi_nprocs))
    interpInfo%allHeaderIndex(:,:,:) = 0
    interpInfo%allLat(:,:,:) = 0.0d0
    interpInfo%allLon(:,:,:) = 0.0d0
    interpInfo%allNumHeaderUsed(:,:) = allNumHeaderUsed(:,:)

    if ( interpInfo%hco%rotated ) then
      allocate(interpInfo%allLatRot(interpInfo%hco%numSubGrid,numHeaderUsedMax,numStep,mpi_nprocs))
      allocate(interpInfo%allLonRot(interpInfo%hco%numSubGrid,numHeaderUsedMax,numStep,mpi_nprocs))
      interpInfo%allLatRot(:,:,:,:) = 0.0d0
      interpInfo%allLonRot(:,:,:,:) = 0.0d0
    end if

    interpInfo%depotIndexBeg(:,:,:,:) = 0
    interpInfo%depotIndexEnd(:,:,:,:) = -1

    ! get observation lat-lon onto all mpi tasks
    step_loop2: do stepIndex = 1, numStep
      numHeaderUsed = 0

      lonVec_r4(:) = 0.0
      latVec_r4(:) = 0.0

      header_loop2: do headerIndex = 1, numHeader

        ! if obs inside window, but zero weight for current stepIndex then skip it
        if ( oti_getTimeInterpWeight(interpInfo%oti, headerIndex, stepIndex) == 0.0d0 ) cycle header_loop2

        numHeaderUsed = numHeaderUsed + 1
        headerIndexVec(numHeaderUsed,stepIndex) = headerIndex

        !- Get LatLon of observation location
        lat_r4 = real(obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex), 4)
        lon_r4 = real(obs_headElem_r(obsSpaceData, OBS_LON, headerIndex), 4)
        if (lon_r4 <  0.0          ) lon_r4 = lon_r4 + 2.0*MPC_PI_R4
        if (lon_r4 >= 2.0*MPC_PI_R4) lon_r4 = lon_r4 - 2.0*MPC_PI_R4

        lonVec_r4(numHeaderUsed) = lon_r4
        latVec_r4(numHeaderUsed) = lat_r4

        ! check for obs outside domain and reject, if requested
        lat_deg_r4 = lat_r4 * MPC_DEGREES_PER_RADIAN_R8
        lon_deg_r4 = lon_r4 * MPC_DEGREES_PER_RADIAN_R8
        ierr = getPositionXY( stateVector%hco%EZscintID,   &
                              xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                              lat_deg_r4, lon_deg_r4, subGridIndex )

        obsOutsideGrid = ( xpos_r4 < 1.0 .or. xpos_r4 > real(niP1) .or.  &
                           ypos_r4 < 1.0 .or. ypos_r4 > real(stateVector%nj) )

        if ( obsOutsideGrid .and. rejectOutsideObs ) then
          ! The observation is outside the domain
          ! With a LAM trial field we must discard this observation
          write(*,*) 's2c_setupInterpInfo: Rejecting OBS outside the stateVector domain, ', headerIndex
          write(*,*) '  position : ', lat_deg_r4, lon_deg_r4, ypos_r4, xpos_r4

          bodyIndexBeg = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
          bodyIndexEnd = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) + bodyIndexBeg -1
          do bodyIndex = bodyIndexBeg, bodyIndexEnd
            call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
          end do
          call obs_headSet_i(obsSpaceData, OBS_ST1, headerIndex,  &
               ibset( obs_headElem_i(obsSpaceData, OBS_ST1, headerIndex), 05))

        end if

      end do header_loop2

      ! gather geographical lat, lon positions of observations from all processors
      call rpn_comm_allgather(real(latVec_r4,8),                numHeaderUsedMax, 'MPI_REAL8', &
                              interpInfo%allLat(:,stepIndex,:), numHeaderUsedMax, 'MPI_REAL8', &
                              'GRID', ierr)
      call rpn_comm_allgather(real(lonVec_r4,8),                numHeaderUsedMax, 'MPI_REAL8', &
                              interpInfo%allLon(:,stepIndex,:), numHeaderUsedMax, 'MPI_REAL8', &
                              'GRID', ierr)

    end do step_loop2

    ! count the total number of grid points for allocation and set up indices
    numGridptTotal = 0
    do stepIndex = 1, numStep
      do procIndex = 1, mpi_nprocs
        do headerIndex = 1, allNumHeaderUsed(stepIndex,procIndex)

          lat_deg_r4 = real(interpInfo%allLat(headerIndex, stepIndex, procIndex) *  &
                       MPC_DEGREES_PER_RADIAN_R8)
          lon_deg_r4 = real(interpInfo%allLon(headerIndex, stepIndex, procIndex) *  &
                       MPC_DEGREES_PER_RADIAN_R8)
          ierr = getPositionXY( stateVector%hco%EZscintID,   &
                                xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                                lat_deg_r4, lon_deg_r4, subGridIndex )

          if ( (subGridIndex == 1) .or. (subGridIndex == 2) ) then
            ! indices for only 1 subgrid, other will have zeros
            interpInfo%depotIndexBeg(subGridIndex, headerIndex, stepIndex, procIndex) = numGridptTotal + 1
            numGridptTotal = numGridptTotal + 4
            interpInfo%depotIndexEnd(subGridIndex, headerIndex, stepIndex, procIndex) = numGridptTotal
          else
            ! locations on both subGrids will be averaged
            interpInfo%depotIndexBeg(1, headerIndex, stepIndex, procIndex) = numGridptTotal + 1
            numGridptTotal = numGridptTotal + 4
            interpInfo%depotIndexEnd(1, headerIndex, stepIndex, procIndex) = numGridptTotal

            interpInfo%depotIndexBeg(2, headerIndex, stepIndex, procIndex) = numGridptTotal + 1
            numGridptTotal = numGridptTotal + 4
            interpInfo%depotIndexEnd(2, headerIndex, stepIndex, procIndex) = numGridptTotal
          end if

        end do ! headerIndex
      end do ! procIndex
    end do

    ! now that we know the size, allocate main arrays for storing interpolation information
    write(*,*) 's2c_setupInterpInfo: numGridptTotal = ', numGridptTotal
    allocate( interpInfo%latIndexDepot(numGridptTotal) )
    allocate( interpInfo%lonIndexDepot(numGridptTotal) )
    allocate( interpInfo%interpWeightDepot(numGridptTotal) )

    step_loop3: do stepIndex = 1, numStep
      do procIndex = 1, mpi_nprocs
        do headerIndex = 1, allNumHeaderUsed(stepIndex,procIndex)

          lat_deg_r4 = real(interpInfo%allLat(headerIndex, stepIndex, procIndex) *  &
                       MPC_DEGREES_PER_RADIAN_R8)
          lon_deg_r4 = real(interpInfo%allLon(headerIndex, stepIndex, procIndex) *  &
                       MPC_DEGREES_PER_RADIAN_R8)
          ierr = getPositionXY( stateVector%hco%EZscintID,   &
                                xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                                lat_deg_r4, lon_deg_r4, subGridIndex )

          if ( xpos_r4 < 1.0 .or. xpos_r4 > real(niP1) .or.  &
               ypos_r4 < 1.0 .or. ypos_r4 > real(stateVector%nj) ) then

            if ( rejectOutsideObs ) then
              ! Assign a realistic lat-lon to this point for rejected obs
              xpos_r4 = real(stateVector%ni)/2.0
              ypos_r4 = real(stateVector%nj)/2.0
              ierr = gdllfxy(stateVector%hco%EZscintID, lat_deg_r4, lon_deg_r4, &
                             xpos_r4, ypos_r4, 1)
              interpInfo%allLon(headerIndex, stepIndex, procIndex) =  &
                   real(lon_deg_r4 * MPC_RADIANS_PER_DEGREE_R4,8)
              interpInfo%allLat(headerIndex, stepIndex, procIndex) =  &
                   real(lat_deg_r4 * MPC_RADIANS_PER_DEGREE_R4,8)
            else
              write(*,*) 's2c_setupInterpInfo: Moving OBS that is outside the stateVector domain, ', headerIndex
              write(*,*) '  position lon, lat = ', lon_deg_r4, lat_deg_r4
              write(*,*) '  position x,   y   = ', xpos_r4, ypos_r4

              ! if obs above or below domain
              if( ypos_r4 < 1.0 ) ypos_r4 = 1.0
              if( ypos_r4 > real(statevector%nj) ) ypos_r4 = real(statevector%nj)

              ! if obs left or right longitude band, move it to the edge of this longitude band
              if( xpos_r4 < 1.0 ) xpos_r4 = 1.0
              if( xpos_r4 > real(statevector%ni) ) xpos_r4 = real(statevector%ni)
              write(*,*) '  new position x, y = ', xpos_r4, ypos_r4

            end if

          end if

          ! Find the lower-left grid point next to the observation
          if ( xpos_r4 /= real(niP1) ) then
            lonIndex = floor(xpos_r4)
          else
            lonIndex = floor(xpos_r4) - 1
          end if
          if ( xpos2_r4 /= real(niP1) ) then
            lonIndex2 = floor(xpos2_r4)
          else
            lonIndex2 = floor(xpos2_r4) - 1
          end if

          if ( ypos_r4 /= real(statevector%nj) ) then
            latIndex = floor(ypos_r4)
          else
            latIndex = floor(ypos_r4) - 1
          end if
          if ( ypos2_r4 /= real(statevector%nj) ) then
            latIndex2 = floor(ypos2_r4)
          else
            latIndex2 = floor(ypos2_r4) - 1
          end if

          if ( stateVector%hco%grtyp == 'U' ) then
            if ( ypos_r4 /= real(stateVector%nj/2) ) then
              latIndex = floor(ypos_r4)
            else
              latIndex = floor(ypos_r4) - 1
            end if
            if ( ypos2_r4 /= real(stateVector%nj/2) ) then
              latIndex2 = floor(ypos2_r4)
            else
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
              write(*,*) 's2c_setupInterpInfo: WARNING, obs position in between Yin and Yang!!!'
              write(*,*) '   xpos, ypos = ', xpos_r4, ypos_r4
            end if
            if ( ypos2_r4 > real(stateVector%nj/2) .and.  &
                 ypos2_r4 < real((stateVector%nj/2)+1) ) then
              write(*,*) 's2c_setupInterpInfo: WARNING, obs position in between Yin and Yang!!!'
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

          do subGridForInterp = 1, numSubGridsForInterp

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

            indexBeg = interpInfo%depotIndexBeg(subGridIndex, headerIndex, stepIndex, procIndex)
            indexEnd = interpInfo%depotIndexEnd(subGridIndex, headerIndex, stepIndex, procIndex)

            if ( allocated(stateVector%hco%mask) ) then

              if ( stateVector%hco%mask(lonIndex,latIndex) == 1 ) then
                interpInfo%interpWeightDepot(indexBeg    ) = (1.d0-dldx) * (1.d0-dldy)
              else
                interpInfo%interpWeightDepot(indexBeg    ) = 0.d0
              end if

              if ( stateVector%hco%mask(lonIndexP1,latIndex) == 1 ) then
                interpInfo%interpWeightDepot(indexBeg + 1) = dldx  * (1.d0-dldy)
              else
                interpInfo%interpWeightDepot(indexBeg + 1) = 0.d0
              end if

              if ( stateVector%hco%mask(lonIndex,latIndex + 1) == 1 ) then
                interpInfo%interpWeightDepot(indexBeg + 2) = (1.d0-dldx) *       dldy
              else
                interpInfo%interpWeightDepot(indexBeg + 2) = 0.d0
              end if

              if ( stateVector%hco%mask(lonIndexP1,latIndex + 1) == 1 ) then
                interpInfo%interpWeightDepot(indexBeg + 3) =       dldx  *       dldy
              else
                interpInfo%interpWeightDepot(indexBeg + 3) = 0.d0
              end if

              weightsSum = sum(interpInfo%interpWeightDepot(indexBeg:indexBeg + 3))
              if ( weightsSum > 0.d0 ) then

                interpInfo%interpWeightDepot(indexBeg:indexBeg + 3) = &
                interpInfo%interpWeightDepot(indexBeg:indexBeg + 3) / weightsSum

              else

                write(*,*) 's2c_setupInterpInfo: Rejecting OBS outside the grid domain, ', headerIndex
                write(*,*) '  position : ', lat_deg_r4, lon_deg_r4, ypos_r4, xpos_r4
                bodyIndexBeg = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
                bodyIndexEnd = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) + bodyIndexBeg -1
                do bodyIndex = bodyIndexBeg, bodyIndexEnd
                  call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
                end do
                call obs_headSet_i(obsSpaceData, OBS_ST1, headerIndex,  &
                     ibset( obs_headElem_i(obsSpaceData, OBS_ST1, headerIndex), 05))

              end if

            else

              interpInfo%interpWeightDepot(indexBeg    ) = (1.d0-dldx) * (1.d0-dldy)
              interpInfo%interpWeightDepot(indexBeg + 1) = dldx  * (1.d0-dldy)
              interpInfo%interpWeightDepot(indexBeg + 2) = (1.d0-dldx) *       dldy
              interpInfo%interpWeightDepot(indexBeg + 3) =       dldx  *       dldy

            end if

            ! divide weight by number of subGrids
            interpInfo%interpWeightDepot(indexBeg:indexEnd) = &
                 interpInfo%interpWeightDepot(indexBeg:indexEnd) / real(numSubGridsForInterp,8)

            interpInfo%latIndexDepot(indexBeg    ) = latIndex
            interpInfo%latIndexDepot(indexBeg + 1) = latIndex
            interpInfo%latIndexDepot(indexBeg + 2) = latIndex + 1
            interpInfo%latIndexDepot(indexBeg + 3) = latIndex + 1

            interpInfo%lonIndexDepot(indexBeg    ) = lonIndex
            interpInfo%lonIndexDepot(indexBeg + 1) = lonIndexP1
            interpInfo%lonIndexDepot(indexBeg + 2) = lonIndex
            interpInfo%lonIndexDepot(indexBeg + 3) = lonIndexP1

            ! compute the rotated lat, lon
            if ( interpInfo%hco%rotated .and.  &
                 (gsv_varExist(varName='UU') .or.  &
                  gsv_varExist(varName='VV')) ) then
              lat = interpInfo%allLat(headerIndex, stepIndex, procIndex)
              lon = interpInfo%allLon(headerIndex, stepIndex, procIndex)
              call uvr_RotateLatLon( interpInfo%uvr,   & ! INOUT
                                     subGridIndex,     & ! IN
                                     latRot, lonRot,   & ! OUT (radians)
                                     lat, lon,         & ! IN  (radians)
                                     'ToLatLonRot')      ! IN
              interpInfo%allLatRot(subGridIndex, headerIndex, stepIndex, procIndex) = latRot
              interpInfo%allLonRot(subGridIndex, headerIndex, stepIndex, procIndex) = lonRot
            end if

          end do ! subGridForInterp


        end do ! headerIndex
      end do ! procIndex

    end do step_loop3

    ! gather the headerIndexVec arrays onto all processors
    call rpn_comm_allgather(headerIndexVec,    numHeaderUsedMax*numStep, 'MPI_INTEGER', &
                            interpInfo%allHeaderIndex, numHeaderUsedMax*numStep, 'MPI_INTEGER', &
                            'GRID',ierr)

    deallocate(headerIndexVec)
    deallocate(latVec_r4)
    deallocate(lonVec_r4)
    deallocate(allNumHeaderUsed)

    interpInfo%initialized = .true.

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 's2c_setupInterpInfo: FINISHED'

  end subroutine s2c_setupInterpInfo

  !---------------------------------------------------------
  ! s2c_tl
  !---------------------------------------------------------
  subroutine s2c_tl( statevector, column, columng, obsSpaceData )
    ! **Purpose:** 
    ! Tangent linear version of the horizontal
    ! interpolation, used for the increment (or perturbations).
    !
    implicit none

    ! arguments
    type(struct_gsv)           :: stateVector, statevector_trial
    type(struct_obs)           :: obsSpaceData
    type(struct_columnData)    :: column, columng

    ! locals
    type(struct_gsv)           :: stateVector_VarsLevs
    integer :: kIndex, kIndex2, kCount, stepIndex, numStep, mykEndExtended
    integer :: headerIndex, numHeader, numHeaderMax, yourNumHeader
    integer :: procIndex, nsize, ierr, headerUsedIndex
    real(8) :: weight
    real(8), pointer     :: allCols_ptr(:,:)
    real(8), pointer     :: ptr4d(:,:,:,:), ptr3d_UV(:,:,:)
    real(8), allocatable :: cols_hint(:,:,:)
    real(8), allocatable :: cols_send(:,:)
    real(8), allocatable :: cols_recv(:,:)
    real(8), allocatable :: cols_send_1proc(:)
    character(len=4)     :: varName
    real(8), pointer     :: onecolumn(:)

    if ( mpi_myid == 0 ) write(*,*) 's2c_tl: Horizontal interpolation StateVector --> ColumnData'
    call tmg_start(167,'S2C_TL')

    call tmg_start(160,'S2C_BARR')
    call rpn_comm_barrier('GRID',ierr)
    call tmg_stop(160)

    if ( .not. stateVector%allocated ) then 
      call utl_abort('s2c_tl_new: stateVector must be allocated')
    end if

    ! calculate delP_T/delP_M on the grid
    call clp_calcPressure_tl(statevector_trial, statevector)

    call vtr_transform( statevector, & ! INOUT
                        'TTHUtoGZ_tl') ! IN

    call gsv_allocate( statevector_VarsLevs, statevector%numstep, &
                       statevector%hco, statevector%vco,          &
                       mpi_local_opt=.true., mpi_distribution_opt='VarsLevs' )
    call gsv_transposeTilesToVarsLevs( statevector, statevector_VarsLevs )

    numStep = stateVector_VarsLevs%numStep
    numHeader = obs_numheader(obsSpaceData)
    call rpn_comm_allreduce(numHeader, numHeaderMax, 1,  &
                            'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)

    if ( .not. interpInfo_tlad%initialized ) then
      call s2c_setupInterpInfo( interpInfo_tlad, obsSpaceData, stateVector_VarsLevs,  &
                                timeInterpType_tlad,  rejectOutsideObs=.false. )
    end if

    ! arrays for interpolated column for 1 level/variable and each time step
    allocate(cols_hint(maxval(interpInfo_tlad%allNumHeaderUsed),numStep,mpi_nprocs))
    cols_hint(:,:,:) = 0.0d0

    ! arrays for sending/receiving time interpolated column for 1 level/variable
    allocate(cols_send(numHeaderMax,mpi_nprocs))
    cols_send(:,:) = 0.0d0

    allocate(cols_recv(numHeaderMax,mpi_nprocs))
    cols_recv(:,:) = 0.0d0

    allocate(cols_send_1proc(numHeaderMax))

    ! set contents of column to zero
    allCols_ptr => col_getAllColumns(column)
    if ( numHeader > 0 ) allCols_ptr(:,:) = 0.0d0

    ptr4d => gsv_getField_r8(stateVector_VarsLevs)

    mykEndExtended = stateVector_VarsLevs%mykBeg + maxval(stateVector_VarsLevs%allkCount(:)) - 1

    kCount = 0
    k_loop: do kIndex = stateVector_VarsLevs%mykBeg, mykEndExtended

      kCount = kCount + 1

      if ( kIndex <= stateVector_VarsLevs%mykEnd ) then
        varName = gsv_getVarNameFromK(statevector,kIndex)

        if ( varName == 'UU' .or. varName == 'VV' ) then
          ptr3d_UV => gsv_getFieldUV_r8(stateVector_VarsLevs,kIndex)
        end if

        call tmg_start(161,'S2CTL_HINTERP')
        !$OMP PARALLEL DO PRIVATE (stepIndex, procIndex, yourNumHeader, headerIndex)
        step_loop: do stepIndex = 1, numStep
          if ( maxval(interpInfo_tlad%allNumHeaderUsed(stepIndex,:)) == 0 ) cycle step_loop

          ! interpolate to the columns destined for all procs for all steps and one lev/var
          do procIndex = 1, mpi_nprocs
            yourNumHeader = interpInfo_tlad%allNumHeaderUsed(stepIndex,procIndex)
            if ( yourNumHeader > 0 ) then
              if ( varName == 'UU' ) then
                call myezuvint_tl( cols_hint(1:yourNumHeader,stepIndex,procIndex), 'UU',  &
                                   ptr4d(:,:,kIndex,stepIndex), ptr3d_UV(:,:,stepIndex),  &
                                   interpInfo_tlad, stepIndex, procIndex )
              else if ( varName == 'VV' ) then
                call myezuvint_tl( cols_hint(1:yourNumHeader,stepIndex,procIndex), 'VV',  &
                                   ptr3d_UV(:,:,stepIndex), ptr4d(:,:,kIndex,stepIndex),  &
                                   interpInfo_tlad, stepIndex, procIndex )
              else
                call myezsint( cols_hint(1:yourNumHeader,stepIndex,procIndex), varName,  &
                               ptr4d(:,:,kIndex,stepIndex), interpInfo_tlad, stepIndex, procIndex )
              end if
            end if
          end do

        end do step_loop
        !$OMP END PARALLEL DO
        call tmg_stop(161)

        ! interpolate in time to the columns destined for all procs and one level/variable
        do procIndex = 1, mpi_nprocs
          cols_send_1proc(:) = 0.0d0
          do stepIndex = 1, numStep
            !$OMP PARALLEL DO PRIVATE (headerUsedIndex, headerIndex, weight)
            do headerUsedIndex = 1, interpInfo_tlad%allNumHeaderUsed(stepIndex, procIndex)
              headerIndex = interpInfo_tlad%allHeaderIndex(headerUsedIndex,stepIndex,procIndex)
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

      call tmg_start(160,'S2C_BARR')
      call rpn_comm_barrier('GRID',ierr)
      call tmg_stop(160)

      call tmg_start(163,'S2CTL_ALLTOALL')
      ! mpi communication: alltoall for one level/variable
      nsize = numHeaderMax
      if(mpi_nprocs > 1) then
        call rpn_comm_alltoall(cols_send, nsize, 'MPI_REAL8',  &
                               cols_recv, nsize, 'MPI_REAL8', 'GRID', ierr)
      else
        cols_recv(:,1) = cols_send(:,1)
      end if
      call tmg_stop(163)

      ! reorganize ensemble of distributed columns
      !$OMP PARALLEL DO PRIVATE (procIndex, kIndex2, headerIndex)
      proc_loop: do procIndex = 1, mpi_nprocs
        kIndex2 = statevector_VarsLevs%allkBeg(procIndex) + kCount - 1
        if ( kIndex2 > stateVector_VarsLevs%allkEnd(procIndex) ) cycle proc_loop
        do headerIndex = 1, numHeader
          allCols_ptr(kIndex2,headerIndex) = cols_recv(headerIndex,procIndex)
        end do
      end do proc_loop
      !$OMP END PARALLEL DO

    end do k_loop

    ! output the interpolated GZ_T/GZ_M for the first header
    write(*,*) 'statevector->Column 1. GZ_T:'
    onecolumn => col_getColumn(column,1,'GZ_T')
    write(*,*) onecolumn
    write(*,*) 'statevector->Column 1. GZ_M:'
    onecolumn => col_getColumn(column,1,'GZ_M')
    write(*,*) onecolumn

    ! Do final preparations of columnData objects (compute GZ increment)
    if ( col_varExist('TT') .and. col_varExist('HU') .and.  &
         col_varExist('P0') ) then
      call tt2phi_tl(column,columng,obsSpaceData)
    end if

    ! output the interpolated GZ_T/GZ_M for the first header
    write(*,*) 'Psfc->Column 1. GZ_T:'
    onecolumn => col_getColumn(column,1,'GZ_T')
    write(*,*) onecolumn
    write(*,*) 'Psfc->Column 1. GZ_M:'
    onecolumn => col_getColumn(column,1,'GZ_M')
    write(*,*) onecolumn
    nullify(onecolumn)

    deallocate(cols_hint)
    deallocate(cols_send)
    deallocate(cols_recv)
    deallocate(cols_send_1proc)

    call gsv_deallocate( statevector_VarsLevs )

    call tmg_stop(167)

  end subroutine s2c_tl

  !---------------------------------------------------------
  ! s2c_ad
  !---------------------------------------------------------
  subroutine s2c_ad( statevector, column, columng, obsSpaceData )
    ! **Purpose:** 
    ! Adjoint version of the horizontal interpolation,
    ! used for the cost function gradient with respect to the increment.
    !
    implicit none

    ! arguments
    type(struct_gsv)           :: stateVector, statevector_trial
    type(struct_obs)           :: obsSpaceData
    type(struct_columnData)    :: column, columng

    ! locals
    type(struct_gsv)           :: stateVector_VarsLevs
    integer :: kIndex, kIndex2, kCount, stepIndex, numStep, mykEndExtended
    integer :: headerIndex, numHeader, numHeaderMax, yourNumHeader
    integer :: procIndex, nsize, ierr, headerUsedIndex
    character(len=4)     :: varName
    real(8) :: weight
    real(8), pointer     :: allCols_ptr(:,:)
    real(8), pointer     :: ptr4d(:,:,:,:), ptr3d_UV(:,:,:)
    real(8), allocatable :: cols_hint(:,:,:)
    real(8), allocatable :: cols_send(:,:)
    real(8), allocatable :: cols_recv(:,:)
    real(8), pointer :: field_out(:,:,:,:)
    character(len=3) :: execOldNew

    execOldNew = 'new'

    if(mpi_myid == 0) write(*,*) 's2c_ad: Adjoint of horizontal interpolation StateVector --> ColumnData'
    call tmg_start(168,'S2C_AD')

    call tmg_start(160,'S2C_BARR')
    call rpn_comm_barrier('GRID',ierr)
    call tmg_stop(160)

    if ( .not. stateVector%allocated ) then 
      call utl_abort('s2c_ad_new: stateVector must be allocated')
    end if

    call gsv_allocate( statevector_VarsLevs, statevector%numstep, &
                       statevector%hco, statevector%vco,          &
                       mpi_local_opt=.true., mpi_distribution_opt='VarsLevs' )
    call gsv_zero( statevector_VarsLevs )

    if ( .not. interpInfo_tlad%initialized ) then
      call s2c_setupInterpInfo( interpInfo_tlad, obsSpaceData, stateVector_VarsLevs,  &
                                timeInterpType_tlad, rejectOutsideObs=.false. )
    end if

    ! Mass fields (TT,PS,HU) to hydrostatic geopotential
    if (gsv_varExist(statevector,'TT') .and.  &
        gsv_varExist(statevector,'HU') .and. gsv_varExist(statevector,'P0') ) then
       call tt2phi_ad(column,columng,obsSpaceData)
    end if

    numStep = stateVector_VarsLevs%numStep
    numHeader = obs_numheader(obsSpaceData)
    call rpn_comm_allreduce(numHeader, numHeaderMax, 1,  &
                            'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)

    ! arrays for interpolated column for 1 level/variable and each time step
    allocate(cols_hint(maxval(interpInfo_tlad%allNumHeaderUsed),numStep,mpi_nprocs))
    cols_hint(:,:,:) = 0.0d0

    ! arrays for sending/receiving time interpolated column for 1 level/variable
    allocate(cols_send(numHeaderMax,mpi_nprocs))
    cols_send(:,:) = 0.0d0

    allocate(cols_recv(numHeaderMax,mpi_nprocs))
    cols_recv(:,:) = 0.0d0

    ! set contents of column to zero
    allCols_ptr => col_getAllColumns(column)

    ptr4d => gsv_getField_r8(stateVector_VarsLevs)
    mykEndExtended = stateVector_VarsLevs%mykBeg + maxval(stateVector_VarsLevs%allkCount(:)) - 1

    kCount = 0
    k_loop: do kIndex = stateVector_VarsLevs%mykBeg, mykEndExtended

      kCount = kCount + 1

      ! reorganize ensemble of distributed columns
      !$OMP PARALLEL DO PRIVATE (procIndex, kIndex2, headerIndex)
      proc_loop: do procIndex = 1, mpi_nprocs
        kIndex2 = statevector_VarsLevs%allkBeg(procIndex) + kCount - 1
        if ( kIndex2 > stateVector_VarsLevs%allkEnd(procIndex) ) cycle proc_loop
        do headerIndex = 1, numHeader
          cols_send(headerIndex,procIndex) = allCols_ptr(kIndex2,headerIndex)
        end do
      end do proc_loop
      !$OMP END PARALLEL DO

      call tmg_start(160,'S2C_BARR')
      call rpn_comm_barrier('GRID',ierr)
      call tmg_stop(160)

      call tmg_start(164,'S2CAD_ALLTOALL')
      ! mpi communication: alltoall for one level/variable
      nsize = numHeaderMax
      if(mpi_nprocs > 1) then
        call rpn_comm_alltoall(cols_send, nsize, 'MPI_REAL8',  &
                               cols_recv, nsize, 'MPI_REAL8', 'GRID', ierr)
      else
        cols_recv(:,1) = cols_send(:,1)
      end if
      call tmg_stop(164)

      if ( kIndex <= stateVector_VarsLevs%mykEnd ) then
        varName = gsv_getVarNameFromK(statevector,kIndex)

        if ( varName == 'UU' .or. varName == 'VV' ) then
          ptr3d_UV => gsv_getFieldUV_r8(stateVector_VarsLevs,kIndex)
        end if

        ! interpolate in time to the columns destined for all procs and one level/variable
        do procIndex = 1, mpi_nprocs
          do stepIndex = 1, numStep
            !$OMP PARALLEL DO PRIVATE (headerIndex, headerUsedIndex, weight)
            do headerIndex = 1, interpInfo_tlad%allNumHeaderUsed(stepIndex, procIndex)

              headerUsedIndex = interpInfo_tlad%allHeaderIndex(headerIndex,stepIndex,procIndex)
              weight = oti_getTimeInterpWeightMpiGlobal(interpInfo_tlad%oti,  &
                                                        headerUsedIndex,stepIndex,procIndex)

              cols_hint(headerIndex,stepIndex,procIndex) =  &
                   weight * cols_recv(headerUsedIndex,procIndex)

            end do
            !$OMP END PARALLEL DO
          end do
        end do

        call tmg_start(162,'S2CAD_HINTERP')
        !$OMP PARALLEL DO PRIVATE (stepIndex, procIndex, yourNumHeader)
        step_loop: do stepIndex = 1, numStep
          if ( maxval(interpInfo_tlad%allNumHeaderUsed(stepIndex,:)) == 0 ) cycle step_loop

          ! interpolate to the columns destined for all procs for all steps and one lev/var
          do procIndex = 1, mpi_nprocs
            yourNumHeader = interpInfo_tlad%allNumHeaderUsed(stepIndex,procIndex)
            if ( yourNumHeader > 0 ) then
              if ( varName == 'UU' ) then
                call myezuvint_ad( cols_hint(1:yourNumHeader,stepIndex,procIndex), 'UU',  &
                                   ptr4d(:,:,kIndex,stepIndex), ptr3d_UV(:,:,stepIndex),  &
                                   interpInfo_tlad, stepIndex, procIndex )
              else if ( varName == 'VV' ) then
                call myezuvint_ad( cols_hint(1:yourNumHeader,stepIndex,procIndex), 'VV',  &
                                   ptr3d_UV(:,:,stepIndex), ptr4d(:,:,kIndex,stepIndex),  &
                                   interpInfo_tlad, stepIndex, procIndex )
              else
                call myezsintad( cols_hint(1:yourNumHeader,stepIndex,procIndex), varName,  &
                                 ptr4d(:,:,kIndex,stepIndex), interpInfo_tlad, stepIndex,  &
                                 procIndex )
              end if
            end if
          end do

        end do step_loop
        !$OMP END PARALLEL DO
        call tmg_stop(162)

      end if ! if kIndex <= mykEnd

    end do k_loop

    deallocate(cols_hint)
    deallocate(cols_send)
    deallocate(cols_recv)

    call tmg_start(160,'S2C_BARR')
    call rpn_comm_barrier('GRID',ierr)
    call tmg_stop(160)

    call gsv_transposeTilesToVarsLevsAd( statevector_VarsLevs, statevector )

    !if ( execOldNew == 'new' ) then
      call vtr_transform( statevector, & ! INOUT
                          'TTHUtoGZ_ad') ! IN
    !end if 

    ! Adjoint of calculate delP_T/delP_M on the grid
    call clp_calcPressure_ad(statevector_trial, statevector)

    write(*,*) 'MAZIAR, s2c_ad for execOldNew=', execOldNew
    write(*,*) 'TT='
    field_out => gsv_getField_r8(statevector,'TT')
    write(*,*) field_out(statevector%myLonBeg,statevector%myLatBeg,:,1)
    write(*,*) 'HU='
    field_out => gsv_getField_r8(statevector,'HU')
    write(*,*) field_out(statevector%myLonBeg,statevector%myLatBeg,:,1)
    write(*,*) 'P0='
    field_out => gsv_getField_r8(statevector,'P0')
    write(*,*) field_out(statevector%myLonBeg,statevector%myLatBeg,:,1)

    call gsv_deallocate( statevector_VarsLevs )

    call tmg_stop(168)

  end subroutine s2c_ad

  !---------------------------------------------------------
  ! s2c_nl
  !---------------------------------------------------------
  subroutine s2c_nl( stateVector, obsSpaceData, column, timeInterpType, varName_opt, &
                     dealloc_opt, moveObsAtPole_opt )
    ! **Purpose:** 
    ! Non-linear version of the horizontal interpolation,
    ! used for a full field (usually the background state when computing
    ! the innovation vector).
    !
    implicit none

    ! arguments
    type(struct_gsv)           :: stateVector
    type(struct_obs)           :: obsSpaceData
    type(struct_columnData)    :: column
    character(len=*)           :: timeInterpType
    character(len=*), optional :: varName_opt
    logical, optional          :: dealloc_opt
    logical, optional          :: moveObsAtPole_opt

    ! locals
    integer :: kIndex, kIndex2, kCount, levIndex, stepIndex, numStep, mykEndExtended
    integer :: headerIndex, numHeader, numHeaderMax, yourNumHeader
    integer :: procIndex, nsize, ierr, iset, headerUsedIndex, varIndex
    integer :: ezdefset, ezqkdef
    integer :: s1, s2, s3, s4, index1, index2, index3, index4
    real(8) :: gzSfc_col
    real(8) :: weight
    character(len=4)     :: varName
    real(8), pointer     :: column_ptr(:), ptr2d_r8(:,:), allCols_ptr(:,:)
    real(4), pointer     :: ptr4d_r4(:,:,:,:), ptr3d_UV_r4(:,:,:)
    real(8), allocatable :: field2d(:,:), field2d_UV(:,:)
    real(8), allocatable :: cols_hint(:,:,:)
    real(8), allocatable :: cols_send(:,:)
    real(8), allocatable :: cols_recv(:,:)
    real(8), allocatable :: cols_send_1proc(:)
    integer, allocatable :: displs(:), nsizes(:)
    logical              :: beSilent, dealloc, moveObsAtPole

    call tmg_start(169,'S2C_NL')

    write(*,*) 's2c_nl: STARTING'

    if ( .not. stateVector%allocated ) then 
      call utl_abort('s2c_nl: stateVector must be allocated')
    end if

    if ( stateVector%mpi_distribution /= 'VarsLevs' ) then 
      call utl_abort('s2c_nl: stateVector must by VarsLevs distributed')
    end if

    if ( present(dealloc_opt) ) then
      dealloc = dealloc_opt
    else
      dealloc = .true.
    end if

    if ( present(moveObsAtPole_opt) ) then
      moveObsAtPole = moveObsAtPole_opt
    else
      moveObsAtPole = .false.
    end if

    numStep = stateVector%numStep
    numHeader = obs_numheader(obsSpaceData)
    call rpn_comm_allreduce(numHeader, numHeaderMax, 1,  &
                            'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)

    if ( .not. interpInfo_nl%initialized ) then
      call tmg_start(165,'S2CNL_SETUPS')
      ! also reject obs outside (LAM) domain and optionally move obs near 
      ! numerical pole to first/last analysis grid latitude
      call s2c_latLonChecks( obsSpaceData, moveObsAtPole )

      ! compute and collect all obs grids onto all mpi tasks
      call s2c_setupInterpInfo( interpInfo_nl, obsSpaceData, stateVector,  &
                                timeInterpType, rejectOutsideObs=.true. )
      if ( mpi_myid == 0 ) then
        do stepIndex = 1, numStep
          write(*,*) 's2c_nl: stepIndex, allNumHeaderUsed = ',  &
                     stepIndex, interpInfo_nl%allNumHeaderUsed(stepIndex,:)
        end do
      end if
      call tmg_stop(165)
    end if

    ! arrays for interpolated column for 1 level/variable and each time step
    allocate(cols_hint(maxval(interpInfo_nl%allNumHeaderUsed),numStep,mpi_nprocs))
    cols_hint(:,:,:) = 0.0d0

    ! arrays for sending/receiving time interpolated column for 1 level/variable
    allocate(cols_send(numHeaderMax,mpi_nprocs))
    cols_send(:,:) = 0.0d0

    allocate(cols_recv(numHeaderMax,mpi_nprocs))
    cols_recv(:,:) = 0.0d0

    allocate(cols_send_1proc(numHeaderMax))

    ! set contents of column to zero (1 variable or all)
    allCols_ptr => col_getAllColumns(column,varName_opt)
    if ( numHeader > 0 ) allCols_ptr(:,:) = 0.0d0

    ptr4d_r4    => gsv_getField_r4(stateVector)

    allocate(field2d(stateVector%ni,stateVector%nj))
    allocate(field2d_UV(stateVector%ni,stateVector%nj))

    mykEndExtended = stateVector%mykBeg + maxval(stateVector%allkCount(:)) - 1

    kCount = 0
    k_loop: do kIndex = stateVector%mykBeg, mykEndExtended
      kCount = kCount + 1

      if ( kIndex <= stateVector%mykEnd ) then
        varName = gsv_getVarNameFromK(statevector,kIndex)

        if ( varName == 'UU' .or. varName == 'VV' ) then
          ptr3d_UV_r4 => gsv_getFieldUV_r4(stateVector,kIndex)
        end if

        call tmg_start(166,'S2CNL_HINTERP')
        !$OMP PARALLEL DO PRIVATE (stepIndex, field2d, field2d_UV, procIndex, yourNumHeader)
        step_loop: do stepIndex = 1, numStep
          if ( maxval(interpInfo_nl%allNumHeaderUsed(stepIndex,:)) == 0 ) cycle step_loop

          ! copy over field
          field2d(:,:) = real(ptr4d_r4(:,:,kIndex,stepIndex),8)
          if ( varName == 'UU' .or. varName == 'VV' ) then
            field2d_UV(:,:) = real(ptr3d_UV_r4(:,:,stepIndex),8)
          end if

          ! interpolate to the columns destined for all procs for all steps and one lev/var
          do procIndex = 1, mpi_nprocs
            yourNumHeader = interpInfo_nl%allNumHeaderUsed(stepIndex,procIndex)
            if ( yourNumHeader > 0 ) then
              if ( varName == 'UU' ) then
                call myezuvint_nl( cols_hint(1:yourNumHeader,stepIndex,procIndex), 'UU',  &
                                   field2d, field2d_UV, interpInfo_nl, stepIndex, procIndex )
              else if ( varName == 'VV' ) then
                call myezuvint_nl( cols_hint(1:yourNumHeader,stepIndex,procIndex), 'VV',  &
                                   field2d_UV, field2d, interpInfo_nl, stepIndex, procIndex )
              else
                call myezsint( cols_hint(1:yourNumHeader,stepIndex,procIndex), varName,  &
                               field2d, interpInfo_nl, stepIndex, procIndex )
              end if
            end if
          end do

        end do step_loop
        !$OMP END PARALLEL DO
        call tmg_stop(166)

        ! interpolate in time to the columns destined for all procs and one level/variable
        do procIndex = 1, mpi_nprocs
          cols_send_1proc(:) = 0.0d0
          do stepIndex = 1, numStep
            !$OMP PARALLEL DO PRIVATE (headerIndex, headerUsedIndex, weight)
            do headerIndex = 1, interpInfo_nl%allNumHeaderUsed(stepIndex, procIndex)
              headerUsedIndex = interpInfo_nl%allHeaderIndex(headerIndex,stepIndex,procIndex)
              weight = oti_getTimeInterpWeightMpiGlobal(interpInfo_nl%oti,  &
                                                        headerUsedIndex,stepIndex,procIndex)
              cols_send_1proc(headerUsedIndex) = cols_send_1proc(headerUsedIndex) &
                            + weight * cols_hint(headerIndex,stepIndex,procIndex)

            end do
            !$OMP END PARALLEL DO
          end do
          cols_send(:,procIndex) = cols_send_1proc(:)
        end do

      else

        ! this value of k does not exist on this mpi task
        cols_send(:,:) = 0.0d0

      end if ! if kIndex <= mykEnd

      call tmg_start(160,'S2C_BARR')
      call rpn_comm_barrier('GRID',ierr)
      call tmg_stop(160)

      ! mpi communication: alltoall for one level/variable
      nsize = numHeaderMax
      if(mpi_nprocs > 1) then
        call rpn_comm_alltoall(cols_send, nsize, 'MPI_REAL8',  &
                               cols_recv, nsize, 'MPI_REAL8', 'GRID', ierr)
      else
        cols_recv(:,1) = cols_send(:,1)
      end if

      ! reorganize ensemble of distributed columns
      !$OMP PARALLEL DO PRIVATE (procIndex, kIndex2, headerIndex)
      proc_loop: do procIndex = 1, mpi_nprocs
        kIndex2 = statevector%allkBeg(procIndex) + kCount - 1
        if ( kIndex2 > stateVector%allkEnd(procIndex) ) cycle proc_loop
        do headerIndex = 1, numHeader
          allCols_ptr(kIndex2,headerIndex) = cols_recv(headerIndex,procIndex)
        end do
      end do proc_loop
      !$OMP END PARALLEL DO

    end do k_loop

    ! Interpolate surface GZ separately, only exists on mpi task 0
    GZsfcPresent: if ( stateVector%GZsfcPresent ) then

      if ( mpi_myid == 0 ) then
        varName = 'GZ'
        step_loop_gz: do stepIndex = 1, numStep

          if ( maxval(interpInfo_nl%allNumHeaderUsed(stepIndex,:)) == 0 ) cycle step_loop_gz

          ! interpolate to the columns destined for all procs for all steps and one lev/var
          !$OMP PARALLEL DO PRIVATE (procIndex, yourNumHeader, ptr2d_r8)
          do procIndex = 1, mpi_nprocs
            yourNumHeader = interpInfo_nl%allNumHeaderUsed(stepIndex,procIndex)
            if ( yourNumHeader > 0 ) then
              ptr2d_r8 => gsv_getGZsfc(stateVector)
              call myezsint( cols_hint(1:yourNumHeader,stepIndex,procIndex), varName,  &
                             ptr2d_r8(:,:), interpInfo_nl, stepIndex, procIndex )

            end if
          end do
          !$OMP END PARALLEL DO

        end do step_loop_gz

        ! interpolate in time to the columns destined for all procs and one level/variable
        do procIndex = 1, mpi_nprocs
          cols_send(:,procIndex) = 0.0d0
          do stepIndex = 1, numStep
            !$OMP PARALLEL DO PRIVATE (headerIndex, headerUsedIndex)
            do headerIndex = 1, interpInfo_nl%allNumHeaderUsed(stepIndex, procIndex)
              headerUsedIndex = interpInfo_nl%allHeaderIndex(headerIndex,stepIndex,procIndex)
              ! just copy, since surface GZ same for all time steps
              cols_send(headerUsedIndex,procIndex) = cols_hint(headerIndex,stepIndex,procIndex)
            end do
            !$OMP END PARALLEL DO
          end do
        end do

      end if

      ! mpi communication: scatter data from task 0
      nsize = numHeaderMax
      if(mpi_nprocs > 1) then
        allocate(displs(mpi_nprocs))
        allocate(nsizes(mpi_nprocs))
        do procIndex = 1, mpi_nprocs
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

      do headerIndex = 1, numHeader
        call col_setGZsfc(column, headerIndex, cols_recv(headerIndex,1))
      end do

    end if GZsfcPresent

    deallocate(field2d_UV)
    deallocate(field2d)
    deallocate(cols_hint)
    deallocate(cols_send)
    deallocate(cols_recv)
    deallocate(cols_send_1proc)

    if( dealloc ) then
      deallocate(interpInfo_nl%interpWeightDepot)
      deallocate(interpInfo_nl%latIndexDepot)
      deallocate(interpInfo_nl%lonIndexDepot)
      if ( interpInfo_nl%hco%rotated ) then
        deallocate(interpInfo_nl%allLonRot)
        deallocate(interpInfo_nl%allLatRot)
      end if
      deallocate(interpInfo_nl%allLon)
      deallocate(interpInfo_nl%allLat)
      deallocate(interpInfo_nl%allHeaderIndex)
      deallocate(interpInfo_nl%depotIndexBeg)
      deallocate(interpInfo_nl%depotIndexEnd)
      deallocate(interpInfo_nl%allNumHeaderUsed)
      call oti_deallocate(interpInfo_nl%oti)

      interpInfo_nl%initialized = .false.
    end if

    ! impose a lower limit on HU
    if(col_varExist('HU')) then
      do headerIndex = 1, numHeader
        column_ptr => col_getColumn(column,headerIndex,'HU')
        column_ptr(:) = max(column_ptr(:),col_rhumin)
      end do
    end if

    write(*,*) 's2c_nl: FINISHED'

    call tmg_stop(169)

  end subroutine s2c_nl

  ! -------------------------------------------------
  ! myezsint: Scalar field horizontal interpolation
  ! -------------------------------------------------
  subroutine myezsint( column_out, varName, field_in, interpInfo, stepIndex, procIndex ) 
    ! **Purpose:** 
    ! Scalar horizontal interpolation, replaces the
    ! ezsint routine from rmnlib.
    !
    implicit none

    ! arguments
    real(8)                 :: column_out(:)
    character(len=*)        :: varName
    real(8)                 :: field_in(:,:)
    type(struct_interpInfo) :: interpInfo
    integer                 :: stepIndex, procIndex

    ! locals
    integer :: lonIndex, latIndex, gridptIndex, headerIndex, subGridIndex, numColumn
    real(8) :: interpValue, weight

    numColumn = size( column_out )

    do headerIndex = 1, numColumn

      ! Interpolate the model state to the obs point
      interpValue = 0.0d0

      do subGridIndex = 1, interpInfo%hco%numSubGrid

        do gridptIndex =  &
             interpInfo%depotIndexBeg(subGridIndex, headerIndex, stepIndex, procIndex), &
             interpInfo%depotIndexEnd(subGridIndex, headerIndex, stepIndex, procIndex)

          lonIndex = interpInfo%lonIndexDepot(gridptIndex)
          latIndex = interpInfo%latIndexDepot(gridptIndex)
          weight = interpInfo%interpWeightDepot(gridptIndex)

          interpValue = interpValue + weight * field_in(lonIndex, latIndex)

        end do

      end do
      column_out(headerIndex) = interpValue

    end do

  end subroutine myezsint

  ! -------------------------------------------------------------
  ! myezsintad: Adjoint of scalar field horizontal interpolation
  ! -------------------------------------------------------------
  subroutine myezsintad( column_in, varName, field_out, interpInfo, stepIndex, procIndex ) 
    ! **Purpose:** 
    ! Adjoint of the scalar horizontal interpolation.
    !
    implicit none

    ! arguments
    real(8)                 :: column_in(:)
    character(len=*)        :: varName
    real(8)                 :: field_out(:,:)
    type(struct_interpInfo) :: interpInfo
    integer                 :: stepIndex, procIndex

    ! locals
    integer :: lonIndex, latIndex, gridptIndex, headerIndex, subGridIndex, numColumn
    real(8) :: weight

    numColumn = size( column_in )

    do headerIndex = 1, numColumn

      ! Interpolate the model state to the obs point

      do subGridIndex = 1, interpInfo%hco%numSubGrid

        do gridptIndex =  &
             interpInfo%depotIndexBeg(subGridIndex, headerIndex, stepIndex, procIndex), &
             interpInfo%depotIndexEnd(subGridIndex, headerIndex, stepIndex, procIndex)

          lonIndex = interpInfo%lonIndexDepot(gridptIndex)
          latIndex = interpInfo%latIndexDepot(gridptIndex)
          weight = interpInfo%interpWeightDepot(gridptIndex)

          field_out(lonIndex, latIndex) = field_out(lonIndex, latIndex) +  &
                                          weight * column_in(headerIndex)

        end do

      end do

    end do

  end subroutine myezsintad

  ! -------------------------------------------------------------
  ! myezuvint_nl: Vector field horizontal interpolation
  ! -------------------------------------------------------------
  subroutine myezuvint_nl( column_out, varName, fieldUU_in, fieldVV_in,  &
                           interpInfo, stepIndex, procIndex ) 
    ! **Purpose:** 
    ! Vector horizontal interpolation, replaces the
    ! ezuvint routine from rmnlib.
    !
    implicit none

    ! arguments
    real(8)                 :: column_out(:)
    character(len=*)        :: varName
    real(8)                 :: fieldUU_in(:,:), fieldVV_in(:,:)
    type(struct_interpInfo) :: interpInfo
    integer                 :: stepIndex, procIndex

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

        indexBeg = interpInfo%depotIndexBeg(subGridIndex, headerIndex, stepIndex, procIndex)
        indexEnd = interpInfo%depotIndexEnd(subGridIndex, headerIndex, stepIndex, procIndex)

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
          lat = interpInfo%allLat(headerIndex, stepIndex, procIndex)
          lon = interpInfo%allLon(headerIndex, stepIndex, procIndex)
          latRot = interpInfo%allLatRot(subGridIndex, headerIndex, stepIndex, procIndex)
          lonRot = interpInfo%allLonRot(subGridIndex, headerIndex, stepIndex, procIndex)

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
                           interpInfo, stepIndex, procIndex ) 
    ! **Purpose:** 
    ! Vector horizontal interpolation, replaces the
    ! ezuvint routine from rmnlib.
    !
    implicit none

    ! arguments
    real(8)                 :: column_out(:)
    character(len=*)        :: varName
    real(8)                 :: fieldUU_in(:,:), fieldVV_in(:,:)
    type(struct_interpInfo) :: interpInfo
    integer                 :: stepIndex, procIndex

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

        indexBeg = interpInfo%depotIndexBeg(subGridIndex, headerIndex, stepIndex, procIndex)
        indexEnd = interpInfo%depotIndexEnd(subGridIndex, headerIndex, stepIndex, procIndex)

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
          lat = interpInfo%allLat(headerIndex, stepIndex, procIndex)
          lon = interpInfo%allLon(headerIndex, stepIndex, procIndex)
          latRot = interpInfo%allLatRot(subGridIndex, headerIndex, stepIndex, procIndex)
          lonRot = interpInfo%allLonRot(subGridIndex, headerIndex, stepIndex, procIndex)

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
                           interpInfo, stepIndex, procIndex ) 
    ! **Purpose:** 
    ! Adjoint of the vector horizontal interpolation.
    !
    implicit none

    ! arguments
    real(8)                 :: column_in(:)
    character(len=*)        :: varName
    real(8)                 :: fieldUU_out(:,:), fieldVV_out(:,:)
    type(struct_interpInfo) :: interpInfo
    integer                 :: stepIndex, procIndex

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

        indexBeg = interpInfo%depotIndexBeg(subGridIndex, headerIndex, stepIndex, procIndex)
        indexEnd = interpInfo%depotIndexEnd(subGridIndex, headerIndex, stepIndex, procIndex)

        if ( indexEnd < IndexBeg ) cycle subGrid_loop

        ! now rotate the wind vector and return the desired component
        if ( interpInfo%hco%rotated ) then
          lat = interpInfo%allLat(headerIndex, stepIndex, procIndex)
          lon = interpInfo%allLon(headerIndex, stepIndex, procIndex)
          latRot = interpInfo%allLatRot(subGridIndex, headerIndex, stepIndex, procIndex)
          lonRot = interpInfo%allLonRot(subGridIndex, headerIndex, stepIndex, procIndex)

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

  !------------------------------------------------------------------
  ! getPositionXY
  !------------------------------------------------------------------
  function getPositionXY( gdid, xpos_r4, ypos_r4, xpos2_r4, ypos2_r4,  &
                          lat_deg_r4, lon_deg_r4, subGridIndex ) result(ierr)
    ! **Purpose:** 
    ! Compute the grid XY position from a lat-lon. This
    ! simply calls the ezsint routine gdxyfll for simple grids. For
    ! Yin-Yan grids it can return locations from both the Yin and Yan
    ! subgrids when in the overlap region, depending on the logical 
    ! variable `useSingleValueOverlap`.
    !
    implicit none

    ! arguments
    integer :: ierr
    integer :: gdid
    integer :: subGridIndex
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    real(4) :: lat_deg_r4, lon_deg_r4

    ! locals
    integer :: numSubGrids
    integer :: ezget_nsubGrids, ezget_subGridids, gdxyfll, ezgprm, gdgaxes
    integer, allocatable :: EZscintIDvec(:)
    character(len=1) :: grtyp
    integer :: ni, nj, ig1, ig2, ig3, ig4, lonIndex, latIndex
    real :: lonrot, latrot
    real, allocatable :: ax_yin(:), ay_yin(:), ax_yan(:), ay_yan(:)

    ! this controls which approach to use for interpolation within the YIN-YAN overlap
    logical :: useSingleValueOverlap = .true.  

    numSubGrids = ezget_nsubGrids(gdid)
    xpos2_r4 = -999.0
    ypos2_r4 = -999.0

    if ( numSubGrids == 1 ) then

      ! Not a Yin-Yang grid, call the standard ezscint routine
      ierr = gdxyfll(gdid, xpos_r4, ypos_r4, lat_deg_r4, lon_deg_r4, 1)
      subGridIndex = 1

    else

      ! This is a Yin-Yang grid, do something different

      allocate(EZscintIDvec(numSubGrids))
      ierr = ezget_subGridids(gdid, EZscintIDvec)   
      ! get ni nj of subGrid, assume same for both YIN and YANG
      ierr = ezgprm(EZscintIDvec(1), grtyp, ni, nj, ig1, ig2, ig3, ig4)

      ! first check YIN
      ierr = gdxyfll(EZscintIDvec(1), xpos_r4, ypos_r4, lat_deg_r4, lon_deg_r4, 1)

      ! compute rotated lon and lat at obs location
      allocate(ax_yin(ni),ay_yin(nj))
      ierr = gdgaxes(EZscintIDvec(1), ax_yin, ay_yin)
      lonIndex = floor(xpos_r4)
      if ( lonIndex >= 1 .and. (lonIndex+1) <= ni ) then
        lonrot = ax_yin(lonIndex) + (ax_yin(lonIndex+1) - ax_yin(lonIndex)) *  &
                 (xpos_r4 - lonIndex)
      else
        lonrot = -999.0
      end if
      latIndex = floor(ypos_r4)
      if ( latIndex >= 1 .and. (latIndex+1) <= nj ) then
        latrot = ay_yin(latIndex) + (ay_yin(latIndex+1) - ay_yin(latIndex)) *  &
                 (ypos_r4 - latIndex)
      else
        latrot = -999.0
      end if
      deallocate(ax_yin,ay_yin)
      subGridIndex = 1

      if ( useSingleValueOverlap ) then

        ! this approach is most similar to how ezsint works, preferentially take YIN

        if ( lonrot < 45.0 .or. lonrot > 315.0 .or. latrot < -45.0 .or. latrot > 45.0 ) then
          ! Outside YIN, therefore use YANG (assume it is inside YANG)
          ierr = gdxyfll(EZscintIDvec(2), xpos_r4, ypos_r4, lat_deg_r4, lon_deg_r4, 1)
          ypos_r4 = ypos_r4 + real(nj) ! shift from YANG position to Supergrid position
          subGridIndex = 2
        else
          subGridIndex = 1
        end if

      else ! not useSingleValueOverlap

        ! this approach returns both the YIN and YAN locations when point is inside both

        if ( lonrot < 45.0 .or. lonrot > 315.0 .or. latrot < -45.0 .or. latrot > 45.0 ) then
          ! Outside YIN, therefore use YANG (assume it is inside YANG)
          ierr = gdxyfll(EZscintIDvec(2), xpos_r4, ypos_r4, lat_deg_r4, lon_deg_r4, 1)
          ypos_r4 = ypos_r4 + real(nj) ! shift from YANG position to Supergrid position
          subGridIndex = 2
        else
          ! inside YIN, check if also inside YANG
          allocate(ax_yan(ni),ay_yan(nj))
          ierr = gdgaxes(EZscintIDvec(2), ax_yan, ay_yan)
          ierr = gdxyfll(EZscintIDvec(2), xpos2_r4, ypos2_r4, lat_deg_r4, lon_deg_r4, 1)
          if ( lonIndex >= 1 .and. (lonIndex+1) <= ni ) then
            lonrot = ax_yan(lonIndex) + (ax_yan(lonIndex+1) - ax_yan(lonIndex)) *  &
                     (xpos2_r4 - lonIndex)
          else
            lonrot = -999.0
          end if
          latIndex = floor(ypos2_r4)
          if ( latIndex >= 1 .and. (latIndex+1) <= nj ) then
            latrot = ay_yan(latIndex) + (ay_yan(latIndex+1) - ay_yan(latIndex)) *  &
                     (ypos2_r4 - latIndex)
          else
            latrot = -999.0
          end if
          deallocate(ax_yan,ay_yan)
          if ( lonrot < 45.0 .or. lonrot > 315.0 .or. latrot < -45.0 .or. latrot > 45.0 ) then
            ! outside YANG, only inside YIN
            xpos2_r4 = -999.0
            ypos2_r4 = -999.0
            subGridIndex = 1
          else
            ! inside both YIN and YANG
            ypos2_r4 = ypos2_r4 + real(nj) ! shift from YANG position to Supergrid position
            subGridIndex = 3
          end if
        end if

      end if

      deallocate(EZscintIDvec)

    end if    

    if ( subGridIndex /= 3 ) then
      ! when only returning 1 position, copy values to pos2
      xpos2_r4 = xpos_r4
      ypos2_r4 = ypos_r4
    end if

  end function getPositionXY

  !---------------------------------------------------------
  ! s2c_bgcheck_bilin
  !---------------------------------------------------------
  subroutine s2c_bgcheck_bilin(column,statevector,obsSpaceData)
    ! **Purpose:**
    ! Special version of s2c_tl used for background check. This should
    ! be replaced by direct call to s2c_tl. It is not general enough to
    ! be used for new analysis variables.
    !
    implicit none

    ! arguments
    type(struct_columnData) :: column
    type(struct_gsv) :: statevector
    type(struct_obs) :: obsSpaceData

    ! locals
    integer :: jlev, jk, jk2, jgl, jlon, headerIndex
    integer :: lonIndex, ila, ierr, subGridIndex
    integer :: extraLongitude
    real(8) :: lat, lon
    real(4) :: lat_r4, lon_r4, lat_deg_r4, lon_deg_r4, xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    real(8) :: dldy, dlw1, dlw2, dlw3, dlw4, dldx, ypos, xpos
    real(8), allocatable ::zgd(:,:,:)
    real(8), pointer :: uu_column(:),vv_column(:),hu_column(:)
    real(8), pointer :: tt_column(:),tr_column(:),ps_column(:),tg_column(:)
    real(8), pointer :: vis_column(:),gust_column(:)
    real(8), pointer :: field_ptr(:,:,:), uu_ptr(:,:,:), vv_ptr(:,:,:)

    ! Note: We assume here the all the obs between the poles and the last grid points
    !       (i.e. outside the grid) have been moved within the grid by suprep

    if (statevector%hco%global) then
      extraLongitude = 1
    else
      extraLongitude = 0
    end if

    allocate(zgd(statevector%ni+extraLongitude,statevector%nj,statevector%nk))
  
    zgd(:,:,:)=0.0d0
    field_ptr => gsv_getField3D_r8(statevector)
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
      ierr = getPositionXY( stateVector % hco % EZscintID,   &
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
      if ( xpos /= real(statevector%ni + extraLongitude,8) ) then
        lonIndex = floor(xpos)
      else
        lonIndex = floor(xpos) - 1
      end if

      if ( ypos /= real(statevector%nj,8) ) then
        ILA = floor(ypos)
      else
        ILA = floor(ypos) - 1
      end if

      !- 2.3 Compute the 4 weights of the bilinear interpolation
      dldx = xpos - real(lonIndex,8)
      dldy = ypos - real(ILA,8)

      dlw1 = (1.d0-dldx) * (1.d0-dldy)
      dlw2 =       dldx  * (1.d0-dldy)
      dlw3 = (1.d0-dldx) *       dldy
      dlw4 =       dldx  *       dldy

      !- 2.4 Interpolate the model state to the obs point
      if(col_varExist('UU'))  uu_column   => col_getColumn(column,headerIndex,'UU')
      if(col_varExist('VV'))  vv_column   => col_getColumn(column,headerIndex,'VV')
      if(col_varExist('HU'))  hu_column   => col_getColumn(column,headerIndex,'HU')
      if(col_varExist('TT'))  tt_column   => col_getColumn(column,headerIndex,'TT')
      if(col_varExist('P0'))  ps_column   => col_getColumn(column,headerIndex,'P0')
      if(col_varExist('TG'))  tg_column   => col_getColumn(column,headerIndex,'TG')
      if(col_varExist('LVIS'))vis_column  => col_getColumn(column,headerIndex,'LVIS')
      if(col_varExist('WGE')) gust_column => col_getColumn(column,headerIndex,'WGE')
     
      do jk = 1, gsv_getNumLev(statevector,'MM')
        if(gsv_varExist(statevector,'UU')) then
          jk2=jk+gsv_getOffsetFromVarName(statevector,'UU')
          uu_column(jk) =   dlw1*zgd(lonIndex  ,ila,jk2)  &
                          + dlw2*zgd(lonIndex+1,ila,jk2)  &
                          + dlw3*zgd(lonIndex  ,ila+1,jk2)  &
                          + dlw4*zgd(lonIndex+1,ila+1,jk2)
        end if
        if(gsv_varExist(statevector,'VV')) then
          jk2=jk+gsv_getOffsetFromVarName(statevector,'VV')
          vv_column(jk) =   dlw1*zgd(lonIndex  ,ila,jk2)  &
                          + dlw2*zgd(lonIndex+1,ila,jk2)  &
                          + dlw3*zgd(lonIndex  ,ila+1,jk2)  &
                          + dlw4*zgd(lonIndex+1,ila+1,jk2)
        end if
      end do
      do jk = 1, gsv_getNumLev(statevector,'TH')
        if(gsv_varExist(statevector,'HU')) then
          jk2=jk+gsv_getOffsetFromVarName(statevector,'HU')
          hu_column(jk) =   dlw1*zgd(lonIndex  ,ila,jk2)  &
                          + dlw2*zgd(lonIndex+1,ila,jk2)  &
                          + dlw3*zgd(lonIndex  ,ila+1,jk2)  &
                          + dlw4*zgd(lonIndex+1,ila+1,jk2)
        end if
        if(gsv_varExist(statevector,'TT')) then
          jk2=jk+gsv_getOffsetFromVarName(statevector,'TT')
          tt_column(jk) =   dlw1*zgd(lonIndex  ,ila,jk2)  &
                          + dlw2*zgd(lonIndex+1,ila,jk2)  &
                          + dlw3*zgd(lonIndex  ,ila+1,jk2)  &
                          + dlw4*zgd(lonIndex+1,ila+1,jk2)
        end if
        if(gsv_varExist(statevector,'LVIS')) then
          jk2=jk+gsv_getOffsetFromVarName(statevector,'LVIS')
          vis_column(jk) =   dlw1*zgd(lonIndex  ,ila,jk2)  &
                           + dlw2*zgd(lonIndex+1,ila,jk2)  &
                           + dlw3*zgd(lonIndex  ,ila+1,jk2)  &
                           + dlw4*zgd(lonIndex+1,ila+1,jk2)
        end if
      end do
      if(gsv_varExist(statevector,'P0')) then
        jk2=1+gsv_getOffsetFromVarName(statevector,'P0')
        ps_column(1) =   dlw1*zgd(lonIndex  ,ila,jk2)  &
                       + dlw2*zgd(lonIndex+1,ila,jk2)  &
                       + dlw3*zgd(lonIndex  ,ila+1,jk2)  &
                       + dlw4*zgd(lonIndex+1,ila+1,jk2)
      end if
      if(gsv_varExist(statevector,'WGE')) then
        jk2=1+gsv_getOffsetFromVarName(statevector,'WGE')
        gust_column(1) = dlw1*zgd(lonIndex  ,ila,jk2)  &
                       + dlw2*zgd(lonIndex+1,ila,jk2)  &
                       + dlw3*zgd(lonIndex  ,ila+1,jk2)  &
                       + dlw4*zgd(lonIndex+1,ila+1,jk2)
      end if
      if(gsv_varExist(statevector,'TG')) then
        jk2=1+gsv_getOffsetFromVarName(statevector,'TG')
        tg_column(1) =   dlw1*zgd(lonIndex  ,ila,jk2)  &
                       + dlw2*zgd(lonIndex+1,ila,jk2)  &
                       + dlw3*zgd(lonIndex  ,ila+1,jk2)  &
                       + dlw4*zgd(lonIndex+1,ila+1,jk2)
      end if
    end do

    deallocate(zgd)

  end subroutine s2c_bgcheck_bilin

  !--------------------------------------------------------------------------
  ! s2c_column_hbilin  
  !--------------------------------------------------------------------------
  subroutine s2c_column_hbilin(field,vlev,nlong,nlat,nlev,xlong,xlat, &
                               plong,plat,vprof,vlevout,nlevout)
    ! **Purpose:** 
    ! Horizontal bilinear interpolation from a 3D field to a profile at (plong,plat).
    ! Assumes vertical interpolation not needed or already done.
    !
    ! This version can be used with fields that are not part of the background state,
    ! such as climatologies.
    !
    ! This version does not depend in column_data and gridstatevector modules.
    !
    ! **Author:** Y. Rochon, Nov 2015 
    !
    implicit none

    ! arguments:
    integer, intent(in) :: nlong            ! number or longitudes
    integer, intent(in) :: nlat             ! number or latitudes
    integer, intent(in) :: nlev             ! number of vertical levels
    integer, intent(in) :: nlevout          ! number of target vertical levels
    real(8), intent(in) :: field(nlong,nlat,nlev) ! 3D field
    real(8), intent(in) :: vlev(nlev)       ! vertical levels of input field (in pressure)
    real(8), intent(in) :: xlong(nlong)     ! longitudes (radians)
    real(8), intent(in) :: xlat(nlat)       ! latitudes (radians)
    real(8), intent(in) :: plong            ! target longitude (radians)
    real(8), intent(in) :: plat             ! target latitude (radian)
    real(8), intent(in) :: vlevout(nlevout) ! target vertical levels (in pressure)
    real(8), intent(out) :: vprof(nlevout)  ! profile at (plong,plat)
    
    ! locals:
    real(8) :: lnvlev(nlev),lnvlevout(nlevout),plong2
    integer :: ilev,lonIndex,latIndex,i,j

    real(8) :: DLDX, DLDY, DLDP, DLW1, DLW2, DLW3, DLW4

    ! Find near lat/long grid points

    plong2 = plong
    if (plong2 < 0.0) plong2 = 2.D0*MPC_PI_R8 + plong2
    do lonIndex = 2, nlong
      if  (xlong(lonIndex-1) < xlong(lonIndex)) then
        if (plong2 >= xlong(lonIndex-1) .and. plong2 <= xlong(lonIndex)) exit
      else 
        ! Assumes this is a transition between 360 to 0 (if it exists). Skip over.
      end if
    end do
    lonIndex = lonIndex-1

    do latIndex = 2, nlat
      if (plat <= xlat(latIndex)) exit
    end do
    latIndex = latIndex-1

    ! Set lat/long interpolation weights

    DLDX = (plong - xlong(lonIndex))/(xlong(lonIndex+1)-xlong(lonIndex))
    DLDY = (plat - xlat(latIndex))/(xlat(latIndex+1)-xlat(latIndex))

    DLW1 = (1.d0-DLDX) * (1.d0-DLDY)
    DLW2 =       DLDX  * (1.d0-DLDY)
    DLW3 = (1.d0-DLDX) *       DLDY
    DLW4 =       DLDX  *       DLDY

    ! Set vertical interpolation weights (assumes pressure vertical coordinate)

    lnvlevout(:) = log(vlevout(:))    
    lnvlev(:) = log(vlev(:))    

    ilev = 1
    do i = 1, nlevout
      do j = ilev, nlev          
        if (lnvlevout(i) < lnvlev(j)) exit ! assumes lnvlevout and lnvlev increase with index
      end do
      ilev = j-1
      if (ilev < 1) then
        ilev = 1
      else if (ilev >= nlev) then
        ilev = nlev-1
      end if

      DLDP = (lnvlev(ilev+1)-lnvlevout(i))/(lnvlev(ilev+1)-lnvlev(ilev))
          
      vprof(i) = DLDP* (DLW1 * field(lonIndex,latIndex,ilev)      &
                      + DLW2 * field(lonIndex+1,latIndex,ilev)    &
                      + DLW3 * field(lonIndex,latIndex+1,ilev)    &
                      + DLW4 * field(lonIndex+1,latIndex+1,ilev)) &
        + (1.d0-DLDP)* (DLW1 * field(lonIndex,latIndex,ilev+1)    &
                      + DLW2 * field(lonIndex+1,latIndex,ilev+1)  &
                      + DLW3 * field(lonIndex,latIndex+1,ilev+1)  &
                      + DLW4 * field(lonIndex+1,latIndex+1,ilev+1))                               
    end do

  end subroutine s2c_column_hbilin   

end module stateToColumn_mod

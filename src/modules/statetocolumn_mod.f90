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
  ! MODULE stateToColumn (prefix='s2c' category='3. High-level transformations')
  !
  ! :Purpose: Non-linear, tangent-linear and adjoint versions of
  !           horizontal-temporal interpolation between a gridStateVector object
  !           and a columnData object.
  !
  use mathPhysConstants_mod
  use mpi, only : mpi_status_size ! this is the mpi library module
  use mpi_mod
  use mpivar_mod 
  use gridstatevector_mod
  use obsSpaceData_mod
  use columnData_mod
  use analysisgrid_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use obsTimeInterp_mod
  use tt2phi_mod
  use windRotation_mod
  use utilities_mod
  use gridVariableTransforms_mod
  use varNameList_mod
  use physicsFunctions_mod
  use timeCoord_mod
  use slantprofilelatlon_mod
  use tovs_nl_mod
  use codtyp_mod

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
    real(8), pointer          :: allLat(:,:,:,:) => null()         ! (headerUsed, proc, step, kIndex)
    real(8), pointer          :: allLon(:,:,:,:) => null()         ! (headerUsed, proc, step, kIndex)
    real(8), pointer          :: allLatRot(:,:,:,:,:) => null()    ! (subGrid, headerUsed, proc, step, kIndex)
    real(8), pointer          :: allLonRot(:,:,:,:,:) => null()    ! (subGrid, headerUsed, proc, step, kIndex)

    ! interpolation weights and lat/lon indices are accessed via the 'depotIndexBeg/End'
    integer, pointer          :: depotIndexBeg(:,:,:,:,:) => null()    ! (subGrid, headerUsed, proc, step, kIndex)
    integer, pointer          :: depotIndexEnd(:,:,:,:,:) => null()    ! (subGrid, headerUsed, proc, step, kIndex)
    real(8), allocatable      :: interpWeightDepot(:)                ! (depotIndex)
    integer, pointer          :: latIndexDepot(:)                    ! (depotIndex)
    integer, pointer          :: lonIndexDepot(:)                    ! (depotIndex)
  end type struct_interpInfo

  type(struct_interpInfo) :: interpInfo_tlad, interpInfo_nl

  real(8), pointer :: allLatOneLev(:,:)
  real(8), pointer :: allLonOneLev(:,:)

  character(len=20), parameter :: timeInterpType_tlad = 'LINEAR' ! hardcoded type of time interpolation for increment

  ! "special" values of the footprint radius
  real(4), parameter :: nearestNeighbourFootprint = -2.0
  real(4), parameter ::             lakeFootprint = -1.0
  real(4), parameter ::         bilinearFootprint =  0.0

  integer, external    :: get_max_rss
  logical, save :: slantPath_nl
  logical, save :: slantPath_tlad
  logical, save :: nmlAlreadyRead = .false.

contains 

  !---------------------------------------------------------
  ! latlonChecksAnlGrid
  !---------------------------------------------------------
  subroutine latlonChecksAnlGrid(obsSpaceData, moveObsAtPole)
    !
    ! :Purpose: Check the lat/lon of observations and modify if necessary
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
    write(*,*) 'latlonChecksAnlGrid: STARTING'
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
      ierr = utl_getPositionXY( hco_anl % EZscintID,  &
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
            write(*,*) 'latlonChecksAnlGrid: Moving OBS inside the GLOBAL ANALYSIS grid, ', headerIndex
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
            write(*,*) 'latlonChecksAnlGrid: OBS outside the GLOBAL ANALYSIS grid, but NOT moved, ', headerIndex
          end if

        else
          ! The observation is outside the domain
          ! In LAM Analysis mode we must discard this observation
          write(*,*) 'latlonChecksAnlGrid: Rejecting OBS outside the LAM ANALYSIS grid domain, ', headerIndex
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

    write(*,*) 'latlonChecksAnlGrid: END'

  end subroutine latlonChecksAnlGrid

  !---------------------------------------------------------
  ! s2c_setupInterpInfo
  !---------------------------------------------------------
  subroutine s2c_setupInterpInfo( interpInfo, obsSpaceData, stateVector,  &
                                  timeInterpType, rejectOutsideObs, &
                                  inputStateVectorType )
    ! :Purpose: Setup all of the information needed to quickly
    !           perform the horizontal interpolation to the observation
    !           locations.
    !
    implicit none

    ! arguments
    type(struct_interpInfo)    :: interpInfo
    type(struct_obs)           :: obsSpaceData
    type(struct_gsv), target   :: stateVector
    logical                    :: rejectOutsideObs
    character(len=*)           :: timeInterpType
    character(len=*)           :: inputStateVectorType

    ! locals
    type(struct_gsv)          :: stateVector_VarsLevs_1Step, stateVector_Tiles_allVar_1Step, stateVector_Tiles_1Step, stateVector_1Step
    type(struct_gsv), pointer :: stateVector_Tiles_ptr
    integer :: numHeader, numHeaderUsedMax, headerIndex, headerUsedIndex, bodyIndex, kIndex, kIndexCount, myKBeg
    integer :: numStep, stepIndex, fnom, fclos, nulnam, ierr
    integer :: bodyIndexBeg, bodyIndexEnd, procIndex, niP1, numGridptTotal, numHeaderUsed
    integer :: subGridIndex, subGridForInterp, numSubGridsForInterp
    real(8) :: latRot, lonRot, lat, lon
    real(4) :: lon_r4, lat_r4, lon_deg_r4, lat_deg_r4
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    real(4) :: footprintRadius_r4 ! (metres)
    integer, allocatable :: numGridpt(:), allNumHeaderUsed(:,:), headerIndexVec(:,:)
    real(8), allocatable :: lat_send_r8(:,:), lat_recv_r8(:,:), lon_send_r8(:,:), lon_recv_r8(:,:)
    real(4), allocatable :: footprintRadiusVec_r4(:), allFootprintRadius_r4(:,:,:)
    integer :: gdllfxy
    logical :: obsOutsideGrid
    character(len=4), pointer :: varNames(:)
    character(len=4)          :: varLevel, varName
    real(8), allocatable :: latColumn(:,:), lonColumn(:,:)
    real(8), allocatable :: latLev_T(:), lonLev_T(:), latLev_M(:), lonLev_M(:)
    real(4), pointer :: height3D_r4_ptr1(:,:,:), height3D_r4_ptr2(:,:,:), height3D_T_r4(:,:,:), height3D_M_r4(:,:,:)
    real(8), pointer :: height3D_r8_ptr1(:,:,:)
    logical :: thisProcIsAsender(mpi_nprocs)
    integer :: sendsizes(mpi_nprocs), recvsizes(mpi_nprocs), senddispls(mpi_nprocs), recvdispls(mpi_nprocs), allkBeg(mpi_nprocs)
    integer :: codeType, nlev_T, nlev_M, levIndex 
    integer :: maxkcount, numkToSend 
    integer :: firstHeaderIndexUsed(mpi_nprocs)
    logical :: doSlantPath, firstHeaderSlantPath 

    namelist /nams2c/ slantPath_nl, slantPath_tlad 

    write(*,*) 's2c_setupInterpInfo: STARTING'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    numStep = stateVector%numStep
    numHeader = obs_numheader(obsSpaceData)

    call oti_setup(interpInfo%oti, obsSpaceData, numStep, timeInterpType, flagObsOutside_opt=.true.)

    if ((stateVector%heightSfcPresent) .and. ( mpi_myid == 0)) then
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

    write(*,*) 's2c_setupInterpInfo: inputStateVectorType=', inputStateVectorType

    if ( .not. nmlAlreadyRead ) then
      nmlAlreadyRead = .true.

      ! default values
      slantPath_nl = .false.
      slantPath_tlad = .false.

      ! reading namelist variables
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam, nml=nams2c, iostat=ierr)
      if ( ierr /= 0 .and. mpi_myid == 0 ) write(*,*) 's2c_setupInterpInfomyWarning: nams2c is missing in the namelist. The default value will be taken.'
      if ( mpi_myid == 0 ) write(*, nml=nams2c)
      ierr = fclos(nulnam)
    end if

    doSlantPath = .false.
    if ( slantPath_nl   .and. inputStateVectorType == 'nl' ) doSlantPath = .true.
    if ( slantPath_tlad .and. inputStateVectorType /= 'nl' ) doSlantPath = .true.
    write(*,*) 's2c_setupInterpInfo: doSlantPath=', doSlantPath

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

    ! allocate arrays that will be returned
    allocate(interpInfo%allNumHeaderUsed(numStep,mpi_nprocs))
    allocate(interpInfo%depotIndexBeg(interpInfo%hco%numSubGrid,numHeaderUsedMax,mpi_nprocs,numStep,mykBeg:stateVector%mykEnd))
    allocate(interpInfo%depotIndexEnd(interpInfo%hco%numSubGrid,numHeaderUsedMax,mpi_nprocs,numStep,mykBeg:stateVector%mykEnd))
    allocate(interpInfo%allHeaderIndex(numHeaderUsedMax,numStep,mpi_nprocs))
    allocate(interpInfo%allLat(numHeaderUsedMax,mpi_nprocs,numStep,mykBeg:stateVector%mykEnd))
    allocate(interpInfo%allLon(numHeaderUsedMax,mpi_nprocs,numStep,mykBeg:stateVector%mykEnd))
    nullify(allLatOneLev)
    nullify(allLonOneLev)
    allocate(allLatOneLev(numHeaderUsedMax,mpi_nprocs))
    allocate(allLonOneLev(numHeaderUsedMax,mpi_nprocs))
    allocate(allFootprintRadius_r4(numHeaderUsedMax,numStep,mpi_nprocs))
    allocate(numGridpt(interpInfo%hco%numSubGrid))
    interpInfo%allHeaderIndex(:,:,:) = 0
    interpInfo%allLat(:,:,:,:) = 0.0d0
    interpInfo%allLon(:,:,:,:) = 0.0d0
    allFootprintRadius_r4(:,:,:) = bilinearFootprint
    interpInfo%allNumHeaderUsed(:,:) = allNumHeaderUsed(:,:)

    if ( interpInfo%hco%rotated ) then
      allocate(interpInfo%allLatRot(interpInfo%hco%numSubGrid,numHeaderUsedMax,mpi_nprocs,numStep,mykBeg:stateVector%mykEnd))
      allocate(interpInfo%allLonRot(interpInfo%hco%numSubGrid,numHeaderUsedMax,mpi_nprocs,numStep,mykBeg:stateVector%mykEnd))
      interpInfo%allLatRot(:,:,:,:,:) = 0.0d0
      interpInfo%allLonRot(:,:,:,:,:) = 0.0d0
    end if

    interpInfo%depotIndexBeg(:,:,:,:,:) = 0
    interpInfo%depotIndexEnd(:,:,:,:,:) = -1

    ! prepare for extracting the 3D height for slant-path calculation
    if ( doSlantPath .and. &
         stateVector%varExistList(vnl_varListIndex('Z_T')) .and. &
         stateVector%varExistList(vnl_varListIndex('Z_M')) ) then

      write(*,*) 's2c_setupInterpInfo: extracting 3D heights for slant-path for ', inputStateVectorType 

      if ( inputStateVectorType == 'nl' ) then
        nullify(varNames)
        call gsv_varNamesList(varNames, stateVector)
        call gsv_allocate( stateVector_VarsLevs_1Step, 1, &
                           stateVector%hco, stateVector%vco, &
                           mpi_local_opt=.true., mpi_distribution_opt='VarsLevs', &
                           dataKind_opt=4, varNames_opt=varNames )

        height3D_r4_ptr1 => gsv_getField3D_r4(stateVector)
        height3D_r4_ptr2 => gsv_getField3D_r4(stateVector_VarsLevs_1Step)
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

        height3D_r4_ptr1 => gsv_getField3D_r4(stateVector_Tiles_allVar_1Step,'Z_T')
        height3D_r4_ptr2 => gsv_getField3D_r4(stateVector_Tiles_1Step,'Z_T')
        height3D_r4_ptr2(:,:,:) = height3D_r4_ptr1(:,:,:)

        height3D_r4_ptr1 => gsv_getField3D_r4(stateVector_Tiles_allVar_1Step,'Z_M')
        height3D_r4_ptr2 => gsv_getField3D_r4(stateVector_Tiles_1Step,'Z_M')
        height3D_r4_ptr2(:,:,:) = height3D_r4_ptr1(:,:,:)

        call gsv_deallocate(stateVector_Tiles_allVar_1Step)

      else
        stateVector_Tiles_ptr => gvt_getStateVectorTrial('height')

        call gsv_allocate( stateVector_Tiles_1Step, 1, &
                           stateVector%hco, stateVector%vco, &
                           mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                           dataKind_opt=4, varNames_opt=(/'Z_M','Z_T'/) )

        height3D_r8_ptr1 => gsv_getField3D_r8(stateVector_Tiles_ptr,'Z_T')
        height3D_r4_ptr2 => gsv_getField3D_r4(stateVector_Tiles_1Step,'Z_T')
        height3D_r4_ptr2(:,:,:) = height3D_r8_ptr1(:,:,:)

        height3D_r8_ptr1 => gsv_getField3D_r8(stateVector_Tiles_ptr,'Z_M')
        height3D_r4_ptr2 => gsv_getField3D_r4(stateVector_Tiles_1Step,'Z_M')
        height3D_r4_ptr2(:,:,:) = height3D_r8_ptr1(:,:,:)

      end if ! inputStateVectorType 

      nlev_T = gsv_getNumLev(stateVector,'TH')
      nlev_M = gsv_getNumLev(stateVector,'MM')
      if ( mpi_myid == 0 ) then
        call gsv_allocate( stateVector_1Step, 1, &
                           stateVector%hco, stateVector%vco, &
                           mpi_local_opt=.false., &
                           dataKind_opt=4, varNames_opt=(/'Z_M','Z_T'/) )

        height3D_T_r4 => gsv_getField3D_r4(stateVector_1Step,'Z_T')
        height3D_M_r4 => gsv_getField3D_r4(stateVector_1Step,'Z_M')

      else
        allocate(height3D_T_r4(stateVector%ni,stateVector%nj,nlev_T))
        allocate(height3D_M_r4(stateVector%ni,stateVector%nj,nlev_M))
      end if

      ! now bring all the heights to processor 0
      call gsv_transposeTilesToStep(stateVector_1Step, stateVector_Tiles_1Step, 1)

      ! broadcast 3D height field (single precision) to all the processors
      call rpn_comm_bcast(height3D_T_r4, size(height3D_T_r4), 'MPI_REAL4', 0, 'GRID', ierr)
      call rpn_comm_bcast(height3D_M_r4, size(height3D_M_r4), 'MPI_REAL4', 0, 'GRID', ierr)

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

      header_loop2: do headerIndex = 1, numHeader

        ! if obs inside window, but zero weight for current stepIndex then skip it
        if ( oti_getTimeInterpWeight(interpInfo%oti, headerIndex, stepIndex) == 0.0d0 ) cycle header_loop2

        numHeaderUsed = numHeaderUsed + 1
        headerIndexVec(numHeaderUsed,stepIndex) = headerIndex

        footprintRadiusVec_r4(numHeaderUsed) = s2c_getFootprintRadius(obsSpaceData, headerIndex)

      end do header_loop2

      call rpn_comm_allgather(footprintRadiusVec_r4,                numHeaderUsedMax, 'MPI_REAL4', &
                              allFootprintRadius_r4(:,stepIndex,:), numHeaderUsedMax, 'MPI_REAL4', &
                              'GRID', ierr)

      allocate(latColumn(numHeaderUsedMax,allkBeg(1):stateVector%nk))
      allocate(lonColumn(numHeaderUsedMax,allkBeg(1):stateVector%nk))
      latColumn(:,:) = 0.0d0
      lonColumn(:,:) = 0.0d0

      if ( doSlantPath .and. &
           stateVector%varExistList(vnl_varListIndex('Z_T')) .and. &
           stateVector%varExistList(vnl_varListIndex('Z_M')) ) then

        allocate(latLev_T(nlev_T))
        allocate(lonLev_T(nlev_T))
        allocate(latLev_M(nlev_M))
        allocate(lonLev_M(nlev_M))
        latLev_T(:) = 0.0d0
        lonLev_T(:) = 0.0d0
        latLev_M(:) = 0.0d0
        lonLev_M(:) = 0.0d0

        firstHeaderSlantPath  = .true.
        header_loop3: do headerUsedIndex = 1, numHeaderUsed
          headerIndex = headerIndexVec(headerUsedIndex,stepIndex)

          !- Get LatLon of observation location
          lat_r4 = real(obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex), 4)
          lon_r4 = real(obs_headElem_r(obsSpaceData, OBS_LON, headerIndex), 4)
          if (lon_r4 <  0.0          ) lon_r4 = lon_r4 + 2.0*MPC_PI_R4
          if (lon_r4 >= 2.0*MPC_PI_R4) lon_r4 = lon_r4 - 2.0*MPC_PI_R4

          codeType = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)

          if ( tvs_isIdBurpTovs(codeType) ) then
            if ( firstHeaderSlantPath ) then
              write(*,'(a,i3,a,i8)') 's2c_setupInterpInfo: start slant-path for TOVS. stepIndex=',stepIndex,' and numHeaderUsed=',numHeaderUsed
              firstHeaderSlantPath = .false.
            end if

            ! calculate lat/lon along the line of sight
            call tmg_start(199,'slp_calcLatLonTovs')
            call slp_calcLatLonTovs( obsSpaceData, stateVector%hco, headerIndex, & ! IN
                                     height3D_T_r4, height3D_M_r4,               & ! IN
                                     latLev_T, lonLev_T,                         & ! OUT
                                     latLev_M, lonLev_M )                          ! OUT
            call tmg_stop(199)

          else
            latLev_T(:) = real(lat_r4,8)
            lonLev_T(:) = real(lon_r4,8)
            latLev_M(:) = real(lat_r4,8)
            lonLev_M(:) = real(lon_r4,8)

          end if !tvs_isIdBurpTovs

          ! check if the slanted lat/lon is inside the domain
          call tmg_start(198,'latlonChecks')
          call latlonChecks ( obsSpaceData, stateVector%hco, & ! IN
                              headerIndex, rejectOutsideObs, & ! IN
                              latLev_T, lonLev_T,            & ! IN/OUT
                              latLev_M, lonLev_M )             ! IN/OUT 
          call tmg_stop(198)

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
              latColumn(headerUsedIndex,kIndex) = real(lat_r4,8)
              lonColumn(headerUsedIndex,kIndex) = real(lon_r4,8)
            end if

          end do

        end do header_loop3

        ! MPI communication for the slant-path lat/lon
        ! all tasks are senders
        do procIndex = 1, mpi_nprocs
          thisProcIsAsender(procIndex) = .true.
        end do

        maxkCount = maxval(stateVector%allkCount(:) + stateVector%allkBeg(:) - allkBeg(:))
        numkToSend = min(mpi_nprocs,stateVector%nk)

        allocate(lat_recv_r8(numHeaderUsedMax,mpi_nprocs))
        lat_recv_r8(:,:) = 0.0d0
        allocate(lat_send_r8(numHeaderUsedMax,mpi_nprocs))
        lat_send_r8(:,:) = 0.0d0
        allocate(lon_recv_r8(numHeaderUsedMax,mpi_nprocs))
        lon_recv_r8(:,:) = 0.0d0
        allocate(lon_send_r8(numHeaderUsedMax,mpi_nprocs))
        lon_send_r8(:,:) = 0.0d0

        ! only send the data from tasks with data, same amount to all
        sendsizes(:) = 0
        do procIndex = 1, numkToSend
          sendsizes(procIndex) = numHeaderUsedMax
        end do
        senddispls(1) = 0
        do procIndex = 2, mpi_nprocs
          senddispls(procIndex) = senddispls(procIndex-1) + sendsizes(procIndex-1)
        end do

        ! all tasks recv only from those with data
        recvsizes(:) = 0
        if ( (1+mpi_myid) <= numkToSend ) then
          do procIndex = 1, mpi_nprocs
            if ( thisProcIsAsender(procIndex) ) then
              recvsizes(procIndex) = numHeaderUsedMax
            end if
          end do
        end if
        recvdispls(1) = 0
        do procIndex = 2, mpi_nprocs
          recvdispls(procIndex) = recvdispls(procIndex-1) + recvsizes(procIndex-1)
        end do

        ! loop to send (at most) 1 level to (at most) all other mpi tasks
        do kIndexCount = 1, maxkCount
          do procIndex = 1, mpi_nprocs
            ! compute kIndex value being sent
            kIndex = kIndexCount + allkBeg(procIndex) - 1
            if ( kIndex <= stateVector%allkEnd(procIndex) ) then
              if( procIndex > numkToSend ) then
                write(*,*) 'procIndex, numkToSend = ', procIndex, numkToSend
                call utl_abort('ERROR: with numkToSend?')
              end if

              lat_send_r8(:,procIndex) = latColumn(:,kIndex)
              lon_send_r8(:,procIndex) = lonColumn(:,kIndex)
            end if
          end do

          call mpi_alltoallv(lat_send_r8, sendsizes, senddispls, mpi_datyp_real8,  &
                             lat_recv_r8, recvsizes, recvdispls, mpi_datyp_real8, mpi_comm_grid, ierr)
          call mpi_alltoallv(lon_send_r8, sendsizes, senddispls, mpi_datyp_real8,  &
                             lon_recv_r8, recvsizes, recvdispls, mpi_datyp_real8, mpi_comm_grid, ierr)

          do procIndex = 1, mpi_nprocs
            ! all tasks copy the received step data into correct slot
            kIndex = kIndexCount + mykBeg - 1
            if ( kIndex <= stateVector%mykEnd ) then
              interpInfo%allLat(:,procIndex,stepIndex,kIndex) = lat_recv_r8(:,procIndex)
              interpInfo%allLon(:,procIndex,stepIndex,kIndex) = lon_recv_r8(:,procIndex)
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
          call tmg_start(198,'latlonChecks')
          call latlonChecks ( obsSpaceData, stateVector%hco, & ! IN
                              headerIndex, rejectOutsideObs, & ! IN
                              latLev_T, lonLev_T,            & ! IN/OUT
                              latLev_M, lonLev_M )             ! IN/OUT 
          call tmg_stop(198)

          latColumn(headerUsedIndex,allkBeg(1)) = latLev_T(1)
          lonColumn(headerUsedIndex,allkBeg(1)) = lonLev_T(1)
        end do

        ! gather geographical lat, lon positions of observations from all processors
        call rpn_comm_allgather(latColumn(:,allkBeg(1)), numHeaderUsedMax, 'MPI_REAL8', &
                                allLatOneLev(:,:), numHeaderUsedMax, 'MPI_REAL8', 'GRID', ierr)
        call rpn_comm_allgather(lonColumn(:,allkBeg(1)), numHeaderUsedMax, 'MPI_REAL8', &
                                allLonOneLev(:,:), numHeaderUsedMax, 'MPI_REAL8', 'GRID', ierr)

        k_loop: do kIndex = mykBeg, statevector%mykEnd
          interpInfo%allLat(:,:,stepIndex,kIndex) = allLatOneLev(:,:)
          interpInfo%allLon(:,:,stepIndex,kIndex) = allLonOneLev(:,:)
        end do k_loop

        deallocate(latLev_T)
        deallocate(lonLev_T)
        deallocate(latLev_M)
        deallocate(lonLev_M)

      end if ! doSlantPath 

      deallocate(lonColumn)
      deallocate(latColumn)

    end do step_loop2

    if ( doSlantPath .and. &
         stateVector%varExistList(vnl_varListIndex('Z_T')) .and. &
         stateVector%varExistList(vnl_varListIndex('Z_M')) ) then
      if ( mpi_myid == 0 ) then
        call gsv_deallocate(stateVector_1Step)
      else
        deallocate(height3D_T_r4)
        deallocate(height3D_M_r4)
      end if
    end if
    deallocate(footprintRadiusVec_r4)

    write(*,*) 's2c_setupInterpInfo: latlonChecks and lat/lon MPI comm finished.'

    ! gather the headerIndexVec arrays onto all processors
    call rpn_comm_allgather(headerIndexVec,            numHeaderUsedMax*numStep, 'MPI_INTEGER', &
                            interpInfo%allHeaderIndex, numHeaderUsedMax*numStep, 'MPI_INTEGER', &
                            'GRID',ierr)

    numGridptTotal = 0
    call tmg_start(173,'S2CNL_SETUPROTLL')
    do kIndex = mykBeg, statevector%mykEnd
      do stepIndex = 1, numStep
        do procIndex = 1, mpi_nprocs
          do headerIndex = 1, allNumHeaderUsed(stepIndex,procIndex)

            ! Compute the rotated lat/lon, needed for the winds

            lat_deg_r4 = real(interpInfo%allLat(headerIndex, procIndex, stepIndex, kIndex) *  &
                         MPC_DEGREES_PER_RADIAN_R8)
            lon_deg_r4 = real(interpInfo%allLon(headerIndex, procIndex, stepIndex, kIndex) *  &
                         MPC_DEGREES_PER_RADIAN_R8)
            ierr = utl_getPositionXY( stateVector%hco%EZscintID,   &
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
                lat = interpInfo%allLat(headerIndex, procIndex, stepIndex, kIndex)
                lon = interpInfo%allLon(headerIndex, procIndex, stepIndex, kIndex)
                call uvr_RotateLatLon( interpInfo%uvr,   & ! INOUT
                                       subGridIndex,     & ! IN
                                       latRot, lonRot,   & ! OUT (radians)
                                       lat, lon,         & ! IN  (radians)
                                       'ToLatLonRot')      ! IN
                interpInfo%allLatRot(subGridIndex, headerIndex, procIndex, stepIndex, kIndex) = latRot
                interpInfo%allLonRot(subGridIndex, headerIndex, procIndex, stepIndex, kIndex) = lonRot
              end if

            end do ! subGridForInterp

            ! Count total number of grid points for allocating interp depot

            footprintRadius_r4 = allFootprintRadius_r4(headerIndex, stepIndex, procIndex)

            call s2c_setupHorizInterp(footprintRadius_r4, interpInfo, obsSpaceData, &
                                      stateVector, headerIndex, kIndex, stepIndex, &
                                      procIndex, numGridpt)

            if ( (subGridIndex == 1) .or. (subGridIndex == 2) ) then
              ! indices for only 1 subgrid, other will have zeros
              interpInfo%depotIndexBeg(subGridIndex, headerIndex, procIndex, stepIndex, kIndex) = numGridptTotal + 1
              numGridptTotal = numGridptTotal + numGridpt(subGridIndex)
              interpInfo%depotIndexEnd(subGridIndex, headerIndex, procIndex, stepIndex, kIndex) = numGridptTotal
            else
              ! locations on both subGrids will be averaged
              interpInfo%depotIndexBeg(1, headerIndex, procIndex, stepIndex, kIndex) = numGridptTotal + 1
              numGridptTotal = numGridptTotal + numGridpt(1)
              interpInfo%depotIndexEnd(1, headerIndex, procIndex, stepIndex, kIndex) = numGridptTotal

              interpInfo%depotIndexBeg(2, headerIndex, procIndex, stepIndex, kIndex) = numGridptTotal + 1
              numGridptTotal = numGridptTotal + numGridpt(2)
              interpInfo%depotIndexEnd(2, headerIndex, procIndex, stepIndex, kIndex) = numGridptTotal
            end if

          end do ! headerIndex
        end do ! kIndex
      end do ! stepIndex
    end do ! procIndex
    call tmg_stop(173)

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! now that we know the size, allocate main arrays for storing interpolation information
    write(*,*) 's2c_setupInterpInfo: numGridptTotal = ', numGridptTotal
    allocate( interpInfo%latIndexDepot(numGridptTotal) )
    allocate( interpInfo%lonIndexDepot(numGridptTotal) )
    allocate( interpInfo%interpWeightDepot(numGridptTotal) )

    call tmg_start(175,'S2CNL_SETUPWEIGHTS')
    !$OMP PARALLEL DO PRIVATE (procIndex, stepIndex, kIndex, headerIndex, footprintRadius_r4, numGridpt)
    do procIndex = 1, mpi_nprocs
      do stepIndex = 1, numStep
        do kIndex = mykBeg, statevector%mykEnd

          do headerIndex = 1, allNumHeaderUsed(stepIndex,procIndex)

            footprintRadius_r4 = allFootprintRadius_r4(headerIndex, stepIndex, procIndex)

            call s2c_setupHorizInterp(footprintRadius_r4, interpInfo, obsSpaceData, &
                                      stateVector, headerIndex, kIndex, stepIndex, &
                                      procIndex, numGridpt)

          end do ! headerIndex

        end do ! kIndex
      end do ! stepIndex
    end do ! procIndex
    !$OMP END PARALLEL DO
    call tmg_stop(175)

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
  subroutine s2c_tl( statevector, column, columng, obsSpaceData )
    !
    ! :Purpose: Tangent linear version of the horizontal
    !           interpolation, used for the increment (or perturbations).
    !
    implicit none

    ! arguments
    type(struct_gsv)           :: stateVector
    type(struct_obs)           :: obsSpaceData
    type(struct_columnData)    :: column
    type(struct_columnData)    :: columng

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
    logical              :: rejectOutsideObs
    character(len=4)     :: varName
    character(len=4), pointer :: varNames(:)

    if ( mpi_myid == 0 ) write(*,*) 's2c_tl: Horizontal interpolation StateVector --> ColumnData'
    call tmg_start(167,'S2C_TL')

    call tmg_start(160,'S2C_BARR')
    call rpn_comm_barrier('GRID',ierr)
    call tmg_stop(160)

    if ( .not. stateVector%allocated ) then 
      call utl_abort('s2c_tl: stateVector must be allocated')
    end if

    ! check the column and statevector have same nk/varNameList
    call checkColumnStatevectorMatch(column,statevector)

    ! calculate delP_T/delP_M on the grid
    if ( statevector%varExistList(vnl_varListIndex('P_T')) .and. &
         statevector%varExistList(vnl_varListIndex('P_M')) ) then
      call gvt_transform( statevector, & ! INOUT
                          'PsfcToP_tl')  ! IN
    end if

    ! calculate del Z_T/Z_M on the grid
    if ( statevector%varExistList(vnl_varListIndex('Z_T')) .and. &
         statevector%varExistList(vnl_varListIndex('Z_M')) ) then
      call gvt_transform( statevector, & ! INOUT
                          'TTHUtoHeight_tl') ! IN
    end if

    nullify(varNames)
    call gsv_varNamesList(varNames, statevector)
    call gsv_allocate( statevector_VarsLevs, statevector%numstep, &
                       statevector%hco, statevector%vco,          &
                       mpi_local_opt=.true., mpi_distribution_opt='VarsLevs', &
                       varNames_opt=varNames )
    deallocate(varNames)
    call gsv_transposeTilesToVarsLevs( statevector, statevector_VarsLevs )

    numStep = stateVector_VarsLevs%numStep
    numHeader = obs_numheader(obsSpaceData)
    call rpn_comm_allreduce(numHeader, numHeaderMax, 1,  &
                            'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)

    if ( .not. interpInfo_tlad%initialized ) then
      rejectOutsideObs = .false.
      call s2c_setupInterpInfo( interpInfo_tlad, obsSpaceData, stateVector_VarsLevs,  &
                                timeInterpType_tlad,  rejectOutsideObs, &
                                inputStateVectorType='tl' )
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
                                   interpInfo_tlad, kIndex, stepIndex, procIndex )
              else if ( varName == 'VV' ) then
                call myezuvint_tl( cols_hint(1:yourNumHeader,stepIndex,procIndex), 'VV',  &
                                   ptr3d_UV(:,:,stepIndex), ptr4d(:,:,kIndex,stepIndex),  &
                                   interpInfo_tlad, kIndex, stepIndex, procIndex )
              else
                call myezsint( cols_hint(1:yourNumHeader,stepIndex,procIndex), varName,  &
                               ptr4d(:,:,kIndex,stepIndex), interpInfo_tlad, kIndex, stepIndex, procIndex )
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
    !
    ! :Purpose: Adjoint version of the horizontal interpolation,
    !           used for the cost function gradient with respect to the increment.
    !
    implicit none

    ! arguments
    type(struct_gsv)           :: stateVector
    type(struct_obs)           :: obsSpaceData
    type(struct_columnData)    :: column
    type(struct_columnData)    :: columng

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
    logical              :: rejectOutsideObs
    character(len=4), pointer :: varNames(:)

    if(mpi_myid == 0) write(*,*) 's2c_ad: Adjoint of horizontal interpolation StateVector --> ColumnData'
    call tmg_start(168,'S2C_AD')

    call tmg_start(160,'S2C_BARR')
    call rpn_comm_barrier('GRID',ierr)
    call tmg_stop(160)

    if ( .not. stateVector%allocated ) then 
      call utl_abort('s2c_ad: stateVector must be allocated')
    end if

    nullify(varNames)
    call gsv_varNamesList(varNames, statevector)
    call gsv_allocate( statevector_VarsLevs, statevector%numstep, &
                       statevector%hco, statevector%vco,          &
                       mpi_local_opt=.true., mpi_distribution_opt='VarsLevs', &
                       varNames_opt=varNames )
    deallocate(varNames)
    call gsv_zero( statevector_VarsLevs )

    if ( .not. interpInfo_tlad%initialized ) then
      rejectOutsideObs = .false.
      call s2c_setupInterpInfo( interpInfo_tlad, obsSpaceData, stateVector_VarsLevs,  &
                                timeInterpType_tlad, rejectOutsideObs, &
                                inputStateVectorType='ad' )
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
            do headerUsedIndex = 1, interpInfo_tlad%allNumHeaderUsed(stepIndex, procIndex)

              headerIndex = interpInfo_tlad%allHeaderIndex(headerUsedIndex,stepIndex,procIndex)
              weight = oti_getTimeInterpWeightMpiGlobal(interpInfo_tlad%oti,  &
                                                        headerIndex,stepIndex,procIndex)

              cols_hint(headerUsedIndex,stepIndex,procIndex) =  &
                   weight * cols_recv(headerIndex,procIndex)

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
                                   interpInfo_tlad, kIndex, stepIndex, procIndex )
              else if ( varName == 'VV' ) then
                call myezuvint_ad( cols_hint(1:yourNumHeader,stepIndex,procIndex), 'VV',  &
                                   ptr3d_UV(:,:,stepIndex), ptr4d(:,:,kIndex,stepIndex),  &
                                   interpInfo_tlad, kIndex, stepIndex, procIndex )
              else
                call myezsintad( cols_hint(1:yourNumHeader,stepIndex,procIndex), varName,  &
                                 ptr4d(:,:,kIndex,stepIndex), interpInfo_tlad, kIndex, stepIndex,  &
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

    if ( statevector%varExistList(vnl_varListIndex('Z_T')) .and. &
         statevector%varExistList(vnl_varListIndex('Z_M')) ) then
      call gvt_transform( statevector, & ! INOUT
                          'TTHUtoHeight_ad') ! IN
    end if

    ! Adjoint of calculate delP_T/delP_M on the grid
    if ( statevector%varExistList(vnl_varListIndex('P_T')) .and. &
         statevector%varExistList(vnl_varListIndex('P_M')) ) then
      call gvt_transform( statevector, & ! INOUT
                          'PsfcToP_ad')  ! IN
    end if

    call gsv_deallocate( statevector_VarsLevs )

    call tmg_stop(168)

  end subroutine s2c_ad

  !---------------------------------------------------------
  ! s2c_nl
  !---------------------------------------------------------
  subroutine s2c_nl( stateVector, obsSpaceData, column, timeInterpType, varName_opt, &
                     dealloc_opt, moveObsAtPole_opt )
    ! :Purpose: Non-linear version of the horizontal interpolation,
    !           used for a full field (usually the background state when computing
    !           the innovation vector).
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
    type(struct_gsv) :: stateVector_VarsLevs 
    integer :: kIndex, kIndex2, kCount, stepIndex, numStep, mykEndExtended
    integer :: headerIndex, numHeader, numHeaderMax, yourNumHeader
    integer :: procIndex, nsize, ierr, headerUsedIndex
    integer :: kIndexHeightSfc
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
    logical              :: dealloc, moveObsAtPole, rejectOutsideObs
    character(len=4), pointer :: varNames(:)

    call tmg_start(169,'S2C_NL')

    write(*,*) 's2c_nl: STARTING'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if ( .not. stateVector%allocated ) then 
      call utl_abort('s2c_nl: stateVector must be allocated')
    end if

    if ( stateVector%mpi_distribution /= 'Tiles' ) then 
      call utl_abort('s2c_nl: stateVector must by Tiles distributed')
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

    ! check the column and statevector have same nk/varNameList
    call checkColumnStatevectorMatch(column,statevector)

    call tmg_start(171,'s2c_nl_calcPandZ')
    ! calculate P_T/P_M on the grid
    if ( statevector%varExistList(vnl_varListIndex('P_T')) .and. &
         statevector%varExistList(vnl_varListIndex('P_M')) ) then
      call gvt_transform( stateVector, & ! INOUT
                          'PsfcToP_nl')  ! IN
    end if

    ! calculate Z_T/Z_M on the grid
    if ( statevector%varExistList(vnl_varListIndex('Z_T')) .and. &
         statevector%varExistList(vnl_varListIndex('Z_M')) ) then
      call gvt_transform( stateVector, & ! INOUT
                          'TTHUtoHeight_nl') ! IN
    end if
    call tmg_stop(171)

    nullify(varNames)
    call gsv_varNamesList(varNames, statevector)
    call gsv_allocate( statevector_VarsLevs, stateVector%numstep, &
                       stateVector%hco, stateVector%vco, mpi_local_opt=.true., &
                       mpi_distribution_opt='VarsLevs', dataKind_opt=4, &
                       allocHeightSfc_opt=.true., varNames_opt=varNames )
    deallocate(varNames)
    call gsv_transposeTilesToVarsLevs( stateVector, stateVector_VarsLevs )

    numStep = stateVector_VarsLevs%numStep
    numHeader = obs_numheader(obsSpaceData)
    call rpn_comm_allreduce(numHeader, numHeaderMax, 1,  &
                            'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)

    if ( .not. interpInfo_nl%initialized ) then
      call tmg_start(165,'S2CNL_SETUPS')
      ! also reject obs outside (LAM) domain and optionally move obs near 
      ! numerical pole to first/last analysis grid latitude
      call latlonChecksAnlGrid( obsSpaceData, moveObsAtPole )

      ! Do not reject obs for global domain
      rejectOutsideObs = .not. stateVector_VarsLevs%hco%global
      write(*,*) 's2c_nl: rejectOutsideObs = ', rejectOutsideObs

      ! compute and collect all obs grids onto all mpi tasks
      call s2c_setupInterpInfo( interpInfo_nl, obsSpaceData, stateVector_VarsLevs,  &
                                timeInterpType, rejectOutsideObs, &
                                inputStateVectorType='nl' )
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

    ptr4d_r4    => gsv_getField_r4(stateVector_VarsLevs)

    allocate(field2d(stateVector_VarsLevs%ni,stateVector_VarsLevs%nj))
    allocate(field2d_UV(stateVector_VarsLevs%ni,stateVector_VarsLevs%nj))

    mykEndExtended = stateVector_VarsLevs%mykBeg + maxval(stateVector_VarsLevs%allkCount(:)) - 1

    kCount = 0
    k_loop: do kIndex = stateVector_VarsLevs%mykBeg, mykEndExtended
      kCount = kCount + 1

      if ( kIndex <= stateVector_VarsLevs%mykEnd ) then
        varName = gsv_getVarNameFromK(stateVector_VarsLevs,kIndex)

        if ( varName == 'UU' .or. varName == 'VV' ) then
          ptr3d_UV_r4 => gsv_getFieldUV_r4(stateVector_VarsLevs,kIndex)
        end if

        call tmg_start(166,'S2CNL_HINTERP')
        step_loop: do stepIndex = 1, numStep
          if ( maxval(interpInfo_nl%allNumHeaderUsed(stepIndex,:)) == 0 ) cycle step_loop

          ! copy over field
          field2d(:,:) = real(ptr4d_r4(:,:,kIndex,stepIndex),8)
          if ( varName == 'UU' .or. varName == 'VV' ) then
            field2d_UV(:,:) = real(ptr3d_UV_r4(:,:,stepIndex),8)
          end if

          ! interpolate to the columns destined for all procs for all steps and one lev/var
          !$OMP PARALLEL DO PRIVATE (procIndex, yourNumHeader)
          do procIndex = 1, mpi_nprocs
            yourNumHeader = interpInfo_nl%allNumHeaderUsed(stepIndex,procIndex)
            if ( yourNumHeader > 0 ) then
              if ( varName == 'UU' ) then
                call myezuvint_nl( cols_hint(1:yourNumHeader,stepIndex,procIndex), 'UU',  &
                                   field2d, field2d_UV, interpInfo_nl, kindex, stepIndex, procIndex )
              else if ( varName == 'VV' ) then
                call myezuvint_nl( cols_hint(1:yourNumHeader,stepIndex,procIndex), 'VV',  &
                                   field2d_UV, field2d, interpInfo_nl, kindex, stepIndex, procIndex )
              else
                call myezsint( cols_hint(1:yourNumHeader,stepIndex,procIndex), varName,  &
                               field2d, interpInfo_nl, kindex, stepIndex, procIndex )
              end if
            end if
          end do
          !$OMP END PARALLEL DO

        end do step_loop
        call tmg_stop(166)

        ! interpolate in time to the columns destined for all procs and one level/variable
        do procIndex = 1, mpi_nprocs
          cols_send_1proc(:) = 0.0d0
          do stepIndex = 1, numStep
            !$OMP PARALLEL DO PRIVATE (headerIndex, headerUsedIndex, weight)
            do headerUsedIndex = 1, interpInfo_nl%allNumHeaderUsed(stepIndex, procIndex)
              headerIndex = interpInfo_nl%allHeaderIndex(headerUsedIndex,stepIndex,procIndex)
              weight = oti_getTimeInterpWeightMpiGlobal(interpInfo_nl%oti,  &
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
        cols_send(:,:) = 0.0d0

      end if ! if kIndex <= mykEnd

      call tmg_start(160,'S2C_BARR')
      call rpn_comm_barrier('GRID',ierr)
      call tmg_stop(160)

      ! mpi communication: alltoall for one level/variable
      call tmg_start(172,'s2c_nl_allToAll')
      nsize = numHeaderMax
      if(mpi_nprocs > 1) then
        call rpn_comm_alltoall(cols_send, nsize, 'MPI_REAL8',  &
                               cols_recv, nsize, 'MPI_REAL8', 'GRID', ierr)
      else
        cols_recv(:,1) = cols_send(:,1)
      end if
      call tmg_stop(172)

      ! reorganize ensemble of distributed columns
      !$OMP PARALLEL DO PRIVATE (procIndex, kIndex2, headerIndex)
      proc_loop: do procIndex = 1, mpi_nprocs
        kIndex2 = stateVector_VarsLevs%allkBeg(procIndex) + kCount - 1
        if ( kIndex2 > stateVector_VarsLevs%allkEnd(procIndex) ) cycle proc_loop
        do headerIndex = 1, numHeader
          allCols_ptr(kIndex2,headerIndex) = cols_recv(headerIndex,procIndex)
        end do
      end do proc_loop
      !$OMP END PARALLEL DO

    end do k_loop

    ! Interpolate surface height separately, only exists on mpi task 0
    HeightSfcPresent: if ( stateVector_VarsLevs%HeightSfcPresent ) then

      if ( mpi_myid == 0 ) then
        varName = 'GZ'
        kIndexHeightSfc = 0    
        step_loop_height: do stepIndex = 1, numStep

          if ( maxval(interpInfo_nl%allNumHeaderUsed(stepIndex,:)) == 0 ) cycle step_loop_height

          ! interpolate to the columns destined for all procs for all steps and one lev/var
          !$OMP PARALLEL DO PRIVATE (procIndex, yourNumHeader, ptr2d_r8)
          do procIndex = 1, mpi_nprocs
            yourNumHeader = interpInfo_nl%allNumHeaderUsed(stepIndex,procIndex)
            if ( yourNumHeader > 0 ) then
              ptr2d_r8 => gsv_getHeightSfc(stateVector_VarsLevs)
              call myezsint( cols_hint(1:yourNumHeader,stepIndex,procIndex), varName,  &
                             ptr2d_r8(:,:), interpInfo_nl, kIndexHeightSfc, stepIndex, procIndex )

            end if
          end do
          !$OMP END PARALLEL DO

        end do step_loop_height

        ! interpolate in time to the columns destined for all procs and one level/variable
        do procIndex = 1, mpi_nprocs
          cols_send(:,procIndex) = 0.0d0
          do stepIndex = 1, numStep
            !$OMP PARALLEL DO PRIVATE (headerIndex, headerUsedIndex)
            do headerUsedIndex = 1, interpInfo_nl%allNumHeaderUsed(stepIndex, procIndex)
              headerIndex = interpInfo_nl%allHeaderIndex(headerUsedIndex,stepIndex,procIndex)
              ! just copy, since surface height same for all time steps
              cols_send(headerIndex,procIndex) = cols_hint(headerUsedIndex,stepIndex,procIndex)
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
        call col_setHeightSfc(column, headerIndex, cols_recv(headerIndex,1))
      end do

    end if HeightSfcPresent

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
    if(col_varExist(column,'HU')) then
      do headerIndex = 1, numHeader
        column_ptr => col_getColumn(column,headerIndex,'HU')
        column_ptr(:) = max(column_ptr(:),col_rhumin)
      end do
    end if

    call gsv_deallocate( statevector_VarsLevs )

    write(*,*) 's2c_nl: FINISHED'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    call tmg_stop(169)

  end subroutine s2c_nl

  ! -------------------------------------------------
  ! myezsint: Scalar field horizontal interpolation
  ! -------------------------------------------------
  subroutine myezsint( column_out, varName, field_in, interpInfo, kIndex, stepIndex, procIndex )
    !
    ! :Purpose: Scalar horizontal interpolation, replaces the
    !           ezsint routine from rmnlib.
    !
    implicit none

    ! arguments
    real(8)                 :: column_out(:)
    character(len=*)        :: varName
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
             interpInfo%depotIndexBeg(subGridIndex, headerIndex, procIndex, stepIndex, kIndex), &
             interpInfo%depotIndexEnd(subGridIndex, headerIndex, procIndex, stepIndex, kIndex)

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
  subroutine myezsintad( column_in, varName, field_out, interpInfo, kIndex, stepIndex, procIndex )
    !
    ! :Purpose: Adjoint of the scalar horizontal interpolation.
    !
    implicit none

    ! Arguments:
    real(8)                 :: column_in(:)
    character(len=*)        :: varName
    real(8)                 :: field_out(:,:)
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
             interpInfo%depotIndexBeg(subGridIndex, headerIndex, procIndex, stepIndex, kIndex), &
             interpInfo%depotIndexEnd(subGridIndex, headerIndex, procIndex, stepIndex, kIndex)

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
                           interpInfo, kIndex, stepIndex, procIndex )
    !
    ! :Purpose: Vector horizontal interpolation, replaces the
    !           ezuvint routine from rmnlib.
    !
    implicit none

    ! arguments
    real(8)                 :: column_out(:)
    character(len=*)        :: varName
    real(8)                 :: fieldUU_in(:,:)
    real(8)                 :: fieldVV_in(:,:)
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

        indexBeg = interpInfo%depotIndexBeg(subGridIndex, headerIndex, procIndex, stepIndex, kIndex)
        indexEnd = interpInfo%depotIndexEnd(subGridIndex, headerIndex, procIndex, stepIndex, kIndex)

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
          lat = interpInfo%allLat(headerIndex, procIndex, stepIndex, kIndex)
          lon = interpInfo%allLon(headerIndex, procIndex, stepIndex, kIndex)
          latRot = interpInfo%allLatRot(subGridIndex, headerIndex, procIndex, stepIndex, kIndex)
          lonRot = interpInfo%allLonRot(subGridIndex, headerIndex, procIndex, stepIndex, kIndex)

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
    real(8)                 :: fieldUU_in(:,:)
    real(8)                 :: fieldVV_in(:,:)
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

        indexBeg = interpInfo%depotIndexBeg(subGridIndex, headerIndex, procIndex, stepIndex, kIndex)
        indexEnd = interpInfo%depotIndexEnd(subGridIndex, headerIndex, procIndex, stepIndex, kIndex)

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
          lat = interpInfo%allLat(headerIndex, procIndex, stepIndex, kIndex)
          lon = interpInfo%allLon(headerIndex, procIndex, stepIndex, kIndex)
          latRot = interpInfo%allLatRot(subGridIndex, headerIndex, procIndex, stepIndex, kIndex)
          lonRot = interpInfo%allLonRot(subGridIndex, headerIndex, procIndex, stepIndex, kIndex)

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
    real(8)                 :: fieldUU_out(:,:)
    real(8)                 :: fieldVV_out(:,:)
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

        indexBeg = interpInfo%depotIndexBeg(subGridIndex, headerIndex, procIndex, stepIndex, kIndex)
        indexEnd = interpInfo%depotIndexEnd(subGridIndex, headerIndex, procIndex, stepIndex, kIndex)

        if ( indexEnd < IndexBeg ) cycle subGrid_loop

        ! now rotate the wind vector and return the desired component
        if ( interpInfo%hco%rotated ) then
          lat = interpInfo%allLat(headerIndex, procIndex, stepIndex, kIndex)
          lon = interpInfo%allLon(headerIndex, procIndex, stepIndex, kIndex)
          latRot = interpInfo%allLatRot(subGridIndex, headerIndex, procIndex, stepIndex, kIndex)
          lonRot = interpInfo%allLonRot(subGridIndex, headerIndex, procIndex, stepIndex, kIndex)

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
      ierr = utl_getPositionXY( stateVector % hco % EZscintID,   &
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

  end subroutine s2c_bgcheck_bilin

  !--------------------------------------------------------------------------
  ! s2c_column_hbilin
  !--------------------------------------------------------------------------
  subroutine s2c_column_hbilin(field,vlev,nlong,nlat,nlev,xlong,xlat, &
                               plong,plat,vprof,vlevout,nlevout)
    !
    ! :Purpose: Horizontal bilinear interpolation from a 3D field to a profile at (plong,plat).
    !           Assumes vertical interpolation not needed or already done.
    !
    !           This version can be used with fields that are not part of the background state,
    !           such as climatologies.
    !
    !           This version does not depend in column_data and gridstatevector modules.
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

  !--------------------------------------------------------------------------
  ! s2c_setupHorizInterp
  !--------------------------------------------------------------------------
  subroutine s2c_setupHorizInterp(footprintRadius_r4, interpInfo, obsSpaceData, &
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
    type(struct_obs)       , intent(inout) :: obsSpaceData
    type(struct_gsv)       , intent(in)    :: stateVector
    integer                , intent(in)    :: headerIndex, kIndex, stepIndex
    integer                , intent(in)    :: procIndex
    integer                , intent(out)   :: numGridpt(interpInfo%hco%numSubGrid)

    if ( footprintRadius_r4 > 0.0 ) then

      call s2c_setupFootprintInterp(footprintRadius_r4, interpInfo, obsSpaceData, stateVector, headerIndex, kIndex, stepIndex, procIndex, numGridpt)

    else if ( footprintRadius_r4 == bilinearFootprint ) then

      call s2c_setupBilinearInterp(interpInfo, obsSpaceData, stateVector, headerIndex, kIndex, stepIndex, procIndex, numGridpt)

    else if ( footprintRadius_r4 == lakeFootprint ) then

      call s2c_setupLakeInterp(interpInfo, obsSpaceData, stateVector, headerIndex, kIndex, stepIndex, procIndex, numGridpt)

    else if ( footprintRadius_r4 == nearestNeighbourFootprint ) then

      call s2c_setupNearestNeighbor(interpInfo, obsSpaceData, stateVector, headerIndex, kIndex, stepIndex, procIndex, numGridpt)

    else

      write(*,*) 'footprint radius = ',footprintRadius_r4
      call utl_abort('s2c_setupHorizInterp: footprint radius not permitted')

    end if

  end subroutine s2c_setupHorizInterp

  !------------------------------------------------------------------
  ! s2c_getFootprintRadius
  !------------------------------------------------------------------
  function s2c_getFootprintRadius( obsSpaceData, headerIndex ) result(fpr)
    !
    !:Purpose: To determine the footprint radius (metres) of the observation.
    !          In the case of bilinear horizontal interpolation,
    !          the returned footprint is zero (default).
    !
    implicit none
    real(4)                       :: fpr

    ! Arguments:
    type(struct_obs), intent(in)  :: obsSpaceData
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

    else

      fpr = bilinearFootprint

    end if

  end function s2c_getFootprintRadius

  !--------------------------------------------------------------------------
  ! s2c_setupBilinearInterp
  !--------------------------------------------------------------------------
  subroutine s2c_setupBilinearInterp(interpInfo, obsSpaceData, stateVector, &
                                     headerIndex, kIndex, stepIndex, procIndex,&
                                     numGridpt)
    !
    !:Purpose: To determine the grid points and their associated weights
    !          for the bilinear horizontal interpolation.
    !
    implicit none

    ! Arguments:
    type(struct_interpInfo), intent(inout) :: interpInfo
    type(struct_obs)       , intent(inout) :: obsSpaceData
    type(struct_gsv)       , intent(in)    :: stateVector
    integer                , intent(in)    :: headerIndex, kIndex, stepIndex
    integer                , intent(in)    :: procIndex
    integer                , intent(out)   :: numGridpt(interpInfo%hco%numSubGrid)

    ! Locals:
    integer :: localHeaderIndex, bodyIndex, depotIndex
    integer :: ierr
    integer :: bodyIndexBeg, bodyIndexEnd, niP1
    integer :: latIndex, lonIndex, latIndex2, lonIndex2, lonIndexP1
    integer :: subGridIndex, subGridForInterp, numSubGridsForInterp
    integer :: ipoint, gridptCount
    integer :: latIndexVec(4), lonIndexVec(4)
    integer :: mask(2,2)
    real(8) :: WeightVec(4)
    real(8) :: dldx, dldy
    real(8) :: weightsSum
    real(4) :: lon_deg_r4, lat_deg_r4
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4

    numGridpt(:) = 0

    lat_deg_r4 = real(interpInfo%allLat(headerIndex, procIndex, stepIndex, kIndex) *  &
                 MPC_DEGREES_PER_RADIAN_R8)
    lon_deg_r4 = real(interpInfo%allLon(headerIndex, procIndex, stepIndex, kIndex) *  &
                 MPC_DEGREES_PER_RADIAN_R8)
    ierr = utl_getPositionXY( stateVector%hco%EZscintID,   &
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

    if ( allocated(stateVector%hco%mask) ) then
      mask(1,1) = stateVector%hco%mask(lonIndex  ,latIndex)
      mask(2,1) = stateVector%hco%mask(lonIndexP1,latIndex)
      mask(1,2) = stateVector%hco%mask(lonIndex  ,latIndex + 1)
      mask(2,2) = stateVector%hco%mask(lonIndexP1,latIndex + 1)
    else
      mask(:,:) = 1
    end if

    do subGridForInterp = 1, numSubGridsForInterp

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

      if ( mask(1,1) == 1 ) then
        gridptCount = gridptCount + 1
        latIndexVec(gridptCount) = latIndex
        lonIndexVec(gridptCount) = lonIndex
        WeightVec(gridptCount) = (1.d0-dldx) * (1.d0-dldy)
      end if

      if ( mask(2,1) == 1 ) then
        gridptCount = gridptCount + 1
        latIndexVec(gridptCount) = latIndex
        lonIndexVec(gridptCount) = lonIndexP1
        WeightVec(gridptCount) =       dldx  * (1.d0-dldy)
      end if

      if ( mask(1,2) == 1 ) then
        gridptCount = gridptCount + 1
        latIndexVec(gridptCount) = latIndex + 1
        lonIndexVec(gridptCount) = lonIndex
        WeightVec(gridptCount) = (1.d0-dldx) *       dldy
      end if

      if ( mask(2,2) == 1 ) then
        gridptCount = gridptCount + 1
        latIndexVec(gridptCount) = latIndex + 1
        lonIndexVec(gridptCount) = lonIndexP1
        WeightVec(gridptCount) =       dldx  *       dldy
      end if

      weightsSum = sum(WeightVec(1:gridptCount))
      if ( weightsSum > 0.d0 ) then

        WeightVec(1:gridptCount) = WeightVec(1:gridptCount) / weightsSum

      else

        if ( procIndex == mpi_myid + 1 ) then

          localHeaderIndex = interpInfo%allHeaderIndex(headerIndex,stepIndex,procIndex)

          write(*,*) 's2c_setupBilinearInterp: Rejecting OBS outside the grid domain, index ', localHeaderIndex
          write(*,*) ' lat-lon (deg) y x : ', lat_deg_r4, lon_deg_r4, ypos_r4, xpos_r4
          bodyIndexBeg = obs_headElem_i(obsSpaceData, OBS_RLN, localHeaderIndex)
          bodyIndexEnd = obs_headElem_i(obsSpaceData, OBS_NLV, localHeaderIndex) + bodyIndexBeg -1
          do bodyIndex = bodyIndexBeg, bodyIndexEnd
            call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
          end do
          call obs_headSet_i(obsSpaceData, OBS_ST1, localHeaderIndex,  &
       ibset( obs_headElem_i(obsSpaceData, OBS_ST1, localHeaderIndex), 05))

        end if

      end if

      ! divide weight by number of subGrids
      WeightVec(1:gridptCount) = WeightVec(1:gridptCount) / real(numSubGridsForInterp,8)

      if ( allocated(interpInfo%interpWeightDepot) ) then

        depotIndex = interpInfo%depotIndexBeg(subGridIndex, headerIndex, procIndex, stepIndex, kIndex)

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
  subroutine s2c_setupFootprintInterp(fpr, interpInfo, obsSpaceData, &
                                      stateVector, headerIndex, kIndex, &
                                      stepIndex, procIndex, numGridpt)
    !
    !:Purpose: To determine the grid points and their associated weights
    !          for the footprint horizontal interpolation.
    !
    implicit none

    ! Arguments:
    real(4)                , intent(in)    :: fpr ! footprint radius (metres)
    type(struct_interpInfo), intent(inout) :: interpInfo
    type(struct_obs)       , intent(inout) :: obsSpaceData
    type(struct_gsv)       , intent(in)    :: stateVector
    integer                , intent(in)    :: headerIndex, kIndex, stepIndex
    integer                , intent(in)    :: procIndex
    integer                , intent(out)   :: numGridpt(interpInfo%hco%numSubGrid)

    ! Locals:
    integer :: localHeaderIndex, bodyIndex, depotIndex
    integer :: ierr
    integer :: bodyIndexBeg, bodyIndexEnd
    integer :: latIndexCentre, lonIndexCentre, latIndexCentre2, lonIndexCentre2
    integer :: subGridIndex, subGridForInterp, numSubGridsForInterp
    real(4) :: lon_deg_r4, lat_deg_r4
    real(8) :: lon_rad, lat_rad
    real(8) :: grid_lon_rad, grid_lat_rad
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    real(4) :: grid_lon_deg_r4, grid_lat_deg_r4
    integer :: ipoint, gridptCount
    integer :: top, bottom, left, right, rectangleCount
    real(8) :: dist
    integer :: lonIndex, latIndex, rectangleIndex, rectangleSize
    integer :: lonIndexVec(statevector%ni*statevector%nj), latIndexVec(statevector%ni*statevector%nj)
    integer :: rectLonIndex(2*(statevector%ni+statevector%nj)-4), rectLatIndex(4*(statevector%ni+statevector%nj)-4)
    logical :: inside, reject

    ! external functions
    integer :: gdllfxy

    reject = .false.

    numGridpt(:) = 0

    ! Determine the grid point nearest the observation.

    lat_rad = interpInfo%allLat(headerIndex, procIndex, stepIndex, kIndex)
    lon_rad = interpInfo%allLon(headerIndex, procIndex, stepIndex, kIndex)
    lat_deg_r4 = real(lat_rad * MPC_DEGREES_PER_RADIAN_R8)
    lon_deg_r4 = real(lon_rad * MPC_DEGREES_PER_RADIAN_R8)
    ierr = utl_getPositionXY( stateVector%hco%EZscintID,   &
                              xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                              lat_deg_r4, lon_deg_r4, subGridIndex )

    lonIndexCentre = nint(xpos_r4)
    latIndexCentre = nint(ypos_r4)
    lonIndexCentre2 = nint(xpos2_r4)
    latIndexCentre2 = nint(ypos2_r4)

    if ( subGridIndex == 3 ) then
      write(*,*) 's2c_setupFootprintInterp: revise code'
      call utl_abort('s2c_setupFootprintInterp: both subGrids involved in interpolation.')
      numSubGridsForInterp = 2
      subGridIndex = 1
    else
      ! only 1 subGrid involved in interpolation
      numSubGridsForInterp = 1
    end if

    do subGridForInterp = 1, numSubGridsForInterp

      gridptCount = 0

      ! If observation is not on the grid, don't use it.
      if ( lonIndexCentre < 1 .or. lonIndexCentre > statevector%ni .or.  &
           latIndexCentre < 1 .or. latIndexCentre > statevector%nj ) reject = .true.

      if ( allocated(stateVector%hco%mask) ) then
        if ( stateVector%hco%mask(lonIndexCentre,latIndexCentre) == 0 ) reject = .true.
      end if

      if ( .not. reject ) then

        !<<<  These four lines assure to use at least the nearrest neighbor
        !     in the case where the footprint size is smaller than the grid spacing.
        gridptCount = 1
        lonIndexVec(gridptCount) = lonIndexCentre
        latIndexVec(gridptCount) = latIndexCentre
        rectangleSize = 1
        !>>>

        inside = .true.
        WHILE_INSIDE: do while(inside)

          inside = .false.

          ! Set up rectangle. We will look in this rectangle
          ! for neighbors.

          top    = latIndexCentre + rectangleSize
          bottom = latIndexCentre - rectangleSize
          left   = lonIndexCentre - rectangleSize
          right  = lonIndexCentre + rectangleSize

          rectangleCount = 0
          latIndex = bottom
          do lonIndex = left, right
            rectangleCount = rectangleCount + 1
            rectLonIndex(rectangleCount) = lonIndex
            rectLatIndex(rectangleCount) = latIndex
          end do
          lonIndex = right
          do latIndex = bottom + 1, top
            rectangleCount = rectangleCount + 1
            rectLonIndex(rectangleCount) = lonIndex
            rectLatIndex(rectangleCount) = latIndex
          end do
          latIndex = top
          do lonIndex = right - 1, left, -1
            rectangleCount = rectangleCount + 1
            rectLonIndex(rectangleCount) = lonIndex
            rectLatIndex(rectangleCount) = latIndex
          end do
          lonIndex = left
          do latIndex = top - 1, bottom + 1, -1
            rectangleCount = rectangleCount + 1
            rectLonIndex(rectangleCount) = lonIndex
            rectLatIndex(rectangleCount) = latIndex
          end do

          do rectangleIndex = 1, rectangleCount

            lonIndex = rectLonIndex(rectangleIndex)
            if (lonIndex >= 1 .and. lonIndex <= statevector%ni) then

              xpos_r4 = real(lonIndex)

              latIndex = rectLatIndex(rectangleIndex)
              if (latIndex >= 1 .and. latIndex <= statevector%nj) then

                ypos_r4 = real(latIndex)

                ierr = gdllfxy(stateVector%hco%EZscintID, grid_lat_deg_r4, grid_lon_deg_r4, &
                       xpos_r4, ypos_r4, 1)

                if(grid_lon_deg_r4 < 0.0) grid_lon_deg_r4 = grid_lon_deg_r4 + 360.0

                grid_lat_rad = real(grid_lat_deg_r4,8)*MPC_RADIANS_PER_DEGREE_R8
                grid_lon_rad = real(grid_lon_deg_r4,8)*MPC_RADIANS_PER_DEGREE_R8

                ! Compute distance between grid point and observation point.
                dist = phf_calcDistance(grid_lat_rad, grid_lon_rad, lat_rad, lon_rad)

                ! If the point is within the footprint, add it to the neighborhood.
                if(dist < fpr) then

                  ! Ignore points that are masked out.
                  if ( allocated(stateVector%hco%mask) ) then
                    if (stateVector%hco%mask(lonIndex, latIndex) == 0) then
                      reject = .true.
                      exit WHILE_INSIDE
                    end if
                  end if

                  if ( .not. reject ) then
                    inside = .true.

                    gridptCount = gridptCount + 1
                    lonIndexVec(gridptCount) = lonIndex
                    latIndexVec(gridptCount) = latIndex

                  end if

                end if

              end if

            end if

          end do ! rectangleIndex

          rectangleSize = rectangleSize + 1

        end do WHILE_INSIDE

      end if ! not reject

      if ( reject ) then

        if ( procIndex == mpi_myid + 1 ) then

          localHeaderIndex = interpInfo%allHeaderIndex(headerIndex,stepIndex,procIndex)

          write(*,*) 's2c_setupFootprintInterp: Rejecting OBS outside the grid domain, index ', localHeaderIndex
          write(*,*) ' lat-lon (deg) : ', lat_deg_r4, lon_deg_r4
          bodyIndexBeg = obs_headElem_i(obsSpaceData, OBS_RLN, localHeaderIndex)
          bodyIndexEnd = obs_headElem_i(obsSpaceData, OBS_NLV, localHeaderIndex) + bodyIndexBeg -1
          do bodyIndex = bodyIndexBeg, bodyIndexEnd
            call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
          end do
          call obs_headSet_i(obsSpaceData, OBS_ST1, localHeaderIndex,  &
                             ibset( obs_headElem_i(obsSpaceData, OBS_ST1, localHeaderIndex), 05))

        end if

      else

        if ( allocated(interpInfo%interpWeightDepot) ) then

          depotIndex = interpInfo%depotIndexBeg(subGridIndex, headerIndex, procIndex, stepIndex, kIndex)

          do ipoint=1,gridptCount

            interpInfo%interpWeightDepot(depotIndex) = 1.0d0 / real(gridptCount,8)
            interpInfo%latIndexDepot(depotIndex)     = latIndexVec(ipoint)
            interpInfo%lonIndexDepot(depotIndex)     = lonIndexVec(ipoint)
            depotIndex = depotIndex + 1

          end do

        end if

        numGridpt(subGridIndex) = gridptCount

      end if

    end do ! subGrid

  end subroutine s2c_setupFootprintInterp

  !--------------------------------------------------------------------------
  ! s2c_setupLakeInterp
  !--------------------------------------------------------------------------
  subroutine s2c_setupLakeInterp(interpInfo, obsSpaceData, stateVector, &
                                 headerIndex, kIndex, stepIndex, procIndex, &
                                 numGridpt)
    !
    !:Purpose: To determine the grid points and their associated weights
    !          for the lake horizontal interpolation.
    !
    implicit none

    ! Arguments:
    type(struct_interpInfo), intent(inout) :: interpInfo
    type(struct_obs)       , intent(inout) :: obsSpaceData
    type(struct_gsv)       , intent(in)    :: stateVector
    integer                , intent(in)    :: headerIndex, kIndex, stepIndex
    integer                , intent(in)    :: procIndex
    integer                , intent(out)   :: numGridpt(interpInfo%hco%numSubGrid)

    ! Locals:
    integer :: localHeaderIndex, bodyIndex, depotIndex
    integer :: ierr
    integer :: bodyIndexBeg, bodyIndexEnd
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

    numGridpt(:) = 0

    reject = .false.

    numGridpt(:) = 0

    ! Determine the grid point nearest the observation.

    lat_rad = interpInfo%allLat(headerIndex, procIndex, stepIndex, kIndex)
    lon_rad = interpInfo%allLon(headerIndex, procIndex, stepIndex, kIndex)
    lat_deg_r4 = real(lat_rad * MPC_DEGREES_PER_RADIAN_R8)
    lon_deg_r4 = real(lon_rad * MPC_DEGREES_PER_RADIAN_R8)
    ierr = utl_getPositionXY( stateVector%hco%EZscintID,   &
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

!     It happens that the lake location is closest to a grid point
!     where MASK(I,J) = 0 while there are other grid points for the
!     same lake where MASK(I,J) = 1. Code needs modifications
!     for this case.

      ! If observation is not on the grid, don't use it.
      if ( lonIndexCentre < 1 .or. lonIndexCentre > statevector%ni .or.  &
           latIndexCentre < 1 .or. latIndexCentre > statevector%nj ) reject = .true.

      if ( allocated(stateVector%hco%mask) ) then
        if ( stateVector%hco%mask(lonIndexCentre,latIndexCentre) == 0 ) reject = .true.
      end if

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
                if(stateVector%hco%mask(lonIndex,latIndex) == 1 .and. .not. lake(lonIndex,latIndex)) then
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

          depotIndex = interpInfo%depotIndexBeg(subGridIndex, headerIndex, procIndex, stepIndex, kIndex)

          do ipoint=1,gridptCount

            interpInfo%interpWeightDepot(depotIndex) = 1.0d0 / real(gridptCount,8)
            interpInfo%latIndexDepot(depotIndex)     = latIndexVec(ipoint)
            interpInfo%lonIndexDepot(depotIndex)     = lonIndexVec(ipoint)
            depotIndex = depotIndex + 1

          end do

        end if

        numGridpt(subGridIndex) = gridptCount

      else

        if ( procIndex == mpi_myid + 1 ) then

          localHeaderIndex = interpInfo%allHeaderIndex(headerIndex,stepIndex,procIndex)

          write(*,*) 's2c_setupLakeInterp: Rejecting OBS outside the grid domain, index ', localHeaderIndex
          write(*,*) ' lat-lon (deg) : ', lat_deg_r4, lon_deg_r4
          write(*,*) ' grid index : ',lonIndexCentre, latIndexCentre
          write(*,*) ' mask : ',stateVector%hco%mask(lonIndexCentre,latIndexCentre)
          bodyIndexBeg = obs_headElem_i(obsSpaceData, OBS_RLN, localHeaderIndex)
          bodyIndexEnd = obs_headElem_i(obsSpaceData, OBS_NLV, localHeaderIndex) + bodyIndexBeg -1
          do bodyIndex = bodyIndexBeg, bodyIndexEnd
            call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
          end do
          call obs_headSet_i(obsSpaceData, OBS_ST1, localHeaderIndex,  &
                             ibset( obs_headElem_i(obsSpaceData, OBS_ST1, localHeaderIndex), 05))

        end if

      end if ! not reject

    end do ! subGrid

  end subroutine s2c_setupLakeInterp

  !--------------------------------------------------------------------------
  ! s2c_setupNearestNeighbor
  !--------------------------------------------------------------------------
  subroutine s2c_setupNearestNeighbor(interpInfo, obsSpaceData, stateVector, &
                                      headerIndex, kIndex, stepIndex, procIndex, numGridpt)
    !
    !:Purpose: Determine the nearest grid points to the observtions location
    !
    implicit none

    ! arguments
    type(struct_interpInfo), intent(inout) :: interpInfo
    type(struct_obs)       , intent(inout) :: obsSpaceData
    type(struct_gsv)       , intent(in)    :: stateVector
    integer                , intent(in)    :: headerIndex, kIndex, stepIndex, procIndex
    integer                , intent(out)   :: numGridpt(interpInfo%hco%numSubGrid)

    ! locals
    integer :: localHeaderIndex, bodyIndex, depotIndex
    integer :: ierr
    integer :: bodyIndexBeg, bodyIndexEnd
    integer :: latIndex, lonIndex
    integer :: subGridIndex
    real(4) :: lon_deg_r4, lat_deg_r4
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4

    if ( stateVector%hco%grtyp == 'U' ) then
      call utl_abort('s2c_setupNearestNeighbor: Yin-Yang grid not supported')
    end if

    numGridpt(:) = 0

    lat_deg_r4 = real(interpInfo%allLat(headerIndex, procIndex, stepIndex, kIndex) *  &
                 MPC_DEGREES_PER_RADIAN_R8)
    lon_deg_r4 = real(interpInfo%allLon(headerIndex, procIndex, stepIndex, kIndex) *  &
                 MPC_DEGREES_PER_RADIAN_R8)

    ierr = utl_getPositionXY( stateVector%hco%EZscintID,   &
                              xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                              lat_deg_r4, lon_deg_r4, subGridIndex )

    latIndex = nint(ypos_r4)
    lonIndex = nint(xpos_r4)

    ! Handle periodicity in longitude
    if ( lonIndex == statevector%ni+1 .and. stateVector%hco%grtyp == 'G' ) lonIndex = 1

    ! Test bounds
    if ( lonIndex < 1 .or. lonIndex > statevector%ni .or. &
         latIndex < 1 .or. latIndex > statevector%nj  ) then

      if ( procIndex == mpi_myid + 1 ) then
        
        localHeaderIndex = interpInfo%allHeaderIndex(headerIndex,stepIndex,procIndex)

        write(*,*) 's2c_setupNearestNeighbor: Rejecting OBS outside the grid domain, index ', localHeaderIndex
        write(*,*) ' lat-lon (deg) y x : ', lat_deg_r4, lon_deg_r4, ypos_r4, xpos_r4
        bodyIndexBeg = obs_headElem_i(obsSpaceData, OBS_RLN, localHeaderIndex)
        bodyIndexEnd = obs_headElem_i(obsSpaceData, OBS_NLV, localHeaderIndex) + bodyIndexBeg -1
        do bodyIndex = bodyIndexBeg, bodyIndexEnd
          call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
        end do
        call obs_headSet_i(obsSpaceData, OBS_ST1, localHeaderIndex,  &
             ibset( obs_headElem_i(obsSpaceData, OBS_ST1, localHeaderIndex), 05))

      end if

    end if

    if ( allocated(interpInfo%interpWeightDepot) ) then
      
      depotIndex = interpInfo%depotIndexBeg(subGridIndex, headerIndex, procIndex, stepIndex, kIndex)

      interpInfo%interpWeightDepot(depotIndex) = 1.d0
      interpInfo%latIndexDepot    (depotIndex) = latIndex
      interpInfo%lonIndexDepot    (depotIndex) = lonIndex

    end if

    numGridpt(subGridIndex) = 1

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
    if ( column%nk /= gsv_getNumK(statevector) ) call utl_abort('checkColumnStatevectorMatch: column%nk /= gsv_getNumK(statevector)')

    ! loop through k and check varNames are same between column/statevector
    do kIndex = 1, column%nk
      if ( gsv_getVarNameFromK(statevector,kIndex) /= col_getVarNameFromK(column,kIndex) ) call utl_abort('checkColumnStatevectorMatch: varname in column and statevector do not match')
    end do

  end subroutine checkColumnStatevectorMatch

  !--------------------------------------------------------------------------
  ! latlonChecks
  !--------------------------------------------------------------------------
  subroutine latlonChecks( obsSpaceData, hco, headerIndex, rejectOutsideObs, latLev_T, lonLev_T, latLev_M, lonLev_M )
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

    ! Locals:
    integer :: ierr
    integer :: bodyIndex, bodyIndexBeg, bodyIndexEnd, niP1, subGridIndex
    integer :: nlev_T, lev_T, nlev_M, lev_M
    real(4) :: lon_r4, lat_r4, lon_deg_r4, lat_deg_r4
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    real(4) :: xposRecalc_r4, yposRecalc_r4, xpos2Recalc_r4, ypos2Recalc_r4
    logical :: latlonOutsideGrid, latlonRecalcOutsideGrid, rejectHeader

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
    ierr = utl_getPositionXY( hco%EZscintID,   &
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
      ierr = utl_getPositionXY( hco%EZscintID,   &
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
      xpos_r4 = real(hco%ni)/2.0
      ypos_r4 = real(hco%nj)/2.0
      ierr = gdllfxy(hco%EZscintID, lat_deg_r4, lon_deg_r4, &
                     xpos_r4, ypos_r4, 1)

      lonLev_T(:) = real(lon_deg_r4 * MPC_RADIANS_PER_DEGREE_R4,8)
      latLev_T(:) = real(lat_deg_r4 * MPC_RADIANS_PER_DEGREE_R4,8)
      lonLev_M(:) = real(lon_deg_r4 * MPC_RADIANS_PER_DEGREE_R4,8)
      latLev_M(:) = real(lat_deg_r4 * MPC_RADIANS_PER_DEGREE_R4,8)

    end if

  end subroutine latlonChecks

end module stateToColumn_mod

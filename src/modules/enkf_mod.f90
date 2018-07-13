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

!--------------------------------------------------------------------------
!! MODULE enkf (prefix="enkf")
!!
!! *Purpose*: Implementation of the EnKF in MIDAS.
!!
!--------------------------------------------------------------------------
MODULE enkf_mod
  use mpi_mod
  use ramDisk_mod
  use gridStateVector_mod
  use mathPhysConstants_mod
  use utilities_mod
  use fileNames_mod
  use varNameList_mod
  use obsTimeInterp_mod
  use tt2phi_mod
  use obsSpaceData_mod
  use columnData_mod
  implicit none
  save
  private

  ! public procedures
  public :: enkf_readMember, enkf_setupColumnsFromEnsemble
  public :: enkf_computeColumnsMean, enkf_computeColumnsPerturbations
  public :: enkf_extractObsRealBodyColumn, enkf_extractObsIntBodyColumn
  public :: enkf_gatherHX

  integer, external :: get_max_rss

contains

  !--------------------------------------------------------------------------
  !! subroutine enkf_readMember
  !! *Purpose*: locally read one member on each MPI process and put in stateVector
  !--------------------------------------------------------------------------
  subroutine enkf_readMember(stateVector, ensPathName, memberIndex)
    implicit none

    ! arguments
    type(struct_gsv) :: stateVector
    character(len=*) :: ensPathName
    integer          :: memberIndex

    ! locals
    real(4), pointer   :: ptr4d_r4(:,:,:,:)
    real(4)            :: multFactor
    integer            :: stepIndex, varIndex, ierr
    character(len=256) :: ensFileName
    character(len=2)   :: typvar
    character(len=12)  :: etiket
    character(len=4)   :: varName

    write(*,*) 'enkf_readMember: starting'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if ( .not. stateVector%allocated ) then
      call utl_abort('enkf_readMember: stateVector object not allocated!')
    end if

    ! Get filename for the requested member
    call fln_ensFileName(ensFileName, ensPathName, memberIndex)

    ! Read the file and remove from ramDisk
    write(*,*) 'enkf_readMember: starting to read member number ', memberIndex
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    typvar = ' '
    etiket = ' '
    call gsv_readFile(stateVector, ensFileName, etiket, typvar, readGZsfc_opt=.true.)
    ierr = ram_remove(ensFileName)

    ! Do unit conversion
    VAR_LOOP: do varIndex = 1, vnl_numvarmax
      varName = vnl_varNameList(varIndex)
      if ( .not. gsv_varExist(stateVector, varName) ) cycle VAR_LOOP
      ptr4d_r4 => gsv_getField_r4(stateVector, varName)

      if ( trim(varName) == 'UU' .or. trim(varName) == 'VV') then
        multFactor = MPC_M_PER_S_PER_KNOT_R4 ! knots -> m/s
      else if ( trim(varName) == 'P0' ) then
        multFactor = MPC_PA_PER_MBAR_R4 ! hPa -> Pa
      else
        multFactor = 1.0 ! no conversion
      end if

      !$OMP PARALLEL DO PRIVATE (stepIndex)
      STEP_LOOP: do stepIndex = 1, stateVector%numStep
        if ( multFactor /= 1.0 ) ptr4d_r4(:,:,:,stepIndex) = multFactor * ptr4d_r4(:,:,:,stepIndex)
        if (trim(varName) == 'TT' ) then
          ptr4d_r4(:,:,:,stepIndex) = ptr4d_r4(:,:,:,stepIndex) + MPC_K_C_DEGREE_OFFSET_R4
        end if
        if (trim(varName) == 'HU' ) then
          ptr4d_r4(:,:,:,stepIndex) = sngl(max(real(ptr4d_r4(:,:,:,stepIndex),8),MPC_MINIMUM_HU_R8))
        end if
      end do STEP_LOOP
      !$OMP END PARALLEL DO

    end do VAR_LOOP

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'enkf_readMember: finished reading the ensemble member'

  end subroutine enkf_readMember


  subroutine enkf_allGatherObsGrid(obsSpaceData, stateVector, allNumObs, allObsGid, allHeaderIndexVec)
    implicit none

    ! arguments
    type(struct_obs) :: obsSpaceData
    type(struct_gsv) :: stateVector
    integer, pointer :: allNumObs(:,:), allObsGid(:,:), allHeaderIndexVec(:,:,:)

    ! locals
    integer :: numHeader, numHeaderMax, headerIndex, bodyIndex, numStep, stepIndex, ierr
    integer :: bodyIndexBeg, bodyIndexEnd, procIndex
    integer :: ig1obs, ig2obs, ig3obs, ig4obs
    real(8) :: zig1, zig2, zig3, zig4, stepObsIndex
    real(4) :: lon_r4, lat_r4, lon_deg_r4, lat_deg_r4, xpos_r4, ypos_r4
    integer, allocatable :: numObs(:), headerIndexVec(:,:)
    real(4), allocatable :: lonVec(:), latVec(:), allLonVec(:,:), allLatVec(:,:)
    integer :: ezgdef, gdxyfll, gdllfxy

    write(*,*) 'enkf_allGatherObsGrid: STARTING'

    numStep = stateVector%numStep
    numHeader = obs_numheader(obsSpaceData)
    call rpn_comm_allreduce(numHeader, numHeaderMax, 1,  &
         'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)

    allObsGid(:,:) = 0
    allNumObs(:,:) = 0
    allHeaderIndexVec(:,:,:) = 0

    ! temporary arrays
    allocate(lonVec(numHeaderMax))
    allocate(latVec(numHeaderMax))
    allocate(allLonVec(numHeaderMax,mpi_nprocs))
    allocate(allLatVec(numHeaderMax,mpi_nprocs))
    allocate(numObs(numStep))
    allocate(headerIndexVec(numHeaderMax,numStep))
    headerIndexVec(:,:) = 0

    numObs(:) = 0

    STEP_LOOP: do stepIndex = 1, numStep

      lonVec(:) = 0.0
      latVec(:) = 0.0

      HEADER_LOOP: do headerIndex = 1, numHeader

        if ( oti_timeInterpWeightAllZero(headerIndex) .and. stepIndex == 1 ) then
          ! obs is outside of assimilation window and stepIndex is 1 (has to go somewhere)

          write(*,*) 'enkf_allGatherObsGrid: Observation time outside assimilation window: ',  &
               obs_headElem_i(obsSpaceData,OBS_DAT,headerIndex),obs_headElem_i(obsSpaceData,OBS_ETM,headerIndex)

          ! flag these observations as out of time domain and turn off its assimilation flag
          bodyIndexBeg = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
          bodyIndexEnd = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) + bodyIndexBeg -1
          do bodyIndex = bodyIndexBeg, bodyIndexEnd
            call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, 0)
          end do
          call obs_headSet_i(obsSpaceData, OBS_ST1, headerIndex,  &
               ibset( obs_headElem_i(obsSpaceData, OBS_ST1, headerIndex), 05))
        else
          ! if obs inside window, but zero weight for current stepIndex then skip it
          if ( oti_getTimeInterpWeight(headerIndex,stepIndex) == 0.0d0 ) cycle HEADER_LOOP

        end if

        numObs(stepIndex) = numObs(stepIndex) + 1
        headerIndexVec(numObs(stepIndex),stepIndex) = headerIndex

        !- Get LatLon of observation location
        lat_r4 = real(obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex), 4)
        lon_r4 = real(obs_headElem_r(obsSpaceData, OBS_LON, headerIndex), 4)
        if (lon_r4 <  0.0          ) lon_r4 = lon_r4 + 2.0*MPC_PI_R4
        if (lon_r4 >= 2.0*MPC_PI_R4) lon_r4 = lon_r4 - 2.0*MPC_PI_R4

        lat_deg_r4 = lat_r4 * MPC_DEGREES_PER_RADIAN_R4 ! Radian To Degree
        lon_deg_r4 = lon_r4 * MPC_DEGREES_PER_RADIAN_R4

        !
        !- Find the position in the trial field grid
        !
        ierr = gdxyfll(stateVector%hco%EZscintID, xpos_r4, ypos_r4, lat_deg_r4, lon_deg_r4, 1)

        if ( xpos_r4 >= 1.0 .and. xpos_r4 <= real(stateVector%ni) .and.  &
             ypos_r4 >= 1.0 .and. ypos_r4 <= real(stateVector%nj) ) then

          lonVec(numObs(stepIndex)) = lon_r4 * MPC_DEGREES_PER_RADIAN_R4
          latVec(numObs(stepIndex)) = lat_r4 * MPC_DEGREES_PER_RADIAN_R4

        else
          ! The observation is outside the domain
          ! With a LAM trial field we must discard this observation
          write(*,*) 'enkf_allGatherObsGrid: Rejecting OBS outside the TRIAL field domain, ', headerIndex
          write(*,*) '  position : ', lat_deg_r4, lon_deg_r4, ypos_r4, xpos_r4

          bodyIndexBeg = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
          bodyIndexEnd = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) + bodyIndexBeg -1
          do bodyIndex = bodyIndexBeg, bodyIndexEnd
            call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, 0)
          end do
          call obs_headSet_i(obsSpaceData, OBS_ST1, headerIndex,  &
               ibset( obs_headElem_i(obsSpaceData, OBS_ST1, headerIndex), 05))

          ! However, we must assigned a realistic lat-lon to this point
          ! to avoid problem later in Hx computation.
          ierr = gdllfxy(stateVector%hco%EZscintID, lat_deg_r4, lon_deg_r4, real(stateVector%ni)/2.0,  &
               real(stateVector%nj)/2.0, 1) ! Middle of the domain
          lonVec(numObs(stepIndex)) = lon_deg_r4
          latVec(numObs(stepIndex)) = lat_deg_r4
        end if

      end do HEADER_LOOP

      ! gather the number of obs over all processors for each timestep
      call rpn_comm_allgather(numObs(stepIndex),      1, 'MPI_INTEGER', &
                              allNumObs(stepIndex,:), 1, 'MPI_INTEGER', &
                              'GRID',ierr)

      ! gather lon-lat of observations from all processors
      call rpn_comm_allgather(lonVec,    numHeaderMax, 'MPI_REAL4', &
                              allLonVec, numHeaderMax, 'MPI_REAL4', &
                              'GRID', ierr)
      call rpn_comm_allgather(latVec,    numHeaderMax, 'MPI_REAL4', &
                              allLatVec, numHeaderMax, 'MPI_REAL4', &
                              'GRID', ierr)

      zig1 = 0.0D0
      zig2 = 0.0D0
      zig3 = 1.0D0
      zig4 = 1.0D0
      call utl_cxgaig('L',ig1obs,ig2obs,ig3obs,ig4obs,zig1,zig2,zig3,zig4)

      do procIndex = 1, mpi_nprocs
        if (allNumObs(stepIndex,procIndex) > 0) then
          allObsGid(stepIndex,procIndex) = ezgdef(allNumObs(stepIndex,procIndex),  &
               1,'Y','L',ig1obs,ig2obs,ig3obs,ig4obs,  &
               allLonVec(1:allNumObs(stepIndex,procIndex),procIndex),  &
               allLatVec(1:allNumObs(stepIndex,procIndex),procIndex))
        else
          allObsGid(stepIndex,procIndex) = -999
        end if
      end do

    end do STEP_LOOP

    ! gather the headerIndexVec arrays onto all processors
    call rpn_comm_allgather(headerIndexVec,    numHeaderMax*numStep, 'MPI_INTEGER', &
                            allHeaderIndexVec, numHeaderMax*numStep, 'MPI_INTEGER', &
                            'GRID',ierr)

    deallocate(lonVec)
    deallocate(latVec)
    deallocate(allLonVec)
    deallocate(allLatVec)
    deallocate(numObs)
    deallocate(headerIndexVec)

    write(*,*) 'enkf_allGatherObsGrid: FINISHED'

  end subroutine enkf_allGatherObsGrid


  subroutine enkf_setupColumnsFromEnsemble(stateVector,obsSpaceData,columns)
    implicit none

    ! arguments
    type(struct_gsv)        :: stateVector
    type(struct_obs)        :: obsSpaceData
    type(struct_columnData) :: columns(:)

    ! locals
    logical :: verbose = .true.
    integer :: varIndex, levIndex, numLev, stepIndex, numStep
    integer :: headerIndex, numHeader, numHeaderMax, yourNumHeader
    integer :: memberIndex, nEns, procIndex, nsize, ierr, iset, headerIndex2
    integer :: ezdefset, ezqkdef
    real(8) :: gzSfc_col
    real(8) :: weight
    character(len=4)     :: varName
    integer, pointer     :: allNumObs(:,:), allObsGid(:,:), allHeaderIndexVec(:,:,:)
    real(8), pointer     :: column_ptr(:), ptr2d_r8(:,:), allCols_ptr(:,:), allCols_VV_ptr(:,:)
    real(4), pointer     :: ptr4d_r4(:,:,:,:), ptrUU4d_r4(:,:,:,:), ptrVV4d_r4(:,:,:,:)
    real(4), allocatable :: cols_hint_r4(:,:,:),cols_hint_VV_r4(:,:,:)
    real(4), allocatable :: cols_send_r4(:,:),cols_send_VV_r4(:,:)
    real(4), allocatable :: cols_recv_r4(:,:),cols_recv_VV_r4(:,:)
    real(8), allocatable :: cols_send_1proc_r8(:), cols_send_VV_1proc_r8(:)
    real(4), allocatable :: gzSfc_r4(:,:)
    logical              :: beSilent

    numStep = stateVector%numStep
    numHeader = obs_numheader(obsSpaceData)
    call rpn_comm_allreduce(numHeader, numHeaderMax, 1,  &
                            'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)
    nEns = size(columns)
    write(*,*) 'enkf_setupColumnsFromEnsemble: nEns =', nEns

    ! compute and collect all obs grids onto all mpi tasks
    allocate(allObsGid(numStep,mpi_nprocs))
    allocate(allNumObs(numStep,mpi_nprocs))
    allocate(allHeaderIndexVec(numHeaderMax,numStep,mpi_nprocs))
    call enkf_allGatherObsGrid(obsSpaceData, stateVector, allNumObs, allObsGid, allHeaderIndexVec)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if ( mpi_myid == 0 ) then
      do stepIndex = 1, numStep
        write(*,*) 'enkf_setupColumnsFromEnsemble: stepIndex, allNumObs = ', stepIndex, allNumObs(stepIndex,:)
      end do
    end if

    ! arrays for interpolated columns for 1 level/variable and each time step
    allocate(cols_hint_r4(maxval(allNumObs),numStep,mpi_nprocs))
    allocate(cols_hint_VV_r4(maxval(allNumObs),numStep,mpi_nprocs))
    cols_hint_r4(:,:,:) = 0.0
    cols_hint_VV_r4(:,:,:) = 0.0

    ! arrays for sending/receiving time interpollated columns for 1 level/variable
    allocate(cols_send_r4(numHeaderMax,mpi_nprocs))
    allocate(cols_send_VV_r4(numHeaderMax,mpi_nprocs))
    cols_send_r4(:,:) = 0.0
    cols_send_VV_r4(:,:) = 0.0

    allocate(cols_recv_r4(numHeaderMax,mpi_nprocs))
    allocate(cols_recv_VV_r4(numHeaderMax,mpi_nprocs))
    cols_recv_r4(:,:) = 0.0
    cols_recv_VV_r4(:,:) = 0.0

    allocate(cols_send_1proc_r8(numHeaderMax))
    allocate(cols_send_VV_1proc_r8(numHeaderMax))

    allocate(gzSfc_r4(stateVector%myLonBeg:stateVector%myLonEnd,  &
                      stateVector%myLatBeg:stateVector%myLatEnd))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    write(*,*) 'before setting to zero' ; call flush(6)
    ! set contents of all columns to zero
    do memberIndex = 1, nEns
      allCols_ptr => col_getAllColumns(columns(memberIndex))
      allCols_ptr(:,:) = 0.0d0
      if ( memberIndex == 1 ) write(*,*) 'shape(allCols_ptr) = ', shape(allCols_ptr)
    end do
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    VAR_LOOP: do varIndex = 1, vnl_numvarmax
      varName = vnl_varNameList(varIndex)
      if ( .not. gsv_varExist(varName=varName) ) cycle
      if ( trim(varName) == 'VV' ) cycle  ! handled together with UU

      write(*,*) 'enkf_setupColumnsFromEnsemble: doing interpolation for varName = ', varName
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      numLev = gsv_getNumLevFromVarName(stateVector,varName)
      LEV_LOOP: do levIndex = 1, numLev

        call tmg_start(140,'ENKF_HINTERP')
        STEP_LOOP: do stepIndex = 1, numStep

          if ( maxval(allNumObs(stepIndex,:)) == 0 ) cycle STEP_LOOP

          ! interpolate to the columns destined for all procs for all steps and one level/variable
          PROC_LOOP: do procIndex = 1, mpi_nprocs
            yourNumHeader = allNumObs(stepIndex,procIndex)
            if ( yourNumHeader > 0 ) then
              iset = ezdefset(allObsGid(stepIndex,procIndex),stateVector%hco%EZscintID)
              if (trim(varName) == 'UU') then
                ptrUU4d_r4 => gsv_getField_r4(stateVector, 'UU')
                ptrVV4d_r4 => gsv_getField_r4(stateVector, 'VV')
                ierr = utl_ezuvint( cols_hint_r4(1:yourNumHeader,stepIndex,procIndex),  &
                                    cols_hint_VV_r4(1:yourNumHeader,stepIndex,procIndex),  &
                                    ptrUU4d_r4(:,:,levIndex,stepIndex), ptrVV4d_r4(:,:,levIndex,stepIndex),  &
                                    interpDegree='LINEAR' )
              else
                ptr4d_r4 => gsv_getField_r4(stateVector, varName)
                ierr = utl_ezsint( cols_hint_r4(1:yourNumHeader,stepIndex,procIndex),  &
                                   ptr4d_r4(:,:,levIndex,stepIndex), &
                                   interpDegree='LINEAR' )
              end if
            end if
          end do PROC_LOOP

        end do STEP_LOOP
        call tmg_stop(140)

        call tmg_start(141,'ENKF_TINTERP')
        ! interpolate in time to the columns destined for all procs and one level/variable
        PROC_LOOP2: do procIndex = 1, mpi_nprocs
          cols_send_1proc_r8(:) = 0.0d0
          if (trim(varName) == 'UU') then
            cols_send_VV_1proc_r8(:) = 0.0d0
          end if
          STEP_LOOP2: do stepIndex = 1, numStep
            !$OMP PARALLEL DO PRIVATE (headerIndex, headerIndex2, weight)
            do headerIndex = 1, allNumObs(stepIndex, procIndex)
              headerIndex2 = allHeaderIndexVec(headerIndex,stepIndex,procIndex)
              weight = oti_getTimeInterpWeightMpiGlobal(headerIndex2,stepIndex,procIndex)
              cols_send_1proc_r8(headerIndex2) = cols_send_1proc_r8(headerIndex2) &
                            + weight * real(cols_hint_r4(headerIndex,stepIndex,procIndex),8)

              if (trim(varName) == 'UU') then
                cols_send_VV_1proc_r8(headerIndex2) = cols_send_VV_1proc_r8(headerIndex2) &
                                 + weight * real(cols_hint_VV_r4(headerIndex,stepIndex,procIndex),8)
              end if

            end do
            !$OMP END PARALLEL DO
          end do STEP_LOOP2
          cols_send_r4(:,procIndex) = real(cols_send_1proc_r8(:),4)
          if (trim(varName) == 'UU') then
            cols_send_VV_r4(:,procIndex) = real(cols_send_VV_1proc_r8(:),4)
          end if
        end do PROC_LOOP2
        call tmg_stop(141)

        call tmg_start(150,'ENKF_BARR')
        call rpn_comm_barrier('GRID',ierr)
        call tmg_stop(150)

        call tmg_start(142,'ENKF_ALLTOALL')
        ! mpi communication: alltoall for one level/variable
        nsize = numHeaderMax
        if(mpi_nprocs > 1) then
          call rpn_comm_alltoall(cols_send_r4, nsize, 'MPI_REAL4',  &
                                 cols_recv_r4, nsize, 'MPI_REAL4', 'GRID', ierr)
        else
          cols_recv_r4(:,1) = cols_send_r4(:,1)
        end if

        if (trim(varName) == 'UU') then
          if(mpi_nprocs > 1) then
            call rpn_comm_alltoall(cols_send_VV_r4, nsize, 'MPI_REAL4',  &
                                   cols_recv_VV_r4, nsize, 'MPI_REAL4', 'GRID', ierr)
          else
            cols_recv_VV_r4(:,1) = cols_send_VV_r4(:,1)
          end if
        end if
        call tmg_stop(142)

        call tmg_start(143,'ENKF_RESHUFFLE')
        ! reorganize ensemble of distributed columns
        do memberIndex = 1, nEns
          allCols_ptr => col_getAllColumns(columns(memberIndex),varName)
          if (trim(varName) == 'UU') then
            allCols_VV_ptr => col_getAllColumns(columns(memberIndex),'VV')
          end if

          !$OMP PARALLEL DO PRIVATE (headerIndex)
          do headerIndex = 1, numHeader
            allCols_ptr(levIndex,headerIndex) = real(cols_recv_r4(headerIndex,memberIndex),8)
            if (trim(varName) == 'UU') then
              allCols_VV_ptr(levIndex,headerIndex) = real(cols_recv_VV_r4(headerIndex,memberIndex),8)
            end if
          end do
          !$OMP END PARALLEL DO

        end do
        call tmg_stop(143)

      end do LEV_LOOP

    end do VAR_LOOP

    ! Interpolate surface GZ separately
    varName = 'GZ'
    write(*,*) 'enkf_setupColumnsFromEnsemble: doing interpolation for varName = ', varName
    STEP_LOOP_GZ: do stepIndex = 1, numStep

      if ( maxval(allNumObs(stepIndex,:)) == 0 ) cycle STEP_LOOP_GZ

      ! interpolate to the columns destined for all procs for all steps and one level/variable
      PROC_LOOP_GZ: do procIndex = 1, mpi_nprocs
        yourNumHeader = allNumObs(stepIndex,procIndex)
        if ( yourNumHeader > 0 ) then
          iset = ezdefset(allObsGid(stepIndex,procIndex),stateVector%hco%EZscintID)
          ptr2d_r8 => gsv_getGZsfc(stateVector)
          gzSfc_r4(:,:) = real(ptr2d_r8(:,:),4)
          ierr = utl_ezsint( cols_hint_r4(1:yourNumHeader,stepIndex,procIndex),  &
                             gzSfc_r4(:,:), interpDegree='LINEAR' )
        end if
      end do PROC_LOOP_GZ

    end do STEP_LOOP_GZ

    ! interpolate in time to the columns destined for all procs and one level/variable
    PROC_LOOP_GZ2: do procIndex = 1, mpi_nprocs
      cols_send_r4(:,procIndex) = 0.0
      STEP_LOOP_GZ2: do stepIndex = 1, numStep
        !$OMP PARALLEL DO PRIVATE (headerIndex, headerIndex2)
        do headerIndex = 1, allNumObs(stepIndex, procIndex)
          headerIndex2 = allHeaderIndexVec(headerIndex,stepIndex,procIndex)
          ! just copy, since surface GZ same for all time steps
          cols_send_r4(headerIndex2,procIndex) = cols_hint_r4(headerIndex,stepIndex,procIndex)
        end do
        !$OMP END PARALLEL DO
      end do STEP_LOOP_GZ2
    end do PROC_LOOP_GZ2

    ! mpi communication: alltoall for one level/variable
    nsize = numHeaderMax
    if(mpi_nprocs > 1) then
      call rpn_comm_alltoall(cols_send_r4, nsize, 'MPI_REAL4',  &
                             cols_recv_r4, nsize, 'MPI_REAL4', 'GRID', ierr)
    else
      cols_recv_r4(:,1) = cols_send_r4(:,1)
    end if

    ! reorganize ensemble of distributed columns
    do memberIndex = 1, nEns
      !$OMP PARALLEL DO PRIVATE (headerIndex, gzSfc_col)
      do headerIndex = 1, numHeader
        gzSfc_col = real(cols_recv_r4(headerIndex,memberIndex),8)
        call col_setGZsfc(columns(memberIndex), headerIndex, gzSfc_col)
      end do
      !$OMP END PARALLEL DO
    end do

    deallocate(allObsGid)
    deallocate(allNumObs)
    deallocate(allHeaderIndexVec)

    ! Free up memory now that ensemble member no longer needed
    call gsv_deallocate(stateVector)

    ! Do final preparations of columnData objects (compute GZ and pressure)
    call tmg_start(144,'ENKF_PRESSURE')
    do memberIndex = 1, nEns
      if (col_varExist('P0')) then
        beSilent = .true.
        if ( memberIndex == 1 ) beSilent = .false.
        call col_calcPressure(columns(memberIndex),beSilent_opt=beSilent)
      end if
    end do
    call tmg_stop(144)

    call tmg_start(145,'ENKF_HEIGHT')
    do memberIndex = 1, nEns
      if (col_varExist('TT') .and. col_varExist('HU') .and. col_varExist('P0') .and. col_getNumLev(columns(1),'MM') > 1) then
        beSilent = .true.
        if ( memberIndex == 1 ) beSilent = .false.
        call tt2phi(columns(memberIndex),beSilent_opt=beSilent)
      end if
    end do
    call tmg_stop(145)

    if ( verbose ) then
      write(*,*) '======================='
      write(*,*) 'Contents of columns(1):'
      column_ptr => col_getColumn(columns(1),1,'UU')
      write(*,*) 'column UU = ', column_ptr(:)
      column_ptr => col_getColumn(columns(1),1,'VV')
      write(*,*) 'column VV = ', column_ptr(:)
      column_ptr => col_getColumn(columns(1),1,'TT')
      write(*,*) 'column TT = ', column_ptr(:)
      column_ptr => col_getColumn(columns(1),1,'HU')
      write(*,*) 'column LQ = ', column_ptr(:)
      column_ptr => col_getColumn(columns(1),1,'P0')
      write(*,*) 'column P0 = ', column_ptr(:)
      column_ptr => col_getColumn(columns(1),1,'TG')
      write(*,*) 'column TG = ', column_ptr(:)
      write(*,*) '======================='

      do levIndex = 1, col_getNumLev(columns(1),'MM')
        write(*,*) 'enkf_setupColumnsFromEnsemble: levIndex, col_getPressure(columns(1),levIndex,1,MM) = ',  &
                   levIndex,col_getPressure(columns(1),levIndex,1,'MM')
      end do
      do levIndex = 1, col_getNumLev(columns(1),'MM')
        write(*,*) 'enkf_setupColumnsFromEnsemble: levIndex, col_getHeight(columns(1),levIndex,1,MM) = ',  &
                   levIndex,col_getHeight(columns(1),levIndex,1,'MM')
      end do
    end if

  end subroutine enkf_setupColumnsFromEnsemble


  subroutine enkf_computeColumnsMean(column_mean, columns)
    implicit none

    ! arguments
    type(struct_columnData) :: column_mean, columns(:)

    ! locals
    logical :: verbose = .true.
    integer :: memberIndex, nEns, levIndex
    real(8) :: multFactor
    real(8), pointer :: column_ptr(:)

    call tmg_start(146,'ENKF_COLSMEAN')

    nEns = size(columns)
    multFactor = 1.0d0 / real(nEns,8)
    write(*,*) 'enkf_computeColumnsMean: nEns =', nEns

    call col_copyLatLon(columns(1),column_mean)

    call col_zero(column_mean)

    do memberIndex = 1, nEns

        column_mean%all(:,:) = column_mean%all(:,:) +  &
                               multFactor * columns(memberIndex)%all(:,:)

        column_mean%gz_T(:,:) = column_mean%gz_T(:,:) +  &
                                multFactor * columns(memberIndex)%gz_T(:,:)
        column_mean%gz_M(:,:) = column_mean%gz_M(:,:) +  &
                                multFactor * columns(memberIndex)%gz_M(:,:)
        column_mean%gz_sfc(:) = column_mean%gz_sfc(:) +  &
                                multFactor * columns(memberIndex)%gz_sfc(:)

        column_mean%pressure_T(:,:) = column_mean%pressure_T(:,:) +  &
                                      multFactor * columns(memberIndex)%pressure_T(:,:)
        column_mean%pressure_M(:,:) = column_mean%pressure_M(:,:) +  &
                                      multFactor * columns(memberIndex)%pressure_M(:,:)

        column_mean%dP_dPsfc_T(:,:) = column_mean%dP_dPsfc_T(:,:) +  &
                                      multFactor * columns(memberIndex)%dP_dPsfc_T(:,:)
        column_mean%dP_dPsfc_M(:,:) = column_mean%dP_dPsfc_M(:,:) +  &
                                      multFactor * columns(memberIndex)%dP_dPsfc_M(:,:)

    end do

    !if (col_varExist('P0')) then
    !  call col_calcPressure(column_mean)
    !end if

    if ( verbose ) then
      write(*,*) '======================='
      write(*,*) 'Contents of column_mean:'
      column_ptr => col_getColumn(column_mean,1,'UU')
      write(*,*) 'column UU = ', column_ptr(:)
      column_ptr => col_getColumn(column_mean,1,'VV')
      write(*,*) 'column VV = ', column_ptr(:)
      column_ptr => col_getColumn(column_mean,1,'TT')
      write(*,*) 'column TT = ', column_ptr(:)
      column_ptr => col_getColumn(column_mean,1,'HU')
      write(*,*) 'column LQ = ', column_ptr(:)
      column_ptr => col_getColumn(column_mean,1,'P0')
      write(*,*) 'column P0 = ', column_ptr(:)
      column_ptr => col_getColumn(column_mean,1,'TG')
      write(*,*) 'column TG = ', column_ptr(:)
      write(*,*) '======================='

      do levIndex = 1, col_getNumLev(column_mean,'MM')
        write(*,*) 'enkf_setupColumnsFromEnsemble: levIndex, col_getPressure(column_mean,levIndex,1,MM) = ',  &
                   levIndex,col_getPressure(column_mean,levIndex,1,'MM')
      end do
      do levIndex = 1, col_getNumLev(column_mean,'MM')
        write(*,*) 'enkf_setupColumnsFromEnsemble: levIndex, col_getHeight(column_mean,levIndex,1,MM) = ',  &
                   levIndex,col_getHeight(column_mean,levIndex,1,'MM')
      end do
    end if

    call tmg_stop(146)

  end subroutine enkf_computeColumnsMean


  subroutine enkf_computeColumnsPerturbations(columns, column_mean)
    implicit none

    ! arguments
    type(struct_columnData) :: columns(:), column_mean

    ! locals
    logical :: verbose = .true.
    integer :: memberIndex, nEns, levIndex
    real(8), pointer :: column_ptr(:)

    call tmg_start(147,'ENKF_COLSPERTS')

    nEns = size(columns)
    write(*,*) 'enkf_computeColumnsPerturbations: nEns =', nEns

    !
    ! Remove ensemble mean from all variables, except: gz_sfc, dP_dPsfc_T/M, oltv
    !
    do memberIndex = 1, nEns

        columns(memberIndex)%all(:,:) = columns(memberIndex)%all(:,:) -  &
                                        column_mean%all(:,:)

        columns(memberIndex)%gz_T(:,:) = columns(memberIndex)%gz_T(:,:) -  &
                                         column_mean%gz_T(:,:)
        columns(memberIndex)%gz_M(:,:) = columns(memberIndex)%gz_M(:,:) -  &
                                         column_mean%gz_M(:,:)

        columns(memberIndex)%pressure_T(:,:) = columns(memberIndex)%pressure_T(:,:) -  &
                                      column_mean%pressure_T(:,:)
        columns(memberIndex)%pressure_M(:,:) = columns(memberIndex)%pressure_M(:,:) -  &
                                      column_mean%pressure_M(:,:)

    end do

    if ( verbose ) then
      write(*,*) '======================='
      write(*,*) 'Contents of columns(1):'
      column_ptr => col_getColumn(columns(1),1,'UU')
      write(*,*) 'column UU = ', column_ptr(:)
      column_ptr => col_getColumn(columns(1),1,'VV')
      write(*,*) 'column VV = ', column_ptr(:)
      column_ptr => col_getColumn(columns(1),1,'TT')
      write(*,*) 'column TT = ', column_ptr(:)
      column_ptr => col_getColumn(columns(1),1,'HU')
      write(*,*) 'column LQ = ', column_ptr(:)
      column_ptr => col_getColumn(columns(1),1,'P0')
      write(*,*) 'column P0 = ', column_ptr(:)
      column_ptr => col_getColumn(columns(1),1,'TG')
      write(*,*) 'column TG = ', column_ptr(:)
      write(*,*) '======================='

      do levIndex = 1, col_getNumLev(columns(1),'MM')
        write(*,*) 'enkf_setupColumnsFromEnsemble: levIndex, col_getPressure(columns(1),levIndex,1,MM) = ',  &
                   levIndex,col_getPressure(columns(1),levIndex,1,'MM')
      end do
      do levIndex = 1, col_getNumLev(columns(1),'MM')
        write(*,*) 'enkf_setupColumnsFromEnsemble: levIndex, col_getHeight(columns(1),levIndex,1,MM) = ',  &
                   levIndex,col_getHeight(columns(1),levIndex,1,'MM')
      end do
    end if

    call tmg_stop(147)

  end subroutine enkf_computeColumnsPerturbations


  subroutine enkf_extractObsRealBodyColumn(outputVector, obsSpaceData, obsColumnIndex)
    implicit none

    ! arguments
    real(8)          :: outputVector(:)
    type(struct_obs) :: obsSpaceData
    integer          :: obsColumnIndex

    ! locals
    integer :: bodyIndex

    call tmg_start(148,'ENKF_EXTRACTBODY')

    do bodyIndex = 1, obs_numBody(obsSpaceData)
      outputVector(bodyIndex) = obs_bodyElem_r(obsSpaceData,obsColumnIndex,bodyIndex)
    end do

    call tmg_stop(148)

  end subroutine enkf_extractObsRealBodyColumn


  subroutine enkf_extractObsIntBodyColumn(outputVector, obsSpaceData, obsColumnIndex)
    implicit none

    ! arguments
    integer          :: outputVector(:)
    type(struct_obs) :: obsSpaceData
    integer          :: obsColumnIndex

    ! locals
    integer :: bodyIndex

    call tmg_start(148,'ENKF_EXTRACTBODY')

    do bodyIndex = 1, obs_numBody(obsSpaceData)
      outputVector(bodyIndex) = obs_bodyElem_i(obsSpaceData,obsColumnIndex,bodyIndex)
    end do

    call tmg_stop(148)

  end subroutine enkf_extractObsIntBodyColumn


  subroutine enkf_gatherHX(HXens,HXensT_mpiglobal)
    implicit none

    ! arguments
    real(8) :: HXens(:,:)
    real(8),pointer :: HXensT_mpiglobal(:,:)

    ! locals
    integer :: ierr, nEns, numBody, procIndex, memberIndex, numBody_mpiglobal
    integer :: allNumBody(mpi_nprocs), displs(mpi_nprocs)

    numBody = size(HXens,1)
    nEns     = size(HXens,2)

    write(*,*) 'enkf_gatherHX: starting'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    call rpn_comm_gather( numBody, 1, 'mpi_integer', allNumBody, 1, 'mpi_integer', &
                          0, 'GRID', ierr )
    if ( mpi_myid == 0 ) then
      displs(1) = 0
      do procIndex = 2, mpi_nprocs
        displs(procIndex) = displs(procIndex-1) + allNumBody(procIndex-1)
      end do
    else
      displs(:) = 0
    end if

    numBody_mpiglobal = sum(allNumBody(:))
    if( mpi_myid == 0 ) then
      allocate(HXensT_mpiglobal(nEns,numBody_mpiglobal))
    else
      allocate(HXensT_mpiglobal(nEns,1))
    end if

    do memberIndex = 1, nEns
      call rpn_comm_gatherv( HXens(:,memberIndex), numBody, 'mpi_double_precision', &
                             HXensT_mpiglobal(memberIndex,:), allNumBody, displs, &
                             'mpi_double_precision', 0, 'GRID', ierr )
    end do

    write(*,*) 'enkf_gatherHX: finished'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine enkf_gatherHX


end module enkf_mod

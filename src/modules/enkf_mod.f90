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
  use timeCoord_mod
  use tt2phi_mod
  use obsSpaceData_mod
  use columnData_mod
  implicit none
  save
  private

  ! public procedures
  public :: enkf_readMember, enkf_setupBackgroundColumnsFromEnsemble

  integer, external :: get_max_rss

contains

  !--------------------------------------------------------------------------
  !! subroutine enkf_readMember
  !! *Purpose*: locally read one member on each MPI process and put in stateVector
  !--------------------------------------------------------------------------
  subroutine enkf_readMember(stateVector, ensPathName, memberIndex, ctrlVarHumidity)
    implicit none

    ! arguments
    type(struct_gsv) :: stateVector
    character(len=*) :: ensPathName
    integer          :: memberIndex
    character(len=*) :: ctrlVarHumidity

    ! locals
    real(4), pointer   :: ptr4d_r4(:,:,:,:)
    real(4)            :: multFactor
    integer            :: stepIndex, varIndex, ierr
    character(len=256) :: ensFileName
    character(len=2)   :: typvar
    character(len=12)  :: etiket
    character(len=4)   :: varName
    logical            :: HUcontainsLQinFile

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
    call gsv_readFile(stateVector, ensFileName, etiket, typvar, HUcontainsLQinFile, readGZsfc_opt=.true.)
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

        if (trim(varName) == 'HU' .and. ctrlVarHumidity == 'LQ' .and. .not.HUcontainsLQinFile ) then
          ptr4d_r4(:,:,:,stepIndex) = sngl(log(max(real(ptr4d_r4(:,:,:,stepIndex),8),MPC_MINIMUM_HU_R8)))
        else if (trim(varName) == 'HU' .and. ctrlVarHumidity == 'HU' .and. .not.HUcontainsLQinFile ) then
          ptr4d_r4(:,:,:,stepIndex) = sngl(max(real(ptr4d_r4(:,:,:,stepIndex),8),MPC_MINIMUM_HU_R8))
        end if
      end do STEP_LOOP
      !$OMP END PARALLEL DO

    end do VAR_LOOP

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'enkf_readMember: finished reading the ensemble member'

  end subroutine enkf_readMember


  subroutine enkf_allGatherObsGrid(obsSpaceData, stateVector, allNumObs, allObsGid, headerIndexVec)
    implicit none

    ! arguments
    type(struct_obs) :: obsSpaceData
    type(struct_gsv) :: stateVector
    integer, pointer :: allNumObs(:,:), allObsGid(:,:), headerIndexVec(:,:)

    ! locals
    integer :: numHeader, numHeaderMax, headerIndex, bodyIndex, numStep, stepIndex, ierr
    integer :: bodyIndexBeg, bodyIndexEnd, procIndex
    integer :: ig1obs, ig2obs, ig3obs, ig4obs
    real(8) :: zig1, zig2, zig3, zig4, stepObsIndex
    real(4) :: lon_r4, lat_r4, lon_deg_r4, lat_deg_r4, xpos_r4, ypos_r4
    integer, allocatable :: numObs(:)
    real(4), allocatable :: lonVec(:), latVec(:), allLonVec(:,:), allLatVec(:,:)
    integer :: ezgdef, ezsetopt, gdxyfll, gdllfxy

    write(*,*) 'enkf_allGatherObsGrid: STARTING'

    ierr = ezsetopt('INTERP_DEGREE', 'LINEAR')

    numStep = stateVector%numStep
    numHeader = obs_numheader(obsSpaceData)
    call rpn_comm_allreduce(numHeader, numHeaderMax, 1,  &
         'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)

    allObsGid(:,:) = 0
    allNumObs(:,:) = 0
    headerIndexVec(:,:) = 0

    ! temporary arrays
    allocate(lonVec(numHeaderMax))
    allocate(latVec(numHeaderMax))
    allocate(allLonVec(numHeaderMax,mpi_nprocs))
    allocate(allLatVec(numHeaderMax,mpi_nprocs))
    allocate(numObs(numStep))

    numObs(:) = 0

    STEP_LOOP: do stepIndex = 1, numStep

      lonVec(:) = 0.0
      latVec(:) = 0.0

      HEADER_LOOP: do headerIndex=1, numHeader

        if ( tim_timeInterpWeightAllZero(headerIndex) .and. stepIndex == 1 ) then
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
          if ( tim_getTimeInterpWeight(headerIndex,stepIndex) == 0.0d0 ) cycle HEADER_LOOP
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

      ! gather and compute the number of obs over all processors for each timestep
      call rpn_comm_allgather(numObs(stepIndex), 1, 'MPI_INTEGER',       &
           allNumObs(stepIndex,:), 1, 'MPI_INTEGER', &
           'GRID',ierr)

      ! gather lon-lat of observations from all processors
      call rpn_comm_allgather(lonVec, numHeaderMax, 'MPI_REAL4',       &
           allLonVec, numHeaderMax, 'MPI_REAL4', &
           'GRID', ierr)
      call rpn_comm_allgather(latVec, numHeaderMax, 'MPI_REAL4',       &
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

    deallocate(lonVec)
    deallocate(latVec)
    deallocate(allLonVec)
    deallocate(allLatVec)
    deallocate(numObs)

    write(*,*) 'enkf_allGatherObsGrid: FINISHED'

  end subroutine enkf_allGatherObsGrid


  subroutine enkf_setupBackgroundColumnsFromEnsemble(stateVector,obsSpaceData,columns)
    implicit none

    ! arguments
    type(struct_gsv)        :: stateVector
    type(struct_obs)        :: obsSpaceData
    type(struct_columnData) :: columns(:)

    ! locals
    logical :: verbose = .true.
    integer :: varIndex, levIndex, numLev, stepIndex, numStep, headerIndex, numObs
    integer :: memberIndex, nEns, procIndex, nsize, ierr, iset, headerIndex2
    integer :: ezdefset, ezuvint, ezsint, ezqkdef
    real(8) :: gzSfc_col
    character(len=4)     :: varName
    integer, pointer     :: allNumObs(:,:), allObsGid(:,:), headerIndexVec(:,:)
    real(8), pointer     :: column_ptr(:), ptr2d_r8(:,:)
    real(4), pointer     :: ptr4d_r4(:,:,:,:), ptrUU4d_r4(:,:,:,:), ptrVV4d_r4(:,:,:,:)
    real(4), allocatable :: cols_send_r4(:,:,:),cols_send_VV_r4(:,:,:)
    real(4), allocatable :: cols_recv_r4(:,:,:),cols_recv_VV_r4(:,:,:)
    real(4), allocatable :: gzSfc_r4(:,:)

    numStep = stateVector%numStep
    nEns = size(columns)
    write(*,*) 'enkf_setupBackgroundColumnsFromEnsemble: nEns =', nEns

    ! compute and collect all obs grids onto all mpi tasks
    allocate(allObsGid(numStep,mpi_nprocs))
    allocate(allNumObs(numStep,mpi_nprocs))
    allocate(headerIndexVec(obs_numheader(obsSpaceData),numStep))
    call enkf_allGatherObsGrid(obsSpaceData, stateVector, allNumObs, allObsGid, headerIndexVec)

    ! arrays for sending/receiving columns for 1 level/variable
    allocate(cols_send_r4(maxval(allNumObs),numStep,mpi_nprocs))
    allocate(cols_send_VV_r4(maxval(allNumObs),numStep,mpi_nprocs))
    cols_send_r4(:,:,:) = 0.0
    cols_send_VV_r4(:,:,:) = 0.0

    allocate(cols_recv_r4(maxval(allNumObs),numStep,mpi_nprocs))
    allocate(cols_recv_VV_r4(maxval(allNumObs),numStep,mpi_nprocs))
    cols_recv_r4(:,:,:) = 0.0
    cols_recv_VV_r4(:,:,:) = 0.0

    allocate(gzSfc_r4(stateVector%myLonBeg:stateVector%myLonEnd,  &
                      stateVector%myLatBeg:stateVector%myLatEnd))

    VAR_LOOP: do varIndex = 1, vnl_numvarmax
      varName = vnl_varNameList(varIndex)
      if ( .not. gsv_varExist(varName=varName) ) cycle
      if ( trim(varName) == 'VV' ) cycle  ! handled together with UU

      numLev = gsv_getNumLevFromVarName(stateVector,varName)
      LEV_LOOP: do levIndex = 1, numLev

        write(*,*) 'doing interpolation for varName, level=', varName, levIndex; call flush(6)
        STEP_LOOP: do stepIndex = 1, numStep

          if ( maxval(allNumObs(stepIndex,:)) == 0 ) cycle STEP_LOOP

          ! interpolate to the columns destined for all procs for all steps and one level/variable
          PROC_LOOP: do procIndex = 1, mpi_nprocs
            numObs = allNumObs(stepIndex,procIndex)
            if ( numObs > 0 ) then
              iset = ezdefset(allObsGid(stepIndex,procIndex),stateVector%hco%EZscintID)
              if (trim(varName).eq.'UU') then
                ptrUU4d_r4 => gsv_getField_r4(stateVector, 'UU')
                ptrVV4d_r4 => gsv_getField_r4(stateVector, 'VV')
                ierr = ezuvint( cols_send_r4(1:numObs,stepIndex,procIndex),  &
                                cols_send_VV_r4(1:numObs,stepIndex,procIndex),  &
                                ptrUU4d_r4(:,:,levIndex,stepIndex), ptrVV4d_r4(:,:,levIndex,stepIndex) )
              else
                ptr4d_r4 => gsv_getField_r4(stateVector, varName)
                ierr = ezsint( cols_send_r4(1:numObs,stepIndex,procIndex),  &
                               ptr4d_r4(:,:,levIndex,stepIndex) )
              end if
            end if
          end do PROC_LOOP

        end do STEP_LOOP

        ! mpi communication: alltoall for one level/variable
        write(*,*) 'doing alltoall for varName, level=', varName, levIndex; call flush(6)
        nsize = maxval(allNumObs) * numStep
        if(mpi_nprocs > 1) then
          call rpn_comm_alltoall(cols_send_r4, nsize, 'MPI_REAL4',  &
                                 cols_recv_r4, nsize, 'MPI_REAL4', 'GRID', ierr)
        else
          cols_recv_r4(:,:,1) = cols_send_r4(:,:,1)
        end if

        if (trim(varName).eq.'UU') then
          if(mpi_nprocs > 1) then
            call rpn_comm_alltoall(cols_send_VV_r4, nsize, 'MPI_REAL4',  &
                                   cols_recv_VV_r4, nsize, 'MPI_REAL4', 'GRID', ierr)
          else
            cols_recv_VV_r4(:,:,1) = cols_send_VV_r4(:,:,1)
          end if
        end if

        ! reorganize ensemble of distributed columns
        write(*,*) 'reorganizing received data into columns'; call flush(6)
        do memberIndex = 1, nEns
          do stepIndex = 1, numStep
            do headerIndex = 1, allNumObs(stepIndex, mpi_myid+1)
              headerIndex2 = headerIndexVec(headerIndex,stepIndex)

              column_ptr => col_getColumn(columns(memberIndex),headerIndex2,varName)
              column_ptr(levIndex) = column_ptr(levIndex) + tim_getTimeInterpWeight(headerIndex2,stepIndex)  &
                                                          * real(cols_recv_r4(headerIndex,stepIndex,memberIndex),8)

              if (trim(varName).eq.'UU') then
                column_ptr => col_getColumn(columns(memberIndex),headerIndex2,'VV')
                column_ptr(levIndex) = column_ptr(levIndex) + tim_getTimeInterpWeight(headerIndex2,stepIndex)  &
                                                            * real(cols_recv_VV_r4(headerIndex,stepIndex,memberIndex),8)
              end if

            end do
          end do
        end do

      end do LEV_LOOP

    end do VAR_LOOP

    ! Interpolate surface GZ separately
    varName = 'GZ'
    write(*,*) 'doing interpolation for varName=', varName; call flush(6)
    STEP_LOOP_GZ: do stepIndex = 1, numStep

      if ( maxval(allNumObs(stepIndex,:)) == 0 ) cycle STEP_LOOP_GZ

      ! interpolate to the columns destined for all procs for all steps and one level/variable
      PROC_LOOP_GZ: do procIndex = 1, mpi_nprocs
        numObs = allNumObs(stepIndex,procIndex)
        if ( numObs > 0 ) then
          iset = ezdefset(allObsGid(stepIndex,procIndex),stateVector%hco%EZscintID)
          ptr2d_r8 => gsv_getGZsfc(stateVector)
          gzSfc_r4(:,:) = ptr2d_r8(:,:)
          ierr = ezsint( cols_send_r4(1:numObs,stepIndex,procIndex),  &
                         gzSfc_r4(:,:) )
        end if
      end do PROC_LOOP_GZ

    end do STEP_LOOP_GZ

    ! mpi communication: alltoall for one level/variable
    write(*,*) 'doing alltoall for varName, level=', varName; call flush(6)
    nsize = maxval(allNumObs) * numStep
    if(mpi_nprocs > 1) then
      call rpn_comm_alltoall(cols_send_r4, nsize, 'MPI_REAL4',  &
                             cols_recv_r4, nsize, 'MPI_REAL4', 'GRID', ierr)
    else
      cols_recv_r4(:,:,1) = cols_send_r4(:,:,1)
    end if

    ! reorganize ensemble of distributed columns
    write(*,*) 'reorganizing received data into columns'; call flush(6)
    do memberIndex = 1, nEns
      do stepIndex = 1, numStep
        do headerIndex = 1, allNumObs(stepIndex, mpi_myid+1)
          headerIndex2 = headerIndexVec(headerIndex,stepIndex)

          gzSfc_col = real(cols_recv_r4(headerIndex,stepIndex,memberIndex),8)
          call col_setGZsfc(columns(memberIndex), headerIndex2, gzSfc_col)
        end do
      end do
    end do

    deallocate(allObsGid)
    deallocate(allNumObs)
    deallocate(headerIndexVec)

    ! Do final preparations of columnData objects (compute GZ and pressure)
    do memberIndex = 1, nEns
      if (col_varExist('P0')) then
        call col_calcPressure(columns(memberIndex))
      end if
      if (col_varExist('TT') .and. col_varExist('HU') .and. col_varExist('P0') .and. col_getNumLev(columns(1),'MM') > 1) then
        call tt2phi(columns(memberIndex))
      end if
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
        write(*,*) 'enkf_setupBackgroundColumnsFromEnsemble: levIndex, col_getPressure(columns(1),levIndex,1,MM) = ',  &
                   levIndex,col_getPressure(columns(1),levIndex,1,'MM')
      end do
      do levIndex = 1, col_getNumLev(columns(1),'MM')
        write(*,*) 'enkf_setupBackgroundColumnsFromEnsemble: levIndex, col_getHeight(columns(1),levIndex,1,MM) = ',  &
                   levIndex,col_getHeight(columns(1),levIndex,1,'MM')
      end do
    end if

  end subroutine enkf_setupBackgroundColumnsFromEnsemble


end module enkf_mod

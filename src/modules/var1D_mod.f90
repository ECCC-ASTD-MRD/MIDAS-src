
module var1D_mod
  ! MODULE var1D_mod (prefix='var1D' category='4. Data Object transformations')
  !
  ! :Purpose: contains all 1Dvar-related methods.
  !
  use columnData_mod
  use gridStatevector_mod
  use horizontalCoord_mod
  use midasMpi_mod 
  use obsSpaceData_mod
  use timeCoord_mod
  use verticalCoord_mod
  use codeprecision_mod
  use mathphysconstants_mod

  implicit none
  save
  private

  ! public procedures
  public :: var1D_Setup, var1D_Finalize
  public :: var1D_transferColumnToYGrid
  ! public variables
  public :: var1D_validHeaderIndex, var1D_validHeaderCount

  logical              :: initialized = .false.
  integer, external    :: get_max_rss
  integer, allocatable :: var1D_validHeaderIndex(:)    ! pointeur vers les colonnes assimilables pour minimiser la taille du vecteur de controle
  integer              :: var1D_validHeaderCount !taille effective de  var1D_validHeaderIndex

contains

  !--------------------------------------------------------------------------
  !  var1D_setup
  !--------------------------------------------------------------------------
  subroutine var1D_setup(obsSpaceData)
    !
    ! :Purpose: to setup var1D module
    !
    implicit none

    ! Arguments:
    type (struct_obs), intent(in) :: obsSpaceData

    ! Locals:
    integer :: countGood, headerIndex
    integer :: bodyStart, bodyEnd, bodyIndex

    if(mmpi_myid == 0) write(*,*) 'var1D_setup: Starting'
    if(mmpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

    !we want to count how many obs are really assimilable to minimize controlvector size
    var1D_validHeaderCount = 0
    allocate(  var1D_validHeaderIndex(obs_numHeader(obsSpaceData)) )
    do headerIndex = 1, obs_numHeader(obsSpaceData)
      bodyStart = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
      bodyEnd = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) + bodyStart - 1
      countGood = 0
      do bodyIndex = bodyStart, bodyEnd
        if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) == obs_assimilated) countGood = countGood + 1
      end do
      if (countGood > 0)  then
        var1D_validHeaderCount = var1D_validHeaderCount + 1
        var1D_validHeaderIndex(var1D_validHeaderCount) = headerIndex
        if (var1D_validHeaderCount == 1) write(*,*) 'first OBS', headerIndex
      end if
    end do
    initialized = .true.
    write(*,*) 'var1D_setup: var1D_validHeaderCount, obs_numHeader(obsSpaceData)', var1D_validHeaderCount, obs_numHeader(obsSpaceData)

  end subroutine var1D_setup

  !--------------------------------------------------------------------------
  ! var1D_Finalize
  !--------------------------------------------------------------------------
  subroutine var1D_Finalize()
    !
    ! :Purpose: to deallocate memory used by internal module structures
    !
    implicit none

    if (initialized) then
       deallocate( var1D_validHeaderIndex )
    end if

  end subroutine var1D_Finalize

  !--------------------------------------------------------------------------
  ! var1D_transferColumnToYGrid
  !--------------------------------------------------------------------------
  subroutine var1D_transferColumnToYGrid( stateVector, obsSpaceData, column, varList)
    !
    ! :Purpose: to transfer content of a columndata object to a statevector object
    !           without interpolation (to be used in 1DVar mode to write increments on Y grid).
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout)        :: stateVector
    type(struct_obs), intent(in)           :: obsSpaceData
    type(struct_columnData), intent(inout) :: column
    character(len=4), intent(in)           :: varList(:)

    ! Locals:
    type(struct_hco), pointer :: hco_yGrid
    integer :: varIndex, globalObsIndex, obsIndex, taskIndex, headerIndex
    integer, allocatable :: var1D_validHeaderCountAllTasks(:), obsOffset(:)
    real(8), pointer :: myColumn(:), myField(:,:,:)
    real(8), allocatable, target :: dummy(:)
    integer :: var1D_validHeaderCountMpiGlobal, var1D_validHeaderCountMax, ierr, status
    real(8) :: lat, lon
    integer :: varDim, tag

    call rpn_comm_barrier("GRID",ierr)
    allocate( obsOffset(0:mmpi_nprocs-1) )
    if (mmpi_myid ==0) then
      allocate( var1D_validHeaderCountAllTasks(mmpi_nprocs) )
    else
      allocate(var1D_validHeaderCountAllTasks(1))
    end if

    call rpn_comm_gather(var1D_validHeaderCount  , 1, 'MPI_INTEGER', var1D_validHeaderCountAllTasks, 1,'MPI_INTEGER', 0, "GRID", ierr )
    if (mmpi_myId ==0) then
      var1D_validHeaderCountMpiGlobal = sum( var1D_validHeaderCountAllTasks(:) )
      var1D_validHeaderCountMax = maxval( var1D_validHeaderCountAllTasks(:) )
      obsOffset(0) = 0
      do taskIndex = 1, mmpi_nprocs - 1
        obsOffset(taskIndex) = obsOffset(taskIndex - 1) + var1D_validHeaderCountAllTasks(taskIndex)
      end do
      write(*,*) 'var1D_transferColumnToYGrid: obsOffset: ', obsOffset(:)
    end if
    call rpn_comm_bcast( obsOffset, mmpi_nprocs, 'MPI_INTEGER', 0,  "GRID",ierr )
    call rpn_comm_bcast( var1D_validHeaderCountMax, 1, 'MPI_INTEGER', 0,  "GRID",ierr )

    call hco_setupYgrid(hco_Ygrid, 1, var1D_validHeaderCountMpiGlobal)
    if (mmpi_myId ==0) then
      call gsv_allocate(stateVector, numstep=tim_nstepobsinc, hco_ptr=hco_Ygrid, vco_ptr=column%vco, &
           datestamp_opt=tim_getDatestamp(), mpi_local_opt=.false., &
           dataKind_opt=pre_incrReal, allocHeight_opt=.false., allocPressure_opt=.false., &
           besilent_opt=.false.)
      write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
    end if

    write(*,*) 'var1D_transferColumnToYGrid: start of lat-lon dissemination'
    do obsIndex = 1, var1D_validHeaderCountMax
      if (obsIndex <= var1D_validHeaderCount ) then
        headerIndex = var1D_validHeaderIndex(obsIndex)      
        lat = obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex)
        lon = obs_headElem_r(obsSpaceData, OBS_LON, headerIndex)
      else
        lat = MPC_missingValue_R8
        lon = MPC_missingValue_R8
      end if
      if (mmpi_myId == 0) then
        if ( obsIndex <= var1D_validHeaderCount ) then
          hco_yGrid%lat2d_4(1, obsIndex) = lat
          hco_yGrid%lon2d_4(1, obsIndex) = lon
        end if
      else
        tag = 2 * mmpi_myID
        call rpn_comm_send( lat, 1, 'mpi_real8', 0, tag,     'GRID', ierr )
        call rpn_comm_send( lon, 1, 'mpi_real8', 0, tag + 1, 'GRID', ierr )
      end if

      if (mmpi_myId == 0) then
        do taskIndex = 1,  mmpi_nprocs - 1
          tag = 2 * taskIndex
          call rpn_comm_recv( lat, 1, 'mpi_real8', taskIndex, tag, 'GRID', status, ierr )
          call rpn_comm_recv( lon, 1, 'mpi_real8', taskIndex, tag+1, 'GRID', status, ierr )
          if (lat /= MPC_missingValue_R8 .and. lon /= MPC_missingValue_R8) then 
            globalObsIndex = obsIndex + obsOffset(taskIndex)
            hco_yGrid%lat2d_4(1, globalObsIndex) = lat
            hco_yGrid%lon2d_4(1, globalObsIndex) = lon
          end if
        end do
      end if
      call rpn_comm_barrier("GRID",ierr)
    end do

    call rpn_comm_barrier("GRID",ierr)
    write(*,*) 'var1D_transferColumnToYGrid: end of lat-lon dissemination'
    
    do varIndex = 1, size(varList)
      write(*,*) 'var1D_transferColumnToYGrid: start of dissemination for ', varList(varIndex)
      if (mmpi_myId == 0 ) then
        call gsv_getField(stateVector, myField, varName_opt=varList(varIndex), stepIndex_opt=1)
        varDim = gsv_getNumLevFromVarName(stateVector, varList(varIndex))
      end if
      call rpn_comm_bcast(varDim, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
      allocate(dummy(varDim))
      dummy(:) = MPC_missingValue_R8
      do obsIndex = 1, var1D_validHeaderCountMax
        if (obsIndex <= var1D_validHeaderCount ) then
          headerIndex = var1D_validHeaderIndex(obsIndex) 
          myColumn => col_getColumn(column, headerIndex, varName_opt=varList(varIndex))
        else
          myColumn => dummy
        end if
        if (mmpi_myId == 0) then
          if ( obsIndex <= var1D_validHeaderCount ) then
            myField(1, obsIndex, :) = myColumn(:)
          end if
        else
          tag = mmpi_myId
          call rpn_comm_send(myColumn , varDim, 'mpi_real8', 0, tag, 'GRID', ierr )
        end if

        if (mmpi_myId == 0) then
          do taskIndex = 1,  mmpi_nprocs - 1
            tag = taskIndex
            call rpn_comm_recv(myColumn,  varDim, 'mpi_real8', taskIndex, tag, 'GRID', status, ierr )
            if (all( myColumn /=  MPC_missingValue_R8)) then
              globalObsIndex = obsIndex + obsOffset(taskIndex)
              myField(1, globalObsIndex, :) = myColumn(:)
            end if
          end do
        end if
      end do

      write(*,*) 'var1D_transferColumnToYGrid: end of dissemination for ', varList(varIndex)
      deallocate(dummy)

    end do

    call rpn_comm_barrier("GRID", ierr)
    deallocate( obsOffset )
    deallocate( var1D_validHeaderCountAllTasks )

  end subroutine var1D_transferColumnToYGrid

end module var1D_mod


module midasMpi_mod
  !
  ! MODULE midasMpi_mod (prefix='mmpi' category='8. Low-level utilities and constants')
  !
  !:Purpose:  Subroutines/functions and public variables related to general aspects of mpi.
  !           Also, subroutine and public variables related to the mpi decomposition
  !           specific to the MIDAS code.
  !
  use mpi_f08 ! this is the Fortran 2008 MPI library module
  !use rpn_comm, only: rpn_comm_mype, rpn_comm_init, rpn_comm_finalize
  use utilities_mod

  implicit none
  save
  private

  ! Public variables
  logical, public, protected :: mmpi_doBarrier = .true.
  integer, public, protected :: mmpi_myid      = 0
  integer, public, protected :: mmpi_myidHost  = 0
  integer, public, protected :: mmpi_nprocs    = 0
  integer, public, protected :: mmpi_myidx     = 0
  integer, public, protected :: mmpi_myidy     = 0
  integer, public, protected :: mmpi_npex      = 0
  integer, public, protected :: mmpi_npey      = 0
  integer, public, protected :: mmpi_numthread = 0
  type(mpi_comm), public, protected :: mmpi_comm_EW, mmpi_comm_NS, mmpi_comm_GRID, mmpi_mpicomm_SHARED
  integer, public, protected :: mmpi_maxTagValue
  integer, public, protected, allocatable :: mmpi_nodeMasters(:)

  ! Public procedures
  public :: mmpi_initialize,mmpi_getptopo
  public :: mmpi_allreduce_sumreal8scalar,mmpi_allgather_string
  public :: mmpi_allreduce_sumR8_1d, mmpi_allreduce_sumR8_2d
  public :: mmpi_reduce_sumR8_1d, mmpi_reduce_sumR8_2d, mmpi_reduce_sumR8_3d
  public :: mmpi_setup_latbands, mmpi_setup_lonbands
  public :: mmpi_setup_m, mmpi_setup_n
  public :: mmpi_setup_levels
  public :: mmpi_setup_varslevels
  public :: mmpi_myidXfromLon, mmpi_myidYfromLat
  public :: mmpi_bcast, mmpi_gather, mmpi_allGather, mmpi_alltoall, mmpi_alltoallv
  public :: mmpi_allReduce, mmpi_gatherv, mmpi_reduce, mmpi_scatterv
  public :: mmpi_send, mmpi_recv, mmpi_sendrecv, mmpi_finalize, mmpi_barrier
  public :: mmpi_stopAndWait4Debug, mmpi_gathervDisplacements

  ! module interfaces
  ! -----------------

  ! general interface for mpi_bcast
  interface mmpi_bcast
    module procedure mmpi_bcast_character
    module procedure mmpi_bcast_logical
    module procedure mmpi_bcast_integer
    module procedure mmpi_bcast_real4
    module procedure mmpi_bcast_real8
  end interface mmpi_bcast

  ! general interface for mpi_gather
  interface mmpi_gather
    module procedure mmpi_gather_logical
    module procedure mmpi_gather_integer
    module procedure mmpi_gather_integer8
    module procedure mmpi_gather_real4
    module procedure mmpi_gather_real8
  end interface mmpi_gather

  ! general interface for mpi_allGather
  interface mmpi_allGather
    module procedure mmpi_allGather_logical
    module procedure mmpi_allGather_integer
    module procedure mmpi_allGather_real4
    module procedure mmpi_allGather_real8
  end interface mmpi_allGather

  ! general interface for mpi_alltoall
  interface mmpi_alltoall
    module procedure mmpi_alltoall_integer
    module procedure mmpi_alltoall_integer8
    module procedure mmpi_alltoall_real4
    module procedure mmpi_alltoall_real8
  end interface mmpi_alltoall

  ! general interface for mpi_alltoallv
  interface mmpi_alltoallv
    module procedure mmpi_alltoallv_real4
    module procedure mmpi_alltoallv_real8
  end interface mmpi_alltoallv

  ! general interface for mpi_allReduce
  interface mmpi_allReduce
    module procedure mmpi_allReduce_logical
    module procedure mmpi_allReduce_integer
    module procedure mmpi_allReduce_integer8
    module procedure mmpi_allReduce_real4
    module procedure mmpi_allReduce_real8
    module procedure mmpi_allReduce_scalar_integer
  end interface mmpi_allReduce

  ! general interface for mpi_gatherv
  interface mmpi_gatherv
    module procedure mmpi_gatherv_logical
    module procedure mmpi_gatherv_logical_displs
    module procedure mmpi_gatherv_integer
    module procedure mmpi_gatherv_integer_displs
    module procedure mmpi_gatherv_real4
    module procedure mmpi_gatherv_real4_displs
    module procedure mmpi_gatherv_real8
    module procedure mmpi_gatherv_real8_displs
  end interface mmpi_gatherv

  ! general interface for mpi_reduce
  interface mmpi_reduce
    module procedure mmpi_reduce_integer
    module procedure mmpi_reduce_real8
  end interface mmpi_reduce

  ! general interface for mpi_scatterv
  interface mmpi_scatterv
    module procedure mmpi_scatterv_real4
    module procedure mmpi_scatterv_real8
  end interface mmpi_scatterv

  ! general interface for mpi_send
  interface mmpi_send
    module procedure mmpi_send_real8
  end interface mmpi_send

  ! general interface for mpi_recv
  interface mmpi_recv
    module procedure mmpi_recv_real8
  end interface mmpi_recv

  ! general interface for mpi_sendrecv
  interface mmpi_sendrecv
    module procedure mmpi_sendrecv_real8
  end interface mmpi_sendrecv

  ! general interface for private function 'handleLength'
  interface handleLength
    module procedure handleLengthOnly
    module procedure handleLengthWithRespectToCommunicator
  end interface handleLength
contains

  !--------------------------------------------------------------------------
  ! mmpi_initialize
  !--------------------------------------------------------------------------
  subroutine mmpi_initialize()
    !
    !:Purpose: Initialize MPI, including special communicators and
    !          several useful public variables.
    !
    implicit none

    ! Locals:
    integer :: mythread, numthread
    integer :: omp_get_thread_num, omp_get_num_threads
    integer :: rpn_comm_mype
    integer :: ierr, numNodeMasters
    integer(kind=MPI_ADDRESS_KIND) :: maxTagValue
    integer, allocatable :: allMyidHost(:)
    logical :: flag

    ! read MPI topology in namelist 'ptopo_nml'
    ! will initialize 'mmpi_npex', 'mmpi_npey' and 'mmpi_nprocs'
    call mmpi_getptopo

    ! Initialize MPI
    !call mpi_init(ierr)
    !call handleMpiError(ierr, 'Error when calling MPI_INIT in ''mmpi_initialize''')

    ! We need to call 'rpn_comm_init' because there is a call to
    ! 'RPN_COMM_xch_halo_8' and 'RPN_COMM_adj_halo8' in
    ! 'lamAnalysisGridTransforms_mod'.
    call rpn_comm_init(mmpi_getptopo, mmpi_myid, mmpi_nprocs, mmpi_npex, mmpi_npey)
    call handleMpiError(ierr, 'Error when calling RPN_COMM_INIT in ''mmpi_initialize''')

    ! get rank as 'mmpi_myid'
    call mpi_comm_rank(mpi_comm_world, mmpi_myid, ierr)
    call handleMpiError(ierr, 'Error when calling MPI_COMM_RANK for global communicator in ''mmpi_initialize''')

    ! this is a special mpi communicator for using shared memory arrays
    call mpi_comm_split_type(mpi_comm_world, mpi_comm_type_shared, 0,  &
                             mpi_info_null,  mmpi_mpicomm_SHARED,  ierr)
    call handleMpiError(ierr, 'Error when calling MPI_COMM_SPLIT_TYPE for shared communicator in ''mmpi_initialize''')

    if(mmpi_nprocs.lt.1) then
      mmpi_nprocs=1
      mmpi_npex=1
      mmpi_npey=1
      mmpi_myid=0
      mmpi_myidHost=0
      mmpi_myidx=0
      mmpi_myidy=0
    else
      !mmpi_myidx = mod(mmpi_myid,mmpi_npex)
      !mmpi_myidy = (mmpi_myid-mmpi_myidx)/mmpi_npey
      ierr = rpn_comm_mype(mmpi_myid,mmpi_myidx,mmpi_myidy) ! this routine always return success

      call mpi_comm_rank(mmpi_mpicomm_shared, mmpi_myidHost, ierr)
      call handleMpiError(ierr, 'Error when calling MPI_COMM_RANK for shared communicator in ''mmpi_initialize''')
    endif

    write(*,*) 'mmpi_initialize: mmpi_myid, mmpi_myidx, mmpi_myidy, mmpi_myidHost = ', &
                                 mmpi_myid, mmpi_myidx, mmpi_myidy, mmpi_myidHost

    ! create some MPI communicators to facilitate

    ! mmpi_comm_GRID = rpn_comm_comm(mmpi_rpn_comm_grid)
    ! Since, we are only considering a single grid, we can assume here that
    mmpi_comm_GRID = mpi_comm_world

    ! Initiliazing the 'EW' communicator, with RPN_COMM, it used to be the command
    !       mmpi_comm_EW  = rpn_comm_comm('EW')
    mmpi_comm_EW = mpi_comm_null
    call mpi_comm_split(mmpi_comm_GRID, mmpi_myidy+1, mmpi_myidx+1, mmpi_comm_EW, ierr)
    call handleMpiError(ierr, 'Error when calling  MPI_COMM_SPLIT for ''EW'' communicator in ''mmpi_initialize''')

    ! Initiliazing the 'NS' communicator, with RPN_COMM, it used to be the command
    !        mmpi_comm_NS = rpn_comm_comm('NS')
    mmpi_comm_NS = mpi_comm_null
    call mpi_comm_split(mmpi_comm_GRID, mmpi_myidx+1, mmpi_myidy+1, mmpi_comm_NS, ierr)
    call handleMpiError(ierr, 'Error when calling  MPI_COMM_SPLIT for ''NS'' communicator in ''mmpi_initialize''')

    ! Determine list of node masters (i.e. first task on each node)
    allocate(allMyidHost(mmpi_nprocs))
    call mpi_allgather(mmpi_myidHost, 1, mpi_integer, &
                       allMyidHost,   1, mpi_integer, mmpi_comm_GRID, ierr)
    numNodeMasters = count(allMyidHost(:) == 0)
    allocate(mmpi_nodeMasters(numNodeMasters))
    mmpi_nodeMasters = utl_findlocs(allMyidHost,0) - 1
    deallocate(allMyidHost)
    write(*,*) 'mmpi_initialize: mmpi_nodeMasters = ', mmpi_nodeMasters(:)

    !$OMP PARALLEL PRIVATE(numthread,mythread)
    mythread=omp_get_thread_num()
    numthread=omp_get_num_threads()
    if(mythread.eq.0) then
      write(*,*) 'mmpi_initialize: NUMBER OF THREADS=',numthread
      mmpi_numthread = numthread
    end if
    !$OMP END PARALLEL

    ! get some other useful values
    call mpi_comm_get_attr(mpi_comm_world, mpi_tag_ub, maxTagValue, flag, ierr)
    if (flag) then
      if (mmpi_myid == 0) write(*,*) 'mmpi_initialize: Maximum mpi tag value = ', maxTagValue
    else
      call utl_abort('mmpi_initialize: Could not obtain maximum tag value')
    end if

    mmpi_maxTagValue = int(maxTagValue)

    write(*,*) ' '
    if(mmpi_doBarrier) then
      write(*,*) 'mmpi_initialize: MPI_BARRIERs will be done to help with interpretation of timings'
    else
      write(*,*) 'mmpi_initialize: no MPI_BARRIERs will be done'
    endif
    write(*,*) ' '

  end subroutine mmpi_initialize

  !--------------------------------------------------------------------------
  ! mmpi_finalize
  !--------------------------------------------------------------------------
  subroutine mmpi_finalize
    !
    !:Purpose: Finalize the MPI mode to terminate the MPI program gracefully.
    !
    implicit none

    ! Locals:
    integer :: ierr

    ! We need to call 'rpn_comm_finalize' because there was a
    ! 'rpn_comm_init' in 'mmpi_initialize'
    ! call mpi_finalize(ierr)
    call rpn_comm_finalize(ierr)

    call handleMpiError(ierr, 'mmpi_finalize')

  end subroutine mmpi_finalize

  !--------------------------------------------------------------------------
  ! mmpi_barrier
  !--------------------------------------------------------------------------
  subroutine mmpi_barrier(communicator_opt)
    !
    !:Purpose: Execute 'mpi_barrier' while catching any error that may be raised.
    !
    implicit none

    ! Arguments:
    type(mpi_comm), optional, intent(in)  :: communicator_opt ! string identifying the MPI communicator

    ! Locals:
    integer :: ierr, nameLength
    type(mpi_comm) :: communicator
    character(len=MPI_MAX_OBJECT_NAME) :: commName

    if (mmpi_doBarrier) then
      communicator = handleCommunicator(communicator_opt)

      call mpi_barrier(communicator, ierr)

      if (ierr /= 0) then
        call mpi_comm_get_name(communicator, commName, nameLength)
        call handleMpiError(ierr, 'mmpi_barrier for communicator ''' // commName // '''')
      end if
    end if

  end subroutine mmpi_barrier

  !--------------------------------------------------------------------------
  ! mmpi_stopAndWait4Debug
  !--------------------------------------------------------------------------
  subroutine mmpi_stopAndWait4Debug(message)
    !
    !:Purpose: Stop the execution for the process reaching a call to the
    !          subroutine, then wait until all MPI processes reached such a
    !          call to mmpi_stopAndWait4Debug.
    !          Intended **for debugging puposes only** since it can cause
    !          unwanted MPI deadlocks - processes waiting infinitely because
    !          not all MPI processes will ever reach a call to
    !          mmpi_stopAndWait4Debug.
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: message ! input message to appear in the listing before waiting

    ! Locals:
    integer :: ierr
    type(mpi_comm), parameter :: communicator = MPI_COMM_WORLD

    write(6,9000) message
9000 format(//,4X,"!!!---ALL STOP---!!!",/,8X,"Debugging message: ",A)
    flush(6)

    call mmpi_barrier(communicator)

    call mpi_abort(communicator, 1, ierr)

  end subroutine mmpi_stopAndWait4Debug

  !--------------------------------------------------------------------------
  ! mmpi_getptopo
  !--------------------------------------------------------------------------
  subroutine mmpi_getptopo
    !
    !:Purpose: Read the file 'ptopo_nml' to get the value of 'npex'
    !          and 'npey' in 'PTOPO' namelist
    !
    implicit none

    ! Locals:
    integer :: ierr
    integer :: nulnam, fnom, fclos

    ! Namelist variables
    integer :: npex  ! number of MPI tasks in 'x' direction (set automatically by launch script)
    integer :: npey  ! number of MPI tasks in 'y' direction (set automatically by launch script)
    namelist /ptopo/ npex, npey

    npex=1
    npey=1

    call utl_tmg_start(181,'low-level--readNML')

    nulnam=0
    ierr=fnom(nulnam,'ptopo_nml','FTN+SEQ+R/O',0)

    if(ierr.ne.0) call utl_abort('mpi_getptopo: Error opening file ptopo_nml')
    read(nulnam,nml=ptopo,iostat=ierr)
    if(ierr.ne.0) call utl_abort('mpi_getptopo: Error reading namelist')

    write(*,nml=ptopo)

    ierr=fclos(nulnam)

    call utl_tmg_stop(181)

    mmpi_npex   = npex
    mmpi_npey   = npey
    mmpi_nprocs = npex*npey

  end subroutine mmpi_getptopo

  !--------------------------------------------------------------------------
  ! mmpi_allreduce_sumreal8scalar
  !--------------------------------------------------------------------------
  subroutine mmpi_allreduce_sumreal8scalar(sendRecvValue)
    !
    !:Purpose: Version of mpi_allReduce that always performs sum in
    !          the same order.
    !
    implicit none

    ! Arguments:
    real(8), intent(inout) :: sendRecvValue ! value to be summed over all mpi tasks

    ! Locals:
    integer :: root
    real(8), allocatable :: allvalues(:)

    ! do a barrier so that timing on reduce operation is accurate
    call utl_tmg_start(171,'low-level--mpi_allreduce_barr')
    call mmpi_barrier
    call utl_tmg_stop(171)

    call utl_tmg_start(170,'low-level--mpi_allreduce_sum8')

    call mmpi_allreduce(mmpi_myid, root, MPI_MIN)

    ! gather values to be added onto 1 processor
    allocate(allvalues(mmpi_nprocs))
    call mmpi_gather(sendRecvValue, allvalues, procID_opt = root)

    ! sum the values on the "root" mpi task and broadcast to group
    if( mmpi_myid == root ) sendRecvValue = sum(allvalues(:))
    deallocate(allvalues)
    call mmpi_bcast(sendRecvValue, procID_opt = root)

    call utl_tmg_stop(170)

  end subroutine mmpi_allreduce_sumreal8scalar

  !--------------------------------------------------------------------------
  ! mmpi_allReduce_sumR8_1d
  !--------------------------------------------------------------------------
  subroutine mmpi_allreduce_sumR8_1d(sendRecvVector)
    !
    ! :Purpose: Perform sum of 1d array over all MPI tasks, guaranteed to
    !           always be in the same order.
    !
    implicit none

    ! Arguments:
    real(8), intent(inout)  :: sendRecvVector(:) ! 1-D vector to be summed over all mpi tasks

    ! Locals:
    integer :: numElements, root
    real(8), allocatable :: all_sendRecvVector(:,:)

    ! do a barrier so that timing on reduce operation is accurate
    call utl_tmg_start(171,'low-level--mpi_allreduce_barr')
    call mmpi_barrier
    call utl_tmg_stop(171)

    call utl_tmg_start(170,'low-level--mpi_allreduce_sum8')

    numElements = size(sendRecvVector)

    call mmpi_allreduce(mmpi_myid, root, MPI_MIN)

    ! gather vectors to be added onto 1 processor
    allocate(all_sendRecvVector(numElements,0:mmpi_nprocs-1))
    call mmpi_gather(sendRecvVector, all_sendRecvVector, procID_opt = root)

    ! sum the values on the "root" mpi task and broadcast to group
    if ( mmpi_myid == root ) sendRecvVector(:) = sum(all_sendRecvVector(:,:),2)
    deallocate(all_sendRecvVector)

    call mmpi_bcast(sendRecvVector, procID_opt = root)

    call utl_tmg_stop(170)

  end subroutine mmpi_allreduce_sumR8_1d

  !--------------------------------------------------------------------------
  ! mmpi_allreduce_sumR8_2d
  !--------------------------------------------------------------------------
  subroutine mmpi_allreduce_sumR8_2d(sendRecvVector)
    !
    ! :Purpose: Perform sum of 2d array over all MPI tasks guaranteed
    !           to always be in the same order.
    !
    implicit none

    ! Arguments:
    real(8), intent(inout)  :: sendRecvVector(:,:) ! 2-D vector to be summed over all mpi tasks

    ! Locals:
    integer :: numElements1, numElements2, root
    real(8), allocatable :: all_sendRecvVector(:,:,:)

    ! do a barrier so that timing on reduce operation is accurate
    call utl_tmg_start(171,'low-level--mpi_allreduce_barr')
    call mmpi_barrier
    call utl_tmg_stop(171)

    call utl_tmg_start(170,'low-level--mpi_allreduce_sum8')

    numElements1 = size(sendRecvVector,1)
    numElements2 = size(sendRecvVector,2)

    call mmpi_allreduce(mmpi_myid, root, MPI_MIN)

    ! gather vectors to be added onto 1 processor
    allocate(all_sendRecvVector(numElements1,numElements2,0:mmpi_nprocs-1))
    call mmpi_gather(sendRecvVector, all_sendRecvVector, procID_opt = root)

    ! sum the values on the "root" mpi task and broadcast to group
    if ( mmpi_myid == root ) sendRecvVector(:,:) = sum(all_sendRecvVector(:,:,:),3)
    deallocate(all_sendRecvVector)

    call mmpi_bcast(sendRecvVector, procID_opt = root)

    call utl_tmg_stop(170)

  end subroutine mmpi_allreduce_sumR8_2d

  !--------------------------------------------------------------------------
  ! mmpi_reduce_sumR8_1d
  !--------------------------------------------------------------------------
  subroutine mmpi_reduce_sumR8_1d(sendVector, recvVector)
    !
    ! :Purpose: Perform sum of 1d array over all MPI tasks guaranteed to
    !           always be in the same order.
    !
    implicit none

    ! Arguments:
    real(8), intent(in)  :: sendVector(:) ! 1-D vector to be summed over all mpi tasks
    real(8), intent(out) :: recvVector(:) ! 1-D vector to be summed over all mpi tasks

    ! Locals:
    integer, parameter :: ROOT = 0
    integer :: numElements
    real(8), allocatable :: all_sendRecvVector(:,:)

    ! do a barrier so that timing on reduce operation is accurate
    call utl_tmg_start(171,'low-level--mpi_allreduce_barr')
    call mmpi_barrier
    call utl_tmg_stop(171)

    call utl_tmg_start(170,'low-level--mpi_allreduce_sum8')

    numElements = size(sendVector)

    ! gather vectors to be added onto 1 processor
    if ( mmpi_myid == ROOT ) then
      allocate(all_sendRecvVector(numElements,0:mmpi_nprocs-1))
    else
      allocate(all_sendRecvVector(1,1))
    end if
    call mmpi_gather(sendVector, all_sendRecvVector, procID_opt = ROOT)

    ! sum the values on the "root" mpi task
    if ( mmpi_myid == ROOT ) recvVector(:) = sum(all_sendRecvVector(:,:),2)
    deallocate(all_sendRecvVector)

    call utl_tmg_stop(170)

  end subroutine mmpi_reduce_sumR8_1d

  !--------------------------------------------------------------------------
  ! mmpi_reduce_sumR8_2d
  !--------------------------------------------------------------------------
  subroutine mmpi_reduce_sumR8_2d(sendVector, recvVector)
    !
    ! :Purpose: Perform sum of 2d array over all MPI tasks guaranteed to
    !           always be in the same order.
    !
    implicit none

    ! Arguments:
    real(8), intent(in)  :: sendVector(:,:) ! 2-D vector to be summed over all mpi tasks
    real(8), intent(out) :: recvVector(:,:) ! 2-D vector to be summed over all mpi tasks

    ! Locals:
    integer, parameter :: ROOT = 0
    integer :: numElements1, numElements2
    real(8), allocatable :: all_sendRecvVector(:,:,:)

    ! do a barrier so that timing on reduce operation is accurate
    call utl_tmg_start(171,'low-level--mpi_allreduce_barr')
    call mmpi_barrier
    call utl_tmg_stop(171)

    call utl_tmg_start(170,'low-level--mpi_allreduce_sum8')

    numElements1 = size(sendVector,1)
    numElements2 = size(sendVector,2)

    ! gather vectors to be added onto 1 processor
    if ( mmpi_myid == ROOT ) then
      allocate(all_sendRecvVector(numElements1,numElements2,0:mmpi_nprocs-1))
    else
      allocate(all_sendRecvVector(1,1,1))
    end if
    call mmpi_gather(sendVector, all_sendRecvVector, procID_opt = ROOT)

    ! sum the values on the "root" mpi task
    if ( mmpi_myid == ROOT ) recvVector(:,:) = sum(all_sendRecvVector(:,:,:),3)
    deallocate(all_sendRecvVector)

    call utl_tmg_stop(170)

  end subroutine mmpi_reduce_sumR8_2d

  !--------------------------------------------------------------------------
  ! mmpi_reduce_sumR8_3d
  !--------------------------------------------------------------------------
  subroutine mmpi_reduce_sumR8_3d(sendVector, recvVector)
    !
    ! :Purpose: Perform sum of 3d array over all MPI tasks guaranteed to
    !           always be in the same order.
    !
    implicit none

    ! Arguments:
    real(8), intent(in)  :: sendVector(:,:,:) ! 3-D vector to be summed over all mpi tasks
    real(8), intent(out) :: recvVector(:,:,:) ! 3-D vector to be summed over all mpi tasks

    ! Locals:
    integer, parameter :: ROOT = 0
    integer :: numElements1, numElements2, numElements3
    real(8), allocatable :: all_sendRecvVector(:,:,:,:)

    ! do a barrier so that timing on reduce operation is accurate
    call utl_tmg_start(171,'low-level--mpi_allreduce_barr')
    call mmpi_barrier
    call utl_tmg_stop(171)

    call utl_tmg_start(170,'low-level--mpi_allreduce_sum8')

    numElements1 = size(sendVector,1)
    numElements2 = size(sendVector,2)
    numElements3 = size(sendVector,3)

    ! gather vectors to be added onto 1 processor
    if ( mmpi_myid == ROOT ) then
      allocate(all_sendRecvVector(numElements1,numElements2,numElements3,0:mmpi_nprocs-1))
    else
      allocate(all_sendRecvVector(1,1,1,1))
    end if
    call mmpi_gather(sendVector, all_sendRecvVector, procID_opt = ROOT)

    ! sum the values on the "root" mpi task
    if ( mmpi_myid == ROOT ) recvVector(:,:,:) = sum(all_sendRecvVector(:,:,:,:),4)
    deallocate(all_sendRecvVector)

    call utl_tmg_stop(170)

  end subroutine mmpi_reduce_sumR8_3d

  !--------------------------------------------------------------------------
  ! mmpi_allgather_string
  !--------------------------------------------------------------------------
  subroutine mmpi_allgather_string(str_list, str_list_all, nlist, nchar)
    !
    ! :Purpose: Performs the MPI 'allgather' routine for an array of strings
    !
    implicit none

    ! Arguments:
    integer,              intent(in)  :: nlist
    integer,              intent(in)  :: nchar
    character(len=nchar), intent(in)  :: str_list(nlist)
    character(len=nchar), intent(out) :: str_list_all(nlist,mmpi_nprocs)

    ! Locals:
    integer :: num_list(nlist*nchar),num_list_all(nlist*nchar,mmpi_nprocs)
    integer :: ierr,ilist,ichar,iproc

    ! Convert strings to integer sequences
    do ilist=1,nlist
       do ichar=1,nchar
          num_list((ilist-1)*nchar+ichar) = iachar(str_list(ilist)(ichar:ichar))
       end do
    end do

    ! Perform allgather with converted integer sequences
    call mpi_allgather(num_list,     nlist*nchar, MPI_INTEGER, &
                       num_list_all, nlist*nchar, MPI_INTEGER, mmpi_comm_GRID, ierr)
    call handleMpiError(ierr, 'mmpi_allgather_string')

    ! Convert integer sequences to stnid character strings
    do iproc=1,mmpi_nprocs
       do ilist=1,nlist
          do ichar=1,nchar
             str_list_all(ilist,iproc)(ichar:ichar) = achar(num_list_all((ilist-1)*nchar+ichar,iproc))
          end do
       end do
    end do

  end subroutine mmpi_allgather_string

  !--------------------------------------------------------------------------
  ! mmpi_setup_latbands
  !--------------------------------------------------------------------------
  subroutine mmpi_setup_latbands(nj, latPerPE, latPerPEmax, myLatBeg, myLatEnd,  &
                                 myLatHalfBeg_opt, myLatHalfEnd_opt, divisible_opt)
    !
    !:Purpose: Compute parameters that define the mpi distribution of
    !          latitudes over tasks in Y direction (npey)
    !
    implicit none

    ! Arguments:
    integer          , intent(in)  :: nj
    integer          , intent(out) :: latPerPE
    integer          , intent(out) :: latPerPEmax
    integer          , intent(out) :: myLatBeg
    integer          , intent(out) :: myLatEnd
    integer, optional, intent(out) :: myLatHalfBeg_opt
    integer, optional, intent(out) :: myLatHalfEnd_opt
    logical, optional, intent(out) :: divisible_opt

    ! Locals:
    integer :: latPerPEmin, njlath, ierr
    logical, save :: firstCall = .true.

    latPerPEmin = floor(real(nj) / real(mmpi_npey))
    myLatBeg = 1 + (mmpi_myidy * latPerPEmin)
    if( mmpi_myidy < (mmpi_npey-1) ) then
      myLatEnd = (1 + mmpi_myidy) * latPerPEmin
    else
      myLatEnd = nj
    end if
    latPerPE = myLatEnd - myLatBeg + 1

    call mpi_allReduce(latPerPE,latPerPEmax,1, mpi_integer, mpi_max, &
                       mmpi_comm_NS, ierr)

    call handleMpiError(ierr, 'mmpi_setup_lonbands')

    if( firstCall ) then
      write(*,'(a,4i8)') 'mmpi_setup_latbands: latPerPE, latPerPEmax, myLatBeg, myLatEnd = ',  &
           latPerPE, latPerPEmax, myLatBeg, myLatEnd
      firstCall = .false.
    end if

    if (present(myLatHalfBeg_opt).and.present(myLatHalfEnd_opt)) then
      njlath = (nj + 1) / 2
      if (myLatBeg <= njlath .and. myLatEnd <= njlath) then
        myLatHalfBeg_opt = myLatBeg
        myLatHalfEnd_opt = myLatEnd
      elseif (myLatBeg >= njlath .and. myLatEnd >= njlath) then
        myLatHalfBeg_opt = 1 + nj - myLatEnd
        myLatHalfEnd_opt = 1 + nj - myLatBeg
      else
        myLatHalfBeg_opt = min(myLatBeg, 1 + nj - myLatEnd)
        myLatHalfEnd_opt = njlath
      end if
    end if

    if( present(divisible_opt) ) then
      divisible_opt = (latPerPEmin * mmpi_npey == nj)
    end if

  end subroutine mmpi_setup_latbands

  !--------------------------------------------------------------------------
  ! mmpi_myidYfromLat
  !--------------------------------------------------------------------------
  function mmpi_myidYfromLat(latIndex, nj) result(IP_y)
    !
    !:Purpose: Use same logic as setup_latbands to compute myidy
    !          corresponding to a latitude grid index
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: latIndex
    integer, intent(in) :: nj
    ! Result:
    integer :: IP_y

    IP_y = (latIndex-1) / floor( real(nj) / real(mmpi_npey) )
    IP_y = min( mmpi_npey-1, IP_y )

  end function mmpi_myidYfromLat

  !--------------------------------------------------------------------------
  ! mmpi_setup_lonbands
  !--------------------------------------------------------------------------
  subroutine mmpi_setup_lonbands(ni, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd, divisible_opt)
    !
    !:Purpose: Compute parameters that define the MPI distribution of
    !          longitudes over tasks in X direction (npex)
    !
    implicit none

    ! Arguments:
    integer          , intent(in)  :: ni
    integer          , intent(out) :: lonPerPE
    integer          , intent(out) :: lonPerPEmax
    integer          , intent(out) :: myLonBeg
    integer          , intent(out) :: myLonEnd
    logical, optional, intent(out) :: divisible_opt

    ! Locals:
    integer :: lonPerPEmin, ierr
    logical, save :: firstCall = .true.

    lonPerPEmin = floor(real(ni) / real(mmpi_npex))
    myLonBeg = 1 + (mmpi_myidx * lonPerPEmin)
    if( mmpi_myidx < (mmpi_npex-1) ) then
      myLonEnd = (1 + mmpi_myidx) * lonPerPEmin
    else
      myLonEnd = ni
    end if
    lonPerPE = myLonEnd - myLonBeg + 1

    call mpi_allReduce(lonPerPE,lonPerPEmax,1, mpi_integer, mpi_max, &
                       mmpi_comm_EW, ierr)

    call handleMpiError(ierr, 'mmpi_setup_lonbands')

    if( firstCall ) then
      write(*,'(a,4i8)') 'mmpi_setup_lonbands: lonPerPE, lonPerPEmax, myLonBeg, myLonEnd = ', &
           lonPerPE, lonPerPEmax, myLonBeg, myLonEnd
      firstCall = .false.
    end if

    if( present(divisible_opt) ) then
      divisible_opt = (lonPerPEmin * mmpi_npex == ni)
    end if

  end subroutine mmpi_setup_lonbands

  !--------------------------------------------------------------------------
  ! mmpi_myidXfromLon
  !--------------------------------------------------------------------------
  function mmpi_myidXfromLon(lonIndex, ni) result(IP_x)
    !
    !:Purpose: Use same logic as setup_lonbands to compute myidx
    !          corresponding to a longitude grid index
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: lonIndex
    integer, intent(in) :: ni
    ! Result:
    integer :: IP_x

    IP_x = (lonIndex-1) / floor( real(ni) / real(mmpi_npex) )
    IP_x = min( mmpi_npex-1, IP_x )

  end function mmpi_myidXfromLon

  !--------------------------------------------------------------------------
  ! mmpi_setup_m
  !--------------------------------------------------------------------------
  subroutine mmpi_setup_m(ntrunc, mymBeg, mymEnd, mymSkip, mymCount)
    !
    !:Purpose: Compute parameters that define the mpi distribution of
    !          wavenumber m over tasks in Y direction (npey)
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: ntrunc
    integer, intent(out) :: mymBeg
    integer, intent(out) :: mymEnd
    integer, intent(out) :: mymSkip
    integer, intent(out) :: mymCount

    ! Locals:
    integer :: jm

    mymBeg = mmpi_myidy
    mymEnd = ntrunc
    mymSkip = mmpi_npey
    mymCount = 0
    do jm = mymBeg, mymEnd, mymSkip
      mymCount = mymCount + 1
    end do

    write(*,'(a,4i8)') 'mmpi_setup_m: mymBeg, mymEnd, mymSkip, mymCount = ', mymBeg, mymEnd, mymSkip, mymCount

  end subroutine mmpi_setup_m

  !--------------------------------------------------------------------------
  ! mmpi_setup_n
  !--------------------------------------------------------------------------
  subroutine mmpi_setup_n(ntrunc, mynBeg, mynEnd, mynSkip, mynCount)
    !
    !:Purpose: Compute parameters that define the mpi distribution of
    !          wavenumber n over tasks in X direction (npex)
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: ntrunc
    integer, intent(out) :: mynBeg
    integer, intent(out) :: mynEnd
    integer, intent(out) :: mynSkip
    integer, intent(out) :: mynCount

    ! Locals:
    integer :: jn

    mynBeg = mmpi_myidx
    mynEnd = ntrunc
    mynSkip = mmpi_npex
    mynCount = 0
    do jn = mynBeg, mynEnd, mynSkip
      mynCount = mynCount + 1
    end do

    write(*,'(a,4i8)') 'mmpi_setup_n: mynBeg, mynEnd, mynSkip, mynCount = ', mynBeg, mynEnd, mynSkip, mynCount

  end subroutine mmpi_setup_n

  !--------------------------------------------------------------------------
  ! mmpi_setup_levels
  !--------------------------------------------------------------------------
  subroutine mmpi_setup_levels(numlevels, myLevBeg, myLevEnd, myLevCount)
    !
    !:Purpose: Compute parameters that define the mpi distribution of
    !          levels over tasks in X direction (npex)
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: numlevels
    integer, intent(out) :: myLevBeg
    integer, intent(out) :: myLevEnd
    integer, intent(out) :: myLevCount

    ! Locals:
    integer :: jlev
    integer :: procIndex
    integer :: factor
    integer :: myLevCounts(mmpi_npex)

    ! when possible, always divide into even number of levels per MPI task
    if(mod(numlevels, 2) /= 0) then
      write(*,*) 'mmpi_setup_levels: total number of levels is not even, now=', numlevels
      write(*,*) '                   therefore, if global grid, may not be able to do '
      write(*,*) '                   transforms of vor/div <-> u/v'
      factor = 1
    else
      factor = 2
    end if

    myLevCounts(:) = 0
    do procIndex = 1, mmpi_npex
      do jlev = procIndex, (numlevels / factor), mmpi_npex
        myLevCounts(procIndex) = myLevCounts(procIndex) + 1
      end do
    end do
    do procIndex = 1, mmpi_npex
      myLevCounts(procIndex) = myLevCounts(procIndex) * factor
    end do

    myLevCount = myLevCounts(mmpi_myidx + 1)

    if (myLevCount > 0) then
      myLevBeg = 1
      do procIndex = 1, mmpi_myidx
        myLevBeg = myLevBeg + myLevCounts(procIndex)
      end do
      myLevEnd = myLevBeg + myLevCount - 1
    else
      myLevBeg = 1
      myLevEnd = 0
    end if

    write(*,'(a,3i8)') 'mmpi_setup_levels: myLevBeg, myLevEnd, myLevCount = ',  &
         myLevBeg, myLevEnd, myLevCount

  end subroutine mmpi_setup_levels

  !--------------------------------------------------------------------------
  ! mmpi_setup_varslevels
  !--------------------------------------------------------------------------
  subroutine mmpi_setup_varslevels(numVarLev, myVarLevBeg, myVarLevEnd, myVarLevCount)
    !
    !:Purpose: Compute parameters that define the mpi distribution of
    !          variables/levels (i.e. 1->nk) over all tasks (nprocs)
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: numVarLev
    integer, intent(out) :: myVarLevBeg
    integer, intent(out) :: myVarLevEnd
    integer, intent(out) :: myVarLevCount

    ! Locals:
    integer :: varLevIndex
    integer :: procIndex
    integer :: myVarLevCounts(mmpi_nprocs)

    myVarLevCounts(:) = 0
    do procIndex = 1, mmpi_nprocs
      do varLevIndex = procIndex, numVarLev, mmpi_nprocs
        myVarLevCounts(procIndex) = myVarLevCounts(procIndex) + 1
      end do
    end do

    myVarLevCount = myVarLevCounts(mmpi_myid + 1)

    myVarLevBeg = 1
    do procIndex = 1, mmpi_myid
      myVarLevBeg = myVarLevBeg + myVarLevCounts(procIndex)
    end do
    myVarLevEnd = myVarLevBeg + myVarLevCount - 1

    write(*,'(a,3i8)') 'mmpi_setup_varslevels: myVarLevBeg, myVarLevEnd, myVarLevCount = ',  &
         myVarLevBeg, myVarLevEnd, myVarLevCount

  end subroutine mmpi_setup_varslevels

  !--------------------------------------------------------------------------
  ! mmpi_bcast_character
  !--------------------------------------------------------------------------
  subroutine mmpi_bcast_character(charData, procID_opt)
    !
    !:Purpose: Calling 'mpi_bcast' for character array
    !
    implicit none

    ! Arguments:
    character(len=*),  intent(inout) :: charData   ! Input character to broadcast
    integer, optional, intent(in)    :: procID_opt ! MPI rank to broadcast from

    ! Locals:
    integer :: ierr
    integer :: procID

    procID = handleProcID(procID_opt)

    call mpi_bcast(charData, len(charData), MPI_CHARACTER, procID, mmpi_comm_GRID, ierr)

    call handleMpiError(ierr, 'mmpi_bcast_character')

  end subroutine mmpi_bcast_character

  !--------------------------------------------------------------------------
  ! mmpi_bcast_logical
  !--------------------------------------------------------------------------
  subroutine mmpi_bcast_logical(logicalData, length_opt, procID_opt)
    !
    !:Purpose: Calling 'mpi_bcast' for a logical value or array
    !
    implicit none

    ! Arguments:
    logical,           intent(inout) :: logicalData(..) ! Input logical data to broadcast
    integer, optional, intent(in)    :: length_opt      ! size of the input data
    integer, optional, intent(in)    :: procID_opt      ! MPI rank to broadcast from

    ! Locals:
    integer :: ierr
    integer :: procID, length

    length = handleLength(logicalData, length_opt)
    procID = handleProcID(procID_opt)

    call mpi_bcast(logicalData, length, MPI_LOGICAL, procID, mmpi_comm_GRID, ierr)

    call handleMpiError(ierr, 'mmpi_bcast_logical')

  end subroutine mmpi_bcast_logical

  !--------------------------------------------------------------------------
  ! mmpi_bcast_integer
  !--------------------------------------------------------------------------
  subroutine mmpi_bcast_integer(integerData, length_opt, procID_opt)
    !
    !:Purpose: Calling 'mpi_bcast' for an integer scalar or array
    !
    implicit none

    ! Arguments:
    integer, contiguous, intent(inout) :: integerData(..) ! Input integer data to broadcast
    integer, optional,   intent(in)    :: length_opt      ! size of the input data
    integer, optional,   intent(in)    :: procID_opt      ! MPI rank to broadcast from

    ! Locals:
    integer :: ierr
    integer :: procID, length

    length = handleLength(integerData, length_opt)
    procID = handleProcID(procID_opt)

    call mpi_bcast(integerData, length, MPI_INTEGER, procID, mmpi_comm_GRID, ierr)

    call handleMpiError(ierr, 'mmpi_bcast_integer')

  end subroutine mmpi_bcast_integer

  !--------------------------------------------------------------------------
  ! mmpi_bcast_real4
  !--------------------------------------------------------------------------
  subroutine mmpi_bcast_real4(real4Data, length_opt, procID_opt)
    !
    !:Purpose: Calling 'mpi_bcast' for a real(4) scalar or array
    !
    implicit none

    ! Arguments:
    real(4), contiguous, intent(inout) :: real4Data(..) ! Input real(4) data to broadcast
    integer, optional,   intent(in)    :: length_opt    ! size of the input data
    integer, optional,   intent(in)    :: procID_opt    ! MPI rank to broadcast from

    ! Locals:
    integer :: ierr
    integer :: procID, length

    length = handleLength(real4Data, length_opt)
    procID = handleProcID(procID_opt)

    call mpi_bcast(real4Data, length, MPI_REAL4, procID, mmpi_comm_GRID, ierr)

    call handleMpiError(ierr, 'mmpi_bcast_real4')

  end subroutine mmpi_bcast_real4

  !--------------------------------------------------------------------------
  ! mmpi_bcast_real8
  !--------------------------------------------------------------------------
  subroutine mmpi_bcast_real8(real8Data, length_opt, procID_opt)
    !
    !:Purpose: Calling 'mpi_bcast' for a real(8) array
    !
    implicit none

    ! Arguments:
    real(8), contiguous, intent(inout) :: real8Data(..) ! Input real(8) data to broadcast
    integer, optional,   intent(in)    :: length_opt    ! size of the input data
    integer, optional,   intent(in)    :: procID_opt    ! MPI rank to broadcast from

    ! Locals:
    integer :: ierr
    integer :: procID, length

    length = handleLength(real8Data, length_opt)
    procID = handleProcID(procID_opt)

    call mpi_bcast(real8Data, length, MPI_REAL8, procID, mmpi_comm_GRID, ierr)

    call handleMpiError(ierr, 'mmpi_bcast_real8')

  end subroutine mmpi_bcast_real8

  !--------------------------------------------------------------------------
  ! mmpi_gather_logical
  !--------------------------------------------------------------------------
  subroutine mmpi_gather_logical(sending, receiving, length_opt, procID_opt)
    !
    !:Purpose: Calling 'mpi_gather' for a logical scalar or array
    !
    implicit none

    ! Arguments:
    logical, contiguous, intent(in)  :: sending(..)     ! logical data sent to 'procID_opt'
    logical, contiguous, intent(out) :: receiving(..,:) ! logical array which stores the data received
    integer, optional,   intent(in)  :: length_opt      ! size of the input data
    integer, optional,   intent(in)  :: procID_opt      ! MPI rank to gather to

    ! Locals:
    integer :: ierr
    integer :: procID, length

    length = handleLength(sending, length_opt)
    procID = handleProcID(procID_opt)

    call mpi_gather(sending,   length, mpi_logical, &
                    receiving, length, mpi_logical, procID, mmpi_comm_GRID, ierr)

    call handleMpiError(ierr, 'mmpi_gather_logical')

  end subroutine mmpi_gather_logical

  !--------------------------------------------------------------------------
  ! mmpi_gather_integer
  !--------------------------------------------------------------------------
  subroutine mmpi_gather_integer(sending, receiving, length_opt, procID_opt)
    !
    !:Purpose: Calling 'mpi_gather' for an integer scalar or array
    !
    implicit none

    ! Arguments:
    integer, contiguous, intent(in)  :: sending(..)     ! integer data sent to 'procID_opt'
    integer, contiguous, intent(out) :: receiving(..,:) ! integer array which stores the data received
    integer, optional,   intent(in)  :: length_opt      ! size of the input data
    integer, optional,   intent(in)  :: procID_opt      ! MPI rank to gather to

    ! Locals:
    integer :: ierr
    integer :: procID, length

    length = handleLength(sending, length_opt)
    procID = handleProcID(procID_opt)

    call mpi_gather(sending,   length, mpi_integer, &
                    receiving, length, mpi_integer, procID, mmpi_comm_GRID, ierr)

    call handleMpiError(ierr, 'mmpi_gather_integer')

  end subroutine mmpi_gather_integer

  !--------------------------------------------------------------------------
  ! mmpi_gather_integer8
  !--------------------------------------------------------------------------
  subroutine mmpi_gather_integer8(sending, receiving, length_opt, procID_opt)
    !
    !:Purpose: Calling 'mpi_gather' for an integer(8) scalar or array
    !
    implicit none

    ! Arguments:
    integer(8), contiguous, intent(in)  :: sending(..)     ! integer(8) data sent to 'procID_opt'
    integer(8), contiguous, intent(out) :: receiving(..,:) ! integer(8) array which stores the data received
    integer, optional,      intent(in)  :: length_opt      ! size of the input data
    integer, optional,      intent(in)  :: procID_opt      ! MPI rank to gather to

    ! Locals:
    integer :: ierr
    integer :: procID, length

    length = handleLength(sending, length_opt)
    procID = handleProcID(procID_opt)

    call mpi_gather(sending,   length, mpi_integer8, &
                    receiving, length, mpi_integer8, procID, mmpi_comm_GRID, ierr)

    call handleMpiError(ierr, 'mmpi_gather_integer8')

  end subroutine mmpi_gather_integer8

  !--------------------------------------------------------------------------
  ! mmpi_gather_real4
  !--------------------------------------------------------------------------
  subroutine mmpi_gather_real4(sending, receiving, length_opt, procID_opt)
    !
    !:Purpose: Calling 'mpi_gather' for an real(4) scalar or array
    !
    implicit none

    ! Arguments:
    real(4), contiguous, intent(in)  :: sending(..)     ! real(4) data sent to 'procID_opt'
    real(4), contiguous, intent(out) :: receiving(..,:) ! real(4) array which stores the data received
    integer, optional,   intent(in)  :: length_opt      ! size of the input data
    integer, optional,   intent(in)  :: procID_opt      ! MPI rank to gather to

    ! Locals:
    integer :: ierr
    integer :: procID, length

    length = handleLength(sending, length_opt)
    procID = handleProcID(procID_opt)

    call mpi_gather(sending,   length, mpi_real4, &
                    receiving, length, mpi_real4, procID, mmpi_comm_GRID, ierr)

    call handleMpiError(ierr, 'mmpi_gather_real4')

  end subroutine mmpi_gather_real4

  !--------------------------------------------------------------------------
  ! mmpi_gather_real8
  !--------------------------------------------------------------------------
  subroutine mmpi_gather_real8(sending, receiving, length_opt, procID_opt)
    !
    !:Purpose: Calling 'mpi_gather' for an real(8) scalar or array
    !
    implicit none

    ! Arguments:
    real(8), contiguous, intent(in)  :: sending(..)     ! real(8) data sent to 'procID_opt'
    real(8), contiguous, intent(out) :: receiving(..,:) ! real(8) array which stores the data received
    integer, optional,   intent(in)  :: length_opt      ! size of the input data
    integer, optional,   intent(in)  :: procID_opt      ! MPI rank to gather to

    ! Locals:
    integer :: ierr
    integer :: procID, length

    length = handleLength(sending, length_opt)
    procID = handleProcID(procID_opt)

    call mpi_gather(sending,   length, mpi_real8, &
                    receiving, length, mpi_real8, procID, mmpi_comm_GRID, ierr)

    call handleMpiError(ierr, 'mmpi_gather_real8')

  end subroutine mmpi_gather_real8

  !--------------------------------------------------------------------------
  ! mmpi_allGather_logical
  !--------------------------------------------------------------------------
  subroutine mmpi_allGather_logical(sending, receiving, length_opt, communicator_opt)
    !
    !:Purpose: Calling 'mpi_allGather' for a logical scalar or array
    !
    implicit none

    ! Arguments:
    logical,        contiguous, intent(in)  :: sending(..)      ! logical data sent to all MPI ranks
    logical,        contiguous, intent(out) :: receiving(..,:)  ! logical array which stores the data received
    integer,        optional,   intent(in)  :: length_opt       ! size of the input data
    type(mpi_comm), optional,   intent(in)  :: communicator_opt ! the MPI communicator

    ! Locals:
    integer :: ierr, length
    type(mpi_comm) :: communicator

    length = handleLength(sending, length_opt)
    communicator = handleCommunicator(communicator_opt)

    call mpi_allGather(sending,   length, mpi_logical,  &
                       receiving, length, mpi_logical, communicator, ierr)

    call handleMpiError(ierr, 'mmpi_allGather_logical')

  end subroutine mmpi_allGather_logical

  !--------------------------------------------------------------------------
  ! mmpi_allGather_integer
  !--------------------------------------------------------------------------
  subroutine mmpi_allGather_integer(sending, receiving, length_opt, communicator_opt)
    !
    !:Purpose: Calling 'mpi_allGather' for a integer scalar or array
    !
    implicit none

    ! Arguments:
    integer,        contiguous, intent(in)  :: sending(..)      ! integer data sent to all MPI ranks
    integer,        contiguous, intent(out) :: receiving(..,:)  ! integer array which stores the data received
    integer,        optional,   intent(in)  :: length_opt       ! size of the input data
    type(mpi_comm), optional,   intent(in)  :: communicator_opt ! the MPI communicator

    ! Locals:
    integer :: ierr, length
    type(mpi_comm) :: communicator

    length = handleLength(sending, length_opt)
    communicator = handleCommunicator(communicator_opt)

    call mpi_allGather(sending,   length, mpi_integer,  &
                       receiving, length, mpi_integer, communicator, ierr)

    call handleMpiError(ierr, 'mmpi_allGather_integer')

  end subroutine mmpi_allGather_integer

  !--------------------------------------------------------------------------
  ! mmpi_allGather_real4
  !--------------------------------------------------------------------------
  subroutine mmpi_allGather_real4(sending, receiving, length_opt, communicator_opt)
    !
    !:Purpose: Calling 'mpi_allGather' for a real(4) scalar or array
    !
    implicit none

    ! Arguments:
    real(4),        contiguous, intent(in)  :: sending(..)     ! real(4) data sent to all MPI ranks
    real(4),        contiguous, intent(out) :: receiving(..,:) ! real(4) array which stores the data received
    integer,        optional,   intent(in)  :: length_opt      ! size of the input data
    type(mpi_comm), optional,   intent(in) :: communicator_opt ! the MPI communicator

    ! Locals:
    integer :: ierr, length
    type(mpi_comm) :: communicator

    length = handleLength(sending, length_opt)
    communicator = handleCommunicator(communicator_opt)

    call mpi_allGather(sending,   length, mpi_real4,  &
                       receiving, length, mpi_real4, communicator, ierr)

    call handleMpiError(ierr, 'mmpi_allGather_real4')

  end subroutine mmpi_allGather_real4

  !--------------------------------------------------------------------------
  ! mmpi_allGather_real8
  !--------------------------------------------------------------------------
  subroutine mmpi_allGather_real8(sending, receiving, length_opt, communicator_opt)
    !
    !:Purpose: Calling 'mpi_allGather' for a real(8) scalar or array
    !
    implicit none

    ! Arguments:
    real(8),        contiguous, intent(in)  :: sending(..)      ! real(8) data sent to all MPI ranks
    real(8),        contiguous, intent(out) :: receiving(..,:)  ! real(8) array which stores the data received
    integer,        optional,   intent(in)  :: length_opt       ! size of the input data
    type(mpi_comm), optional,   intent(in)  :: communicator_opt ! the MPI communicator

    ! Locals:
    integer :: ierr, length
    type(mpi_comm) :: communicator

    length = handleLength(sending, length_opt)
    communicator = handleCommunicator(communicator_opt)

    call mpi_allGather(sending,   length, mpi_real8,  &
                       receiving, length, mpi_real8, communicator, ierr)

    call handleMpiError(ierr, 'mmpi_allGather_real8')

  end subroutine mmpi_allGather_real8

  !--------------------------------------------------------------------------
  ! mmpi_alltoall_integer
  !--------------------------------------------------------------------------
  subroutine mmpi_alltoall_integer(sending, receiving, length_opt, communicator_opt)
    !
    !:Purpose: Calling 'mpi_alltoall' for a integer scalar or array
    !
    implicit none

    ! Arguments:
    integer,        contiguous, intent(in)  :: sending(..)      ! integer data sent to all MPI ranks
    integer,        contiguous, intent(out) :: receiving(..)    ! integer array which stores the data received
    integer,        optional,   intent(in)  :: length_opt       ! size of the input data
    type(mpi_comm), optional,   intent(in)  :: communicator_opt ! the MPI communicator

    ! Locals:
    integer :: ierr, length
    type(mpi_comm) :: communicator

    communicator = handleCommunicator(communicator_opt)
    length = handleLength(sending, communicator, length_opt)

    call mpi_alltoall(sending,   length, mpi_integer,  &
                      receiving, length, mpi_integer, communicator, ierr)

    call handleMpiError(ierr, 'mmpi_alltoall_integer')

  end subroutine mmpi_alltoall_integer

  !--------------------------------------------------------------------------
  ! mmpi_alltoall_integer8
  !--------------------------------------------------------------------------
  subroutine mmpi_alltoall_integer8(sending, receiving, length_opt, communicator_opt)
    !
    !:Purpose: Calling 'mpi_alltoall' for a integer(8) scalar or array
    !
    implicit none

    ! Arguments:
    integer(8),     contiguous, intent(in)  :: sending(..)      ! integer(8) data sent to all MPI ranks
    integer(8),     contiguous, intent(out) :: receiving(..)    ! integer(8) array which stores the data received
    integer,        optional,   intent(in)  :: length_opt       ! size of the input data
    type(mpi_comm), optional,   intent(in)  :: communicator_opt ! the MPI communicator

    ! Locals:
    integer :: ierr, length
    type(mpi_comm) :: communicator

    communicator = handleCommunicator(communicator_opt)
    length = handleLength(sending, communicator, length_opt)

    call mpi_alltoall(sending,   length, mpi_integer8,  &
                      receiving, length, mpi_integer8, communicator, ierr)

    call handleMpiError(ierr, 'mmpi_alltoall_integer8')

  end subroutine mmpi_alltoall_integer8

  !--------------------------------------------------------------------------
  ! mmpi_alltoall_real4
  !--------------------------------------------------------------------------
  subroutine mmpi_alltoall_real4(sending, receiving, length_opt, communicator_opt)
    !
    !:Purpose: Calling 'mpi_alltoall' for a real(4) scalar or array
    !
    implicit none

    ! Arguments:
    real(4),        contiguous, intent(in)  :: sending(..)      ! real(4) data sent to all MPI ranks
    real(4),        contiguous, intent(out) :: receiving(..)    ! real(4) array which stores the data received
    integer,        optional,   intent(in)  :: length_opt       ! size of the input data
    type(mpi_comm), optional,   intent(in)  :: communicator_opt ! the MPI communicator

    ! Locals:
    integer :: ierr, length
    type(mpi_comm) :: communicator

    communicator = handleCommunicator(communicator_opt)
    length = handleLength(sending, communicator, length_opt)

    call mpi_alltoall(sending,   length, mpi_real4,  &
                      receiving, length, mpi_real4, communicator, ierr)

    call handleMpiError(ierr, 'mmpi_alltoall_real4')

  end subroutine mmpi_alltoall_real4

  !--------------------------------------------------------------------------
  ! mmpi_alltoall_real8
  !--------------------------------------------------------------------------
  subroutine mmpi_alltoall_real8(sending, receiving, length_opt, communicator_opt)
    !
    !:Purpose: Calling 'mpi_alltoall' for a real(8) scalar or array
    !
    implicit none

    ! Arguments:
    real(8),        contiguous, intent(in)  :: sending(..)      ! real(8) data sent to all MPI ranks
    real(8),        contiguous, intent(out) :: receiving(..)    ! real(8) array which stores the data received
    integer,        optional,   intent(in)  :: length_opt       ! size of the input data
    type(mpi_comm), optional,   intent(in)  :: communicator_opt ! the MPI communicator

    ! Locals:
    integer :: ierr, length
    type(mpi_comm) :: communicator

    communicator = handleCommunicator(communicator_opt)
    length = handleLength(sending, communicator, length_opt)

    call mpi_alltoall(sending,   length, mpi_real8,  &
                      receiving, length, mpi_real8, communicator, ierr)

    call handleMpiError(ierr, 'mmpi_alltoall_real8')

  end subroutine mmpi_alltoall_real8

  !--------------------------------------------------------------------------
  ! mmpi_alltoallv_real4
  !--------------------------------------------------------------------------
  subroutine mmpi_alltoallv_real4(sending,   sendsizes, senddispls, &
                                  receiving, recvsizes, recvdispls, &
                                  communicator_opt)
    !
    !:Purpose: Calling 'mpi_alltoallv' for a real(4) scalar or array
    !
    implicit none

    ! Arguments:
    real(4),        contiguous, intent(in)  :: sending(..)      ! real(4) data sent to all MPI ranks
    integer,                    intent(in)  :: sendsizes(:)     ! array containing the size of each array to be sent
    integer,                    intent(in)  :: senddispls(:)    ! displacement offsets in the input array
    real(4),        contiguous, intent(out) :: receiving(..)    ! real(4) array which stores the data received
    integer,                    intent(in)  :: recvsizes(:)     ! array containing the size of each array to be received
    integer,                    intent(in)  :: recvdispls(:)    ! displacement offsets in the output array
    type(mpi_comm), optional,   intent(in)  :: communicator_opt ! the MPI communicator

    ! Locals:
    integer :: ierr
    type(mpi_comm) :: communicator

    communicator = handleCommunicator(communicator_opt)

    call mpi_alltoallv(sending,   sendsizes, senddispls, mpi_real4, &
                       receiving, recvsizes, recvdispls, mpi_real4, &
                       communicator, ierr)

    call handleMpiError(ierr, 'mmpi_alltoallv_real4')

  end subroutine mmpi_alltoallv_real4

  !--------------------------------------------------------------------------
  ! mmpi_alltoallv_real8
  !--------------------------------------------------------------------------
  subroutine mmpi_alltoallv_real8(sending,   sendsizes, senddispls, &
                                  receiving, recvsizes, recvdispls, &
                                  communicator_opt)
    !
    !:Purpose: Calling 'mpi_alltoallv' for a real(8) scalar or array
    !
    implicit none

    real(8),        contiguous, intent(in)  :: sending(..)      ! real(8) data sent to all MPI ranks
    integer,                    intent(in)  :: sendsizes(:)     ! array containing the size of each array to be sent
    integer,                    intent(in)  :: senddispls(:)    ! displacement offsets in the input array
    real(8),        contiguous, intent(out) :: receiving(..)    ! real(8) array which stores the data received
    integer,                    intent(in)  :: recvsizes(:)     ! array containing the size of each array to be received
    integer,                    intent(in)  :: recvdispls(:)    ! displacement offsets in the output array
    type(mpi_comm), optional,   intent(in)  :: communicator_opt ! the MPI communicator

    ! Locals:
    integer :: ierr
    type(mpi_comm) :: communicator

    communicator = handleCommunicator(communicator_opt)

    call mpi_alltoallv(sending,   sendsizes, senddispls, mpi_real8, &
                       receiving, recvsizes, recvdispls, mpi_real8, &
                       communicator, ierr)

    call handleMpiError(ierr, 'mmpi_alltoallv_real8')

  end subroutine mmpi_alltoallv_real8

  !--------------------------------------------------------------------------
  ! mmpi_allReduce_logical
  !--------------------------------------------------------------------------
  subroutine mmpi_allReduce_logical(sending, receiving, operation, length_opt)
    !
    !:Purpose: Calling 'mpi_allReduce' for a logical scalar or array
    !
    implicit none

    ! Arguments:
    logical, contiguous, intent(in)  :: sending(..)   ! logical data sent to all MPI ranks
    logical, contiguous, intent(out) :: receiving(..) ! logical array which stores the result of the reduce operation
    type(mpi_op),        intent(in)  :: operation     ! operation used to reduce the data
    integer, optional,   intent(in)  :: length_opt    ! size of the input data

    ! Locals:
    integer :: ierr, length

    length = handleLength(sending, length_opt)

    call mpi_allReduce(sending, receiving, length, mpi_logical, operation, &
                       mmpi_comm_grid, ierr)

    call handleMpiError(ierr, 'mmpi_allReduce_logical')

  end subroutine mmpi_allReduce_logical

  !--------------------------------------------------------------------------
  ! mmpi_allReduce_integer
  !--------------------------------------------------------------------------
  subroutine mmpi_allReduce_integer(sending, receiving, operation, length_opt)
    !
    !:Purpose: Calling 'mpi_allReduce' for a integer scalar or array
    !
    implicit none

    ! Arguments:
    integer, contiguous, intent(in)  :: sending(..)   ! integer data sent to all MPI ranks
    integer, contiguous, intent(out) :: receiving(..) ! integer array which stores the result of the reduce operation
    type(mpi_op),        intent(in)  :: operation     ! operation used to reduce the data
    integer, optional,   intent(in)  :: length_opt    ! size of the input data

    ! Locals:
    integer :: ierr, length

    length = handleLength(sending, length_opt)

    call mpi_allReduce(sending, receiving, length, mpi_integer, operation, &
                       mmpi_comm_grid, ierr)

    call handleMpiError(ierr, 'mmpi_allReduce_integer')

  end subroutine mmpi_allReduce_integer

  !--------------------------------------------------------------------------
  ! mmpi_allReduce_integer8
  !--------------------------------------------------------------------------
  subroutine mmpi_allReduce_integer8(sending, receiving, operation, length_opt)
    !
    !:Purpose: Calling 'mpi_allReduce' for a integer(8) scalar or array
    !
    implicit none

    ! Arguments:
    integer(8), contiguous, intent(in)  :: sending(..)   ! integer(8) data sent to all MPI ranks
    integer(8), contiguous, intent(out) :: receiving(..) ! integer(8) array which stores the result of the reduce operation
    type(mpi_op),           intent(in)  :: operation     ! operation used to reduce the data
    integer, optional,      intent(in)  :: length_opt    ! size of the input data

    ! Locals:
    integer :: ierr, length

    length = handleLength(sending, length_opt)

    call mpi_allReduce(sending, receiving, length, mpi_integer8, operation, &
                       mmpi_comm_grid, ierr)

    call handleMpiError(ierr, 'mmpi_allReduce_integer8')

  end subroutine mmpi_allReduce_integer8

  !--------------------------------------------------------------------------
  ! mmpi_allReduce_real4
  !--------------------------------------------------------------------------
  subroutine mmpi_allReduce_real4(sending, receiving, operation, length_opt)
    !
    !:Purpose: Calling 'mpi_allReduce' for a real(4) scalar or array
    !
    implicit none

    ! Arguments:
    real(4), contiguous, intent(in)  :: sending(..)   ! real(4) data sent to all MPI ranks
    real(4), contiguous, intent(out) :: receiving(..) ! real(4) array which stores the result of the reduce operation
    type(mpi_op),        intent(in)  :: operation     ! operation used to reduce the data
    integer, optional,   intent(in)  :: length_opt    ! size of the input data

    ! Locals:
    integer :: ierr, length

    length = handleLength(sending, length_opt)

    call mpi_allReduce(sending, receiving, length, mpi_real4, operation, &
                       mmpi_comm_grid, ierr)

    call handleMpiError(ierr, 'mmpi_allReduce_real4')

  end subroutine mmpi_allReduce_real4

  !--------------------------------------------------------------------------
  ! mmpi_allReduce_real8
  !--------------------------------------------------------------------------
  subroutine mmpi_allReduce_real8(sending, receiving, operation, length_opt)
    !
    !:Purpose: Calling 'mpi_allReduce' for a real(8) scalar or array
    !
    implicit none

    ! Arguments:
    real(8), contiguous, intent(in)  :: sending(..)   ! real(8) data sent to all MPI ranks
    real(8), contiguous, intent(out) :: receiving(..) ! real(8) array which stores the result of the reduce operation
    type(mpi_op),        intent(in)  :: operation     ! operation used to reduce the data
    integer, optional,   intent(in)  :: length_opt    ! size of the input data

    ! Locals:
    integer :: ierr, length

    length = handleLength(sending, length_opt)

    call mpi_allReduce(sending, receiving, length, mpi_real8, operation, &
                       mmpi_comm_grid, ierr)

    call handleMpiError(ierr, 'mmpi_allReduce_real8')

  end subroutine mmpi_allReduce_real8

  !--------------------------------------------------------------------------
  ! mmpi_allReduce_scalar_integer
  !--------------------------------------------------------------------------
  subroutine mmpi_allReduce_scalar_integer(localGlobalValue)
    !
    !:Purpose: Perform 'mmpi_allReduce' to sum integer values over all
    !          mpi tasks and copy result back to same variable.
    !
    implicit none

    ! Arguments:
    integer, intent(inout) :: localGlobalValue ! input value used in the sum and sum is put back in that variable

    ! Locals:
    integer :: localValue, globalValue

    localValue = localGlobalValue
    call mmpi_allReduce(localValue, globalValue, mpi_sum)
    localGlobalValue = globalValue

  end subroutine mmpi_allReduce_scalar_integer


  !--------------------------------------------------------------------------
  ! mmpi_gathervDisplacements
  !--------------------------------------------------------------------------
  subroutine mmpi_gathervDisplacements(length, allLengths, displacements)
    !
    !:Purpose: Compute the displacements offsets for a future 'gatherv' call
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: length           ! size of the data on this MPI rank
    integer, intent(out) :: allLengths(:)    ! array containing the size of the data for each MPI rank
    integer, intent(out) :: displacements(:) ! offsets in the array to look for the data for each MPI rank

    ! Locals:
    integer :: procIndex

    call mmpi_allGather(length, allLengths)

    if ( mmpi_myid == 0 ) then
      displacements(1) = 0
      do procIndex = 2, mmpi_nprocs
        displacements(procIndex) = displacements(procIndex-1) + allLengths(procIndex-1)
      end do
    else
      displacements(:) = 0
    end if

  end subroutine mmpi_gathervDisplacements

  !--------------------------------------------------------------------------
  ! mmpi_gatherv_logical
  !--------------------------------------------------------------------------
  subroutine mmpi_gatherv_logical(sending, receiving, length_opt)
    !
    !:Purpose: Calling 'mpi_gatherv' for a logical scalar or array
    !          It computes 'allLengths' and 'displacements' locally.
    !
    implicit none

    ! Arguments:
    logical, contiguous, intent(in)  :: sending(..)   ! logical data sent to all MPI ranks
    logical, contiguous, intent(out) :: receiving(..) ! logical array which stores the data received
    integer, optional,   intent(in)  :: length_opt    ! size of the input array

    ! Locals:
    integer :: length
    integer :: allLengths(mmpi_nprocs)
    integer :: displacements(mmpi_nprocs)

    length = handleLength(sending, length_opt)

    call mmpi_gathervDisplacements(length, allLengths, displacements)

    call mmpi_gatherv(sending, receiving, allLengths, displacements, length)

  end subroutine mmpi_gatherv_logical

  !--------------------------------------------------------------------------
  ! mmpi_gatherv_logical_displs
  !--------------------------------------------------------------------------
  subroutine mmpi_gatherv_logical_displs(sending, receiving, allLengths, displacements, length_opt)
    !
    !:Purpose: Calling 'mpi_gatherv' for a logical scalar or
    !          array when 'allLengths' and 'displacements' are both
    !          provided.
    !
    implicit none

    ! Arguments:
    logical, contiguous, intent(in)  :: sending(..)      ! logical data sent to all MPI ranks
    logical, contiguous, intent(out) :: receiving(..)    ! logical array which stores the data received
    integer,             intent(in)  :: allLengths(:)    ! array containing the size of the data for each MPI rank
    integer,             intent(in)  :: displacements(:) ! offsets in the array to look for the data for each MPI rank
    integer, optional,   intent(in)  :: length_opt       ! size of the input array

    ! Locals:
    integer, parameter :: procID = 0
    integer :: ierr, length

    length = handleLength(sending, length_opt)

    call mpi_gatherv(sending,   length,                    mpi_logical, &
                     receiving, allLengths, displacements, mpi_logical, &
                     procID,    mmpi_comm_GRID, ierr )

    call handleMpiError(ierr, 'mmpi_gatherv_logical_displs')

  end subroutine mmpi_gatherv_logical_displs

  !--------------------------------------------------------------------------
  ! mmpi_gatherv_integer
  !--------------------------------------------------------------------------
  subroutine mmpi_gatherv_integer(sending, receiving, length_opt)
    !
    !:Purpose: Calling 'mpi_gatherv' for a integer scalar or array
    !          It computes 'allLengths' and 'displacements' locally.
    !
    implicit none

    ! Arguments:
    integer, contiguous, intent(in)  :: sending(..)   ! integer data sent to all MPI ranks
    integer, contiguous, intent(out) :: receiving(..) ! integer array which stores the data received
    integer, optional,   intent(in)  :: length_opt    ! size of the input array

    ! Locals:
    integer :: length
    integer :: allLengths(mmpi_nprocs)
    integer :: displacements(mmpi_nprocs)

    length = handleLength(sending, length_opt)

    call mmpi_gathervDisplacements(length, allLengths, displacements)

    call mmpi_gatherv(sending, receiving, allLengths, displacements, length)

  end subroutine mmpi_gatherv_integer

  !--------------------------------------------------------------------------
  ! mmpi_gatherv_integer_displs
  !--------------------------------------------------------------------------
  subroutine mmpi_gatherv_integer_displs(sending, receiving, allLengths, displacements, length_opt)
    !
    !:Purpose: Calling 'mpi_gatherv' for a logical scalar or
    !          array when 'allLengths' and 'displacements' are both
    !          provided.
    !
    implicit none

    ! Arguments:
    integer, contiguous, intent(in)  :: sending(..)      ! integer data sent to all MPI ranks
    integer, contiguous, intent(out) :: receiving(..)    ! integer array which stores the data received
    integer,             intent(in)  :: allLengths(:)    ! array containing the size of the data for each MPI rank
    integer,             intent(in)  :: displacements(:) ! offsets in the array to look for the data for each MPI rank
    integer, optional,   intent(in)  :: length_opt       ! size of the input array

    ! Locals:
    integer, parameter :: procID = 0
    integer :: ierr, length

    length = handleLength(sending, length_opt)

    call mpi_gatherv(sending,   length,                    mpi_integer, &
                     receiving, allLengths, displacements, mpi_integer, &
                     procID,    mmpi_comm_GRID, ierr )

    call handleMpiError(ierr, 'mmpi_gatherv_integer_displs')

  end subroutine mmpi_gatherv_integer_displs

  !--------------------------------------------------------------------------
  ! mmpi_gatherv_real4
  !--------------------------------------------------------------------------
  subroutine mmpi_gatherv_real4(sending, receiving, length_opt)
    !
    !:Purpose: Calling 'mpi_gatherv' for a real(4) scalar or array
    !          It computes 'allLengths' and 'displacements' locally.
    !
    implicit none

    ! Arguments:
    real(4), contiguous, intent(in)  :: sending(..)   ! real(4) data sent to all MPI ranks
    real(4), contiguous, intent(out) :: receiving(..) ! real(4) array which stores the data received
    integer, optional,   intent(in)  :: length_opt    ! size of the input array

    ! Locals:
    integer, parameter :: procID = 0
    integer :: length
    integer :: allLengths(mmpi_nprocs)
    integer :: displacements(mmpi_nprocs)

    length = handleLength(sending, length_opt)

    call mmpi_gathervDisplacements(length, allLengths, displacements)

    call mmpi_gatherv(sending, receiving, allLengths, displacements, length)

  end subroutine mmpi_gatherv_real4

  !--------------------------------------------------------------------------
  ! mmpi_gatherv_real4_displs
  !--------------------------------------------------------------------------
  subroutine mmpi_gatherv_real4_displs(sending, receiving, allLengths, displacements, length_opt)
    !
    !:Purpose: Calling 'mpi_gatherv' for a real(4) scalar or array
    !          when 'allLengths' and 'displacements' are both
    !          provided.
    !
    implicit none

    ! Arguments:
    real(4), contiguous, intent(in)  :: sending(..)      ! real(4) data sent to all MPI ranks
    real(4), contiguous, intent(out) :: receiving(..)    ! real(4) array which stores the data received
    integer,             intent(in)  :: allLengths(:)    ! array containing the size of the data for each MPI rank
    integer,             intent(in)  :: displacements(:) ! offsets in the array to look for the data for each MPI rank
    integer, optional,   intent(in)  :: length_opt       ! size of the input array

    ! Locals:
    integer, parameter :: procID = 0
    integer :: ierr, length

    length = handleLength(sending, length_opt)

    call mpi_gatherv(sending,   length,                    mpi_real4, &
                     receiving, allLengths, displacements, mpi_real4, &
                     procID,    mmpi_comm_GRID, ierr )

    call handleMpiError(ierr, 'mmpi_gatherv_real4_displs')

  end subroutine mmpi_gatherv_real4_displs

  !--------------------------------------------------------------------------
  ! mmpi_gatherv_real8
  !--------------------------------------------------------------------------
  subroutine mmpi_gatherv_real8(sending, receiving, length_opt)
    !
    !:Purpose: Calling 'mpi_gatherv' for a real8 scalar or array
    !          It computes 'allLengths' and 'displacements' locally.
    !
    implicit none

    ! Arguments:
    real(8), contiguous, intent(in)  :: sending(..)   ! real(8) data sent to all MPI ranks
    real(8), contiguous, intent(out) :: receiving(..) ! real(8) array which stores the data received
    integer, optional,   intent(in)  :: length_opt    ! size of the input array

    ! Locals:
    integer, parameter :: procID = 0
    integer :: length
    integer :: allLengths(mmpi_nprocs)
    integer :: displacements(mmpi_nprocs)

    length = handleLength(sending, length_opt)

    call mmpi_gathervDisplacements(length, allLengths, displacements)

    call mmpi_gatherv(sending, receiving, allLengths, displacements, length)

  end subroutine mmpi_gatherv_real8

  !--------------------------------------------------------------------------
  ! mmpi_gatherv_real8_displs
  !--------------------------------------------------------------------------
  subroutine mmpi_gatherv_real8_displs(sending, receiving, allLengths, displacements, length_opt)
    !
    !:Purpose: Calling 'mpi_gatherv' for a real8 scalar or array
    !          when 'allLengths' and 'displacements' are both
    !          provided.
    !
    implicit none

    ! Arguments:
    real(8), contiguous, intent(in)  :: sending(..)      ! real(8) data sent to all MPI ranks
    real(8), contiguous, intent(out) :: receiving(..)    ! real(8) array which stores the data received
    integer,             intent(in)  :: allLengths(:)    ! array containing the size of the data for each MPI rank
    integer,             intent(in)  :: displacements(:) ! offsets in the array to look for the data for each MPI rank
    integer, optional,   intent(in)  :: length_opt       ! size of the input array

    ! Locals:
    integer, parameter :: procID = 0
    integer :: ierr, length

    length = handleLength(sending, length_opt)

    call mpi_gatherv(sending,   length,                    mpi_real8, &
                     receiving, allLengths, displacements, mpi_real8, &
                     procID,    mmpi_comm_GRID, ierr )

    call handleMpiError(ierr, 'mmpi_gatherv_real8_displs')

  end subroutine mmpi_gatherv_real8_displs

  !--------------------------------------------------------------------------
  ! mmpi_reduce_integer
  !--------------------------------------------------------------------------
  subroutine mmpi_reduce_integer(sending, receiving, operation, length_opt, procID_opt)
    !
    !:Purpose: Calling 'mpi_reduce' for a integer scalar or array
    !
    implicit none

    ! Arguments:
    integer, contiguous, intent(in)  :: sending(..)   ! integer data sent to all MPI ranks
    integer, contiguous, intent(out) :: receiving(..) ! integer array which stores the result of the reduce operation
    type(mpi_op),        intent(in)  :: operation     ! operation used to reduce the data
    integer, optional,   intent(in)  :: length_opt    ! size of the input data
    integer, optional,   intent(in)  :: procID_opt    ! MPI rank which received the result of the reduce operation

    ! Locals:
    integer :: ierr, length, procID

    length = handleLength(sending, length_opt)
    procID = handleProcID(procID_opt)

    call mpi_reduce(sending, receiving, length, mpi_integer, operation, &
                    procID, mmpi_comm_GRID, ierr)

    call handleMpiError(ierr, 'mmpi_reduce_integer')

  end subroutine mmpi_reduce_integer

  !--------------------------------------------------------------------------
  ! mmpi_reduce_real8
  !--------------------------------------------------------------------------
  subroutine mmpi_reduce_real8(sending, receiving, operation, length_opt, procID_opt)
    !
    !:Purpose: Calling 'mpi_reduce' for a real(8) scalar or array
    !
    implicit none

    ! Arguments:
    real(8), contiguous, intent(in)  :: sending(..)   ! real(8) data sent to all MPI ranks
    real(8), contiguous, intent(out) :: receiving(..) ! real(8) array which stores the result of the reduce operation
    type(mpi_op),        intent(in)  :: operation     ! operation used to reduce the data
    integer, optional,   intent(in)  :: length_opt    ! size of the input data
    integer, optional,   intent(in)  :: procID_opt    ! MPI rank which received the result of the reduce operation

    ! Locals:
    integer :: ierr, length, procID

    length = handleLength(sending, length_opt)
    procID = handleProcID(procID_opt)

    call mpi_reduce(sending, receiving, length, mpi_real8, operation, &
                    procID, mmpi_comm_GRID, ierr)

    call handleMpiError(ierr, 'mmpi_reduce_real8')

  end subroutine mmpi_reduce_real8

  !--------------------------------------------------------------------------
  ! mmpi_scatterv_real4
  !--------------------------------------------------------------------------
  subroutine mmpi_scatterv_real4(sending, receiving, allLengths, displacements, length_opt)
    !
    !:Purpose: Calling 'mpi_scatterv' for a real4 scalar or array
    !
    implicit none

    ! Arguments:
    real(4), contiguous, intent(in)  :: sending(..)      ! real(4) data sent to all MPI ranks
    real(4), contiguous, intent(out) :: receiving(..)    ! real(4) array which stores the data received
    integer,             intent(in)  :: allLengths(:)    ! array containing the size of the data for each MPI rank
    integer,             intent(in)  :: displacements(:) ! offsets in the array to look for the data for each MPI rank
    integer, optional,   intent(in)  :: length_opt       ! size of the input array

    ! Locals:
    integer :: ierr, length
    integer, parameter :: procID = 0

    length = handleLength(receiving, length_opt)

    call mpi_scatterv(sending,   allLengths, displacements, mpi_real4, &
                      receiving, length,                    mpi_real4, &
                      procID, mmpi_comm_GRID, ierr)

    call handleMpiError(ierr, 'mmpi_scatterv_real4')

  end subroutine mmpi_scatterv_real4

  !--------------------------------------------------------------------------
  ! mmpi_scatterv_real8
  !--------------------------------------------------------------------------
  subroutine mmpi_scatterv_real8(sending, receiving, allLengths, displacements, length_opt)
    !
    !:Purpose: Calling 'mpi_scatterv' for a real8 scalar or array
    !
    implicit none

    ! Arguments:
    real(8), contiguous, intent(in)  :: sending(..)      ! real(8) data sent to all MPI ranks
    real(8), contiguous, intent(out) :: receiving(..)    ! real(8) array which stores the data received
    integer,             intent(in)  :: allLengths(:)    ! array containing the size of the data for each MPI rank
    integer,             intent(in)  :: displacements(:) ! offsets in the array to look for the data for each MPI rank
    integer, optional,   intent(in)  :: length_opt       ! size of the input array

    ! Locals:
    integer :: ierr, length
    integer, parameter :: procID = 0

    length = handleLength(receiving, length_opt)

    call mpi_scatterv(sending,   allLengths, displacements, mpi_real8, &
                      receiving, length,                    mpi_real8, &
                      procID, mmpi_comm_GRID, ierr)

    call handleMpiError(ierr, 'mmpi_scatterv_real8')

  end subroutine mmpi_scatterv_real8

  !--------------------------------------------------------------------------
  ! mmpi_send_real8
  !--------------------------------------------------------------------------
  subroutine mmpi_send_real8(data, tag, procID, length_opt, communicator_opt)
    !
    !:Purpose: Calling 'mpi_send' for a real(8) scalar or array
    !
    implicit none

    ! Arguments:
    real(8),      contiguous, intent(in) :: data(..)         ! real(8) data sent to all MPI ranks
    integer,                  intent(in) :: tag              ! MPI rank which sends the data
    integer,        optional, intent(in) :: procID           ! MPI rank which sends the data
    integer,        optional, intent(in) :: length_opt       ! size of the input array
    type(mpi_comm), optional, intent(in) :: communicator_opt ! the MPI communicator

    ! Locals:
    integer        :: ierr, length
    type(mpi_comm) :: communicator

    length = handleLength(data, length_opt)
    communicator = handleCommunicator(communicator_opt)

    call mpi_send(data, length, mpi_real8, procID, tag, communicator, ierr)

    call handleMpiError(ierr, 'mmpi_send_real8')

  end subroutine mmpi_send_real8

  !--------------------------------------------------------------------------
  ! mmpi_recv_real8
  !--------------------------------------------------------------------------
  subroutine mmpi_recv_real8(data, tag, procID, length_opt, communicator_opt)
    !
    !:Purpose: Calling 'mpi_recv' for a real(8) scalar or array
    !
    implicit none

    ! Arguments:
    real(8),     contiguous,  intent(out) :: data(..)         ! real(8) array which stores the data received
    integer,                  intent(in)  :: tag              ! tag which identified the data received
    integer,                  intent(in)  :: procID           ! MPI rank which receives the data
    integer,        optional, intent(in)  :: length_opt       ! size of the input array
    type(mpi_comm), optional, intent(in)  :: communicator_opt ! the MPI communicator

    ! Locals:
    integer :: ierr, length
    type(mpi_comm) :: communicator

    length = handleLength(data, length_opt)
    communicator = handleCommunicator(communicator_opt)

    call mpi_recv(data, length, mpi_real8, procID, tag, communicator, &
                  MPI_STATUS_IGNORE, ierr)

    call handleMpiError(ierr, 'mmpi_recv_real8')

  end subroutine mmpi_recv_real8

  !--------------------------------------------------------------------------
  ! mmpi_sendrecv_real8
  !--------------------------------------------------------------------------
  subroutine mmpi_sendrecv_real8(sending,   sendProcID, sendTag, &
                                 receiving, recvProcID, recvTag, &
                                 length_opt, communicator_opt)
    !
    !:Purpose: Calling 'mpi_sendrecv' for a real(8) scalar or array
    !
    implicit none

    ! Arguments:
    real(8),      contiguous, intent(in)  :: sending(..)      ! real(8) data sent to all MPI ranks
    integer,                  intent(in)  :: sendProcID       ! MPI rank which sends the data
    integer,                  intent(in)  :: sendTag          ! tag which identified the data sent
    real(8),      contiguous, intent(out) :: receiving(..)    ! real(8) array which stores the data received
    integer,                  intent(in)  :: recvProcID       ! MPI rank which receives the data
    integer,                  intent(in)  :: recvTag          ! tag which identified the data received
    integer,        optional, intent(in)  :: length_opt       ! size of the input array
    type(mpi_comm), optional, intent(in)  :: communicator_opt ! the MPI communicator

    ! Locals:
    integer :: ierr, length
    type(mpi_comm) :: communicator

    length = handleLength(sending, length_opt)
    communicator = handleCommunicator(communicator_opt)

    call mpi_sendrecv(sending,   length, mpi_real8, sendProcID, sendTag, &
                      receiving, length, mpi_real8, recvProcID, recvTag, &
                      communicator, MPI_STATUS_IGNORE, ierr)

    call handleMpiError(ierr, 'mmpi_sendrecv_real8')

  end subroutine mmpi_sendrecv_real8

  !--------------------------------------------------------------------------
  ! handleCommunicator (private)
  !--------------------------------------------------------------------------
  function handleCommunicator(communicator_opt) result(communicator)
    !
    !:Purpose: Process 'communicator_opt' optional argument
    !
    implicit none

    ! Arguments:
    type(mpi_comm), optional, intent(in)  :: communicator_opt ! MPI communicator identifier
    ! Result:
    type(mpi_comm) :: communicator ! MPI communicator identifier

    if ( present(communicator_opt) ) then
      communicator = communicator_opt
    else
      communicator = mmpi_comm_GRID
    end if

  end function handleCommunicator

  !--------------------------------------------------------------------------
  ! handleProcID (private)
  !--------------------------------------------------------------------------
  function handleProcID(procID_opt) result(procID)
    !
    !:Purpose: Process 'procID_opt' optional argument
    !
    implicit none

    ! Arguments:
    integer, optional, intent(in)  :: procID_opt ! optional MPI rank
    ! Result:
    integer :: procID ! default procID if 'procID_opt' is not provided

    if (present(procID_opt)) then
      procID = procID_opt
    else
      procID = 0
    end if

  end function handleProcID

  !--------------------------------------------------------------------------
  ! handleLengthOnly (private)
  !--------------------------------------------------------------------------
  function handleLengthOnly(inputData, length_opt) result(length)
    !
    !:Purpose: Process 'length_opt' optional argument
    !
    implicit none

    ! Arguments:
    type(*),           intent(in)  :: inputData(..) ! input array (rank and size will be used to find the length)
    integer, optional, intent(in)  :: length_opt    ! optional length
    ! Result:
    integer :: length ! length if 'length_opt' is not provided

    if ( present(length_opt) ) then
      length = length_opt
    else
      if ( rank(inputData) == 0 ) then
        length = 1
      else
        length = size(inputData)
      end if ! 'else' of 'if ( rank(inputData) == 0 ) then'
    end if ! 'else' of 'if ( present(length_opt) ) then'

  end function handleLengthOnly

  !--------------------------------------------------------------------------
  ! handleLengthWithRespectToCommunicator (private)
  !--------------------------------------------------------------------------
  function handleLengthWithRespectToCommunicator(inputData, communicator, length_opt) result(length)
    !
    !:Purpose: Process 'length_opt' optional argument when a
    !          communicator has to be taken into account.
    !
    implicit none

    ! Arguments:
    type(*),           intent(in) :: inputData(..) ! input array (rank and size will be used to find the length)
    type(mpi_comm),    intent(in) :: communicator  ! string identifying the MPI communicator
    integer, optional, intent(in) :: length_opt    ! optional length
    ! Result:
    integer :: length ! length if 'length_opt' is not provided
    integer :: ierr, nameLength
    character(len=MPI_MAX_OBJECT_NAME) :: commName

    length = handleLength(inputData, length_opt)

    if ( .not. present(length_opt) ) then
      ! If 'communicator_opt' is provided, then we will find the
      ! size to the data to exchange by taking into account that
      ! the input data is allocated with an added dimension to
      ! store the data for each MPI rank.

      ! The input array is of the form 'inputData(..,mmpi_{nprocs,npex,npey})'
      ! and the size expected by 'mpi_alltoall' does not include
      ! the last dimension.
      if ( communicator == mmpi_comm_GRID ) then
        length = length/mmpi_nprocs
      else if ( communicator == mmpi_comm_NS ) then
        length = length/mmpi_npey
      else if ( communicator == mmpi_comm_EW ) then
        length = length/mmpi_npex
      else
        call mpi_comm_get_name(communicator, commName, nameLength, ierr)
        call utl_abort('handleLengthWithRespectToCommunicator: cannot guess the size of the data for communicator='''// commName // '''')
      end if
    end if ! 'if ( .not. present(length_opt) ) then'

  end function handleLengthWithRespectToCommunicator

  !--------------------------------------------------------------------------
  ! handleMpiError (private)
  !--------------------------------------------------------------------------
  subroutine handleMpiError(errCode, context)
    !
    !:Purpose: Handle gracefully the MPI error code
    !
    implicit none

    ! Arguments:
    integer,          intent(in) :: errCode ! error code from MPI routine
    character(len=*), intent(in) :: context ! string containing the context if an error is raised

    ! Locals:
    character(len=MPI_MAX_ERROR_STRING) :: errorMsg
    integer :: resultlen, ierror

    if (errCode /= MPI_SUCCESS) then
      call MPI_Error_string(errcode, errorMsg, resultlen, ierror)
      call utl_abort('MPI error found in ' // context // ' : ' // trim(errorMsg))
    end if
  end subroutine handleMpiError

end module midasMpi_mod

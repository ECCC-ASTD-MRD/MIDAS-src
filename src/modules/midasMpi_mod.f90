
module midasMpi_mod
  !
  ! MODULE midasMpi_mod (prefix='mmpi' category='8. Low-level utilities and constants')
  !
  !:Purpose:  Subroutines/functions and public variables related to general aspects of mpi.
  !           Also, subroutine and public variables related to the mpi decomposition
  !           specific to the MIDAS code.
  !
  use mpi
  use utilities_mod

  implicit none
  save
  private

  ! Public variables
  logical, public, protected :: mmpi_doBarrier = .true.
  integer, public, protected :: mmpi_myid   = 0
  integer, public, protected :: mmpi_myidHost = 0
  integer, public, protected :: mmpi_nprocs = 0
  integer, public, protected :: mmpi_myidx  = 0
  integer, public, protected :: mmpi_myidy  = 0
  integer, public, protected :: mmpi_npex   = 0
  integer, public, protected :: mmpi_npey   = 0
  integer, public, protected :: mmpi_numthread = 0
  integer, public, protected :: mmpi_comm_EW, mmpi_comm_NS, mmpi_comm_GRID, mmpi_mpicomm_SHARED
  integer, public, protected :: mmpi_datyp_real4, mmpi_datyp_real8, mmpi_datyp_int
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
  public :: mmpi_bcast, mmpi_gather, mmpi_allGather

  ! module interfaces
  ! -----------------

  ! general interface for rpn_comm_bcast and rpn_comm_bcastc
  interface mmpi_bcast
    module procedure mmpi_bcast_character
    module procedure mmpi_bcast_logical
    module procedure mmpi_bcast_integer
    module procedure mmpi_bcast_real4
    module procedure mmpi_bcast_real8
  end interface mmpi_bcast

  ! general interface for rpn_comm_gather
  interface mmpi_gather
    module procedure mmpi_gather_logical
    module procedure mmpi_gather_integer
    module procedure mmpi_gather_integer8
    module procedure mmpi_gather_real4
    module procedure mmpi_gather_real8
  end interface mmpi_gather

  ! general interface for rpn_comm_allGather
  interface mmpi_allGather
    module procedure mmpi_allGather_logical
    !module procedure mmpi_allGather_integer
    !module procedure mmpi_allGather_integer8
    module procedure mmpi_allGather_real4
    module procedure mmpi_allGather_real8
  end interface mmpi_allGather

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
    integer :: mythread,numthread,omp_get_thread_num,omp_get_num_threads,rpn_comm_mype
    integer :: ierr, numNodeMasters
    integer :: rpn_comm_comm, rpn_comm_datyp
    integer(kind=MPI_ADDRESS_KIND) :: maxTagValue
    integer, allocatable :: allMyidHost(:)
    logical :: flag

    ! Namelist variables
    integer :: npex  ! number of MPI tasks in 'x' direction (set automatically by launch script)
    integer :: npey  ! number of MPI tasks in 'y' direction (set automatically by launch script)

    ! Initilize MPI
    npex=0
    npey=0
    call rpn_comm_init(mmpi_getptopo,mmpi_myid,mmpi_nprocs,npex,npey)

    ! this is a special mpi communicator (not rpn_comm) for using shared memory arrays
    call mpi_comm_split_type(mpi_comm_world, mpi_comm_type_shared, 0,  &
                             mpi_info_null, mmpi_mpicomm_SHARED,ierr)

    if(mmpi_nprocs.lt.1) then
      mmpi_nprocs=1
      mmpi_npex=1
      mmpi_npey=1
      mmpi_myid=0
      mmpi_myidHost=0
      mmpi_myidx=0
      mmpi_myidy=0
    else
      ierr = rpn_comm_mype(mmpi_myid,mmpi_myidx,mmpi_myidy)
      mmpi_npex=npex
      mmpi_npey=npey
      call mpi_comm_rank(mmpi_mpicomm_shared, mmpi_myidHost, ierr)
    endif

    write(*,*) 'mmpi_initialize: mmpi_myid, mmpi_myidx, mmpi_myidy, mmpi_myidHost = ', &
                                 mmpi_myid, mmpi_myidx, mmpi_myidy, mmpi_myidHost

    ! Determine list of node masters (i.e. first task on each node)
    allocate(allMyidHost(mmpi_nprocs))
    call rpn_comm_allgather(mmpi_myidHost,    1, 'mpi_integer',  &
                            allMyidHost, 1, 'mpi_integer', &
                            'GRID', ierr)
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

    ! create standard mpi handles to rpn_comm mpi communicators to facilitate 
    ! use of standard mpi routines
    mmpi_comm_EW = rpn_comm_comm('EW')
    mmpi_comm_NS = rpn_comm_comm('NS')
    mmpi_comm_GRID = rpn_comm_comm('GRID')

    mmpi_datyp_real4 = rpn_comm_datyp('MPI_REAL4')
    mmpi_datyp_real8 = rpn_comm_datyp('MPI_REAL8')
    mmpi_datyp_int = rpn_comm_datyp('MPI_INTEGER')

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
  ! mmpi_getptopo
  !--------------------------------------------------------------------------
  subroutine mmpi_getptopo( npex, npey )
    !
    !:Purpose: Subroutine called by the rpn_comm MPI initializing
    !          subroutine rpn_comm_init.
    !
    implicit none

    ! Arguments:
    integer, intent(out) :: npex
    integer, intent(out) :: npey

    ! Locals:
    integer :: ierr
    integer :: nulnam,fnom,fclos
    namelist /ptopo/npex,npey

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

  end subroutine mmpi_getptopo 

  !--------------------------------------------------------------------------
  ! mmpi_allreduce_sumreal8scalar
  !--------------------------------------------------------------------------
  subroutine mmpi_allreduce_sumreal8scalar( sendRecvValue, comm )
    !
    !:Purpose: Version of mpi_allReduce that always performs sum in
    !          the same order.
    !
    implicit none

    ! Arguments:
    real(8),          intent(inout) :: sendRecvValue ! value to be summed over all mpi tasks
    character(len=*), intent(in)    :: comm          ! rpn_comm communicator

    ! Locals:
    integer :: nsize, ierr, root, rank
    real(8), allocatable :: allvalues(:)

    ! do a barrier so that timing on reduce operation is accurate
    call utl_tmg_start(171,'low-level--mpi_allreduce_barr')
    if(mmpi_doBarrier) call rpn_comm_barrier(comm,ierr)
    call utl_tmg_stop(171)

    call utl_tmg_start(170,'low-level--mpi_allreduce_sum8')

    ! determine number of processors in the communicating group
    call rpn_comm_size(comm,nsize,ierr)

    ! determine where to gather the values: first task in group
    call rpn_comm_rank(comm,rank,ierr)
    call rpn_comm_allreduce(rank,root,1,"MPI_INTEGER","MPI_MIN",comm,ierr)

    ! gather values to be added onto 1 processor
    allocate(allvalues(nsize))
    call rpn_comm_gather(sendRecvValue, 1, "MPI_DOUBLE_PRECISION", allvalues, 1, "MPI_DOUBLE_PRECISION", root, comm, ierr)

    ! sum the values on the "root" mpi task and broadcast to group
    if(rank.eq.root) sendRecvValue = sum(allvalues(:))
    deallocate(allvalues)
    call rpn_comm_bcast(sendRecvValue, 1, "MPI_DOUBLE_PRECISION", root, comm, ierr)

    call utl_tmg_stop(170)

  end subroutine mmpi_allreduce_sumreal8scalar

  !--------------------------------------------------------------------------
  ! mmpi_allReduce_sumR8_1d
  !--------------------------------------------------------------------------
  subroutine mmpi_allreduce_sumR8_1d( sendRecvVector, comm )
    !
    ! :Purpose: Perform sum of 1d array over all MPI tasks, guaranteed to
    !           always be in the same order.
    !
    implicit none

    ! Arguments:
    real(8)         , intent(inout)  :: sendRecvVector(:) ! 1-D vector to be summed over all mpi tasks
    character(len=*), intent(in)     :: comm              ! rpn_comm communicator

    ! Locals:
    integer :: nprocs_mpi, numElements, ierr, root, rank
    real(8), allocatable :: all_sendRecvVector(:,:)

    ! do a barrier so that timing on reduce operation is accurate
    call utl_tmg_start(171,'low-level--mpi_allreduce_barr')
    if ( mmpi_doBarrier ) call rpn_comm_barrier(comm,ierr)
    call utl_tmg_stop(171)

    call utl_tmg_start(170,'low-level--mpi_allreduce_sum8')

    numElements = size(sendRecvVector)

    ! determine number of processors in the communicating group
    call rpn_comm_size(comm,nprocs_mpi,ierr)

    ! determine where to gather the values: first task in group
    call rpn_comm_rank(comm,rank,ierr)
    call rpn_comm_allreduce(rank,root,1,"mpi_integer","mpi_min",comm,ierr)

    ! gather vectors to be added onto 1 processor
    allocate(all_sendRecvVector(numElements,0:nprocs_mpi-1))
    call rpn_comm_gather(sendRecvVector    , numElements, "mpi_double_precision", &
                         all_sendRecvVector, numElements, "mpi_double_precision", &
                         root, comm, ierr)

    ! sum the values on the "root" mpi task and broadcast to group
    if ( rank == root ) sendRecvVector(:) = sum(all_sendRecvVector(:,:),2)
    deallocate(all_sendRecvVector)
    call rpn_comm_bcast(sendRecvVector, numElements, "mpi_double_precision", &
                        root, comm, ierr)

    call utl_tmg_stop(170)

  end subroutine mmpi_allreduce_sumR8_1d

  !--------------------------------------------------------------------------
  ! mmpi_allreduce_sumR8_2d
  !--------------------------------------------------------------------------
  subroutine mmpi_allreduce_sumR8_2d( sendRecvVector, comm )
    !
    ! :Purpose: Perform sum of 2d array over all MPI tasks guaranteed
    !           to always be in the same order.
    !
    implicit none

    ! Arguments:
    real(8)         , intent(inout)  :: sendRecvVector(:,:) ! 2-D vector to be summed over all mpi tasks
    character(len=*), intent(in)     :: comm                ! rpn_comm communicator

    ! Locals:
    integer :: nprocs_mpi, numElements1, numElements2, ierr, root, rank
    real(8), allocatable :: all_sendRecvVector(:,:,:)

    ! do a barrier so that timing on reduce operation is accurate
    call utl_tmg_start(171,'low-level--mpi_allreduce_barr')
    if ( mmpi_doBarrier ) call rpn_comm_barrier(comm,ierr)
    call utl_tmg_stop(171)

    call utl_tmg_start(170,'low-level--mpi_allreduce_sum8')

    numElements1 = size(sendRecvVector,1)
    numElements2 = size(sendRecvVector,2)

    ! determine number of processors in the communicating group
    call rpn_comm_size(comm,nprocs_mpi,ierr)

    ! determine where to gather the values: first task in group
    call rpn_comm_rank(comm,rank,ierr)
    call rpn_comm_allreduce(rank,root,1,"mpi_integer","mpi_min",comm,ierr)

    ! gather vectors to be added onto 1 processor
    allocate(all_sendRecvVector(numElements1,numElements2,0:nprocs_mpi-1))
    call rpn_comm_gather(sendRecvVector    , numElements1*numElements2, "mpi_double_precision", &
                         all_sendRecvVector, numElements1*numElements2, "mpi_double_precision", &
                         root, comm, ierr)

    ! sum the values on the "root" mpi task and broadcast to group
    if ( rank == root ) sendRecvVector(:,:) = sum(all_sendRecvVector(:,:,:),3)
    deallocate(all_sendRecvVector)
    call rpn_comm_bcast(sendRecvVector, numElements1*numElements2, "mpi_double_precision", &
                        root, comm, ierr)

    call utl_tmg_stop(170)

  end subroutine mmpi_allreduce_sumR8_2d
   
  !--------------------------------------------------------------------------
  ! mmpi_reduce_sumR8_1d
  !--------------------------------------------------------------------------
  subroutine mmpi_reduce_sumR8_1d( sendVector, recvVector, root, comm )
    !
    ! :Purpose: Perform sum of 1d array over all MPI tasks guaranteed to
    !           always be in the same order.
    !
    implicit none

    ! Arguments:
    real(8)         , intent(in)  :: sendVector(:) ! 1-D vector to be summed over all mpi tasks
    real(8)         , intent(out) :: recvVector(:) ! 1-D vector to be summed over all mpi tasks
    integer         , intent(in)  :: root          ! mpi task id where data is put
    character(len=*), intent(in)  :: comm          ! rpn_comm communicator

    ! Locals:
    integer :: nprocs_mpi, numElements, ierr, rank
    real(8), allocatable :: all_sendRecvVector(:,:)

    ! do a barrier so that timing on reduce operation is accurate
    call utl_tmg_start(171,'low-level--mpi_allreduce_barr')
    if ( mmpi_doBarrier ) call rpn_comm_barrier(comm,ierr)
    call utl_tmg_stop(171)

    call utl_tmg_start(170,'low-level--mpi_allreduce_sum8')

    numElements = size(sendVector)

    ! determine number of processors in the communicating group
    call rpn_comm_size(comm,nprocs_mpi,ierr)

    ! determine rank of group
    call rpn_comm_rank(comm,rank,ierr)

    ! gather vectors to be added onto 1 processor
    if ( rank == root ) then
      allocate(all_sendRecvVector(numElements,0:nprocs_mpi-1))
    else
      allocate(all_sendRecvVector(1,1))
    end if
    call rpn_comm_gather(sendVector        , numElements, "mpi_double_precision", &
                         all_sendRecvVector, numElements, "mpi_double_precision", &
                         root, comm, ierr)

    ! sum the values on the "root" mpi task
    if ( rank == root ) recvVector(:) = sum(all_sendRecvVector(:,:),2)
    deallocate(all_sendRecvVector)

    call utl_tmg_stop(170)

  end subroutine mmpi_reduce_sumR8_1d
   
  !--------------------------------------------------------------------------
  ! mmpi_reduce_sumR8_2d
  !--------------------------------------------------------------------------
  subroutine mmpi_reduce_sumR8_2d( sendVector, recvVector, root, comm )
    !
    ! :Purpose: Perform sum of 2d array over all MPI tasks guaranteed to
    !           always be in the same order.
    !
    implicit none

    ! Arguments:
    real(8)         , intent(in)  :: sendVector(:,:) ! 2-D vector to be summed over all mpi tasks
    real(8)         , intent(out) :: recvVector(:,:) ! 2-D vector to be summed over all mpi tasks
    integer         , intent(in)  :: root            ! mpi task id where data will be put
    character(len=*), intent(in)  :: comm            ! rpn_comm communicator

    ! Locals:
    integer :: nprocs_mpi, numElements1, numElements2, ierr, rank
    real(8), allocatable :: all_sendRecvVector(:,:,:)

    ! do a barrier so that timing on reduce operation is accurate
    call utl_tmg_start(171,'low-level--mpi_allreduce_barr')
    if ( mmpi_doBarrier ) call rpn_comm_barrier(comm,ierr)
    call utl_tmg_stop(171)

    call utl_tmg_start(170,'low-level--mpi_allreduce_sum8')

    numElements1 = size(sendVector,1)
    numElements2 = size(sendVector,2)

    ! determine number of processors in the communicating group
    call rpn_comm_size(comm,nprocs_mpi,ierr)

    ! determine rank of group
    call rpn_comm_rank(comm,rank,ierr)

    ! gather vectors to be added onto 1 processor
    if ( rank == root ) then
      allocate(all_sendRecvVector(numElements1,numElements2,0:nprocs_mpi-1))
    else
      allocate(all_sendRecvVector(1,1,1))
    end if
    call rpn_comm_gather(sendVector        , numElements1*numElements2, "mpi_double_precision", &
                         all_sendRecvVector, numElements1*numElements2, "mpi_double_precision", &
                         root, comm, ierr)

    ! sum the values on the "root" mpi task
    if ( rank == root ) recvVector(:,:) = sum(all_sendRecvVector(:,:,:),3)
    deallocate(all_sendRecvVector)

    call utl_tmg_stop(170)

  end subroutine mmpi_reduce_sumR8_2d

  !--------------------------------------------------------------------------
  ! mmpi_reduce_sumR8_3d
  !--------------------------------------------------------------------------
  subroutine mmpi_reduce_sumR8_3d( sendVector, recvVector, root, comm )
    !
    ! :Purpose: Perform sum of 3d array over all MPI tasks guaranteed to
    !           always be in the same order.
    !
    implicit none

    ! Arguments:
    real(8)         , intent(in)  :: sendVector(:,:,:) ! 3-D vector to be summed over all mpi tasks
    real(8)         , intent(out) :: recvVector(:,:,:) ! 3-D vector to be summed over all mpi tasks
    integer         , intent(in)  :: root              ! mpi task id where data is put
    character(len=*), intent(in)  :: comm              ! rpn_comm communicator

    ! Locals:
    integer :: nprocs_mpi, numElements1, numElements2, numElements3, ierr, rank
    real(8), allocatable :: all_sendRecvVector(:,:,:,:)

    ! do a barrier so that timing on reduce operation is accurate
    call utl_tmg_start(171,'low-level--mpi_allreduce_barr')
    if ( mmpi_doBarrier ) call rpn_comm_barrier(comm,ierr)
    call utl_tmg_stop(171)

    call utl_tmg_start(170,'low-level--mpi_allreduce_sum8')

    numElements1 = size(sendVector,1)
    numElements2 = size(sendVector,2)
    numElements3 = size(sendVector,3)

    ! determine number of processors in the communicating group
    call rpn_comm_size(comm,nprocs_mpi,ierr)

    ! determine rank of group
    call rpn_comm_rank(comm,rank,ierr)

    ! gather vectors to be added onto 1 processor
    if ( rank == root ) then
      allocate(all_sendRecvVector(numElements1,numElements2,numElements3,0:nprocs_mpi-1))
    else
      allocate(all_sendRecvVector(1,1,1,1))
    end if
    call rpn_comm_gather(sendVector        , numElements1*numElements2*numElements3, "mpi_double_precision", &
                         all_sendRecvVector, numElements1*numElements2*numElements3, "mpi_double_precision", &
                         root, comm, ierr)

    ! sum the values on the "root" mpi task
    if ( rank == root ) recvVector(:,:,:) = sum(all_sendRecvVector(:,:,:,:),4)
    deallocate(all_sendRecvVector)

    call utl_tmg_stop(170)

  end subroutine mmpi_reduce_sumR8_3d
  
  !--------------------------------------------------------------------------
  ! mmpi_allgather_string
  !--------------------------------------------------------------------------
  subroutine mmpi_allgather_string( str_list, str_list_all, nlist, nchar, nproc, comm, ierr )
    ! 
    ! :Purpose: Performs the MPI 'allgather' routine for an array of strings
    !
    implicit none

    ! Arguments:
    integer             , intent(in) :: nlist
    integer             , intent(in) :: nchar
    character(len=nchar), intent(in) :: str_list(nlist)
    character(len=*)    , intent(in) :: comm
    integer             , intent(in) :: nproc
    character(len=nchar), intent(out) :: str_list_all(nlist,nproc)
    integer             , intent(out) :: ierr

    ! Locals:
    integer :: num_list(nlist*nchar),num_list_all(nlist*nchar,nproc)
    integer :: ilist,ichar,iproc
              
    ! Convert strings to integer sequences

    do ilist=1,nlist
       do ichar=1,nchar
          num_list((ilist-1)*nchar+ichar) = iachar(str_list(ilist)(ichar:ichar))
       end do
    end do

    ! Perform allgather with converted integer sequences

    call rpn_comm_allgather(num_list,nlist*nchar,"MPI_INTEGER",num_list_all,nlist*nchar,"MPI_INTEGER",comm,ierr)
       
    ! Convert integer sequences to stnid character strings
          
    do iproc=1,nproc
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
    call rpn_comm_allreduce(latPerPE,latPerPEmax,1,'MPI_INTEGER','MPI_MAX','NS',ierr)

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
    !:Purpose: Compute parameters that define the mpi distribution of
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
    call rpn_comm_allreduce(lonPerPE,lonPerPEmax,1,'MPI_INTEGER','MPI_MAX','EW',ierr)

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
    !:Purpose: Calling 'rpn_comm_bcastc' for character array
    !
    implicit none

    ! Arguments:
    character(len=*),  intent(inout) :: charData
    integer, optional, intent(in)    :: procID_opt

    ! Locals:
    integer :: ierr
    integer :: procID

    if (present(procID_opt)) then
      procID = procID_opt
    else
      procID = 0
    end if

    call rpn_comm_bcastc(charData, len(charData), 'MPI_CHARACTER', procID, 'GRID', ierr)
    call handleMpiError(ierr, 'mmpi_bcast_character')

  end subroutine mmpi_bcast_character

  !--------------------------------------------------------------------------
  ! mmpi_bcast_logical
  !--------------------------------------------------------------------------
  subroutine mmpi_bcast_logical(logicalData, length_opt, procID_opt)
    !
    !:Purpose: Calling 'rpn_comm_bcast' for a logical value or array
    !
    implicit none

    ! Arguments:
    logical,           intent(inout) :: logicalData(..)
    integer, optional, intent(in)    :: length_opt
    integer, optional, intent(in)    :: procID_opt

    ! Locals:
    integer :: ierr
    integer :: procID, length

    call handleLengthProcID(length_opt, procID_opt, length, procID, &
                            rank(logicalData), size(logicalData))

    call rpn_comm_bcast(logicalData, length, 'MPI_LOGICAL', procID, 'GRID', ierr)
    call handleMpiError(ierr, 'mmpi_bcast_logical')

  end subroutine mmpi_bcast_logical

  !--------------------------------------------------------------------------
  ! mmpi_bcast_integer
  !--------------------------------------------------------------------------
  subroutine mmpi_bcast_integer(integerData, length_opt, procID_opt)
    !
    !:Purpose: Calling 'rpn_comm_bcast' for an integer scalar or array
    !
    implicit none

    ! Arguments:
    integer, contiguous, intent(inout) :: integerData(..)
    integer, optional,   intent(in)    :: length_opt
    integer, optional,   intent(in)    :: procID_opt

    ! Locals:
    integer :: ierr
    integer :: procID, length

    call handleLengthProcID(length_opt, procID_opt, length, procID, &
                            rank(integerData), size(integerData))

    call rpn_comm_bcast(integerData, length, 'MPI_INTEGER', procID, 'GRID', ierr)
    call handleMpiError(ierr, 'mmpi_bcast_integer')

  end subroutine mmpi_bcast_integer

  !--------------------------------------------------------------------------
  ! mmpi_bcast_real4
  !--------------------------------------------------------------------------
  subroutine mmpi_bcast_real4(real4Data, length_opt, procID_opt)
    !
    !:Purpose: Calling 'rpn_comm_bcast' for a real(4) scalar or array
    !
    implicit none

    ! Arguments:
    real(4), contiguous, intent(inout) :: real4Data(..)
    integer, optional,   intent(in)    :: length_opt
    integer, optional,   intent(in)    :: procID_opt

    ! Locals:
    integer :: ierr
    integer :: procID, length

    call handleLengthProcID(length_opt, procID_opt, length, procID, rank(real4Data), size(real4Data))

    call rpn_comm_bcast(real4Data, length, 'MPI_REAL4', procID, 'GRID', ierr)
    call handleMpiError(ierr, 'mmpi_bcast_real4')

  end subroutine mmpi_bcast_real4

  !--------------------------------------------------------------------------
  ! mmpi_bcast_real8
  !--------------------------------------------------------------------------
  subroutine mmpi_bcast_real8(real8Data, length_opt, procID_opt)
    !
    !:Purpose: Calling 'rpn_comm_bcast' for a real(8) array
    !
    implicit none

    ! Arguments:
    real(8), contiguous, intent(inout) :: real8Data(..)
    integer, optional,   intent(in)    :: length_opt
    integer, optional,   intent(in)    :: procID_opt

    ! Locals:
    integer :: ierr
    integer :: procID, length

    call handleLengthProcID(length_opt, procID_opt, length, procID, rank(real8Data), size(real8Data))

    call rpn_comm_bcast(real8Data, length, 'MPI_REAL8', procID, 'GRID', ierr)
    call handleMpiError(ierr, 'mmpi_bcast_real8')

  end subroutine mmpi_bcast_real8

  !--------------------------------------------------------------------------
  ! mmpi_gather_logical
  !--------------------------------------------------------------------------
  subroutine mmpi_gather_logical(sending, receiving, length_opt, procID_opt)
    !
    !:Purpose: Calling 'rpn_comm_gather' for a logical scalar or array
    !
    implicit none

    ! Arguments:
    logical, contiguous, intent(in)  :: sending(..)
    logical, contiguous, intent(out) :: receiving(..,:)
    integer, optional,   intent(in)  :: length_opt
    integer, optional,   intent(in)  :: procID_opt

    ! Locals:
    integer :: ierr
    integer :: procID, length

    call handleLengthProcID(length_opt, procID_opt, length, procID, rank(sending), size(sending))

    call rpn_comm_gather(sending,   length, 'mpi_logical',  &
                         receiving, length, 'mpi_logical', procID, 'grid', ierr)

    call handleMpiError(ierr, 'mmpi_gather_logical')

  end subroutine mmpi_gather_logical

  !--------------------------------------------------------------------------
  ! mmpi_gather_integer
  !--------------------------------------------------------------------------
  subroutine mmpi_gather_integer(sending, receiving, length_opt, procID_opt)
    !
    !:Purpose: Calling 'rpn_comm_gather' for an integer scalar or array
    !
    implicit none

    ! Arguments:
    integer, contiguous, intent(in)  :: sending(..)
    integer, contiguous, intent(out) :: receiving(..,:)
    integer, optional,   intent(in)  :: length_opt
    integer, optional,   intent(in)  :: procID_opt

    ! Locals:
    integer :: ierr
    integer :: procID, length

    call handleLengthProcID(length_opt, procID_opt, length, procID, rank(sending), size(sending))

    call rpn_comm_gather(sending,   length, 'mpi_integer',  &
                         receiving, length, 'mpi_integer', procID, 'grid', ierr)

    call handleMpiError(ierr, 'mmpi_gather_integer')

  end subroutine mmpi_gather_integer

  !--------------------------------------------------------------------------
  ! mmpi_gather_integer8
  !--------------------------------------------------------------------------
  subroutine mmpi_gather_integer8(sending, receiving, length_opt, procID_opt)
    !
    !:Purpose: Calling 'rpn_comm_gather' for an integer8 scalar or array
    !
    implicit none

    ! Arguments:
    integer(8), contiguous, intent(in)  :: sending(..)
    integer(8), contiguous, intent(out) :: receiving(..,:)
    integer, optional,      intent(in)  :: length_opt
    integer, optional,      intent(in)  :: procID_opt

    ! Locals:
    integer :: ierr
    integer :: procID, length

    call handleLengthProcID(length_opt, procID_opt, length, procID, rank(sending), size(sending))

    call rpn_comm_gather(sending,   length, 'mpi_integer8',  &
                         receiving, length, 'mpi_integer8', procID, 'grid', ierr)

    call handleMpiError(ierr, 'mmpi_gather_integer8')

  end subroutine mmpi_gather_integer8

  !--------------------------------------------------------------------------
  ! mmpi_gather_real4
  !--------------------------------------------------------------------------
  subroutine mmpi_gather_real4(sending, receiving, length_opt, procID_opt)
    !
    !:Purpose: Calling 'rpn_comm_gather' for an real(4) scalar or array
    !
    implicit none

    ! Arguments:
    real(4), contiguous, intent(in)  :: sending(..)
    real(4), contiguous, intent(out) :: receiving(..,:)
    integer, optional,   intent(in)  :: length_opt
    integer, optional,   intent(in)  :: procID_opt

    ! Locals:
    integer :: ierr
    integer :: procID, length

    call handleLengthProcID(length_opt, procID_opt, length, procID, rank(sending), size(sending))

    call rpn_comm_gather(sending,   length, 'mpi_real4',  &
                         receiving, length, 'mpi_real4', procID, 'grid', ierr)

    call handleMpiError(ierr, 'mmpi_gather_real4')

  end subroutine mmpi_gather_real4

  !--------------------------------------------------------------------------
  ! mmpi_gather_real8
  !--------------------------------------------------------------------------
  subroutine mmpi_gather_real8(sending, receiving, length_opt, procID_opt)
    !
    !:Purpose: Calling 'rpn_comm_gather' for an real(8) scalar or array
    !
    implicit none

    ! Arguments:
    real(8), contiguous, intent(in)  :: sending(..)
    real(8), contiguous, intent(out) :: receiving(..,:)
    integer, optional,   intent(in)  :: length_opt
    integer, optional,   intent(in)  :: procID_opt

    ! Locals:
    integer :: ierr
    integer :: procID, length

    call handleLengthProcID(length_opt, procID_opt, length, procID, rank(sending), size(sending))

    call rpn_comm_gather(sending,   length, 'mpi_real8',  &
                         receiving, length, 'mpi_real8', procID, 'grid', ierr)

    call handleMpiError(ierr, 'mmpi_gather_real8')

  end subroutine mmpi_gather_real8

  !--------------------------------------------------------------------------
  ! mmpi_allGather_logical
  !--------------------------------------------------------------------------
  subroutine mmpi_allGather_logical(sending, receiving, length_opt)
    !
    !:Purpose: Calling 'rpn_comm_allGather' for a logical scalar or array
    !
    implicit none

    ! Arguments:
    logical, contiguous, intent(in)  :: sending(..)
    logical, contiguous, intent(out) :: receiving(..,:)
    integer, optional,   intent(in)  :: length_opt

    ! Locals:
    integer :: ierr, length

    call handleLength(length_opt, length, rank(sending), size(sending))

    call rpn_comm_allGather(sending,   length, 'mpi_logical',  &
                            receiving, length, 'mpi_logical', 'grid', ierr)

    call handleMpiError(ierr, 'mmpi_allGather_logical')

  end subroutine mmpi_allGather_logical

  !--------------------------------------------------------------------------
  ! mmpi_allGather_real4
  !--------------------------------------------------------------------------
  subroutine mmpi_allGather_real4(sending, receiving, length_opt)
    !
    !:Purpose: Calling 'rpn_comm_allGather' for a real(4) scalar or array
    !
    implicit none

    ! Arguments:
    real(4), contiguous, intent(in)  :: sending(..)
    real(4), contiguous, intent(out) :: receiving(..,:)
    integer, optional,   intent(in)  :: length_opt

    ! Locals:
    integer :: ierr, length

    call handleLength(length_opt, length, rank(sending), size(sending))

    call rpn_comm_allGather(sending,   length, 'mpi_real4',  &
                            receiving, length, 'mpi_real4', 'grid', ierr)

    call handleMpiError(ierr, 'mmpi_allGather_real4')

  end subroutine mmpi_allGather_real4

  !--------------------------------------------------------------------------
  ! mmpi_allGather_real8
  !--------------------------------------------------------------------------
  subroutine mmpi_allGather_real8(sending, receiving, length_opt)
    !
    !:Purpose: Calling 'rpn_comm_allGather' for a real(8) scalar or array
    !
    implicit none

    ! Arguments:
    real(8), contiguous, intent(in)  :: sending(..)
    real(8), contiguous, intent(out) :: receiving(..,:)
    integer, optional,   intent(in)  :: length_opt

    ! Locals:
    integer :: ierr, length

    call handleLength(length_opt, length, rank(sending), size(sending))

    call rpn_comm_allGather(sending,   length, 'mpi_real8',  &
                            receiving, length, 'mpi_real8', 'grid', ierr)

    call handleMpiError(ierr, 'mmpi_allGather_real8')

  end subroutine mmpi_allGather_real8

  !--------------------------------------------------------------------------
  ! handleLengthProcID
  !--------------------------------------------------------------------------
  subroutine handleLengthProcID(length_opt, procID_opt, length, procID, rankSend, sizeSend)
    !
    !:Purpose: Process 'length_opt' and 'procID_opt' optional arguments
    !
    implicit none

    ! Arguments:
    integer, optional, intent(in)  :: length_opt
    integer, optional, intent(in)  :: procID_opt
    integer,           intent(out) :: length
    integer,           intent(out) :: procID
    integer,           intent(in)  :: rankSend
    integer,           intent(in)  :: sizeSend

    call handleLength(length_opt, length, rankSend, sizeSend)
    call handleProcID(procID_opt, procID)

  end subroutine handleLengthProcID

  !--------------------------------------------------------------------------
  ! handleProcID
  !--------------------------------------------------------------------------
  subroutine handleProcID(procID_opt, procID)
    !
    !:Purpose: Process 'procID_opt' optional argument
    !
    implicit none

    ! Arguments:
    integer, optional, intent(in)  :: procID_opt
    integer,           intent(out) :: procID

    if (present(procID_opt)) then
      procID = procID_opt
    else
      procID = 0
    end if

  end subroutine handleProcID

  !--------------------------------------------------------------------------
  ! handleLength
  !--------------------------------------------------------------------------
  subroutine handleLength(length_opt, length, rankSend, sizeSend)
    !
    !:Purpose: Process 'length_opt' optional argument
    !
    implicit none

    ! Arguments:
    integer, optional, intent(in)  :: length_opt
    integer,           intent(out) :: length
    integer,           intent(in)  :: rankSend
    integer,           intent(in)  :: sizeSend

    if (present(length_opt)) then
      length = length_opt
    else
      if ( rankSend == 0 ) then
        length = 1
      else
        length = sizeSend
      end if
    end if

  end subroutine handleLength

  !--------------------------------------------------------------------------
  ! handleMpiError
  !--------------------------------------------------------------------------
  subroutine handleMpiError(errCode, context)
    !
    !:Purpose: Handle gracefully the MPI error code
    !
    implicit none

    ! Arguments:
    integer,          intent(in) :: errCode
    character(len=*), intent(in) :: context

    ! Locals:
    character(len=MPI_MAX_ERROR_STRING) :: errorMsg
    integer :: resultlen, ierror

    if (errCode /= MPI_SUCCESS) then
      call MPI_Error_string(errcode, errorMsg, resultlen, ierror)
      call utl_abort('MPI error found in ' // context // ' : ' // trim(errorMsg))
    end if
  end subroutine handleMpiError

end module midasMpi_mod

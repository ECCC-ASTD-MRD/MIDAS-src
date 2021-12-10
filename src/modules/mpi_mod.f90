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

module mpi_mod
  !
  ! MODULE mpi_mod (prefix='mpi' category='7. Low-level data objects and utilities')
  !
  ! :Purpose: Subroutine and public variables related to general aspects of mpi.
  !
  use utilities_mod
  implicit none
  save
  private

  ! public variables
  public :: mpi_myid,mpi_nprocs,mpi_myidx,mpi_myidy,mpi_npex,mpi_npey
  public :: mpi_comm_EW, mpi_comm_NS, mpi_comm_GRID, mpi_doBarrier
  public :: mpi_datyp_real4, mpi_datyp_real8, mpi_datyp_int
  ! public procedures
  public :: mpi_initialize,mpi_getptopo
  public :: mpi_allreduce_sumreal8scalar,mpi_allgather_string
  public :: mpi_allreduce_sumR8_1d

  integer :: mpi_myid   = 0
  integer :: mpi_nprocs = 0
  integer :: mpi_myidx  = 0
  integer :: mpi_myidy  = 0
  integer :: mpi_npex   = 0
  integer :: mpi_npey   = 0

  integer :: mpi_comm_EW, mpi_comm_NS, mpi_comm_GRID
  integer :: mpi_datyp_real4, mpi_datyp_real8, mpi_datyp_int

  logical :: mpi_doBarrier = .true.

  contains

  subroutine mpi_initialize()
    implicit none
    integer :: mythread,numthread,omp_get_thread_num,omp_get_num_threads,rpn_comm_mype
    integer :: npex,npey,ierr
    integer :: rpn_comm_comm, rpn_comm_datyp

    ! Initilize MPI
    npex=0
    npey=0
    call rpn_comm_init(mpi_getptopo,mpi_myid,mpi_nprocs,npex,npey)
    if(mpi_nprocs.lt.1) then
      mpi_nprocs=1
      mpi_npex=1
      mpi_npey=1
      mpi_myid=0
      mpi_myidx=0
      mpi_myidy=0
    else
      ierr = rpn_comm_mype(mpi_myid,mpi_myidx,mpi_myidy)
      mpi_npex=npex
      mpi_npey=npey
    endif

    write(*,*) 'mpi_initialize: mpi_myid, mpi_myidx, mpi_myidy = ', mpi_myid, mpi_myidx, mpi_myidy

    !$OMP PARALLEL PRIVATE(numthread,mythread)
    mythread=omp_get_thread_num()
    numthread=omp_get_num_threads()
    if(mythread.eq.0) write(*,*) 'mpi_initialize: NUMBER OF THREADS=',numthread     
    !$OMP END PARALLEL

    ! create standard mpi handles to rpn_comm mpi communicators to facilitate 
    ! use of standard mpi routines
    mpi_comm_EW = rpn_comm_comm('EW')
    mpi_comm_NS = rpn_comm_comm('NS')
    mpi_comm_GRID = rpn_comm_comm('GRID')

    mpi_datyp_real4 = rpn_comm_datyp('MPI_REAL4')
    mpi_datyp_real8 = rpn_comm_datyp('MPI_REAL8')
    mpi_datyp_int = rpn_comm_datyp('MPI_INTEGER')

    write(*,*) ' '
    if(mpi_doBarrier) then
      write(*,*) 'mpi_initialize: MPI_BARRIERs will be done to help with interpretation of timings'
    else
      write(*,*) 'mpi_initialize: no MPI_BARRIERs will be done'
    endif
    write(*,*) ' '

  end subroutine mpi_initialize


  subroutine mpi_getptopo( npex, npey )

    implicit none

    integer, intent(out) :: npex
    integer, intent(out) :: npey

    integer :: ierr
    namelist /ptopo/npex,npey
    integer :: nulnam,fnom,fclos

    npex=1
    npey=1

    nulnam=0
    ierr=fnom(nulnam,'ptopo_nml','FTN+SEQ+R/O',0)
    if(ierr.ne.0) call utl_abort('mpi_getptopo: Error opening file ptopo_nml')
    read(nulnam,nml=ptopo,iostat=ierr)
    if(ierr.ne.0) call utl_abort('mpi_getptopo: Error reading namelist')
    write(*,nml=ptopo)
    ierr=fclos(nulnam)

  end subroutine mpi_getptopo 


  subroutine mpi_allreduce_sumreal8scalar( sendRecvValue, comm )

    implicit none

    real(8), intent(inout)       :: sendRecvValue ! value to be summed over all mpi tasks
    character(len=*), intent(in) :: comm          ! rpn_comm communicator

    integer :: nsize, ierr, root, rank
    real(8), allocatable :: allvalues(:)

    ! do a barrier so that timing on reduce operation is accurate
    call tmg_start(109,'MPI_SUMSCA_BARR')
    if(mpi_doBarrier) call rpn_comm_barrier(comm,ierr)
    call tmg_stop(109)

    call tmg_start(16,'allreduce_sum8')

    ! determine number of processors in the communicating group
    call rpn_comm_size(comm,nsize,ierr)

    ! determine where to gather the values: first task in group
    call rpn_comm_rank(comm,rank,ierr)
    call rpn_comm_allreduce(rank,root,1,"MPI_INTEGER","MPI_MIN",comm,ierr)

    ! gather values to be added onto 1 processor
    allocate(allvalues(nsize))
    call rpn_comm_gather(sendRecvValue, 1, "MPI_DOUBLE_PRECISION", allvalues, 1, "MPI_DOUBLE_PRECISION", root, comm, ierr)

    ! sum the values and broadcast to group
    if(rank.eq.root) sendRecvValue = sum(allvalues(:))
    deallocate(allvalues)
    call rpn_comm_bcast(sendRecvValue, 1, "MPI_DOUBLE_PRECISION", root, comm, ierr)

    call tmg_stop(16)

  end subroutine mpi_allreduce_sumreal8scalar


  subroutine mpi_allreduce_sumR8_1d( sendRecvVector, comm )
    !
    ! :Purpose: Perform sum of 1d array over all MPI tasks.
    !
    implicit none

    real(8)         , intent(inout)  :: sendRecvVector(:) ! 1-D vector to be summed over all mpi tasks
    character(len=*), intent(in)     :: comm              ! rpn_comm communicator

    integer :: nprocs_mpi, numElements, ierr, root, rank
    real(8), allocatable :: all_sendRecvVector(:,:)

    ! do a barrier so that timing on reduce operation is accurate
    call tmg_start(109,'MPI_SUMSCA_BARR')
    if ( mpi_doBarrier ) call rpn_comm_barrier(comm,ierr)
    call tmg_stop(109)

    call tmg_start(16,'allreduce_sum8')

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

    ! sum the values and broadcast to group
    if ( rank == root ) sendRecvVector(:) = sum(all_sendRecvVector(:,:),2)
    deallocate(all_sendRecvVector)
    call rpn_comm_bcast(sendRecvVector, numElements, "mpi_double_precision", &
                        root, comm, ierr)

    call tmg_stop(16)

  end subroutine mpi_allreduce_sumR8_1d
 
  
  subroutine mpi_allgather_string( str_list, str_list_all, nlist, nchar, nproc, comm, ierr )
    ! 
    ! :Purpose: Performs the MPI 'allgather' routine for an array of strings
    !
    implicit none

    character(len=nchar), intent(in) :: str_list(nlist)
    character(len=*)    , intent(in) :: comm
    integer             , intent(in) :: nlist
    integer             , intent(in) :: nchar
    integer             , intent(in) :: nproc
    character(len=nchar), intent(out) :: str_list_all(nlist,nproc)
    integer, intent(out) :: ierr

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

  end subroutine mpi_allgather_string

end module mpi_mod

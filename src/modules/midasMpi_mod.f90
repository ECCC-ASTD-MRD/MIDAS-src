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

module midasMpi_mod
  !
  ! MODULE midasMpi_mod (prefix='mmpi' category='8. Low-level utilities and constants')
  !
  ! :Purpose: Subroutine and public variables related to general aspects of mpi.
  !           Also, subroutine and public variables related to the mpi decomposition
  !           specific to the MIDAS code.
  !
  use utilities_mod
  implicit none
  save
  private

  ! public variables
  public :: mmpi_myid,mmpi_nprocs,mmpi_myidx,mmpi_myidy,mmpi_npex,mmpi_npey,mmpi_numthread
  public :: mmpi_comm_EW, mmpi_comm_NS, mmpi_comm_GRID, mmpi_doBarrier
  public :: mmpi_datyp_real4, mmpi_datyp_real8, mmpi_datyp_int
  ! public procedures
  public :: mmpi_initialize,mmpi_getptopo
  public :: mmpi_allreduce_sumreal8scalar,mmpi_allgather_string
  public :: mmpi_allreduce_sumR8_1d, mmpi_allreduce_sumR8_2d
  public :: mmpi_reduce_sumR8_1d, mmpi_reduce_sumR8_2d, mmpi_reduce_sumR8_3d
  public :: mmpi_setup_latbands, mmpi_setup_lonbands
  public :: mmpi_setup_m, mmpi_setup_n
  public :: mmpi_setup_levels
  public :: mmpi_setup_varslevels
  public :: mmpi_myidXfromLon, mmpi_myidYfromLat

  integer :: mmpi_myid   = 0
  integer :: mmpi_nprocs = 0
  integer :: mmpi_myidx  = 0
  integer :: mmpi_myidy  = 0
  integer :: mmpi_npex   = 0
  integer :: mmpi_npey   = 0
  integer :: mmpi_numthread = 0

  integer :: mmpi_comm_EW, mmpi_comm_NS, mmpi_comm_GRID
  integer :: mmpi_datyp_real4, mmpi_datyp_real8, mmpi_datyp_int

  logical :: mmpi_doBarrier = .true.

  contains

  subroutine mmpi_initialize()
    implicit none
    integer :: mythread,numthread,omp_get_thread_num,omp_get_num_threads,rpn_comm_mype
    integer :: ierr
    integer :: rpn_comm_comm, rpn_comm_datyp

    ! Namelist variables
    integer :: npex
    integer :: npey

    ! Initilize MPI
    npex=0
    npey=0
    call rpn_comm_init(mmpi_getptopo,mmpi_myid,mmpi_nprocs,npex,npey)
    if(mmpi_nprocs.lt.1) then
      mmpi_nprocs=1
      mmpi_npex=1
      mmpi_npey=1
      mmpi_myid=0
      mmpi_myidx=0
      mmpi_myidy=0
    else
      ierr = rpn_comm_mype(mmpi_myid,mmpi_myidx,mmpi_myidy)
      mmpi_npex=npex
      mmpi_npey=npey
    endif

    write(*,*) 'mmpi_initialize: mmpi_myid, mmpi_myidx, mmpi_myidy = ', mmpi_myid, mmpi_myidx, mmpi_myidy

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

    write(*,*) ' '
    if(mmpi_doBarrier) then
      write(*,*) 'mmpi_initialize: MPI_BARRIERs will be done to help with interpretation of timings'
    else
      write(*,*) 'mmpi_initialize: no MPI_BARRIERs will be done'
    endif
    write(*,*) ' '

  end subroutine mmpi_initialize


  subroutine mmpi_getptopo( npex, npey )

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

  end subroutine mmpi_getptopo 


  subroutine mmpi_allreduce_sumreal8scalar( sendRecvValue, comm )

    implicit none

    real(8), intent(inout)       :: sendRecvValue ! value to be summed over all mpi tasks
    character(len=*), intent(in) :: comm          ! rpn_comm communicator

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


  subroutine mmpi_allreduce_sumR8_1d( sendRecvVector, comm )
    !
    ! :Purpose: Perform sum of 1d array over all MPI tasks.
    !
    implicit none

    real(8)         , intent(inout)  :: sendRecvVector(:) ! 1-D vector to be summed over all mpi tasks
    character(len=*), intent(in)     :: comm              ! rpn_comm communicator

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


  subroutine mmpi_allreduce_sumR8_2d( sendRecvVector, comm )
    !
    ! :Purpose: Perform sum of 2d array over all MPI tasks.
    !
    implicit none

    real(8)         , intent(inout)  :: sendRecvVector(:,:) ! 2-D vector to be summed over all mpi tasks
    character(len=*), intent(in)     :: comm                ! rpn_comm communicator

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
 
  
  subroutine mmpi_reduce_sumR8_1d( sendVector, recvVector, root, comm )
    !
    ! :Purpose: Perform sum of 1d array over all MPI tasks.
    !
    implicit none

    real(8)         , intent(in)  :: sendVector(:) ! 1-D vector to be summed over all mpi tasks
    real(8)         , intent(out) :: recvVector(:) ! 1-D vector to be summed over all mpi tasks
    integer         , intent(in)  :: root          ! mpi task id where data is put
    character(len=*), intent(in)  :: comm          ! rpn_comm communicator

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
 
  
  subroutine mmpi_reduce_sumR8_2d( sendVector, recvVector, root, comm )
    !
    ! :Purpose: Perform sum of 2d array over all MPI tasks.
    !
    implicit none

    real(8)         , intent(in)  :: sendVector(:,:) ! 2-D vector to be summed over all mpi tasks
    real(8)         , intent(out) :: recvVector(:,:) ! 2-D vector to be summed over all mpi tasks
    integer         , intent(in)  :: root            ! mpi task id where data will be put
    character(len=*), intent(in)  :: comm            ! rpn_comm communicator

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


  subroutine mmpi_reduce_sumR8_3d( sendVector, recvVector, root, comm )
    !
    ! :Purpose: Perform sum of 3d array over all MPI tasks.
    !
    implicit none

    real(8)         , intent(in)  :: sendVector(:,:,:) ! 3-D vector to be summed over all mpi tasks
    real(8)         , intent(out) :: recvVector(:,:,:) ! 3-D vector to be summed over all mpi tasks
    integer         , intent(in)  :: root              ! mpi task id where data is put
    character(len=*), intent(in)  :: comm              ! rpn_comm communicator

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
 
  
  subroutine mmpi_allgather_string( str_list, str_list_all, nlist, nchar, nproc, comm, ierr )
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

  end subroutine mmpi_allgather_string


  subroutine mmpi_setup_latbands(nj, latPerPE, latPerPEmax, myLatBeg, myLatEnd,  &
                                   myLatHalfBeg_opt, myLatHalfEnd_opt, divisible_opt)
    ! :Purpose: compute parameters that define the mpi distribution of
    !          latitudes over tasks in Y direction (npey)
    implicit none
    integer           :: nj
    integer           :: latPerPE
    integer           :: latPerPEmax
    integer           :: myLatBeg
    integer           :: myLatEnd
    integer           :: njlath
    integer, optional :: myLatHalfBeg_opt
    integer, optional :: myLatHalfEnd_opt
    logical, optional :: divisible_opt

    integer :: latPerPEmin, ierr
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


  function mmpi_myidYfromLat(latIndex, nj) result(IP_y)
    ! :Purpose: use same logic as setup_latbands to compute myidy
    !          corresponding to a latitude grid index
    implicit none
    integer :: latIndex
    integer :: nj
    integer :: IP_y

    IP_y = (latIndex-1) / floor( real(nj) / real(mmpi_npey) )
    IP_y = min( mmpi_npey-1, IP_y )

  end function mmpi_myidYfromLat


  subroutine mmpi_setup_lonbands(ni, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd, divisible_opt)
    ! :Purpose: compute parameters that define the mpi distribution of
    !          longitudes over tasks in X direction (npex)
    implicit none

    integer          :: ni
    integer          :: lonPerPE
    integer          :: lonPerPEmax
    integer          :: myLonBeg
    integer          :: myLonEnd
    logical, optional :: divisible_opt

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


  function mmpi_myidXfromLon(lonIndex, ni) result(IP_x)
    ! :Purpose: use same logic as setup_lonbands to compute myidx
    !          corresponding to a longitude grid index
    implicit none

    integer :: lonIndex
    integer :: ni
    integer :: IP_x

    IP_x = (lonIndex-1) / floor( real(ni) / real(mmpi_npex) )
    IP_x = min( mmpi_npex-1, IP_x )

  end function mmpi_myidXfromLon


  subroutine mmpi_setup_m(ntrunc, mymBeg, mymEnd, mymSkip, mymCount)
    ! :Purpose: compute parameters that define the mpi distribution of
    !          wavenumber m over tasks in Y direction (npey)
    implicit none
    integer :: ntrunc
    integer :: mymBeg
    integer :: mymEnd
    integer :: mymSkip
    integer :: mymCount
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

 
  subroutine mmpi_setup_n(ntrunc, mynBeg, mynEnd, mynSkip, mynCount)
    ! :Purpose: compute parameters that define the mpi distribution of
    !          wavenumber n over tasks in X direction (npex)
    implicit none

    integer :: ntrunc
    integer :: mynBeg
    integer :: mynEnd
    integer :: mynSkip
    integer :: mynCount
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


  subroutine mmpi_setup_levels(numlevels, myLevBeg, myLevEnd, myLevCount)
    ! :Purpose: compute parameters that define the mpi distribution of
    !          levels over tasks in X direction (npex)
    implicit none

    integer :: numlevels
    integer :: myLevBeg
    integer :: myLevEnd
    integer :: myLevCount
    integer :: jlev
    integer :: jproc
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
    do jproc = 1, mmpi_npex
      do jlev = jproc, (numlevels / factor), mmpi_npex
        myLevCounts(jproc) = myLevCounts(jproc) + 1
      end do
    end do
    do jproc = 1, mmpi_npex
      myLevCounts(jproc) = myLevCounts(jproc) * factor
    end do

    myLevCount = myLevCounts(mmpi_myidx + 1)

    if (myLevCount > 0) then
      myLevBeg = 1
      do jproc = 1, mmpi_myidx
        myLevBeg = myLevBeg + myLevCounts(jproc)
      end do
      myLevEnd = myLevBeg + myLevCount - 1
    else
      myLevBeg = 1
      myLevEnd = 0
    end if

    write(*,'(a,3i8)') 'mmpi_setup_levels: myLevBeg, myLevEnd, myLevCount = ',  &
         myLevBeg, myLevEnd, myLevCount

  end subroutine mmpi_setup_levels


  subroutine mmpi_setup_varslevels(numk, mykBeg, mykEnd, mykCount)
    ! :Purpose: compute parameters that define the mpi distribution of
    !          variables/levels (i.e. 1->nk) over all tasks (nprocs)
    implicit none

    integer :: numk
    integer :: mykBeg
    integer :: mykEnd
    integer :: mykCount
    integer :: jk
    integer :: jproc
    integer :: mykCounts(mmpi_nprocs)

    mykCounts(:) = 0
    do jproc = 1, mmpi_nprocs
      do jk = jproc, numk, mmpi_nprocs
        mykCounts(jproc) = mykCounts(jproc) + 1
      end do
    end do

    mykCount = mykCounts(mmpi_myid + 1)

    mykBeg = 1
    do jproc = 1, mmpi_myid
      mykBeg = mykBeg + mykCounts(jproc)
    end do
    mykEnd = mykBeg + mykCount - 1

  end subroutine mmpi_setup_varslevels

end module midasMpi_mod

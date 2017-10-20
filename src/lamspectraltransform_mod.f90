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
!! MODULE LamSpectralTransform (prefix="lst")
!! 
!! *Purpose*: Bi-Fourier spectral transform for limited area applications.
!!            Depends on ffft8 and setfft8 routines in ARMNLIB.
!!
!--------------------------------------------------------------------------
module LamSpectralTransform_mod
  use mpi
  use mpivar_mod
  use MathPhysConstants_mod, only: MPC_RADIANS_PER_DEGREE_R8, MPC_PI_R8
  use earthconstants_mod,    only: RA
  use utilities_mod
  implicit none
  save
  private

  ! public derived type
  public :: struct_lst
  ! public procedures
  public :: lst_Setup, lst_Laplacian, lst_VarTransform
  public :: lst_ReshapeTrunc ! only for standalone tests

  type :: struct_lst
     integer              :: id              ! Transform ID number
     integer              :: nla             ! First dimension of VAR spectral array
     integer              :: nphase          ! Second dimension of VAR spectral array
     real(8), allocatable :: k_r8(:)         ! (real) Total Wavenumber associated with each
                                             !  nla spectral coefficient
     integer, allocatable :: k(:)            ! Total Wavenumber bin associated with each
                                             !  nla spectral coefficient
     integer, allocatable :: m(:)            ! Wavenumber in x associated with each
                                             !  nla spectral coefficient
     integer, allocatable :: n(:)            ! Wavenumber in y associated with each
                                             !  nla spectral coefficient
     real(8), allocatable :: Weight(:)       ! Weight associated with each
                                             !  nla spectral coefficient
     integer, allocatable :: nePerK(:)       ! Number of spectral element in each
                                             !  total wavenumber bands
     integer, allocatable :: ilaFromEK(:,:)  ! ila index associated to each spectral element
                                             !  of total wavenumber band
     real(8), allocatable :: NormFactor(:,:)
     real(8), allocatable :: NormFactorAd(:,:)
     integer              :: nlaGlobal       ! First dimension of VAR global spectral array
     integer, allocatable :: ilaGlobal(:)    ! Position of the local vector element in the global
                                             !  vector
  end type struct_lst

  type :: struct_lst_local
     integer                         :: ni
     integer                         :: nj
     integer                         :: mmax, nmax
     integer                         :: ktrunc
     integer                         :: nla
     integer                         :: maxnla
     integer                         :: nlaGlobal
     integer                         :: mymBeg, mymEnd, mymSkip, mymCount, mymActiveCount,maxmActiveCount
     integer                         :: mynBeg, mynEnd, mynSkip, mynCount
     integer, allocatable            :: nla_Index(:,:)
     real(8), allocatable            :: lapxy(:)
     real(8), allocatable            :: ilapxy(:)
     integer, allocatable            :: allmBeg(:), allmEnd(:), allmSkip(:)
     integer, allocatable            :: mymIndex(:)
     integer, allocatable            :: allnBeg(:), allnEnd(:), allnSkip(:)
     integer, allocatable            :: mynIndex(:)
     logical                         :: allocated = .false.
     integer                         :: latPerPE, myLatBeg, myLatEnd
     integer                         :: lonPerPE, myLonBeg, myLonEnd
     integer                         :: myLevBeg, myLevEnd, myLevCount, maxLevCount
     integer, allocatable            :: allLatBeg(:), allLatEnd(:)
     integer, allocatable            :: allLonBeg(:), allLonEnd(:)
     integer, allocatable            :: allLevBeg(:), allLevEnd(:)
     integer, allocatable            :: KfromMNglb(:,:)
     character(len=10)               :: MpiMode
     character(len=3)                :: gridDataOrder ! Ordering the gridded data: 'ijk' or 'kij'
     integer                         :: sendType_LevToLon, recvType_LevToLon
     integer                         :: sendType_LonToLev, recvType_LonToLev
  end type struct_lst_local

  integer,parameter      :: nMaxLst = 10
  integer                :: nLstAlreadyAllocated = 0
  type(struct_lst_local) :: lst(nMaxLst)

  character(len=*), parameter     :: TransformType = 'SinCos'
  integer, parameter              :: nip = 2     ! Padding
  integer, parameter              :: njp = 2     ! Padding
  integer, parameter              :: nphase = 4  ! For Sin&Cos we have a) sin(m)sin(n)
                                                 !                     b) cos(m)cos(n)
                                                 !                     c) sin(m)cos(n)
                                                 !                     d) cos(m)sin(n)

  logical, parameter :: verbose = .false.

contains

!--------------------------------------------------------------------------
! LST_SETUP
!--------------------------------------------------------------------------
  subroutine lst_Setup( lst_out, ni_in, nj_in, dlon_in, ktrunc_in,    &
                        MpiMode, maxlevels_in, gridDataOrder )
    implicit none

    integer,          intent(in)    :: ni_in, nj_in    
                                     ! Global Grid point data horizontal dimensions
    character(len=*), intent(in)    :: MpiMode
                                     ! MPI Strategy
    real(8),          intent(in)    :: dlon_in
                                     ! Grid Spacing in Radians
    integer,          intent(in)    :: ktrunc_in
                                     ! Spectral Truncation (global)
    integer, intent(in), optional   :: maxlevels_in
                                     ! Number of levels; Only needed when MpiMode = LatLev
    character(len=*), intent(in), optional :: gridDataOrder
                                     ! 'ijk' or 'kij'
    type(struct_lst), intent(out)   :: lst_out
                                     ! Parameters available to the outside world

    real(8), allocatable            :: Kr8fromMN(:,:)

    integer, allocatable            :: KfromMN(:,:)
    integer, allocatable            :: my_KfromMNglb(:,:)
    integer                         :: kref, mref, nref, id
    integer                         :: m, n, k, ila, nfact
    integer                         :: ier, ilaglb, i, j, p
    real(8)                         :: a, b, r
    real(8)                         :: dlon, dx2, fac, ca, cp, cb, cq
    real(8)                         :: NormFactor1, NormFactor2, NormFactor3
    real(8)                         :: NormFactorAd1, NormFactorAd2, NormFactorAd3
    real(8)                         :: factor, factorAd
    character(len=60)               :: kreftype

    integer(kind=MPI_ADDRESS_KIND) :: lowerBound, extent
    integer :: realSize, sendType, recvType, ierr

    !
    !- 1.  Set variables needed by the LAM Spectral Transform in VAR
    !
    if (verbose) write(*,*) 'Entering lst_Setup'

    kreftype = 'MAX' ! hardwired

    nLstAlreadyAllocated = nLstAlreadyAllocated + 1
    if (nLstAlreadyAllocated <= nMaxLst) then
       id = nLstAlreadyAllocated
       lst_out%id = id
       write(*,*)
       write(*,*) "lst_setup: Setting transform id = ",id
    else
       call utl_abort('lst_setup: Too many Spectral Transforms!!!')
    end if

    lst(id)%ni = ni_in
    lst(id)%nj = nj_in
    lst(id)%ktrunc = ktrunc_in

    !  1.1 Check grid dimensions and set padding for the RPN DFT routines

    ! We need to padd the input array such as ...                              
    ! O O O O O O O O O
    ! O O O O O O O O O
    ! X X X X X X X O O
    ! X X X X X X X O O
    ! X X X X X X X O O
    ! X X X X X X X O O

    if (mod(lst(id)%ni,2) /= 0 .or. mod(lst(id)%nj,2) /= 0) then
       write(6,*) ' The regular Sin & Cos Transform requires that', &
                  ' dimensions be EVEN. Fields MUST be periodic' , &
                  ' but the last colum and row MUST NOT BE a '   , &
                  ' repetition of the first colum and row. '
       call utl_abort('lst_setup')
    end if

    nfact = lst(id)%ni
    call ngfft( nfact ) ! INOUT
    if ( nfact /= lst(id)%ni ) then
       write(6,*) 'Error: A fast transform cannot be used in X'
       write(6,6130) lst(id)%ni, nfact
       call utl_abort('lst_setup')
    end if

    nfact = lst(id)%nj
    call ngfft( nfact ) ! INOUT
    if ( nfact /= lst(id)%nj ) then
       write(6,*) 'Error: A fast transform cannot be used in Y'
       write(6,6140) lst(id)%nj, nfact
       call utl_abort('lst_setup')
    end if

6130 FORMAT('N = ni = ', I4,' the nearest factorizable N = ',I4)
6140 FORMAT('N = nj = ', I4,' the nearest factorizable N = ',I4)

    !  1.2 Set the number of phases
    lst_out%nphase = nphase

    !- 1.3 Maximum of integer wavenumbers in x and y directions
    lst(id)%mmax = lst(id)%ni/2
    lst(id)%nmax = lst(id)%nj/2

    write(*,'(A,f8.1)') ' lst_Setup: Your grid spacing (in km) = ', RA*dlon_in/1000.0
    write(*,*) '           Max wavenumbers in x-axis = ', lst(id)%mmax            
    write(*,*) '           Max wavenumbers in y-axis = ', lst(id)%nmax
    if      ( lst(id)%ktrunc > nint(sqrt(real(lst(id)%mmax)**2+real(lst(id)%nmax)**2)) ) then
      write(*,*)
      write(*,*) 'lst_Setup: Warning: Truncation is larger than sqrt(mmax^2+nmax^2)'
      write(*,*) '           NO TRUNCATION will be applied'
    else if ( lst(id)%ktrunc > min(lst(id)%mmax,lst(id)%nmax) ) then
      write(*,*)
      if ( lst(id)%ktrunc > lst(id)%mmax ) then
         write(*,*) 'lst_Setup: Warning: Truncation is larger than mmax only'
         write(*,*) '           TRUNCATION will be applied only above nmax'
         write(*,'(A,f8.1)') '          i.e., for wavelenght (in km) in y-axis smaller than ',&
                     (lst(id)%nj*RA*dlon_in/1000.0)/lst(id)%ktrunc
      else
         write(*,*) 'lst_Setup: Warning: Truncation is larger than nmax only'
         write(*,*) '           TRUNCATION will be applied only above mmax'
         write(*,'(A,f8.1)') '          i.e., for wavelenght (in km) in x-axis smaller than ',&
                     (lst(id)%ni*RA*dlon_in/1000.0)/lst(id)%ktrunc
      end if
    else
      write(*,*)
      write(*,*) 'lst_Setup: TRUNCATION will be applied above k = ',lst(id)%ktrunc
      write(*,'(A,f8.1)') '          i.e., for wavelenght (in km) in x-axis smaller than ',&
                     (lst(id)%ni*RA*dlon_in/1000.0)/lst(id)%ktrunc
      write(*,'(A,f8.1)') '          i.e., for wavelenght (in km) in y-axis smaller than ',&
                     (lst(id)%nj*RA*dlon_in/1000.0)/lst(id)%ktrunc
    end if

    !- 1.4 MPI Strategy

    lst(id)%MpiMode = MpiMode
    select case ( trim(lst(id)%MpiMode) )
    case ('NoMpi')
       !- 1.4.1 No MPI

       ! range of LONS handled by this processor (ALL) in GRIDPOINT SPACE
       lst(id)%lonPerPE   = lst(id)%ni
       lst(id)%myLonBeg   = 1
       lst(id)%myLonEnd   = lst(id)%ni

       ! range of LATS handled by this processor (ALL) in GRIDPOINT SPACE
       lst(id)%latPerPE   = lst(id)%nj
       lst(id)%myLatBeg   = 1
       lst(id)%myLatEnd   = lst(id)%nj

       ! range of M handled by this processor (ALL) in SPECTRAL SPACE
       lst(id)%mymBeg     = 0
       lst(id)%mymEnd     = lst(id)%mmax
       lst(id)%mymSkip    = 1
       lst(id)%mymCount   = lst(id)%mmax + 1

       ! range of N handled by this processor (ALL) in SPECTRAL SPACE
       lst(id)%mynBeg     = 0
       lst(id)%mynEnd     = lst(id)%nmax
       lst(id)%mynSkip    = 1
       lst(id)%mynCount   = lst(id)%nmax + 1

       ! set a dummy range of LEVELS handled by this processor
       lst(id)%myLevBeg   = -1
       lst(id)%myLevEnd   = -1
       lst(id)%myLevCount = 0

    case ('LatLonMN')
       !- 1.4.2 MPI 2D: Distribution of lon/lat tiles (gridpoint space) and n/m (spectral space)

       ! range of LONS handled by this processor in GRIDPOINT SPACE
       call mpivar_setup_lonbands( lst(id)%ni,                                       & ! IN
                                   lst(id)%lonPerPE,lst(id)%myLonBeg,lst(id)%myLonEnd) ! OUT

       ! range of LATS handled by this processor in GRIDPOINT SPACE
       call mpivar_setup_latbands( lst(id)%nj,                                       & ! IN
                                   lst(id)%latPerPE,lst(id)%myLatBeg,lst(id)%myLatEnd) ! OUT

       ! range of M handled by this processor in SPECTRAL SPACE
       call mpivar_setup_m( lst(id)%mmax,                                                     & ! IN
                            lst(id)%mymBeg, lst(id)%mymEnd, lst(id)%mymSkip, lst(id)%mymCount ) ! OUT

       ! range of N handled by this processor in SPECTRAL SPACE
       call mpivar_setup_n( lst(id)%nmax,                                                     & ! IN
                            lst(id)%mynBeg, lst(id)%mynEnd, lst(id)%mynSkip, lst(id)%mynCount ) ! OUT

       ! range of LEVELS TEMPORARILY handled by this processor DURING THE SPECTRAL TRANSFORM
       if ( .not.present(maxlevels_in) ) then
          call utl_abort('lst_setup: ERROR, number of levels must be specified with MpiMode LatLonMN')
       end if
       ! 2D MPI decomposition: split levels across npex
       call mpivar_setup_levels_npex( maxlevels_in,                                        & ! IN
                                      lst(id)%myLevBeg,lst(id)%myLevEnd,lst(id)%myLevCount ) ! OUT

    case default
       write(*,*)
       write(*,*) 'Error: MpiMode Unknown ', trim(MpiMode)
       call utl_abort('lst_setup')
    end select

    write(*,*)
    write(*,*) ' I am processor ', mpi_myid+1, ' on a total of ', mpi_nprocs
    write(*,*) '          mband info = ', lst(id)%mymBeg, lst(id)%mymEnd, lst(id)%mymSkip, lst(id)%mymCount
    write(*,*) '          nband info = ', lst(id)%mynBeg, lst(id)%mynEnd, lst(id)%mynSkip, lst(id)%mynCount
    write(*,*) '          level info = ', lst(id)%myLevBeg, lst(id)%myLevEnd, lst(id)%myLevCount

    ! Set M index
    allocate(lst(id)%mymIndex(lst(id)%mymBeg:lst(id)%mymEnd))
    lst(id)%mymIndex(:)=0
    do m = lst(id)%mymBeg, lst(id)%mymEnd, lst(id)%mymSkip
      if (m == lst(id)%mymBeg ) then
        lst(id)%mymIndex(m) = 1
      else
        lst(id)%mymIndex(m) = lst(id)%mymIndex(m-lst(id)%mymSkip) + 1
      end if
    end do
    
    ! Set N index
    allocate(lst(id)%mynIndex(lst(id)%mynBeg:lst(id)%mynEnd))
    lst(id)%mynIndex(:)=0
    do n = lst(id)%mynBeg, lst(id)%mynEnd, lst(id)%mynSkip
      if (n == lst(id)%mynBeg ) then
        lst(id)%mynIndex(n) = 1
      else
        lst(id)%mynIndex(n) = lst(id)%mynIndex(n-lst(id)%mynSkip) + 1
      end if
    end do

    if ( trim(lst(id)%MpiMode) /= 'NoMpi') then

      ! Gathering with respect to Longitude
      allocate(lst(id)%allLonBeg(mpi_npex))
      call rpn_comm_allgather(lst(id)%myLonBeg ,1,"mpi_integer",       &
                              lst(id)%allLonBeg,1,"mpi_integer","EW",ier)
      if ( mpi_myid == 0 ) write(*,*) 'AllLonBeg =', lst(id)%allLonBeg(:)

      allocate(lst(id)%allLonEnd(mpi_npex))
      call rpn_comm_allgather(lst(id)%myLonEnd ,1,"mpi_integer",       &
                              lst(id)%allLonEnd,1,"mpi_integer","EW",ier)
      if ( mpi_myid == 0 ) write(*,*) 'AllLonEnd =', lst(id)%allLonEnd(:)

      ! Gathering with respect to Latitude
      allocate(lst(id)%allLatBeg(mpi_npey))
      call rpn_comm_allgather(lst(id)%myLatBeg ,1,"mpi_integer",       &
                              lst(id)%allLatBeg,1,"mpi_integer","NS",ier)
      if ( mpi_myid == 0 ) write(*,*) 'AllLatBeg =', lst(id)%allLatBeg(:)

      allocate(lst(id)%allLatEnd(mpi_npey))
      call rpn_comm_allgather(lst(id)%myLatEnd ,1,"mpi_integer",       &
                              lst(id)%allLatEnd,1,"mpi_integer","NS",ier)
      if ( mpi_myid == 0 ) write(*,*) 'AllLatEnd =', lst(id)%allLatEnd(:)

      ! Gathering with respect to M
      allocate(lst(id)%allmBeg(mpi_npey))
      call rpn_comm_allgather(lst(id)%mymBeg ,1,"mpi_integer",       &
                              lst(id)%allmBeg,1,"mpi_integer","NS",ier)
      if ( mpi_myid == 0 ) write(*,*) 'AllmBeg =', lst(id)%allmBeg(:)

      allocate(lst(id)%allmEnd(mpi_npey))
      call rpn_comm_allgather(lst(id)%mymEnd ,1,"mpi_integer",       &
                              lst(id)%allmEnd,1,"mpi_integer","NS",ier)
      if ( mpi_myid == 0 ) write(*,*) 'allmEnd =', lst(id)%allmEnd(:)
    
      allocate(lst(id)%allmSkip(mpi_npey))
      call rpn_comm_allgather(lst(id)%mymSkip ,1,"mpi_integer",       &
                              lst(id)%allmSkip,1,"mpi_integer","NS",ier)
      if ( mpi_myid == 0 ) write(*,*) 'allmSkip = ', lst(id)%allmSkip(:)

      allocate(lst(id)%allnBeg(mpi_npex))
      call rpn_comm_allgather(lst(id)%mynBeg ,1,"mpi_integer",       &
                              lst(id)%allnBeg,1,"mpi_integer","EW",ier)
      if ( mpi_myid == 0 ) write(*,*) 'AllnBeg =', lst(id)%allnBeg(:)

      allocate(lst(id)%allnEnd(mpi_npex))
      call rpn_comm_allgather(lst(id)%mynEnd ,1,"mpi_integer",       &
                              lst(id)%allnEnd,1,"mpi_integer","EW",ier)
      if ( mpi_myid == 0 ) write(*,*) 'AllnEnd =', lst(id)%allnEnd(:)
    
      allocate(lst(id)%allnSkip(mpi_npex))
      call rpn_comm_allgather(lst(id)%mynSkip ,1,"mpi_integer",       &
                              lst(id)%allnSkip,1,"mpi_integer","EW",ier)
      if ( mpi_myid == 0 ) write(*,*) 'AllnSkip = ', lst(id)%allnSkip(:)

      ! Gathering with respect to levels
      call rpn_comm_allreduce(lst(id)%myLevCount,lst(id)%maxLevCount, &
                                1,"MPI_INTEGER","MPI_MAX","GRID",ier)
      if ( mpi_myid == 0 ) write(*,*) 'MaxLevCount =',lst(id)%maxLevCount

      allocate(lst(id)%allLevBeg(mpi_npex))
      call rpn_comm_allgather(lst(id)%myLevBeg ,1,"mpi_integer",       &
                              lst(id)%allLevBeg,1,"mpi_integer","EW",ier)
      if ( mpi_myid == 0 ) write(*,*) 'AllLevBeg =', lst(id)%allLevBeg(:)

      allocate(lst(id)%allLevEnd(mpi_npex))
      call rpn_comm_allgather(lst(id)%myLevEnd ,1,"mpi_integer",       &
                              lst(id)%allLevEnd,1,"mpi_integer","EW",ier)
      if ( mpi_myid == 0 ) write(*,*) 'AllLevEnd =', lst(id)%allLevEnd(:)

      ! Setup mpi derived types used in transposes
      ! ... mpi_type_vector(count, blocklength, stride, ...)
      ! ... mpi_type_create_resized(oldtype, lowerbound, extent(in bytes), newtype, ierr)
   
      call mpi_type_size(MPI_REAL8, realSize, ierr)
      lowerBound = 0

      ! create the send type for LevToLon
      extent = lst(id)%maxLevCount * lst(id)%lonPerPE * realSize
      call mpi_type_vector(lst(id)%latPerPE, lst(id)%maxLevCount * lst(id)%lonPerPE,  &
           lst(id)%maxLevCount * lst(id)%ni, MPI_REAL8, sendtype, ierr)
      call mpi_type_create_resized(sendtype, lowerBound , extent, lst(id)%sendType_LevToLon, ierr);
      call mpi_type_commit(lst(id)%sendType_LevToLon,ierr)

      ! create the receive type for LevToLon
      extent = lst(id)%maxLevCount * realSize
      call mpi_type_vector(lst(id)%lonPerPE * lst(id)%latPerPE , lst(id)%maxLevCount,  &
           maxlevels_in, MPI_REAL8, recvtype, ierr);
      call mpi_type_create_resized(recvtype, lowerBound, extent, lst(id)%recvType_LevToLon, ierr);
      call mpi_type_commit(lst(id)%recvType_LevToLon, ierr)

      ! create the send type for LonToLev
      extent = lst(id)%maxLevCount * realSize
      call mpi_type_vector(lst(id)%lonPerPE * lst(id)%latPerPE , lst(id)%maxLevCount,  &
           maxlevels_in, MPI_REAL8, sendtype, ierr);
      call mpi_type_create_resized(sendtype, lowerBound, extent, lst(id)%sendType_LonToLev, ierr);
      call mpi_type_commit(lst(id)%sendType_LonToLev, ierr)
      
      ! create the recv type for LonToLev
      extent = lst(id)%maxLevCount * lst(id)%lonPerPE * realSize
      call mpi_type_vector(lst(id)%latPerPE, lst(id)%maxLevCount * lst(id)%lonPerPE,  &
           lst(id)%maxLevCount * lst(id)%ni, MPI_REAL8, recvtype, ierr)
      call mpi_type_create_resized(recvtype, lowerBound , extent, lst(id)%recvType_LonToLev, ierr);
      call mpi_type_commit(lst(id)%recvType_LonToLev,ierr)
      
     end if

    !- 1.5 Compute the Total Wavenumber associated with weach m,n pairs and
    !      the number of spectral element in the VAR array (nla) 
    !      FOR THE LOCAL PROCESSOR
    allocate( Kr8fromMN(0:lst(id)%mmax,0:lst(id)%nmax) )
    Kr8FromMN(:,:) = -1.d0

    allocate( KfromMN(0:lst(id)%mmax,0:lst(id)%nmax) )
    KFromMN(:,:) = -1

    ! Denis et al., MWR, 2002
    !mref = lst(id)%mmax
    !nref = lst(id)%nmax

    ! old var code (L. Fillion)
    mref = lst(id)%ni-1
    nref = lst(id)%nj-1

    select case(trim(kreftype))
    case ('MAX','max')
       ! old var code (L. Fillion)
       write(*,*)
       write(*,*) 'lst_Setup: Using legacy (old var code) total wavenumber normalization'
       kref = max(mref,nref)
    case ('MIN','min')
       ! Denis et al., MWR, 2002
       write(*,*)
       write(*,*) 'lst_Setup: Using Denis et al. total wavenumber normalization'
       kref = min(mref,nref)
    case default
       write(*,*)
       write(*,*) 'Unknown KREFTYPE in lst_setup : ', trim(kreftype)
       call utl_abort('lst_setup')
    end select

    ila = 0
    do n = lst(id)%mynBeg, lst(id)%mynEnd, lst(id)%mynSkip
      do m = lst(id)%mymBeg, lst(id)%mymEnd, lst(id)%mymSkip
         a = real(m,8)/real(mref,8)
         b = real(n,8)/real(nref,8)
         r = real(kref,8) * sqrt( (a**2) + (b**2) ) ! Ellipse Shape if nref /= mref
         k = nint(r) ! or ceiling(r) as in Denis et al. ?
         if ( k <= lst(id)%ktrunc ) then
            ila = ila +1
            KFromMN(m,n) = k
            Kr8FromMN(m,n) = r
         end if
      end do
    end do

    if ( ila == 0 ) then
       write(*,*)
       write(*,*) 'There are no spectral elements associated to this mpi task!'
    !   write(*,*) 'Your options: decrease the number of mpi tasks and/or increase the truncation'
    !   call utl_abort('lst_setup')
    end if

    lst(id)%nla = ila
    lst_out%nla = lst(id)%nla ! Number of spectral element per phase in the VAR array
    if ( trim(lst(id)%MpiMode) /= 'NoMpi') then
       call rpn_comm_allreduce(lst(id)%nla,lst(id)%maxnla, &
            1,"MPI_INTEGER","MPI_MAX","GRID",ier)
       if ( mpi_myid == 0 ) write(*,*) 'MaxNLA =',lst(id)%maxnla
    end if

    allocate( lst(id)%KfromMNglb(0:lst(id)%mmax,0:lst(id)%nmax) )
    allocate( my_KfromMNglb(0:lst(id)%mmax,0:lst(id)%nmax) )
    my_KfromMNglb = 0
    my_KFromMNglb(:,:) = KFromMN(:,:)
    if ( trim(lst(id)%MpiMode) /= 'NoMpi') then
      call rpn_comm_allreduce(my_KFromMNglb,lst(id)%KFromMNglb, &
                                (lst(id)%mmax+1)*(lst(id)%nmax+1),"MPI_INTEGER","MPI_MAX","GRID",ier)
    end if
    deallocate(my_KfromMNglb) 

    lst(id)%mymActiveCount=0
    do m = lst(id)%mymBeg, lst(id)%mymEnd, lst(id)%mymSkip
      if ( KfromMN(m,0) /= -1 ) lst(id)%mymActiveCount = lst(id)%mymActiveCount + 1
    end do
    if ( trim(lst(id)%MpiMode) /= 'NoMpi') then
       call rpn_comm_allreduce(lst(id)%mymActiveCount,lst(id)%maxmActiveCount, &
                               1,"MPI_INTEGER","MPI_MAX","GRID",ier)
       if ( mpi_myid == 0 ) write(*,*) 'MaxmActiveCount =',lst(id)%maxmActiveCount
    end if

    !- 1.6 VAR spectral element ordering &
    !      Total Wavenumbers and Weights associated with each spectral element
    !      FOR THE LOCAL PROCESSOR
    allocate( lst(id)%nla_Index(0:lst(id)%mmax,0:lst(id)%nmax) )

    allocate( lst_out%k_r8(1:lst_out%nla) )
    allocate( lst_out%k(1:lst_out%nla) )
    allocate( lst_out%m(1:lst_out%nla) )
    allocate( lst_out%n(1:lst_out%nla) )
    allocate( lst_out%Weight(1:lst_out%nla) )
    allocate( lst_out%nePerK(0:lst(id)%ktrunc))
    allocate( lst_out%ilaFromEK(1:lst_out%nla,0:lst(id)%ktrunc))
    allocate( lst_out%ilaGlobal(1:lst_out%nla) )

    lst(id)%nla_Index(:,:) = -1
    lst_out%ilaFromEK(:,:) = -1
    lst_out%NEPerK(:)      =  0

    ila    = 0
    ilaglb = 0
    do n = 0, lst(id)%nmax
       do m = 0, lst(id)%mmax
        k    = KfromMN(m,n)

        if ( lst(id)%KfromMNglb(m,n) /= -1 ) ilaglb = ilaglb + 1 ! Global Index

        if ( k /= -1 ) then
          ila = ila+1

          ! Internal index
          lst(id)%nla_Index(m,n) = ila

          ! Outgoing (public) variables
          lst_out%nePerK(k) = lst_out%nePerK(k) + 1
          lst_out%ilaFromEK(lst_out%nePerK(k),k) = ila
          lst_out%k_r8(ila) = Kr8fromMN(m,n)
          lst_out%k(ila) = k
          lst_out%m(ila) = m
          lst_out%n(ila) = n
          lst_out%ilaGlobal(ila) = ilaglb

          ! Spectral coefficient weight associated with this index
          if ( m == 0 .and. n == 0) then
            lst_out%Weight(ila) = 1.0d0
          else if (m /= 0 .and. n /= 0) then
            lst_out%Weight(ila) = 4.0d0
          else
            lst_out%Weight(ila) = 2.0d0
          end if

       end if

      end do
    end do

    lst_out%nlaGlobal = ilaglb ! Number of spectral element per phase in the VAR mpi global array

    deallocate( Kr8fromMN )
    deallocate( KfromMN )

    !- 1.7 Gridded data ordering (input/output)
    if (present(gridDataOrder)) then
       lst(id)%gridDataOrder = trim(gridDataOrder)
    else
       lst(id)%gridDataOrder = 'ijk' ! default value
    end if

    select case ( trim(lst(id)%gridDataOrder) )
    case ('ijk')
       write(*,*) 'lst_setup: gridded data ordering = IJK' 
    case ('kij')
       write(*,*) 'lst_setup: gridded data ordering = KIJ'
    case default
       write(*,*)
       write(*,*) 'Error: gridDataOrder Unknown ', trim(gridDataOrder)
       call utl_abort('lst_setup')
    end select

    !
    !- 2.  Set factors for parseval identity
    !
    allocate( lst_out%NormFactor  (lst(id)%nla,nphase))
    allocate( lst_out%NormFactorAd(lst(id)%nla,nphase))

    Normfactor1   = 1.0d0
    Normfactor2   = 0.5d0 * sqrt(2.0d0)
    Normfactor3   = 0.5d0
    NormfactorAd1 =      1.0d0  * real((lst(id)%ni * lst(id)%nj),8)
    NormfactorAd2 = sqrt(2.0d0) * real((lst(id)%ni * lst(id)%nj),8)
    NormfactorAd3 =      2.0d0  * real((lst(id)%ni * lst(id)%nj),8)

    do ila = 1,lst(id)%nla

      m = lst_out%m(ila)
      n = lst_out%n(ila)

      do p = 1, nphase
         if      ( p == 1) then
            i = 2*m+1
            j = 2*n+1
         else if ( p == 2) then
            i = 2*m+1
            j = 2*n+2
         else if ( p == 3) then
            i = 2*m+2
            j = 2*n+1
         else if ( p == 4) then
            i = 2*m+2
            j = 2*n+2
         else
            call utl_abort('lst_Setup: Error in NormFactor')
         end if

         if ( i == 1 .or. j == 1) then  
            if ( i == 1 .and. j == 1) then
               factor   = Normfactor1
               factorAd = NormfactorAd1
            else
               factor   = Normfactor2
               factorAd = NormfactorAd2
            end if
         else
            factor   = Normfactor3
            factorAd = NormfactorAd3
         end if
         
         lst_out%NormFactor  (ila,p) = factor
         lst_out%NormFactorAd(ila,p) = factorAd
      end do
      
   end do

    !
    !- 3.  Set variables needed by Forward and Inverse Laplacian
    !
    allocate(lst(id)%lapxy (lst(id)%nla))
    allocate(lst(id)%ilapxy(lst(id)%nla))

    dlon = dlon_in
    dx2  = (RA*dlon)**2
    fac  = 2.d0/dx2

    do ila = 1,lst(id)%nla
      ca = 2.d0*MPC_PI_R8 * lst_out%m(ila)
      cp = cos(ca/lst(id)%ni)
      cb = 2.d0*MPC_PI_R8 * lst_out%n(ila)
      cq = cos(cb/lst(id)%nj)

      lst(id)%lapxy(ila) = fac * (cp + cq - 2.d0)
      if ( lst(id)%lapxy(ila) /= 0.d0 ) then
         lst(id)%ilapxy(ila) = 1.d0 / lst(id)%lapxy(ila)
      else
         lst(id)%ilapxy(ila) = 0.d0
      end if
    end do

    !
    !- 4. Finalized
    !
    lst(id)%allocated = .true.

  end subroutine lst_Setup

!--------------------------------------------------------------------------
! LST_VARTRANSFORM
!--------------------------------------------------------------------------
  subroutine lst_VarTransform( id, SpectralStateVar, GridState,        &
                               TransformDirection, nk)
    implicit none

    integer,          intent(in)    :: id
                                     ! LST ID
    integer,          intent(in)    :: nk
                                     ! Grid point data dimensions
    character(len=*), intent(in)    :: TransformDirection
                                     ! SpectralToGridPoint or
                                     ! GridPointToSpectral
    real(8),          intent(inout) :: GridState(:,:,:)
                                     ! 3D field in grid point space
    real(8),          intent(inout) :: SpectralStateVar(:,:,:)
                                     ! 3D spectral coefficients

    !
    !- 0. Some tests...
    !
    if (verbose) write(*,*) 'Entering lst_varTransform'

    call idcheck(id)

    !
    !- 1. Call the appropriate transform
    !
    if ( trim(lst(id)%gridDataOrder) == 'ijk' ) then
       call lst_VarTransform_ijk( id, SpectralStateVar, GridState, &
                                  TransformDirection, nk)
    else
       call lst_VarTransform_kij( id, SpectralStateVar, GridState, &
                                  TransformDirection, nk)
    end if

  end subroutine lst_VarTransform

!--------------------------------------------------------------------------
! LST_VARTRANSFORM_IJK
!--------------------------------------------------------------------------
  subroutine lst_VarTransform_ijk( id, SpectralStateVar, GridState,        &
                                   TransformDirection, nk)
    implicit none

    integer,          intent(in)    :: id
                                     ! LST ID
    integer,          intent(in)    :: nk
                                     ! Grid point data dimensions
    character(len=*), intent(in)    :: TransformDirection
                                     ! SpectralToGridPoint or
                                     ! GridPointToSpectral
    real(8),          intent(inout) :: GridState(lst(id)%myLonBeg:lst(id)%myLonEnd,lst(id)%myLatBeg:lst(id)%myLatEnd,nk)
                                     ! 3D field in grid point space
    real(8),          intent(inout) :: SpectralStateVar(:,:,:)
                                     ! 3D spectral coefficients

    integer                         :: m, n, k, ni_l, nj_l, nip_l, njp_l
    integer                         :: iStart, iEnd, jStart, jEnd, kStart, kEnd

    integer                         :: i, j

    real(8), allocatable            :: Step0(:,:,:)
    real(8), allocatable            :: Step1(:,:,:)
    real(8), allocatable            :: Step2(:,:,:)
    real(8), allocatable            :: Step3(:,:,:)

    character(len=1)                :: TransformAxe

    !
    !- 1. Data pre-processing (Input -> Step0)
    !
    if (verbose) write(*,*) 'Entering lst_varTransform_ijk'

    !- 1.1 Settings and Data Selection

    if ( trim(TransformDirection) == 'GridPointToSpectral' ) then
       iStart = 1 
       iEnd   = lst(id)%ni
       jStart = lst(id)%myLatBeg
       jEnd   = lst(id)%myLatEnd
       if ( trim(lst(id)%MpiMode) == 'NoMpi') then
         kStart= 1
         kEnd  = nk
       else
         kStart= lst(id)%myLevBeg
         kEnd  = lst(id)%myLevEnd
       end if
       allocate( Step0(iStart:iEnd,jStart:jEnd,kStart:kEnd) )
       if ( trim(lst(id)%MpiMode) == 'NoMpi') then
         Step0(:,:,:) = GridState(:,:,:)
       else
         call transpose2d_LonToLev( Step0,            & ! OUT
                                    GridState, nk, id ) ! IN
       end if
    end if

    !
    !- 1.  First pass (Step0 -> Step1)
    !

    !- 1.1 Settings and Data Selection
    select case ( trim(TransformDirection) )
    case ('GridPointToSpectral')
       TransformAxe = 'i'
       ni_l  = lst(id)%ni
       nip_l = nip
       nj_l  = lst(id)%latPerPE
       njp_l = 0
       if ( trim(lst(id)%MpiMode) == 'NoMpi') then
         kStart= 1
         kEnd  = nk
       else
         kStart= lst(id)%myLevBeg
         kEnd  = lst(id)%myLevEnd
       end if
    case ('SpectralToGridPoint')
       TransformAxe = 'j'
       ni_l  = 2*lst(id)%mymCount
       nip_l = 0
       nj_l  = lst(id)%nj
       njp_l = njp
       if ( trim(lst(id)%MpiMode) == 'NoMpi') then
         kStart= 1
         kEnd  = nk
       else
         kStart= lst(id)%myLevBeg
         kEnd  = lst(id)%myLevEnd
       end if
    case default
       write(*,*)
       write(*,*) 'Error: TranformDirection Unknown ', trim(TransformDirection)
       call utl_abort('lst_VarTransform')
    end select

    allocate( Step1(ni_l+nip_l,nj_l+njp_l,kStart:kEnd) )

    !- 1.2 Spectral transform
    if ( trim(TransformDirection) == 'SpectralToGridPoint' ) then
       if ( trim(lst(id)%MpiMode) == 'NoMpi' ) then
         call lst_ReshapeTrunc( Step1,                       & ! OUT
                                SpectralStateVar,            & ! IN
                                'ToRPN', kStart, kEnd, id )    ! IN
       else
         call transpose2d_NToLev( Step1,                   & ! OUT
                                  SpectralStateVar, nk, id ) ! IN
       end if
    else
       Step1(1:lst(id)%ni,1:lst(id)%latPerPE,:) = Step0(1:lst(id)%ni,lst(id)%myLatBeg:lst(id)%myLatEnd,:)
       deallocate( Step0 )
    end if

    call lst_transform1d( Step1,                       & ! INOUT
                          TransformDirection,          & ! IN
                          TransformAxe,                & ! IN
                          ni_l, nj_l, nip_l, njp_l,    & ! IN
                          kStart, kEnd )                 ! IN

    !
    !- 2.0 Second pass (Step1 -> Step2)
    !   

    !- 2.1 Settings
    if ( trim(TransformDirection) == 'GridPointToSpectral' ) then
       TransformAxe = 'j'
       ni_l  = 2*lst(id)%mymCount
       nip_l = 0
       nj_l  = lst(id)%nj
       njp_l = njp
       if ( trim(lst(id)%MpiMode) == 'NoMpi') then
         kStart= 1
         kEnd  = nk
       else
         kStart= lst(id)%myLevBeg
         kEnd  = lst(id)%myLevEnd
       end if
    else
       TransformAxe = 'i'
       ni_l  = lst(id)%ni
       nip_l = nip
       nj_l  = lst(id)%latPerPE
       njp_l = 0
       if ( trim(lst(id)%MpiMode) == 'NoMpi') then
         kStart= 1
         kEnd  = nk
       else
         kStart= lst(id)%myLevBeg
         kEnd  = lst(id)%myLevEnd
       end if
    end if

    allocate( Step2(ni_l+nip_l,nj_l+njp_l,kStart:kEnd) )

    !- 2.2 Communication between processors
    
    if ( trim(TransformDirection) == 'GridPointToSpectral' ) then
       if      ( trim(lst(id)%MpiMode) == 'NoMpi') then
         Step2(:,1:lst(id)%nj,:) = Step1(:,1:lst(id)%nj,:)
       else
         call transpose2d_LatToM( Step2,    & ! OUT
                                  Step1, id ) ! IN
       end if
    else
       if      ( trim(lst(id)%MpiMode) == 'NoMpi') then
         Step2(:,1:lst(id)%nj,:) = Step1(:,1:lst(id)%nj,:)
       else
          call transpose2d_MToLat( Step2,    & ! OUT
                                   Step1, id ) ! IN
       end if
    end if

    deallocate( Step1 )

    !- 2.3 Spectral Transform
    call lst_transform1d( Step2,                      & ! INOUT
                          TransformDirection,         & ! IN
                          TransformAxe,               & ! IN
                          ni_l, nj_l, nip_l, njp_l,   & ! IN
                          kStart, kEnd )                ! IN

    !
    !- 3.0 Post-processing (Step2 -> Step3 -> Output)
    ! 

    select case ( trim(TransformDirection) )
    case ('GridPointToSpectral')
       iStart = 1 
       iEnd   = 2*lst(id)%mymCount
       jStart = 1
       jEnd   = 2*lst(id)%mynCount
       kStart= 1
       kEnd  = nk
    case ('SpectralToGridPoint')
       iStart = 1
       iEnd   = lst(id)%lonPerPE
       jStart = 1
       jEnd   = lst(id)%latPerPE
       kStart = 1
       kEnd   = nk
    case default
       write(*,*)
       write(*,*) 'Error: TranformDirection Unknown ', trim(TransformDirection)
       call utl_abort('lst_VarTransform')
    end select

    if ( trim(TransformDirection) == 'GridPointToSpectral' ) then

       ! Communication between processors (Truncation (if applicable) will occur in this step)
       if ( trim(lst(id)%MpiMode) == 'NoMpi' ) then
         call lst_ReshapeTrunc( Step2,             &  ! IN
                                SpectralStateVar,  &  ! OUT
                                'ToVAR', kStart, kEnd, id )     ! IN
       else
         call transpose2d_LevToN(SpectralStateVar, & ! OUT
                                 Step2, nk, id )     ! IN
       end if

    else

       allocate( Step3(iStart:iEnd,jStart:jEnd,kStart:kEnd) )

       ! Communication between processors
       if ( trim(lst(id)%MpiMode) == 'NoMpi') then
         Step3(:,:,:) = Step2(1:lst(id)%ni,:,:)
       else
         call transpose2d_LevToLon( Step3,                          & ! OUT
                                    Step2(1:lst(id)%ni,:,:), nk, id ) ! IN
       end if

       GridState(lst(id)%myLonBeg:lst(id)%myLonEnd,lst(id)%myLatBeg:lst(id)%myLatEnd,:) = Step3(1:lst(id)%lonPerPE,1:lst(id)%latPerPE,:)

       deallocate(Step3)

    end if

    deallocate(Step2)

  end subroutine lst_VarTransform_ijk

!--------------------------------------------------------------------------
! LST_VARTRANSFORM_KIJ
!--------------------------------------------------------------------------
  subroutine lst_VarTransform_kij( id, SpectralStateVar, GridState,        &
                                   TransformDirection, nk)
    implicit none

    integer,          intent(in)    :: id
                                     ! LST ID
    integer,          intent(in)    :: nk
                                     ! Grid point data dimensions
    character(len=*), intent(in)    :: TransformDirection
                                     ! SpectralToGridPoint or
                                     ! GridPointToSpectral
    real(8),          intent(inout) :: GridState(nk,lst(id)%myLonBeg:lst(id)%myLonEnd,lst(id)%myLatBeg:lst(id)%myLatEnd)
                                     ! 3D field in grid point space
    real(8),          intent(inout) :: SpectralStateVar(:,:,:)
                                     ! 3D spectral coefficients

    integer                         :: m, n, k, ni_l, nj_l, nip_l, njp_l
    integer                         :: iStart, iEnd, jStart, jEnd, kStart, kEnd

    integer                         :: i, j

    real(8), allocatable            :: Step0(:,:,:)
    real(8), allocatable            :: Step1(:,:,:)
    real(8), allocatable            :: Step2(:,:,:)
    real(8), allocatable            :: Step3(:,:,:)

    character(len=1)                :: TransformAxe

    !
    !- 1. Data pre-processing (Input -> Step0)
    !
    if (verbose) write(*,*) 'Entering lst_varTransform_kij'

    !- 1.1 Settings and Data Selection

    if ( trim(TransformDirection) == 'GridPointToSpectral' ) then
       iStart = 1 
       iEnd   = lst(id)%ni
       jStart = lst(id)%myLatBeg
       jEnd   = lst(id)%myLatEnd
       if ( trim(lst(id)%MpiMode) == 'NoMpi') then
         kStart= 1
         kEnd  = nk
       else
         kStart= lst(id)%myLevBeg
         kEnd  = lst(id)%myLevEnd
       end if
       allocate( Step0(kStart:kEnd,iStart:iEnd,jStart:jEnd) )
       if ( trim(lst(id)%MpiMode) == 'NoMpi') then
         Step0(:,:,:) = GridState(:,:,:)
       else
         call transpose2d_LonToLev_kij( Step0,            & ! OUT
                                        GridState, nk, id ) ! IN
       end if
    end if

    !
    !- 1.  First pass (Step0 -> Step1)
    !

    !- 1.1 Settings and Data Selection
    select case ( trim(TransformDirection) )
    case ('GridPointToSpectral')
       TransformAxe = 'i'
       ni_l  = lst(id)%ni
       nip_l = nip
       nj_l  = lst(id)%latPerPE
       njp_l = 0
       if ( trim(lst(id)%MpiMode) == 'NoMpi') then
         kStart= 1
         kEnd  = nk
       else
         kStart= lst(id)%myLevBeg
         kEnd  = lst(id)%myLevEnd
       end if
    case ('SpectralToGridPoint')
       TransformAxe = 'j'
       ni_l  = 2*lst(id)%mymCount
       nip_l = 0
       nj_l  = lst(id)%nj
       njp_l = njp
       if ( trim(lst(id)%MpiMode) == 'NoMpi') then
         kStart= 1
         kEnd  = nk
       else
         kStart= lst(id)%myLevBeg
         kEnd  = lst(id)%myLevEnd
       end if
    case default
       write(*,*)
       write(*,*) 'Error: TranformDirection Unknown ', trim(TransformDirection)
       call utl_abort('lst_VarTransform')
    end select

    allocate( Step1(kStart:kEnd,ni_l+nip_l,nj_l+njp_l) )

    !- 1.2 Spectral transform
    if ( trim(TransformDirection) == 'SpectralToGridPoint' ) then
       if ( trim(lst(id)%MpiMode) == 'NoMpi' ) then
         call lst_ReshapeTrunc_kij( Step1,                       & ! OUT
                                    SpectralStateVar,            & ! IN
                                    'ToRPN', kStart, kEnd, id )    ! IN
       else
         call transpose2d_NToLev_kij( Step1,                   & ! OUT
                                      SpectralStateVar, nk, id ) ! IN
       end if
    else
       Step1(:,1:lst(id)%ni,1:lst(id)%latPerPE) = Step0(:,1:lst(id)%ni,lst(id)%myLatBeg:lst(id)%myLatEnd)
       deallocate( Step0 )
    end if

    call lst_transform1d_kij( Step1,                       & ! INOUT
                              TransformDirection,          & ! IN
                              TransformAxe,                & ! IN
                              ni_l, nj_l, nip_l, njp_l,    & ! IN
                              kStart, kEnd )                 ! IN

    !
    !- 2.0 Second pass (Step1 -> Step2)
    !   

    !- 2.1 Settings
    if ( trim(TransformDirection) == 'GridPointToSpectral' ) then
       TransformAxe = 'j'
       ni_l  = 2*lst(id)%mymCount
       nip_l = 0
       nj_l  = lst(id)%nj
       njp_l = njp
       if ( trim(lst(id)%MpiMode) == 'NoMpi') then
         kStart= 1
         kEnd  = nk
       else
         kStart= lst(id)%myLevBeg
         kEnd  = lst(id)%myLevEnd
       end if
    else
       TransformAxe = 'i'
       ni_l  = lst(id)%ni
       nip_l = nip
       nj_l  = lst(id)%latPerPE
       njp_l = 0
       if ( trim(lst(id)%MpiMode) == 'NoMpi') then
         kStart= 1
         kEnd  = nk
       else
         kStart= lst(id)%myLevBeg
         kEnd  = lst(id)%myLevEnd
       end if
    end if

    allocate( Step2(kStart:kEnd,ni_l+nip_l,nj_l+njp_l) )

    !- 2.2 Communication between processors
    
    if ( trim(TransformDirection) == 'GridPointToSpectral' ) then
       if      ( trim(lst(id)%MpiMode) == 'NoMpi') then
         Step2(:,1:lst(id)%nj,:) = Step1(:,1:lst(id)%nj,:)
       else
         call transpose2d_LatToM_kij( Step2,    & ! OUT
                                      Step1, id ) ! IN
       end if
    else
       if      ( trim(lst(id)%MpiMode) == 'NoMpi') then
         Step2(:,1:lst(id)%nj,:) = Step1(:,1:lst(id)%nj,:)
       else
          call transpose2d_MToLat_kij( Step2,    & ! OUT
                                       Step1, id ) ! IN
       end if
    end if

    deallocate( Step1 )

    !- 2.3 Spectral Transform
    call lst_transform1d_kij( Step2,                      & ! INOUT
                              TransformDirection,         & ! IN
                              TransformAxe,               & ! IN
                              ni_l, nj_l, nip_l, njp_l,   & ! IN
                              kStart, kEnd )                ! IN

    !
    !- 3.0 Post-processing (Step2 -> Step3 -> Output)
    ! 

    select case ( trim(TransformDirection) )
    case ('GridPointToSpectral')
       iStart = 1 
       iEnd   = 2*lst(id)%mymCount
       jStart = 1
       jEnd   = 2*lst(id)%mynCount
       kStart= 1
       kEnd  = nk
    case ('SpectralToGridPoint')
       iStart = 1
       iEnd   = lst(id)%lonPerPE
       jStart = 1
       jEnd   = lst(id)%latPerPE
       kStart = 1
       kEnd   = nk
    case default
       write(*,*)
       write(*,*) 'Error: TranformDirection Unknown ', trim(TransformDirection)
       call utl_abort('lst_VarTransform')
    end select

    if ( trim(TransformDirection) == 'GridPointToSpectral' ) then

       ! Communication between processors (Truncation (if applicable) will occur in this step)
       if ( trim(lst(id)%MpiMode) == 'NoMpi' ) then
         call lst_ReshapeTrunc_kij( Step2,             &  ! IN
                                    SpectralStateVar,  &  ! OUT
                                    'ToVAR', kStart, kEnd, id )     ! IN
       else
         call transpose2d_LevToN_kij(SpectralStateVar, & ! OUT
                                     Step2, nk, id )     ! IN
       end if

    else

       allocate( Step3(kStart:kEnd,iStart:iEnd,jStart:jEnd) )

       ! Communication between processors
       if ( trim(lst(id)%MpiMode) == 'NoMpi') then
         Step3(:,:,:) = Step2(:,1:lst(id)%ni,:)
       else
         call transpose2d_LevToLon_kij( Step3,                          & ! OUT
                                        Step2(:,1:lst(id)%ni,:), nk, id ) ! IN
       end if

       GridState(:,lst(id)%myLonBeg:lst(id)%myLonEnd,lst(id)%myLatBeg:lst(id)%myLatEnd) = Step3(:,1:lst(id)%lonPerPE,1:lst(id)%latPerPE)

       deallocate(Step3)

    end if

    deallocate(Step2)

  end subroutine lst_VarTransform_kij

!--------------------------------------------------------------------------
! LST_TRANSFORM1D
!--------------------------------------------------------------------------
  subroutine lst_transform1d( Field3d,            &
                              TransformDirection, &
                              TransformAxe,       &
                              ni_l, nj_l, nip_l, njp_l, kStart, kEnd)
    implicit none

    integer,          intent(in)        :: ni_l, nj_l, kStart, kEnd
                                         ! Grid point data dimensions
    integer,          intent(in)        :: nip_l, njp_l
                                         ! Extra point in spectral space
    character(len=*), intent(in)        :: TransformDirection
                                         ! SpectralToGridPoint or
                                         ! GridPointToSpectral
    character(len=*), intent(in)        :: TransformAxe
                                         ! 'i' or 'j'
    real(8),          intent(inout)     :: Field3d(1:ni_l+nip_l,1:nj_l+njp_l,kStart:kEnd)  
                                         ! InOut 3D field

    integer :: nit, njt
    integer :: way, type
    integer :: maxsize
    integer :: axe, n, nlot, nfact, np, lot, nk

    !
    !- 1.  Set some options
    !
    if (verbose) write(*,*) 'Entering lst_transform1d'
    call tmg_start(24,'LST_FFT')

    !- 1.1 Transform Direction
    select case ( trim(TransformDirection) )
    case ('GridPointToSpectral')
       way = -1
    case ('SpectralToGridPoint')
       way = +1
    case default
       write(*,*)
       write(*,*) 'Error: TranformDirection Unknown ', trim(TransformDirection)
       call utl_abort('lst_VarTransform')
    end select

    nk = kEnd - kStart + 1

    !
    !- 2.  Do the transforms in one direction
    !
    select case ( trim(TransformAxe) )
    case ('i')
       !- 2.1 First pass  --> Along INDEX "I"
       axe = 0
       n   = ni_l
       np  = nip_l
       nlot= nj_l
    case ('j')
       !- 2.2 Second pass --> Along INDEX "J"
       axe = 1
       n   = nj_l
       np  = njp_l
       nlot= ni_l
    case default
       write(*,*)
       write(*,*) 'Error: TranformAxe Unknown ', trim(TransformAxe)
       call utl_abort('lst_VarTransform')
    end select

    !- 1.2 Fast or Slow Fourier Transform ?
    nfact = n
    call ngfft( nfact ) ! INOUT

    if (nfact == n ) then
       call setfft8( n ) ! IN
    else
       call utl_abort('lst_VarTransform: This module can only handle fast sin&cos FFT')
    end if

    select case ( trim(TransformAxe) )
    case ('i')
!$OMP PARALLEL DO PRIVATE(lot)
       do lot = 1, nlot
          call ffft8(Field3d(:,lot,:), 1, n+np, nk, way)
       end do
!$OMP END PARALLEL DO
    case ('j')
!$OMP PARALLEL DO PRIVATE(lot)
       do lot = 1, nlot
          call ffft8(Field3d(lot,:,:), 1, n+np, nk, way)
       end do
!$OMP END PARALLEL DO
    end select

    !*     subroutine ffft8( a, inc, jump, lot, isign )
    !*     a      is the array containing input & output data
    !*     inc    is the increment within each data 'vector'
    !*            (e.g. inc=1 for consecutively stored data)
    !*     jump   is the increment between the start of each data vector
    !*     lot    is the number of data vectors
    !*     isign  = +1 for transform from spectral to gridpoint
    !*            = -1 for transform from gridpoint to spectral

    call tmg_stop(24)

  end subroutine lst_transform1d

!--------------------------------------------------------------------------
! LST_TRANSFORM1D_KIJ
!--------------------------------------------------------------------------
  subroutine lst_transform1d_kij( Field3d,            &
                                  TransformDirection, &
                                  TransformAxe,       &
                                  ni_l, nj_l, nip_l, njp_l, kStart, kEnd)
    implicit none

    integer,          intent(in)        :: ni_l, nj_l, kStart, kEnd
                                         ! Grid point data dimensions
    integer,          intent(in)        :: nip_l, njp_l
                                         ! Extra point in spectral space
    character(len=*), intent(in)        :: TransformDirection
                                         ! SpectralToGridPoint or
                                         ! GridPointToSpectral
    character(len=*), intent(in)        :: TransformAxe
                                         ! 'i' or 'j'
    real(8),          intent(inout)     :: Field3d(kStart:kEnd,1:ni_l+nip_l,1:nj_l+njp_l)  
                                         ! InOut 3D field

    integer :: nit, njt
    integer :: way, type
    integer :: maxsize
    integer :: axe, n, nlot, nfact, np, lot, nk, k

    !
    !- 1.  Set some options
    !
    if (verbose) write(*,*) 'Entering lst_transform1d_kij'
    call tmg_start(24,'LST_FFT')

    !- 1.1 Transform Direction
    select case ( trim(TransformDirection) )
    case ('GridPointToSpectral')
       way = -1
    case ('SpectralToGridPoint')
       way = +1
    case default
       write(*,*)
       write(*,*) 'Error: TranformDirection Unknown ', trim(TransformDirection)
       call utl_abort('lst_VarTransform')
    end select

    nk = kEnd - kStart + 1

    !
    !- 2.  Do the transforms in one direction
    !
    select case ( trim(TransformAxe) )
    case ('i')
       !- 2.1 First pass  --> Along INDEX "I"
       axe = 0
       n   = ni_l
       np  = nip_l
       nlot= nj_l
    case ('j')
       !- 2.2 Second pass --> Along INDEX "J"
       axe = 1
       n   = nj_l
       np  = njp_l
       nlot= ni_l
    case default
       write(*,*)
       write(*,*) 'Error: TranformAxe Unknown ', trim(TransformAxe)
       call utl_abort('lst_VarTransform')
    end select

    !- 1.2 Fast or Slow Fourier Transform ?
    nfact = n
    call ngfft( nfact ) ! INOUT

    if (nfact == n ) then
       call setfft8( n ) ! IN
    else
       call utl_abort('lst_VarTransform: This module can only handle fast sin&cos FFT')
    end if

    select case ( trim(TransformAxe) )
    case ('i')
!$OMP PARALLEL DO PRIVATE(lot,k)
       do lot = 1, nlot
          call ffft8(Field3d(:,:,lot), nk, 1, nk, way)
       end do
!$OMP END PARALLEL DO
    case ('j')
!$OMP PARALLEL DO PRIVATE(lot,k)
       do lot = 1, nlot
          call ffft8(Field3d(:,lot,:), nk, 1, nk, way)
       end do
!$OMP END PARALLEL DO
    end select

    !*     subroutine ffft8( a, inc, jump, lot, isign )
    !*     a      is the array containing input & output data
    !*     inc    is the increment within each data 'vector'
    !*            (e.g. inc=1 for consecutively stored data)
    !*     jump   is the increment between the start of each data vector
    !*     lot    is the number of data vectors
    !*     isign  = +1 for transform from spectral to gridpoint
    !*            = -1 for transform from gridpoint to spectral

    call tmg_stop(24)

  end subroutine lst_transform1d_kij

!--------------------------------------------------------------------------
! LST_Transpose2d_LonToLev
!--------------------------------------------------------------------------
  SUBROUTINE transpose2d_LonToLev(gd_out, gd_in, nk, id)
    implicit none

    integer, intent(in) :: nk, id

    real(8), intent(in) :: gd_in(lst(id)%myLonBeg:lst(id)%myLonEnd, lst(id)%myLatBeg:lst(id)%myLatEnd, nk)
    real(8), intent(out):: gd_out(lst(id)%ni, lst(id)%myLatBeg:lst(id)%myLatEnd, lst(id)%myLevBeg:lst(id)%myLevEnd)

    real(8) :: gd_send(lst(id)%lonPerPE, lst(id)%latPerPE, lst(id)%maxLevCount, mpi_npex)
    real(8) :: gd_recv(lst(id)%lonPerPE, lst(id)%latPerPE, lst(id)%maxLevCount, mpi_npex)

    integer :: yourid, nsize, ierr, levIndex, levIndex2, latIndex, latIndex2, lonIndex, lonIndex2

    if (verbose) write(*,*) 'Entering transpose2d_LonToLev'
    call rpn_comm_barrier("GRID",ierr)

    call tmg_start(28,'TRANSP_2D_LEVtoLON')

!$OMP PARALLEL DO PRIVATE(yourid,levIndex,levIndex2)
    do yourid = 0, (mpi_npex-1)
       do levIndex = lst(id)%allLevBeg(yourid+1), lst(id)%allLevEnd(yourid+1)
          levIndex2 = levIndex-lst(id)%allLevBeg(yourid+1)+1
          gd_send(:,:,levIndex2,yourid+1) = gd_in(:,:,levIndex)
       end do
    end do
!$OMP END PARALLEL DO

    nsize = lst(id)%lonPerPE*lst(id)%maxLevCount*lst(id)%latPerPE
    if ( mpi_npex > 1 ) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","EW",ierr)
    else
      gd_recv(:,:,:,1) = gd_send(:,:,:,1)
    end if

!$OMP PARALLEL DO PRIVATE(yourid,levIndex,levIndex2)
    do levIndex = lst(id)%myLevBeg, lst(id)%myLevEnd
       levIndex2 = levIndex - lst(id)%myLevBeg + 1
       do yourid = 0, (mpi_npex-1)
          gd_out(lst(id)%allLonBeg(yourid+1):lst(id)%allLonEnd(yourid+1),:,levIndex) = gd_recv(:,:,levIndex2,yourid+1)
       end do
    end do
!$OMP END PARALLEL DO

    call tmg_stop(28)

  END SUBROUTINE transpose2d_LonToLev

!--------------------------------------------------------------------------
! LST_Transpose2d_LonToLev_kij
!--------------------------------------------------------------------------
  SUBROUTINE transpose2d_LonToLev_kij(gd_out, gd_in, nk, id)
    implicit none

    integer, intent(in) :: nk, id

    real(8), intent(in) :: gd_in (nk,lst(id)%myLonBeg:lst(id)%myLonEnd, lst(id)%myLatBeg:lst(id)%myLatEnd)
    real(8), intent(out):: gd_out(lst(id)%myLevBeg:lst(id)%myLevEnd,lst(id)%ni, lst(id)%myLatBeg:lst(id)%myLatEnd)

    !real(8) :: gd_send(lst(id)%maxLevCount,lst(id)%lonPerPE, lst(id)%latPerPE, mpi_npex)
    !real(8) :: gd_recv(lst(id)%maxLevCount,lst(id)%lonPerPE, lst(id)%latPerPE, mpi_npex)

    integer :: yourid, nsize, ierr, levIndex, levIndex2, latIndex, latIndex2, lonIndex, lonIndex2

    if (verbose) write(*,*) 'Entering transpose2d_LonToLev_kij'
    call rpn_comm_barrier("GRID",ierr)

    call tmg_start(28,'TRANSP_2D_LEVtoLON')

!!$!$OMP PARALLEL DO PRIVATE(yourid,latIndex,latIndex2,levIndex,levIndex2,lonIndex,lonIndex2)
!!$    do yourid = 0, (mpi_npex-1)
!!$       do latIndex = lst(id)%myLatBeg, lst(id)%myLatEnd
!!$          latIndex2 = latIndex - lst(id)%myLatBeg + 1
!!$          do lonIndex = lst(id)%myLonBeg, lst(id)%myLonEnd
!!$             lonIndex2 = lonIndex - lst(id)%myLonBeg + 1
!!$             do levIndex = lst(id)%allLevBeg(yourid+1), lst(id)%allLevEnd(yourid+1)
!!$                levIndex2 = levIndex-lst(id)%allLevBeg(yourid+1)+1
!!$                gd_send(levIndex2,lonIndex2,latIndex2,yourid+1) = gd_in(levIndex,lonIndex,latIndex)
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$!$OMP END PARALLEL DO

    nsize = lst(id)%lonPerPE*lst(id)%maxLevCount*lst(id)%latPerPE
    if ( mpi_npex > 1 ) then
      !call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
      !                       gd_recv,nsize,"mpi_double_precision","EW",ierr)
      call mpi_alltoall(gd_in,      1, lst(id)%sendType_LonToLev,  &
                        gd_out,     1, lst(id)%recvType_LonToLev, mpi_comm_EW, ierr)
    else
      !gd_recv(:,:,:,1) = gd_send(:,:,:,1)
       gd_out(:,:,:) = gd_in(:,:,:)
    end if

!!$!$OMP PARALLEL DO PRIVATE(yourid,levIndex,levIndex2,lonIndex,lonIndex2,latIndex,latIndex2)
!!$    do yourid = 0, (mpi_npex-1)
!!$       do latIndex = lst(id)%myLatBeg, lst(id)%myLatEnd
!!$          latIndex2 = latIndex - lst(id)%myLatBeg + 1
!!$          do lonIndex = lst(id)%allLonBeg(yourid+1), lst(id)%allLonEnd(yourid+1)
!!$             lonIndex2 = lonIndex - lst(id)%allLonBeg(yourid+1) + 1
!!$             do levIndex = lst(id)%myLevBeg, lst(id)%myLevEnd
!!$                levIndex2 = levIndex - lst(id)%myLevBeg + 1
!!$                gd_out(levIndex,lonIndex,latIndex) = gd_recv(levIndex2,lonIndex2,latIndex2,yourid+1)
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$!$OMP END PARALLEL DO

    call tmg_stop(28)

  END SUBROUTINE transpose2d_LonToLev_kij

!--------------------------------------------------------------------------
! LST_Transpose2d_LevToLon
!--------------------------------------------------------------------------
  SUBROUTINE transpose2d_LevToLon(gd_out,gd_in,nk,id)
    implicit none

    integer, intent(in) :: nk, id

    real(8), intent(out):: gd_out(lst(id)%myLonBeg:lst(id)%myLonEnd, lst(id)%myLatBeg:lst(id)%myLatEnd, nk)
    real(8), intent(in) :: gd_in(lst(id)%ni, lst(id)%myLatBeg:lst(id)%myLatEnd, lst(id)%myLevBeg:lst(id)%myLevEnd)

    real(8) :: gd_send(lst(id)%lonPerPE, lst(id)%latPerPE, lst(id)%maxLevCount, mpi_npex)
    real(8) :: gd_recv(lst(id)%lonPerPE, lst(id)%latPerPE, lst(id)%maxLevCount, mpi_npex)

    integer :: yourid, nsize, ierr, levIndex, levIndex2, latIndex, latIndex2, lonIndex, lonIndex2

    if (verbose) write(*,*) 'Entering transpose2d_LevToLon'
    call rpn_comm_barrier("GRID",ierr)

    call tmg_start(28,'TRANSP_2D_LEVtoLON')

!$OMP PARALLEL DO PRIVATE(yourid,levIndex,levIndex2)
    do levIndex = lst(id)%myLevBeg, lst(id)%myLevEnd
       levIndex2 = levIndex - lst(id)%myLevBeg + 1
       do yourid = 0, (mpi_npex-1)
          gd_send(:,:,levIndex2,yourid+1) = gd_in(lst(id)%allLonBeg(yourid+1):lst(id)%allLonEnd(yourid+1),:,levIndex)
        end do
    end do
!$OMP END PARALLEL DO

    nsize = lst(id)%lonPerPE*lst(id)%maxLevCount*lst(id)%latPerPE
    if ( mpi_npex > 1 ) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","EW",ierr)
    else
      gd_recv(:,:,:,1) = gd_send(:,:,:,1)
    end if

!$OMP PARALLEL DO PRIVATE(yourid,levIndex,levIndex2)
    do yourid = 0, (mpi_npex-1)
       do levIndex=lst(id)%allLevBeg(yourid+1),lst(id)%allLevEnd(yourid+1)
          levIndex2=levIndex-lst(id)%allLevBeg(yourid+1)+1
          gd_out(:,:,levIndex) = gd_recv(:,:,levIndex2,yourid+1)
       end do
    end do
!$OMP END PARALLEL DO

    call tmg_stop(28)

  END SUBROUTINE transpose2d_LevtoLon

!--------------------------------------------------------------------------
! LST_Transpose2d_LevToLon_kij
!--------------------------------------------------------------------------
  SUBROUTINE transpose2d_LevToLon_kij(gd_out,gd_in,nk,id)
    implicit none

    integer, intent(in) :: nk, id

    real(8), intent(out):: gd_out(nk,lst(id)%myLonBeg:lst(id)%myLonEnd, lst(id)%myLatBeg:lst(id)%myLatEnd)
    real(8), intent(in) :: gd_in(lst(id)%myLevBeg:lst(id)%myLevEnd,lst(id)%ni, lst(id)%myLatBeg:lst(id)%myLatEnd)

  !  real(8) :: gd_send(lst(id)%maxLevCount, lst(id)%lonPerPE, lst(id)%latPerPE, mpi_npex)
  !  real(8) :: gd_recv(lst(id)%maxLevCount, lst(id)%lonPerPE, lst(id)%latPerPE, mpi_npex)

    integer :: yourid, nsize, ierr, levIndex, levIndex2, latIndex, latIndex2, lonIndex, lonIndex2

    if (verbose) write(*,*) 'Entering transpose2d_LevToLon_kij'
    call rpn_comm_barrier("GRID",ierr)

    call tmg_start(28,'TRANSP_2D_LEVtoLON')

!!$!$OMP PARALLEL DO PRIVATE(yourid,levIndex,levIndex2,lonIndex,lonIndex2,latIndex,latIndex2)
!!$    do yourid = 0, (mpi_npex-1)
!!$       do latIndex = lst(id)%myLatBeg, lst(id)%myLatEnd
!!$          latIndex2 = latIndex - lst(id)%myLatBeg + 1
!!$          do lonIndex = lst(id)%allLonBeg(yourid+1), lst(id)%allLonEnd(yourid+1)
!!$             lonIndex2 = lonIndex - lst(id)%allLonBeg(yourid+1) + 1
!!$             do levIndex = lst(id)%myLevBeg, lst(id)%myLevEnd
!!$                levIndex2 = levIndex - lst(id)%myLevBeg + 1
!!$                gd_send(levIndex2,lonIndex2,latIndex2,yourid+1) = gd_in(levIndex,lonIndex,latIndex)
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$!$OMP END PARALLEL DO

    nsize = lst(id)%lonPerPE*lst(id)%maxLevCount*lst(id)%latPerPE
    if ( mpi_npex > 1 ) then
      !call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
      !                       gd_recv,nsize,"mpi_double_precision","EW",ierr)
       call mpi_alltoall(gd_in,      1, lst(id)%sendType_LevToLon,  &
                         gd_out,     1, lst(id)%recvType_LevToLon, mpi_comm_EW, ierr) 
    else
      !gd_recv(:,:,:,1) = gd_send(:,:,:,1)
      gd_out(:,:,:) = gd_in(:,:,:)
    end if

!!$!$OMP PARALLEL DO PRIVATE(yourid,latIndex,latIndex2,levIndex,levIndex2,lonIndex,lonIndex2)
!!$    do yourid = 0, (mpi_npex-1)
!!$       do latIndex = lst(id)%myLatBeg, lst(id)%myLatEnd
!!$          latIndex2 = latIndex - lst(id)%myLatBeg + 1
!!$          do lonIndex = lst(id)%myLonBeg, lst(id)%myLonEnd
!!$             lonIndex2 = lonIndex - lst(id)%myLonBeg + 1
!!$             do levIndex=lst(id)%allLevBeg(yourid+1),lst(id)%allLevEnd(yourid+1)
!!$                levIndex2=levIndex-lst(id)%allLevBeg(yourid+1)+1
!!$                gd_out(levIndex,lonIndex,latIndex) = gd_recv(levIndex2,lonIndex2,latIndex2,yourid+1)
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$!$OMP END PARALLEL DO

    call tmg_stop(28)

  END SUBROUTINE transpose2d_LevtoLon_kij

!--------------------------------------------------------------------------
! LST_Transpose2d_LatToM
!--------------------------------------------------------------------------
  SUBROUTINE transpose2d_LatToM(gd_out, gd_in, id)
    implicit none

    integer, intent(in)  :: id

    real(8), intent(out) :: gd_out(2*lst(id)%mymCount,lst(id)%nj+njp  ,lst(id)%myLevBeg:lst(id)%myLevEnd)
    real(8), intent(in)  :: gd_in (lst(id)%ni+nip    ,lst(id)%latPerPE,lst(id)%myLevBeg:lst(id)%myLevEnd)

    real(8) :: gd_recv(lst(id)%maxmActiveCount,2,lst(id)%latPerPE, lst(id)%maxLevCount, mpi_npey)
    real(8) :: gd_send(lst(id)%maxmActiveCount,2,lst(id)%latPerPE, lst(id)%maxLevCount, mpi_npey)

    integer :: yourid, mIndex, icount, nsize, ierr, levIndex, levIndex2, latIndex, latIndex2

    if (verbose) write(*,*) 'Entering transpose2d_LatToM'
    call rpn_comm_barrier("GRID",ierr)

    call tmg_start(27,'TRANSP_2D_MtoLAT')

!$OMP PARALLEL DO PRIVATE(yourid,latIndex,levIndex,levIndex2,icount,mIndex)
    do yourid = 0, (mpi_npey-1)
       do levIndex = lst(id)%myLevBeg, lst(id)%myLevEnd
          levIndex2 = levIndex - lst(id)%myLevBeg + 1 
          do latIndex = 1, lst(id)%latPerPE
             gd_send(:,:,latIndex,levIndex2,yourid+1) = 0.d0
             icount = 0
             do mIndex = lst(id)%allmBeg(yourid+1), lst(id)%allmEnd(yourid+1), lst(id)%allmSkip(yourid+1)
                if ( lst(id)%KfromMNglb(mIndex,0) /= -1 ) then
                   icount = icount + 1
                   gd_send(icount,1,latIndex,levIndex2,yourid+1) = gd_in(2*mIndex+1, latIndex, levIndex)
                   gd_send(icount,2,latIndex,levIndex2,yourid+1) = gd_in(2*mIndex+2, latIndex, levIndex)
                end if
             end do
          end do
       end do
    end do
!$OMP END PARALLEL DO

    nsize = lst(id)%maxmActiveCount*2*lst(id)%maxLevCount*lst(id)%latPerPE
    if ( mpi_npey > 1 ) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","NS",ierr)
    else
      gd_recv(:,:,:,:,1) = gd_send(:,:,:,:,1)
    end if

!$OMP PARALLEL DO PRIVATE(yourid,latIndex,latIndex2,levIndex,levIndex2,icount,mIndex)
    do yourid = 0, (mpi_npey-1)
       do levIndex = lst(id)%myLevBeg, lst(id)%myLevEnd
          levIndex2 = levIndex - lst(id)%myLevBeg + 1
          do latIndex = lst(id)%allLatBeg(yourid+1), lst(id)%allLatEnd(yourid+1)
             latIndex2 = latIndex - lst(id)%allLatBeg(yourid+1) + 1
             icount = 0
             do mIndex = lst(id)%mymBeg, lst(id)%mymEnd, lst(id)%mymSkip
                if ( lst(id)%KfromMNglb(mIndex,0) /= -1 ) then
                   icount = icount + 1
                   gd_out(2*lst(id)%mymIndex(mIndex)-1,latIndex,levIndex) = gd_recv(icount,1,latIndex2,levIndex2,yourid+1)
                   gd_out(2*lst(id)%mymIndex(mIndex)  ,latIndex,levIndex) = gd_recv(icount,2,latIndex2,levIndex2,yourid+1)
                else
                   gd_out(2*lst(id)%mymIndex(mIndex)-1,latIndex,levIndex) = 0.d0
                   gd_out(2*lst(id)%mymIndex(mIndex)  ,latIndex,levIndex) = 0.d0
                end if
             end do
          end do
       end do
    end do
!$OMP END PARALLEL DO

    call tmg_stop(27)

  END SUBROUTINE transpose2d_LatToM

!--------------------------------------------------------------------------
! LST_Transpose2d_LatToM_kij
!--------------------------------------------------------------------------
  SUBROUTINE transpose2d_LatToM_kij(gd_out, gd_in, id)
    implicit none

    integer, intent(in)  :: id

    real(8), intent(out) :: gd_out(lst(id)%myLevBeg:lst(id)%myLevEnd,2*lst(id)%mymCount,lst(id)%nj+njp  )
    real(8), intent(in)  :: gd_in (lst(id)%myLevBeg:lst(id)%myLevEnd,lst(id)%ni+nip    ,lst(id)%latPerPE)

    real(8) :: gd_recv(lst(id)%maxLevCount,lst(id)%maxmActiveCount,2,lst(id)%latPerPE, mpi_npey)
    real(8) :: gd_send(lst(id)%maxLevCount,lst(id)%maxmActiveCount,2,lst(id)%latPerPE, mpi_npey)

    integer :: yourid, mIndex, icount, nsize, ierr, levIndex, levIndex2, latIndex, latIndex2

    if (verbose) write(*,*) 'Entering transpose2d_LatToM_kij'
    call rpn_comm_barrier("GRID",ierr)

    call tmg_start(27,'TRANSP_2D_MtoLAT')

!$OMP PARALLEL DO PRIVATE(yourid,latIndex,levIndex,levIndex2,icount,mIndex)
    do yourid = 0, (mpi_npey-1)
       do latIndex = 1, lst(id)%latPerPE
          gd_send(:,:,:,latIndex,yourid+1) = 0.d0
          icount = 0
          do mIndex = lst(id)%allmBeg(yourid+1), lst(id)%allmEnd(yourid+1), lst(id)%allmSkip(yourid+1)
             if ( lst(id)%KfromMNglb(mIndex,0) /= -1 ) then
                icount = icount + 1
                do levIndex = lst(id)%myLevBeg, lst(id)%myLevEnd
                   levIndex2 = levIndex - lst(id)%myLevBeg + 1
                   gd_send(levIndex2,icount,1,latIndex,yourid+1) = gd_in(levIndex,2*mIndex+1, latIndex)
                   gd_send(levIndex2,icount,2,latIndex,yourid+1) = gd_in(levIndex,2*mIndex+2, latIndex)
                end do
             end if
          end do
       end do
    end do
!$OMP END PARALLEL DO

    nsize = lst(id)%maxmActiveCount*2*lst(id)%maxLevCount*lst(id)%latPerPE
    if ( mpi_npey > 1 ) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","NS",ierr)
    else
      gd_recv(:,:,:,:,1) = gd_send(:,:,:,:,1)
    end if

!$OMP PARALLEL DO PRIVATE(yourid,latIndex,latIndex2,levIndex,levIndex2,icount,mIndex)
    do yourid = 0, (mpi_npey-1)
       do latIndex = lst(id)%allLatBeg(yourid+1), lst(id)%allLatEnd(yourid+1)
          latIndex2 = latIndex - lst(id)%allLatBeg(yourid+1) + 1
          icount = 0
          do mIndex = lst(id)%mymBeg, lst(id)%mymEnd, lst(id)%mymSkip
             if ( lst(id)%KfromMNglb(mIndex,0) /= -1 ) then
                icount = icount + 1
                do levIndex = lst(id)%myLevBeg, lst(id)%myLevEnd
                   levIndex2 = levIndex - lst(id)%myLevBeg + 1
                   gd_out(levIndex,2*lst(id)%mymIndex(mIndex)-1,latIndex) = gd_recv(levIndex2,icount,1,latIndex2,yourid+1)
                   gd_out(levIndex,2*lst(id)%mymIndex(mIndex)  ,latIndex) = gd_recv(levIndex2,icount,2,latIndex2,yourid+1)
                end do
             else
                gd_out(:,2*lst(id)%mymIndex(mIndex)-1,latIndex) = 0.d0
                gd_out(:,2*lst(id)%mymIndex(mIndex)  ,latIndex) = 0.d0
             end if
          end do
       end do
    end do
!$OMP END PARALLEL DO

    call tmg_stop(27)

  END SUBROUTINE transpose2d_LatToM_kij

!--------------------------------------------------------------------------
! Transpose2d_MToLat
!--------------------------------------------------------------------------
  SUBROUTINE transpose2d_MtoLat(gd_out, gd_in, id)
    implicit none

    integer, intent(in)  :: id

    real(8), intent(in)   :: gd_in (2*lst(id)%mymCount,lst(id)%nj+njp  ,lst(id)%myLevBeg:lst(id)%myLevEnd)
    real(8), intent(out)  :: gd_out(lst(id)%ni+nip    ,lst(id)%latPerPE,lst(id)%myLevBeg:lst(id)%myLevEnd)

    real(8) :: gd_recv(lst(id)%maxmActiveCount,2,lst(id)%latPerPE,lst(id)%maxLevCount, mpi_npey)
    real(8) :: gd_send(lst(id)%maxmActiveCount,2,lst(id)%latPerPE,lst(id)%maxLevCount, mpi_npey)

    integer :: yourid, mIndex, icount, nsize, ierr, levIndex, levIndex2, latIndex, latIndex2

    if (verbose) write(*,*) 'Entering transpose2d_MToLat'
    call rpn_comm_barrier("GRID",ierr)

    call tmg_start(27,'TRANSP_2D_MtoLAT')

!$OMP PARALLEL DO PRIVATE(yourid,latIndex,latIndex2,levIndex,levIndex2,icount,mIndex)
    do yourid = 0, (mpi_npey-1)
       do levIndex = lst(id)%myLevBeg, lst(id)%myLevEnd
          levIndex2 = levIndex - lst(id)%myLevBeg + 1
          do latIndex = lst(id)%allLatBeg(yourid+1), lst(id)%allLatEnd(yourid+1)
             latIndex2 = latIndex - lst(id)%allLatBeg(yourid+1) + 1
             gd_send(:,:,latIndex2,levIndex2,yourid+1) = 0.d0
             icount = 0
             do mIndex = lst(id)%mymBeg, lst(id)%mymEnd, lst(id)%mymSkip
                if ( lst(id)%KfromMNglb(mIndex,0) /= -1 ) then
                   icount = icount+1
                   gd_send(icount,1,latIndex2,levIndex2,yourid+1) = gd_in(2*lst(id)%mymIndex(mIndex)-1,latIndex,levIndex)
                   gd_send(icount,2,latIndex2,levIndex2,yourid+1) = gd_in(2*lst(id)%mymIndex(mIndex)  ,latIndex,levIndex)
                end if
             end do
          end do
       end do
    end do
!$OMP END PARALLEL DO

    nsize = lst(id)%maxmActiveCount*2*lst(id)%maxLevCount*lst(id)%latPerPE
    if ( mpi_npey > 1 ) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","NS",ierr)
    else
      gd_recv(:,:,:,:,1) = gd_send(:,:,:,:,1)
    end if

!$OMP PARALLEL DO PRIVATE(yourid,latIndex,levIndex,levIndex2,icount,mIndex)
    do yourid = 0, (mpi_npey-1)
       do levIndex = lst(id)%myLevBeg, lst(id)%myLevEnd
          levIndex2 = levIndex - lst(id)%myLevBeg + 1
          do latIndex = 1, lst(id)%latPerPE
             icount = 0
             do mIndex = lst(id)%allmBeg(yourid+1), lst(id)%allmEnd(yourid+1), lst(id)%allmSkip(yourid+1)
                if ( lst(id)%KfromMNglb(mIndex,0) /= -1 ) then
                   icount = icount+1
                   gd_out(2*mIndex+1,latIndex,levIndex) = gd_recv(icount,1,latIndex,levIndex2,yourid+1)
                   gd_out(2*mIndex+2,latIndex,levIndex) = gd_recv(icount,2,latIndex,levIndex2,yourid+1)
                else
                   gd_out(2*mIndex+1,latIndex,levIndex) = 0.d0
                   gd_out(2*mIndex+2,latIndex,levIndex) = 0.d0
                end if
             end do
          end do
       end do
    end do
!$OMP END PARALLEL DO

    call tmg_stop(27)

  END SUBROUTINE transpose2d_MtoLat

!--------------------------------------------------------------------------
! Transpose2d_MToLat_kij
!--------------------------------------------------------------------------
  SUBROUTINE transpose2d_MtoLat_kij(gd_out, gd_in, id)
    implicit none

    integer, intent(in)  :: id

    real(8), intent(in)   :: gd_in (lst(id)%myLevBeg:lst(id)%myLevEnd,2*lst(id)%mymCount,lst(id)%nj+njp  )
    real(8), intent(out)  :: gd_out(lst(id)%myLevBeg:lst(id)%myLevEnd,lst(id)%ni+nip    ,lst(id)%latPerPE)

    real(8) :: gd_recv(lst(id)%maxLevCount,lst(id)%maxmActiveCount,2,lst(id)%latPerPE, mpi_npey)
    real(8) :: gd_send(lst(id)%maxLevCount,lst(id)%maxmActiveCount,2,lst(id)%latPerPE, mpi_npey)

    integer :: yourid, mIndex, icount, nsize, ierr, levIndex, levIndex2, latIndex, latIndex2

    if (verbose) write(*,*) 'Entering transpose2d_MToLat_kij'
    call rpn_comm_barrier("GRID",ierr)

    call tmg_start(27,'TRANSP_2D_MtoLAT')

!$OMP PARALLEL DO PRIVATE(yourid,latIndex,latIndex2,levIndex,levIndex2,icount,mIndex)
    do yourid = 0, (mpi_npey-1)
       do latIndex = lst(id)%allLatBeg(yourid+1), lst(id)%allLatEnd(yourid+1)
          latIndex2 = latIndex - lst(id)%allLatBeg(yourid+1) + 1
          gd_send(:,:,:,latIndex2,yourid+1) = 0.d0
          icount = 0
          do mIndex = lst(id)%mymBeg, lst(id)%mymEnd, lst(id)%mymSkip
             if ( lst(id)%KfromMNglb(mIndex,0) /= -1 ) then
                icount = icount+1
                do levIndex = lst(id)%myLevBeg, lst(id)%myLevEnd
                   levIndex2 = levIndex - lst(id)%myLevBeg + 1
                   gd_send(levIndex2,icount,1,latIndex2,yourid+1) = gd_in(levIndex,2*lst(id)%mymIndex(mIndex)-1,latIndex)
                   gd_send(levIndex2,icount,2,latIndex2,yourid+1) = gd_in(levIndex,2*lst(id)%mymIndex(mIndex)  ,latIndex)
                end do
             end if
          end do
       end do
    end do
!$OMP END PARALLEL DO

    nsize = lst(id)%maxmActiveCount*2*lst(id)%maxLevCount*lst(id)%latPerPE
    if ( mpi_npey > 1 ) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","NS",ierr)
    else
      gd_recv(:,:,:,:,1) = gd_send(:,:,:,:,1)
    end if

!$OMP PARALLEL DO PRIVATE(yourid,latIndex,levIndex,levIndex2,icount,mIndex)
    do yourid = 0, (mpi_npey-1)
       do latIndex = 1, lst(id)%latPerPE
          icount = 0
          do mIndex = lst(id)%allmBeg(yourid+1), lst(id)%allmEnd(yourid+1), lst(id)%allmSkip(yourid+1)
             if ( lst(id)%KfromMNglb(mIndex,0) /= -1 ) then
                icount = icount+1
                do levIndex = lst(id)%myLevBeg, lst(id)%myLevEnd
                   levIndex2 = levIndex - lst(id)%myLevBeg + 1
                   gd_out(levIndex,2*mIndex+1,latIndex) = gd_recv(levIndex2,icount,1,latIndex,yourid+1)
                   gd_out(levIndex,2*mIndex+2,latIndex) = gd_recv(levIndex2,icount,2,latIndex,yourid+1)
                end do
             else
                gd_out(:,2*mIndex+1,latIndex) = 0.d0
                gd_out(:,2*mIndex+2,latIndex) = 0.d0
             end if
          end do
       end do
    end do
!$OMP END PARALLEL DO

    call tmg_stop(27)

  END SUBROUTINE transpose2d_MtoLat_kij

!--------------------------------------------------------------------------
! LST_Transpose2d_LevToN
!--------------------------------------------------------------------------
  SUBROUTINE transpose2d_LevToN(SpectralStateVar, gd_in, nk, id)
    implicit none

    integer, intent(in) :: nk, id
    real(8), intent(out):: SpectralStateVar(lst(id)%nla,nphase,nk)
    real(8), intent(in) :: gd_in (2*lst(id)%mymCount, lst(id)%nj+njp, lst(id)%myLevBeg:lst(id)%myLevEnd)

    real(8) :: gd_send(lst(id)%maxnla, 4, lst(id)%maxLevCount, mpi_npex)
    real(8) :: gd_recv(lst(id)%maxnla, 4, lst(id)%maxLevCount, mpi_npex)

    integer :: yourid, nsize, ierr, levIndex, levIndex2, nIndex, mIndex, icount, ila

    if (verbose) write(*,*) 'Entering transpose2d_LevToN'
    call rpn_comm_barrier("GRID",ierr)

    call tmg_start(26,'TRANSP_2D_LEVtoN')

!$OMP PARALLEL DO PRIVATE(yourid,mIndex,levIndex,levIndex2,nIndex,icount)
    do yourid = 0, (mpi_npex-1)
       do levIndex = lst(id)%myLevBeg, lst(id)%myLevEnd
          levIndex2 = levIndex - lst(id)%myLevBeg + 1
          gd_send(:,:,levIndex2,yourid+1) = 0.d0
          icount = 0
          do nIndex = lst(id)%allnBeg(yourid+1), lst(id)%allnEnd(yourid+1), lst(id)%allnSkip(yourid+1)
             do mIndex = lst(id)%mymBeg, lst(id)%mymEnd, lst(id)%mymSkip 
                if ( lst(id)%KfromMNglb(mIndex,nIndex) /= -1 ) then
                   icount = icount + 1
                   gd_send(icount,1,levIndex2,yourid+1) = gd_in(2*lst(id)%mymIndex(mIndex)-1,2*nIndex+1,levIndex)
                   gd_send(icount,2,levIndex2,yourid+1) = gd_in(2*lst(id)%mymIndex(mIndex)-1,2*nIndex+2,levIndex)
                   gd_send(icount,3,levIndex2,yourid+1) = gd_in(2*lst(id)%mymIndex(mIndex)  ,2*nIndex+1,levIndex)
                   gd_send(icount,4,levIndex2,yourid+1) = gd_in(2*lst(id)%mymIndex(mIndex)  ,2*nIndex+2,levIndex)
                end if
             end do
          end do
       end do
    end do
!$OMP END PARALLEL DO

    nsize = lst(id)%maxnla*4*lst(id)%maxLevCount
    if ( mpi_npex > 1 ) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","EW",ierr)
    else
      gd_recv(:,:,:,1) = gd_send(:,:,:,1)
    end if

!$OMP PARALLEL DO PRIVATE(yourid,levIndex,levIndex2)
    do yourid = 0, (mpi_npex-1)
       do levIndex = lst(id)%allLevBeg(yourid+1), lst(id)%allLevEnd(yourid+1)
          levIndex2 = levIndex - lst(id)%allLevBeg(yourid+1) + 1
             SpectralStateVar(:,:,levIndex) = gd_recv(1:lst(id)%nla,:,levIndex2,yourid+1)
       end do
    end do
!$OMP END PARALLEL DO

    call tmg_stop(26)

  END SUBROUTINE transpose2d_LevToN

!--------------------------------------------------------------------------
! LST_Transpose2d_LevToN_kij
!--------------------------------------------------------------------------
  SUBROUTINE transpose2d_LevToN_kij(SpectralStateVar, gd_in, nk, id)
    implicit none

    integer, intent(in) :: nk, id
    real(8), intent(out):: SpectralStateVar(lst(id)%nla,nphase,nk)
    real(8), intent(in) :: gd_in (lst(id)%myLevBeg:lst(id)%myLevEnd,2*lst(id)%mymCount, lst(id)%nj+njp)

    real(8) :: gd_send(lst(id)%maxnla, 4, lst(id)%maxLevCount, mpi_npex)
    real(8) :: gd_recv(lst(id)%maxnla, 4, lst(id)%maxLevCount, mpi_npex)

    integer :: yourid, nsize, ierr, levIndex, levIndex2, nIndex, mIndex, icount, ila

    if (verbose) write(*,*) 'Entering transpose2d_LevToN_kij'
    call rpn_comm_barrier("GRID",ierr)

    call tmg_start(26,'TRANSP_2D_LEVtoN')

!$OMP PARALLEL DO PRIVATE(yourid,mIndex,levIndex,levIndex2,nIndex,icount)
    do yourid = 0, (mpi_npex-1)
       do levIndex = lst(id)%myLevBeg, lst(id)%myLevEnd
          levIndex2 = levIndex - lst(id)%myLevBeg + 1
          gd_send(:,:,levIndex2,yourid+1) = 0.d0
          icount = 0
          do nIndex = lst(id)%allnBeg(yourid+1), lst(id)%allnEnd(yourid+1), lst(id)%allnSkip(yourid+1)
             do mIndex = lst(id)%mymBeg, lst(id)%mymEnd, lst(id)%mymSkip 
                if ( lst(id)%KfromMNglb(mIndex,nIndex) /= -1 ) then
                   icount = icount + 1
                   gd_send(icount,1,levIndex2,yourid+1) = gd_in(levIndex,2*lst(id)%mymIndex(mIndex)-1,2*nIndex+1)
                   gd_send(icount,2,levIndex2,yourid+1) = gd_in(levIndex,2*lst(id)%mymIndex(mIndex)-1,2*nIndex+2)
                   gd_send(icount,3,levIndex2,yourid+1) = gd_in(levIndex,2*lst(id)%mymIndex(mIndex)  ,2*nIndex+1)
                   gd_send(icount,4,levIndex2,yourid+1) = gd_in(levIndex,2*lst(id)%mymIndex(mIndex)  ,2*nIndex+2)
                end if
             end do
          end do
       end do
    end do
!$OMP END PARALLEL DO

    nsize = lst(id)%maxnla*4*lst(id)%maxLevCount
    if ( mpi_npex > 1 ) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","EW",ierr)
    else
      gd_recv(:,:,:,1) = gd_send(:,:,:,1)
    end if

!$OMP PARALLEL DO PRIVATE(yourid,levIndex,levIndex2)
    do yourid = 0, (mpi_npex-1)
       do levIndex = lst(id)%allLevBeg(yourid+1), lst(id)%allLevEnd(yourid+1)
          levIndex2 = levIndex - lst(id)%allLevBeg(yourid+1) + 1
             SpectralStateVar(:,:,levIndex) = gd_recv(1:lst(id)%nla,:,levIndex2,yourid+1)
       end do
    end do
!$OMP END PARALLEL DO

    call tmg_stop(26)

  END SUBROUTINE transpose2d_LevToN_kij

!--------------------------------------------------------------------------
! LST_Transpose2d_NToLev
!--------------------------------------------------------------------------
  SUBROUTINE transpose2d_NToLev(gd_out, SpectralStateVar, nk, id)
    implicit none

    integer, intent(in) :: nk, id

    real(8), intent(in):: SpectralStateVar(lst(id)%nla,nphase,nk)
    real(8), intent(out):: gd_out(2*lst(id)%mymCount, lst(id)%nj+njp    , lst(id)%myLevBeg:lst(id)%myLevEnd)

    real(8) :: gd_send(lst(id)%maxnla, 4, lst(id)%maxLevCount, mpi_npex)
    real(8) :: gd_recv(lst(id)%maxnla, 4, lst(id)%maxLevCount, mpi_npex)

    integer :: yourid, nsize, ierr, levIndex, levIndex2, nIndex, mIndex, icount,ila

    if (verbose) write(*,*) 'Entering transpose2d_NToLev'
    call rpn_comm_barrier("GRID",ierr)

    call tmg_start(26,'TRANSP_2D_LEVtoN')

!$OMP PARALLEL DO PRIVATE(yourid,levIndex,levIndex2)
    do yourid = 0, (mpi_npex-1)
       do levIndex = lst(id)%allLevBeg(yourid+1), lst(id)%allLevEnd(yourid+1)
          levIndex2 = levIndex - lst(id)%allLevBeg(yourid+1) + 1
             gd_send(1:lst(id)%nla,:,levIndex2,yourid+1) = SpectralStateVar(:,:,levIndex)
             if (lst(id)%nla < lst(id)%maxnla) gd_send(lst(id)%nla+1:lst(id)%maxnla,:,levIndex2,yourid+1) = 0.d0
       end do
    end do
!$OMP END PARALLEL DO

    nsize = lst(id)%maxnla*4*lst(id)%maxLevCount
    if ( mpi_npex > 1 ) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","EW",ierr)
    else
      gd_recv(:,:,:,1) = gd_send(:,:,:,1)
    end if

!$OMP PARALLEL DO PRIVATE(yourid,nIndex,levIndex,levIndex2,mIndex,icount)
    do yourid = 0, (mpi_npex-1)
       do levIndex = lst(id)%myLevBeg, lst(id)%myLevEnd
          levIndex2 = levIndex - lst(id)%myLevBeg + 1
          icount = 0
          do nIndex = lst(id)%allnBeg(yourid+1), lst(id)%allnEnd(yourid+1), lst(id)%allnSkip(yourid+1)
             do mIndex = lst(id)%mymBeg, lst(id)%mymEnd, lst(id)%mymSkip
                if ( lst(id)%KfromMNglb(mIndex,nIndex) /= -1 ) then
                   icount = icount + 1
                   gd_out(2*lst(id)%mymIndex(mIndex)-1,2*nIndex+1,levIndex) = gd_recv(icount,1,levIndex2,yourid+1)
                   gd_out(2*lst(id)%mymIndex(mIndex)-1,2*nIndex+2,levIndex) = gd_recv(icount,2,levIndex2,yourid+1)
                   gd_out(2*lst(id)%mymIndex(mIndex)  ,2*nIndex+1,levIndex) = gd_recv(icount,3,levIndex2,yourid+1)
                   gd_out(2*lst(id)%mymIndex(mIndex)  ,2*nIndex+2,levIndex) = gd_recv(icount,4,levIndex2,yourid+1)
                else
                   gd_out(2*lst(id)%mymIndex(mIndex)-1,2*nIndex+1,levIndex) = 0.d0
                   gd_out(2*lst(id)%mymIndex(mIndex)-1,2*nIndex+2,levIndex) = 0.d0
                   gd_out(2*lst(id)%mymIndex(mIndex)  ,2*nIndex+1,levIndex) = 0.d0
                   gd_out(2*lst(id)%mymIndex(mIndex)  ,2*nIndex+2,levIndex) = 0.d0
                end if
             end do
          end do
       end do
    end do
!$OMP END PARALLEL DO

    call tmg_stop(26)

  END SUBROUTINE transpose2d_NToLev

!--------------------------------------------------------------------------
! LST_Transpose2d_NToLev_kij
!--------------------------------------------------------------------------
  SUBROUTINE transpose2d_NToLev_kij(gd_out, SpectralStateVar, nk, id)
    implicit none

    integer, intent(in) :: nk, id

    real(8), intent(in):: SpectralStateVar(lst(id)%nla,nphase,nk)
    real(8), intent(out):: gd_out(lst(id)%myLevBeg:lst(id)%myLevEnd,2*lst(id)%mymCount, lst(id)%nj+njp)

    real(8) :: gd_send(lst(id)%maxnla, 4, lst(id)%maxLevCount, mpi_npex)
    real(8) :: gd_recv(lst(id)%maxnla, 4, lst(id)%maxLevCount, mpi_npex)

    integer :: yourid, nsize, ierr, levIndex, levIndex2, nIndex, mIndex, icount,ila

    if (verbose) write(*,*) 'Entering transpose2d_NToLev_kij'
    call rpn_comm_barrier("GRID",ierr)

    call tmg_start(26,'TRANSP_2D_LEVtoN')

!$OMP PARALLEL DO PRIVATE(yourid,levIndex,levIndex2)
    do yourid = 0, (mpi_npex-1)
       do levIndex = lst(id)%allLevBeg(yourid+1), lst(id)%allLevEnd(yourid+1)
          levIndex2 = levIndex - lst(id)%allLevBeg(yourid+1) + 1
             gd_send(1:lst(id)%nla,:,levIndex2,yourid+1) = SpectralStateVar(:,:,levIndex)
             if (lst(id)%nla < lst(id)%maxnla) gd_send(lst(id)%nla+1:lst(id)%maxnla,:,levIndex2,yourid+1) = 0.d0
       end do
    end do
!$OMP END PARALLEL DO

    nsize = lst(id)%maxnla*4*lst(id)%maxLevCount
    if ( mpi_npex > 1 ) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","EW",ierr)
    else
      gd_recv(:,:,:,1) = gd_send(:,:,:,1)
    end if

!$OMP PARALLEL DO PRIVATE(yourid,nIndex,levIndex,levIndex2,mIndex,icount)
    do yourid = 0, (mpi_npex-1)
       do levIndex = lst(id)%myLevBeg, lst(id)%myLevEnd
          levIndex2 = levIndex - lst(id)%myLevBeg + 1
          icount = 0
          do nIndex = lst(id)%allnBeg(yourid+1), lst(id)%allnEnd(yourid+1), lst(id)%allnSkip(yourid+1)
             do mIndex = lst(id)%mymBeg, lst(id)%mymEnd, lst(id)%mymSkip
                if ( lst(id)%KfromMNglb(mIndex,nIndex) /= -1 ) then
                   icount = icount + 1
                   gd_out(levIndex,2*lst(id)%mymIndex(mIndex)-1,2*nIndex+1) = gd_recv(icount,1,levIndex2,yourid+1)
                   gd_out(levIndex,2*lst(id)%mymIndex(mIndex)-1,2*nIndex+2) = gd_recv(icount,2,levIndex2,yourid+1)
                   gd_out(levIndex,2*lst(id)%mymIndex(mIndex)  ,2*nIndex+1) = gd_recv(icount,3,levIndex2,yourid+1)
                   gd_out(levIndex,2*lst(id)%mymIndex(mIndex)  ,2*nIndex+2) = gd_recv(icount,4,levIndex2,yourid+1)
                else
                   gd_out(levIndex,2*lst(id)%mymIndex(mIndex)-1,2*nIndex+1) = 0.d0
                   gd_out(levIndex,2*lst(id)%mymIndex(mIndex)-1,2*nIndex+2) = 0.d0
                   gd_out(levIndex,2*lst(id)%mymIndex(mIndex)  ,2*nIndex+1) = 0.d0
                   gd_out(levIndex,2*lst(id)%mymIndex(mIndex)  ,2*nIndex+2) = 0.d0
                end if
             end do
          end do
       end do
    end do
!$OMP END PARALLEL DO

    call tmg_stop(26)

  END SUBROUTINE transpose2d_NToLev_kij

!--------------------------------------------------------------------------
! LST_ReshapeTrunc
!--------------------------------------------------------------------------
  subroutine lst_ReshapeTrunc( SpectralStateRpn, SpectralStateVar,  &
                               Direction, kStart, kEnd, id )
    implicit none

    integer,          intent(in)    :: id, kStart, kEnd
    character(len=*), intent(in)    :: Direction ! ToVAR or ToRPN
    real(8),          intent(inout) :: SpectralStateRpn(2*lst(id)%mymCount,2*lst(id)%mynCount,kStart:kEnd)
    real(8),          intent(inout) :: SpectralStateVar(lst(id)%nla     ,nphase,kStart:kEnd)

    integer k, m, n, ila

    if (verbose) write(*,*) 'Entering lst_ReshapeTrunc'

    select case ( trim(Direction) )
    case ('ToVAR')
      ! Truncation (if applicable) will be applied here
!$OMP PARALLEL DO PRIVATE (n,m,ila,k) 
      do n = lst(id)%mynBeg, lst(id)%mynEnd, lst(id)%mynSkip
        do m = lst(id)%mymBeg, lst(id)%mymEnd, lst(id)%mymSkip
          ila = lst(id)%nla_Index(m,n)
          if ( ila /= -1 ) then
            do k = kStart, kEnd 
              SpectralStateVar(ila,1,k) = SpectralStateRpn(2*lst(id)%mymIndex(m)-1,2*lst(id)%mynIndex(n)-1,k)
              SpectralStateVar(ila,2,k) = SpectralStateRpn(2*lst(id)%mymIndex(m)-1,2*lst(id)%mynIndex(n)  ,k)
              SpectralStateVar(ila,3,k) = SpectralStateRpn(2*lst(id)%mymIndex(m)  ,2*lst(id)%mynIndex(n)-1,k)
              SpectralStateVar(ila,4,k) = SpectralStateRpn(2*lst(id)%mymIndex(m)  ,2*lst(id)%mynIndex(n),  k)
            end do
          end if
        end do
      end do
!$OMP END PARALLEL DO

    case ('ToRPN')
      SpectralStateRpn(:,:,:) = 0.0d0
!$OMP PARALLEL DO PRIVATE (n,m,ila,k)
      do n = lst(id)%mynBeg, lst(id)%mynEnd, lst(id)%mynSkip
        do m = lst(id)%mymBeg, lst(id)%mymEnd, lst(id)%mymSkip
          ila = lst(id)%nla_Index(m,n)
          if ( ila /= -1 ) then
            do k = kStart, kEnd
              SpectralStateRpn(2*lst(id)%mymIndex(m)-1,2*lst(id)%mynIndex(n)-1,k) = SpectralStateVar(ila,1,k)
              SpectralStateRpn(2*lst(id)%mymIndex(m)-1,2*lst(id)%mynIndex(n)  ,k) = SpectralStateVar(ila,2,k)
              SpectralStateRpn(2*lst(id)%mymIndex(m)  ,2*lst(id)%mynIndex(n)-1,k) = SpectralStateVar(ila,3,k)
              SpectralStateRpn(2*lst(id)%mymIndex(m)  ,2*lst(id)%mynIndex(n)  ,k) = SpectralStateVar(ila,4,k)
            end do
          end if
        end do
      end do
!$OMP END PARALLEL DO

    case default
      write(6,*)
      write(6,*) 'lst_ReshapeTrunc: Unknown Direction', trim(Direction)
      call utl_abort('lst_ReshapeTrunc')
    end select

  end subroutine lst_ReshapeTrunc

!--------------------------------------------------------------------------
! LST_ReshapeTrunc_kij
!--------------------------------------------------------------------------
  subroutine lst_ReshapeTrunc_kij( SpectralStateRpn, SpectralStateVar,  &
                                   Direction, kStart, kEnd, id )
    implicit none

    integer,          intent(in)    :: id, kStart, kEnd
    character(len=*), intent(in)    :: Direction ! ToVAR or ToRPN
    real(8),          intent(inout) :: SpectralStateRpn(kStart:kEnd,2*lst(id)%mymCount,2*lst(id)%mynCount)
    real(8),          intent(inout) :: SpectralStateVar(lst(id)%nla     ,nphase,kStart:kEnd)

    integer k, m, n, ila

    if (verbose) write(*,*) 'Entering lst_ReshapeTrunc_kij'

    select case ( trim(Direction) )
    case ('ToVAR')
      ! Truncation (if applicable) will be applied here
!$OMP PARALLEL DO PRIVATE (n,m,ila,k) 
      do n = lst(id)%mynBeg, lst(id)%mynEnd, lst(id)%mynSkip
        do m = lst(id)%mymBeg, lst(id)%mymEnd, lst(id)%mymSkip
          ila = lst(id)%nla_Index(m,n)
          if ( ila /= -1 ) then
            do k = kStart, kEnd 
              SpectralStateVar(ila,1,k) = SpectralStateRpn(k,2*lst(id)%mymIndex(m)-1,2*lst(id)%mynIndex(n)-1)
              SpectralStateVar(ila,2,k) = SpectralStateRpn(k,2*lst(id)%mymIndex(m)-1,2*lst(id)%mynIndex(n)  )
              SpectralStateVar(ila,3,k) = SpectralStateRpn(k,2*lst(id)%mymIndex(m)  ,2*lst(id)%mynIndex(n)-1)
              SpectralStateVar(ila,4,k) = SpectralStateRpn(k,2*lst(id)%mymIndex(m)  ,2*lst(id)%mynIndex(n)  )
            end do
          end if
        end do
      end do
!$OMP END PARALLEL DO

    case ('ToRPN')
      SpectralStateRpn(:,:,:) = 0.0d0
!$OMP PARALLEL DO PRIVATE (n,m,ila,k)
      do n = lst(id)%mynBeg, lst(id)%mynEnd, lst(id)%mynSkip
        do m = lst(id)%mymBeg, lst(id)%mymEnd, lst(id)%mymSkip
          ila = lst(id)%nla_Index(m,n)
          if ( ila /= -1 ) then
            do k = kStart, kEnd
              SpectralStateRpn(k,2*lst(id)%mymIndex(m)-1,2*lst(id)%mynIndex(n)-1) = SpectralStateVar(ila,1,k)
              SpectralStateRpn(k,2*lst(id)%mymIndex(m)-1,2*lst(id)%mynIndex(n)  ) = SpectralStateVar(ila,2,k)
              SpectralStateRpn(k,2*lst(id)%mymIndex(m)  ,2*lst(id)%mynIndex(n)-1) = SpectralStateVar(ila,3,k)
              SpectralStateRpn(k,2*lst(id)%mymIndex(m)  ,2*lst(id)%mynIndex(n)  ) = SpectralStateVar(ila,4,k)
            end do
          end if
        end do
      end do
!$OMP END PARALLEL DO

    case default
      write(6,*)
      write(6,*) 'lst_ReshapeTrunc_kij: Unknown Direction', trim(Direction)
      call utl_abort('lst_ReshapeTrunc_kij')
    end select

  end subroutine lst_ReshapeTrunc_kij

!--------------------------------------------------------------------------
! LST_Laplacian
!--------------------------------------------------------------------------
  subroutine lst_Laplacian( id, GridState, Mode, nk)
    implicit none

    integer,          intent(in)    :: id
                                     ! LST ID
    integer,          intent(in)    :: nk
                                     ! Grid point data dimensions
    real(8),          intent(inout) :: GridState(lst(id)%myLonBeg:lst(id)%myLonEnd,lst(id)%myLatBeg:lst(id)%myLatEnd,nk)  
                                     ! 3D field in grid point space
    character(len=*), intent(in)    :: Mode
                                     ! Forward or Inverse

    real(8), allocatable            :: SpectralStateVar(:,:,:)
    real(8), allocatable            :: factor(:)

    integer :: k, ila, p

    character(len=24)   :: kind

    if (verbose) write(*,*) 'Entering lst_Laplacian'

    allocate( SpectralStateVar(lst(id)%nla,nphase,nk) )
    allocate( factor(lst(id)%nla) )

    call idcheck(id)

    !
    !- 1.  Set Mode-dependent factors
    !
    select case ( trim(Mode) )
    case ('Forward')
      factor(:) = lst(id)%lapxy(:)
    case ('Inverse')
      factor(:) = lst(id)%ilapxy(:)
    case default
      write(6,*)
      write(6,*) 'lst_Laplacian: Error: Mode Unknown ', trim(Mode)
      call utl_abort('lst_Laplacian')
    end select

    !
    !- 2. Grid Point Space -> Spectral Space
    !
    kind = 'GridPointToSpectral'
    call lst_VarTransform( id,                    & ! IN
                           SpectralStateVar,      & ! OUT
                           GridState,             & ! IN
                           kind, nk     )           ! IN    

    !
    !- 3. Laplacian (forward or inverse) Transform
    !
!$OMP PARALLEL DO PRIVATE (k,ila,p)
    do k = 1, nk
      do ila = 1, lst(id)%nla
        do p = 1, nphase
          SpectralStateVar(ila,p,k) = factor(ila) * SpectralStateVar(ila,p,k)
        end do
      end do
    end do
!$OMP END PARALLEL DO

    !
    !- 4. Spectral Space -> Grid Point Space
    !
    kind = 'SpectralToGridPoint'
    call lst_VarTransform( id,                    & ! IN
                           SpectralStateVar,      & ! IN
                           GridState,             & ! OUT
                           kind, nk     )           ! IN

    deallocate( SpectralStateVar )
    deallocate( factor )

  end subroutine lst_Laplacian

!--------------------------------------------------------------------------
!   IDCHECK
!--------------------------------------------------------------------------
  subroutine idcheck(id)
    implicit none

    integer, intent(in) :: id

    if ( .not. lst(id)%allocated) then
       write(*,*)
       write(*,*) "transform ID ", id
       call utl_abort('lst_IDCHECK: Unknown transform ID')
    end if

  end subroutine idcheck

!--------------------------------------------------------------------------
!   NGFFT
!--------------------------------------------------------------------------
  subroutine ngfft(n)
    implicit none

    integer, intent(inout) :: n ! le plus petit entier >= n qui factorise

    integer, parameter :: l = 3
    integer :: k(l) , m
    data m , k / 8 , 2 , 3 , 5 /

    integer :: i, j

    if ( n <= m ) n = m + 1
    n = n - 1
1   n = n + 1
    i = n
2   do j = 1, l
       if (mod(i,k(j)) == 0 ) go to 4
    end do
    go to 1
4   i = i/k(j)
    if ( i /= 1 ) go to 2

  end subroutine ngfft

end module LamSpectralTransform_mod
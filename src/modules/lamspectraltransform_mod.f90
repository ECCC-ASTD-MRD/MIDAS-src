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

module lamSpectralTransform_mod
  ! MODULE lamSpectralTransform_mod (prefix='lst' category='3. High-level transformations')
  ! 
  ! :Purpose: Bi-Fourier spectral transform for limited-area applications.
  !           Depends on ffft8 and setfft8 routines in ARMNLIB.
  !
  use mpi
  use mpi_mod
  use mpivar_mod
  use MathPhysConstants_mod
  use earthConstants_mod
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
     integer, allocatable :: nePerKglobal(:) ! Number of spectral element in each
                                             !  total wavenumber bands over ALL PROCESSORS
     integer, allocatable :: ilaFromEK(:,:)  ! ila index associated to each spectral element
                                             !  of total wavenumber band
     real(8), allocatable :: NormFactor(:,:)
     real(8), allocatable :: NormFactorAd(:,:)
     integer              :: nlaGlobal       ! First dimension of VAR global spectral array
     integer, allocatable :: ilaGlobal(:)    ! Position of the local vector element in the global
                                             !  vector
     integer              :: nphase          ! Second dimension of VAR spectral array
     integer              :: ni
     integer              :: nj
     integer              :: nip
     integer              :: njp
     integer              :: mmax, nmax
     integer              :: ktrunc
     integer              :: nla             ! First dimension of VAR spectral array
     integer              :: maxnla
     integer              :: mymBeg, mymEnd, mymSkip, mymCount, mymActiveCount,maxmActiveCount
     integer              :: mynBeg, mynEnd, mynSkip, mynCount
     integer, allocatable :: nla_Index(:,:)
     real(8), allocatable :: lapxy(:)
     real(8), allocatable :: ilapxy(:)
     integer, allocatable :: allmBeg(:), allmEnd(:), allmSkip(:)
     integer, allocatable :: mymIndex(:)
     integer, allocatable :: allnBeg(:), allnEnd(:), allnSkip(:)
     integer, allocatable :: mynIndex(:)
     logical              :: allocated = .false.
     integer              :: latPerPE, latPerPEmax, myLatBeg, myLatEnd
     integer              :: lonPerPE, lonPerPEmax, myLonBeg, myLonEnd
     integer              :: myLevBeg, myLevEnd, myLevCount, maxLevCount
     integer, allocatable :: allLatBeg(:), allLatEnd(:), allLatPerPE(:)
     integer, allocatable :: allLonBeg(:), allLonEnd(:), allLonPerPE(:)
     integer, allocatable :: allLevBeg(:), allLevEnd(:)
     integer, allocatable :: KfromMNglb(:,:)
     character(len=10)    :: MpiMode
     character(len=3)     :: gridDataOrder ! Ordering the gridded data: 'ijk' or 'kij'
     integer              :: sendType_LevToLon, recvType_LevToLon
     integer              :: sendType_LonToLev, recvType_LonToLev
     logical              :: lonLatDivisible
   end type struct_lst

  ! TransformType = 'SinCos'
  integer, parameter              :: nip = 2     ! Padding
  integer, parameter              :: njp = 2     ! Padding
  integer, parameter              :: nphase = 4  ! For Sin&Cos we have a) sin(m)sin(n)
                                                 !                     b) cos(m)cos(n)
                                                 !                     c) sin(m)cos(n)
                                                 !                     d) cos(m)sin(n)

  logical, parameter :: verbose = .false.

contains

  !--------------------------------------------------------------------------
  ! lst_setup
  !--------------------------------------------------------------------------
  subroutine lst_Setup(lst, ni_in, nj_in, dlon_in, ktrunc_in,    &
                       MpiMode, maxlevels_opt, gridDataOrder_opt)
    implicit none

    type(struct_lst)                :: lst
                                     ! Parameters available to the outside world

    integer,          intent(in)    :: ni_in, nj_in    
                                     ! Global Grid point data horizontal dimensions
    character(len=*), intent(in)    :: MpiMode
                                     ! MPI Strategy
    real(8),          intent(in)    :: dlon_in
                                     ! Grid Spacing in Radians
    integer,          intent(in)    :: ktrunc_in
                                     ! Spectral Truncation (global)
    integer, intent(in), optional   :: maxlevels_opt
                                     ! Number of levels; Only needed when MpiMode = LatLev
    character(len=*), intent(in), optional :: gridDataOrder_opt
                                     ! 'ijk' or 'kij'

    real(8), allocatable            :: Kr8fromMN(:,:)

    integer, allocatable            :: KfromMN(:,:)
    integer, allocatable            :: my_KfromMNglb(:,:)
    integer                         :: kref, mref, nref
    integer                         :: m, n, k, kMax, ila, nfact_lon, nfact_lat
    integer                         :: ier, ilaglb, i, j, p
    real(8)                         :: dlon, dx2, fac, ca, cp, cb, cq, r
    real(8)                         :: NormFactor1, NormFactor2, NormFactor3
    real(8)                         :: NormFactorAd1, NormFactorAd2, NormFactorAd3
    real(8)                         :: factor, factorAd
    character(len=60)               :: kreftype
    logical                         :: divisibleLon, divisibleLat
    integer(kind=MPI_ADDRESS_KIND)  :: lowerBound, extent
    integer :: realSize, sendType, recvType, ierr

    !
    !- 1.  Set variables needed by the LAM Spectral Transform in VAR
    !
    if (verbose) write(*,*) 'Entering lst_Setup'

    if (lst%allocated) then
       call utl_abort('lst_setup: this structure is already allocated')
    end if

    kreftype = 'MAX' ! hardwired

    lst%ni     = ni_in
    lst%nj     = nj_in
    lst%nphase = nphase
    lst%nip    = nip
    lst%njp    = njp

    !  1.1 Check grid dimensions and set padding for the RPN DFT routines

    ! We need to padd the input array such as ...                              
    ! O O O O O O O O O
    ! O O O O O O O O O
    ! X X X X X X X O O
    ! X X X X X X X O O
    ! X X X X X X X O O
    ! X X X X X X X O O

    if (mod(lst%ni,2) /= 0 .or. mod(lst%nj,2) /= 0) then
       write(*,*) ' The regular Sin & Cos Transform requires that', &
                  ' dimensions be EVEN. Fields MUST be periodic' , &
                  ' but the last colum and row MUST NOT BE a '   , &
                  ' repetition of the first colum and row. '
       call utl_abort('lst_setup')
    end if

    nfact_lon = lst%ni
    call ngfft(nfact_lon) ! INOUT
    nfact_lat = lst%nj
    call ngfft(nfact_lat) ! INOUT

    if (nfact_lon /= lst%ni .or. nfact_lat /= lst%nj) then
      if (nfact_lon /= lst%ni) then
        write(*,*) 'Error: A fast transform cannot be used in X'
        write(6,6130) lst%ni, nfact_lon
      end if
      if (nfact_lat /= lst%nj) then
        write(*,*) 'Error: A fast transform cannot be used in Y'
        write(6,6140) lst%nj, nfact_lat
      end if
      call utl_abort('lst_setup')
    end if

6130 FORMAT('N = ni = ', I4,' the nearest factorizable N = ',I4)
6140 FORMAT('N = nj = ', I4,' the nearest factorizable N = ',I4)

    !- 1.2 Maximum of integer wavenumbers in x and y directions
    lst%mmax = lst%ni/2
    lst%nmax = lst%nj/2

    write(*,'(A,f8.1)') ' lst_Setup: Your grid spacing (in km) = ', ec_ra*dlon_in/1000.0
    write(*,*) '           Max wavenumbers in x-axis = ', lst%mmax            
    write(*,*) '           Max wavenumbers in y-axis = ', lst%nmax

    !- 1.3 MPI Strategy

    lst%MpiMode = MpiMode
    select case (trim(lst%MpiMode))
    case ('NoMpi')
       !- 1.3.1 No MPI

       ! range of LONS handled by this processor (ALL) in GRIDPOINT SPACE
       lst%lonPerPE   = lst%ni
       lst%myLonBeg   = 1
       lst%myLonEnd   = lst%ni

       ! range of LATS handled by this processor (ALL) in GRIDPOINT SPACE
       lst%latPerPE   = lst%nj
       lst%myLatBeg   = 1
       lst%myLatEnd   = lst%nj

       ! range of M handled by this processor (ALL) in SPECTRAL SPACE
       lst%mymBeg     = 0
       lst%mymEnd     = lst%mmax
       lst%mymSkip    = 1
       lst%mymCount   = lst%mmax + 1

       ! range of N handled by this processor (ALL) in SPECTRAL SPACE
       lst%mynBeg     = 0
       lst%mynEnd     = lst%nmax
       lst%mynSkip    = 1
       lst%mynCount   = lst%nmax + 1

       ! set a dummy range of LEVELS handled by this processor
       lst%myLevBeg   = -1
       lst%myLevEnd   = -1
       lst%myLevCount = 0

    case ('LatLonMN')
       !- 1.3.2 MPI 2D: Distribution of lon/lat tiles (gridpoint space) and n/m (spectral space)

       ! range of LONS handled by this processor in GRIDPOINT SPACE
       call mpivar_setup_lonbands(lst%ni,                        & ! IN
                                  lst%lonPerPE, lst%lonPerPEmax, & ! OUT
                                  lst%myLonBeg, lst%myLonEnd,    & ! OUT
                                  divisible_opt=divisibleLon)      ! OUT

       ! range of LATS handled by this processor in GRIDPOINT SPACE
       call mpivar_setup_latbands(lst%nj,                        & ! IN
                                  lst%latPerPE, lst%latPerPEmax, & ! OUT
                                  lst%myLatBeg, lst%myLatEnd,    & ! OUT
                                  divisible_opt=divisibleLat)      ! OUT

       lst%lonLatDivisible = (divisibleLon .and. divisibleLat)
       if(mpi_myid == 0) write(*,*) 'lst_setup: lonLatDivisible = ', lst%lonLatDivisible

       ! range of M handled by this processor in SPECTRAL SPACE
       call mpivar_setup_m(lst%mmax,                                        & ! IN
                           lst%mymBeg, lst%mymEnd, lst%mymSkip, lst%mymCount) ! OUT

       ! range of N handled by this processor in SPECTRAL SPACE
       call mpivar_setup_n(lst%nmax,                                        & ! IN
                           lst%mynBeg, lst%mynEnd, lst%mynSkip, lst%mynCount) ! OUT

       ! range of LEVELS TEMPORARILY handled by this processor DURING THE SPECTRAL TRANSFORM
       if (.not.present(maxlevels_opt)) then
          call utl_abort('lst_setup: ERROR, number of levels must be specified with MpiMode LatLonMN')
       end if
       ! 2D MPI decomposition: split levels across npex
       call mpivar_setup_levels(maxlevels_opt,                          & ! IN
                                lst%myLevBeg,lst%myLevEnd,lst%myLevCount) ! OUT

    case default
       write(*,*)
       write(*,*) 'Error: MpiMode Unknown ', trim(MpiMode)
       call utl_abort('lst_setup')
    end select

    write(*,*)
    write(*,*) ' I am processor ', mpi_myid+1, ' on a total of ', mpi_nprocs
    write(*,*) '          mband info = ', lst%mymBeg, lst%mymEnd, lst%mymSkip, lst%mymCount
    write(*,*) '          nband info = ', lst%mynBeg, lst%mynEnd, lst%mynSkip, lst%mynCount
    write(*,*) '          level info = ', lst%myLevBeg, lst%myLevEnd, lst%myLevCount

    ! Set M index
    allocate(lst%mymIndex(lst%mymBeg:lst%mymEnd))
    lst%mymIndex(:)=0
    do m = lst%mymBeg, lst%mymEnd, lst%mymSkip
      if (m == lst%mymBeg) then
        lst%mymIndex(m) = 1
      else
        lst%mymIndex(m) = lst%mymIndex(m-lst%mymSkip) + 1
      end if
    end do
    
    ! Set N index
    allocate(lst%mynIndex(lst%mynBeg:lst%mynEnd))
    lst%mynIndex(:)=0
    do n = lst%mynBeg, lst%mynEnd, lst%mynSkip
      if (n == lst%mynBeg) then
        lst%mynIndex(n) = 1
      else
        lst%mynIndex(n) = lst%mynIndex(n-lst%mynSkip) + 1
      end if
    end do

    if (trim(lst%MpiMode) /= 'NoMpi') then

      ! Gathering with respect to Longitude
      allocate(lst%allLonBeg(mpi_npex))
      call rpn_comm_allgather(lst%myLonBeg ,1,"mpi_integer",       &
                              lst%allLonBeg,1,"mpi_integer","EW",ier)
      if (mpi_myid == 0) write(*,*) 'AllLonBeg =', lst%allLonBeg(:)

      allocate(lst%allLonEnd(mpi_npex))
      call rpn_comm_allgather(lst%myLonEnd ,1,"mpi_integer",       &
                              lst%allLonEnd,1,"mpi_integer","EW",ier)
      if (mpi_myid == 0) write(*,*) 'AllLonEnd =', lst%allLonEnd(:)

      allocate(lst%allLonPerPE(mpi_npex))
      call rpn_comm_allgather(lst%lonPerPE ,1,"mpi_integer",       &
                              lst%allLonPerPE,1,"mpi_integer","EW",ier)
      if (mpi_myid == 0) write(*,*) 'AllLonPerPE =', lst%allLonPerPE(:)

      ! Gathering with respect to Latitude
      allocate(lst%allLatBeg(mpi_npey))
      call rpn_comm_allgather(lst%myLatBeg ,1,"mpi_integer",       &
                              lst%allLatBeg,1,"mpi_integer","NS",ier)
      if (mpi_myid == 0) write(*,*) 'AllLatBeg =', lst%allLatBeg(:)

      allocate(lst%allLatEnd(mpi_npey))
      call rpn_comm_allgather(lst%myLatEnd ,1,"mpi_integer",       &
                              lst%allLatEnd,1,"mpi_integer","NS",ier)
      if (mpi_myid == 0) write(*,*) 'AllLatEnd =', lst%allLatEnd(:)

      allocate(lst%allLatPerPE(mpi_npey))
      call rpn_comm_allgather(lst%latPerPE ,1,"mpi_integer",       &
                              lst%allLatPerPE,1,"mpi_integer","NS",ier)
      if (mpi_myid == 0) write(*,*) 'AllLatPerPE =', lst%allLatPerPE(:)

      ! Gathering with respect to M
      allocate(lst%allmBeg(mpi_npey))
      call rpn_comm_allgather(lst%mymBeg ,1,"mpi_integer",       &
                              lst%allmBeg,1,"mpi_integer","NS",ier)
      if (mpi_myid == 0) write(*,*) 'AllmBeg =', lst%allmBeg(:)

      allocate(lst%allmEnd(mpi_npey))
      call rpn_comm_allgather(lst%mymEnd ,1,"mpi_integer",       &
                              lst%allmEnd,1,"mpi_integer","NS",ier)
      if (mpi_myid == 0) write(*,*) 'allmEnd =', lst%allmEnd(:)
    
      allocate(lst%allmSkip(mpi_npey))
      call rpn_comm_allgather(lst%mymSkip ,1,"mpi_integer",       &
                              lst%allmSkip,1,"mpi_integer","NS",ier)
      if (mpi_myid == 0) write(*,*) 'allmSkip = ', lst%allmSkip(:)

      allocate(lst%allnBeg(mpi_npex))
      call rpn_comm_allgather(lst%mynBeg ,1,"mpi_integer",       &
                              lst%allnBeg,1,"mpi_integer","EW",ier)
      if (mpi_myid == 0) write(*,*) 'AllnBeg =', lst%allnBeg(:)

      allocate(lst%allnEnd(mpi_npex))
      call rpn_comm_allgather(lst%mynEnd ,1,"mpi_integer",       &
                              lst%allnEnd,1,"mpi_integer","EW",ier)
      if (mpi_myid == 0) write(*,*) 'AllnEnd =', lst%allnEnd(:)
    
      allocate(lst%allnSkip(mpi_npex))
      call rpn_comm_allgather(lst%mynSkip ,1,"mpi_integer",       &
                              lst%allnSkip,1,"mpi_integer","EW",ier)
      if (mpi_myid == 0) write(*,*) 'AllnSkip = ', lst%allnSkip(:)

      ! Gathering with respect to levels
      call rpn_comm_allreduce(lst%myLevCount,lst%maxLevCount, &
                                1,"MPI_INTEGER","MPI_MAX","GRID",ier)
      if (mpi_myid == 0) write(*,*) 'MaxLevCount =',lst%maxLevCount

      allocate(lst%allLevBeg(mpi_npex))
      call rpn_comm_allgather(lst%myLevBeg ,1,"mpi_integer",       &
                              lst%allLevBeg,1,"mpi_integer","EW",ier)
      if (mpi_myid == 0) write(*,*) 'AllLevBeg =', lst%allLevBeg(:)

      allocate(lst%allLevEnd(mpi_npex))
      call rpn_comm_allgather(lst%myLevEnd ,1,"mpi_integer",       &
                              lst%allLevEnd,1,"mpi_integer","EW",ier)
      if (mpi_myid == 0) write(*,*) 'AllLevEnd =', lst%allLevEnd(:)

      ! Setup mpi derived types used in transposes (only used when grid is divisible)
      ! ... mpi_type_vector(count, blocklength, stride, ...)
      ! ... mpi_type_create_resized(oldtype, lowerbound, extent(in bytes), newtype, ierr)
   
      call mpi_type_size(MPI_REAL8, realSize, ierr)
      lowerBound = 0

      ! create the send type for LevToLon
      extent = lst%maxLevCount * lst%lonPerPE * realSize
      call mpi_type_vector(lst%latPerPE, lst%maxLevCount * lst%lonPerPE,  &
           lst%maxLevCount * lst%ni, MPI_REAL8, sendtype, ierr)
      call mpi_type_create_resized(sendtype, lowerBound , extent, lst%sendType_LevToLon, ierr);
      call mpi_type_commit(lst%sendType_LevToLon,ierr)

      ! create the receive type for LevToLon
      extent = lst%maxLevCount * realSize
      call mpi_type_vector(lst%lonPerPE * lst%latPerPE , lst%maxLevCount,  &
           maxlevels_opt, MPI_REAL8, recvtype, ierr);
      call mpi_type_create_resized(recvtype, lowerBound, extent, lst%recvType_LevToLon, ierr);
      call mpi_type_commit(lst%recvType_LevToLon, ierr)

      ! create the send type for LonToLev
      extent = lst%maxLevCount * realSize
      call mpi_type_vector(lst%lonPerPE * lst%latPerPE , lst%maxLevCount,  &
           maxlevels_opt, MPI_REAL8, sendtype, ierr);
      call mpi_type_create_resized(sendtype, lowerBound, extent, lst%sendType_LonToLev, ierr);
      call mpi_type_commit(lst%sendType_LonToLev, ierr)
      
      ! create the recv type for LonToLev
      extent = lst%maxLevCount * lst%lonPerPE * realSize
      call mpi_type_vector(lst%latPerPE, lst%maxLevCount * lst%lonPerPE,  &
           lst%maxLevCount * lst%ni, MPI_REAL8, recvtype, ierr)
      call mpi_type_create_resized(recvtype, lowerBound , extent, lst%recvType_LonToLev, ierr);
      call mpi_type_commit(lst%recvType_LonToLev,ierr)
      
     end if

    !- 1.4 Compute the Total Wavenumber associated with weach m,n pairs and
    !      the number of spectral element in the VAR array (nla) 
    !      FOR THE LOCAL PROCESSOR
    allocate(Kr8fromMN(0:lst%mmax,0:lst%nmax))
    Kr8FromMN(:,:) = -1.d0

    allocate(KfromMN(0:lst%mmax,0:lst%nmax))
    KFromMN(:,:) = -1

    ! Denis et al., MWR, 2002
    !mref = lst%mmax
    !nref = lst%nmax

    ! old var code (L. Fillion)
    mref = lst%ni-1
    nref = lst%nj-1

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

    kMax = nint(lst_totalWaveNumber(lst%mmax,lst%nmax,mref,nref,kref))

    if (ktrunc_in == -1) then ! no truncation case
      lst%ktrunc = kMax
    else
      if (ktrunc_in > 0) then
        lst%ktrunc = ktrunc_in
      else
        call utl_abort('lst_setup: invalid truncation')
      end if
    end if

    if      (lst%ktrunc >= kMax) then
      write(*,*)
      write(*,*) 'lst_Setup: Warning: Truncation is larger than kMax'
      write(*,*) '           NO TRUNCATION will be applied'
    else if (lst%ktrunc > min(lst%mmax,lst%nmax)) then
      write(*,*)
      if (lst%ktrunc > lst%mmax) then
         write(*,*) 'lst_Setup: Warning: Truncation is larger than mmax only'
         write(*,*) '           TRUNCATION will be applied only above nmax'
         write(*,'(A,f8.1)') '          i.e., for wavelenght (in km) in y-axis smaller than ',&
                     (lst%nj*ec_ra*dlon_in/1000.0)/lst%ktrunc
      else
         write(*,*) 'lst_Setup: Warning: Truncation is larger than nmax only'
         write(*,*) '           TRUNCATION will be applied only above mmax'
         write(*,'(A,f8.1)') '          i.e., for wavelenght (in km) in x-axis smaller than ',&
                     (lst%ni*ec_ra*dlon_in/1000.0)/lst%ktrunc
      end if
    else
      write(*,*)
      write(*,*) 'lst_Setup: TRUNCATION will be applied above k = ',lst%ktrunc
      write(*,'(A,f8.1)') '          i.e., for wavelenght (in km) in x-axis smaller than ',&
                     (lst%ni*ec_ra*dlon_in/1000.0)/lst%ktrunc
      write(*,'(A,f8.1)') '          i.e., for wavelenght (in km) in y-axis smaller than ',&
                     (lst%nj*ec_ra*dlon_in/1000.0)/lst%ktrunc
    end if

    ila = 0
    do n = lst%mynBeg, lst%mynEnd, lst%mynSkip
      do m = lst%mymBeg, lst%mymEnd, lst%mymSkip
         r = lst_totalWaveNumber(m,n,mref,nref,kref)
         k = nint(r) ! or ceiling(r) as in Denis et al. ?
         if (k <= lst%ktrunc) then
            ila = ila +1
            KFromMN(m,n) = k
            Kr8FromMN(m,n) = r
         end if
      end do
    end do

    if (ila == 0) then
       write(*,*)
       write(*,*) 'There are no spectral elements associated to this mpi task!'
    end if

    lst%nla = ila     ! Number of spectral element per phase in the VAR array
    if (trim(lst%MpiMode) /= 'NoMpi') then
       call rpn_comm_allreduce(lst%nla,lst%maxnla, &
            1,"MPI_INTEGER","MPI_MAX","GRID",ier)
       if (mpi_myid == 0) write(*,*) 'MaxNLA =',lst%maxnla
    end if

    allocate(lst%KfromMNglb(0:lst%mmax,0:lst%nmax))
    allocate(my_KfromMNglb(0:lst%mmax,0:lst%nmax))
    my_KfromMNglb = 0
    my_KFromMNglb(:,:) = KFromMN(:,:)
    if (trim(lst%MpiMode) /= 'NoMpi') then
      call rpn_comm_allreduce(my_KFromMNglb,lst%KFromMNglb, &
                                (lst%mmax+1)*(lst%nmax+1),"MPI_INTEGER","MPI_MAX","GRID",ier)
    end if
    deallocate(my_KfromMNglb) 

    lst%mymActiveCount=0
    do m = lst%mymBeg, lst%mymEnd, lst%mymSkip
      if (KfromMN(m,0) /= -1) lst%mymActiveCount = lst%mymActiveCount + 1
    end do
    if (trim(lst%MpiMode) /= 'NoMpi') then
       call rpn_comm_allreduce(lst%mymActiveCount,lst%maxmActiveCount, &
                               1,"MPI_INTEGER","MPI_MAX","GRID",ier)
       if (mpi_myid == 0) write(*,*) 'MaxmActiveCount =',lst%maxmActiveCount
    end if

    !- 1.5 VAR spectral element ordering &
    !      Total Wavenumbers and Weights associated with each spectral element
    !      FOR THE LOCAL PROCESSOR
    allocate(lst%nla_Index(0:lst%mmax,0:lst%nmax))

    allocate(lst%k_r8(1:lst%nla))
    allocate(lst%k(1:lst%nla))
    allocate(lst%m(1:lst%nla))
    allocate(lst%n(1:lst%nla))
    allocate(lst%Weight(1:lst%nla))
    allocate(lst%nePerK(0:lst%ktrunc))
    allocate(lst%nePerKglobal(0:lst%ktrunc))
    allocate(lst%ilaFromEK(1:lst%nla,0:lst%ktrunc))
    allocate(lst%ilaGlobal(1:lst%nla))

    lst%nla_Index(:,:) = -1
    lst%ilaFromEK(:,:) = -1
    lst%NEPerK(:)      =  0
    lst%NEPerKglobal(:)=  0

    ila    = 0
    ilaglb = 0
    do n = 0, lst%nmax
       do m = 0, lst%mmax
        k    = KfromMN(m,n)

        if (lst%KfromMNglb(m,n) /= -1) ilaglb = ilaglb + 1 ! Global Index

        if (k /= -1) then
          ila = ila+1

          ! Internal index
          lst%nla_Index(m,n) = ila

          ! Outgoing (public) variables
          lst%nePerK(k) = lst%nePerK(k) + 1
          lst%ilaFromEK(lst%nePerK(k),k) = ila
          lst%k_r8(ila) = Kr8fromMN(m,n)
          lst%k(ila) = k
          lst%m(ila) = m
          lst%n(ila) = n
          lst%ilaGlobal(ila) = ilaglb

          ! Spectral coefficient weight associated with this index
          if (m == 0 .and. n == 0) then
            lst%Weight(ila) = 1.0d0
          else if (m /= 0 .and. n /= 0) then
            lst%Weight(ila) = 4.0d0
          else
            lst%Weight(ila) = 2.0d0
          end if

       end if

      end do
    end do

    lst%nlaGlobal = ilaglb ! Number of spectral element per phase in the VAR mpi global array

    if (trim(lst%MpiMode) /= 'NoMpi') then
      call rpn_comm_allreduce(lst%nePerK, lst%nePerKglobal, lst%ktrunc+1, &
                              "mpi_integer", "mpi_sum", "GRID", ierr)
    end if

    deallocate(Kr8fromMN)
    deallocate(KfromMN)

    !- 1.7 Gridded data ordering (input/output)
    if (present(gridDataOrder_opt)) then
       lst%gridDataOrder = trim(gridDataOrder_opt)
    else
       lst%gridDataOrder = 'ijk' ! default value
    end if

    select case (trim(lst%gridDataOrder))
    case ('ijk')
       write(*,*) 'lst_setup: gridded data ordering = IJK' 
    case ('kij')
       write(*,*) 'lst_setup: gridded data ordering = KIJ'
    case default
       write(*,*)
       write(*,*) 'Error: gridDataOrder Unknown ', trim(gridDataOrder_opt)
       call utl_abort('lst_setup')
    end select

    !
    !- 2.  Set factors for parseval identity
    !
    allocate(lst%NormFactor  (lst%nla,lst%nphase))
    allocate(lst%NormFactorAd(lst%nla,lst%nphase))

    Normfactor1   = 1.0d0
    Normfactor2   = 0.5d0 * sqrt(2.0d0)
    Normfactor3   = 0.5d0
    NormfactorAd1 =      1.0d0  * real((lst%ni * lst%nj),8)
    NormfactorAd2 = sqrt(2.0d0) * real((lst%ni * lst%nj),8)
    NormfactorAd3 =      2.0d0  * real((lst%ni * lst%nj),8)

    do ila = 1,lst%nla

      m = lst%m(ila)
      n = lst%n(ila)

      do p = 1, lst%nphase
         if      (p == 1) then
            i = 2*m+1
            j = 2*n+1
         else if (p == 2) then
            i = 2*m+1
            j = 2*n+2
         else if (p == 3) then
            i = 2*m+2
            j = 2*n+1
         else if (p == 4) then
            i = 2*m+2
            j = 2*n+2
         else
            call utl_abort('lst_Setup: Error in NormFactor')
         end if

         if (i == 1 .or. j == 1) then  
            if (i == 1 .and. j == 1) then
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
         
         lst%NormFactor  (ila,p) = factor
         lst%NormFactorAd(ila,p) = factorAd
      end do
      
   end do

    !
    !- 3.  Set variables needed by Forward and Inverse Laplacian
    !
    allocate(lst%lapxy (lst%nla))
    allocate(lst%ilapxy(lst%nla))

    dlon = dlon_in
    dx2  = (ec_ra*dlon)**2
    fac  = 2.d0/dx2

    do ila = 1,lst%nla
      ca = 2.d0*MPC_PI_R8 * lst%m(ila)
      cp = cos(ca/lst%ni)
      cb = 2.d0*MPC_PI_R8 * lst%n(ila)
      cq = cos(cb/lst%nj)

      lst%lapxy(ila) = fac * (cp + cq - 2.d0)
      if (lst%lapxy(ila) /= 0.d0) then
         lst%ilapxy(ila) = 1.d0 / lst%lapxy(ila)
      else
         lst%ilapxy(ila) = 0.d0
      end if
    end do

    !
    !- 4. Finalized
    !
    lst%allocated = .true.

  end subroutine lst_Setup

  !--------------------------------------------------------------------------
  ! lst_totalWaveNumber
  !--------------------------------------------------------------------------
  function lst_totalWaveNumber(m,n,mref,nref,kref) result(r)
    integer :: m, n, mref, nref, kref
    real(8) :: a, b, r

    a = real(m,8)/real(mref,8)
    b = real(n,8)/real(nref,8)
    r = real(kref,8) * sqrt((a**2) + (b**2)) ! Ellipse Shape if nref /= mref

  end function lst_totalWaveNumber

  !--------------------------------------------------------------------------
  ! lst_VARTRANSFORM
  !--------------------------------------------------------------------------
  subroutine lst_VarTransform(lst, SpectralStateVar, GridState, &
                              TransformDirection, nk)
    implicit none

    type(struct_lst), intent(in)    :: lst

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

    !
    !- 1. Call the appropriate transform
    !
    if (trim(lst%gridDataOrder) == 'ijk') then
       call lst_VarTransform_ijk(lst, SpectralStateVar, GridState, &
                                  TransformDirection, nk)
    else
       call lst_VarTransform_kij(lst, SpectralStateVar, GridState, &
                                  TransformDirection, nk)
    end if

  end subroutine lst_VarTransform

  !--------------------------------------------------------------------------
  ! lst_VARTRANSFORM_IJK
  !--------------------------------------------------------------------------
  subroutine lst_VarTransform_ijk(lst, SpectralStateVar, GridState, &
                                  TransformDirection, nk)
    implicit none

    type(struct_lst), intent(in)    :: lst

    integer,          intent(in)    :: nk
                                     ! Grid point data dimensions
    character(len=*), intent(in)    :: TransformDirection
                                     ! SpectralToGridPoint or
                                     ! GridPointToSpectral
    real(8),          intent(inout) :: GridState(lst%myLonBeg:lst%myLonEnd,lst%myLatBeg:lst%myLatEnd,nk)
                                     ! 3D field in grid point space
    real(8),          intent(inout) :: SpectralStateVar(:,:,:)
                                     ! 3D spectral coefficients

    integer                         :: ni_l, nj_l, nip_l, njp_l
    integer                         :: iStart, iEnd, jStart, jEnd, kStart, kEnd

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

    if (trim(TransformDirection) == 'GridPointToSpectral') then
       iStart = 1 
       iEnd   = lst%ni
       jStart = lst%myLatBeg
       jEnd   = lst%myLatEnd
       if (trim(lst%MpiMode) == 'NoMpi') then
         kStart= 1
         kEnd  = nk
       else
         kStart= lst%myLevBeg
         kEnd  = lst%myLevEnd
       end if
       allocate(Step0(iStart:iEnd,jStart:jEnd,kStart:kEnd))
       if (trim(lst%MpiMode) == 'NoMpi') then
         Step0(:,:,:) = GridState(:,:,:)
       else
         call transpose2d_LonToLev(Step0,            & ! OUT
                                   GridState, nk, lst) ! IN
       end if
    end if

    !
    !- 1.  First pass (Step0 -> Step1)
    !

    !- 1.1 Settings and Data Selection
    select case (trim(TransformDirection))
    case ('GridPointToSpectral')
       TransformAxe = 'i'
       ni_l  = lst%ni
       nip_l = lst%nip
       nj_l  = lst%latPerPE
       njp_l = 0
       if (trim(lst%MpiMode) == 'NoMpi') then
         kStart= 1
         kEnd  = nk
       else
         kStart= lst%myLevBeg
         kEnd  = lst%myLevEnd
       end if
    case ('SpectralToGridPoint')
       TransformAxe = 'j'
       ni_l  = 2*lst%mymCount
       nip_l = 0
       nj_l  = lst%nj
       njp_l = lst%njp
       if (trim(lst%MpiMode) == 'NoMpi') then
         kStart= 1
         kEnd  = nk
       else
         kStart= lst%myLevBeg
         kEnd  = lst%myLevEnd
       end if
    case default
       write(*,*)
       write(*,*) 'Error: TranformDirection Unknown ', trim(TransformDirection)
       call utl_abort('lst_VarTransform')
    end select

    allocate(Step1(ni_l+nip_l,nj_l+njp_l,kStart:kEnd))

    !- 1.2 Spectral transform
    if (trim(TransformDirection) == 'SpectralToGridPoint') then
       if (trim(lst%MpiMode) == 'NoMpi') then
         call lst_ReshapeTrunc(Step1,                       & ! OUT
                               SpectralStateVar,            & ! IN
                               'ToRPN', kStart, kEnd, lst)    ! IN
       else
         call transpose2d_NToLev(Step1,                   & ! OUT
                                 SpectralStateVar, nk, lst) ! IN
       end if
    else
       Step1(1:lst%ni,1:lst%latPerPE,:) = Step0(1:lst%ni,lst%myLatBeg:lst%myLatEnd,:)
       deallocate(Step0)
    end if

    call lst_transform1d(Step1,                       & ! INOUT
                         TransformDirection,          & ! IN
                         TransformAxe,                & ! IN
                         ni_l, nj_l, nip_l, njp_l,    & ! IN
                         kStart, kEnd)                  ! IN

    !
    !- 2.0 Second pass (Step1 -> Step2)
    !   

    !- 2.1 Settings
    if (trim(TransformDirection) == 'GridPointToSpectral') then
       TransformAxe = 'j'
       ni_l  = 2*lst%mymCount
       nip_l = 0
       nj_l  = lst%nj
       njp_l = lst%njp
       if (trim(lst%MpiMode) == 'NoMpi') then
         kStart= 1
         kEnd  = nk
       else
         kStart= lst%myLevBeg
         kEnd  = lst%myLevEnd
       end if
    else
       TransformAxe = 'i'
       ni_l  = lst%ni
       nip_l = lst%nip
       nj_l  = lst%latPerPE
       njp_l = 0
       if (trim(lst%MpiMode) == 'NoMpi') then
         kStart= 1
         kEnd  = nk
       else
         kStart= lst%myLevBeg
         kEnd  = lst%myLevEnd
       end if
    end if

    allocate(Step2(ni_l+nip_l,nj_l+njp_l,kStart:kEnd))

    !- 2.2 Communication between processors
    
    if (trim(TransformDirection) == 'GridPointToSpectral') then
       if      (trim(lst%MpiMode) == 'NoMpi') then
         Step2(:,1:lst%nj,:) = Step1(:,1:lst%nj,:)
       else
         call transpose2d_LatToM(Step2,    & ! OUT
                                 Step1, lst) ! IN
       end if
    else
       if      (trim(lst%MpiMode) == 'NoMpi') then
         Step2(:,1:lst%nj,:) = Step1(:,1:lst%nj,:)
       else
          call transpose2d_MToLat(Step2,    & ! OUT
                                  Step1, lst) ! IN
       end if
    end if

    deallocate(Step1)

    !- 2.3 Spectral Transform
    call lst_transform1d(Step2,                      & ! INOUT
                         TransformDirection,         & ! IN
                         TransformAxe,               & ! IN
                         ni_l, nj_l, nip_l, njp_l,   & ! IN
                         kStart, kEnd)                 ! IN

    !
    !- 3.0 Post-processing (Step2 -> Step3 -> Output)
    ! 

    select case (trim(TransformDirection))
    case ('GridPointToSpectral')
       iStart = 1 
       iEnd   = 2*lst%mymCount
       jStart = 1
       jEnd   = 2*lst%mynCount
       kStart= 1
       kEnd  = nk
    case ('SpectralToGridPoint')
       iStart = 1
       iEnd   = lst%lonPerPE
       jStart = 1
       jEnd   = lst%latPerPE
       kStart = 1
       kEnd   = nk
    case default
       write(*,*)
       write(*,*) 'Error: TranformDirection Unknown ', trim(TransformDirection)
       call utl_abort('lst_VarTransform')
    end select

    if (trim(TransformDirection) == 'GridPointToSpectral') then

       ! Communication between processors (Truncation (if applicable) will occur in this step)
       if (trim(lst%MpiMode) == 'NoMpi') then
         call lst_ReshapeTrunc(Step2,                     &  ! IN
                               SpectralStateVar,          &  ! OUT
                               'ToVAR', kStart, kEnd, lst)   ! IN
       else
         call transpose2d_LevToN(SpectralStateVar, & ! OUT
                                 Step2, nk, lst)     ! IN
       end if

    else

       allocate(Step3(iStart:iEnd,jStart:jEnd,kStart:kEnd))

       ! Communication between processors
       if (trim(lst%MpiMode) == 'NoMpi') then
         Step3(:,:,:) = Step2(1:lst%ni,:,:)
       else
         call transpose2d_LevToLon(Step3,                      & ! OUT
                                   Step2(1:lst%ni,:,:), nk, lst) ! IN
       end if

       GridState(lst%myLonBeg:lst%myLonEnd,lst%myLatBeg:lst%myLatEnd,:) = Step3(1:lst%lonPerPE,1:lst%latPerPE,:)

       deallocate(Step3)

    end if

    deallocate(Step2)

  end subroutine lst_VarTransform_ijk

  !--------------------------------------------------------------------------
  ! lst_VARTRANSFORM_KIJ
  !--------------------------------------------------------------------------
  subroutine lst_VarTransform_kij(lst, SpectralStateVar, GridState, &
                                  TransformDirection, nk)
    implicit none

    type(struct_lst), intent(in)    :: lst

    integer,          intent(in)    :: nk
                                     ! Grid point data dimensions
    character(len=*), intent(in)    :: TransformDirection
                                     ! SpectralToGridPoint or
                                     ! GridPointToSpectral
    real(8),          intent(inout) :: GridState(nk,lst%myLonBeg:lst%myLonEnd,lst%myLatBeg:lst%myLatEnd)
                                     ! 3D field in grid point space
    real(8),          intent(inout) :: SpectralStateVar(:,:,:)
                                     ! 3D spectral coefficients

    integer                         :: ni_l, nj_l, nip_l, njp_l
    integer                         :: iStart, iEnd, jStart, jEnd, kStart, kEnd

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

    if (trim(TransformDirection) == 'GridPointToSpectral') then
       iStart = 1 
       iEnd   = lst%ni
       jStart = lst%myLatBeg
       jEnd   = lst%myLatEnd
       if (trim(lst%MpiMode) == 'NoMpi') then
         kStart= 1
         kEnd  = nk
       else
         kStart= lst%myLevBeg
         kEnd  = lst%myLevEnd
       end if
       allocate(Step0(kStart:kEnd,iStart:iEnd,jStart:jEnd))
       if (trim(lst%MpiMode) == 'NoMpi') then
         Step0(:,:,:) = GridState(:,:,:)
       else
         if(lst%lonLatDivisible) then
           call transpose2d_LonToLev_kij_mpitypes(Step0,            & ! OUT
                                                  GridState, nk, lst) ! IN
         else
           call transpose2d_LonToLev_kij(Step0,            & ! OUT
                                         GridState, nk, lst) ! IN
         end if
       end if
    end if

    !
    !- 1.  First pass (Step0 -> Step1)
    !

    !- 1.1 Settings and Data Selection
    select case (trim(TransformDirection))
    case ('GridPointToSpectral')
       TransformAxe = 'i'
       ni_l  = lst%ni
       nip_l = lst%nip
       nj_l  = lst%latPerPE
       njp_l = 0
       if (trim(lst%MpiMode) == 'NoMpi') then
         kStart= 1
         kEnd  = nk
       else
         kStart= lst%myLevBeg
         kEnd  = lst%myLevEnd
       end if
    case ('SpectralToGridPoint')
       TransformAxe = 'j'
       ni_l  = 2*lst%mymCount
       nip_l = 0
       nj_l  = lst%nj
       njp_l = lst%njp
       if (trim(lst%MpiMode) == 'NoMpi') then
         kStart= 1
         kEnd  = nk
       else
         kStart= lst%myLevBeg
         kEnd  = lst%myLevEnd
       end if
    case default
       write(*,*)
       write(*,*) 'Error: TranformDirection Unknown ', trim(TransformDirection)
       call utl_abort('lst_VarTransform')
    end select

    allocate(Step1(kStart:kEnd,ni_l+nip_l,nj_l+njp_l))

    !- 1.2 Spectral transform
    if (trim(TransformDirection) == 'SpectralToGridPoint') then
       if (trim(lst%MpiMode) == 'NoMpi') then
         call lst_ReshapeTrunc_kij(Step1,                       & ! OUT
                                   SpectralStateVar,            & ! IN
                                   'ToRPN', kStart, kEnd, lst)    ! IN
       else
         call transpose2d_NToLev_kij(Step1,                   & ! OUT
                                     SpectralStateVar, nk, lst) ! IN
       end if
    else
       Step1(:,1:lst%ni,1:lst%latPerPE) = Step0(:,1:lst%ni,lst%myLatBeg:lst%myLatEnd)
       deallocate(Step0)
    end if

    call lst_transform1d_kij(Step1,                        & ! INOUT
                              TransformDirection,          & ! IN
                              TransformAxe,                & ! IN
                              ni_l, nj_l, nip_l, njp_l,    & ! IN
                              kStart, kEnd)                  ! IN

    !
    !- 2.0 Second pass (Step1 -> Step2)
    !   

    !- 2.1 Settings
    if (trim(TransformDirection) == 'GridPointToSpectral') then
       TransformAxe = 'j'
       ni_l  = 2*lst%mymCount
       nip_l = 0
       nj_l  = lst%nj
       njp_l = lst%njp
       if (trim(lst%MpiMode) == 'NoMpi') then
         kStart= 1
         kEnd  = nk
       else
         kStart= lst%myLevBeg
         kEnd  = lst%myLevEnd
       end if
    else
       TransformAxe = 'i'
       ni_l  = lst%ni
       nip_l = lst%nip
       nj_l  = lst%latPerPE
       njp_l = 0
       if (trim(lst%MpiMode) == 'NoMpi') then
         kStart= 1
         kEnd  = nk
       else
         kStart= lst%myLevBeg
         kEnd  = lst%myLevEnd
       end if
    end if

    allocate(Step2(kStart:kEnd,ni_l+nip_l,nj_l+njp_l))

    !- 2.2 Communication between processors
    
    if (trim(TransformDirection) == 'GridPointToSpectral') then
       if      (trim(lst%MpiMode) == 'NoMpi') then
         Step2(:,1:lst%nj,:) = Step1(:,1:lst%nj,:)
       else
         call transpose2d_LatToM_kij(Step2,    & ! OUT
                                     Step1, lst) ! IN
       end if
    else
       if      (trim(lst%MpiMode) == 'NoMpi') then
         Step2(:,1:lst%nj,:) = Step1(:,1:lst%nj,:)
       else
          call transpose2d_MToLat_kij(Step2,    & ! OUT
                                      Step1, lst) ! IN
       end if
    end if

    deallocate(Step1)

    !- 2.3 Spectral Transform
    call lst_transform1d_kij(Step2,                      & ! INOUT
                             TransformDirection,         & ! IN
                             TransformAxe,               & ! IN
                             ni_l, nj_l, nip_l, njp_l,   & ! IN
                             kStart, kEnd)                 ! IN

    !
    !- 3.0 Post-processing (Step2 -> Step3 -> Output)
    ! 

    select case (trim(TransformDirection))
    case ('GridPointToSpectral')
       iStart = 1 
       iEnd   = 2*lst%mymCount
       jStart = 1
       jEnd   = 2*lst%mynCount
       kStart= 1
       kEnd  = nk
    case ('SpectralToGridPoint')
       iStart = 1
       iEnd   = lst%lonPerPE
       jStart = 1
       jEnd   = lst%latPerPE
       kStart = 1
       kEnd   = nk
    case default
       write(*,*)
       write(*,*) 'Error: TranformDirection Unknown ', trim(TransformDirection)
       call utl_abort('lst_VarTransform')
    end select

    if (trim(TransformDirection) == 'GridPointToSpectral') then

       ! Communication between processors (Truncation (if applicable) will occur in this step)
       if (trim(lst%MpiMode) == 'NoMpi') then
         call lst_ReshapeTrunc_kij(Step2,                    &  ! IN
                                   SpectralStateVar,         &  ! OUT
                                   'ToVAR', kStart, kEnd, lst)  ! IN
       else
         call transpose2d_LevToN_kij(SpectralStateVar, & ! OUT
                                     Step2, nk, lst)     ! IN
       end if

    else

       allocate(Step3(kStart:kEnd,iStart:iEnd,jStart:jEnd))

       ! Communication between processors
       if (trim(lst%MpiMode) == 'NoMpi') then
         Step3(:,:,:) = Step2(:,1:lst%ni,:)
       else
         if(lst%lonLatDivisible) then
           call transpose2d_LevToLon_kij_mpitypes(Step3,                      & ! OUT
                                                  Step2(:,1:lst%ni,:), nk, lst) ! IN
         else
           call transpose2d_LevToLon_kij(Step3,                      & ! OUT
                                         Step2(:,1:lst%ni,:), nk, lst) ! IN
         end if
       end if

       GridState(:,lst%myLonBeg:lst%myLonEnd,lst%myLatBeg:lst%myLatEnd) = Step3(:,1:lst%lonPerPE,1:lst%latPerPE)

       deallocate(Step3)

    end if

    deallocate(Step2)

  end subroutine lst_VarTransform_kij

  !--------------------------------------------------------------------------
  ! lst_TRANSFORM1D
  !--------------------------------------------------------------------------
  subroutine lst_transform1d(Field3d,            &
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

    integer :: way
    integer :: axe, n, nlot, nfact, np, lot, nk

    !
    !- 1.  Set some options
    !
    if (verbose) write(*,*) 'Entering lst_transform1d'
    call utl_tmg_start(150,'low-level--lst_fft')

    !- 1.1 Transform Direction
    select case (trim(TransformDirection))
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
    select case (trim(TransformAxe))
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
    call ngfft(nfact) ! INOUT

    if (nfact == n) then
       call setfft8(n) ! IN
    else
       call utl_abort('lst_VarTransform: This module can only handle fast sin&cos FFT')
    end if

    select case (trim(TransformAxe))
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

    call tmg_stop(150)

  end subroutine lst_transform1d

  !--------------------------------------------------------------------------
  ! lst_TRANSFORM1D_KIJ
  !--------------------------------------------------------------------------
  subroutine lst_transform1d_kij(Field3d,            &
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

    integer :: way
    integer :: axe, n, nlot, nfact, np, lot, nk, k

    !
    !- 1.  Set some options
    !
    if (verbose) write(*,*) 'Entering lst_transform1d_kij'
    call utl_tmg_start(151,'low-level--lst_fft')

    !- 1.1 Transform Direction
    select case (trim(TransformDirection))
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
    select case (trim(TransformAxe))
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
    call ngfft(nfact) ! INOUT

    if (nfact == n) then
       call setfft8(n) ! IN
    else
       call utl_abort('lst_VarTransform: This module can only handle fast sin&cos FFT')
    end if

    select case (trim(TransformAxe))
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

    call tmg_stop(151)

  end subroutine lst_transform1d_kij

  !--------------------------------------------------------------------------
  ! lst_Transpose2d_LonToLev
  !--------------------------------------------------------------------------
  subroutine transpose2d_LonToLev(gd_out, gd_in, nk, lst)
    implicit none

    type(struct_lst)    :: lst

    integer, intent(in) :: nk
    real(8), intent(in) :: gd_in(lst%myLonBeg:lst%myLonEnd, lst%myLatBeg:lst%myLatEnd, nk)
    real(8), intent(out):: gd_out(lst%ni, lst%myLatBeg:lst%myLatEnd, lst%myLevBeg:lst%myLevEnd)

    real(8) :: gd_send(lst%lonPerPEmax, lst%latPerPEmax, lst%maxLevCount, mpi_npex)
    real(8) :: gd_recv(lst%lonPerPEmax, lst%latPerPEmax, lst%maxLevCount, mpi_npex)
    integer :: yourid, nsize, ierr, levIndex, levIndex2

    if (verbose) write(*,*) 'Entering transpose2d_LonToLev'
    call rpn_comm_barrier("GRID",ierr)

    call utl_tmg_start(155,'low-level--lst_transpose_LEVtoLON')

    !$OMP PARALLEL DO PRIVATE(yourid,levIndex,levIndex2)
    do yourid = 0, (mpi_npex-1)
      gd_send(:,:,:,yourid+1) = 0.0d0
      do levIndex = lst%allLevBeg(yourid+1), lst%allLevEnd(yourid+1)
        levIndex2 = levIndex-lst%allLevBeg(yourid+1)+1
        gd_send(1:lst%lonPerPE,1:lst%latPerPE,levIndex2,yourid+1) = gd_in(:,:,levIndex)
      end do
    end do
    !$OMP END PARALLEL DO

    nsize = lst%lonPerPEmax * lst%maxLevCount * lst%latPerPEmax
    if (mpi_npex > 1) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","EW",ierr)
    else
      gd_recv(:,:,:,1) = gd_send(:,:,:,1)
    end if

    !$OMP PARALLEL DO PRIVATE(yourid,levIndex,levIndex2)
    do levIndex = lst%myLevBeg, lst%myLevEnd
      levIndex2 = levIndex - lst%myLevBeg + 1
      do yourid = 0, (mpi_npex-1)
        gd_out(lst%allLonBeg(yourid+1):lst%allLonEnd(yourid+1),:,levIndex) =  &
            gd_recv(1:lst%allLonPerPE(yourid+1),1:lst%latPerPE,levIndex2,yourid+1)
      end do
    end do
    !$OMP END PARALLEL DO

    call tmg_stop(155)

  end subroutine transpose2d_LonToLev

  !--------------------------------------------------------------------------
  ! lst_Transpose2d_LonToLev_kij_mpitypes
  !--------------------------------------------------------------------------
  subroutine transpose2d_LonToLev_kij_mpitypes(gd_out, gd_in, nk, lst)
    implicit none

    type(struct_lst)    :: lst

    integer, intent(in) :: nk
    real(8), intent(in) :: gd_in (nk,lst%myLonBeg:lst%myLonEnd, lst%myLatBeg:lst%myLatEnd)
    real(8), intent(out):: gd_out(lst%myLevBeg:lst%myLevEnd,lst%ni, lst%myLatBeg:lst%myLatEnd)

    integer :: nsize, ierr

    if (verbose) write(*,*) 'Entering transpose2d_LonToLev_kij'
    call rpn_comm_barrier("GRID",ierr)

    call utl_tmg_start(155,'low-level--lst_transpose_LEVtoLON')

    nsize = lst%lonPerPE * lst%maxLevCount * lst%latPerPE
    if (mpi_npex > 1) then
      call mpi_alltoall(gd_in,      1, lst%sendType_LonToLev,  &
                        gd_out,     1, lst%recvType_LonToLev, mpi_comm_EW, ierr)
    else
       gd_out(:,:,:) = gd_in(:,:,:)
    end if

    call tmg_stop(155)

  end subroutine transpose2d_LonToLev_kij_mpitypes

  !--------------------------------------------------------------------------
  ! lst_Transpose2d_LonToLev_kij
  !--------------------------------------------------------------------------
  subroutine transpose2d_LonToLev_kij(gd_out, gd_in, nk, lst)
    implicit none

    type(struct_lst)    :: lst

    integer, intent(in) :: nk
    real(8), intent(in) :: gd_in (nk,lst%myLonBeg:lst%myLonEnd, lst%myLatBeg:lst%myLatEnd)
    real(8), intent(out):: gd_out(lst%myLevBeg:lst%myLevEnd,lst%ni, lst%myLatBeg:lst%myLatEnd)

    real(8) :: gd_send(lst%maxLevCount,lst%lonPerPEmax, lst%latPerPEmax, mpi_npex)
    real(8) :: gd_recv(lst%maxLevCount,lst%lonPerPEmax, lst%latPerPEmax, mpi_npex)
    integer :: yourid, nsize, ierr, levIndex, levIndex2, latIndex, latIndex2, lonIndex, lonIndex2

    if (verbose) write(*,*) 'Entering transpose2d_LonToLev_kij'
    call rpn_comm_barrier("GRID",ierr)

    call utl_tmg_start(155,'low-level--lst_transpose_LEVtoLON')

    !$OMP PARALLEL DO PRIVATE(yourid,latIndex,latIndex2,levIndex,levIndex2,lonIndex,lonIndex2)
    do yourid = 0, (mpi_npex-1)
      gd_send(:,:,:,yourid+1) = 0.0d0
      do latIndex = lst%myLatBeg, lst%myLatEnd
        latIndex2 = latIndex - lst%myLatBeg + 1
        do lonIndex = lst%myLonBeg, lst%myLonEnd
          lonIndex2 = lonIndex - lst%myLonBeg + 1
          do levIndex = lst%allLevBeg(yourid+1), lst%allLevEnd(yourid+1)
            levIndex2 = levIndex-lst%allLevBeg(yourid+1)+1
            gd_send(levIndex2,lonIndex2,latIndex2,yourid+1) = gd_in(levIndex,lonIndex,latIndex)
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    nsize = lst%lonPerPEmax * lst%maxLevCount * lst%latPerPEmax
    if (mpi_npex > 1) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","EW",ierr)
    else
      gd_recv(:,:,:,1) = gd_send(:,:,:,1)
    end if

    !$OMP PARALLEL DO PRIVATE(yourid,levIndex,levIndex2,lonIndex,lonIndex2,latIndex,latIndex2)
    do yourid = 0, (mpi_npex-1)
      do latIndex = lst%myLatBeg, lst%myLatEnd
        latIndex2 = latIndex - lst%myLatBeg + 1
        do lonIndex = lst%allLonBeg(yourid+1), lst%allLonEnd(yourid+1)
          lonIndex2 = lonIndex - lst%allLonBeg(yourid+1) + 1
          do levIndex = lst%myLevBeg, lst%myLevEnd
            levIndex2 = levIndex - lst%myLevBeg + 1
            gd_out(levIndex,lonIndex,latIndex) = gd_recv(levIndex2,lonIndex2,latIndex2,yourid+1)
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call tmg_stop(155)

  end subroutine transpose2d_LonToLev_kij

  !--------------------------------------------------------------------------
  ! lst_Transpose2d_LevToLon
  !--------------------------------------------------------------------------
  subroutine transpose2d_LevToLon(gd_out,gd_in,nk,lst)
    implicit none

    type(struct_lst)    :: lst

    integer, intent(in) :: nk
    real(8), intent(out):: gd_out(lst%myLonBeg:lst%myLonEnd, lst%myLatBeg:lst%myLatEnd, nk)
    real(8), intent(in) :: gd_in(lst%ni, lst%myLatBeg:lst%myLatEnd, lst%myLevBeg:lst%myLevEnd)

    real(8) :: gd_send(lst%lonPerPEmax, lst%latPerPEmax, lst%maxLevCount, mpi_npex)
    real(8) :: gd_recv(lst%lonPerPEmax, lst%latPerPEmax, lst%maxLevCount, mpi_npex)
    integer :: yourid, nsize, ierr, levIndex, levIndex2

    if (verbose) write(*,*) 'Entering transpose2d_LevToLon'
    call rpn_comm_barrier("GRID",ierr)

    call utl_tmg_start(155,'low-level--lst_transpose_LEVtoLON')

    !$OMP PARALLEL DO PRIVATE(yourid,levIndex,levIndex2)
    do levIndex = lst%myLevBeg, lst%myLevEnd
      levIndex2 = levIndex - lst%myLevBeg + 1
      do yourid = 0, (mpi_npex-1)
        gd_send(:,:,levIndex2,yourid+1) = 0.0d0
        gd_send(1:lst%allLonPerPE(yourid+1),1:lst%latPerPE,levIndex2,yourid+1) =  &
             gd_in(lst%allLonBeg(yourid+1):lst%allLonEnd(yourid+1),:,levIndex)
      end do
    end do
    !$OMP END PARALLEL DO
    
    nsize = lst%lonPerPEmax * lst%maxLevCount * lst%latPerPEmax
    if (mpi_npex > 1) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","EW",ierr)
    else
      gd_recv(:,:,:,1) = gd_send(:,:,:,1)
    end if

    !$OMP PARALLEL DO PRIVATE(yourid,levIndex,levIndex2)
    do yourid = 0, (mpi_npex-1)
      do levIndex=lst%allLevBeg(yourid+1),lst%allLevEnd(yourid+1)
        levIndex2=levIndex-lst%allLevBeg(yourid+1)+1
        gd_out(:,:,levIndex) = gd_recv(1:lst%lonPerPE,1:lst%latPerPE,levIndex2,yourid+1)
      end do
    end do
    !$OMP END PARALLEL DO

    call tmg_stop(155)

  end subroutine transpose2d_LevtoLon

  !--------------------------------------------------------------------------
  ! lst_Transpose2d_LevToLon_kij_mpitypes
  !--------------------------------------------------------------------------
  subroutine transpose2d_LevToLon_kij_mpitypes(gd_out,gd_in,nk,lst)
    implicit none

    type(struct_lst)    :: lst

    integer, intent(in) :: nk
    real(8), intent(out):: gd_out(nk,lst%myLonBeg:lst%myLonEnd, lst%myLatBeg:lst%myLatEnd)
    real(8), intent(in) :: gd_in(lst%myLevBeg:lst%myLevEnd,lst%ni, lst%myLatBeg:lst%myLatEnd)

    integer :: nsize, ierr

    if (verbose) write(*,*) 'Entering transpose2d_LevToLon_kij'
    call rpn_comm_barrier("GRID",ierr)

    call utl_tmg_start(155,'low-level--lst_transpose_LEVtoLON')

    nsize = lst%lonPerPE*lst%maxLevCount*lst%latPerPE
    if (mpi_npex > 1) then
      call mpi_alltoall(gd_in,      1, lst%sendType_LevToLon,  &
                        gd_out,     1, lst%recvType_LevToLon, mpi_comm_EW, ierr) 
    else
      gd_out(:,:,:) = gd_in(:,:,:)
    end if

    call tmg_stop(155)

  end subroutine transpose2d_LevtoLon_kij_mpitypes

  !--------------------------------------------------------------------------
  ! lst_Transpose2d_LevToLon_kij
  !--------------------------------------------------------------------------
  subroutine transpose2d_LevToLon_kij(gd_out,gd_in,nk,lst)
    implicit none

    type(struct_lst)    :: lst

    integer, intent(in) :: nk
    real(8), intent(out):: gd_out(nk,lst%myLonBeg:lst%myLonEnd, lst%myLatBeg:lst%myLatEnd)
    real(8), intent(in) :: gd_in(lst%myLevBeg:lst%myLevEnd,lst%ni, lst%myLatBeg:lst%myLatEnd)

    real(8) :: gd_send(lst%maxLevCount, lst%lonPerPEmax, lst%latPerPEmax, mpi_npex)
    real(8) :: gd_recv(lst%maxLevCount, lst%lonPerPEmax, lst%latPerPEmax, mpi_npex)
    integer :: yourid, nsize, ierr, levIndex, levIndex2, latIndex, latIndex2, lonIndex, lonIndex2

    if (verbose) write(*,*) 'Entering transpose2d_LevToLon_kij'
    call rpn_comm_barrier("GRID",ierr)

    call utl_tmg_start(155,'low-level--lst_transpose_LEVtoLON')
    
    !$OMP PARALLEL DO PRIVATE(yourid,levIndex,levIndex2,lonIndex,lonIndex2,latIndex,latIndex2)
    do yourid = 0, (mpi_npex-1)
      gd_send(:,:,:,yourid+1) = 0.0d0
      do latIndex = lst%myLatBeg, lst%myLatEnd
        latIndex2 = latIndex - lst%myLatBeg + 1
        do lonIndex = lst%allLonBeg(yourid+1), lst%allLonEnd(yourid+1)
          lonIndex2 = lonIndex - lst%allLonBeg(yourid+1) + 1
          do levIndex = lst%myLevBeg, lst%myLevEnd
            levIndex2 = levIndex - lst%myLevBeg + 1
            gd_send(levIndex2,lonIndex2,latIndex2,yourid+1) = gd_in(levIndex,lonIndex,latIndex)
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    nsize = lst%lonPerPEmax * lst%maxLevCount * lst%latPerPEmax
    if (mpi_npex > 1) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","EW",ierr)
    else
      gd_recv(:,:,:,1) = gd_send(:,:,:,1)
    end if

    !$OMP PARALLEL DO PRIVATE(yourid,latIndex,latIndex2,levIndex,levIndex2,lonIndex,lonIndex2)
    do yourid = 0, (mpi_npex-1)
      do latIndex = lst%myLatBeg, lst%myLatEnd
        latIndex2 = latIndex - lst%myLatBeg + 1
        do lonIndex = lst%myLonBeg, lst%myLonEnd
          lonIndex2 = lonIndex - lst%myLonBeg + 1
          do levIndex=lst%allLevBeg(yourid+1),lst%allLevEnd(yourid+1)
            levIndex2=levIndex-lst%allLevBeg(yourid+1)+1
            gd_out(levIndex,lonIndex,latIndex) = gd_recv(levIndex2,lonIndex2,latIndex2,yourid+1)
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call tmg_stop(155)

  end subroutine transpose2d_LevtoLon_kij

  !--------------------------------------------------------------------------
  ! lst_Transpose2d_LatToM
  !--------------------------------------------------------------------------
  subroutine transpose2d_LatToM(gd_out, gd_in, lst)
    implicit none

    type(struct_lst)     :: lst
    
    real(8), intent(out) :: gd_out(2*lst%mymCount,lst%nj+lst%njp,lst%myLevBeg:lst%myLevEnd)
    real(8), intent(in)  :: gd_in (lst%ni+lst%nip,lst%latPerPE,lst%myLevBeg:lst%myLevEnd)

    real(8) :: gd_recv(lst%maxmActiveCount,2,lst%latPerPEmax, lst%maxLevCount, mpi_npey)
    real(8) :: gd_send(lst%maxmActiveCount,2,lst%latPerPEmax, lst%maxLevCount, mpi_npey)
    integer :: yourid, mIndex, icount, nsize, ierr, levIndex, levIndex2, latIndex, latIndex2

    if (verbose) write(*,*) 'Entering transpose2d_LatToM'
    call rpn_comm_barrier("GRID",ierr)

    call utl_tmg_start(154,'low-level--lst_transpose_MtoLAT')

    !$OMP PARALLEL DO PRIVATE(yourid,latIndex,levIndex,levIndex2,icount,mIndex)
    do yourid = 0, (mpi_npey-1)
      do levIndex = lst%myLevBeg, lst%myLevEnd
        levIndex2 = levIndex - lst%myLevBeg + 1 
        gd_send(:,:,:,levIndex2,yourid+1) = 0.d0
        do latIndex = 1, lst%latPerPE
          icount = 0
          do mIndex = lst%allmBeg(yourid+1), lst%allmEnd(yourid+1), lst%allmSkip(yourid+1)
            if (lst%KfromMNglb(mIndex,0) /= -1) then
              icount = icount + 1
              gd_send(icount,1,latIndex,levIndex2,yourid+1) = gd_in(2*mIndex+1, latIndex, levIndex)
              gd_send(icount,2,latIndex,levIndex2,yourid+1) = gd_in(2*mIndex+2, latIndex, levIndex)
            end if
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    nsize = lst%maxmActiveCount * 2 * lst%maxLevCount * lst%latPerPEmax
    if (mpi_npey > 1) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","NS",ierr)
    else
      gd_recv(:,:,:,:,1) = gd_send(:,:,:,:,1)
    end if

    !$OMP PARALLEL DO PRIVATE(yourid,latIndex,latIndex2,levIndex,levIndex2,icount,mIndex)
    do yourid = 0, (mpi_npey-1)
      do levIndex = lst%myLevBeg, lst%myLevEnd
        levIndex2 = levIndex - lst%myLevBeg + 1
        do latIndex = lst%allLatBeg(yourid+1), lst%allLatEnd(yourid+1)
          latIndex2 = latIndex - lst%allLatBeg(yourid+1) + 1
          icount = 0
          do mIndex = lst%mymBeg, lst%mymEnd, lst%mymSkip
            if (lst%KfromMNglb(mIndex,0) /= -1) then
              icount = icount + 1
              gd_out(2*lst%mymIndex(mIndex)-1,latIndex,levIndex) = gd_recv(icount,1,latIndex2,levIndex2,yourid+1)
              gd_out(2*lst%mymIndex(mIndex)  ,latIndex,levIndex) = gd_recv(icount,2,latIndex2,levIndex2,yourid+1)
            else
              gd_out(2*lst%mymIndex(mIndex)-1,latIndex,levIndex) = 0.d0
              gd_out(2*lst%mymIndex(mIndex)  ,latIndex,levIndex) = 0.d0
            end if
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    
    call tmg_stop(154)

  end subroutine transpose2d_LatToM

  !--------------------------------------------------------------------------
  ! lst_Transpose2d_LatToM_kij
  !--------------------------------------------------------------------------
  subroutine transpose2d_LatToM_kij(gd_out, gd_in, lst)
    implicit none

    type(struct_lst)     :: lst

    real(8), intent(out) :: gd_out(lst%myLevBeg:lst%myLevEnd,2*lst%mymCount,lst%nj+lst%njp )
    real(8), intent(in)  :: gd_in (lst%myLevBeg:lst%myLevEnd,lst%ni+lst%nip,lst%latPerPE)

    real(8) :: gd_recv(lst%maxLevCount,lst%maxmActiveCount,2,lst%latPerPEmax, mpi_npey)
    real(8) :: gd_send(lst%maxLevCount,lst%maxmActiveCount,2,lst%latPerPEmax, mpi_npey)
    integer :: yourid, mIndex, icount, nsize, ierr, levIndex, levIndex2, latIndex, latIndex2

    if (verbose) write(*,*) 'Entering transpose2d_LatToM_kij'
    call rpn_comm_barrier("GRID",ierr)

    call utl_tmg_start(154,'low-level--lst_transpose_MtoLAT')

    !$OMP PARALLEL DO PRIVATE(yourid,latIndex,levIndex,levIndex2,icount,mIndex)
    do yourid = 0, (mpi_npey-1)
      gd_send(:,:,:,:,yourid+1) = 0.d0
      do latIndex = 1, lst%latPerPE
        icount = 0
        do mIndex = lst%allmBeg(yourid+1), lst%allmEnd(yourid+1), lst%allmSkip(yourid+1)
          if (lst%KfromMNglb(mIndex,0) /= -1) then
            icount = icount + 1
            do levIndex = lst%myLevBeg, lst%myLevEnd
              levIndex2 = levIndex - lst%myLevBeg + 1
              gd_send(levIndex2,icount,1,latIndex,yourid+1) = gd_in(levIndex,2*mIndex+1, latIndex)
              gd_send(levIndex2,icount,2,latIndex,yourid+1) = gd_in(levIndex,2*mIndex+2, latIndex)
            end do
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    nsize = lst%maxmActiveCount * 2 * lst%maxLevCount * lst%latPerPEmax
    if (mpi_npey > 1) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","NS",ierr)
    else
      gd_recv(:,:,:,:,1) = gd_send(:,:,:,:,1)
    end if

    !$OMP PARALLEL DO PRIVATE(yourid,latIndex,latIndex2,levIndex,levIndex2,icount,mIndex)
    do yourid = 0, (mpi_npey-1)
      do latIndex = lst%allLatBeg(yourid+1), lst%allLatEnd(yourid+1)
        latIndex2 = latIndex - lst%allLatBeg(yourid+1) + 1
        icount = 0
        do mIndex = lst%mymBeg, lst%mymEnd, lst%mymSkip
          if (lst%KfromMNglb(mIndex,0) /= -1) then
            icount = icount + 1
            do levIndex = lst%myLevBeg, lst%myLevEnd
              levIndex2 = levIndex - lst%myLevBeg + 1
              gd_out(levIndex,2*lst%mymIndex(mIndex)-1,latIndex) = gd_recv(levIndex2,icount,1,latIndex2,yourid+1)
              gd_out(levIndex,2*lst%mymIndex(mIndex)  ,latIndex) = gd_recv(levIndex2,icount,2,latIndex2,yourid+1)
            end do
          else
            gd_out(:,2*lst%mymIndex(mIndex)-1,latIndex) = 0.d0
            gd_out(:,2*lst%mymIndex(mIndex)  ,latIndex) = 0.d0
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call tmg_stop(154)

  end subroutine transpose2d_LatToM_kij

  !--------------------------------------------------------------------------
  ! Transpose2d_MToLat
  !--------------------------------------------------------------------------
  subroutine transpose2d_MtoLat(gd_out, gd_in, lst)
    implicit none

    type(struct_lst)      :: lst

    real(8), intent(in)   :: gd_in (2*lst%mymCount,lst%nj+lst%njp,lst%myLevBeg:lst%myLevEnd)
    real(8), intent(out)  :: gd_out(lst%ni+lst%nip,lst%latPerPE,lst%myLevBeg:lst%myLevEnd)

    real(8) :: gd_recv(lst%maxmActiveCount,2,lst%latPerPEmax,lst%maxLevCount, mpi_npey)
    real(8) :: gd_send(lst%maxmActiveCount,2,lst%latPerPEmax,lst%maxLevCount, mpi_npey)
    integer :: yourid, mIndex, icount, nsize, ierr, levIndex, levIndex2, latIndex, latIndex2

    if (verbose) write(*,*) 'Entering transpose2d_MToLat'
    call rpn_comm_barrier("GRID",ierr)

    call utl_tmg_start(154,'low-level--lst_transpose_MtoLAT')

    !$OMP PARALLEL DO PRIVATE(yourid,latIndex,latIndex2,levIndex,levIndex2,icount,mIndex)
    do yourid = 0, (mpi_npey-1)
      do levIndex = lst%myLevBeg, lst%myLevEnd
        levIndex2 = levIndex - lst%myLevBeg + 1
        gd_send(:,:,:,levIndex2,yourid+1) = 0.d0
        do latIndex = lst%allLatBeg(yourid+1), lst%allLatEnd(yourid+1)
          latIndex2 = latIndex - lst%allLatBeg(yourid+1) + 1
          icount = 0
          do mIndex = lst%mymBeg, lst%mymEnd, lst%mymSkip
            if (lst%KfromMNglb(mIndex,0) /= -1) then
              icount = icount+1
              gd_send(icount,1,latIndex2,levIndex2,yourid+1) = gd_in(2*lst%mymIndex(mIndex)-1,latIndex,levIndex)
              gd_send(icount,2,latIndex2,levIndex2,yourid+1) = gd_in(2*lst%mymIndex(mIndex)  ,latIndex,levIndex)
            end if
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    nsize = lst%maxmActiveCount * 2 * lst%maxLevCount * lst%latPerPEmax
    if (mpi_npey > 1) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","NS",ierr)
    else
      gd_recv(:,:,:,:,1) = gd_send(:,:,:,:,1)
    end if

    !$OMP PARALLEL DO PRIVATE(yourid,latIndex,levIndex,levIndex2,icount,mIndex)
    do yourid = 0, (mpi_npey-1)
      do levIndex = lst%myLevBeg, lst%myLevEnd
        levIndex2 = levIndex - lst%myLevBeg + 1
        do latIndex = 1, lst%latPerPE
          icount = 0
          do mIndex = lst%allmBeg(yourid+1), lst%allmEnd(yourid+1), lst%allmSkip(yourid+1)
            if (lst%KfromMNglb(mIndex,0) /= -1) then
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

    call tmg_stop(154)
    
  end subroutine transpose2d_MtoLat

  !--------------------------------------------------------------------------
  ! Transpose2d_MToLat_kij
  !--------------------------------------------------------------------------
  subroutine transpose2d_MtoLat_kij(gd_out, gd_in, lst)
    implicit none

    type(struct_lst)      :: lst

    real(8), intent(in)   :: gd_in (lst%myLevBeg:lst%myLevEnd,2*lst%mymCount,lst%nj+lst%njp)
    real(8), intent(out)  :: gd_out(lst%myLevBeg:lst%myLevEnd,lst%ni+lst%nip,lst%latPerPE)

    real(8) :: gd_recv(lst%maxLevCount,lst%maxmActiveCount,2,lst%latPerPEmax, mpi_npey)
    real(8) :: gd_send(lst%maxLevCount,lst%maxmActiveCount,2,lst%latPerPEmax, mpi_npey)
    integer :: yourid, mIndex, icount, nsize, ierr, levIndex, levIndex2, latIndex, latIndex2

    if (verbose) write(*,*) 'Entering transpose2d_MToLat_kij'
    call rpn_comm_barrier("GRID",ierr)

    call utl_tmg_start(154,'low-level--lst_transpose_MtoLAT')

    !$OMP PARALLEL DO PRIVATE(yourid,latIndex,latIndex2,levIndex,levIndex2,icount,mIndex)
    do yourid = 0, (mpi_npey-1)
      gd_send(:,:,:,:,yourid+1) = 0.d0
      do latIndex = lst%allLatBeg(yourid+1), lst%allLatEnd(yourid+1)
        latIndex2 = latIndex - lst%allLatBeg(yourid+1) + 1
        icount = 0
        do mIndex = lst%mymBeg, lst%mymEnd, lst%mymSkip
          if (lst%KfromMNglb(mIndex,0) /= -1) then
            icount = icount+1
            do levIndex = lst%myLevBeg, lst%myLevEnd
              levIndex2 = levIndex - lst%myLevBeg + 1
              gd_send(levIndex2,icount,1,latIndex2,yourid+1) = gd_in(levIndex,2*lst%mymIndex(mIndex)-1,latIndex)
              gd_send(levIndex2,icount,2,latIndex2,yourid+1) = gd_in(levIndex,2*lst%mymIndex(mIndex)  ,latIndex)
            end do
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    nsize = lst%maxmActiveCount * 2 * lst%maxLevCount * lst%latPerPEmax
    if (mpi_npey > 1) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","NS",ierr)
    else
      gd_recv(:,:,:,:,1) = gd_send(:,:,:,:,1)
    end if

    !$OMP PARALLEL DO PRIVATE(yourid,latIndex,levIndex,levIndex2,icount,mIndex)
    do yourid = 0, (mpi_npey-1)
      do latIndex = 1, lst%latPerPE
        icount = 0
        do mIndex = lst%allmBeg(yourid+1), lst%allmEnd(yourid+1), lst%allmSkip(yourid+1)
          if (lst%KfromMNglb(mIndex,0) /= -1) then
            icount = icount+1
            do levIndex = lst%myLevBeg, lst%myLevEnd
              levIndex2 = levIndex - lst%myLevBeg + 1
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

    call tmg_stop(154)

  end subroutine transpose2d_MtoLat_kij

  !--------------------------------------------------------------------------
  ! lst_Transpose2d_LevToN
  !--------------------------------------------------------------------------
  subroutine transpose2d_LevToN(SpectralStateVar, gd_in, nk, lst)
    implicit none

    type(struct_lst)    :: lst

    integer, intent(in) :: nk
    real(8), intent(out):: SpectralStateVar(lst%nla,lst%nphase,nk)
    real(8), intent(in) :: gd_in (2*lst%mymCount, lst%nj+lst%njp, lst%myLevBeg:lst%myLevEnd)

    real(8) :: gd_send(lst%maxnla, 4, lst%maxLevCount, mpi_npex)
    real(8) :: gd_recv(lst%maxnla, 4, lst%maxLevCount, mpi_npex)
    integer :: yourid, nsize, ierr, levIndex, levIndex2, nIndex, mIndex, icount

    if (verbose) write(*,*) 'Entering transpose2d_LevToN'
    call rpn_comm_barrier("GRID",ierr)

    call utl_tmg_start(153,'low-level--lst_transpose_NtoLEV')

    !$OMP PARALLEL DO PRIVATE(yourid,mIndex,levIndex,levIndex2,nIndex,icount)
    do yourid = 0, (mpi_npex-1)
      do levIndex = lst%myLevBeg, lst%myLevEnd
        levIndex2 = levIndex - lst%myLevBeg + 1
        gd_send(:,:,levIndex2,yourid+1) = 0.d0
        icount = 0
        do nIndex = lst%allnBeg(yourid+1), lst%allnEnd(yourid+1), lst%allnSkip(yourid+1)
          do mIndex = lst%mymBeg, lst%mymEnd, lst%mymSkip 
            if (lst%KfromMNglb(mIndex,nIndex) /= -1) then
              icount = icount + 1
              gd_send(icount,1,levIndex2,yourid+1) = gd_in(2*lst%mymIndex(mIndex)-1,2*nIndex+1,levIndex)
              gd_send(icount,2,levIndex2,yourid+1) = gd_in(2*lst%mymIndex(mIndex)-1,2*nIndex+2,levIndex)
              gd_send(icount,3,levIndex2,yourid+1) = gd_in(2*lst%mymIndex(mIndex)  ,2*nIndex+1,levIndex)
              gd_send(icount,4,levIndex2,yourid+1) = gd_in(2*lst%mymIndex(mIndex)  ,2*nIndex+2,levIndex)
            end if
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    
    nsize = lst%maxnla * 4 * lst%maxLevCount
    if (mpi_npex > 1) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","EW",ierr)
    else
      gd_recv(:,:,:,1) = gd_send(:,:,:,1)
    end if

    !$OMP PARALLEL DO PRIVATE(yourid,levIndex,levIndex2)
    do yourid = 0, (mpi_npex-1)
      do levIndex = lst%allLevBeg(yourid+1), lst%allLevEnd(yourid+1)
        levIndex2 = levIndex - lst%allLevBeg(yourid+1) + 1
        SpectralStateVar(:,:,levIndex) = gd_recv(1:lst%nla,:,levIndex2,yourid+1)
      end do
    end do
    !$OMP END PARALLEL DO
    
    call tmg_stop(153)

  end subroutine transpose2d_LevToN

  !--------------------------------------------------------------------------
  ! lst_Transpose2d_LevToN_kij
  !--------------------------------------------------------------------------
  subroutine transpose2d_LevToN_kij(SpectralStateVar, gd_in, nk, lst)
    implicit none

    type(struct_lst)    :: lst

    integer, intent(in) :: nk
    real(8), intent(out):: SpectralStateVar(lst%nla,lst%nphase,nk)
    real(8), intent(in) :: gd_in (lst%myLevBeg:lst%myLevEnd,2*lst%mymCount, lst%nj+lst%njp)

    real(8) :: gd_send(lst%maxnla, 4, lst%maxLevCount, mpi_npex)
    real(8) :: gd_recv(lst%maxnla, 4, lst%maxLevCount, mpi_npex)
    integer :: yourid, nsize, ierr, levIndex, levIndex2, nIndex, mIndex, icount

    if (verbose) write(*,*) 'Entering transpose2d_LevToN_kij'
    call rpn_comm_barrier("GRID",ierr)

    call utl_tmg_start(153,'low-level--lst_transpose_NtoLEV')

    !$OMP PARALLEL DO PRIVATE(yourid,mIndex,levIndex,levIndex2,nIndex,icount)
    do yourid = 0, (mpi_npex-1)
      do levIndex = lst%myLevBeg, lst%myLevEnd
        levIndex2 = levIndex - lst%myLevBeg + 1
        gd_send(:,:,levIndex2,yourid+1) = 0.d0
        icount = 0
        do nIndex = lst%allnBeg(yourid+1), lst%allnEnd(yourid+1), lst%allnSkip(yourid+1)
          do mIndex = lst%mymBeg, lst%mymEnd, lst%mymSkip 
            if (lst%KfromMNglb(mIndex,nIndex) /= -1) then
              icount = icount + 1
              gd_send(icount,1,levIndex2,yourid+1) = gd_in(levIndex,2*lst%mymIndex(mIndex)-1,2*nIndex+1)
              gd_send(icount,2,levIndex2,yourid+1) = gd_in(levIndex,2*lst%mymIndex(mIndex)-1,2*nIndex+2)
              gd_send(icount,3,levIndex2,yourid+1) = gd_in(levIndex,2*lst%mymIndex(mIndex)  ,2*nIndex+1)
              gd_send(icount,4,levIndex2,yourid+1) = gd_in(levIndex,2*lst%mymIndex(mIndex)  ,2*nIndex+2)
            end if
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    nsize = lst%maxnla * 4 * lst%maxLevCount
    if (mpi_npex > 1) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","EW",ierr)
    else
      gd_recv(:,:,:,1) = gd_send(:,:,:,1)
    end if

    !$OMP PARALLEL DO PRIVATE(yourid,levIndex,levIndex2)
    do yourid = 0, (mpi_npex-1)
      do levIndex = lst%allLevBeg(yourid+1), lst%allLevEnd(yourid+1)
        levIndex2 = levIndex - lst%allLevBeg(yourid+1) + 1
        SpectralStateVar(:,:,levIndex) = gd_recv(1:lst%nla,:,levIndex2,yourid+1)
      end do
    end do
    !$OMP END PARALLEL DO
    
    call tmg_stop(153)

  end subroutine transpose2d_LevToN_kij

  !--------------------------------------------------------------------------
  ! lst_Transpose2d_NToLev
  !--------------------------------------------------------------------------
  subroutine transpose2d_NToLev(gd_out, SpectralStateVar, nk, lst)
    implicit none

    type(struct_lst)    :: lst

    integer, intent(in) :: nk
    real(8), intent(in):: SpectralStateVar(lst%nla,lst%nphase,nk)
    real(8), intent(out):: gd_out(2*lst%mymCount, lst%nj+lst%njp, lst%myLevBeg:lst%myLevEnd)

    real(8) :: gd_send(lst%maxnla, 4, lst%maxLevCount, mpi_npex)
    real(8) :: gd_recv(lst%maxnla, 4, lst%maxLevCount, mpi_npex)
    integer :: yourid, nsize, ierr, levIndex, levIndex2, nIndex, mIndex, icount

    if (verbose) write(*,*) 'Entering transpose2d_NToLev'
    call rpn_comm_barrier("GRID",ierr)

    call utl_tmg_start(153,'low-level--lst_transpose_NtoLEV')

    !$OMP PARALLEL DO PRIVATE(yourid,levIndex,levIndex2)
    do yourid = 0, (mpi_npex-1)
      do levIndex = lst%allLevBeg(yourid+1), lst%allLevEnd(yourid+1)
        levIndex2 = levIndex - lst%allLevBeg(yourid+1) + 1
        gd_send(1:lst%nla,:,levIndex2,yourid+1) = SpectralStateVar(:,:,levIndex)
        if (lst%nla < lst%maxnla) gd_send(lst%nla+1:lst%maxnla,:,levIndex2,yourid+1) = 0.d0
      end do
    end do
    !$OMP END PARALLEL DO

    nsize = lst%maxnla * 4* lst%maxLevCount
    if (mpi_npex > 1) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","EW",ierr)
    else
      gd_recv(:,:,:,1) = gd_send(:,:,:,1)
    end if

    !$OMP PARALLEL DO PRIVATE(yourid,nIndex,levIndex,levIndex2,mIndex,icount)
    do yourid = 0, (mpi_npex-1)
      do levIndex = lst%myLevBeg, lst%myLevEnd
        levIndex2 = levIndex - lst%myLevBeg + 1
        icount = 0
        do nIndex = lst%allnBeg(yourid+1), lst%allnEnd(yourid+1), lst%allnSkip(yourid+1)
          do mIndex = lst%mymBeg, lst%mymEnd, lst%mymSkip
            if (lst%KfromMNglb(mIndex,nIndex) /= -1) then
              icount = icount + 1
              gd_out(2*lst%mymIndex(mIndex)-1,2*nIndex+1,levIndex) = gd_recv(icount,1,levIndex2,yourid+1)
              gd_out(2*lst%mymIndex(mIndex)-1,2*nIndex+2,levIndex) = gd_recv(icount,2,levIndex2,yourid+1)
              gd_out(2*lst%mymIndex(mIndex)  ,2*nIndex+1,levIndex) = gd_recv(icount,3,levIndex2,yourid+1)
              gd_out(2*lst%mymIndex(mIndex)  ,2*nIndex+2,levIndex) = gd_recv(icount,4,levIndex2,yourid+1)
            else
              gd_out(2*lst%mymIndex(mIndex)-1,2*nIndex+1,levIndex) = 0.d0
              gd_out(2*lst%mymIndex(mIndex)-1,2*nIndex+2,levIndex) = 0.d0
              gd_out(2*lst%mymIndex(mIndex)  ,2*nIndex+1,levIndex) = 0.d0
              gd_out(2*lst%mymIndex(mIndex)  ,2*nIndex+2,levIndex) = 0.d0
            end if
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call tmg_stop(153)
    
  end subroutine transpose2d_NToLev

  !--------------------------------------------------------------------------
  ! lst_Transpose2d_NToLev_kij
  !--------------------------------------------------------------------------
  subroutine transpose2d_NToLev_kij(gd_out, SpectralStateVar, nk, lst)
    implicit none

    type(struct_lst)    :: lst

    integer, intent(in) :: nk
    real(8), intent(in) :: SpectralStateVar(lst%nla,lst%nphase,nk)
    real(8), intent(out):: gd_out(lst%myLevBeg:lst%myLevEnd,2*lst%mymCount,lst%nj+lst%njp)

    real(8) :: gd_send(lst%maxnla, 4, lst%maxLevCount, mpi_npex)
    real(8) :: gd_recv(lst%maxnla, 4, lst%maxLevCount, mpi_npex)
    integer :: yourid, nsize, ierr, levIndex, levIndex2, nIndex, mIndex, icount

    if (verbose) write(*,*) 'Entering transpose2d_NToLev_kij'
    call rpn_comm_barrier("GRID",ierr)

    call utl_tmg_start(153,'low-level--lst_transpose_NtoLEV')

    !$OMP PARALLEL DO PRIVATE(yourid,levIndex,levIndex2)
    do yourid = 0, (mpi_npex-1)
      do levIndex = lst%allLevBeg(yourid+1), lst%allLevEnd(yourid+1)
        levIndex2 = levIndex - lst%allLevBeg(yourid+1) + 1
        gd_send(1:lst%nla,:,levIndex2,yourid+1) = SpectralStateVar(:,:,levIndex)
        if (lst%nla < lst%maxnla) gd_send(lst%nla+1:lst%maxnla,:,levIndex2,yourid+1) = 0.d0
      end do
    end do
    !$OMP END PARALLEL DO
    
    nsize = lst%maxnla * 4 * lst%maxLevCount
    if (mpi_npex > 1) then
      call rpn_comm_alltoall(gd_send,nsize,"mpi_double_precision",  &
                             gd_recv,nsize,"mpi_double_precision","EW",ierr)
    else
      gd_recv(:,:,:,1) = gd_send(:,:,:,1)
    end if

    !$OMP PARALLEL DO PRIVATE(yourid,nIndex,levIndex,levIndex2,mIndex,icount)
    do yourid = 0, (mpi_npex-1)
      do levIndex = lst%myLevBeg, lst%myLevEnd
        levIndex2 = levIndex - lst%myLevBeg + 1
        icount = 0
        do nIndex = lst%allnBeg(yourid+1), lst%allnEnd(yourid+1), lst%allnSkip(yourid+1)
          do mIndex = lst%mymBeg, lst%mymEnd, lst%mymSkip
            if (lst%KfromMNglb(mIndex,nIndex) /= -1) then
              icount = icount + 1
              gd_out(levIndex,2*lst%mymIndex(mIndex)-1,2*nIndex+1) = gd_recv(icount,1,levIndex2,yourid+1)
              gd_out(levIndex,2*lst%mymIndex(mIndex)-1,2*nIndex+2) = gd_recv(icount,2,levIndex2,yourid+1)
              gd_out(levIndex,2*lst%mymIndex(mIndex)  ,2*nIndex+1) = gd_recv(icount,3,levIndex2,yourid+1)
              gd_out(levIndex,2*lst%mymIndex(mIndex)  ,2*nIndex+2) = gd_recv(icount,4,levIndex2,yourid+1)
            else
              gd_out(levIndex,2*lst%mymIndex(mIndex)-1,2*nIndex+1) = 0.d0
              gd_out(levIndex,2*lst%mymIndex(mIndex)-1,2*nIndex+2) = 0.d0
              gd_out(levIndex,2*lst%mymIndex(mIndex)  ,2*nIndex+1) = 0.d0
              gd_out(levIndex,2*lst%mymIndex(mIndex)  ,2*nIndex+2) = 0.d0
            end if
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call tmg_stop(153)

  end subroutine transpose2d_NToLev_kij

  !--------------------------------------------------------------------------
  ! lst_ReshapeTrunc
  !--------------------------------------------------------------------------
  subroutine lst_ReshapeTrunc(SpectralStateRpn, SpectralStateVar,  &
                               Direction, kStart, kEnd, lst)
    implicit none

    type(struct_lst)                :: lst

    integer,          intent(in)    :: kStart, kEnd
    character(len=*), intent(in)    :: Direction ! ToVAR or ToRPN
    real(8),          intent(inout) :: SpectralStateRpn(2*lst%mymCount,2*lst%mynCount,kStart:kEnd)
    real(8),          intent(inout) :: SpectralStateVar(lst%nla     ,lst%nphase,kStart:kEnd)

    integer k, m, n, ila

    if (verbose) write(*,*) 'Entering lst_ReshapeTrunc'

    select case (trim(Direction))
    case ('ToVAR')
      ! Truncation (if applicable) will be applied here
      !$OMP PARALLEL DO PRIVATE (n,m,ila,k) 
      do n = lst%mynBeg, lst%mynEnd, lst%mynSkip
        do m = lst%mymBeg, lst%mymEnd, lst%mymSkip
          ila = lst%nla_Index(m,n)
          if (ila /= -1) then
            do k = kStart, kEnd 
              SpectralStateVar(ila,1,k) = SpectralStateRpn(2*lst%mymIndex(m)-1,2*lst%mynIndex(n)-1,k)
              SpectralStateVar(ila,2,k) = SpectralStateRpn(2*lst%mymIndex(m)-1,2*lst%mynIndex(n)  ,k)
              SpectralStateVar(ila,3,k) = SpectralStateRpn(2*lst%mymIndex(m)  ,2*lst%mynIndex(n)-1,k)
              SpectralStateVar(ila,4,k) = SpectralStateRpn(2*lst%mymIndex(m)  ,2*lst%mynIndex(n),  k)
            end do
          end if
        end do
      end do
      !$OMP END PARALLEL DO

    case ('ToRPN')
      SpectralStateRpn(:,:,:) = 0.0d0
      !$OMP PARALLEL DO PRIVATE (n,m,ila,k)
      do n = lst%mynBeg, lst%mynEnd, lst%mynSkip
        do m = lst%mymBeg, lst%mymEnd, lst%mymSkip
          ila = lst%nla_Index(m,n)
          if (ila /= -1) then
            do k = kStart, kEnd
              SpectralStateRpn(2*lst%mymIndex(m)-1,2*lst%mynIndex(n)-1,k) = SpectralStateVar(ila,1,k)
              SpectralStateRpn(2*lst%mymIndex(m)-1,2*lst%mynIndex(n)  ,k) = SpectralStateVar(ila,2,k)
              SpectralStateRpn(2*lst%mymIndex(m)  ,2*lst%mynIndex(n)-1,k) = SpectralStateVar(ila,3,k)
              SpectralStateRpn(2*lst%mymIndex(m)  ,2*lst%mynIndex(n)  ,k) = SpectralStateVar(ila,4,k)
            end do
          end if
        end do
      end do
      !$OMP END PARALLEL DO

    case default
      write(*,*)
      write(*,*) 'lst_ReshapeTrunc: Unknown Direction', trim(Direction)
      call utl_abort('lst_ReshapeTrunc')
    end select

  end subroutine lst_ReshapeTrunc

  !--------------------------------------------------------------------------
  ! lst_ReshapeTrunc_kij
  !--------------------------------------------------------------------------
  subroutine lst_ReshapeTrunc_kij(SpectralStateRpn, SpectralStateVar,  &
                                   Direction, kStart, kEnd, lst)
    implicit none

    type(struct_lst)                :: lst

    integer,          intent(in)    :: kStart, kEnd
    character(len=*), intent(in)    :: Direction ! ToVAR or ToRPN
    real(8),          intent(inout) :: SpectralStateRpn(kStart:kEnd,2*lst%mymCount,2*lst%mynCount)
    real(8),          intent(inout) :: SpectralStateVar(lst%nla     ,lst%nphase,kStart:kEnd)

    integer k, m, n, ila

    if (verbose) write(*,*) 'Entering lst_ReshapeTrunc_kij'

    select case (trim(Direction))
    case ('ToVAR')
      ! Truncation (if applicable) will be applied here
      !$OMP PARALLEL DO PRIVATE (n,m,ila,k) 
      do n = lst%mynBeg, lst%mynEnd, lst%mynSkip
        do m = lst%mymBeg, lst%mymEnd, lst%mymSkip
          ila = lst%nla_Index(m,n)
          if (ila /= -1) then
            do k = kStart, kEnd 
              SpectralStateVar(ila,1,k) = SpectralStateRpn(k,2*lst%mymIndex(m)-1,2*lst%mynIndex(n)-1)
              SpectralStateVar(ila,2,k) = SpectralStateRpn(k,2*lst%mymIndex(m)-1,2*lst%mynIndex(n) )
              SpectralStateVar(ila,3,k) = SpectralStateRpn(k,2*lst%mymIndex(m)  ,2*lst%mynIndex(n)-1)
              SpectralStateVar(ila,4,k) = SpectralStateRpn(k,2*lst%mymIndex(m)  ,2*lst%mynIndex(n) )
            end do
          end if
        end do
      end do
      !$OMP END PARALLEL DO

    case ('ToRPN')
      SpectralStateRpn(:,:,:) = 0.0d0
      !$OMP PARALLEL DO PRIVATE (n,m,ila,k)
      do n = lst%mynBeg, lst%mynEnd, lst%mynSkip
        do m = lst%mymBeg, lst%mymEnd, lst%mymSkip
          ila = lst%nla_Index(m,n)
          if (ila /= -1) then
            do k = kStart, kEnd
              SpectralStateRpn(k,2*lst%mymIndex(m)-1,2*lst%mynIndex(n)-1) = SpectralStateVar(ila,1,k)
              SpectralStateRpn(k,2*lst%mymIndex(m)-1,2*lst%mynIndex(n) ) = SpectralStateVar(ila,2,k)
              SpectralStateRpn(k,2*lst%mymIndex(m)  ,2*lst%mynIndex(n)-1) = SpectralStateVar(ila,3,k)
              SpectralStateRpn(k,2*lst%mymIndex(m)  ,2*lst%mynIndex(n) ) = SpectralStateVar(ila,4,k)
            end do
          end if
        end do
      end do
      !$OMP END PARALLEL DO

    case default
      write(*,*)
      write(*,*) 'lst_ReshapeTrunc_kij: Unknown Direction', trim(Direction)
      call utl_abort('lst_ReshapeTrunc_kij')
    end select

  end subroutine lst_ReshapeTrunc_kij

  !--------------------------------------------------------------------------
  ! lst_Laplacian
  !--------------------------------------------------------------------------
  subroutine lst_Laplacian(lst, GridState, Mode, nk)
    implicit none

    type(struct_lst), intent(in)    :: lst

    integer,          intent(in)    :: nk
                                     ! Grid point data dimensions
    real(8),          intent(inout) :: GridState(lst%myLonBeg:lst%myLonEnd,lst%myLatBeg:lst%myLatEnd,nk)  
                                     ! 3D field in grid point space
    character(len=*), intent(in)    :: Mode
                                     ! Forward or Inverse

    real(8), allocatable            :: SpectralStateVar(:,:,:)
    real(8), allocatable            :: factor(:)

    integer :: k, ila, p

    character(len=24)   :: kind

    if (verbose) write(*,*) 'Entering lst_Laplacian'

    allocate(SpectralStateVar(lst%nla,lst%nphase,nk))
    allocate(factor(lst%nla))

    !
    !- 1.  Set Mode-dependent factors
    !
    select case (trim(Mode))
    case ('Forward')
      factor(:) = lst%lapxy(:)
    case ('Inverse')
      factor(:) = lst%ilapxy(:)
    case default
      write(*,*)
      write(*,*) 'lst_Laplacian: Error: Mode Unknown ', trim(Mode)
      call utl_abort('lst_Laplacian')
    end select

    !
    !- 2. Grid Point Space -> Spectral Space
    !
    kind = 'GridPointToSpectral'
    call lst_VarTransform(lst,                   & ! IN
                          SpectralStateVar,      & ! OUT
                          GridState,             & ! IN
                          kind, nk)                ! IN    

    !
    !- 3. Laplacian (forward or inverse) Transform
    !
    !$OMP PARALLEL DO PRIVATE (k,ila,p)
    do k = 1, nk
      do ila = 1, lst%nla
        do p = 1, lst%nphase
          SpectralStateVar(ila,p,k) = factor(ila) * SpectralStateVar(ila,p,k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    !
    !- 4. Spectral Space -> Grid Point Space
    !
    kind = 'SpectralToGridPoint'
    call lst_VarTransform(lst,                   & ! IN
                          SpectralStateVar,      & ! IN
                          GridState,             & ! OUT
                          kind, nk)                ! IN

    deallocate(SpectralStateVar)
    deallocate(factor)

  end subroutine lst_Laplacian

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

    if (n <= m) n = m + 1
    n = n - 1
1   n = n + 1
    i = n
2   do j = 1, l
       if (mod(i,k(j)) == 0) go to 4
    end do
    go to 1
4   i = i/k(j)
    if (i /= 1) go to 2

  end subroutine ngfft

end module lamSpectralTransform_mod

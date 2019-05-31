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
!CANADA, H9P 1J3; or send e-mail to ec.service.rpn.ec@canada.ca
!-------------------------------------- LICENCE END --------------------------------------

module diffusion_mod
  ! MODULE diffusion_mod (prefix='diff' category='3. High-level transformations')
  !
  ! :Purpose: Storage for data required by the diffusion operator used to model
  !           background-error horizontal correlations.
  !
  ! :Reference: Weaver, A. T., and P. Courtier, 2001: Correlation modelling on
  !             the sphere using a generalized diffusion equation.
  !             Q. J. R. Meteorol. Soc., 127, 1815-1846.
  !
  ! :Basic equations: Lcorr^2 = 2*k*dt*numt   (1)
  !                   stab    = k*dt/dx^2     (2)
  !
  ! :Origin: diffus_mod.f from the envar experiments by Anna Shlyaeva
  !

  implicit none
  private

  ! Public subroutines and functions
  public :: diff_setup, diff_finalize, diff_Csqrt, diff_Csqrtadj

  type struct_diff
     ! number of grid points in x and y directions (includes perimeter of land points)
     integer :: ni, nj
     integer :: numt
     real(8) :: dt
     real(8), allocatable :: cosyhalf(:), cosyinv(:), cosyinvsq(:)
     real(8), allocatable :: mhalfx(:,:), mhalfy(:,:)
     real(8), allocatable :: khalfx(:,:), khalfy(:,:)
     real(8) :: dlon, dlat        ! grid spacing in radians
     ! These are used in subroutines diff_Csqrt and diff_Csqrtadj
     real(8), allocatable :: Winv(:,:), Wsqrt(:,:), Winvsqrt(:,:)
     real(8), allocatable :: Lambda(:,:)
     real(8), allocatable :: diff1x_ap(:,:),diff1x_bp_inv(:,:)
     real(8), allocatable :: diff1x_c(:,:)
     real(8), allocatable :: diff1y_ap(:,:),diff1y_bp_inv(:,:)
     real(8), allocatable :: diff1y_c(:,:)
     logical :: limplicit
  end type struct_diff

  integer, parameter :: nMaxDiff = 10
  integer            :: nDiffAlreadyAllocated = 0
  type(struct_diff)  :: diff(nMaxDiff)

!*************************************************************************

contains

  integer function diff_setup(hco, corr_len, stab, nsamp, limplicit)

    use horizontalCoord_mod
    use EarthConstants_mod, only : rayt
    use randomNumber_mod
    use utilities_mod

    implicit none

    ! Arguments:
    type(struct_hco), pointer :: hco
    real, intent(in) :: corr_len     ! Horizontal correlation length scale (km)
    real, intent(in) :: stab         ! Stability criteria (definitely < 0.5)
    integer, intent(in) :: nsamp     ! Number of samples in the estimation of the normalization factors by randomization.
    logical, intent(in) :: limplicit ! Indicate to use the implicit formulation of the diffusion operator (.true.) or the explicit version (.false.).

    ! Locals:

    ! latitudes on the analysis rotated grid, in radians
    real(8), allocatable :: latr(:)

    real, allocatable :: buf2d(:,:)

    integer :: lonIndex, latIndex, isamp, k, l, t

    real(8) :: mindxy, maxL

    real(8) :: a,b

    real(8), allocatable :: Lcorr(:,:)
    real(8), allocatable :: kappa(:,:)
    real(8), allocatable :: W(:,:)
    real(8), allocatable :: m(:,:)
    real(8), allocatable :: xin(:,:)

    ! diff_norm_fact is the name of the RPN format file for the normalization factors.
    character(len=*), parameter :: diff_norm_fact = './diffusmod.std'

    ! Variables and functions required to write to RPN Standard files.

    integer  :: nmax

    integer  :: std_unit, ierr, key, npak, datyp
    integer  :: nii, njj, nkk
    integer  :: dateo
    integer  :: deet, npas
    integer  :: ip1, ip2, ip3
    integer  :: ig1, ig2, ig3, ig4
    character(len=1) :: grtyp
    character(len=2) :: typvar
    character(len=12) :: etiket
    logical  :: rewrit, file_exist
    real     :: dumwrk(1)

    integer  :: fnom,fstouv,fstecr,fstlir,fstfrm,fclos
    external :: fnom,fstouv,fstecr,fstlir,fstfrm,fclos

    integer  :: diffID

    integer :: NI, NJ


    if(nDiffAlreadyAllocated == nMaxDiff) then
      write(*,*) 'diff_setup: The maximum number of diffusion operators have already been allocated! ',nMaxDiff
      call utl_abort('diff_setup')
    endif

    nDiffAlreadyAllocated = nDiffAlreadyAllocated + 1
    diffID = nDiffAlreadyAllocated

    NI = hco%ni
    NJ = hco%nj
    diff(diffID)%dlon = hco%dlon
    diff(diffID)%dlat = hco%dlat
    diff(diffID)%limplicit = limplicit

    allocate(diff(diffID)%cosyhalf(NJ), diff(diffID)%cosyinv(NJ), diff(diffID)%cosyinvsq(NJ))
    allocate(diff(diffID)%Winv(NI,NJ), diff(diffID)%Wsqrt(NI,NJ), diff(diffID)%Winvsqrt(NI,NJ))
    allocate(diff(diffID)%khalfx(NI-1,NJ), diff(diffID)%khalfy(NI,NJ-1))
    allocate(diff(diffID)%mhalfx(NI-1,NJ), diff(diffID)%mhalfy(NI,NJ-1))
    allocate(diff(diffID)%Lambda(NI,NJ))

    allocate(latr(NJ))
    allocate(Lcorr(NI,NJ))
    allocate(kappa(NI,NJ))
    allocate(W(NI,NJ))
    allocate(m(NI,NJ))
    allocate(xin(NI,NJ))

    allocate(diff(diffID)%diff1x_ap(NI,NJ))
    allocate(diff(diffID)%diff1x_bp_inv(NI,NJ))
    allocate(diff(diffID)%diff1x_c(NI,NJ))
    allocate(diff(diffID)%diff1y_ap(NJ,NI))
    allocate(diff(diffID)%diff1y_bp_inv(NJ,NI))
    allocate(diff(diffID)%diff1y_c(NJ,NI))

    latr(:) = hco%lat(:)

    diff(diffID)%cosyinv(:) = 1.0d0/cos(latr(:))
    diff(diffID)%cosyinvsq(:) = diff(diffID)%cosyinv(:)*diff(diffID)%cosyinv(:)

    ! cosinus of latitudes on staggered grid
    diff(diffID)%cosyhalf(:) = cos(latr(:) + 0.5d0*diff(diffID)%dlat)

    ! minimum grid spacing over the domain (radians, but with cos(latr))
    mindxy = min(minval(cos(latr(:)))*diff(diffID)%dlon,diff(diffID)%dlat)

    Lcorr(:,:) = corr_len/(rayt/1000.0) ! lengthscale in radians
    maxL = maxval(Lcorr(2:NI-1,2:NJ-1)) ! maximum lengthscale over domain

    ! set main parameters for diffusion operator
    kappa(:,:) = Lcorr(:,:)**2          ! arbitrarily set k to L^2 (in radians)
    diff(diffID)%dt = stab*(mindxy**2)/(maxL**2) ! determine dt from stability criteria (2)
!    diff(diffID)%numt = 1.0d0/(2.0d0*diff(diffID)%dt)       ! determine number of timesteps from (1)
    diff(diffID)%numt = ceiling(1.0d0/(4.0d0*diff(diffID)%dt))*2; ! make sure it is an even integer
    diff(diffID)%dt = 1.0d0/(2.0d0*dble(diff(diffID)%numt)); ! recompute dt

    ! interpolate diffusion coefficient onto 2 staggered lat-lon grids
    do latIndex = 1, NJ
       do lonIndex = 1, NI-1
          diff(diffID)%khalfx(lonIndex,latIndex) = (kappa(lonIndex,latIndex) + kappa(lonIndex+1,latIndex))/2.0d0
       end do
    end do
    do latIndex = 1, NJ-1
       do lonIndex = 1, NI
          diff(diffID)%khalfy(lonIndex,latIndex) = (kappa(lonIndex,latIndex) + kappa(lonIndex,latIndex+1))/2.0d0
       end do
    end do

    ! print this stuff in listing file for user information:
    if ( .not. limplicit) then
       write(*,*)
       write(*,*) 'Number of timesteps = ',diff(diffID)%numt 
       write(*,*) 'Stability           = ',maxval(kappa)*diff(diffID)%dt/(mindxy**2)
       write(*,*)
    end if

    ! this is the matrix necessary for defining the inner product: for lat-lon grid, only cos(y)
    ! Actually, only the diagonal of the matrix is stored in the array W, since the matrix is diagonal.
    W(1,:) = cos(latr(:))
    do lonIndex = 2, NI
       W(lonIndex,:) = W(1,:)
    end do
    diff(diffID)%Winv(:,:) = 1.0d0/W(:,:)
    diff(diffID)%Wsqrt(:,:) = sqrt(W(:,:))
    diff(diffID)%Winvsqrt(:,:) = 1.0d0/diff(diffID)%Wsqrt(:,:)

    ! Get mask from hco

    ! land mask (1=water, 0=land)
    do latIndex = 1, NJ
       do lonIndex = 1, NI
          if(hco%mask(lonIndex,latIndex) == 1) then
             m(lonIndex,latIndex) = 1.0d0
          else
             m(lonIndex,latIndex) = 0.0d0
          end if
       end do
    end do
    m(:,1) = 0.0d0
    m(1,:) = 0.0d0
    m(:,NJ) = 0.0d0
    m(NI,:) = 0.0d0

    ! define mask on staggered grids
    do latIndex = 1, NJ
       do lonIndex = 1, NI-1
          if(sum(m(lonIndex:lonIndex+1,latIndex)) < 2.0d0) then
             diff(diffID)%mhalfx(lonIndex,latIndex) = 0.0d0
          else
             diff(diffID)%mhalfx(lonIndex,latIndex) = 1.0d0
          end if
       end do
    end do

    do latIndex = 1, NJ-1
       do lonIndex = 1, NI
          if(sum(m(lonIndex,latIndex:latIndex+1)) < 2.0d0) then
             diff(diffID)%mhalfy(lonIndex,latIndex) = 0.0d0
          else
             diff(diffID)%mhalfy(lonIndex,latIndex) = 1.0d0
          end if
       end do
    end do

    diff(diffID)%ni = ni
    diff(diffID)%nj = nj

    ! specify number of timesteps and timestep length for implicit 1D diffusion
    if (limplicit) then
       diff(diffID)%numt = 5
       diff(diffID)%dt = 1.0d0/(2.0d0*dble(2*diff(diffID)%numt)-3.0d0)
    end if

    ! compute the LU decomposition for the implicit 1D diffusion
    diff(diffID)%diff1x_ap(:,:) = 0.0d0
    diff(diffID)%diff1x_bp_inv(:,:) = 0.0d0      
    diff(diffID)%diff1x_c(:,:) = 0.0d0
    diff(diffID)%diff1y_ap(:,:) = 0.0d0
    diff(diffID)%diff1y_bp_inv(:,:) = 0.0d0      
    diff(diffID)%diff1y_c(:,:) = 0.0d0

    !$OMP PARALLEL DO PRIVATE(latIndex,lonIndex,a,b)
    do latIndex = 2, NJ-1
       lonIndex = 2
       diff(diffID)%diff1x_bp_inv(lonIndex,latIndex)= 1.0d0/(1.0d0 + &
            diff(diffID)%dt*diff(diffID)%cosyinvsq(latIndex)*(diff(diffID)%mhalfx(lonIndex,latIndex)*diff(diffID)%khalfx(lonIndex,latIndex) + &
            diff(diffID)%mhalfx(lonIndex-1,latIndex)*diff(diffID)%khalfx(lonIndex-1,latIndex))/(diff(diffID)%dlon*diff(diffID)%dlon))
       do lonIndex = 3, NI-1
          ! elements of the tri-diagonal coefficient matrix
          a = - diff(diffID)%dt*diff(diffID)%cosyinvsq(latIndex)*diff(diffID)%mhalfx(lonIndex-1,latIndex)*diff(diffID)%khalfx(lonIndex-1,latIndex)/(diff(diffID)%dlon*diff(diffID)%dlon)
          b = 1 + diff(diffID)%dt*diff(diffID)%cosyinvsq(latIndex)*(diff(diffID)%mhalfx(lonIndex,latIndex)*diff(diffID)%khalfx(lonIndex,latIndex) + &
               diff(diffID)%mhalfx(lonIndex-1,latIndex)*diff(diffID)%khalfx(lonIndex-1,latIndex))/(diff(diffID)%dlon*diff(diffID)%dlon)
          diff(diffID)%diff1x_c(lonIndex,latIndex) = -diff(diffID)%dt*diff(diffID)%cosyinvsq(latIndex)*diff(diffID)%mhalfx(lonIndex,latIndex)*diff(diffID)%khalfx(lonIndex,latIndex)/(diff(diffID)%dlon*diff(diffID)%dlon)
          diff(diffID)%diff1x_ap(lonIndex,latIndex) = a*diff(diffID)%diff1x_bp_inv(lonIndex-1,latIndex)
          diff(diffID)%diff1x_bp_inv(lonIndex,latIndex) = 1.0d0/(b-a*diff(diffID)%diff1x_c(lonIndex-1,latIndex) &
               *diff(diffID)%diff1x_bp_inv(lonIndex-1,latIndex))
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(latIndex,lonIndex,a,b)
    do lonIndex = 2, NI-1
       latIndex = 2
       diff(diffID)%diff1y_bp_inv(latIndex,lonIndex) = 1.0d0/(1.0d0 + &
            diff(diffID)%dt*diff(diffID)%cosyinvsq(latIndex)*(diff(diffID)%cosyhalf(latIndex)*diff(diffID)%mhalfy(lonIndex,latIndex) &
            *diff(diffID)%khalfy(lonIndex,latIndex) + &
            diff(diffID)%cosyhalf(latIndex-1)*diff(diffID)%mhalfy(lonIndex,latIndex-1)*diff(diffID)%khalfy(lonIndex,latIndex-1))/(diff(diffID)%dlat*diff(diffID)%dlat))
       do latIndex = 3, NJ-1
          ! elements of the tri-diagonal coefficient matrix
          a = - diff(diffID)%dt*diff(diffID)%cosyinv(latIndex)*diff(diffID)%cosyhalf(latIndex-1)* &
               diff(diffID)%mhalfy(lonIndex,latIndex-1)*diff(diffID)%khalfy(lonIndex,latIndex-1)/(diff(diffID)%dlat*diff(diffID)%dlat)
          b = 1 + diff(diffID)%dt*diff(diffID)%cosyinv(latIndex)*(diff(diffID)%cosyhalf(latIndex)*diff(diffID)%mhalfy(lonIndex,latIndex)* &
               diff(diffID)%khalfy(lonIndex,latIndex) + &
               diff(diffID)%cosyhalf(latIndex-1)*diff(diffID)%mhalfy(lonIndex,latIndex-1)*diff(diffID)%khalfy(lonIndex,latIndex-1)) &
               /(diff(diffID)%dlat*diff(diffID)%dlat)
          diff(diffID)%diff1y_c(latIndex,lonIndex) = - diff(diffID)%dt*diff(diffID)%cosyinv(latIndex)*diff(diffID)%cosyhalf(latIndex)* &
               diff(diffID)%mhalfy(lonIndex,latIndex)*diff(diffID)%khalfy(lonIndex,latIndex)/(diff(diffID)%dlat*diff(diffID)%dlat)
          diff(diffID)%diff1y_ap(latIndex,lonIndex) = a*diff(diffID)%diff1y_bp_inv(latIndex-1,lonIndex)
          diff(diffID)%diff1y_bp_inv(latIndex,lonIndex) = 1.0d0/(b-a*diff(diffID)%diff1y_c(latIndex-1,lonIndex) &
               *diff(diffID)%diff1y_bp_inv(latIndex-1,lonIndex))
       end do
    end do
    !$OMP END PARALLEL DO

    if (limplicit) then
       write (etiket, FMT='(''KM'',i3.3,''IMPLICI'')') int(corr_len)
    else
       write (etiket, FMT='(''KM'',i3.3,''STAB'',f3.1)') int(corr_len), stab
    end if

    std_unit = 0

    inquire (file=diff_norm_fact, exist=file_exist)

    if(file_exist) then

       ierr = fnom(std_unit, diff_norm_fact, 'RND+R/O', 0)
       nmax = fstouv(std_unit, 'RND')

       allocate(buf2d(NI,NJ))

       key = fstlir(buf2d, std_unit, nii, njj, nkk,    &
            -1, etiket, -1, -1, -1, ' ', 'LAMB')

       ierr = fstfrm(std_unit)
       ierr = fclos(std_unit)

       diff(diffID)%Lambda = dble(buf2d)

       deallocate(buf2d)

    end if

    if(nii /= ni .or. njj /= nj .or. key <= 0 .or.     &
         (.not. file_exist)) then

       call rng_setup(1)

       if(nsamp < NI*NJ) then

          ! compute normalization:  Lambda = inverse stddev of (Diffuse * W^-1/2)
          write(*,*) 'Randomization estimation of the '//   &
               'normalization for diffusion...'
          write(*,*) 'Will use ',nsamp,' samples.'
          call flush(6)
          diff(diffID)%Lambda = 0.0d0
          do isamp = 1, nsamp

             do latIndex = 1, NJ
                do lonIndex = 1, NI
                   xin(lonIndex,latIndex) = diff(diffID)%Winvsqrt(lonIndex,latIndex)*rng_gaussian()
                end do
             end do
             if(limplicit) then
                do t = 1, diff(diffID)%numt
                   call diffusion1x_implicit(diffID,xin,xin)
                   call diffusion1y_implicit(diffID,xin,xin)
                end do
             else
                call diffusion_explicit(diffID,xin,xin)
             end if
             diff(diffID)%Lambda = diff(diffID)%Lambda + xin*xin

          end do

          do latIndex = 1, NJ
             do lonIndex = 1, NI
                diff(diffID)%Lambda(lonIndex,latIndex) = sqrt(diff(diffID)%Lambda(lonIndex,latIndex)/dble(nsamp-1)) ! normalization: inverse of rms of ens
                if(diff(diffID)%Lambda(lonIndex,latIndex) > 0.0d0) then
                   diff(diffID)%Lambda(lonIndex,latIndex) = 1.0d0/diff(diffID)%Lambda(lonIndex,latIndex);
                end if
             end do
          end do

       else

          write(*,*) 'Exact calculation of the '//      &
               'normalization for diffusion...'
          call flush(6)

          do latIndex = 1, NJ
             write(*,*) 'Doing row latIndex = ',latIndex,' of ',NJ
             call flush(6)
             do lonIndex = 1, NI

                xin = 0.0d0
                xin(lonIndex,latIndex) = 1.0d0
                if(limplicit) then
                   do t = 1, diff(diffID)%numt
                      call diffusion1x_implicit(diffID,xin,xin)
                      call diffusion1y_implicit(diffID,xin,xin)
                   end do
                else
                   call diffusion_explicit(diffID,xin,xin)
                end if
                xin = diff(diffID)%Winvsqrt*xin
                diff(diffID)%Lambda(lonIndex,latIndex) = 0.0d0
                do l = 1, NJ
                   do k = 1, NI
                      diff(diffID)%Lambda(lonIndex,latIndex) = diff(diffID)%Lambda(lonIndex,latIndex) + xin(k,l)*xin(k,l)
                   end do
                end do
                diff(diffID)%Lambda(lonIndex,latIndex) = sqrt(diff(diffID)%Lambda(lonIndex,latIndex))
                if(diff(diffID)%Lambda(lonIndex,latIndex) > 0.0d0) then
                   diff(diffID)%Lambda(lonIndex,latIndex) = 1.0d0/diff(diffID)%Lambda(lonIndex,latIndex);
                end if

             end do
          end do

       end if

       npak = 0
       dateo = 0
       deet = 0
       npas = 0
       ip1 = 0
       ip2 = 0
       ip3 = 0
       typvar = 'X'
       grtyp = 'X'
       ig1 = 0
       ig2 = 0
       ig3 = 0
       ig4 = 0
       datyp = 1
       rewrit = .FALSE.

       ierr = fnom(std_unit, diff_norm_fact, 'RND', 0)
       nmax = fstouv(std_unit, 'RND')

!         ierr = fstecr(real(khalfx), dumwrk, npak, std_unit,
!        X        dateo, deet, npas,
!        X        NI-1, NJ, 1, ip1, ip2, ip3,
!        X        typvar, 'KHX', etiket,
!        X        grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)
!   
!         ierr = fstecr(real(khalfy), dumwrk, npak, std_unit,
!        X        dateo, deet, npas,
!        X        NI, NJ-1, 1, ip1, ip2, ip3,
!        X        typvar, 'KHY', etiket,
!        X        grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)
!   
!         ierr = fstecr(real(mhalfx), dumwrk, npak, std_unit,
!        X        dateo, deet, npas,
!        X        NI-1, NJ, 1, ip1, ip2, ip3,
!        X        typvar, 'MHX', etiket,
!        X        grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)
!   
!         ierr = fstecr(real(mhalfy), dumwrk, npak, std_unit,
!        X        dateo, deet, npas,
!        X        NI, NJ-1, 1, ip1, ip2, ip3,
!        X        typvar, 'MHY', etiket,
!        X        grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)
!   
       ierr = fstecr(real(diff(diffID)%Lambda), dumwrk, npak, std_unit,  &
            dateo, deet, npas,                              &
            NI, NJ, 1, ip1, ip2, ip3,                       &
            typvar, 'LAMB', etiket,                         &
            grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)

!         ierr = fstecr(real(Winv), dumwrk, npak, std_unit,
!        X        dateo, deet, npas,
!        X        NI, NJ, 1, ip1, ip2, ip3,
!        X        typvar, 'WINV', etiket,
!        X        grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)
!   
!         ierr = fstecr(real(Wsqrt), dumwrk, npak, std_unit,
!        X        dateo, deet, npas,
!        X        NI, NJ, 1, ip1, ip2, ip3,
!        X        typvar, 'WSQR', etiket,
!        X        grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)
!   
!         ierr = fstecr(real(Winvsqrt), dumwrk, npak, std_unit,
!        X        dateo, deet, npas,
!        X        NI, NJ, 1, ip1, ip2, ip3,
!        X        typvar, 'WISQ', etiket,
!        X        grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)

       ierr = fstfrm(std_unit)
       ierr = fclos(std_unit)

    end if

    deallocate(xin)
    deallocate(m)
    deallocate(W)
    deallocate(kappa)
    deallocate(Lcorr)
    deallocate(latr)

    diff_setup = diffID

  end function diff_setup


  subroutine diff_finalize(diffID)
    !
    !:Purpose: To finalize the diffusion operator module, and to free up memory.
    implicit none

    ! Arguments:
    integer, intent(in)  :: diffID

    deallocate(diff(diffID)%diff1y_c)
    deallocate(diff(diffID)%diff1y_bp_inv)
    deallocate(diff(diffID)%diff1y_ap)
    deallocate(diff(diffID)%diff1x_c)
    deallocate(diff(diffID)%diff1x_bp_inv)
    deallocate(diff(diffID)%diff1x_ap)

    deallocate(diff(diffID)%Lambda)
    deallocate(diff(diffID)%mhalfy, diff(diffID)%mhalfx)
    deallocate(diff(diffID)%khalfy, diff(diffID)%khalfx)
    deallocate(diff(diffID)%Winvsqrt, diff(diffID)%Wsqrt, diff(diffID)%Winv)
    deallocate(diff(diffID)%cosyinvsq, diff(diffID)%cosyinv, diff(diffID)%cosyhalf)

  end subroutine diff_finalize


  subroutine diffusion_explicit(diffID, xin, xout)
    !
    !:Purpose: To compute Lsqrt*xin (diffusion over numt/2 timesteps), and to
    !          specify initial conditions
    implicit none

    ! Arguments:
    integer, intent(in)  :: diffID
    real(8), intent(in)  :: xin(diff(diffID)%ni,diff(diffID)%nj)
    real(8), intent(out) :: xout(diff(diffID)%ni,diff(diffID)%nj)

    ! Locals:
    integer :: t, latIndex, lonIndex

    real(8) :: xlast(diff(diffID)%ni,diff(diffID)%nj)

    call tmg_start(192,'diffusion_explicit')

    xlast(:,:) = xin(:,:)
    ! iterate difference equations
    do t = 1, diff(diffID)%numt/2
       !$OMP PARALLEL DO PRIVATE(latIndex,lonIndex)
       do latIndex = 2, diff(diffID)%nj-1
          do lonIndex = 2, diff(diffID)%ni-1
             xout(lonIndex,latIndex) = xlast(lonIndex,latIndex) + diff(diffID)%dt*(                           &
                  diff(diffID)%cosyinvsq(latIndex) * ( diff(diffID)%mhalfx(lonIndex  ,latIndex) * diff(diffID)%khalfx(lonIndex  ,latIndex)*    &
                  ( xlast(lonIndex+1,latIndex) - xlast(lonIndex  ,latIndex) ) / diff(diffID)%dlon -           &
                  diff(diffID)%mhalfx(lonIndex-1,latIndex) * diff(diffID)%khalfx(lonIndex-1,latIndex)*                     &
                  ( xlast(lonIndex  ,latIndex) - xlast(lonIndex-1,latIndex) ) / diff(diffID)%dlon ) / diff(diffID)%dlon    &
                  + diff(diffID)%cosyinv(latIndex) *                                     &
                  ( diff(diffID)%cosyhalf(latIndex  ) * diff(diffID)%mhalfy(lonIndex,latIndex  ) * diff(diffID)%khalfy(lonIndex,latIndex  )*   &
                  ( xlast(lonIndex,latIndex+1) - xlast(lonIndex,latIndex  ) ) / diff(diffID)%dlat -           &
                  diff(diffID)%cosyhalf(latIndex-1) * diff(diffID)%mhalfy(lonIndex,latIndex-1) * diff(diffID)%khalfy(lonIndex,latIndex-1)*   &
                  ( xlast(lonIndex,latIndex  ) - xlast(lonIndex,latIndex-1) ) / diff(diffID)%dlat             &
                  ) / diff(diffID)%dlat                                    &
                  )
          end do
       end do
       !$OMP END PARALLEL DO
       xlast(:,:) = xout(:,:)
!            xlast(diff(diffID)%ni,:) = 0.0d0
!            xlast(:,diff(diffID)%nj) = 0.0d0
    end do
    xout(:,:) = xlast(:,:)

    call tmg_stop(192)

  end subroutine diffusion_explicit


  subroutine diff_Csqrt(diffID, xin, xout)

    implicit none

    ! Arguments
    integer, intent(in)  :: diffID
    real(8), intent(in)  :: xin(diff(diffID)%ni,diff(diffID)%nj)
    real(8), intent(out) :: xout(diff(diffID)%ni,diff(diffID)%nj)

    ! Locals:
    integer :: t

    ! compute Csqrt

    ! this is the C^1/2 required for the forward model: Csqrt = Lambda * Diffuse * W^-1/2
    xout(:,:) = diff(diffID)%Winvsqrt(:,:) * xin(:,:)
    if (diff(diffID)%limplicit) then
       do t = 1, diff(diffID)%numt
          call diffusion1x_implicit(diffID,xout,xout)
          call diffusion1y_implicit(diffID,xout,xout)
       end do
    else
       call diffusion_explicit(diffID,xout,xout)
    end if
    xout(:,:) = diff(diffID)%Lambda(:,:) * xout(:,:)

  end subroutine diff_Csqrt


  subroutine diff_Csqrtadj(diffID, xin, xout)

    implicit none

    ! Arguments:
    integer, intent(in)  :: diffID
    real(8), intent(in)  :: xin(diff(diffID)%ni,diff(diffID)%nj)
    real(8), intent(out) :: xout(diff(diffID)%ni,diff(diffID)%nj)

    ! Locals:
    integer :: t

    ! compute Csqrtadj

    ! this is the (C^1/2)^T required for the adjoint: Csqrt^T = W^1/2 * Diffuse * W^-1 * Lambda
    xout(:,:) = diff(diffID)%Lambda(:,:) * xin(:,:)
    xout(:,:) = diff(diffID)%Winv(:,:) * xout(:,:)
    if(diff(diffID)%limplicit) then
       do t = 1, diff(diffID)%numt
          call diffusion1y_implicit(diffID,xout,xout)
          call diffusion1x_implicit(diffID,xout,xout)
       end do
    else
       call diffusion_explicit(diffID,xout,xout)
    end if
    xout(:,:) = diff(diffID)%Wsqrt(:,:) * xout(:,:)

  end subroutine diff_Csqrtadj


  subroutine diffusion1x_implicit(diffID, xin, xout)
    !
    !:Purpose: To compute Lsqrt*xin (diffusion over 1 timestep, loop over
    !          timesteps is external to the subroutine), and to specify initial
    !          conditions
    implicit none

    ! Arguments:
    integer, intent(in)  :: diffID
    real(8), intent(in)  :: xin(diff(diffID)%ni,diff(diffID)%nj)
    real(8), intent(out) :: xout(diff(diffID)%ni,diff(diffID)%nj)

    ! Locals:
    integer :: t, latIndex, lonIndex

    real(8) :: xlast(diff(diffID)%ni,diff(diffID)%nj)
    real(8) :: dp(diff(diffID)%ni)

    !$OMP PARALLEL DO PRIVATE(latIndex,lonIndex)
    do latIndex = 1, diff(diffID)%nj
       do lonIndex = 1, diff(diffID)%ni
          xlast(lonIndex,latIndex) = xin(lonIndex,latIndex)
       end do
    end do
    !$OMP END PARALLEL DO

    do t = 1, 1
       !$OMP PARALLEL DO PRIVATE(latIndex,lonIndex,dp)
       do latIndex = 2, diff(diffID)%nj-1
          lonIndex = 2
          dp(lonIndex) = xlast(lonIndex,latIndex)
          do lonIndex = 3,diff(diffID)%ni-1
             dp(lonIndex) = xlast(lonIndex,latIndex) - diff(diffID)%diff1x_ap(lonIndex,latIndex)*dp(lonIndex-1)
          end do
          lonIndex = diff(diffID)%ni-1
          xout(lonIndex,latIndex) = dp(lonIndex)*diff(diffID)%diff1x_bp_inv(lonIndex,latIndex)
          do lonIndex = diff(diffID)%ni-2,2,-1
             xout(lonIndex,latIndex) = (dp(lonIndex)-diff(diffID)%diff1x_c(lonIndex,latIndex)*xout(lonIndex+1,latIndex)) &
                  *diff(diffID)%diff1x_bp_inv(lonIndex,latIndex)
          end do
       end do
       !$OMP END PARALLEL DO
    end do  ! do t

    do latIndex = 1, diff(diffID)%nj
       xout(1,latIndex) = xin(1,latIndex)
       xout(diff(diffID)%ni,latIndex) = xin(diff(diffID)%ni,latIndex)
    end do

  end subroutine diffusion1x_implicit


  subroutine diffusion1y_implicit(diffID, xin, xout)
    !
    !:Purpose: To compute Lsqrt*xin (diffusion over 1 timestep, loop over 
    !          timesteps is external to the subroutine), and to specify initial
    !          conditions
    implicit none

    ! Arguments:
    integer, intent(in)  :: diffID
    real(8), intent(in)  :: xin(diff(diffID)%ni,diff(diffID)%nj)
    real(8), intent(out) :: xout(diff(diffID)%ni,diff(diffID)%nj)

    ! Locals:
    integer :: t, latIndex, lonIndex

    real(8) :: xlast(diff(diffID)%nj,diff(diffID)%ni)
    real(8) :: dp(diff(diffID)%nj)

! NOTE:for improved efficiency, the 2D fields used internally are !
!      ordered (diff(diffID)%nj,diff(diffID)%ni) and NOT (diff(diffID)%ni,diff(diffID)%nj) as in the rest of the code !

    !$OMP PARALLEL DO PRIVATE(latIndex,lonIndex)
    do latIndex = 1, diff(diffID)%nj
       do lonIndex = 1, diff(diffID)%ni
          xlast(latIndex,lonIndex) = xin(lonIndex,latIndex)
       end do
    end do
    !$OMP END PARALLEL DO

    do t = 1, 1
       !$OMP PARALLEL DO PRIVATE(latIndex,lonIndex,dp)
       do lonIndex = 2, diff(diffID)%ni-1
          latIndex = 2
          dp(latIndex) = xlast(latIndex,lonIndex)
          do latIndex = 3, diff(diffID)%nj-1
             dp(latIndex) = xlast(latIndex,lonIndex) - diff(diffID)%diff1y_ap(latIndex,lonIndex)*dp(latIndex-1)
          end do
          latIndex = diff(diffID)%nj-1
          xout(lonIndex,latIndex) = dp(latIndex)*diff(diffID)%diff1y_bp_inv(latIndex,lonIndex)
          do latIndex = diff(diffID)%nj-2, 2, -1
             xout(lonIndex,latIndex) = (dp(latIndex)-diff(diffID)%diff1y_c(latIndex,lonIndex)*xout(lonIndex,latIndex+1)) &
                  *diff(diffID)%diff1y_bp_inv(latIndex,lonIndex)
          end do
       end do
       !$OMP END PARALLEL DO
    end do  ! do t

    do lonIndex = 1, diff(diffID)%ni
       xout(lonIndex,1) = xin(lonIndex,1)
       xout(lonIndex,diff(diffID)%nj) = xin(lonIndex,diff(diffID)%nj)
    end do

  end subroutine diffusion1y_implicit


  subroutine diffusion(diffID, xin, xout)

    implicit none

    ! Arguments:
    integer, intent(in)  :: diffID
    real(8), intent(in)  :: xin(diff(diffID)%ni,diff(diffID)%nj)
    real(8), intent(out) :: xout(diff(diffID)%ni,diff(diffID)%nj)

    ! Locals:
    integer :: t

    xout(:,:) = xin(:,:)
    if (diff(diffID)%limplicit) then
       do t = 1, diff(diffID)%numt
          call diffusion1x_implicit(diffID,xout,xout)
          call diffusion1y_implicit(diffID,xout,xout)
       end do
    else
       call diffusion_explicit(diffID,xout,xout)
    end if

  end subroutine diffusion

end module diffusion_mod

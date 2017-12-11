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

!--------------------------------------------------------------------------
!! MODULE diffusion (prefix="diff")
!!
!! *Purpose*: 
!! Storage for data required by the diffusion operator used to model
!! background-error horizontal correlations.
!!
!! Reference:
!! Weaver, A. T., and P. Courtier, 2001: Correlation modelling on the
!! sphere using a generalized diffusion equation.
!! Q. J. R. Meteorol. Soc., 127, 1815-1846.
!!
!! Basic equations:
!! Lcorr^2 = 2*k*dt*numt   (1)
!! stab    = k*dt/dx^2     (2)
!!
!! Author: Alain Caya (alain.caya@canada.ca)
!!
!! *Origin*: diffus_mod.f from the 3D-Var sea ice analysis system
!!
!--------------------------------------------------------------------------
module diffusion_mod

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
  end type struct_diff

  integer, parameter :: nMaxDiff = 10
  integer            :: nDiffAlreadyAllocated = 0
  type(struct_diff)  :: diff(nMaxDiff)

!*************************************************************************

contains

  integer function diff_setup(hco, corr_len, stab, nsamp)

    use horizontalCoord_mod
    use EarthConstants_mod, only : rayt
    use randomNumber_mod
    use utilities_mod

    implicit none

    type(struct_hco), pointer :: hco

    ! Horizontal correlation length scale (km)
    real, intent(in) :: corr_len

    ! Stability criteria (definitely < 0.5)
    real, intent(in) :: stab

    ! Number of samples in the estimation of the normalization factors by randomization.
    integer, intent(in) :: nsamp


    ! Local Variables

    ! latitudes on the analysis rotated grid, in radians
    real(8), allocatable :: latr(:)

    real, allocatable :: buf2d(:,:)
    integer, allocatable :: mask(:,:)

    integer :: i, j, isamp, k, l

    real(8) :: mindxy, maxL

    real(8), allocatable :: Lcorr(:,:)
    real(8), allocatable :: kappa(:,:)
    real(8), allocatable :: W(:,:)
    real(8), allocatable :: m(:,:)
    real(8), allocatable :: xin(:,:)

    ! diff_norm_fact is the name of the RPN format file for the normalization factors.
    character(len=*), parameter :: diff_norm_fact = './diffusmod.std'

    character(len=*), parameter :: mask_file = './analysisgrid'

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

    integer :: ni, nj


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
    do j=1,NJ
       do i=1,NI-1
          diff(diffID)%khalfx(i,j) = (kappa(i,j) + kappa(i+1,j))/2.0d0
       end do
    end do
    do j=1,NJ-1
       do i=1,NI
          diff(diffID)%khalfy(i,j) = (kappa(i,j) + kappa(i,j+1))/2.0d0
       end do
    end do

    ! print this stuff in listing file for user information:
    write(*,*)
    write(*,*) 'Number of timesteps = ',diff(diffID)%numt 
    write(*,*) 'Stability           = ',maxval(kappa)*diff(diffID)%dt/(mindxy**2)
    write(*,*)

    ! this is the matrix necessary for defining the inner product: for lat-lon grid, only cos(y)
    ! Actually, only the diagonal of the matrix is stored in the array W, since the matrix is diagonal.
    W(1,:) = cos(latr(:))
    do i=2,NI
       W(i,:) = W(1,:)
    end do
    diff(diffID)%Winv(:,:) = 1.0d0/W(:,:)
    diff(diffID)%Wsqrt(:,:) = sqrt(W(:,:))
    diff(diffID)%Winvsqrt(:,:) = 1.0d0/diff(diffID)%Wsqrt(:,:)

    ! Read mask from file

    std_unit = 0

    inquire (file=mask_file, exist=file_exist)

    if(file_exist) then

       ierr = fnom(std_unit, mask_file, 'RND+R/O', 0)
       nmax = fstouv(std_unit, 'RND')

       allocate(mask(NI,NJ))

       key = fstlir(mask, std_unit, nii, njj, nkk,    &
            -1, ' ', -1, -1, -1, '@@', ' ')

       ierr = fstfrm(std_unit)
       ierr = fclos(std_unit)

    else

       write(*,*) 'diff_setup: File: ',trim(mask_file)
       write(*,*) 'is missing and required for the mask'
       call utl_abort('diff_setup')

    end if

    ! land mask (1=water, 0=land)
    do j=1,NJ
       do i=1,NI
          if(mask(i,j) == 1) then
             m(i,j) = 1.0d0
          else
             m(i,j) = 0.0d0
          end if
       end do
    end do
    m(:,1) = 0.0d0
    m(1,:) = 0.0d0
    m(:,NJ) = 0.0d0
    m(NI,:) = 0.0d0

    ! define mask on staggered grids
    do j=1,NJ
       do i=1,NI-1
          if(sum(m(i:i+1,j)) < 2.0d0) then
             diff(diffID)%mhalfx(i,j) = 0.0d0
          else
             diff(diffID)%mhalfx(i,j) = 1.0d0
          end if
       end do
    end do

    do j=1,NJ-1
       do i=1,NI
          if(sum(m(i,j:j+1)) < 2.0d0) then
             diff(diffID)%mhalfy(i,j) = 0.0d0
          else
             diff(diffID)%mhalfy(i,j) = 1.0d0
          end if
       end do
    end do

    diff(diffID)%ni = ni
    diff(diffID)%nj = nj

    write (etiket, FMT='(''KM'',i3.3,''STAB'',f3.1)') int(corr_len), stab

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
          do isamp=1,nsamp

             do j=1,NJ
                do i=1,NI
                   xin(i,j) = diff(diffID)%Winvsqrt(i,j)*rng_gaussian()
                end do
             end do
             call diffusion(diffID,xin,xin)
             diff(diffID)%Lambda = diff(diffID)%Lambda + xin*xin

             write(*,*) 'Sample ',isamp,' done.'

          end do

          do j=1,NJ
             do i=1,NI
                diff(diffID)%Lambda(i,j) = sqrt(diff(diffID)%Lambda(i,j)/dble(nsamp-1)) ! normalization: inverse of rms of ens
                if(diff(diffID)%Lambda(i,j) > 0.0d0) then
                   diff(diffID)%Lambda(i,j) = 1.0d0/diff(diffID)%Lambda(i,j);
                end if
             end do
          end do

       else

          write(*,*) 'Exact calculation of the '//      &
               'normalization for diffusion...'
          call flush(6)

          do j=1,NJ
             write(*,*) 'Doing row j = ',j,' of ',NJ
             call flush(6)
             do i=1,NI

                xin = 0.0d0
                xin(i,j) = 1.0d0
                call diffusion(diffID,xin,xin)
                xin = diff(diffID)%Winvsqrt*xin
                diff(diffID)%Lambda(i,j) = 0.0d0
                do l=1,NJ
                   do k=1,NI
                      diff(diffID)%Lambda(i,j) = diff(diffID)%Lambda(i,j) + xin(k,l)*xin(k,l)
                   end do
                end do
                diff(diffID)%Lambda(i,j) = sqrt(diff(diffID)%Lambda(i,j))
                if(diff(diffID)%Lambda(i,j) > 0.0d0) then
                   diff(diffID)%Lambda(i,j) = 1.0d0/diff(diffID)%Lambda(i,j);
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
       rewrit=.FALSE.

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

    deallocate(mask)
    deallocate(xin)
    deallocate(m)
    deallocate(W)
    deallocate(kappa)
    deallocate(Lcorr)
    deallocate(latr)

    diff_setup = diffID

  end function diff_setup

!**********************************************************

  subroutine diff_finalize(diffID)

!-----------------------------------------------------------------------
!
! Purpose:
! Finalize the diffusion operator module.
! Free up memory.
!
! Author: Alain Caya (alain.caya@canada.ca)
!
!-----------------------------------------------------------------------

    implicit none

    integer, intent(in)  :: diffID

    deallocate(diff(diffID)%Lambda)
    deallocate(diff(diffID)%mhalfy, diff(diffID)%mhalfx)
    deallocate(diff(diffID)%khalfy, diff(diffID)%khalfx)
    deallocate(diff(diffID)%Winvsqrt, diff(diffID)%Wsqrt, diff(diffID)%Winv)
    deallocate(diff(diffID)%cosyinvsq, diff(diffID)%cosyinv, diff(diffID)%cosyhalf)

  end subroutine diff_finalize

!**********************************************************

  subroutine diffusion(diffID, xin, xout)

!-----------------------------------------------------------------------
!
! Purpose:
! compute Lsqrt*xin (diffusion over numt/2 timesteps)
! specify initial conditions
!
! Author: Alain Caya (alain.caya@canada.ca)
!
!-----------------------------------------------------------------------

    implicit none

    integer, intent(in)  :: diffID
    real(8), intent(in)  :: xin(diff(diffID)%ni,diff(diffID)%nj)
    real(8), intent(out) :: xout(diff(diffID)%ni,diff(diffID)%nj)

    integer :: t, j, i

    real(8) :: xlast(diff(diffID)%ni,diff(diffID)%nj)

    call tmg_start(2,'diffusion')

    xlast(:,:) = xin(:,:)
    ! iterate difference equations
    do t=1,diff(diffID)%numt/2
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(j,i)
       do j=2,diff(diffID)%nj-1
          do i=2,diff(diffID)%ni-1
             xout(i,j) = xlast(i,j) + diff(diffID)%dt*(                           &
                  diff(diffID)%cosyinvsq(j) * ( diff(diffID)%mhalfx(i  ,j) * diff(diffID)%khalfx(i  ,j)*    &
                  ( xlast(i+1,j) - xlast(i  ,j) ) / diff(diffID)%dlon -           &
                  diff(diffID)%mhalfx(i-1,j) * diff(diffID)%khalfx(i-1,j)*                     &
                  ( xlast(i  ,j) - xlast(i-1,j) ) / diff(diffID)%dlon ) / diff(diffID)%dlon    &
                  + diff(diffID)%cosyinv(j) *                                     &
                  ( diff(diffID)%cosyhalf(j  ) * diff(diffID)%mhalfy(i,j  ) * diff(diffID)%khalfy(i,j  )*   &
                  ( xlast(i,j+1) - xlast(i,j  ) ) / diff(diffID)%dlat -           &
                  diff(diffID)%cosyhalf(j-1) * diff(diffID)%mhalfy(i,j-1) * diff(diffID)%khalfy(i,j-1)*   &
                  ( xlast(i,j  ) - xlast(i,j-1) ) / diff(diffID)%dlat             &
                  ) / diff(diffID)%dlat                                    &
                  )
          end do
       end do
!$OMP END DO
!$OMP END PARALLEL
       xlast(:,:) = xout(:,:)
!            xlast(diff(diffID)%ni,:) = 0.0d0
!            xlast(:,diff(diffID)%nj) = 0.0d0
    end do
    xout(:,:) = xlast(:,:)

    call tmg_stop(2)

  end subroutine diffusion

!**********************************************************

  subroutine diff_Csqrt(diffID, xin, xout)

    implicit none

    integer, intent(in)  :: diffID
    real(8), intent(in)  :: xin(diff(diffID)%ni,diff(diffID)%nj)
    real(8), intent(out) :: xout(diff(diffID)%ni,diff(diffID)%nj)

    ! compute Csqrt

    ! this is the C^1/2 required for the forward model: Csqrt = Lambda * Diffuse * W^-1/2
    xout(:,:) = diff(diffID)%Winvsqrt(:,:) * xin(:,:)
    call diffusion(diffID,xout,xout)
    xout(:,:) = diff(diffID)%Lambda(:,:) * xout(:,:)

  end subroutine diff_Csqrt

!**********************************************************

  subroutine diff_Csqrtadj(diffID, xin, xout)

    implicit none

    integer, intent(in)  :: diffID
    real(8), intent(in)  :: xin(diff(diffID)%ni,diff(diffID)%nj)
    real(8), intent(out) :: xout(diff(diffID)%ni,diff(diffID)%nj)

    ! compute Csqrtadj

    ! this is the (C^1/2)^T required for the adjoint: Csqrt^T = W^1/2 * Diffuse * W^-1 * Lambda
    xout(:,:) = diff(diffID)%Lambda(:,:) * xin(:,:)
    xout(:,:) = diff(diffID)%Winv(:,:) * xout(:,:)
    call diffusion(diffID,xout,xout)
    xout(:,:) = diff(diffID)%Wsqrt(:,:) * xout(:,:)

  end subroutine diff_Csqrtadj

end module diffusion_mod

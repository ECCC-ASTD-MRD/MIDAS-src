      module diffus_mod

! <next 7 lines automatically updated by CVS, do not edit>
! $RCSfile: diffus_mod.f,v $
! $Revision: 115 $
! $Date: 2012-10-26 14:31:56 -0400 (Fri, 26 Oct 2012) $
! $Author: armagr5 $
! $State: Exp $
! $Source: /users/dor/arma/gr5/.ocm/.cvs/I1/src/3dvar/diffus_mod.f,v $
! $Name:  $

!-----------------------------------------------------------------------
!
! Purpose:
! Storage for data required by the diffusion operator used to model
! background-error horizontal correlations.
!
! Reference:
! Weaver, A. T., and P. Courtier, 2001: Correlation modelling on the
! sphere using a generalized diffusion equation.
! Q. J. R. Meteorol. Soc., 127, 1815-1846.
!
! Basic equations:
! Lcorr^2 = 2*k*dt*numt   (1)
! stab    = k*dt/dx^2     (2)
!
! Author: Alain Caya (alain.caya@ec.gc.ca)
!
!-----------------------------------------------------------------------

      use analysis_mod, only : NI, NJ

      implicit none
      private

      ! Public subroutines and functions

      public :: sudiff, fin_diff, Csqrt, Csqrtadj, corr_diff, 
     &          diffus_data, loc_diff, test_diff

      double precision :: dlon, dlat        ! grid spacing in radians

      double precision, allocatable :: cosyhalf(:), cosyinv(:),
     $                                 cosyinvsq(:), cosy(:)
      double precision, allocatable :: Winv(:,:), Wsqrt(:,:),
     $                                 Winvsqrt(:,:)
      double precision, allocatable :: mhalfx(:,:), mhalfy(:,:)
      real*8, allocatable :: diff1x_ap(:,:),diff1x_bp_inv(:,:)
      real*8, allocatable :: diff1x_c(:,:)
      real*8, allocatable :: diff1y_ap(:,:),diff1y_bp_inv(:,:)
      real*8, allocatable :: diff1y_c(:,:)

      integer, parameter :: niter_imp = 4

      type diffus_data
        logical :: limplicit
        real :: corr_len
        integer :: numt
        double precision :: dt
        double precision, allocatable :: khalfx(:,:), khalfy(:,:)
        double precision, allocatable :: Lambda(:,:)
      end type diffus_data
!     NI          ! number of grid points in x (includes perimeter of land points)
!     NJ          ! number of grid points in y (includes perimeter of land points)
!     grd_dx      ! grid spacing in degrees
!     grd_dy
!     latr        ! latitudes on the rotated analysis grid, in radians

!     Namelist with default values

!     Horizontal correlation length scale (km)
!      real :: corr_len = 0.0
!     Stability criteria (definitely < 0.5)
      real :: stab     = 0.3
!     Number of samples in the estimation of the normalization factors by randomization.
      integer :: nsamp = 23241

      type(diffus_data) :: corr_diff, loc_diff

      namelist /diffusion_nml/ stab, nsamp

!*************************************************************************

      contains

      subroutine sudiff(paramlist, diff_norm_fact, corr_len, limplicit, diff)

      use            const_mod, only : earth_radius, deg2rad
      use         analysis_mod, only : grd_dx, grd_dy,
     $                                 Grd_jref, Grd_latr, mask
!     $     , latr

      implicit none

      character(len=500), intent(in) :: paramlist, diff_norm_fact
!     diff_norm_fact is the name of the RPN format file for the normalization factors.
      real, intent(in) :: corr_len
      type(diffus_data), intent(out) :: diff

!     Local Variables

      integer :: i, j, isamp, k, l
      integer :: iunit

      double precision :: mindxy, maxL

      double precision :: Lcorr(NI,NJ)
      double precision :: kappa(NI,NJ)
      double precision :: latr(NJ)
      double precision :: W(NI,NJ)
      double precision :: m(NI,NJ)
      double precision :: xin(NI,NJ)

!     Function
      double precision :: gasdev

      character*10 :: date,time

!     Variables and functions required to write to RPN Standard files.

      integer  :: nmax

      integer  :: std_unit, ierr, key, npak, datyp
      integer  :: nii, njj, nkk
      integer  :: dateo
      integer  :: deet,npas
      integer  :: ip1,ip2,ip3
      integer  :: ig1, ig2, ig3, ig4
      character*1 :: grtyp
      character*2 :: typvar
      character*12 :: etiket
      logical  :: rewrit, file_exist
      real     :: dumwrk(1)

      real, allocatable :: buf2d(:,:)

      integer  :: fnom,fstouv,fstecr,fstlir,fstfrm,fclos
      external :: fnom,fstouv,fstecr,fstlir,fstfrm,fclos

      write(*,*)
      write(*,*)
     $     '***********************************************************'
      write(*,*)
     $     '** Module diffus_mod setting up...                       **'
      write(*,*)
     $     '** $Revision: 115 $ last vers. in which the file changed **'
      write(*,*)
     $     '** $Date: 2012-10-26 14:31:56 -0400 (Fri, 26 Oct 2012) $ **'
      write(*,*)
     $     '***********************************************************'
      write(*,*)

!     Begin by reading the namelist inputs

      diff%limplicit = limplicit
      iunit = 0
      ierr = fnom(iunit, paramlist, 'OLD+SEQ+R/O', 0)
      read(iunit, nml = diffusion_nml)
      ierr = fclos(iunit)

      write(*, nml = diffusion_nml)

      diff%corr_len = corr_len
!     Export the length scale factor

      allocate(cosyhalf(NJ), cosyinv(NJ), cosyinvsq(NJ))
      allocate(Winv(NI,NJ), Wsqrt(NI,NJ), Winvsqrt(NI,NJ))
      allocate(diff%khalfx(NI-1,NJ), diff%khalfy(NI,NJ-1))
      allocate(mhalfx(NI-1,NJ), mhalfy(NI,NJ-1))
      allocate(diff%Lambda(NI,NJ))

      allocate(diff1x_ap(NI,NJ))
      allocate(diff1x_bp_inv(NI,NJ))
      allocate(diff1x_c(NI,NJ))
      allocate(diff1y_ap(NJ,NI))
      allocate(diff1y_bp_inv(NJ,NI))
      allocate(diff1y_c(NJ,NI))

      dlon = grd_dx*deg2rad
      dlat = grd_dy*deg2rad

      ! Could also use latr directly from analysis_mod.
      ! That would avoid the calculation and the need of Grd_latr and Grd_jref
      do j=1,NJ
         latr(j) = (Grd_latr + dble(j-Grd_jref)*grd_dy)*deg2rad
      enddo

      cosyinv(:) = 1.0d0/cos(latr(:))
      cosyinvsq(:) = cosyinv(:)*cosyinv(:)

!     cosinus of latitudes on staggered grid
      cosyhalf(:) = cos(latr(:) + 0.5d0*dlat)

!     minimum grid spacing over the domain (radians, but with cos(latr))
      mindxy = min(minval(cos(latr(:)))*dlon,dlat)

      Lcorr(:,:) = diff%corr_len/earth_radius  ! lengthscale in radians
      maxL = maxval(Lcorr(2:NI-1,2:NJ-1)) ! maximum lengthscale over domain

!     set main parameters for diffusion operator
      kappa = Lcorr**2          ! arbitrarily set k to L^2 (in radians)
      diff%dt = stab*(mindxy**2)/(maxL**2) ! determine dt from stability criteria (2)
c$$$      numt = 1.0d0/(2.0d0*dt)       ! determine number of timesteps from (1)
      diff%numt = ceiling(1.0d0/(4.0d0*diff%dt))*2; ! make sure it is an even integer
      diff%dt = 1.0d0/(2.0d0*dble(diff%numt)); ! recompute dt

      diff%numt = 10
      diff%dt = 0.5 /dble(diff%numt * niter_imp)

!     interpolate diffusion coefficient onto 2 staggered lat-lon grids
      do j=1,NJ
         do i=1,NI-1
            diff%khalfx(i,j) = (kappa(i,j) + kappa(i+1,j))/2.0d0
         enddo
      enddo
      do j=1,NJ-1
         do i=1,NI
            diff%khalfy(i,j) = (kappa(i,j) + kappa(i,j+1))/2.0d0
         enddo
      enddo

!     print this stuff in listing file for user information:
      write(*,*)
      write(*,*) 'Number of timesteps = ',diff%numt 
      write(*,*) 'Stability           = ',maxval(kappa)*
     &                                    diff%dt/(mindxy**2)
      write(*,*)

!     this is the matrix necessary for defining the inner product: for lat-lon grid, only cos(y)
!     Actually, only the diagonal of the matrix is stored in the array W, since the matrix is diagonal.
      W(1,:) = cos(latr(:))
      do i=2,NI
         W(i,:) = W(1,:)
      enddo
      Winv = 1.0d0/W
      Wsqrt = sqrt(W)
      Winvsqrt = 1.0d0/Wsqrt

! land mask (1=water, 0=land)
      do j=1,NJ
         do i=1,NI
            if(mask(i,j) == 0.0) then
               m(i,j) = 1.0d0
            else
               m(i,j) = 0.0d0
            endif
         enddo
      enddo
      m(:,1) = 0.0d0
      m(1,:) = 0.0d0
      m(:,NJ) = 0.0d0
      m(NI,:) = 0.0d0

! define mask on staggered grids
      do j=1,NJ
         do i=1,NI-1
            if(sum(m(i:i+1,j)) < 2.0d0) then
               mhalfx(i,j) = 0.0d0
            else
               mhalfx(i,j) = 1.0d0
            endif
         enddo
      enddo

      do j=1,NJ-1
         do i=1,NI
            if(sum(m(i,j:j+1)) < 2.0d0) then
               mhalfy(i,j) = 0.0d0
            else
               mhalfy(i,j) = 1.0d0
            endif
         enddo
      enddo

! specify number of timesteps and timestep length for implicit 1D diffusion
      if (diff%limplicit) then
        diff%numt=5
        diff%dt=1.0d0/(2.0d0*dble(2*diff%numt)-3.0d0)
      endif

! compute the LU decomposition for the implicit 1D diffusion
      diff1x_ap(:,:)=0.0d0
      diff1x_bp_inv(:,:)=0.0d0      
      diff1x_c(:,:)=0.0d0
      diff1y_ap(:,:)=0.0d0
      diff1y_bp_inv(:,:)=0.0d0      
      diff1y_c(:,:)=0.0d0

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(j,i,a,b)
      do j=2,NJ-1
         i=2
         diff1x_bp_inv(i,j)= 1.0d0/(1.0d0 + 
     +        dt2*cosyinvsq(j)*(mhalfx(i,j)*diff%khalfx(i,j) + 
     +        mhalfx(i-1,j)*diff%khalfx(i-1,j))/(dlon*dlon))
         do i=3,NI-1
            ! elements of the tri-diagonal coefficient matrix
            a=   - dt2*cosyinvsq(j)*mhalfx(i-1,j)*diff%khalfx(i-1,j)
     +                /(dlon*dlon)
            b=    1 + dt2*cosyinvsq(j)*(mhalfx(i,j)*diff%khalfx(i,j) + 
     +                mhalfx(i-1,j)*diff%khalfx(i-1,j))/(dlon*dlon)
            diff1x_c(i,j)= - dt2*cosyinvsq(j)*mhalfx(i,j)*diff%khalfx(i,j)
     +                /(dlon*dlon)
            diff1x_ap(i,j)=a*diff1x_bp_inv(i-1,j)
            diff1x_bp_inv(i,j)=1.0d0/(b-a*diff1x_c(i-1,j)
     +                         *diff1x_bp_inv(i-1,j))
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(j,i,a,b)
      do i=2,NI-1
         j=2
         diff1y_bp_inv(j,i)= 1.0d0/(1.0d0 + 
     +     dt2*cosyinvsq(j)*(cosyhalf(j)*mhalfy(i,j)*diff%khalfy(i,j) + 
     +     cosyhalf(j-1)*mhalfy(i,j-1)*diff%khalfy(i,j-1))/(dlat*dlat))
         do j=3,NJ-1
            ! elements of the tri-diagonal coefficient matrix
            a= - dt2*cosyinv(j)*cosyhalf(j-1)*
     +           mhalfy(i,j-1)*diff%khalfy(i,j-1)/(dlat*dlat)
            b= 1 + dt2*cosyinv(j)*(cosyhalf(j)*mhalfy(i,j)*diff%khalfy(i,j) + 
     +         cosyhalf(j-1)*mhalfy(i,j-1)*diff%khalfy(i,j-1))/(dlat*dlat)
            diff1y_c(j,i)= - dt2*cosyinv(j)*cosyhalf(j)*
     +           mhalfy(i,j)*diff%khalfy(i,j)/(dlat*dlat)
            diff1y_ap(j,i)=a*diff1y_bp_inv(j-1,i)
            diff1y_bp_inv(j,i)=1.0d0/(b-a*diff1y_c(j-1,i)
     +                         *diff1y_bp_inv(j-1,i))
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL


      write (etiket, FMT='(''KM'',i3.3,''STAB'',f3.1)')
     $     int(diff%corr_len), stab

      std_unit = 0

      inquire (file=diff_norm_fact, exist=file_exist)

      if(file_exist) then

         ierr = fnom(std_unit, diff_norm_fact, 'RND+R/O', 0)
         nmax = fstouv(std_unit, 'RND')

         allocate(buf2d(NI,NJ))

         key = fstlir(buf2d, std_unit, nii, njj, nkk,
     X        -1, etiket, -1, -1, -1, ' ', 'LAMB')

         ierr = fstfrm(std_unit)
         ierr = fclos(std_unit)

         diff%Lambda = dble(buf2d)

         deallocate(buf2d)

      endif

      if(nii /= ni .or. njj /= nj .or. key <= 0 .or.
     $     (.not. file_exist)) then

         call DATE_AND_TIME(date, time)
         write(*,*) 'date is: ',date
         write(*,*) 'time is: ',time

         if(nsamp < NI*NJ) then

! compute normalization:  Lambda = inverse stddev of (Diffuse * W^-1/2)
            write(*,*) 'Randomization estimation of the '//
     $           'normalization for diffusion...'
            write(*,*) 'Will use ',nsamp,' samples.'
            call flush(6)
            diff%Lambda = 0.0d0
            do isamp=1,nsamp

               do j=1,NJ
                  do i=1,NI
                     xin(i,j) = Winvsqrt(i,j)*gasdev(1)
                  enddo
               enddo
               if(diff%limplicit) then
                  do t=1,numt2
                     call diffusion1x_implicit(xin,xin,diff)
                     call diffusion1y_implicit(xin,xin,diff)
                  enddo
               else
                 call diffusion(xin,xin,diff)
               endif
               diff%Lambda = diff%Lambda + xin*xin

               call DATE_AND_TIME(date, time)
               write(*,*) 'Sample ',isamp,' done. Time is: ',time
               call flush(6)

            enddo

            do j=1,NJ
               do i=1,NI
                  diff%Lambda(i,j) = sqrt(diff%Lambda(i,j)/
     &                               dble(nsamp-1)) ! normalization: inverse of rms of ens
                  if(diff%Lambda(i,j) > 0.0d0) then
                     diff%Lambda(i,j) = 1.0d0/diff%Lambda(i,j);
                  endif
               enddo
            enddo

         else

            write(*,*) 'Exact calculation of the '//
     $           'normalization for diffusion...'
            call flush(6)

            do j=1,NJ
               write(*,*) 'Doing row j = ',j,' of ',NJ
               call flush(6)
               do i=1,NI

                  xin = 0.0d0
                  xin(i,j) = 1.0d0
                  if(diff%limplicit) then
                     do t=1,numt2
                        call diffusion1x_implicit(xin,xin,diff)
                        call diffusion1y_implicit(xin,xin,diff)
                     enddo
                  else
                    call diffusion(xin,xin,diff)
                  endif
                  xin = Winvsqrt*xin
                  diff%Lambda(i,j) = 0.0d0
                  do l=1,NJ
                     do k=1,NI
                        diff%Lambda(i,j) = diff%Lambda(i,j) + 
     &                                     xin(k,l)*xin(k,l)
                     enddo
                  enddo
                  diff%Lambda(i,j) = sqrt(diff%Lambda(i,j))
                  if(diff%Lambda(i,j) > 0.0d0) then
                     diff%Lambda(i,j) = 1.0d0/diff%Lambda(i,j);
                  endif

               enddo
               call DATE_AND_TIME(date, time)
               write(*,*) 'time is: ',time
               call flush(6)
            enddo

         endif

         call DATE_AND_TIME(date, time)
         write(*,*) 'time is: ',time
         call flush(6)

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

c$$$      ierr = fstecr(real(khalfx), dumwrk, npak, std_unit,
c$$$     X        dateo, deet, npas,
c$$$     X        NI-1, NJ, 1, ip1, ip2, ip3,
c$$$     X        typvar, 'KHX', etiket,
c$$$     X        grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)
c$$$
c$$$      ierr = fstecr(real(khalfy), dumwrk, npak, std_unit,
c$$$     X        dateo, deet, npas,
c$$$     X        NI, NJ-1, 1, ip1, ip2, ip3,
c$$$     X        typvar, 'KHY', etiket,
c$$$     X        grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)
c$$$
c$$$      ierr = fstecr(real(mhalfx), dumwrk, npak, std_unit,
c$$$     X        dateo, deet, npas,
c$$$     X        NI-1, NJ, 1, ip1, ip2, ip3,
c$$$     X        typvar, 'MHX', etiket,
c$$$     X        grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)
c$$$
c$$$      ierr = fstecr(real(mhalfy), dumwrk, npak, std_unit,
c$$$     X        dateo, deet, npas,
c$$$     X        NI, NJ-1, 1, ip1, ip2, ip3,
c$$$     X        typvar, 'MHY', etiket,
c$$$     X        grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)
c$$$
         ierr = fstecr(real(diff%Lambda), dumwrk, npak, std_unit,
     X        dateo, deet, npas,
     X        NI, NJ, 1, ip1, ip2, ip3,
     X        typvar, 'LAMB', etiket,
     X        grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)

c$$$      ierr = fstecr(real(Winv), dumwrk, npak, std_unit,
c$$$     X        dateo, deet, npas,
c$$$     X        NI, NJ, 1, ip1, ip2, ip3,
c$$$     X        typvar, 'WINV', etiket,
c$$$     X        grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)
c$$$
c$$$      ierr = fstecr(real(Wsqrt), dumwrk, npak, std_unit,
c$$$     X        dateo, deet, npas,
c$$$     X        NI, NJ, 1, ip1, ip2, ip3,
c$$$     X        typvar, 'WSQR', etiket,
c$$$     X        grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)
c$$$
c$$$      ierr = fstecr(real(Winvsqrt), dumwrk, npak, std_unit,
c$$$     X        dateo, deet, npas,
c$$$     X        NI, NJ, 1, ip1, ip2, ip3,
c$$$     X        typvar, 'WISQ', etiket,
c$$$     X        grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)

         ierr = fstfrm(std_unit)
         ierr = fclos(std_unit)

      endif

      end subroutine sudiff

!**********************************************************

      subroutine fin_diff

!-----------------------------------------------------------------------
!
! Purpose:
! Finalize the diffusion operator module.
! Free up memory.
!
! Author: Alain Caya (alain.caya@ec.gc.ca)
!
!-----------------------------------------------------------------------

      implicit none

      deallocate(cosyhalf, cosyinv, cosyinvsq)
      deallocate(Winv, Wsqrt, Winvsqrt)
      deallocate(mhalfx, mhalfy)

      deallocate(diff1x_ap)
      deallocate(diff1x_bp_inv)
      deallocate(diff1x_c)
      deallocate(diff1y_ap)
      deallocate(diff1y_bp_inv)
      deallocate(diff1y_c)

      if (allocated(corr_diff%khalfx)) deallocate(corr_diff%khalfx)
      if (allocated(corr_diff%khalfy)) deallocate(corr_diff%khalfy)
      if (allocated(corr_diff%Lambda)) deallocate(corr_diff%Lambda)
      if (allocated(loc_diff%khalfx)) deallocate(loc_diff%khalfx)
      if (allocated(loc_diff%khalfy)) deallocate(loc_diff%khalfy)
      if (allocated(loc_diff%Lambda)) deallocate(loc_diff%Lambda)

      end subroutine fin_diff

!**********************************************************

      subroutine diffusion(xout,xin,diff)

!-----------------------------------------------------------------------
!
! Purpose:
! compute Lsqrt*xin (diffusion over numt/2 timesteps)
! specify initial conditions
!
! Author: Alain Caya (alain.caya@ec.gc.ca)
!
!-----------------------------------------------------------------------

      implicit none

      double precision, intent(in)  :: xin(NI,NJ)
      double precision, intent(out) :: xout(NI,NJ)
      type(diffus_data), intent(in) :: diff
      integer :: t, j, i

      double precision :: xlast(NI,NJ)

      xlast = xin
!     iterate difference equations
      do t=1,diff%numt/2
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(j,i)
         do j=2,NJ-1
            do i=2,NI-1
               xout(i,j) = xlast(i,j) + diff%dt*(
     &              cosyinvsq(j) * ( mhalfx(i  ,j) * diff%khalfx(i  ,j)*
     $              ( xlast(i+1,j) - xlast(i  ,j) ) / dlon -
     &              mhalfx(i-1,j) * diff%khalfx(i-1,j)*
     $              ( xlast(i  ,j) - xlast(i-1,j) ) / dlon ) / dlon
     &              + cosyinv(j) *
     $              ( cosyhalf(j  ) * mhalfy(i,j  ) * diff%khalfy(i,j )*
     $              ( xlast(i,j+1) - xlast(i,j  ) ) / dlat -
     &                cosyhalf(j-1) *mhalfy(i,j-1) * diff%khalfy(i,j-1)*
     $              ( xlast(i,j  ) - xlast(i,j-1) ) / dlat
     &                     ) / dlat
     $              )
            enddo
         enddo
!$OMP end DO
!$OMP end PARALLEL
         xlast = xout
c$$$         xlast(NI,:) = 0.0d0
c$$$         xlast(:,NJ) = 0.0d0
      enddo
      xout = xlast

      end subroutine diffusion

!**********************************************************

      subroutine diffusion1x_implicit(xout,xin,diff)

!-----------------------------------------------------------------------
!
! Purpose:
! compute Lsqrt*xin (diffusion over numt2 timesteps)
! specify initial conditions
!
! Author: Mark Buehner (mark.buehner@ec.gc.ca)
!
!-----------------------------------------------------------------------

      implicit none

      double precision, intent(in)  :: xin(NI,NJ)
      double precision, intent(out) :: xout(NI,NJ)

      integer :: t, j, i

      double precision :: xlast(NI,NJ)
      double precision :: dp(NI)

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(j,i)
      do j=1,NJ
         do i=1,NI
           xlast(i,j) = xin(i,j)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

!      do t=1,numt2
      do t=1,1
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(j,i,dp)
         do j=2,NJ-1
            i=2
            dp(i)=xlast(i,j)
            do i=3,NI-1
               dp(i)=xlast(i,j) - diff1x_ap(i,j)*dp(i-1)
            enddo
            i=NI-1
            xout(i,j)=dp(i)*diff1x_bp_inv(i,j)
            do i=NI-2,2,-1
               xout(i,j)=(dp(i)-diff1x_c(i,j)*xout(i+1,j))
     +                   *diff1x_bp_inv(i,j)
            enddo
!            if(t.lt.numt2) then
!               do i=1,NI
!                  xlast(i,j) = xout(i,j)
!               enddo
!            endif
         enddo
!$OMP END DO
!$OMP END PARALLEL
      enddo  ! do t

      do j=1,NJ
        xout(1,j) = xin(1,j)
        xout(NI,j) = xin(NI,j)
      enddo

      end subroutine diffusion1x_implicit
!**********************************************************

      subroutine diffusion1y_implicit(xout,xin,diff)

!-----------------------------------------------------------------------
!
! Purpose:
! compute Lsqrt*xin (diffusion over numt2 timesteps)
! specify initial conditions
!
! Author: Mark Buehner (mark.buehner@ec.gc.ca)
!
!-----------------------------------------------------------------------

      implicit none

      double precision, intent(in)  :: xin(NI,NJ)
      double precision, intent(out) :: xout(NI,NJ)

      integer :: t, j, i

      double precision :: xlast(NJ,NI)
      double precision :: dp(NJ)

! NOTE:for improved efficiency, the 2D fields used internally are !
!      ordered (NJ,NI) and NOT (NI,NJ) as in the rest of the code !

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(j,i)
      do j=1,NJ
         do i=1,NI
            xlast(j,i) = xin(i,j)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

!      do t=1,numt2
      do t=1,1
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(j,i,dp)
         do i=2,NI-1
            j=2
            dp(j)=xlast(j,i)
            do j=3,NJ-1
               dp(j)=xlast(j,i) - diff1y_ap(j,i)*dp(j-1)
            enddo
            j=NJ-1
            xout(i,j)=dp(j)*diff1y_bp_inv(j,i)
            do j=NJ-2,2,-1
               xout(i,j)=(dp(j)-diff1y_c(j,i)*xout(i,j+1))
     +                   *diff1y_bp_inv(j,i)
            enddo 
!            if(t.lt.numt2) then
!               do j=1,NJ
!                  xlast(j,i) = xout(i,j)
!               enddo
!            endif
         enddo
!$OMP END DO
!$OMP END PARALLEL
      enddo  ! do t

      do i=1,NI
        xout(i,1) = xin(i,1)
        xout(i,NJ) = xin(i,NJ)
      enddo

      end subroutine diffusion1y_implicit


!**********************************************************

      subroutine Csqrt(xout,xin,diff)

      implicit none

      double precision, intent(in)  :: xin(NI,NJ)
      double precision, intent(out) :: xout(NI,NJ)
      type(diffus_data), intent(in) :: diff

!     compute Csqrt

!     this is the C^1/2 required for the forward model: Csqrt = Lambda * Diffuse * W^-1/2
      xout = Winvsqrt * xin
      if (diff%limplicit) then
         do t=1,diff%numt
            call diffusion1x_implicit(xout,xout,diff)
            call diffusion1y_implicit(xout,xout,diff)
         enddo
      else
         call diffusion(xout,xout,diff)
      endif
      xout = diff%Lambda * xout

      end subroutine Csqrt

!**********************************************************

      subroutine Csqrtadj(xout,xin,diff)

      implicit none

      double precision, intent(in)  :: xin(NI,NJ)
      double precision, intent(out) :: xout(NI,NJ)
      type(diffus_data), intent(in) :: diff

!     compute Csqrtadj

!     this is the (C^1/2)^T required for the adjoint: Csqrt^T = W^1/2 * Diffuse * W^-1 * Lambda
      xout = diff%Lambda * xin
      xout = Winv * xout
      if(diff%limplicit) then
         do t=1,numt2
            call diffusion1y_implicit(xout,xout,diff)
            call diffusion1x_implicit(xout,xout,diff)
         enddo
      else
         call diffusion(xout,xout,diff)
      endif
      xout = Wsqrt * xout

      end subroutine Csqrtadj

!**********************************************************

      subroutine test_diff(xout)
        double precision :: xin(NI,NJ)
        double precision, intent(out) :: xout(NI,NJ)

        xin = 0.
        xin(200,200) = 1.
        call Csqrt(xout, xin, corr_diff)
        
      end subroutine test_diff

!**********************************************************

      end module diffus_mod

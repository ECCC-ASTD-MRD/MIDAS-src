module hydrostatic_mod
  use mpi_mod
  use columnData_mod
  use utilities_mod

  implicit none
  private

!modgps00base

  ! Maximum number of gps levels:
  integer, parameter :: ngpssize  = 100

  ! Maximum number of gps extra fictitious low levels:
  integer, parameter :: ngpsxlow  = 20

  ! Associated maximum number of control variables:
  integer, parameter :: ngpscvmx  = 3*ngpssize+1


!modgps01ctphys

  ! Avogadro constant:
  real(8), parameter           :: p_Avog  = 6.02214129e23D0        ! From CODATA

  ! Boltzmann constant:
  real(8), parameter           :: p_Boltz = 1.3806488D-23        ! From CODATA

  ! Air properties:
  real(8), parameter           :: p_md    = 28.965516D0            ! From Aparicio(2011)
  real(8), parameter           :: p_mw    = 18.015254D0            ! From Aparicio(2011)
  real(8), parameter           :: p_wa    = p_md/p_mw
  real(8), parameter           :: p_wb    = (p_md-p_mw)/p_mw

  ! Gas constants:
  real(8), parameter           :: p_R     = p_Avog*p_Boltz          ! per mol
  real(8), parameter           :: p_Rd    = p_Avog*p_Boltz/(1.D-3*p_md)   ! per air mass

  ! Units and scales:
  real(8), parameter           :: p_TC    = 273.15D0
  real(8), parameter           :: p_knot  = 0.514444D0

  ! Standard GEM gravity:
  real(8), parameter           :: p_g_GEM = 9.80616D0              ! m/s2


!modgps02wgs84const

  ! Semimajor axis (a) (m)                             [*Defining constant*]
  real(8), parameter :: WGS_a          = 6378137.0D0
  
  ! Flattening (f)                                     [*Defining constant*]
  real(8), parameter :: WGS_f          = 1.D0 / 298.257223563D0
  
  ! Earth's Gravitational Constant (GM) (m3/s2).       [*Defining constant*]
  ! This is the GM value with Earth's atmosphere included.
  real(8), parameter :: WGS_GM         = 3986004.418D8
  
  ! Angular velocity of the Earth (omega) (radians/s). [*Defining constant*]
  ! Standard Earth, rotating with a constant angular velocity.
  real(8), parameter :: WGS_Omega      = 7292115.D-11


  !  ***********************************************************
  !  Parameter values for special applications:
  !  ***********************************************************
  
  ! Earth's Gravitational Constant (GM_KEP) (m3/s2).
  !
  ! This is the old GM value with Earth's atmosphere included,
  ! but still used for consistency in the transformations
  ! instantaneous keplerian <-> cartesian state vector.
  real(8), parameter :: WGS_GM_KEP     = 3986005.0D8
  
  ! Earth's Atmosphere Gravitational Constant (GMA) (m3/s2).
  ! This is the GM value of the Earth's atmosphere alone.
  real(8), parameter :: WGS_GMA        = 3.5D8
  
  ! Earth's Gravitational Constant (GMPrime) (m3/s2).
  ! This is the GM value with Earth's atmosphere excluded.
  real(8), parameter :: WGS_GMPrime    = 3986000.9D8
  
  ! Angular velocity of the Earth (omegaPrime) (radians/s).
  ! Standard Earth, rotating with a constant angular velocity (IAU, GRS67).
  real(8), parameter :: WGS_OmegaPrime = 7292115.1467D-11
  
  ! Angular velocity of the Earth (omegaStar) (radians/s).
  ! Standard Earth, in a precessing frame (TU: Julian centuries since J2000.0)
  ! OmegaStar = OmegaStar0 + OmegaStar1 * TU
  real(8), parameter :: WGS_OmegaStar0 = 7292115.8553D-11
  real(8), parameter :: WGS_OmegaStar1 = 4.3D-15


  !  ***********************************************************
  !  Derived geometric constants:
  !  ***********************************************************
  
  ! Inverse flattening: 1/f = WGS_1f
  real(8), parameter :: WGS_1f         = 298.257223563D0
  
  ! Second degree zonal coefficient:
  real(8), parameter :: WGS_C20        = -0.484166774985D-3
  
  ! Semiminor axis:
  real(8), parameter :: WGS_b          = 6356752.3142D0
  
  ! First eccentricity:
  real(8), parameter :: WGS_e          = 8.1819190842622D-2
  
  ! First eccentricity squared:
  real(8), parameter :: WGS_e2         = 6.69437999014D-3
  
  ! Second eccentricity:
  real(8), parameter :: WGS_ePrime     = 8.2094437949696D-2
  
  ! Second eccentricity squared:
  real(8), parameter :: WGS_ePrime2    = 6.73949674228D-3
  
  ! Linear eccentricity:
  real(8), parameter :: WGS_ELinear    = 5.2185400842339D5
  
  ! Polar radius of curvature:
  real(8), parameter :: WGS_c          = 6399593.6258D0
  
  ! Focal length:
  real(8), parameter :: WGS_EFocal     = 521854.00897D0
  
  ! Axis ratio:
  real(8), parameter :: WGS_ba         = 0.996647189335D0
  
  ! Mean radius of semiaxes:
  real(8), parameter :: WGS_R1         = 6371008.7714D0
  
  ! Radius of sphere of equal area:
  real(8), parameter :: WGS_R2         = 6371007.1809D0
  
  ! Radius of sphere of equal volume:
  real(8), parameter :: WGS_R3         = 6371000.7900D0


  !  ***********************************************************
  !  Derived physical constants:
  !  ***********************************************************
  
  ! Theoretical (Normal) Gravity potential of the ellipsoid (m2/s2):
  real(8), parameter :: WGS_U0         = 62636860.8497D0
  
  ! Theoretical (Normal) Gravity at the equator (m/s2):
  real(8), parameter :: WGS_GammaE     = 9.7803253359D0
  
  ! Theoretical (Normal) Gravity at the pole (m/s2):
  real(8), parameter :: WGS_GammaP     = 9.8321849378D0
  
  ! Mean Value of the Theoretical (Normal) Gravity (m/s2):
  real(8), parameter :: WGS_GammaM     = 9.7976432222D0
  
  ! Theoretical (Normal) Gravity Formula Constant:
  real(8), parameter :: WGS_TNGk       = 0.00193185265241D0
  
  ! Mass of the Earth (Atmosphere Included):
  real(8), parameter :: WGS_Mass       = 5.9733328e24D0
  
  ! m = omega^2 a^2 b / GM
  real(8), parameter :: WGS_m          = 0.00344978650684D0
  
  ! Dynamical Ellipticity H:
  real(8), parameter :: WGS_H          = 1.D0 / 305.4413D0
  
  ! Universal Constant of Gravitation (value used in WGS) (m3/kg*s2):
  real(8), parameter :: WGS_G          = 6.673D-11
  
  ! Earth's principal moments of inertia (A, B, C) (kg m2):
  real(8), parameter :: WGS_PMI_A      = 8.0091029D37
  real(8), parameter :: WGS_PMI_B      = 8.0092559D37
  real(8), parameter :: WGS_PMI_C      = 8.0354872D37


!modgps03diff

  type gps_diff
     real(8)           :: Var
     real(8)           :: DVar(ngpscvmx)
  end type gps_diff


  contains

  subroutine hydrostat
    implicit none
    ! 
    ! Purpose: transfer the hydrostatic calculation outside gps_mod.f90
    !

!modgps04profile

    real(8), intent(in)  :: rLat
    real(8), intent(in)  :: rLon
    real(8), intent(in)  :: rAzm
    real(8), intent(in)  :: rMT
    real(8), intent(in)  :: Rad
    real(8), intent(in)  :: geoid
    real(8), intent(in)  :: rP0
    real(8), intent(in)  :: rPP (ngpssize)
    real(8), intent(in)  :: rDP (ngpssize)
    real(8), intent(in)  :: rTT (ngpssize)
    real(8), intent(in)  :: rHU (ngpssize)
    real(8), intent(in)  :: rUU (ngpssize)
    real(8), intent(in)  :: rVV (ngpssize)

    integer :: lev_T
    integer :: nlev_M,nlev_T

    real(8), parameter :: delta = 0.6077686814144D0

    real(8) :: h0,dh,Rgh,Eot,Eot2, sLat, cLat
    real(8) :: p, t, q, x
    real(8) :: tr, z
    real(8) :: mold, dd, dw, dx, n0, nd1, nw1, tvm
    real(8) :: xi(ngpssize), tv(ngpssize)

    nlev_T = col_getNumLev(columng,'TH')
    nlev_M = col_getNumLev(columng,'MM')

!!! LOOP OVER columnIndex

    do lev_T = 1, nlev_T
      p = col_getPressure(column,lev_T,columnIndex,'TH')
      t = col_getElem(column,lev_T,columnIndex,'TT')
      q = col_getElem(column,lev_T,columnIndex,'HU')
      !
      ! Log(P)
      !
      xi(lev_T) = log(p)
      !
      ! Virtual temperature (K) (corrected of compressibility)
      !
      tv(lev_T) = (1.D0 + delta * q) * t * gpscompressibility(p,t,q)
    enddo

    sLat = sin(rLat)
    cLat = cos(rLat)
    dx  = xi(nlev_T) - log(prf%P0)
    Rgh = gpsgravitysrf(sLat)
    z   = (-p_Rd / Rgh) * tv(nlev_T) * dx
!!! calculate rMT
    prf%gst(nlev_T) = rMT + z
    do i = nlev_T-1, 1, -1
      dx = xi(i) - xi(i+1)
      tvm = 0.5D0 * (tv(i) + tv(i+1))
      !
      ! Gravity acceleration (includes 2nd-order Eotvos effect)
      !
      h0  = prf%gst(i+1)
      Eot = 2 * WGS_OmegaPrime * cLat * p_knot * rUU(i)
      Eot2= ((p_knot*rUU(i)) ** 2 + (p_knot * rVV(i)) ** 2) / WGS_a
      Rgh = gpsgravityalt(sLat, h0) - Eot - Eot2
      dh  = (-p_Rd / Rgh) * tvm * dx
      Rgh = gpsgravityalt(sLat,h0+0.5D0*dh) - Eot - Eot2
      !
      ! Height increment
      !
      z   = (-p_Rd / Rgh) * tvm * dx
      prf%gst(i) = prf%gst(i+1) + z
    enddo

  end subroutine hydrostat


  function gpscompressibility(p,t,q)
    type(gps_diff), intent(in)  :: p,t,q
    type(gps_diff)              :: gpscompressibility

    real(8), parameter   :: a0= 1.58123D-6
    real(8), parameter   :: a1=-2.9331D-8
    real(8), parameter   :: a2= 1.1043D-10
    real(8), parameter   :: b0= 5.707D-6
    real(8), parameter   :: b1=-2.051D-8
    real(8), parameter   :: c0= 1.9898D-4
    real(8), parameter   :: c1=-2.376D-6
    real(8), parameter   :: d = 1.83D-11
    real(8), parameter   :: e =-0.765D-8

    type(gps_diff)         :: x,tc,pt,tc2,x2

    x  = p_wa*q/(1.D0+p_wb*q)
    ! Estimate, from CIPM, Picard (2008)
    tc = t-p_TC
    pt = 1.D2*p/t
    tc2= tc*tc
    x2 = x*x
    gpscompressibility = 1.D0-pt*(a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2)+pt*pt*(d+e*x2)
  end function gpscompressibility


  pure function gpsgravitysrf(sLat)
    !  Normal gravity on ellipsoidal surface:
    !  Input:  Latitude
    !          sin(Latitude)
    !
    !  Output: Normal gravity
    !          gpsgravitysrf         : m/s2
    !
    real(8), intent(in)  :: sLat
    real(8)              :: gpsgravitysrf
    
    real(8)              :: ks2
    real(8)              :: e2s

    ks2 = WGS_TNGk * sLat*sLat
    e2s = 1.D0 - WGS_e2 * sLat*sLat
    gpsgravitysrf = WGS_GammaE * (1.D0 + ks2) / sqrt(e2s)
  end function gpsgravitysrf


  pure function gpsgravityalt(sLat, Altitude)
    ! Normal gravity above the ellipsoidal surface:
    ! Input:  Latitude, altitude
    !         sin(Latitude)
    !         Altitude               : m
    !
    ! Output: Normal gravity
    !         gpsgravityalt          : m/s2
    !
    real(8), intent(in)  :: sLat
    real(8), intent(in)  :: Altitude
    real(8)              :: gpsgravityalt

    real(8)              :: C1
    real(8)              :: C2

    C1 =-2.D0/WGS_a*(1.D0+WGS_f+WGS_m-2*WGS_f*sLat*sLat)
    C2 = 3.D0/WGS_a**2
    gpsgravityalt = gpsgravitysrf(sLat)*                                   &
         (1.D0 + C1 * Altitude + C2 * Altitude**2)
  end function gpsgravityalt


end module hydrostatic_mod

MODULE earthConstants_mod
  ! MODULE earthConstants_mod (prefix='ec' category='8. Low-level utilities and constants')
  !
  !:Purpose: Define public constants related to the Earth.
  !
  !          Prefixes:
  !           * ec_ (Earth constants), for miscellaneous values from diverse sources or
  !           * ec_wgs_ (World Geodetic Syatem) when they are explicitly from WGS84 specification.
  !
  use mathPhysConstants_mod
  !
  ! The following constants should ultimately be taken from module
  ! modgps02wgs84const OR module modgps06gravity.  They have been placed here
  ! as an intermediate step.
  real(8), parameter :: ec_Omega  = 7292000.D-11               ! Approximation, see for instance ec_wgs_Omega below.
  real(8), parameter :: ec_rg     = 9.80616000000000000000D+00 ! ...616 ? This is Lambert IGF at 45 lat (g45). The WMO standard for conversions is ...665.
  real(8), parameter :: ec_ra     = 6371229.00000000000000D+00 ! ECMWF topographical mean radius (Authalic radius, eg ec_wgs_R2,+average topography)
  real(8), parameter :: ec_r1sa   = 1.D0/ec_ra                 ! 1.56955588945241177170D-07          ! Modified to explicitly 1/ec_ra
  !real(8), parameter :: RV       = 4.61524993308387870310D+02 ! WV gas constant? Seems unused.
  real(8), parameter :: ec_rayt   = 6371220.00000000000000D+00 ! Another topographical mean (NOAA)
  !
  ! Further values through the code:
  !
  real(8), parameter :: ec_RDay    = 86400.D0                                   ! Standard day
  real(8), parameter :: ec_JYear   = 365.25D0                                   ! Standard (Julian) year
  real(8), parameter :: ec_SidYear = ec_JYear*ec_RDAY*2.D0*MPC_PI_R8/6.283076D0 ! Not really sidereal year, seems rather anomalistic year (odd).
  real(8), parameter :: ec_SidDay  = ec_RDAY/(1.D0+ec_RDAY/ec_SidYear)          ! Approx sidereal day (wrt perihelion?), odd expression.
  real(8), parameter :: ec_ROmega  = 2.D0*MPC_PI_R8/ec_SidDay                   ! Angular rotation (approx., again)
  !                                                                             ! ec_SidDay ~= 86164.0996
  !                                                                             ! ec_ROmega ~= 7292115.D-11
  !                                                                             ! First seven digits as ec_wgs_Omega

  ! WGS (World Geodetic System, 1984) constants

  ! Semimajor axis (a) (m)                             [*Defining constant*]
  real(8), parameter :: ec_wgs_a = 6378137.0D0

  ! Flattening (f)                                     [*Defining constant*]
  real(8), parameter :: ec_wgs_f = 1.D0 / 298.257223563D0

  ! m = omega^2 a^2 b / GM
  real(8), parameter :: ec_wgs_m = 0.00344978650684D0

  ! Theoretical (Normal) Gravity Formula Constant:
  real(8), parameter :: ec_wgs_TNGk = 0.00193185265241D0

  ! First eccentricity squared:
  real(8), parameter :: ec_wgs_e2 = 6.69437999014D-3

  ! Theoretical (Normal) Gravity at the equator (m/s2):
  real(8), parameter :: ec_wgs_GammaE = 9.7803253359D0

  ! Earth's Gravitational Constant (GM) (m3/s2).       [*Defining constant*]
  ! This is the GM value with Earth's atmosphere included.
  real(8), parameter :: ec_wgs_GM         = 3986004.418D8

  ! Angular velocity of the Earth (omega) (radians/s). [*Defining constant*]
  ! Standard Earth, rotating with a constant angular velocity.
  real(8), parameter :: ec_wgs_Omega      = 7292115.D-11

  !  ***********************************************************
  !  Parameter values for special applications:
  !  ***********************************************************

  ! Earth's Gravitational Constant (GM_KEP) (m3/s2).
  !
  ! This is the old GM value with Earth's atmosphere included,
  ! but still used for consistency in the transformations
  ! instantaneous keplerian <-> cartesian state vector.
  real(8), parameter :: ec_wgs_GM_KEP     = 3986005.0D8

  ! Earth's Atmosphere Gravitational Constant (GMA) (m3/s2).
  ! This is the GM value of the Earth's atmosphere alone.
  real(8), parameter :: ec_wgs_GMA        = 3.5D8

  ! Earth's Gravitational Constant (GMPrime) (m3/s2).
  ! This is the GM value with Earth's atmosphere excluded.
  real(8), parameter :: ec_wgs_GMPrime    = 3986000.9D8

  ! Angular velocity of the Earth (omegaPrime) (radians/s).
  ! Standard Earth, rotating with a constant angular velocity (IAU, GRS67).
  real(8), parameter :: ec_wgs_OmegaPrime = 7292115.1467D-11

  ! Angular velocity of the Earth (omegaStar) (radians/s).
  ! Standard Earth, in a precessing frame (TU: Julian centuries since J2000.0)
  ! OmegaStar = OmegaStar0 + OmegaStar1 * TU
  real(8), parameter :: ec_wgs_OmegaStar0 = 7292115.8553D-11
  real(8), parameter :: ec_wgs_OmegaStar1 = 4.3D-15

  !  ***********************************************************
  !  Derived geometric constants:
  !  ***********************************************************

  ! Inverse flattening: 1/f = ec_wgs_1f
  real(8), parameter :: ec_wgs_1f         = 298.257223563D0

  ! Second degree zonal coefficient:
  real(8), parameter :: ec_wgs_C20        = -0.484166774985D-3

  ! Semiminor axis:
  real(8), parameter :: ec_wgs_b          = 6356752.3142D0

  ! First eccentricity:
  real(8), parameter :: ec_wgs_e          = 8.1819190842622D-2

  ! Second eccentricity:
  real(8), parameter :: ec_wgs_ePrime     = 8.2094437949696D-2

  ! Second eccentricity squared:
  real(8), parameter :: ec_wgs_ePrime2    = 6.73949674228D-3

  ! Linear eccentricity:
  real(8), parameter :: ec_wgs_ELinear    = 5.2185400842339D5

  ! Polar radius of curvature:
  real(8), parameter :: ec_wgs_c          = 6399593.6258D0

  ! Focal length:
  real(8), parameter :: ec_wgs_EFocal     = 521854.00897D0

  ! Axis ratio:
  real(8), parameter :: ec_wgs_ba         = 0.996647189335D0

  ! Mean radius of semiaxes:
  real(8), parameter :: ec_wgs_R1         = 6371008.7714D0

  ! Radius of sphere of equal area:
  real(8), parameter :: ec_wgs_R2         = 6371007.1809D0

  ! Radius of sphere of equal volume:
  real(8), parameter :: ec_wgs_R3         = 6371000.7900D0

  !  ***********************************************************
  !  Derived physical constants:
  !  ***********************************************************

  ! Theoretical (Normal) Gravity potential of the ellipsoid (m2/s2):
  real(8), parameter :: ec_wgs_U0         = 62636860.8497D0

  ! Theoretical (Normal) Gravity at the pole (m/s2):
  real(8), parameter :: ec_wgs_GammaP     = 9.8321849378D0

  ! Mean Value of the Theoretical (Normal) Gravity (m/s2):
  real(8), parameter :: ec_wgs_GammaM     = 9.7976432222D0

  ! Mass of the Earth (Atmosphere Included):
  real(8), parameter :: ec_wgs_Mass       = 5.9733328D24

  ! Dynamical Ellipticity H:
  real(8), parameter :: ec_wgs_H          = 1.D0 / 305.4413D0

  ! Universal Constant of Gravitation (m3/kg*s2):
  ! (value to be used only within WGS, for internal consistency)
  real(8), parameter :: ec_wgs_G          = 6.673D-11

  ! Earth's principal moments of inertia (A, B, C) (kg m2):
  real(8), parameter :: ec_wgs_PMI_A      = 8.0091029D37
  real(8), parameter :: ec_wgs_PMI_B      = 8.0092559D37
  real(8), parameter :: ec_wgs_PMI_C      = 8.0354872D37

end MODULE earthConstants_mod

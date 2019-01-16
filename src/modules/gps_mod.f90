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
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------

!--------------------------------------------------------------------------
!! MODULE gps (prefix='gps' category='4. Observation operators')
!!
!! Note: prefix not used for all public variables
!!
!! *Purpose*: Code related to GPS-RO and ground-based GPS observation operators. 
!!
!--------------------------------------------------------------------------
module gps_mod
  use mpi_mod
  use utilities_mod
  implicit none
  save
  private

  ! public types
  public :: gps_profile, gps_profilezd, gps_diff

  ! public variables
  public :: gps_numROProfiles, gps_vRO_IndexPrf, gps_vRO_Jacobian, gps_vRO_lJac
  public :: LEVELGPSRO, GPSRO_MAXPRFSIZE, NUMGPSSATS, IGPSSAT, SURFMIN, HSFMIN, HTPMAX, BGCKBAND, WGPS
  public :: gpsgravitysrf, p_tc, p_knot, max_gps_data, vgpsztd_jacobian, vgpsztd_ljac, dzmin
  public :: ltestop, llblmet, lbevis, irefopt, iztdop, lassmet, l1obs, yzderrwgt, numgpsztd
  public :: vgpsztd_index, ngpscvmx, dzmax, yztderr, ysferrwgt

  ! public procedures
  public :: gps_setupro, gps_iprofile_from_index
  public :: gps_setupgb, gps_i_from_index
  public :: gps_struct1sw, gps_bndopv1, gps_refopv, gps_structztd, gps_ztdopv, gps_pw


!modgps00base
  
  ! 32-bit integers
  integer, parameter     :: i4 = selected_int_kind(9)
  
  ! Short floats
  integer, parameter     :: sp = selected_real_kind(6)
  
  ! Long floats
  integer, parameter     :: dp = selected_real_kind(12)

  ! Maximum number of gps levels:
  integer(i4), parameter :: ngpssize  = 100

  ! Maximum number of gps extra fictitious low levels:
  integer(i4), parameter :: ngpsxlow  = 20

  ! Associated maximum number of control variables:
  integer(i4), parameter :: ngpscvmx  = 2*ngpssize+1

  
!modgps01ctphys
  
  ! Avogadro constant:
  real(dp), parameter           :: p_Avog  = 6.02214129e23_dp        ! From CODATA

  ! Boltzmann constant:
  real(dp), parameter           :: p_Boltz = 1.3806488e-23_dp        ! From CODATA

  ! Air properties:
  real(dp), parameter           :: p_md    = 28.965516_dp            ! From Aparicio(2011)
  real(dp), parameter           :: p_mw    = 18.015254_dp            ! From Aparicio(2011)
  real(dp), parameter           :: p_wa    = p_md/p_mw
  real(dp), parameter           :: p_wb    = (p_md-p_mw)/p_mw

  ! Gas constants:
  real(dp), parameter           :: p_R     = p_Avog*p_Boltz          ! per mol
  real(dp), parameter           :: p_Rd    = p_Avog*p_Boltz/(1.e-3_dp*p_md)   ! per air mass

  ! Units and scales:
  real(dp), parameter           :: p_TC    = 273.15_dp
  real(dp), parameter           :: p_knot  = 0.514444_dp

  ! Standard GEM gravity:
  real(dp), parameter           :: p_g_GEM = 9.80616_dp              ! m/s2


!modgps02wgs84const

  ! Semimajor axis (a) (m)                             [*Defining constant*]
  real(dp), parameter :: WGS_a          = 6378137.0_dp
  
  ! Flattening (f)                                     [*Defining constant*]
  real(dp), parameter :: WGS_f          = 1._dp / 298.257223563_dp
  
  ! Earth's Gravitational Constant (GM) (m3/s2).       [*Defining constant*]
  ! This is the GM value with Earth's atmosphere included.
  real(dp), parameter :: WGS_GM         = 3986004.418e8_dp
  
  ! Angular velocity of the Earth (omega) (radians/s). [*Defining constant*]
  ! Standard Earth, rotating with a constant angular velocity.
  real(dp), parameter :: WGS_Omega      = 7292115.e-11_dp
  
  !  ***********************************************************
  !  Parameter values for special applications:
  !  ***********************************************************
  
  ! Earth's Gravitational Constant (GM_KEP) (m3/s2).
  !
  ! This is the old GM value with Earth's atmosphere included,
  ! but still used for consistency in the transformations
  ! instantaneous keplerian <-> cartesian state vector.
  real(dp), parameter :: WGS_GM_KEP     = 3986005.0e8_dp
  
  ! Earth's Atmosphere Gravitational Constant (GMA) (m3/s2).
  ! This is the GM value of the Earth's atmosphere alone.
  real(dp), parameter :: WGS_GMA        = 3.5e8_dp
  
  ! Earth's Gravitational Constant (GMPrime) (m3/s2).
  ! This is the GM value with Earth's atmosphere excluded.
  real(dp), parameter :: WGS_GMPrime    = 3986000.9e8_dp
  
  ! Angular velocity of the Earth (omegaPrime) (radians/s).
  ! Standard Earth, rotating with a constant angular velocity (IAU, GRS67).
  real(dp), parameter :: WGS_OmegaPrime = 7292115.1467e-11_dp
  
  ! Angular velocity of the Earth (omegaStar) (radians/s).
  ! Standard Earth, in a precessing frame (TU: Julian centuries since J2000.0)
  ! OmegaStar = OmegaStar0 + OmegaStar1 * TU
  real(dp), parameter :: WGS_OmegaStar0 = 7292115.8553e-11_dp
  real(dp), parameter :: WGS_OmegaStar1 = 4.3e-15_dp
  
  
  !  ***********************************************************
  !  Derived geometric constants:
  !  ***********************************************************
  
  ! Inverse flattening: 1/f = WGS_1f
  real(dp), parameter :: WGS_1f         = 298.257223563_dp
  
  ! Second degree zonal coefficient:
  real(dp), parameter :: WGS_C20        = -0.484166774985e-3_dp
  
  ! Semiminor axis:
  real(dp), parameter :: WGS_b          = 6356752.3142_dp
  
  ! First eccentricity:
  real(dp), parameter :: WGS_e          = 8.1819190842622e-2_dp
  
  ! First eccentricity squared:
  real(dp), parameter :: WGS_e2         = 6.69437999014e-3_dp
  
  ! Second eccentricity:
  real(dp), parameter :: WGS_ePrime     = 8.2094437949696e-2_dp
  
  ! Second eccentricity squared:
  real(dp), parameter :: WGS_ePrime2    = 6.73949674228e-3_dp
  
  ! Linear eccentricity:
  real(dp), parameter :: WGS_ELinear    = 5.2185400842339e5_dp
  
  ! Polar radius of curvature:
  real(dp), parameter :: WGS_c          = 6399593.6258_dp
  
  ! Focal length:
  real(dp), parameter :: WGS_EFocal     = 521854.00897_dp
  
  ! Axis ratio:
  real(dp), parameter :: WGS_ba         = 0.996647189335_dp
  
  ! Mean radius of semiaxes:
  real(dp), parameter :: WGS_R1         = 6371008.7714_dp
  
  ! Radius of sphere of equal area:
  real(dp), parameter :: WGS_R2         = 6371007.1809_dp
  
  ! Radius of sphere of equal volume:
  real(dp), parameter :: WGS_R3         = 6371000.7900_dp
  
  
  !  ***********************************************************
  !  Derived physical constants:
  !  ***********************************************************
  
  ! Theoretical (Normal) Gravity potential of the ellipsoid (m2/s2):
  real(dp), parameter :: WGS_U0         = 62636860.8497_dp
  
  ! Theoretical (Normal) Gravity at the equator (m/s2):
  real(dp), parameter :: WGS_GammaE     = 9.7803253359_dp
  
  ! Theoretical (Normal) Gravity at the pole (m/s2):
  real(dp), parameter :: WGS_GammaP     = 9.8321849378_dp
  
  ! Mean Value of the Theoretical (Normal) Gravity (m/s2):
  real(dp), parameter :: WGS_GammaM     = 9.7976432222_dp
  
  ! Theoretical (Normal) Gravity Formula Constant:
  real(dp), parameter :: WGS_TNGk       = 0.00193185265241_dp
  
  ! Mass of the Earth (Atmosphere Included):
  real(dp), parameter :: WGS_Mass       = 5.9733328e24_dp
  
  ! m = omega^2 a^2 b / GM
  real(dp), parameter :: WGS_m          = 0.00344978650684_dp
  
  ! Dynamical Ellipticity H:
  real(dp), parameter :: WGS_H          = 1._dp / 305.4413_dp
  
  ! Universal Constant of Gravitation (value used in WGS) (m3/kg*s2):
  real(dp), parameter :: WGS_G          = 6.673e-11_dp
  
  ! Earth's principal moments of inertia (A, B, C) (kg m2):
  real(dp), parameter :: WGS_PMI_A      = 8.0091029e37_dp
  real(dp), parameter :: WGS_PMI_B      = 8.0092559e37_dp
  real(dp), parameter :: WGS_PMI_C      = 8.0354872e37_dp


!modgps03diff

  type gps_diff
     real(dp)           :: Var
     real(dp)           :: DVar(ngpscvmx)
  end type gps_diff
  
  interface assignment(=)
     module procedure gpsdiffasfd, gpsdiffasff
  end interface
  
  interface operator(+)
     module procedure gpsdiffsmfd, gpsdiffsmdf, gpsdiffsmfi, gpsdiffsmif, gpsdiffsmff
  end interface
  
  interface operator(-)
     module procedure gpsdiffsbfd, gpsdiffsbdf, gpsdiffsbfi, gpsdiffsbif, gpsdiffsbff
  end interface
  
  interface operator(*)
     module procedure gpsdiffmlfd, gpsdiffmldf, gpsdiffmlfi, gpsdiffmlif, gpsdiffmlff
  end interface
  
  interface operator(/)
     module procedure gpsdiffdvfd, gpsdiffdvdf, gpsdiffdvfi, gpsdiffdvif, gpsdiffdvff
  end interface
  
  interface operator(**)
     module procedure gpsdiffpwfd, gpsdiffpwdf, gpsdiffpwfi, gpsdiffpwif, gpsdiffpwff
  end interface

  interface sqrt
     module procedure gpsdiffsqr
  end interface

  interface exp
     module procedure gpsdiffexp
  end interface
  
  interface log
     module procedure gpsdifflog
  end interface

  interface cos
     module procedure gpsdiffcos
  end interface

  interface tan
     module procedure gpsdifftan
  end interface

  interface acos
     module procedure gpsdiffacos
  end interface

  interface atan
     module procedure gpsdiffatan
  end interface

  interface erf
     module procedure gpsdifferf
  end interface

!modgps04profile

  type gps_profile
     integer(i4)                                     :: ngpslev
     real(dp)                                        :: rLat
     real(dp)                                        :: rLon
     real(dp)                                        :: rAzm
     real(dp)                                        :: rMT
     real(dp)                                        :: Rad
     real(dp)                                        :: geoid
     real(dp)                                        :: RadN
     real(dp)                                        :: RadM

     type(gps_diff)                                   :: P0

     type(gps_diff)    , dimension(ngpssize)          :: pst
     type(gps_diff)    , dimension(ngpssize)          :: tst
     type(gps_diff)    , dimension(ngpssize)          :: qst
     type(gps_diff)    , dimension(ngpssize)          :: rst
     type(gps_diff)    , dimension(ngpssize)          :: gst

     logical                                         :: bbst
     type(gps_diff)    , dimension(ngpssize)          :: dst
     type(gps_diff)    , dimension(ngpssize+ngpsxlow) :: ast
     type(gps_diff)    , dimension(ngpssize+ngpsxlow) :: bst
  end type gps_profile


!modgps04profilezd
  
  type gps_profilezd
     integer(i4)                                     :: ngpslev
     real(dp)                                        :: rLat
     real(dp)                                        :: rLon
     real(dp)                                        :: rMT

     type(gps_diff)                                   :: P0
     
     type(gps_diff)    , dimension(ngpssize)          :: pst
     type(gps_diff)    , dimension(ngpssize)          :: tst
     type(gps_diff)    , dimension(ngpssize)          :: qst
     type(gps_diff)    , dimension(ngpssize)          :: rst
     type(gps_diff)    , dimension(ngpssize)          :: gst
     type(gps_diff)    , dimension(ngpssize)          :: ztd
     logical                                         :: bpst
  end type gps_profilezd


!modgpsro_mod

!
! Values determined by input data:
!
  integer                                :: gps_numROProfiles
  integer         , allocatable          :: gps_vRO_IndexPrf(:)    ! index for each profile
  real*8          , allocatable          :: gps_vRO_Jacobian(:,:,:)
  logical         , allocatable          :: gps_vRO_lJac(:)

!     Contents of previous comdeck comgpsro
!     -------------------------------------
!*    Control variables for GPSRO observations - constant within job
!
!     LEVELGPSRO: Data level to use (1 for bending angle, 2 for refractivity)
!     GPSRO_MAXPRFSIZE: Maximal number of data that is expected from a profile (default 300)
!     SURFMIN:  Minimum allowed distance to the model surface (default 1000 m)
!     HSFMIN:   Minimum allowed MSL height of an obs          (default 4000 m)
!     HTPMAX:   Maximum allowed MSL height of an obs          (default 40000 m)
!     BGCKBAND: Maximum allowed deviation abs(O-P)/P          (default 0.05)
!
!     J.M. Aparicio, Apr 2008
!          
  INTEGER LEVELGPSRO, GPSRO_MAXPRFSIZE,NUMGPSSATS,IGPSSAT(50)
  REAL*8  SURFMIN, HSFMIN, HTPMAX, BGCKBAND, WGPS(50)

  NAMELIST /NAMGPSRO/ LEVELGPSRO,GPSRO_MAXPRFSIZE,SURFMIN,HSFMIN,HTPMAX,BGCKBAND,NUMGPSSATS,IGPSSAT,WGPS


!modgpsztd_mod

  integer, parameter      ::  max_gps_sites = 1200
  integer, parameter      ::  max_gps_data  = max_gps_sites*24     ! (max_gps_sites) * (max_num_obs in 6h)

  integer                 :: numGPSZTD                ! number of ZTD data to be assimilated
  integer , allocatable   :: vGPSZTD_Index (:)        ! INDEX_HEADER in CMA (ObsSpace) for each ZTD observation
  real*8  , allocatable   :: vGPSZTD_Jacobian (:,:)   ! Jacobian for each ZTD observation (numGPSZTD,ncv)
                                                              ! ncv = 2*nlev+1 = 161  (TTx80, LQx80, P0)
  logical , allocatable   :: vGPSZTD_lJac (:)         ! logical = true once Jacobian computed/stored

!*    Namelist variables for Ground-based GPS (ZTD)
!
!     DZMIN:      Minimum DZ = Zobs-Zmod (m) for which DZ adjustment to ZTD 
!                 will be made.
!     YSFERRWGT:  Weighting factor multiplier for GPS surface met errors (to 
!                 account for time series observations with error correlations)
!     DZMAX:      Maximum DZ (m) over which the ZTD data are rejected
!                 due to topography (used in SOBSSFC when LTOPOFILT = .TRUE.)
!     YZTDERR:    If < 0 then read ZTD errors from data blocks in input
!                 files (i.e. the formal errors). 
!                 If > 0 then use value as a constant error (m) for all ZTD
!                 observations.
!                 If = 0 then compute error as a function of ZWD.
!     LASSMET:    Flag to assimilate GPS Met surface P, T, T-Td
!     LLBLMET:    Flag to indicate that surface met data have been blacklisted
!                 for GPS sites close to surface weather stations.
!     YZDERRWGT:  Weighting factor multiplier for GPS ZTD errors (to account
!                 for time series observations with error correlations)
!     LBEVIS:     .true.  = use Bevis(1994)  refractivity (k1,k2,k3) constants
!                 .false. = use Rueger(2002) refractivity (k1,k2,k3) constants
!     IREFOPT:    1 = conventional expression for refractivity N using k1,k2,k3
!                 2 = Aparicio & Laroche refractivity N (incl. compressibility)
!     L1OBS       Flag to select a single ZTD observation using criteria in
!                 subroutine DOBSGPSGB
!     LTESTOP     Flag to test ZTD observation operator (Omp and Bgck modes only)
!                 Runs subroutine SETFGEGPS to do the test.
!     IZTDOP      1 = normal mode: use stored ZTD profiles to get ZTDmod
!                 2 = Vedel & Huang ZTD formulation: ZTDmod = ZHD(Pobs) + ZWD
!
  REAL*8  DZMIN, YZTDERR, YSFERRWGT, YZDERRWGT
  REAL(8) :: DZMAX = 1000.0D0 ! need to give it a default value here in case setup not called
  LOGICAL LASSMET, LLBLMET, LBEVIS, L1OBS, LTESTOP
  INTEGER IREFOPT, IZTDOP

  NAMELIST /NAMGPSGB/ DZMIN, DZMAX, YZTDERR, LASSMET, YSFERRWGT,  &
       LLBLMET, YZDERRWGT, LBEVIS, L1OBS, LTESTOP, IREFOPT, IZTDOP


contains


!modgps02wgs84grav

  !  Normal gravity on ellipsoidal surface:
  !  Input:  Latitude
  !          sin(Latitude)
  !
  !  Output: Normal gravity
  !          gpsgravitysrf         : m/s2
  !
  pure function gpsgravitysrf(sLat)
    real(dp), intent(in)  :: sLat
    real(dp)              :: gpsgravitysrf
    
    real(dp)              :: ks2
    real(dp)              :: e2s

    ks2 = WGS_TNGk * sLat*sLat
    e2s = 1._dp - WGS_e2 * sLat*sLat
    gpsgravitysrf = WGS_GammaE * (1._dp + ks2) / sqrt(e2s)
  end function gpsgravitysrf

  ! Normal gravity above the ellipsoidal surface:
  ! Input:  Latitude, altitude
  !         sin(Latitude)
  !         Altitude               : m
  !
  ! Output: Normal gravity
  !         gpsgravityalt          : m/s2
  !
  pure function gpsgravityalt(sLat, Altitude)
    real(dp), intent(in)  :: sLat
    real(dp), intent(in)  :: Altitude
    real(dp)              :: gpsgravityalt

    real(dp)              :: C1
    real(dp)              :: C2

    C1 =-2._dp/WGS_a*(1._dp+WGS_f+WGS_m-2*WGS_f*sLat*sLat)
    C2 = 3._dp/WGS_a**2
    gpsgravityalt = gpsgravitysrf(sLat)*                                   &
         (1._dp + C1 * Altitude + C2 * Altitude**2)
  end function gpsgravityalt

  ! Geopotential energy at a given point.
  ! Result is based on the WGS84 approximate expression for the
  ! gravity acceleration as a function of latitude and altitude,
  ! integrated with the trapezoidal rule.
  ! Input:  Latitude, altitude
  !         Latitude               : rad
  !         Altitude               : m
  !
  ! Output: Geopotential
  !         gpsgeopotential                              : m2/s2
  pure function gpsgeopotential(Latitude, Altitude)
    real(dp), intent(in)  :: Latitude
    real(dp), intent(in)  :: Altitude
    real(dp)              :: gpsgeopotential

    real(dp)              :: dh, sLat
    integer               :: n, i
    real(dp), allocatable :: hi(:)
    real(dp), allocatable :: gi(:)
    
    dh = 500._dp
    n = 1 + int(Altitude/dh)

    allocate(hi(0:n))
    allocate(gi(0:n))

    sLat=sin(Latitude)

    do i = 0, n-1
       hi(i) = i * dh
       gi(i) = gpsgravityalt(sLat, hi(i))
    enddo
    hi(n) = Altitude
    gi(n) = gpsgravityalt(sLat, hi(n))

    gpsgeopotential = 0._dp
    do i = 1, n
       gpsgeopotential = gpsgeopotential + 0.5_dp * (gi(i)+gi(i-1)) * (hi(i)-hi(i-1))
    enddo

    deallocate(hi)
    deallocate(gi)
  end function gpsgeopotential

  subroutine gpsRadii(Latitude, RadN, RadM)
    real(dp), intent(in)  :: Latitude
    real(dp), intent(out) :: RadN, RadM
    real(dp)              :: sLat, e2s

    sLat = sin(Latitude)
    e2s = 1._dp - WGS_e2 * sLat * sLat
    RadN = WGS_a / sqrt(e2s)
    RadM = WGS_a * (1._dp - WGS_e2) / (e2s*sqrt(e2s))
  end subroutine gpsRadii


!modgps03diff
  pure subroutine gpsdiffasfd(gd1, d2)
    type(gps_diff), intent(out) :: gd1
    real(dp)     , intent(in)  :: d2
    
    gd1%Var  = d2
    gd1%DVar = 0._dp
  end subroutine gpsdiffasfd

  pure subroutine gpsdiffasff(gd1, gd2)
    type(gps_diff), intent(out) :: gd1
    type(gps_diff), intent(in)  :: gd2
    
    gd1%Var  = gd2%Var
    gd1%DVar = gd2%DVar
  end subroutine gpsdiffasff

  pure function gpsdiffsmfd(gd1, d2)
    type(gps_diff), intent(in)  :: gd1
    real(dp)     , intent(in)  :: d2
    type(gps_diff)              :: gpsdiffsmfd
    
    gpsdiffsmfd%Var  = gd1%Var  + d2
    gpsdiffsmfd%DVar = gd1%DVar
  end function gpsdiffsmfd

  pure function gpsdiffsmdf(d1, gd2)
    real(dp)     , intent(in)  :: d1
    type(gps_diff), intent(in)  :: gd2
    type(gps_diff)              :: gpsdiffsmdf
    
    gpsdiffsmdf%Var  = d1 + gd2%Var
    gpsdiffsmdf%DVar =      gd2%DVar
  end function gpsdiffsmdf

  pure function gpsdiffsmfi(gd1, i2)
    type(gps_diff), intent(in)  :: gd1
    integer(i4)  , intent(in)  :: i2
    type(gps_diff)              :: gpsdiffsmfi
    
    gpsdiffsmfi%Var  = gd1%Var  + i2
    gpsdiffsmfi%DVar = gd1%DVar
  end function gpsdiffsmfi

  pure function gpsdiffsmif(i1, gd2)
    integer(i4)  , intent(in)  :: i1
    type(gps_diff), intent(in)  :: gd2
    type(gps_diff)              :: gpsdiffsmif
    
    gpsdiffsmif%Var  = i1 + gd2%Var
    gpsdiffsmif%DVar =      gd2%DVar
  end function gpsdiffsmif

  pure function gpsdiffsmff(gd1, gd2)
    type(gps_diff), intent(in)  :: gd1
    type(gps_diff), intent(in)  :: gd2
    type(gps_diff)              :: gpsdiffsmff
    
    gpsdiffsmff%Var  = gd1%Var  + gd2%Var
    gpsdiffsmff%DVar = gd1%DVar + gd2%DVar
  end function gpsdiffsmff
  
  pure function gpsdiffsbfd(gd1, d2)
    type(gps_diff), intent(in)  :: gd1
    real(dp)     , intent(in)  :: d2
    type(gps_diff)              :: gpsdiffsbfd
    
    gpsdiffsbfd%Var  = gd1%Var  - d2
    gpsdiffsbfd%DVar = gd1%DVar
  end function gpsdiffsbfd

  pure function gpsdiffsbdf(d1, gd2)
    real(dp)     , intent(in)  :: d1
    type(gps_diff), intent(in)  :: gd2
    type(gps_diff)              :: gpsdiffsbdf
    
    gpsdiffsbdf%Var  = d1 - gd2%Var
    gpsdiffsbdf%DVar =    - gd2%DVar
  end function gpsdiffsbdf

  pure function gpsdiffsbfi(gd1, i2)
    type(gps_diff), intent(in)  :: gd1
    integer(i4)  , intent(in)  :: i2
    type(gps_diff)              :: gpsdiffsbfi
    
    gpsdiffsbfi%Var  = gd1%Var  - i2
    gpsdiffsbfi%DVar = gd1%DVar
  end function gpsdiffsbfi

  pure function gpsdiffsbif(i1, gd2)
    integer(i4)  , intent(in)  :: i1
    type(gps_diff), intent(in)  :: gd2
    type(gps_diff)              :: gpsdiffsbif
    
    gpsdiffsbif%Var  = i1 - gd2%Var
    gpsdiffsbif%DVar =    - gd2%DVar
  end function gpsdiffsbif

  pure function gpsdiffsbff(gd1, gd2)
    type(gps_diff), intent(in)  :: gd1
    type(gps_diff), intent(in)  :: gd2
    type(gps_diff)              :: gpsdiffsbff
    
    gpsdiffsbff%Var  = gd1%Var  - gd2%Var
    gpsdiffsbff%DVar = gd1%DVar - gd2%DVar
  end function gpsdiffsbff

  pure function gpsdiffmlfd(gd1, d2)
    type(gps_diff), intent(in)  :: gd1
    real(dp)     , intent(in)  :: d2
    type(gps_diff)              :: gpsdiffmlfd
    
    gpsdiffmlfd%Var  = d2 * gd1%Var
    gpsdiffmlfd%DVar = d2 * gd1%DVar
  end function gpsdiffmlfd

  pure function gpsdiffmldf(d1, gd2)
    real(dp)     , intent(in)  :: d1
    type(gps_diff), intent(in)  :: gd2
    type(gps_diff)              :: gpsdiffmldf
    
    gpsdiffmldf%Var  = d1 * gd2%Var
    gpsdiffmldf%DVar = d1 * gd2%DVar
  end function gpsdiffmldf

  pure function gpsdiffmlfi(gd1, i2)
    type(gps_diff), intent(in)  :: gd1
    integer(i4)  , intent(in)  :: i2
    type(gps_diff)              :: gpsdiffmlfi
    
    gpsdiffmlfi%Var  = i2 * gd1%Var
    gpsdiffmlfi%DVar = i2 * gd1%DVar
  end function gpsdiffmlfi

  pure function gpsdiffmlif(i1, gd2)
    integer(i4)  , intent(in)  :: i1
    type(gps_diff), intent(in)  :: gd2
    type(gps_diff)              :: gpsdiffmlif
    
    gpsdiffmlif%Var  = i1 * gd2%Var
    gpsdiffmlif%DVar = i1 * gd2%DVar
  end function gpsdiffmlif

  pure function gpsdiffmlff(gd1, gd2)
    type(gps_diff), intent(in)  :: gd1
    type(gps_diff), intent(in)  :: gd2
    type(gps_diff)              :: gpsdiffmlff
    
    gpsdiffmlff%Var  = gd1%Var * gd2%Var
    gpsdiffmlff%DVar = (gd2%Var * gd1%DVar) + (gd1%Var * gd2%DVar)
  end function gpsdiffmlff

  pure function gpsdiffdvfd(gd1, d2)
    type(gps_diff), intent(in)  :: gd1
    real(dp)     , intent(in)  :: d2
    type(gps_diff)              :: gpsdiffdvfd
    
    gpsdiffdvfd%Var  = gd1%Var  / d2
    gpsdiffdvfd%DVar = gd1%DVar / d2
  end function gpsdiffdvfd

  pure function gpsdiffdvdf(d1, gd2)
    real(dp)     , intent(in)  :: d1
    type(gps_diff), intent(in)  :: gd2
    type(gps_diff)              :: gpsdiffdvdf
    
    gpsdiffdvdf%Var  =  d1 / gd2%Var
    gpsdiffdvdf%DVar = (-d1 / gd2%Var**2) * gd2%DVar
  end function gpsdiffdvdf

  pure function gpsdiffdvfi(gd1, i2)
    type(gps_diff), intent(in)  :: gd1
    integer(i4)  , intent(in)  :: i2
    type(gps_diff)              :: gpsdiffdvfi
    
    gpsdiffdvfi%Var  = gd1%Var  / i2
    gpsdiffdvfi%DVar = gd1%DVar / i2
  end function gpsdiffdvfi

  pure function gpsdiffdvif(i1, gd2)
    integer(i4)  , intent(in)  :: i1
    type(gps_diff), intent(in)  :: gd2
    type(gps_diff)              :: gpsdiffdvif
    
    gpsdiffdvif%Var  = i1 / gd2%Var
    gpsdiffdvif%DVar = (-i1 / gd2%Var**2) * gd2%DVar
  end function gpsdiffdvif

  pure function gpsdiffdvff(gd1, gd2)
    type(gps_diff), intent(in)  :: gd1
    type(gps_diff), intent(in)  :: gd2
    type(gps_diff)              :: gpsdiffdvff
    real(dp)                   :: onegd2
    
    onegd2 = 1._dp / gd2%Var
    gpsdiffdvff%Var  = gd1%Var * onegd2
    gpsdiffdvff%DVar = onegd2 * gd1%DVar - (gd1%Var*onegd2*onegd2) * gd2%DVar
  end function gpsdiffdvff

  pure function gpsdiffpwfd(gd1, d2)
    type(gps_diff), intent(in)  :: gd1
    real(dp)     , intent(in)  :: d2
    type(gps_diff)              :: gpsdiffpwfd
    
    gpsdiffpwfd%Var  = gd1%Var  ** d2
    gpsdiffpwfd%DVar = (d2*(gd1%Var**(d2-1._dp))) * gd1%DVar
  end function gpsdiffpwfd

  pure function gpsdiffpwdf(d1, gd2)
    real(dp)     , intent(in)  :: d1
    type(gps_diff), intent(in)  :: gd2
    type(gps_diff)              :: gpsdiffpwdf
    
    gpsdiffpwdf%Var  =  d1 ** gd2%Var
    gpsdiffpwdf%DVar = (log(d1)*d1**gd2%Var) * gd2%DVar
  end function gpsdiffpwdf

  pure function gpsdiffpwfi(gd1, i2)
    type(gps_diff), intent(in)  :: gd1
    integer(i4)  , intent(in)  :: i2
    type(gps_diff)              :: gpsdiffpwfi
    
    gpsdiffpwfi%Var  = gd1%Var  ** i2
    gpsdiffpwfi%DVar = (i2*(gd1%Var**(i2-1))) * gd1%DVar
  end function gpsdiffpwfi

  pure function gpsdiffpwif(i1, gd2)
    integer(i4)  , intent(in)  :: i1
    type(gps_diff), intent(in)  :: gd2
    type(gps_diff)              :: gpsdiffpwif
    
    gpsdiffpwif%Var  = i1 ** gd2%Var
    gpsdiffpwif%DVar = (log(1._dp*i1)*i1**gd2%Var) * gd2%DVar
  end function gpsdiffpwif

  pure function gpsdiffpwff(gd1, gd2)
    type(gps_diff), intent(in)  :: gd1
    type(gps_diff), intent(in)  :: gd2
    type(gps_diff)              :: gpsdiffpwff
    
    gpsdiffpwff%Var  = gd1%Var ** gd2%Var
    gpsdiffpwff%DVar = ( gd2%Var * ( gd1%Var**(gd2%Var-1) ) ) * gd1%DVar +    &
         (log(gd1%Var)*(gd1%Var**gd2%Var))*gd2%DVar
  end function gpsdiffpwff

  pure function gpsdiffsqr(gd1)
    type(gps_diff), intent(in)  :: gd1
    type(gps_diff)              :: gpsdiffsqr
    
    gpsdiffsqr%Var  = sqrt( gd1%Var )
    gpsdiffsqr%DVar = (0.5_dp / sqrt( gd1%Var )) * gd1%DVar
  end function gpsdiffsqr

  pure function gpsdiffexp(gd1)
    type(gps_diff), intent(in)  :: gd1
    type(gps_diff)              :: gpsdiffexp
    
    gpsdiffexp%Var  = exp(gd1%Var)
    gpsdiffexp%DVar = gd1%DVar * exp(gd1%Var)
  end function gpsdiffexp
  
  pure function gpsdifflog(gd1)
    type(gps_diff), intent(in)  :: gd1
    type(gps_diff)              :: gpsdifflog
    
    gpsdifflog%Var  = log(gd1%Var)
    gpsdifflog%DVar = gd1%DVar / gd1%Var
  end function gpsdifflog

  pure function gpsdiffcos(gd1)
    type(gps_diff), intent(in)  :: gd1
    type(gps_diff)              :: gpsdiffcos
    
    gpsdiffcos%Var  = cos(gd1%Var)
    gpsdiffcos%DVar = gd1%DVar * (-1._dp*sin(gd1%Var))
  end function gpsdiffcos

  pure function gpsdifftan(gd1)
    type(gps_diff), intent(in)  :: gd1
    type(gps_diff)              :: gpsdifftan
    
    gpsdifftan%Var  = tan(gd1%Var)
    gpsdifftan%DVar = (1._dp/cos(gd1%Var)**2) * gd1%DVar
  end function gpsdifftan

  pure function gpsdiffacos(gd1)
    type(gps_diff), intent(in)  :: gd1
    type(gps_diff)              :: gpsdiffacos
    
    gpsdiffacos%Var  = acos(gd1%Var)
    gpsdiffacos%DVar = gd1%DVar * (-1._dp/(1._dp-gd1%Var*gd1%Var))
  end function gpsdiffacos

  pure function gpsdiffatan(gd1)
    type(gps_diff), intent(in)  :: gd1
    type(gps_diff)              :: gpsdiffatan
    
    gpsdiffatan%Var  = atan(gd1%Var)
    gpsdiffatan%DVar = (1._dp/(1._dp+gd1%Var**2)) * gd1%DVar
  end function gpsdiffatan

  pure function gpsdifferf(gd1)
    type(gps_diff), intent(in)  :: gd1
    type(gps_diff)              :: gpsdifferf
    real(dp) , parameter :: pi = 4*ATAN(1.0) 
   ! real(dp)                   ::m_sqrtpi
    gpsdifferf%Var  = erf(gd1%Var)
    gpsdifferf%DVar = ((2._dp/sqrt(pi)) * exp(-gd1%Var**2)) * gd1%DVar
  end function gpsdifferf

!modgps04profile

  subroutine gps_struct1sw(ngpslev,rLat,rLon,rAzm,rMT,Rad,geoid,    &
       rP0,rPP,rDP,rTT,rHU,rUU,rVV,prf)
    integer(i4)     , intent(in)  :: ngpslev
    real(dp)        , intent(in)  :: rLat
    real(dp)        , intent(in)  :: rLon
    real(dp)        , intent(in)  :: rAzm
    real(dp)        , intent(in)  :: rMT
    real(dp)        , intent(in)  :: Rad
    real(dp)        , intent(in)  :: geoid
    real(dp)        , intent(in)  :: rP0
    real(dp)        , intent(in)  :: rPP (ngpssize)
    real(dp)        , intent(in)  :: rDP (ngpssize)
    real(dp)        , intent(in)  :: rTT (ngpssize)
    real(dp)        , intent(in)  :: rHU (ngpssize)
    real(dp)        , intent(in)  :: rUU (ngpssize)
    real(dp)        , intent(in)  :: rVV (ngpssize)

    type(gps_profile), intent(out) :: prf

    integer(i4)                   :: i


    real(dp), parameter           :: delta = 0.6077686814144_dp

    type(gps_diff)                 :: cmp(ngpssize)
    real(dp)                      :: h0,dh,Rgh,Eot,Eot2, sLat, cLat
    type(gps_diff)                 :: p, t, q, x
    type(gps_diff)                 :: tr, z
    type(gps_diff)                 :: mold, dd, dw, dx, n0, nd1, nw1, tvm
    type(gps_diff)                 :: xi(ngpssize), tv(ngpssize)

    prf%ngpslev = ngpslev
    prf%rLat    = rLat
    prf%rLon    = rLon
    prf%rAzm    = rAzm
    prf%rMT     = rMT
    prf%Rad     = Rad
    prf%geoid   = geoid
    call gpsRadii(rLat, prf%RadN, prf%RadM)

    !
    ! Fill pressure placeholders:
    !
    prf%P0%Var               = 0.01_dp*rP0
    prf%P0%DVar              = 0._dp
    prf%P0%DVar(2*ngpslev+1) = 0.01_dp
    do i=1,ngpslev
       prf%pst(i)%Var               = 0.01_dp*rPP(i)
       prf%pst(i)%DVar              = 0._dp
       prf%pst(i)%DVar(2*ngpslev+1) = 0.01_dp*rDP(i)
    enddo

    !
    ! Fill temperature placeholders:
    !
    do i = 1, ngpslev
       prf%tst(i)%Var               = rTT(i)+p_TC
       prf%tst(i)%DVar              = 0._dp
       prf%tst(i)%DVar(i)           = 1._dp
    enddo

    !
    ! Fill moisture placeholders:
    !
    do i = 1, ngpslev
       prf%qst(i)%Var               = rHU(i)
       prf%qst(i)%DVar              = 0._dp
       prf%qst(i)%DVar(ngpslev+i)   = 1._dp
    enddo

    ! Compressibility:
    do i = 1, ngpslev
       cmp(i)= gpscompressibility(prf%pst(i),prf%tst(i),prf%qst(i))
    enddo

    ! Refractivity:
    do i = 1, ngpslev
       p  = prf%pst(i)
       t  = prf%tst(i)
       q  = prf%qst(i)
       x  = p_wa*q/(1._dp+p_wb*q)

       ! Densities (molar, total, dry, water vapor):
       mold  = p/t * (100._dp/(p_R*cmp(i)))               ! note that p is in hPa
       dd = mold * (1._dp-x) * (p_md/1000._dp)
       dw = mold * x         * (p_mw/1000._dp)
       ! Aparicio (2011) expression
       tr = p_TC/t-1._dp
       nd1= ( 222.682_dp+   0.069_dp*tr) * dd
       nw1= (6701.605_dp+6385.886_dp*tr) * dw
       n0 = (nd1+nw1)
       prf%rst(i) = n0*(1._dp+(1.e-6_dp/6._dp)*n0)
    enddo

    !
    ! Hydrostatic equation
    !
    do i = 1, ngpslev
       p = prf%pst(i)
       t = prf%tst(i)
       q = prf%qst(i)
       !
       ! Log(P)
       !
       xi(i) = log(p)
       !
       ! Virtual temperature (K) (corrected of compressibility)
       !
       tv(i) = (1._dp+delta*q) * t * cmp(i)
    enddo

    sLat=sin(rLat)
    cLat=cos(rLat)
    dx  = xi(ngpslev)-log(prf%P0)
    Rgh = gpsgravitysrf(sLat)
    z   = (-p_Rd/Rgh) * tv(ngpslev) * dx
    prf%gst(ngpslev) = rMT + z
    do i=ngpslev-1,1,-1
       dx = xi(i)-xi(i+1)
       tvm = 0.5_dp*(tv(i)+tv(i+1))
       !
       ! Gravity acceleration (includes 2nd-order Eotvos effect)
       !
       h0  = prf%gst(i+1)%Var
       Eot = 2*WGS_OmegaPrime*cLat*p_knot*rUU(i)
       Eot2= ((p_knot*rUU(i))**2+(p_knot*rVV(i))**2)/WGS_a
       Rgh = gpsgravityalt(sLat, h0)-Eot-Eot2
       dh  = (-p_Rd/Rgh) * tvm%Var * dx%Var
       Rgh = gpsgravityalt(sLat, h0+0.5_dp*dh)-Eot-Eot2
       !
       ! Height increment
       !
       z   = (-p_Rd/Rgh) * tvm * dx
       prf%gst(i) = prf%gst(i+1) + z
    enddo

    prf%bbst=.false.
  end subroutine gps_struct1sw

  function gpscompressibility(p,t,q)
    type(gps_diff), intent(in)  :: p,t,q
    type(gps_diff)              :: gpscompressibility

    real(dp), parameter   :: a0= 1.58123e-6_dp
    real(dp), parameter   :: a1=-2.9331e-8_dp
    real(dp), parameter   :: a2= 1.1043e-10_dp
    real(dp), parameter   :: b0= 5.707e-6_dp
    real(dp), parameter   :: b1=-2.051e-8_dp
    real(dp), parameter   :: c0= 1.9898e-4_dp
    real(dp), parameter   :: c1=-2.376e-6_dp
    real(dp), parameter   :: d = 1.83e-11_dp
    real(dp), parameter   :: e =-0.765e-8_dp

    type(gps_diff)         :: x,tc,pt,tc2,x2

    x  = p_wa*q/(1._dp+p_wb*q)
    ! Estimate, from CIPM, Picard (2008)
    tc = t-p_TC
    pt = 1.e2_dp*p/t
    tc2= tc*tc
    x2 = x*x
    gpscompressibility = 1._dp-pt*(a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2)+pt*pt*(d+e*x2)
  end function gpscompressibility


!modgps04profilezd

  subroutine gps_structztd(ngpslev,rLat,rLon,rMT,rP0,rPP,rDP,rTT,rHU,lbevis,refopt,prf)
!
! This subroutine fills GPS profiles of type gps_profilezd (for ZTD operator)
!
    integer(i4)       , intent(in)  :: ngpslev          ! number of profile levels
    real(dp)          , intent(in)  :: rLat             ! radians
    real(dp)          , intent(in)  :: rLon             ! radians
    real(dp)          , intent(in)  :: rMT              ! height (ASL) of model surface (m)
    real(dp)          , intent(in)  :: rP0              ! surface pressure (Pa)
    real(dp)          , intent(in)  :: rPP (ngpssize)   ! pressure P at each level (Pa)
    real(dp)          , intent(in)  :: rDP (ngpssize)   ! dP/dP0 at each level (Pa/Pa)
    real(dp)          , intent(in)  :: rTT (ngpssize)   ! temperature T at each level (C)
    real(dp)          , intent(in)  :: rHU (ngpssize)   ! q at each level
    
!   lbevis determines which set of refractivity constants to use (Bevis or Rueger)
    logical           , intent(in)  :: lbevis

!   refopt=1 --> use conventional expression for refractivity N
!   refopt=2 --> use new Aparicio & Laroche refractivity N
    integer           , intent(in)  :: refopt

    type(gps_profilezd), intent(out) :: prf

!!        ******** PARAMETERS *************

    real(dp), parameter :: delta = 0.6077686814144_dp
    real(dp), parameter :: eps   = 0.6219800221014_dp

! Reuger (2002) refractivity constants (MKS units)
    real(dp), parameter :: k1r = 0.776890_dp
    real(dp), parameter :: k2r = 0.712952_dp
    real(dp), parameter :: k3r = 3754.63_dp
! Bevis (1994) refractivity constants (MKS units)
    real(dp), parameter :: k1b = 0.776000_dp
    real(dp), parameter :: k2b = 0.704000_dp
    real(dp), parameter :: k3b = 3739.000_dp

!    real(dp), parameter :: Avog  = 6.02214e26_dp
!    real(dp), parameter :: Boltz = 1.38065e-23_dp
!    real(dp), parameter :: mwDAir= 28.966_dp
!    real(dp), parameter :: Rd = Avog*Boltz/mwDAir    --> p_Rd
!    real(dp), parameter :: Rd = 287.05_dp            --> p_Rd
!    real(dp), parameter :: Rg = 9.80616_dp           --> p_g_GEM
    
!    real(dp), parameter :: R  = 8.314472_dp          --> p_R
!    real(dp), parameter :: md = 28.965516_dp         --> p_md
!    real(dp), parameter :: mw = 18.015254_dp         --> p_mw
!    real(dp), parameter :: wa = md/mw                --> p_wa
!    real(dp), parameter :: wb = (md-mw)/mw           --> p_wb

!!       ******** VARIABLES *************

    real(dp)            :: a0,a1,a2,b0,b1,c0,c1,d,e
    type(gps_diff)       :: tc, pt, tc2, x2, tr
    type(gps_diff)       :: mold, dd, dw, dx, n0, nd1, nw1
    integer(i4)         :: i
    real(dp)            :: k1, k2, k3, k2p
    real(dp)            :: h0, dh, Rgh, sLat, ptop
    type(gps_diff)       :: p, t, q, x, na, tvm, z
    type(gps_diff)       :: xi(ngpssize), tv(ngpssize), cmp(ngpssize), N(ngpssize) 

    prf%ngpslev = ngpslev
    prf%rLat    = rLat
    prf%rLon    = rLon
    prf%rMT     = rMT
    prf%bpst    = .false.
    !
    ! Fill pressure (P) placeholders (Pa):
    !
    prf%P0%Var               = rP0
    prf%P0%DVar              = 0._dp
    prf%P0%DVar(2*ngpslev+1) = 1._dp
    do i = 1, ngpslev
       prf%pst(i)%Var               = rPP(i)
       prf%pst(i)%DVar              = 0._dp
       prf%pst(i)%DVar(2*ngpslev+1) = rDP(i)
    enddo
    ! Pressure at model top (Pa)
    ptop = rPP(1)
    prf%bpst = .true.
    !
    ! Fill temperature (T) placeholders (C--> K):
    !
    do i = 1, ngpslev
       prf%tst(i)%Var               = rTT(i)+p_TC
       prf%tst(i)%DVar              = 0._dp
       prf%tst(i)%DVar(i)           = 1._dp
    enddo

    !
    ! Fill moisture (Q) placeholders (kg/kg):
    !
    do i = 1, ngpslev
       prf%qst(i)%Var               = rHU(i)
       prf%qst(i)%DVar              = 0._dp
       prf%qst(i)%DVar(ngpslev+i)   = 1._dp
    enddo

    if ( refopt == 2 ) then  ! use Aparicio & Laroche refractivity
    ! This code is copied from modgps04profile.cdk90
    !
    ! Compressibility:
    !
      a0 = 1.58123e-6_dp
      a1 = -2.9331e-8_dp
      a2 = 1.1043e-10_dp
      b0 = 5.707e-6_dp
      b1 = -2.051e-8_dp
      c0 = 1.9898e-4_dp
      c1 = -2.376e-6_dp
      d = 1.83e-11_dp
      e = -0.765e-8_dp
      do i = 1, ngpslev
        p  = prf%pst(i)
        t  = prf%tst(i)
        q  = prf%qst(i)
        x  = p_wa*q/(1._dp+p_wb*q)
        ! Estimate, from CIPM, Piccard (2008)
        tc = t-p_TC
        pt = p/t
        tc2 = tc*tc
        x2 = x*x
        cmp(i) = 1._dp-pt*(a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2)+pt*pt*(d+e*x2)
      enddo

    ! Refractivity:
      do i = 1, ngpslev
        p  = prf%pst(i)
        t  = prf%tst(i)
        q  = prf%qst(i)
        x  = p_wa*q/(1._dp+p_wb*q)

        ! Densities (molar, total, dry, water vapor):
        mold  = p/(p_R*t*cmp(i))
        dd = mold * (1._dp-x) * (p_md/1000._dp)
        dw = mold * x         * (p_mw/1000._dp)
        ! Aparicio (2011) expression
        tr = p_TC/t-1._dp
        nd1= ( 222.682_dp+   0.069_dp*tr) * dd
        nw1= (6701.605_dp+6385.886_dp*tr) * dw
        n0 = (nd1+nw1)
        na = n0*(1._dp+1.e-6_dp*n0/6._dp)
        N(i) = na
      enddo

    endif

    ! Refractivity constants
    if ( lbevis ) then
      k1 = k1b
      k2 = k2b
      k3 = k3b
    else
      k1 = k1r
      k2 = k2r
      k3 = k3r
    endif
    k2p = k2-(eps*k1)

    ! Virtual temperature Tv and log(P) profiles
    !
    do i = 1, ngpslev
       p = prf%pst(i)
       t = prf%tst(i)
       q = prf%qst(i)
       xi(i) = log(p)
       tv(i) = (1._dp+delta*q) * t
    enddo
    
    ! Geometric height (m) profile from lowest model level to top  --> prf%gst
    sLat = sin(rLat)
    dx  = xi(ngpslev)-log(prf%P0)
    Rgh = gpsgravitysrf(sLat)
    z   = (-p_Rd/Rgh) * tv(ngpslev) * dx
    prf%gst(ngpslev) = rMT + z
    do i = ngpslev-1, 1, -1
       dx = xi(i)-xi(i+1)
       tvm = 0.5_dp*(tv(i)+tv(i+1))
       !
       ! Gravity acceleration
       !
       h0  = prf%gst(i+1)%Var
       Rgh = gpsgravityalt(sLat, h0)
       dh  = (-p_Rd/Rgh) * tvm%Var * dx%Var
       Rgh = gpsgravityalt(sLat, h0+0.5_dp*dh)
       !
       ! Height increment (m)
       !
       z   = (-p_Rd/Rgh) * tvm * dx
       prf%gst(i) = prf%gst(i+1) + z
    enddo

    ! Profile of dZTD/dp --> prf%rst
    do i = 1, ngpslev
       p  = prf%pst(i)
       t  = prf%tst(i)
       q  = prf%qst(i)
       if ( refopt == 1 ) then
         na = (k1/tv(i)) + (k2p*(q/(eps*t))) + (k3*(q/(eps*t**2)))
       else
         na = N(i) / p
       endif
       prf%rst(i) = 1.e-6_dp * na * (p_Rd*tv(i))/gpsgravityalt(sLat, prf%gst(i)%Var)
    enddo

    ! ZTD (m) profile from model top down to lowest model level --> prf%ztd
    prf%ztd(1) = 1.e-6_dp * ((k1*p_Rd*ptop)/(gpsgravityalt(sLat, prf%gst(1)%Var)))
    do i = 2, ngpslev
      !
      ! ZTD increment = Avg(dZTD/dP) * delta_P
      !
      z = ((prf%rst(i-1) + prf%rst(i))/2._dp) * (prf%pst(i)-prf%pst(i-1))
      prf%ztd(i) = prf%ztd(i-1) + z
    enddo

  end subroutine gps_structztd

  subroutine gpsdpress(nlev,rHYB,rP0,rPT,rPR,rCF,rDP)
!  
! Computes dP/dP0 for HYBRID or ETA vertical grids
!
    integer           , intent(in)    :: nlev
    real(dp)          , intent(in)    :: rP0
    real(dp)          , intent(in)    :: rPT
    real(dp)          , intent(in)    :: rPR
    real(dp)          , intent(in)    :: rCF  ! = 1.0 for eta level grid
    real(dp)          , intent(in)    :: rHYB (ngpssize)
    real(dp)          , intent(out)   :: rDP  (ngpssize)

    integer(i4)      :: i, ngpslev
    real(dp)         :: pr1

    ngpslev = nlev
    if ( abs(rCF-1._dp) .lt. 0.01_dp ) then ! eta
        do i = 1, ngpslev
          rDP(i) = rHYB(i)
        enddo
    else                                    ! hybrid
        rDP(1) = 0._dp
        pr1 = 1._dp/(1._dp - rPT/rPR)
        do i = 2, ngpslev
          rDP(i) = ( (rHYB(i) - rPT/rPR)*pr1 )**rCF
        enddo
    endif

  end subroutine gpsdpress


!modgps05refstruct

  subroutine gpscmp(prf, cmp)
    type(gps_profile) :: prf
    type(gps_diff)   , intent(out):: cmp(:)

    integer(i4)      :: i, ngpslev
    type(gps_diff)               :: p, t, q
    !type(gps_diff)               :: Zd,Zn,Zo,Za,Zw,Zt
    type(gps_diff)               :: x,tc,pt,tc2,x2,ZtC
    real(dp)                    :: a0,a1,a2,b0,b1,c0,c1,d,e
    real(dp), parameter         :: md=28.965516_dp
    real(dp), parameter         :: mw=18.015254_dp
    real(dp), parameter         :: wa=md/mw
    real(dp), parameter         :: wb=(md-mw)/mw
    !
    a0=1.58123e-6_dp
    a1=-2.9331e-8_dp
    a2=1.1043e-10_dp
    b0=5.707e-6_dp
    b1=-2.051e-8_dp
    c0=1.9898e-4_dp
    c1=-2.376e-6_dp
    d =1.83e-11_dp
    e =-0.765e-8_dp
    !
    ngpslev = prf%ngpslev
    do i = 1, ngpslev
       p  = prf%pst(i)
       t  = prf%tst(i)
       q  = prf%qst(i)
       x  = wa*q/(1._dp+wb*q)
       ! First implementation (2007)
       !Zn=1._dp+(0.03913_dp-1.408_dp/(0.08314472_dp*t))*p/(83.14472_dp*t)
       !Zo=1._dp+(0.03183_dp-1.378_dp/(0.08314472_dp*t))*p/(83.14472_dp*t)
       !Za=1._dp+(0.03219_dp-1.363_dp/(0.08314472_dp*t))*p/(83.14472_dp*t)
       !Zw=1._dp+(0.03049_dp-5.536_dp/(0.08314472_dp*t))*p/(83.14472_dp*t)
       !Zd=0.78_dp*Zn+0.21_dp*Zo+0.01_dp*Za
       !Zt=(1._dp-q)*Zd+q*Zw
       ! Better estimate, from CIPM, Piccard (2008)
       tc = t-p_TC
       pt = 1.e2_dp*p/t
       tc2= tc*tc
       x2 = x*x
       ZtC= 1._dp-pt*(a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2)       !+pt*pt*(d+e*x2)
       ! Either choose Zt (First implementation) or ZtC (CIPM, better)
       cmp(i)=ZtC
    enddo
  end subroutine gpscmp

  subroutine gpsden(prf, den)
    type(gps_profile)              :: prf
    type(gps_diff)   , intent(out) :: den(:)
    type(gps_diff)                 :: mold, dd, dw, cmp(ngpssize)
    integer(i4)      :: i, ngpslev
    real(dp), parameter         :: R=8.314472_dp
    real(dp), parameter         :: md=28.965516_dp
    real(dp), parameter         :: mw=18.015254_dp
    real(dp), parameter         :: wa=md/mw
    real(dp), parameter         :: wb=(md-mw)/mw
    type(gps_diff)               :: p, t, q, x

    call gpscmp(prf, cmp)
    ngpslev = prf%ngpslev
    do i = 1, ngpslev
       p  = prf%pst(i)
       t  = prf%tst(i)
       q  = prf%qst(i)
       x  = wa*q/(1._dp+wb*q)

       ! Densities (molar, total, dry, water vapor):
       mold  = 100._dp*p/(R*t*cmp(i))               ! note that p is in hPa
       dd = mold * (1._dp-x) * (md/1000._dp)
       dw = mold * x         * (mw/1000._dp)
       den(i)=dd+dw
    enddo
  end subroutine gpsden


!modgps07geostruct
  
  subroutine gpsbvf(prf, bvf)
    type(gps_profile)              :: prf
    type(gps_diff)   , intent(out) :: bvf(ngpssize)

    type(gps_diff)                 :: den(ngpssize), dddz(ngpssize)
    integer(i4)                   :: i, ngpslev, im, ip
    real(dp)                      :: g, sLat

    call gpsden(prf, den)

    ngpslev = prf%ngpslev
    sLat=sin(prf%rLat)
    do i = 1, ngpslev
       ip=i+1
       im=i-1
       if (i==1)       im=1
       if (i==ngpslev) ip=ngpslev
       dddz(i)=den(i)*(log(den(ip))-log(den(im)))/(prf%gst(ip)-prf%gst(im))
       g=gpsgravityalt(sLat, prf%gst(i)%Var)
       bvf(i)=sqrt((-g)/den(i)*dddz(i))
    enddo
  end subroutine gpsbvf


!modgps08refop

  !
  ! GPSRO Refractivity operator
  ! On input:
  ! -hv       an array of height values
  ! -prf      local profile
  ! On output:
  ! -refopv   an array of refractivity values (with derivatives)
  !
  pure subroutine gps_refopv(hv, nval, prf, refopv)
    real(dp)             , intent(in) :: hv(:)
    integer(i4)          , intent(in) :: nval
    type(gps_profile)     , intent(in) :: prf
    type(gps_diff)        , intent(out):: refopv(:)
    
    integer(i4)                       :: iSize, i, ngpslev
    integer(i4)                       :: j, jloc
    real(dp)                          :: h
    
    type(gps_diff)                     :: dz

    type(gps_diff)                     :: dzm
    type(gps_diff)                     :: dzp
    
    ngpslev=prf%ngpslev
    iSize = size(hv)
    if (nval < iSize) iSize=nval
    !
    ! Given a height
    !
    do i = 1, iSize
       h = hv(i)
       !
       ! Search where it is located
       !
       if (h > prf%gst(1)%Var) then
          jloc = 1
       endif
       
       do j=1, ngpslev-1
          if ((h <= prf%gst(j)%Var) .and. (h > prf%gst(j+1)%Var)) then
             jloc = j
             exit
          endif
       enddo
       
       if (h <= prf%gst(ngpslev)%Var) then
          jloc = ngpslev-1
       endif
       !
       ! Linear-log interpolation
       !
       dz  = prf%gst(jloc) - prf%gst(jloc+1)
       
       dzm = h - prf%gst(jloc+1)
       dzp = prf%gst(jloc) - h
       
       refopv(i) = exp( (dzm * log(prf%rst(jloc)) + dzp * log(prf%rst(jloc+1))) / dz )
    enddo
  end subroutine gps_refopv

  subroutine gpshgtopv(pr, prf, hgtopv)
    real(dp)             , intent(in) :: pr
    type(gps_profile)                  :: prf
    type(gps_diff)        , intent(out):: hgtopv
    
    integer(i4)                       :: j, jloc, ngpslev
    real(dp)                          :: p
    type(gps_diff)                     :: vpm
    type(gps_diff)                     :: vpp
    
    type(gps_diff)                     :: dpr
    
    type(gps_diff)                     :: dxm
    type(gps_diff)                     :: dxp
    
    type(gps_diff)                     :: Hm
    type(gps_diff)                     :: Hp
    
    type(gps_diff)                     :: H

    ngpslev=prf%ngpslev
    !
    ! Given a pressure
    !
    p = pr
    !
    ! Search where it is located
    !
    if (p < prf%pst(1)%Var) then
       jloc = 1
    endif
    
    do j=1, ngpslev-1
       if ((p >= prf%pst(j)%Var) .and. (p < prf%pst(j+1)%Var)) then
          jloc = j
          exit
       endif
    enddo
    
    if (p >= prf%pst(ngpslev)%Var) then
       jloc = ngpslev-1
    endif
    !
    ! Find properties in that band
    !
    vpm = log(prf%pst(jloc))
    vpp = log(prf%pst(jloc+1))
    
    dpr  = vpp-vpm
    
    dxm = (vpp-log(p)) / dpr
    dxp = (log(p)-vpm) / dpr
    
    Hm  = prf%gst(jloc)
    Hp  = prf%gst(jloc+1)
    
    H   = dxm * Hm + dxp * Hp
    
    hgtopv = H
  end subroutine gpshgtopv

  subroutine gpstemopv(pr, nval, prf, temopv)
    real(dp)             , intent(in) :: pr(:)
    integer(i4)          , intent(in) :: nval
    type(gps_profile)                  :: prf
    type(gps_diff)        , intent(out):: temopv(:)
    
    integer                           :: iSize, ngpslev
    integer(i4)                       :: i, j, jloc
    real(dp)                          :: p
    type(gps_diff)                     :: vpm
    type(gps_diff)                     :: vpp
    
    type(gps_diff)                     :: dpr
    
    type(gps_diff)                     :: dxm
    type(gps_diff)                     :: dxp
    
    type(gps_diff)                     :: Tm
    type(gps_diff)                     :: Tp
    
    type(gps_diff)                     :: T

    ngpslev=prf%ngpslev
    iSize = size(pr)
    if (nval < iSize) iSize=nval
    do i = 1, iSize
       !
       ! Given a pressure
       !
       p = pr(i)
       !
       ! Search where it is located
       !
       if (p < prf%pst(1)%Var) then
          jloc = 1
       endif
    
       do j=1, ngpslev-1
          if ((p >= prf%pst(j)%Var) .and. (p < prf%pst(j+1)%Var)) then
             jloc = j
             exit
          endif
       enddo
    
       if (p >= prf%pst(ngpslev)%Var) then
          jloc = ngpslev-1
       endif
       !
       ! Find properties in that band
       !
       vpm = log(prf%pst(jloc))
       vpp = log(prf%pst(jloc+1))
       
       dpr  = vpp-vpm
    
       dxm = (vpp-log(p)) / dpr
       dxp = (log(p)-vpm) / dpr
    
       Tm  = prf%tst(jloc)
       Tp  = prf%tst(jloc+1)
    
       T   = dxm * Tm + dxp * Tp
       
       temopv(i) = T
    enddo
  end subroutine gpstemopv

  subroutine gpswmropv(pr, prf, wmropv)
    real(dp)             , intent(in) :: pr(:)
    type(gps_profile)                  :: prf
    type(gps_diff)        , intent(out):: wmropv(:)

    integer                           :: iSize, ngpslev
    integer(i4)                       :: i, j, jloc
    real(dp)                          :: p
    type(gps_diff)                     :: vpm
    type(gps_diff)                     :: vpp
    
    type(gps_diff)                     :: dpr
    
    type(gps_diff)                     :: dxm
    type(gps_diff)                     :: dxp
    
    type(gps_diff)                     :: Rm
    type(gps_diff)                     :: Rp
    
    type(gps_diff)                     :: R

    ngpslev=prf%ngpslev
    iSize = size(pr)
    do i = 1, iSize
       !
       ! Given a pressure
       !
       p = pr(i)
       !
       ! Search where it is located
       !
       if (p < prf%pst(1)%Var) then
          jloc = 1
       endif
    
       do j=1, ngpslev-1
          if ((p >= prf%pst(j)%Var) .and. (p < prf%pst(j+1)%Var)) then
             jloc = j
             exit
          endif
       enddo
    
       if (p >= prf%pst(ngpslev)%Var) then
          jloc = ngpslev-1
       endif
       !
       ! Find properties in that band
       !
       vpm = log(prf%pst(jloc))
       vpp = log(prf%pst(jloc+1))
       
       dpr  = vpp-vpm
    
       dxm = (vpp-log(p)) / dpr
       dxp = (log(p)-vpm) / dpr

       Rm  = prf%qst(jloc)
       Rp  = prf%qst(jloc+1)
    
       R   = dxm * Rm + dxp * Rp
       
       wmropv(i) = R * 28.97_dp / 18.01528_dp
    enddo
  end subroutine gpswmropv

  subroutine gpsbvfopv(hv, nval, prf, bvfopv)
    real(dp)             , intent(in) :: hv(:)
    integer(i4)          , intent(in) :: nval
    type(gps_profile)                  :: prf
    type(gps_diff)        , intent(out):: bvfopv(:)
    
    integer(i4)                       :: iSize, i, ngpslev
    integer(i4)                       :: j, jloc
    real(dp)                          :: h

    type(gps_diff)                     :: bvf(ngpssize)

    type(gps_diff)                     :: gpm
    type(gps_diff)                     :: gpp
    
    type(gps_diff)                     :: dz
    
    type(gps_diff)                     :: dxm
    type(gps_diff)                     :: dxp
    
    type(gps_diff)                     :: BVm
    type(gps_diff)                     :: BVp
    
    call gpsbvf(prf,bvf)

    ngpslev=prf%ngpslev
    iSize = size(hv)
    if (nval < iSize) iSize=nval
    !
    ! Given a height
    !
    do i = 1, iSize
       h = hv(i)
       !
       ! Search where it is located
       !
       if (h > prf%gst(1)%Var) then
          jloc = 1
       endif
       
       do j=1, ngpslev-1
          if ((h <= prf%gst(j)%Var) .and. (h > prf%gst(j+1)%Var)) then
             jloc = j
             exit
          endif
       enddo
       
       if (h <= prf%gst(ngpslev)%Var) then
          jloc = ngpslev-1
       endif
       !
       ! Find properties in that band
       !
       gpm = prf%gst(jloc)
       gpp = prf%gst(jloc+1)
       
       dz  = gpm - gpp
       
       dxm = (h-gpp) / dz
       dxp = (gpm-h) / dz
       
       BVm = bvf(jloc)
       BVp = bvf(jloc+1)
       
       bvfopv (i) = dxm * BVm + dxp * BVp
    enddo
  end subroutine gpsbvfopv


!modgps08ztdop
  !
  ! GB-GPS ZTD operator
  ! On input:
  !   -hv       height of ZTD observation Zobs (m)
  !   -prf      local model profile (type gps_profilezd)
  !   -lbevis   true/false --> use Bevis instead of Rueger k values
  !   -dzmin    Minimum DZ = Zobs-Zmod (m) for which DZ adjustment to ZTD will be made
  !             when Zobs < Zmod.
  !   -mode     1 = normal mode: use stored ZTD profiles
  !             2 = Vedel & Huang ZTD formulation: ZTD = ZHD(Pobs) + ZWD
  !                 Pobs computed from P0 using CMC hydrostatic extrapolation.
  !
  ! On output:
  !   -ZTDopv   ZTD (m) at height of observation (with derivatives)
  !   -rPobs    Pressure (Pa) at height of observation
  !
  pure subroutine gps_ztdopv(hv, prf, lbevis, dzmin, ZTDopv, rPobs, mode)
    real(dp)             , intent(in) :: hv
    type(gps_profilezd)   , intent(in) :: prf
    logical              , intent(in) :: lbevis
    real(dp)             , intent(in) :: dzmin
    type(gps_diff)        , intent(out):: ZTDopv
    real(dp)             , intent(out):: rPobs
    integer              , intent(in) :: mode
    
    integer(i4)                       :: ngpslev
    integer(i4)                       :: j, jloc
    real(dp)                          :: h, x, lat, sLat, dh
    real(dp)                          :: k1, k2, k3, k2p
    real(dp)                          :: zcon, zcon1, zconh, zfph, zconw
    
    type(gps_diff)                     :: dz, tvsfc, tobs, qobs, tvobs, naobs, Pobs
    type(gps_diff)                     :: dztddp, dztddpm
    
    type(gps_diff)                     :: zhd, tbar, qbar, qtterm, zsum, ztmobs, zqmobs
    type(gps_diff)                     :: zpbar, ztbar, zqbar, zrmean, zwd
    
    type(gps_diff)                     :: dzm, dzp

    real(dp), parameter :: delta = 0.6077686814144_dp
    real(dp), parameter :: eps   = 0.6219800221014_dp
    real(dp), parameter :: kappa = (1.0_dp/eps)-1.0_dp
    real(dp), parameter :: gamma = 0.0065_dp    ! -dT/dz (K/m)
!    real(dp), parameter :: Rg  = 9.80616_dp   --> p_g_GEM
    real(dp), parameter :: Rgm = 9.784_dp
!    real(dp), parameter :: Rd = 287.05_dp     --> p_Rd
    real(dp), parameter :: dzmax = 100.0
! Reuger (2002) refractivity constants (MKS units)
    real(dp), parameter :: k1r = 0.776890_dp
    real(dp), parameter :: k2r = 0.712952_dp
    real(dp), parameter :: k3r = 3754.63_dp
! Bevis (1994) refractivity constants (MKS units)
    real(dp), parameter :: k1b = 0.776000_dp
    real(dp), parameter :: k2b = 0.704000_dp
    real(dp), parameter :: k3b = 3739.000_dp

    ! Refractivity constants to use
    if ( lbevis ) then
      k1 = k1b
      k2 = k2b
      k3 = k3b
    else
      k1 = k1r
      k2 = k2r
      k3 = k3r
    endif
    k2p = k2-(eps*k1)

    ngpslev = prf%ngpslev
    lat     = prf%rLat
    sLat    = sin(lat)
    !
    ! Given obs height hv
    !
    h  = hv
    dh = h - prf%gst(ngpslev)%Var   ! dh = Zgps-Zmod
    !
    ! Search where it is located
    !
    do j = 1, ngpslev-1
      if ((h <= prf%gst(j)%Var) .and. (h > prf%gst(j+1)%Var)) then
        jloc = j   ! the model level above the observation
        exit
      endif
    enddo

    if (h <= prf%gst(ngpslev)%Var) then  ! obs is at or below model lowest level
      jloc = ngpslev
    endif
    
    if ( mode == 2 ) then
    
!  Compute ZTD the Vedel and Huang (2004) way: (as in old s/r gpsztdop.ftn)

       zcon  = 1.0e-06_dp*p_Rd
       zcon1 = zcon*k1
       zconw = zcon/eps
       zconh = zcon1/Rgm
       zfph = (1.0_dp - 2.66e-03_dp*cos(2.0*lat) - 2.8e-07_dp*h)

!  Pressure at obs height (CMC hydrostatic extrapolation from Psfc)
       x      = p_g_GEM/(p_Rd*gamma)
       tvsfc  = prf%tst(ngpslev)*(1._dp+delta*prf%qst(ngpslev))
       Pobs   = prf%pst(ngpslev)*(((tvsfc-gamma*dh)/tvsfc)**x)
!  Dry delay ZHD (m) at obs height
       zhd    = (zconh/zfph) * Pobs

! Integrate column q/T on pressure levels to get model ZWD
       do j = 1, ngpslev-1
         tbar = (prf%tst(j) + prf%tst(j+1))*0.5_dp
         qbar = (prf%qst(j) + prf%qst(j+1))*0.5_dp
         qtterm = ((qbar + kappa*qbar**2 )/gpsgravityalt(sLat,prf%gst(j)%Var))*(k2p + k3/tbar)
         if ( j == 1 ) then
           zsum = qtterm*(prf%pst(j+1)-prf%pst(j))
         else
           zsum = zsum + qtterm*(prf%pst(j+1)-prf%pst(j))
         endif
       enddo
       
! Compute ZWD at obs height using Higgins method (HU constant over dh layer)
       ztmobs = prf%tst(ngpslev) - (gamma * dh)
       zqmobs = prf%qst(ngpslev)
       zpbar  = (Pobs + prf%pst(ngpslev)) * 0.5_dp
       ztbar  = (ztmobs + prf%tst(ngpslev)) * 0.5_dp
       zqbar  = (zqmobs + prf%qst(ngpslev)) * 0.5_dp
!      Mean (wet) refractivity of dz layer
       zrmean = 1.0e-06_dp*(k2p*((zpbar*zqbar)/(eps*ztbar)) + k3*((zpbar*zqbar)/(eps*ztbar**2)))
       
!      Make sure adjusted ZWD >= 0
       if ( (zsum%Var*zconw)-(zrmean%Var*dh) > 0._dp ) then
         zwd = (zsum*zconw) - (zrmean*dh)
       else
         zwd = (zsum*zconw)
       endif
       
! Compute ZTD as sum of ZHD and ZWD      
       ZTDopv = zhd + zwd

    
    else   !  mode = 1: Compute ZTD using stored ZTD profile
    

      if ( jloc /= ngpslev ) then
        !
        ! Linear-log interpolation in height between levels when obs above lowest level
        !
        dz  = prf%gst(jloc) - prf%gst(jloc+1)

        dzm = h - prf%gst(jloc+1)
        dzp = prf%gst(jloc) - h

        ZTDopv = exp( (dzm*log(prf%ztd(jloc)) + dzp*log(prf%ztd(jloc+1))) / dz )
        Pobs   = exp( (dzm*log(prf%pst(jloc)) + dzp*log(prf%pst(jloc+1))) / dz )

      else   ! jloc = ngpslev ; obs is at or below model lowest level
        !
        if ( abs(dh) <= dzmin ) then  ! take lowest level values when obs is close to sfc
          ZTDopv = prf%ztd(jloc)
          Pobs   = prf%pst(jloc)
        else ! otherwise do extrapolation from lowest level values
          x      = p_g_GEM/(p_Rd*gamma)
          tvsfc  = prf%tst(jloc)*(1._dp+delta*prf%qst(jloc))
          Pobs   = prf%pst(jloc)*(((tvsfc-gamma*dh)/tvsfc)**x)
          if ( abs(dh) <= dzmax ) then
            dztddpm = prf%rst(jloc)   ! lowest level value of dZTD/dp
          else
            tobs   = prf%tst(jloc)-gamma*dh
            qobs   = prf%qst(jloc)
            tvobs  = tvsfc-gamma*dh
            naobs  = (k1/tvobs) + (k2p*(qobs/(eps*tobs))) + (k3*(qobs/(eps*tobs**2)))
            dztddp = 1.e-6_dp * naobs * (p_Rd*tvobs)/gpsgravityalt(sLat, h)
            dztddpm = (dztddp + prf%rst(jloc))/2._dp  ! mean value of dZTD/dp over dh layer
          endif
          ZTDopv = prf%ztd(jloc) + dztddpm*(Pobs-prf%pst(jloc))
        endif

      endif
    
    endif
    
    rPobs = Pobs%Var

  end subroutine gps_ztdopv

  subroutine gps_pw(prf, PW)
  !
  !  Subroutine to compute lowest level PW (kg/m2) using layer mean Q and layer delta_p (Pa)
  !
  !   Author:  S. Macpherson,  2010-2012
  !
    type(gps_profilezd)     , intent(in)  :: prf
    real(dp)               , intent(out) :: PW

    integer(i4)                       :: i, ngpslev
    real(dp)                          :: qbar, gt, gb, g, lat, sLat
    real(dp)                          :: pt, pb

    ngpslev = prf%ngpslev
    lat     = prf%rLat
    sLat    = sin(lat)

    PW = 0.0_dp

    do i = 1, ngpslev-1
      qbar = 0.5_dp * (prf%qst(i+1)%Var + prf%qst(i)%Var)
      gt  = gpsgravityalt(sLat, prf%gst(i)%Var)
      gb  = gpsgravityalt(sLat, prf%gst(i+1)%Var)
      pt  = prf%pst(i)%Var
      pb  = prf%pst(i+1)%Var
      g   = 0.5_dp * (gt + gb)
      PW = PW + (qbar/g)*(pb-pt)
    enddo

  end subroutine gps_pw


!modgps09bend

  subroutine gpsbend(prf)
    type(gps_profile)     :: prf

    type(gps_diff)                     :: sum,ta,tb,tm,trap,simp,boole,num,fa,fb,fm,nm,alpha_B
    type(gps_diff)                     :: sa,sb,sm,ra,rm,rb,dlnndra,dlnndrb,dlnndrm
    type(gps_diff)                     :: s1,s2,s3,s4,s5,r1,r2,r3,r4,r5
    type(gps_diff)                     :: nu1,nu2,nu3,nu4,nu5,n1,n2,n3,n4,n5
    type(gps_diff)                     :: t1,t2,t3,t4,t5,dlnndr1,dlnndr2,dlnndr3,dlnndr4,dlnndr5
    type(gps_diff)                     :: f1,f2,f3,f4,f5
    type(gps_diff)                     :: r  (ngpssize)
    type(gps_diff)                     :: ref(ngpssize)
    type(gps_diff)                     :: nu (ngpssize)
    type(gps_diff)                     :: lnu(ngpssize)
    type(gps_diff)                     :: n  (ngpssize)
    type(gps_diff)                     :: dlgnudr(ngpssize-1)
    type(gps_diff)                     :: rsq(ngpssize)
    type(gps_diff)                     :: nsq(ngpssize)
    type(gps_diff)                     :: x  (-ngpsxlow+1:ngpssize)
    type(gps_diff)                     :: xsq(-ngpsxlow+1:ngpssize)
    type(gps_diff)                     :: s(ngpssize),t(ngpssize)
    integer                           :: i,j,ngpslev
    logical                           :: lok

    if (.not. prf%bbst) then
       ngpslev=prf%ngpslev

       ! Radial distances and impact parameters:
       do i=1,ngpslev
          prf%dst(i)= (prf%Rad+prf%geoid+prf%gst(i))
          prf%ast(i)= prf%dst(i) * (1._dp+1.e-6_dp*prf%rst(i))
       enddo
       ! Extended lower levels:
       do i=ngpslev+1,ngpslev+ngpsxlow
          prf%ast(i)= prf%ast(i-1)-50._dp
       enddo

       ! Standard levels:
       do i=1,ngpslev
          r  (i)=prf%dst(ngpslev-i+1)
          ref(i)=prf%rst(ngpslev-i+1)
          !ref(i)=300._dp*exp((-1._dp/7000._dp)*(r(i)%Var-prf%Rad))
       enddo
       ! Extended upper levels:
       do i=ngpslev+1,ngpssize
          r  (i)=r  (i-1)+1000._dp
          ref(i)=ref(i-1)*exp(-1000._dp/7000_dp)
       enddo

       ! log n and x:
       do i=1,ngpssize
          nu(i)=1.e-6_dp*ref(i)
          lnu(i)=log(nu(i))
          n (i)=1._dp+nu(i)
          x (i)=n(i)*r(i)
          rsq(i)=r(i)**2
          nsq(i)=n(i)**2
          xsq(i)=x(i)**2
       enddo
       do i=0,-ngpsxlow+1,-1
          x  (i)=x(i+1)-50._dp
          xsq(i)=x(i)**2
       enddo

       ! Radial derivatives of log refractivity.
       ! Refractivity will be assumed exponential within each shell.
       ! We store the derivative of log(nu).
       ! dn/dr = nu * dlgnudr
       do i=1,ngpssize-1
          dlgnudr(i)=(lnu(i+1)-lnu(i))/(r(i+1)-r(i))
       enddo

       ! Evaluation of complete bending for ray tangent at r(i):
       do i=1,ngpslev
          ! Check that ray is not trapped
          lok=.true.
          do j = i+1,ngpssize
             lok= lok .and. (x(j)%Var .gt. x(i)%Var)
          enddo
          if (lok) then
             s(i)=0._dp
             t(i)=1._dp
             do j=i+1,ngpssize
                s(j)=sqrt(nsq(i)*rsq(j)-xsq(i))
                t(j)=s(j)/sqrt(xsq(j)-xsq(i))
             enddo

             ! Trapezoid integration:
             sum=0._dp
             do j=i, ngpssize-1
                sa=s(j)
                sb=s(j+1)
                ta=t(j)
                tb=t(j+1)
                dlnndra=dlgnudr(j)*nu(j  )/n(j  )
                dlnndrb=dlgnudr(j)*nu(j+1)/n(j+1)
                fa=dlnndra*ta/sqrt(xsq(i)+sa*sa)
                fb=dlnndrb*tb/sqrt(xsq(i)+sb*sb)
                sum=sum+(1._dp/2._dp)*(fa+fb)*(sb-sa)
             enddo
             trap=(-2)*r(i)*sum

             ! Simpson 1/3 integration:
             sum=0._dp
             do j=i, ngpssize-1
                sa=s(j)
                sb=s(j+1)
                sm=0.5_dp*(sa+sb)
                !
                ra=r(j)
                rb=r(j+1)
                rm=sqrt(xsq(i)+sm*sm)/n(i)
                !
                num=nu(j)*exp(dlgnudr(j)*(rm-ra))
                nm=(1._dp+num)
                !
                ta=t(j)
                tb=t(j+1)
                tm=sm/sqrt(nm*nm*rm*rm-xsq(i))
                !
                dlnndra=dlgnudr(j)*nu(j  )/n(j  )
                dlnndrb=dlgnudr(j)*nu(j+1)/n(j+1)
                dlnndrm=dlgnudr(j)*num    /nm
                !
                fa=dlnndra*ta/sqrt(xsq(i)+sa*sa)
                fb=dlnndrb*tb/sqrt(xsq(i)+sb*sb)
                fm=dlnndrm*tm/sqrt(xsq(i)+sm*sm)
                !
                sum=sum+(1._dp/6._dp)*(fa+4*fm+fb)*(sb-sa)
             enddo
             simp=(-2)*r(i)*sum

             ! Boole 2/45 integration:
             sum=0._dp
             do j=i, ngpssize-1
                s1=s(j)
                s5=s(j+1)
                s2=0.75_dp*s1+0.25_dp*s5
                s3=0.50_dp*s1+0.50_dp*s5
                s4=0.25_dp*s1+0.75_dp*s5
                !
                r1=r(j)
                r5=r(j+1)
                r2=sqrt(xsq(i)+s2*s2)/n(i)
                r3=sqrt(xsq(i)+s3*s3)/n(i)
                r4=sqrt(xsq(i)+s4*s4)/n(i)
                !
                nu1=nu(j)
                nu2=nu(j)*exp(dlgnudr(j)*(r2-r1))
                nu3=nu(j)*exp(dlgnudr(j)*(r3-r1))
                nu4=nu(j)*exp(dlgnudr(j)*(r4-r1))
                nu5=nu(j+1)
                n1=n(j)
                n2=(1._dp+nu2)
                n3=(1._dp+nu3)
                n4=(1._dp+nu4)
                n5=n(j+1)
                !
                t1=t(j)
                t2=s2/sqrt(n2*n2*r2*r2-xsq(i))
                t3=s3/sqrt(n3*n3*r3*r3-xsq(i))
                t4=s4/sqrt(n4*n4*r4*r4-xsq(i))
                t5=t(j+1)
                !
                dlnndr1=dlgnudr(j)*nu(j  )/n(j  )
                dlnndr5=dlgnudr(j)*nu(j+1)/n(j+1)
                dlnndr2=dlgnudr(j)*nu2    /n2
                dlnndr3=dlgnudr(j)*nu3    /n3
                dlnndr4=dlgnudr(j)*nu4    /n4
                !
                f1=dlnndr1*t1/sqrt(xsq(i)+s1*s1)
                f2=dlnndr2*t2/sqrt(xsq(i)+s2*s2)
                f3=dlnndr3*t3/sqrt(xsq(i)+s3*s3)
                f4=dlnndr4*t4/sqrt(xsq(i)+s4*s4)
                f5=dlnndr5*t5/sqrt(xsq(i)+s5*s5)
                !
                sum=sum+(1._dp/90._dp)*(7*f1+32*f2+12*f3+32*f4+7*f5)*(s5-s1)
             enddo
             boole=(-2)*r(i)*sum

             prf%bst(ngpslev-i+1)=boole
          else
             prf%bst(ngpslev-i+1)=-10._dp
          endif
       enddo

       ! Extended low levels:
       do i=0,-ngpsxlow+1,-1
          lok=.true.
          do j = 1,ngpssize
             lok= lok .and. (x(j)%Var .gt. x(i)%Var)
          enddo
          if (lok) then
             do j=1,ngpssize
                s(j)=sqrt(nsq(1)*rsq(j)-xsq(i))
                t(j)=s(j)/sqrt(xsq(j)-xsq(i))
             enddo

             ! Simpson integration:
             sum=0._dp
             do j=1, ngpssize-1
                sa=s(j)
                sb=s(j+1)
                sm=0.5_dp*(sa+sb)
                !
                ra=r(j)
                rb=r(j+1)
                rm=sqrt(xsq(i)+sm*sm)/n(1)
                !
                num=nu(j)*exp(dlgnudr(j)*(rm-ra))
                nm=(1._dp+num)
                !
                ta=t(j)
                tb=t(j+1)
                tm=sm/sqrt(nm*nm*rm*rm-xsq(i))
                !
                dlnndra=dlgnudr(j)*nu(j  )/n(j  )
                dlnndrb=dlgnudr(j)*nu(j+1)/n(j+1)
                dlnndrm=dlgnudr(j)*num/nm
                !
                fa=dlnndra*ta/sqrt(xsq(i)+sa*sa)
                fb=dlnndrb*tb/sqrt(xsq(i)+sb*sb)
                fm=dlnndrm*tm/sqrt(xsq(i)+sm*sm)
                !
                sum=sum+(1._dp/6._dp)*(fa+4*fm+fb)*(sb-sa)
             enddo
             simp=(-2)*(x(i)/n(1))*sum
             alpha_B=acos(x(i)/x(1))
             prf%bst(ngpslev-i+1)=simp-2*alpha_B
          else
             prf%bst(ngpslev-i+1)=-10._dp
          endif
       enddo

       prf%bbst=.true.
    endif
  end subroutine gpsbend

  subroutine gpsbend1(prf)
    type(gps_profile)     :: prf

    type(gps_diff)                     :: r  (ngpssize)
    type(gps_diff)                     :: ref(ngpssize)
    type(gps_diff)                     :: nu (ngpssize)
    type(gps_diff)                     :: lnu(ngpssize)
    type(gps_diff)                     :: n  (ngpssize)
    type(gps_diff)                     :: dlgnudr(ngpssize-1)
    type(gps_diff)                     :: x  (-ngpsxlow+1:ngpssize)

    type(gps_diff)                     :: angle0,angle,angleB,bend,nu0,th,sum,nexp
    real(dp)                          :: dxn
    integer                           :: ngpslev,i,j,jmin
    logical                           :: lok, lok2

    if (.not. prf%bbst) then
       ngpslev=prf%ngpslev

       ! Radial distances and impact parameters:
       do i=1,ngpslev
          prf%dst(i)= (prf%Rad+prf%geoid+prf%gst(i))
          prf%ast(i)= prf%dst(i) * (1._dp+1.e-6_dp*prf%rst(i))
       enddo
       ! Extended lower levels:
       do i=ngpslev+1,ngpslev+ngpsxlow
          prf%ast(i)= prf%ast(i-1)-50._dp
       enddo

       ! Standard levels:
       do i=1,ngpslev
          r  (i)=prf%dst(ngpslev-i+1)
          ref(i)=prf%rst(ngpslev-i+1)
       enddo
       ! Extended upper levels:
       do i=ngpslev+1,ngpssize
          r  (i)=r  (i-1)+1000._dp
          ref(i)=ref(i-1)*exp(-1000._dp/7000_dp)
       enddo

       ! log n and x:
       do i=1,ngpssize
          nu(i) = 1.e-6_dp*ref(i)
          lnu(i)= log(nu(i))
          n (i) = 1._dp+nu(i)
          x (i) = n(i)*r(i)
       enddo
       dxn=20._dp
       do i=0,-ngpsxlow+1,-1
          x (i) = x(i+1)-dxn
       enddo

       ! Radial derivatives of log refractivity.
       ! Refractivity will be assumed exponential within each shell.
       ! We store the derivative of log(nu).
       ! dn/dr = nu * dlgnudr
       do i=1,ngpssize-1
          dlgnudr(i)=(lnu(i+1)-lnu(i))/(r(i+1)-r(i))
       enddo

       ! Evaluation of complete bending for ray tangent at r(i):
       do i=-ngpsxlow+1,ngpslev
          lok=.true.
          lok2=.false.
          ! Check that the ray is not trapped
          ! For low impact (reflected, jmin<1) rays, begin at the surface
          jmin = i
          if (jmin < 1) jmin=1
          do j = jmin+1,ngpssize
             lok= lok .and. (x(j)%Var .gt. x(i)%Var)
          enddo
          if (lok) then
             ! Integration:
             sum=0._dp
             if (i.ge.1) then
                ! Direct rays
                angleB=0._dp
             else
                ! Reflected
                angleB=sqrt(2*(-i+1)*dxn/x(1))
             endif
             angle0=angleB
             do j=jmin, ngpssize-1
                th=r(j+1)-r(j)
                nu0=nu(j)
                nexp=dlgnudr(j)
                call gpsbendlayer(r(j), th, nu0, nexp, angle0, angle, bend, lok2)
                sum=sum+bend
                angle0=angle
             enddo
          endif
          if (lok2) then
             prf%bst(ngpslev-i+1)=(-2)*(sum+angleB)
          else
             prf%bst(ngpslev-i+1)=-10._dp
          endif
       enddo
    endif
  end subroutine gpsbend1

  subroutine gpsbendlayer(ra, th, nu0, nexp, angle0, angle, bend, lok)
    type(gps_diff)        , intent(in) :: ra, th     ! Radius of inner shell (ra) and shell thickness (th)   (m)
    type(gps_diff)        , intent(in) :: nu0, nexp  ! Refraction index coefs: n=1+nu0*exp(nexp*(r-ra)); nexp in 1/m
    type(gps_diff)        , intent(in) :: angle0     ! Ray angle above horizon at ra
    type(gps_diff)        , intent(out):: angle      ! Ray angle above horizon at rb
    type(gps_diff)        , intent(out):: bend       ! Accumulated bending over the layer
    logical              , intent(out):: lok

    type(gps_diff) :: rb,angle0i,dh,hi,rai,nu0i,anglei
    integer :: i,numunits
    
    lok=.false.
    if (th%Var.lt.0._dp) return

    ! Radius of the outer shell:
    rb = ra + th

    ! Divide layer in smaller layers:
    numunits=10
    dh =th/(1._dp*numunits)
    angle0i=angle0
    bend   =0._dp
    do i = 1, numunits
       hi =(i-1)*dh
       rai=ra+(i-1)*dh
       nu0i=nu0*exp(nexp*hi)
       call gpsbendunit(rai, dh, nu0i, nexp, angle0i, anglei, bend, lok)
       angle0i=anglei
       if (.not.lok) return
    enddo
    angle=anglei
  end subroutine gpsbendlayer

  subroutine gpsbendunit(ra, th, nu0, nexp, angle0, angle, bend, lok)
    type(gps_diff)        , intent(in) :: ra, th     ! Radius of inner shell (ra) and shell thickness (th)  (m)
    type(gps_diff)        , intent(in) :: nu0, nexp  ! Refraction index coefs: n=1+nu0*exp(nexp*(r-ra)); nexp in 1/m
    type(gps_diff)        , intent(in) :: angle0     ! Ray angle above horizon at ra
    type(gps_diff)        , intent(out):: angle      ! Ray angle above horizon at rb
    type(gps_diff)        , intent(inout):: bend       ! Accumulated bending over the layer
    logical              , intent(out):: lok

    type(gps_diff) :: rb, nu, dlnndh, g0,g1,g2,f0,f1,f2,x,a,b,c,disc,ds,bendi,g1av

    lok=.false.
    if (th%Var.lt.0._dp) return

    ! Radius of the outer shell:
    rb = ra + th

    ! Excess refraction index:
    !    at ra:    nu0
    !    at rb:    nu0*exp(nexp*h)
    nu     = nu0*exp(0.5_dp*nexp*th)
    dlnndh = nexp*nu/(1._dp+nu)

    ! Geometric trajectory:
    ! g(x) = g0+g1*x+g2*x^2
    !
    g0=0._dp
    g1=tan(angle0)
    g2=0.5_dp*dlnndh*cos(angle0)

    ! Outer circle:
    ! f(x) = f0+f1*x+f2*x^2
    !
    f0=th
    f1=0._dp
    f2=(-0.5_dp)/rb

    ! Difference:
    a=f2-g2
    b=f1-g1
    c=f0-g0

    ! Discriminant:
    disc=b*b-4*a*c
    if (disc%Var.lt.0._dp) then
       lok=.false.
       return
    else
       x =((-1)*b-sqrt(disc))/(2*a)
       g1av=g1+g2*x
       ds=x*(1._dp+(g2*x)**2)
       bendi = 2 * g2 * ds
       angle = angle0+atan(x/rb)+bendi
       bend  = bend + bendi
       if (angle%Var .gt. 0) lok=.true.
    endif
  end subroutine gpsbendunit
  subroutine j_point(a,z,ngpslev,j2)
    integer                , intent(in) :: ngpslev
    type(gps_diff)          , intent(in) :: z  (:)
    real(dp)               , intent(in) :: a  
    integer                , intent(out):: j2  
    integer                             :: j 
    j2=0
    do j=2,ngpslev
       if ((z(j-1)%Var>a) .and. (a>z(j)%Var)) then !    
          j2=j
          exit  
       endif
    enddo
  end subroutine j_point
  subroutine gpsbndopv(impv, azmv, nval, prf, bstv)
    real(dp)             , intent(in) :: impv(:), azmv(:)
    integer(i4)          , intent(in) :: nval
    type(gps_profile)                  :: prf
    type(gps_diff)        , intent(out):: bstv(:)
    
    integer                           :: iSize, i, j, ngpslev, jlocm, jlocp
    real(dp)                          :: imp1,azm1,rad, rad0
    real(dp)                          :: imp(ngpssize+ngpsxlow)
    type(gps_diff)                     :: am, ap, da, dam, dap

    call gpsbend(prf)
    ngpslev=prf%ngpslev
    iSize = size(impv)
    if (nval < iSize) iSize=nval
    rad0=prf%rad
    !
    ! Given an impact
    !
    do i = 1, iSize
       imp1 = impv(i)
       azm1 = azmv(i)
       rad=1._dp/(cos(azm1)**2/prf%radM+sin(azm1)**2/prf%radN)
       do j=1, ngpslev+ngpsxlow
          imp(j)=prf%ast(j)%Var
       enddo
       !
       ! Search where it is located
       !
       jlocm = -1000
       jlocp = -1000
       if (imp1 > imp(1)) then
          jlocm = 1
          jlocp = 2
       endif

       do j=1, ngpslev+ngpsxlow-1
          if ((imp1 <= imp(j)) .and. (abs(prf%bst(j)%Var) < 1._dp)) then
             jlocm = j
          endif
       enddo

       do j=jlocm+1, ngpslev+ngpsxlow
          if ((imp1 >  imp(j)) .and. (abs(prf%bst(j)%Var) < 1._dp)) then
             jlocp = j
             exit
          endif
       enddo
      
       if (jlocm == -1000) jlocm = ngpslev+ngpsxlow-1
       if (jlocp == -1000) jlocp = ngpslev+ngpsxlow

       !
       ! Find properties in that band
       !
       am = prf%ast(jlocm)
       ap = prf%ast(jlocp)
       
       da = am - ap
       dam = (imp1-ap) / da
       dap = (am-imp1) / da
       
       ! Use loglinear interpolation for most data (notably direct rays)
       if (prf%bst(jlocm)%Var > 1.e-6_dp .and. prf%bst(jlocp)%Var > 1.e-6_dp) then
          bstv(i)=exp(dam*log(prf%bst(jlocm))+dap*log(prf%bst(jlocp)))*(rad/rad0)
       else
          ! Use linear interpolation for near-zero or negative bending (most reflected rays)
          bstv(i)=(dam*prf%bst(jlocm)+dap*prf%bst(jlocp))*(rad/rad0)
       endif
    enddo
  end subroutine gpsbndopv

  subroutine gps_bndopv1(impv, azmv, nval, prf, bstv)
    real(dp)             , intent(in) :: impv(:), azmv(:)
    integer(i4)          , intent(in) :: nval
    type(gps_profile)                  :: prf
    type(gps_diff)        , intent(out):: bstv(:)
    
    integer                           :: iSize, i, j,j2,j3,j4, ngpslev, jlocm, jlocp
    real(dp)                          :: imp1,azm1,rad, rad0
    real(dp)                          :: imp(ngpssize+ngpsxlow)
    type(gps_diff)                     :: am, ap, da, dam, dap
    type(gps_diff)                    :: h(ngpssize), nu(ngpssize), lnu(ngpssize), n(ngpssize), z(ngpssize)
    type(gps_diff)                    :: N0a, N1a, ka, NAa,Aa, Ba, Ba2, Ba3, delta_alpha, delta_alpha_top, z_0, h_0
    real(dp)                          :: a2, a, gz(ngpssize), cazm, sazm
    real(dp)               ,parameter :: pi = 4*ATAN(1.0)


    call gpsbend1(prf)
    ngpslev=prf%ngpslev
    do i=1,ngpslev 
       h(i)  = prf%geoid+prf%gst(i)
       nu(i) = prf%rst(i)
       lnu(i)=log(nu(i)) 
       n(i)  = 1._dp+nu(i)*1e-6_dp
       z(i)  = n(i)*(prf%Rad+prf%geoid+prf%gst(i))
    enddo
    do i=1,ngpslev-1
        gz(i) = (lnu(i+1)%Var - lnu(i)%Var) / (h(i+1)%Var-  h(i)%Var)
    enddo
    
    iSize = size(impv)
    if (nval < iSize) iSize=nval
    rad0=prf%rad
    !
    ! Given an impact
    !
    do i = 1, iSize
       a2 = impv(i)*impv(i)
       a  = impv(i)
       call j_point(a,z,ngpslev,j2)
       if  (j2/=0) then
           h_0   = h(j2)+(((a-z(j2))/(z(j2-1)-z(j2)))*( h(j2-1) -h(j2)))  
           cazm  = cos(azmv(j2))
           sazm  = sin(azmv(j2)) 
           N0a   =  nu(j2)
           NAa   =  nu(j2)*exp( gz(j2)*(h_0-h(j2)))
           delta_alpha = 0.d0
           delta_alpha_top = 0.d0
           j3 = j2-1
           z_0 = a
           do while (j3>=1) 
              N1a = nu(j3)
              ka = log(N0a/N1a)/(z(j3) - z(j3+1))
              j4=j3
              do while  (ka%Var<0.or.ka%Var>0.1) 
                 j3=j3-1
                 N1a = nu(j3) 
                 ka = log(N0a/N1a)/(z(j3) - z(j4)) 
              enddo 
              Aa = 1e-6_dp* sqrt((pi/2.d0*a)/ka )*NAa* exp(ka*(z_0-a))  
              if (z_0%Var==a) then
                 Ba = erf(sqrt(ka*(z(j3)-a)))
              else  
                 Ba2 = erf(sqrt(ka*(z(j3)-a)))
                 Ba3 = erf(sqrt((ka*(z_0-a))))
                 Ba = Ba2 - Ba3
               endif 
               delta_alpha = delta_alpha + 2*ka*Ba*Aa 
               N0a = N1a
               NAa = N1a 
               z_0 = z(j3)
               j3=j3-1
            enddo
            Ba = erf(1-erf(sqrt((ka*(z_0-a)))))
            delta_alpha_top = 2*Aa*ka*Ba 
            if ((abs(delta_alpha%Var +delta_alpha_top%Var))>1) then
               j2=0
            endif
            bstv(i)   = delta_alpha +delta_alpha_top
       endif
       if (j2==0) then
          imp1 = impv(i)
          azm1 = azmv(i)
          rad=1._dp/(cos(azm1)**2/prf%radM+sin(azm1)**2/prf%radN)
          do j=1, ngpslev+ngpsxlow
             imp(j)=prf%ast(j)%Var
          enddo
          !
          ! Search where it is located
          !
          jlocm = -1000
          jlocp = -1000
          if (imp1 > imp(1)) then
             jlocm = 1
             jlocp = 2
          endif

          do j=1, ngpslev+ngpsxlow-1
             if ((imp1 <= imp(j)) .and. (abs(prf%bst(j)%Var) < 1._dp)) then
                jlocm = j
             endif
          enddo

          do j=jlocm+1, ngpslev+ngpsxlow
             if ((imp1 >  imp(j)) .and. (abs(prf%bst(j)%Var) < 1._dp)) then
                jlocp = j
             exit
             endif
          enddo
       
          if (jlocm == -1000) jlocm = ngpslev+ngpsxlow-1
          if (jlocp == -1000) jlocp = ngpslev+ngpsxlow

          !
          ! Find properties in that band
          !
          am = prf%ast(jlocm)
          ap = prf%ast(jlocp)
       
          da = am - ap
          dam = (imp1-ap) / da
          dap = (am-imp1) / da
       
          ! Use loglinear interpolation for most data (notably direct rays)
          if (prf%bst(jlocm)%Var > 1.e-6_dp .and. prf%bst(jlocp)%Var > 1.e-6_dp) then
             bstv(i)=exp(dam*log(prf%bst(jlocm))+dap*log(prf%bst(jlocp)))*(rad/rad0)
          else
          ! Use linear interpolation for near-zero or negative bending (most reflected rays)
             bstv(i)=(dam*prf%bst(jlocm)+dap*prf%bst(jlocp))*(rad/rad0)
          endif
       endif
    enddo
  end subroutine gps_bndopv1


!modgpsro_mod

  subroutine gps_setupro
    implicit none
    integer nulnam,ierr,fnom,fclos,j
!
!   Define default values:
!
    LEVELGPSRO = 2
    GPSRO_MAXPRFSIZE = 300
    SURFMIN    = 0.d0
    HSFMIN     = 0.d0
    HTPMAX     = 70000.d0
    BGCKBAND   = 0.05d0
    NUMGPSSATS = 0
!
!   Override with NML values:
!     
    nulnam=0
    ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=NAMGPSRO,iostat=ierr)
    if(ierr.ne.0) call utl_abort('gps_setupro: Error reading namelist')
    if(mpi_myid.eq.0) write(*,nml=NAMGPSRO)
    ierr=fclos(nulnam)
    if(mpi_myid.eq.0) write(*,*)'NAMGPSRO',LEVELGPSRO,GPSRO_MAXPRFSIZE,SURFMIN,HSFMIN,HTPMAX,BGCKBAND,NUMGPSSATS
!
!   Force a min/max values for the effective Fresnel widths per satellite:
!
    DO J=1,NUMGPSSATS
       IF (WGPS(J).LT. 100.d0) WGPS(J)= 100.d0
       IF (WGPS(J).GT.1500.d0) WGPS(J)=1500.d0
    ENDDO

  end subroutine gps_setupro

  integer function gps_iprofile_from_index(index)
    implicit none
    integer, intent(in) :: index
    integer i

    gps_iprofile_from_index=-1
    do i=1,size(gps_vRO_IndexPrf)
       if (index.eq.gps_vRO_IndexPrf(i)) then
          gps_iprofile_from_index=i
          return
       endif
    enddo
    return
  end function gps_iprofile_from_index


!modgpsztd_mod

  SUBROUTINE GPS_SETUPGB
!
!**s/r GPS_SETUPGB : Initialisation of ground-based GPS
!
!Author  : Stephen Macpherson *ARMA/MRD August 2008
!Revsions:
!     Stephen Macpherson  December 2012
!        -- modifcation to GB-GPS namelist parameters
!
!    -------------------
!*    Purpose: to read and initialize GB-GPS namelist parameters and print information
!*             on options selected.
!
    IMPLICIT NONE
    INTEGER J
    integer :: nulnam,ierr,fnom,fclos

!*  .  1.1 Default values
!   .      --------------

    DZMIN  = 2.0D0
    DZMAX  = 1000.0D0
    YZTDERR = 0.012D0
    LASSMET = .TRUE.
    YSFERRWGT = 1.0D0
    LLBLMET = .FALSE.
    YZDERRWGT = 1.0D0
    LBEVIS = .TRUE.
    IREFOPT = 1
    L1OBS = .FALSE.
    LTESTOP = .FALSE.
    IZTDOP = 1

    nulnam=0
    ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=NAMGPSGB,iostat=ierr)
    if(ierr.ne.0) call utl_abort('gps_setupgb: Error reading namelist')
    if(mpi_myid.eq.0) write(*,nml=NAMGPSGB)
    ierr=fclos(nulnam)

    IF (L1OBS.and.mpi_myid.eq.0) THEN
      write(*,*)' '
      write(*,*)' ******************************************'
      write(*,*)' *        GB-GPS OBSERVATIONS             *'
      write(*,*)' *                                        *'
      write(*,*)' *        ONE OBSERVATION MODE            *'
      write(*,*)' *                                        *'
      write(*,*)' ******************************************'
      write(*,*)' '
    ENDIF

!   Options to fix/adjust model ZTD to observation height and
!   assimilate GPS met data

    if(mpi_myid.eq.0) then
      write(*,*)' '
      write(*,*)' ******************************************'
      write(*,*)' *        GB-GPS OBSERVATIONS             *'
      write(*,*)' * DZ ADJUSTMENT IN gps_ztdopv IF DZ>DZMIN *'
      write(*,*)' * ZTD NOT ASSIM. IF DZ > DZMAX           *'
      write(*,*)' *                                        *'
      write(*,*)' ******************************************'
      write(*,*) ' '
      write(*,*) 'DZMIN, DZMAX = ', DZMIN, DZMAX
      write(*,*) ' '

      IF (LASSMET) THEN
        IF ( LLBLMET ) THEN
        write(*,*)' '
        write(*,*)' *****************************************'
        write(*,*)' *          GB-GPS OBSERVATIONS          *'
        write(*,*)' *     GPS MET DATA ARE ASSIMILATED      *'
        write(*,*)' *     BUT BLACKLISTED NEAR SYNO STNS    *'
        write(*,*)' *                                       *'
        write(*,*)' *****************************************'
        write(*,*) 'YSFERRWGT = ', YSFERRWGT
        write(*,*) 'YZDERRWGT = ', YZDERRWGT
        write(*,*) ' '        
        ELSE
        write(*,*)' '
        write(*,*)' *****************************************'
        write(*,*)' *          GB-GPS OBSERVATIONS          *'
        write(*,*)' *     GPS MET DATA ARE ASSIMILATED      *'
        write(*,*)' *                                       *'
        write(*,*)' *****************************************'
        write(*,*) 'YSFERRWGT = ', YSFERRWGT
        write(*,*) 'YZDERRWGT = ', YZDERRWGT
        write(*,*) ' '
        ENDIF
      ELSE
        write(*,*)' '
        write(*,*)' *****************************************'
        write(*,*)' *          GB-GPS OBSERVATIONS          *'
        write(*,*)' *   GPS MET DATA ARE NOT ASSIMILATED    *'
        write(*,*)' *                                       *'
        write(*,*)' *****************************************'
        write(*,*) 'YZDERRWGT = ', YZDERRWGT
        write(*,*) ' '
      ENDIF

      IF (YZTDERR .LT. 0.0D0) THEN
        write(*,*)' '
        write(*,*)' *****************************************'
        write(*,*)' *          GB-GPS OBSERVATIONS          *'
        write(*,*)' *    ZTD OBSERVATION ERROR FROM FERR    *'
        write(*,*)' *                                       *'
        write(*,*)' *****************************************'
      ELSE IF (YZTDERR .GT. 0.0D0) THEN
        write(*,*)' '
        write(*,*)' *****************************************'
        write(*,*)' *          GB-GPS OBSERVATIONS          *'
        write(*,*)' *     ZTD OBSERVATION ERROR IS FIXED    *'
        write(*,*)' *                                       *'
        write(*,*)' *****************************************'
        write(*,*)' '
        write(*,*)'YZTDERR (mm) = ', YZTDERR*1000.D0
      ELSE
        write(*,*)' '
        write(*,*)' *****************************************'
        write(*,*)' *          GB-GPS OBSERVATIONS          *'
        write(*,*)' *   ZTD OBSERVATION ERROR IS FROM ZWD   *'
        write(*,*)' *   USING SD(O-P) STATS (REGRESSION)    *'
        write(*,*)' *                                       *'
        write(*,*)' *****************************************'
        write(*,*)' '
      ENDIF

      IF (IREFOPT .EQ. 1) THEN
        IF (LBEVIS) THEN
        write(*,*)' '
        write(*,*)' *****************************************'
        write(*,*)' *          GB-GPS OBSERVATIONS          *'
        write(*,*)' *                                       *'
        write(*,*)' *  CONVENTIONAL REFACTIVITY N USING     *'
        write(*,*)' *  BEVIS 92 K1, K2, K3 TO COMPUTE ZTD   *'
        write(*,*)' *****************************************'
        write(*,*)' ' 
        ELSE
        write(*,*)' '
        write(*,*)' *****************************************'
        write(*,*)' *          GB-GPS OBSERVATIONS          *'
        write(*,*)' *                                       *'
        write(*,*)' *  CONVENTIONAL REFACTIVITY N USING     *'
        write(*,*)' *  RUEGER 02 K1, K2, K3 TO COMPUTE ZTD  *'
        write(*,*)' *****************************************'
        write(*,*)' ' 
        ENDIF
        IF (IZTDOP .EQ. 1) THEN
        write(*,*)' '
        write(*,*)' *****************************************'
        write(*,*)' *          GB-GPS OBSERVATIONS          *'
        write(*,*)' *                                       *'
        write(*,*)' *   NORMAL ZTD OPERATOR -- ZTD COMPUTED *'
        write(*,*)' *           FROM ZTD(K) PROFILE         *'
        write(*,*)' *****************************************'
        write(*,*)' ' 
        ELSE
        write(*,*)' '
        write(*,*)' *****************************************'
        write(*,*)' *          GB-GPS OBSERVATIONS          *'
        write(*,*)' *                                       *'
        write(*,*)' *   ORIGINAL OPERATOR -- ZTD = ZHD+ZWD  *'
        write(*,*)' *        VEDEL AND HUANG (2004)         *'
        write(*,*)' *****************************************'
        write(*,*)' ' 
        ENDIF        
      ELSE
        write(*,*)' '
        write(*,*)' *****************************************'
        write(*,*)' *          GB-GPS OBSERVATIONS          *'
        write(*,*)' *                                       *'
        write(*,*)' *  APARICIO & LAROCHE REFRACTIVITY N    *'
        write(*,*)' *         USED TO COMPUTE ZTD           *'
        write(*,*)' *****************************************'
        write(*,*)' '       
      ENDIF

    endif

  END subroutine gps_setupgb

  integer function gps_i_from_index(index)
    integer, intent(in) :: index
    integer i

    gps_i_from_index = -1
    do i = 1, size(vGPSZTD_Index)
       if (index .eq. vGPSZTD_Index(i)) then
          gps_i_from_index = i
          return
       endif
    enddo
    return
  end function gps_i_from_index


end module gps_mod

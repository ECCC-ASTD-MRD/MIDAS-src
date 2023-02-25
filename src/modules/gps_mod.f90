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

module gps_mod
  ! MODULE gps_mod (prefix='gps' category='5. Observation operators')
  !
  ! :Purpose: Code related to GPS-RO and ground-based GPS observation operators.
  !
  ! :Note: prefix not used for all public variables
  !
  use midasMpi_mod
  use utilities_mod
  use mathPhysConstants_mod
  use earthConstants_mod
  implicit none
  save
  private

  ! public types
  public :: gps_profile, gps_profilezd, gps_diff

  ! public variables
  public :: gps_numROProfiles, gps_vRO_IndexPrf
  public :: gps_Level_RO, gps_RO_MAXPRFSIZE, gps_SURFMIN, gps_HSFMIN, gps_HTPMAX, gps_HTPMAXER, gps_BGCKBAND, gps_WGPS
  public :: gps_roError, gps_roBNorm
  public :: gps_gravitysrf, gps_gb_maxdata, gps_gb_dzmin
  public :: gps_gb_ltestop, gps_gb_llblmet, gps_gb_lbevis, gps_gb_irefopt, gps_gb_iztdop, gps_gb_lassmet, gps_gb_l1obs, gps_gb_yzderrwgt, gps_gb_numztd
  public :: gps_ZTD_Index, gps_ncvmx, gps_gb_dzmax, gps_gb_yztderr, gps_gb_ysferrwgt

  ! public procedures
  public :: gps_setupro, gps_iprofile_from_index
  public :: gps_setupgb, gps_iztd_from_index
  public :: gps_struct1sw, gps_struct1sw_v2, gps_bndopv1, gps_refopv, gps_structztd_v2, gps_ztdopv, gps_pw
  public :: gps_geopotential

  ! public constants
  integer, parameter, public :: gps_Level_RO_Bnd       = 1
  integer, parameter, public :: gps_Level_RO_Ref       = 2
  integer, parameter, public :: gps_Level_RO_BndandRef = 3
  public :: gps_p_md, gps_p_mw, gps_p_wa, gps_p_wb

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
  integer(i4), parameter :: gps_ncvmx  = 4*ngpssize

!modgps01ctphys
  
  ! Avogadro constant:
  real(dp), parameter           :: p_Avog  = 6.02214129e23_dp        ! From CODATA

  ! Boltzmann constant:
  real(dp), parameter           :: p_Boltz = 1.3806488e-23_dp        ! From CODATA

  ! Air properties (public):
  real(dp), parameter           :: gps_p_md    = 28.965516_dp            ! From Aparicio(2011)
  real(dp), parameter           :: gps_p_mw    = 18.015254_dp            ! From Aparicio(2011)
  real(dp), parameter           :: gps_p_wa    = gps_p_md/gps_p_mw
  real(dp), parameter           :: gps_p_wb    = (gps_p_md-gps_p_mw)/gps_p_mw

  ! Gas constants:
  real(dp), parameter           :: p_R     = p_Avog*p_Boltz          ! per mol
  real(dp), parameter           :: p_Rd    = p_Avog*p_Boltz/(1.e-3_dp*gps_p_md)   ! per air mass

!modgps03diff

  type gps_diff
     real(dp)           :: Var
     real(dp)           :: DVar(gps_ncvmx)
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
     integer(i4)                                  :: ngpslev
     real(dp)                                     :: rLat
     real(dp)                                     :: rLon
     real(dp)                                     :: rAzm
     real(dp)                                     :: rMT
     real(dp)                                     :: Rad
     real(dp)                                     :: geoid
     real(dp)                                     :: RadN
     real(dp)                                     :: RadM

     type(gps_diff)                               :: P0

     type(gps_diff), dimension(ngpssize)          :: pst
     type(gps_diff), dimension(ngpssize)          :: tst
     type(gps_diff), dimension(ngpssize)          :: qst
     type(gps_diff), dimension(ngpssize)          :: rst
     type(gps_diff), dimension(ngpssize)          :: gst

     logical                                      :: bbst
     type(gps_diff), dimension(ngpssize)          :: dst
     type(gps_diff), dimension(ngpssize+ngpsxlow) :: ast
     type(gps_diff), dimension(ngpssize+ngpsxlow) :: bst
  end type gps_profile

!modgps04profilezd
  
  type gps_profilezd
     integer(i4)                                  :: ngpslev
     real(dp)                                     :: rLat
     real(dp)                                     :: rLon
     real(dp)                                     :: rMT

     type(gps_diff)                               :: P0
     
     type(gps_diff), dimension(ngpssize)          :: pst
     type(gps_diff), dimension(ngpssize)          :: tst
     type(gps_diff), dimension(ngpssize)          :: qst
     type(gps_diff), dimension(ngpssize)          :: rst
     type(gps_diff), dimension(ngpssize)          :: gst
     type(gps_diff), dimension(ngpssize)          :: ztd
     logical                                      :: bpst
  end type gps_profilezd

!modgpsro_mod

!
! Values determined by input data:
!
  integer                                :: gps_numROProfiles
  integer         , allocatable          :: gps_vRO_IndexPrf(:,:)   ! index for each profile

  ! Public versions of namelist variables
  INTEGER gps_Level_RO, gps_RO_MAXPRFSIZE
  REAL*8  gps_SurfMin, gps_HsfMin, gps_HtpMax, gps_BgckBand, gps_HtpMaxEr
  REAL*4  gps_Wgps(0:1023,4)
  character(len=20) :: gps_roError
  LOGICAL :: gps_roBNorm, gps_gpsroEotvos


!modgpsztd_mod

  integer, parameter      ::  max_gps_sites = 1200
  integer, parameter      ::  gps_gb_maxdata  = max_gps_sites*24     ! (max_gps_sites) * (max_num_obs in 6h)

  integer                 :: gps_gb_numZTD            ! number of ZTD data to be assimilated
  integer , allocatable   :: gps_ZTD_Index (:)        ! INDEX_HEADER in CMA (ObsSpace) for each ZTD observation

  ! Public versions of namelist variables
  REAL*8 gps_gb_DZMIN, gps_gb_YZTDERR, gps_gb_YSFERRWGT, gps_gb_YZDERRWGT
  REAL(8) :: gps_gb_DZMAX = 1000.d0 ! need to give it a default value here in case setup not called
  INTEGER gps_gb_IREFOPT, gps_gb_IZTDOP
  LOGICAL gps_gb_LASSMET, gps_gb_LLBLMET, gps_gb_LBEVIS, gps_gb_L1OBS, gps_gb_LTESTOP

contains

!modgps02wgs84grav

  pure function gps_gravitysrf(sLat)
    !
    !:Purpose: Normal gravity on ellipsoidal surface
    !
    implicit none

    real(dp)              :: gps_gravitysrf ! Normal gravity (m/s2)

    ! Arguments:
    real(dp), intent(in)  :: sLat ! sin(Latitude)

    ! Locals:
    real(dp)              :: ks2
    real(dp)              :: e2s

    ks2 = ec_wgs_TNGk * sLat*sLat
    e2s = 1._dp - ec_wgs_e2 * sLat*sLat
    gps_gravitysrf = ec_wgs_GammaE * (1._dp + ks2) / sqrt(e2s)
  end function gps_gravitysrf

  pure function gps_gravityalt(sLat, Altitude)
    !
    !:Purpose: Normal gravity above the ellipsoidal surface
    !
    implicit none

    real(dp)              :: gps_gravityalt ! Normal gravity (m/s2)

    ! Arguments:
    real(dp), intent(in)  :: sLat     ! sin(Latitude)
    real(dp), intent(in)  :: Altitude ! Altitude (m)

    real(dp)              :: C1
    real(dp)              :: C2

    C1 =-2._dp/ec_wgs_a*(1._dp+ec_wgs_f+ec_wgs_m-2*ec_wgs_f*sLat*sLat)
    C2 = 3._dp/ec_wgs_a**2
    gps_gravityalt = gps_gravitysrf(sLat)*                                   &
         (1._dp + C1 * Altitude + C2 * Altitude**2)
  end function gps_gravityalt

  pure function gps_geopotential(Latitude, Altitude)
    !
    !:Purpose: Geopotential energy at a given point.
    !          Result is based on the WGS84 approximate expression for the
    !          gravity acceleration as a function of latitude and altitude,
    !          integrated with the trapezoidal rule.
    !
    implicit none

    real(dp)              :: gps_geopotential ! Geopotential (m2/s2)

    ! Arguments:
    real(dp), intent(in)  :: Latitude ! (rad)
    real(dp), intent(in)  :: Altitude ! (m)

    ! Locals:
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
       gi(i) = gps_gravityalt(sLat, hi(i))
    end do
    hi(n) = Altitude
    gi(n) = gps_gravityalt(sLat, hi(n))

    gps_geopotential = 0._dp
    do i = 1, n
       gps_geopotential = gps_geopotential + 0.5_dp * (gi(i)+gi(i-1)) * (hi(i)-hi(i-1))
    end do

    deallocate(hi)
    deallocate(gi)
  end function gps_geopotential

  subroutine gpsRadii(Latitude, RadN, RadM)
    implicit none

    ! Arguments:
    real(dp), intent(in)  :: Latitude
    real(dp), intent(out) :: RadN, RadM

    ! Locals:
    real(dp)              :: sLat, e2s

    sLat = sin(Latitude)
    e2s = 1._dp - ec_wgs_e2 * sLat * sLat
    RadN = ec_wgs_a / sqrt(e2s)
    RadM = ec_wgs_a * (1._dp - ec_wgs_e2) / (e2s*sqrt(e2s))
  end subroutine gpsRadii


!modgps03diff
  pure subroutine gpsdiffasfd(gd1, d2)
    implicit none
    
    ! Arguments:
    type(gps_diff), intent(out) :: gd1
    real(dp)     , intent(in)  :: d2
    
    gd1%Var  = d2
    gd1%DVar = 0._dp
  end subroutine gpsdiffasfd

  pure subroutine gpsdiffasff(gd1, gd2)
    implicit none

    ! Arguments:
    type(gps_diff), intent(out) :: gd1
    type(gps_diff), intent(in)  :: gd2
    
    gd1%Var  = gd2%Var
    gd1%DVar = gd2%DVar
  end subroutine gpsdiffasff

  pure function gpsdiffsmfd(gd1, d2)
    implicit none

    type(gps_diff)              :: gpsdiffsmfd

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    real(dp)      , intent(in)  :: d2
    
    gpsdiffsmfd%Var  = gd1%Var  + d2
    gpsdiffsmfd%DVar = gd1%DVar
  end function gpsdiffsmfd

  pure function gpsdiffsmdf(d1, gd2)
    implicit none

    type(gps_diff)              :: gpsdiffsmdf

    ! Arguments:
    real(dp)      , intent(in)  :: d1
    type(gps_diff), intent(in)  :: gd2
    
    gpsdiffsmdf%Var  = d1 + gd2%Var
    gpsdiffsmdf%DVar =      gd2%DVar
  end function gpsdiffsmdf

  pure function gpsdiffsmfi(gd1, i2)
    implicit none

    type(gps_diff)              :: gpsdiffsmfi

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    integer(i4)   , intent(in)  :: i2
    
    gpsdiffsmfi%Var  = gd1%Var  + i2
    gpsdiffsmfi%DVar = gd1%DVar
  end function gpsdiffsmfi

  pure function gpsdiffsmif(i1, gd2)
    implicit none
   
    type(gps_diff)              :: gpsdiffsmif

    ! Arguments:
    integer(i4)   , intent(in)  :: i1
    type(gps_diff), intent(in)  :: gd2
    
    gpsdiffsmif%Var  = i1 + gd2%Var
    gpsdiffsmif%DVar =      gd2%DVar
  end function gpsdiffsmif

  pure function gpsdiffsmff(gd1, gd2)
    implicit none
   
    type(gps_diff)              :: gpsdiffsmff

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    type(gps_diff), intent(in)  :: gd2
    
    gpsdiffsmff%Var  = gd1%Var  + gd2%Var
    gpsdiffsmff%DVar = gd1%DVar + gd2%DVar
  end function gpsdiffsmff
  
  pure function gpsdiffsbfd(gd1, d2)
    implicit none 

    type(gps_diff)              :: gpsdiffsbfd

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    real(dp)      , intent(in)  :: d2
    
    gpsdiffsbfd%Var  = gd1%Var  - d2
    gpsdiffsbfd%DVar = gd1%DVar
  end function gpsdiffsbfd

  pure function gpsdiffsbdf(d1, gd2)
    implicit none

    type(gps_diff)              :: gpsdiffsbdf

    ! Arguments:
    real(dp)      , intent(in)  :: d1
    type(gps_diff), intent(in)  :: gd2
    
    gpsdiffsbdf%Var  = d1 - gd2%Var
    gpsdiffsbdf%DVar =    - gd2%DVar
  end function gpsdiffsbdf

  pure function gpsdiffsbfi(gd1, i2)
    implicit none

    type(gps_diff)              :: gpsdiffsbfi

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    integer(i4)   , intent(in)  :: i2
    
    gpsdiffsbfi%Var  = gd1%Var  - i2
    gpsdiffsbfi%DVar = gd1%DVar
  end function gpsdiffsbfi

  pure function gpsdiffsbif(i1, gd2)
    implicit none

    type(gps_diff)              :: gpsdiffsbif

    ! Arguments:
    integer(i4)   , intent(in)  :: i1
    type(gps_diff), intent(in)  :: gd2
    
    gpsdiffsbif%Var  = i1 - gd2%Var
    gpsdiffsbif%DVar =    - gd2%DVar
  end function gpsdiffsbif

  pure function gpsdiffsbff(gd1, gd2)
    implicit none

    type(gps_diff)              :: gpsdiffsbff

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    type(gps_diff), intent(in)  :: gd2
    
    gpsdiffsbff%Var  = gd1%Var  - gd2%Var
    gpsdiffsbff%DVar = gd1%DVar - gd2%DVar
  end function gpsdiffsbff

  pure function gpsdiffmlfd(gd1, d2)
    implicit none

    type(gps_diff)              :: gpsdiffmlfd

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    real(dp)      , intent(in)  :: d2
    
    gpsdiffmlfd%Var  = d2 * gd1%Var
    gpsdiffmlfd%DVar = d2 * gd1%DVar
  end function gpsdiffmlfd

  pure function gpsdiffmldf(d1, gd2)
    implicit none

    type(gps_diff)              :: gpsdiffmldf

    ! Arguments:
    real(dp)      , intent(in)  :: d1
    type(gps_diff), intent(in)  :: gd2
    
    gpsdiffmldf%Var  = d1 * gd2%Var
    gpsdiffmldf%DVar = d1 * gd2%DVar
  end function gpsdiffmldf

  pure function gpsdiffmlfi(gd1, i2)
    implicit none

    type(gps_diff)              :: gpsdiffmlfi

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    integer(i4)   , intent(in)  :: i2
    
    gpsdiffmlfi%Var  = i2 * gd1%Var
    gpsdiffmlfi%DVar = i2 * gd1%DVar
  end function gpsdiffmlfi

  pure function gpsdiffmlif(i1, gd2)
    implicit none

    type(gps_diff)              :: gpsdiffmlif

    ! Arguments:
    integer(i4)   , intent(in)  :: i1
    type(gps_diff), intent(in)  :: gd2
    
    gpsdiffmlif%Var  = i1 * gd2%Var
    gpsdiffmlif%DVar = i1 * gd2%DVar
  end function gpsdiffmlif

  pure function gpsdiffmlff(gd1, gd2)
    implicit none

    type(gps_diff)              :: gpsdiffmlff

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    type(gps_diff), intent(in)  :: gd2
    
    gpsdiffmlff%Var  = gd1%Var * gd2%Var
    gpsdiffmlff%DVar = (gd2%Var * gd1%DVar) + (gd1%Var * gd2%DVar)
  end function gpsdiffmlff

  pure function gpsdiffdvfd(gd1, d2)
    implicit none

    type(gps_diff)              :: gpsdiffdvfd

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    real(dp)      , intent(in)  :: d2
    
    gpsdiffdvfd%Var  = gd1%Var  / d2
    gpsdiffdvfd%DVar = gd1%DVar / d2
  end function gpsdiffdvfd

  pure function gpsdiffdvdf(d1, gd2)
    implicit none

    type(gps_diff)              :: gpsdiffdvdf

    ! Arguments:
    real(dp)      , intent(in)  :: d1
    type(gps_diff), intent(in)  :: gd2
    
    gpsdiffdvdf%Var  =  d1 / gd2%Var
    gpsdiffdvdf%DVar = (-d1 / gd2%Var**2) * gd2%DVar
  end function gpsdiffdvdf

  pure function gpsdiffdvfi(gd1, i2)
    implicit none

    type(gps_diff)              :: gpsdiffdvfi

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    integer(i4)   , intent(in)  :: i2
    
    gpsdiffdvfi%Var  = gd1%Var  / i2
    gpsdiffdvfi%DVar = gd1%DVar / i2
  end function gpsdiffdvfi

  pure function gpsdiffdvif(i1, gd2)
    implicit none

    type(gps_diff)              :: gpsdiffdvif

    ! Arguments:
    integer(i4)   , intent(in)  :: i1
    type(gps_diff), intent(in)  :: gd2
    
    gpsdiffdvif%Var  = i1 / gd2%Var
    gpsdiffdvif%DVar = (-i1 / gd2%Var**2) * gd2%DVar
  end function gpsdiffdvif

  pure function gpsdiffdvff(gd1, gd2)
    implicit none

    type(gps_diff)              :: gpsdiffdvff

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    type(gps_diff), intent(in)  :: gd2

    ! Locals:
    real(dp)                    :: onegd2
    
    onegd2 = 1._dp / gd2%Var
    gpsdiffdvff%Var  = gd1%Var * onegd2
    gpsdiffdvff%DVar = onegd2 * gd1%DVar - (gd1%Var*onegd2*onegd2) * gd2%DVar
  end function gpsdiffdvff

  pure function gpsdiffpwfd(gd1, d2)
    implicit none

    type(gps_diff)              :: gpsdiffpwfd

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    real(dp)      , intent(in)  :: d2
    
    gpsdiffpwfd%Var  = gd1%Var  ** d2
    gpsdiffpwfd%DVar = (d2*(gd1%Var**(d2-1._dp))) * gd1%DVar
  end function gpsdiffpwfd

  pure function gpsdiffpwdf(d1, gd2)
    implicit none
    
    type(gps_diff)              :: gpsdiffpwdf

    ! Arguments:
    real(dp)      , intent(in)  :: d1
    type(gps_diff), intent(in)  :: gd2
    
    gpsdiffpwdf%Var  =  d1 ** gd2%Var
    gpsdiffpwdf%DVar = (log(d1)*d1**gd2%Var) * gd2%DVar
  end function gpsdiffpwdf

  pure function gpsdiffpwfi(gd1, i2)
    implicit none

    type(gps_diff)              :: gpsdiffpwfi

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    integer(i4)   , intent(in)  :: i2
    
    gpsdiffpwfi%Var  = gd1%Var  ** i2
    gpsdiffpwfi%DVar = (i2*(gd1%Var**(i2-1))) * gd1%DVar
  end function gpsdiffpwfi

  pure function gpsdiffpwif(i1, gd2)
    implicit none

    type(gps_diff)              :: gpsdiffpwif

    ! Arguments:
    integer(i4)   , intent(in)  :: i1
    type(gps_diff), intent(in)  :: gd2
    
    gpsdiffpwif%Var  = i1 ** gd2%Var
    gpsdiffpwif%DVar = (log(1._dp*i1)*i1**gd2%Var) * gd2%DVar
  end function gpsdiffpwif

  pure function gpsdiffpwff(gd1, gd2)
    implicit none

    type(gps_diff)              :: gpsdiffpwff

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    type(gps_diff), intent(in)  :: gd2
    
    gpsdiffpwff%Var  = gd1%Var ** gd2%Var
    gpsdiffpwff%DVar = ( gd2%Var * ( gd1%Var**(gd2%Var-1) ) ) * gd1%DVar +    &
         (log(gd1%Var)*(gd1%Var**gd2%Var))*gd2%DVar
  end function gpsdiffpwff

  pure function gpsdiffsqr(gd1)
    implicit none

    type(gps_diff)              :: gpsdiffsqr

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    
    gpsdiffsqr%Var  = sqrt( gd1%Var )
    gpsdiffsqr%DVar = (0.5_dp / sqrt( gd1%Var )) * gd1%DVar
  end function gpsdiffsqr

  pure function gpsdiffexp(gd1)
    implicit none

    type(gps_diff)              :: gpsdiffexp

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    
    gpsdiffexp%Var  = exp(gd1%Var)
    gpsdiffexp%DVar = gd1%DVar * exp(gd1%Var)
  end function gpsdiffexp
  
  pure function gpsdifflog(gd1)
    implicit none

    type(gps_diff)              :: gpsdifflog

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    
    gpsdifflog%Var  = log(gd1%Var)
    gpsdifflog%DVar = gd1%DVar / gd1%Var
  end function gpsdifflog

  pure function gpsdiffcos(gd1)
    implicit none

    type(gps_diff)              :: gpsdiffcos

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    
    gpsdiffcos%Var  = cos(gd1%Var)
    gpsdiffcos%DVar = gd1%DVar * (-1._dp*sin(gd1%Var))
  end function gpsdiffcos

  pure function gpsdifftan(gd1)
    implicit none

    type(gps_diff)              :: gpsdifftan

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    
    gpsdifftan%Var  = tan(gd1%Var)
    gpsdifftan%DVar = (1._dp/cos(gd1%Var)**2) * gd1%DVar
  end function gpsdifftan

  pure function gpsdiffacos(gd1)
    implicit none

    type(gps_diff)              :: gpsdiffacos

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    
    gpsdiffacos%Var  = acos(gd1%Var)
    gpsdiffacos%DVar = gd1%DVar * (-1._dp/(1._dp-gd1%Var*gd1%Var))
  end function gpsdiffacos

  pure function gpsdiffatan(gd1)
    implicit none

    type(gps_diff)              :: gpsdiffatan

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1
    
    gpsdiffatan%Var  = atan(gd1%Var)
    gpsdiffatan%DVar = (1._dp/(1._dp+gd1%Var**2)) * gd1%DVar
  end function gpsdiffatan

  pure function gpsdifferf(gd1)
    implicit none

    type(gps_diff)              :: gpsdifferf

    ! Arguments:
    type(gps_diff), intent(in)  :: gd1

    ! Locals:
    real(dp) , parameter :: pi = MPC_PI_R8
    ! real(dp)                   ::m_sqrtpi
    gpsdifferf%Var  = erf(gd1%Var)
    gpsdifferf%DVar = ((2._dp/sqrt(pi)) * exp(-gd1%Var**2)) * gd1%DVar
  end function gpsdifferf

!modgps04profile

  subroutine gps_struct1sw(ngpslev,rLat,rLon,rAzm,rMT,Rad,geoid,    &
       rP0,rPP,rDP,rTT,rHU,rUU,rVV,prf,printHeight)
    implicit none
 
    ! Arguments:
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
    logical         , optional    :: printHeight

    type(gps_profile), intent(out) :: prf

    ! Locals:
    integer(i4)                    :: i

    real(dp) , parameter           :: delta = 0.6077686814144_dp

    type(gps_diff)                 :: cmp(ngpssize)
    real(dp)                       :: h0,dh,Rgh,Eot,Eot2, sLat, cLat
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
    end do

    !
    ! Fill temperature placeholders:
    !
    do i = 1, ngpslev
       prf%tst(i)%Var               = rTT(i)+MPC_K_C_DEGREE_OFFSET_R8
       prf%tst(i)%DVar              = 0._dp
       prf%tst(i)%DVar(i)           = 1._dp
    end do

    !
    ! Fill moisture placeholders:
    !
    do i = 1, ngpslev
       prf%qst(i)%Var               = rHU(i)
       prf%qst(i)%DVar              = 0._dp
       prf%qst(i)%DVar(ngpslev+i)   = 1._dp
    end do

    ! Compressibility:
    do i = 1, ngpslev
       cmp(i)= gpscompressibility(prf%pst(i),prf%tst(i),prf%qst(i))
    end do

    ! Refractivity:
    do i = 1, ngpslev
       p  = prf%pst(i)
       t  = prf%tst(i)
       q  = prf%qst(i)
       x  = gps_p_wa*q/(1._dp+gps_p_wb*q)

       ! Densities (molar, total, dry, water vapor):
       mold  = p/t * (100._dp/(p_R*cmp(i)))               ! note that p is in hPa
       dd = mold * (1._dp-x) * (gps_p_md/1000._dp)
       dw = mold * x         * (gps_p_mw/1000._dp)
       ! Aparicio (2011) expression
       tr = MPC_K_C_DEGREE_OFFSET_R8/t-1._dp
       nd1= ( 222.682_dp+   0.069_dp*tr) * dd
       nw1= (6701.605_dp+6385.886_dp*tr) * dw
       n0 = (nd1+nw1)
       prf%rst(i) = n0*(1._dp+(1.e-6_dp/6._dp)*n0)
    end do

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
    end do

    sLat=sin(rLat)
    cLat=cos(rLat)
    dx  = xi(ngpslev)-log(prf%P0)
    Rgh = gps_gravitysrf(sLat)
    z   = (-p_Rd/Rgh) * tv(ngpslev) * dx
    prf%gst(ngpslev) = rMT + z
    do i=ngpslev-1,1,-1
       dx = xi(i)-xi(i+1)
       tvm = 0.5_dp*(tv(i)+tv(i+1))
       !
       ! Gravity acceleration (includes 2nd-order Eotvos effect)
       !
       h0  = prf%gst(i+1)%Var
       Eot = 2*ec_wgs_OmegaPrime*cLat*rUU(i)
       Eot2= (rUU(i)**2+rVV(i)**2)/ec_wgs_a
       Rgh = gps_gravityalt(sLat, h0)-Eot-Eot2
       dh  = (-p_Rd/Rgh) * tvm%Var * dx%Var
       Rgh = gps_gravityalt(sLat, h0+0.5_dp*dh)-Eot-Eot2
       !
       ! Height increment
       !
       z   = (-p_Rd/Rgh) * tvm * dx
       prf%gst(i) = prf%gst(i+1) + z
    end do

    if ( present(printHeight) ) then
      if ( printHeight ) then
        write(*,*) 'gps_struct1sw, height='
        write(*,*) prf%gst(1:ngpslev)%Var

        printHeight = .false.
      end if
    end if

    prf%bbst=.false.
  end subroutine gps_struct1sw

  subroutine gps_struct1sw_v2(ngpslev,rLat,rLon,rAzm,rMT,Rad,geoid,    &
       rP0,rPP,rTT,rHU,rUU,rVV,rALT,prf)
    implicit none

    ! Arguments:
    integer(i4)     , intent(in)  :: ngpslev
    real(dp)        , intent(in)  :: rLat
    real(dp)        , intent(in)  :: rLon
    real(dp)        , intent(in)  :: rAzm
    real(dp)        , intent(in)  :: rMT
    real(dp)        , intent(in)  :: Rad
    real(dp)        , intent(in)  :: geoid
    real(dp)        , intent(in)  :: rP0
    real(dp)        , intent(in)  :: rPP (ngpssize)
    real(dp)        , intent(in)  :: rTT (ngpssize)
    real(dp)        , intent(in)  :: rHU (ngpssize)
    real(dp)        , intent(in)  :: rUU (ngpssize)
    real(dp)        , intent(in)  :: rVV (ngpssize)
    real(dp)        , intent(in)  :: rALT (ngpssize)

    type(gps_profile), intent(out) :: prf

    ! Locals
    integer(i4)                   :: i
    real(dp)                      :: rALT_E(ngpssize)


    real(dp) , parameter           :: delta = 0.6077686814144_dp

    type(gps_diff)                 :: cmp(ngpssize)
    type(gps_diff)                 :: p, t, q, x
    type(gps_diff)                 :: tr
    type(gps_diff)                 :: mold, dd, dw, n0, nd1, nw1

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
    prf%P0%DVar(4*ngpslev)   = 0.01_dp
    do i=1,ngpslev
       prf%pst(i)%Var               = 0.01_dp*rPP(i)
       prf%pst(i)%DVar              = 0._dp
       prf%pst(i)%DVar(3*ngpslev+i) = 0.01_dp
    end do

    !
    ! Fill temperature placeholders:
    !
    do i = 1, ngpslev
       prf%tst(i)%Var               = rTT(i)+MPC_K_C_DEGREE_OFFSET_R8
       prf%tst(i)%DVar              = 0._dp
       prf%tst(i)%DVar(i)           = 1._dp
    end do

    !
    ! Fill moisture placeholders:
    !
    do i = 1, ngpslev
       prf%qst(i)%Var               = rHU(i)
       prf%qst(i)%DVar              = 0._dp
       prf%qst(i)%DVar(ngpslev+i)   = 1._dp
    end do

    !
    ! Fill altitude placeholders:
    !
    if (gps_gpsroEotvos) then
      call gpsro_Eotvos_dH(ngpslev, rLat, rALT, rUU, rVV, rALT_E)
    else
      rALT_E(1:ngpslev) = rALT(1:ngpslev)
    end if
    do i = 1, ngpslev
       prf%gst(i)%Var                 = rALT_E(i)
       prf%gst(i)%DVar                = 0._dp
       prf%gst(i)%DVar(2*ngpslev+i)   = 1._dp
    end do

    ! Compressibility:
    do i = 1, ngpslev
       cmp(i)= gpscompressibility(prf%pst(i),prf%tst(i),prf%qst(i))
    end do

    ! Refractivity:
    do i = 1, ngpslev
       p  = prf%pst(i)
       t  = prf%tst(i)
       q  = prf%qst(i)
       x  = gps_p_wa*q/(1._dp+gps_p_wb*q)

       ! Densities (molar, total, dry, water vapor):
       mold  = p/t * (100._dp/(p_R*cmp(i)))               ! note that p is in hPa
       dd = mold * (1._dp-x) * (gps_p_md/1000._dp)
       dw = mold * x         * (gps_p_mw/1000._dp)
       ! Aparicio (2011) expression
       tr = MPC_K_C_DEGREE_OFFSET_R8/t-1._dp
       nd1= ( 222.682_dp+   0.069_dp*tr) * dd
       nw1= (6701.605_dp+6385.886_dp*tr) * dw
       n0 = (nd1+nw1)
       prf%rst(i) = n0*(1._dp+(1.e-6_dp/6._dp)*n0)
    end do

    prf%bbst=.false.
  end subroutine gps_struct1sw_v2

  function gpscompressibility(p,t,q)
    implicit none

    type(gps_diff)              :: gpscompressibility

    ! Arguments:
    type(gps_diff), intent(in)  :: p,t,q

    ! Locals:
    real(dp) , parameter   :: a0= 1.58123e-6_dp
    real(dp) , parameter   :: a1=-2.9331e-8_dp
    real(dp) , parameter   :: a2= 1.1043e-10_dp
    real(dp) , parameter   :: b0= 5.707e-6_dp
    real(dp) , parameter   :: b1=-2.051e-8_dp
    real(dp) , parameter   :: c0= 1.9898e-4_dp
    real(dp) , parameter   :: c1=-2.376e-6_dp
    real(dp) , parameter   :: d = 1.83e-11_dp
    real(dp) , parameter   :: e =-0.765e-8_dp

    type(gps_diff)         :: x,tc,pt,tc2,x2

    x  = gps_p_wa*q/(1._dp+gps_p_wb*q)
    ! Estimate, from CIPM, Picard (2008)
    tc = t-MPC_K_C_DEGREE_OFFSET_R8
    pt = 1.e2_dp*p/t
    tc2= tc*tc
    x2 = x*x
    gpscompressibility = 1._dp-pt*(a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2)+pt*pt*(d+e*x2)
  end function gpscompressibility

  subroutine gpsro_Eotvos_dH(ngpslev, rLat, rALT, rUU, rVV, rALT_E)
    implicit none

    ! Arguments:
    integer,  intent(in)   :: ngpslev
    real(dp), intent(in)   :: rLat
    real(dp), intent(in)   :: rALT(ngpslev)
    real(dp), intent(in)   :: rUU (ngpslev)
    real(dp), intent(in)   :: rVV (ngpslev)
    real(dp), intent(out)  :: rALT_E(ngpslev)

    ! Locals
    integer                :: i
    real(dp)               :: cLat, dALT, Eot, Eot2, dALTE, ddAL, acc

    cLat=cos(rLat)
    rALT_E(ngpslev) = rALT(ngpslev)
    acc = 0.d0
    do i = ngpslev-1, 1, -1
      dALT = rALT(i) - rALT(i+1)
      Eot = 2*ec_wgs_OmegaPrime*cLat*rUU(i)
      Eot2= (rUU(i)**2+rVV(i)**2)/ec_wgs_a
      dALTE = dALT*(1.d0+(Eot+Eot2)/ec_rg)
      ddAL = dALTE - dALT
      acc = acc + ddAL
      rALT_E(i) = rALT(i) + acc
      !write(*,'(A15,I4,8F15.8)')'EOTVOS shift', i, rALT(i), rALT_E(i), dALT, Eot, Eot2, ddAL, acc
    end do
  end subroutine gpsro_Eotvos_dH

!modgps04profilezd

  subroutine gps_structztd(ngpslev,rLat,rLon,rMT,rP0,rPP,rDP,rTT,rHU,lbevis,&
                           refopt,prf)
    !
    !:Purpose: This subroutine fills GPS profiles of type gps_profilezd (for ZTD
    !          operator)
    !
    !:Arguments:
    !     :refopt:
    !               =1 --> use conventional expression for refractivity N
    !
    !               =2 --> use new Aparicio & Laroche refractivity N
    implicit none
    
    ! Arguments:
    integer(i4)       , intent(in)  :: ngpslev          ! number of profile levels
    real(dp)          , intent(in)  :: rLat             ! radians
    real(dp)          , intent(in)  :: rLon             ! radians
    real(dp)          , intent(in)  :: rMT              ! height (ASL) of model surface (m)
    real(dp)          , intent(in)  :: rP0              ! surface pressure (Pa)
    real(dp)          , intent(in)  :: rPP (ngpssize)   ! pressure P at each level (Pa)
    real(dp)          , intent(in)  :: rDP (ngpssize)   ! dP/dP0 at each level (Pa/Pa)
    real(dp)          , intent(in)  :: rTT (ngpssize)   ! temperature T at each level (C)
    real(dp)          , intent(in)  :: rHU (ngpssize)   ! q at each level
    logical           , intent(in)  :: lbevis ! determines which set of refractivity constants to use (Bevis or Rueger)
    integer           , intent(in)  :: refopt

    type(gps_profilezd), intent(out) :: prf


    ! Locals:

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
!    real(dp), parameter :: md = 28.965516_dp         --> gps_p_md
!    real(dp), parameter :: mw = 18.015254_dp         --> gps_p_mw
!    real(dp), parameter :: wa = md/mw                --> gps_p_wa
!    real(dp), parameter :: wb = (md-mw)/mw           --> gps_p_wb

!!       ******** VARIABLES *************

    real(dp)             :: a0,a1,a2,b0,b1,c0,c1,d,e
    type(gps_diff)       :: tc, pt, tc2, x2, tr
    type(gps_diff)       :: mold, dd, dw, dx, n0, nd1, nw1
    integer(i4)         :: i
    real(dp)             :: k1, k2, k3, k2p
    real(dp)             :: h0, dh, Rgh, sLat, ptop
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
    end do
    ! Pressure at model top (Pa)
    ptop = rPP(1)
    prf%bpst = .true.
    !
    ! Fill temperature (T) placeholders (C--> K):
    !
    do i = 1, ngpslev
       prf%tst(i)%Var               = rTT(i)+MPC_K_C_DEGREE_OFFSET_R8
       prf%tst(i)%DVar              = 0._dp
       prf%tst(i)%DVar(i)           = 1._dp
    end do

    !
    ! Fill moisture (Q) placeholders (kg/kg):
    !
    do i = 1, ngpslev
       prf%qst(i)%Var               = rHU(i)
       prf%qst(i)%DVar              = 0._dp
       prf%qst(i)%DVar(ngpslev+i)   = 1._dp
    end do

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
        x  = gps_p_wa*q/(1._dp+gps_p_wb*q)
        ! Estimate, from CIPM, Piccard (2008)
        tc = t-MPC_K_C_DEGREE_OFFSET_R8
        pt = p/t
        tc2 = tc*tc
        x2 = x*x
        cmp(i) = 1._dp-pt*(a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2)+pt*pt*(d+e*x2)
      end do

    ! Refractivity:
      do i = 1, ngpslev
        p  = prf%pst(i)
        t  = prf%tst(i)
        q  = prf%qst(i)
        x  = gps_p_wa*q/(1._dp+gps_p_wb*q)

        ! Densities (molar, total, dry, water vapor):
        mold  = p/(p_R*t*cmp(i))
        dd = mold * (1._dp-x) * (gps_p_md/1000._dp)
        dw = mold * x         * (gps_p_mw/1000._dp)
        ! Aparicio (2011) expression
        tr = MPC_K_C_DEGREE_OFFSET_R8/t-1._dp
        nd1= ( 222.682_dp+   0.069_dp*tr) * dd
        nw1= (6701.605_dp+6385.886_dp*tr) * dw
        n0 = (nd1+nw1)
        na = n0*(1._dp+1.e-6_dp*n0/6._dp)
        N(i) = na
      end do

    end if

    ! Refractivity constants
    if ( lbevis ) then
      k1 = k1b
      k2 = k2b
      k3 = k3b
    else
      k1 = k1r
      k2 = k2r
      k3 = k3r
    end if
    k2p = k2-(eps*k1)

    ! Virtual temperature Tv and log(P) profiles
    !
    do i = 1, ngpslev
       p = prf%pst(i)
       t = prf%tst(i)
       q = prf%qst(i)
       xi(i) = log(p)
       tv(i) = (1._dp+delta*q) * t
    end do
    
    ! Geometric height (m) profile from lowest model level to top  --> prf%gst
    sLat = sin(rLat)
    dx  = xi(ngpslev)-log(prf%P0)
    Rgh = gps_gravitysrf(sLat)
    z   = (-p_Rd/Rgh) * tv(ngpslev) * dx
    prf%gst(ngpslev) = rMT + z
    do i = ngpslev-1, 1, -1
       dx = xi(i)-xi(i+1)
       tvm = 0.5_dp*(tv(i)+tv(i+1))
       !
       ! Gravity acceleration
       !
       h0  = prf%gst(i+1)%Var
       Rgh = gps_gravityalt(sLat, h0)
       dh  = (-p_Rd/Rgh) * tvm%Var * dx%Var
       Rgh = gps_gravityalt(sLat, h0+0.5_dp*dh)
       !
       ! Height increment (m)
       !
       z   = (-p_Rd/Rgh) * tvm * dx
       prf%gst(i) = prf%gst(i+1) + z
    end do

    ! Profile of dZTD/dp --> prf%rst
    do i = 1, ngpslev
       p  = prf%pst(i)
       t  = prf%tst(i)
       q  = prf%qst(i)
       if ( refopt == 1 ) then
         na = (k1/tv(i)) + (k2p*(q/(eps*t))) + (k3*(q/(eps*t**2)))
       else
         na = N(i) / p
       end if
       prf%rst(i) = 1.e-6_dp * na * (p_Rd*tv(i))/gps_gravityalt(sLat, prf%gst(i)%Var)
    end do

    ! ZTD (m) profile from model top down to lowest model level --> prf%ztd
    prf%ztd(1) = 1.e-6_dp * ((k1*p_Rd*ptop)/(gps_gravityalt(sLat, prf%gst(1)%Var)))
    do i = 2, ngpslev
      !
      ! ZTD increment = Avg(dZTD/dP) * delta_P
      !
      z = ((prf%rst(i-1) + prf%rst(i))/2._dp) * (prf%pst(i)-prf%pst(i-1))
      prf%ztd(i) = prf%ztd(i-1) + z
    end do

  end subroutine gps_structztd

  subroutine gps_structztd_v2(ngpslev,rLat,rLon,rMT,rP0,rPP,rTT,rHU,rALT,&
                              lbevis,refopt,prf)
    !
    !:Purpose: This subroutine fills GPS profiles of type gps_profilezd (for ZTD
    !          operator)
    !
    !:Arguments:
    !     :refopt:  =1 --> use conventional expression for refractivity N
    !
    !               =2 --> use new Aparicio & Laroche refractivity N
    implicit none

    ! Arguments:
    integer(i4)       , intent(in)  :: ngpslev          ! number of profile levels
    real(dp)          , intent(in)  :: rLat             ! radians
    real(dp)          , intent(in)  :: rLon             ! radians
    real(dp)          , intent(in)  :: rMT              ! height (ASL) of model surface (m)
    real(dp)          , intent(in)  :: rP0              ! surface pressure (Pa)
    real(dp)          , intent(in)  :: rPP (ngpssize)   ! pressure P at each level (Pa)
    real(dp)          , intent(in)  :: rTT (ngpssize)   ! temperature T at each level (C)
    real(dp)          , intent(in)  :: rHU (ngpssize)   ! q at each level
    real(dp)          , intent(in)  :: rALT (ngpssize)   ! altitude at each level
    logical           , intent(in)  :: lbevis ! determines which set of refractivity constants to use (Bevis or Rueger)
    integer           , intent(in)  :: refopt

    type(gps_profilezd), intent(out) :: prf


    ! Locals:
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
!    real(dp), parameter :: md = 28.965516_dp         --> gps_p_md
!    real(dp), parameter :: mw = 18.015254_dp         --> gps_p_mw
!    real(dp), parameter :: wa = md/mw                --> gps_p_wa
!    real(dp), parameter :: wb = (md-mw)/mw           --> gps_p_wb

!!       ******** VARIABLES *************

    real(dp)             :: a0,a1,a2,b0,b1,c0,c1,d,e
    type(gps_diff)       :: tc, pt, tc2, x2, tr
    type(gps_diff)       :: mold, dd, dw, n0, nd1, nw1
    integer(i4)         :: i
    real(dp)             :: k1, k2, k3, k2p
    real(dp)             :: sLat, ptop
    type(gps_diff)       :: p, t, q, x, na, z
    type(gps_diff)       :: tv(ngpssize), cmp(ngpssize), N(ngpssize) 

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
    prf%P0%DVar(4*ngpslev)   = 1._dp
    do i = 1, ngpslev
       prf%pst(i)%Var               = rPP(i)
       prf%pst(i)%DVar              = 0._dp
       prf%pst(i)%DVar(3*ngpslev+i) = 1._dp
    end do
    ! Pressure at model top (Pa)
    ptop = rPP(1)
    prf%bpst = .true.
    !
    ! Fill temperature (T) placeholders (C--> K):
    !
    do i = 1, ngpslev
       prf%tst(i)%Var               = rTT(i)+MPC_K_C_DEGREE_OFFSET_R8
       prf%tst(i)%DVar              = 0._dp
       prf%tst(i)%DVar(i)           = 1._dp
    end do

    !
    ! Fill moisture (Q) placeholders (kg/kg):
    !
    do i = 1, ngpslev
       prf%qst(i)%Var               = rHU(i)
       prf%qst(i)%DVar              = 0._dp
       prf%qst(i)%DVar(ngpslev+i)   = 1._dp
    end do

    !
    ! Fill altitude (AL) placeholders (m):
    !
    do i = 1, ngpslev
       prf%gst(i)%Var               = rALT(i)
       prf%gst(i)%DVar              = 0._dp
       prf%gst(i)%DVar(2*ngpslev+i) = 1._dp
    end do

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
        x  = gps_p_wa*q/(1._dp+gps_p_wb*q)
        ! Estimate, from CIPM, Piccard (2008)
        tc = t-MPC_K_C_DEGREE_OFFSET_R8
        pt = p/t
        tc2 = tc*tc
        x2 = x*x
        cmp(i) = 1._dp-pt*(a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2)+pt*pt*(d+e*x2)
      end do

    ! Refractivity:
      do i = 1, ngpslev
        p  = prf%pst(i)
        t  = prf%tst(i)
        q  = prf%qst(i)
        x  = gps_p_wa*q/(1._dp+gps_p_wb*q)

        ! Densities (molar, total, dry, water vapor):
        mold  = p/(p_R*t*cmp(i))
        dd = mold * (1._dp-x) * (gps_p_md/1000._dp)
        dw = mold * x         * (gps_p_mw/1000._dp)
        ! Aparicio (2011) expression
        tr = MPC_K_C_DEGREE_OFFSET_R8/t-1._dp
        nd1= ( 222.682_dp+   0.069_dp*tr) * dd
        nw1= (6701.605_dp+6385.886_dp*tr) * dw
        n0 = (nd1+nw1)
        na = n0*(1._dp+1.e-6_dp*n0/6._dp)
        N(i) = na
      end do

    end if

    ! Refractivity constants
    if ( lbevis ) then
      k1 = k1b
      k2 = k2b
      k3 = k3b
    else
      k1 = k1r
      k2 = k2r
      k3 = k3r
    end if
    k2p = k2-(eps*k1)

    ! Virtual temperature Tv and log(P) profiles
    !
    do i = 1, ngpslev
       p = prf%pst(i)
       t = prf%tst(i)
       q = prf%qst(i)
       tv(i) = (1._dp+delta*q) * t
    end do

    sLat = sin(rLat)

    ! Profile of dZTD/dp --> prf%rst
    do i = 1, ngpslev
       p  = prf%pst(i)
       t  = prf%tst(i)
       q  = prf%qst(i)
       if ( refopt == 1 ) then
         na = (k1/tv(i)) + (k2p*(q/(eps*t))) + (k3*(q/(eps*t**2)))
       else
         na = N(i) / p
       end if
       prf%rst(i) = 1.e-6_dp * na * (p_Rd*tv(i))/gps_gravityalt(sLat, prf%gst(i)%Var)
    end do

    ! ZTD (m) profile from model top down to lowest model level --> prf%ztd
    prf%ztd(1) = 1.e-6_dp * ((k1*p_Rd*ptop)/(gps_gravityalt(sLat, prf%gst(1)%Var)))
    do i = 2, ngpslev
      !
      ! ZTD increment = Avg(dZTD/dP) * delta_P
      !
      z = ((prf%rst(i-1) + prf%rst(i))/2._dp) * (prf%pst(i)-prf%pst(i-1))
      prf%ztd(i) = prf%ztd(i-1) + z
    end do

  end subroutine gps_structztd_v2

!modgps05refstruct

  subroutine gpscmp(prf, cmp)
    implicit none

    ! Arguments:
    type(gps_profile) :: prf
    type(gps_diff)   , intent(out):: cmp(:)

    ! Locals:
    integer(i4)      :: i, ngpslev
    type(gps_diff)               :: p, t, q
    !type(gps_diff)               :: Zd,Zn,Zo,Za,Zw,Zt
    type(gps_diff)               :: x,tc,pt,tc2,x2,ZtC
    real(dp)                     :: a0,a1,a2,b0,b1,c0,c1,d,e
    real(dp) , parameter         :: md=28.965516_dp
    real(dp) , parameter         :: mw=18.015254_dp
    real(dp) , parameter         :: wa=md/mw
    real(dp) , parameter         :: wb=(md-mw)/mw
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
       tc = t-MPC_K_C_DEGREE_OFFSET_R8
       pt = 1.e2_dp*p/t
       tc2= tc*tc
       x2 = x*x
       ZtC= 1._dp-pt*(a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2)       !+pt*pt*(d+e*x2)
       ! Either choose Zt (First implementation) or ZtC (CIPM, better)
       cmp(i)=ZtC
    end do
  end subroutine gpscmp

  subroutine gpsden(prf, den)
    implicit none

    ! Arguments:
    type(gps_profile)              :: prf
    type(gps_diff)   , intent(out) :: den(:)

    ! Locals:
    type(gps_diff)                 :: mold, dd, dw, cmp(ngpssize)
    integer(i4)      :: i, ngpslev
    real(dp) , parameter         :: R=8.314472_dp
    real(dp) , parameter         :: md=28.965516_dp
    real(dp) , parameter         :: mw=18.015254_dp
    real(dp) , parameter         :: wa=md/mw
    real(dp) , parameter         :: wb=(md-mw)/mw
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
    end do
  end subroutine gpsden


!modgps07geostruct
  
  subroutine gpsbvf(prf, bvf)
    implicit none
 
    ! Arguments:
    type(gps_profile)              :: prf
    type(gps_diff)   , intent(out) :: bvf(ngpssize)

    ! Locals:
    type(gps_diff)                 :: den(ngpssize), dddz(ngpssize)
    integer(i4)                    :: i, ngpslev, im, ip
    real(dp)                       :: g, sLat

    call gpsden(prf, den)

    ngpslev = prf%ngpslev
    sLat=sin(prf%rLat)
    do i = 1, ngpslev
       ip=i+1
       im=i-1
       if (i==1)       im=1
       if (i==ngpslev) ip=ngpslev
       dddz(i)=den(i)*(log(den(ip))-log(den(im)))/(prf%gst(ip)-prf%gst(im))
       g=gps_gravityalt(sLat, prf%gst(i)%Var)
       bvf(i)=sqrt((-g)/den(i)*dddz(i))
    end do
  end subroutine gpsbvf


!modgps08refop

  pure subroutine gps_refopv(hv, nval, prf, refopv)
    !
    !:Purpose: GPSRO Refractivity operator
    !
    implicit none

    ! Arguments:
    real(dp)              , intent(in) :: hv(:) ! an array of height values
    integer(i4)           , intent(in) :: nval
    type(gps_profile)     , intent(in) :: prf   ! local profile
    type(gps_diff)        , intent(out):: refopv(:) ! an array of refractivity values (with derivatives)

    ! Locals:
    integer(i4)                       :: iSize, i, ngpslev
    integer(i4)                       :: j, jloc
    real(dp)                           :: h
    
    type(gps_diff)                     :: dz

    type(gps_diff)                     :: dzm
    type(gps_diff)                     :: dzp
    
    ngpslev=prf%ngpslev
    iSize = size(hv)
    if (nval < iSize) iSize=nval
    !
    ! Given a height
    !
    jloc = 1
    do i = 1, iSize
       h = hv(i)
       !
       ! Search where it is located
       !
       if (h > prf%gst(1)%Var) then
          jloc = 1
       end if
       
       do j=1, ngpslev-1
          if ((h <= prf%gst(j)%Var) .and. (h > prf%gst(j+1)%Var)) then
             jloc = j
             exit
          end if
       end do
       
       if (h <= prf%gst(ngpslev)%Var) then
          jloc = ngpslev-1
       end if
       !
       ! Interpolation/extrapolation
       !
       if (h >= prf%gst(ngpslev)%Var) then
          !
          ! Either linear-log interpolation
          !
          dz  = prf%gst(jloc) - prf%gst(jloc+1)
       
          dzm = h - prf%gst(jloc+1)
          dzp = prf%gst(jloc) - h
       
          refopv(i) = exp( (dzm * log(prf%rst(jloc)) + dzp * log(prf%rst(jloc+1))) / dz )
       else
          !
          ! Or exp extrapolation at the lower edge
          ! (better standard exp profile than linear-log, which may be unstable)
          !
          dzm = h - prf%gst(jloc+1)
          refopv(i) = prf%rst(jloc+1) * exp((-1._dp/6500._dp)*dzm)
       end if
    end do
  end subroutine gps_refopv

  subroutine gpshgtopv(pr, prf, hgtopv)
    implicit none
 
    ! Arguments:
    real(dp)              , intent(in) :: pr
    type(gps_profile)                  :: prf
    type(gps_diff)        , intent(out):: hgtopv

    ! Locals:
    integer(i4)                       :: j, jloc, ngpslev
    real(dp)                           :: p
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
    end if
    
    do j=1, ngpslev-1
       if ((p >= prf%pst(j)%Var) .and. (p < prf%pst(j+1)%Var)) then
          jloc = j
          exit
       end if
    end do
    
    if (p >= prf%pst(ngpslev)%Var) then
       jloc = ngpslev-1
    end if
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
    implicit none

    ! Arguments:
    real(dp)              , intent(in) :: pr(:)
    integer(i4)           , intent(in) :: nval
    type(gps_profile)                  :: prf
    type(gps_diff)        , intent(out):: temopv(:)

    ! Locals:
    integer                           :: iSize, ngpslev
    integer(i4)                        :: i, j, jloc
    real(dp)                           :: p
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
       end if
    
       do j=1, ngpslev-1
          if ((p >= prf%pst(j)%Var) .and. (p < prf%pst(j+1)%Var)) then
             jloc = j
             exit
          end if
       end do
    
       if (p >= prf%pst(ngpslev)%Var) then
          jloc = ngpslev-1
       end if
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
    end do
  end subroutine gpstemopv

  subroutine gpswmropv(pr, prf, wmropv)
    implicit none

    ! Arguments:
    real(dp)              , intent(in) :: pr(:)
    type(gps_profile)                  :: prf
    type(gps_diff)        , intent(out):: wmropv(:)

    ! Locals:
    integer                           :: iSize, ngpslev
    integer(i4)                        :: i, j, jloc
    real(dp)                           :: p
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
       end if
    
       do j=1, ngpslev-1
          if ((p >= prf%pst(j)%Var) .and. (p < prf%pst(j+1)%Var)) then
             jloc = j
             exit
          end if
       end do
    
       if (p >= prf%pst(ngpslev)%Var) then
          jloc = ngpslev-1
       end if
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
    end do
  end subroutine gpswmropv

  subroutine gpsbvfopv(hv, nval, prf, bvfopv)
    implicit none

    ! Arguments:
    real(dp)              , intent(in) :: hv(:)
    integer(i4)           , intent(in) :: nval
    type(gps_profile)                  :: prf
    type(gps_diff)        , intent(out):: bvfopv(:)

    ! Locals:
    integer(i4)                        :: iSize, i, ngpslev
    integer(i4)                        :: j, jloc
    real(dp)                           :: h

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
       end if
       
       do j=1, ngpslev-1
          if ((h <= prf%gst(j)%Var) .and. (h > prf%gst(j+1)%Var)) then
             jloc = j
             exit
          end if
       end do
       
       if (h <= prf%gst(ngpslev)%Var) then
          jloc = ngpslev-1
       end if
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
    end do
  end subroutine gpsbvfopv


!modgps08ztdop

  pure subroutine gps_ztdopv(hv, prf, lbevis, dzmin, ZTDopv, rPobs, mode)
    !
    !:Purpose: GB-GPS ZTD operator
    !
    !:Arguments:
    !   :dzmin:   Minimum DZ = Zobs-Zmod (m) for which DZ adjustment to ZTD will
    !             be made when Zobs < Zmod.
    !   :mode:    1 = normal mode: use stored ZTD profiles
    !
    !             2 = Vedel & Huang ZTD formulation: ZTD = ZHD(Pobs) + ZWD
    !                 Pobs computed from P0 using CMC hydrostatic extrapolation.
    !
    implicit none

    ! Arguments:
    real(dp)              , intent(in) :: hv    ! height of ZTD observation Zobs (m)
    type(gps_profilezd)   , intent(in) :: prf   ! local model profile (type gps_profilezd)
    logical               , intent(in) :: lbevis! true/false --> use Bevis instead of Rueger k values
    real(dp)              , intent(in) :: dzmin
    type(gps_diff)        , intent(out):: ZTDopv! ZTD (m) at height of observation (with derivatives)
    real(dp)              , intent(out):: rPobs ! Pressure (Pa) at height of observation
    integer               , intent(in) :: mode

    ! Locals:
    integer(i4)                        :: ngpslev
    integer(i4)                        :: j, jloc
    real(dp)                           :: h, x, lat, sLat, dh
    real(dp)                           :: k1, k2, k3, k2p
    real(dp)                           :: zcon, zcon1, zconh, zfph, zconw
    
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
    end if
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
      end if
    end do

    if (h <= prf%gst(ngpslev)%Var) then  ! obs is at or below model lowest level
      jloc = ngpslev
    end if
    
    if ( mode == 2 ) then
    
!  Compute ZTD the Vedel and Huang (2004) way: (as in old s/r gpsztdop.ftn)

       zcon  = 1.0e-06_dp*p_Rd
       zcon1 = zcon*k1
       zconw = zcon/eps
       zconh = zcon1/Rgm
       zfph = (1.0_dp - 2.66e-03_dp*cos(2.0*lat) - 2.8e-07_dp*h)

!  Pressure at obs height (CMC hydrostatic extrapolation from Psfc)
       x      = ec_rg/(p_Rd*gamma)
       tvsfc  = prf%tst(ngpslev)*(1._dp+delta*prf%qst(ngpslev))
       Pobs   = prf%pst(ngpslev)*(((tvsfc-gamma*dh)/tvsfc)**x)
!  Dry delay ZHD (m) at obs height
       zhd    = (zconh/zfph) * Pobs

! Integrate column q/T on pressure levels to get model ZWD
       do j = 1, ngpslev-1
         tbar = (prf%tst(j) + prf%tst(j+1))*0.5_dp
         qbar = (prf%qst(j) + prf%qst(j+1))*0.5_dp
         qtterm = ((qbar + kappa*qbar**2 )/gps_gravityalt(sLat,prf%gst(j)%Var))*(k2p + k3/tbar)
         if ( j == 1 ) then
           zsum = qtterm*(prf%pst(j+1)-prf%pst(j))
         else
           zsum = zsum + qtterm*(prf%pst(j+1)-prf%pst(j))
         end if
       end do
       
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
       end if
       
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
          x      = ec_rg/(p_Rd*gamma)
          tvsfc  = prf%tst(jloc)*(1._dp+delta*prf%qst(jloc))
          Pobs   = prf%pst(jloc)*(((tvsfc-gamma*dh)/tvsfc)**x)
          if ( abs(dh) <= dzmax ) then
            dztddpm = prf%rst(jloc)   ! lowest level value of dZTD/dp
          else
            tobs   = prf%tst(jloc)-gamma*dh
            qobs   = prf%qst(jloc)
            tvobs  = tvsfc-gamma*dh
            naobs  = (k1/tvobs) + (k2p*(qobs/(eps*tobs))) + (k3*(qobs/(eps*tobs**2)))
            dztddp = 1.e-6_dp * naobs * (p_Rd*tvobs)/gps_gravityalt(sLat, h)
            dztddpm = (dztddp + prf%rst(jloc))/2._dp  ! mean value of dZTD/dp over dh layer
          end if
          ZTDopv = prf%ztd(jloc) + dztddpm*(Pobs-prf%pst(jloc))
        end if

      end if
    
    end if
    
    rPobs = Pobs%Var

  end subroutine gps_ztdopv

  subroutine gps_pw(prf, PW)
  !
  !:Purpose: To compute lowest level PW (kg/m2) using layer mean Q and layer
  !          delta_p (Pa)
  !
    implicit none

    ! Arguments:
    type(gps_profilezd)     , intent(in)  :: prf
    real(dp)                , intent(out) :: PW

    ! Locals:
    integer(i4)                       :: i, ngpslev
    real(dp)                          :: qbar, gt, gb, g, lat, sLat
    real(dp)                          :: pt, pb

    ngpslev = prf%ngpslev
    lat     = prf%rLat
    sLat    = sin(lat)

    PW = 0.0_dp

    do i = 1, ngpslev-1
      qbar = 0.5_dp * (prf%qst(i+1)%Var + prf%qst(i)%Var)
      gt  = gps_gravityalt(sLat, prf%gst(i)%Var)
      gb  = gps_gravityalt(sLat, prf%gst(i+1)%Var)
      pt  = prf%pst(i)%Var
      pb  = prf%pst(i+1)%Var
      g   = 0.5_dp * (gt + gb)
      PW = PW + (qbar/g)*(pb-pt)
    end do

  end subroutine gps_pw


!modgps09bend

  subroutine gpsbend(prf)
    implicit none
 
    ! Arguments:
    type(gps_profile)     :: prf

    ! Locals:
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
    integer                            :: i,j,ngpslev
    logical                            :: lok

    if (.not. prf%bbst) then
       ngpslev=prf%ngpslev

       ! Radial distances and impact parameters:
       do i=1,ngpslev
          prf%dst(i)= (prf%Rad+prf%geoid+prf%gst(i))
          prf%ast(i)= prf%dst(i) * (1._dp+1.e-6_dp*prf%rst(i))
       end do
       ! Extended lower levels:
       do i=ngpslev+1,ngpslev+ngpsxlow
          prf%ast(i)= prf%ast(i-1)-50._dp
       end do

       ! Standard levels:
       do i=1,ngpslev
          r  (i)=prf%dst(ngpslev-i+1)
          ref(i)=prf%rst(ngpslev-i+1)
          !ref(i)=300._dp*exp((-1._dp/7000._dp)*(r(i)%Var-prf%Rad))
       end do
       ! Extended upper levels:
       do i=ngpslev+1,ngpssize
          r  (i)=r  (i-1)+1000._dp
          ref(i)=ref(i-1)*exp(-1000._dp/7000_dp)
       end do

       ! log n and x:
       do i=1,ngpssize
          nu(i)=1.e-6_dp*ref(i)
          lnu(i)=log(nu(i))
          n (i)=1._dp+nu(i)
          x (i)=n(i)*r(i)
          rsq(i)=r(i)**2
          nsq(i)=n(i)**2
          xsq(i)=x(i)**2
       end do
       do i=0,-ngpsxlow+1,-1
          x  (i)=x(i+1)-50._dp
          xsq(i)=x(i)**2
       end do

       ! Radial derivatives of log refractivity.
       ! Refractivity will be assumed exponential within each shell.
       ! We store the derivative of log(nu).
       ! dn/dr = nu * dlgnudr
       do i=1,ngpssize-1
          dlgnudr(i)=(lnu(i+1)-lnu(i))/(r(i+1)-r(i))
       end do

       ! Evaluation of complete bending for ray tangent at r(i):
       do i=1,ngpslev
          ! Check that ray is not trapped
          lok=.true.
          do j = i+1,ngpssize
             lok= lok .and. (x(j)%Var .gt. x(i)%Var)
          end do
          if (lok) then
             s(i)=0._dp
             t(i)=1._dp
             do j=i+1,ngpssize
                s(j)=sqrt(nsq(i)*rsq(j)-xsq(i))
                t(j)=s(j)/sqrt(xsq(j)-xsq(i))
             end do

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
             end do
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
             end do
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
             end do
             boole=(-2)*r(i)*sum

             prf%bst(ngpslev-i+1)=boole
          else
             prf%bst(ngpslev-i+1)=-10._dp
          end if
       end do

       ! Extended low levels:
       do i=0,-ngpsxlow+1,-1
          lok=.true.
          do j = 1,ngpssize
             lok= lok .and. (x(j)%Var .gt. x(i)%Var)
          end do
          if (lok) then
             do j=1,ngpssize
                s(j)=sqrt(nsq(1)*rsq(j)-xsq(i))
                t(j)=s(j)/sqrt(xsq(j)-xsq(i))
             end do

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
             end do
             simp=(-2)*(x(i)/n(1))*sum
             alpha_B=acos(x(i)/x(1))
             prf%bst(ngpslev-i+1)=simp-2*alpha_B
          else
             prf%bst(ngpslev-i+1)=-10._dp
          end if
       end do

       prf%bbst=.true.
    end if
  end subroutine gpsbend

  subroutine gpsbend1(prf)
    implicit none

    ! Arguments:
    type(gps_profile)     :: prf

    ! Locals:
    type(gps_diff)                     :: r  (ngpssize)
    type(gps_diff)                     :: ref(ngpssize)
    type(gps_diff)                     :: nu (ngpssize)
    type(gps_diff)                     :: lnu(ngpssize)
    type(gps_diff)                     :: n  (ngpssize)
    type(gps_diff)                     :: dlgnudr(ngpssize-1)
    type(gps_diff)                     :: x  (-ngpsxlow+1:ngpssize)

    type(gps_diff)                     :: angle0,angle,angleB,bend,nu0,th,sum,nexp
    real(dp)                           :: dxn
    integer                            :: ngpslev,i,j,jmin
    logical                            :: lok, lok2

    if (.not. prf%bbst) then
       ngpslev=prf%ngpslev

       ! Radial distances and impact parameters:
       do i=1,ngpslev
          prf%dst(i)= (prf%Rad+prf%geoid+prf%gst(i))
          prf%ast(i)= prf%dst(i) * (1._dp+1.e-6_dp*prf%rst(i))
       end do
       ! Extended lower levels:
       do i=ngpslev+1,ngpslev+ngpsxlow
          prf%ast(i)= prf%ast(i-1)-50._dp
       end do

       ! Standard levels:
       do i=1,ngpslev
          r  (i)=prf%dst(ngpslev-i+1)
          ref(i)=prf%rst(ngpslev-i+1)
       end do
       ! Extended upper levels:
       do i=ngpslev+1,ngpssize
          r  (i)=r  (i-1)+1000._dp
          ref(i)=ref(i-1)*exp(-1000._dp/7000_dp)
       end do

       ! log n and x:
       do i=1,ngpssize
          nu(i) = 1.e-6_dp*ref(i)
          lnu(i)= log(nu(i))
          n (i) = 1._dp+nu(i)
          x (i) = n(i)*r(i)
       end do
       dxn=20._dp
       do i=0,-ngpsxlow+1,-1
          x (i) = x(i+1)-dxn
       end do

       ! Radial derivatives of log refractivity.
       ! Refractivity will be assumed exponential within each shell.
       ! We store the derivative of log(nu).
       ! dn/dr = nu * dlgnudr
       do i=1,ngpssize-1
          dlgnudr(i)=(lnu(i+1)-lnu(i))/(r(i+1)-r(i))
       end do

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
          end do
          if (lok) then
             ! Integration:
             sum=0._dp
             if (i.ge.1) then
                ! Direct rays
                angleB=0._dp
             else
                ! Reflected
                angleB=sqrt(2*(-i+1)*dxn/x(1))
             end if
             angle0=angleB
             do j=jmin, ngpssize-1
                th=r(j+1)-r(j)
                nu0=nu(j)
                nexp=dlgnudr(j)
                call gpsbendlayer(r(j), th, nu0, nexp, angle0, angle, bend, lok2)
                sum=sum+bend
                angle0=angle
             end do
          end if
          if (lok2) then
             prf%bst(ngpslev-i+1)=(-2)*(sum+angleB)
          else
             prf%bst(ngpslev-i+1)=-10._dp
          end if
       end do
    end if
  end subroutine gpsbend1

  subroutine gpsbendlayer(ra, th, nu0, nexp, angle0, angle, bend, lok)
    !
    !:Arguments:
    !     :nu0, nexp:  Refraction index coefs: n=1+nu0*exp(nexp*(r-ra));
    !                  nexp in in 1/m
    implicit none

    ! Arguments:
    type(gps_diff), intent(in) :: ra    ! Radius of inner shell (m)
    type(gps_diff), intent(in) :: th    ! Shell thickness (m)
    type(gps_diff), intent(in) :: nu0   
    type(gps_diff), intent(in) :: nexp
    type(gps_diff), intent(in) :: angle0 ! Ray angle above horizon at ra
    type(gps_diff), intent(out):: angle  ! Ray angle above horizon at rb
    type(gps_diff), intent(out):: bend   ! Accumulated bending over the layer
    logical       , intent(out):: lok

    ! Locals:
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
    end do
    angle=anglei
  end subroutine gpsbendlayer

  subroutine gpsbendunit(ra, th, nu0, nexp, angle0, angle, bend, lok)
    !
    !:Arguments:
    !    :nu0, nexp:  Refraction index coefs: n=1+nu0*exp(nexp*(r-ra));
    !                 nexp in 1/m
    implicit none

    ! Arguments:
    type(gps_diff), intent(in) :: ra ! Radius of inner shell (m)
    type(gps_diff), intent(in) :: th ! Shell thickness (m)
    type(gps_diff), intent(in) :: nu0, nexp
    type(gps_diff), intent(in) :: angle0 ! Ray angle above horizon at ra
    type(gps_diff), intent(out):: angle  ! Ray angle above horizon at rb
    type(gps_diff), intent(inout):: bend ! Accumulated bending over the layer
    logical       , intent(out):: lok

    ! Locals:
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
    end if
  end subroutine gpsbendunit

  subroutine gpsbndopv(impv, azmv, nval, prf, bstv)
    implicit none

    ! Arguments:
    real(dp)              , intent(in) :: impv(:), azmv(:)
    integer(i4)           , intent(in) :: nval
    type(gps_profile)                  :: prf
    type(gps_diff)        , intent(out):: bstv(:)

    ! Locals:
    integer                            :: iSize, i, j, ngpslev, jlocm, jlocp
    real(dp)                           :: imp1,azm1,rad, rad0
    real(dp)                           :: imp(ngpssize+ngpsxlow)
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
       end do
       !
       ! Search where it is located
       !
       jlocm = -1000
       jlocp = -1000
       if (imp1 > imp(1)) then
          jlocm = 1
          jlocp = 2
       end if

       do j=1, ngpslev+ngpsxlow-1
          if ((imp1 <= imp(j)) .and. (abs(prf%bst(j)%Var) < 1._dp)) then
             jlocm = j
          end if
       end do

       do j=jlocm+1, ngpslev+ngpsxlow
          if ((imp1 >  imp(j)) .and. (abs(prf%bst(j)%Var) < 1._dp)) then
             jlocp = j
             exit
          end if
       end do
      
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
       end if
    end do
  end subroutine gpsbndopv

  subroutine gps_bndopv1(impv, azmv, nval, prf, bstv)
    ! :Purpose: Computation of the observation operator for Bending 
    !
    ! :Note: The Operator is based  from Assimilation experiments withCHAMP GPS radio occultation measurements
    !         (S. B. HEALY and J.-N. THEPAUT, 2005)
    implicit none

    ! Arguments:
    real(dp)              , intent(in) :: impv(:), azmv(:)
    integer(i4)           , intent(in) :: nval
    type(gps_profile)                  :: prf
    type(gps_diff)        , intent(out):: bstv(:)

    ! Locals:
    integer                            :: levIndexObs, ngpslev,last_levIndexObs,numLevels,levelHigh,levelLow,levIndexAnl 
    type(gps_diff)                     :: h(ngpssize), nu(ngpssize), lnu(ngpssize), n(ngpssize), z(ngpssize)
    type(gps_diff)                     :: N0a, N1a, ka, NAa,Aa, Ba, Ba2, Ba3, delta_alpha, delta_alpha_top, z_0, h_0
    real(dp)                           :: a2, a, gz, cazm, sazm, last_a

    ! model levels 
    ngpslev=prf%ngpslev
    do levIndexAnl = 1,ngpslev 
      h(levIndexAnl)  = prf%geoid+prf%gst(levIndexAnl)
      nu(levIndexAnl) = prf%rst(levIndexAnl)
      lnu(levIndexAnl)=log(nu(levIndexAnl)) 
      n(levIndexAnl)  = 1._dp+nu(levIndexAnl)*1e-6_dp
      z(levIndexAnl)  = n(levIndexAnl)*(prf%Rad+prf%geoid+prf%gst(levIndexAnl))
    end do
    ! number of observed levels in the profile
    numLevels  = size(impv)
    if (nval < numLevels) numLevels=nval
    
    do levIndexObs =  numLevels,1,-1
      a2 = impv(levIndexObs)*impv(levIndexObs)
      a  = impv(levIndexObs)
      cazm = cos(azmv(levIndexObs))
      sazm = sin(azmv(levIndexObs))
      !find model levels that bracket the observation
      !   note to self:   like in GEM, level=1 is the highest level
      do levIndexAnl = 1, ngpslev-1
        levelHigh = levIndexAnl - 1
        levelLow  = levIndexAnl
        if (z(levIndexAnl)%VaR< a) exit 
          levelLow  = 0
      end do
    
      if  (levelLow/=0) then
        h_0  = h(levelLow)+(((a-z(levelLow))/(z(levelHigh)-z(levelLow)))*(h(levelHigh)-h(levelLow)))  
        N0a  = nu(levelLow)
        gz   = (lnu(levelLow+1)%Var-lnu(levelLow)%Var)/(h(levelLow+1)%Var-h(levelLow)%Var)
        NAa  = nu(levelLow)*exp(gz*(h_0-h(levelLow)))
        delta_alpha = 0.d0
        delta_alpha_top = 0.d0
        z_0 = a
        do while ((levelHigh)>=1) 
          N1a = nu(levelHigh)
          ka  = log(N0a/N1a)/(z(levelHigh) - z(levelLow))
          ! Test of Reflected
          do while (ka%Var<0.or.ka%Var>0.1) 
            levelHigh = levelHigh - 1
            N1a = nu(levelHigh)
            ka = log(N0a/N1a)/(z(levelHigh) - z(levelLow))
          end do
          Aa = 1e-6_dp* sqrt((MPC_PI_R8/2.d0*a)/ka )*NAa* exp(ka*(z_0-a))  
          if (z_0%Var==a) then
            Ba  = erf(sqrt(ka*(z(levelHigh)-a)))
          else  
            Ba2 = erf(sqrt(ka*(z(levelHigh)-a)))
            Ba3 = erf(sqrt((ka*(z_0-a))))
            Ba  = Ba2 - Ba3
          end if   
          delta_alpha = delta_alpha+2*ka*Ba*Aa 
          N0a = N1a
          NAa = N1a 
          z_0 = z(levelHigh)

          levelLow  = levelHigh 
          levelHigh = levelLow-1
        end do
        Ba = erf(1-erf(sqrt((ka*(z_0-a)))))
        delta_alpha_top = 2*Aa*ka*Ba 
        last_a = a
        last_levIndexObs = levIndexObs
        bstv(levIndexObs)= delta_alpha +delta_alpha_top
      else  ! (levelLow==0) 
        ! Use loglinear extrapolation for most data (notably direct rays)
        if (a>(1._dp+prf%rst(ngpslev)%Var*1e-6_dp)*prf%Rad) then
          bstv(levIndexObs)=bstv(last_levIndexObs)*exp((-1._dp/6500._dp)*(a-last_a))
        else
        ! Use linear extrapolation (most reflected rays) from Information content in reflected
        ! signals during GPS Radio Occultation observations (Josep Aparicio et al.,2017) 
          bstv(levIndexObs)=bstv(last_levIndexObs)*exp((-1._dp/6500._dp)*(a-last_a))-2*acos(a/((1._dp+prf%rst(ngpslev)*1e-6_dp)*prf%Rad))
        end if

      end if
    end do
  end subroutine gps_bndopv1


!modgpsro_mod

  subroutine gps_setupro
    implicit none

    ! Locals:
    integer :: nulnam,ierr,fnom,fclos,SatID
    
    ! Namelist variables for GPS-RO
    INTEGER :: LEVELGPSRO       ! Data level to use (1 for bending angle, 2 for refractivity)
    INTEGER :: GPSRO_MAXPRFSIZE ! Maximal number of data that is expected from a profile (default 300)
    REAL(8) :: SURFMIN          ! Minimum allowed distance to the model surface (default 0 m)
    REAL(8) :: HSFMIN           ! Minimum allowed MSL height of an obs          (default 0 m)
    REAL(8) :: HTPMAX           ! Maximum allowed MSL height of an obs          (default 70000 m)
    REAL(8) :: BGCKBAND         ! Maximum allowed deviation abs(O-P)/P          (default 0.05)
    REAL(8) :: HTPMAXER         ! Maximum MSL height to evaluate the obs error  (default to HTPMAX)
    REAL(4) :: WGPS(0:1023,4)   ! WGPS values for each satellite sensor
    character(len=20) :: gpsroError ! key for using dynamic/static refractivity error estimation (default 'DYNAMIC')
    LOGICAL :: gpsroBNorm       ! Choose to normalize based on B=H(x) (default=.True.), or approximate exponential reference
    LOGICAL :: gpsroEotvos      ! Add an operator-only Eotvos correction to local gravity (shift of altitudes, default False)

    NAMELIST /NAMGPSRO/ LEVELGPSRO, GPSRO_MAXPRFSIZE, SURFMIN, HSFMIN, HTPMAX, HTPMAXER, &
                        BGCKBAND, WGPS, gpsroError, gpsroBNorm, gpsroEotvos


!
!   Define default values:
!
    LEVELGPSRO = gps_Level_RO_Ref
    GPSRO_MAXPRFSIZE = 300
    SURFMIN    = 0.d0
    HSFMIN     = 0.d0
    HTPMAX     = 70000.d0
    HTPMAXER   = -1.d0
    BGCKBAND   = 0.05d0
    gpsroError = 'DYNAMIC'
    gpsroBNorm = .True.
    gpsroEotvos= .False.
!
!   Force a pre-NML default for the effective data weight of all
!   GPSRO satellites. This array has rows 0-1023 (following BUFR element
!   SATID), and 4 cols. The 4 parameters for each SATID are used to
!   represent data correlation, a combined property of the satellite
!   hardware and provider postprocessing.
!   The default assumes no correlation. 
!
    WGPS = 0.
    WGPS(:,1) = 1.
!
!   Override with NML values:
!     
    nulnam=0
    ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=NAMGPSRO,iostat=ierr)
    if(ierr.ne.0) call utl_abort('gps_setupro: Error reading namelist')
    !if(mmpi_myid.eq.0) write(*,nml=NAMGPSRO)
    ierr=fclos(nulnam)
    if (HTPMAXER < 0.0D0) HTPMAXER = HTPMAX
    gps_Level_RO      = LEVELGPSRO
    gps_RO_MAXPRFSIZE = GPSRO_MAXPRFSIZE
    gps_SurfMin       = SURFMIN
    gps_HsfMin        = HSFMIN
    gps_HtpMax        = HTPMAX
    gps_HtpMaxEr      = HTPMAXER
    gps_BgckBand      = BGCKBAND
    gps_roError       = gpsroError
    gps_roBNorm       = gpsroBNorm
    gps_WGPS = WGPS
    gps_gpsroEotvos = gpsroEotvos

    if(mmpi_myid.eq.0) then
      write(*,*)'NAMGPSRO',gps_Level_RO, gps_RO_MAXPRFSIZE, gps_SurfMin, gps_HsfMin, &
           gps_HtpMax, gps_HtpMaxEr, gps_BgckBand, trim(gps_roError), gps_roBNorm, gpsroEotvos
      do SatID = 0, 1023
        if (WGPS(SatID,2) /= 0.) then
          write(*,*)'WGPS', SatID, gps_WGPS(SatID, 1:4)
        end if
      end do
    end if
  end subroutine gps_setupro

  integer function gps_iprofile_from_index(index)
    implicit none

    ! Arguments:
    integer, intent(in) :: index

    ! Locals:
    integer i

    gps_iprofile_from_index=-1
    do i=1,gps_numROProfiles
       if (index.eq.gps_vRO_IndexPrf(i, 1)) then
          gps_iprofile_from_index=i
          return
       end if
    end do
    return
  end function gps_iprofile_from_index


!modgpsztd_mod

  subroutine gps_setupgb
    !
    !:Purpose: Initialisation of ground-based GPS - to read and to initialize
    !          GB-GPS namelist parameters and print information on options
    !          selected.
    !
    implicit none

    ! Locals:
    integer :: nulnam,ierr,fnom,fclos

    ! Namelist variables for Ground-based GPS (ZTD)
    REAL(8) :: DZMIN            ! Minimum DZ = Zobs-Zmod (m) for which DZ adjustment to ZTD will be made
    REAL(8) :: DZMAX = 1000.0D0 ! Maximum DZ (m) over which ZTD rejected due to topography (when LTOPOFILT = .TRUE.)
    REAL(8) :: YZTDERR          ! If < 0 use errors in input files; if > 0 use value as constant error (m); if 0 compute error as f(ZWD)
    REAL(8) :: YSFERRWGT        ! Scale factor for GPS surface met errors (account for time series obs with error correlations)
    REAL(8) :: YZDERRWGT        ! Scale factor for GPS ZTD errors (account for time series obs with error correlations)
    LOGICAL :: LASSMET          ! Choose to assimilate GPS Met surface P, T, T-Td
    LOGICAL :: LLBLMET          ! Indicate that surface met data blacklisted for GPS sites close to surface weather stations.
    LOGICAL :: LBEVIS           ! If .true. use Bevis(1994); if .false. use Rueger(2002) refractivity (k1,k2,k3) constants
    LOGICAL :: L1OBS            ! Choose to select a single ZTD observation
    LOGICAL :: LTESTOP          ! Choose to test ZTD observation operator (Omp and Bgck modes only)
    INTEGER :: IREFOPT          ! 1 = conventional expression for N using k1,k2,k3; 2 = Aparicio & Laroche N (incl. compressibility)
    INTEGER :: IZTDOP           ! 1 = use stored ZTD profiles to get ZTDmod; 2 = Vedel & Huang ZTD formulation: ZTDmod = ZHD(Pobs) + ZWD

    NAMELIST /NAMGPSGB/ DZMIN, DZMAX, YZTDERR, LASSMET, YSFERRWGT,  &
         LLBLMET, YZDERRWGT, LBEVIS, L1OBS, LTESTOP, IREFOPT, IZTDOP

!*  .  1.1 Default values
!!  .      --------------

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
    gps_gb_DZMIN     = DZMIN
    gps_gb_DZMAX     = DZMAX
    gps_gb_YZTDERR   = YZTDERR
    gps_gb_LASSMET   = LASSMET
    gps_gb_YSFERRWGT = YSFERRWGT
    gps_gb_LLBLMET   = LLBLMET
    gps_gb_YZDERRWGT = YZDERRWGT
    gps_gb_LBEVIS    = LBEVIS
    gps_gb_IREFOPT   = IREFOPT
    gps_gb_L1OBS     = L1OBS
    gps_gb_LTESTOP   = LTESTOP
    gps_gb_IZTDOP    = IZTDOP
    if(mmpi_myid.eq.0) write(*,nml=NAMGPSGB)
    ierr=fclos(nulnam)

    IF (L1OBS.and.mmpi_myid.eq.0) THEN
      write(*,*)' '
      write(*,*)' ******************************************'
      write(*,*)' *        GB-GPS OBSERVATIONS             *'
      write(*,*)' *                                        *'
      write(*,*)' *        ONE OBSERVATION MODE            *'
      write(*,*)' *                                        *'
      write(*,*)' ******************************************'
      write(*,*)' '
    END IF

!   Options to fix/adjust model ZTD to observation height and
!   assimilate GPS met data

    if(mmpi_myid.eq.0) then
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
        END IF
      ELSE
        write(*,*)' '
        write(*,*)' *****************************************'
        write(*,*)' *          GB-GPS OBSERVATIONS          *'
        write(*,*)' *   GPS MET DATA ARE NOT ASSIMILATED    *'
        write(*,*)' *                                       *'
        write(*,*)' *****************************************'
        write(*,*) 'YZDERRWGT = ', YZDERRWGT
        write(*,*) ' '
      END IF

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
      END IF

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
        END IF
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
        END IF       
      ELSE
        write(*,*)' '
        write(*,*)' *****************************************'
        write(*,*)' *          GB-GPS OBSERVATIONS          *'
        write(*,*)' *                                       *'
        write(*,*)' *  APARICIO & LAROCHE REFRACTIVITY N    *'
        write(*,*)' *         USED TO COMPUTE ZTD           *'
        write(*,*)' *****************************************'
        write(*,*)' '       
      END IF

    end if

  end subroutine gps_setupgb

  integer function gps_iztd_from_index(index)
    implicit none
 
    ! Arguments:
    integer, intent(in) :: index

    ! Locals:
    integer i

    gps_iztd_from_index = -1
    do i = 1, size(gps_ZTD_Index)
       if (index .eq. gps_ZTD_Index(i)) then
          gps_iztd_from_index = i
          return
       end if
    end do
    return
  end function gps_iztd_from_index

end module gps_mod

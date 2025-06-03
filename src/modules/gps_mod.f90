
module gps_mod
  ! MODULE gps_mod (prefix='gps' category='5. Observation operators')
  !
  ! :Purpose: Code related to GPS-RO and ground-based GPS observation operators.
  !
  use midasMpi_mod
  use utilities_mod
  use mathPhysConstants_mod
  use earthConstants_mod
  use physicsFunctions_mod
  use timeCoord_mod

  implicit none
  save
  private

  ! Public derived types
  public :: gps_diff          ! Multidual real numbers (they keep track of the derivatives)
  public :: gps_profile       ! Encapsulate vertical profiles for GPS RO
  public :: gps_profilezd     ! Encapsulate vertical profiles for GPS GB

  ! Public variables
  ! RO
  integer      , public, protected, allocatable :: gps_vRO_IndexPrf(:, :) ! header index for each profile
  real(8)      , public, protected, allocatable :: gps_vRO_dR(:, :)       ! vertical asphericity shift for each datum
  integer      , public, protected :: gps_numROProfiles                   ! number of RO Profiles handled
  integer      , public, protected :: gps_level_RO
  integer      , public, protected :: gps_RO_maxPrfSize
  real(8)      , public, protected :: gps_surfMin
  real(8)      , public, protected :: gps_hSfMin
  real(8)      , public, protected :: gps_hTpMax
  real(8)      , public, protected :: gps_bgckBand
  real(8)      , public, protected :: gps_hTpMaxEr
  real(8)      , public, protected :: gps_roNsigma
  real(4)      , public, protected :: gps_wGPS(0:1023, 4)
  character(20), public, protected :: gps_roError
  logical      , public, protected :: gps_roBNorm
  logical      , public, protected :: gps_roEotvos
  logical      , public, protected :: gps_roCurvAnisot
  logical      , public, protected :: gps_roNCurv
  integer      , public, protected :: gps_roNFlavour
  real(8)      , public, protected :: gps_roYear
  ! ZTD
  integer      , public, protected, allocatable :: gps_ZTD_Index (:)      ! header index in CMA for each ZTD observation
  integer      , public, protected :: gps_gb_numZTD                       ! number of ZTD data to be assimilated
  real(8)      , public, protected :: gps_gb_dzMin
  real(8)      , public, protected :: gps_gb_dzMax
  real(8)      , public, protected :: gps_gb_yZTDErr
  real(8)      , public, protected :: gps_gb_ySFErrWgt
  real(8)      , public, protected :: gps_gb_yZDErrWgt
  integer      , public, protected :: gps_gb_iRefOpt
  integer      , public, protected :: gps_gb_iZTDOp
  logical      , public, protected :: gps_gb_lAssMet
  logical      , public, protected :: gps_gb_lLBLMet
  logical      , public, protected :: gps_gb_lBevis
  logical      , public, protected :: gps_gb_l1Obs
  logical      , public, protected :: gps_gb_lTestOp

  ! Public variables (parameters)
  public :: gps_ncvmx, gps_gb_maxdata
  public :: gps_p_md, gps_p_mw, gps_p_wa, gps_p_wb
  public :: gps_level_RO_Bnd, gps_level_RO_Ref, gps_level_RO_BndAndRef

  ! Public procedures
  public :: gps_setupro, gps_iprofile_from_index
  public :: gps_setNumROProfiles, gps_setROIndexPrf
  public :: gps_setupgb, gps_iztd_from_index
  public :: gps_setNumZTD, gps_setZTDIndex
  public :: gps_struct1sw, gps_struct1sw_v2, gps_bndopv, gps_refopv, gps_structztd_v2, gps_ztdopv, gps_pw, gps_pwopv

  integer    , parameter :: gps_level_RO_Bnd       = 1
  integer    , parameter :: gps_level_RO_Ref       = 2
  integer    , parameter :: gps_level_RO_BndAndRef = 3

  ! 32-bit integers
  integer    , parameter :: i4 = selected_int_kind(9)

  ! Short floats
  integer    , parameter :: sp = selected_real_kind(6)

  ! Long floats
  integer    , parameter :: dp = selected_real_kind(12)

  ! Maximum number of gps levels:
  integer(i4), parameter :: ngpssize  = 100

  ! Associated maximum number of control variables:
  integer(i4), parameter :: gps_ncvmx = 4*ngpssize

  ! Avogadro constant:
  real(dp)   , parameter :: p_Avog    = 6.02214129e23_dp                    ! From CODATA

  ! Boltzmann constant:
  real(dp)   , parameter :: p_Boltz   = 1.3806488e-23_dp                    ! From CODATA

  ! Air properties (public):
  real(dp)   , parameter :: gps_p_md  = 28.965516_dp                        ! From Aparicio(2011)
  real(dp)   , parameter :: gps_p_mw  = 18.015254_dp                        ! From Aparicio(2011)
  real(dp)   , parameter :: gps_p_wa  = gps_p_md/gps_p_mw
  real(dp)   , parameter :: gps_p_wb  = (gps_p_md-gps_p_mw)/gps_p_mw

  ! Gas constants:
  real(dp)   , parameter :: p_R       = p_Avog*p_Boltz                      ! per mol
  real(dp)   , parameter :: p_Rd      = p_Avog*p_Boltz/(1.e-3_dp*gps_p_md)  ! per air mass

  ! Multiduals of real numbers, with gps_ncvmx dual units
  type gps_diff
    real(dp) :: Var
    real(dp) :: DVar(gps_ncvmx)
  end type gps_diff

  ! Overloaded multidual operations =, +, -, *, /, **, sqrt, exp, log, cos, tan, acos, atan, erf

  ! Assignment onto a diff multidual (of a double, and of a diff)
  interface assignment(=)
    module procedure gpsdiffasfd, gpsdiffasff
  end interface assignment(=)

  ! Addition involving diff multiduals (diff+double, double+diff, diff+int, int+diff, diff+diff)
  interface operator(+)
    module procedure gpsdiffsmfd, gpsdiffsmdf, gpsdiffsmfi, gpsdiffsmif, gpsdiffsmff
  end interface operator(+)

  ! Subtraction involving diff multiduals (diff-double, double-diff, diff-int, int-diff, diff-diff)
  interface operator(-)
    module procedure gpsdiffsbfd, gpsdiffsbdf, gpsdiffsbfi, gpsdiffsbif, gpsdiffsbff
  end interface operator(-)

  ! Multiplication involving diff multiduals (diff*double, double*diff, diff*int, int*diff, diff*diff)
  interface operator(*)
    module procedure gpsdiffmlfd, gpsdiffmldf, gpsdiffmlfi, gpsdiffmlif, gpsdiffmlff
  end interface operator(*)

  ! Division involving diff multiduals (diff/double, double/diff, diff/int, int/diff, diff/diff)
  interface operator(/)
    module procedure gpsdiffdvfd, gpsdiffdvdf, gpsdiffdvfi, gpsdiffdvif, gpsdiffdvff
  end interface operator(/)

  ! Power involving diff multiduals (diff^double, double^diff, diff^int, int^diff, diff^diff)
  interface operator(**)
    module procedure gpsdiffpwfd, gpsdiffpwdf, gpsdiffpwfi, gpsdiffpwif, gpsdiffpwff
  end interface operator(**)

  ! Square root of a diff multidual
  interface sqrt
    module procedure gpsdiffsqr
  end interface sqrt

  ! Exp of a diff multidual
  interface exp
    module procedure gpsdiffexp
  end interface exp

  ! Logarithm of a diff multidual
  interface log
    module procedure gpsdifflog
  end interface log

  ! Cosine of a diff multidual
  interface cos
    module procedure gpsdiffcos
  end interface cos

  ! Tangent of a diff multidual
  interface tan
    module procedure gpsdifftan
  end interface tan

  ! Arc Cosine of a diff multidual
  interface acos
    module procedure gpsdiffacos
  end interface acos

  ! Arc Tangent of a diff multidual
  interface atan
    module procedure gpsdiffatan
  end interface atan

  ! Error function of a diff multidual
  interface erf
    module procedure gpsdifferf
  end interface erf

  type gps_profile
    integer(i4)                         :: ngpslev
    real(dp)                            :: rLat
    real(dp)                            :: rLon
    real(dp)                            :: rAzm
    real(dp)                            :: rMT
    real(dp)                            :: Rad
    real(dp)                            :: geoid
    real(dp)                            :: RadN
    real(dp)                            :: RadM

    type(gps_diff)                      :: P0

    type(gps_diff)                      :: pst (ngpssize)
    type(gps_diff)                      :: tst (ngpssize)
    type(gps_diff)                      :: qst (ngpssize)
    type(gps_diff)                      :: rst (ngpssize)
    type(gps_diff)                      :: gst (ngpssize)
    type(gps_diff)                      :: vst (ngpssize)
    type(gps_diff)                      :: lrmd(ngpssize)
  end type gps_profile

  type gps_profilezd
    integer(i4)                         :: ngpslev
    real(dp)                            :: rLat
    real(dp)                            :: rLon
    real(dp)                            :: rMT

    type(gps_diff)                      :: P0

    type(gps_diff)                      :: pst (ngpssize)
    type(gps_diff)                      :: tst (ngpssize)
    type(gps_diff)                      :: qst (ngpssize)
    type(gps_diff)                      :: rst (ngpssize)
    type(gps_diff)                      :: gst (ngpssize)
    type(gps_diff)                      :: ztd (ngpssize)
    logical                             :: bpst
  end type gps_profilezd

  integer    , parameter :: max_gps_sites  = 1200
  integer    , parameter :: gps_gb_maxdata = max_gps_sites*24               ! (max_gps_sites) * (max_num_obs in 6h)

  contains

  !--------------------------------------------------------------------------
  ! gpsdiffasfd
  !--------------------------------------------------------------------------
  pure subroutine gpsdiffasfd(gd1, d2)
    ! :Purpose: Overloaded =
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(out)    :: gd1
    real(dp)           , intent(in)     :: d2

    gd1%Var     = d2
    gd1%DVar(:) = 0._dp
  end subroutine gpsdiffasfd

  !--------------------------------------------------------------------------
  ! gpsdiffasff
  !--------------------------------------------------------------------------
  pure subroutine gpsdiffasff(gd1, gd2)
    ! :Purpose: Overloaded =
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(out)    :: gd1
    type(gps_diff)     , intent(in)     :: gd2

    gd1%Var     = gd2%Var
    gd1%DVar(:) = gd2%DVar(:)
  end subroutine gpsdiffasff

  !--------------------------------------------------------------------------
  ! gpsdiffsmfd
  !--------------------------------------------------------------------------
  pure function gpsdiffsmfd(gd1, d2)
    ! :Purpose: Overloaded +
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    real(dp)           , intent(in)     :: d2
    ! Result:
    type(gps_diff)                      :: gpsdiffsmfd

    gpsdiffsmfd%Var     = gd1%Var  + d2
    gpsdiffsmfd%DVar(:) = gd1%DVar(:)
  end function gpsdiffsmfd

  !--------------------------------------------------------------------------
  ! gpsdiffsmdd
  !--------------------------------------------------------------------------
  pure function gpsdiffsmdf(d1, gd2)
    ! :Purpose: Overloaded +
    implicit none

    ! Arguments:
    real(dp)           , intent(in)     :: d1
    type(gps_diff)     , intent(in)     :: gd2
    ! Result:
    type(gps_diff)                      :: gpsdiffsmdf

    gpsdiffsmdf%Var     = d1 + gd2%Var
    gpsdiffsmdf%DVar(:) =      gd2%DVar(:)
  end function gpsdiffsmdf

  !--------------------------------------------------------------------------
  ! gpsdiffsmfi
  !--------------------------------------------------------------------------
  pure function gpsdiffsmfi(gd1, i2)
    ! :Purpose: Overloaded +
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    integer(i4)        , intent(in)     :: i2
    ! Result:
    type(gps_diff)                      :: gpsdiffsmfi

    gpsdiffsmfi%Var     = gd1%Var  + i2
    gpsdiffsmfi%DVar(:) = gd1%DVar(:)
  end function gpsdiffsmfi

  !--------------------------------------------------------------------------
  ! gpsdiffsmif
  !--------------------------------------------------------------------------
  pure function gpsdiffsmif(i1, gd2)
    ! :Purpose: Overloaded +
    implicit none

    ! Arguments:
    integer(i4)        , intent(in)     :: i1
    type(gps_diff)     , intent(in)     :: gd2
    ! Result:
    type(gps_diff)                      :: gpsdiffsmif

    gpsdiffsmif%Var     = i1 + gd2%Var
    gpsdiffsmif%DVar(:) =      gd2%DVar(:)
  end function gpsdiffsmif

  !--------------------------------------------------------------------------
  ! gpsdiffsmff
  !--------------------------------------------------------------------------
  pure function gpsdiffsmff(gd1, gd2)
    ! :Purpose: Overloaded +
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    type(gps_diff)     , intent(in)     :: gd2
    ! Result:
    type(gps_diff)                      :: gpsdiffsmff

    gpsdiffsmff%Var     = gd1%Var     + gd2%Var
    gpsdiffsmff%DVar(:) = gd1%DVar(:) + gd2%DVar(:)
  end function gpsdiffsmff

  !--------------------------------------------------------------------------
  ! gpsdiffsbfd
  !--------------------------------------------------------------------------
  pure function gpsdiffsbfd(gd1, d2)
    ! :Purpose: Overloaded -
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    real(dp)           , intent(in)     :: d2
    ! Result:
    type(gps_diff)                      :: gpsdiffsbfd

    gpsdiffsbfd%Var     = gd1%Var  - d2
    gpsdiffsbfd%DVar(:) = gd1%DVar(:)
  end function gpsdiffsbfd

  !--------------------------------------------------------------------------
  ! gpsdiffsbdf
  !--------------------------------------------------------------------------
  pure function gpsdiffsbdf(d1, gd2)
    ! :Purpose: Overloaded -
    implicit none

    ! Arguments:
    real(dp)           , intent(in)     :: d1
    type(gps_diff)     , intent(in)     :: gd2
    ! Result:
    type(gps_diff)                      :: gpsdiffsbdf

    gpsdiffsbdf%Var     = d1 - gd2%Var
    gpsdiffsbdf%DVar(:) =    - gd2%DVar(:)
  end function gpsdiffsbdf

  !--------------------------------------------------------------------------
  ! gpsdiffsbfi
  !--------------------------------------------------------------------------
  pure function gpsdiffsbfi(gd1, i2)
    ! :Purpose: Overloaded -
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    integer(i4)        , intent(in)     :: i2
    ! Result:
    type(gps_diff)                      :: gpsdiffsbfi

    gpsdiffsbfi%Var     = gd1%Var  - i2
    gpsdiffsbfi%DVar(:) = gd1%DVar(:)
  end function gpsdiffsbfi

  !--------------------------------------------------------------------------
  ! gpsdiffsbif
  !--------------------------------------------------------------------------
  pure function gpsdiffsbif(i1, gd2)
    ! :Purpose: Overloaded -
    implicit none

    ! Arguments:
    integer(i4)        , intent(in)     :: i1
    type(gps_diff)     , intent(in)     :: gd2
    ! Result:
    type(gps_diff)                      :: gpsdiffsbif

    gpsdiffsbif%Var     = i1 - gd2%Var
    gpsdiffsbif%DVar(:) =    - gd2%DVar(:)
  end function gpsdiffsbif

  !--------------------------------------------------------------------------
  ! gpsdiffsbff
  !--------------------------------------------------------------------------
  pure function gpsdiffsbff(gd1, gd2)
    ! :Purpose: Overloaded -
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    type(gps_diff)     , intent(in)     :: gd2
    ! Result:
    type(gps_diff)                      :: gpsdiffsbff

    gpsdiffsbff%Var     = gd1%Var     - gd2%Var
    gpsdiffsbff%DVar(:) = gd1%DVar(:) - gd2%DVar(:)
  end function gpsdiffsbff

  !--------------------------------------------------------------------------
  ! gpsdiffmlfd
  !--------------------------------------------------------------------------
  pure function gpsdiffmlfd(gd1, d2)
    ! :Purpose: Overloaded *
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    real(dp)           , intent(in)     :: d2
    ! Result:
    type(gps_diff)                      :: gpsdiffmlfd

    gpsdiffmlfd%Var     = d2 * gd1%Var
    gpsdiffmlfd%DVar(:) = d2 * gd1%DVar(:)
  end function gpsdiffmlfd

  !--------------------------------------------------------------------------
  ! gpsdiffmldf
  !--------------------------------------------------------------------------
  pure function gpsdiffmldf(d1, gd2)
    ! :Purpose: Overloaded *
    implicit none

    ! Arguments:
    real(dp)           , intent(in)     :: d1
    type(gps_diff)     , intent(in)     :: gd2
    ! Result:
    type(gps_diff)                      :: gpsdiffmldf

    gpsdiffmldf%Var  = d1 * gd2%Var
    gpsdiffmldf%DVar = d1 * gd2%DVar
  end function gpsdiffmldf

  !--------------------------------------------------------------------------
  ! gpsdiffmlfi
  !--------------------------------------------------------------------------
  pure function gpsdiffmlfi(gd1, i2)
    ! :Purpose: Overloaded *
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    integer(i4)        , intent(in)     :: i2
    ! Result:
    type(gps_diff)                      :: gpsdiffmlfi

    gpsdiffmlfi%Var     = i2 * gd1%Var
    gpsdiffmlfi%DVar(:) = i2 * gd1%DVar(:)
  end function gpsdiffmlfi

  !--------------------------------------------------------------------------
  ! gpsdiffmlif
  !--------------------------------------------------------------------------
  pure function gpsdiffmlif(i1, gd2)
    ! :Purpose: Overloaded *
    implicit none

    ! Arguments:
    integer(i4)        , intent(in)     :: i1
    type(gps_diff)     , intent(in)     :: gd2
    ! Result:
    type(gps_diff)                      :: gpsdiffmlif

    gpsdiffmlif%Var  = i1 * gd2%Var
    gpsdiffmlif%DVar = i1 * gd2%DVar
  end function gpsdiffmlif

  !--------------------------------------------------------------------------
  ! gpsdiffmlff
  !--------------------------------------------------------------------------
  pure function gpsdiffmlff(gd1, gd2)
    ! :Purpose: Overloaded *
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    type(gps_diff)     , intent(in)     :: gd2
    ! Result:
    type(gps_diff)                      :: gpsdiffmlff

    gpsdiffmlff%Var     =  gd1%Var * gd2%Var
    gpsdiffmlff%DVar(:) = (gd2%Var * gd1%DVar(:)) + (gd1%Var * gd2%DVar(:))
  end function gpsdiffmlff

  !--------------------------------------------------------------------------
  ! gpsdiffdvfd
  !--------------------------------------------------------------------------
  pure function gpsdiffdvfd(gd1, d2)
    ! :Purpose: Overloaded /
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    real(dp)           , intent(in)     :: d2
    ! Result:
    type(gps_diff)                      :: gpsdiffdvfd

    gpsdiffdvfd%Var     = gd1%Var     / d2
    gpsdiffdvfd%DVar(:) = gd1%DVar(:) / d2
  end function gpsdiffdvfd

  !--------------------------------------------------------------------------
  ! gpsdiffdvdf
  !--------------------------------------------------------------------------
  pure function gpsdiffdvdf(d1, gd2)
    ! :Purpose: Overloaded /
    implicit none

    ! Arguments:
    real(dp)           , intent(in)     :: d1
    type(gps_diff)     , intent(in)     :: gd2
    ! Result:
    type(gps_diff)                      :: gpsdiffdvdf

    gpsdiffdvdf%Var     =   d1 / gd2%Var
    gpsdiffdvdf%DVar(:) = (-d1 / gd2%Var**2) * gd2%DVar(:)
  end function gpsdiffdvdf

  !--------------------------------------------------------------------------
  ! gpsdiffdvfi
  !--------------------------------------------------------------------------
  pure function gpsdiffdvfi(gd1, i2)
    ! :Purpose: Overloaded /
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    integer(i4)        , intent(in)     :: i2
    ! Result:
    type(gps_diff)                      :: gpsdiffdvfi

    gpsdiffdvfi%Var     = gd1%Var     / i2
    gpsdiffdvfi%DVar(:) = gd1%DVar(:) / i2
  end function gpsdiffdvfi

  !--------------------------------------------------------------------------
  ! gpsdiffdvif
  !--------------------------------------------------------------------------
  pure function gpsdiffdvif(i1, gd2)
    ! :Purpose: Overloaded /
    implicit none

    ! Arguments:
    integer(i4)        , intent(in)     :: i1
    type(gps_diff)     , intent(in)     :: gd2
    ! Result:
    type(gps_diff)                      :: gpsdiffdvif

    gpsdiffdvif%Var     =   i1 / gd2%Var
    gpsdiffdvif%DVar(:) = (-i1 / gd2%Var**2) * gd2%DVar(:)
  end function gpsdiffdvif

  !--------------------------------------------------------------------------
  ! gpsdiffdvff
  !--------------------------------------------------------------------------
  pure function gpsdiffdvff(gd1, gd2)
    ! :Purpose: Overloaded /
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    type(gps_diff)     , intent(in)     :: gd2
    ! Result:
    type(gps_diff)                      :: gpsdiffdvff

    ! Locals:
    real(dp)                            :: onegd2

    onegd2 = 1._dp / gd2%Var
    gpsdiffdvff%Var     = gd1%Var * onegd2
    gpsdiffdvff%DVar(:) = onegd2 * gd1%DVar(:) - (gd1%Var*onegd2*onegd2) * gd2%DVar(:)
  end function gpsdiffdvff

  !--------------------------------------------------------------------------
  ! gpsdiffpwfd
  !--------------------------------------------------------------------------
  pure function gpsdiffpwfd(gd1, d2)
    ! :Purpose: Overloaded **
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    real(dp)           , intent(in)     :: d2
    ! Result:
    type(gps_diff)                      :: gpsdiffpwfd

    gpsdiffpwfd%Var     =      gd1%Var** d2
    gpsdiffpwfd%DVar(:) = (d2*(gd1%Var**(d2-1._dp))) * gd1%DVar(:)
  end function gpsdiffpwfd

  !--------------------------------------------------------------------------
  ! gpsdiffpwdf
  !--------------------------------------------------------------------------
  pure function gpsdiffpwdf(d1, gd2)
    ! :Purpose: Overloaded **
    implicit none

    ! Arguments:
    real(dp)           , intent(in)     :: d1
    type(gps_diff)     , intent(in)     :: gd2
    ! Result:
    type(gps_diff)                      :: gpsdiffpwdf

    gpsdiffpwdf%Var     =          d1**gd2%Var
    gpsdiffpwdf%DVar(:) = (log(d1)*d1**gd2%Var) * gd2%DVar(:)
  end function gpsdiffpwdf

  !--------------------------------------------------------------------------
  ! gpsdiffpwfi
  !--------------------------------------------------------------------------
  pure function gpsdiffpwfi(gd1, i2)
    ! :Purpose: Overloaded **
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    integer(i4)        , intent(in)     :: i2
    ! Result:
    type(gps_diff)                      :: gpsdiffpwfi

    gpsdiffpwfi%Var     =      gd1%Var** i2
    gpsdiffpwfi%DVar(:) = (i2*(gd1%Var**(i2-1))) * gd1%DVar(:)
  end function gpsdiffpwfi

  !--------------------------------------------------------------------------
  ! gpsdiffpwif
  !--------------------------------------------------------------------------
  pure function gpsdiffpwif(i1, gd2)
    ! :Purpose: Overloaded **
    implicit none

    ! Arguments:
    integer(i4)        , intent(in)     :: i1
    type(gps_diff)     , intent(in)     :: gd2
    ! Result:
    type(gps_diff)                      :: gpsdiffpwif

    gpsdiffpwif%Var     =                i1**gd2%Var
    gpsdiffpwif%DVar(:) = (log(1._dp*i1)*i1**gd2%Var) * gd2%DVar(:)
  end function gpsdiffpwif

  !--------------------------------------------------------------------------
  ! gpsdiffpwff
  !--------------------------------------------------------------------------
  pure function gpsdiffpwff(gd1, gd2)
    ! :Purpose: Overloaded **
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    type(gps_diff)     , intent(in)     :: gd2
    ! Result:
    type(gps_diff)                      :: gpsdiffpwff

    gpsdiffpwff%Var     =               gd1%Var** gd2%Var
    gpsdiffpwff%DVar(:) = ( gd2%Var * ( gd1%Var**(gd2%Var-1) ) ) * gd1%DVar(:) + &
         (log(gd1%Var)*(gd1%Var**gd2%Var))*gd2%DVar(:)
  end function gpsdiffpwff

  !--------------------------------------------------------------------------
  ! gpsdiffsqr
  !--------------------------------------------------------------------------
  pure function gpsdiffsqr(gd1)
    ! :Purpose: Overloaded sqrt
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    ! Result:
    type(gps_diff)                      :: gpsdiffsqr

    gpsdiffsqr%Var     =           sqrt( gd1%Var )
    gpsdiffsqr%DVar(:) = (0.5_dp / sqrt( gd1%Var )) * gd1%DVar(:)
  end function gpsdiffsqr

  !--------------------------------------------------------------------------
  ! gpsdiffexp
  !--------------------------------------------------------------------------
  pure function gpsdiffexp(gd1)
    ! :Purpose: Overloaded exp
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    ! Result:
    type(gps_diff)                      :: gpsdiffexp

    gpsdiffexp%Var     = exp(gd1%Var)
    gpsdiffexp%DVar(:) = exp(gd1%Var) * gd1%DVar(:)
  end function gpsdiffexp

  !--------------------------------------------------------------------------
  ! gpsdifflog
  !--------------------------------------------------------------------------
  pure function gpsdifflog(gd1)
    ! :Purpose: Overloaded log
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    ! Result:
    type(gps_diff)                      :: gpsdifflog

    gpsdifflog%Var     =    log(gd1%Var)
    gpsdifflog%DVar(:) = (1._dp/gd1%Var) * gd1%DVar(:)
  end function gpsdifflog

  !--------------------------------------------------------------------------
  ! gpsdiffcos
  !--------------------------------------------------------------------------
  pure function gpsdiffcos(gd1)
    ! :Purpose: Overloaded cos
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    ! Result:
    type(gps_diff)                      :: gpsdiffcos

    gpsdiffcos%Var     =         cos(gd1%Var)
    gpsdiffcos%DVar(:) = (-1._dp*sin(gd1%Var)) * gd1%DVar(:)
  end function gpsdiffcos

  !--------------------------------------------------------------------------
  ! gpsdifftan
  !--------------------------------------------------------------------------
  pure function gpsdifftan(gd1)
    ! :Purpose: Overloaded tan
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    ! Result:
    type(gps_diff)                      :: gpsdifftan

    gpsdifftan%Var     =        tan(gd1%Var)
    gpsdifftan%DVar(:) = (1._dp/cos(gd1%Var)**2) * gd1%DVar(:)
  end function gpsdifftan

  !--------------------------------------------------------------------------
  ! gpsdiffacos
  !--------------------------------------------------------------------------
  pure function gpsdiffacos(gd1)
    ! :Purpose: Overloaded acos
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    ! Result:
    type(gps_diff)                      :: gpsdiffacos

    gpsdiffacos%Var     = acos(gd1%Var)
    gpsdiffacos%DVar(:) = (-1._dp/(1._dp-gd1%Var*gd1%Var)) * gd1%DVar(:)
  end function gpsdiffacos

  !--------------------------------------------------------------------------
  ! gpsdiffatan
  !--------------------------------------------------------------------------
  pure function gpsdiffatan(gd1)
    ! :Purpose: Overloaded atan
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    ! Result:
    type(gps_diff)                      :: gpsdiffatan

    gpsdiffatan%Var     = atan(gd1%Var)
    gpsdiffatan%DVar(:) = (1._dp/(1._dp+gd1%Var**2)) * gd1%DVar(:)
  end function gpsdiffatan

  !--------------------------------------------------------------------------
  ! gpsdifferf
  !--------------------------------------------------------------------------
  pure function gpsdifferf(gd1)
    ! :Purpose: Overloaded erf
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: gd1
    ! Result:
    type(gps_diff)                      :: gpsdifferf

    ! Locals:
    real(dp)           , parameter      :: pi = MPC_PI_R8

    gpsdifferf%Var     =                     erf( gd1%Var   )
    gpsdifferf%DVar(:) = ((2._dp/sqrt(pi)) * exp(-gd1%Var**2)) * gd1%DVar(:)
  end function gpsdifferf

  !--------------------------------------------------------------------------
  ! gps_struct1sw
  !--------------------------------------------------------------------------
  pure subroutine gps_struct1sw(ngpslev, rLat, rLon, rAzm, rMT, Rad, geoid,    &
       rP0, rPP, rDP, rTT, rHU, rUU, rVV, prf)
    ! :Purpose: Create a profile for GNSSRO from local model state
    !           Version where level altitude is derived from pressure state
    implicit none

    ! Arguments:
    integer(i4)        , intent(in)     :: ngpslev
    real(dp)           , intent(in)     :: rLat
    real(dp)           , intent(in)     :: rLon
    real(dp)           , intent(in)     :: rAzm
    real(dp)           , intent(in)     :: rMT
    real(dp)           , intent(in)     :: Rad
    real(dp)           , intent(in)     :: geoid
    real(dp)           , intent(in)     :: rP0
    real(dp)           , intent(in)     :: rPP (ngpssize)
    real(dp)           , intent(in)     :: rDP (ngpssize)
    real(dp)           , intent(in)     :: rTT (ngpssize)
    real(dp)           , intent(in)     :: rHU (ngpssize)
    real(dp)           , intent(in)     :: rUU (ngpssize)
    real(dp)           , intent(in)     :: rVV (ngpssize)
    type(gps_profile)  , intent(out)    :: prf

    ! Locals:
    integer(i4)                         :: levIndex
    real(dp)                            :: delta
    real(dp)                            :: h0, dh, Rgh, Eot, Eot2, sLat, cLat, q1, q2, q3, q4, md, mw, ym2k, wa, wb
    type(gps_diff)                      :: p, t, q, x
    type(gps_diff)                      :: tr, z
    type(gps_diff)                      :: mold, dd, dw, dx, n0, nd1, nw1, tvm
    type(gps_diff)                      :: cmp(ngpssize), xi(ngpssize), tv(ngpssize)

    if (gps_roNFlavour == 0) then
      md = gps_p_md
      mw = gps_p_mw
      wa = gps_p_wa
      wb = gps_p_wb
      delta = 0.6077686814144_dp
      q1 =  222.682_dp
      q2 =    0.069_dp
      q3 = 6701.605_dp
      q4 = 6385.886_dp
    else
      ym2k = gps_roYear - 2000._dp
      md = 28.96496_dp +    1.30d-5*ym2k + 4.41d-8 * ym2k*ym2k
      mw = 18.01525_dp
      wa = md/mw
      wb = (md-mw)/mw
      delta = wb
      q1 =  222.654_dp + 0.000259_dp*ym2k + 2.24d-6 * ym2k*ym2k
      q2 =    0.097_dp
      q3 = 6703.497_dp
      q4 = 6393.484_dp
    end if

    prf%ngpslev = ngpslev
    prf%rLat    = rLat
    prf%rLon    = rLon
    prf%rAzm    = rAzm
    prf%rMT     = rMT
    prf%Rad     = Rad
    prf%geoid   = geoid
    call phf_gpsRadii(rLat, prf%RadN, prf%RadM)

    !
    ! Fill pressure placeholders:
    !
    prf%P0%Var               = 0.01_dp*rP0
    prf%P0%DVar              = 0._dp
    prf%P0%DVar(2*ngpslev+1) = 0.01_dp
    do levIndex = 1, ngpslev
      prf%pst(levIndex)%Var                    = 0.01_dp*rPP(levIndex)
      prf%pst(levIndex)%DVar                   = 0._dp
      prf%pst(levIndex)%DVar(2*ngpslev+1)      = 0.01_dp*rDP(levIndex)
    end do

    !
    ! Fill temperature placeholders:
    !
    do levIndex = 1, ngpslev
      prf%tst(levIndex)%Var                    = rTT(levIndex)+MPC_K_C_DEGREE_OFFSET_R8
      prf%tst(levIndex)%DVar                   = 0._dp
      prf%tst(levIndex)%DVar(levIndex)         = 1._dp
    end do

    !
    ! Fill moisture placeholders:
    !
    do levIndex = 1, ngpslev
      prf%qst(levIndex)%Var                    = rHU(levIndex)
      prf%qst(levIndex)%DVar                   = 0._dp
      prf%qst(levIndex)%DVar(ngpslev+levIndex) = 1._dp
    end do

    ! Compressibility:
    do levIndex = 1, ngpslev
      cmp(levIndex)= gps_compressibility(prf%pst(levIndex), prf%tst(levIndex), prf%qst(levIndex))
    end do

    ! Refractivity:
    do levIndex = 1, ngpslev
      p  = prf%pst(levIndex)
      t  = prf%tst(levIndex)
      q  = prf%qst(levIndex)
      x  = wa*q/(1._dp+wb*q)

      ! Densities (molar, total, dry, water vapor):
      mold  = p/t * (100._dp/(p_R*cmp(levIndex)))                           ! p in hPa
      dd = mold * (1._dp-x) * (md/1000._dp)
      dw = mold * x         * (mw/1000._dp)
      ! Aparicio (2011) expression
      tr = MPC_K_C_DEGREE_OFFSET_R8/t-1._dp
      nd1= (q1+q2*tr) * dd
      nw1= (q3+q4*tr) * dw
      n0 = (nd1+nw1)
      prf%rst(levIndex) = n0*(1._dp+(1.e-6_dp/6._dp)*n0)
    end do

    ! Midpoint refractivity between layers (2 to ngpslev)
    do levIndex = 2, ngpslev
      p  = exp( 0.5_dp*( log(prf%pst(levIndex-1))+log(prf%pst(levIndex)) ) )
      t  =      0.5_dp*(     prf%tst(levIndex-1) +    prf%tst(levIndex))
      q  = exp( 0.5_dp*( log(prf%qst(levIndex-1))+log(prf%qst(levIndex)) ) )
      x  = wa*q/(1._dp+wb*q)

      ! Densities (molar, total, dry, water vapor):
      mold  = p/t * (100._dp/(p_R*0.5_dp*(cmp(levIndex-1)+cmp(levIndex))))  ! p in hPa
      dd = mold * (1._dp-x) * (md/1000._dp)
      dw = mold * x         * (mw/1000._dp)
      ! Aparicio (2011) expression
      tr = MPC_K_C_DEGREE_OFFSET_R8/t-1._dp
      nd1= (q1+q2*tr) * dd
      nw1= (q3+q4*tr) * dw
      n0 = (nd1+nw1)
      prf%lrmd(levIndex) = log(n0*(1._dp+(1.e-6_dp/6._dp)*n0))
    end do
    prf%lrmd(1) = prf%lrmd(2)+log(prf%pst(1))-log(prf%pst(2))

    !
    ! Hydrostatic equation
    !
    do levIndex = 1, ngpslev
      p = prf%pst(levIndex)
      t = prf%tst(levIndex)
      q = prf%qst(levIndex)
      !
      ! Log(P)
      !
      xi(levIndex) = log(p)
      !
      ! Virtual temperature (K) (corrected of compressibility)
      !
      tv(levIndex) = (1._dp+delta*q) * t * cmp(levIndex)
      prf%vst(levIndex) = tv(levIndex)
    end do

    sLat=sin(rLat)
    cLat=cos(rLat)
    dx  = xi(ngpslev)-log(prf%P0)
    Rgh = phf_gravitysrf(sLat)
    z   = (-p_Rd/Rgh) * tv(ngpslev) * dx
    prf%gst(ngpslev) = rMT + z
    do levIndex = ngpslev-1, 1, -1
      dx = xi(levIndex)-xi(levIndex+1)
      tvm = 0.5_dp*(tv(levIndex)+tv(levIndex+1))
      !
      ! Gravity acceleration (includes 2nd-order Eotvos effect)
      !
      h0  = prf%gst(levIndex+1)%Var
      Eot = 2*ec_wgs_OmegaPrime*cLat*rUU(levIndex)
      Eot2= (rUU(levIndex)**2+rVV(levIndex)**2)/ec_wgs_a
      Rgh = phf_gravityalt(sLat, h0)-Eot-Eot2
      dh  = (-p_Rd/Rgh) * tvm%Var * dx%Var
      Rgh = phf_gravityalt(sLat, h0+0.5_dp*dh)-Eot-Eot2
      !
      ! Height increment
      !
      z   = (-p_Rd/Rgh) * tvm * dx
      prf%gst(levIndex) = prf%gst(levIndex+1) + z
    end do

  end subroutine gps_struct1sw

  !--------------------------------------------------------------------------
  ! gps_struct1sw_v2
  !--------------------------------------------------------------------------
  pure subroutine gps_struct1sw_v2(ngpslev, rLat, rLon, rAzm, rMT, Rad, geoid,    &
       rP0, rPP, rTT, rHU, rUU, rVV, rALT, prf)
    ! :Purpose: Create a profile for GNSSRO from local model state
    !           Version where level altitude is a state variable
    implicit none

    ! Arguments:
    integer(i4)        , intent(in)     :: ngpslev
    real(dp)           , intent(in)     :: rLat
    real(dp)           , intent(in)     :: rLon
    real(dp)           , intent(in)     :: rAzm
    real(dp)           , intent(in)     :: rMT
    real(dp)           , intent(in)     :: Rad
    real(dp)           , intent(in)     :: geoid
    real(dp)           , intent(in)     :: rP0
    real(dp)           , intent(in)     :: rPP (ngpssize)
    real(dp)           , intent(in)     :: rTT (ngpssize)
    real(dp)           , intent(in)     :: rHU (ngpssize)
    real(dp)           , intent(in)     :: rUU (ngpssize)
    real(dp)           , intent(in)     :: rVV (ngpssize)
    real(dp)           , intent(in)     :: rALT (ngpssize)
    type(gps_profile)  , intent(out)    :: prf

    ! Locals:
    integer(i4)                         :: levIndex
    real(dp)                            :: rALT_E(ngpssize), q1, q2, q3, q4, md, mw, wa, wb, ym2k
    real(dp)                            :: delta
    type(gps_diff)                      :: cmp(ngpssize)
    type(gps_diff)                      :: p, t, q, x
    type(gps_diff)                      :: tr, tv
    type(gps_diff)                      :: mold, dd, dw, n0, nd1, nw1

    if (gps_roNFlavour == 0) then
      md = gps_p_md
      mw = gps_p_mw
      wa = gps_p_wa
      wb = gps_p_wb
      delta = 0.6077686814144_dp
      q1 =  222.682_dp
      q2 =    0.069_dp
      q3 = 6701.605_dp
      q4 = 6385.886_dp
    else
      ym2k = gps_roYear - 2000._dp
      md = 28.96496_dp +    1.30d-5*ym2k + 4.41d-8 * ym2k*ym2k
      mw = 18.01525_dp
      wa = md/mw
      wb = (md-mw)/mw
      delta = wb
      q1 =  222.654_dp + 0.000259_dp*ym2k + 2.24d-6 * ym2k*ym2k
      q2 =    0.097_dp
      q3 = 6703.497_dp
      q4 = 6393.484_dp
    end if

    prf%ngpslev = ngpslev
    prf%rLat    = rLat
    prf%rLon    = rLon
    prf%rAzm    = rAzm
    prf%rMT     = rMT
    prf%Rad     = Rad
    prf%geoid   = geoid
    call phf_gpsRadii(rLat, prf%RadN, prf%RadM)

    !
    ! Fill pressure placeholders:
    !
    prf%P0%Var                                   = 0.01_dp*rP0
    prf%P0%DVar                                  = 0._dp
    prf%P0%DVar(4*ngpslev)                       = 0.01_dp
    do levIndex = 1, ngpslev
      prf%pst(levIndex)%Var                      = 0.01_dp*rPP(levIndex)
      prf%pst(levIndex)%DVar                     = 0._dp
      prf%pst(levIndex)%DVar(3*ngpslev+levIndex) = 0.01_dp
    end do

    !
    ! Fill temperature placeholders:
    !
    do levIndex = 1, ngpslev
      prf%tst(levIndex)%Var                      = rTT(levIndex)+MPC_K_C_DEGREE_OFFSET_R8
      prf%tst(levIndex)%DVar                     = 0._dp
      prf%tst(levIndex)%DVar(levIndex)           = 1._dp
    end do

    !
    ! Fill moisture placeholders:
    !
    do levIndex = 1, ngpslev
      prf%qst(levIndex)%Var                      = rHU(levIndex)
      prf%qst(levIndex)%DVar                     = 0._dp
      prf%qst(levIndex)%DVar(ngpslev+levIndex)   = 1._dp
    end do

    !
    ! Fill altitude placeholders:
    !
    if (gps_roEotvos) then
      call gpsro_Eotvos_dH(ngpslev, rLat, rALT, rUU, rVV, rALT_E)
    else
      rALT_E(1:ngpslev) = rALT(1:ngpslev)
    end if
    do levIndex = 1, ngpslev
      prf%gst(levIndex)%Var                      = rALT_E(levIndex)
      prf%gst(levIndex)%DVar                     = 0._dp
      prf%gst(levIndex)%DVar(2*ngpslev+levIndex) = 1._dp
    end do

    ! Compressibility:
    do levIndex = 1, ngpslev
      cmp(levIndex)= gps_compressibility(prf%pst(levIndex), prf%tst(levIndex), prf%qst(levIndex))
    end do

    ! Refractivity:
    do levIndex = 1, ngpslev
      p  = prf%pst(levIndex)
      t  = prf%tst(levIndex)
      q  = prf%qst(levIndex)
      x  = wa*q/(1._dp+wb*q)
      ! Densities (molar, total, dry, water vapor):
      mold  = p/t * (100._dp/(p_R*cmp(levIndex)))                            ! p in hPa
      dd = mold * (1._dp-x) * (md/1000._dp)
      dw = mold * x         * (mw/1000._dp)
      ! Aparicio (2011) expression
      tr = MPC_K_C_DEGREE_OFFSET_R8/t-1._dp
      nd1= (q1+q2*tr) * dd
      nw1= (q3+q4*tr) * dw
      n0 = (nd1+nw1)
      prf%rst(levIndex) = n0*(1._dp+(1.e-6_dp/6._dp)*n0)
    end do

    ! Midpoint refractivity between layers (2 to ngpslev)
    do levIndex = 2, ngpslev
      p  = exp( 0.5_dp*( log(prf%pst(levIndex-1))+log(prf%pst(levIndex)) ) )
      t  =      0.5_dp*(     prf%tst(levIndex-1) +    prf%tst(levIndex))
      q  = exp( 0.5_dp*( log(prf%qst(levIndex-1))+log(prf%qst(levIndex)) ) )
      x  = wa*q/(1._dp+wb*q)

      ! Densities (molar, total, dry, water vapor):
      mold  = p/t * (100._dp/(p_R*0.5_dp*(cmp(levIndex-1)+cmp(levIndex))))   ! p in hPa
      dd = mold * (1._dp-x) * (md/1000._dp)
      dw = mold * x         * (mw/1000._dp)
      ! Aparicio (2011) expression
      tr = MPC_K_C_DEGREE_OFFSET_R8/t-1._dp
      nd1= (q1+q2*tr) * dd
      nw1= (q3+q4*tr) * dw
      n0 = (nd1+nw1)
      prf%lrmd(levIndex) = log(n0*(1._dp+(1.e-6_dp/6._dp)*n0))
    end do
    prf%lrmd(1) = prf%lrmd(2)+log(prf%pst(1))-log(prf%pst(2))

    !
    ! Virtual temperature
    !
    do levIndex = 1, ngpslev
      t = prf%tst(levIndex)
      q = prf%qst(levIndex)
      !
      ! Virtual temperature (K) corrected for compressibility
      !
      tv = (1._dp+delta*q) * t * cmp(levIndex)
      prf%vst(levIndex) = tv
    end do
  end subroutine gps_struct1sw_v2

  !--------------------------------------------------------------------------
  ! gps_compressibility
  !--------------------------------------------------------------------------
  pure function gps_compressibility(p, t, q)
    ! :Purpose: Air compressibility factor, as a function of p, t, q
    implicit none

    ! Arguments:
    type(gps_diff)     , intent(in)     :: p
    type(gps_diff)     , intent(in)     :: t
    type(gps_diff)     , intent(in)     :: q
    ! Result:
    type(gps_diff)                      :: gps_compressibility

    ! Locals:
    real(dp)           , parameter      :: a0 = 1.58123e-6_dp
    real(dp)           , parameter      :: a1 =-2.9331e-8_dp
    real(dp)           , parameter      :: a2 = 1.1043e-10_dp
    real(dp)           , parameter      :: b0 = 5.707e-6_dp
    real(dp)           , parameter      :: b1 =-2.051e-8_dp
    real(dp)           , parameter      :: c0 = 1.9898e-4_dp
    real(dp)           , parameter      :: c1 =-2.376e-6_dp
    real(dp)           , parameter      :: d  = 1.83e-11_dp
    real(dp)           , parameter      :: e  =-0.765e-8_dp
    type(gps_diff)                      :: x, tc, pt, tc2, x2

    x  = gps_p_wa*q/(1._dp+gps_p_wb*q)
    ! Estimate, from CIPM, Picard (2008)
    tc = t-MPC_K_C_DEGREE_OFFSET_R8
    pt = 1.e2_dp*p/t
    tc2= tc*tc
    x2 = x*x
    gps_compressibility = 1._dp-pt*(a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2)+pt*pt*(d+e*x2)
  end function gps_compressibility

  !--------------------------------------------------------------------------
  ! gpsro_Eotvos_dH
  !--------------------------------------------------------------------------
  pure subroutine gpsro_Eotvos_dH(ngpslev, rLat, rALT, rUU, rVV, rALT_E)
    ! :Purpose: Evaluate Eotvos vertical shift
    implicit none

    ! Arguments:
    integer(i4)        , intent(in)     :: ngpslev
    real(dp)           , intent(in)     :: rLat
    real(dp)           , intent(in)     :: rALT(ngpslev)
    real(dp)           , intent(in)     :: rUU (ngpslev)
    real(dp)           , intent(in)     :: rVV (ngpslev)
    real(dp)           , intent(out)    :: rALT_E(ngpslev)

    ! Locals:
    integer(i4)                         :: levIndex
    real(dp)                            :: cLat, dALT, Eot, Eot2, dALTE, ddAL, acc

    cLat=cos(rLat)
    rALT_E(ngpslev) = rALT(ngpslev)
    acc = 0._dp
    do levIndex = ngpslev-1, 1, -1
      dALT = rALT(levIndex) - rALT(levIndex+1)
      Eot = 2*ec_wgs_OmegaPrime*cLat*rUU(levIndex)
      Eot2= (rUU(levIndex)**2+rVV(levIndex)**2)/ec_wgs_a
      dALTE = dALT*(1._dp+(Eot+Eot2)/ec_rg)
      ddAL = dALTE - dALT
      acc = acc + ddAL
      rALT_E(levIndex) = rALT(levIndex) + acc
    end do
  end subroutine gpsro_Eotvos_dH

  !--------------------------------------------------------------------------
  ! gps_structztd
  !--------------------------------------------------------------------------
  subroutine gps_structztd(ngpslev, rLat, rLon, rMT, rP0, rPP, rDP, rTT, rHU, lBevis, &
                           refopt, prf)
    !
    ! :Purpose: This subroutine fills GPS profiles of type gps_profilezd (for ZTD
    !           operator)
    !
    ! :Arguments:
    !     :refopt:
    !               =1 --> use conventional expression for refractivity N
    !
    !               =2 --> use new Aparicio & Laroche refractivity N
    implicit none

    ! Arguments:
    integer(i4)        , intent(in)     :: ngpslev          ! number of profile levels
    real(dp)           , intent(in)     :: rLat             ! radians
    real(dp)           , intent(in)     :: rLon             ! radians
    real(dp)           , intent(in)     :: rMT              ! height (ASL) of model surface (m)
    real(dp)           , intent(in)     :: rP0              ! surface pressure (Pa)
    real(dp)           , intent(in)     :: rPP (ngpssize)   ! pressure P at each level (Pa)
    real(dp)           , intent(in)     :: rDP (ngpssize)   ! dP/dP0 at each level (Pa/Pa)
    real(dp)           , intent(in)     :: rTT (ngpssize)   ! temperature T at each level (C)
    real(dp)           , intent(in)     :: rHU (ngpssize)   ! q at each level
    logical            , intent(in)     :: lBevis           ! determines refractivity constants to use (Bevis or Rueger)
    integer(i4)        , intent(in)     :: refopt
    type(gps_profilezd), intent(out)    :: prf

    ! Locals:
    ! ******** PARAMETERS *************
    real(dp)           , parameter      :: delta = 0.6077686814144_dp
    real(dp)           , parameter      :: eps   = 0.6219800221014_dp
    ! Reuger (2002) refractivity constants (MKS units)
    real(dp)           , parameter      :: k1r = 0.776890_dp
    real(dp)           , parameter      :: k2r = 0.712952_dp
    real(dp)           , parameter      :: k3r = 3754.63_dp
    ! Bevis (1994) refractivity constants (MKS units)
    real(dp)           , parameter      :: k1b = 0.776000_dp
    real(dp)           , parameter      :: k2b = 0.704000_dp
    real(dp)           , parameter      :: k3b = 3739.000_dp
    ! ******** VARIABLES *************
    type(gps_diff)                      :: tr
    type(gps_diff)                      :: mold, dd, dw, dx, n0, nd1, nw1
    integer(i4)                         :: levIndex
    real(dp)                            :: k1, k2, k3, k2p
    real(dp)                            :: h0, dh, Rgh, sLat, ptop
    type(gps_diff)                      :: p, t, q, x, na, tvm, z
    type(gps_diff)                      :: xi(ngpssize), tv(ngpssize), cmp(ngpssize), N(ngpssize)

    prf%ngpslev = ngpslev
    prf%rLat    = rLat
    prf%rLon    = rLon
    prf%rMT     = rMT
    prf%bpst    = .false.

    ! Fill pressure (P) placeholders (Pa):
    prf%P0%Var               = rP0
    prf%P0%DVar              = 0._dp
    prf%P0%DVar(2*ngpslev+1) = 1._dp
    do levIndex = 1, ngpslev
      prf%pst(levIndex)%Var                    = rPP(levIndex)
      prf%pst(levIndex)%DVar                   = 0._dp
      prf%pst(levIndex)%DVar(2*ngpslev+1)      = rDP(levIndex)
    end do
    ! Pressure at model top (Pa)
    ptop = rPP(1)
    prf%bpst = .true.

    ! Fill temperature (T) placeholders (C--> K):
    do levIndex = 1, ngpslev
      prf%tst(levIndex)%Var                    = rTT(levIndex)+MPC_K_C_DEGREE_OFFSET_R8
      prf%tst(levIndex)%DVar                   = 0._dp
      prf%tst(levIndex)%DVar(levIndex)         = 1._dp
    end do

    ! Fill moisture (Q) placeholders (kg/kg):
    do levIndex = 1, ngpslev
      prf%qst(levIndex)%Var                    = rHU(levIndex)
      prf%qst(levIndex)%DVar                   = 0._dp
      prf%qst(levIndex)%DVar(ngpslev+levIndex) = 1._dp
    end do

    if ( refopt == 2 ) then  ! use Aparicio & Laroche refractivity
      ! Compressibility:
      do levIndex = 1, ngpslev
        cmp(levIndex)= gps_compressibility(prf%pst(levIndex), prf%tst(levIndex), prf%qst(levIndex))
      end do

      ! Refractivity:
      do levIndex = 1, ngpslev
        p  = prf%pst(levIndex)
        t  = prf%tst(levIndex)
        q  = prf%qst(levIndex)
        x  = gps_p_wa*q/(1._dp+gps_p_wb*q)

        ! Densities (molar, total, dry, water vapor):
        mold  = p/(p_R*t*cmp(levIndex))
        dd = mold * (1._dp-x) * (gps_p_md/1000._dp)
        dw = mold * x         * (gps_p_mw/1000._dp)

        ! Aparicio (2011) expression
        tr = MPC_K_C_DEGREE_OFFSET_R8/t-1._dp
        nd1= ( 222.682_dp+   0.069_dp*tr) * dd
        nw1= (6701.605_dp+6385.886_dp*tr) * dw
        n0 = (nd1+nw1)
        na = n0*(1._dp+1.e-6_dp*n0/6._dp)
        N(levIndex) = na
      end do
    end if

    ! Refractivity constants
    if ( lBevis ) then
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
    do levIndex = 1, ngpslev
      p = prf%pst(levIndex)
      t = prf%tst(levIndex)
      q = prf%qst(levIndex)
      xi(levIndex) = log(p)
      tv(levIndex) = (1._dp+delta*q) * t
    end do

    ! Geometric height (m) profile from lowest model level to top  --> prf%gst
    sLat = sin(rLat)
    dx  = xi(ngpslev)-log(prf%P0)
    Rgh = phf_gravitysrf(sLat)
    z   = (-p_Rd/Rgh) * tv(ngpslev) * dx
    prf%gst(ngpslev) = rMT + z
    do levIndex = ngpslev-1, 1, -1
      dx = xi(levIndex)-xi(levIndex+1)
      tvm = 0.5_dp*(tv(levIndex)+tv(levIndex+1))

      ! Gravity acceleration
      h0  = prf%gst(levIndex+1)%Var
      Rgh = phf_gravityalt(sLat, h0)
      dh  = (-p_Rd/Rgh) * tvm%Var * dx%Var
      Rgh = phf_gravityalt(sLat, h0+0.5_dp*dh)

      ! Height increment (m)
      z   = (-p_Rd/Rgh) * tvm * dx
      prf%gst(levIndex) = prf%gst(levIndex+1) + z
    end do

    ! Profile of dZTD/dp --> prf%rst
    do levIndex = 1, ngpslev
      p  = prf%pst(levIndex)
      t  = prf%tst(levIndex)
      q  = prf%qst(levIndex)
      if ( refopt == 1 ) then
        na = (k1/tv(levIndex)) + (k2p*(q/(eps*t))) + (k3*(q/(eps*t**2)))
      else
        na = N(levIndex) / p
      end if
      prf%rst(levIndex) = 1.e-6_dp * na * (p_Rd*tv(levIndex))/phf_gravityalt(sLat, prf%gst(levIndex)%Var)
    end do

    ! ZTD (m) profile from model top down to lowest model level --> prf%ztd
    prf%ztd(1) = 1.e-6_dp * ((k1*p_Rd*ptop)/(phf_gravityalt(sLat, prf%gst(1)%Var)))
    do levIndex = 2, ngpslev
      ! ZTD increment = Avg(dZTD/dP) * delta_P
      z = ((prf%rst(levIndex-1) + prf%rst(levIndex))/2._dp) * (prf%pst(levIndex)-prf%pst(levIndex-1))
      prf%ztd(levIndex) = prf%ztd(levIndex-1) + z
    end do
  end subroutine gps_structztd

  !--------------------------------------------------------------------------
  ! gps_structztd_v2
  !--------------------------------------------------------------------------
  subroutine gps_structztd_v2(ngpslev, rLat, rLon, rMT, rP0, rPP, rTT, rHU, rALT,   &
                              lBevis, refopt, prf)
    !
    ! :Purpose: This subroutine fills GPS profiles of type gps_profilezd (for ZTD
    !           operator)
    !
    ! :Arguments:
    !     :refopt:  =1 --> use conventional expression for refractivity N
    !
    !               =2 --> use new Aparicio & Laroche refractivity N
    implicit none

    ! Arguments:
    integer(i4)        , intent(in)     :: ngpslev          ! number of profile levels
    real(dp)           , intent(in)     :: rLat             ! radians
    real(dp)           , intent(in)     :: rLon             ! radians
    real(dp)           , intent(in)     :: rMT              ! height (ASL) of model surface (m)
    real(dp)           , intent(in)     :: rP0              ! surface pressure (Pa)
    real(dp)           , intent(in)     :: rPP (ngpssize)   ! pressure P at each level (Pa)
    real(dp)           , intent(in)     :: rTT (ngpssize)   ! temperature T at each level (C)
    real(dp)           , intent(in)     :: rHU (ngpssize)   ! q at each level
    real(dp)           , intent(in)     :: rALT (ngpssize)  ! altitude at each level
    logical            , intent(in)     :: lBevis           ! determines refractivity constants to use (Bevis or Rueger)
    integer(i4)        , intent(in)     :: refopt
    type(gps_profilezd), intent(out)    :: prf

    ! Locals:
    ! ******** PARAMETERS *************
    real(dp)           , parameter      :: delta = 0.6077686814144_dp
    real(dp)           , parameter      :: eps   = 0.6219800221014_dp
    ! Reuger (2002) refractivity constants (MKS units)
    real(dp)           , parameter      :: k1r = 0.776890_dp
    real(dp)           , parameter      :: k2r = 0.712952_dp
    real(dp)           , parameter      :: k3r = 3754.63_dp
    ! Bevis (1994) refractivity constants (MKS units)
    real(dp)           , parameter      :: k1b = 0.776000_dp
    real(dp)           , parameter      :: k2b = 0.704000_dp
    real(dp)           , parameter      :: k3b = 3739.000_dp
    ! ******** VARIABLES *************
    type(gps_diff)                      :: tr
    type(gps_diff)                      :: mold, dd, dw, n0, nd1, nw1
    integer(i4)                         :: levIndex
    real(dp)                            :: k1, k2, k3, k2p
    real(dp)                            :: sLat, ptop
    type(gps_diff)                      :: p, t, q, x, na, z
    type(gps_diff)                      :: tv(ngpssize), cmp(ngpssize), N(ngpssize)

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
    do levIndex = 1, ngpslev
      prf%pst(levIndex)%Var                      = rPP(levIndex)
      prf%pst(levIndex)%DVar                     = 0._dp
      prf%pst(levIndex)%DVar(3*ngpslev+levIndex) = 1._dp
    end do
    ! Pressure at model top (Pa)
    ptop = rPP(1)
    prf%bpst = .true.
    !
    ! Fill temperature (T) placeholders (C--> K):
    !
    do levIndex = 1, ngpslev
      prf%tst(levIndex)%Var                      = rTT(levIndex)+MPC_K_C_DEGREE_OFFSET_R8
      prf%tst(levIndex)%DVar                     = 0._dp
      prf%tst(levIndex)%DVar(levIndex)           = 1._dp
    end do

    !
    ! Fill moisture (Q) placeholders (kg/kg):
    !
    do levIndex = 1, ngpslev
      prf%qst(levIndex)%Var                      = rHU(levIndex)
      prf%qst(levIndex)%DVar                     = 0._dp
      prf%qst(levIndex)%DVar(ngpslev+levIndex)   = 1._dp
    end do

    !
    ! Fill altitude (AL) placeholders (m):
    !
    do levIndex = 1, ngpslev
      prf%gst(levIndex)%Var                      = rALT(levIndex)
      prf%gst(levIndex)%DVar                     = 0._dp
      prf%gst(levIndex)%DVar(2*ngpslev+levIndex) = 1._dp
    end do

    if ( refopt == 2 ) then  ! use Aparicio & Laroche refractivity
      ! Compressibility:
      do levIndex = 1, ngpslev
        cmp(levIndex)= gps_compressibility(prf%pst(levIndex), prf%tst(levIndex), prf%qst(levIndex))
      end do

      ! Refractivity:
      do levIndex = 1, ngpslev
        p  = prf%pst(levIndex)
        t  = prf%tst(levIndex)
        q  = prf%qst(levIndex)
        x  = gps_p_wa*q/(1._dp+gps_p_wb*q)

        ! Densities (molar, total, dry, water vapor):
        mold  = p/(p_R*t*cmp(levIndex))
        dd = mold * (1._dp-x) * (gps_p_md/1000._dp)
        dw = mold * x         * (gps_p_mw/1000._dp)
        ! Aparicio (2011) expression
        tr = MPC_K_C_DEGREE_OFFSET_R8/t-1._dp
        nd1= ( 222.682_dp+   0.069_dp*tr) * dd
        nw1= (6701.605_dp+6385.886_dp*tr) * dw
        n0 = (nd1+nw1)
        na = n0*(1._dp+1.e-6_dp*n0/6._dp)
        N(levIndex) = na
      end do
    end if

    ! Refractivity constants
    if ( lBevis ) then
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
    do levIndex = 1, ngpslev
      p = prf%pst(levIndex)
      t = prf%tst(levIndex)
      q = prf%qst(levIndex)
      tv(levIndex) = (1._dp+delta*q) * t
    end do

    sLat = sin(rLat)

    ! Profile of dZTD/dp --> prf%rst
    do levIndex = 1, ngpslev
      p  = prf%pst(levIndex)
      t  = prf%tst(levIndex)
      q  = prf%qst(levIndex)
      if ( refopt == 1 ) then
        na = (k1/tv(levIndex)) + (k2p*(q/(eps*t))) + (k3*(q/(eps*t**2)))
      else
        na = N(levIndex) / p
      end if
      prf%rst(levIndex) = 1.e-6_dp * na * (p_Rd*tv(levIndex))/phf_gravityalt(sLat, prf%gst(levIndex)%Var)
    end do

    ! ZTD (m) profile from model top down to lowest model level --> prf%ztd
    prf%ztd(1) = 1.e-6_dp * ((k1*p_Rd*ptop)/(phf_gravityalt(sLat, prf%gst(1)%Var)))
    do levIndex = 2, ngpslev
      !
      ! ZTD increment = Avg(dZTD/dP) * delta_P
      !
      z = ((prf%rst(levIndex-1) + prf%rst(levIndex))/2._dp) * (prf%pst(levIndex)-prf%pst(levIndex-1))
      prf%ztd(levIndex) = prf%ztd(levIndex-1) + z
    end do
  end subroutine gps_structztd_v2

  !--------------------------------------------------------------------------
  ! gps_refopv
  !--------------------------------------------------------------------------
  pure subroutine gps_refopv(hv, nval, prf, refopv)
    !
    ! :Purpose: GPSRO Refractivity operator
    !
    implicit none

    ! Arguments:
    real(dp)           , intent(in)     :: hv(:)     ! an array of height values
    integer(i4)        , intent(in)     :: nval      ! number of values in the array
    type(gps_profile)  , intent(in)     :: prf       ! local profile
    type(gps_diff)     , intent(out)    :: refopv(:) ! an array of refractivity values (with derivatives)

    ! Locals:
    integer(i4)                         :: iSize, iObsIndex, ngpslev
    integer(i4)                         :: levIndex, jloc
    real(dp)                            :: h
    type(gps_diff)                      :: dz
    type(gps_diff)                      :: dzm
    type(gps_diff)                      :: dzp
    type(gps_diff)                      :: m, dt, tav2

    ngpslev=prf%ngpslev
    iSize = size(hv)
    if (nval < iSize) iSize=nval
    !
    ! Given a height
    !
    jloc = 1
    do iObsIndex = 1, iSize
      h = hv(iObsIndex)
      !
      ! Search where it is located
      !
      if (h > prf%gst(1)%Var) then
        jloc = 1
      end if

      do levIndex = 1, ngpslev-1
        if ((h <= prf%gst(levIndex)%Var) .and. (h > prf%gst(levIndex+1)%Var)) then
          jloc = levIndex
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
        ! Either quasi-linear-log interpolation
        !
        dz  = prf%gst(jloc) - prf%gst(jloc+1)

        dzm = h - prf%gst(jloc+1)
        dzp = prf%gst(jloc) - h

        if (.not. gps_roNCurv) then
          ! Perfect linear-log (zero curvature)
          refopv(iObsIndex) = exp( (dzm * log(prf%rst(jloc)) + dzp * log(prf%rst(jloc+1))) / dz )
        else
          ! Quasi-linear-log, with fixed curvature m (depends on dT/dz)
          dt = prf%tst(jloc) - prf%tst(jloc+1)
          tav2 = prf%tst(jloc) * prf%tst(jloc+1)
          m = ec_wgs_GammaM * dt/dz /(2*MPC_RGAS_DRY_AIR_R8*tav2)
          refopv(iObsIndex) = exp( (dzm * log(prf%rst(jloc)) + dzp * log(prf%rst(jloc+1))) / dz + m*dzp*dzm)
        end if
      else
        !
        ! Or exp extrapolation at the lower edge
        ! (better standard exp profile than linear-log, which may be unstable)
        !
        dzm = h - prf%gst(jloc+1)
        refopv(iObsIndex) = prf%rst(jloc+1) * exp((-1._dp/6500._dp)*dzm)
      end if
    end do
  end subroutine gps_refopv

  !--------------------------------------------------------------------------
  ! gps_pwopv
  !--------------------------------------------------------------------------
  subroutine gps_pwopv(prf, pwst)
    !
    ! :Purpose: Precipitable WV operator
    !
    implicit none

    ! Arguments:
    type(gps_profile)  , intent(in)     :: prf       ! local profile
    type(gps_diff)     , intent(out)    :: pwst      ! the PWV (with derivatives)

    ! Locals:
    integer(i4)                         :: levIndex, ngpslev
    real(dp)                            :: mw, wa, wb
    type(gps_diff)                      :: p, t, q, cmp, x, mold, dw(ngpssize)
    type(gps_diff)                      :: tcwv1, tcwv2, gp, gm, dwp, dwm, k

    mw = gps_p_mw
    wa = gps_p_wa
    wb = gps_p_wb
    ngpslev=prf%ngpslev

    ! Vertical integral
    ! PW density
    do levIndex = 1, ngpslev
      p  = prf%pst(levIndex)
      t  = prf%tst(levIndex)
      q  = prf%qst(levIndex)
      ! Densities (molar, water vapor):
      cmp= gps_compressibility(p, t, q)
      x  = wa*q/(1._dp+wb*q)
      mold  = p/t * (100._dp/(p_R*cmp))               ! p in hPa
      dw(levIndex) = mold * x         * (mw/1000._dp)
    end do

    ! Vertical integral
    tcwv1 = 0._dp
    tcwv2 = 0._dp
    do levIndex = 1, ngpslev-1
      gp  = prf%gst(levIndex)
      gm  = prf%gst(levIndex+1)
      dwp = dw(levIndex)
      dwm = dw(levIndex+1)
      k = (log(dwp)-log(dwm))/(gp-gm)
      tcwv1 = tcwv1 + 0.5_dp*(dwm+dwp)*(gp-gm)
      tcwv2 = tcwv2 + (dwp-dwm)/k
   end do
   write(*,*)'gps_pwopv', 'TCWV', levIndex, tcwv1%Var, tcwv2%Var
   pwst = tcwv1
  end subroutine gps_pwopv

  !--------------------------------------------------------------------------
  ! gps_ztdopv
  !--------------------------------------------------------------------------
  pure subroutine gps_ztdopv(hv, prf, lBevis, dzmin, ZTDopv, rPobs, mode)
    !
    ! :Purpose: GB-GPS ZTD operator
    !
    ! :Arguments:
    !   :dzmin:   Minimum DZ = Zobs-Zmod (m) for which DZ adjustment to ZTD will
    !             be made when Zobs < Zmod.
    !   :mode:    1 = normal mode: use stored ZTD profiles
    !
    !             2 = Vedel & Huang ZTD formulation: ZTD = ZHD(Pobs) + ZWD
    !                 Pobs computed from P0 using CMC hydrostatic extrapolation.
    !
    implicit none

    ! Arguments:
    real(dp)           , intent(in)     :: hv    ! height of ZTD observation Zobs (m)
    type(gps_profilezd), intent(in)     :: prf   ! local model profile (type gps_profilezd)
    logical            , intent(in)     :: lBevis! true/false --> use Bevis instead of Rueger k values
    real(dp)           , intent(in)     :: dzmin
    type(gps_diff)     , intent(out)    :: ZTDopv! ZTD (m) at height of observation (with derivatives)
    real(dp)           , intent(out)    :: rPobs ! Pressure (Pa) at height of observation
    integer(i4)        , intent(in)     :: mode

    ! Locals:
    integer(i4)                         :: ngpslev
    integer(i4)                         :: levIndex, jloc
    real(dp)                            :: h, x, lat, sLat, dh
    real(dp)                            :: k1, k2, k3, k2p
    real(dp)                            :: zcon, zcon1, zconh, zfph, zconw
    type(gps_diff)                      :: dz, tvsfc, tobs, qobs, tvobs, naobs, Pobs
    type(gps_diff)                      :: dztddp, dztddpm
    type(gps_diff)                      :: zhd, tbar, qbar, qtterm, zsum, ztmobs, zqmobs
    type(gps_diff)                      :: zpbar, ztbar, zqbar, zrmean, zwd
    type(gps_diff)                      :: dzm, dzp
    real(dp)           , parameter      :: delta = 0.6077686814144_dp
    real(dp)           , parameter      :: eps   = 0.6219800221014_dp
    real(dp)           , parameter      :: kappa = (1.0_dp/eps)-1.0_dp
    real(dp)           , parameter      :: gamma = 0.0065_dp    ! -dT/dz (K/m)
    real(dp)           , parameter      :: Rgm = 9.784_dp
    real(dp)           , parameter      :: dzmax = 100.0
    ! Reuger (2002) refractivity constants (MKS units)
    real(dp)           , parameter      :: k1r = 0.776890_dp
    real(dp)           , parameter      :: k2r = 0.712952_dp
    real(dp)           , parameter      :: k3r = 3754.63_dp
    ! Bevis (1994) refractivity constants (MKS units)
    real(dp)           , parameter      :: k1b = 0.776000_dp
    real(dp)           , parameter      :: k2b = 0.704000_dp
    real(dp)           , parameter      :: k3b = 3739.000_dp

    ! Refractivity constants to use
    if ( lBevis ) then
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
    do levIndex = 1, ngpslev-1
      jloc = MPC_missingValue_INT
      if ((h <= prf%gst(levIndex)%Var) .and. (h > prf%gst(levIndex+1)%Var)) then
        jloc = levIndex   ! the model level above the observation
        exit
      end if
    end do

    if (h <= prf%gst(ngpslev)%Var) then  ! obs is at or below model lowest level
      jloc = ngpslev
    end if

    if ( mode == 2 ) then
      ! Compute ZTD the Vedel and Huang (2004) way: (as in old s/r gpsztdop.ftn)
      zcon  = 1.0e-06_dp*p_Rd
      zcon1 = zcon*k1
      zconw = zcon/eps
      zconh = zcon1/Rgm
      zfph = (1.0_dp - 2.66e-03_dp*cos(2.0*lat) - 2.8e-07_dp*h)

      ! Pressure at obs height (CMC hydrostatic extrapolation from Psfc)
      x      = ec_rg/(p_Rd*gamma)
      tvsfc  = prf%tst(ngpslev)*(1._dp+delta*prf%qst(ngpslev))
      Pobs   = prf%pst(ngpslev)*(((tvsfc-gamma*dh)/tvsfc)**x)
      ! Dry delay ZHD (m) at obs height
      zhd    = (zconh/zfph) * Pobs

      ! Integrate column q/T on pressure levels to get model ZWD
      do levIndex = 1, ngpslev-1
        tbar = (prf%tst(levIndex) + prf%tst(levIndex+1))*0.5_dp
        qbar = (prf%qst(levIndex) + prf%qst(levIndex+1))*0.5_dp
        qtterm = ((qbar + kappa*qbar**2 )/phf_gravityalt(sLat,prf%gst(levIndex)%Var))*(k2p + k3/tbar)
        if ( levIndex == 1 ) then
          zsum = qtterm*(prf%pst(levIndex+1)-prf%pst(levIndex))
        else
          zsum = zsum + qtterm*(prf%pst(levIndex+1)-prf%pst(levIndex))
        end if
      end do

      ! Compute ZWD at obs height using Higgins method (HU constant over dh layer)
      ztmobs = prf%tst(ngpslev) - (gamma * dh)
      zqmobs = prf%qst(ngpslev)
      zpbar  = (Pobs + prf%pst(ngpslev)) * 0.5_dp
      ztbar  = (ztmobs + prf%tst(ngpslev)) * 0.5_dp
      zqbar  = (zqmobs + prf%qst(ngpslev)) * 0.5_dp
      ! Mean (wet) refractivity of dz layer
      zrmean = 1.0e-06_dp*(k2p*((zpbar*zqbar)/(eps*ztbar)) + k3*((zpbar*zqbar)/(eps*ztbar**2)))

      ! Make sure adjusted ZWD >= 0
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
            dztddp = 1.e-6_dp * naobs * (p_Rd*tvobs)/phf_gravityalt(sLat, h)
            dztddpm = (dztddp + prf%rst(jloc))/2._dp  ! mean value of dZTD/dp over dh layer
          end if
          ZTDopv = prf%ztd(jloc) + dztddpm*(Pobs-prf%pst(jloc))
        end if
      end if
    end if
    rPobs = Pobs%Var
  end subroutine gps_ztdopv

  !--------------------------------------------------------------------------
  ! gps_pw
  !--------------------------------------------------------------------------
  pure subroutine gps_pw(prf, PW)
    !
    ! :Purpose: To compute lowest level PW (kg/m2) using layer mean Q and layer
    !           delta_p (Pa)
    !
    implicit none

    ! Arguments:
    type(gps_profilezd), intent(in)     :: prf
    real(dp)           , intent(out)    :: PW

    ! Locals:
    integer(i4)                         :: levIndex, ngpslev
    real(dp)                            :: qbar, gt, gb, g, lat, sLat
    real(dp)                            :: pt, pb

    ngpslev = prf%ngpslev
    lat     = prf%rLat
    sLat    = sin(lat)

    PW = 0.0_dp

    do levIndex = 1, ngpslev-1
      qbar = 0.5_dp * (prf%qst(levIndex+1)%Var + prf%qst(levIndex)%Var)
      gt  = phf_gravityalt(sLat, prf%gst(levIndex)%Var)
      gb  = phf_gravityalt(sLat, prf%gst(levIndex+1)%Var)
      pt  = prf%pst(levIndex)%Var
      pb  = prf%pst(levIndex+1)%Var
      g   = 0.5_dp * (gt + gb)
      PW  = PW + (qbar/g)*(pb-pt)
    end do
  end subroutine gps_pw

  !--------------------------------------------------------------------------
  ! gps_bndopv
  !--------------------------------------------------------------------------
  pure subroutine gps_bndopv(impv, azmv, nval, prf, bstv)
    ! :Purpose: GPSRO Bending angle operator
    !
    ! :Note: The operator is loosely based on Fjeldbo 1971,
    !        but adds some elements from Healy 2001, Burrows et al 2014 and Aparicio et al 2018
    implicit none

    ! Arguments:
    real(dp)           , intent(in)     :: impv(:)
    real(dp)           , intent(in)     :: azmv(:)
    integer(i4)        , intent(in)     :: nval
    type(gps_profile)  , intent(in)     :: prf
    type(gps_diff)     , intent(out)    :: bstv(:)

    ! Locals:
    integer(i4)                         :: levIndexObs, ngpslev, numLevels, levelHgh, levelLow, levIndexAnl
    type(gps_diff)                      :: nu(ngpssize), lnu(ngpssize), n(ngpssize), z(ngpssize), dlnudz(ngpssize)
    type(gps_diff)                      :: v(ngpssize), dtvdz(ngpssize), dzip, dsum, mvec(ngpssize), mvec1(ngpssize)
    type(gps_diff)                      :: Ntp, alpha_top, kat, zk, lrmd0, drmd, kvec1(ngpssize)
    real(dp)                            :: a, cazm, sazm, hi, bref, R
    integer(i4)                         :: imin, imax
    type(gps_diff)                      :: sum4, zimd, zipd, Dim, num, nm, k, m, zmin, zmax, zma, D0, D1, D2, k1, m1
    type(gps_diff)                      :: aI00, aI10, aI01, aI02
    logical                             :: lTooHigh, lReflection

    ! Model levels
    ngpslev=prf%ngpslev
    R = prf%Rad
    imin = -1
    imax = -1
    zmin = 1.e7_dp
    zmax =-1.e7_dp
    do levIndexAnl = 1, ngpslev
      nu (levIndexAnl) = prf%rst(levIndexAnl)
      lnu(levIndexAnl) = log(nu(levIndexAnl))
      n  (levIndexAnl) = 1._dp+nu(levIndexAnl)*1e-6_dp
      z  (levIndexAnl) = n(levIndexAnl)*((R+prf%geoid)+prf%gst(levIndexAnl))
      v  (levIndexAnl) = prf%vst(levIndexAnl)
      ! Record min and max z:
      if (z(levIndexAnl)%Var < zmin%Var) then
        zmin = z(levIndexAnl)
        imin = levIndexAnl
      end if
      if (z(levIndexAnl)%Var > zmax%Var) then
        zmax = z(levIndexAnl)
        imax = levIndexAnl
      end if
    end do
    ! Vertical gradients
    do levIndexAnl = 2, ngpslev
      dzip = z(levIndexAnl-1)-z(levIndexAnl)
      dlnudz(levIndexAnl) = (lnu(levIndexAnl-1)-lnu(levIndexAnl))/dzip
      dtvdz (levIndexAnl) = (v  (levIndexAnl-1)-v  (levIndexAnl))/dzip
      lrmd0= 0.5_dp*        (lnu(levIndexAnl-1)+lnu(levIndexAnl))
      drmd = prf%lrmd(levIndexAnl)-lrmd0
      mvec(levIndexAnl) = 4 * drmd / (dzip*dzip)
      k1 = (-9.8_dp)/(p_Rd*v(levIndexAnl))
      m1 = (-1)*k1*dtvdz(levIndexAnl)/(2*v(levIndexAnl))
      kvec1(levIndexAnl) = k1
      mvec1(levIndexAnl) = m1
    end do
    dlnudz(1) = dlnudz(2)
    dtvdz (1) = dtvdz (2)
    mvec  (1) = mvec  (2)
    ! Number of observed levels in the profile
    numLevels  = size(impv)
    if (nval < numLevels) numLevels=nval

    do levIndexObs =  1, numLevels
      ! Each observation is either:
      ! 1- Above the highest layer (lTooHigh, above zmax, nearly always level 1)
      ! 2- Between the highest and the lowest layer. We then search i such that z(i-1) <= a <= z(i).
      ! 3- Below the lowest layer (lReflection, below zmin, often at level ngpslev, but not necessarily
      a  = impv(levIndexObs)
      hi = a-R
      bref = 0.03_dp*exp(-hi/6500._dp)
      cazm = cos(azmv(levIndexObs))
      sazm = sin(azmv(levIndexObs))
      !find model levels that bracket the observation (1 is the highest level)
      levelHgh = 1
      levelLow = 2
      lTooHigh = .False.
      lReflection = .False.
      if (a >= zmax%Var) then
        ! Obs above lid
        lTooHigh = .True.
        levelHgh = imax
        levelLow = imax+1
      else if (a >= zmin%Var) then
        ! Normal case
        do levIndexAnl = 2, ngpslev
          levelHgh = levIndexAnl - 1
          levelLow = levIndexAnl
          if (z(levelLow)%Var <= a .and. a <= z(levelHgh)%Var) exit
        end do
      else
        ! Obs below lowest level, assume reflection
        lReflection = .True.
        levelHgh = imin-1
        levelLow = imin
      end if
      sum4 = 0._dp
      if (.not.lTooHigh) then
        do levIndexAnl = levelLow, 2, -1
          zimd = z(levIndexAnl)
          zipd = z(levIndexAnl-1)
          dzip = zipd-zimd
          Dim = nu(levIndexAnl)
          num = 1.e-6_dp*Dim
          nm = 1._dp+num
          k = dlnudz(levIndexAnl)
          m = mvec(levIndexAnl)
          D0 = num*(k-dzip*m)/nm
          D1 = num*(m*(dzip**2*m+2*nm)-2*dzip*k*m+k**2)/nm**2
          D2 = (num*(k-dzip*m)*(m*(num*(dzip**2*m-6) - dzip**2*m - 6) - 2*(num - 1)*dzip*k*m + (num - 1)*k**2))/(2*nm**3)
          if (levIndexAnl == levelLow .and. .not. lReflection) then
            aI00 = Integral_I00_a(a, zipd)
            aI10 = Integral_I10_a(a, zipd)
            aI01 = Integral_I01_a(a, zipd)
            aI02 = Integral_I02_a(a, zipd)
          else
            aI00 = Integral_I00_b(a, zimd, zipd)
            aI10 = Integral_I10_b(a, zimd, zipd)
            aI01 = Integral_I01_b(a, zimd, zipd)
            aI02 = Integral_I02_b(a, zimd, zipd)
          end if
          dsum = D0*(aI00 - 1._dp/(4*a) * aI10) + D1 * aI01 + D2 * aI02
          sum4 = sum4 + dsum
        end do
        sum4 = (-sqrt(2*a))*sum4
        if (lReflection) then
          sum4 = sum4 - acos(a/zmin)
        end if
        zma = zmax-a
      else
        zma = 0._dp
      end if
      !
      ! Top
      !
      kat = (-1._dp)*dlnudz(1)
      ! Prevent negative scale height. Should be approx 1 / 6500-8000m. Force limits between 1/10000m and 1/3000m.
      if (kat%Var < 1._dp/10000._dp) then
        kat = 1._dp/10000._dp
      end if
      if (kat%Var > 1._dp/3000._dp) then
        kat = 1._dp/3000._dp
      end if
      Ntp = nu(imax)
      if (.not.lTooHigh) then
        zk = kat*zma
        alpha_top = 1.e-6_dp*sqrt(2*MPC_PI_R8*a*kat)*Ntp*exp(zk)*(1._dp-erf(sqrt(zk)))
      else
        alpha_top = 1.e-6_dp*sqrt(2*MPC_PI_R8*a*kat)*Ntp
      end if
      sum4 = sum4 + alpha_top
      bstv(levIndexObs) = sum4
    end do
  end subroutine gps_bndopv

  !--------------------------------------------------------------------------
  ! Integral_I00_a
  !--------------------------------------------------------------------------
  pure function Integral_I00_a(a, z2)
    ! :Purpose: Integral for the bending evaluation
    implicit none
    real(dp)           , intent(in)     :: a
    type(gps_diff)     , intent(in)     :: z2

    type(gps_diff)                      :: Integral_I00_a

    Integral_I00_a = 2 * sqrt(z2-a)
  end function Integral_I00_a

  !--------------------------------------------------------------------------
  ! Integral_I00_b
  !--------------------------------------------------------------------------
  pure function Integral_I00_b(a, z1, z2)
    ! :Purpose: Integral for the bending evaluation
    implicit none
    real(dp)           , intent(in)     :: a
    type(gps_diff)     , intent(in)     :: z1, z2

    type(gps_diff)                      :: Integral_I00_b

    Integral_I00_b = 2 * (sqrt(z2-a) - sqrt(z1-a))
  end function Integral_I00_b

  !--------------------------------------------------------------------------
  ! Integral_I10_a
  !--------------------------------------------------------------------------
  pure function Integral_I10_a(a, z2)
    ! :Purpose: Integral for the bending evaluation
    implicit none
    real(dp)           , intent(in)     :: a
    type(gps_diff)     , intent(in)     :: z2

    type(gps_diff)                      :: Integral_I10_a

    Integral_I10_a = (2._dp/3._dp) * (z2-a)**1.5_dp
  end function Integral_I10_a

  !--------------------------------------------------------------------------
  ! Integral_I10_b
  !--------------------------------------------------------------------------
  pure function Integral_I10_b(a, z1, z2)
    ! :Purpose: Integral for the bending evaluation
    implicit none
    real(dp)           , intent(in)     :: a
    type(gps_diff)     , intent(in)     :: z1, z2

    type(gps_diff)                      :: Integral_I10_b

    Integral_I10_b = (2._dp/3._dp) * ((z2-a)**1.5_dp - (z1-a)**1.5_dp)
  end function Integral_I10_b

  !--------------------------------------------------------------------------
  ! Integral_I01_a
  !--------------------------------------------------------------------------
  pure function Integral_I01_a(a, z2)
    implicit none
    real(dp)           , intent(in)     :: a
    type(gps_diff)     , intent(in)     :: z2

    type(gps_diff)                      :: Integral_I01_a

    Integral_I01_a = (2._dp/3._dp) * (z2-a)**1.5_dp
  end function Integral_I01_a

  !--------------------------------------------------------------------------
  ! Integral_I01_b
  !--------------------------------------------------------------------------
  pure function Integral_I01_b(a, z1, z2)
    ! :Purpose: Integral for the bending evaluation
    implicit none
    real(dp)           , intent(in)     :: a
    type(gps_diff)     , intent(in)     :: z1, z2

    type(gps_diff)                      :: Integral_I01_b

    Integral_I01_b = (2._dp/3._dp) * ( sqrt(z2-a) * (2*a+z2-3*z1) - sqrt(z1-a) * (2*a-2*z1) )
  end function Integral_I01_b

  !--------------------------------------------------------------------------
  ! Integral_I02_a
  !--------------------------------------------------------------------------
  pure function Integral_I02_a(a, z2)
    ! :Purpose: Integral for the bending evaluation
    implicit none
    real(dp)           , intent(in)     :: a
    type(gps_diff)     , intent(in)     :: z2

    type(gps_diff)                      :: Integral_I02_a

    Integral_I02_a = (2._dp/5._dp) * (z2-a)**2.5_dp
  end function Integral_I02_a

  !--------------------------------------------------------------------------
  ! Integral_I02_b
  !--------------------------------------------------------------------------
  pure function Integral_I02_b(a, z1, z2)
    ! :Purpose: Integral for the bending evaluation
    implicit none
    real(dp)           , intent(in)     :: a
    type(gps_diff)     , intent(in)     :: z1, z2

    type(gps_diff)                      :: Integral_I02_b

    Integral_I02_b = (2._dp/15._dp) * ( sqrt(z2-a) * (8*a*a+4*a*(z2-5*z1)+3*z2*z2-10*z1*z2+15*z1*z1) - sqrt(z1-a) * (8*a*a-16*a*z1+8*z1*z1) )
  end function Integral_I02_b

  !--------------------------------------------------------------------------
  ! gps_setupro
  !--------------------------------------------------------------------------
  subroutine gps_setupro
    ! :Purpose: Initialisation of Radio Occultations - to read and to initialize
    !           GPSRO namelist parameters and print information on options
    !           selected.
    implicit none

    ! Locals:
    integer(i4) :: ierr, satID

    ! Namelist variables for GPS-RO
    integer(i4) :: levelGPSRO       ! Data level to use (1 for bending angle, 2 for refractivity)
    integer(i4) :: gpsro_maxPrfSize ! Maximal number of data that is expected from a profile (default 300)
    real(dp)    :: surfMin          ! Minimum allowed distance to the model surface (default 0 m)
    real(dp)    :: hSfMin           ! Minimum allowed MSL height of an obs          (default 0 m)
    real(dp)    :: hTpMax           ! Maximum allowed MSL height of an obs          (default 70000 m)
    real(dp)    :: bgckBand         ! Maximum allowed deviation abs(O-P)/P          (default 0.05)
    real(dp)    :: hTpMaxEr         ! Maximum MSL height to evaluate the obs error  (default to hTpMax)
    real(sp)    :: wGPS(0:1023, 4)  ! wGPS values for each satellite sensor
    character(20) :: gpsroError     ! key for using dynamic/static refractivity error estimation (default 'DYNAMIC')
    logical     :: gpsroBNorm       ! Normalize based on B=H(x) (default=.True.), otherwise, exponential reference
    logical     :: gpsroEotvos      ! Add an operator-only Eotvos corr to local gravity (shift of alt, default False)
    real(dp)    :: gpsroNsigma      ! Factor applied to obs err for bckg check when gpsroBNorm is .true. (default 1.d6)
    logical     :: gpsroCurvAnisot  ! Apply vert shift to account for asphericity of ellipsoid wrt ref tangent sphere
    logical     :: gpsroNCurv       ! Add small curvature in the log-linear vertical interpolation of N
    integer(i4) :: gpsroNFlavour    ! Choice of refractivity constants: 0 for AL11, 1 for the 2015-25 update.
    integer(i4) :: dateStamp, hour, day, month, ndays, yyyy, jd
    NAMELIST /NAMGPSRO/ levelGPSRO, gpsro_maxPrfSize, surfMin, hSfMin, hTpMax, hTpMaxER, &
        bgckBand, wGPS, gpsroError, gpsroBNorm, gpsroEotvos, gpsroNsigma, gpsroCurvAnisot, gpsroNCurv, gpsroNFlavour

    !
    !   Define default values:
    !
    levelGPSRO        = gps_Level_RO_Ref
    gpsro_maxPrfSize  = 300
    surfMin           = 0._dp
    hSfMin            = 0._dp
    hTpMax            = 70000._dp
    hTpMaxEr          = -1._dp
    bgckBand          = 0.05_dp
    gpsroError        = 'DYNAMIC'
    gpsroBNorm        = .True.
    gpsroEotvos       = .False.
    gpsroNsigma       = 1000000._dp
    gpsroCurvAnisot   = .False.
    gpsroNCurv        = .False.
    gpsroNFlavour     = 0
    !
    !   Force a pre-NML default for the effective data weight of all
    !   GPSRO satellites. This array has rows 0-1023 (following BUFR element
    !   SATID), and 4 cols. The 4 parameters for each SATID are used to
    !   represent data correlation, a combined property of the satellite
    !   hardware and provider postprocessing.
    !   The default assumes no correlation.
    !
    wGPS              = 0._sp
    wGPS(:,1)         = 1._sp
    !
    !   Override with NML values:
    !
    call utl_tmg_start(181, 'low-level--readNML')
    read(utl_flnml, nml=NAMGPSRO, iostat=ierr)
    if (ierr /= 0) call utl_abort('gps_setupro: Error reading namelist')
    call utl_tmg_stop(181)
    !
    if (hTpMaxEr < 0._dp) hTpMaxEr = hTpMax
    gps_level_RO      = levelGPSRO
    gps_RO_maxPrfSize = gpsro_maxPrfSize
    gps_surfMin       = surfMin
    gps_hSfMin        = hSfMin
    gps_hTpMax        = hTpMax
    gps_hTpMaxEr      = hTpMaxEr
    gps_bgckBand      = bgckBand
    gps_roError       = gpsroError
    gps_roBNorm       = gpsroBNorm
    gps_wGPS          = wGPS
    gps_roEotvos      = gpsroEotvos
    gps_roNsigma      = gpsroNsigma
    gps_roCurvAnisot  = gpsroCurvAnisot
    gps_roNCurv       = gpsroNCurv
    gps_roNFlavour    = gpsroNFlavour

    dateStamp = tim_getDatestamp()
    call tim_dateStampToYYYYMMDDHH(dateStamp, hour, day, month, ndays, yyyy, .False.)

    ! Determine Julian day number from date (RMNLIB function)
    call jDatec(jd, yyyy, month, day)

    ! Continuous Julian date, in Julian years (offset to read 2000 at Jan 1st 2000):
    gps_roYear = 2000._dp+(jd-2451545._dp)/365.25_dp

    ! Safety limits:
    if (gps_roYear < 1990._dp) gps_roYear = 1990._dp
    if (gps_roYear > 2100._dp) gps_roYear = 2100._dp

    if (mmpi_myid == 0) then
      write(*,*)'NAMGPSRO',gps_level_RO, gps_RO_maxPrfSize, gps_surfMin, gps_hSfMin, &
           gps_hTpMax, gps_hTpMaxEr, gps_bgckBand, trim(gps_roError), gps_roBNorm, gps_roEotvos, &
           gps_roNsigma, gps_roCurvAnisot, gps_roNCurv, gps_roNFlavour
      do satID = 0, 1023
        if ( .not. utl_isEqual(wGPS(satID, 2), 0.) ) then
          write(*,*)'wGPS', satID, gps_wGPS(satID, 1:4)
        end if
      end do
      write(*,*)'gps_setupRO: Epoch:', yyyy, ndays, month, day, hour, gps_roYear
    end if
  end subroutine gps_setupro

  !--------------------------------------------------------------------------
  ! gps_setNumROProfiles
  !--------------------------------------------------------------------------
  subroutine gps_setNumROProfiles(numROProfiles)
    ! :Purpose: Create RO buffer space for numROProfiles locations
    implicit none

    ! Arguments:
    integer(i4)        , intent(in)     :: numROProfiles

    gps_numROProfiles = numROProfiles
    if (gps_numROProfiles > 0) then
      if (.not.allocated(gps_vRO_IndexPrf)) allocate(gps_vRO_IndexPrf(gps_numROProfiles, 10))
      if (.not.allocated(gps_vRO_dR)      ) allocate(gps_vRO_dR(      gps_numROProfiles, gps_RO_maxPrfSize))
    end if
  end subroutine gps_setNumROProfiles

  !--------------------------------------------------------------------------
  ! gps_setROIndexPrf
  !--------------------------------------------------------------------------
  subroutine gps_setROIndexPrf(profileIndex, headerIndex, varNum, iSat, iDsc, dR)
    ! :Purpose: Fill RO buffer space of one profile
    implicit none

    ! Arguments:
    integer(i4)        , intent(in)     :: profileIndex    ! index within RO Profiles
    integer(i4)        , intent(in)     :: headerIndex     ! index within all observations
    integer(i4)        , intent(in)     :: varNum          ! observation variable kind (15036, ref; or 15037, bnd)
    integer(i4)        , intent(in)     :: iSat            ! satellite ID (receiver)
    integer(i4)        , intent(in)     :: iDsc            ! descendin/ascending profile
    real(dp)           , intent(in)     :: dR(:)           ! individual vertical shift associated to asphericity

    gps_vRO_IndexPrf(profileIndex, 1) = headerIndex
    gps_vRO_IndexPrf(profileIndex, 2) = varNum
    gps_vRO_IndexPrf(profileIndex, 3) = iSat
    gps_vRO_IndexPrf(profileIndex, 4) = iDsc
    gps_vRO_dR(      profileIndex, :) = dR(:)
  end subroutine gps_setROIndexPrf

  !--------------------------------------------------------------------------
  ! gps_iprofile_from_index
  !--------------------------------------------------------------------------
  integer function gps_iprofile_from_index(headerIndex)
    ! :Purpose: Find an obs profile stored in RO buffer space
    !           They are uniquely identified with the CMA index.
    implicit none

    ! Arguments:
    integer(i4)        , intent(in)     :: headerIndex

    ! Locals:
    integer(i4)                         :: profileIndex

    gps_iprofile_from_index = -1
    do profileIndex = 1, gps_numROProfiles
      if (headerIndex == gps_vRO_IndexPrf(profileIndex, 1)) then
        gps_iprofile_from_index = profileIndex
        return
      end if
    end do
    return
  end function gps_iprofile_from_index

  !--------------------------------------------------------------------------
  ! gps_setupgb
  !--------------------------------------------------------------------------
  subroutine gps_setupgb
    !
    ! :Purpose: Initialisation of ground-based GPS - to read and to initialize
    !           GB-GPS namelist parameters and print information on options
    !           selected.
    !
    implicit none

    ! Locals:
    integer(i4)                         :: ierr

    ! Namelist variables for Ground-based GPS (ZTD)
    real(dp)    :: dzMin     ! Minimum DZ = Zobs-Zmod (m) for which DZ adjustment to ZTD will be made
    real(dp)    :: dzMax     ! Maximum DZ (m) over which ZTD rejected due to topography (when LTOPOFILT = .TRUE.)
    real(dp)    :: yZTDErr   ! If < 0 err from inp files; if > 0 use value as const err (m); if =0 error=f(ZWD)
    real(dp)    :: ySFErrWgt ! Scale factor for GPS surf met err (account for time series obs with err corr)
    real(dp)    :: yZDErrWgt ! Scale factor for GPS ZTD errors (account for time series obs with err corr)
    logical     :: lAssMet   ! Choose to assimilate GPS Met surface P, T, T-Td
    logical     :: lLBLMet   ! Indicate that surf met data blacklisted for GPS sites close to surface weather stations
    logical     :: lBevis    ! If .true. use Bevis(1994); if .false. use Rueger(2002) refractivity (k1,k2,k3) constants
    logical     :: l1Obs     ! Choose to select a single ZTD observation
    logical     :: lTestOp   ! Choose to test ZTD observation operator (Omp and Bgck modes only)
    integer(i4) :: iRefOpt   ! 1 = conventional expression for N using k1,k2,k3; 2 = Aparicio&Laroche (incl. compress)
    integer(i4) :: iZTDOp    ! 1 = use stored ZTD profiles to get ZTDmod; 2 = Vedel & Huang ZTD: ZTDmod=ZHD(Pobs)+ZWD

    NAMELIST /NAMGPSGB/ dzMin, dzMax, yZTDErr, lAssMet, ySFErrWgt,  &
         lLBLMet, yZDErrWgt, lBevis, l1Obs, lTestOp, iRefOpt, iZTDOp

    !*  .  1.1 Default values
    !!  .      --------------

    dzMin             = 2._dp
    dzMax             = 1000._dp
    yZTDErr           = 0.012_dp
    lAssMet           = .TRUE.
    ySFErrWgt         = 1._dp
    lLBLMet           = .FALSE.
    yZDErrWgt         = 1._dp
    lBevis            = .TRUE.
    iRefOpt           = 1
    l1Obs             = .FALSE.
    lTestOp           = .FALSE.
    iZTDOp            = 1

    call utl_tmg_start(181, 'low-level--readNML')
    read(utl_flnml, nml=NAMGPSGB, iostat=ierr)
    if (ierr /= 0) call utl_abort('gps_setupgb: Error reading namelist')
    gps_gb_dzMin      = dzMin
    gps_gb_dzMax      = dzMax
    gps_gb_yZTDErr    = yZTDErr
    gps_gb_lAssMet    = lAssMet
    gps_gb_ySFErrWgt  = ySFErrWgt
    gps_gb_lLBLMet    = lLBLMet
    gps_gb_yZDErrWgt  = yZDErrWgt
    gps_gb_lBevis     = lBevis
    gps_gb_iRefOpt    = iRefOpt
    gps_gb_l1Obs      = l1Obs
    gps_gb_lTestOp    = lTestOp
    gps_gb_iZTDOp     = iZTDOp
    if (mmpi_myid == 0) write(*, nml=NAMGPSGB)
    call utl_tmg_stop(181)

    if (l1Obs .and. mmpi_myid == 0) then
      write(*,*)    ' ******************************************'
      write(*,*)    ' * GB-GPS: ONE OBSERVATION MODE           *'
      write(*,*)    ' ******************************************'
    end if

    !   Options to fix/adjust model ZTD to observation height and
    !   assimilate GPS met data

    if (mmpi_myid == 0) then
      write(*,*)    ' ******************************************'
      write(*,*)    ' *        GB-GPS OBSERVATIONS             *'
      write(*,*)    ' * DZ ADJUSTMENT IN gps_ztdopv if dz>dzMin*'
      write(*,*)    ' * ZTD NOT ASSIM. if dz > dzMax           *'
      write(*,*)    ' ******************************************'
      write(*,*)    'dzMin, dzMax = ', dzMin, dzMax
      if (lAssMet) then
        if ( lLBLMet ) then
          write(*,*)' *****************************************'
          write(*,*)' *          GB-GPS OBSERVATIONS          *'
          write(*,*)' *     GPS MET DATA ARE ASSIMILATED      *'
          write(*,*)' *     BUT BLACKLISTED NEAR SYNO STNS    *'
          write(*,*)' *****************************************'
          write(*,*)'ySFErrWgt = ', ySFErrWgt
          write(*,*)'yZDErrWgt = ', yZDErrWgt
        else
          write(*,*)' *****************************************'
          write(*,*)' *          GB-GPS OBSERVATIONS          *'
          write(*,*)' *     GPS MET DATA ARE ASSIMILATED      *'
          write(*,*)' *****************************************'
          write(*,*)'ySFErrWgt = ', ySFErrWgt
          write(*,*)'yZDErrWgt = ', yZDErrWgt
        end if
      else
        write(*,*)  ' *****************************************'
        write(*,*)  ' *          GB-GPS OBSERVATIONS          *'
        write(*,*)  ' *   GPS MET DATA ARE NOT ASSIMILATED    *'
        write(*,*)  ' *****************************************'
        write(*,*)  'yZDErrWgt = ', yZDErrWgt
      end if

      if (yZTDErr < 0._dp) then
        write(*,*)  ' *****************************************'
        write(*,*)  ' *          GB-GPS OBSERVATIONS          *'
        write(*,*)  ' *    ZTD OBSERVATION ERROR FROM FERR    *'
        write(*,*)  ' *****************************************'
      else if (yZTDErr > 0._dp) then
        write(*,*)  ' *****************************************'
        write(*,*)  ' *          GB-GPS OBSERVATIONS          *'
        write(*,*)  ' *     ZTD OBSERVATION ERROR IS FIXED    *'
        write(*,*)  ' *****************************************'
        write(*,*)  'yZTDErr (mm) = ', yZTDErr*1000.D0
      else
        write(*,*)  ' *****************************************'
        write(*,*)  ' *          GB-GPS OBSERVATIONS          *'
        write(*,*)  ' *   ZTD OBSERVATION ERROR IS FROM ZWD   *'
        write(*,*)  ' *   USING SD(O-P) STATS (REGRESSION)    *'
        write(*,*)  ' *****************************************'
      end if

      if (iRefOpt == 1) then
        if (lBevis) then
          write(*,*)' *****************************************'
          write(*,*)' *          GB-GPS OBSERVATIONS          *'
          write(*,*)' *  CONVENTIONAL REFACTIVITY N USING     *'
          write(*,*)' *  BEVIS 92 K1, K2, K3 TO COMPUTE ZTD   *'
          write(*,*)' *****************************************'
        else
          write(*,*)' *****************************************'
          write(*,*)' *          GB-GPS OBSERVATIONS          *'
          write(*,*)' *  CONVENTIONAL REFACTIVITY N USING     *'
          write(*,*)' *  RUEGER 02 K1, K2, K3 TO COMPUTE ZTD  *'
          write(*,*)' *****************************************'
        end if

        if (iZTDOp == 1) then
          write(*,*)' *****************************************'
          write(*,*)' *          GB-GPS OBSERVATIONS          *'
          write(*,*)' *   NORMAL ZTD OPERATOR -- ZTD COMPUTED *'
          write(*,*)' *           FROM ZTD(K) PROFILE         *'
          write(*,*)' *****************************************'
        else
          write(*,*)' *****************************************'
          write(*,*)' *          GB-GPS OBSERVATIONS          *'
          write(*,*)' *   ORIGINAL OPERATOR -- ZTD = ZHD+ZWD  *'
          write(*,*)' *        VEDEL AND HUANG (2004)         *'
          write(*,*)' *****************************************'
        end if
      else
        write(*,*)  ' *****************************************'
        write(*,*)  ' *          GB-GPS OBSERVATIONS          *'
        write(*,*)  ' *  APARICIO & LAROCHE REFRACTIVITY N    *'
        write(*,*)  ' *         USED TO COMPUTE ZTD           *'
        write(*,*)  ' *****************************************'
      end if
    end if
  end subroutine gps_setupgb

  !--------------------------------------------------------------------------
  ! gps_setNumZTD
  !--------------------------------------------------------------------------
  subroutine gps_setNumZTD(numZTD)
    ! :Purpose: Create ZTD buffer space for numZTD places
    implicit none

    ! Arguments:
    integer(i4)        , intent(in)     :: numZTD

    gps_gb_numZTD = numZTD

    if (allocated(gps_ZTD_Index)) deallocate(gps_ZTD_Index)
    if (gps_gb_numztd > 0) then
      allocate(gps_ZTD_Index(gps_gb_numztd))
    end if
  end subroutine gps_setNumZTD

  !--------------------------------------------------------------------------
  ! gps_setZTDIndex
  !--------------------------------------------------------------------------
  subroutine gps_setZTDIndex(iztd, headerIndex)
    ! :Purpose: Fill ZTD buffer space of an obs
    implicit none

    ! Arguments:
    integer(i4)        , intent(in)     :: iztd
    integer(i4)        , intent(in)     :: headerIndex

    gps_ZTD_Index(iztd) = headerIndex
  end subroutine gps_setZTDIndex

  !--------------------------------------------------------------------------
  ! gps_iztd_from_index
  !--------------------------------------------------------------------------
  integer function gps_iztd_from_index(headerIndex)
    ! :Purpose: Find a GPSGB obs in stored ZTD buffer space.
    !           They are uniquely identified with the CMA header index.
    implicit none

    ! Arguments:
    integer(i4)        , intent(in)     :: headerIndex

    ! Locals:
    integer(i4)                         :: iObs

    gps_iztd_from_index = -1
    do iObs = 1, size(gps_ZTD_Index)
      if (headerIndex == gps_ZTD_Index(iObs)) then
        gps_iztd_from_index = iObs
        return
      end if
    end do
    return
  end function gps_iztd_from_index

end module gps_mod

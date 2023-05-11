
module regions_mod
  ! MODULE regions_mod (prefix='reg' category='7. Low-level data objects')
  !
  ! :Purpose: This is a subset of the EnKF module that deals with 
  !           operations involving the generation of a set of regions. 
  !           Each region consists of a small number of observations 
  !           which is to be assimilated simultaneously in the 
  !           sequential algorithm. 
  !
  use mathPhysConstants_mod
  use midasMpi_mod
  implicit none

  private
  public :: struct_reg
  public :: reg_getlatitude, reg_getblock, reg_locatestn, reg_init_struct

  ! REG_MXOBS : maximum number of observations in a region
  ! REG_MXOBS must be at least as large as the number of observations
  ! in an individual radiosonde (121 currently)
  integer, parameter :: REG_MXOBS = 900

  type :: struct_reg
   ! horizontal radius in km of the area in which observations are selected
   real(8) :: r0_km
   ! horizontal radius in km of the area in which an observation can have an impact
   real(8) :: r1_km
   ! the vertical extend in log(P) of the area in which an observation can have an 
   ! impact
   real(8) :: rz_logp 
   ! as r0km but in radians and adjusted considering that we want an integer
   ! number of latitude bands
   real(8) :: r0_rad 
   real(8) :: r1_rad ! equivalent of r1 in radians
   ! the horizontal length scale for the hadamard product (0.5 r1)
   real(8) :: scale_hor
   ! the vertical length scale for the hadamard product (in log (p))
   real(8) :: scale_ver
   ! number of latitude bands
   integer :: nlatband
   ! estimated maximum number of regions that could possibly be formed
   integer :: mxreg
  end type struct_reg
      
contains

  subroutine reg_getblock(nlatband, r0, latmin, latmax, nlonblock)
    !
    ! :Purpose: We have nlatband latitude bands (bounded by latitudes
    !           latmin and latmax). For each band determine how many 
    !           longitidunal blocks of at most r0*root(2) radians it 
    !           contains. Output is in nlonblock.
    !
    implicit none

    integer :: i, nlatband
    real(8) :: bar, dpie, lat, latmin(nlatband), latmax(nlatband), r0
    integer :: nlonblock(nlatband)

    dpie=MPC_PI_R8
      
    nlonblock(1)=1
    nlonblock(nlatband)=1
      
    do i = 2, nlatband-1 
     if (latmin(i)*latmax(i).gt.0) then
      lat = min(abs(latmin(i)),abs(latmax(i)))
     else
       ! at the equator
       lat = 0
     end if
     bar = (2**0.5)*r0/cos(lat)
     nlonblock(i) = int((2.0D0*dpie)/bar)+1
   end do 

  end subroutine reg_getblock


  subroutine reg_getlatitude(r0, nlatband, latmin, latcenter, latmax)
    !
    ! :Purpose: the circle is covered with nlatband latitude bands. 
    !           the polar caps (with radius r0 radians) form the first and 
    !           last band. Intermediate bands are of width r0*(2**0.5)
    !           For reach band we have the southern most latitude latmin,
    !           the central latitude latcenter (at the poles for the extreme 
    !           bands) and the northern most latitude latcenter
    !      
    implicit none
    integer :: i, nlatband
    real(8) :: dpie, r0, latmin(nlatband), latcenter(nlatband)
    real(8) :: latmax(nlatband)

    dpie = MPC_PI_R8

    latmin(1) = -0.5D0*dpie
    latcenter(1) = -0.5D0*dpie
    latmax(1) = -0.5D0*dpie+r0
 
    if (mmpi_myid == 0) write(*,*) 'get latitude r0 and nlatband: ', r0, nlatband
    if (mmpi_myid == 0) write(*,*) 'at ', 1, latmin(1), latcenter(1), latmax(1)

    latmin(nlatband) = 0.5D0*dpie-r0
    latcenter(nlatband) = 0.5D0*dpie
    latmax(nlatband) = 0.5D0*dpie

    if (mmpi_myid == 0) write(*,*) 'at ', nlatband, latmin(nlatband), latcenter(nlatband), latmax(nlatband)
    do i = 2, nlatband-1
      latmin(i) = latmax(i-1)
      latmax(i) = latmin(i)+r0*(2**0.5)
      latcenter(i) = latmin(i)+0.5D0*r0*(2**0.5)
      if (mmpi_myid == 0) write(*,*) 'at ', i, latmin(i), latcenter(i), latmax(i)
    end do
 
  end subroutine reg_getlatitude


  subroutine getr0(r0km, r0, nlatband)
    !
    ! :Purpose: The user specifies a radius of r0km (in km) within which 
    !           stations can be taken together for simultaneous analysis
    !           in a single batch. For simplicity of the search procedures
    !           we - instead - use somewhat smaller squares with sides of 
    !           at most r0km*root(2). We do have circles at the two poles.
    !           Here we compute the radius r0 in radians such that we arrive 
    !           exactly at nlatband latitude bands. The equator separates 
    !           two latitude bands.
    !         
    implicit none
    real(8), intent(in)  :: r0km
    real(8), intent(out) :: r0
    integer, intent(out) :: nlatband
    real(8)              :: dpie, eps, difference, r0new, dminpole, nbandr
    integer              :: nbandi

    dpie = MPC_PI_R8
      
    ! figure out if r0km*root(2) is a divisor of 10000-r0km

    eps = 0.01
    ! we put a circle with radius r0km around each pole.
    dminpole = 10000.-r0km
    ! potential fractional number of latitude bands.
    nbandr = dminpole/(r0km*(2**0.5))
    nbandi = nint(nbandr)
    difference = abs(dble(nbandi)-nbandr)
    if (difference.gt.eps) then
      ! increasing number of bands to nbandi
      nbandi = nbandi+1
    end if
    ! decrease r0km to r0new
    r0new = 10000./(dble(nbandi)*2**0.5 + 1.)
    ! compute the corresponding r0 in radians       
    r0 = (r0new/20000.)*dpie
    ! two hemispheres with a pole
    nlatband = 2*(nbandi+1)
      
  end subroutine getr0

  subroutine reg_init_struct(lsc, r0_km, r1_km, rz_logp)
    !
    ! :Purpose: store and derive parameters related to the localization
    !           in a structure.
    !
    implicit none
    type(struct_reg), intent(out) :: lsc
    real(8), intent(in) :: r0_km
    real(8), intent(in) :: r1_km
    real(8), intent(in) :: rz_logp
    real(8)  :: dpie
    integer :: ione
  
    lsc%r0_km   = r0_km
    lsc%r1_km   = r1_km
    lsc%rz_logp = rz_logp

    call getr0(lsc%r0_km, lsc%r0_rad, lsc%nlatband)
    lsc%r1_rad = lsc%r1_km/6370.
    lsc%scale_hor = 0.5*lsc%r1_rad
    lsc%scale_ver = 0.5*lsc%rz_logp

    dpie = MPC_PI_R8
    lsc%mxreg = floor(dpie/(lsc%r1_rad+(2**0.5)*lsc%r0_rad)**2)

    ione = 1
    lsc%mxreg = max(lsc%mxreg, ione)

  end subroutine reg_init_struct


  subroutine reg_locatestn(r0, lat, lon, nlatband, nlonblock, &
                           nblockoffset, iblock)
    !
    ! :Purpose: locate the lat-lon block in which position (lat,lon) falls
    !
    ! :Arguments:
    !       r0:        radius of a region (in radians)
    !       (lat,lon): latitude and longitude in radians
    !       nlatband:  number of latitude bands
    !       nlonblock: number of longitudinal blocks at each latitude band
    !       nblockoffset: offset (the number of blocks south of each
    !                  latitude band)
    !       iblock:   block number
    !
    implicit none
    integer :: iblock, ilat, ilon, nlatband
    integer :: nlonblock(nlatband), nblockoffset(nlatband)
    real(4) :: lat, lon
    real(8) :: dpie, r0

    dpie = MPC_PI_R8

    if (lat.le.(-0.5D0*dpie+r0)) then
      ilat = 1
    else if (lat.gt.(0.5D0*dpie-r0)) then
      ilat = nlatband
    else
      lat = lat+0.5D0*dpie-r0
      ilat = 1+ceiling(lat/(r0*(2**0.5)))
    end if
    ilat = min(max(ilat,1),nlatband)
    ilon = ceiling((lon/(2.0D0*dpie))*nlonblock(ilat))
    ! added for points at the Greenwich meridian
    ilon = min(max(ilon,1),nlonblock(ilat))
    iblock = nblockoffset(ilat)+ilon

  end subroutine reg_locatestn
      
end module regions_mod

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

!--------------------------------------------------------------------------
!! MODULE WindRotation (prefix="uvr")
!!
!! *Purpose*: To transform winds FROM the rotated spherical 
!!            coordinate system TO the non-rotated spherical coordinate system.
!!
!--------------------------------------------------------------------------
module WindRotation_mod
  use MathPhysConstants_mod, only: MPC_RADIANS_PER_DEGREE_R8,MPC_DEGREES_PER_RADIAN_R8, MPC_PI_R8
  use HorizontalCoord_mod
  implicit none
  save
  private

  ! Public Subroutines
  public :: uvr_Setup, uvr_RotateWind, uvr_RotateLatLon
  public :: uvr_RotateWindAdj

  integer, parameter :: msize = 3
  real(8)            :: grd_rot_8   (msize,msize) ! Real Earth to Rotated Earth
  real(8)            :: grd_rotinv_8(msize,msize) ! Rotated Earth to Real Earth

  logical :: initialized = .false.

  contains

!--------------------------------------------------------------------------
! uvr_Setup
!--------------------------------------------------------------------------
  subroutine uvr_Setup(hco_in)
    implicit none

    type(struct_hco), intent(in) :: hco_in

    !
    !-  Compute the rotation matrices (grd_rot_8 and grd_rotinv_8)
    !
    write(*,*)
    write(*,*) 'uvr_Setup: Starting...'

    print*,hco_in % xlon1
    print*,hco_in % xlat1
    print*,hco_in % xlon2
    print*,hco_in % xlat2

    call sugrdpar(hco_in % xlon1, hco_in % xlat1, hco_in % xlon2, hco_in % xlat2 )

    initialized = .true.

    write(*,*)
    write(*,*) 'uvr_Setup: Done!'

  end subroutine uvr_Setup

!--------------------------------------------------------------------------
! SUGRDPAR
!--------------------------------------------------------------------------
  SUBROUTINE sugrdpar(grd_xlon1, grd_xlat1, grd_xlon2, grd_xlat2 )
    !
    !- Compute the rotation matrix (r_8) that allows transformation
    !  from the non-rotated to the rotated spherical coordinate system.
    IMPLICIT NONE

    real(8), intent(in) :: grd_xlon1, grd_xlat1, grd_xlon2, grd_xlat2 

    integer :: ierr, ji, jj, j1, j2
    integer :: ii0, ij0, Idum, Imargin

    real(8) :: zxlon1_8,zxlat1_8,zxlon2_8,zxlat2_8

    real(8) :: a_8, b_8, c_8, d_8, xyz1_8(msize), xyz2_8(msize)
    real(8) :: zrot_t(msize,msize), zunit(msize,msize)

    zxlon1_8 = grd_xlon1
    zxlat1_8 = grd_xlat1
    zxlon2_8 = grd_xlon2
    zxlat2_8 = grd_xlat2

    if ( zxlon1_8 == zxlon2_8 .and. zxlat1_8 == zxlat2_8 ) then
      !
      !- 1. Non rotated grid (for test case only)
      !
      write(*,*)
      write(*,*) 'uvr_sugrdpar: Warning: This grid is not rotated !!!'
      grd_rot_8(1,1) = 1.0d0       
      grd_rot_8(1,2) = 0.0d0
      grd_rot_8(1,3) = 0.0d0
      grd_rot_8(2,1) = 0.0d0
      grd_rot_8(2,2) = 1.0d0
      grd_rot_8(2,3) = 0.0d0
      grd_rot_8(3,1) = 0.0d0
      grd_rot_8(3,2) = 0.0d0
      grd_rot_8(3,3) = 1.0d0

    else
      !
      !- 2. Rotated grid
      !
      call vllacar( xyz1_8, zxlon1_8, zxlat1_8 )
      call vllacar( xyz2_8, zxlon2_8, zxlat2_8 )

      !- 2.1 Compute a = cos(alpha) & b = sin(alpha)
      a_8 = (xyz1_8(1)*xyz2_8(1)) + (xyz1_8(2)*xyz2_8(2)) &
                                  + (xyz1_8(3)*xyz2_8(3))

      b_8 = sqrt (((xyz1_8(2)*xyz2_8(3)) - (xyz2_8(2)*xyz1_8(3)))**2 &
               +  ((xyz2_8(1)*xyz1_8(3)) - (xyz1_8(1)*xyz2_8(3)))**2 &
               +  ((xyz1_8(1)*xyz2_8(2)) - (xyz2_8(1)*xyz1_8(2)))**2)

      !- 2.2 Compute c = norm(-r1) & d = norm(r4)
      c_8 = sqrt ( xyz1_8(1)**2 + xyz1_8(2)**2 + xyz1_8(3)**2 )
    
      d_8 = sqrt ( ( ( (a_8*xyz1_8(1)) - xyz2_8(1) ) / b_8 )**2 + &
                 ( ( (  a_8*xyz1_8(2)) - xyz2_8(2) ) / b_8 )**2 + &
                 ( ( (  a_8*xyz1_8(3)) - xyz2_8(3) ) / b_8 )**2  )

      !- 2.3 Compute the forward rotation matrix
      grd_rot_8(1,1) = -xyz1_8(1) / c_8
      grd_rot_8(1,2) = -xyz1_8(2) / c_8
      grd_rot_8(1,3) = -xyz1_8(3) / c_8
      grd_rot_8(2,1) = ( ((a_8*xyz1_8(1)) - xyz2_8(1)) / b_8 ) / d_8
      grd_rot_8(2,2) = ( ((a_8*xyz1_8(2)) - xyz2_8(2)) / b_8 ) / d_8
      grd_rot_8(2,3) = ( ((a_8*xyz1_8(3)) - xyz2_8(3)) / b_8 ) / d_8
      grd_rot_8(3,1) = ( (xyz1_8(2)*xyz2_8(3)) - (xyz2_8(2)*xyz1_8(3))) / b_8
      grd_rot_8(3,2) = ( (xyz2_8(1)*xyz1_8(3)) - (xyz1_8(1)*xyz2_8(3))) / b_8
      grd_rot_8(3,3) = ( (xyz1_8(1)*xyz2_8(2)) - (xyz2_8(1)*xyz1_8(2))) / b_8

    end if

    do j1 = 1, msize
      do j2 = 1, msize
        write(*,*) 'sugrdpar: grd_rot_8(j1,j2) =',j1,j2,grd_rot_8(j1,j2)
      end do
    end do

    !
    !- 3. Compute the inverse matrix (transpose is the inverse)
    !
    do j1 = 1, msize
      do j2 = 1, msize
        grd_rotinv_8(j1,j2) = grd_rot_8(j2,j1)
      end do
    end do

    zunit(:,:) = 0.d0
    call mxma8x(zunit,grd_rotinv_8,grd_rot_8,msize,msize,msize)
    do j1 = 1, msize
      do j2 = 1, msize
        write(*,*) 'sugrdpar: unit = ', j1, j2, zunit(j1,j2)
      end do
    end do

  end SUBROUTINE sugrdpar

!--------------------------------------------------------------------------
! VLLACAR
!--------------------------------------------------------------------------
  subroutine vllacar(F_xyz_8, F_lon, F_lat)
    IMPLICIT NONE

    real(8), intent(out) :: F_xyz_8(msize)
    real(8), intent(in)  :: F_lon, F_lat ! In degrees

    F_xyz_8(1) = cos(MPC_RADIANS_PER_DEGREE_R8*F_lat) * cos(MPC_RADIANS_PER_DEGREE_R8*F_lon)
    F_xyz_8(2) = cos(MPC_RADIANS_PER_DEGREE_R8*F_lat) * sin(MPC_RADIANS_PER_DEGREE_R8*F_lon)
    F_xyz_8(3) = sin(MPC_RADIANS_PER_DEGREE_R8*F_lat)

  end subroutine vllacar

!--------------------------------------------------------------------------
! MXMA8X
!--------------------------------------------------------------------------
  subroutine mxma8x(pmat3,pmat1,pmat2,kdimi1,kdimj1,kdimj2)
    !
    !- Matrix times matrix.
    !
    IMPLICIT NONE
    integer kdimi1,kdimj1,kdimj2
    real*8 pmat3(kdimi1,kdimj2)
    real*8 pmat1(kdimi1,kdimj1)
    real*8 pmat2(kdimj1,kdimj2)

    integer ji1,jj2,jj

    pmat3(:,:) = 0.d0

    do jj2 = 1, kdimj2
      do jj = 1, kdimj1
        do ji1 = 1, kdimi1
          pmat3(ji1,jj2) = pmat3(ji1,jj2) + pmat1(ji1,jj) * pmat2(jj,jj2)
        end do
      end do
    end do
 
  end subroutine mxma8x

!--------------------------------------------------------------------------
! UVR_RotateWind
!--------------------------------------------------------------------------
  subroutine uvr_RotateWind( uwind, vwind, Lat, Lon, LatRot, LonRot, mode, &
                             nlev )
    !
    !- Go from tangential wind components from one sphere to another
    !  (same origin!). Fast version used by Variational analysis. 
    !
    IMPLICIT NONE

    integer, intent(in)      :: nlev
    real(8), intent(in)      :: Lat, Lon       ! In radians
    real(8), intent(in)      :: LatRot, LonRot ! In radians
    real(8), intent(inout)   :: uwind(nlev), vwind(nlev)
    character(*), intent(in) :: mode ! ToMetWind or ToRotWind

    real(8) :: coslatr, sinlatr, coslonr, sinlonr, coslon, sinlon, rsinlat

    real(8), allocatable :: xyz(:,:)
    real(8), allocatable :: uvcart(:,:)

    if ( .not. initialized ) then
      write(*,*)
      write(*,*) 'uvr_RotateWind: WindRotation module is not initialize'
      stop
    endif

    allocate(xyz   (msize,nlev))
    allocate(uvcart(msize,nlev))

    coslatr = cos(LatRot)
    sinlatr = sin(LatRot)
    coslonr = cos(LonRot)
    sinlonr = sin(LonRot)
    coslon  = cos(Lon)
    sinlon  = sin(Lon)
    rsinlat = 1.0d0/sin(Lat)
    
    if ( trim(mode) == 'ToMetWind' ) then 

      call fuvacart( xyz,                            & ! OUT
                     uwind, vwind, sinlonr, coslonr, & ! IN
                     sinlatr, coslatr, nlev )          ! IN
    
      call mxm8( grd_rotinv_8, msize, xyz, msize, & ! IN
                 uvcart,                          & ! OUT
                 nlev )                             ! IN

      call fcartauv( uwind, vwind,                         & ! OUT
                     uvcart, coslon, sinlon, rsinlat, nlev ) ! IN

    else if ( trim(mode) == 'ToRotWind' ) then
      write(*,*) 
      write(*,*) 'uvr_RotateWind: mode ToRotWind is not available yet'
      stop
    else
      write(*,*) 
      write(*,*) 'uvr_RotateWind: Unknown transform name: ', trim(mode)
      write(*,*) '                mode = ToMetWind or ToRotWind'
      stop
    end if

    deallocate(xyz)
    deallocate(uvcart)

  end subroutine uvr_RotateWind

!--------------------------------------------------------------------------
! UVR_RotateWindAdj
!--------------------------------------------------------------------------
  subroutine uvr_RotateWindAdj( uwind, vwind, Lat, Lon, LatRot, LonRot, mode, &
                                nlev )
    !
    !- Adjoint of : Go from tangential wind components from one sphere to another
    !  (same origin!). Fast version used by Variational analysis. 
    !
    IMPLICIT NONE

    integer, intent(in)      :: nlev
    real(8), intent(in)      :: Lat, Lon       ! In radians
    real(8), intent(in)      :: LatRot, LonRot ! In radians
    real(8), intent(inout)   :: uwind(nlev), vwind(nlev)
    character(*), intent(in) :: mode ! ToMetWind or ToRotWind

    real(8) :: coslatr, sinlatr, coslonr, sinlonr, coslon, sinlon, rsinlat

    real(8), allocatable :: xyz(:,:)
    real(8), allocatable :: uvcart(:,:)

    if ( .not. initialized ) then
      write(*,*)
      write(*,*) 'uvr_RotateWindAdj: WindRotation module is not initialize'
      stop
    endif

    allocate(xyz   (msize,nlev))
    allocate(uvcart(msize,nlev))

    coslatr = cos(LatRot)
    sinlatr = sin(LatRot)
    coslonr = cos(LonRot)
    sinlonr = sin(LonRot)
    coslon  = cos(Lon)
    sinlon  = sin(Lon)
    rsinlat = 1.0d0/sin(Lat)

    if ( trim(mode) == 'ToMetWind' ) then 

      call fcartauvAdj( uvcart,                                     & ! OUT
                        uwind, vwind, coslon, sinlon, rsinlat, nlev ) ! IN

      call mxm8( grd_rot_8, msize, uvcart, msize, & ! IN
                 xyz,                             & ! OUT
                 nlev )                             ! IN

      call fuvacartAdj( uwind, vwind,             & ! OUT
                        xyz, sinlonr, coslonr,    & ! IN
                        sinlatr, coslatr, nlev )    ! IN

    else if ( trim(mode) == 'ToRotWind' ) then
      write(*,*) 
      write(*,*) 'uvr_RotateWindAdj: mode ToRotWind is not available yet'
      stop
    else
      write(*,*) 
      write(*,*) 'uvr_RotateWindAdj: Unknown transform name: ', trim(mode)
      write(*,*) '                   mode = ToMetWind or ToRotWind'
      stop
    end if

    deallocate(xyz)
    deallocate(uvcart)

  end subroutine uvr_RotateWindAdj

!--------------------------------------------------------------------------
! FUVACART
!--------------------------------------------------------------------------
  subroutine fuvacart(pcart,pu,pv,psinlon,pcoslon,psinlat,pcoslat,knk)
    !
    !- Compute the winds in the cartesian space from
    !  the tangential wind vector components.
    !  Fast version of uvacart special to Var analysis.
    !
    !arguments
    !    out    pvcart - Cartesian components of the tangential wind vector
    !    IN     pu     - Zonal component of the wind on the tangential plane.
    !           pv     - Meridional component of the wind on the tangential plane.
    !
    IMPLICIT NONE

    integer, intent(in)  :: knk
    real(8), intent(in)  :: psinlon,pcoslon,psinlat,pcoslat
    real(8), intent(out) :: pcart(msize,knk)
    real(8), intent(in)  :: pu(knk),pv(knk)

    integer :: ji,jj,jk
    real(8) :: zfac1,zfac2

    zfac1 = pcoslon*psinlat
    zfac2 = psinlon*psinlat

    do jk = 1, knk
       pcart(1,jk) = -(pu(jk)*psinlon) - (pv(jk)*zfac1)
       pcart(2,jk) =  (pu(jk)*pcoslon) - (pv(jk)*zfac2)
       pcart(3,jk) =   pv(jk)*pcoslat
    end do

  end subroutine fuvacart

!--------------------------------------------------------------------------
! FUCACARTAdj
!--------------------------------------------------------------------------
  subroutine fuvacartAdj(pu,pv,pcart,psinlonr,pcoslonr,psinlatr,pcoslatr,knk)
    !
    !- Adjoint of fuvacart.
    !
    IMPLICIT NONE

    integer, intent(in)  ::  knk
    real(8), intent(in)  :: psinlonr,pcoslonr,psinlatr,pcoslatr
    real(8), intent(inout)  :: pcart(msize,knk)
    real(8), intent(out) :: pu(knk), pv(knk)
    
    integer :: jk
    real(8) :: zfac1,zfac2

    do jk = 1, knk
      pu(jk)=0.0d0  ! dont use sub. zero.ftn here... *4 versus *8 diff...
      pv(jk)=0.0d0
    end do

    zfac1=psinlonr*psinlatr
    zfac2=pcoslonr*psinlatr

    do jk =1, knk
      pu(jk) = pu(jk)+pcart(2,jk)*pcoslonr
      pv(jk) = pv(jk)-pcart(2,jk)*zfac1
      pcart(2,jk) = 0.0d0
      pu(jk) = pu(jk)-pcart(1,jk)*psinlonr
      pv(jk) = pv(jk)-pcart(1,jk)*zfac2
      pcart(1,jk) = 0.0d0
      pv(jk) = pv(jk)+pcart(3,jk)*pcoslatr
      pcart(3,jk) = 0.0d0
    end do

  end subroutine fuvacartAdj

!--------------------------------------------------------------------------
! FCARTAUV
!--------------------------------------------------------------------------
  subroutine fcartauv(pu,pv,pcart,pcoslonob,psinlonob,prsinlatob,knk)
    !
    !- Compute the components of the winds in the rotated system of
    !  coordinates from the winds in the rotated cartesian space
    !  Fast version of cartauv used by Variational analysis.
    !
    IMPLICIT NONE

    integer, intent(in)  :: knk

    real(8), intent(in)  :: pcart(msize,knk)
    real(8), intent(out) :: pu(knk), pv(knk)
    real(8), intent(in)  :: pcoslonob, psinlonob, prsinlatob

    integer :: jk
    real(8) :: zfac1, zfac2

    zfac1 = -prsinlatob*pcoslonob
    zfac2 = -prsinlatob*psinlonob

    do jk=1,knk
      pu(jk) = -psinlonob*pcart(1,jk) + pcoslonob*pcart(2,jk)
      pv(jk) = zfac1*pcart(1,jk) + zfac2*pcart(2,jk)
    end do

  end subroutine fcartauv

!--------------------------------------------------------------------------
! FCARTAUVAdj
!--------------------------------------------------------------------------
  subroutine fcartauvAdj(pcart,pu,pv,pcoslonob,psinlonob,prsinlatob,knk)
    !
    !- Adjoint of cartauv.
    !
    IMPLICIT NONE
    integer, intent(in)    :: knk
    real(8), intent(in)    :: pcoslonob, psinlonob, prsinlatob
    real(8), intent(inout) :: pu(knk), pv(knk)
    real(8), intent(out)   :: pcart(msize,knk)

    integer :: jk
    real(8) ::  zfac1,zfac2

    do jk=1,knk
      pcart(1,jk)=0.0d0
      pcart(2,jk)=0.0d0
      pcart(3,jk)=0.0d0
    end do

    zfac1=prsinlatob*pcoslonob
    zfac2=prsinlatob*psinlonob

    do jk = 1, knk
      pcart(1,jk) = pcart(1,jk)-pv(jk)*zfac1
      pcart(2,jk) = pcart(2,jk)-pv(jk)*zfac2
      pv(jk) = 0.0d0
      pcart(2,jk) = pcart(2,jk)+pu(jk)*pcoslonob
      pcart(1,jk) = pcart(1,jk)-pu(jk)*psinlonob
      pu(jk) = 0.0d0
    end do

  end subroutine fcartauvAdj

!--------------------------------------------------------------------------
! MXM8
!--------------------------------------------------------------------------
  subroutine mxm8(a,nar,b,nac,c,nbc)
    implicit none

    integer, intent(in)  :: nar, nac, nbc
    real(8), intent(in)  :: a(nar,nac), b(nac,nbc)
    real(8), intent(out) :: c(nar,nbc)

    integer :: i,j,k

    do j = 1, nbc
       do i = 1, nar
          c(i,j) = 0.0d0
          do k = 1, nac
             c(i,j) = c(i,j) + a(i,k)*b(k,j)
          end do
       end do
    end do

  end subroutine mxm8

!--------------------------------------------------------------------------
! uvr_RotateLatLon
!--------------------------------------------------------------------------
  subroutine uvr_RotateLatLon(LatOut,LonOut,LatIn,LonIn,mode)
    !
    !- Go from (lat,lon) of one Cartesian frame to (lat,lon) of another
    !  Cartesian frame given the rotation matrix.
    !
    IMPLICIT NONE

    real(8), intent(in)          :: LatIn , LonIn   ! In radians
    real(8), intent(out)         :: LatOut, LonOut  ! In radians

    character(len=*), intent(in) :: mode

    real(8) :: CartIn(msize),CartOut(msize)
    real(8) :: rLon,rLat

    if ( .not. initialized ) then
      write(*,*)
      write(*,*) 'uvr_RotateLatLon: WindRotation module is not initialize'
      stop
    endif

    rLon = LonIn * MPC_DEGREES_PER_RADIAN_R8 ! To degress
    rLat = LatIn * MPC_DEGREES_PER_RADIAN_R8 ! To degrees

    call vllacar( CartIn,      & ! OUT
                  rLon, rLat )   ! IN

    if ( trim(mode) == 'ToLatLonRot' ) then
      call mxv( CartOut,              & ! OUT
                grd_rot_8, CartIn,    & ! IN
                msize, msize)           ! IN
    else if ( trim(mode) == 'ToLatLon' ) then 
      call mxv( CartOut,              & ! OUT
                grd_rotinv_8, CartIn, & ! IN
                msize, msize)           ! IN
    else
      write(*,*) 
      write(*,*) 'uvr_RotateLatLon: Unknown transform name: ', trim(mode)
      write(*,*) '                  mode = ToLatLonRot or ToLatLon'
      stop
    end if

    call carall( LonOut, LatOut, & ! OUT
                 CartOut )         ! IN

    LonOut = LonOut * MPC_RADIANS_PER_DEGREE_R8 ! To radians
    LatOut = LatOut * MPC_RADIANS_PER_DEGREE_R8 ! To radians

  end subroutine uvr_RotateLatLon

!--------------------------------------------------------------------------
! CARALL
!--------------------------------------------------------------------------
  subroutine carall(plon,plat,pcart)
    !
    !- Returns (lat,lon) (degrees) of an input Cartesian position vector 
    !  on the unit sphere.
    !
    IMPLICIT NONE

    real(8), intent(out) :: plat, plon
    real(8), intent(in)  :: pcart(msize)

    plat = asin(pcart(3))
    plat = plat * MPC_DEGREES_PER_RADIAN_R8
    
    if ( pcart(1) == 0.d0 ) then
       if ( pcart(2) == 0.d0 ) then  ! point is located at the pole
          plon = 0.0d0  ! can be any longitude... set it to zero simply.
       else if ( pcart(2) > 0.0d0 ) then
          plon = 90.0d0
       else if ( pcart(2) < 0.0d0 ) then
          plon = 270.0d0
       end if
    else
       plon = atan2(pcart(2),pcart(1))
       if (plon < 0.0d0) plon = plon + 2.d0 * MPC_PI_R8
       plon = plon * MPC_DEGREES_PER_RADIAN_R8
    end if

  end subroutine carall

!--------------------------------------------------------------------------
! MXV
!--------------------------------------------------------------------------
  subroutine mxv(pvec2,pmat,pvec1,kdimi,kdimj)
    !
    !- Matrix times vector.
    !
    IMPLICIT NONE

    integer :: kdimi, kdimj
    real(8) :: pvec2(kdimi), pmat(kdimi,kdimj), pvec1(kdimj)

    integer :: ji,jj

    pvec2(:) = 0.0d0
    do jj = 1,kdimj
      do ji = 1,kdimi
        pvec2(ji) = pvec2(ji) + pmat(ji,jj) * pvec1(jj)
      end do
    end do

  end subroutine mxv

end module WindRotation_mod
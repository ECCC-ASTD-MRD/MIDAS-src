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
!! MODULE windRotation (prefix="uvr")
!!
!! *Purpose*: To transform winds FROM the rotated spherical 
!!            coordinate system TO the non-rotated spherical coordinate system.
!!
!--------------------------------------------------------------------------
module windRotation_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use utilities_mod
  implicit none
  save
  private

  ! Public derived type definition
  public :: struct_uvr

  ! Public Subroutines
  public :: uvr_Setup, uvr_RotateWind, uvr_RotateWindAdj, uvr_RotateLatLon

  integer, parameter :: msize = 3
  integer, parameter :: maxNumSubGrid = 2

  type :: struct_uvr
    real(8) :: grd_rot_8   (msize,msize,maxNumSubGrid) ! Real Earth to Rotated Earth
    real(8) :: grd_rotinv_8(msize,msize,maxNumSubGrid) ! Rotated Earth to Real Earth
    logical :: initialized = .false.
  end type struct_uvr

  contains

  !--------------------------------------------------------------------------
  ! uvr_Setup
  !--------------------------------------------------------------------------
  subroutine uvr_Setup( uvr, hco_in )
    implicit none

    ! arguments
    type(struct_uvr), pointer    :: uvr
    type(struct_hco), intent(in) :: hco_in

    !
    !-  Compute the rotation matrices (grd_rot_8 and grd_rotinv_8)
    !
    write(*,*)
    write(*,*) 'uvr_Setup: Starting...  for grid type = ', hco_in%grtyp

    if ( associated(uvr) ) then
      call utl_abort('uvr_setup: the supplied uvr pointer is not null!')
    endif

    allocate(uvr)

    write(*,*) hco_in % xlon1
    write(*,*) hco_in % xlat1
    write(*,*) hco_in % xlon2
    write(*,*) hco_in % xlat2

    call sugrdpar(uvr, 1, hco_in % xlon1, hco_in % xlat1, hco_in % xlon2, hco_in % xlat2 )

    if ( hco_in%grtyp == 'U' ) then

      write(*,*) 'uvr_setup: doing setup for YAN grid'
      write(*,*) hco_in % xlon1_yan
      write(*,*) hco_in % xlat1_yan
      write(*,*) hco_in % xlon2_yan
      write(*,*) hco_in % xlat2_yan

      call sugrdpar(uvr, 2, hco_in % xlon1_yan, hco_in % xlat1_yan, hco_in % xlon2_yan, hco_in % xlat2_yan )

    end if

    uvr%initialized = .true.

    write(*,*)
    write(*,*) 'uvr_Setup: Done!'

  end subroutine uvr_Setup

  !--------------------------------------------------------------------------
  ! SUGRDPAR
  !--------------------------------------------------------------------------
  subroutine sugrdpar( uvr, subGridIndex, grd_xlon1, grd_xlat1, grd_xlon2, grd_xlat2 )
    !
    ! Compute the rotation matrix (r_8) that allows transformation
    ! from the non-rotated to the rotated spherical coordinate system.
    implicit none

    ! arguments
    type(struct_uvr), pointer :: uvr
    integer                   :: subGridIndex
    real(8), intent(in)       :: grd_xlon1, grd_xlat1, grd_xlon2, grd_xlat2 

    ! locals
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
      uvr%grd_rot_8(1,1,subGridIndex) = 1.0d0       
      uvr%grd_rot_8(1,2,subGridIndex) = 0.0d0
      uvr%grd_rot_8(1,3,subGridIndex) = 0.0d0
      uvr%grd_rot_8(2,1,subGridIndex) = 0.0d0
      uvr%grd_rot_8(2,2,subGridIndex) = 1.0d0
      uvr%grd_rot_8(2,3,subGridIndex) = 0.0d0
      uvr%grd_rot_8(3,1,subGridIndex) = 0.0d0
      uvr%grd_rot_8(3,2,subGridIndex) = 0.0d0
      uvr%grd_rot_8(3,3,subGridIndex) = 1.0d0

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
      uvr%grd_rot_8(1,1,subGridIndex) = -xyz1_8(1) / c_8
      uvr%grd_rot_8(1,2,subGridIndex) = -xyz1_8(2) / c_8
      uvr%grd_rot_8(1,3,subGridIndex) = -xyz1_8(3) / c_8
      uvr%grd_rot_8(2,1,subGridIndex) = ( ((a_8*xyz1_8(1)) - xyz2_8(1)) / b_8 ) / d_8
      uvr%grd_rot_8(2,2,subGridIndex) = ( ((a_8*xyz1_8(2)) - xyz2_8(2)) / b_8 ) / d_8
      uvr%grd_rot_8(2,3,subGridIndex) = ( ((a_8*xyz1_8(3)) - xyz2_8(3)) / b_8 ) / d_8
      uvr%grd_rot_8(3,1,subGridIndex) = ( (xyz1_8(2)*xyz2_8(3)) - (xyz2_8(2)*xyz1_8(3))) / b_8
      uvr%grd_rot_8(3,2,subGridIndex) = ( (xyz2_8(1)*xyz1_8(3)) - (xyz1_8(1)*xyz2_8(3))) / b_8
      uvr%grd_rot_8(3,3,subGridIndex) = ( (xyz1_8(1)*xyz2_8(2)) - (xyz2_8(1)*xyz1_8(2))) / b_8

    end if

    do j1 = 1, msize
      do j2 = 1, msize
        write(*,*) 'sugrdpar: grd_rot_8(j1,j2) =',j1,j2,uvr%grd_rot_8(j1,j2,subGridIndex)
      end do
    end do

    !
    !- 3. Compute the inverse matrix (transpose is the inverse)
    !
    do j1 = 1, msize
      do j2 = 1, msize
        uvr%grd_rotinv_8(j1,j2,subGridIndex) = uvr%grd_rot_8(j2,j1,subGridIndex)
      end do
    end do

    zunit(:,:) = 0.d0
    call mxma8x(zunit,uvr%grd_rotinv_8(:,:,subGridIndex),uvr%grd_rot_8(:,:,subGridIndex),msize,msize,msize)
    do j1 = 1, msize
      do j2 = 1, msize
        write(*,*) 'sugrdpar: unit = ', j1, j2, zunit(j1,j2)
      end do
    end do

  end subroutine sugrdpar

  !--------------------------------------------------------------------------
  ! vllacar
  !--------------------------------------------------------------------------
  subroutine vllacar( F_xyz_8, F_lon, F_lat )
    implicit none

    real(8), intent(out) :: F_xyz_8(msize)
    real(8), intent(in)  :: F_lon, F_lat ! In degrees

    F_xyz_8(1) = cos(MPC_RADIANS_PER_DEGREE_R8*F_lat) * cos(MPC_RADIANS_PER_DEGREE_R8*F_lon)
    F_xyz_8(2) = cos(MPC_RADIANS_PER_DEGREE_R8*F_lat) * sin(MPC_RADIANS_PER_DEGREE_R8*F_lon)
    F_xyz_8(3) = sin(MPC_RADIANS_PER_DEGREE_R8*F_lat)

  end subroutine vllacar

  !--------------------------------------------------------------------------
  ! mxma8x
  !--------------------------------------------------------------------------
  subroutine mxma8x( pmat3, pmat1, pmat2, kdimi1, kdimj1, kdimj2 )
    !
    ! Matrix times matrix.
    !
    implicit none

    ! arguments
    integer :: kdimi1,kdimj1,kdimj2
    real(8) :: pmat3(kdimi1,kdimj2)
    real(8) :: pmat1(kdimi1,kdimj1)
    real(8) :: pmat2(kdimj1,kdimj2)

    ! locals
    integer :: ji1,jj2,jj

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
  ! uvr_rotateWind
  !--------------------------------------------------------------------------
  subroutine uvr_rotateWind( uvr, subGridIndex, uwind, vwind, Lat, Lon, LatRot, LonRot, mode )
    !
    ! Go from tangential wind components from one sphere to another
    ! (same origin!). Fast version used by Variational analysis. 
    !
    implicit none

    ! arguments
    type(struct_uvr), pointer :: uvr
    integer, intent(in)       :: subGridIndex
    real(8), intent(in)       :: Lat, Lon       ! In radians
    real(8), intent(in)       :: LatRot, LonRot ! In radians
    real(8), intent(inout)    :: uwind, vwind
    character(*), intent(in)  :: mode ! ToMetWind or ToRotWind

    ! locals
    integer :: index1, index2
    real(8) :: coslatr, sinlatr, coslonr, sinlonr, coslon, sinlon, rsinlat
    real(8) :: xyz(msize), uvcart(msize)

    if ( .not. uvr%initialized ) then
      write(*,*)
      call utl_abort('uvr_RotateWind: WindRotation module is not initialize')
    endif

    coslatr = cos(LatRot)
    sinlatr = sin(LatRot)
    coslonr = cos(LonRot)
    sinlonr = sin(LonRot)
    coslon  = cos(Lon)
    sinlon  = sin(Lon)
    if ( Lat /= 0.0d0 ) then
      rsinlat = 1.0d0/sin(Lat)
    else
      !write(*,*) 'rotateWind: latitude is zero!', Lat
      rsinlat = 1.0d0/sin(1.0d-8)
    end if
    
    if ( trim(mode) == 'ToMetWind' ) then 

      xyz(1) = -uwind*sinlonr - vwind*coslonr*sinlatr
      xyz(2) =  uwind*coslonr - vwind*sinlonr*sinlatr
      xyz(3) =                  vwind*coslatr
    
      uvcart(:) = 0.0d0
      do index2 = 1, msize
        do index1 = 1, msize
          uvcart(index1) = uvcart(index1) +   &
               uvr%grd_rotinv_8(index1,index2,subGridIndex)*xyz(index2)
        end do
      end do

      uwind = -uvcart(1)*sinlon         + uvcart(2)*coslon
      vwind = -uvcart(1)*rsinlat*coslon - uvcart(2)*rsinlat*sinlon

    else if ( trim(mode) == 'ToRotWind' ) then
      write(*,*) 
      call utl_abort('uvr_RotateWind: mode ToRotWind is not available yet')
    else
      write(*,*) 
      write(*,*) 'uvr_RotateWind: Unknown transform name: ', trim(mode)
      write(*,*) '                mode = ToMetWind or ToRotWind'
      call utl_abort('uvr_RotateWind')
    end if

  end subroutine uvr_RotateWind

!--------------------------------------------------------------------------
! uvr_rotateWindAdj
!--------------------------------------------------------------------------
  subroutine uvr_rotateWindAdj( uvr, subGridIndex, uwind, vwind, Lat, Lon, LatRot, LonRot, mode )
    !
    ! Adjoint of : Go from tangential wind components from one sphere to another
    ! (same origin!). Fast version used by Variational analysis. 
    !
    implicit none

    ! arguments
    type(struct_uvr), pointer :: uvr
    integer, intent(in)       :: subGridIndex
    real(8), intent(in)       :: Lat, Lon       ! In radians
    real(8), intent(in)       :: LatRot, LonRot ! In radians
    real(8), intent(inout)    :: uwind, vwind
    character(*), intent(in)  :: mode ! ToMetWind or ToRotWind

    ! locals
    integer :: index1, index2
    real(8) :: coslatr, sinlatr, coslonr, sinlonr, coslon, sinlon, rsinlat
    real(8) :: xyz(msize), uvcart(msize)

    if ( .not. uvr%initialized ) then
      write(*,*)
      call utl_abort('uvr_RotateWindAdj: WindRotation module is not initialize')
    endif

    coslatr = cos(LatRot)
    sinlatr = sin(LatRot)
    coslonr = cos(LonRot)
    sinlonr = sin(LonRot)
    coslon  = cos(Lon)
    sinlon  = sin(Lon)
    if ( Lat /= 0.0d0 ) then
      rsinlat = 1.0d0/sin(Lat)
    else
      !write(*,*) 'rotateWindAdj: latitude is zero!', Lat
      rsinlat = 1.0d0/sin(1.0d-8)
    end if

    if ( trim(mode) == 'ToMetWind' ) then 

      uvcart(1) = -uwind*sinlon - vwind*rsinlat*coslon
      uvcart(2) =  uwind*coslon - vwind*rsinlat*sinlon
      uvcart(3) = 0.0d0

      xyz(:) = 0.0d0
      do index2 = 1, msize
        do index1 = 1, msize
          xyz(index1) = xyz(index1) + uvr%grd_rot_8(index1,index2,subGridIndex)*uvcart(index2)
        end do
      end do

      uwind = -xyz(1)*sinlonr         + xyz(2)*coslonr
      vwind = -xyz(1)*coslonr*sinlatr - xyz(2)*sinlonr*sinlatr + xyz(3)*coslatr

    else if ( trim(mode) == 'ToRotWind' ) then
      write(*,*) 
      call utl_abort('uvr_RotateWindAdj: mode ToRotWind is not available yet')
    else
      write(*,*) 
      write(*,*) 'uvr_RotateWindAdj: Unknown transform name: ', trim(mode)
      write(*,*) '                   mode = ToMetWind or ToRotWind'
      call utl_abort('uvr_RotateWindAdj')
    end if

  end subroutine uvr_RotateWindAdj

  !--------------------------------------------------------------------------
  ! uvr_rotateLatLon
  !--------------------------------------------------------------------------
  subroutine uvr_rotateLatLon( uvr, subGridIndex, LatOut, LonOut, LatIn, LonIn, mode )
    !
    ! Go from (lat,lon) of one Cartesian frame to (lat,lon) of another
    ! Cartesian frame given the rotation matrix.
    !
    implicit none

    ! arguments
    type(struct_uvr), pointer :: uvr
    integer, intent(in)       :: subGridIndex
    real(8), intent(in)       :: LatIn , LonIn   ! In radians
    real(8), intent(out)      :: LatOut, LonOut  ! In radians
    character(len=*), intent(in) :: mode

    ! locals
    real(8) :: CartIn(msize),CartOut(msize)
    real(8) :: rLon,rLat

    if ( .not. uvr%initialized ) then
      write(*,*)
      call utl_abort('uvr_RotateLatLon: WindRotation module is not initialize')
    endif

    rLon = LonIn * MPC_DEGREES_PER_RADIAN_R8 ! To degress
    rLat = LatIn * MPC_DEGREES_PER_RADIAN_R8 ! To degrees

    call vllacar( CartIn,      & ! OUT
                  rLon, rLat )   ! IN

    if ( trim(mode) == 'ToLatLonRot' ) then
      call mxv( CartOut,              & ! OUT
                uvr%grd_rot_8(:,:,subGridIndex), CartIn,    & ! IN
                msize, msize)           ! IN
    else if ( trim(mode) == 'ToLatLon' ) then 
      call mxv( CartOut,              & ! OUT
                uvr%grd_rotinv_8(:,:,subGridIndex), CartIn, & ! IN
                msize, msize)           ! IN
    else
      write(*,*) 
      write(*,*) 'uvr_RotateLatLon: Unknown transform name: ', trim(mode)
      write(*,*) '                  mode = ToLatLonRot or ToLatLon'
      call utl_abort('uvr_RotateLatLon')
    end if

    call carall( LonOut, LatOut, & ! OUT
                 CartOut )         ! IN

    LonOut = LonOut * MPC_RADIANS_PER_DEGREE_R8 ! To radians
    LatOut = LatOut * MPC_RADIANS_PER_DEGREE_R8 ! To radians

  end subroutine uvr_RotateLatLon

  !--------------------------------------------------------------------------
  ! carall
  !--------------------------------------------------------------------------
  subroutine carall( plon, plat, pcart )
    !
    ! Returns (lat,lon) (degrees) of an input Cartesian position vector 
    ! on the unit sphere.
    !
    implicit none

    ! arguments
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
  ! mxv
  !--------------------------------------------------------------------------
  subroutine mxv( pvec2, pmat, pvec1, kdimi, kdimj )
    !
    ! Matrix times vector.
    !
    implicit none

    ! arguments
    integer :: kdimi, kdimj
    real(8) :: pvec2(kdimi), pmat(kdimi,kdimj), pvec1(kdimj)

    ! locals
    integer :: ji,jj

    pvec2(:) = 0.0d0
    do jj = 1,kdimj
      do ji = 1,kdimi
        pvec2(ji) = pvec2(ji) + pmat(ji,jj) * pvec1(jj)
      end do
    end do

  end subroutine mxv

end module windRotation_mod

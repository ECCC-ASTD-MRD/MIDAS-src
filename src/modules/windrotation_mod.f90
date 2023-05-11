
module windRotation_mod
  ! MODULE windRotation (prefix='uvr' category='4. Data Object transformations')
  !
  ! :Purpose: To transform winds FROM the rotated spherical coordinate system
  !           TO the non-rotated spherical coordinate system.
  !
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use utilities_mod
  use midasMpi_mod
  implicit none
  save
  private

  public :: struct_uvr ! Public derived type definition

  ! Public Subroutines
  public :: uvr_setup, uvr_rotateWind_nl, uvr_rotateWind_tl, uvr_rotateWind_ad, uvr_rotateLatLon

  integer, parameter :: msize = 3
  integer, parameter :: maxNumSubGrid = 2

  type :: struct_uvr
    real(8) :: grd_rot_8   (msize,msize,maxNumSubGrid) ! Real Earth to Rotated Earth
    real(8) :: grd_rotinv_8(msize,msize,maxNumSubGrid) ! Rotated Earth to Real Earth
    logical :: initialized = .false.
  end type struct_uvr

  ! Minimum value for latitude used in TL/AD wind rotation routines
  real(8), parameter :: lat_minValue = 1D-10

  contains

  subroutine uvr_setup( uvr, hco_in )
    !
    ! :Purpose:    Setup the information for wind rotation
    !
    implicit none

    ! Arguments:
    type(struct_uvr), pointer    :: uvr    ! Wind rotation object
    type(struct_hco), intent(in) :: hco_in ! Horizontal grid object

    if ( associated(uvr) ) then
      if ( uvr%initialized ) then
        if ( mmpi_myid == 0 ) write(*,*) 'uvr_setup: already initialized, returning'
        return
      else
        call utl_abort('uvr_setup: the supplied non-null uvr pointer is not initialized!')
      end if
    end if

    !
    !-  Compute the rotation matrices (grd_rot_8 and grd_rotinv_8)
    !
    if ( mmpi_myid == 0 ) then
      write(*,*)
      write(*,*) 'uvr_setup: Starting...  for grid type = ', hco_in%grtyp
    end if

    allocate(uvr)

    if ( mmpi_myid == 0 ) then
      write(*,*) hco_in % xlon1
      write(*,*) hco_in % xlat1
      write(*,*) hco_in % xlon2
      write(*,*) hco_in % xlat2
    end if

    call sugrdpar(uvr, 1, hco_in % xlon1, hco_in % xlat1, hco_in % xlon2, hco_in % xlat2 )

    if ( hco_in%grtyp == 'U' ) then

      if ( mmpi_myid == 0 ) then
        write(*,*) 'uvr_setup: doing setup for YAN grid'
        write(*,*) hco_in % xlon1_yan
        write(*,*) hco_in % xlat1_yan
        write(*,*) hco_in % xlon2_yan
        write(*,*) hco_in % xlat2_yan
      end if

      call sugrdpar(uvr, 2, hco_in % xlon1_yan, hco_in % xlat1_yan, hco_in % xlon2_yan, hco_in % xlat2_yan )

    end if

    uvr%initialized = .true.

    if ( mmpi_myid == 0 ) then
      write(*,*)
      write(*,*) 'uvr_setup: Done!'
    end if

  end subroutine uvr_setup

  subroutine sugrdpar( uvr, subGridIndex, grd_xlon1, grd_xlat1, grd_xlon2, grd_xlat2 )
    !
    ! :Purpose: Compute the rotation matrix (r_8) that allows transformation
    !           from the non-rotated to the rotated spherical coordinate system.
    !
    implicit none

    ! Arguments:
    type(struct_uvr), pointer :: uvr          ! Wind rotation object
    integer                   :: subGridIndex ! Horizontal subGrid index = 2
    real(8), intent(in)       :: grd_xlon1    ! Horizontal grid xlon1_yan
    real(8), intent(in)       :: grd_xlat1    ! Horizontal grid xlat1_yan
    real(8), intent(in)       :: grd_xlon2    ! Horizontal grid xlon2_yan
    real(8), intent(in)       :: grd_xlat2    ! Horizontal grid xlat2_yan

    ! locals
    integer :: j1, j2
    real(8) :: zxlon1_8,zxlat1_8,zxlon2_8,zxlat2_8
    real(8) :: a_8, b_8, c_8, d_8, xyz1_8(msize), xyz2_8(msize)
    real(8) :: zunit(msize,msize)

    zxlon1_8 = grd_xlon1
    zxlat1_8 = grd_xlat1
    zxlon2_8 = grd_xlon2
    zxlat2_8 = grd_xlat2

    if ( zxlon1_8 == zxlon2_8 .and. zxlat1_8 == zxlat2_8 ) then
      !
      !- 1. Non rotated grid (for test case only)
      !
      if ( mmpi_myid == 0 ) then
        write(*,*)
        write(*,*) 'uvr_sugrdpar: Warning: This grid is not rotated !!!'
      end if
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
        if ( mmpi_myid == 0 ) then
          write(*,*) 'sugrdpar: grd_rot_8(j1,j2) =',j1,j2,uvr%grd_rot_8(j1,j2,subGridIndex)
        end if
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
        if ( mmpi_myid == 0 ) then
          write(*,*) 'sugrdpar: unit = ', j1, j2, zunit(j1,j2)
        end if
      end do
    end do

  end subroutine sugrdpar

  subroutine vllacar( F_xyz_8, F_lon, F_lat )
    !
    ! :Purpose:    Compute parameters of rotated grid
    ! 
    implicit none

    real(8), intent(out) :: F_xyz_8(msize) ! output
    real(8), intent(in)  :: F_lon          ! Input in degrees  
    real(8), intent(in)  :: F_lat          ! Input in degrees

    F_xyz_8(1) = cos(MPC_RADIANS_PER_DEGREE_R8*F_lat) * cos(MPC_RADIANS_PER_DEGREE_R8*F_lon)
    F_xyz_8(2) = cos(MPC_RADIANS_PER_DEGREE_R8*F_lat) * sin(MPC_RADIANS_PER_DEGREE_R8*F_lon)
    F_xyz_8(3) = sin(MPC_RADIANS_PER_DEGREE_R8*F_lat)

  end subroutine vllacar

  subroutine mxma8x( pmat3, pmat1, pmat2, kdimi1, kdimj1, kdimj2 )
    !
    ! :Purpose:    Compute a product of two matrices.
    !
    implicit none

    ! Arguments:
    real(8) :: pmat3(kdimi1,kdimj2)  ! output
    real(8) :: pmat1(kdimi1,kdimj1)  ! input matrix one
    real(8) :: pmat2(kdimj1,kdimj2)  ! input matrix two
    integer :: kdimi1                ! first  dimension of the first  matrix
    integer :: kdimj1                ! second dimension of the first  matrix 
    integer :: kdimj2                ! second dimension of the second matrix

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

  subroutine uvr_rotateWind_nl( uvr, subGridIndex, uwind, vwind, Lat, Lon, LatRot, LonRot, mode )
    !
    ! :Purpose: Go from tangential wind components from one sphere to another
    !           (same origin!). Original ezsint version used for computing innovation.
    !
    implicit none

    ! Arguments:
    type(struct_uvr), pointer :: uvr          ! Wind rotation object
    integer, intent(in)       :: subGridIndex ! Current subgrid index
    real(8), intent(inout)    :: uwind        ! interpUU
    real(8), intent(inout)    :: vwind        ! interpVV
    real(8), intent(in)       :: Lat          ! Latitude in radians
    real(8), intent(in)       :: Lon          ! Longitude in radians
    real(8), intent(in)       :: LatRot       ! Rotated latitude in radians
    real(8), intent(in)       :: LonRot       ! Rotated longitude in radians 
    character(*), intent(in)  :: mode         ! ToMetWind or ToRotWind

    ! locals
    integer :: index1, index2
    real(8) :: coslatr, sinlatr, coslonr, sinlonr, coslat, sinlat, coslon, sinlon, ezCoeff_C, ezCoeff_D
    real(8) :: xyz(msize), uvcart(msize)

    if ( .not. uvr%initialized ) then
      write(*,*)
      call utl_abort('uvr_rotateWind_nl: WindRotation module is not initialize')
    endif

    coslatr = cos(LatRot)
    sinlatr = sin(LatRot)
    coslonr = cos(LonRot)
    sinlonr = sin(LonRot)
    coslat  = cos(Lat)
    sinlat  = sin(Lat)
    coslon  = cos(Lon)
    sinlon  = sin(Lon)
    
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

      uwind  = uvcart(2)*coslon - uvcart(1)*sinlon
      ezCoeff_C = uvcart(1)*coslon + uvcart(2)*sinlon
      ezCoeff_D = sqrt( ezCoeff_C**2 + uvcart(3)**2 )
      vwind     = sign( ezCoeff_D, uvcart(3)*coslat - ezCoeff_C*sinlat )

    else if ( trim(mode) == 'ToRotWind' ) then
      write(*,*) 
      call utl_abort('uvr_rotateWind_nl: mode ToRotWind is not available yet')
    else
      write(*,*) 
      write(*,*) 'uvr_rotateWind_nl: Unknown transform name: ', trim(mode)
      write(*,*) '                mode = ToMetWind or ToRotWind'
      call utl_abort('uvr_rotateWind_nl')
    end if

  end subroutine uvr_rotateWind_nl

  subroutine uvr_rotateWind_tl( uvr, subGridIndex, uwind, vwind, Lat_in, Lon_in, LatRot_in, LonRot_in, mode )
    !
    ! :Purpose: Go from tangential wind components from one sphere to another
    !           (same origin!). Fast version used by Variational analysis. 
    !
    implicit none

    ! Arguments:
    type(struct_uvr), pointer :: uvr          ! Wind rotation object
    integer, intent(in)       :: subGridIndex ! Current subgrid index
    real(8), intent(inout)    :: uwind        ! interpUU
    real(8), intent(inout)    :: vwind        ! interpVV
    real(8), intent(in)       :: Lat_in       ! Latitude in radians
    real(8), intent(in)       :: Lon_in       ! Longitude in radians
    real(8), intent(in)       :: LatRot_in    ! Rotated latitude in radians
    real(8), intent(in)       :: LonRot_in    ! Rotated longitude in radians 
    character(*), intent(in)  :: mode         ! ToMetWind or ToRotWind

    ! locals
    integer :: index1, index2
    real(8) :: coslatr, sinlatr, coslonr, sinlonr, coslon, sinlon, rsinlat
    real(8) :: lat, lon, latRot, lonRot
    real(8) :: xyz(msize), uvcart(msize)

    if ( .not. uvr%initialized ) then
      write(*,*)
      call utl_abort('uvr_rotateWind_tl: WindRotation module is not initialize')
    endif

    ! Special case for the equator
    if (abs(lat_in) < lat_minValue) then
      lat = lat_minValue
      lon = lon_in
      call uvr_RotateLatLon( uvr,   & ! INOUT
                             subGridIndex,     & ! IN
                             latRot, lonRot,   & ! OUT (radians)
                             lat, lon,         & ! IN  (radians)
                             'ToLatLonRot')      ! IN
    else
      lat = lat_in
      lon = lon_in
      latRot = latRot_in
      lonRot = lonRot_in
    end if

    coslatr = cos(LatRot)
    sinlatr = sin(LatRot)
    coslonr = cos(LonRot)
    sinlonr = sin(LonRot)
    coslon  = cos(Lon)
    sinlon  = sin(Lon)
    if ( Lat /= 0.0d0 ) then
      rsinlat = 1.0d0/sin(Lat)
    else
      call utl_abort('uvr_rotateWind_tl: cannot be used for points on the equator')
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
      call utl_abort('uvr_rotateWind_tl: mode ToRotWind is not available yet')
    else
      write(*,*) 
      write(*,*) 'uvr_rotateWind_tl: Unknown transform name: ', trim(mode)
      write(*,*) '                mode = ToMetWind or ToRotWind'
      call utl_abort('uvr_rotateWind_tl')
    end if

  end subroutine uvr_rotateWind_tl

  subroutine uvr_rotateWind_ad( uvr, subGridIndex, uwind, vwind, Lat_in, Lon_in, LatRot_in, LonRot_in, mode )
    !
    ! :Purpose: Adjoint of : Go from tangential wind components from one sphere to another
    !           (same origin!). Fast version used by Variational analysis. 
    !
    implicit none

    ! Arguments:
    type(struct_uvr), pointer :: uvr          ! Wind rotation object
    integer, intent(in)       :: subGridIndex ! Current subgrid index
    real(8), intent(inout)    :: uwind        ! interpUU
    real(8), intent(inout)    :: vwind        ! interpVV
    real(8), intent(in)       :: Lat_in       ! Latitude in radians
    real(8), intent(in)       :: Lon_in       ! Longitude in radians
    real(8), intent(in)       :: LatRot_in    ! Rotated latitude in radians
    real(8), intent(in)       :: LonRot_in    ! Rotated longitude in radians 
    character(*), intent(in)  :: mode         ! ToMetWind or ToRotWind

    ! locals
    integer :: index1, index2
    real(8) :: coslatr, sinlatr, coslonr, sinlonr, coslon, sinlon, rsinlat
    real(8) :: lat, lon, latRot, lonRot
    real(8) :: xyz(msize), uvcart(msize)

    if ( .not. uvr%initialized ) then
      write(*,*)
      call utl_abort('uvr_rotateWind_ad: WindRotation module is not initialize')
    endif

    ! Special case for the equator
    if (abs(lat_in) < lat_minValue) then
      lat = lat_minValue
      lon = lon_in
      call uvr_RotateLatLon( uvr,   & ! INOUT
                             subGridIndex,     & ! IN
                             latRot, lonRot,   & ! OUT (radians)
                             lat, lon,         & ! IN  (radians)
                             'ToLatLonRot')      ! IN
    else
      lat = lat_in
      lon = lon_in
      latRot = latRot_in
      lonRot = lonRot_in
    end if

    coslatr = cos(LatRot)
    sinlatr = sin(LatRot)
    coslonr = cos(LonRot)
    sinlonr = sin(LonRot)
    coslon  = cos(Lon)
    sinlon  = sin(Lon)
    if ( Lat /= 0.0d0 ) then
      rsinlat = 1.0d0/sin(Lat)
    else
      call utl_abort('uvr_rotateWind_ad: cannot be used for points on the equator')
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
      call utl_abort('uvr_rotateWind_ad: mode ToRotWind is not available yet')
    else
      write(*,*) 
      write(*,*) 'uvr_rotateWind_ad: Unknown transform name: ', trim(mode)
      write(*,*) '                   mode = ToMetWind or ToRotWind'
      call utl_abort('uvr_rotateWind_ad')
    end if

  end subroutine uvr_rotateWind_ad

  subroutine uvr_rotateLatLon( uvr, subGridIndex, LatOut, LonOut, LatIn, LonIn, mode )
    !
    ! :Purpose: Go from (lat,lon) of one Cartesian frame to (lat,lon) of another
    !           Cartesian frame given the rotation matrix.
    !
    implicit none

    ! Arguments:
    type(struct_uvr), pointer :: uvr          ! Wind rotation object
    integer, intent(in)       :: subGridIndex ! Current subgrid index
    real(8), intent(in)       :: LatIn        ! Input latitude in radians
    real(8), intent(in)       :: LonIn        ! Input longitude in radians
    real(8), intent(out)      :: LatOut       ! Output latitude in radians
    real(8), intent(out)      :: LonOut       ! Output longitude in radians 
    character(*), intent(in)  :: mode         ! ToLatLonRot or ToLatLon

    ! locals
    real(8) :: CartIn(msize),CartOut(msize)
    real(8) :: rLon,rLat

    if ( .not. uvr%initialized ) then
      write(*,*)
      call utl_abort('uvr_rotateLatLon: WindRotation module is not initialize')
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
      write(*,*) 'uvr_rotateLatLon: Unknown transform name: ', trim(mode)
      write(*,*) '                  mode = ToLatLonRot or ToLatLon'
      call utl_abort('uvr_rotateLatLon')
    end if

    call carall( LonOut, LatOut, & ! OUT
                 CartOut )         ! IN

    LonOut = LonOut * MPC_RADIANS_PER_DEGREE_R8 ! To radians
    LatOut = LatOut * MPC_RADIANS_PER_DEGREE_R8 ! To radians

  end subroutine uvr_rotateLatLon


  subroutine carall( plon, plat, pcart )
    !
    ! :Purpose: Returns (lat,lon) (degrees) of an input Cartesian position vector 
    !           on the unit sphere.
    !
    implicit none

    ! Arguments:
    real(8), intent(out) :: plat         ! output latitude 
    real(8), intent(out) :: plon         ! output longitude
    real(8), intent(in)  :: pcart(msize) ! input Cartesian vector

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

  subroutine mxv( pvec2, pmat, pvec1, kdimi, kdimj )
    !
    ! :Purpose: Compute a product : matrix times vector.
    !
    implicit none

    ! Arguments:
    real(8) :: pvec2(kdimi)      ! output vector
    real(8) :: pmat(kdimi,kdimj) ! input matrix
    real(8) :: pvec1(kdimj)      ! input vector
    integer :: kdimi             ! first dimension
    integer :: kdimj             ! second dimension
    

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

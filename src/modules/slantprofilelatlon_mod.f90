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

module slantprofilelatlon_mod
  ! MODULE slantprofilelatlon_mod (prefix='slp' category='7. Low-level data objects and utilities')
  !
  ! :Purpose: To calculate latitudes/longitudes on slant-path based on
  !           ColumnData.
  !
  use mpi_mod
  use earthConstants_mod
  use mathPhysConstants_mod
  use utilities_mod
  use obsSpaceData_mod
  use horizontalCoord_mod

  implicit none
  save
  private

  ! public procedures
  public :: slp_calcLatLonTovs

  ! private module variables and derived types
  real(4), save :: toleranceHeightDiff
  integer, save :: maxNumIteration
  logical, save :: nmlAlreadyRead = .false.


contains 

  subroutine slp_calcLatLonTovs(obsSpaceData, hco, headerIndex, height3D_T_r4, height3D_M_r4, latSlantLev_T, lonSlantLev_T, latSlantLev_M, lonSlantLev_M )
    !
    ! :Purpose: call the computation of lat/lon on the slant path for radiance 
    !           observations, iteratively. To replace the vertical columns with 
    !           line-of-sight slanted columns.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(in) :: obsSpaceData
    type(struct_hco), intent(in) :: hco
    integer, intent(in)  :: headerIndex
    real(4), intent(in)  :: height3D_T_r4(:,:,:)
    real(4), intent(in)  :: height3D_M_r4(:,:,:)
    real(8), intent(out)  :: latSlantLev_T(:)
    real(8), intent(out)  :: lonSlantLev_T(:)
    real(8), intent(out)  :: latSlantLev_M(:)
    real(8), intent(out)  :: lonSlantLev_M(:)

    ! Locals:
    real(4) :: heightInterp, heightIntersect, heightDiff 
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    real(8) :: lat, lon, latSlant, lonSlant
    integer :: subGridIndex, lonIndex, latIndex
    integer :: ierr, fnom, fclos, nulnam
    integer :: nlev_T, lev_T, nlev_M, lev_M
    integer :: numIteration
    logical :: doIteration

    namelist /namSlantPath/ toleranceHeightDiff, maxNumIteration

    if ( .not. nmlAlreadyRead ) then
      nmlAlreadyRead = .true.

      ! default values
      toleranceHeightDiff = 10.0
      maxNumIteration = 1

      ! reading namelist variables
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam, nml=namSlantPath, iostat=ierr)
      if (ierr /= 0) write(*,*) 'slp_calcLatLonTovs: namSlantPath is missing in the namelist. The default value will be taken.'
      if (mpi_myid == 0) write(*, nml=namSlantPath)
      ierr = fclos(nulnam)
    end if

    nlev_M = size(height3D_M_r4,3)
    nlev_T = size(height3D_T_r4,3)

    lat = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)
    lon = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex)
    if (lon <  0.0d0          ) lon = lon + 2.0d0*MPC_PI_R8
    if (lon >= 2.0d0*MPC_PI_R8) lon = lon - 2.0d0*MPC_PI_R8

    ! put the last lat/lon at the surface
    latSlantLev_T(nlev_T) = lat
    lonSlantLev_T(nlev_T) = lon

    ! loop through rest of thermo levels
    do lev_T = 1, nlev_T-1

      ! find the interpolated height 
      call tmg_start(197,'heightBilinearInterp')
      call heightBilinearInterp(lat, lon, hco, height3D_T_r4(:,:,lev_T), heightInterp)
      call tmg_stop(197)

      doIteration = .true.
      numIteration = 0
      while_doIteration: do while (doIteration)

        call tmg_start(196,'findIntersectLatlon')
        call findIntersectLatlon(obsSpaceData, headerIndex, heightInterp, latSlant, lonSlant)
        call tmg_stop(196)

        ! find the interpolated height 
        call tmg_start(197,'heightBilinearInterp')
        call heightBilinearInterp(latSlant, lonSlant, hco, height3D_T_r4(:,:,lev_T), heightIntersect)
        call tmg_stop(197)

        heightDiff = abs(heightInterp-heightIntersect)
        if ( heightDiff > toleranceHeightDiff .or. &
             numIteration >= maxNumIteration ) doIteration = .false.

        numIteration = numIteration + 1

      end do while_doIteration

      latSlantLev_T(lev_T) = latSlant
      lonSlantLev_T(lev_T) = lonSlant
    end do

    ! put the last lat/lon at the surface
    latSlantLev_M(nlev_M) = lat
    lonSlantLev_M(nlev_M) = lon

    ! loop through rest of momentum levels
    do lev_M = 1, nlev_M-1

      ! find the interpolated height 
      call tmg_start(197,'heightBilinearInterp')
      call heightBilinearInterp(lat, lon, hco, height3D_M_r4(:,:,lev_M), heightInterp)
      call tmg_stop(197)

      doIteration = .true.
      numIteration = 0
      while_doIteration2: do while (doIteration)

        call tmg_start(196,'findIntersectLatlon')
        call findIntersectLatlon(obsSpaceData, headerIndex, heightInterp, latSlant, lonSlant)
        call tmg_stop(196)

        ! find the interpolated height 
        call tmg_start(197,'heightBilinearInterp')
        call heightBilinearInterp(latSlant, lonSlant, hco, height3D_M_r4(:,:,lev_M), heightIntersect)
        call tmg_stop(197)

        heightDiff = abs(heightInterp-heightIntersect)
        if ( heightDiff > toleranceHeightDiff .or. &
             numIteration >= maxNumIteration ) doIteration = .false.

        numIteration = numIteration + 1

      end do while_doIteration2

      latSlantLev_M(lev_M) = latSlant
      lonSlantLev_M(lev_M) = lonSlant
    end do

  end subroutine slp_calcLatLonTovs


  subroutine findIntersectLatlon(obsSpaceData, headerIndex, height_r4, latSlant, lonSlant)
    !
    !:Purpose: Computation of lat/lon of the intersection between model level 
    !          and the slant line-of-sight for radiance observations.
    !
    implicit none
    ! Arguments:
    type(struct_obs), intent(in)  :: obsSpaceData
    real(4), intent(in)  :: height_r4
    integer, intent(in)  :: headerIndex
    real(8), intent(out) :: latSlant
    real(8), intent(out) :: lonSlant 

    ! Locals:
    real(8) :: lat, lon, geometricHeight
    real(8) :: zenithAngle, zenithAngle_rad, azimuthAngle, azimuthAngle_rad, elevationAngle_rad, distAlongPath
    real(8) :: obsCordGlb(3), slantPathCordGlb(3), unitx(3), unity(3), unitz(3), unitSatLoc(3), unitSatGlb(3)

    ! read lat/lon/angles from obsSpaceData
    lat = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)
    lon = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex)
    if (lon <  0.0d0          ) lon = lon + 2.0d0*MPC_PI_R8
    if (lon >= 2.0d0*MPC_PI_R8) lon = lon - 2.0d0*MPC_PI_R8
    azimuthAngle = obs_headElem_r(obsSpaceData,OBS_AZA,headerIndex)
    zenithAngle = obs_headElem_r(obsSpaceData,OBS_SZA,headerIndex)

    ! convert angles to radian unit
    azimuthAngle_rad = azimuthAngle * MPC_RADIANS_PER_DEGREE_R8
    zenithAngle_rad = zenithAngle * MPC_RADIANS_PER_DEGREE_R8
    elevationAngle_rad = 0.5d0 * MPC_PI_R8 - zenithAngle_rad

    obsCordGlb  = RA * (/ cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat) /)
    unitz = (/  cos(lat)*cos(lon), cos(lat)*sin(lon) , sin(lat) /)
    unitx = (/ -sin(lon)         , cos(lon)          , 0.d0     /)
    unity = (/ -sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat) /)

    ! unit vector towards satellite in local coordinate
    unitSatLoc = (/ cos(elevationAngle_rad)*sin(azimuthAngle_rad) , &
                    cos(elevationAngle_rad)*cos(azimuthAngle_rad) , &
                    sin(elevationAngle_rad) /)
    ! unit vector towards satellite in global coordinate
    unitSatGlb = unitSatLoc(1) * unitx + unitSatLoc(2) * unity + unitSatLoc(3) * unitz

    ! Geometric altitude
    geometricHeight = real(height_r4,8)

    ! distance along line of sight
    distAlongPath = geometricHeight / cos(zenithAngle_rad)

    slantPathCordGlb(:) = obsCordGlb(:) + distAlongPath * unitSatGlb(:) 

    latSlant = atan(slantPathCordGlb(3)/sqrt(slantPathCordGlb(1)**2+slantPathCordGlb(2)**2))
    lonSlant = atan2(slantPathCordGlb(2),slantPathCordGlb(1))
    if (lonSlant <  0.0d0          ) lonSlant = lonSlant + 2.0d0*MPC_PI_R8
    if (lonSlant >= 2.0d0*MPC_PI_R8) lonSlant = lonSlant - 2.0d0*MPC_PI_R8

  end subroutine findIntersectLatlon


  subroutine heightBilinearInterp(lat, lon, hco, height_r4, heightInterp_r4)
    !
    !:Purpose: To interpolate the 2D height field to the obs location
    !
    implicit none

    ! Arguments:
    real(8), intent(in)          :: lat
    real(8), intent(in)          :: lon
    type(struct_hco), intent(in) :: hco
    real(4), intent(in)          :: height_r4(:,:)
    real(4), intent(out)         :: heightInterp_r4

    ! Locals:
    integer :: ierr, niP1
    integer :: latIndex, lonIndex, latIndex2, lonIndex2, lonIndexP1
    integer :: subGridIndex, subGridForInterp, numSubGridsForInterp
    integer :: ipoint, gridptCount
    integer :: latIndexVec(8), lonIndexVec(8)
    real(8) :: WeightVec(8)
    real(8) :: dldx, dldy
    real(8) :: weightsSum
    real(4) :: lon_deg_r4, lat_deg_r4
    real(4) :: xpos_r4, ypos_r4, xpos2_r4, ypos2_r4
    logical :: latlonOutsideGrid

    lat_deg_r4 = real(lat * MPC_DEGREES_PER_RADIAN_R8)
    lon_deg_r4 = real(lon * MPC_DEGREES_PER_RADIAN_R8)
    ierr = utl_getPositionXY( hco%EZscintID,   &
                              xpos_r4, ypos_r4, xpos2_r4, ypos2_r4, &
                              lat_deg_r4, lon_deg_r4, subGridIndex )

    ! Allow for periodicity in Longitude for global Gaussian grid
    if ( hco%grtyp == 'G' .or. (hco%grtyp == 'Z' .and. hco%global) ) then
      niP1 = hco%ni + 1
    else
      niP1 = hco%ni
    end if

    ! Check if the interpolated lat/lon is outside the hco domain
    latlonOutsideGrid = ( xpos_r4 < 1.0 .or. &
                          xpos_r4 > real(niP1) .or. &
                          ypos_r4 < 1.0 .or. &
                          ypos_r4 > real(hco%nj) )

    if ( latlonOutsideGrid ) then
      write(*,*) 'heightBilinearInterp: interpolated lat/lon outside the hco domain.'
      write(*,*) '  position   lon,       lat = ', lon_deg_r4, lat_deg_r4
      write(*,*) '  position     x,       y   = ', xpos_r4, ypos_r4

      ! if above or below domain
      if( ypos_r4 < 1.0 ) ypos_r4 = 1.0
      if( ypos_r4 > real(hco%nj) ) ypos_r4 = real(hco%nj)

      ! if on the left or right longitude band, move it to the edge of this longitude band
      if( xpos_r4 < 1.0 ) xpos_r4 = 1.0
      if( xpos_r4 > real(niP1) ) xpos_r4 = real(niP1)
      write(*,*) '  new position x,       y   = ', xpos_r4, ypos_r4

    end if

    ! Find the lower-left grid point next to the observation
    if ( xpos_r4 == real(niP1) ) then
      lonIndex = floor(xpos_r4) - 1
    else
      lonIndex = floor(xpos_r4)
    end if
    if ( xpos2_r4 == real(niP1) ) then
      lonIndex2 = floor(xpos2_r4) - 1
    else
      lonIndex2 = floor(xpos2_r4)
    end if

    if ( ypos_r4 == real(hco%nj) ) then
      latIndex = floor(ypos_r4) - 1
    else
      latIndex = floor(ypos_r4)
    end if
    if ( ypos2_r4 == real(hco%nj) ) then
      latIndex2 = floor(ypos2_r4) - 1
    else
      latIndex2 = floor(ypos2_r4)
    end if

    if ( hco%grtyp == 'U' ) then
      if ( ypos_r4 == real(hco%nj/2) ) then
        latIndex = floor(ypos_r4) - 1
      else
        latIndex = floor(ypos_r4)
      end if
      if ( ypos2_r4 == real(hco%nj/2) ) then
        latIndex2 = floor(ypos2_r4) - 1
      else
        latIndex2 = floor(ypos2_r4)
      end if
    end if

    ! Handle periodicity in longitude
    lonIndexP1 = lonIndex + 1
    if ( lonIndexP1 == hco%ni + 1 ) lonIndexP1 = 1

    ! Check if location is in between Yin and Yang (should not happen)
    if ( hco%grtyp == 'U' ) then
      if ( ypos_r4 > real(hco%nj/2) .and.  &
           ypos_r4 < real((hco%nj/2)+1) ) then
        write(*,*) 'heightBilinearInterp: WARNING, obs position in between Yin and Yang!!!'
        write(*,*) '   xpos, ypos = ', xpos_r4, ypos_r4
      end if
      if ( ypos2_r4 > real(hco%nj/2) .and.  &
           ypos2_r4 < real((hco%nj/2)+1) ) then
        write(*,*) 'heightBilinearInterp: WARNING, obs position in between Yin and Yang!!!'
        write(*,*) '   xpos2, ypos2 = ', xpos2_r4, ypos2_r4
      end if
    end if

    if ( subGridIndex == 3 ) then
      ! both subGrids involved in interpolation, so first treat subGrid 1
      numSubGridsForInterp = 2
      subGridIndex = 1
    else
      ! only 1 subGrid involved in interpolation
      numSubGridsForInterp = 1
    end if

    gridptCount = 0

    do subGridForInterp = 1, numSubGridsForInterp

      ! Compute the 4 weights of the bilinear interpolation
      if ( subGridForInterp == 1 ) then
        ! when only 1 subGrid involved, subGridIndex can be 1 or 2
        dldx = real(xpos_r4,8) - real(lonIndex,8)
        dldy = real(ypos_r4,8) - real(latIndex,8)
      else
        ! when 2 subGrids, subGridIndex is set to 1 for 1st iteration, 2 for second
        subGridIndex = 2
        lonIndex = lonIndex2
        latIndex = latIndex2
        lonIndexP1 = lonIndex2 + 1
        dldx = real(xpos2_r4,8) - real(lonIndex,8)
        dldy = real(ypos2_r4,8) - real(latIndex,8)
      end if

      gridptCount = gridptCount + 1
      latIndexVec(gridptCount) = latIndex
      lonIndexVec(gridptCount) = lonIndex
      WeightVec(gridptCount) = (1.d0-dldx) * (1.d0-dldy)

      gridptCount = gridptCount + 1
      latIndexVec(gridptCount) = latIndex
      lonIndexVec(gridptCount) = lonIndexP1
      WeightVec(gridptCount) =       dldx  * (1.d0-dldy)

      gridptCount = gridptCount + 1
      latIndexVec(gridptCount) = latIndex + 1
      lonIndexVec(gridptCount) = lonIndex
      WeightVec(gridptCount) = (1.d0-dldx) *       dldy

      gridptCount = gridptCount + 1
      latIndexVec(gridptCount) = latIndex + 1
      lonIndexVec(gridptCount) = lonIndexP1
      WeightVec(gridptCount) =       dldx  *       dldy

    end do ! subGrid

    weightsSum = sum(WeightVec(1:gridptCount))
    if ( weightsSum > 0.d0 ) then

      WeightVec(1:gridptCount) = WeightVec(1:gridptCount) / weightsSum

    else

      call utl_abort('heightBilinearInterp: weightsSum smaller than 0.')

    end if

    ! perform the bi-linear interpolation
    heightInterp_r4 = 0.0
    do ipoint = 1, gridptCount

      if ( latIndexVec(ipoint) > hco%nj ) then
        write(*,*) 'heightBilinearInterp: latIndexVec(ipoint) > hco%nj: latIndex=',latIndexVec(ipoint),', lonIndex=',lonIndexVec(ipoint)
        write(*,*) 'lat_deg_r4=',lat_deg_r4,', lon_deg_r4=',lon_deg_r4
        write(*,*) 'ypos_r4=',ypos_r4,', xpos_r4=',xpos_r4
      end if
      if ( latIndexVec(ipoint) < 1 ) then
        write(*,*) 'heightBilinearInterp: latIndexVec(ipoint) < 1: latIndex=',latIndexVec(ipoint),', lonIndex=',lonIndexVec(ipoint)
        write(*,*) 'lat_deg_r4=',lat_deg_r4,', lon_deg_r4=',lon_deg_r4
        write(*,*) 'ypos_r4=',ypos_r4,', xpos_r4=',xpos_r4
      end if

      if ( lonIndexVec(ipoint) > hco%ni ) then
        write(*,*) 'heightBilinearInterp: lonIndexVec(ipoint) > hco%ni: latIndex=',latIndexVec(ipoint),', lonIndex=',lonIndexVec(ipoint)
        write(*,*) 'lat_deg_r4=',lat_deg_r4,', lon_deg_r4=',lon_deg_r4
        write(*,*) 'ypos_r4=',ypos_r4,', xpos_r4=',xpos_r4
      end if
      if ( lonIndexVec(ipoint) < 1 ) then
        write(*,*) 'heightBilinearInterp: lonIndexVec(ipoint) < 1: latIndex=',latIndexVec(ipoint),', lonIndex=',lonIndexVec(ipoint)
        write(*,*) 'lat_deg_r4=',lat_deg_r4,', lon_deg_r4=',lon_deg_r4
        write(*,*) 'ypos_r4=',ypos_r4,', xpos_r4=',xpos_r4
      end if

      latlonOutsideGrid = ( latIndexVec(ipoint) < 1 .or. &
                            latIndexVec(ipoint) > hco%nj .or. &
                            lonIndexVec(ipoint) < 1 .or. &
                            lonIndexVec(ipoint) > hco%ni )

      if ( latlonOutsideGrid ) &
        call utl_abort('heightBilinearInterp: lat/lon outside the domain.')

      heightInterp_r4 = heightInterp_r4 + &
                    real(WeightVec(ipoint),4) * &
                    height_r4(lonIndexVec(ipoint), latIndexVec(ipoint))
    end do

  end subroutine heightBilinearInterp

end module slantprofilelatlon_mod

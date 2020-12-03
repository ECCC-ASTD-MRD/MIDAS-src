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
  use tovs_nl_mod
  use codtyp_mod
  use getGridPosition_mod

  implicit none
  save
  private

  ! public procedures
  public :: slp_calcLatLonTovs, slp_calcLatLonRO, slp_calcLatLonRadar, slp_radar_getlatlonHRfromRange, slp_radar_getRangefromH, slp_radar_getHfromRange

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
    real(4) :: heightInterp_r4, heightIntersect_r4, heightDiff_r4 
    real(8) :: lat, lon, latSlant, lonSlant, azimuthAngle
    integer :: idatyp
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
      read(nulnam, nml = namSlantPath, iostat = ierr)
      if (ierr /= 0) write(*,*) 'slp_calcLatLonTovs: namSlantPath is missing in the namelist. The default value will be taken.'
      if (mpi_myid == 0) write(*, nml = namSlantPath)
      ierr = fclos(nulnam)
    end if

    nlev_M = size(height3D_M_r4,3)
    nlev_T = size(height3D_T_r4,3)

    lat = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)
    lon = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex)
    if (lon <  0.0d0          ) lon = lon + 2.0d0*MPC_PI_R8
    if (lon >= 2.0d0*MPC_PI_R8) lon = lon - 2.0d0*MPC_PI_R8
    azimuthAngle = tvs_getCorrectedSatelliteAzimuth(obsSpaceData, headerIndex)

    !SSMIS special case
    idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
    if ( idatyp == codtyp_get_codtyp('ssmis') ) then
      azimuthAngle = obs_missingValue_R
    end if

    ! loop through all thermo levels
    do lev_T = 1, nlev_T

      ! find the interpolated height 
      call tmg_start(197,'heightBilinearInterp')
      call heightBilinearInterp(lat, lon, hco, height3D_T_r4(:,:,lev_T), heightInterp_r4)
      call tmg_stop(197)

      doIteration = .true.
      numIteration = 0
      while_doIteration: do while (doIteration)

        numIteration = numIteration + 1

        call tmg_start(196,'findIntersectLatlon')
        call findIntersectLatlon(obsSpaceData, headerIndex, heightInterp_r4, azimuthAngle, latSlant, lonSlant)
        call tmg_stop(196)

        ! find the interpolated height 
        call tmg_start(197,'heightBilinearInterp')
        call heightBilinearInterp(latSlant, lonSlant, hco, height3D_T_r4(:,:,lev_T), heightIntersect_r4)
        call tmg_stop(197)

        heightDiff_r4 = abs(heightInterp_r4-heightIntersect_r4)
        if ( heightDiff_r4 < toleranceHeightDiff .or. &
             numIteration >= maxNumIteration ) doIteration = .false.

      end do while_doIteration

      latSlantLev_T(lev_T) = latSlant
      lonSlantLev_T(lev_T) = lonSlant
    end do

    ! loop through all momentum levels
    do lev_M = 1, nlev_M

      ! find the interpolated height 
      call tmg_start(197,'heightBilinearInterp')
      call heightBilinearInterp(lat, lon, hco, height3D_M_r4(:,:,lev_M), heightInterp_r4)
      call tmg_stop(197)

      doIteration = .true.
      numIteration = 0
      while_doIteration2: do while (doIteration)

        numIteration = numIteration + 1

        call tmg_start(196,'findIntersectLatlon')
        call findIntersectLatlon(obsSpaceData, headerIndex, heightInterp_r4, azimuthAngle, latSlant, lonSlant)
        call tmg_stop(196)

        ! find the interpolated height 
        call tmg_start(197,'heightBilinearInterp')
        call heightBilinearInterp(latSlant, lonSlant, hco, height3D_M_r4(:,:,lev_M), heightIntersect_r4)
        call tmg_stop(197)

        heightDiff_r4 = abs(heightInterp_r4-heightIntersect_r4)
        if ( heightDiff_r4 < toleranceHeightDiff .or. &
             numIteration >= maxNumIteration ) doIteration = .false.

      end do while_doIteration2

      latSlantLev_M(lev_M) = latSlant
      lonSlantLev_M(lev_M) = lonSlant
    end do

  end subroutine slp_calcLatLonTovs


  subroutine findIntersectLatlon(obsSpaceData, headerIndex, height_r4, azimuthAngle, latSlant, lonSlant)
    !
    !:Purpose: Computation of lat/lon of the intersection between model level 
    !          and the slant line-of-sight for radiance observations.
    !
    implicit none
    ! Arguments:
    type(struct_obs), intent(in)  :: obsSpaceData
    real(4), intent(in)  :: height_r4
    integer, intent(in)  :: headerIndex
    real(8), intent(in)  :: azimuthAngle
    real(8), intent(out) :: latSlant
    real(8), intent(out) :: lonSlant 

    ! Locals:
    real(8) :: lat, lon, geometricHeight
    real(8) :: zenithAngle, zenithAngle_rad, azimuthAngle_rad, elevationAngle_rad, distAlongPath
    real(8) :: obsCordGlb(3), slantPathCordGlb(3), unitx(3), unity(3), unitz(3), unitSatLoc(3), unitSatGlb(3)

    ! read lat/lon/angles from obsSpaceData
    lat = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)
    lon = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex)
    if (lon <  0.0d0          ) lon = lon + 2.0d0*MPC_PI_R8
    if (lon >= 2.0d0*MPC_PI_R8) lon = lon - 2.0d0*MPC_PI_R8

    if ( azimuthAngle /= obs_missingValue_R) then

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

    else

      latSlant = lat
      lonSlant = lon

    end if

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
    ierr = gpos_getPositionXY( hco%EZscintID,   &
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

  subroutine slp_calcLatLonRO(obsSpaceData, hco, headerIndex, &
                              height3D_T_r4, height3D_M_r4, &
                              latSlantLev_T, lonSlantLev_T, latSlantLev_M, lonSlantLev_M )
    !
    ! :Purpose: call the computation of lat/lon on the slant path for GPSRO
    !           observations, iteratively. To replace the vertical columns with 
    !           slanted columns.
    !
    implicit none

    ! Arguments:
    type(struct_obs)     :: obsSpaceData
    type(struct_hco)     :: hco 
    integer, intent(in)  :: headerIndex
    real(4), intent(in)  :: height3D_T_r4(:,:,:)
    real(4), intent(in)  :: height3D_M_r4(:,:,:)
    real(8), intent(out) :: latSlantLev_T(:)
    real(8), intent(out) :: lonSlantLev_T(:)
    real(8), intent(out) :: latSlantLev_M(:)
    real(8), intent(out) :: lonSlantLev_M(:)

    ! Locals:
    real(8) :: latr, lonr, height, rad, dH, hmin
    integer :: bodyIndex, imin
    real(8), allocatable :: Lat_Obs(:), Lon_Obs(:), Hgt_Obs(:)
    real(4), allocatable :: H_M(:,:), H_T(:,:)
    integer :: nObs, iObs
    integer :: levIndex, nlev_M, nlev_T
    
    nlev_M = size(height3D_M_r4,3)
    nlev_T = size(height3D_T_r4,3)

    ! Header, and Obs counting:
    call obs_set_current_header_list(obsSpaceData,'RO')
    call obs_set_current_body_list(obsSpaceData, headerIndex)
    Rad = obs_headElem_r(obsSpaceData,OBS_TRAD,headerIndex)
    nObs = 0
    BODY: do
      bodyIndex = obs_getBodyIndex(obsSpaceData)
      if (bodyIndex < 0) exit BODY
      nObs = nObs + 1
    end do BODY

    call obs_set_current_body_list(obsSpaceData, headerIndex)

    allocate(Hgt_Obs(nObs), Lat_Obs(nObs), Lon_Obs(nObs))
    allocate(H_M(nObs,nlev_M), H_T(nObs,nlev_T))

    do iObs = 1,nObs
      bodyIndex = obs_getBodyIndex(obsSpaceData)
      if (bodyIndex < 0) exit

      height = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
      ! If the vertical coordinate is an impact parameter (6e6<himp<7e6), subtract radius:
      if (6.e6 < height .and. height < 7.e6) height = height-rad
      Hgt_Obs(iObs) = height
      latr = obs_bodyElem_r(obsSpaceData,OBS_ROLA,bodyIndex)
      lonr = obs_bodyElem_r(obsSpaceData,OBS_ROLO,bodyIndex)
      if (lonr <  0.d0          ) lonr = lonr + 2.d0*MPC_PI_R8
      if (lonr >= 2.d0*MPC_PI_R8) lonr = lonr - 2.d0*MPC_PI_R8
      Lat_Obs(iObs) = latr
      Lon_Obs(iObs) = lonr
      do levIndex = 1, nlev_M
        call heightBilinearInterp(latr, lonr, hco, height3D_M_r4(:,:,levIndex), &
                                  H_M(iObs,levIndex))
      end do
      do levIndex = 1, nlev_T
        call heightBilinearInterp(latr, lonr, hco, height3D_T_r4(:,:,levIndex), &
                                  H_T(iObs,levIndex))
      end do
    end do

    do levIndex = 1, nlev_M
      hmin = 1.e30
      imin = -1
      do iObs = 1, nObs
        dH = abs(H_M(iObs,levIndex)-Hgt_Obs(iObs))
        if (dH < hmin) then
          hmin = dH
          imin = iObs
        end if
      end do
      latSlantLev_M(levIndex) = Lat_Obs(imin)
      lonSlantLev_M(levIndex) = Lon_Obs(imin)
    end do

    do levIndex = 1, nlev_T
      hmin = 1.e30
      imin = -1
      do iObs = 1, nObs
        dH = abs(H_T(iObs,levIndex)-Hgt_Obs(iObs))
        if (dH < hmin) then
          hmin = dH
          imin = iObs
        end if
      end do
      latSlantLev_T(levIndex) = Lat_Obs(imin)
      lonSlantLev_T(levIndex) = Lon_Obs(imin)
    end do

    deallocate(H_T, H_M)
    deallocate(Lon_Obs, Lat_Obs, Hgt_Obs)

  end subroutine slp_calcLatLonRO
 
  subroutine slp_calcLatLonRadar(obsSpaceData, hco, headerIndex, height3D_T_r4, height3D_M_r4, latSlantLev_T, lonSlantLev_T, latSlantLev_M, lonSlantLev_M)
    !
    ! :Purpose: call the computation of lat/lon on the slant path for radar 
    !           observations, iteratively. To replace the vertical columns with 
    !           radar beam columns .
    !
    implicit none
    ! Arguments:
    type(struct_obs), intent(in) :: obsSpaceData
    type(struct_hco), intent(in) :: hco 
    integer, intent(in)  :: headerIndex
    real(4), intent(in)  :: height3D_T_r4(:,:,:)
    real(4), intent(in)  :: height3D_M_r4(:,:,:)
    real(8), intent(out) :: latSlantLev_T(:)
    real(8), intent(out) :: lonSlantLev_T(:)
    real(8), intent(out) :: latSlantLev_M(:)
    real(8), intent(out) :: lonSlantLev_M(:)
    ! Locals:
    real(8) :: lat, lon, latSlant, lonSlant, rele, rzam, ralt, rans, rane
    integer :: nlev_M,lev_M, nlev_T,lev_T

    nlev_M = size(height3D_M_r4,3)
    nlev_T = size(height3D_T_r4,3)

    lat  = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)
    lon  = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex)
    if (lon <  0.0d0          ) lon = lon + 2.0d0*MPC_PI_R8
    if (lon >= 2.0d0*MPC_PI_R8) lon = lon - 2.0d0*MPC_PI_R8

    rele = obs_headElem_r(obsSpaceData, OBS_RELE, headerIndex) *  MPC_RADIANS_PER_DEGREE_R8
    rzam = obs_headElem_r(obsSpaceData, OBS_RZAM, headerIndex) *  MPC_RADIANS_PER_DEGREE_R8
    ralt = obs_headElem_r(obsSpaceData, OBS_ALT,  headerIndex)
    rans = obs_headElem_r(obsSpaceData, OBS_RANS, headerIndex)
    rane = obs_headElem_r(obsSpaceData, OBS_RANE, headerIndex)
    
    ! put the last lat/lon at the surface of thermo level
    latSlantLev_T(nlev_T) = lat
    lonSlantLev_T(nlev_T) = lon
    ! Loop through rest of thermo level
    do lev_T = 1, nlev_T-1
      call findIntersectLatlonRadar(lat, lon, rele, rzam,  ralt, rans, rane, hco, height3D_T_r4(:,:,lev_T), latSlant, lonSlant)
      latSlantLev_T(lev_T) = latSlant
      lonSlantLev_T(lev_T) = lonSlant
    end do
    
    ! put the last lat/lon at the surface of momentum levels
    latSlantLev_M(nlev_M) = lat
    lonSlantLev_M(nlev_M) = lon
    ! loop through rest of momentum levels
    do lev_M = 1, nlev_M-1
      ! find the intersection 
      call findIntersectLatlonRadar(lat, lon, rele, rzam,  ralt, rans, rane, hco, height3D_M_r4(:,:,lev_M), latSlant, lonSlant)
      latSlantLev_M(lev_M) = latSlant
      lonSlantLev_M(lev_M) = lonSlant
    end do

  end subroutine slp_calcLatLonRadar

  subroutine findIntersectLatlonRadar(lat, lon, rele, rzam,  ralt, rans, rane,  hco, field2d_height, latSlant, lonSlant)
    !
    ! :Purpose: Computation of lat lon  of the intersection 
    !             of radar beam with the levels model
    !
    ! :NOTE: Bisection method  
    !
    implicit none
    !Argument
    type(struct_hco), intent(in) :: hco
    real(8), intent(in)          :: lat, lon, rele, rzam, ralt, rans, rane
    real(4), intent(in)          :: field2d_height(hco%ni,hco%nj)
    real(8)                      :: upper_bound, lower_bound, tolerance
    real(8), intent(out)         :: latSlant,lonSlant
    ! Local
    real(8)                      :: h_radar, mid, difference_heights, d_radar
    integer                      :: iteration,maximum_iteration
    real(4)                      :: heightInterp_r4

    lower_bound  = rans ! m max useable range for weather radars
    upper_bound  = rane ! m minimum range
    tolerance = 1       ! tolerance of the bisection method 
    iteration = 0
    maximum_iteration = 30
    ! bisection method 
    mid=0.
    do while ((abs((upper_bound - lower_bound)/2.)>tolerance).or.(iteration>maximum_iteration))
      mid = (upper_bound + lower_bound)/2
      call slp_radar_getlatlonHRfromRange(lat, lon, rele, rzam, ralt, mid, latSlant, lonSlant, h_radar, d_radar)
      call heightBilinearInterp(latSlant, lonSlant, hco, field2d_height, heightInterp_r4)
      difference_heights = h_radar-heightInterp_r4
      if (difference_heights>0.) then 
        upper_bound = mid
      else
        lower_bound = mid
      end if 
      iteration = iteration+1
    end do
   
    if (lonSlant <  0.0d0          ) lonSlant = lonSlant + 2.0d0*MPC_PI_R8
    if (lonSlant >= 2.0d0*MPC_PI_R8) lonSlant = lonSlant - 2.0d0*MPC_PI_R8
  
  end subroutine findIntersectLatlonRadar
 
  subroutine slp_radar_getlatlonHRfromRange(lat, lon, rele, rzam, ralt, r_radar, latSlant,lonSlant, h_radar,d_radar)
    !
    ! :Purpose: Computation of  lat-lon , height  of the trajectory
    !            along the radar beam from range of the radar beam
    !
    implicit none
    ! Argument 
    real(8), intent(in)  :: lat, lon, rele, rzam, ralt
    real(8), intent(in)  :: r_radar
    real(8), intent(out) :: latSlant, lonSlant, h_radar,d_radar
    ! Local
    real(8)              :: Re 
    
    ! Radius of sphere of equal area from earthconstants_mod.f90
    ! earth_r2 = 6371007.1809
    ! efective radius of the earth
    Re = earth_r2 * (4./3.)
    ! height radar beam  from range at rele and ralt 
    h_radar = sqrt(r_radar**2.+(earth_r2+ralt)**2.+(2.*r_radar*(earth_r2+ralt)*sin(rele)))-(earth_r2) 
    ! distance radar beam  from range at rele and ralt 
    d_radar = atan(r_radar*cos(rele)/(r_radar*sin(rele)+Re+ralt))*Re
    ! second version -Book:Doppler Radar and Wheather Observation (Doviak and Zrnic)
    !             d_radar = asin(r_range*cos(ele)/(Radius+alt))*Re    
    ! lat lon of the path along the radar beam  
    latSlant = asin( sin(lat)*cos(d_radar/earth_r2) + cos(lat)*sin(d_radar/earth_r2)*cos(rzam))
    lonSlant = lon + atan2(sin(rzam)*sin(d_radar/earth_r2)*cos(lat), cos(d_radar/earth_r2)-sin(lat)*sin(latSlant))

  end subroutine slp_radar_getlatlonHRfromRange

  subroutine slp_radar_getHfromRange(r_radar, ralt, rele, h_radar)
    !
    ! :Purpose: Computation of height of the radar beam
    !            from range of the radar beam
    !
    implicit none
    ! Argument 
    real(8) , intent(in)  :: r_radar, ralt, rele
    real(8) , intent(out) :: h_radar
    ! Local
    real(8)               :: Re 
   
    ! Radius of sphere of equal area from earthconstants_mod.f90
    ! earth_r2 = 6371007.1809
    ! efective radius of the earth
    Re = earth_r2*(4./3.)
    ! height of radar beam  from range at rele and ralt 
    h_radar = sqrt(r_radar**2.+(earth_r2+ralt)**2.+2.*r_radar*(earth_r2+ralt)*sin(rele))-(earth_r2)

  end subroutine slp_radar_getHfromRange

  subroutine slp_radar_getRangefromH(h_radar, ralt, rele, r_radar)
    !
    ! :Purpose: Computation of range of the radar beam from height of the radar beam
    !
    implicit none
    ! Argument
    real(8) , intent(in)  :: h_radar,  ralt, rele
    real(8) , intent(out) :: r_radar 
    ! Local
    real(8)               :: a, b, c 
  
    ! Radius of sphere of equal area from earthconstants_mod.f90
    ! earth_r2 = 6371007.1809
    a = 1.
    b = 2.*(earth_r2+ralt)*sin(rele)
    c = -(earth_r2+h_radar)**2+(earth_r2+ralt)**2.
    ! range of radar beam  from height and elevation of the radar beam 
    r_radar  = (-b+sqrt(b**2-4*a*c))/2*a
     
  end subroutine slp_radar_getRangefromH

end module slantprofilelatlon_mod


module getGridPosition_mod
  ! MODULE getGridPosition_mod (prefix='gpos' category='8. Low-level utilities and constants')
  !
  !:Purpose:  A place to collect numerous interpolation related routines.
  !           The main task of the module is to compute the grid XY position from a lat-lon.
  !           This simply calls the ezsint routine gdxyfll for simple grids. For
  !           Yin-Yan grids it calls the function gpos_xyfll_yinYanGrid
  !           (see below in this module). There is also support for
  !           RPN Y grids, in which case it calls the subroutine gpos_xyfll_unstructGrid.
  !
  use kdTree2_mod
  use mathPhysConstants_mod
  use message_mod
  use physicsFunctions_mod
  use utilities_mod

  implicit none
  save
  private

  ! Public procedures
  public :: gpos_getPositionXY, gpos_gridIsOrca

  interface gpos_getPositionXY
    module procedure gpos_getPosXY_vector
    module procedure gpos_getPosXY_scalar
  end interface gpos_getPositionXY

  integer, parameter :: maxNumLocalGridPointsSearch = 3000

  type(kdtree2), pointer :: tree => null()

contains
  !---------------------------------------------------------
  ! gpos_getPosXY_vector
  !---------------------------------------------------------
  function gpos_getPosXY_vector(gdid, xpos_r4, ypos_r4, xpos2_r4, ypos2_r4,  &
                                lat_deg_r4, lon_deg_r4, subGridIndex) result(ierr)
    !
    ! :Purpose: Compute the grid XY position from a lat-lon. This
    !           simply calls the ezsint routine gdxyfll for simple grids. For
    !           Yin-Yan grids it calls the function gpos_xyfll_yinYanGrid
    !           (see below in this module). There is also support for
    !           RPN Y grids, in which case it calls the subroutine gpos_xyfll_unstructGrid.
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: gdid
    real(4), intent(out) :: xpos_r4(:)
    real(4), intent(out) :: ypos_r4(:)
    real(4), intent(out) :: xpos2_r4(:)
    real(4), intent(out) :: ypos2_r4(:)
    real(4), intent(in)  :: lat_deg_r4(:)
    real(4), intent(in)  :: lon_deg_r4(:)
    integer, intent(out) :: subGridIndex(:)
    ! Result:
    integer :: ierr  ! returned value of function

    ! Locals:
    integer :: numSubGrids, numPoints, pointIndex
    integer :: ezget_nsubGrids, gdxyfll, ezgprm
    character(len=1) :: grtyp
    integer :: ni, nj, ig1, ig2, ig3, ig4

    numPoints = size(xpos_r4)
    numSubGrids = ezget_nsubGrids(gdid)
    xpos2_r4 = mpc_missingValue_R4
    ypos2_r4 = mpc_missingValue_R4

    if (numSubGrids == 1) then

      ierr = ezgprm(gdid, grtyp, ni, nj, ig1, ig2, ig3, ig4)

      if (grtyp == 'Y') then

        !$omp critical
        do pointIndex = 1, numPoints
          ierr = gpos_xyfll_unstructGrid(gdid, xpos_r4(pointIndex), ypos_r4(pointIndex),  &
                                         lat_deg_r4(pointIndex), lon_deg_r4(pointIndex))
        end do
        !$omp end critical

      else

        ! Not an unstructured grid, nor a Yin-Yan grid, call the standard ezscint routine
        ierr = gdxyfll(gdid, xpos_r4, ypos_r4, lat_deg_r4, lon_deg_r4, numPoints)

      end if

      subGridIndex = 1

    else

      ! This is a Yin-Yan grid, do something different
      !$omp critical
      ierr = gpos_xyfll_yinYanGrid(gdid, xpos_r4, ypos_r4, &
                                   xpos2_r4, ypos2_r4, &
                                   lat_deg_r4, lon_deg_r4, subGridIndex)
      !$omp end critical

    end if

    if ( all(subGridIndex(:) /= 3) ) then
      ! when only returning 1 position, copy values to pos2
      xpos2_r4(:) = xpos_r4(:)
      ypos2_r4(:) = ypos_r4(:)
    end if

  end function gpos_getPosXY_vector

  !---------------------------------------------------------
  ! gpos_getPosXY_scalar
  !---------------------------------------------------------
  function gpos_getPosXY_scalar(gdid, xpos_r4, ypos_r4, xpos2_r4, ypos2_r4,  &
                                lat_deg_r4, lon_deg_r4, subGridIndex, maxGridSpacing_opt) result(ierr)
    !
    ! :Purpose: Compute the grid XY position from a lat-lon. This
    !           simply calls the ezsint routine gdxyfll for simple grids. For
    !           Yin-Yan grids it calls the function gpos_xyfll_yinYanGrid
    !           (see below in this module). There is also support for
    !           RPN Y grids, in which case it calls the subroutine gpos_xyfll_unstructGrid.
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: gdid
    real(4), intent(out) :: xpos_r4
    real(4), intent(out) :: ypos_r4
    real(4), intent(out) :: xpos2_r4
    real(4), intent(out) :: ypos2_r4
    real(4), intent(in)  :: lat_deg_r4
    real(4), intent(in)  :: lon_deg_r4
    integer, intent(out) :: subGridIndex
    real(8), optional, intent(in) :: maxGridSpacing_opt
    ! Result:
    integer :: ierr  ! returned value of function

    ! Locals:
    integer :: numSubGrids
    integer :: ezget_nsubGrids, gdxyfll, ezgprm
    character(len=1) :: grtyp
    integer :: ni, nj, ig1, ig2, ig3, ig4
    real(4) :: xpos_r4_vec(1), ypos_r4_vec(1), xpos2_r4_vec(1), ypos2_r4_vec(1)
    real(4) :: lat_deg_r4_vec(1),lon_deg_r4_vec(1)
    integer :: subGridIndex_vec(1)

    numSubGrids = ezget_nsubGrids(gdid)
    xpos2_r4 = mpc_missingValue_R4
    ypos2_r4 = mpc_missingValue_R4

    if (numSubGrids == 1) then

      ierr = ezgprm(gdid, grtyp, ni, nj, ig1, ig2, ig3, ig4)

      if (grtyp == 'Y') then

        !$omp critical
        if (present(maxGridSpacing_opt)) then
          ierr = gpos_xyfll_unstructGrid(gdid, xpos_r4, ypos_r4, lat_deg_r4, lon_deg_r4, maxGridSpacing_opt)
        else
          ierr = gpos_xyfll_unstructGrid(gdid, xpos_r4, ypos_r4, lat_deg_r4, lon_deg_r4)
        end if
        !$omp end critical

      else

        lat_deg_r4_vec(1) = lat_deg_r4
        lon_deg_r4_vec(1) = lon_deg_r4
        ! Not an unstructured grid, nor a Yin-Yan grid, call the standard ezscint routine
        ierr = gdxyfll(gdid, xpos_r4_vec, ypos_r4_vec, lat_deg_r4_vec, lon_deg_r4_vec, 1)
        xpos_r4 = xpos_r4_vec(1)
        ypos_r4 = ypos_r4_vec(1)
      end if

      subGridIndex = 1

    else

      ! This is a Yin-Yan grid, do something different
      !$omp critical
      lat_deg_r4_vec(1) = lat_deg_r4
      lon_deg_r4_vec(1) = lon_deg_r4
      ierr = gpos_xyfll_yinYanGrid(gdid, xpos_r4_vec, ypos_r4_vec, xpos2_r4_vec, ypos2_r4_vec, &
                                   lat_deg_r4_vec, lon_deg_r4_vec, subGridIndex_vec)
      xpos_r4 = xpos_r4_vec(1)
      ypos_r4 = ypos_r4_vec(1)
      xpos2_r4 = xpos2_r4_vec(1)
      ypos2_r4 = ypos2_r4_vec(1)
      subGridIndex = subGridIndex_vec(1)
      !$omp end critical

    end if

    if ( subGridIndex /= 3 ) then
      ! when only returning 1 position, copy values to pos2
      xpos2_r4 = xpos_r4
      ypos2_r4 = ypos_r4
    end if

  end function gpos_getPosXY_scalar

  !---------------------------------------------------------
  ! gpos_xyfll_yinYanGrid
  !---------------------------------------------------------
  function gpos_xyfll_yinYanGrid(gdid, xpos_r4, ypos_r4, xpos2_r4, ypos2_r4,  &
                                 lat_deg_r4, lon_deg_r4, subGridIndex) result(ierr)
    !
    ! :Purpose: Compute the grid XY position from a lat-lon for a Yin-Yan grid.
    !           It returns locations from both the Yin and Yan
    !           subgrids when in the overlap region.
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: gdid
    integer, intent(out) :: subGridIndex(:)
    real(4), intent(in)  :: lat_deg_r4(:)
    real(4), intent(in)  :: lon_deg_r4(:)
    real(4), intent(out) :: xpos_r4(:)
    real(4), intent(out) :: ypos_r4(:)
    real(4), intent(out) :: xpos2_r4(:)
    real(4), intent(out) :: ypos2_r4(:)
    ! Result:
    integer :: ierr  ! returned value of function

    ! Locals:
    integer :: ezget_subGridids, gdgaxes, gdxyfll, ezgprm
    integer :: EZscintIDvec(2)
    integer :: lonIndex, latIndex, numPoints, pointIndex
    integer :: ig1, ig2, ig3, ig4
    logical :: axesDifferent
    character(len=1) :: grtyp
    real(4), allocatable :: lonrot(:), latrot(:)
    real(4), allocatable :: xposYan_r4(:), yposYan_r4(:)
    ! Save some variables to reduce execution time
    integer, save :: EZscintIDvec1_old = MPC_missingValue_INT
    integer, save :: ni, nj
    real(4), allocatable, save :: ax_yin(:), ay_yin(:)

    numPoints = size(lat_deg_r4)
    ierr = ezget_subGridids(gdid, EZscintIDvec)
    allocate(xposYan_r4(numPoints))
    allocate(yposYan_r4(numPoints))
    allocate(latrot(numPoints))
    allocate(lonrot(numPoints))

    ! Compute positions on YIN
    ierr = gdxyfll(EZscintIDvec(1), xpos_r4, ypos_r4, lat_deg_r4, lon_deg_r4, numPoints)

    ! compute rotated lon and lat at obs location
    axesDifferent = (EZscintIDvec1_old /= EZscintIDvec(1))
    if (axesDifferent) then
      write(*,*) 'gpos_getPositionXY: axesDifferent, compute needed parameters'
      ! get ni nj of subGrid, assume same for both YIN and YAN
      ierr = ezgprm(EZscintIDvec(1), grtyp, ni, nj, ig1, ig2, ig3, ig4)
      if (allocated(ax_yin)) deallocate(ax_yin,ay_yin)
      allocate(ax_yin(ni),ay_yin(nj))
      ierr = gdgaxes(EZscintIDvec(1), ax_yin, ay_yin)
      EZscintIDvec1_old = EZscintIDvec(1)
    end if

    pointLoop: do pointIndex = 1, numPoints

      lonIndex = floor(xpos_r4(pointIndex))
      if ( lonIndex >= 1 .and. (lonIndex+1) <= ni ) then
        lonrot(pointIndex) = ax_yin(lonIndex) + (ax_yin(lonIndex+1) - ax_yin(lonIndex)) *  &
             (xpos_r4(pointIndex) - lonIndex)
      else
        lonrot(pointIndex) = MPC_missingValue_R4
      end if

      latIndex = floor(ypos_r4(pointIndex))
      if ( latIndex >= 1 .and. (latIndex+1) <= nj ) then
        latrot(pointIndex) = ay_yin(latIndex) + (ay_yin(latIndex+1) - ay_yin(latIndex)) *  &
             (ypos_r4(pointIndex) - latIndex)
      else
        latrot(pointIndex) = MPC_missingValue_R4
      end if
      subGridIndex(pointIndex) = 1

    end do pointLoop

    if ( any( lonrot(:) < 45.0 .or. lonrot(:) > 315.0 .or.  &
              latrot(:) < -45.0 .or. latrot(:) > 45.0 ) ) then

      ierr = gdxyfll(EZscintIDvec(2), xposYan_r4, yposYan_r4, lat_deg_r4, lon_deg_r4, numPoints)

      pointLoop2: do pointIndex = 1, numPoints

        if ( lonrot(pointIndex) < 45.0 .or. lonrot(pointIndex) > 315.0 .or.  &
             latrot(pointIndex) < -45.0 .or. latrot(pointIndex) > 45.0 ) then
          ! this approach is most similar to how ezsint works, preferentially take YIN
          ! Outside YIN, therefore use YAN (assume it is inside YAN)
          xpos_r4(pointIndex) = xposYan_r4(pointIndex)
          ypos_r4(pointIndex) = yposYan_r4(pointIndex) + real(nj) ! shift from YAN position to Supergrid position
          subGridIndex(pointIndex) = 2
        else
          subGridIndex(pointIndex) = 1
        end if

      end do pointLoop2

    end if

    ! These arguments should probably be removed
    xpos2_r4(:) = xpos_r4(:)
    ypos2_r4(:) = ypos_r4(:)

    deallocate(xposYan_r4)
    deallocate(yposYan_r4)
    deallocate(latrot)
    deallocate(lonrot)

  end function gpos_xyfll_yinYanGrid

  !---------------------------------------------------------
  ! gpos_xyfll_unstructGrid
  !---------------------------------------------------------
  function gpos_xyfll_unstructGrid(gdid, xpos_r4, ypos_r4, lat_deg_r4, lon_deg_r4, maxGridSpacing_opt) result(ierr)

    ! :Purpose: This function is used to interpolate from a RPN grid to a location
    !           specified by a latitude and longitude.
    !           It has special treatment for the ORCA12 and ORCA025 tri-polar grids
    !           used for sea ice and ocean models.
    !           The kdtree2 module is used to efficiently
    !           perform this task. The kdtree itself is constructed on the first call.
    implicit none

    ! Arguments:
    integer, intent(in)  :: gdid
    real(4), intent(in)  :: lat_deg_r4
    real(4), intent(in)  :: lon_deg_r4
    real(4), intent(out) :: xpos_r4    ! the x-position in the grid units, in[1,ni]
    real(4), intent(out) :: ypos_r4    ! the y-position in the grid units, in[1,nj]
    real(8), optional, intent(in) :: maxGridSpacing_opt
    ! Result:
    integer :: ierr  ! returned value of function

    ! Locals:
    real(8) :: lon_rad_r8, lat_rad_r8, pertGridLonRad, pertGridLatRad
    real(8), parameter :: deltaGrid = 0.00001d0
    real(8) :: pertPosition(3)
    real(8) :: gridSpacing
    real(8) :: gridSpacingSquared, lowerLeftCornerDistSquared
    real(8) :: lowerRightCornerDistSquared, upperLeftCornerDistSquared
    integer, save :: gdidOld = MPC_missingValue_INT
    integer :: nx, ny
    integer, save :: ni, nj
    real(8), allocatable, save :: grid_lon_rad(:,:), grid_lat_rad(:,:)
    real(4), save :: maxGridSpacing
    integer, save :: startXIndex, startYIndex, endXIndex, endYIndex
    character(len=1) :: grtyp
    integer :: ig1, ig2, ig3, ig4
    real(4), allocatable :: grid_lat_deg_r4(:,:), grid_lon_deg_r4(:,:)
    integer :: closePointsIndex
    integer :: gridIndex
    integer :: xIndex, yIndex, xIndexMin, xIndexMax, yIndexMin, yIndexMax
    integer :: ezgprm, gdll ! Functions
    integer                   :: numLocalGridPointsFound
    real(kdkind), allocatable :: positionArray(:,:)
    type(kdtree2_result)      :: searchResults(maxNumLocalGridPointsSearch)
    real(kdkind)              :: maxRadiusSquared
    real(kdkind)              :: refPosition(3)

    integer :: indexFound

    if (gdid /= gdidOld .and. gdidOld > 0) then
      write(*,*) 'gpos_xyfll_unstructGrid: gdid gdidOld = ', gdid, gdidOld
      call utl_abort('gpos_xyfll_unstructGrid: only one grid expected. Change code!')
    end if

    ! create the kdtree on the first call
    if (.not. associated(tree)) then
      write(*,*) 'gpos_xyfll_unstructGrid: start creating kdtree'
      call msg_memUsage('gpos_xyfll_unstructGrid')
      ierr = ezgprm(gdid, grtyp, ni, nj, ig1, ig2, ig3, ig4)

      startXIndex = 1
      startYIndex = 1
      endXIndex = ni
      endYIndex = nj
      if (present(maxGridSpacing_opt)) then
        maxGridSpacing = real(maxGridSpacing_opt,4)
      else
        maxGridSpacing = 10000.0
      end if

      if (allocated(grid_lat_rad)) deallocate(grid_lat_rad, grid_lon_rad)
      allocate(grid_lat_rad(ni, nj), grid_lon_rad(ni, nj))
      allocate(grid_lat_deg_r4(ni, nj), grid_lon_deg_r4(ni, nj))
      ierr = gdll(gdid, grid_lat_deg_r4, grid_lon_deg_r4)
      grid_lat_rad(:,:) = real(grid_lat_deg_r4(:,:), 8) * MPC_RADIANS_PER_DEGREE_R8
      grid_lon_rad(:,:) = real(grid_lon_deg_r4(:,:), 8) * MPC_RADIANS_PER_DEGREE_R8
      deallocate(grid_lat_deg_r4, grid_lon_deg_r4)
      gdidOld = gdid

      if( gpos_gridIsOrca(ni, nj, grid_lat_rad, grid_lon_rad) ) then
        ! The last 2 columns (i=ni-1 and ni) are repetitions of the first 2 columns (i=1 and 2)
        ! so both ends are removed from the search area.
        startXIndex = 2
        endXIndex = ni - 1
        ! The last line (j=nj) is a repetition in reverse order of the line (j=nj-2)
        ! and is thus eliminated from the search domain.
        endYIndex = nj - 1
      end if
      if(ni == 1442 .and. nj == 1021) then
        ! This is orca025
        maxGridSpacing = 40000.0
      end if

      nx = endXIndex - startXIndex + 1
      ny = endYIndex - startYIndex + 1

      allocate(positionArray(3, nx * ny))

      gridIndex = 0
      do yIndex = startYIndex, endYIndex
        do xIndex = startXIndex, endXIndex
          gridIndex = gridIndex + 1
          positionArray(:, gridIndex) = kdtree2_3dPosition(grid_lon_rad(xIndex, yIndex), &
                                                           grid_lat_rad(xIndex, yIndex))
        end do
      end do
      tree => kdtree2_create(positionArray, sort = .true., rearrange = .true.)
      write(*,*) 'gpos_xyfll_unstructGrid: done creating kdtree'
      call msg_memUsage('gpos_xyfll_unstructGrid')

    else

      ierr = 0
      nx = endXIndex - startXIndex + 1

    end if

    ! do the search
    maxRadiusSquared = maxGridSpacing**2
    lon_rad_r8 = real(lon_deg_r4, 8) * MPC_RADIANS_PER_DEGREE_R8
    lat_rad_r8 = real(lat_deg_r4, 8) * MPC_RADIANS_PER_DEGREE_R8
    refPosition(:) = kdtree2_3dPosition(lon_rad_r8, lat_rad_r8)

    call kdtree2_r_nearest(tp = tree, qv = refPosition, r2 = maxRadiusSquared, &
                           nfound = numLocalGridPointsFound, &
                           nalloc = maxNumLocalGridPointsSearch, results = searchResults)
    if (numLocalGridPointsFound > maxNumLocalGridPointsSearch) then
      write(*,*) 'gpos_xyfll_unstructGrid: ***where*** nfound > nalloc. obs lon, lat = ', lon_deg_r4, lat_deg_r4
      do indexFound = 1, maxNumLocalGridPointsSearch
        gridIndex = searchResults(indexFound)%idx
        yIndex = (gridIndex - 1) / nx + startYIndex
        xIndex = gridIndex - (yIndex - startYIndex) * nx + startXIndex - 1
        write(*,*) indexFound, sqrt(searchResults(indexFound)%dis) * 0.001, &
                   searchResults(indexFound)%idx, xIndex, yIndex, &
                   grid_lon_rad(xIndex, yIndex) * MPC_DEGREES_PER_RADIAN_R8, &
                   grid_lat_rad(xIndex, yIndex) * MPC_DEGREES_PER_RADIAN_R8
      end do
      call utl_abort('gpos_xyfll_unstructGrid: the parameter maxNumLocalGridPointsSearch must be increased')
    end if

    if (numLocalGridPointsFound < 4) then
      write(*,*) 'gpos_xyfll_unstructGrid: numLocalGridPointsFound = ', numLocalGridPointsFound
      write(*,*) 'gpos_xyfll_unstructGrid: search radius           = ', maxGridSpacing, ' meters'
      write(*,*) 'gpos_xyfll_unstructGrid: obs lon, lat = ', lon_deg_r4, lat_deg_r4
      write(*,*) 'gpos_xyfll_unstructGrid: ORCA025 WARNING: the search did not find 4 close points.'
      write(*,*) 'gpos_xyfll_unstructGrid: the x- and y-positions are set to -999.'
      xpos_r4 = MPC_missingValue_R4
      ypos_r4 = MPC_missingValue_R4
      return
    end if

    if (searchResults(1)%dis > maxRadiusSquared) then
      write(*,*) 'gpos_xyfll_unstructGrid: No grid point found within ', maxGridSpacing, ' meters'
      write(*,*) 'of the reference location lat-lon (degrees): ', lat_deg_r4, lon_deg_r4
      write(*,*) 'gpos_xyfll_unstructGrid: the x- and y-positions are set to -999.'
      xpos_r4 = MPC_missingValue_R4
      ypos_r4 = MPC_missingValue_R4
      return
    end if

    ! We found the closest grid point to the reference point.
    ! Now we need to determine in which of the 4 quadrants around the grid point the reference point lies.

    closePointsIndex = 1

    gridIndex = searchResults(closePointsIndex)%idx
    yIndex = (gridIndex - 1) / nx + startYIndex
    xIndex = gridIndex - (yIndex - startYIndex) * nx + startXIndex - 1

    ! Perturbing the closest grid point in the x direction

    if (xIndex < ni) then

      if (grid_lon_rad(xIndex + 1, yIndex) - grid_lon_rad(xIndex, yIndex) > 5.0d0) then
        pertGridLonRad = grid_lon_rad(xIndex, yIndex) + (grid_lon_rad(xIndex + 1, yIndex) - &
                         (grid_lon_rad(xIndex, yIndex) + 2.d0 * MPC_PI_R8)) * deltaGrid
      else if (grid_lon_rad(xIndex + 1, yIndex) - grid_lon_rad(xIndex, yIndex) < -5.0d0) then
        pertGridLonRad = grid_lon_rad(xIndex, yIndex) + ((grid_lon_rad(xIndex + 1, yIndex) + &
                         2.d0 * MPC_PI_R8) - grid_lon_rad(xIndex, yIndex)) * deltaGrid
      else
        pertGridLonRad = grid_lon_rad(xIndex, yIndex) + (grid_lon_rad(xIndex + 1, yIndex) - &
                         grid_lon_rad(xIndex, yIndex)) * deltaGrid
      end if
      pertGridLatRad = grid_lat_rad(xIndex, yIndex) + (grid_lat_rad(xIndex + 1, yIndex) - &
                       grid_lat_rad(xIndex, yIndex)) * deltaGrid

      ! Calculate the distance between the reference point and the perturbed grid point.
      pertPosition(:) = kdtree2_3dPosition(pertGridLonRad, pertGridLatRad)

      gridSpacingSquared = sum((pertPosition(:) - refPosition(:))**2)

      if (gridSpacingSquared < searchResults(closePointsIndex)%dis) then
        xIndexMin = xIndex
        xIndexMax = xIndex + 1
      else
        if (xIndex > 1) then
          xIndexMin = xIndex - 1
          xIndexMax = xIndex
        else
          write(*,*) 'xIndex yIndex ', xIndex, yIndex
          write(*,*) 'Closest distance = ', searchResults(closePointsIndex)%dis
          write(*,*) 'Perturbed distance = ', gridSpacingSquared
          write(*,*) 'gpos_xyfll_unstructGrid: of the reference location lat-lon (degrees): ', &
                     lat_deg_r4, lon_deg_r4
          write(*,*) 'gpos_xyfll_unstructGrid: ORCA025 WARNING (1): the reference point is outside the grid.'
          write(*,*) 'gpos_xyfll_unstructGrid: the x- and y-positions are set to -999.'
          xpos_r4 = MPC_missingValue_R4
          ypos_r4 = MPC_missingValue_R4
          return
        end if
      end if

    else

      if (grid_lon_rad(xIndex - 1, yIndex) - grid_lon_rad(xIndex, yIndex) > 5.0d0) then
        pertGridLonRad = grid_lon_rad(xIndex, yIndex) + (grid_lon_rad(xIndex - 1, yIndex) - &
                         (grid_lon_rad(xIndex, yIndex) + 2.d0 * MPC_PI_R8)) * deltaGrid
      else if (grid_lon_rad(xIndex - 1, yIndex) - grid_lon_rad(xIndex, yIndex) < -5.0d0) then
        pertGridLonRad = grid_lon_rad(xIndex, yIndex) + ((grid_lon_rad(xIndex - 1, yIndex) + &
                         2.d0 * MPC_PI_R8) - grid_lon_rad(xIndex, yIndex)) * deltaGrid
      else
        pertGridLonRad = grid_lon_rad(xIndex, yIndex) + (grid_lon_rad(xIndex - 1, yIndex) - &
                         grid_lon_rad(xIndex, yIndex)) * deltaGrid
      end if
      pertGridLatRad = grid_lat_rad(xIndex, yIndex) + (grid_lat_rad(xIndex - 1, yIndex) - &
                       grid_lat_rad(xIndex, yIndex)) * deltaGrid

      ! Calculate the distance between the reference point and the perturbed grid point.
      pertPosition(:) = kdtree2_3dPosition(pertGridLonRad, pertGridLatRad)

      gridSpacingSquared = sum((pertPosition(:) - refPosition(:))**2)

      if (gridSpacingSquared < searchResults(closePointsIndex)%dis) then
        xIndexMin = xIndex - 1
        xIndexMax = xIndex
      else
        write(*,*) 'xIndex yIndex ', xIndex, yIndex
        write(*,*) 'Closest distance = ', searchResults(closePointsIndex)%dis
        write(*,*) 'Perturbed distance = ', gridSpacingSquared
        write(*,*) 'gpos_xyfll_unstructGrid: ORCA025 WARNING (2): of the reference location lat-lon (degrees): ', &
                   lat_deg_r4, lon_deg_r4
        write(*,*) 'gpos_xyfll_unstructGrid: the x- and y-positions are set to -999.'
        xpos_r4 = MPC_missingValue_R4
        ypos_r4 = MPC_missingValue_R4
        return
      end if

    end if

    ! Perturbing the closest grid point in the y direction

    if (yIndex < nj) then

      if (grid_lon_rad(xIndex, yIndex + 1) - grid_lon_rad(xIndex, yIndex) > 5.0d0) then
        pertGridLonRad = grid_lon_rad(xIndex, yIndex) + (grid_lon_rad(xIndex, yIndex + 1) - &
                        (grid_lon_rad(xIndex, yIndex) + 2.d0 * MPC_PI_R8)) * deltaGrid
      else if (grid_lon_rad(xIndex, yIndex + 1) - grid_lon_rad(xIndex, yIndex) < -5.0d0) then
        pertGridLonRad = grid_lon_rad(xIndex, yIndex) + ((grid_lon_rad(xIndex, yIndex + 1) + &
                         2.d0 * MPC_PI_R8) - grid_lon_rad(xIndex,yIndex)) * deltaGrid
      else
        pertGridLonRad = grid_lon_rad(xIndex, yIndex) + (grid_lon_rad(xIndex, yIndex + 1) - &
                         grid_lon_rad(xIndex, yIndex)) * deltaGrid
      end if
      pertGridLatRad = grid_lat_rad(xIndex, yIndex) + (grid_lat_rad(xIndex, yIndex + 1) - &
                       grid_lat_rad(xIndex, yIndex)) * deltaGrid

      ! Calculate the distance between the reference point and the perturbed grid point.
      pertPosition(:) = kdtree2_3dPosition(pertGridLonRad, pertGridLatRad)

      gridSpacingSquared = sum((pertPosition(:) - refPosition(:))**2)

      if (gridSpacingSquared < searchResults(closePointsIndex)%dis) then
        yIndexMin = yIndex
        yIndexMax = yIndex + 1
      else
        if (yIndex > 1) then
          yIndexMin = yIndex - 1
          yIndexMax = yIndex
        else
          write(*,*) 'xIndex yIndex ', xIndex, yIndex
          write(*,*) 'Closest distance = ', searchResults(closePointsIndex)%dis
          write(*,*) 'Perturbed distance = ', gridSpacingSquared
          write(*,*) 'gpos_xyfll_unstructGrid: of the reference location lat-lon (degrees): ', lat_deg_r4, lon_deg_r4
          write(*,*) 'gpos_xyfll_unstructGrid: ORCA025 WARNING (3): the reference point is outside the grid.'
          write(*,*) 'gpos_xyfll_unstructGrid: the x- and y-positions are set to -999.'
          xpos_r4 = MPC_missingValue_R4
          ypos_r4 = MPC_missingValue_R4
          return
        end if
      end if

    else

      if (grid_lon_rad(xIndex, yIndex - 1) - grid_lon_rad(xIndex, yIndex) > 5.0d0) then
        pertGridLonRad = grid_lon_rad(xIndex, yIndex) + (grid_lon_rad(xIndex, yIndex - 1) - &
                         (grid_lon_rad(xIndex, yIndex) + 2.d0 * MPC_PI_R8)) * deltaGrid
      else if (grid_lon_rad(xIndex, yIndex - 1) - grid_lon_rad(xIndex, yIndex) < -5.0d0) then
        pertGridLonRad = grid_lon_rad(xIndex, yIndex) + ((grid_lon_rad(xIndex, yIndex - 1) + &
                         2.d0 * MPC_PI_R8) - grid_lon_rad(xIndex, yIndex)) * deltaGrid
      else
        pertGridLonRad = grid_lon_rad(xIndex, yIndex) + (grid_lon_rad(xIndex, yIndex - 1) - &
                         grid_lon_rad(xIndex, yIndex)) * deltaGrid
      end if
      pertGridLatRad = grid_lat_rad(xIndex, yIndex) + (grid_lat_rad(xIndex, yIndex - 1) - &
                       grid_lat_rad(xIndex, yIndex)) * deltaGrid

      ! Calculate the distance between the reference point and the perturbed grid point.
      pertPosition(:) = kdtree2_3dPosition(pertGridLonRad, pertGridLatRad)

      gridSpacingSquared = sum((pertPosition(:) - refPosition(:))**2)

      if (gridSpacingSquared < searchResults(closePointsIndex)%dis) then
        yIndexMin = yIndex - 1
        yIndexMax = yIndex
      else
        write(*,*) 'xIndex yIndex ', xIndex, yIndex
        write(*,*) 'Closest distance = ', searchResults(closePointsIndex)%dis
        write(*,*) 'Perturbed distance = ', gridSpacingSquared
        write(*,*) 'gpos_xyfll_unstructGrid: of the reference location lat-lon (degrees): ', lat_deg_r4, lon_deg_r4
        write(*,*) 'ORCA025 WARNING 4: the reference point is outside the grid. We will set the positions to -999.!!!'
        xpos_r4 = MPC_missingValue_R4
        ypos_r4 = MPC_missingValue_R4
        return
      end if

    end if

    ! Calculate real x and y position in the grid.
    ! kdtree returns distance in square meters !

    gridSpacingSquared = phf_calcDistance(grid_lat_rad(xIndexMin, yIndexMin), grid_lon_rad(xIndexMin, yIndexMin), &
           grid_lat_rad(xIndexMax, yIndexMin), grid_lon_rad(xIndexMax, yIndexMin))**2

    lowerLeftCornerDistSquared = phf_calcDistance(grid_lat_rad(xIndexMin, yIndexMin), grid_lon_rad(xIndexMin, yIndexMin), &
           lat_rad_r8, lon_rad_r8)**2

    lowerRightCornerDistSquared = phf_calcDistance(grid_lat_rad(xIndexMax, yIndexMin), grid_lon_rad(xIndexMax, yIndexMin), &
           lat_rad_r8, lon_rad_r8)**2

    xpos_r4 = real(real(xIndexMin,8) + (lowerLeftCornerDistSquared + gridSpacingSquared - &
                                        lowerRightCornerDistSquared) / (2.0 * (gridSpacingSquared)), 4)

    gridSpacingSquared = phf_calcDistance(grid_lat_rad(xIndexMin, yIndexMin), grid_lon_rad(xIndexMin, yIndexMin), &
           grid_lat_rad(xIndexMin, yIndexMax), grid_lon_rad(xIndexMin, yIndexMax))**2

    upperLeftCornerDistSquared = phf_calcDistance(grid_lat_rad(xIndexMin, yIndexMax), grid_lon_rad(xIndexMin, yIndexMax), &
           lat_rad_r8, lon_rad_r8)**2

    ypos_r4 = real(real(yIndexMin,8) + (lowerLeftCornerDistSquared + gridSpacingSquared  - &
                                        upperLeftCornerDistSquared) / (2.0 * (gridSpacingSquared)), 4)

    if ( abs(ypos_r4 - yIndexMin) > 2.0 .or. abs(xpos_r4 - xIndexMin) > 4.3 ) then
      write(*,*) 'xpos_r4 = ',xpos_r4
      write(*,*) 'ypos_r4 = ',ypos_r4
      write(*,*) 'xIndexMin = ',xIndexMin
      write(*,*) 'yIndexMin = ',yIndexMin
      do closePointsIndex = 1,min(numLocalGridPointsFound,5)
        gridIndex = searchResults(closePointsIndex)%idx
        yIndex = (gridIndex - 1) / nx + startYIndex
        xIndex = gridIndex - (yIndex - startYIndex) * nx
        write(*,*) 'closePointsIndex: ', closePointsIndex
        write(*,*) 'xIndex yIndex ', xIndex, yIndex
        write(*,*) 'idx: ', searchResults(closePointsIndex)%idx
        write(*,*) 'distance (meters): ', sqrt(searchResults(closePointsIndex)%dis)
      end do
      write(*,*) 'lowerLeftCornerDist = ', sqrt(lowerLeftCornerDistSquared)
      write(*,*) 'lowerRightCornerDist = ', sqrt(lowerRightCornerDistSquared)
      write(*,*) 'upperLeftCornerDist = ', sqrt(upperLeftCornerDistSquared)
      write(*,*) 'y gridSpacing = ', sqrt(gridSpacingSquared)
      gridSpacing = phf_calcDistance(grid_lat_rad(xIndexMin, yIndexMin), grid_lon_rad(xIndexMin, yIndexMin), &
           grid_lat_rad(xIndexMax, yIndexMin), grid_lon_rad(xIndexMax, yIndexMin))
      write(*,*) 'x gridSpacing = ', gridSpacing
      write(*,*) 'lon-lat (deg):'
      write(*,*) 'reference point: ', lon_deg_r4, lat_deg_r4
      write(*,*) xIndexMin, yIndexMin, grid_lon_rad(xIndexMin, yIndexMin) * MPC_DEGREES_PER_RADIAN_R8, &
                 grid_lat_rad(xIndexMin,yIndexMin) * MPC_DEGREES_PER_RADIAN_R8
      write(*,*) xIndexMax, yIndexMin, grid_lon_rad(xIndexMax, yIndexMin) * MPC_DEGREES_PER_RADIAN_R8, &
                 grid_lat_rad(xIndexMax, yIndexMin) * MPC_DEGREES_PER_RADIAN_R8
      write(*,*) xIndexMax+1, yIndexMin, grid_lon_rad(xIndexMax + 1, yIndexMin) * MPC_DEGREES_PER_RADIAN_R8,&
                 grid_lat_rad(xIndexMax + 1, yIndexMin) * MPC_DEGREES_PER_RADIAN_R8
      write(*,*) xIndexMin, yIndexMax, grid_lon_rad(xIndexMin, yIndexMax) * MPC_DEGREES_PER_RADIAN_R8, &
                 grid_lat_rad(xIndexMin, yIndexMax) * MPC_DEGREES_PER_RADIAN_R8
      write(*,*) xIndexMax, yIndexMax, grid_lon_rad(xIndexMax, yIndexMax) * MPC_DEGREES_PER_RADIAN_R8, &
                 grid_lat_rad(xIndexMax, yIndexMax) * MPC_DEGREES_PER_RADIAN_R8
      write(*,*) xIndexMax+1, yIndexMax, grid_lon_rad(xIndexMax + 1, yIndexMax) * MPC_DEGREES_PER_RADIAN_R8, &
                 grid_lat_rad(xIndexMax + 1, yIndexMax) * MPC_DEGREES_PER_RADIAN_R8
      write(*,*) xIndexMin, yIndexMax + 1, grid_lon_rad(xIndexMin, yIndexMax + 1) * MPC_DEGREES_PER_RADIAN_R8, &
                 grid_lat_rad(xIndexMin, yIndexMax + 1) * MPC_DEGREES_PER_RADIAN_R8
      write(*,*) xIndexMax, yIndexMax + 1, grid_lon_rad(xIndexMax, yIndexMax + 1) * MPC_DEGREES_PER_RADIAN_R8, &
                 grid_lat_rad(xIndexMax, yIndexMax + 1) * MPC_DEGREES_PER_RADIAN_R8
      write(*,*) xIndexMax + 1, yIndexMax + 1, grid_lon_rad(xIndexMax + 1, yIndexMax + 1) * MPC_DEGREES_PER_RADIAN_R8, &
                 grid_lat_rad(xIndexMax + 1, yIndexMax + 1) * MPC_DEGREES_PER_RADIAN_R8
    end if

  end function gpos_xyfll_unstructGrid

  !---------------------------------------------------------
  ! gpos_gridIsOrca
  !---------------------------------------------------------
  function gpos_gridIsOrca( ni, nj, latitude, longitude ) result(gridIsOrca)

    ! :Purpose: Check if the grid is of the orca tri-polar family
    implicit none

   ! Arguments:
    integer, intent(in) :: ni                ! first  dimension of the grid
    integer, intent(in) :: nj                ! second dimension of the grid
    real(8), intent(in) :: longitude(ni,nj)  ! longitude (degrees or radians)
    real(8), intent(in) :: latitude(ni,nj)   ! latitude  (degrees or radians)
    ! Result:
    logical :: gridIsOrca  ! returned value of function

    ! Locals:
    integer :: xIndex, yIndex

    gridIsOrca = .true.

    ! The last 2 columns (ni-1 and ni) are repetitions of the first 2 columns (i=1 and 2)
    ! Check that is true
    do yIndex = 1, nj
      if(.not. utl_isEqual(latitude(1,yIndex), latitude(ni-1,yIndex))   .or. &
         .not. utl_isEqual(latitude(2,yIndex), latitude(ni,  yIndex))   .or. &
         .not. utl_isEqual(longitude(1,yIndex), longitude(ni-1,yIndex)) .or. &
         .not. utl_isEqual(longitude(2,yIndex), longitude(ni,  yIndex))) then
        gridIsOrca = .false.
        write(*,*) '1. Not an orca grid',latitude(1,yIndex),latitude(ni-1,yIndex),latitude(2,yIndex),latitude(ni,  yIndex),longitude(1,yIndex),longitude(ni-1,yIndex),longitude(2,yIndex),longitude(ni,  yIndex)
        return
      end if
    end do

    ! The last line (j=nj) is a repetition in reverse order of the line (j=nj-2)
    ! Check that is true (does not work at (ni/2 + 1) for orca025 lats as with the letkf ocean unit test)
    do xIndex = 2, ni
      if( .not. utl_isEqual(longitude(xIndex,nj), longitude(ni-xIndex+2,nj-2)) ) then
        gridIsOrca = .false.
        write(*,*) '2. Not an orca grid',xIndex,longitude(xIndex,nj),longitude(ni-xIndex+2,nj-2)
        return
      end if
    end do
    do xIndex = 2, ni/2
      if( .not. utl_isEqual(latitude(xIndex,nj), latitude(ni-xIndex+2,nj-2)) ) then
        gridIsOrca = .false.
        write(*,*) '3. Not an orca grid',xIndex,latitude(xIndex,nj),latitude(ni-xIndex+2,nj-2)
        return
      end if
    end do
    do xIndex = ni/2 + 2, ni
      if( .not. utl_isEqual(latitude(xIndex,nj), latitude(ni-xIndex+2,nj-2)) ) then
        gridIsOrca = .false.
        write(*,*) '3. Not an orca grid',xIndex,latitude(xIndex,nj),latitude(ni-xIndex+2,nj-2)
        return
      end if
    end do

  end function gpos_gridIsOrca

end module getGridPosition_mod

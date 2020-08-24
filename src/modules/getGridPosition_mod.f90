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

module getGridPosition_mod
  ! MODULE getGridPosition_mod (prefix='gpos' category='7. Low-level data objects and utilities')
  !
  ! :Purpose: A place to collect numerous interpolation related routines
  !
  use kdtree2_mod
  use mathPhysConstants_mod
  use utilities_mod
  use earthconstants_mod

  implicit none
  save
  private

  ! public procedures
  public :: gpos_getPositionXY

  integer, parameter :: maxNumLocalGridPointsSearch = 100
  integer, external  :: get_max_rss

  type(kdtree2), pointer :: tree => null()

contains


  function gpos_getPositionXY( gdid, xpos_r4, ypos_r4, xpos2_r4, ypos2_r4,  &
                          lat_deg_r4, lon_deg_r4, subGridIndex ) result(ierr)
    !
    ! :Purpose: Compute the grid XY position from a lat-lon. This
    !           simply calls the ezsint routine gdxyfll for simple grids. For
    !           Yin-Yan grids it can return locations from both the Yin and Yan
    !           subgrids when in the overlap region, depending on the logical 
    !           variable `useSingleValueOverlap`. There is also support for
    !           RPN Y grids, in which case it calls the subroutine gpos_xyfll.
    !
    implicit none

    ! arguments
    integer :: ierr  ! returned value of function
    integer, intent(in) :: gdid
    integer, intent(out) :: subGridIndex
    real(4), intent(out) :: xpos_r4
    real(4), intent(out) :: ypos_r4
    real(4), intent(out) :: xpos2_r4
    real(4), intent(out) :: ypos2_r4
    real(4), intent(in) :: lat_deg_r4
    real(4), intent(in) :: lon_deg_r4

    ! locals
    integer :: numSubGrids
    integer :: ezget_nsubGrids, ezget_subGridids, gdxyfll, ezgprm, gdgaxes
    integer :: EZscintIDvec(2)
    integer, save :: EZscintIDvec1_old = -999
    character(len=1) :: grtyp
    integer :: ni, nj, ig1, ig2, ig3, ig4, lonIndex, latIndex
    real :: lonrot, latrot
    real, allocatable, save :: ax_yin(:), ay_yin(:), ax_yan(:), ay_yan(:)
    logical :: axesDifferent

    ! this controls which approach to use for interpolation within the YIN-YAN overlap
    logical :: useSingleValueOverlap = .true.  

    numSubGrids = ezget_nsubGrids(gdid)
    xpos2_r4 = -999.0
    ypos2_r4 = -999.0

    if ( numSubGrids == 1 ) then

      ierr = ezgprm(gdid, grtyp, ni, nj, ig1, ig2, ig3, ig4)

      if (grtyp == 'Y') then

         ierr = gpos_xyfll(gdid, xpos_r4, ypos_r4, lat_deg_r4, lon_deg_r4)

      else

        ! Not a Yin-Yang grid, call the standard ezscint routine
         ierr = gdxyfll(gdid, xpos_r4, ypos_r4, lat_deg_r4, lon_deg_r4, 1)

      end if

      subGridIndex = 1

    else

      ! This is a Yin-Yang grid, do something different

      ierr = ezget_subGridids(gdid, EZscintIDvec)
      ! get ni nj of subGrid, assume same for both YIN and YANG
      ierr = ezgprm(EZscintIDvec(1), grtyp, ni, nj, ig1, ig2, ig3, ig4)

      ! first check YIN
      ierr = gdxyfll(EZscintIDvec(1), xpos_r4, ypos_r4, lat_deg_r4, lon_deg_r4, 1)

      ! compute rotated lon and lat at obs location
      axesDifferent = (EZscintIDvec1_old /= EZscintIDvec(1))
      if (axesDifferent) then
        write(*,*) 'gpos_getPositionXY: axesDifferent, compute needed parameters'
        if (allocated(ax_yin)) deallocate(ax_yin,ay_yin)
        allocate(ax_yin(ni),ay_yin(nj))
        ierr = gdgaxes(EZscintIDvec(1), ax_yin, ay_yin)
        EZscintIDvec1_old = EZscintIDvec(1)
      end if
      lonIndex = floor(xpos_r4)
      if ( lonIndex >= 1 .and. (lonIndex+1) <= ni ) then
        lonrot = ax_yin(lonIndex) + (ax_yin(lonIndex+1) - ax_yin(lonIndex)) *  &
                 (xpos_r4 - lonIndex)
      else
        lonrot = -999.0
      end if
      latIndex = floor(ypos_r4)
      if ( latIndex >= 1 .and. (latIndex+1) <= nj ) then
        latrot = ay_yin(latIndex) + (ay_yin(latIndex+1) - ay_yin(latIndex)) *  &
                 (ypos_r4 - latIndex)
      else
        latrot = -999.0
      end if
      subGridIndex = 1

      if ( useSingleValueOverlap ) then

        ! this approach is most similar to how ezsint works, preferentially take YIN

        if ( lonrot < 45.0 .or. lonrot > 315.0 .or. latrot < -45.0 .or. latrot > 45.0 ) then
          ! Outside YIN, therefore use YANG (assume it is inside YANG)
          ierr = gdxyfll(EZscintIDvec(2), xpos_r4, ypos_r4, lat_deg_r4, lon_deg_r4, 1)
          ypos_r4 = ypos_r4 + real(nj) ! shift from YANG position to Supergrid position
          subGridIndex = 2
        else
          subGridIndex = 1
        end if

      else ! not useSingleValueOverlap

        ! this approach returns both the YIN and YAN locations when point is inside both

        if ( lonrot < 45.0 .or. lonrot > 315.0 .or. latrot < -45.0 .or. latrot > 45.0 ) then
          ! Outside YIN, therefore use YANG (assume it is inside YANG)
          ierr = gdxyfll(EZscintIDvec(2), xpos_r4, ypos_r4, lat_deg_r4, lon_deg_r4, 1)
          ypos_r4 = ypos_r4 + real(nj) ! shift from YANG position to Supergrid position
          subGridIndex = 2
        else
          ! inside YIN, check if also inside YANG
          allocate(ax_yan(ni),ay_yan(nj))
          ierr = gdgaxes(EZscintIDvec(2), ax_yan, ay_yan)
          ierr = gdxyfll(EZscintIDvec(2), xpos2_r4, ypos2_r4, lat_deg_r4, lon_deg_r4, 1)
          if ( lonIndex >= 1 .and. (lonIndex+1) <= ni ) then
            lonrot = ax_yan(lonIndex) + (ax_yan(lonIndex+1) - ax_yan(lonIndex)) *  &
                     (xpos2_r4 - lonIndex)
          else
            lonrot = -999.0
          end if
          latIndex = floor(ypos2_r4)
          if ( latIndex >= 1 .and. (latIndex+1) <= nj ) then
            latrot = ay_yan(latIndex) + (ay_yan(latIndex+1) - ay_yan(latIndex)) *  &
                     (ypos2_r4 - latIndex)
          else
            latrot = -999.0
          end if
          deallocate(ax_yan,ay_yan)
          if ( lonrot < 45.0 .or. lonrot > 315.0 .or. latrot < -45.0 .or. latrot > 45.0 ) then
            ! outside YANG, only inside YIN
            xpos2_r4 = -999.0
            ypos2_r4 = -999.0
            subGridIndex = 1
          else
            ! inside both YIN and YANG
            ypos2_r4 = ypos2_r4 + real(nj) ! shift from YANG position to Supergrid position
            subGridIndex = 3
          end if
        end if

      end if

    end if    

    if ( subGridIndex /= 3 ) then
      ! when only returning 1 position, copy values to pos2
      xpos2_r4 = xpos_r4
      ypos2_r4 = ypos_r4
    end if

  end function gpos_getPositionXY

  function gpos_xyfll( gdid, xpos_r4, ypos_r4, lat_deg_r4, lon_deg_r4 ) result(ierr)

    ! :Purpose: This function is used to interpolate from a RPN grid to a location
    !           specified by a latitude and longitude.
    !           It has special treatment for the ORCA12 and ORCA025 tri-polar grids
    !           used for sea ice and ocean models.
    !           The kdtree2 module is used to efficiently
    !           perform this task. The kdtree itself is constructed on the first call.
    implicit none

    ! arguments
    integer :: ierr  ! returned value of function
    integer, intent(in)  :: gdid
    real(4), intent(in)  :: lat_deg_r4, lon_deg_r4
    real(4), intent(out) :: xpos_r4, ypos_r4

    ! Local Variables

    real(8) :: lon_rad_r8, lat_rad_r8

    integer, save :: gdid_old = -999
    integer, save :: nx
    integer :: ni, nj, ny
    real(8), allocatable, save  :: grid_lon_rad(:,:), grid_lat_rad(:,:)
    logical :: global
    real(4), save :: max_grid_spacing

    character(len=1) :: grtyp
    integer :: ig1, ig2, ig3, ig4
    real(4), allocatable :: grid_lat_deg_r4(:,:), grid_lon_deg_r4(:,:)
    integer :: xIndex, yIndex, gridIndex

    ! Functions
    integer :: ezgprm, gdll

    integer                   :: numLocalGridPointsFound
    real(kdkind), allocatable :: positionArray(:,:)
    type(kdtree2_result)      :: searchResults(maxNumLocalGridPointsSearch)
    real(kdkind)              :: maxRadius
    real(kdkind)              :: refPosition(3)

    if ( gdid /= gdid_old .and. gdid_old > 0) then
      write(*,*) 'gpos_xyfll: gdid gdid_old = ',gdid,gdid_old
      call utl_abort('gpos_xyfll: only one grid expected. Change code !')
    end if

    ! create the kdtree on the first call
    if (.not. associated(tree)) then
      write(*,*) 'gpos_xyfll: start creating kdtree'
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
      ierr = ezgprm(gdid, grtyp, ni, nj, ig1, ig2, ig3, ig4)

      global = .false.
      max_grid_spacing = 10000.0

      if(ni == 4322 .and. nj == 3059) then
        ! This is orca12
        global = .true.
      end if
      if(ni == 1442 .and. nj == 1021) then
        ! This is orca025
        global = .true.
        max_grid_spacing = 27000.0
      end if

      if(global) then
        ! The last 2 columns (i=4321 and 4322) are repetitions of the first 2 columns (i=1 and 2)
        ! and are thus eliminated from the search domain.
        ! The last line (j=3059) is a repetition in reverse order of the line j=3057
        ! and is thus eliminated from the search domain.
        nx = ni - 2
        ny = nj - 1
      else
        nx = ni
        ny = nj
      end if

      allocate(positionArray(3,nx*ny))
      if (allocated(grid_lat_rad)) deallocate(grid_lat_rad, grid_lon_rad)
      allocate(grid_lat_rad(ni,nj), grid_lon_rad(ni,nj))
      allocate(grid_lat_deg_r4(ni,nj), grid_lon_deg_r4(ni,nj))
      ierr = gdll(gdid, grid_lat_deg_r4, grid_lon_deg_r4)
      grid_lat_rad(:,:) = real(grid_lat_deg_r4(:,:),8)*MPC_RADIANS_PER_DEGREE_R8
      grid_lon_rad(:,:) = real(grid_lon_deg_r4(:,:),8)*MPC_RADIANS_PER_DEGREE_R8
      deallocate(grid_lat_deg_r4, grid_lon_deg_r4)
      gdid_old = gdid

      gridIndex = 0
      do yIndex = 1, ny
        do xIndex = 1, nx
          gridIndex = gridIndex + 1
          positionArray(1,gridIndex) = RA * sin(grid_lon_rad(xIndex,yIndex)) * cos(grid_lat_rad(xIndex,yIndex))
          positionArray(2,gridIndex) = RA * cos(grid_lon_rad(xIndex,yIndex)) * cos(grid_lat_rad(xIndex,yIndex))
          positionArray(3,gridIndex) = RA * sin(grid_lat_rad(xIndex,yIndex))
        end do
      end do
      tree => kdtree2_create(positionArray, sort=.true., rearrange=.true.) 
      write(*,*) 'gpos_xyfll: done creating kdtree'
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    else

      ierr = 0

    end if

    ! do the search
    maxRadius = max_grid_spacing**2
    lon_rad_r8 = real(lon_deg_r4,8)*MPC_RADIANS_PER_DEGREE_R8
    lat_rad_r8 = real(lat_deg_r4,8)*MPC_RADIANS_PER_DEGREE_R8
    refPosition(1) = RA * sin(lon_rad_r8) * cos(lat_rad_r8)
    refPosition(2) = RA * cos(lon_rad_r8) * cos(lat_rad_r8)
    refPosition(3) = RA * sin(lat_rad_r8)

    call kdtree2_r_nearest(tp=tree, qv=refPosition, r2=maxRadius, nfound=numLocalGridPointsFound,&
                           nalloc=maxNumLocalGridPointsSearch, results=searchResults)
    if (numLocalGridPointsFound > maxNumLocalGridPointsSearch) then
      call utl_abort('gpos_xyfll: the parameter maxNumLocalGridPointsSearch must be increased')
    end if

    if (numLocalGridPointsFound < 1) then
      write(*,*) 'gpos_xyfll: numLocalGridPointsFound = ',numLocalGridPointsFound
      call utl_abort('gpos_xyfll: the search did not found close points.')
    end if

    ! Get closest grid point from the observation
    gridIndex = searchResults(1)%idx
    yIndex = (gridIndex-1)/nx + 1
    xIndex = gridIndex - (yIndex-1)*nx
    xpos_r4 = real(xIndex)
    ypos_r4 = real(yIndex)

  end function gpos_xyfll

end module getGridPosition_mod

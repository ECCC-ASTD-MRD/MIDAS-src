
module oceanMask_mod
  ! MODULE oceanMask_mod (prefix='ocm' category='7. Low-level data objects')
  !
  !:Purpose:  The horizontal mask indicating land (.false.) and water (.true.) grid points.
  !           This mask is either:
  !                 * In the case of variables on ocean depth levels, it varies with vertical level.
  !                 * In other cases it is a single 2D field used for all variables.
  !
  use rmn_fst98
  use midasMpi_mod
  use kdTree2_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use varNameList_mod
  use utilities_mod
  use mathPhysConstants_mod
  use earthConstants_mod
  use message_mod

  implicit none
  save
  private

  ! Public derived type
  public :: struct_ocm

  ! Public subroutines and functions
  public :: ocm_readMaskFromFile, ocm_deallocate
  public :: ocm_copyMask, ocm_communicateMask
  public :: ocm_farFromLand
  public :: ocm_copyToInt, ocm_copyFromInt
  public :: ocm_computeMinGridSpacing

  type struct_ocm
    ! This is the derived type of the ocean mask object
    integer                   :: nLev
    logical, pointer          :: mask(:,:,:) => null()
    logical                   :: maskPresent = .false.
    type(struct_hco), pointer :: hco         => null()
  end type struct_ocm

  contains

  !--------------------------------------------------------------------------
  ! ocm_readMaskFromFile
  !--------------------------------------------------------------------------
  subroutine ocm_readMaskFromFile(oceanMask, hco, vco, inputFileName)
    !
    ! :Purpose: Check if any mask fields exist for surface or ocean depth levels.
    !
    implicit none

    ! Arguments:
    type(struct_ocm),          intent(inout) :: oceanMask
    type(struct_hco), pointer, intent(in)    :: hco
    type(struct_vco),          intent(in)    :: vco
    character(len=*),          intent(in)    :: inputFileName

    ! Constants:
    character(len=*), parameter :: netcdfFileExtention = '_fst'
    ! Locals:
    integer :: nulfile, ierr, ni_file, nj_file, nk_file
    integer :: ikey, levIndex
    integer, allocatable :: mask(:,:)
    integer :: maxkeys
    integer :: ip1
    logical :: fileExist
    character(len=len_trim(inputFileName)+len_trim(netcdfFileExtention)) :: fileName
    ! external definitions
    integer, external :: fnom, fclos

    if (trim(utl_fileType(inputFileName)) == 'NetCDF') then
      fileName = trim(inputFileName) // netcdfFileExtention
      inquire(file = trim(fileName), exist = fileExist)
      if (.not. fileExist) then
        call utl_abort('ocm_readMaskFromFile: mandatory FST file does not exist: '//&
                       trim(fileName))
      end if
    else
      fileName = trim(inputFileName)
    end if
    call msg('ocm_readMaskFromFile', 'File name = '//trim(fileName))

    ! Check if any mask is present in file, return if not
    if ( .not. vnl_varNamePresentInFile(' ', fileName = trim(fileName), typvar_opt='@@')) then
      return
    end if

    !- Open input file
    nulfile = 0
    ierr = fnom(nulfile, trim(fileName), 'RND+OLD+R/O', 0)

    if (ierr >= 0) then
      maxkeys = fstouv(nulfile, 'RND+OLD')
    else
      call utl_abort('ocm_readMaskFromFile problem opening input file')
    end if

    if (nulfile == 0) then
      call utl_abort('ocm_readMaskFromFile: unit number for input file not valid')
    end if

    ! Read mask for all fields
    if ( vco%nLev_depth > 0 ) then
      lev_loop: do levIndex = 1, vco%nLev_depth

        ip1 = vco%ip1_depth(levIndex)

        ! Make sure that the mask for this variable has the same grid size as hco
        ikey = fstinf(nulfile, ni_file, nj_file, nk_file, &
                      -1, ' ', ip1, -1, -1, '@@', ' ')

        if (ikey < 0) then
          call msg('ocm_readMaskFromFile', 'Searched for mask with ip1: ' // str(ip1))
          call utl_abort('ocm_readMaskFromFile: cannot find mask for this ip1 in file ' // trim(fileName))
        end if

        do while (ni_file /= hco%ni .or. nj_file /= hco%nj)
          ikey = fstsui(nulfile, ni_file, nj_file, nk_file)
        end do

        if (ni_file == hco%ni .and. nj_file == hco%nj) then
          call msg('ocm_readMaskFromFile', 'Read mask for ip1: ' // str(ip1))
          if (.not. associated(oceanMask%mask)) then
            call ocm_allocate(oceanMask, hco, vco%nLev_depth)
          end if
          if (.not. allocated(mask)) allocate(mask(hco%ni, hco%nj))
          ierr = fstluk(mask(:,:), ikey, ni_file, nj_file, nk_file)
          call ocm_copyFromInt(oceanMask, mask, levIndex)
          if (ierr < 0) then
            call utl_abort('ocm_readMaskFromFile: error when reading mask record')
          end if
        else
          ! Special cases for variables that are on a different horizontal grid in LAM (e.g. TG)
          write(*,*)
          call msg('ocm_readMaskFromFile', 'Mask is on a different horizontal grid')
          call msg('ocm_readMaskFromFile', 'ni = ' // str(ni_file) //','// str(hco%ni)// &
                                           ', nj = ' // str(nj_file)//','// str(hco%nj))
          call utl_abort('ocm_readMaskFromFile: This is not allowed at the moment')
        end if

      end do lev_loop

    else
      ! ***No depth levels, so just look for any mask field***

      ! Make sure that the mask for this variable has the same grid size as hco
      ikey = fstinf(nulfile, ni_file, nj_file, nk_file, &
                    -1, ' ', -1, -1, -1, '@@', ' ')

      if (ikey < 0) then
        call utl_abort('ocm_readMaskFromFile: cannot find any mask in file ' // trim(fileName))
      end if

      do while (ni_file /= hco%ni .or. nj_file /= hco%nj)
        ikey = fstsui(nulfile, ni_file, nj_file, nk_file)
      end do

      if (ni_file == hco%ni .and. nj_file == hco%nj) then
        call msg('ocm_readMaskFromFile', 'Reading mask')
        if (.not. associated(oceanMask%mask)) then
          call ocm_allocate(oceanMask, hco, 1)
          if (.not. allocated(mask)) allocate(mask(hco%ni, hco%nj))
          ierr = fstluk(mask(:,:), ikey, ni_file, nj_file, nk_file)
          call ocm_copyFromInt(oceanMask, mask, 1)
          if (ierr < 0) then
            call utl_abort('ocm_readMaskFromFile: error when reading mask record')
          end if
        end if
      else
        ! Special cases for variables that are on a different horizontal grid in LAM (e.g. TG)
        write(*,*)
        call msg('ocm_readMaskFromFile', 'Mask is on a different horizontal grid')
        call msg('ocm_readMaskFromFile', 'ni = ' // str(ni_file)//','//str(hco%ni)//&
                                         ', nj = '// str(nj_file)//','//str(hco%nj))
        call utl_abort('ocm_readMaskFromFile: This is not allowed at the moment')
      end if

    end if

    ierr = fstfrm(nulfile)
    ierr = fclos(nulfile)

    if (allocated(mask)) deallocate(mask)

  end subroutine ocm_readMaskFromFile

  !--------------------------------------------------------------------------
  ! ocm_farFromLand
  !--------------------------------------------------------------------------
  function ocm_farFromLand(oceanMask, levIndex, lon, lat, distanceToLand) result(farFromLand)
    !
    ! :Purpose: Determine if the supplied lat/lat location is far from the
    !           nearest land based on the supplied 'distanceToLand'.
    !
    implicit none

    ! Arguments:
    type(struct_ocm), intent(in) :: oceanMask
    integer,          intent(in) :: levIndex
    real(8),          intent(in) :: lon
    real(8),          intent(in) :: lat
    real(8),          intent(in) :: distanceToLand
    ! Result:
    logical                      :: farFromLand

    ! Locals:
    integer, parameter           :: maxNumLocalGridPointsSearch = 200000
    type(kdtree2), save, pointer :: tree => null()
    integer                      :: ni, nj, xIndex, yIndex, gridIndex
    integer                      :: numTotalLandPoints, numLocalGridPointsFound
    real(kdkind), allocatable    :: positionArray(:,:)
    type(kdtree2_result)         :: searchResults(maxNumLocalGridPointsSearch)
    real(kdkind)                 :: searchRadiusSquared
    real(kdkind)                 :: refPosition(3)

    ! do some basic checks
    if (.not.oceanMask%maskPresent .or. .not.associated(oceanMask%mask)) then
      call utl_abort('ocm_farFromLand: mask is not allocated')
    end if
    if (levIndex < 0 .or. levIndex > oceanMask%nLev) then
      call utl_abort('ocm_farFromLand: specified levIndex not valid')
    end if

    ni = oceanMask%hco%ni
    nj = oceanMask%hco%nj
    searchRadiusSquared = (1.1D0 * distanceToLand)**2

    ! create the kdtree on the first call
    if (.not. associated(tree)) then
      write(*,*) 'ocm_farFromLand: start creating kdtree'
      call msg_memUsage('ocm_farFromLand')

      numTotalLandPoints = count(.not. oceanMask%mask(:,:,levIndex))
      allocate(positionArray(3,numTotalLandPoints))

      gridIndex = 0
      do xIndex = 1, ni
        do yIndex = 1, nj
          if (.not. oceanMask%mask(xIndex,yIndex,levIndex)) then
            gridIndex = gridIndex + 1
            positionArray(:,gridIndex) = &
                 kdtree2_3dPosition(real(oceanMask%hco%lon2d_4(xIndex,yIndex),8), &
                                    real(oceanMask%hco%lat2d_4(xIndex,yIndex),8))
          end if
        end do
      end do
      tree => kdtree2_create(positionArray, sort=.true., rearrange=.true.)
      write(*,*) 'ocm_farFromLand: done creating kdtree'
      call msg_memUsage('ocm_farFromLand')

    end if

    ! do the search
    refPosition(:) = kdtree2_3dPosition(lon, lat)

    call kdtree2_r_nearest(tp=tree, qv=refPosition, r2=searchRadiusSquared, &
                           nfound=numLocalGridPointsFound, &
                           nalloc=maxNumLocalGridPointsSearch, results=searchResults)
    if (numLocalGridPointsFound > maxNumLocalGridPointsSearch) then
      call utl_abort('ocm_farFromLand: the parameter maxNumLocalGridPointsSearch must be increased')
    end if
    if (numLocalGridPointsFound == 0) then
      farFromLand = .true.
    else
      farFromLand = ( sqrt(searchResults(1)%dis) > distanceToLand )
    end if

  end function ocm_farFromLand

  !--------------------------------------------------------------------------
  ! ocm_copyMask
  !--------------------------------------------------------------------------
  subroutine ocm_copyMask(oceanMask_in, oceanMask_out, beSilent_opt)
    !
    ! :Purpose: Copy the mask data from one instance of oceanMask to
    !           another. If the destination instance is not already
    !           allocated, then this will also be done.
    !
    implicit none

    ! Arguments:
    type(struct_ocm),  intent(in)    :: oceanMask_in
    logical, optional, intent(in)    :: beSilent_opt
    type(struct_ocm),  intent(inout) :: oceanMask_out

    ! Locals:
    logical :: beSilent

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .false.
    end if

    if (.not.oceanMask_in%maskPresent .or. .not.associated(oceanMask_in%mask)) then
      if ( .not. beSilent ) write(*,*) 'ocm_copyMask: no input mask, do nothing'
      return
    end if

    if (.not.associated(oceanMask_out%mask)) then
      call ocm_allocate(oceanMask_out, oceanMask_in%hco, oceanMask_in%nLev)
    end if

    if ( .not. beSilent ) write(*,*) 'ocm_copyMask: copying over the horizontal mask'
    oceanMask_out%mask(:,:,:) = oceanMask_in%mask(:,:,:)
    oceanMask_out%maskPresent = .true.

  end subroutine ocm_copyMask

  !--------------------------------------------------------------------------
  ! ocm_communicateMask
  !--------------------------------------------------------------------------
  subroutine ocm_communicateMask(oceanMask)
    !
    ! :Purpose: Copy mask fields from task 0 to all others
    !
    implicit none

    ! Arguments:
    type(struct_ocm), intent(inout) :: oceanMask

    ! Locals:
    type(struct_hco), pointer :: hco_temp

    write(*,*) 'ocm_communicateMask: starting'

    call mmpi_bcast(oceanMask%maskPresent)
    if (.not.oceanMask%maskPresent) then
      write(*,*) 'ocm_communicateMask: mask not present, return'
      return
    end if

    call mmpi_bcast(oceanMask%nLev)

    ! special treatment of hco object since EZscintID not properly communicated
    nullify(hco_temp)
    if (mmpi_myid > 0 .and. associated(oceanMask%hco)) then
      hco_temp => oceanMask%hco
      nullify(oceanMask%hco)
    end if
    call hco_mpiBcast(oceanMask%hco)

    if (associated(hco_temp)) then
      call hco_deallocate(oceanMask%hco)
      oceanMask%hco => hco_temp
    end if

    if (.not.associated(oceanMask%mask)) then
      call ocm_allocate(oceanMask,oceanMask%hco,oceanMask%nLev)
    end if
    ! We only need to broadcast the first level 'oceanMask%mask(:,:,1)'
    call mmpi_bcast(oceanMask%mask(:,:,1))

    write(*,*) 'ocm_communicateMask: finished'

  end subroutine ocm_communicateMask

  !--------------------------------------------------------------------------
  ! ocm_allocate
  !--------------------------------------------------------------------------
  subroutine ocm_allocate(oceanMask,hco,nLev)
    !
    ! :Purpose: Allocate the object, if it isn't already.
    !
    implicit none

    ! Arguments:
    type(struct_ocm),          intent(inout) :: oceanMask
    type(struct_hco), pointer, intent(in)    :: hco
    integer,                   intent(in)    :: nLev

    if (.not.associated(oceanMask%mask)) then
      allocate(oceanMask%mask(hco%ni,hco%nj,nLev))
      oceanMask%maskPresent = .true.
      oceanMask%hco         => hco
      oceanMask%nLev        = nLev
    end if

  end subroutine ocm_allocate

  !--------------------------------------------------------------------------
  ! ocm_deallocate
  !--------------------------------------------------------------------------
  subroutine ocm_deallocate(oceanMask)
    !
    ! :Purpose: Deallocate object.
    !
    implicit none

    ! Arguments:
    type(struct_ocm), intent(inout) :: oceanMask

    if (associated(oceanMask%mask)) then
      deallocate(oceanMask%mask)
      nullify(oceanMask%mask)
      nullify(oceanMask%hco)
      oceanMask%maskPresent = .false.
      oceanMask%nLev        = 0
    end if

  end subroutine ocm_deallocate

  !--------------------------------------------------------------------------
  ! ocm_logicalToInt
  !--------------------------------------------------------------------------
  subroutine ocm_copyToInt(oceanMask, intArray, maskLev)
    !
    ! :Purpose: Convert the selected level of the logical oceanMask
    !           object into integer values where true is 1 and false is 0.
    !
    implicit none

    ! Arguments:
    type(struct_ocm), intent(inout) :: oceanMask
    integer,          intent(out)   :: intArray(:,:)
    integer,          intent(in)    :: maskLev

    intArray(:,:) = 0
    where(oceanMask%mask(:,:,maskLev)) intArray(:,:) = 1

  end subroutine ocm_copyToInt

  !--------------------------------------------------------------------------
  ! ocm_logicalToInt
  !--------------------------------------------------------------------------
  subroutine ocm_copyFromInt(oceanMask, intArray, maskLev)
    !
    ! :Purpose: Convert a 2D integer array into logical values
    !           where true is 1 and false is 0 and copy into
    !           the selected level of the oceanMask object.
    !
    implicit none

    ! Arguments:
    type(struct_ocm), intent(inout) :: oceanMask
    integer,          intent(in)    :: intArray(:,:)
    integer,          intent(in)    :: maskLev

    oceanMask%mask(:,:,maskLev) = .false.
    where(intArray(:,:) == 1) oceanMask%mask(:,:,maskLev) = .true.

  end subroutine ocm_copyFromInt

  !--------------------------------------------------------------------------
  ! ocm_computeMinGridSpacing
  !--------------------------------------------------------------------------
  subroutine ocm_computeMinGridSpacing(oceanMask, hco, minGridSpacing)
    !
    ! :Purpose: Compute minGridSpacing taking into account ocean mask.
    !           hco_setupFromFile compute it for all points, land and ocean resulting in
    !           minGridSpacing = 0, as the ORCA025 grid has the numerical pole on land
    !           in Canada around (107W, 66N).
    !
    implicit none

    ! Arguments:
    type(struct_ocm),          intent(inout) :: oceanMask
    type(struct_hco), pointer, intent(in)    :: hco
    real(8)                  , intent(out)   :: minGridSpacing ! horisontal min grid spacing in meters

    ! Locals:
    real(4), parameter :: absMaxLat = 85. ! abs of latitude threshold where to compute waterMinGridSpacing
    real(8) :: minDeltaLat, minDeltaLon
    real(8) :: deltaLon, deltaLon1, deltaLon2, deltaLon3
    real(8) :: deltaLat, deltaLat1, deltaLat2, deltaLat3
    integer :: latIndex, lonIndex

    call msg('ocm_computeMinGridSpacing', 'Computing minGridSpacing using ocean-land mask')

    minDeltaLat = 1.0d6
    do lonIndex = 1, hco%ni - 1
      do latIndex = 1, hco%nj - 1

        if (oceanMask%mask(lonIndex, latIndex, 1)) then

          deltaLat1 = abs(hco%lat2d_4(lonIndex, latIndex) - hco%lat2d_4(lonIndex    , latIndex + 1))
          deltaLat2 = abs(hco%lat2d_4(lonIndex, latIndex) - hco%lat2d_4(lonIndex + 1, latIndex    ))
          deltaLat3 = abs(hco%lat2d_4(lonIndex, latIndex) - hco%lat2d_4(lonIndex + 1, latIndex + 1))

          deltaLat = max(deltaLat1, deltaLat2, deltaLat3)
          if (deltaLat < minDeltaLat) minDeltaLat = deltaLat

        end if

      end do
    end do

    minDeltaLon = 1.0d6
    do lonIndex = 1, hco%ni - 1
      do latIndex = 1, hco%nj - 1

        if(abs(hco%lat2d_4(lonIndex, latIndex)) * MPC_DEGREES_PER_RADIAN_R8 < absMaxLat) then

          if (oceanMask%mask(lonIndex, latIndex, 1)) then

            deltaLon1 = abs(hco%lon2d_4(lonIndex, latIndex) - hco%lon2d_4(lonIndex    , latIndex + 1))
            deltaLon2 = abs(hco%lon2d_4(lonIndex, latIndex) - hco%lon2d_4(lonIndex + 1, latIndex    ))
            deltaLon3 = abs(hco%lon2d_4(lonIndex, latIndex) - hco%lon2d_4(lonIndex + 1, latIndex + 1))

            if (deltaLon1 > MPC_PI_R8) deltaLon1 = deltaLon1 - 2.0d0 * MPC_PI_R8
            deltaLon1 = abs(deltaLon1 * cos(hco%lat2d_4(lonIndex,latIndex)))
            if (deltaLon2 > MPC_PI_R8) deltaLon2 = deltaLon2 - 2.0d0 * MPC_PI_R8
            deltaLon2 = abs(deltaLon2 * cos(hco%lat2d_4(lonIndex,latIndex)))
            if (deltaLon3 > MPC_PI_R8) deltaLon3 = deltaLon3 - 2.0d0 * MPC_PI_R8
            deltaLon3 = abs(deltaLon3 * cos(hco%lat2d_4(lonIndex,latIndex)))

            deltaLon = max(deltaLon1, deltaLon2, deltaLon3)
            if (deltaLon < minDeltaLon) minDeltaLon = deltaLon

          end if

        end if

      end do
    end do

    minGridSpacing = ec_ra * sqrt(2.0d0) * min(minDeltaLon, minDeltaLat)

    if ( .not. utl_isEqual(minGridSpacing, hco%minGridSpacing) ) then
      call msg('ocm_computeMinGridSpacing', &
               'minDeltaLat= '//str(minDeltaLat * MPC_DEGREES_PER_RADIAN_R8)//' deg')
      call msg('ocm_computeMinGridSpacing', &
               'minDeltaLon= '//str(minDeltaLon * MPC_DEGREES_PER_RADIAN_R8)//' deg')
      call msg('ocm_computeMinGridSpacing', &
               'Updated ocean water points minGridSpacing: '//str(hco%minGridSpacing)//' m')
    end if

    if (minGridSpacing > 1.0d6) then
      call utl_abort('ocm_computeMinGridSpacing: minGridSpacing is greater than 1000 km.')
    end if

  end subroutine ocm_computeMinGridSpacing

end module oceanMask_mod

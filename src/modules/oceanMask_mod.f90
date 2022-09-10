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

module oceanMask_mod
  ! MODULE oceanMask_mod (prefix='ocm' category='7. Low-level data objects')
  !
  ! :Purpose: The horizontal mask indicating land (=.false.) and water (=.true.) grid points.
  !           This mask is either:
  !                 * In the case of variables on ocean depth levels, it varies with vertical level.
  !                 * In other cases it is a single 2D field used for all variables.
  !
  use midasMpi_mod
  use kdtree2_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use utilities_mod
  implicit none
  save
  private

  ! public structure definition
  public :: struct_ocm

  ! public subroutines and functions
  public :: ocm_readMaskFromFile, ocm_deallocate
  public :: ocm_copyMask, ocm_communicateMask
  public :: ocm_farFromLand
  public :: ocm_copyToInt, ocm_copyFromInt

  type struct_ocm
    ! This is the derived type of the ocean mask object
    integer                   :: nLev
    logical, pointer          :: mask(:,:,:) => null()
    logical                   :: maskPresent = .false.
    type(struct_hco), pointer :: hco
  end type struct_ocm

  integer, external  :: get_max_rss

  contains

  !--------------------------------------------------------------------------
  ! ocm_readMaskFromFile
  !--------------------------------------------------------------------------
  subroutine ocm_readMaskFromFile(oceanMask, hco, vco, filename)
    !
    ! :Purpose: Check if any mask fields exist for surface or ocean depth levels.
    !
    !
    implicit none

    ! arguments:
    type(struct_ocm),          intent(inout) :: oceanMask
    type(struct_hco), pointer, intent(inout) :: hco
    type(struct_vco),          intent(in)    :: vco
    character(len=*),          intent(in)    :: fileName

    ! locals:
    integer :: nulfile, ierr, ip1, ni_file, nj_file, nk_file
    integer :: ikey, levIndex
    integer :: fnom, fstouv, fclos, fstfrm, fstluk, fstinf, fstsui
    integer, allocatable :: mask(:,:)
    integer :: maxkeys

    ! Check if any mask is present in file, return if not
    if ( .not. utl_varNamePresentInFile(' ',fileName_opt=trim(fileName),typvar_opt='@@') ) then
      return
    end if

    !- Open input file
    nulfile = 0
    write(*,*) 'ocm_readMaskFromFile: file name = ',trim(fileName)
    ierr = fnom(nulfile,trim(fileName),'RND+OLD+R/O',0)

    if ( ierr >= 0 ) then
      maxkeys = fstouv(nulfile,'RND+OLD')
    else
      call utl_abort('ocm_readMaskFromFile: problem opening input file')
    end if

    if (nulfile == 0 ) then
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
          write(*,*) 'ocm_readMaskFromFile: searched for mask with ip1 = ', ip1
          call utl_abort('ocm_readMaskFromFile: cannot find mask for this ip1 in file ' // trim(fileName))
        end if

        do while (ni_file /= hco%ni .or. nj_file /= hco%nj)
          ikey = fstsui(nulfile, ni_file, nj_file, nk_file)
        end do

        if (ni_file == hco%ni .and. nj_file == hco%nj) then
          write(*,*) 'ocm_readMaskFromFile: read mask for ip1 = ', ip1
          if (.not. associated(oceanMask%mask)) then
            call ocm_allocate(oceanMask,hco,vco%nLev_depth)
          end if
          if (.not. allocated(mask)) allocate(mask(hco%ni,hco%nj))
          ierr = fstluk(mask(:,:), ikey, ni_file, nj_file, nk_file)
          call ocm_copyFromInt(oceanMask,mask,levIndex)
          if (ierr < 0) then
            call utl_abort('ocm_readMaskFromFile: error when reading mask record')
          end if
        else
          ! Special cases for variables that are on a different horizontal grid in LAM (e.g. TG)
          write(*,*)
          write(*,*) 'ocm_readMaskFromFile: mask is on a different horizontal grid'
          write(*,*) 'ni = ', ni_file, hco%ni, ', nj = ', nj_file, hco%nj
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
        write(*,*) 'ocm_readMaskFromFile: reading mask'
        if (.not. associated(oceanMask%mask)) then
          call ocm_allocate(oceanMask,hco,1)
          if (.not. allocated(mask)) allocate(mask(hco%ni,hco%nj))
          ierr = fstluk(mask(:,:), ikey, ni_file, nj_file, nk_file)
          call ocm_copyFromInt(oceanMask,mask,1)
          if (ierr < 0) then
            call utl_abort('ocm_readMaskFromFile: error when reading mask record')
          end if
        end if
      else
        ! Special cases for variables that are on a different horizontal grid in LAM (e.g. TG)
        write(*,*)
        write(*,*) 'ocm_readMaskFromFile: mask is on a different horizontal grid'
        write(*,*) 'ni = ', ni_file, hco%ni, ', nj = ', nj_file, hco%nj
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

    ! arguments:
    type(struct_ocm), intent(in) :: oceanMask
    integer,          intent(in) :: levIndex
    real(8),          intent(in) :: lon
    real(8),          intent(in) :: lat
    real(8),          intent(in) :: distanceToLand
    logical                      :: farFromLand

    ! locals:
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
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

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
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

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
  subroutine ocm_copyMask(oceanMask_in,oceanMask_out)
    !
    ! :Purpose: Copy the mask data from one instance of oceanMask to
    !           another. If the destination instance is not already
    !           allocated, then this will also be done.
    !
    implicit none

    ! arguments:
    type(struct_ocm), intent(in)    :: oceanMask_in
    type(struct_ocm), intent(inout) :: oceanMask_out

    if (.not.oceanMask_in%maskPresent .or. .not.associated(oceanMask_in%mask)) then
      write(*,*) 'ocm_copyMask: no input mask, do nothing'
      return
    end if

    if (.not.associated(oceanMask_out%mask)) then
      call ocm_allocate(oceanMask_out, oceanMask_in%hco, oceanMask_in%nLev)
    end if

    write(*,*) 'ocm_copyMask: copying over the horizontal mask'
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

    ! arguments:
    type(struct_ocm), intent(inout) :: oceanMask

    ! locals:
    integer                   :: ierr
    type(struct_hco), pointer :: hco_temp

    write(*,*) 'ocm_communicateMask: starting'

    call rpn_comm_bcast(oceanMask%maskPresent, 1,  &
                        'MPI_LOGICAL', 0, 'GRID', ierr)
    if (.not.oceanMask%maskPresent) then
      write(*,*) 'ocm_communicateMask: mask not present, return'
      return
    end if
    
    call rpn_comm_bcast(oceanMask%nLev, 1,  &
                        'MPI_INTEGER', 0, 'GRID', ierr)

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
    call rpn_comm_bcast(oceanMask%mask, oceanMask%hco%ni*oceanMask%hco%nj*1,  &
                        'MPI_LOGICAL', 0, 'GRID', ierr)

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

    ! arguments:
    type(struct_ocm),          intent(inout) :: oceanMask
    type(struct_hco), pointer, intent(inout) :: hco
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

    ! arguments:
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

    ! arguments:
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

    ! arguments:
    type(struct_ocm), intent(inout) :: oceanMask
    integer,          intent(in)    :: intArray(:,:)
    integer,          intent(in)    :: maskLev

    oceanMask%mask(:,:,maskLev) = .false.
    where(intArray(:,:) == 1) oceanMask%mask(:,:,maskLev) = .true.
    
  end subroutine ocm_copyFromInt

end module oceanMask_mod

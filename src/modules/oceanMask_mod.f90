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
  ! MODULE oceanMask_mod (prefix='ocm' category='7. Low-level data objects and utilities')
  !
  ! :Purpose: The horizontal mask indicating land (=.false.) and water (=.true.) grid points.
  !           This mask is either:
  !                 * In the case of variables on ocean depth levels, it varies with vertical level.
  !                 * In other cases it is a single 2D field used for all variables.
  !
  use horizontalCoord_mod
  use verticalCoord_mod
  use utilities_mod
  implicit none
  save
  private

  ! public structure definition
  public :: struct_ocm

  ! public subroutines and functions
  public :: ocm_readMaskFromFile, ocm_allocate, ocm_deallocate
  public :: ocm_copyMask, ocm_communicateMask
  public :: ocm_copyToInt, ocm_copyFromInt

  type struct_ocm
    ! This is the derived type of the ocean mask object
    integer             :: ni
    integer             :: nj
    integer             :: nLev
    logical, pointer    :: mask(:,:,:) => null()
    logical             :: maskPresent      = .false.
  end type struct_ocm

  contains

  !--------------------------------------------------------------------------
  ! ocm_readMaskFromFile
  !--------------------------------------------------------------------------
  subroutine ocm_readMaskFromFile(oceanMask, hco, vco, filename)
    !
    ! :Purpose: Check if any mask fields exist for surface or ocean depth levels.
    !
    ! :Note:    This is a temporary version of the subroutine that only
    !           reads the first mask found. Eventually we will need to
    !           store masks separately for each level.
    !
    implicit none

    ! arguments
    type(struct_ocm)              :: oceanMask
    type(struct_hco)              :: hco
    type(struct_vco)              :: vco
    character(len=*), intent(in)  :: fileName

    ! locals
    integer :: nulfile, ierr, ip1, ni_file, nj_file, nk_file
    integer :: ikey, levIndex
    integer :: fnom, fstouv, fclos, fstfrm, fstlir, fstinf
    integer, allocatable :: mask(:,:)

    ! Check if any mask is present in file, return if not
    if ( .not. utl_varNamePresentInFile(' ',fileName_opt=trim(fileName),typvar_opt='@@') ) then
      return
    end if

    !- Open input field
    nulfile = 0
    write(*,*) 'ocm_readMaskFromFile: file name = ',trim(fileName)
    ierr = fnom(nulfile,trim(fileName),'RND+OLD+R/O',0)
       
    if ( ierr >= 0 ) then
      ierr  =  fstouv(nulfile,'RND+OLD')
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

        if (ni_file == hco%ni .and. nj_file == hco%nj) then
          write(*,*) 'ocm_readMaskFromFile: read mask for ip1 = ', ip1
          if (.not. associated(oceanMask%mask)) then
            call ocm_allocate(oceanMask,hco%ni,hco%nj,vco%nLev_depth)
          end if
          if (.not.allocated(mask)) allocate(mask(hco%ni,hco%nj))
          ierr = fstlir(mask(:,:), nulfile, ni_file, nj_file, nk_file,  &
                        -1, ' ', ip1, -1, -1, '@@', ' ')
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

      if (ni_file == hco%ni .and. nj_file == hco%nj) then
        write(*,*) 'ocm_readMaskFromFile: reading mask'
        if (.not. associated(oceanMask%mask)) then
          call ocm_allocate(oceanMask,hco%ni,hco%nj,1)
          if (.not.allocated(mask)) allocate(mask(hco%ni,hco%nj))
          ierr = fstlir(mask(:,:), nulfile, ni_file, nj_file, nk_file,  &
                        -1, ' ', -1, -1, -1, '@@', ' ')
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
  ! ocm_copyMask
  !--------------------------------------------------------------------------
  subroutine ocm_copyMask(oceanMask_in,oceanMask_out)
    implicit none
    ! arguments
    type(struct_ocm)  :: oceanMask_in, oceanMask_out

    if (.not.oceanMask_in%maskPresent .or. .not.associated(oceanMask_in%mask)) then
      write(*,*) 'ocm_copyMask: no input mask, do nothing'
      return
    end if

    if (.not.associated(oceanMask_out%mask)) then
      call ocm_allocate(oceanMask_out, oceanMask_in%ni, oceanMask_in%nj, oceanMask_in%nLev)
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

    ! arguments
    type(struct_ocm) :: oceanMask
    ! locals
    integer :: ierr

    write(*,*) 'ocm_communicateMask: starting'

    call rpn_comm_bcast(oceanMask%maskPresent, 1,  &
                        'MPI_LOGICAL', 0, 'GRID', ierr)
    call rpn_comm_bcast(oceanMask%ni,   1,  &
                        'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(oceanMask%nj,   1,  &
                        'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(oceanMask%nLev, 1,  &
                        'MPI_INTEGER', 0, 'GRID', ierr)

    if (.not.oceanMask%maskPresent) then
      write(*,*) 'ocm_communicateMask: mask not present, return'
      return
    end if
    
    if (.not.associated(oceanMask%mask)) then
      call ocm_allocate(oceanMask,oceanMask%ni,oceanMask%nj,oceanMask%nLev)
    end if
    call rpn_comm_bcast(oceanMask%mask, oceanMask%ni*oceanMask%nj*1,  &
                        'MPI_LOGICAL', 0, 'GRID', ierr)

    write(*,*) 'ocm_communicateMask: finished'

  end subroutine ocm_communicateMask

  !--------------------------------------------------------------------------
  ! ocm_allocate
  !--------------------------------------------------------------------------
  subroutine ocm_allocate(oceanMask,ni,nj,nLev)
    ! :Purpose: Allocate object.

    implicit none

    ! Arguments:
    type(struct_ocm)     :: oceanMask
    integer              :: ni
    integer              :: nj
    integer              :: nLev

    if (.not.associated(oceanMask%mask)) then
      allocate(oceanMask%mask(ni,nj,nLev))
      oceanMask%maskPresent = .true.
      oceanMask%ni          = ni
      oceanMask%nj          = nj
      oceanMask%nLev        = nLev
    end if
    
  end subroutine ocm_allocate

  !--------------------------------------------------------------------------
  ! ocm_deallocate
  !--------------------------------------------------------------------------
  subroutine ocm_deallocate(oceanMask)
    ! :Purpose: Deallocate object.

    implicit none

    ! Arguments:
    type(struct_ocm)     :: oceanMask

    if (associated(oceanMask%mask)) then
      deallocate(oceanMask%mask)
      nullify(oceanMask%mask)
      oceanMask%maskPresent = .false.
      oceanMask%ni          = 0
      oceanMask%nj          = 0
      oceanMask%nLev        = 0
    end if
    
  end subroutine ocm_deallocate

  !--------------------------------------------------------------------------
  ! ocm_logicalToInt
  !--------------------------------------------------------------------------
  subroutine ocm_copyToInt(oceanMask, intArray, maskLev)
    ! :Purpose: Convert a 2D logical array into integer values
    !           where true is 1 and false is 0.

    implicit none

    ! Arguments:
    type(struct_ocm)     :: oceanMask
    integer, intent(out) :: intArray(:,:)
    integer, intent(in)  :: maskLev

    intArray(:,:) = 0
    where(oceanMask%mask(:,:,maskLev)) intArray(:,:) = 1
    
  end subroutine ocm_copyToInt

  !--------------------------------------------------------------------------
  ! ocm_logicalToInt
  !--------------------------------------------------------------------------
  subroutine ocm_copyFromInt(oceanMask, intArray, maskLev)
    ! :Purpose: Convert a 2D integer array into logical values
    !           where true is 1 and false is 0.

    implicit none

    ! Arguments:
    type(struct_ocm)    :: oceanMask
    integer, intent(in) :: intArray(:,:)
    integer, intent(in) :: maskLev

    oceanMask%mask(:,:,maskLev) = .false.
    where(intArray(:,:) == 1) oceanMask%mask(:,:,maskLev) = .true.
    
  end subroutine ocm_copyFromInt

end module oceanMask_mod

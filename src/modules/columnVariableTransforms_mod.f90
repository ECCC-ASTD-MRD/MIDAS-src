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

module columnVariableTransforms_mod
  ! MODULE columnVariableTransforms (prefix='cvt' category='3. High-level transformations')
  !
  ! :Purpose: To store various functions for variable transforms using inputs
  !           from column(s). Outputs are also placed in a column.
  !
  use mpivar_mod
  use mathPhysConstants_mod
  use earthConstants_mod
  use timeCoord_mod
  use columnData_mod
  use verticalCoord_mod
  use utilities_mod
  use varNameList_mod
  use tt2phi_mod
  use utilities_mod
  
  implicit none
  save
  private

  ! public procedures
  public :: cvt_transform

CONTAINS

  !--------------------------------------------------------------------------
  ! cvt_transform
  !--------------------------------------------------------------------------
  subroutine cvt_transform(column, columnRef, transform)
    implicit none
   
    ! Arguments
    type(struct_columnData)      :: column
    type(struct_columnData)      :: columnRef
    character(len=*), intent(in) :: transform

    select case(trim(transform))

    case ('TTHUtoHeight_tl')
      if ( .not. col_varExist(column,'TT')  ) then
        call utl_abort('cvt_transform: for TTHUtoHeight_tl, variable TT must be allocated in column')
      end if
      if ( .not. col_varExist(column,'HU')  ) then
        call utl_abort('cvt_transform: for TTHUtoHeight_tl, variable HU must be allocated in column')
      end if
      if ( .not. col_varExist(column,'P0')  ) then
        call utl_abort('cvt_transform: for TTHUtoHeight_tl, variable P0 must be allocated in column')
      end if
      if ( .not. col_varExist(column,'Z_T')  ) then
        call utl_abort('cvt_transform: for TTHUtoHeight_tl, variable Z_T must be allocated in column')
      end if
      if ( .not. col_varExist(column,'Z_M')  ) then
        call utl_abort('cvt_transform: for TTHUtoHeight_tl, variable Z_M must be allocated in column')
      end if
      call TTHUtoHeight_tl(column,columnRef)

    case ('TTHUtoHeight_ad')
      if ( .not. col_varExist(column,'TT')  ) then
        call utl_abort('cvt_transform: for TTHUtoHeight_ad, variable TT must be allocated in column')
      end if
      if ( .not. col_varExist(column,'HU')  ) then
        call utl_abort('cvt_transform: for TTHUtoHeight_ad, variable HU must be allocated in column')
      end if
      if ( .not. col_varExist(column,'P0')  ) then
        call utl_abort('cvt_transform: for TTHUtoHeight_ad, variable P0 must be allocated in column')
      end if
      if ( .not. col_varExist(column,'Z_T')  ) then
        call utl_abort('cvt_transform: for TTHUtoHeight_ad, variable Z_T must be allocated in column')
      end if
      if ( .not. col_varExist(column,'Z_M')  ) then
        call utl_abort('cvt_transform: for TTHUtoHeight_ad, variable Z_M must be allocated in column')
      end if
      call TTHUtoHeight_ad(column,columnRef)

    case ('PsfcToP_tl')
      if ( .not. col_varExist(column,'P_T')  ) then
        call utl_abort('cvt_transform: for PsfcToP_tl, variable P_T must be allocated in column')
      end if
      if ( .not. col_varExist(column,'P_M')  ) then
        call utl_abort('cvt_transform: for PsfcToP_tl, variable P_M must be allocated in column')
      end if
      if ( .not. col_varExist(column,'P0')  ) then
        call utl_abort('cvt_transform: for PsfcToP_tl, variable P0 must be allocated in column')
      end if
      call PsfcToP_tl(column,columnRef)

    case ('PsfcToP_ad')
      if ( .not. col_varExist(column,'P_T')  ) then
        call utl_abort('cvt_transform: for PsfcToP_ad, variable P_T must be allocated in column')
      end if
      if ( .not. col_varExist(column,'P_M')  ) then
        call utl_abort('cvt_transform: for PsfcToP_ad, variable P_M must be allocated in column')
      end if
      if ( .not. col_varExist(column,'P0')  ) then
        call utl_abort('cvt_transform: for PsfcToP_ad, variable P0 must be allocated in column')
      end if
      call PsfcToP_ad(column,columnRef)

    case default
      write(*,*)
      write(*,*) 'Unsupported function : ', trim(transform)
      call utl_abort('cvt_transform')
    end select

  end subroutine cvt_transform

  !--------------------------------------------------------------------------
  ! TTHUtoHeight_tl
  !--------------------------------------------------------------------------
  subroutine TTHUtoHeight_tl(column,columnRef)
    implicit none

    type(struct_columnData)    :: column, columnRef

    call tt2phi_tl(column, columnRef)

  end subroutine TTHUtoHeight_tl

  !--------------------------------------------------------------------------
  ! TTHUtoHeight_ad
  !--------------------------------------------------------------------------
  subroutine TTHUtoHeight_ad(column,columnRef)
    implicit none

    type(struct_columnData)    :: column, columnRef

    call tt2phi_ad(column,columnRef)

  end subroutine TTHUtoHeight_ad

  !--------------------------------------------------------------------------
  ! PsfcToP_tl
  !--------------------------------------------------------------------------
  subroutine PsfcToP_tl(column,columnRef)
    implicit none

    type(struct_columnData)    :: column, columnRef

    call cvt_calcPressure_tl(column,columnRef)

  end subroutine PsfcToP_tl

  !--------------------------------------------------------------------------
  ! PsfcToP_ad
  !--------------------------------------------------------------------------
  subroutine PsfcToP_ad(column,columnRef)
    implicit none

    type(struct_columnData)    :: column, columnRef

    call cvt_calcPressure_ad(column,columnRef)

  end subroutine PsfcToP_ad

  !--------------------------------------------------------------------------
  ! cvt_calcPressure_tl
  !--------------------------------------------------------------------------
  subroutine cvt_calcPressure_tl(column, columnRef, beSilent_opt)
    !
    !:Purpose: calculation of the pressure increment on the grid.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(inout) :: column    ! column that will contain the P_T/P_M increments
    type(struct_columnData), intent(in)    :: columnRef ! column that has the Psfc
    logical, optional                      :: beSilent_opt

    ! Locals:
    real(8)          :: Psfc
    real(8), pointer :: delPsfc(:,:), PsfcRef(:,:)
    real(8), pointer :: delP_T(:,:), delP_M(:,:)
    real(8), pointer :: dP_dPsfc_T(:), dP_dPsfc_M(:)
    integer          :: status, colIndex
    integer          :: lev_M, lev_T, nlev_T, nlev_M, numColumns
    logical          :: beSilent

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .true.
    end if

    if ( .not. beSilent ) write(*,*) 'cvt_calcPressure_tl: computing delP_T/delP_M on the column'

    nullify(dP_dPsfc_T)
    nullify(dP_dPsfc_M)
    nullify(delPsfc)
    nullify(delP_T)
    nullify(delP_M)

    delP_M  => col_getAllColumns(column,'P_M')
    delP_T  => col_getAllColumns(column,'P_T')
    delPsfc => col_getAllColumns(column,'P0')
    PsfcRef => col_getAllColumns(columnRef,'P0')

    nlev_T = col_getNumLev(column,'TH')
    nlev_M = col_getNumLev(column,'MM')
    numColumns = col_getNumCol(column)

    do colIndex = 1, numColumns

      Psfc = PsfcRef(1,colIndex)

      ! dP_dPsfc_M
      nullify(dP_dPsfc_M)
      status = vgd_dpidpis(column%vco%vgrid, &
                           column%vco%ip1_M, &
                           dP_dPsfc_M, &
                           Psfc)
      if( status .ne. VGD_OK ) call utl_abort('cvt_calcPressure_tl: ERROR with vgd_dpidpis')
      ! calculate delP_M
      do lev_M = 1, nlev_M
        delP_M(lev_M,colIndex) = dP_dPsfc_M(lev_M) * delPsfc(1,colIndex)
      end do
      deallocate(dP_dPsfc_M)

      ! dP_dPsfc_T
      nullify(dP_dPsfc_T)
      status = vgd_dpidpis(column%vco%vgrid, &
                           column%vco%ip1_T, &
                           dP_dPsfc_T, &
                           Psfc)
      if( status .ne. VGD_OK ) call utl_abort('cvt_calcPressure_tl: ERROR with vgd_dpidpis')
      ! calculate delP_T
      do lev_T = 1, nlev_T
        delP_T(lev_T,colIndex) = dP_dPsfc_T(lev_T) * delPsfc(1,colIndex)
      end do
      deallocate(dP_dPsfc_T)

    end do

  end subroutine cvt_calcPressure_tl

  !--------------------------------------------------------------------------
  ! cvt_calcPressure_ad
  !--------------------------------------------------------------------------
  subroutine cvt_calcPressure_ad(column, columnRef, beSilent_opt)
    !
    !:Purpose: adjoint of calculation of the pressure on the grid.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(inout) :: column    ! column that will contain ncrement of Psfc.
    type(struct_columnData), intent(in)    :: columnRef ! column that has the Psfc
    logical, optional                      :: beSilent_opt

    ! Locals:
    real(8)          :: Psfc
    real(8), pointer :: delPsfc(:,:), PsfcRef(:,:)
    real(8), pointer :: delP_T(:,:), delP_M(:,:)
    real(8), pointer :: dP_dPsfc_T(:), dP_dPsfc_M(:)
    integer          :: status, colIndex
    integer          :: lev_M, lev_T, nlev_T, nlev_M, numColumns
    logical          :: beSilent

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .true.
    end if

    if ( .not. beSilent ) write(*,*) 'cvt_calcPressure_ad: computing delP_T/delP_M on the column'

    nullify(delPsfc)
    nullify(PsfcRef)
    nullify(delP_T)
    nullify(delP_M)
    nullify(dP_dPsfc_T)
    nullify(dP_dPsfc_M)

    delP_M  => col_getAllColumns(column,'P_M')
    delP_T  => col_getAllColumns(column,'P_T')
    delPsfc => col_getAllColumns(column,'P0')
    PsfcRef => col_getAllColumns(columnRef,'P0')

    nlev_T = col_getNumLev(column,'TH')
    nlev_M = col_getNumLev(column,'MM')
    numColumns = col_getNumCol(column)

    do colIndex = 1, numColumns

      Psfc = PsfcRef(1,colIndex)

      ! dP_dPsfc_M
      nullify(dP_dPsfc_M)
      status = vgd_dpidpis(column%vco%vgrid, &
                           column%vco%ip1_M, &
                           dP_dPsfc_M, &
                           Psfc)
      if( status .ne. VGD_OK ) call utl_abort('cvt_calcPressure_ad: ERROR with vgd_dpidpis')
      ! calculate delP_M
      do lev_M = 1, nlev_M
        delPsfc(1,colIndex) = delPsfc(1,colIndex) + &
             dP_dPsfc_M(lev_M) * delP_M(lev_M,colIndex)
      end do
      deallocate(dP_dPsfc_M)

      ! dP_dPsfc_T
      nullify(dP_dPsfc_T)
      status = vgd_dpidpis(column%vco%vgrid, &
                           column%vco%ip1_T, &
                           dP_dPsfc_T, &
                           Psfc)
      if( status .ne. VGD_OK ) call utl_abort('cvt_calcPressure_ad: ERROR with vgd_dpidpis')
      ! calculate delP_T
      do lev_T = 1, nlev_T
        delPsfc(1,colIndex) = delPsfc(1,colIndex) + dP_dPsfc_T(lev_T) * delP_T(lev_T,colIndex)
      end do
      deallocate(dP_dPsfc_T)

    end do

  end subroutine cvt_calcPressure_ad

end module columnVariableTransforms_mod

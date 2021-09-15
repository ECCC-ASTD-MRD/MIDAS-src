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
  use calcHeightAndPressure_mod
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
  subroutine cvt_transform(columnInc, columnRefOnIncLev, transform)
    implicit none
   
    ! Arguments
    type(struct_columnData)      :: columnInc
    type(struct_columnData)      :: columnRefOnIncLev
    character(len=*), intent(in) :: transform

    select case(trim(transform))

    case ('TTHUtoHeight_tl')
      if ( .not. col_varExist(columnInc,'TT')  ) then
        call utl_abort('cvt_transform: for TTHUtoHeight_tl, variable TT must be allocated in columnInc')
      end if
      if ( .not. col_varExist(columnInc,'HU')  ) then
        call utl_abort('cvt_transform: for TTHUtoHeight_tl, variable HU must be allocated in columnInc')
      end if
      if ( .not. col_varExist(columnInc,'P0')  ) then
        call utl_abort('cvt_transform: for TTHUtoHeight_tl, variable P0 must be allocated in columnInc')
      end if
      if ( .not. col_varExist(columnInc,'Z_T')  ) then
        call utl_abort('cvt_transform: for TTHUtoHeight_tl, variable Z_T must be allocated in columnInc')
      end if
      if ( .not. col_varExist(columnInc,'Z_M')  ) then
        call utl_abort('cvt_transform: for TTHUtoHeight_tl, variable Z_M must be allocated in columnInc')
      end if
      call TTHUtoHeight_tl(columnInc,columnRefOnIncLev)

    case ('TTHUtoHeight_ad')
      if ( .not. col_varExist(columnInc,'TT')  ) then
        call utl_abort('cvt_transform: for TTHUtoHeight_ad, variable TT must be allocated in columnInc')
      end if
      if ( .not. col_varExist(columnInc,'HU')  ) then
        call utl_abort('cvt_transform: for TTHUtoHeight_ad, variable HU must be allocated in columnInc')
      end if
      if ( .not. col_varExist(columnInc,'P0')  ) then
        call utl_abort('cvt_transform: for TTHUtoHeight_ad, variable P0 must be allocated in columnInc')
      end if
      if ( .not. col_varExist(columnInc,'Z_T')  ) then
        call utl_abort('cvt_transform: for TTHUtoHeight_ad, variable Z_T must be allocated in columnInc')
      end if
      if ( .not. col_varExist(columnInc,'Z_M')  ) then
        call utl_abort('cvt_transform: for TTHUtoHeight_ad, variable Z_M must be allocated in columnInc')
      end if
      call TTHUtoHeight_ad(columnInc,columnRefOnIncLev)

    case ('PsfcToP_tl')
      if ( .not. col_varExist(columnInc,'P_T')  ) then
        call utl_abort('cvt_transform: for PsfcToP_tl, variable P_T must be allocated in columnInc')
      end if
      if ( .not. col_varExist(columnInc,'P_M')  ) then
        call utl_abort('cvt_transform: for PsfcToP_tl, variable P_M must be allocated in columnInc')
      end if
      if ( .not. col_varExist(columnInc,'P0')  ) then
        call utl_abort('cvt_transform: for PsfcToP_tl, variable P0 must be allocated in columnInc')
      end if
      call PsfcToP_tl(columnInc,columnRefOnIncLev)

    case ('PsfcToP_ad')
      if ( .not. col_varExist(columnInc,'P_T')  ) then
        call utl_abort('cvt_transform: for PsfcToP_ad, variable P_T must be allocated in columnInc')
      end if
      if ( .not. col_varExist(columnInc,'P_M')  ) then
        call utl_abort('cvt_transform: for PsfcToP_ad, variable P_M must be allocated in columnInc')
      end if
      if ( .not. col_varExist(columnInc,'P0')  ) then
        call utl_abort('cvt_transform: for PsfcToP_ad, variable P0 must be allocated in columnInc')
      end if
      call PsfcToP_ad(columnInc,columnRefOnIncLev)

    case default
      write(*,*)
      write(*,*) 'Unsupported function : ', trim(transform)
      call utl_abort('cvt_transform')
    end select

  end subroutine cvt_transform

  !--------------------------------------------------------------------------
  ! TTHUtoHeight_tl
  !--------------------------------------------------------------------------
  subroutine TTHUtoHeight_tl(columnInc,columnRefOnIncLev)
    implicit none

    type(struct_columnData)    :: columnInc, columnRefOnIncLev

    call czp_tt2phi_tl(columnInc, columnRefOnIncLev)

  end subroutine TTHUtoHeight_tl

  !--------------------------------------------------------------------------
  ! TTHUtoHeight_ad
  !--------------------------------------------------------------------------
  subroutine TTHUtoHeight_ad(columnInc,columnRefOnIncLev)
    implicit none

    type(struct_columnData)    :: columnInc, columnRefOnIncLev

    call czp_tt2phi_ad(columnInc,columnRefOnIncLev)

  end subroutine TTHUtoHeight_ad

  !--------------------------------------------------------------------------
  ! PsfcToP_tl
  !--------------------------------------------------------------------------
  subroutine PsfcToP_tl(columnInc,columnRefOnIncLev)
    implicit none

    type(struct_columnData)    :: columnInc, columnRefOnIncLev

    call czp_calcColumnPressure_tl(columnInc,columnRefOnIncLev)

  end subroutine PsfcToP_tl

  !--------------------------------------------------------------------------
  ! PsfcToP_ad
  !--------------------------------------------------------------------------
  subroutine PsfcToP_ad(columnInc,columnRefOnIncLev)
    implicit none

    type(struct_columnData)    :: columnInc, columnRefOnIncLev

    call czp_calcColumnPressure_ad(columnInc,columnRefOnIncLev)

  end subroutine PsfcToP_ad

end module columnVariableTransforms_mod

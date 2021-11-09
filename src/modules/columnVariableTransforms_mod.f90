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
  subroutine cvt_transform(columnInc, transform, columnRefOnIncLev_opt)
    implicit none
   
    ! Arguments
    type(struct_columnData), intent(inout)        :: columnInc
    character(len=*), intent(in)                  :: transform
    type(struct_columnData), optional, intent(in) :: columnRefOnIncLev_opt
    
    select case(trim(transform))

    case ('ZandP_tl')
      if (.not. present(columnRefOnIncLev_opt)) then
        call utl_abort('cvt_transform: columnRefOnIncLev_opt required')
      end if
      call czp_calcZandP_tl(columnInc, columnRefOnIncLev_opt)

    case ('ZandP_ad')
      if (.not. present(columnRefOnIncLev_opt)) then
        call utl_abort('cvt_transform: columnRefOnIncLev_opt required')
      end if
      call czp_calcZandP_ad(columnInc, columnRefOnIncLev_opt)

    case ('PsfcToP_nl')
      call czp_calcPressure_nl(columnInc)

    case ('PsfcToP_tl')
      if (.not. present(columnRefOnIncLev_opt)) then
        call utl_abort('cvt_transform: columnRefOnIncLev_opt required')
      end if
      if ( .not. col_varExist(columnInc,'P_T')  ) then
        call utl_abort('cvt_transform: for PsfcToP_tl, variable P_T must be allocated in columnInc')
      end if
      if ( .not. col_varExist(columnInc,'P_M')  ) then
        call utl_abort('cvt_transform: for PsfcToP_tl, variable P_M must be allocated in columnInc')
      end if
      if ( .not. col_varExist(columnInc,'P0')  ) then
        call utl_abort('cvt_transform: for PsfcToP_tl, variable P0 must be allocated in columnInc')
      end if
      call PsfcToP_tl(columnInc,columnRefOnIncLev_opt)

    case ('PsfcToP_ad')
      if (.not. present(columnRefOnIncLev_opt)) then
        call utl_abort('cvt_transform: columnRefOnIncLev_opt required')
      end if
      if ( .not. col_varExist(columnInc,'P_T')  ) then
        call utl_abort('cvt_transform: for PsfcToP_ad, variable P_T must be allocated in columnInc')
      end if
      if ( .not. col_varExist(columnInc,'P_M')  ) then
        call utl_abort('cvt_transform: for PsfcToP_ad, variable P_M must be allocated in columnInc')
      end if
      if ( .not. col_varExist(columnInc,'P0')  ) then
        call utl_abort('cvt_transform: for PsfcToP_ad, variable P0 must be allocated in columnInc')
      end if
      call PsfcToP_ad(columnInc,columnRefOnIncLev_opt)

    case default
      write(*,*)
      write(*,*) 'Unsupported function : ', trim(transform)
      call utl_abort('cvt_transform')
    end select

  end subroutine cvt_transform

  !--------------------------------------------------------------------------
  ! PsfcToP_tl
  !--------------------------------------------------------------------------
  subroutine PsfcToP_tl(columnInc,columnRefOnIncLev)
    implicit none

    type(struct_columnData)    :: columnInc, columnRefOnIncLev

    call czp_calcPressure_tl(columnInc,columnRefOnIncLev)

  end subroutine PsfcToP_tl

  !--------------------------------------------------------------------------
  ! PsfcToP_ad
  !--------------------------------------------------------------------------
  subroutine PsfcToP_ad(columnInc,columnRefOnIncLev)
    implicit none

    type(struct_columnData)    :: columnInc, columnRefOnIncLev

    call czp_calcPressure_ad(columnInc,columnRefOnIncLev)

  end subroutine PsfcToP_ad

end module columnVariableTransforms_mod

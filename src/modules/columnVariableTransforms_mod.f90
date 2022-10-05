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
  ! MODULE columnVariableTransforms (prefix='cvt' category='4. Data Object transformations')
  !
  ! :Purpose: To store various functions for variable transforms using inputs
  !           from column(s). Outputs are also placed in a column.
  !
  use columnData_mod
  use utilities_mod
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

    case ('ZandP_nl')
      call czp_calcZandP_nl(columnInc)

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

    case default
      write(*,*)
      write(*,*) 'Unsupported function : ', trim(transform)
      call utl_abort('cvt_transform')
    end select

  end subroutine cvt_transform

end module columnVariableTransforms_mod

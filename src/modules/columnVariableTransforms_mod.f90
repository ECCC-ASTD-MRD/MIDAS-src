
module columnVariableTransforms_mod
  ! MODULE columnVariableTransforms_mod (prefix='cvt' category='4. Data Object transformations')
  !
  !:Purpose:  To store various functions for variable transforms using inputs
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
   
    ! Arguments:
    type(struct_columnData),           intent(inout) :: columnInc
    character(len=*),                  intent(in)    :: transform
    type(struct_columnData), optional, intent(in)    :: columnRefOnIncLev_opt
    
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

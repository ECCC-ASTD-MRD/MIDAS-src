
module runtimeInfo_mod
  !
  ! MODULE runtimeInfo_mod (prefix='rti' category='8. Low-level utilities and constants')
  !
  !:Purpose: Store and print the MIDAS version number for display in the listing.
  !
  use mpi_f08

  implicit none

  save
  private

  ! Public routines
  public :: rti_abort

contains

  subroutine rti_abort(message)
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: message

    ! Locals:
    integer :: ierr

    write(6,9000) message
9000 format(//,4X,"!!!---ABORT---!!!",/,8X,"MIDAS stopped in ",A)
    flush(6)

    call mpi_abort(mpi_comm_world, 1, ierr)

  end subroutine rti_abort

end module runtimeInfo_mod

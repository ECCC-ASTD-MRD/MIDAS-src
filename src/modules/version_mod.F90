
module version_mod
  !
  ! MODULE version_mod (prefix='ver' category='8. Low-level utilities and constants')
  !
  !:Purpose: Store and print the MIDAS version number for display in the listing.
  !
  implicit none
  save
  private

  ! Public routines
  public :: ver_printNameAndVersion

  ! Acquire the 'VERSION' variable
#include <midas_build_info.h>

contains

  subroutine ver_printNameAndVersion(progName, progDescription)
    !
    !:Purpose: Print the program name, description and version to listing
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: progName
    character(len=*), intent(in) :: progDescription

    write(*,*) " --------------------------------------------"
    write(*,*) " ---  START OF MAIN PROGRAM midas-", trim(progName), " ---"
    write(*,*) " ---  ", trim(progDescription), " ---"
    write(*,*) " ---  Revision: ", trim(GIT_VERSION)
    write(*,*) " --------------------------------------------"

  end subroutine ver_printNameAndVersion

end module version_mod

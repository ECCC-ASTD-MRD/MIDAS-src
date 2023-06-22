
module version_mod
  !
  ! MODULE version_mod (prefix='ver' category='8. Low-level utilities and constants')
  !
  ! :Purpose: Store and print the MIDAS version number for display in the listing.
  !
  implicit none
  save
  private

  ! public routines
  public :: ver_printNameAndVersion

contains

  subroutine ver_printNameAndVersion(progName, progDescription)
    !
    !:Purpose: Print the program name, description and version to listing
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: progName
    character(len=*), intent(in) :: progDescription

    ! Locals:
    character(len=100), parameter :: ver_version = "GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE"

    write(*,*) " --------------------------------------------"
    write(*,*) " ---  START OF MAIN PROGRAM midas-", trim(progName), " ---"
    write(*,*) " ---  ", trim(progDescription), " ---"
    write(*,*) " ---  Revision: ", trim(ver_version)
    write(*,*) " --------------------------------------------"

  end subroutine ver_printNameAndVersion

end module version_mod

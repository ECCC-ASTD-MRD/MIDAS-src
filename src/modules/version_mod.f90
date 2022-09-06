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
    ! :Purpose: Print the program name, description and version to listing
    implicit none

    ! Arguments:
    character(len=*) :: progName
    character(len=*) :: progDescription

    ! Locals:
    character(len=100), parameter :: ver_version = "GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE"

    write(*,*) " --------------------------------------------"
    write(*,*) " ---  START OF MAIN PROGRAM midas-", trim(progName), " ---"
    write(*,*) " ---  ", trim(progDescription), " ---"
    write(*,*) " ---  Revision: ", trim(ver_version)
    write(*,*) " --------------------------------------------"

  end subroutine ver_printNameAndVersion

end module version_mod

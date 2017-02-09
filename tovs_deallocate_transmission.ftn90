!-------------------------------------- LICENCE BEGIN ------------------------------------
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

!--------------------------------------------------------------------------
!! *Purpose*: Memory deallocation for derived type transmission
!!
!! @author j. halle *cmda/aes  june 7, 2005
!!
!! revision 001  : s. heilliette arma
!!    changes for rttov 9.0 and 10.0
!!
!! arguments
!!     -ptransmission: array of derived type transmission
!!
!--------------------------------------------------------------------------
subroutine tovs_deallocate_transmission(ptransmission)
  use tovs_nl_mod,only : transmission_type
  use utilities_mod
  implicit none

  integer :: istat(2)

  type(transmission_type) :: ptransmission

  istat(:) = 0
  deallocate( ptransmission % tau_total    ,stat= istat(1))
  deallocate( ptransmission % tau_levels   ,stat= istat(2))

  if( any(istat /= 0) ) then
     write(*,*) ' tovs_deallocate_transmission: istat = ', istat(:)
     write(*,'("  tovs_deallocate_transmission: memory allocation error")')
     call utl_abort('tovs_deallocate_transmission')
  end if

  return

end subroutine tovs_deallocate_transmission

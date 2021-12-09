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
!-------------------------------------- LICENCE end --------------------------------------

module bgckOcean_mod
  ! MODULE bgckOcean_mod (prefix='sstbg' category='1. High-level functionality')
  !
  ! :Purpose: to perform ocean data background Check
  !
  use mpi_mod
  use MathPhysConstants_mod
  use utilities_mod
  use obsSpaceData_mod
  use obsErrors_mod
  use columnData_mod

  implicit none
  save
  private

  ! Public functions/subroutines
  public :: sstbg_bgCheck
  
  
  contains

  !----------------------------------------------------------------------------------------
  ! sstbg_bgCheck
  !----------------------------------------------------------------------------------------
  subroutine sstbg_bgCheck( column, obsData )
    !
    !: Purpose: to compute SST data background Check  
    !           
    
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)    :: column  ! column object on trial levels
    type(struct_obs)       , intent(inout) :: obsData ! obsSpaceData object

    ! Locals:
    character(len=*), parameter :: myName = 'ose_compute_hbht_bdiff'
    
    write(*,*) myName//': coucou'
  
  end subroutine sstbg_bgCheck
  
end module bgckOcean_mod  

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

MODULE biasCorrectionConv_mod
  ! MODULE biasCorrectionConv_mod (prefix="bcc" category='1. High-level functionality')
  !
  ! :Purpose: Performs bias correction for conventional observations
  !
  use utilities_mod
  use ramDisk_mod
  use MathPhysConstants_mod
  use obsSpaceData_mod
  use controlVector_mod
  use mpi_mod
  use mpivar_mod
  use timeCoord_mod
  use columnData_mod
  use codePrecision_mod
  use localizationFunction_mod
  use HorizontalCoord_mod
  use verticalCoord_mod
  use gridStateVector_mod
  use stateToColumn_mod
  use codtyp_mod
  use timeCoord_mod
  use clib_interfaces_mod

  implicit none
  save
  private

  public               :: bcc_readConfig
 

  
  integer, external            :: fnom, fclos
 
  namelist /nambiasconv/ 
  
CONTAINS
 
  !-----------------------------------------------------------------------
  ! bcc_readConfig
  !-----------------------------------------------------------------------
  subroutine bcc_readConfig()
    !
    ! :Purpose: Read nambiasconv namelist section
    !
    implicit none
    !Locals:
    integer  :: ierr,nulnam
  
    ! set default values for namelist variables
   
    ! read in the namelist NAMBIASCONV
    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=nambiasconv,iostat=ierr)
    if ( ierr /= 0 .and. mpi_myid == 0 )  &
         write(*,*) 'bcs_readConfig: WARNING: Error reading namelist, assume it will not be used!'
    if ( mpi_myid == 0 ) write(*,nml=nambiasconv)
    ierr = fclos(nulnam)
    

  end subroutine bcc_readConfig


end MODULE biasCorrectionSat_mod


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

module cmaFiles_mod
  ! MODULE cmaFiles_mod (prefix='cma' category='3. Observation input/output')
  !
  ! :Purpose: Read/Write "cma" format observation files, as used by the EnKF
  !
  use obsSpaceData_mod
  implicit none
  save
  private

  ! public procedures
  public :: cma_readFiles, cma_writeFiles

contains

  subroutine cma_readFiles(obsSpaceData)
    implicit none
    type(struct_obs) :: obsSpaceData

    character(len=20)   :: fileNameObsHdr, fileNameObsBdy
    integer :: unitObsHdr, unitObsBdy, unitHX
    integer :: fnom, fclos, ierr
    real(8) :: HX(1,1)

    ! open the input files
    fileNameObsHdr = './obs/cmaheader'
    unitObsHdr = 0
    ierr = fnom(unitObsHdr, fileNameObsHdr, 'FTN+SEQ+UNF+R/O', 0)

    fileNameObsBdy = './obs/cmabdy'
    unitObsBdy = 0
    ierr = fnom(unitObsBdy, fileNameObsBdy, 'FTN+SEQ+UNF+R/O', 0)

    unitHX = -1

    ! write out obsSpaceData and HX files
    call obs_read(obsSpaceData, HX, unitObsHdr, unitObsBdy, unitHX)

    ierr = fclos(unitObsHdr)
    ierr = fclos(unitObsBdy)

  end subroutine cma_readFiles


  subroutine cma_writeFiles(obsSpaceData,HXens_mpiglobal)
    implicit none

    ! arguments
    type(struct_obs)  :: obsSpaceData
    real(8)           :: HXens_mpiglobal(:,:)

    ! locals
    character(len=20) :: fileNameObsHdr, fileNameObsBdy, fileNameHX, fileNameDim
    integer           :: ierr, unitObsHdr, unitObsBdy, unitHX, unitDim, nEns
    integer           :: fnom, fclos

    nEns     = size(HXens_mpiglobal,1)

    ! open the output files
    fileNameObsHdr = 'cmaheaderout'
    unitObsHdr = 0
    ierr = fnom(unitObsHdr, fileNameObsHdr, 'FTN+SEQ+UNF+R/W', 0)

    fileNameObsBdy = 'cmabdyout'
    unitObsBdy = 0
    ierr = fnom(unitObsBdy, fileNameObsBdy, 'FTN+SEQ+UNF+R/W', 0)

    fileNameHX     = 'cmahxout'
    unitHX = 0
    ierr = fnom(unitHX, fileNameHX, 'FTN+SEQ+UNF+R/W', 0)

    fileNameDim    = 'cmadimout'
    unitDim = 0
    ierr = fnom(unitDim, fileNameDim, 'FTN+SEQ+R/W', 0)

    ! write out obsSpaceData and HX files
    call obs_write(obsSpaceData, HXens_mpiglobal, nEns, unitObsHdr, unitObsBdy, unitHX, unitDim)

    ierr = fclos(unitObsHdr)
    ierr = fclos(unitObsBdy)
    ierr = fclos(unitHX)
    ierr = fclos(unitDim)

  end subroutine cma_writeFiles

end module cmaFiles_mod

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

program midas_ensDiagnostics
  ! :Purpose: Compute diagnostics related to imbalance and spin-up in a data assimilation cycle     
  use version_mod
  use mpi_mod
  use timeCoord_mod
  use utilities_mod
  implicit none

  character(len=256) :: ensPathName,ensFileName
  character(len=4) :: charmem
  integer, allocatable :: dateStampList(:)
  integer :: ierr, nulnam
  integer, external :: fnom, fclos
  integer :: iEns,nEns ! ensemble size
  character(len=12) :: prefix ! first part of input filenames. e.g. '2019061300'
  character(len=12) :: postfix ! final part of input filenames. e.g. 'balance'

  NAMELIST /namEnsDiagnostics/nEns,prefix,postfix

  call ver_printNameAndVersion('ensDiagnostics','Program to estimate imbalance in a model integration')
  call mpi_initialize
  write(*,*) 'hello from mpi-process: ',mpi_myid
  !- Read the namelist
  nulnam = 0
  ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namEnsDiagnostics, iostat=ierr)
  if ( ierr /= 0) call utl_abort('midas-ensDiagnostics: Error reading namelist')
  if (mpi_myid == 0) then
    write(*,nml=namEnsDiagnostics)      
    write(*,*) 'ensemble size: ',nEns
  endif
  iens=1
  write(charmem,'(I4.4)') iens
  ensPathName = '../input/inputs/'
  ensFileName = trim(ensPathName)//trim(prefix)//charmem//trim(postfix)
  write(*,*) 'full input filename:',ensFileName 
  ierr = fclos(nulnam)
 
  !- 1. Initialize date/time-related info
  call tim_setup()
  allocate(dateStampList(tim_nstepobsinc))
  call tim_getstamplist(dateStampList,tim_nstepobsinc,tim_getDatestamp())
  write(*,*) 'dateStamp of first time of trajectory: ',dateStampList(1)
  write(*,*) 'dateStamp of last time of trajectory:  ',dateStampList(tim_nstepobsinc) 
  call rpn_comm_finalize(ierr)

end program midas_ensDiagnostics      


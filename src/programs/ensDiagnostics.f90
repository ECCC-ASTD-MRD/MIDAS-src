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
  use ensembleStateVector_mod
  use gridStateVector_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use analysisGrid_mod
  use timeCoord_mod
  use utilities_mod
  implicit none

  type(struct_ens), pointer :: ensembleTrl
  type(struct_vco), pointer :: vco_ens => null()
  type(struct_hco), pointer :: hco_ens => null()
  character(len=256) :: ensPathName,ensFileName
  character(len=4) :: charmem
  integer, allocatable :: dateStampList(:)
  integer :: ierr, nulnam
  integer, external :: fnom, fclos, fstopc
  integer :: iEns,nEns ! ensemble size
  character(len=12) :: prefix ! first part of input filenames. e.g. '2019061300'

  NAMELIST /namEnsDiagnostics/nEns,prefix

  call ver_printNameAndVersion('ensDiagnostics','Program to estimate imbalance in a model integration')
  call mpi_initialize
  write(*,*) 'hello from mpi-process: ',mpi_myid
  ! Avoid printing lots of stuff to listing for std file I/O
  ierr = fstopc('MSGLVL','ERRORS',0)
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
  ensFileName = trim(ensPathName)//trim(prefix)//charmem
  write(*,*) 'full input filename:',ensFileName 
  ierr = fclos(nulnam)
 
  !- 1. Initialize date/time-related info
  call tim_setup()
  allocate(dateStampList(tim_nstepobsinc))
  call tim_getstamplist(dateStampList,tim_nstepobsinc,tim_getDatestamp())
  write(*,*) 'dateStamp of first time of trajectory: ',dateStampList(1)
  write(*,*) 'dateStamp of last time of trajectory:  ',dateStampList(tim_nstepobsinc)

  !- 2. Initialize variables and grids
  call gsv_setup

  call hco_SetupFromFile(hco_ens, ensFileName, ' ', 'ENSFILEGRID')
  ! Note that the input file only contains P0 and PR fields, it does not have a vertical
  ! grid descriptor.
  call vco_setupFromFile(vco_ens, ensFileName)
  ! PLH: not sure what agd_Setup does and if this is needed
  call agd_SetupFromHCO(hco_ens)
  ! Allocate ensembles, read the Trl ensemble
  allocate(ensembleTrl)
  ! PLH: need to make sure no interpolation will be done
  call ens_allocate(ensembleTrl, nEns, tim_nstepobsinc, hco_ens, vco_ens, &
                    dateStampList)
  write(*,*) 'proper exit from ens_allocate'
  call ens_readEnsemble(ensembleTrl, ensPathName, biPeriodic=.false.)
  write(*,*) 'expecting issues with filenames prior to this'  
  call rpn_comm_finalize(ierr)

end program midas_ensDiagnostics      


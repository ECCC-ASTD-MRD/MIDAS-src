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
  integer :: i,idum, istep, j,ni,ierr, nulnam, num, numK, numStep
  integer :: myLonBeg, myLonEnd, myLatBeg, myLatEnd
  integer, external :: fnom, fclos, fstopc
  integer :: iEns,nEns ! ensemble size
  real*8, dimension(:,:), allocatable :: weight
  real*8, dimension(:), allocatable  :: MeanValMem, RainRate, Imbalance
  real*8, dimension(:,:,:,:), allocatable :: dp0dt2
  real*8  :: MeanVal, MeanValPrev
  logical :: debug
  character(len=12) :: prefix ! first part of input filenames. e.g. '2019061300'
  real(4), pointer :: onevar(:,:,:,:)

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
  iEns=1
  write(charmem,'(I4.4)') iEns
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
  write(*,*) 'ni, nj: ',hco_ens%ni, hco_ens%nj
  allocate(weight(hco_ens%ni,hco_ens%nj))
  call hco_weight(hco_ens, weight)

  call vco_setupFromFile(vco_ens, ensFileName)
  ! PLH: not sure what agd_Setup does and if this is needed
  call agd_SetupFromHCO(hco_ens)
  ! Allocate ensembles, read the Trl ensemble
  allocate(ensembleTrl)
  call ens_allocate(ensembleTrl, nEns, tim_nstepobsinc, hco_ens, vco_ens, &
                    dateStampList)
  call ens_readEnsemble(ensembleTrl, ensPathName, biPeriodic=.false.)
  num = ens_getNumMembers(ensembleTrl)
  numStep = ens_getNumStep(ensembleTrl)
  allocate(RainRate(numStep))
  allocate(Imbalance(numStep))
  numK = ens_getNumK(ensembleTrl)
  call ens_getLatLonBounds(ensembleTrl, myLonBeg, myLonEnd, myLatBeg, myLatEnd)
  if (mpi_myid == 0) then
    write(*,*) 'number of members that is stored: ',num
    write(*,*) 'number of time steps that is stored: ',numStep
    write(*,*) 'number of variables: ',numK
    write(*,*) 'precision: ',ens_getDataKind(ensembleTrl)
    idum = 1
    write(*,*) 'k for P0: ',ens_getKFromLevVarName(ensembleTrl,idum,'P0')
    write(*,*) 'k for PR: ',ens_getKFromLevVarName(ensembleTrl,idum,'PR')
    write(*,*) 'longitude bounds: ',myLonBeg, myLonEnd
    write(*,*) 'latitude bounds: ',myLatBeg, myLatEnd
  endif
  allocate(dp0dt2(nEns,numStep,myLonBeg:myLonEnd,myLatBeg:myLatEnd))
  onevar => ens_getOneLev_r4(ensembleTrl,1)
  allocate(MeanValMem(nEns))
  do j=myLatBeg,myLatEnd
    do i=myLonBeg,myLonEnd
      do istep=2,numStep-1
        do iEns=1,nEns !  estimate the second time derivative (not normalized by timestep)
          dp0dt2(iEns,istep,i,j) = onevar(iEns,istep+1,i,j) + &
                  onevar(iEns,istep-1,i,j) - 2.0D0*onevar(iEns,istep,i,j)
        enddo
      enddo  
    enddo
  enddo

  Imbalance=0.0D0
  do istep=2,numStep-1
    MeanValMem=0.0D0
    do j=myLatBeg,myLatEnd
      do i=myLonBeg,myLonEnd
        do iEns=1,nEns ! the mean square value is of interest
          MeanValMem(iens) = MeanValMem(iEns) + (dp0dt2(iEns,istep,i,j)**2)*weight(i,j)
        enddo
      enddo
    enddo
    do iEns=1,nEns ! for each member average over the horizontal grid
      call mpi_allreduce_sumreal8scalar(MeanValMem(iEns),'GRID')
    enddo
    MeanVal=0.0D0
    do iEns=1,nEns ! for each member, we move to the rms 
      MeanVal=MeanVal+MeanValMem(iEns)**0.5
    enddo
    MeanVal=MeanVal/dble(nEns)
    write(*,*) 'second derivative of P0: ',istep,MeanVal
    imbalance(istep)=MeanVal
  enddo
  if (mpi_myid == 0) then
    ierr = fnom(17,'imbalance.dat','FTN+SQN+R/W',0)
    do istep=2,numStep-1
      write(17,'(I4,x,E12.5)') istep,imbalance(istep)
    enddo
    ierr = fclos(17)
  endif    
  onevar => ens_getOneLev_r4(ensembleTrl,2)
  RainRate = 0.0D0
  do istep=1,numStep
    MeanValMem = 0.0D0
    do j=myLatBeg,myLatEnd
      do i=myLonBeg,myLonEnd
        do iEns=1,nEns
          MeanValMem(iens) = MeanValMem(iEns) + onevar(iEns,istep,i,j)*weight(i,j)
        enddo  
      enddo
    enddo
    do iEns=1,nEns ! for each member average over the horizontal grid
      call mpi_allreduce_sumreal8scalar(MeanValMem(iEns),'GRID')
    enddo
    MeanVal=0.0D0
    do iEns=1,nEns
      MeanVal=MeanVal+MeanValMem(iEns)
    enddo
    MeanVal=MeanVal/dble(nEns)
    if (istep > 1) then
      if (istep == 2) then ! at time zero accumulated precip was zero
        RainRate(istep) = MeanVal      
      else if (MeanVal < MeanValPrev) then 
        ! the PR variable is reset to zero at the central analysis time
        RainRate(istep) = MeanVal
      else      
        RainRate(istep) = MeanVal - MeanValPrev
      endif      
      write(*,*) 'step and rain rate: ',istep,RainRate(istep)
    endif
    MeanValPrev = MeanVal
  enddo
  deallocate(MeanValMem)
  if (mpi_myid == 0) then
    ierr = fnom(17,'rainrate.dat','FTN+SQN+R/W',0)
    do istep=2,numStep
      write(17,'(I4,x,E12.5)') istep,RainRate(istep)
    enddo
    ierr= fclos(17)
  endif    
  call rpn_comm_finalize(ierr)

end program midas_ensDiagnostics      


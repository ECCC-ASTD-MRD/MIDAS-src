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
  use mathPhysConstants_mod
  use ensembleStateVector_mod
  use gridStateVector_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use timeCoord_mod
  use utilities_mod
  implicit none

  type(struct_ens), pointer :: ensembleTrl
  type(struct_vco), pointer :: vco_ens => null()
  type(struct_hco), pointer :: hco_ens => null()
  character(len=256) :: ensPathName,ensFileName
  character(len=4) :: charmem
  integer, allocatable :: dateStampList(:)
  integer :: lonIndex,IndexP0,IndexPR, surfaceIndex, memberIndex, stepIndex, latIndex
  integer :: ierr, nulnam, numK, numStep
  integer :: myLonBeg, myLonEnd, myLatBeg, myLatEnd, unitNum
  integer, external :: fnom, fclos, fstopc
  real(8), allocatable :: weight(:,:)
  real(8), allocatable :: MeanValMem(:), RainRate(:), Imbalance(:)
  real(8), allocatable :: dp0dt2(:,:,:,:)
  real(8)  :: MeanVal, MeanValPrev
  real(4), pointer :: onevar(:,:,:,:)
  ! namelist variables
  integer :: nEns ! ensemble size
  character(len=256) :: pathName ! directory with input files
  character(len=12) :: prefix ! first part of input filenames. e.g. '2019061300'

  NAMELIST /namEnsDiagnostics/nEns,pathName,prefix

  call ver_printNameAndVersion('ensDiagnostics','Program to estimate imbalance in a model integration')
  call mpi_initialize
  call tmg_init(mpi_myid, 'TMG_INFO')
  call utl_tmg_start(0,'MAIN')
  ! Avoid printing lots of stuff to listing for std file I/O
  ierr = fstopc('MSGLVL','ERRORS',0)

  !- Read the namelist
  !- all values should be provided.
  nEns = mpc_missingValue_int
  pathName = 'UNDEFINED'
  prefix = 'UNDEFINED'

  nulnam = 0
  ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namEnsDiagnostics, iostat=ierr)
  if ( ierr /= 0) call utl_abort('midas-ensDiagnostics: Error reading namelist')
  if (nEns == mpc_missingValue_int) call utl_abort('midas-ensDiagnostics: set namelist value nEns')
  if (pathName == 'UNDEFINED') call utl_abort('midas-ensDiagnostics: set namelist value pathName')
  if (prefix == 'UNDEFINED') call utl_abort('midas-ensDiagnostics: set namelist value prefix')
  ierr = fclos(nulnam)
  if (mpi_myid == 0) then
    write(*,nml=namEnsDiagnostics)      
    write(*,*) 'ensemble size: ',nEns
    write(*,*) 'pathname: ',pathname
  end if

  memberIndex = 1
  write(charmem,'(I4.4)') memberIndex
  ensPathName = pathname
  ensFileName = trim(ensPathName)//'/'//trim(prefix)//charmem
  write(*,*) 'full input filename:',ensFileName 
 
  !- 1. Initialize date/time-related info
  call tim_setup
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

  ! Allocate ensembles, read the Trl ensemble
  allocate(ensembleTrl)
  call ens_allocate(ensembleTrl, nEns, tim_nstepobsinc, hco_ens, vco_ens, &
                    dateStampList)
  call ens_readEnsemble(ensembleTrl, ensPathName, biPeriodic = .false.)
  numStep = ens_getNumStep(ensembleTrl)
  allocate(RainRate(numStep))
  allocate(Imbalance(numStep))
  numK = ens_getNumK(ensembleTrl)
  call ens_getLatLonBounds(ensembleTrl, myLonBeg, myLonEnd, myLatBeg, myLatEnd)
  if (mpi_myid == 0) then
    write(*,*) 'number of time steps that is stored: ',numStep
    write(*,*) 'number of variables: ',numK
    write(*,*) 'precision: ',ens_getDataKind(ensembleTrl)
    write(*,*) 'longitude bounds: ',myLonBeg, myLonEnd
    write(*,*) 'latitude bounds: ',myLatBeg, myLatEnd
  end if
  surfaceIndex = 1 ! "surface" variables like P0 and PR have level index 1
  indexP0 = ens_getKFromLevVarName(ensembleTrl, surfaceIndex, 'P0')
  indexPR = ens_getKFromLevVarName(ensembleTrl, surfaceIndex, 'PR')

  allocate(dp0dt2(nEns,numStep,myLonBeg:myLonEnd,myLatBeg:myLatEnd))
  onevar => ens_getOneLev_r4(ensembleTrl,IndexP0)
  allocate(MeanValMem(nEns))
  do latIndex = myLatBeg,myLatEnd
    do lonIndex = myLonBeg,myLonEnd
      do stepIndex = 2,numStep-1
        do memberIndex = 1,nEns !  estimate the second time derivative (not normalized by timestep)
          dp0dt2(memberIndex,stepIndex,lonIndex,latIndex) = &
                  onevar(memberIndex,stepIndex+1,lonIndex,latIndex) + &
                  onevar(memberIndex,stepIndex-1,lonIndex,latIndex) - &
                  2.0D0*onevar(memberIndex,stepIndex,lonIndex,latIndex)
        end do
      end do  
    end do
  end do

  Imbalance = 0.0D0
  do stepIndex = 2,numStep-1
    MeanValMem = 0.0D0
    do latIndex = myLatBeg,myLatEnd
      do lonIndex = myLonBeg,myLonEnd
        do memberIndex = 1,nEns ! the mean square value is of interest
          MeanValMem(memberIndex) = MeanValMem(memberIndex) + &
                  (dp0dt2(memberIndex,stepIndex,lonIndex,latIndex)**2)* &
                  weight(lonIndex,latIndex)
        end do
      end do
    end do
    do memberIndex = 1,nEns ! for each member average over the horizontal grid
      call mpi_allreduce_sumreal8scalar(MeanValMem(memberIndex),'GRID')
    end do
    MeanVal = 0.0D0
    do memberIndex = 1,nEns ! for each member, we move to the rms 
      MeanVal = MeanVal+MeanValMem(memberIndex)**0.5
    end do
    MeanVal = MeanVal/dble(nEns)
    write(*,*) 'second derivative of P0: ',stepIndex,MeanVal
    imbalance(stepIndex) = MeanVal
  end do
  if (mpi_myid == 0) then
    unitNum = 0      
    ierr = fnom(unitNum,'imbalance.dat','FTN+SQN+R/W',0)
    do stepIndex = 2,numStep-1
      write(unitNum,'(I4,x,E12.5)') stepIndex,imbalance(stepIndex)
    end do
    ierr = fclos(unitNum)
  end if

  ! compute the rainrate per timestep (i.e. the increment in the 
  ! cumulative variable that is seen during after one timestep).  
  onevar => ens_getOneLev_r4(ensembleTrl,IndexPR)
  RainRate = 0.0D0
  do stepIndex = 1,numStep
    MeanValMem = 0.0D0
    do latIndex = myLatBeg,myLatEnd
      do lonIndex = myLonBeg,myLonEnd
        do memberIndex = 1,nEns
          MeanValMem(memberIndex) = MeanValMem(memberIndex) + & 
                  onevar(memberIndex,stepIndex,lonIndex,latIndex)* &
                  weight(lonIndex,latIndex)
        end do  
      end do
    end do
    do memberIndex = 1,nEns ! for each member average over the horizontal grid
      call mpi_allreduce_sumreal8scalar(MeanValMem(memberIndex),'GRID')
    end do
    MeanVal = 0.0D0
    do memberIndex = 1,nEns
      MeanVal = MeanVal+MeanValMem(memberIndex)
    end do
    MeanVal = MeanVal/dble(nEns)
    if (stepIndex > 1) then
      if (stepIndex == 2) then ! at time zero accumulated precip was zero
        RainRate(stepIndex) = MeanVal      
      else if (MeanVal < MeanValPrev) then 
        ! the PR variable is reset to zero at the central analysis time
        RainRate(stepIndex) = MeanVal
      else      
        RainRate(stepIndex) = MeanVal - MeanValPrev
      end if      
      write(*,*) 'step and rain rate: ',stepIndex,RainRate(stepIndex)
    end if
    MeanValPrev = MeanVal
  end do
  deallocate(MeanValMem)
  if (mpi_myid == 0) then
    unitNum=0      
    ierr = fnom(unitNum,'rainrate.dat','FTN+SQN+R/W',0)
    do stepIndex = 2,numStep
      write(unitNum,'(I4,x,E12.5)') stepIndex,RainRate(stepIndex)
    enddo
    ierr= fclos(unitNum)
  end if

  call tmg_stop(0)
  call tmg_terminate(mpi_myid, 'TMG_INFO')
  call rpn_comm_finalize(ierr)

end program midas_ensDiagnostics      


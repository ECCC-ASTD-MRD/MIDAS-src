
program main_midasprogramtemplate
  !
  !:Purpose: Include here a short description of the purpose of the program

  use rmn_fnom
  use rmn_date
  use midasMpi_mod
  use version_mod
  use message_mod
  use mathPhysConstants_mod
  use runtimeInfo_mod
  use controlVector_mod
  use gridStateVector_mod
  use bmatrix_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use timeCoord_mod
  use randomNumber_mod
  use lamAnalysisGridTransforms_mod

  implicit none

  type(struct_gsv) :: statevector

  type(struct_vco), pointer :: vco_anl => null()
  type(struct_hco), pointer :: hco_anl => null()
  type(struct_hco), pointer :: hco_core => null()

  integer :: nstamp
  integer :: ierr,status
  integer :: nkgdim
  integer :: idate,itime,ndate,nulnam
  integer :: latPerPE, latPerPEmax, myLatBeg, myLatEnd
  integer :: lonPerPE, lonPerPEmax, myLonBeg, myLonEnd

  real(8), allocatable :: controlVector(:)
  real(8), allocatable :: gdmean(:,:,:)
  real(4), allocatable :: ensemble_r4(:,:,:,:)

  character(len=10) :: cldate
  character(len=12) :: etiket

  ! namelist variables
  integer :: nens ! Ensemble size
  integer :: seed ! Seed for random number generator
  integer :: date ! Date for output standard file

  namelist /NAMENKF/nens, seed, date

  call ver_printNameAndVersion('midas-programtemplate','Example program to show how to use the MIDAS library')

  !
  !- 0. MPI, tmg initialization
  !
  call mmpi_initialize
  call tmg_init(mmpi_myid, 'TMG_MIDASPROGRAMTEMPLATE')

  call msg_memUsage('midas-programtemplate')

  !
  !- 1. Set/Read values for the namelist NAMENKF
  !

  !- 1.1 Setting default values
  nens = 10
  seed = 1
  date = 1900120100

  !- 1.2 Read the namelist
  nulnam = 0
  ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namenkf, iostat=ierr)
  if(ierr /= 0) call rti_abort('main_midasprogramtemplate: Error reading namelist')
  write(*,nml=namenkf)
  ierr = fclos(nulnam)

  ndate = date
  write(cldate,'(I10)') ndate

  call msg_memUsage('midas-programtemplate')

  !
  !- 2.  Initialization
  !

  !- 2.1 Decompose ndate(yyyymmddhh) into date(YYYYMMDD) time(HHMMSShh)
  !      calculate date-time stamp for postproc.ftn
  idate = ndate/100
  itime = (ndate-idate*100)*1000000
  nstamp = tim_yyyymmddhhToDatestamp(idate, itime)
  write(*,*)' idate= ',idate,' time= ',itime
  write(*,*)' date= ',ndate,' stamp= ',nstamp

  !- 2.2 Initialize variables of the model states
  call gsv_setup

  !
  !- Initialize the Temporal grid
  !
  call tim_setup
  call tim_setDatestamp(nstamp)

  !- 2.3 Initialize the Analysis grid
  if (mmpi_myid == 0) write(*,*)''
  if (mmpi_myid == 0) write(*,*)' preproc: Set hco parameters for analysis grid'
  call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis')

  if (hco_anl%global) then
    call lgt_SetupFromHCO(hco_anl)
  else
    !- Initialized the core (Non-Extended) analysis grid
    call hco_SetupFromFile(hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore')
    !- Setup the LAM analysis grid metrics
    call lgt_setupFromHCO(hco_anl, hco_core)
  end if

  call mmpi_setup_latbands(hco_anl%nj,                  & ! IN
                           latPerPE, latPerPEmax, myLatBeg, myLatEnd ) ! OUT
  call mmpi_setup_lonbands(hco_anl%ni,                  & ! IN
                           lonPerPE, lonPerPEmax, myLonBeg, myLonEnd ) ! OUT

  !- 2.4 Initialize the vertical coordinate from the statistics file
  if (hco_anl % global) then
    etiket = 'BGCK_STDDEV'
  else
    etiket = 'STDDEV'
  end if
  call vco_SetupFromFile(vco_anl, './bgcov', etiket)

  !- 2.5 Initialize the B_hi matrix
  call bmat_setup(hco_anl, hco_core, vco_anl)
  call msg_memUsage('midas-programtemplate')

  !
  !- 3. Memory allocations
  !

  !- 3.1 Allocate the statevector
  call gsv_allocate(statevector, 1, hco_anl, vco_anl, &
                    dateStamp_opt=nstamp, mpi_local_opt=.true.)
  nkgdim = statevector%numVarLev
  allocate(ensemble_r4(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim,nEns))

  !- 3.2 Allocate auxillary variables
  allocate(gdmean(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim), STAT=status)
  if (status /= 0) then
    call rti_abort('midas-programtemplate: PROBLEM WITH ALLOCATING OF GDMEAN')
  end if

  allocate(controlVector(cvm_nvadim))

  call msg_memUsage('midas-programtemplate')

  !
  !- 4. Compute an ensemble of random perturbations
  !

  write(*,*) '******************'
  write(*,*) 'COMPUTE the mean of the random perturbations' &
              ,' of all the members'

  call rng_setup(abs(seed+mmpi_myid))

!!!
!!! CODE REMOVED FOR THIS DEMONSTRATION

!!!

  !
  !- 5.  Memory deallocations
  !
  deallocate(gdmean,STAT=status)
  if (status /= 0) then
    call rti_abort('midas-programtemplate: PROBLEM WITH DEALLOCATE OF GDMEAN')
  end if

  deallocate(ensemble_r4)
  deallocate(controlVector)

  call msg_memUsage('midas-programtemplate')

  !
  !- 6.  MPI, tmg finalize
  !
  call tmg_terminate(mmpi_myid, 'TMG_MIDASPROGRAMTEMPLATE' )
  call mmpi_finalize

  call msg_memUsage('midas-programtemplate')

  !
  !- 7.  Ending
  !
  write(*,*) ' --------------------------------'
  write(*,*) ' MIDAS-PROGRAMTEMPLATE ENDS'
  write(*,*) ' --------------------------------'

end program main_midasprogramtemplate

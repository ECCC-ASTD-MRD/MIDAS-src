
program main_randomPert
  !
  !:Purpose: Main program for generating an ensemble of random perturbations
  !          based on the B matrix (can be homogeneous/isotropic or ensemble-based).

  use topLevelControl_mod
  use mpi_mod
  use mpivar_mod
  use mathPhysConstants_mod
  use controlVector_mod
  use gridStateVector_mod
  use bmatrix_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use analysisGrid_mod
  use timeCoord_mod
  use randomNumber_mod
  use utilities_mod
  implicit none

  type(struct_gsv) :: statevector

  type(struct_vco), pointer :: vco_anl => null()
  type(struct_hco), pointer :: hco_anl => null()
  type(struct_hco), pointer :: hco_core => null()

  integer :: fclos,fnom,newdate,nstamp
  integer :: ierr,status
  integer :: jmem,i,j,k,nkgdim
  integer :: idate,itime,ndate,nulnam,cvDim_local
  integer :: get_max_rss
  integer :: LatPerPE, myLatBeg, myLatEnd
  integer :: LonPerPE, myLonBeg, myLonEnd

  real(8), allocatable :: controlVector(:)
  real(8), allocatable :: gdmean(:,:,:)
  real(4), allocatable :: ensemble_r4(:,:,:,:)
  real(8), allocatable :: avg_pturb_var(:)
  real(8), allocatable :: avg_pturb_var_glb(:)
  real(8), allocatable :: pturb_var(:,:,:)

  character(len=10) :: cldate
  character(len=3)  :: clmember
  character(len=25) :: clfiname
  character(len=12) :: etiket
  
  ! namelist variables
  integer :: nens ! Ensemble size
  integer :: seed ! Seed for random number generator
  integer :: date ! Date for output standard file

  namelist /NAMENKF/nens, seed, date

  write(*,'(/,' //                                            &
            '3(" *****************"),/,' //                   &
            '14x,"-- START OF MAIN MAIN_RANDOMPERT --",/,' //   &
            '14x,"-- Generation of random perturbations --",/, ' //  &
            '14x,"-- VAR Revision number ",a," --",/,' //  &
            '3(" *****************"))') top_crevision

  !
  !- 0. MPI, tmg initialization
  !
  call mpi_initialize
  call tmg_init(mpi_myid, 'TMG_RANDOMPERT')

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

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
  if(ierr /= 0) call utl_abort('main_randomPert: Error reading namelist')
  write(*,nml=namenkf)
  ierr = fclos(nulnam)

  ndate = date
  write(cldate,'(I10)') ndate

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- 2.  Initialization
  !

  !- 2.1 Decompose ndate(yyyymmddhh) into date(YYYYMMDD) time(HHMMSShh)
  !      calculate date-time stamp for postproc.ftn 
  idate = ndate/100
  itime = (ndate-idate*100)*1000000
  ierr = newdate(nstamp, idate, itime, 3)
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
  if (mpi_myid == 0) write(*,*)''
  if (mpi_myid == 0) write(*,*)' preproc: Set hco parameters for analysis grid'
  call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis')

  if (hco_anl%global) then
    call agd_SetupFromHCO(hco_anl)
  else
    !- Initialized the core (Non-Extended) analysis grid
    call hco_SetupFromFile(hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore')
    !- Setup the LAM analysis grid metrics
    call agd_SetupFromHCO(hco_anl, hco_core)
  end if

  call mpivar_setup_latbands(hco_anl % nj,                & ! IN
                             latPerPE, myLatBeg, myLatEnd ) ! OUT
  call mpivar_setup_lonbands(hco_anl % ni,                & ! IN
                             lonPerPE, myLonBeg, myLonEnd ) ! OUT

  !- 2.4 Initialize the vertical coordinate from the statistics file
  if (hco_anl % global) then
    etiket = 'BGCK_STDDEV'
  else
    etiket = 'STDDEV'
  end if
  call vco_SetupFromFile(vco_anl, './bgcov', etiket)
 
  !- 2.5 Initialize the B_hi matrix
  call bmat_setup(hco_anl, vco_anl)

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- 3. Memory allocations
  !

  !- 3.1 Allocate the statevector
  call gsv_allocate(statevector, 1, hco_anl, vco_anl, &
                    dateStamp=nstamp, mpi_local=.true.)
  nkgdim = statevector%nk
  allocate(ensemble_r4(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim,nEns))

  !- 3.2 Allocate auxillary variables
  allocate(gdmean(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim), STAT=status)
  if (status /= 0) then
    call utl_abort('main_randomPert: PROBLEM WITH ALLOCATING OF GDMEAN')
  end if

  allocate(controlVector(cvm_nvadim))

  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- 4. Compute an ensemble of random perturbations
  !

  write(*,*) '******************'
  write(*,*) 'COMPUTE the mean of the random perturbations' &
              ,' of all the members'

  call rng_setup(abs(seed+mpi_myid))

!!!
!!! CODE REMOVED FOR THIS DEMONSTRATION

!!!

  !
  !- 5.  Memory deallocations
  !
  deallocate(gdmean,STAT=status)
  if (status /= 0) then
    call utl_abort('main_randomPert: PROBLEM WITH DEALLOCATE OF GDMEAN')
  end if

  deallocate(ensemble_r4)
  deallocate(controlVector)  

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- 6.  MPI, tmg finalize
  !  
  call tmg_terminate(mpi_myid, 'TMG_RANDOMPERT' )
  call rpn_comm_finalize(ierr) 

  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- 7.  Ending
  !
  write(*,*) ' --------------------------------'
  write(*,*) ' MAIN_RANDOMPERT ENDS'
  write(*,*) ' --------------------------------'

end program main_randomPert

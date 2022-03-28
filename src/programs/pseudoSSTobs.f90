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

program midas_pseudoSSTobs
  !
  ! :Purpose: Main program to compute SST bias estimate
  !
  use version_mod
  use ramDisk_mod
  use utilities_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use analysisGrid_mod
  use mpi_mod
  use oceanObservations_mod
  use gridStateVector_mod
  use obsSpaceData_mod
   
  implicit none

  integer, external :: exdb, exfin, fnom, fclos, get_max_rss
  integer :: ierr, istamp, nulnam

  type(struct_hco), pointer   :: hco_anl => null()
  type(struct_vco), pointer   :: vco_anl => null()
  character(len=48),parameter :: obsMpiStrategy = 'LIKESPLITFILES', &
                                 varMode        = 'analysis'
  real(4)                     :: iceFractionThreshold    ! consider no ice condition below this threshold
  real(4)                     :: outputSST               ! output SST value for pseudo observations
  real(4)                     :: outputFreshWaterST      ! output fresh water surface temperature for pseudo obs.  
  integer                     :: seaiceThinning          ! generate pseudo obs in every 'seaiceThinning' points 
  character(len=20)           :: outputFileName
  character(len=20)           :: etiket
  real(4)                     :: seaWaterThreshold       ! to distinguish inland water from sea water  

  namelist /pseudoSSTobs/ iceFractionThreshold, outputSST, outputFreshWaterST, seaiceThinning, outputFileName, etiket, seaWaterThreshold
  
  istamp = exdb('pseudoSSTobs','DEBUT','NON')

  call ver_printNameAndVersion('pseudoSSTobs','Generation of pseudo SST observations')

  ! MPI initialization
  call mpi_initialize

  call tmg_init(mpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'MAIN')
 
  ! 1. Top level setup

  call ram_setup()
 
  ! Do initial set up
  call utl_tmg_start(2,'--SETUP')
  call pseudoSSTobs_setup()
  call tmg_stop(2)

  call oobs_pseudoSST(hco_anl, vco_anl, iceFractionThreshold, outputSST, outputFreshWaterST, &
                      seaiceThinning, outputFileName, etiket, seaWaterThreshold)

  ! 3. Job termination

  istamp = exfin('pseudoSSTobs','FIN','NON')

  call tmg_stop(0)

  call tmg_terminate(mpi_myid, 'TMG_INFO')

  call rpn_comm_finalize(ierr) 

  contains

  subroutine pseudoSSTobs_setup()
    !
    ! :Purpose:  Control of the preprocessing of pseudo SST obs
    !

    implicit none
    
    !Locals:	
    type(struct_hco), pointer   :: hco_core => null()
    character(len=*), parameter :: myName = 'pseudoSSTobs_setup'
    character(len=*), parameter :: gridFile = './analysisgrid'
    
    write(*,*) ''
    write(*,*) '-------------------------------------------------'
    write(*,*) '-- Starting subroutine '//myName//' --'
    write(*,*) '-------------------------------------------------'

    ! Setting default namelist variable values
    iceFractionThreshold   = 0.05 
    outputSST              = 271.4
    outputFreshWaterST     = 271.4
    outputFileName         = ''
    etiket                 = ''
    seaiceThinning         = 5
    seaWaterThreshold      = 0.0
    
    ! Read the namelist
    nulnam = 0
    ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    read(nulnam, nml = pseudoSSTobs, iostat = ierr)
    if (ierr /= 0) call utl_abort(myName//': Error reading namelist')
    if (mpi_myid == 0) write(*, nml = pseudoSSTobs)
    ierr = fclos(nulnam)

    write(*,*)''
    write(*,*) myName//': output SST globally: ', outputSST
    write(*,*) myName//': output fresh water surface temperature  globally: ', outputFreshWaterST
    write(*,*) myName//': sea-ice fraction threshold: ', iceFractionThreshold
    write(*,*) myName//': sea water fraction threshold: ', seaWaterThreshold
    write(*,*) myName//': pseudo SST obs will be generated in every ', seaiceThinning, ' points of the sea-ice field'    
    write(*,*) myName//': output file name: ', outputFileName
    write(*,*) myName//': etiket for output SQLite files: ', etiket   
    !
    !- Initialize the Analysis grid
    !
    if(mpi_myid == 0) write(*,*)''
    if(mpi_myid == 0) write(*,*) myName//': Set hco parameters for analysis grid'
    call hco_SetupFromFile(hco_anl, gridFile, 'ANALYSIS') ! IN

    if (hco_anl % global) then
      call agd_SetupFromHCO(hco_anl) ! IN
    else
      !- Initialize the core (Non-Extended) analysis grid
      if(mpi_myid == 0) write(*,*) myName//': Set hco parameters for core grid'
      call hco_SetupFromFile(hco_core, gridFile, 'COREGRID', 'AnalysisCore') ! IN
      !- Setup the LAM analysis grid metrics
      call agd_SetupFromHCO(hco_anl, hco_core) ! IN
    end if

    !     
    !- Initialisation of the analysis grid vertical coordinate from analysisgrid file
    !
    call vco_SetupFromFile(vco_anl, & ! OUT
                           gridFile) ! IN

    !- Initialize variables of the model states
    !
    call gsv_setup

    call obs_class_initialize('VAR')

    if(mpi_myid == 0) write(*,*) myName//': done.'
    
  end subroutine pseudoSSTobs_setup

end program midas_pseudoSSTobs

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

program midas_sstBias
  !
  ! :Purpose: Main program to compute SST bias correction
  !
  use version_mod
  use ramDisk_mod
  use utilities_mod
  use mpi_mod
  use mpivar_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use obsSpaceData_mod
  use columnData_mod  
  use gridStateVector_mod
  use obsSpaceDiag_mod
  use obsFiles_mod
  use obsFilter_mod  
  use innovation_mod
  use obsErrors_mod
  use statetocolumn_mod
  use biasCorrectionSat_mod
  use increment_mod
  use stateToColumn_mod
  use backgroundCheck_mod
  use analysisGrid_mod
  use SSTbiasEstimation_mod
  
  implicit none

  integer, external :: exdb, exfin, fnom, fclos, get_max_rss
  integer :: ierr, istamp, nulnam

  type(struct_obs),        target :: obsSpaceData
  type(struct_columnData), target :: trlColumnOnAnlLev
  type(struct_hco), pointer       :: hco_anl => null()
  type(struct_vco), pointer       :: vco_anl => null()
  real(8)                         :: horizontalSearchRadius ! horizontal search radius, in km, around input grid points
                                                            ! where to search satellite observations
  character(len=48),parameter     :: obsMpiStrategy = 'LIKESPLITFILES', &
                                     varMode        = 'analysis'

  istamp = exdb('SSTBIASESTIMATION','DEBUT','NON')

  call ver_printNameAndVersion('SSTbiasEstimation','SST Bias Estimation')

  ! MPI initialization
  call mpi_initialize

  call tmg_init(mpi_myid, 'TMG_SSTbiasEstimation' )

  call tmg_start(1,'MAIN')
 
  ! 1. Top level setup

  call ram_setup()
 
  ! Do initial set up
  call tmg_start(2,'SETUP')
  call SSTbiasEstimation_setup('VAR') ! obsColumnMode
  call tmg_stop(2)
  
  call sstb_computeGriddedObservations( obsSpaceData, hco_anl, horizontalSearchRadius )

  ! Deallocate copied obsSpaceData
  call obs_finalize(obsSpaceData)
  
  ! 3. Job termination

  istamp = exfin('SSTBIASESTIMATION','FIN','NON')

  call tmg_stop(1)

  call tmg_terminate(mpi_myid, 'TMG_SSTbiasEstimation' )

  call rpn_comm_finalize(ierr) 

  contains

  subroutine SSTbiasEstimation_setup(obsColumnMode)
    !
    ! :Purpose:  Control of the preprocessing of bias estimation
    !

    implicit none
    !Arguments:
    character(len=*), intent(in) :: obsColumnMode
    !Locals:	
    integer :: datestamp
    type(struct_hco), pointer :: hco_core => null()
    character(len=*), parameter :: myName = 'SSTbiasEstimation_setup'
    namelist /namSSTbiasEstimate/ horizontalSearchRadius
    
    write(*,*) ''
    write(*,*) '-------------------------------------------------'
    write(*,*) '-- Starting subroutine '//myName//' --'
    write(*,*) '-------------------------------------------------'

    !
    !- Initialize the Temporal grid
    !
    call tim_setup

    !     
    !- Initialize observation file names and set datestamp
    !
    call obsf_setup( dateStamp, varMode )

    dateStamp = tim_getDatestampFromFile("./trlm_01")

    if ( dateStamp > 0 ) then
      call tim_setDatestamp(datestamp)     ! IN
    else
      call utl_abort( myName//': Problem getting dateStamp from observation first trial field')
    end if

    ! Setting default namelist variable values
    horizontalSearchRadius = 100.  ! horizontal search radius, in km, around input grid points
                                   ! where to search satellite observations
    ! Read the namelist
    nulnam = 0
    ierr = fnom( nulnam, './flnml', 'FTN+SEQ+R/O', 0 )
    read( nulnam, nml = namSSTbiasEstimate, iostat = ierr )
    if ( ierr /= 0) call utl_abort( myName//': Error reading namelist')
    if ( mpi_myid == 0 ) write(*, nml = namSSTbiasEstimate )
    ierr = fclos( nulnam )

    !
    !- Initialize constants
    !
    if(mpi_myid == 0) call mpc_printConstants(6)

    !
    !- Initialize variables of the model states
    !
    call gsv_setup
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Initialize the Analysis grid
    !
    if(mpi_myid == 0) write(*,*)''
    if(mpi_myid == 0) write(*,*) myName//': Set hco parameters for analysis grid'
    call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

    if ( hco_anl % global ) then
      call agd_SetupFromHCO( hco_anl ) ! IN
    else
      !- Initialize the core (Non-Extended) analysis grid
      if(mpi_myid == 0) write(*,*) myName//': Set hco parameters for core grid'
      call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore' ) ! IN
      !- Setup the LAM analysis grid metrics
      call agd_SetupFromHCO( hco_anl, hco_core ) ! IN
    end if

    !     
    !- Initialisation of the analysis grid vertical coordinate from analysisgrid file
    !
    call vco_SetupFromFile( vco_anl,        & ! OUT
                            './analysisgrid') ! IN

    call col_setVco(trlColumnOnAnlLev,vco_anl)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Setup and read observations
    !
    call inn_setupObs(obsSpaceData, obsColumnMode, obsMpiStrategy, varMode) ! IN
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Basic setup of columnData module
    !
    call col_setup
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if(mpi_myid == 0) write(*,*) myName//': done.'
    
  end subroutine SSTbiasEstimation_setup

end program midas_sstBias

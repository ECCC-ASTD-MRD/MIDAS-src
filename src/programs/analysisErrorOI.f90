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

program midas_analysisErrorOI
  !
  ! :Purpose: Calculate analysis-error standard deviation given
  !           new assimilated observations. It only works for sea ice variables and
  !           uses a simple OI approach.
  !
  use version_mod
  use codePrecision_mod
  use ramDisk_mod
  use utilities_mod
  use mpi_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use obsSpaceData_mod
  use columnData_mod
  use gridStateVector_mod
  use obsSpaceDiag_mod
  use controlVector_mod
  use obsFiles_mod
  use innovation_mod
  use analysisGrid_mod
  use obsErrors_mod
  use analysisError_mod

  implicit none

  integer :: istamp, exdb, exfin
  integer :: ierr, dateStamp
  integer :: get_max_rss
  character(len=48) :: obsMpiStrategy, varMode
  character(len=20) :: trlmFileName
  character(len=15), parameter :: myName = 'analysisErrorOI'

  type(struct_obs)       , target :: obsSpaceData
  type(struct_columnData), target :: trlColumnOnAnlLev
  type(struct_hco)      , pointer :: hco_anl => null()
  type(struct_vco)      , pointer :: vco_anl => null()
  type(struct_hco)      , pointer :: hco_core => null()

  istamp = exdb('ANALYSISERROROI','DEBUT','NON')

  call ver_printNameAndVersion(myName,'Program to calculate the analysis-error standard deviation for sea ice using OI.')

  ! MPI initialization
  call mpi_initialize

  if( mpi_nprocs > 1 ) then
    write(*,*) 'mpi_nprocs = ',mpi_nprocs
    call utl_abort(myName//': this version of the code should only be used with one mpi task.')
  end if

  call tmg_init(mpi_myid, 'TMG_ANALYSISERROROI' )

  call tmg_start(1,'MAIN')

  varMode='analysis'

  ! Setup the ram disk
  call ram_setup

  obsMpiStrategy = 'LIKESPLITFILES'

  !
  !- Initialize the Temporal grid
  !
  call tim_setup

  if (tim_nstepobs > 1 .or. tim_nstepobsinc > 1) then
    call utl_abort(myName//': The program assumes only one time step.')
  end if

  !
  !- Initialize observation file names and set datestamp
  !
  call obsf_setup( dateStamp, varMode )
  if ( dateStamp > 0 ) then
    call tim_setDatestamp(datestamp)     ! IN
  else
    call utl_abort(myName//': Problem getting dateStamp from observation file')
  end if

  !
  !- Initialize constants
  !
  if (mpi_myid == 0) call mpc_printConstants(6)

  trlmFileName = './trlm_01'

  !
  !- Initialize the Analysis grid
  !
  if (mpi_myid == 0) write(*,*)
  if (mpi_myid == 0) write(*,*) 'var_setup: Set hco parameters for analysis grid'
  call hco_SetupFromFile(hco_anl, trlmFileName, ' ') ! IN

  if ( hco_anl % global ) then
    call agd_SetupFromHCO( hco_anl ) ! IN
  else
    !- Initialize the core (Non-Extended) analysis grid
    if (mpi_myid == 0) write(*,*)'var_setup: Set hco parameters for core grid'
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
  call inn_setupObs(obsSpaceData, hco_anl, 'VAR', obsMpiStrategy, varMode) ! IN
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- Basic setup of columnData module
  !
  call col_setup
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- Memory allocation for background column data
  !
  call col_allocate(trlColumnOnAnlLev,obs_numheader(obsSpaceData),mpiLocal_opt=.true.)

  !
  !- Initialize the observation error covariances
  !
  call oer_setObsErrors(obsSpaceData, varMode) ! IN
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  ! Initialize list of analyzed variables.
  !
  call gsv_setup
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Compute the analysis-error
  call aer_analysisError(obsSpaceData, hco_anl, vco_anl, trlmFileName)

  ! Now write out the observation data files
  if ( .not. obsf_filesSplit() ) then 
    write(*,*) 'We read/write global observation files'
    call obs_expandToMpiGlobal(obsSpaceData)
    if (mpi_myid == 0) call obsf_writeFiles(obsSpaceData)
  else
    ! redistribute obs data to how it was just after reading the files
    call obs_MpiRedistribute(obsSpaceData,OBS_IPF)
    call obsf_writeFiles(obsSpaceData)
  end if

  ! Deallocate copied obsSpaceData
  call obs_finalize(obsSpaceData)

  !
  ! 3. Job termination
  !
  istamp = exfin('ANALYSISERROROI','FIN','NON')

  call tmg_stop(1)

  call tmg_terminate(mpi_myid, 'TMG_ANALYSISERROROI' )

  call rpn_comm_finalize(ierr) 

end program midas_analysisErrorOI

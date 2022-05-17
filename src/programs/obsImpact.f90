!------------------------------------- LICENCE BEGIN -------------------------------------
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

program midas_obsimpact
  !
  ! :Purpose: Main program for Observation Impact computation (FSOI)
  !
  use version_mod
  use codePrecision_mod
  use ramDisk_mod
  use utilities_mod
  use mpi_mod
  use mpivar_mod
  use MathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use columnData_mod
  use obsSpaceData_mod
  use controlVector_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use bmatrix_mod
  use bmatrixensemble_mod
  use stateToColumn_mod
  use obsOperators_mod
  use costFunction_mod
  use quasinewton_mod
  use innovation_mod
  use obsFiles_mod
  use obsFilter_mod
  use obsErrors_mod
  use gridVariableTransforms_mod
  use rttov_const, only :inst_name, platform_name
  use tovs_nl_mod
  use rMatrix_mod
  use fsoi_mod

  implicit none

  integer :: istamp,exdb,exfin,ierr, dateStamp

  type(struct_columnData),target :: columnTrlOnTrlLev
  type(struct_gsv)               :: stateVectorTrialHighRes

  character(len=48) :: obsMpiStrategy
  character(len=3)  :: obsColumnMode

  type(struct_vco),        pointer :: vco_anl => null()
  type(struct_hco),        pointer :: hco_anl => null()
  type(struct_hco),        pointer :: hco_trl => null()
  type(struct_vco),        pointer :: vco_trl => null()
  type(struct_hco),        pointer :: hco_core => null()
  integer,external    :: get_max_rss
  logical             :: allocHeightSfc

  istamp = exdb('OBSIMPACT','DEBUT','NON')

  call ver_printNameAndVersion('obsImpact','Calculation of observation impact')

  ! MPI initilization
  call mpi_initialize

  call tmg_init(mpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')

  if (mpi_myid == 0) then
    call utl_writeStatus('VAR3D_BEG')
  end if

  call ram_setup

  !
  !- 1. Settings 
  !
  obsColumnMode  = 'VAR'
  obsMpiStrategy = 'LIKESPLITFILES'

  !- Do initial set up
  call fso_setup

  !
  !- Initialize the Temporal grid
  !
  call tim_setup
  !     
  !- Initialize burp file names and set datestamp
  !
  call obsf_setup( dateStamp, 'FSO' )
  if ( dateStamp > 0 ) then
    call tim_setDatestamp(datestamp)     ! IN
  else
    call utl_abort('obsImpact: Problem getting dateStamp from observation file')
  end if
  !
  !- Initialize constants
  !
  if ( mpi_myid == 0 ) then
    call mpc_printConstants(6)
    call pre_printPrecisions
  end if

  !
  !- Initialize variables of the model states
  !
  call gsv_setup
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
  !
  !- Initialize the Analysis grid
  !
  if (mpi_myid.eq.0) write(*,*)''
  if (mpi_myid.eq.0) write(*,*)'obsImpact: Set hco parameters for analysis grid'
  call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

  if ( hco_anl % global ) then
    hco_core => hco_anl
  else
    !- Iniatilized the core (Non-Exteded) analysis grid
    call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore' ) ! IN
  end if

  !     
  !- Initialisation of the analysis grid vertical coordinate from analysisgrid file !
  call vco_SetupFromFile( vco_anl,        & ! OUT
                          './analysisgrid') ! IN

  call col_setVco(columnTrlOnAnlIncLev,vco_anl)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- Setup and read observations
  !
  call inn_setupObs(obsSpaceData, hco_anl, obsColumnMode, obsMpiStrategy, 'FSO') ! IN
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- Basic setup of columnData module
  !
  call col_setup
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- Memory allocation for background column data
  !
  call col_allocate(columnTrlOnAnlIncLev,obs_numheader(obsSpaceData),mpiLocal_opt=.true.)

  !
  !- Initialize the observation error covariances
  !
  call oer_setObsErrors(obsSpaceData, 'FSO') ! IN
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Reading trials
  call inn_getHcoVcoFromTrlmFile( hco_trl, vco_trl )
  allocHeightSfc = ( vco_trl%Vcode /= 0 )

  call gsv_allocate( stateVectorTrialHighRes, tim_nstepobs, hco_trl, vco_trl,  &
                     dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                     mpi_distribution_opt='Tiles', dataKind_opt=4,  &
                     allocHeightSfc_opt=allocHeightSfc, hInterpolateDegree_opt='LINEAR', &
                     beSilent_opt=.false. )
  call gsv_zero( stateVectorTrialHighRes )
  call gio_readTrials( stateVectorTrialHighRes )
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Horizontally interpolate trials to trial columns
  call inn_setupColumnsOnTrlLev( columnTrlOnTrlLev, obsSpaceData, hco_core, &
                                   stateVectorTrialHighRes )
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- Initialize the background-error covariance, also sets up control vector module (cvm)
  !
  call bmat_setup(hco_anl,hco_core,vco_anl)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  ! - Initialize the gridded variable transform module
  !
  call gvt_setup(hco_anl,hco_core,vco_anl)
  call gvt_setupRefFromTrialFiles('HU')
  call gvt_setupRefFromTrialFiles('height')

  !
  !- 2. Do the actual job
  !

  ! Interpolate trial columns to analysis levels and setup for linearized H
  call inn_setupColumnsOnAnlIncLev(columnTrlOnTrlLev,columnTrlOnAnlIncLev)

  ! Compute observation innovations and prepare obsSpaceData for minimization
  call inn_computeInnovation(columnTrlOnTrlLev,obsSpaceData)

  ! Perform forecast sensitivity to observation calculation using ensemble approach 
  call fso_ensemble(columnTrlOnAnlIncLev,obsSpaceData)

  ! Deallocate memory related to B matrices
  call bmat_finalize()

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
  !- 3. Job termination
  !
  istamp = exfin('OBSIMPACT','FIN','NON')

  if (mpi_myid == 0) then
    call utl_writeStatus('VAR3D_END')
  endif

  call utl_tmg_stop(0)

  call tmg_terminate(mpi_myid, 'TMG_INFO')

  call rpn_comm_finalize(ierr)

end program midas_obsimpact 

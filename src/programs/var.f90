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

program midas_var
  !
  ! :Purpose: Main program for variational minimization or background check 
  !           (depending on the mode selected in the namelist).
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
  use minimization_mod
  use innovation_mod
  use bmatrix_mod
  use obsErrors_mod
  use gridVariableTransforms_mod
  use increment_mod
  use biasCorrectionSat_mod

  implicit none

  integer :: istamp, exdb, exfin
  integer :: ierr, dateStamp
  integer :: get_max_rss
  character(len=9)  :: clmsg
  character(len=48) :: obsMpiStrategy, varMode
  real(8), allocatable :: controlVectorIncr(:)

  type(struct_obs)       , target :: obsSpaceData
  type(struct_columnData), target :: columnTrlOnAnlIncLev
  type(struct_columnData), target :: columnTrlOnTrlLev
  type(struct_gsv)                :: stateVectorIncr
  type(struct_gsv)                :: stateVectorUpdateHighRes
  type(struct_gsv)                :: stateVectorTrial
  type(struct_gsv)                :: statevectorPsfc
  type(struct_gsv)                :: stateVectorLowResTime
  type(struct_gsv)                :: stateVectorLowResTimeSpace
  type(struct_gsv)                :: stateVectorAnalHighRes
  type(struct_gsv)       , target :: stateVectorRefHU
  type(struct_hco)      , pointer :: hco_anl => null()
  type(struct_vco)      , pointer :: vco_anl => null()
  type(struct_hco)      , pointer :: hco_core => null()

  character(len=4), pointer :: varNames(:)
  logical :: allocHeightSfc

  istamp = exdb('VAR','DEBUT','NON')

  call ver_printNameAndVersion('var','Variational Assimilation')

  ! MPI initialization
  call mpi_initialize

  call tmg_init(mpi_myid, 'TMG_VAR' )

  call tmg_start(1,'MAIN')

  if (mpi_myid == 0) then
    clmsg = 'VAR3D_BEG'
    call utl_writeStatus(clmsg)
  end if 

  write(*,*)
  write(*,*) 'Real Kind used for computing the increment =', pre_incrReal
  write(*,*)

  varMode='analysis'

  ! Setup the ram disk
  call ram_setup

  ! Do initial set up
  call tmg_start(2,'PREMIN')

  obsMpiStrategy = 'LIKESPLITFILES'

  !
  !- Initialize the Temporal grid
  !
  call tim_setup

  !     
  !- Initialize observation file names and set datestamp
  !
  call obsf_setup( dateStamp, varMode )
  if ( dateStamp > 0 ) then
    call tim_setDatestamp(datestamp)     ! IN
  else
    call utl_abort('var_setup: Problem getting dateStamp from observation file')
  end if

  !
  !- Initialize constants
  !
  if (mpi_myid == 0) call mpc_printConstants(6)

  !
  !- Initialize the Analysis grid
  !
  if (mpi_myid == 0) write(*,*)''
  if (mpi_myid == 0) write(*,*)'var_setup: Set hco parameters for analysis grid'
  call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

  if ( hco_anl % global ) then
    hco_core => hco_anl
  else
    !- Initialize the core (Non-Extended) analysis grid
    if (mpi_myid == 0) write(*,*)'var_setup: Set hco parameters for core grid'
    call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore' ) ! IN
  end if

  !     
  !- Initialisation of the analysis grid vertical coordinate from analysisgrid file
  !
  call vco_SetupFromFile( vco_anl,        & ! OUT
                          './analysisgrid') ! IN

  call col_setVco(columnTrlOnAnlIncLev,vco_anl)
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
  call col_allocate(columnTrlOnAnlIncLev,obs_numheader(obsSpaceData),mpiLocal_opt=.true.)

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

  !
  !- Initialize the background-error covariance, also sets up control vector module (cvm)
  !
  call bmat_setup(hco_anl,hco_core,vco_anl)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  ! - Initialize the gridded variable transform module
  !
  call gvt_setup(hco_anl,hco_core,vco_anl)

  !
  !- Set up the minimization module, now that the required parameters are known
  !  NOTE: some global variables remain in minimization_mod that must be initialized before
  !        inn_setupColumnsOnTrialLev
  !
  call min_setup( cvm_nvadim, hco_anl ) ! IN
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Reading 15-min trials
  call gsv_readTrialsHighRes( stateVectorUpdateHighRes, &
                              stateVectorTrialOut_opt=stateVectorTrial )
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  if ( stateVectorUpdateHighRes%vco%Vcode == 0 ) then
    allocHeightSfc = .false.
  else
    allocHeightSfc = .true.
  end if

  ! Horizontally interpolate high-resolution stateVectorUpdate to trial columns
  call inn_setupColumnsOnTrialLev( columnTrlOnTrlLev, obsSpaceData, hco_core, &
                                   stateVectorUpdateHighRes )
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Interpolate trial columns to analysis levels and setup for linearized H
  call inn_setupColumnsOnAnlLev( columnTrlOnTrlLev, columnTrlOnAnlIncLev )
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Compute observation innovations and prepare obsSpaceData for minimization
  call inn_computeInnovation( columnTrlOnTrlLev, obsSpaceData )
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Initialize stateVectorRefHU for doing variable transformation of the increments.
  if ( gsv_varExist(stateVectorUpdateHighRes,'HU') ) then
    call gsv_allocate(stateVectorRefHU, tim_nstepobsinc, hco_anl, vco_anl,   &
                      dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                      allocHeightSfc_opt=.true., hInterpolateDegree_opt='LINEAR', &
                      varNames_opt=(/'HU','P0'/) )

    ! First degrade the time steps
    call gsv_allocate( stateVectorLowResTime, tim_nstepobsinc, &
                       stateVectorUpdateHighRes%hco, stateVectorUpdateHighRes%vco,  &
                       dataKind_opt=pre_incrReal, &
                       dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                       allocHeightSfc_opt=allocHeightSfc, hInterpolateDegree_opt='LINEAR', &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_copy( stateVectorUpdateHighRes, stateVectorLowResTime,  &
                   allowTimeMismatch_opt=.true., allowVarMismatch_opt=.true. )

    ! Second interpolate to the low-resolution spatial grid.
    nullify(varNames)
    call gsv_varNamesList(varNames, stateVectorLowResTime)
    call gsv_allocate(stateVectorLowResTimeSpace, tim_nstepobsinc, hco_anl, vco_anl,   &
                      dataKind_opt=pre_incrReal, &
                      dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                      allocHeightSfc_opt=.true., hInterpolateDegree_opt='LINEAR', &
                      varNames_opt=varNames)
    call gsv_interpolate(stateVectorLowResTime, stateVectorLowResTimeSpace)        

    ! Now copy only P0 and HU.
    call gsv_copy( stateVectorLowResTimeSpace, stateVectorRefHU, &
                   allowTimeMismatch_opt=.false., allowVarMismatch_opt=.true. )
    call gsv_deallocate(stateVectorLowResTimeSpace)
    call gsv_deallocate(stateVectorLowResTime)

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
  end if
  call tmg_stop(2)

  allocate(controlVectorIncr(cvm_nvadim),stat=ierr)
  if (ierr /= 0) then
    write(*,*) 'var: Problem allocating memory for ''controlVectorIncr''',ierr
    call utl_abort('aborting in VAR')
  end if

  ! Do minimization of cost function
  call min_minimize( columnTrlOnAnlIncLev, obsSpaceData, controlVectorIncr, &
                     stateVectorRef_opt=stateVectorRefHU )
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Compute satellite bias correction increment and write to file
  call bcs_writebias(controlVectorIncr)

  call gsv_allocate(stateVectorIncr, tim_nstepobsinc, hco_anl, vco_anl, &
       datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true., &
       dataKind_opt=pre_incrReal, allocHeight_opt=.false., allocPressure_opt=.false.)

  ! get final increment with mask if it exists
  call inc_getIncrement(controlVectorIncr, stateVectorIncr, cvm_nvadim, &
                        statevectorRef_opt=stateVectorRefHU)
  call gsv_readMaskFromFile(stateVectorIncr,'./analysisgrid')
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Compute high-resolution analysis on trial grid
  call inc_computeHighResAnalysis(stateVectorIncr, stateVectorUpdateHighRes, &
                                  statevectorPsfc, stateVectorAnalHighRes)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! output the analysis increment
  call tmg_start(6,'WRITEINCR')
  call inc_writeIncrement(stateVectorIncr) ! IN
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
  call tmg_stop(6)

  ! Interpolate analysis stateVector from 7 to 25 timeStep
  call gsv_tInterpolate( stateVectorTrial, stateVectorUpdateHighRes)

  ! Conduct obs-space post-processing diagnostic tasks (some diagnostic
  ! computations controlled by NAMOSD namelist in flnml)
  call osd_ObsSpaceDiag( obsSpaceData, columnTrlOnAnlIncLev, hco_anl )

  ! Deallocate memory related to B matrices
  call bmat_finalize()

  ! Deallocate the high-resolution statevector 
  call gsv_deallocate(stateVectorUpdateHighRes)

  ! compute and write the analysis (as well as the increment on the trial grid)
  call tmg_start(18,'ADDINCREMENT')
  call inc_writeIncrementHighRes(stateVectorIncr, stateVectorTrial, &
                                 statevectorPsfc, stateVectorAnalHighRes)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
  call tmg_stop(18)

  if (mpi_myid == 0) then
    clmsg = 'REBM_DONE'
    call utl_writeStatus(clmsg)
  end if

  call gsv_deallocate(stateVectorIncr)
  if ( stateVectorRefHU%allocated ) call gsv_deallocate(stateVectorRefHU)

  ! write the Hessian
  call min_writeHessian(controlVectorIncr)
  deallocate(controlVectorIncr)

  ! Deallocate memory related to variational bias correction
  call bcs_finalize()

  ! Now write out the observation data files
  if (min_niter > 0) then
    if ( .not. obsf_filesSplit() ) then 
      write(*,*) 'We read/write global observation files'
      call obs_expandToMpiGlobal(obsSpaceData)
      if (mpi_myid == 0) call obsf_writeFiles(obsSpaceData)
    else
      ! redistribute obs data to how it was just after reading the files
      call obs_MpiRedistribute(obsSpaceData,OBS_IPF)
      call obsf_writeFiles(obsSpaceData)
    end if
  end if

  ! Deallocate copied obsSpaceData
  call obs_finalize(obsSpaceData)

  !
  ! 3. Job termination
  !
  istamp = exfin('VAR','FIN','NON')

  if (mpi_myid == 0) then
    clmsg = 'VAR3D_END'
    call utl_writeStatus(clmsg)
  end if

  call tmg_stop(1)

  call tmg_terminate(mpi_myid, 'TMG_VAR' )

  call rpn_comm_finalize(ierr) 

end program midas_var

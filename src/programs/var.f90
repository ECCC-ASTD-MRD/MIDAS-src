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
  use humidityLimits_mod
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
  use rmatrix_mod
  use obsErrors_mod
  use gridVariableTransforms_mod
  use increment_mod
  use biasCorrectionSat_mod
  use varqc_mod
  use tovs_nl_mod
  use stateToColumn_mod

  implicit none

  integer :: istamp, exdb, exfin
  integer :: ierr, dateStamp, nulnam
  integer :: get_max_rss, fclos, fnom
  character(len=9)  :: clmsg
  character(len=48) :: obsMpiStrategy, varMode
  real(8), allocatable :: controlVectorIncr(:)
  real(8), allocatable :: controlVectorIncrSum(:)

  type(struct_obs)       , target :: obsSpaceData
  type(struct_columnData), target :: columnTrlOnAnlIncLev
  type(struct_columnData), target :: columnTrlOnTrlLev
  type(struct_gsv)                :: stateVectorIncr
  type(struct_gsv)                :: stateVectorIncrSum
  type(struct_gsv)                :: stateVectorUpdateHighRes
  type(struct_gsv)                :: stateVectorTrial
  type(struct_gsv)                :: stateVectorPsfcHighRes
  type(struct_gsv)                :: stateVectorPsfc
  type(struct_gsv)                :: stateVectorAnal
  type(struct_hco)      , pointer :: hco_anl => null()
  type(struct_vco)      , pointer :: vco_anl => null()
  type(struct_hco)      , pointer :: hco_trl => null()
  type(struct_vco)      , pointer :: vco_trl => null()
  type(struct_hco)      , pointer :: hco_core => null()

  integer :: outerLoopIndex, ip3ForWriteToFile
  integer :: numIterWithoutVarqc, numInnerLoopIterDone
  integer :: numIterMaxInnerLoopUsed

  logical :: allocHeightSfc, applyLimitOnHU
  logical :: deallocHessian, isMinimizationFinalCall
  logical :: varqcActive, applyVarqcOnNlJo
  logical :: filterObsAndInitOer
  logical :: deallocInterpInfoNL

  integer, parameter :: maxNumberOfOuterLoopIterations = 15

  ! namelist variables
  integer :: numOuterLoopIterations, numIterMaxInnerLoop(maxNumberOfOuterLoopIterations)
  logical :: limitHuInOuterLoop, computeFinalNlJo
  NAMELIST /NAMVAR/ numOuterLoopIterations, numIterMaxInnerLoop, limitHuInOuterLoop
  NAMELIST /NAMVAR/ computeFinalNlJo

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

  varMode='analysis'

  ! Setup the ram disk
  call ram_setup

  ! Do initial set up
  call tmg_start(2,'PREMIN')

  ! Set/Read values for the namelist NAMVAR
  ! Setting default namelist variable values
  numOuterLoopIterations = 1
  limitHuInOuterLoop = .false.
  numIterMaxInnerLoop(:) = 0
  computeFinalNlJo = .false.

  if ( .not. utl_isNamelistPresent('NAMVAR','./flnml') ) then
    if ( mpi_myid == 0 ) then
      write(*,*) 'midas-var: namvar is missing in the namelist.'
      write(*,*) '           The default values will be taken.'
    end if

  else
    ! read in the namelist NAMVAR
    nulnam = 0
    ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    read(nulnam, nml=namvar, iostat=ierr)
    if( ierr /= 0) call utl_abort('midas-var: Error reading namelist')
    ierr = fclos(nulnam)
  end if
  if ( mpi_myid == 0 ) write(*,nml=namvar)

  if ( numOuterLoopIterations > maxNumberOfOuterLoopIterations ) then
    call utl_abort('midas-var: numOuterLoopIterations is greater than max value')
  end if

  if ( numOuterLoopIterations > 1 .and. &
       .not. all(numIterMaxInnerLoop(1:numOuterLoopIterations) > 0) ) then
    call utl_abort('midas-var: some numIterMaxInnerLoop(:) in namelist are negative or zero')
  end if

  obsMpiStrategy = 'LIKESPLITFILES'

  ! Initialize the Temporal grid
  call tim_setup

  ! Initialize observation file names and set datestamp
  call obsf_setup( dateStamp, varMode )
  if ( dateStamp > 0 ) then
    call tim_setDatestamp(datestamp)     ! IN
  else
    call utl_abort('var_setup: Problem getting dateStamp from observation file')
  end if

  ! Initialize constants
  if ( mpi_myid == 0 ) then
    call mpc_printConstants(6)
    call pre_printPrecisions
  end if

  ! Initialize the Analysis grid
  if (mpi_myid == 0) write(*,*)''
  if (mpi_myid == 0) write(*,*)'var_setup: Set hco parameters for analysis grid'
  call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

  if ( hco_anl % global ) then
    hco_core => hco_anl
  else
    ! Initialize the core (Non-Extended) analysis grid
    if (mpi_myid == 0) write(*,*)'var_setup: Set hco parameters for core grid'
    call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore' ) ! IN
  end if

  ! Initialisation of the analysis grid vertical coordinate from analysisgrid file
  call vco_SetupFromFile( vco_anl,        & ! OUT
                          './analysisgrid') ! IN

  call col_setVco(columnTrlOnAnlIncLev,vco_anl)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Setup and read observations
  call inn_setupObs(obsSpaceData, hco_anl, 'VAR', obsMpiStrategy, varMode) ! IN
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Basic setup of columnData module
  call col_setup
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !- Memory allocation for background column data
  call col_allocate(columnTrlOnAnlIncLev,obs_numheader(obsSpaceData),mpiLocal_opt=.true.)

  ! Initialize the observation error covariances
  call oer_setObsErrors(obsSpaceData, varMode) ! IN
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Initialize list of analyzed variables.
  call gsv_setup
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Reading trials
  call inn_getHcoVcoFromTrlmFile( hco_trl, vco_trl )
  allocHeightSfc = ( vco_trl%Vcode /= 0 )

  call gsv_allocate( stateVectorUpdateHighRes, tim_nstepobs, hco_trl, vco_trl,  &
                     dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                     mpi_distribution_opt='Tiles', dataKind_opt=pre_incrReal,  &
                     allocHeightSfc_opt=allocHeightSfc, hInterpolateDegree_opt='LINEAR', &
                     beSilent_opt=.false. )
  call gsv_zero( stateVectorUpdateHighRes )
  call gsv_readTrials( stateVectorUpdateHighRes )
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Initialize the background-error covariance, also sets up control vector module (cvm)
  call bmat_setup(hco_anl,hco_core,vco_anl)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Initialize the gridded variable transform module
  call gvt_setup(hco_anl,hco_core,vco_anl)

  ! Set up the minimization module, now that the required parameters are known
  ! NOTE: some global variables remain in minimization_mod that must be initialized before
  !       inn_setupColumnsOnTrlLev
  call min_setup( cvm_nvadim, hco_anl,                                   & ! IN
                  varqc_opt=varqcActive, nwoqcv_opt=numIterWithoutVarqc )  ! OUT
  allocate(controlVectorIncr(cvm_nvadim),stat=ierr)
  if (ierr /= 0) then
    write(*,*) 'var: Problem allocating memory for ''controlVectorIncr''',ierr
    call utl_abort('aborting in VAR')
  end if
  call utl_reallocate(controlVectorIncrSum,cvm_nvadim)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
  call tmg_stop(2)

  numInnerLoopIterDone = 0

  ! Enter outer-loop
  outer_loop: do outerLoopIndex = 1, numOuterLoopIterations
    write(*,*) 'var: start of outer-loop index=', outerLoopIndex

    ! Initialize stateVectorRefHeight for transforming TT/HU/P0 increments to
    ! height/pressure increments.
    if ( (gsv_varExist(stateVectorUpdateHighRes,'P_T') .and. &
          gsv_varExist(stateVectorUpdateHighRes,'P_M')) .or. &
         (gsv_varExist(stateVectorUpdateHighRes,'Z_T') .and. &
          gsv_varExist(stateVectorUpdateHighRes,'Z_M')) ) then

      call gvt_setupRefFromStateVector( stateVectorUpdateHighRes, 'height' )

      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    end if

    ! Horizontally interpolate high-resolution stateVectorUpdate to trial columns
    deallocInterpInfoNL = ( .not. numOuterLoopIterations > 1)
    call inn_setupColumnsOnTrlLev( columnTrlOnTrlLev, obsSpaceData, hco_core, &
                                   stateVectorUpdateHighRes, &
                                   deallocInterpInfoNL_opt=deallocInterpInfoNL )
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! Interpolate trial columns to analysis levels and setup for linearized H
    call inn_setupColumnsOnAnlIncLev( columnTrlOnTrlLev, columnTrlOnAnlIncLev )
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! Determine if to apply varqc to Jo of non-linear operator
    applyVarqcOnNlJo = ( varqcActive .and. numInnerLoopIterDone > numIterWithoutVarqc )
    if ( mpi_myid == 0 .and. applyVarqcOnNlJo ) write(*,*) 'applying varqc to non-linear Jo'

    ! Compute observation innovations and prepare obsSpaceData for minimization
    filterObsAndInitOer = ( outerLoopIndex == 1 )
    call inn_computeInnovation( columnTrlOnTrlLev, obsSpaceData, &
                                filterObsAndInitOer_opt=filterObsAndInitOer, &
                                applyVarqcOnNlJo_opt=applyVarqcOnNlJo )
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! Initialize stateVectorRefHU for doing variable transformation of the increments.
    if ( gsv_varExist(stateVectorUpdateHighRes,'HU') ) then
      applyLimitOnHU = ( limitHuInOuterLoop .and. outerLoopIndex > 1 )

      call gvt_setupRefFromStateVector( stateVectorUpdateHighRes, 'HU', &
                                        applyLimitOnHU_opt=applyLimitOnHU )

      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    end if

    ! Do minimization of cost function. Use numIterMaxInnerLoop from NAMVAR, instead of
    ! nitermax from NAMMIN, when numOuterLoopIterations > 1
    controlVectorIncr(:) = 0.0d0
    deallocHessian = ( numOuterLoopIterations == 1 )
    isMinimizationFinalCall = ( outerLoopIndex == numOuterLoopIterations )
    call min_minimize( outerLoopIndex, columnTrlOnAnlIncLev, obsSpaceData, controlVectorIncrSum, &
                       controlVectorIncr, numIterMaxInnerLoop(outerLoopIndex), &
                       deallocHessian_opt=deallocHessian, &
                       isMinimizationFinalCall_opt=isMinimizationFinalCall, &
                       numIterMaxInnerLoopUsed_opt=numIterMaxInnerLoopUsed )
    numInnerLoopIterDone = numInnerLoopIterDone + numIterMaxInnerLoopUsed
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! Accumulate control vector increments of all the previous iterations
    controlVectorIncrSum(:) = controlVectorIncrSum(:) + controlVectorIncr(:)

    ! Compute satellite bias correction increment and write to file on last outer-loop 
    ! iteration
    if ( outerLoopIndex == numOuterLoopIterations ) then
      call bcs_writebias(controlVectorIncr)
    end if

    call tvs_deallocateProfilesNlTlAd

    call gsv_allocate(stateVectorIncr, tim_nstepobsinc, hco_anl, vco_anl, &
         datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true., &
         dataKind_opt=pre_incrReal, allocHeight_opt=.false., allocPressure_opt=.false.)

    ! get final increment with mask if it exists
    call inc_getIncrement( controlVectorIncr, stateVectorIncr, cvm_nvadim )
    call gsv_readMaskFromFile( stateVectorIncr, './analysisgrid' )
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! Compute high-resolution analysis on trial grid
    call inc_computeHighResAnalysis( stateVectorIncr,                                  & ! IN
                                     stateVectorUpdateHighRes, stateVectorPsfcHighRes )  ! OUT
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! Impose limits on stateVectorUpdateHighRes only when outer loop is used.
    if ( limitHuInOuterLoop ) then
      write(*,*) 'var: impose limits on stateVectorUpdateHighRes'
      call qlim_saturationLimit( stateVectorUpdateHighRes )
      call qlim_rttovLimit( stateVectorUpdateHighRes )
    end if

    ! prepare to write incremnt when no outer-loop, or sum of increments at last
    ! outer-loop iteration.
    call tmg_start(6,'WRITEINCR')
    if ( numOuterLoopIterations > 1 .and. &
         outerLoopIndex == numOuterLoopIterations ) then

      call gsv_allocate(stateVectorIncrSum, tim_nstepobsinc, hco_anl, vco_anl, &
           datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true., &
           dataKind_opt=pre_incrReal, allocHeight_opt=.false., allocPressure_opt=.false.)
      call inc_getIncrement( controlVectorIncrSum, stateVectorIncrSum, cvm_nvadim )
      call gsv_readMaskFromFile( stateVectorIncrSum, './analysisgrid' )
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      ip3ForWriteToFile = numOuterLoopIterations
      call inc_writeIncrement( stateVectorIncrSum, &                     ! IN
                               ip3ForWriteToFile_opt=ip3ForWriteToFile ) ! IN
    else if ( numOuterLoopIterations == 1 ) then
      ip3ForWriteToFile = 0
      call inc_writeIncrement( stateVectorIncr, &                        ! IN
                               ip3ForWriteToFile_opt=ip3ForWriteToFile ) ! IN
    end if
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call tmg_stop(6)

    call gsv_deallocate( stateVectorIncr )

    write(*,*) 'var: end of outer-loop index=', outerLoopIndex
  end do outer_loop

  if ( gsv_isAllocated(stateVectorIncrSum) ) call gsv_deallocate( stateVectorIncrSum )

  ! Set the QC flags to be consistent with VAR-QC if control analysis
  if ( varqcActive ) call vqc_listrej(obsSpaceData)

  if ( computeFinalNlJo ) then
    ! Horizontally interpolate high-resolution stateVectorUpdate to trial columns
    call inn_setupColumnsOnTrlLev( columnTrlOnTrlLev, obsSpaceData, hco_core, &
                                   stateVectorUpdateHighRes )
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! Determine if to apply varqc to Jo of non-linear operator
    applyVarqcOnNlJo = ( varqcActive .and. numInnerLoopIterDone > numIterWithoutVarqc )
    if ( mpi_myid == 0 .and. applyVarqcOnNlJo ) write(*,*) 'applying varqc to non-linear Jo'

    ! Compute observation innovations and prepare obsSpaceData for minimization
    filterObsAndInitOer = .false.
    call inn_computeInnovation( columnTrlOnTrlLev, obsSpaceData, &
                                destObsColumn_opt=OBS_OMA, &
                                filterObsAndInitOer_opt=filterObsAndInitOer, &
                                applyVarqcOnNlJo_opt=applyVarqcOnNlJo )
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
  end if

  ! Memory deallocations for non diagonal R matrices for radiances
  call rmat_cleanup()

  ! Conduct obs-space post-processing diagnostic tasks (some diagnostic
  ! computations controlled by NAMOSD namelist in flnml)
  call osd_ObsSpaceDiag( obsSpaceData, columnTrlOnAnlIncLev, hco_anl )

  ! Deallocate memory related to B matrices and update stateVector
  call bmat_finalize()

  ! Deallocate structures needed for interpolation
  call s2c_deallocInterpInfo( inputStateVectorType='nl' )
  call s2c_deallocInterpInfo( inputStateVectorType='tlad' )

  ! Post processing of analyis before writing (variable transform+humidity clipping)
  call inc_analPostProcessing( stateVectorPsfcHighRes, stateVectorUpdateHighRes, &  ! IN 
                               stateVectorTrial, stateVectorPsfc, stateVectorAnal ) ! OUT
  call gsv_deallocate( stateVectorUpdateHighRes )

  ! compute and write the analysis (as well as the increment on the trial grid)
  call tmg_start(18,'ADDINCREMENT')
  call inc_writeIncAndAnalHighRes( stateVectorTrial, stateVectorPsfc, &
                                   stateVectorAnal )
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
  call tmg_stop(18)

  if (mpi_myid == 0) then
    clmsg = 'REBM_DONE'
    call utl_writeStatus(clmsg)
  end if

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

  ! Job termination
  istamp = exfin('VAR','FIN','NON')

  if (mpi_myid == 0) then
    clmsg = 'VAR3D_END'
    call utl_writeStatus(clmsg)
  end if

  call tmg_stop(1)

  call tmg_terminate(mpi_myid, 'TMG_VAR' )

  call rpn_comm_finalize(ierr) 

end program midas_var

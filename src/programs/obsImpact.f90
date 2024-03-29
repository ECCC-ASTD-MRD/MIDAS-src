program midas_obsImpact
  !
  !:Purpose: Main program for computing the Forecast Sensitivity to Observation Impact (FSOI)
  !
  !           ---
  !
  !:Algorithm: FSOI partitions the forecast error reduction within the current system from 
  !            assimilating the observations. The total forecast error reduction defined as:
  !
  !             :math:`(e_{t}^{fa})^{T}*C*(e_{t}^{fa})-(e_{t}^{fb})^{T}*C*(e_{t}^{fb})`
  !
  !             where :math:`e_{t}^{fa}=M(x_{0}^{a})-x_{t}^{a}`
  !
  !             and :math:`e_{t}^{fb}=M(x_{0}^{b})-x_{t}^{a}`
  !
  !            and C is the energy norm. 
  !             
  !            --
  !
  !            In this program, there are the hybrid FSOI approach (HFSOI) and the ensemble FSOI (EFSOI) approaches.
  !            The hybrid approach is appropriate for computing the impact of observations assimilated with 4D-EnVar.
  !            It combines the ensemble approach for propagating the sensitivities from forecast to analysis time and
  !            the variational approach for the adjoint of the assimilation procedure. The ensemble approach is 
  !            appropriate for computing the impact of observations assimilated with the LETKF. It relies purely on
  !            analysis and forecast ensemble for the calculation. 
  !
  !
  !            More details on HFSOI can be found in the paper: `HFSOI approach <https://doi.org/10.1175/MWR-D-17-0252.1>`_
  !
  !
  !            More details on EFSOI can be found in the paper:`EFSOI approach <http://doi.org/10.3402/tellusa.v65i0.20038>`_
  !
  !            --
  !
  !:File I/O: The required input files and produced output files are listed as follows.
  !
  !           --
  !
  !============================================== ====================================================================
  ! Input and Output Files                         Description of file
  !============================================== ====================================================================
  ! ``flnml``                                      In - Main namelist file with parameters user may modify
  ! ``trlm_$NN`` (e.g. ``trlm_01``)                In - Background state (a.k.a. trial) files for each timestep
  ! ``analysisgrid``                               In - File defining grid for computing the analysis increment
  ! ``bgcov``                                      In - Static (i.e. NMC) B matrix file for NWP fields
  ! ``ensemble/$YYYYMMDDHH_006_$NNNN``             In - Ensemble member files defining ensemble B matrix
  ! ``ensemble/$YYYYMMDDHH_foradv``                In - advection files
  ! ``obsfiles_$FAM/obs$FAM_$NNNN_$NNNN``          In - Observation file for each "family" and MPI task
  ! ``obserr``                                     In - Observation error statistics
  ! ``forecasts/forecast_*``                       In - Global deterministic forecast from backgroud and analysis
  ! ``forecasts/analysis``                         In - Verifying analysis
  ! ``diafiles_$FAM.updated/obs$FAM_$NNNN_$NNNN``  Out - Updated obs file for each "family" and MPI task in SQLite 
  ! Remainder are files related to radiance obs:
  ! ``stats_$SENSOR_assim``                        In - Satellite radiance observation errors of difference sensors
  ! ``stats_tovs``                                 In - Satellite radiance observation errors 
  ! ``rtcoef_$PLATFORM_$SENSOR.**``                In - RTTOV coefficient files 
  ! ``rttov_h2o_limits.dat``                       In - Min/Max humidity limits 
  ! ``stats_tovs_symmetricObsErr``                 In - user-defined symmetric TOVS errors for all sky
  ! ``Cmat_$PLATFORM_$SENSOR.dat``                 In - Inter-channel observation-error correlations
  !============================================== ====================================================================
  !
  !           --
  !
  !:Synopsis: Below is a summary of the ``obsImpact`` program calling sequence:
  !
  !           - **Initial setups:**
  !
  !             - Initialize FSOI module, read the namelist from ``fso_setup``
  !
  !             - Setup horizontal and vertical grid objects for "analysis
  !               grid" from ``analysisgrid`` file and for "trial grid" from
  !               first trial file: ``trlm_01``.
  !
  !             - Setup ``obsSpaceData`` object and read observations from
  !               files: ``inn_setupObs``.
  !
  !             - Setup ``columnData`` and ``gridStateVector`` modules (read
  !               list of analysis variables from namelist) and allocate column
  !               object for storing trial on analysis levels.
  ! 
  !             - Setup the observation error statistics in ``obsSpaceData``
  !               object: ``oer_setObsErrors``.
  !
  !             - Allocate a stateVector object on the trial grid and then
  !               read the trials: ``gio_readTrials``.
  !
  !             - Setup the B matrices: ``bmat_setup``.
  !
  !             - Setup the ``gridVariableTransforms``. 
  !
  !           - **Computation:** 
  !
  !             - ``incdatr``: Setup dateStamp(one extra timestamp: leadTime) for FSOI.
  !
  !             - ``inn_computeInnovation``: Compute innovation 
  !
  !             - ``calcFcstError``: Read the forcasts from background and analysis, and 
  !                                  the verifying analysis as well, calculate the forecast
  !                                  error norm defined as: 
  !                                  :math:`(C*(e_{t}^{fa}+e_{t}^{fb}))`.
  !
  !             - ``bmat_sqrtBT``: Compute the variable :math:`\hat{v}` for minimization in HFSO mode.             
  !                                       :math:`\hat{v}=B_{t}^{T/2}*C*(e_{t}^{fa}+e_{t}^{fb})` 
  !
  !             - ``minimize``: Do the minimization to apply the adjoint of the 4D-EnVar assimilation  
  !                             to :math:`\hat{v}` with the result being :math:`\hat{a}`.
  !
  !             - ``bmat_sqrt``: Compute :math:`B^{1/2}*\hat{a}` .
  !
  !             - ``s2c_tl``, ``oop_Htl``: Apply the observation operators :math:`H*B^{1/2}*\hat{a}` .
  !
  !             - ``rmat_RsqrtInverseAllObs``: Multiply by the inverse of the observation error variances
  !                                                    :math:`R^{-1}*H*B^{1/2}*\hat{a}` 
  !
  !             - ``obs_bodySet_r``: Multiply the resulting sensitivity value by the innovation and put  
  !                                  the result in the ``obsSpaceDate`` column ``OBS_FSO`` so it can be stored in
  !                                  the observation files 
  !
  !             - ``sumFSO``: Print out the FSOI value, including total and the one from each obs family.
  !
  !           - **Final steps:** 
  !
  !             - ``obsf_writeFiles``: Write the final FSOI value in SQLite file.
  !
  !           --
  !
  !:Options: `List of namelist blocks <../namelists_in_each_program.html#obsimpact>`_
  !          that can affect the ``obsImpact`` program.
  !
  !          * The use of ``obsImpact`` program is controlled mainly by the namelist block
  !            ``&NAMFSO`` read by the ``fso_mod`` module. 
  ! 
  !          * Some of the other relevant namelist blocks used to configure FSOI
  !            are listed in the following table:
  !
  !   
  !========================= ====================== =============================================================
  ! Module                   Namelist               Description of what is controlled
  !========================= ====================== =============================================================
  ! ``timeCoord_mod``        ``NAMTIME``            assimilation time window length, temporal resolution of
  !                                                 the background state and increment
  ! ``bMatrixEnsemble_mod``  ``NAMBEN``             weight and other parameters for ensemble-based B matrix
  !                                                 component
  !========================= ====================== =============================================================
  !
  !
  use version_mod
  use codePrecision_mod
  use ramDisk_mod
  use utilities_mod
  use midasMpi_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use columnData_mod
  use obsSpaceData_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use bMatrix_mod
  use innovation_mod
  use obsFiles_mod
  use obsErrors_mod
  use gridVariableTransforms_mod
  use fsoi_mod

  implicit none

  integer :: istamp,exdb,exfin,ierr, dateStampFromObs

  type(struct_obs),       target :: obsSpaceData
  type(struct_columnData),target :: columnTrlOnAnlIncLev
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
  call mmpi_initialize

  call tmg_init(mmpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')

  if (mmpi_myid == 0) then
    call utl_writeStatus('VAR3D_BEG')
  end if

  call ram_setup

  !
  !- 1. Settings 
  !
  obsColumnMode  = 'VAR'
  obsMpiStrategy = 'LIKESPLITFILES'

  !
  !- Initialize the Temporal grid and set dateStamp from env variable
  !
  call tim_setup()
  !     
  !- Initialize burp file names and set datestamp if not already
  !
  call obsf_setup(dateStampFromObs, 'FSO')
  if (tim_getDateStamp() == 0) then
    if (dateStampFromObs > 0) then
      call tim_setDatestamp(dateStampFromObs)
    else
      call utl_abort('obsImpact: DateStamp was not set')
    end if
  end if
  !
  !- Initialize constants
  !
  if ( mmpi_myid == 0 ) then
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
  if (mmpi_myid.eq.0) write(*,*)''
  if (mmpi_myid.eq.0) write(*,*)'obsImpact: Set hco parameters for analysis grid'
  call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

  !- Do FSO module set up
  call fso_setup(hco_anl)

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
  call col_allocate(columnTrlOnAnlIncLev,obs_numheader(obsSpaceData))

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
    if (mmpi_myid == 0) call obsf_writeFiles(obsSpaceData)
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

  if (mmpi_myid == 0) then
    call utl_writeStatus('VAR3D_END')
  endif

  call utl_tmg_stop(0)

  call tmg_terminate(mmpi_myid, 'TMG_INFO')

  call rpn_comm_finalize(ierr)

end program midas_obsImpact 

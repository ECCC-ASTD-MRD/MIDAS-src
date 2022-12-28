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
  !:Purpose: Main program for Observation Impact computation (FSOI)
  !
  !           ---
  !
  !:Algorithm: FSOI is to partition, with respect to arbitrary subsets 
  !            of observations, the forecast error reduction, defined as:
  !
  !             :math:`(e_{t}^{fa})^{T}*C*(e_{t}^{fa})-(e_{t}^{fb})^{T}*C*(e_{t}^{fb})`
  !
  !            (within the current system) from assimilating these observations.
  !             
  !            --
  !
  !            The hybrid FSOI approach (HFSO) is developed to suit 4D-EnVar 
  !            without requiring an adjoint of the forecast model. It is 
  !            similar to the standard adjoint approach as it uses an iterative
  !            approach to represent the adjoint of the data assimilation
  !            procedure, but relies on an ensemble of forecasts to propagate 
  !            the impact between forecast and analysis times. In this way, the
  !            impact of assimilated observations on 24h forecasts quality is
  !            estimated by combining the impact of observations on the 
  !            analysis increment (e.g. via the 4D-EnVar) and the impact of
  !            analysis increments on reducing forecast errors.
  !
  !            --
  !
  !            More details on HFSO can be found in the paper: `HFSO approach <https://doi.org/10.1175/MWR-D-17-0252.1>`_
  !
  !            --
  !
  !            Ensemble Forecast Sensitivity to Observation (EFSO) is an approach
  !            for estimating the observation impacts in the ensemble-based data 
  !            assimilation systeman. This approach is convenient and no need of 
  !            adjoint of the forecast model. 
  !            
  !            --
  ! 
  !            More details on EFSO can be found in the paper:`EFSO approach <http://doi.org/10.3402/tellusa.v65i0.20038>`_
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
  !             - ``bmat_sqrtBT``: Compute the variable vhat for minimization in HFSO mode.             
  !                                       :math:`vhat=B_{t}^{T/2}*C*(e_{t}^{fa}+e_{t}^{fb})` 
  !
  !             - ``minimize``: Do the minimization of vhat for HFSO mode.
  !
  !             - ``bmat_sqrt``: Compute :math:`B^{1/2}*ahat` .
  !
  !             - ``s2c_tl``, ``oop_Htl``: Put in columndata and get :math:`H*B^{1/2}*ahat` .
  !
  !             - ``rmat_RsqrtInverseAllObs``: Apply observation error variances and get
  !                                                    :math:`R^{-1}*H*B^{1/2}*ahat` 
  !
  !             - ``sumFSO``: Print out the FSOI value, including total and the one from each obs family.
  !
  !           - **Final steps:** 
  !
  !             - ``obsf_writeFiles``: Write the final FSOI value in SQLite file.
  !
  !           --
  !
  !:Options: `List of namelist blocks <../namelists_in_each_program.html#obsImpact>`_
  !          that can affect the ``obsImpact`` program.
  !
  !          * The use of ``obsImpact`` program is controlled mainly by the namelist block
  !            ``&NAMFSO`` read by the ``fso_mod`` module. Here are more details about the variables in ``&NAMFSO`` 
  !
  !          * ``LEADTIME``: the forecast leading time 
  !
  !          * ``NVAMAJ``, ``NITERMAX``, ``NSIMMAX``, ``REPSG``: variables for minimization 
  !
  !          * ``latMinNorm``, ``latMaxNorm``, ``lonMinNorm``, ``lonMaxNorm`` : define the domain of forecast errors
  !
  !          * ``FORCECASTPATH``: the path contains the forecasts from background and analysis and verifying analysis
  !
  !          * ``FSOMODE``:  use HFSO (hybrid FSOI) or EFSO (ensemble FSOI), default is HFSO
  ! 
  !          * ``includeHUnorm``: use wet norm or not, default value is dry norm 
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
  !                                                 the background state and increment
  ! ``bMatrixEnsemble_mod``  ``NAMBEN``             weight and other parameters for ensemble-based B matrix
  !                                                 component
  ! ``bMatrixHI_mod``        ``NAMBHI``             weight and other parameters for the climatological B matrix
  !                                                 component based on homogeneous-isotropic covariances
  !                                                 represented in spectral space
  ! ``obsFiles_mod``         ``NAMWRITEDIAG``       write out in SQLite file
  !========================= ====================== =============================================================
  !
  !
  use version_mod
  use codePrecision_mod
  use ramDisk_mod
  use utilities_mod
  use midasMpi_mod
  use MathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use columnData_mod
  use obsSpaceData_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use bmatrix_mod
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

end program midas_obsimpact 

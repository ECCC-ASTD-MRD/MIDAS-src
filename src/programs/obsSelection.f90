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

program midas_obsSelection
  !
  !:Purpose: Main program for background check, and thinning of all observation types.
  !
  !          ---
  !
  !:Algorithm: For background check of each observation type, several quality control tests 
  !            specific to that observation type are performed on each available observation
  !            to determine if the observation meets the standard to be assimilated later.
  !            Certain bit numbers of the observation flag are modified depending on whether
  !            the observation is assimilable or if not, the reason for which the observation 
  !            was flagged for rejection. The thinning is performed in next step to reduce 
  !            the number of observation data for assimilation. This is to reduce the 
  !            1-computational cost, 2-observation error correlation during 
  !            assimilation stage. The background-checked thinned observations are written
  !            to new files, ready for assimilation.
  !
  !            --
  !
  !            The computed bias correction values are applied to the observation before the 
  !            background check step for certain observation types (e.g. radiances). One of quality
  !            control tests during background check is the ``rogue check`` which is the allowed
  !            distance of the observation from the background state. Innovation vector  
  !            ``y-H(xb)`` is required to measure this distance. The innovations are 
  !            computed at the beginning of the program before the background check starts.
  !
  !            --
  !
  !============================================== ==============================================================
  ! Input and Output Files                        Description of file
  !============================================== ==============================================================
  ! ``flnml``                                      In - Main namelist file with parameters user may modify
  ! ``flnml_static``                               In - The "static" namelist that should not be modified
  ! ``trlm_$NN`` (e.g. ``trlm_01``)                In - Background state (a.k.a. trial) files for each timestep
  ! ``analysisgrid``                               In - File defining grid for computing the analysis increment
  ! ``obsfiles_$FAM/obs$FAM_$NNNN_$NNNN``          In - Observation file for each "family" and MPI task
  ! ``obscov``                                     In - Observation error statistics
  ! ``obserr``                                     In - Observation error statistics
  ! ``obsfiles_$FAM.updated/obs$FAM_$NNNN_$NNNN``  Out - Updated obs file for each "family" and MPI task
  ! Remainder are files related to radiance obs:
  ! ``stats_$SENSOR_assim``                        In - Satellite radiance observation errors of difference sensors
  ! ``stats_tovs``                                 In - Satellite radiance observation errors
  ! ``stats_tovs_symmetricObsErr``                 In - User-defined symmetric TOVS errors for all sky
  ! ``Cmat_$PLATFORM_$SENSOR.dat``                 In - Inter-channel observation-error correlations
  ! ``dynbcor.coeffs.$SENSOR.*.coeffs_$SENSOR``    In - Dynamic bias correction file
  ! ``ceres_global.std``                           In - Surface emmissivity and type?
  ! ``champ_fd_181x91``                            In - NOT USED?
  ! ``rtcoef_$PLATFORM_$SENSOR.dat``               In - RTTOV coefficient files
  ! ``rttov_h2o_limits.dat``                       In - Min/max humidity limits applied to analysis
  ! ``ozoneclim98``                                In - Ozone climatology
  !============================================== ==============================================================
  !
  !           --
  !
  !:Synopsis: Below is a summary of the ``obsSelection`` program calling sequence:
  !
  !             - **Initial setups:**
  !
  !               - Read the NAMOBSSELECTION namelist and check/modify some values.
  !
  !               - Various modules are setup: ``obsFiles_mod``, ``timeCoord_mod``.
  !
  !               - Setup horizontal and vertical grid objects for "analysis
  !                 grid" from ``analysisgrid`` file.
  !
  !               - Setup ``obsSpaceData`` object and read observations from
  !                 files: ``inn_setupObs``.
  !
  !               - Compute and update the stored surface type for some satellite
  !                 radiance instruments. 
  !
  !               - Setup ``columnData`` module (read list of analysis variables 
  !                 from namelist) and allocate column object for storing trial 
  !                 on analysis levels.
  !
  !               - Setup the observation error statistics in ``obsSpaceData``
  !                 object: ``oer_setObsErrors``
  !
  !               - Setup ``gridStateVector`` module.
  !
  !               - Applying optional bias corrections to some observation types.
  !
  !               - Setup horizontal and vertical grid objects for "trial grid" from
  !                 first trial file: ``trlm_01``.
  !
  !               - Allocate a stateVector object on the trial grid and then
  !                 read the trials: ``gio_readTrials``.
  !
  !             - **Computation**
  !
  !               - Compute ``columnTrlOnTrlLev`` and ``columnTrlOnAnlIncLev`` from
  !                 background state: ``inn_setupColumnsOnTrlLev``,
  !                 ``inn_setupColumnsOnAnlIncLev``.
  !
  !               - Compute innovation from background state: ``inn_computeInnovation``.
  !
  !               - Do background check for conventional observation: 
  !                 ``bgck_bgCheck_conv``.
  !
  !               - Update radiance bias correction in ``obsSpaceData`` and apply 
  !                 the bias corrections to the observations and innovations for
  !                 radiances: ``bcs_calcBias``, ``bcs_applyBiasCorrection``.
  !
  !               - Perform background check for multiple observation types.
  !             
  !               - If thinning was requested: 
  !
  !                 - add some cloud parameters and set missing observation flags in 
  !                   observation files (specific for radiances) and 
  !
  !                 - perform thinning for different observation types.
  !
  !               - Write the final background-checked bias corrected results (either
  !                 thinned or not thinned) into the observation file.
  ! 
  !               - If thinning was requested, remove observations which were flagged 
  !                 not to be assimilated from the observation file.
  !
  !             --
  !
  !:Options: `List of namelist blocks <../namelists_in_each_program.html#obsSelection>`_
  !          that can affect the ``obsSelection`` program.
  !
  !          * The use of ``obsSelection`` program is controlled by the namelist block
  !           ``&NAMOBSSELECTION`` read by the ``obsSelection`` program.
  !
  !          * Some of the other relevant namelist blocks used to configure the
  !            ``obsSelection`` are listed in the following table:
  ! 
  !=========================== ========================= =========================================
  ! Module                      Namelist                  Description of what is controlled
  !=========================== ========================= =========================================
  ! ``midas_obsSelection``      ``NAMOBSSELECTION``       if or not thinning is done.
  ! ``biasCorrectionConv_mod``  ``NAMBIASCONV``           variables to performs bias correction 
  !                                                       for conventional observations.
  ! ``biasCorrectionConv_mod``  ``NAMSONDETYPES``         additional variables to performs bias 
  !                                                       correction for radiosondes conventional 
  !                                                       observations.
  ! ``backgroundCheck_mod``     ``NAMBGCKCONV``           variables to perform background check 
  !                                                       for conventional observations.
  ! ``SSTbias_mod``             ``NAMSSTBIASESTIMATE``    variables to perform bias correction 
  !                                                       for satellite SST.
  ! ``biasCorrectionSat_mod``   ``NAMBIASSAT``            variables to perform bias correction 
  !                                                       for satellite radiances.
  ! ``multi_ir_bgck_mod``       ``NAMBGCKIR``             Variables to perform background check 
  !                                                       for hyperspectral infrared radiances.
  ! ``bgckmicrowave_mod``       ``NAMBGCK``               Variables to perform background check 
  !                                                       for microwave radiances.
  ! ``bgckcsr_mod``             ``NAMCSR``                Variables To perform background check 
  !                                                       for CSR radiances.
  ! ``bgckssmis_mod``           ``NAMBGCK``               Variables to perform background check 
  !                                                       for SSMIS radiances.
  ! ``bgckOcean_mod``           ``NAMOCEANBGCHECK``       Variables to perform background check 
  !                                                       for ocean data.
  ! ``bgckOcean_mod``           ``NAMICEBGCHECK``         Variables to perform background check 
  !                                                       for the SeaIce data.
  ! ``burpread_mod``            ``NAMADDTOBURP``          element IDs to add to the BURP file
  ! ``thinning_mod``            ``THIN_HYPER``            variables to perform thinning on 
  !                                                       hyperspectral infrared radiances.
  ! ``thinning_mod``            ``THIN_TOVS``             variables to perform thinning on 
  !                                                       microwave radiances.
  ! ``thinning_mod``            ``thin_csr``              variables to perform thinning on 
  !                                                       CSR radiances.
  ! ``thinning_mod``            ``thin_raobs``            variables to perform thinning on 
  !                                                       radiosonde observations.
  ! ``thinning_mod``            ``thin_scat``             variables to perform thinning on 
  !                                                       scatterometer wind observations.
  ! ``thinning_mod``            ``thin_aircraft``         variables to perform thinning on 
  !                                                       aircraft observations.
  ! ``thinning_mod``            ``thin_surface``          variables to perform thinning on 
  !                                                       surface observations.
  ! ``thinning_mod``            ``thin_gbgps``            variables to perform thinning on 
  !                                                       ground-based GPS observations.
  ! ``thinning_mod``            ``thin_gpsro``            variables to perform thinning on 
  !                                                       GPS radio-occultation observations.
  ! ``thinning_mod``            ``thin_aladin``           variables to perform thinning on 
  !                                                       aladin wind observations.
  ! ``thinning_mod``            ``thin_aladin``           variables to perform thinning on 
  !                                                       aladin wind observations.
  ! ``timeCoord_mod``           ``NAMTIME``               assimilation time window length, 
  !                                                       temporal resolution of the background 
  !                                                       state.
  !=========================== ========================= =========================================
  !
  use version_mod
  use codePrecision_mod
  use ramDisk_mod
  use midasMpi_mod
  use utilities_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use backgroundCheck_mod
  use multi_ir_bgck_mod
  use innovation_mod
  use obsSpaceData_mod
  use columnData_mod
  use obsFiles_mod
  use obsErrors_mod
  use biasCorrectionSat_mod
  use biasCorrectionConv_mod
  use thinning_mod
  use bgckmicrowave_mod
  use bgckssmis_mod
  use bgckcsr_mod
  use bgckOcean_mod 
  use SSTbias_mod
   
  implicit none

  integer :: dateStampFromObs, headerIndex, ierr, nulnam
  type(struct_columnData),target :: columnTrlOnAnlIncLev
  type(struct_columnData),target :: columnTrlOnTrlLev
  type(struct_obs),       target :: obsSpaceData
  type(struct_gsv)               :: stateVectorTrialHighRes
  type(struct_hco),      pointer :: hco_anl => null()
  type(struct_vco),      pointer :: vco_anl => null()
  type(struct_hco),      pointer :: hco_trl => null()
  type(struct_vco),      pointer :: vco_trl => null()
  type(struct_hco),      pointer :: hco_core => null()

  logical :: allocHeightSfc
  integer :: get_max_rss, fnom, fclos

  ! Namelist variables
  logical                        :: doThinning  ! Control whether or not thinning is done

  namelist /namObsSelection/ doThinning

  call ver_printNameAndVersion('obsSelection','Obs Quality Control and Thinning')

  !- 1.0 mpi
  call mmpi_initialize

  !- 1.1 timings
  call tmg_init(mmpi_myid, 'TMG_INFO')
  call utl_tmg_start(0,'Main')

  !- 1.2 Read the namelist for obsSelection program (if it exists)
  ! set default values for namelist variables
  doThinning = .false.
  if (utl_isNamelistPresent('namObsSelection', './flnml')) then
    nulnam = 0
    ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    if (ierr /= 0) call utl_abort('midas-obsSelection: Error opening file flnml')
    read(nulnam, nml = namObsSelection, iostat = ierr)
    if (ierr /= 0) call utl_abort('midas-obsSelection: Error reading namelist namObsSelection')
    if (mmpi_myid == 0) write(*,nml = namObsSelection)
    ierr = fclos(nulnam)
  else
    write(*,*)
    write(*,*) 'midas-obsSelection: Namelist block namObsSelection is missing in the namelist.'
    write(*,*) '                    The default value will be taken.'
    if (mmpi_myid == 0) write(*, nml = namObsSelection)
  end if

  !
  !- Initialize the ram disk
  !
  call ram_setup

  !     
  !- Initialize observation file names, but don't use datestamp
  !
  call obsf_setup(dateStampFromObs, 'bgck')

  !
  !- Initialize the Temporal grid and dateStamp from trial file
  !
  call tim_setup(fileNameForDate_opt = './trlm_01')

  !
  !- Initialize constants
  !
  if (mmpi_myid == 0) then
    call mpc_printConstants(6)
    call pre_printPrecisions
  end if

  !
  !- Initialize the Analysis grid
  !
  if(mmpi_myid == 0) write(*,*)
  if(mmpi_myid == 0) write(*,*) 'midas-obsSelection: Set hco parameters for analysis grid'
  call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

  if ( hco_anl % global ) then
    hco_core => hco_anl
  else
    !- Initialize the core (Non-Extended) analysis grid
    if( mmpi_myid == 0) write(*,*) 'midas-obsSelection: Set hco parameters for core grid'
    call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore' ) ! IN
  end if

  !     
  !- Initialisation of the analysis grid vertical coordinate from analysisgrid file
  !
  call vco_SetupFromFile(vco_anl, './analysisgrid')

  call col_setVco(columnTrlOnAnlIncLev, vco_anl)
  write(*,*) 'Memory Used: ', get_max_rss()/1024,'Mb'

  !
  !- Setup and read observations
  !
  call inn_setupObs(obsSpaceData, hco_anl, 'ALL', 'LIKESPLITFILES', 'bgck')
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! if ssmis, compute the surface type ele and update obspacedata

  call ssbg_computeSsmisSurfaceType(obsSpaceData)

  ! Only for MWHS2 data and if modLSQ option is set to .true. in nambgck namelist, set values for land qualifier
  ! indice and terrain type based on calculations
  if (obs_famExist(obsSpaceData, 'TO')) then
    call mwbg_computeMwhs2SurfaceType(obsSpaceData)
  end if

  !
  !- Basic setup of columnData module
  !
  call col_setup
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- Memory allocation for background column data
  !
  call col_allocate(columnTrlOnAnlIncLev, obs_numheader(obsSpaceData), mpiLocal_opt = .true.)

  !
  !- Initialize the observation error covariances
  !
  call oer_setObsErrors(obsSpaceData, 'bgck')

  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  ! Initialize list of analyzed variables.
  !
  call gsv_setup

  ! Apply optional bias corrections
  if (obs_famExist(obsSpaceData, 'TM')) then
    call sstb_applySatelliteSSTBiasCorrection(obsSpaceData, hco_anl, &
                                              vco_anl, columnTrlOnAnlIncLev)
  end if  

  if (obs_famExist(obsSpaceData, 'AI')) call bcc_applyAIBcor(obsSpaceData)    
  if (obs_famExist(obsSpaceData, 'GP')) call bcc_applyGPBcor(obsSpaceData)
  if (obs_famExist(obsSpaceData, 'UA')) call bcc_applyUABcor(obsSpaceData)
    
  ! Reading trials
  call inn_getHcoVcoFromTrlmFile(hco_trl, vco_trl)
  allocHeightSfc = (vco_trl%Vcode /= 0)

  call gsv_allocate(stateVectorTrialHighRes, tim_nstepobs, hco_trl, vco_trl,  &
                    dateStamp_opt = tim_getDateStamp(), mpi_local_opt = .true., &
                    mpi_distribution_opt = 'Tiles', dataKind_opt = 4, &
                    allocHeightSfc_opt = allocHeightSfc, &
                    hInterpolateDegree_opt = 'LINEAR', beSilent_opt=.false.)
  call gsv_zero(stateVectorTrialHighRes)
  call gio_readTrials(stateVectorTrialHighRes)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Horizontally interpolate trials to trial columns
  call inn_setupColumnsOnTrlLev(columnTrlOnTrlLev, obsSpaceData, hco_core, &
                                stateVectorTrialHighRes)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Interpolate trial columns to analysis levels and setup for linearized H
  call inn_setupColumnsOnAnlIncLev(columnTrlOnTrlLev, columnTrlOnAnlIncLev)

  ! Compute observation innovations and prepare obsSpaceData for minimization
  call inn_computeInnovation(columnTrlOnTrlLev, obsSpaceData, analysisMode_opt = .false.)

  ! 2.2 Perform the background check
  !     The routine also calls compute_HBHT and writes to listings & obsSpaceData

  ! Do the conventional data background check
  call bgck_bgCheck_conv(columnTrlOnAnlIncLev, columnTrlOnTrlLev, hco_anl, obsSpaceData)

  if (obs_famExist(obsSpaceData, 'TO')) then

    ! Satellite radiance bias correction
    call bcs_calcBias(obsSpaceData, columnTrlOnTrlLev)        ! Fill in OBS_BCOR obsSpaceData column 
                                                              ! with computed bias correction
    call bcs_applyBiasCorrection(obsSpaceData, obs_var, 'TO') ! Apply bias correction to OBS
    call bcs_applyBiasCorrection(obsSpaceData, obs_omp, 'TO') ! Apply bias correction to O-F

    ! Do the TO background check
    call irbg_bgCheckIR(columnTrlOnTrlLev, obsSpaceData)
    call mwbg_bgCheckMW(obsSpaceData)
    call csrbg_bgCheckCSR(obsSpaceData)
    call ssbg_bgCheckSsmis(obsSpaceData)

  end if

  ! Do the ocean data background check
  if (obs_famExist(obsSpaceData, 'TM')) call ocebg_bgCheckSST(obsSpaceData, tim_getDateStamp(), &
                                                              columnTrlOnTrlLev, hco_trl)

  ! Do the sea ice data gross background check
  if (obs_famExist(obsSpaceData, 'GL')) call ocebg_bgCheckSeaIce(obsSpaceData)

  if (doThinning) then

    ! Copy original obs files into another directory
    call obsf_copyObsDirectory('./obsOriginal', direction = 'TO')

    ! 2.3 Write obs files after background check, but before thinning
    call obsf_writeFiles(obsSpaceData, writeDiagFiles_opt = .false.)
    
    ! Add cloud parameter data to burp files (AIRS,IASI,CrIS,ATMS,AMSUA,...)
    if (obs_famExist(obsSpaceData, 'TO')) then
      call obsf_updateMissingObsFlags(obsSpaceData)
      call obsf_addCloudParametersAndEmissivity(obsSpaceData)
    end if

    ! Copy the pre-thinning files into another directory
    call obsf_copyObsDirectory('./obsBeforeThinning', direction = 'TO')

    ! Copy original obs files back into usual directory
    call obsf_copyObsDirectory('./obsOriginal', direction = 'FROM')

    ! 2.4 Thinning:  Set bit 11 of flag, one observation type at a time
    call thn_thinHyper(obsSpaceData)
    call thn_thinTovs(obsSpaceData)
    call thn_thinCSR(obsSpaceData)
    call thn_thinRaobs(obsSpaceData)
    call thn_thinScat(obsSpaceData)
    call thn_thinSatWinds(obsSpaceData)
    call thn_thinAircraft(obsSpaceData)
    call thn_thinSurface(obsSpaceData, 'SF') ! surface data thinning
    if (obs_famExist(obsSpaceData, 'TM')) then
      call thn_thinSurface(obsSpaceData, 'TM') ! SST thinning
      call thn_thinSatSST(obsSpaceData)        ! satellite SST thinning
    end if      
    call thn_thinGbGps(obsSpaceData)
    call thn_thinGpsRo(obsSpaceData)
    call thn_thinAladin(obsSpaceData)
    ! if requested, dump the thinned predictors and coefficients to sqlite
    call bcs_dumpBiasToSqliteAfterThinning(obsSpaceData)

  end if

  
  ! 3 Write the final results

  ! 3.1 Into the listings
  write(*,*)
  write(*,*) '> midas-obsSelection: printing the FIRST header and body'
  do headerIndex = 1, min(1,obs_numHeader(obsSpaceData))
    call obs_prnthdr(obsSpaceData, headerIndex)
    call obs_prntbdy(obsSpaceData, headerIndex)
  end do

  ! 3.2 Into the observation files
  write(*,*)
  write(*,*) '> midas-obsSelection: writing to file'
  call obsf_writeFiles(obsSpaceData)
  
  ! Add cloud parameter data to burp files (AIRS,IASI,CrIS,...)
  if (obs_famExist(obsSpaceData, 'TO')) then
    call obsf_updateMissingObsFlags(obsSpaceData)
    call obsf_addCloudParametersAndEmissivity(obsSpaceData)
  end if

  ! cleaning the observation files
  if (doThinning) call obsf_cleanObsFiles()

  !
  ! 4.  Ending
  !
  write(*,*)
  write(*,*) '> midas-obsSelection: Ending'
  call obs_finalize(obsSpaceData) ! deallocate obsSpaceData

  call rpn_comm_finalize(ierr)

  call utl_tmg_stop(0)
  call tmg_terminate(mmpi_myid, 'TMG_INFO')

end program midas_obsSelection

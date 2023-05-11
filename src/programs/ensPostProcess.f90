program midas_ensPostProcess
  !
  !:Purpose: Post-processing program for the local ensemble transform Kalman filter (LETKF).
  !          Many aspects of this program are controlled throught the namelist block
  !          namEnsPostProc defined in epp_postProcess.
  !
  !          ---
  !
  !:Algorithm: The ``ensPostProcess`` program performs several post-processing tasks of
  !            LETK ensemble analyses and the ensemble trials. The following tasks are
  !            performed based on the namelist options.
  !
  !              - **Mean and std**: Ensemble mean and standard deviation are computed
  !                upon reading the ensemble members and after the processes which may
  !                alter the mean or the standard deviation.
  !
  !              - **RTPP**: Relaxation to the prior perturbation will relax the analysis
  !                perturbation to the trial perturbation with a given amount (:math:`\alpha`) in the namelist.
  !                It can be formulated as:
  !                :math:`x^a_i \leftarrow \overline{x^a} + (1-\alpha)(\overline{x^a}-x^a_i)+\alpha(\overline{x^b}-x^b_i)`
  !
  !                More details for the RTPP can be found in paper:
  !                `Impacts of Initial Estimate and Observation Availability on Convective-Scale
  !                Data Assimilation with an Ensemble Kalman Filter
  !                <https://journals.ametsoc.org/view/journals/mwre/132/5/1520-0493_2004_132_1238_ioieao_2.0.co_2.xml>`_
  !
  !              - **RTPS**: Relaxation to the prior spread will relax the analysis ensemble spread (i.e. standard
  !                deviation) to the trial ensemble spread with a given amount (:math:`\alpha`) in the namelist.
  !                With the ensemble spread of analysis (:math:`\sigma^a`)
  !                and trial (:math:`\sigma^b`), it can be formulated as:
  !                :math:`x^a_i \leftarrow \overline{x^a} + (\overline{x^a}-x^a_i)(1+\alpha(\sigma^b-\sigma^a)/\sigma^a)`
  !
  !                More details for the RTPS can be found in paper:
  !                `Evaluating Methods to Account for System Errors in Ensemble Data Assimilation
  !                <https://journals.ametsoc.org/view/journals/mwre/140/9/mwr-d-11-00276.1.xml>`_
  !
  !              - **Humidity limits**: Impose the saturation limits and RTTOV limits
  !                on the humidity for each member.
  !
  !              - **Recenter**: The ensemble analyses from the LETKF are recentered
  !                around the EnVar analysis based on the input coefficent.
  !                It is also called hybrid gain algorithm and more details can be found in
  !                paper:
  !                `Using the hybrid gain algorithm to sample data assimilation uncertainty
  !                <https://rmets.onlinelibrary.wiley.com/doi/10.1002/qj.3426>`_
  !
  !              - **Subsample**: Select 20 ensemble members for the medium range forecasts.
  !
  !              - **Random perturbation**: Generate homogeneous and isotropic random perturbations
  !                and add them to the ensemble analyses.
  !
  !              - **Analysis increments**: After all the processes above, the analysis incrememnts are
  !                computed at the central time for all the members including the control member and
  !                the subsampled 20 members.
  !
  !              - **Mask analysis increments**: For LAM grid, the analysis increments in the blending
  !                area are masked out.
  !
  !            --
  !
  !============================================= ==============================================================
  ! Input and Output Files                        Description of file
  !============================================= ==============================================================
  ! ``flnml``                                     In - Main namelist file with parameters user may modify
  ! ``ensemble_trial/$YYYYMMDDHH_006_$NNNN``      In - Background ensemble member files
  ! ``ensemble_anal/$YYYYMMDDHH_000_$NNNN``       In - Analysis ensemble member files
  ! ``bgcov``                                     In - Bnmc matrix for the random perturbations
  ! ``analysis_grid``                             In - Horizontal grid file on which the random perturbations
  !                                               are generated
  ! ``$YYYYMMDDHH_recentering_analysis``          In - Analysis from EnVar for recentering the ensemble analyses
  ! ``analinc_mask``                              In - Mask file for masking the analysis increments for LAM grid
  ! ``$YYYYMMDDHH_000_inc_$NNNN``                 Out - Ensemble analysis increments for trials
  ! ``$YYYYMMDDHH_000_$NNNN``                     Out - Ensemble analyses recentered and perturbed
  ! ``$YYYYMMDDHH_000_$NNNN_raw``                 Out - Ensemble analyses from LETKF
  ! ``$YYYYMMDDHH_000_analmean``                  Out - Mean ensemble analysis without perturbation
  ! ``$YYYYMMDDHH_000_analrms``                   Out - Analysis ensemble spread (i.e. standard deviation) without perturbation
  ! ``$YYYYMMDDHH_000_analrms_ascii``             Out - Global average of  ``$YYYYMMDDHH_000_analrms``
  ! ``$YYYYMMDDHH_000_analpertmean``              Out - Mean ensemble analysis with perturbation and recentering
  ! ``$YYYYMMDDHH_000_analpertrms``               Out - Analysis ensemble spread with perturbation and recentering
  ! ``$YYYYMMDDHH_000_analpertrms_ascii``         Out - Global average of  ``$YYYYMMDDHH_000_analpertrms``
  ! ``$YYYYMMDDHH_006_trialmean``                 Out - Mean ensemble trial for all time levels
  ! ``$YYYYMMDDHH_006_trialrms``                  Out - Trial ensemble spread at 6H
  ! ``$YYYYMMDDHH_006_trialrms_ascii``            Out - Global average of ``$YYYYMMDDHH_006_trialrms``
  ! ``subspace/$YYYYMMDDHH_000_inc_$NNNN``        Out - Subsampled ensemble analysis increments for progs
  ! ``subspace/$YYYYMMDDHH_000_$NNNN``            Out - Subsampled ensemble analyses recentered and perturbed
  ! ``subspace_unpert/$YYYYMMDDHH_000_$NNNN``     Out - Unperturbed subsampled ensemble analyses recentered
  !============================================= ==============================================================
  !
  !            --
  !
  !:Synopsis: Below is a summary of the ``ensPostProcess`` program calling sequence:
  !
  !             - **Initial setups:**
  !
  !               - Read the NAMLETKF namelist and check/modify some values.
  !
  !               - Various modules are setup: ``gridStateVector_mod``, ``timeCoord_mod`` (and set up dates
  !                 and ``dateStampList`` variables for both trials and increments/analyses).
  !
  !               - Setup horizontal and vertical grid objects from first
  !                 ensemble member file.
  !
  !               - Allocate varaibles and read ensemble analyses and trials
  !
  !             - **Ensemble postprocess:**
  !
  !               - Compute ensemble mean and spread (i.e. standard deviation) for analysis and trial
  !
  !               - Perform RTPP and RTPS and recompute analysis spread
  !
  !               - Recenter ensemble mean and recompute analysis mean and spread
  !
  !               - Apply humidity saturation limit and RTTOV HU limit and recompute analysis mean
  !
  !               - Subsample ensemble analysis
  !
  !               - Add random perturbation to the analyses and compute perturbed analysis spread
  !
  !               - Add random perturbation to the subsampled analyses and recenter subsampled mean to the
  !                 global mean.
  !
  !               - Compute analysis increments for all analyses and subsampled analyses including the ensemble mean
  !
  !               - Mask LAM analysis increments are recompute the ensemble analyses by adding the increments
  !                 to the trials
  !
  !               - Writing the outputs.
  !
  !:Options: `List of namelist blocks <../namelists_in_each_program.html#ensPostProcess>`_
  !          that can affect the ``ensPostProcess`` program.
  !
  !========================== ========================== ==============================================================
  ! Program/Module             Namelist                   Description of what is controlled
  !========================== ========================== ==============================================================
  ! ``midas_ensPostProcess``   ``NAMENSPOSTPROC``         Number of ensemble members and read, write control of trials
  !                                                       and analyses and horizontal interpolation degree
  ! ``enspostprocess_mod``     ``NAMENSPOSTPROCMODULE``   Control the postprocess algorithms
  ! ``timeCoord_mod``          ``NAMTIME``                Temporal resolution of the background state and the analysis
  !========================== ========================== ==============================================================
  !

  use version_mod
  use midasMpi_mod
  use fileNames_mod
  use ensembleStateVector_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use timeCoord_mod
  use utilities_mod
  use ramDisk_mod
  use ensPostProcess_mod
  implicit none

  type(struct_ens)          :: ensembleTrl
  type(struct_ens)          :: ensembleAnl
  type(struct_vco), pointer :: vco_ens => null()
  type(struct_hco), pointer :: hco_ens => null()
  type(struct_gsv)          :: stateVectorHeightSfc, stateVectorCtrlTrl

  character(len=256) :: ensPathNameAnl = 'ensemble_anal'
  character(len=256) :: ensPathNameTrl = 'ensemble_trial'
  character(len=256) :: ensFileName, gridFileName, ctrlFileName
  integer, allocatable :: dateStampList(:)
  integer :: ierr, nulnam, stepIndex
  logical :: targetGridFileExists
  integer, external :: get_max_rss, fstopc, fnom, fclos

  integer :: nEns             ! ensemble size
  logical :: readTrlEnsemble  ! activate reading of trial ensemble
  logical :: readAnlEnsemble  ! activate reading of analysis ensemble
  logical :: writeTrlEnsemble ! activate writing of the trial ensemble (useful when it's interpolated)
  character(len=12) :: hInterpolationDegree ! select degree of horizontal interpolation (if needed)
  NAMELIST /namEnsPostProc/nEns, readTrlEnsemble, readAnlEnsemble, &
                           writeTrlEnsemble, hInterpolationDegree

  call ver_printNameAndVersion('ensPostProcess','Program for post-processing of LETKF analysis ensemble')

  !
  !- 0. MPI, TMG and misc. initialization
  !
  call mmpi_initialize
  call tmg_init(mmpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Avoid printing lots of stuff to listing for std file I/O
  ierr = fstopc('MSGLVL','ERRORS',0)

  ! Setup the ramdisk directory (if supplied)
  call ram_setup

  !- Setting default namelist variable values
  nEns = 256
  readTrlEnsemble  = .true.
  readAnlEnsemble  = .true.
  writeTrlEnsemble = .false.
  hInterpolationDegree = 'LINEAR' ! or 'CUBIC' or 'NEAREST'

  !- Read the namelist
  nulnam = 0
  ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namEnsPostProc, iostat=ierr)
  if ( ierr /= 0) call utl_abort('midas-ensPostProcess: Error reading namelist')
  if ( mmpi_myid == 0 ) write(*,nml=namEnsPostProc)
  ierr = fclos(nulnam)

  if (.not.readTrlEnsemble .and. .not.readAnlEnsemble) then
    call utl_abort('midas-ensPostProcess: must read either Trial or Analysis ensemble')
  end if

  if (writeTrlEnsemble .and. .not.readTrlEnsemble) then
    call utl_abort('midas-ensPostProcess: cannot write Trial ensemble if it is not read')
  end if

  if (writeTrlEnsemble .and. readAnlEnsemble) then
    call utl_abort('midas-ensPostProcess: cannot write Trial ensemble when Analysis ensemble is read')
  end if

  !- 1. Initialize date/time-related info

  ! Setup timeCoord module, set dateStamp with value from trial or analysis ensemble member 1
  if (readTrlEnsemble) then
    call fln_ensFileName(ensFileName, ensPathNameTrl, memberIndex_opt=1, &
                         copyToRamDisk_opt=.false.)
  else
    call fln_ensFileName(ensFileName, ensPathNameAnl, memberIndex_opt=1, &
                         copyToRamDisk_opt=.false.)
  end if
  call tim_setup(fileNameForDate_opt=ensFileName)
  allocate(dateStampList(tim_nstepobsinc))
  call tim_getstamplist(dateStampList,tim_nstepobsinc,tim_getDatestamp())

  write(*,*) 'midas-ensPostProcess: analysis dateStamp = ', tim_getDatestamp()

  !- 2. Initialize variables and grids

  !- Initialize variables of the model states
  call gsv_setup

  !- Initialize the grid from targetgrid file or from trial or analysis ensemble member 1
  if (mmpi_myid == 0) write(*,*) ''
  if (mmpi_myid == 0) write(*,*) 'midas-ensPostProcess: Set hco and vco parameters for ensemble grid'
  inquire(file='targetgrid', exist=targetGridFileExists)
  if (targetGridFileExists) then
    gridFileName = 'targetgrid'
  else if (readTrlEnsemble) then
    call fln_ensFileName(gridFileName, ensPathNameTrl, memberIndex_opt=1, &
                         copyToRamDisk_opt=.false.)
  else
    call fln_ensFileName(gridFileName, ensPathNameAnl, memberIndex_opt=1, &
                         copyToRamDisk_opt=.false.)
  end if
  if (mmpi_myid == 0) then
    write(*,*) 'midas-ensPostProcess: file use to define grid = ', trim(gridFileName)
  end if
  call hco_SetupFromFile(hco_ens, gridFileName, ' ', 'ENSFILEGRID')
  call vco_setupFromFile(vco_ens, gridFileName)

  !- Read the sfc height from trial ensemble member 1 - only if we are doing NWP!
  if (vco_getNumLev(vco_ens,'TH') > 0 .or. vco_getNumLev(vco_ens,'MM') > 0) then
    if (readTrlEnsemble) then
      call fln_ensFileName(ensFileName, ensPathNameTrl, memberIndex_opt=1, &
                           copyToRamDisk_opt=.false.)
    else
      call fln_ensFileName(ensFileName, ensPathNameAnl, memberIndex_opt=1, &
                           copyToRamDisk_opt=.false.)
    end if
    call gsv_allocate(stateVectorHeightSfc, 1, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                      mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                      hInterpolateDegree_opt=hInterpolationDegree, &
                      dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0','TT'/))
    call gio_readFromFile(stateVectorHeightSfc, ensFileName, ' ', ' ',  &
                          containsFullField_opt=.true., readHeightSfc_opt=.true.)
  end if

  !- 3. Allocate and read ensembles

  call utl_tmg_start(2,'--ReadEnsemble')

  !- Allocate ensembles, read the Anl ensemble
  if (readAnlEnsemble) then
    call fln_ensFileName(ensFileName, ensPathNameAnl, resetFileInfo_opt=.true.)
    call ens_allocate(ensembleAnl, nEns, tim_nstepobsinc, hco_ens, vco_ens, &
                      dateStampList, hInterpolateDegree_opt=hInterpolationDegree)
    call ens_readEnsemble(ensembleAnl, ensPathNameAnl, biPeriodic=.false.)
  end if

  !- Allocate ensembles, read the Trl ensemble
  if (readTrlEnsemble) then
    call fln_ensFileName(ensFileName, ensPathNameAnl, resetFileInfo_opt=.true.)
    call ens_allocate(ensembleTrl, nEns, tim_nstepobsinc, hco_ens, vco_ens, &
                      dateStampList, hInterpolateDegree_opt=hInterpolationDegree)
    call ens_readEnsemble(ensembleTrl, ensPathNameTrl, biPeriodic=.false.)

  end if

  call utl_tmg_stop(2)

  !- Allocate and read the Trl control member
  if (readTrlEnsemble .and. readAnlEnsemble) then
    !- Allocate and read control member Trl
    call gsv_allocate( stateVectorCtrlTrl, tim_nstepobsinc, hco_ens, vco_ens, &
                       dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., &
                       hInterpolateDegree_opt=hInterpolationDegree, &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call fln_ensFileName(ctrlFileName, ensPathNameTrl, memberIndex_opt=0, &
                         copyToRamDisk_opt=.false.)
    do stepIndex = 1, tim_nstepobsinc
      call gio_readFromFile( stateVectorCtrlTrl, ctrlFileName, ' ', ' ',  &
                             stepIndex_opt=stepIndex, containsFullField_opt=.true., &
                             readHeightSfc_opt=.false. )
    end do
  end if

  !- 4. Post processing of the analysis results (if desired) and write everything to files
  call epp_postProcess(ensembleTrl, ensembleAnl, &
                       stateVectorHeightSfc, stateVectorCtrlTrl, &
                       writeTrlEnsemble)

  !
  !- 5. MPI, tmg finalize
  !  
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  call utl_tmg_stop(0)

  call tmg_terminate(mmpi_myid, 'TMG_INFO')
  call rpn_comm_finalize(ierr) 

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  if ( mmpi_myid == 0 ) write(*,*) ' --------------------------------'
  if ( mmpi_myid == 0 ) write(*,*) ' MIDAS-ensPostProcess ENDS'
  if ( mmpi_myid == 0 ) write(*,*) ' --------------------------------'

end program midas_ensPostProcess

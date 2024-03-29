program midas_letkf
  !
  !:Purpose: Main program for the local ensemble transform Kalman filter
  !          (LETKF). Several different variations of the LETKF algorithm have
  !          been implemented. Note that the actual calculation of the analyses
  !          is in the subroutine ``enkf_LETKFanalyses``. Many aspects of this
  !          program are controlled through the namelist block NAMLETKF.
  !
  !          ---
  !
  !:Algorithm: The ``letkf`` programs implements several variations of the LETKF
  !            ensemble data assimilation algorithm. The following variations
  !            can be chosen through the namelist:
  !
  !              - **LETKF**: standard LETKF
  !
  !              - **CVLETKF**: LETKF with deterministic approach to cross validation using
  !                the "gain form" (currently operational)
  !
  !              - **CVLETKF-PERTOBS**: LETKF with stochastic approach to cross validation
  !
  !              - **LETKF-Gain**: standard LETKF, but using the "gain form"
  !
  !              - **LETKF-Gain-ME**: standard LETKF, but using the "modulated ensemble" and "gain
  !                form" to account for vertical localization
  !
  !              - **CVLETKF-ME**: LETKF with deterministic approach to cross validation using
  !                the "modulated ensemble" and "gain form" to account for
  !                vertical localization
  !
  !            More detail on the LETKF algorithms with and without cross
  !            validation (deterministic and stochastic) can be found in the
  !            paper:
  !            `LETKF with cross validation <https://journals.ametsoc.org/view/journals/mwre/148/6/MWR-D-19-0402.1.xml>`_
  !
  !            --
  !
  !============================================== ==============================================================
  ! Input and Output Files (NWP applicaton)        Description of file
  !============================================== ==============================================================
  ! ``flnml``                                      In - Main namelist file with parameters user may modify
  ! ``flnml_static``                               In - The "static" namelist that should not be modified
  ! ``ensemble/$YYYYMMDDHH_006_$NNNN``             In - Background ensemble member files
  ! ``obsfiles_$FAM/obs$FAM_$NNNN_$NNNN``          In - Observation file for each "family" and MPI task
  ! ``obserr``                                     In - Observation error statistics
  ! ``$YYYYMMDDHH_000_$NNNN``                      Out - Analysis ensemble member files
  ! ``obsfiles_$FAM.updated/obs$FAM_$NNNN_$NNNN``  Out - Updated obs file for each "family" and MPI task
  ! Remainder are files related to radiance obs:
  ! ``stats_$SENSOR_assim``                        In - Satellite radiance observation errors of difference sensors
  ! ``stats_tovs``                                 In - Satellite radiance observation errors
  ! ``stats_tovs_symmetricObsErr``                 In - User-defined symmetric TOVS errors for all sky
  ! ``Cmat_$PLATFORM_$SENSOR.dat``                 In - Inter-channel observation-error correlations
  ! ``ceres_global.std``                           In - High-res surface type and water fraction for radiance obs
  ! ``rtcoef_$PLATFORM_$SENSOR.dat``               In - RTTOV coefficient files
  ! ``rttov_h2o_limits.dat``                       In - Min/max humidity limits applied to analysis
  ! ``ozoneclim98``                                In - Ozone climatology
  !============================================== ==============================================================
  !
  !           --
  !
  !:Synopsis: Below is a summary of the ``letkf`` program calling sequence:
  !
  !             - **Initial setups:**
  !
  !               - Read the NAMLETKF namelist and check/modify some values.
  !
  !               - Various modules are setup: ``obsFiles_mod``,
  !                 ``gridStateVector_mod``, ``timeCoord_mod`` (and set up dates
  !                 and ``dateStampList`` variables for both trials and
  !                 increments/analyses).
  !
  !               - Setup horizontal and vertical grid objects from first
  !                 ensemble member file and determine if this is an NWP or
  !                 ocean application.
  !
  !               - Setup ``obsSpaceData`` object and read observations from
  !                 files: ``inn_setupObs``.
  !
  !               - Setup the observation error statistics in ``obsSpaceData``
  !                 object: ``oer_setObsErrors``.
  !
  !               - Allocate and some setup of objects for
  !                 ``ensembleObservations_mod``.
  !
  !               - Allocate objects for ``column_mod`` and
  !                 ``gridStateVector_mod``.
  !
  !               - Allocate ensemble object and read trial ensemble:
  !                 ``ens_readEnsemble``.
  !
  !               - Optionally, read a deterministic state for recentering the
  !                 trial ensemble before the analysis (for special applications
  !                 using LETKF for deterministic analysis).
  !
  !               - Compute ensemble mean: ``ens_computeMean``.
  !
  !             - **LETKF computations:**
  !
  !               - Loop over trial members, computing innovation for each, with
  !                 resulting :math:`H(xb)` being stored in ``ensObs`` objects both for
  !                 original ensemble and, optionally, for the modulated
  !                 ensemble members.
  !
  !               - Compute some additional quantities in ``ensObs`` objects based
  !                 on :math:`H(xb)` values.
  !
  !               - Compute :math:`y-H(xb)` for trial ensemble mean.
  !
  !               - Additional observation quality control procedures based on
  !                 quantities computed from trial ensemble.
  !
  !               - Gather ``ensObs`` quantities from all MPI tasks onto all MPI
  !                 tasks.
  !
  !               - Allocate and prepare objects for ``ensembleStateVector_mod``
  !                 used to store the trial and analysis ensembles with
  !                 temporal resolution of analysis (can be same as temporal
  !                 resolution for innovation calculation or only a single time
  !                 step).
  !
  !               - Setup information for interpolating weights from coarse grid
  !                 to the full model grid: ``enkf_setupInterpInfo``.
  !
  !               - Perform LETKF analysis to compute the analysis ensemble:
  !                 ``enkf_LETKFanalyses``
  !
  !             - **Final steps:**
  !
  !               - Optionally, compute :math:`H(xa)` for each member and put values in
  !                 ``ensObs`` object for later output to diag SQLite files.
  !
  !               - Write out (update) the observation files: ``obsf_writeFiles``.
  !
  !               - Compute :math:`y-H(xa)` for ensemble mean analysis, stored in
  !                 ``obsSpaceData``.
  !
  !               - Optionally, do post-processing of the analysis ensemble (same
  !                 functionality as the ``ensPostProcess`` program), or just
  !                 write out the raw analysis ensemble (for later processing by
  !                 ``ensPostProcess``).
  !
  !           --
  !
  !:Options: `List of namelist blocks <../namelists_in_each_program.html#letkf>`_
  !          that can affect the ``letkf`` program.
  !
  !          * Some of the other relevant namelist blocks used to configure the
  !            letkf analysis are listed in the following table:
  ! 
  !======================== ============== ==============================================================
  ! Program/Module           Namelist       Description of what is controlled
  !======================== ============== ==============================================================
  ! ``midas_letkf``          ``NAMLETKF``   LETKF algorithm, number of ensemble members and
  !                                         additional parameters for controlling the LETKF analysis
  ! ``timeCoord_mod``        ``NAMTIME``    assimilation time window length, temporal resolution of
  !                                         the background state and the analysis
  !======================== ============== ==============================================================
  !
  use version_mod
  use midasMpi_mod
  use mathPhysConstants_mod
  use fileNames_mod
  use ensembleObservations_mod
  use ensembleStateVector_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use columnData_mod
  use tovsNL_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use oceanMask_mod
  use timeCoord_mod
  use obsTimeInterp_mod
  use utilities_mod
  use ramDisk_mod
  use stateToColumn_mod
  use obsFiles_mod
  use obsSpaceData_mod
  use obsErrors_mod
  use obsFilter_mod
  use innovation_mod
  use enkf_mod
  use ensPostProcess_mod
  implicit none

  type(struct_obs), target  :: obsSpaceData
  type(struct_ens), target  :: ensembleTrl4D
  type(struct_ens), pointer :: ensembleTrl
  type(struct_ens)          :: ensembleAnl
  type(struct_gsv)          :: stateVectorMeanTrl4D
  type(struct_gsv)          :: stateVectorMemberAnl
  type(struct_gsv)          :: stateVectorMeanAnl
  type(struct_gsv)          :: stateVector4D
  type(struct_gsv)          :: stateVector4Dmod
  type(struct_gsv)          :: stateVectorWithZandP4D
  type(struct_gsv)          :: stateVectorHeightSfc
  type(struct_gsv)          :: stateVectorCtrlTrl
  type(struct_gsv)          :: stateVectorRecenter
  type(struct_columnData)   :: column

  type(struct_eob), target  :: ensObs, ensObs_mpiglobal
  type(struct_eob), pointer :: ensObsGain, ensObsGain_mpiglobal

  type(struct_vco), pointer :: vco_ens => null()
  type(struct_hco), pointer :: hco_ens => null()
  type(struct_ocm)          :: oceanMask

  integer :: memberIndex, middleStepIndex, stepIndex, randomSeedObs
  integer :: nulnam, dateStampFromObs, ierr
  integer :: get_max_rss, fclos, fnom, fstopc
  integer :: nEnsGain, eigenVectorIndex, memberIndexInEnsObs
  integer, allocatable :: dateStampList(:), dateStampListInc(:)

  character(len=256) :: ensFileName, ctrlFileName, recenterFileName
  character(len=9)   :: obsColumnMode
  character(len=48)  :: obsMpiStrategy
  character(len=48)  :: midasMode

  logical :: nwpFields   ! indicates if fields are on momentum and thermo levels
  logical :: oceanFields ! indicates if fields are on depth levels
  logical :: useModulatedEns ! using modulated ensembles is requested by setting numRetainedEigen

  ! interpolation information for weights (in enkf_mod)
  type(struct_enkfInterpInfo) :: wInterpInfo

  real(4), pointer :: field_Psfc(:,:,:,:)

  ! namelist variables
  character(len=20)  :: algorithm  ! name of the chosen LETKF algorithm: 'LETKF', 'CVLETKF'
  logical            :: ensPostProcessing ! do all post-processing of analysis ensemble
  logical            :: recenterInputEns  ! read a deterministic state to recenter ensemble
  integer            :: numSubEns  ! number of sub-ensembles to split the full ensemble
  character(len=256) :: ensPathName ! absolute or relative path to ensemble directory
  integer  :: nEns                 ! ensemble size
  logical  :: randomShuffleSubEns  ! choose to randomly shuffle members into subensembles 
  logical  :: writeLocalEnsObsToFile ! Controls writing the ensObs to file.
  integer  :: maxNumLocalObs       ! maximum number of obs in each local volume to assimilate
  integer  :: weightLatLonStep     ! separation of lat-lon grid points for weight calculation
  integer  :: numRetainedEigen     ! number of retained eigenValues/Vectors of vertical localization matrix
  logical  :: modifyAmsubObsError  ! reduce AMSU-B obs error stddev in tropics
  logical  :: backgroundCheck      ! apply additional background check using ensemble spread
  logical  :: huberize             ! apply huber norm quality control procedure
  logical  :: rejectHighLatIR      ! reject all IR observations at high latitudes
  logical  :: rejectRadNearSfc     ! reject radiance observations near the surface
  logical  :: ignoreEnsDate        ! when reading ensemble, ignore the date
  logical  :: outputOnlyEnsMean    ! when writing ensemble, can choose to only write member zero
  logical  :: outputEnsObs         ! to write trial and analysis ensemble members in observation space to sqlite 
  logical  :: debug                ! debug option to print values to the listings.
  logical  :: readEnsObsFromFile   ! instead of computing innovations, read ensObs%Yb from file.
  real(8)  :: hLocalize(4)         ! horizontal localization radius (in km)
  real(8)  :: hLocalizePressure(3) ! pressures where horizontal localization changes (in hPa)
  real(8)  :: vLocalize            ! vertical localization radius (units: ln(Pressure in Pa) or meters)
  real(8)  :: minDistanceToLand    ! for ice/ocean DA: minimum distance to land for assimilating obs
  character(len=20) :: obsTimeInterpType ! type of time interpolation to obs time
  character(len=20) :: mpiDistribution   ! type of mpiDistribution for weight calculation ('ROUNDROBIN' or 'TILES')
  character(len=12) :: etiket_anl        ! etiket for output files
  
  NAMELIST /NAMLETKF/algorithm, ensPostProcessing, recenterInputEns, nEns, numSubEns, &
                     ensPathName, randomShuffleSubEns,  &
                     hLocalize, hLocalizePressure, vLocalize, minDistanceToLand,  &
                     maxNumLocalObs, weightLatLonStep,  &
                     modifyAmsubObsError, backgroundCheck, huberize, rejectHighLatIR, rejectRadNearSfc,  &
                     ignoreEnsDate, outputOnlyEnsMean, outputEnsObs,  & 
                     obsTimeInterpType, mpiDistribution, etiket_anl, &
                     readEnsObsFromFile, writeLocalEnsObsToFile, &
                     numRetainedEigen, debug

  ! Some high-level configuration settings
  midasMode = 'analysis'
  obsColumnMode = 'ENKFMIDAS'
  obsMpiStrategy = 'LIKESPLITFILES'

  call ver_printNameAndVersion('letkf','Program for Local Ensemble Transform Kalman Filter')

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

  !
  !- 1. Set/Read values for the namelist NAMLETKF
  !

  !- 1.1 Setting default namelist variable values
  algorithm                = 'LETKF'
  ensPostProcessing        = .false.
  recenterInputEns         = .false.
  ensPathName              = 'ensemble'
  nEns                     = 10
  numSubEns                = 2
  randomShuffleSubEns      = .false.
  maxNumLocalObs           = 1000
  weightLatLonStep         = 1
  modifyAmsubObsError      = .false.
  backgroundCheck          = .false.
  huberize                 = .false.
  rejectHighLatIR          = .false.
  rejectRadNearSfc         = .false.
  ignoreEnsDate            = .false.
  outputOnlyEnsMean        = .false.
  outputEnsObs             = .false.
  hLocalize(:)             = -1.0D0
  hLocalizePressure        = (/14.0D0, 140.0D0, 400.0D0/)
  vLocalize                = -1.0D0
  minDistanceToLand        = -1.0D0
  obsTimeInterpType        = 'LINEAR'
  mpiDistribution          = 'ROUNDROBIN'
  etiket_anl               = 'ENS_ANL'
  readEnsObsFromFile       = .false.
  writeLocalEnsObsToFile   = .false.
  numRetainedEigen         = 0
  debug                    = .false.
  
  !- 1.2 Read the namelist
  nulnam = 0
  ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namletkf, iostat=ierr)
  if ( ierr /= 0) call utl_abort('midas-letkf: Error reading namelist')
  if ( mmpi_myid == 0 ) write(*,nml=namletkf)
  ierr = fclos(nulnam)

  !- 1.3 Some minor modifications of namelist values
  if (hLocalize(1) > 0.0D0 .and. hLocalize(2) < 0.0D0) then
    ! if only 1 value given for hLocalize, use it for entire column
    hLocalize(2:4) = hLocalize(1)
    if ( mmpi_myid == 0 ) write(*,*) 'midas-letkf: hLocalize(2:4) are modified after reading namelist. ' // &
                                    'hLocalize(2:4)=', hLocalize(1)
  else if ( hLocalize(1) < 0.0D0 ) then
    call utl_abort('midas-letkf: hLocalize(1) < 0.0D0')
  end if
  hLocalize(:) = hLocalize(:) * 1000.0D0 ! convert from km to m
  hLocalizePressure(:) = log(hLocalizePressure(:) * MPC_PA_PER_MBAR_R8)

  if (minDistanceToLand > 0.0D0) then
    minDistanceToLand = minDistanceToLand * 1000.0D0 ! convert from km to m
  end if

  if ( trim(algorithm) /= 'LETKF'           .and. &
       trim(algorithm) /= 'CVLETKF'         .and. &
       trim(algorithm) /= 'CVLETKF-PERTOBS' .and. &
       trim(algorithm) /= 'LETKF-Gain'      .and. &
       trim(algorithm) /= 'LETKF-Gain-ME'   .and. &
       trim(algorithm) /= 'CVLETKF-ME' ) then
    call utl_abort('midas-letkf: unknown LETKF algorithm: ' // trim(algorithm))
  end if

  if ( numRetainedEigen < 0 ) call utl_abort('midas-letkf: numRetainedEigen should be ' // &
    'equal or greater than zero')

  useModulatedEns = (numRetainedEigen > 0)

  if ( trim(algorithm) == 'LETKF-Gain-ME' .or. trim(algorithm) == 'CVLETKF-ME' ) then
    if ( .not. useModulatedEns ) call utl_abort('midas-letkf: numRetainedEigen should be ' // &
    'equal or greater than one for LETKF algorithm: ' // &
    trim(algorithm))
  else
    if ( useModulatedEns ) call utl_abort('midas-letkf: numRetainedEigen should be ' // &
    'equal to zero for LETKF algorithm: ' // &
    trim(algorithm))
  end if
  
  ! check for NO varying horizontal localization lengthscale in letkf with modulated ensembles.
  if ( .not. all(hLocalize(2:4) == hLocalize(1)) .and. useModulatedEns ) then
    call utl_abort('midas-letkf: Varying horizontal localization lengthscales is NOT allowed in ' // &
    'letkf with modulated ensembles')
  end if

  !
  !- 2.  Initialization
  !

  !- 2.1 Setup the observation file names and get dateStamp from obs
  call obsf_setup(dateStampFromObs, midasMode)

  !- 2.2 Initialize date/time-related info

  ! Setup timeCoord module, set dateStamp from env variable
  call tim_setup()
  if (tim_getDateStamp() == 0) then
    if (dateStampFromObs > 0) then
      ! use dateStamp from obs if not already set by env variable
      call tim_setDateStamp(dateStampFromObs)
    else
      call utl_abort('midas-letkf: DateStamp was not set')
    end if
  end if
  if (tim_nstepobsinc /= 1 .and. tim_nstepobsinc /= tim_nstepobs) then
    call utl_abort('midas-letkf: invalid value for namelist variable DSTEPOBSINC. ' // &
                   'Increments can be either 3D or have same number of time steps as trials')
  end if
  allocate(dateStampList(tim_nstepobs))
  call tim_getstamplist(dateStampList,tim_nstepobs,tim_getDatestamp())
  allocate(dateStampListInc(tim_nstepobsinc))
  call tim_getstamplist(dateStampListInc,tim_nstepobsinc,tim_getDatestamp())

  write(*,*) 'midas-letkf: analysis dateStamp = ',tim_getDatestamp()

  !- 2.3 Initialize variables of the model states
  call gsv_setup

  !- 2.4 Initialize the Ensemble grid
  if (mmpi_myid == 0) write(*,*) ''
  if (mmpi_myid == 0) write(*,*) 'midas-letkf: Set hco and vco parameters for ensemble grid'
  call fln_ensFileName( ensFileName, ensPathName, memberIndex_opt=1, &
                        copyToRamDisk_opt=.false. )
  call hco_SetupFromFile( hco_ens, ensFileName, ' ', 'ENSFILEGRID')
  call vco_setupFromFile( vco_ens, ensFileName )
  if (vco_getNumLev(vco_ens, 'MM') /= vco_getNumLev(vco_ens, 'TH')) then
    call utl_abort('midas-letkf: nLev_M /= nLev_T - currently not supported')
  end if
  nwpFields   = (vco_getNumLev(vco_ens,'TH') > 0 .or. vco_getNumLev(vco_ens,'MM') > 0)
  oceanFields = (vco_getNumLev(vco_ens,'DP') > 0)
  if (.not.nwpFields .and. .not.oceanFields) then
    call utl_abort('midas-letkf: vertical coordinate does not contain nwp nor ocean fields')
  end if

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !- 2.5 Read in the observations and other obs-related set up

  ! Read the observations
  call inn_setupObs( obsSpaceData, hco_ens, obsColumnMode, obsMpiStrategy, midasMode,  &
                     obsClean_opt = .false. )

  ! Initialize obs error covariances and set flag using 'util' column of stats_tovs
  call oer_setObsErrors(obsSpaceData, midasMode, useTovsUtil_opt=.true.) ! IN

  ! Call suprep again to filter out channels according to 'util' column of stats_tovs
  call filt_suprep(obsSpaceData)

  ! Allocate vectors for storing HX values
  call eob_allocate(ensObs, nEns, obs_numBody(obsSpaceData), obsSpaceData)
  if ( outputEnsObs ) allocate(ensObs%Ya_r4(ensObs%numMembers,ensObs%numObs))
  call eob_zero(ensObs)
  if ( useModulatedEns ) then
    nEnsGain = nEns * numRetainedEigen
    allocate(ensObsGain)
    call eob_allocate(ensObsGain, nEnsGain, obs_numBody(obsSpaceData), obsSpaceData)
    call eob_zero(ensObsGain)
  else
    ensObsGain => ensObs
  end if

  ! Set lat, lon, obs values in ensObs
  call eob_setLatLonObs(ensObs)
  if ( useModulatedEns ) call eob_setLatLonObs(ensObsGain)

  !- 2.6 Initialize a single columnData object
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  call col_setup
  call col_setVco(column, vco_ens)
  call col_allocate(column, obs_numheader(obsSpaceData),  &
                    setToZero_opt=.true.)
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !- 2.7 Read the sfc height from ensemble member 1 - only if we are doing NWP
  if ( nwpFields ) then
    call gsv_allocate( stateVectorHeightSfc, 1, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0','TT'/) )
    call gio_readFromFile( stateVectorHeightSfc, ensFileName, ' ', ' ',  &
                           containsFullField_opt=.true., readHeightSfc_opt=.true. )
  end if

  !- 2.8 Allocate statevector related to ensemble mean
  call gsv_allocate( stateVectorMeanTrl4D, tim_nstepobs, hco_ens, vco_ens, &
                     dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeightSfc_opt=.true., &
                     allocHeight_opt=.false., allocPressure_opt=.false. )
  call gsv_zero(stateVectorMeanTrl4D)
  call gsv_allocate( stateVectorMeanAnl, tim_nstepobsinc, hco_ens, vco_ens, &
                     dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeightSfc_opt=.true., &
                     allocHeight_opt=.false., allocPressure_opt=.false. )
  call gsv_zero(stateVectorMeanAnl)
  call gsv_allocate( stateVector4D, tim_nstepobs, hco_ens, vco_ens, &
                     dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeightSfc_opt=.true., &
                     allocHeight_opt=.false., allocPressure_opt=.false. )
  call gsv_zero(stateVector4D)
  if ( useModulatedEns ) then
    ! same as stateVector4D
    call gsv_allocate( stateVector4Dmod, tim_nstepobs, hco_ens, vco_ens, &
                       dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_zero(stateVector4Dmod)
  end if

  !- 2.9 Allocate statevector related to an analysis ensemble member  
  call gsv_allocate( stateVectorMemberAnl, tim_nstepobsinc, hco_ens, vco_ens, &
                     dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeightSfc_opt=.true. )
  call gsv_zero(stateVectorMemberAnl)
  
  !- 2.10 Allocate statevector for storing state with heights and pressures allocated (for s2c_nl)
  call gsv_allocate( stateVectorWithZandP4D, tim_nstepobs, hco_ens, vco_ens, &
                     dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeightSfc_opt=.true. )
  call gsv_zero(stateVectorWithZandP4D)

  !- 2.11 Allocate ensembles, read the Trl ensemble
  call utl_tmg_start(2,'--ReadEnsemble')
  call ens_allocate(ensembleTrl4D, nEns, tim_nstepobs, hco_ens, vco_ens, dateStampList)
  call ens_readEnsemble(ensembleTrl4D, ensPathName, biPeriodic=.false., &
                        ignoreDate_opt=ignoreEnsDate)
  call utl_tmg_stop(2)

  !- 2.12 If desired, read a deterministic state for recentering the ensemble
  if (recenterInputEns) then
    call gsv_allocate( stateVectorRecenter, tim_nstepobs, hco_ens, vco_ens, &
                       dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.false., &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_zero(stateVectorRecenter)
    call fln_ensTrlFileName( recenterFileName, './', tim_getDateStamp() )
    do stepIndex = 1, tim_nstepobs
      call gio_readFromFile( stateVectorRecenter, recenterFileName, ' ', ' ',  &
                             stepIndex_opt=stepIndex, containsFullField_opt=.true., &
                             readHeightSfc_opt=.false. )
    end do
    call ens_recenter( ensembleTrl4D, stateVectorRecenter, recenteringCoeffScalar_opt=1.0d0 )
    call gsv_deallocate( stateVectorRecenter )
  end if
  
  !- 2.13 Compute ensemble mean and copy to meanTrl and meanAnl stateVectors
  call ens_computeMean(ensembleTrl4D)
  call ens_copyEnsMean(ensembleTrl4D, stateVectorMeanTrl4D)
  if (tim_nstepobsinc < tim_nstepobs) then
    call gsv_copy4Dto3D(stateVectorMeanTrl4D, stateVectorMeanAnl)
  else
    call gsv_copy(stateVectorMeanTrl4D, stateVectorMeanAnl)
  end if
  
  !
  !- 3. Compute HX values with results in ensObs/ensObsGain
  !

  !- 3.1 Loop over all members and compute HX for each
  if ( readEnsObsFromFile ) then
    call eob_readFromFiles(ensObs, nEns, inputFilenamePrefix='eob_HX', readObsInfo=.true.)
    if ( useModulatedEns ) then
      call enkf_setupModulationFactor(stateVectorMeanTrl4D%vco, numRetainedEigen, nEns, vLocalize, &
                                      beSilent=.true.)

      ! refresh assimilation flag before reading the files
      call eob_setAssFlag(ensObsGain)
      call eob_readFromFiles(ensObsGain, nEnsGain, inputFilenamePrefix='eobGain_HX', readObsInfo=.false.)
    end if
  else
    do memberIndex = 1, nEns
  
      write(*,*) ''
      write(*,*) 'midas-letkf: apply nonlinear H to ensemble member ', memberIndex
      write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

      ! copy 1 member to a stateVector
      call ens_copyMember(ensembleTrl4D, stateVector4D, memberIndex)
      call gsv_copy(stateVector4D, stateVectorWithZandP4D, allowVarMismatch_opt=.true., &
                    beSilent_opt=.true.)

      ! copy the surface height field
      if (nwpFields) then
        call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorWithZandP4D)
      end if

      ! Compute and set Yb in ensObs
      call s2c_nl( stateVectorWithZandP4D, obsSpaceData, column, hco_ens, &
                  timeInterpType=obsTimeInterpType, dealloc_opt=.false., &
                  beSilent_opt=.true. )

      ! Compute Y-H(X) in OBS_OMP
      call inn_computeInnovation(column, obsSpaceData, beSilent_opt=.true.)

      ! Copy to ensObs: Y-HX for this member
      call eob_setYb(ensObs, memberIndex)

      ! Compute and set Yb in ensObsGain
      do eigenVectorIndex = 1, numRetainedEigen
        if ( mmpi_myid == 0 .and. debug ) then
          write(*,*) 'midas-letkf: apply nonlinear H to modulated member ', &
                    eigenVectorIndex, '/', numRetainedEigen
        end if

        ! modulate the member with eigenvectors of vertical localization matrix
        call enkf_getModulatedState( stateVector4D, stateVectorMeanTrl4D, &
                                    vLocalize, numRetainedEigen, nEns, &
                                    eigenVectorIndex, stateVector4Dmod, &
                                    beSilent=.true. )
        if ( debug ) then
          call gsv_getField(stateVector4Dmod,field_Psfc,'P0')
          write(*,*) 'midas-letkf: max(Psfc)=', maxval(field_Psfc), &
                    ', min(Psfc)=', minval(field_Psfc)
        end if

        call gsv_copy(stateVector4Dmod, stateVectorWithZandP4D, allowVarMismatch_opt=.true., &
                      beSilent_opt=.true.)
        if (nwpFields) then
          call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorWithZandP4D)
        end if
        call s2c_nl( stateVectorWithZandP4D, obsSpaceData, column, hco_ens, &
                    timeInterpType=obsTimeInterpType, dealloc_opt=.false., &
                    beSilent_opt=.true. )

        ! Compute Y-H(X) in OBS_OMP
        call inn_computeInnovation( column, obsSpaceData, filterObsAndInitOer_opt=.false., &
                                    beSilent_opt=.true. )

        ! Copy to ensObsGain: Y-HX for this member
        memberIndexInEnsObs = (eigenVectorIndex - 1) * nEns + memberIndex
        call eob_setYb(ensObsGain, memberIndexInEnsObs)
      end do ! eigenVectorIndex

    end do ! memberIndex
  end if
  if ( gsv_isAllocated(stateVector4Dmod) ) call gsv_deallocate(stateVector4Dmod)

  ! write local ensObs to file
  if (writeLocalEnsObsToFile) then
    call eob_writeToFiles(ensObs, outputFilenamePrefix='eob_HX', writeObsInfo=.true.)
    if (useModulatedEns) then
      call eob_writeToFiles(ensObsGain, outputFilenamePrefix='eobGain_HX', writeObsInfo=.false., &
                            numGroupsToDivideMembers_opt=numRetainedEigen, &
                            maxNumMembersPerGroup_opt=nEns)
    end if
  end if

  !- 3.2 Set some additional information in ensObs/ensObsGain and additional quality
  !      control before finally communicating ensObs/ensObsGain globally

  ! Compute and remove the mean of Yb
  call eob_calcAndRemoveMeanYb(ensObs)
  if ( useModulatedEns ) call eob_calcAndRemoveMeanYb(ensObsGain)

  ! Put HPHT in OBS_HPHT, for writing to obs files
  call eob_setHPHT(ensObs)
  if ( useModulatedEns ) call eob_setHPHT(ensObsGain)

  ! Compute random observation perturbations
  if (trim(algorithm) == 'CVLETKF-PERTOBS') then
    randomSeedObs = 1 + mmpi_myid
    call eob_calcRandPert(ensObs, randomSeedObs)
  end if

  ! Apply obs operators to ensemble mean background for several purposes
  write(*,*) ''
  write(*,*) 'midas-letkf: apply nonlinear H to ensemble mean background'
  write(*,*) ''
  call gsv_copy(stateVectorMeanTrl4D, stateVectorWithZandP4D, allowVarMismatch_opt=.true.)
  if (nwpFields) then
    call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorWithZandP4D)
  end if
  call s2c_nl( stateVectorWithZandP4D, obsSpaceData, column, hco_ens, &
               timeInterpType=obsTimeInterpType, dealloc_opt=.false. )
  call tvs_allocTransmission(col_getNumLev(column,'TH')) ! this will cause radiative transmission profiles to be stored for use in eob_setVertLocation
  call inn_computeInnovation(column, obsSpaceData, beSilent_opt=.false.)

  ! Put y-mean(H(X)) in OBS_OMP for writing to obs files (overwrites y-H(mean(X)))
  call eob_setMeanOMP(ensObs)
  if ( useModulatedEns ) call eob_setMeanOMP(ensObsGain)

  ! Set y for family of interest to mean(H(x)) if doing simulated observations
  call eob_setSimObsVal(ensObs)

  ! Set vertical location for all obs for vertical localization (based on ensemble mean pressure and height)
  if (vLocalize > 0.0d0) then
    if (nwpFields) then
      call eob_setTypeVertCoord(ensObs,'logPressure')
      if ( useModulatedEns ) call eob_setTypeVertCoord(ensObsGain,'logPressure')
    else if (oceanFields) then
      call eob_setTypeVertCoord(ensObs,'depth')
      if ( useModulatedEns ) call eob_setTypeVertCoord(ensObsGain,'depth')
    end if
    call eob_setVertLocation(ensObs, column)
    if ( useModulatedEns ) call eob_setVertLocation(ensObsGain, column)
  end if

  ! Modify the obs error stddev for AMSUB in the tropics
  if (modifyAmsubObsError) call enkf_modifyAmsubObsError(obsSpaceData)

  ! Apply a background check (reject limit is set in the routine)
  if (backgroundCheck) call eob_backgroundCheck(ensObs)

  ! For ice/ocean DA: remove obs that are too close to land
  if (minDistanceToLand > 0.0D0) then
    call ens_getMask(ensembleTrl4D,oceanMask)
    call eob_removeObsNearLand(ensObs, oceanMask, minDistanceToLand)
  end if

  ! Set values of obs_sigi and obs_sigo before hubernorm modifies obs_oer
  call eob_setSigiSigo(ensObs)

  ! Apply huber norm quality control procedure (modifies obs_oer)
  if (huberize) call eob_huberNorm(ensObs)

  !- Reject all IR radiance observation in arctic and antarctic (.i.e |lat|>60. )
  if (rejectHighLatIR) call enkf_rejectHighLatIR(obsSpaceData)

  ! Reject radiance observations too close to the surface
  if (rejectRadNearSfc) call eob_rejectRadNearSfc(ensObs)

  ! Compute inverse of obs error variance (done here to use dynamic GPS-RO, GB-GPS based on mean O-P)
  call eob_setObsErrInv(ensObs)
  if ( useModulatedEns ) call eob_setObsErrInv(ensObsGain)

  call utl_tmg_start(141,'----Barr')
  call rpn_comm_barrier('GRID',ierr)
  call utl_tmg_stop(141)

  ! Clean and globally communicate obs-related data to all mpi tasks
  call eob_allGather(ensObs,ensObs_mpiglobal)
  if ( useModulatedEns ) then
    allocate(ensObsGain_mpiglobal)
    call eob_allGather(ensObsGain,ensObsGain_mpiglobal)
  else
    ensObsGain_mpiglobal => ensObs_mpiglobal
  end if

  ! Print number of assimilated obs per family to the listing
  write(*,*) 'oti_timeBinning: After extra filtering done in midas-letkf'
  call oti_timeBinning(obsSpaceData,tim_nstepobs)

  !
  !- 4. Final preparations for computing analyses
  !

  !- 4.1 Copy trial ensemble to nstepobsinc time steps
  if (tim_nstepobsinc < tim_nstepobs) then
    allocate(ensembleTrl)
    call ens_allocate(ensembleTrl, nEns, tim_nstepobsinc, hco_ens, vco_ens, dateStampListInc)
    call ens_copy4Dto3D(ensembleTrl4D,ensembleTrl)
    call ens_deallocate(ensembleTrl4D)
  else
    ! number of timesteps is the same, so just point to it
    ensembleTrl => ensembleTrl4D
  end if

  !- 4.2 Copy trl ensemble to anl ensemble
  call ens_allocate(ensembleAnl, nEns, tim_nstepobsinc, hco_ens, vco_ens, dateStampListInc)
  call ens_copy(ensembleTrl,ensembleAnl)

  !- 4.3 Setup for interpolating weights from coarse to full resolution
  call enkf_setupInterpInfo(wInterpInfo, stateVectorMeanAnl%hco, weightLatLonStep,  &
                            stateVectorMeanAnl%myLonBeg, stateVectorMeanAnl%myLonEnd,  &
                            stateVectorMeanAnl%myLatBeg, stateVectorMeanAnl%myLatEnd)

  !
  !- 5. Main calculation of ensemble analyses
  !

  !- 5.1 Call to perform LETKF
  call enkf_LETKFanalyses(algorithm, numSubEns, randomShuffleSubEns,  &
                          ensembleAnl, ensembleTrl, &
                          ensObs_mpiglobal, ensObsGain_mpiglobal, &
                          stateVectorMeanAnl, &
                          wInterpInfo, maxNumLocalObs,  &
                          hLocalize, hLocalizePressure, vLocalize, &
                          mpiDistribution, numRetainedEigen)

  !- 5.2 Loop over all analysis members and compute H(Xa_member) (if output is desired) 
  if ( outputEnsObs ) then
  
    do memberIndex = 1, nEns

      write(*,*) ''
      write(*,*) 'midas-letkf: apply nonlinear H to analysis ensemble member ', memberIndex
      write(*,*) ''

      ! copy 1 member to a stateVector
      call ens_copyMember(ensembleAnl, stateVectorMemberAnl, memberIndex)
      
      if (tim_nstepobsinc < tim_nstepobs) then
        ! ensembleAnl is only 3D, so need to make 4D for s2c_nl
        middleStepIndex = (tim_nstepobs + 1) / 2
        call gsv_copy(stateVectorMemberAnl, stateVectorWithZandP4D, allowVarMismatch_opt=.true., stepIndexOut_opt=middleStepIndex)
        call gsv_3dto4d(stateVectorWithZandP4D)
      else
        call gsv_copy(stateVectorMemberAnl, stateVectorWithZandP4D, allowVarMismatch_opt=.true.)
      end if
      
      ! copy the surface height field
      if (nwpFields) then
        call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorWithZandP4D)
      end if  
      
      call s2c_nl( stateVectorWithZandP4D, obsSpaceData, column, hco_ens, &
                   timeInterpType=obsTimeInterpType, dealloc_opt=.false. )

      ! Compute Y-H(Xa) in OBS_OMAM (used instead of OBS_OMA so that obsSpaceData isn't unintentionally modified) 
      call inn_computeInnovation(column, obsSpaceData, destObsColumn_opt=OBS_OMAM, beSilent_opt=.true., &
                                 callFiltTopo_opt=.false., callSetErrGpsgb_opt=.false., analysisMode_opt=.false.)

      ! Copy to ensObs: Y-H(Xa_member) for this member
      call eob_setYa(ensObs, memberIndex, OBS_OMAM)

    end do
  
  end if
  
  !- 6. Output obs files with mean OMP and (unrecentered) OMA

  ! Compute Y-H(Xa_mean) in OBS_OMA
  write(*,*) ''
  write(*,*) 'midas-letkf: apply nonlinear H to ensemble mean analysis'
  write(*,*) ''
  if (tim_nstepobsinc < tim_nstepobs) then
    ! meanAnl is only 3D, so need to make 4D for s2c_nl
    middleStepIndex = (tim_nstepobs + 1) / 2
    call gsv_copy(stateVectorMeanAnl, stateVectorWithZandP4D, allowVarMismatch_opt=.true., stepIndexOut_opt=middleStepIndex)
    call gsv_3dto4d(stateVectorWithZandP4D)
  else
    call gsv_copy(stateVectorMeanAnl, stateVectorWithZandP4D, allowVarMismatch_opt=.true.)
  end if
  if (nwpFields) then
    call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorWithZandP4D)
  end if
  call s2c_nl( stateVectorWithZandP4D, obsSpaceData, column, hco_ens, &
               timeInterpType=obsTimeInterpType )
  call inn_computeInnovation(column, obsSpaceData, destObsColumn_opt=OBS_OMA, beSilent_opt=.false.)
  
  ! Write (update) observation files. 
  if (outputEnsObs) then
    call obsf_writeFiles( obsSpaceData, ensObs_opt=ensObs)
  else
    call obsf_writeFiles( obsSpaceData)
  end if
  
  !- 7. Post processing of the analysis results (if desired) and write everything to files
  if (ensPostProcessing) then
    !- Allocate and read the Trl control member (used to compute control member increment for IAU)
    call gsv_allocate( stateVectorCtrlTrl, tim_nstepobsinc, hco_ens, vco_ens, &
                       dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    if (recenterInputEns) then
      !- Use the deterministic trial, if we are recentering the input ensemble
      ctrlFileName = trim(recenterFileName)
    else
      !- Otherwise, use member 0000
      call fln_ensFileName(ctrlFileName, ensPathName, memberIndex_opt=0, &
                           copyToRamDisk_opt=.false.)
    end if
    do stepIndex = 1, tim_nstepobsinc
      call gio_readFromFile( stateVectorCtrlTrl, ctrlFileName, ' ', ' ',  &
                             stepIndex_opt=stepIndex, containsFullField_opt=.true., &
                             readHeightSfc_opt=.false. )
    end do

    call epp_postProcess(ensembleTrl, ensembleAnl, stateVectorHeightSfc, stateVectorCtrlTrl, &
                         writeTrlEnsemble=.false., outputOnlyEnsMean_opt=outputOnlyEnsMean)
  else
    !- Just write the raw analysis ensemble to files
    if (mmpi_myid == 0) then
      write(*,*) 'midas-letkf: No ensemble post-processing requested, so just write the raw analysis ensemble'
    end if
    call utl_tmg_start(3,'--WriteEnsemble')
    if (.not. outputOnlyEnsMean) then
      call ens_writeEnsemble(ensembleAnl, '.', '', etiket_anl, 'A',  &
                             numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                             containsFullField_opt=.true.)
    end if
    call utl_tmg_stop(3)

  end if

  !
  !- 8. MPI, tmg finalize
  !  
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  call utl_tmg_stop(0)

  call tmg_terminate(mmpi_myid, 'TMG_INFO')
  call rpn_comm_finalize(ierr) 

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  if ( mmpi_myid == 0 ) write(*,*) ' --------------------------------'
  if ( mmpi_myid == 0 ) write(*,*) ' MIDAS-LETKF ENDS'
  if ( mmpi_myid == 0 ) write(*,*) ' --------------------------------'

end program midas_letkf

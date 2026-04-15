program midas_sstBias
  !
  !:Purpose: Main program to compute Sea Surface Temperature (SST)
  !          satellite observations bias estimate.
  !
  !          ---
  !
  !:Algorithm: The bias estimation of SST satellite observations is computed
  !            with respect to insitu observations that are considered unbiased.
  !            The bias estimation is produced for each sensor separately
  !            for day and night time.
  !
  !            --
  !
  !            First, each dataset is put on a regular grid using a small
  !            search radius (of ~25 km). It is currently a 1800x900 Gaussian grid.
  !            Second, the bias estimation at every gridpoint is computed as
  !            an average difference between satellite and insitu observations
  !            between all collocated valid satellite and insitu observations
  !            within a larger search radius (of ~1500km).
  !
  !            --
  !
  !            The resulting bias estimation :math:`B_{a}(k)` at point :math:`k`
  !            is computed as follows:
  !            :math:`B_{a}(k) = (1 - w(k)) * B_{b}(k) * \beta + w(k) * B_{a}(k)`,
  !            where :math:`B_{b}(k)` is a background state of the bias estimation
  !            computed on the previous day,
  !            :math:`\beta` is a background term for zero bias in unobserved areas,
  !            and :math:`w(k)` is a weight which is defined as:
  !            :math:`w(k) = N_{a}(k) / (N_{a}(k) + N_{b})`,
  !            where :math:`N_{a}(k)` is the number of observations involved in the
  !            computation of the current bias estimate :math:`B_{a}(k)` at point :math:`k`
  !            and :math:`N_{b}` is a parameter for corresponding number
  !            of observations used to compute the background state :math:`B_{b}(k)`.
  !
  !=========================================================== ========================================================
  ! Input and Output Files                                     Description of file
  !=========================================================== ========================================================
  ! ``analysisgrid``                                           In - File containing the grid where the bias is computed
  ! ``seaice_analysis``                                        In - File containing ``LG`` and ``VF`` fields
  ! ``obsfiles_$FAM/obs$FAM_$NNNN_$NNNN``                      In - Observation file for each "family" and MPI task
  ! ``searchRadius``                                           In - 'Large' search radius field to compute biases
  ! ``trlm_01``                                                In - Background state of the bias estimation
  ! ``satellite_bias.fst``                                     Out - Bias estimations
  ! ``auxOutput.fst``                                          Out - Auxiliary output (optional):
  !                                                            number of observations and weight fields
  !=========================================================== ========================================================
  !
  !           --
  !
  !:Synopsis: Below is a summary of the ``SSTbias`` program calling sequence:
  !
  !           - **Initial setups:**
  !
  !             - Setup horizontal and vertical grid objects for "analysis
  !               grid" from ``analysisgrid``.
  !
  !             - Setup ``obsSpaceData`` object and read observations from
  !               files: ``inn_setupObs``.
  !
  !             - Setup ``columnData`` and ``gridStateVector`` modules.
  !
  !           - **Computation**
  !
  !             - ``ocm_readMaskFromFile`` get the land-ocean mask
  !
  !             - ``oobs_computeObsData`` compute pseudo observation values
  !                and their coordinates and save them in SQLite files.
  !
  !             - ``sstb_getGriddedObs`` get all datasets on a regular grid
  !
  !             - ``sstb_getGriddedBias`` compute bias estimation for each sensor,
  !                for day time or night time on a regular grid,
  !                and save the results into an output standard file.
  !           --
  !
  !:Options: `List of namelist blocks <../namelists_in_each_program.html#sstbias>`_
  !          that can affect the ``SSTbias`` program.
  !
  !          * The use of ``SSTbias`` program is controlled by the namelist block
  !            ``&namSSTbiasEstimate`` read by the ``SSTbias`` program.
  !
  !              * ``iceFractionThreshold`` the sea-ice fraction threshold to define
  !                the presence of ice
  !
  !              * ``searchRadius`` horizontal search radius for observation gridding
  !
  !              * ``maxBias`` max allowed insitu-satellite difference in degrees
  !
  !              * ``numberPointsBG`` :math:`N_{b}`, number of points to compute
  !                the background bias estimation
  !
  !              * ``sensorList`` name of sensor
  !
  !              *  ``weightMin`` minimum value of weight
  !
  !              *  ``weightMax`` maximum value of weight
  !
  !              *  ``saveAuxFields`` to store or not auxiliary fields: nobs and weight
  !
  !              *  ``bgTermZeroBias`` background term to zero bias
  !
  !           --
  !
  use midasMpi_mod
  use version_mod
  use ramDisk_mod
  use utilities_mod
  use runtimeInfo_mod
  use message_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use obsSpaceData_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use obsFiles_mod
  use innovation_mod
  use sstBias_mod
  use columnData_mod

  implicit none

  integer                     :: dateStampFromObs
  type(struct_obs), target    :: obsSpaceData
  type(struct_hco), pointer   :: hco_anl => null()
  type(struct_vco), pointer   :: vco_anl => null()
  type(struct_columnData)     :: column
  character(len=48),parameter :: obsMpiStrategy = 'LIKESPLITFILES'
  character(len=48),parameter :: varMode        = 'analysis'

  call ver_printNameAndVersion('SSTbias','SST Bias Estimation')

  ! MPI initialization
  call mmpi_initialize

  call tmg_init(mmpi_myid, 'TMG_INFO')

  call rti_tmg_start(0,'Main')
  call rti_printTime()

  ! 1. Top level setup

  ! Read the namelists
  call utl_readNml()

  ! Setup the ramdisk directory (if supplied)
  call ram_setup()

  ! Setup the format of the output RPN standard files to 'XDF' or 'RSF'
  call gio_setup

  ! Do initial set up
  call SSTbias_setup('VAR') ! obsColumnMode

  call sstb_computeBias(obsSpaceData, hco_anl, vco_anl)

  ! Deallocate copied obsSpaceData
  call obs_finalize(obsSpaceData)
  call col_deallocate(column)

  ! 3. Job termination

  call rti_tmg_stop(0)
  call rti_printTime()

  call tmg_terminate(mmpi_myid, 'TMG_INFO')

  call mmpi_finalize

  contains

  subroutine SSTbias_setup(obsColumnMode)
    !
    ! :Purpose:  Control of the preprocessing of bias estimation
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in)  :: obsColumnMode

    ! Locals:
    character(len=*), parameter :: gridFile = './analysisgrid'

    write(*,*) ''
    write(*,*) '-------------------------------------------------'
    write(*,*) '-- Starting subroutine SSTbias_setup --'
    write(*,*) '-------------------------------------------------'

    !
    !- Initialize the Temporal grid and set dateStamp from env variable
    !
    call tim_setup()

    !
    !- Initialize observation file names and set dateStamp
    !
    call obsf_setup(dateStampFromObs, varMode)
    if (tim_getDateStamp() == 0) then
      if (dateStampFromObs > 0) then
        ! use dateStamp from obs if not set by env variable
        call tim_setDateStamp(dateStampFromObs)
      else
        call rti_abort('SSTbias_setup: DateStamp was not set')
      end if
    end if

    !
    !- Initialize constants
    !
    if(mmpi_myid == 0) call mpc_printConstants(6)

    !
    !- Initialize variables of the model states
    !
    call gsv_setup
    call msg_memUsage('midas-sstBias')

    !
    !- Initialize the Analysis grid
    !
    if(mmpi_myid == 0) write(*,*)''
    if(mmpi_myid == 0) write(*,*) 'SSTbias_setup: Set hco parameters for analysis grid'
    call hco_SetupFromFile(hco_anl, gridFile, 'GRID' ) ! IN

    !
    !- Initialisation of the analysis grid vertical coordinate from analysisgrid file
    !
    call vco_SetupFromFile(vco_anl, & ! OUT
                           gridFile ) ! IN

    call col_setVco(column, vco_anl)
    !
    !- Setup and read observations
    !
    call inn_setupObs(obsSpaceData, hco_anl, obsColumnMode, obsMpiStrategy, varMode) ! IN
    call msg_memUsage('midas-sstBias')

    !- Basic setup of columnData module
    call col_setup

    !- Memory allocation for background column data
    call col_allocate(column, obs_numHeader(obsSpaceData))

    if(mmpi_myid == 0) write(*,*) 'SSTbias_setup: done.'

  end subroutine SSTbias_setup

end program midas_sstBias

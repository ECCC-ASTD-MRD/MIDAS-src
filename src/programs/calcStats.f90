program midas_calcStats
  !
  !:Purpose: Swiss-knife type program originally used only for computing background-error covariances
  !          (**B**) based on homogeneous and isotropic (HI) correlations that was extended with time
  !          to compute various statistics and diagnostics from an ensemble of background-error estimates
  !
  !          ---
  !
  !:Algorithm: The ``calcStats`` program performs various types of statistics and diagnostics based on
  !            two modes defined through namelist options.
  !
  !              - **BHI**: Compute HI **B**. The approach for  limited-area applications is based on a bi-Fourier
  !                         spectral representation and the CVT described in <https://doi.org/10.1175/2010WAF2222401.1>.
  !                         For global applications, two formulations based on sperical-harmonics spectral representation
  !                         are available and controlled controlled by ``NAMCOMPUTEBHI``: 1) the legacy CVT formulation
  !                         described in <https://doi.org/10.1175/MWR-D-11-00097.1> and refereces therein,
  !                         and 2) the latbands CVT formulation which can also be run globally and, contrary to the
  !                         legacy formulation, does not impose balance operators. (2) is not restricted to the standard
  !                         weather fields and provides additional options controlled by ``NAMCOMPUTEBHILATBANDS``.
  !              - **TOOLBOX**: The swiss-knife component of this program controlled by ``NAMTOOLBOX`` from the
  !                global and LAM calcstats-related modules. Compute various statistics and diagnostics from
  !                an ensemble of background-error estimates in model-variable and/or control-variable space,
  !                such as:
  !                - HI vertical correlations
  !                - HI horizontal correlation functions
  !                - Variances
  !                - Local correlations
  !                - Optimal covariance localization radii
  !                - Power spectra
  !                Note that the above options are not all available in both global and limited-area
  !                applications.
  !            ---
  !
  !            Vertical coordinates input for the atmospheric case (i.e., when !! is provided): When the ensemble
  !            files do not include TT and UU, then TH and MM records are required at least in the ensemble member
  !            file with suffix $NNNN=0001 in order to set the vertical coordinates. Alternatively, if stats
  !            are not to be generated for TT an HU, TT and UU can be provided as substitutes for TH and MM
  !            in the file with $NNNN=0001.
  !
  !============================================= ==============================================================
  ! Input and Output Files                        Description of file
  !============================================= ==============================================================
  ! ``flnml``                                     In - Main namelist file with parameters user may modify
  ! ``ensemble/$YYYYMMDDHH_$HHH_$NNNN``           In - Background-error estimates
  ! ``various``                                   Out - Too many to be listed here
  !============================================= ==============================================================
  !
  !            --
  !
  !:Synopsis: Below is a summary of the ``calcStats`` program calling sequence:
  !
  !             - **Initial setups:**
  !
  !               - Read the NAMCONF and NAMENS namelist
  !
  !               - Setup horizontal and vertical grid objects using the first ensemble member
  !
  !               - Related modules are setup: ``gridStateVector_mod`` and ``timeCoord_mod``
  !
  !             - **Statistics and Diagnostics:**
  !
  !               - Intialize the module for global or limited-area applications
  !
  !               - Select the appropriate master routine based on the chosen mode
  !
  !               - The I/O and number crunching is performed within the global and limited-area modules
  !
  !
  !======================== ============ ==============================================================
  ! Module                   Namelist     Description of what is controlled
  !======================== ============ ==============================================================
  ! ``timeCoord_mod``        ``NAMTIME``  date of validity of the ensemble of background-error estimates
  !                                       the background state and increment
  ! ``calcstatslglb_mod``    ``mode-dependent``  Too many to be listed here
  ! ``calcstatslam_mod``     ``mode-dependent``  Too many to be listed here
  !======================== ============ ==============================================================
  !
  use midasMpi_mod
  use version_mod
  use fileNames_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use calcStatsGlb_mod
  use calcStatsLam_mod
  use utilities_mod
  use runtimeInfo_mod
  use ramDisk_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use timeCoord_mod

  implicit none

  type(struct_vco), pointer :: vco_ens => null()
  type(struct_hco), pointer :: hco_ens => null()

  character(len=256), parameter :: enspathname = 'ensemble'

  character(len=256) :: ensFileName
  character(len=4), pointer :: anlVar(:)

  ! namelist variables
  character(len=60) :: mode   ! can be 'BHI', 'TOOLBOX', 'STDDEV or 'POWERSPEC'
  integer           :: nens   ! Ensemble size
  integer           :: ip2    ! Ensemble lead time (hour) selected within the file

  NAMELIST /NAMCONF/mode
  NAMELIST /NAMENS/nens,ip2

  call ver_printNameAndVersion('calcStats','Compute the homogeneous-isotropic stats')

  !
  !- 1.  Initilization
  !

  !- 1.1 MPI and TMG
  call mmpi_initialize
  call tmg_init(mmpi_myid, 'TMG_INFO')

  call rti_tmg_start(0,'Main')
  call rti_printTime()

  ! Read the namelists
  call utl_readNml()

  ! Setup the ramdisk directory (if supplied)
  call ram_setup

  ! Setup time
  call tim_setup

  ! Setup the format of the output RPN standard files to 'XDF' or 'RSF'
  call gio_setup

  !- 1.2 Read NAMENS namelist
  nens              = 96                ! default value
  ip2               = -1                ! default value

  call rti_tmg_start(181,'low-level--readNML')
  read (utl_flnml, nml=namens)
  write(*, nml=namens)
  call rti_tmg_stop(181)

  !- 2.3 Initialize variables of the model states
  call gsv_setup

  call fln_ensfileName(ensFileName, ensPathName, memberIndex_opt=1)

  !- 1.3 Initialize the horizontal grid
  nullify(anlVar)
  call gsv_varNamesList(anlVar)
  call hco_SetupFromFile(hco_ens, ensFileName, ' ', 'Ensemble', varName_opt=anlVar(1)) ! IN

  !- 1.4 Initialize the vertical grid
  call vco_SetupFromFile( vco_ens,        & ! OUT
                          ensFileName, ' ') ! IN

  !- 1.5 Read NAMCONF namelist to find the mode
  mode  = 'BHI'  ! default value

  call rti_tmg_start(181,'low-level--readNML')
  read (utl_flnml, nml=namconf)
  write(*, nml=namconf)
  call rti_tmg_stop(181)

  !
  !- 2. Select and launch the appropriate mode
  !

  !- 2.1 Module initialization
  if (hco_ens % global) then
     call csg_setup( nens, hco_ens, vco_ens) ! IN
  else
     call csl_setup( nens, hco_ens, vco_ens, ip2) ! IN
  end if

  !- 2.2 Mode selection
  select case(trim(mode))
  case ('BHI')
     if (hco_ens % global) then
        call csg_computeBhi
     else
        call csl_computeBhi
     end if
  case ('TOOLBOX')
     if (hco_ens % global) then
        call csg_toolbox
     else
        call csl_toolbox
     end if
  case default
     write(*,*)
     write(*,*) 'Unknown value of MODE in global mode: ',mode
     call rti_abort('midas-calcstats')
  end select

  !
  !- 3.  Ending...
  !

  write(*,*)
  write(*,*) '---------------------'
  write(*,*) '> ENDING CALCBMATRIX '
  write(*,*) '---------------------'

  !
  !- 4.  MPI, tmg finalize
  !
  call rti_tmg_stop(0)
  call rti_printTime()

  call tmg_terminate(mmpi_myid, 'TMG_INFO')
  call mmpi_finalize

end program midas_calcStats

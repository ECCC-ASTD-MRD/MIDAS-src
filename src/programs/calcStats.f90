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
  !                         and 2) an experimental-only formulation that mimics the CVT model used in the limited-area mode.
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
  !               - Various modules are setup: ``gridStateVector_mod``, ``timeCoord_mod`` and ``bmatrix_mod``
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
 
  use version_mod
  use midasMpi_mod
  use fileNames_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use calcStatsGlb_mod
  use calcStatsLam_mod
  use utilities_mod
  use ramDisk_mod
  use gridStateVector_mod
  use timeCoord_mod
  implicit none

  type(struct_vco), pointer :: vco_ens => null()
  type(struct_hco), pointer :: hco_ens => null()

  character(len=256), parameter :: enspathname = 'ensemble'

  integer           :: fstopc
  integer           :: nulnam, ierr, fnom, fclos
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
  ierr = fstopc('MSGLVL','ERRORS',0)

  !- 1.1 MPI and TMG
  call mmpi_initialize
  call tmg_init(mmpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')

  ! Setup the ramdisk directory (if supplied)
  call ram_setup

  ! Setup time
  call tim_setup

  !- 1.2 Read NAMENS namelist
  nens              = 96                ! default value
  ip2               = -1                ! default value

  nulnam = 0
  ierr   = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read (nulnam,nml=namens)
  write(*     ,nml=namens)
  ierr=fclos(nulnam)

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

  nulnam = 0
  ierr   = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read (nulnam,nml=namconf)
  write(*     ,nml=namconf)
  ierr   = fclos(nulnam)

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
     call utl_abort('midas-calcstats')
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
  call utl_tmg_stop(0)

  call tmg_terminate(mmpi_myid, 'TMG_INFO')
  call rpn_comm_finalize(ierr) 

end program midas_calcStats

program midas_pseudoSSTobs
  !
  !:Purpose: Main program to produce pseudo SST observations 
  !          in ice-covered areas. Pseudo SST observations are needed
  !          to prevent the propagation of analysis increments to
  !          the ice-covered areas, that may result in undesirable sea-ice melting.
  !
  !          ---
  !
  !:Algorithm: Pseudo SST observations are assigned to the ice-covered 
  !            water points.
  !            First, a global sea-ice analysis is read.
  !            The sea-ice analysis file contains a mandatory sea-water
  !            fraction field.
  !            The grid and land-ocean mask are read 
  !            from the ``analysisgrid`` file.
  !
  !            --
  !
  !            Second, the number of ice-covered water points, including 
  !            concerned inland water points, are computed.
  !            If the number of ice-covered water points is zero,
  !            an empty observation SQLite file is created.
  !            If not, the computation of pseudo observations starts.
  !
  !            --
  !              
  !            First, the index array of ice-covered water points are 
  !            randomly shuffled to prevent the insertion of pseudo 
  !            observations at the same locations
  !            that would lead to spatial correlation of observations.
  !            Second, the pseudo observations of sea surface temperature :math:`T`
  !            are inserted at every ice-covered inland water point :math:`k`, 
  !            where the value of observations is computed as follows:
  !            :math:`T(k)=(1 - w(k)) * T_{fw} + w(k) * T_{s}`,
  !            where :math:`w(k)` is the sea-water fraction at the point :math:`k`,
  !            :math:`T_{fw}` is the temperature of fresh water below the ice,
  !            :math:`T_{s}` is a temperature of the sea water below the ice.
  !            The pseudo observations are inserted into every :math:`N`-th point 
  !            of sea water ice-covered points, 
  !            where the value of observation is defined as
  !            :math:`T_{s}`.
  !             
  !            --
  !  
  !            The computed observation values along with the corresponding 
  !            coordinates are put into ``obsSpaceData``. 
  !            Finally, output SQLite files are created.
  !       
  !            --
  !
  !=========================================================== ======================================================
  ! Input and Output Files                                     Description of file
  !=========================================================== ======================================================
  ! ``analysisgrid``                                           In - File defining sea-ice global grid
  ! ``seaice_analysis``                                        In - File containing ``LG`` and ``VF`` fields
  ! ``obsfiles_sst_pseudo.updated/obssst_pseudo_$NNNN_$NNNN``  Out - Pseudo obs file for each MPI task
  !=========================================================== ======================================================
  !
  !           --
  !
  !:Synopsis: Below is a summary of the ``pseudoSSTobs`` program calling sequence:
  !
  !           - **Initial setups:**
  !
  !             - Setup horizontal and vertical grid objects for "analysis
  !               grid" from ``analysisgrid``.
  !
  !             - Setup ``obsSpaceData`` object.
  !
  !             - Setup ``gridStateVector`` module.
  !
  !           - **Computation**
  !
  !             - ``utl_randomOrderInt`` random shuffle the ice-covered point indices
  !
  !             - ``oobs_computeObsData`` compute pseudo observation locations and values
  !                                       and save them in SQLite files.
  !
  !           --
  !
  !:Options: `List of namelist blocks <../namelists_in_each_program.html#pseudoSSTobs>`_
  !          that can affect the ``pseudoSSTobs`` program.
  !
  !          * The use of ``pseudoSSTobs`` program is controlled by the namelist block
  !            ``&pseudoSSTobs`` read by the ``pseudoSSTobs`` program.
  !
  !          * ``iceFractionThreshold`` the sea-ice fraction threshold to define 
  !                                      the presence of ice at each particular point
  !
  !          * ``outputSST`` the value of :math:`T_{s}` in K; 
  !
  !          * ``outputFreshWaterST`` the value of :math:`T_{fw}` in K;
  !
  !          * ``seaiceThinning`` pseudo observations are inserted into each :math:`N`-th point,
  !                                this parameter controls the observation thinning
  !
  !          * ``outputFileName`` controls the output file names
  !
  !          * ``etiket`` etiket to put into the table "resume" of output SQLite file
  ! 
  !          *  ``seaWaterThreshold`` a threshold to distinguish sea and fresh water
  !
  !          --
  !
  !========================== ================ ==============================================================
  ! Module                    Namelist         Description of what is controlled
  !========================== ================ ==============================================================
  ! ``oceanObservations_mod`` ``pseudoSSTobs`` parameters of pseudo SST observations
  !========================== ================ ==============================================================
  !
  use version_mod
  use ramDisk_mod
  use utilities_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use midasMpi_mod
  use oceanObservations_mod
  use gridStateVector_mod
  use obsSpaceData_mod
  use gridStateVectorFileIO_mod

  implicit none

  integer, external :: exdb, exfin
  integer :: ierr, istamp

  type(struct_hco), pointer   :: hco_anl => null()
  type(struct_vco), pointer   :: vco_anl => null()
  character(len=48),parameter :: obsMpiStrategy = 'LIKESPLITFILES'
  character(len=48),parameter :: varMode        = 'analysis'

  ! namelist variables
  real(8)                     :: iceFractionThreshold    ! consider no ice condition below this threshold
  real(8)                     :: outputSST               ! output SST value for pseudo observations
  real(8)                     :: outputFreshWaterST      ! output fresh water surface temperature for pseudo obs.  
  integer                     :: seaiceThinning          ! generate pseudo obs in every 'seaiceThinning' points 
  character(len=100)          :: outputFileName          ! name of the file containing the generated observations
  real(8)                     :: seaWaterThreshold       ! to distinguish inland water from sea water
  logical                     :: useSalinity             ! to use or not NEMO salinity field to compute freezing point temperature

  namelist /pseudoSSTobs/ iceFractionThreshold, outputSST, outputFreshWaterST, seaiceThinning, &
                          outputFileName, seaWaterThreshold, useSalinity
  
  istamp = exdb('pseudoSSTobs','DEBUT','NON')

  call ver_printNameAndVersion('pseudoSSTobs','Generation of pseudo SST observations')

  ! MPI initialization
  call mmpi_initialize

  call tmg_init(mmpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')
  call utl_printTime()

  ! 1. Top level setup

  ! Read the namelists
  call utl_readNml()

  !- RAM disk usage
  call ram_setup()
 
  ! Setup the format of the output RPN standard files to 'XDF' or 'RSF'
  call gio_setup

  ! Do initial set up
  call pseudoSSTobs_setup()

  call oobs_pseudoSST(hco_anl, vco_anl, iceFractionThreshold, outputSST, outputFreshWaterST, &
                      seaiceThinning, outputFileName, seaWaterThreshold, useSalinity)

  ! 3. Job termination

  istamp = exfin('pseudoSSTobs','FIN','NON')

  call utl_printTime()
  call utl_tmg_stop(0)

  call tmg_terminate(mmpi_myid, 'TMG_INFO')

  call mmpi_finalize

  contains

  subroutine pseudoSSTobs_setup()
    !
    ! :Purpose:  Control of the preprocessing of pseudo SST obs
    !
    implicit none
    
    ! Locals:	
    character(len=*), parameter :: gridFile = './analysisgrid'
    
    write(*,*) ''
    write(*,*) '-------------------------------------------------'
    write(*,*) '-- Starting subroutine pseudoSSTobs_setup --'
    write(*,*) '-------------------------------------------------'

    ! Setting default namelist variable values
    iceFractionThreshold   = 0.05d0
    outputSST              = 271.4d0
    outputFreshWaterST     = 271.4d0
    outputFileName         = ''
    seaiceThinning         = 5
    seaWaterThreshold      = 0.0d0
    useSalinity            = .false.

    ! Read the namelist
    call utl_tmg_start(181,'low-level--readNML')
    read(utl_flnml, nml = pseudoSSTobs, iostat = ierr)
    if (ierr /= 0) call utl_abort('pseudoSSTobs_setup: Error reading namelist')
    if (mmpi_myid == 0) write(*, nml = pseudoSSTobs)
    call utl_tmg_stop(181)

    write(*,*)''
    write(*,*) 'pseudoSSTobs_setup: output SST globally: ', outputSST
    write(*,*) 'pseudoSSTobs_setup: output fresh water surface temperature  globally: ', outputFreshWaterST
    write(*,*) 'pseudoSSTobs_setup: sea-ice fraction threshold: ', iceFractionThreshold
    write(*,*) 'pseudoSSTobs_setup: sea water fraction threshold: ', seaWaterThreshold
    if (useSalinity) then
      write(*,*) 'pseudoSSTobs_setup: Surface salinity field will be used to compute ocean temperature freezing point.'
      write(*,*) 'pseudoSSTobs_setup: WARNING: pseudo SST obs will be generatated in every grid point!'
      write(*,*) 'pseudoSSTobs_setup:          seaiceThinning value is put to 1.'
      seaiceThinning = 1
    else
      write(*,*) 'pseudoSSTobs_setup: pseudo SST obs will be generated in every ', seaiceThinning, ' points of the sea-ice field'    
    end if
    write(*,*) 'pseudoSSTobs_setup: output file name: ', outputFileName
    !
    !- Initialize the Analysis grid
    !
    if(mmpi_myid == 0) write(*,*)''
    if(mmpi_myid == 0) write(*,*) 'pseudoSSTobs_setup: Set hco parameters for analysis grid'
    call hco_SetupFromFile(hco_anl, gridFile, 'ANALYSIS') ! IN

    !     
    !- Initialisation of the analysis grid vertical coordinate from analysisgrid file
    !
    call vco_SetupFromFile(vco_anl, & ! OUT
                           gridFile) ! IN

    !- Initialize variables of the model states
    !
    call gsv_setup

    call obs_class_initialize('VAR')

    if(mmpi_myid == 0) write(*,*) 'pseudoSSTobs_setup: done.'
    
  end subroutine pseudoSSTobs_setup

end program midas_pseudoSSTobs

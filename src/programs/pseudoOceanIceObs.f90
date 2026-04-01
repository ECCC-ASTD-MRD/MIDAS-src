program midas_pseudoOceanIceObs
  !
  !:Purpose: Main program to produce pseudo SST observations
  !          in ice-covered areas and pseudo SIC observations
  !          around the ice edges and in the polynias/gaps in the ice.
  !          Pseudo SST observations are needed
  !          to prevent the propagation of analysis increments to
  !          the ice-covered areas, that may result in undesirable sea-ice melting.
  !          Pseudo SIC observations are needed to preserve a proper ice edge
  !          in coupled ocean-sea-ice data assimilation experiments.
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
  use midasMpi_mod
  use version_mod
  use ramDisk_mod
  use utilities_mod
  use runtimeInfo_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use oceanObservations_mod
  use gridStateVector_mod
  use obsSpaceData_mod
  use gridStateVectorFileIO_mod

  implicit none

  integer, external :: exdb, exfin
  integer :: ierr, istamp

  type(struct_hco), pointer   :: hco_anl => null()
  type(struct_vco), pointer   :: vco_anl => null()

  ! SST namelist variables
  logical                     :: computeSSTobs           ! a flag to switch the computation of SST pseudo obs on and off
  real(8)                     :: iceFractionThresholdSST ! consider no ice condition below this threshold
  real(8)                     :: outputSST               ! output SST value for pseudo observations
  real(8)                     :: outputFreshWaterST      ! output fresh water surface temperature for pseudo obs.
  integer                     :: seaiceThinning          ! generate pseudo obs in every 'seaiceThinning' points
  character(len=100)          :: outputFileNameSST       ! name of the file containing the generated SST observations
  real(8)                     :: seaWaterThreshold       ! to distinguish inland water from sea water
  logical                     :: useSalinity             ! to use or not NEMO salinity field to compute freezing point temperature
  ! SIC namelist variables
  logical                     :: computeSICobs           ! a flag to switch the computation of SIC pseudo obs on and off
  character(len=100)          :: outputFileNameSIC       ! name of the file containing the generated SIC observations
  real(8)                     :: iceFractionThresholdSIC ! consider no ice condition below this threshold
  real(8)                     :: seaIceBand              ! band in km around the ice edge where to put pseudo SIC observations

  namelist /pseudoSSTobs/ iceFractionThresholdSST, outputSST, outputFreshWaterST, seaiceThinning, &
                          outputFileNameSST, seaWaterThreshold, useSalinity, computeSSTobs
  namelist /pseudoSICobs/ outputFileNameSIC, computeSICobs, iceFractionThresholdSIC, seaIceBand

  istamp = exdb('pseudoOceanIceObs','DEBUT','NON')

  call ver_printNameAndVersion('pseudoOceanIceObs','Generation of pseudo SST observations')

  ! MPI initialization
  call mmpi_initialize

  call tmg_init(mmpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')
  call rti_printTime()

  ! 1. Top level setup

  ! Read the namelists
  call utl_readNml()

  !- RAM disk usage
  call ram_setup()

  ! Setup the format of the output RPN standard files to 'XDF' or 'RSF'
  call gio_setup

  ! Do initial set up
  call pseudoOceanIceObs_setup()

  if (computeSSTobs) then
    call oobs_pseudoSST(hco_anl, vco_anl, iceFractionThresholdSST, outputSST, outputFreshWaterST, &
                        seaiceThinning, outputFileNameSST, seaWaterThreshold, useSalinity)
  end if

  if (computeSICobs) then
    call oobs_pseudoSIC(hco_anl, vco_anl, iceFractionThresholdSIC, outputFileNameSIC, seaIceBand)
  end if

  ! 3. Job termination

  istamp = exfin('pseudoOceanIceObs','FIN','NON')

  call rti_printTime()
  call utl_tmg_stop(0)

  call tmg_terminate(mmpi_myid, 'TMG_INFO')

  call mmpi_finalize

  contains

  subroutine pseudoOceanIceObs_setup()
    !
    ! :Purpose:  Control of the preprocessing of pseudo SST and SIC obs
    !
    implicit none

    ! Locals:
    character(len=*), parameter :: gridFile = './analysisgrid'

    write(*,*) ''
    write(*,*) '-------------------------------------------------'
    write(*,*) '-- Starting subroutine pseudoOceanIceObs_setup --'
    write(*,*) '-------------------------------------------------'

    ! Setting default pseudoSSTobs namelist variable values
    computeSSTobs           = .true.
    iceFractionThresholdSST = 0.05d0
    outputSST               = 271.4d0
    outputFreshWaterST      = 271.4d0
    outputFileNameSST       = ''
    seaiceThinning          = 5
    seaWaterThreshold       = 0.0d0
    useSalinity             = .false.

    ! Setting default pseudoSICobs namelist variable values
    computeSICobs           = .false.
    outputFileNameSIC       = ''
    iceFractionThresholdSIC = 0.05d0
    seaIceBand              = 25.0d0

    ! Read the SST namelist
    call utl_tmg_start(181,'low-level--readNML')
    read(utl_flnml, nml = pseudoSSTobs, iostat = ierr)
    if (ierr /= 0) call rti_abort('pseudoOceanIceObs_setup: Error reading SST namelist')
    if (mmpi_myid == 0) write(*, nml = pseudoSSTobs)
    call utl_tmg_stop(181)

    ! Read the SIC namelist
    call utl_tmg_start(181,'low-level--readNML')
    read(utl_flnml, nml = pseudoSICobs, iostat = ierr)
    if (ierr /= 0) call rti_abort('pseudoOceanIceObs_setup: Error reading SIC namelist')
    if (mmpi_myid == 0) write(*, nml = pseudoSICobs)
    call utl_tmg_stop(181)

    if (computeSSTobs) then
      write(*,*) ''
      write(*,*) 'pseudoOceanIceObs_setup: computing SST pseudo observations is demanded.'
      write(*,*) 'pseudoOceanIceObs_setup: output SST globally: ', outputSST
      write(*,*) 'pseudoOceanIceObs_setup: output fresh water surface temperature  globally: ', outputFreshWaterST
      write(*,*) 'pseudoOceanIceObs_setup: sea-ice fraction threshold: ', iceFractionThresholdSST
      write(*,*) 'pseudoOceanIceObs_setup: sea water fraction threshold: ', seaWaterThreshold
      if (useSalinity) then
        write(*,*) 'pseudoOceanIceObs_setup: Surface salinity field will be used to compute ocean temperature freezing point.'
        write(*,*) 'pseudoOceanIceObs_setup: WARNING: pseudo SST obs will be generatated in every grid point!'
        write(*,*) 'pseudoOceanIceObs_setup:          seaiceThinning value is put to 1.'
        seaiceThinning = 1
      else
        write(*,*) 'pseudoOceanIceObs_setup: pseudo SST obs will be generated in every ', seaiceThinning, ' points of the sea-ice field'
      end if
      write(*,*) 'pseudoOceanIceObs_setup: SST output file name: ', outputFileNameSST
    end if

    if (computeSICobs) then
      write(*,*) ''
      write(*,*) 'pseudoOceanIceObs_setup: computing SIC pseudo observations is demanded.'
      write(*,*) 'pseudoOceanIceObs_setup: SIC output file name: ', outputFileNameSIC
      write(*,*) 'pseudoOceanIceObs_setup: sea-ice fraction threshold: ', iceFractionThresholdSIC
      write(*,*) 'pseudoOceanIceObs_setup: sea-ice band where to put pseudo SIC obs: ', seaIceBand
    end if

    !
    !- Initialize the Analysis grid
    !
    if(mmpi_myid == 0) write(*,*)''
    if(mmpi_myid == 0) write(*,*) 'pseudoOceanIceObs_setup: Set hco parameters for analysis grid'
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

    if(mmpi_myid == 0) write(*,*) 'pseudoOceanIceObs_setup: done.'

  end subroutine pseudoOceanIceObs_setup

end program midas_pseudoOceanIceObs

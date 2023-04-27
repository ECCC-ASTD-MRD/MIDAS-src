!--------------------------------------- LICENCE BEGIN -----------------------------------
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
  use analysisGrid_mod
  use midasMpi_mod
  use oceanObservations_mod
  use gridStateVector_mod
  use obsSpaceData_mod
   
  implicit none

  integer, external :: exdb, exfin, fnom, fclos, get_max_rss
  integer :: ierr, istamp, nulnam

  type(struct_hco), pointer   :: hco_anl => null()
  type(struct_vco), pointer   :: vco_anl => null()
  character(len=48),parameter :: obsMpiStrategy = 'LIKESPLITFILES'
  character(len=48),parameter :: varMode        = 'analysis'

  ! namelist variables
  real(4)                     :: iceFractionThreshold    ! consider no ice condition below this threshold
  real(4)                     :: outputSST               ! output SST value for pseudo observations
  real(4)                     :: outputFreshWaterST      ! output fresh water surface temperature for pseudo obs.  
  integer                     :: seaiceThinning          ! generate pseudo obs in every 'seaiceThinning' points 
  character(len=100)          :: outputFileName          ! name of the file containing the generated observations
  character(len=20)           :: etiket                  ! name to write in 'run' column of 'resume' table
  real(4)                     :: seaWaterThreshold       ! to distinguish inland water from sea water

  namelist /pseudoSSTobs/ iceFractionThreshold, outputSST, outputFreshWaterST, seaiceThinning, &
                          outputFileName, etiket, seaWaterThreshold
  
  istamp = exdb('pseudoSSTobs','DEBUT','NON')

  call ver_printNameAndVersion('pseudoSSTobs','Generation of pseudo SST observations')

  ! MPI initialization
  call mmpi_initialize

  call tmg_init(mmpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')
 
  ! 1. Top level setup

  call ram_setup()
 
  ! Do initial set up
  call pseudoSSTobs_setup()

  call oobs_pseudoSST(hco_anl, vco_anl, iceFractionThreshold, outputSST, outputFreshWaterST, &
                      seaiceThinning, outputFileName, etiket, seaWaterThreshold)

  ! 3. Job termination

  istamp = exfin('pseudoSSTobs','FIN','NON')

  call utl_tmg_stop(0)

  call tmg_terminate(mmpi_myid, 'TMG_INFO')

  call rpn_comm_finalize(ierr) 

  contains

  subroutine pseudoSSTobs_setup()
    !
    ! :Purpose:  Control of the preprocessing of pseudo SST obs
    !

    implicit none
    
    !Locals:	
    type(struct_hco), pointer   :: hco_core => null()
    character(len=*), parameter :: myName = 'pseudoSSTobs_setup'
    character(len=*), parameter :: gridFile = './analysisgrid'
    
    write(*,*) ''
    write(*,*) '-------------------------------------------------'
    write(*,*) '-- Starting subroutine '//myName//' --'
    write(*,*) '-------------------------------------------------'

    ! Setting default namelist variable values
    iceFractionThreshold   = 0.05 
    outputSST              = 271.4
    outputFreshWaterST     = 271.4
    outputFileName         = ''
    etiket                 = ''
    seaiceThinning         = 5
    seaWaterThreshold      = 0.0
    
    ! Read the namelist
    nulnam = 0
    ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    read(nulnam, nml = pseudoSSTobs, iostat = ierr)
    if (ierr /= 0) call utl_abort(myName//': Error reading namelist')
    if (mmpi_myid == 0) write(*, nml = pseudoSSTobs)
    ierr = fclos(nulnam)

    write(*,*)''
    write(*,*) myName//': output SST globally: ', outputSST
    write(*,*) myName//': output fresh water surface temperature  globally: ', outputFreshWaterST
    write(*,*) myName//': sea-ice fraction threshold: ', iceFractionThreshold
    write(*,*) myName//': sea water fraction threshold: ', seaWaterThreshold
    write(*,*) myName//': pseudo SST obs will be generated in every ', seaiceThinning, ' points of the sea-ice field'    
    write(*,*) myName//': output file name: ', outputFileName
    write(*,*) myName//': etiket for output SQLite files: ', etiket   
    !
    !- Initialize the Analysis grid
    !
    if(mmpi_myid == 0) write(*,*)''
    if(mmpi_myid == 0) write(*,*) myName//': Set hco parameters for analysis grid'
    call hco_SetupFromFile(hco_anl, gridFile, 'ANALYSIS') ! IN

    if (hco_anl % global) then
      call agd_SetupFromHCO(hco_anl) ! IN
    else
      !- Initialize the core (Non-Extended) analysis grid
      if(mmpi_myid == 0) write(*,*) myName//': Set hco parameters for core grid'
      call hco_SetupFromFile(hco_core, gridFile, 'COREGRID', 'AnalysisCore') ! IN
      !- Setup the LAM analysis grid metrics
      call agd_SetupFromHCO(hco_anl, hco_core) ! IN
    end if

    !     
    !- Initialisation of the analysis grid vertical coordinate from analysisgrid file
    !
    call vco_SetupFromFile(vco_anl, & ! OUT
                           gridFile) ! IN

    !- Initialize variables of the model states
    !
    call gsv_setup

    call obs_class_initialize('VAR')

    if(mmpi_myid == 0) write(*,*) myName//': done.'
    
  end subroutine pseudoSSTobs_setup

end program midas_pseudoSSTobs

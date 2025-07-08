program midas_sstTrial
  !
  !:Purpose: Main program to compute SST background state.
  !
  !          ---
  !
  !:Algorithm: The background state of SST is computed from the previous SST analysis field
  !            using a relaxation towards climatology as follows:
  !            :math:`X_{b}(t) = (X_{a}(t-1) - X_{clim}(t-1)) * \alpha + X_{clim}(t)`,
  !            where :math:`X_{b}(t)` is a resulting background state at time :math:`t`,
  !            :math:`X_{a}(t-1)` is analysis state at time :math:`t-1`,
  !            :math:`\alpha` is a relaxation coefficient, 
  !            :math:`X_{clim}(t)` is a climatology state at time :math:`t`,
  !            computed as an interpolation in time 
  !            between climatological field of the current and the following month
  !            for the current day of the month. 
  !           
  !            --
  ! 
  !=========================================================== ======================================================
  ! Input and Output Files                                     Description of file
  !=========================================================== ======================================================
  ! ``analysis``                                               In - SST analysis field
  ! ``climatology``                                            In - 12 monthly climatological SST fields 
  ! ``trial``                                                  Out - SST background (trial) field
  !=========================================================== ======================================================
  !
  !            --
  !
  !:Synopsis: Below is a summary of the ``SSTtrial`` program calling sequence:
  !
  !            - **Initial setups:**
  !
  !             - Setup horizontal and vertical grid objects for "analysis
  !               grid" from ``analysis``.
  !
  !             - Setup ``gridStateVector`` modules.
  !
  !             - Initialize a time objects ``tim_setup``
  !
  !             - Initialize datastamps from climatology file
  !
  !            - **Computation**
  !
  !             - ``obgd_getClimatology`` to read SST climatological fields from a standard file,
  !               to interpolate the field in time fot the current day :math:`t` in current month :math:`m` as follows   
  !               :math:`X_{clim}(t) = X_{clim}(m) + (t - 1) /(N - 1) * (X_{clim}(m+1) - X_{clim}(m))`,
  !               where :math:`N` is a number of days in the current month         
  !
  !             - ``obgd_computeSSTrial`` to compute the background field and save it into a standard file
  !
  !            --
  !
  !:Options: `List of namelist blocks <../namelists_in_each_program.html#ssttrial>`_
  !          that can affect the ``SSTtrial`` program.
  !
  !            * The use of ``SSTtrial`` program is controlled by the namelist block
  !              ``&NAMSSTTRIAL`` read by the ``SSTtrial`` program.
  !                * ``etiketAnalysis`` etiket to put into output standard file
  !
  !                * ``datestampClim`` list of twelve datestamps of climatological fields
  !                  in the ``climatology`` file
  !
  !                * ``alphaClim`` a parameter defining the relaxation towards climatology.
  !
  !            --
  !
  use midasMpi_mod
  use version_mod
  use ramDisk_mod
  use utilities_mod
  use message_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use oceanBackground_mod
  
  implicit none

  integer, external :: exdb, exfin
  integer :: ierr, istamp

  type(struct_hco), pointer :: hco_anl => null()
  type(struct_vco), pointer :: vco_anl => null()
  integer                   :: trialDateStamp    ! output: datestamp for the trial field
  integer                   :: analysisDateStamp !  input: datastamp of the latest analysis
  integer, parameter        :: nmonthsClim = 12  ! number of climatological fields

  ! namelist variables
  character(len=10) :: etiketAnalysis    ! etiket in the analysis file for grid setup
  integer           :: datestampClim(12) ! datestamps of input climatology fields 
  real(8)           :: alphaClim         ! scaling factor to relax towards climatology
  
  istamp = exdb('SSTTRIAL','DEBUT','NON')

  call ver_printNameAndVersion('SSTtrial','SST trial preparation')

  ! MPI initialization
  call mmpi_initialize

  call tmg_init(mmpi_myid, 'TMG_INFO')
  call utl_tmg_start(0,'Main')
  call utl_printTime()

  ! 1. Top level setup

  ! Read the namelists
  call utl_readNml()

  ! Setup the ramdisk directory (if supplied)
  call ram_setup()
 
  ! Setup the format of the output RPN standard files to 'XDF' or 'RSF'
  call gio_setup

  ! Do initial set up
  call SSTtrial_setup(trialDateStamp, analysisDateStamp)
  
  call obgd_computeSSTrial(hco_anl, vco_anl, trialDateStamp, analysisDateStamp, &
                           nmonthsClim, datestampClim, alphaClim, etiketAnalysis)

  ! 3. Job termination

  istamp = exfin('SSTTRIAL','FIN','NON')
  call utl_tmg_stop(0)
  call utl_printTime()
  call tmg_terminate(mmpi_myid, 'TMG_INFO')
  call mmpi_finalize

  contains
  
  !----------------------------------------------------------------------------------------
  ! SSTtrial_setup
  !----------------------------------------------------------------------------------------
  subroutine SSTtrial_setup(trialDateStamp, analysisDateStamp)
    !
    ! :Purpose:  Control of the preprocessing of trial
    !
    implicit none
    
    ! Arguments:
    integer, intent(out) :: trialDateStamp
    integer, intent(out) :: analysisDateStamp
    
    ! Locals:
    character(len=*), parameter :: gridFile = './analysis'
    integer                     :: prntdate, prnttime, imode, newdate, indexMonth
    namelist /namSSTtrial/ etiketAnalysis, datestampClim, alphaClim
        
    write(*,*) ''
    write(*,*) '-------------------------------------------------'
    write(*,*) '-- Starting subroutine SSTtrial_setup --'
    write(*,*) '-------------------------------------------------'
 
    ! namelist variables default values
    etiketAnalysis = ''
    datestampClim(:) = 0
    alphaClim = 0.983d0
    
    ! Read the namelist
    call utl_tmg_start(181,'low-level--readNML')
    read( utl_flnml, nml = namSSTtrial, iostat = ierr )
    if ( ierr /= 0) call utl_abort( 'SSTtrial_setup: Error reading namelist')
    if ( mmpi_myid == 0 ) write(*, nml = namSSTtrial )
    call utl_tmg_stop(181)

    if(mmpi_myid == 0) then
      write(*,'(1X,"***********************************")')
      write(*,'(1X," datestamps for climatology file:",/)')
      write(*,'(1X,"***********************************")')
      do indexMonth =1, nmonthsClim
        write(*,*) indexMonth, datestampClim(indexMonth)
      end do
      write(*,'(1X,"***********************************")')
    end if  

    !
    !- Initialize the Temporal grid and set dateStamp from gridFile
    !
    call tim_setup(fileNameForDate_opt = gridFile)
    analysisDateStamp = tim_getDateStamp()
    write(*,*) 'SSTtrial_setup: analysis datestamp  = ', analysisDateStamp
    write(*,*) 'SSTtrial_setup:          windowsize = ', tim_windowsize
    
    call incdatr(trialDateStamp, analysisDateStamp, tim_windowsize)
    write(*,*) 'SSTtrial_setup:    trial datestamp  = ', trialDateStamp
    
    imode = -3 ! stamp to printable
    ierr = newdate(trialDateStamp, prntdate, prnttime, imode)
    write(*,*) 'SSTtrial_setup: trial date = ', prntdate
    write(*,*) 'SSTtrial_setup: trial time = ', prnttime
    !
    !- Initialize variables of the model states
    !
    call gsv_setup
    call msg_memUsage('midas-sstTrial')

    !
    !- Initialize the Analysis grid
    !
    if(mmpi_myid == 0) write(*,*)''
    if(mmpi_myid == 0) write(*,*) 'SSTtrial_setup: Set hco parameters for analysis grid'
    call hco_SetupFromFile(hco_anl, gridFile, trim(etiketAnalysis)) ! IN

    !     
    !- Initialisation of the analysis grid vertical coordinate from analysisgrid file
    !
    call vco_SetupFromFile( vco_anl, & ! OUT
                            gridFile ) ! IN

    if(mmpi_myid == 0) write(*,*) 'SSTtrial_setup: done.'

  end subroutine SSTtrial_setup

end program midas_sstTrial

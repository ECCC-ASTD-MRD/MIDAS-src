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

program midas_sstTrial
  !
  ! :Purpose: Main program to compute SST background state
  !           xb(t) = (xa(t-1) - xclim(t-1))*alpha + xclim(t)
  !
  use version_mod
  use ramDisk_mod
  use utilities_mod
  use mpi_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use gridStateVector_mod
  use analysisGrid_mod
  use oceanBackground_mod
  
  implicit none

  integer, external :: exdb, exfin, fnom, fclos, get_max_rss
  integer :: ierr, istamp, nulnam

  type(struct_hco), pointer :: hco_anl => null()
  type(struct_vco), pointer :: vco_anl => null()
  integer                   :: trialDateStamp    ! output: datestamp for the trial field
  integer                   :: analysisDateStamp !  input: datastamp of the latest analysis
  integer, parameter        :: nmonthsClim = 12  ! number of climatological fields

  ! namelist variables
  character(len=10) :: etiketAnalysis    ! etiket in the analysis file for grid setup
  integer           :: datestampClim(12) ! datestamps of input climatology fields 
  real(4)           :: alphaClim         ! scalling factor to relax towards climatology
  
  istamp = exdb('SSTTRIAL','DEBUT','NON')

  call ver_printNameAndVersion('SSTtrial','SST trial preparation')

  ! MPI initialization
  call mpi_initialize

  call tmg_init(mpi_myid, 'TMG_INFO')
  call tmg_start(0,'Main')
  
  ! 1. Top level setup

  call ram_setup()
 
  ! Do initial set up
  call tmg_start(1,'Setup')
  call SSTtrial_setup(trialDateStamp, analysisDateStamp)
  call tmg_stop(1)
  
  call obgd_computeSSTrial(hco_anl, vco_anl, trialDateStamp, analysisDateStamp, &
                           nmonthsClim, datestampClim, alphaClim, etiketAnalysis)
			 
  ! 3. Job termination

  istamp = exfin('SSTTRIAL','FIN','NON')
  call tmg_stop(0)
  call tmg_terminate(mpi_myid, 'TMG_INFO')
  call rpn_comm_finalize(ierr) 

  contains
  
  !----------------------------------------------------------------------------------------
  ! SSTtrial_setup
  !----------------------------------------------------------------------------------------
  subroutine SSTtrial_setup(trialDateStamp, analysisDateStamp)
    !
    ! :Purpose:  Control of the preprocessing of trial
    !

    implicit none
    
    !Arguments:
    integer, intent(out) :: trialDateStamp
    integer, intent(out) :: analysisDateStamp
    
    !Locals:	
    type(struct_hco), pointer   :: hco_core => null()
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
    alphaClim = 0.983
    
    ! Read the namelist
    nulnam = 0
    ierr = fnom( nulnam, './flnml', 'FTN+SEQ+R/O', 0 )
    read( nulnam, nml = namSSTtrial, iostat = ierr )
    if ( ierr /= 0) call utl_abort( 'SSTtrial_setup: Error reading namelist')
    if ( mpi_myid == 0 ) write(*, nml = namSSTtrial )
    ierr = fclos(nulnam)

    if(mpi_myid == 0) then
      write(*,'(1X,"***********************************")')
      write(*,'(1X," datestamps for climatology file:",/)')
      write(*,'(1X,"***********************************")')
      do indexMonth =1, nmonthsClim
        write(*,*) indexMonth, datestampClim(indexMonth)
      end do
      write(*,'(1X,"***********************************")')
    end if  

    !
    !- Initialize the Temporal grid
    !
    call tim_setup(gridFile)
    analysisDateStamp = tim_getDatestampFromFile(gridFile)
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
    write(*,*) 'SSTtrial_setup: memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Initialize the Analysis grid
    !
    if(mpi_myid == 0) write(*,*)''
    if(mpi_myid == 0) write(*,*) 'SSTtrial_setup: Set hco parameters for analysis grid'
    call hco_SetupFromFile(hco_anl, gridFile, trim(etiketAnalysis)) ! IN

    if ( hco_anl % global ) then
      call agd_SetupFromHCO( hco_anl ) ! IN
    else
      !- Initialize the core (Non-Extended) analysis grid
      if(mpi_myid == 0) write(*,*) 'SSTtrial_setup: Set hco parameters for core grid'
      call hco_SetupFromFile( hco_core, gridFile, 'COREGRID', 'AnalysisCore' ) ! IN
      !- Setup the LAM analysis grid metrics
      call agd_SetupFromHCO( hco_anl, hco_core ) ! IN
    end if

    !     
    !- Initialisation of the analysis grid vertical coordinate from analysisgrid file
    !
    call vco_SetupFromFile( vco_anl, & ! OUT
                            gridFile ) ! IN

    if(mpi_myid == 0) write(*,*) 'SSTtrial_setup: done.'

  end subroutine SSTtrial_setup

end program midas_sstTrial

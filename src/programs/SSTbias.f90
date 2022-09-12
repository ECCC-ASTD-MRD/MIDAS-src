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

program midas_sstBias
  !
  ! :Purpose: Main program to compute SST bias estimate
  !
  use version_mod
  use ramDisk_mod
  use utilities_mod
  use midasMpi_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use obsSpaceData_mod
  use gridStateVector_mod
  use obsFiles_mod
  use innovation_mod
  use analysisGrid_mod
  use SSTbias_mod
  use columnData_mod
  
  implicit none

  integer, external :: exdb, exfin, fnom, fclos, get_max_rss
  integer :: ierr, istamp, nulnam

  type(struct_obs), target    :: obsSpaceData
  type(struct_hco), pointer   :: hco_anl => null()
  type(struct_vco), pointer   :: vco_anl => null()
  real(8)                     :: searchRadius            ! horizontal search radius, in km, for obs gridding
  character(len=48),parameter :: obsMpiStrategy = 'LIKESPLITFILES', &
                                 varMode        = 'analysis'
  real(4)                     :: iceFractionThreshold    ! consider no ice condition below this threshold
  real(4)                     :: maxBias                 ! max acceptable difference (insitu - satellite) 
  integer                     :: numberSensors           ! number of sensors to treat
  integer                     :: numberPointsBG          ! parameter, number of matchups of the background bias estimation
  character(len=10)           :: sensorList( 10 )        ! list of sensors
  integer                     :: dateStamp
  type(struct_columnData)     :: column                  ! column data 
  character(len=20)           :: timeInterpType_nl       ! 'NEAREST' or 'LINEAR'
  integer                     :: numObsBatches           ! number of batches for calling interp setup
   
  istamp = exdb('SSTBIASESTIMATION','DEBUT','NON')

  call ver_printNameAndVersion('SSTbias','SST Bias Estimation')

  ! MPI initialization
  call mmpi_initialize

  call tmg_init(mmpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')
 
  ! 1. Top level setup

  call ram_setup()
 
  ! Do initial set up
  call SSTbias_setup('VAR') ! obsColumnMode
  
  call sstb_computeBias(obsSpaceData, hco_anl, vco_anl, iceFractionThreshold, searchRadius, &
                        numberSensors, sensorList, maxBias, numberPointsBG, dateStamp)
			 
  ! Deallocate copied obsSpaceData
  call obs_finalize(obsSpaceData)
  call col_deallocate(column)

  ! 3. Job termination

  istamp = exfin('SSTBIAS','FIN','NON')

  call utl_tmg_stop(0)

  call tmg_terminate(mmpi_myid, 'TMG_INFO')

  call rpn_comm_finalize(ierr) 

  contains

  subroutine SSTbias_setup(obsColumnMode)
    !
    ! :Purpose:  Control of the preprocessing of bias estimation
    !

    implicit none
    !Arguments:
    character(len=*), intent(in)  :: obsColumnMode
    
    !Locals:	
    type(struct_hco), pointer   :: hco_core => null()
    character(len=*), parameter :: gridFile = './analysisgrid'
    integer                     :: sensorIndex
    namelist /namSSTbiasEstimate/ searchRadius, maxBias, iceFractionThreshold, numberPointsBG, &
                                  timeInterpType_nl, numObsBatches, numberSensors, sensorList
    
    write(*,*) ''
    write(*,*) '-------------------------------------------------'
    write(*,*) '-- Starting subroutine SSTbias_setup --'
    write(*,*) '-------------------------------------------------'

    !
    !- Initialize the Temporal grid
    !
    call tim_setup

    !     
    !- Initialize observation file names and set dateStamp
    !
    call obsf_setup(dateStamp, varMode)
    if (dateStamp <= 0) then
      call utl_abort('SSTbias_setup: Problem getting dateStamp from observation file')
    end if

    ! Setting default namelist variable values
    searchRadius = 10.            
    maxBias = 1.                  
    iceFractionThreshold   = 0.05 
    numberSensors = 0             
    numberPointsBG = 0            
    timeInterpType_nl = 'NEAREST'
    numObsBatches = 20
    sensorList(:) = ''
    
    ! Read the namelist
    nulnam = 0
    ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    read(nulnam, nml = namSSTbiasEstimate, iostat = ierr )
    if (ierr /= 0) call utl_abort('SSTbias_setup: Error reading namelist')
    if (mmpi_myid == 0) write(*, nml = namSSTbiasEstimate)
    ierr = fclos( nulnam )

    if (numberSensors == 0) call utl_abort('SSTbias_setup: Number of sensors to treat is not defined!')
    write(*,*)''
    write(*,*) 'SSTbias_setup: sensors to treat: '
    do sensorIndex = 1, numberSensors
      write(*,*) 'SSTbias_setup: sensor index: ', sensorIndex, ', sensor: ', sensorList( sensorIndex )
    end do
       
    !
    !- Initialize constants
    !
    if(mmpi_myid == 0) call mpc_printConstants(6)

    !
    !- Initialize variables of the model states
    !
    call gsv_setup
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Initialize the Analysis grid
    !
    if(mmpi_myid == 0) write(*,*)''
    if(mmpi_myid == 0) write(*,*) 'SSTbias_setup: Set hco parameters for analysis grid'
    call hco_SetupFromFile(hco_anl, gridFile, 'GRID' ) ! IN

    if ( hco_anl % global ) then
      call agd_SetupFromHCO( hco_anl ) ! IN
    else
      !- Initialize the core (Non-Extended) analysis grid
      if(mmpi_myid == 0) write(*,*) 'SSTbias_setup: Set hco parameters for core grid'
      call hco_SetupFromFile( hco_core, gridFile, 'COREGRID', 'AnalysisCore' ) ! IN
      !- Setup the LAM analysis grid metrics
      call agd_SetupFromHCO( hco_anl, hco_core ) ! IN
    end if

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
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !- Basic setup of columnData module
    call col_setup

    !- Memory allocation for background column data
    call col_allocate(column, obs_numHeader(obsSpaceData), mpiLocal_opt = .true.)

    if(mmpi_myid == 0) write(*,*) 'SSTbias_setup: done.'
    
  end subroutine SSTbias_setup

end program midas_sstBias

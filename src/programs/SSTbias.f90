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
  integer :: ierr, istamp

  type(struct_obs), target    :: obsSpaceData
  type(struct_hco), pointer   :: hco_anl => null()
  type(struct_vco), pointer   :: vco_anl => null()
  character(len=48),parameter :: obsMpiStrategy = 'LIKESPLITFILES', &
                                 varMode        = 'analysis'
  type(struct_columnData)     :: column                  ! column data
  integer                     :: dateStampFromObs
  
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
  
  call sstb_computeBias(obsSpaceData, hco_anl, vco_anl)
			 
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
    ! Arguments:
    character(len=*), intent(in)  :: obsColumnMode
    
    ! Locals:	
    type(struct_hco), pointer   :: hco_core => null()
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
        call utl_abort('SSTbias_setup: DateStamp was not set')
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

!-------------------------------------- LICENCE BEGIN ------------------------------------
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

!--------------------------------------------------------------------------
!!
!! *Purpose*: Main program for obervation background check and thinning
!!
!--------------------------------------------------------------------------
program midas_screening
  use ramDisk_mod
  use utilities_mod
  use mpivar_mod
  use MathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use obsSpaceData_mod
  use columnData_mod  
  use gridStateVector_mod
  use obsSpaceDiag_mod
  use controlVector_mod
  use burpFiles_mod
  use obsFilter_mod
  use innovation_mod
  use WindRotation_mod
  use analysisGrid_mod
  use tovs_nl_mod
  use obsErrors_mod
  use obsOperators_mod
  use multi_ir_bgck_mod
  implicit none

  integer :: istamp,exdb,exfin
  integer :: ierr,nconf

  type(struct_obs),       target :: obsSpaceData
  type(struct_columnData),target :: trlColumnOnAnlLev
  type(struct_columnData),target :: trlColumnOnTrlLev

  character(len=9) :: clmsg
  character(len=48) :: obsMpiStrategy, varMode

  NAMELIST /NAMCT0/NCONF
  integer nulnam, fnom, fclos 

  istamp = exdb('VAR','DEBUT','NON')

  write(*,'(/,' //                                                &
            '3(" *****************"),/,' //                       &
            '14x,"-- START OF MAIN PROGRAM MIDAS-SCREENING: --",/,' //   &
            '14x,"-- BACKGROUND CHECK AND SCREENING         --",/, ' //&
            '14x,"-- SCREENING Revision number   ",a," --",/,' //       &
            '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  ! MPI initilization
  call mpi_initialize  

  call tmg_init(mpi_myid, 'TMG_SCREENING' )

  call tmg_start(1,'MAIN')

  if(mpi_myid == 0) then
    clmsg = 'VAR3D_BEG'
    call utl_writeStatus(clmsg)
  endif 

  ! 1. Top level setup
  nconf             = 111

  nulnam=0
  ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  if(ierr.ne.0) call utl_abort('midas-screening: Error opening file flnml')
  read(nulnam,nml=namct0,iostat=ierr)
  if(ierr.ne.0) call utl_abort('midas-screening: Error reading namelist')
  write(*,nml=namct0)
  ierr=fclos(nulnam)

  select case(nconf)
  case (111)
    write(*,*)
    write(*,*) 'midas-screening: Background check for IR sat. data mode selected'
    varMode='bgckIR'
  case (101)
    write(*,*)
    write(*,*) 'midas-screening: Background check for conventional obs mode selected'
    varMode='bgckConv'
  case default
    write(*,*)
    write(*,*) 'midas-screening: Unknown mode ', nconf
    call utl_abort('midas-screening')
  end select

  call ram_setup

  ! 2. Decide on configuration of job

  ! ---BGCHECK (conventional obs)--- !
  if ( trim(varMode) == 'bgckConv' ) then
    if(mpi_myid == 0) write(*,*) 'MIDAS-SCREENING: CONVENTIONNAL BGCHECK MODE'

    ! Do initial set up
    call tmg_start(2,'PREMIN')

    obsMpiStrategy = 'LIKESPLITFILES'

    call screening_setup('ALL') ! obsColumnMode   

    ! Reading, horizontal interpolation and unit conversions of the 3D trial fields
    call inn_setupBackgroundColumns(trlColumnOnTrlLev,obsSpaceData)

    ! Interpolate trial columns to analysis levels and setup for linearized H
    call inn_setupBackgroundColumnsAnl(trlColumnOnTrlLev,trlColumnOnAnlLev)

    ! Compute observation innovations and prepare obsSpaceData for minimization
    call inn_computeInnovation(trlColumnOnTrlLev,obsSpaceData)
    call tmg_stop(2)

    ! Do the background check and output the observation data files
    call bgcheck_conv(trlColumnOnAnlLev,trlColumnOnTrlLev,obsSpaceData)

  ! ---BGCHECK (AIRS, IASI, CrIS)--- !
  else if ( trim(varMode) == 'bgckIR' ) then
    if(mpi_myid == 0) write(*,*) 'MIDAS-SCREENING: HYPERSPECTRAL IR BGCHECK MODE'

    ! Do initial set up
    call tmg_start(2,'PREMIN')

    obsMpiStrategy = 'LIKESPLITFILES'

    call screening_setup('ALL') ! obsColumnMode   

    ! Reading, horizontal interpolation and unit conversions of the 3D trial fields
    call inn_setupBackgroundColumns(trlColumnOnTrlLev,obsSpaceData)

    ! Compute observation innovations and prepare obsSpaceData for minimization
    call inn_computeInnovation(trlColumnOnTrlLev,obsSpaceData)
    call tmg_stop(2)

    ! Do the background check and output the observation data files
    call irbg_bgCheckIR(trlColumnOnTrlLev,obsSpaceData)

  else

    write(*,*) ' MIDAS-SCREENING: ERROR, UNKNOWN NCONF SPECIFIED'

  endif

  ! 3. Job termination

  istamp = exfin('VAR','FIN','NON')

  if(mpi_myid == 0) then
    clmsg = 'VAR3D_END'
    call utl_writeStatus(clmsg)
  endif

  call tmg_stop(1)

  call tmg_terminate(mpi_myid, 'TMG_SCREENING' )

  call rpn_comm_finalize(ierr) 

contains

  subroutine screening_setup(obsColumnMode)
    implicit none

    character (len=*) :: obsColumnMode
    integer :: datestamp
    type(struct_vco),pointer :: vco_anl => null()
    type(struct_vco),pointer :: vco_trl => null()
    type(struct_hco),pointer :: hco_anl => null()
    type(struct_hco),pointer :: hco_core => null()

    integer :: get_max_rss

    write(*,*) ''
    write(*,*) '-----------------------------------'
    write(*,*) '-- Starting subroutine var_setup --'
    write(*,*) '-----------------------------------'

    !
    !- Initialize the Temporal grid
    !
    call tim_setup

    !     
    !- Initialize burp file names and set datestamp
    !
    call burp_setupFiles (datestamp, varMode) ! IN
    call tim_setDatestamp(datestamp)            ! IN

    !
    !- Initialize constants
    !
    if(mpi_myid.eq.0) call mpc_printConstants(6)

    !
    !- Set vertical coordinate parameters from !! record in trial file
    !
    if(mpi_myid.eq.0) write(*,*)''
    if(mpi_myid.eq.0) write(*,*)' preproc: Set vcoord parameters for trial grid'
    call vco_SetupFromFile( vco_trl,     & ! OUT
                            './trlm_01')   ! IN
    call col_setVco(trlColumnOnTrlLev,vco_trl)

    !
    !- Initialize variables of the model states
    !
    call gsv_setup
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Initialize the Analysis grid
    !
    if(mpi_myid.eq.0) write(*,*)''
    if(mpi_myid.eq.0) write(*,*)' preproc: Set hco parameters for analysis grid'
    call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

    if ( hco_anl % global ) then
       call agd_SetupFromHCO( hco_anl ) ! IN
    else
       !- Iniatilized the core (Non-Exteded) analysis grid
       call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore' ) ! IN
       !- Setup the LAM analysis grid metrics
       call agd_SetupFromHCO( hco_anl, hco_core ) ! IN
    end if

    if ( hco_anl % rotated ) then
       call uvr_Setup(hco_anl) ! IN 
    end if

    !     
    !- Initialisation of the analysis grid vertical coordinate from analysisgrid file
    !
    call vco_SetupFromFile( vco_anl,        & ! OUT
                            './analysisgrid') ! IN

    call col_setVco(trlColumnOnAnlLev,vco_anl)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Setup and read observations
    !
    call inn_setupObs(obsSpaceData, obsColumnMode, obsMpiStrategy, varMode) ! IN
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Setup observation operators
    !
    call oop_setup(varMode) ! IN

    !
    !- Basic setup of columnData module
    !
    call col_setup
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Memory allocation for background column data
    !
    call col_allocate(trlColumnOnAnlLev,obs_numheader(obsSpaceData),mpi_local=.true.)
    call col_allocate(trlColumnOnTrlLev,obs_numheader(obsSpaceData),mpi_local=.true.)

    !
    !- Initialize the observation error covariances
    !
    call oer_setObsErrors(obsSpaceData, varMode) ! IN
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine screening_setup

end program midas_screening

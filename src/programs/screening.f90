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

  type(struct_obs),       target  :: obsSpaceData
  type(struct_columnData),target  :: trlColumnOnAnlLev
  type(struct_columnData),target  :: trlColumnOnTrlLev
  type(struct_vco),       pointer :: vco_anl  => null()
  type(struct_vco),       pointer :: vco_trl  => null()
  type(struct_hco),       pointer :: hco_anl  => null()
  type(struct_hco),       pointer :: hco_core => null()

  character(len=48) :: obsMpiStrategy, varMode
  character(len=3)  :: obsColumnMode

  integer :: ierr
  integer :: datestamp
  integer :: get_max_rss
  integer nulnam, fnom, fclos 

  ! Namelist
  integer :: nconf
  NAMELIST /NAMCT0/NCONF

  write(*,'(/,' //                                                &
       '3(" *****************"),/,' //                       &
       '14x,"-- START OF MAIN PROGRAM MIDAS-SCREENING: --",/,' //   &
       '14x,"-- BACKGROUND CHECK AND SCREENING         --",/, ' //&
       '14x,"-- SCREENING Revision number   ",a," --",/,' //       &
       '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  if ( mpi_myid == 0 ) then
    call utl_writeStatus('VAR3D_BEG')
  endif

  !
  !- 1.  Settings and module initializations
  !
  write(*,*)
  write(*,*) '> midas-screening: setup - START'

  !- 1.0 Namelist
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

  !- 1.1 mpi
  call mpi_initialize  

  !- 1.2 timings
  call tmg_init(mpi_myid, 'TMG_SCREENING' )
  call tmg_start(1,'MAIN')
  
  !- 1.3 RAM disk usage
  call ram_setup

  !- 1.4 Temporal grid
  call tim_setup

  !- 1.5 Initialize burp file names and set datestamp
  call burp_setupFiles (datestamp, varMode) ! IN
  call tim_setDatestamp(datestamp)            ! IN

  !- 1.6 Constants
  if ( mpi_myid == 0 ) call mpc_printConstants(6)

  !- 1.7 Setup a column vector following the background vertical grid
  if(mpi_myid.eq.0) write(*,*)''
  if(mpi_myid.eq.0) write(*,*)' preproc: Set vcoord parameters for trial grid'
  call vco_SetupFromFile( vco_trl,     & ! OUT
       './trlm_01')   ! IN
  call col_setVco(trlColumnOnTrlLev,vco_trl)

  !- 1.8 Variables of the model states
  call gsv_setup
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !- 1.9 Set the horizontal domain
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

  !- 1.10 Setup a column vector following the analysis vertical grid
  call vco_SetupFromFile( vco_anl,        & ! OUT
       './analysisgrid') ! IN

  call col_setVco(trlColumnOnAnlLev,vco_anl)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !- 1.11 Setup and read observations
  obsMpiStrategy = 'LIKESPLITFILES'
  obsColumnMode  = 'ALL'
  call inn_setupObs(obsSpaceData, obsColumnMode, obsMpiStrategy, varMode) ! IN
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !- 1.12 Setup observation operators
  call oop_setup(varMode) ! IN

  !- 1.13 Basic setup of columnData module
  call col_setup
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !- 1.14 Memory allocation for background column data
  call col_allocate(trlColumnOnAnlLev,obs_numheader(obsSpaceData),mpi_local=.true.)
  call col_allocate(trlColumnOnTrlLev,obs_numheader(obsSpaceData),mpi_local=.true.)

  !- 1.15 Initialize the observation error covariances
  call oer_setObsErrors(obsSpaceData, varMode) ! IN
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !- 1.16 Reading, horizontal interpolation and unit conversions of the 3D background fields
  call inn_setupBackgroundColumns(trlColumnOnTrlLev,obsSpaceData)

  !
  !- 2.  O-P computation
  !
  write(*,*)
  write(*,*) '> midas-screening: compute innovation'
  call inn_computeInnovation(trlColumnOnTrlLev,obsSpaceData)

  !
  !- 3.  Background Check
  !

  !-  3.1 Conventional data
  if ( trim(varMode) == 'bgckConv' ) then
    write(*,*)
    write(*,*) 'midas-screening: CONVENTIONNAL background check'

    ! Interpolate trial columns to analysis levels and setup for linearized H
    call inn_setupBackgroundColumnsAnl(trlColumnOnTrlLev,trlColumnOnAnlLev)

    ! Do the background check and output the observation data files
    call bgcheck_conv(trlColumnOnAnlLev,trlColumnOnTrlLev,obsSpaceData)

  !-  3.2 Hyperspectral IR (AIRS, IASI & CrIS)
  else if ( trim(varMode) == 'bgckIR' ) then
    write(*,*)
    write(*,*) 'midas-screening: HYPERSPECTRAL IR background check'

    ! Do the background check and output the observation data files
    call irbg_bgCheckIR(trlColumnOnTrlLev,obsSpaceData)

  endif

  !
  !- 4.  Ending
  !
  write(*,*)
  write(*,*) '> midas-screening: Ending'
  call obs_finalize(obsSpaceData) ! deallocate obsSpaceData

  call tmg_stop(1)
  call tmg_terminate(mpi_myid, 'TMG_SCREENING' )

  call rpn_comm_finalize(ierr) 

  if ( mpi_myid == 0 ) then
    call utl_writeStatus('VAR3D_END')
  endif

end program midas_screening

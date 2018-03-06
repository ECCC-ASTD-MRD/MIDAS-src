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
!! *Purpose*: Main program for Observation minus Forecast (O-F) computation
!!
!--------------------------------------------------------------------------
program midas_ominusf
  !
  ! **Purpose**: Main program for Observation minus Forecast (O-F) computation
  !
  use ramDisk_mod
  use utilities_mod
  use mpi_mod
  use mpiVar_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use obsSpaceData_mod
  use columnData_mod  
  use stateToColumn_mod
  use gridStateVector_mod
  use obsFiles_mod
  use obsFilter_mod  
  use innovation_mod
  use analysisGrid_mod
  use tovs_nl_mod
  use obsErrors_mod
  use obsOperators_mod
  use computeHBHT_mod
  implicit none

  type(struct_obs),       target  :: obsSpaceData
  type(struct_columnData),target  :: trlColumnOnAnlLev
  type(struct_columnData),target  :: trlColumnOnTrlLev
  type(struct_gsv)                :: stateVector
  type(struct_vco),       pointer :: vco_anl  => null()
  type(struct_vco),       pointer :: vco_trl  => null()
  type(struct_hco),       pointer :: hco_anl  => null()
  type(struct_hco),       pointer :: hco_core => null()
  type(struct_hco),       pointer :: hco_trl  => null()

  character(len=48) :: obsMpiStrategy
  character(len=3)  :: obsColumnMode
  character(len=10) :: trialFileName

  integer :: datestamp, get_max_rss, headerIndex, ierr
  integer :: fnom, fclos, nulnam

  character(len=20) :: timeInterpType_nl  ! 'NEAREST' or 'LINEAR'

  ! Namelist
  logical :: addHBHT  
  logical :: addSigmaO
  NAMELIST /NAMOMF/addHBHT,addSigmaO,timeInterpType_nl

  write(*,*) " --------------------------------------------"
  write(*,*) " ---  START OF MAIN PROGRAM midas-oMinusF ---"
  write(*,*) " ---  Computation of the innovation       ---"
  write(*,*) " ---  Revision: GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE "
  write(*,*) " --------------------------------------------"

  if ( mpi_myid == 0 ) then
    call utl_writeStatus('VAR3D_BEG')
  endif

  !
  !- 1.  Settings and module initializations
  !
  write(*,*)
  write(*,*) '> midas-OminusF: setup - START'

  obsMpiStrategy = 'LIKESPLITFILES'
  obsColumnMode  = 'VAR'
  trialFileName = './trlm_01'

  !- 1.0 Namelist
  addHBHT   = .false. ! default value
  addSigmaO = .false. 
  timeInterpType_nl='NEAREST'

  nulnam = 0
  ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read(nulnam,nml=namomf,iostat=ierr)
  if (ierr /= 0) call utl_abort('midas-OminusF: Error reading namelist')
  if (mpi_myid == 0) write(*,nml=namomf)
  ierr = fclos(nulnam)

  !- 1.1 mpi
  call mpi_initialize  

  !- 1.2 timings
  call tmg_init(mpi_myid, 'TMG_OMINUSF' )
  call tmg_start(1,'MAIN')

  !- 1.3 RAM disk usage
  call ram_setup

  !- 1.4 Temporal grid
  call tim_setup( fileNameForDate_opt=trim(trialFileName) )

  !- 1.5 Observation file names and get datestamp, but do not use it
  call obsf_setup(dateStamp, 'OminusF' )

  !- 1.6 Constants
  if (mpi_myid == 0) call mpc_printConstants(6)

  !- 1.7 Setup a column vector following the background vertical grid
  if (mpi_myid == 0) write(*,*)''
  if (mpi_myid == 0) write(*,*)' set vcoord parameters for background grid'
  call vco_SetupFromFile( vco_trl,             & ! OUT
                          trim(trialFileName) )  ! IN
  call col_setVco(trlColumnOnTrlLev,vco_trl)

  !- 1.8 Variables of the model states
  call gsv_setup

  !- 1.9 Set the horizontal domain
  if ( addHBHT ) then
    call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS') ! IN
    if ( hco_anl % global ) then
      call agd_SetupFromHCO( hco_anl ) ! IN
    else
      !- Iniatilized the core (Non-Exteded) analysis grid
      call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID') ! IN
      !- Setup the LAM analysis grid metrics
      call agd_SetupFromHCO( hco_anl, hco_core ) ! IN
    end if
  else
    call hco_SetupFromFile(hco_anl, trim(trialFileName), ' ') ! IN
    call agd_SetupFromHCO( hco_anl ) ! IN
  end if

  ! 1.10 Setup a column vector following the analysis vertical grid
  if ( addHBHT ) then
    call vco_SetupFromFile(vco_anl,        & ! OUT
                           './analysisgrid') ! IN
    call col_setVco(trlColumnOnAnlLev,vco_anl)
  end if

  !- 1.11 Setup and read observations
  call inn_setupObs(obsSpaceData, obsColumnMode, obsMpiStrategy, 'OminusF') ! IN

  !- 1.12 Setup observation operators
  call oop_setup('OminusF') ! IN

  !- 1.13 Basic setup of columnData module
  call col_setup

  !- 1.14 Memory allocation for background column data
  call col_allocate(trlColumnOnTrlLev,obs_numheader(obsSpaceData),mpiLocal_opt=.true.)
  if ( addHBHT ) then
    call col_allocate(trlColumnOnAnlLev, obs_numheader(obsSpaceData),mpiLocal_opt=.true.)
  end if

  if ( addSigmaO ) then
    !- 1.15 Initialize the observation error covariances
    write(*,*)
    write(*,*) '> midas-OminusF: Adding sigma_O'
    call oer_setObsErrors(obsSpaceData, 'OminusF')
  end if

  !- 1.16 Reading, horizontal interpolation and unit conversions of the 3D background fields
  call hco_SetupFromFile( hco_trl, trim(trialFileName), ' ', 'Trial' )
  call gsv_allocate( stateVector, tim_nstepobs, hco_trl, vco_trl,  &
                     dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                     mpi_distribution_opt='VarsLevs', dataKind_opt=4, &
                     allocGZsfc_opt=.true., hInterpolateDegree_opt='LINEAR' )
  call tmg_start(9,'readTrials')
  call gsv_readTrials( stateVector )
  call tmg_stop(9)
  call tmg_start(8,'s2c_nl')
  call s2c_nl( stateVector, obsSpaceData, trlColumnOnTrlLev,  &
               moveObsAtPole_opt=.true., timeInterpType_opt=timeInterpType_nl )
  call tmg_stop(8)
  call gsv_deallocate(stateVector)

  write(*,*)
  write(*,*) '> midas-OminusF: setup - END'
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- 2.  O-P computation
  !
 
  !- 2.1 Compute observation innovations
  write(*,*)
  write(*,*) '> midas-OminusF: compute innovation'
  call inn_computeInnovation(trlColumnOnTrlLev,obsSpaceData)

  if ( addHBHT ) then
    write(*,*)
    write(*,*) '> midas-OminusF: Adding HBH^T'
    !- 2.2 Interpolate background columns to analysis levels and setup for linearized H
    call inn_setupBackgroundColumnsAnl(trlColumnOnTrlLev,trlColumnOnAnlLev)

    ! 2.3 Compute the background errors in observation space
    call hbht_compute(trlColumnOnAnlLev,trlColumnOnTrlLev,obsSpaceData)
  end if

  ! 2.4 Write the results

  ! 2.4.1 Into the listings
  write(*,*)
  write(*,*) '> midas-OminusF: printing the FIRST header and body'
  do headerIndex = 1, min(1,obs_numHeader(obsSpaceData))
    call obs_prnthdr(obsSpaceData,headerIndex)
    call obs_prntbdy(obsSpaceData,headerIndex)
  end do
  ! 2.4.2 Into the observation files
  write(*,*)
  write(*,*) '> midas-OminusF: writing to file'
  call obsf_writeFiles(obsSpaceData)

  !
  !- 3.  Ending
  !
  write(*,*)
  write(*,*) '> midas-OminusF: Ending'
  call obs_finalize(obsSpaceData) ! deallocate obsSpaceData

  call tmg_stop(1)
  call tmg_terminate(mpi_myid, 'TMG_OMINUSF' )

  call rpn_comm_finalize(ierr) 

  if ( mpi_myid == 0 ) then
    call utl_writeStatus('VAR3D_END')
  endif

end program midas_ominusf

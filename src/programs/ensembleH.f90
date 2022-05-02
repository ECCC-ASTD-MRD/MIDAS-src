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

program midas_ensembleH
  !
  ! :Purpose: Main program for applying the observation operator to an ensemble
  !           of states as the first step for most EnKF algorithms.
  !
  use version_mod
  use mpi_mod
  use mpivar_mod
  use mathPhysConstants_mod
  use fileNames_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use columnData_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use timeCoord_mod
  use utilities_mod
  use ramDisk_mod
  use statetocolumn_mod
  use obsFiles_mod
  use obsSpaceData_mod
  use obsErrors_mod
  use innovation_mod
  use ensembleObservations_mod
  implicit none

  type(struct_obs), target :: obsSpaceData
  type(struct_gsv)         :: stateVector, statevector_tiles
  type(struct_columnData)  :: column
  type(struct_eob)         :: ensObs, ensObs_mpiglobal

  type(struct_vco), pointer :: vco_ens => null()
  type(struct_hco), pointer :: hco_ens => null()

  integer              :: fclos, fnom, fstopc, ierr
  integer              :: memberIndex, numStep, numBody
  integer              :: nulnam, dateStamp
  integer              :: get_max_rss

  character(len=256)  :: ensFileName
  character(len=9)    :: obsColumnMode
  character(len=48)   :: obsMpiStrategy
  character(len=48)   :: midasMode
  character(len=10)   :: obsFileType

  logical :: beSilent, dealloc

  ! namelist variables
  character(len=256) :: ensPathName
  integer  :: nEns
  NAMELIST /NAMENSEMBLEH/nEns, ensPathName

  midasMode = 'analysis'
  obsColumnMode = 'ENKFMIDAS'
  obsMpiStrategy = 'LIKESPLITFILES'

  call ver_printNameAndVersion('ensembleH','Program for applying H to ensemble')

  !
  !- 0. MPI, TMG initialization
  !
  call mpi_initialize
  call tmg_init(mpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Avoid printing lots of stuff to listing for std file I/O
  ierr = fstopc('MSGLVL','ERRORS',0)

  ! Setup the ramdisk directory (if supplied)
  call ram_setup

  ! Setting default namelist variable values
  nEns             = 10
  ensPathName      = 'ensemble'

  ! Read the namelist
  nulnam = 0
  ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namensembleh, iostat=ierr)
  if ( ierr /= 0) call utl_abort('midas-ensembleH: Error reading namelist')
  if ( mpi_myid == 0 ) write(*,nml=namensembleh)
  ierr = fclos(nulnam)

  ! Read the observations
  call obsf_setup( dateStamp, midasMode, obsFileType_opt = obsFileType )
  if ( obsFileType /= 'BURP' .and. obsFileType /= 'SQLITE' ) then
    call utl_abort('midas-ensembleH: only BURP and SQLITE are valid obs file formats')
  end if

  ! Use the first ensemble member to initialize datestamp and grid
  call fln_ensFileName( ensFileName, ensPathName, memberIndex_opt=1 )

  ! Setup timeCoord module, get datestamp from ensemble member
  call tim_setup( fileNameForDate_opt = ensFileName )
  numStep = tim_nstepobs

  !- Initialize variables of the model states
  call gsv_setup

  !- Initialize the Ensemble grid
  if (mpi_myid == 0) write(*,*) ''
  if (mpi_myid == 0) write(*,*) 'midas-ensembleH: Set hco and vco parameters for ensemble grid'
  call hco_SetupFromFile( hco_ens, ensFileName, ' ', 'ENSFILEGRID')
  call vco_setupFromFile( vco_ens, ensFileName )

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  ! read in the observations
  call inn_setupObs( obsSpaceData, hco_ens, obsColumnMode, obsMpiStrategy, midasMode )

  ! Initialize obs error covariances
  call oer_setObsErrors(obsSpaceData, midasMode)

  call col_setup
  call col_setVco(column, vco_ens)
  call col_allocate(column, obs_numheader(obsSpaceData))
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  ! Initialize the observation error covariances
  call oer_setObsErrors(obsSpaceData, midasMode) ! IN

  ! Allocate statevector to store an ensemble member (keep distribution as members on native grid)
  call gsv_allocate( stateVector, numStep, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='VarsLevs', &
                     dataKind_opt=4, allocHeightSfc_opt=.true. )

  call gsv_allocate( statevector_tiles, numStep, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(), &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeightSfc_opt=.true. )

  ! Allocate and initialize eob object for storing HX values
  numBody = obs_numBody(obsSpaceData)
  call eob_allocate(ensObs, nEns, numBody, obsSpaceData)
  call eob_zero(ensObs)

  ! Set lat, lon, obs values in ensObs
  call eob_setLatLonObs(ensObs)

  do memberIndex = 1, nEns
    write(*,*) ''
    write(*,*) 'midas-ensembleH: read member ', memberIndex
    call utl_tmg_start(2,'--ReadEnsemble')
    call fln_ensFileName( ensFileName, ensPathName, memberIndex_opt=memberIndex, copyToRamDisk_opt=.false.  )
    call gio_readFile( stateVector, ensFileName, ' ', ' ', containsFullField=.true., &
                       readHeightSfc_opt=.true. )
    call gio_fileUnitsToStateUnits( stateVector, containsFullField=.true. )

    call gsv_transposeVarsLevsToTiles(statevector, statevector_tiles)
    call utl_tmg_stop(2)

    write(*,*) ''
    write(*,*) 'midas-ensembleH: call s2c_nl for member ', memberIndex
    write(*,*) ''
    dealloc = .false.
    if ( memberIndex == nEns ) dealloc = .true.
    call s2c_nl( stateVector_tiles, obsSpaceData, column, hco_ens, timeInterpType='LINEAR', dealloc_opt=dealloc )

    write(*,*) ''
    write(*,*) 'midas-ensembleH: apply nonlinear H to member ', memberIndex
    write(*,*) ''
    beSilent = .true.
    if ( memberIndex == 1 ) beSilent = .false.
    call inn_computeInnovation(column, obsSpaceData, beSilent_opt=beSilent)

    ! Copy to ensObs: Y-HX for this member
    call eob_setYb(ensObs, memberIndex)

  end do

  call gsv_deallocate( stateVector_tiles )
  call gsv_deallocate( stateVector )

  ! Clean and globally communicate obs-related data, then write to files
  call eob_allGather(ensObs,ensObs_mpiglobal)
  call eob_writeToFiles(ensObs_mpiglobal)

  !
  !- MPI, tmg finalize
  !  
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  call utl_tmg_stop(0)

  call tmg_terminate(mpi_myid, 'TMG_INFO')
  call rpn_comm_finalize(ierr) 

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- 7.  Ending
  !
  if ( mpi_myid == 0 ) write(*,*) ' --------------------------------'
  if ( mpi_myid == 0 ) write(*,*) ' MIDAS-ENSEMBLEH ENDS'
  if ( mpi_myid == 0 ) write(*,*) ' --------------------------------'

end program midas_ensembleH

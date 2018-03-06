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
!! *Purpose*: Main program for applying the observation operator to an ensemble
!!            of states as the first step for most EnKF algorithms.
!!
!--------------------------------------------------------------------------
program midas_ensembleH
  !
  ! **Purpose**: Main program for applying the observation operator to an ensemble
  ! of states as the first step for most EnKF algorithms.
  !
  use mpi_mod
  use mpivar_mod
  use mathPhysConstants_mod
  use fileNames_mod
  use gridStateVector_mod
  use columnData_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use analysisGrid_mod
  use timeCoord_mod
  use utilities_mod
  use ramDisk_mod
  use enkf_mod
  use statetocolumn_mod
  use obsFiles_mod
  use obsSpaceData_mod
  use obsErrors_mod
  use obsOperators_mod
  use innovation_mod
  implicit none

  type(struct_obs), target             :: obsSpaceData
  type(struct_gsv)                     :: stateVector
  type(struct_columnData), allocatable :: columns(:)
  type(struct_columnData)              :: column_mean

  real(8), allocatable :: HXmean(:)
  real(8), allocatable :: HXens(:,:)
  real(8), pointer     :: HXensT_mpiglobal(:,:)
  real(8), allocatable :: obsVal(:)

  type(struct_vco), pointer :: vco_ens => null()
  type(struct_hco), pointer :: hco_ens => null()
  type(struct_hco), pointer :: hco_ens_core => null()

  integer              :: fclos, fnom, fstopc, ierr
  integer              :: memberIndex, stepIndex, numStep, procIndex, numBody, bodyIndex
  integer              :: nulnam, dateStamp
  integer              :: get_max_rss

  character(len=256)  :: ensFileName
  character(len=9)    :: obsColumnMode
  character(len=48)   :: obsMpiStrategy
  character(len=48)   :: midasMode
  character(len=10)   :: obsFileType

  logical :: beSilent, dealloc

  real(8), pointer    :: column_ptr(:)

  ! namelist variables
  character(len=256) :: ensPathName, ensFileBaseName
  logical  :: useTlmH, obsClean, asciDumpObs
  integer  :: nEns
  NAMELIST /NAMENSEMBLEH/nEns, ensPathName, ensFileBaseName, useTlmH, &
                         obsClean, asciDumpObs

  write(*,'(/,' //  &
        '3(" *****************"),/,' //                   &
        '14x,"-- START OF MIDAS-ENSEMBLEH             --",/,' //   &
        '14x,"-- Program for applying H to ensemble --",/, ' //  &
        '14x,"-- Revision number ",a," --",/,' //  &
        '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  midasMode = 'analysis'

  !
  !- 0. MPI, TMG initialization
  !
  call mpi_initialize
  call tmg_init(mpi_myid, 'TMG_ENSEMBLEH' )

  call tmg_start(1,'MAIN')
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Avoid printing lots of stuff to listing for std file I/O
  ierr = fstopc('MSGLVL','ERRORS',0)

  ! Setup the ramdisk directory (if supplied)
  call ram_setup

  ! Setting default namelist variable values
  nEns            = 10
  ensPathName     = 'ensemble'
  ensFileBaseName = ''
  useTlmH         = .false.
  obsClean        = .false.
  asciDumpObs     = .false.

  ! Read the namelist
  nulnam = 0
  ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namensembleh, iostat=ierr)
  if ( ierr /= 0) call utl_abort('midas-ensembleH: Error reading namelist')
  if ( mpi_myid == 0 ) write(*,nml=namensembleh)
  ierr = fclos(nulnam)

  if ( useTlmH ) call utl_abort('midas-ensembleH: WARNING use of TL of H not tested recently')

  

  ! Read the observations
  call obsf_setup( dateStamp, midasMode, obsFileType_opt = obsFileType )

  ! Use the first ensemble member to initialize datestamp and grid
  call fln_ensFileName( ensFileName, ensPathName, 1 )

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

  if ( hco_ens % global ) then
    call agd_SetupFromHCO( hco_ens ) ! IN
  else
    !- Setup the LAM analysis grid metrics
    hco_ens_core => hco_ens
    call agd_SetupFromHCO( hco_ens, hco_ens_core ) ! IN
  end if

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  call tmg_start(4,'SETUPOBS')
  obsColumnMode = 'ENKFMIDAS'
  ! determine the mpi strategy for observations, based on file type
  if ( obsFileType == 'BURP' .or. obsFileType == 'SQLITE' ) then
    obsMpiStrategy = 'LIKESPLITFILES'
  else
    obsMpiStrategy = 'ROUNDROBIN'
  end if
  ! read in the observations
  call inn_setupObs( obsSpaceData, obsColumnMode, obsMpiStrategy, midasMode,  &
                     obsClean_opt = obsClean )
  ! set up the observation operators
  call oop_setup(midasMode)
  call tmg_stop(4)

  call tmg_start(5,'ALLOC_COLS')
  write(*,*) ''
  write(*,*) 'midas-ensembleH: allocate an ensemble of columnData objects'
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  write(*,*) ''
  call col_setup
  allocate(columns(nEns))
  do memberIndex = 1, nEns
    beSilent = .true.
    if ( memberIndex == 1 ) beSilent = .false.
    call col_setVco(columns(memberIndex), vco_ens)
    call col_allocate(columns(memberIndex), obs_numheader(obsSpaceData),  &
                      mpiLocal_opt=.true., beSilent_opt=beSilent, setToZero_opt=.false.)
  end do
  if ( useTlmH ) then
    write(*,*) 'midas-ensembleH: allocating column_mean'
    call col_setVco(column_mean, vco_ens)
    call col_allocate(column_mean, obs_numheader(obsSpaceData), mpiLocal_opt=.true.)
  end if
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  call tmg_stop(5)

  ! Allocate statevector to store an ensemble member (keep distribution as members on native grid)
  call gsv_allocate( stateVector, numStep, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='VarsLevs', dataKind_opt=4, allocGZsfc_opt=.true. )

  do memberIndex = 1, nEns
    write(*,*) ''
    write(*,*) 'midas-ensembleH: read member ', memberIndex
    call fln_ensFileName( ensFileName, ensPathName, memberIndex )
    call tmg_start(3,'READ_ENSEMBLE')
    call gsv_readFile( stateVector, ensFileName, ' ', ' ', readGZsfc_opt=.true. )
    call gsv_fileUnitsToStateUnits( stateVector, containsFullField=.true. )
    ierr = ram_remove(ensFileName)
    call tmg_stop(3)

    write(*,*) ''
    write(*,*) 'midas-ensembleH: call s2c_nl for member ', memberIndex
    call tmg_start(6,'SETUP_COLS')
    if ( memberIndex == nEns ) then
      dealloc = .true.  ! some deallocation when doing the last member
    else
      dealloc = .false.
    end if
    call s2c_nl( stateVector, obsSpaceData, columns(memberIndex), dealloc_opt=dealloc )
    call tmg_stop(6)
  end do
  call gsv_deallocate( stateVector )

  ! Initialize the observation error covariances
  call oer_setObsErrors(obsSpaceData, midasMode) ! IN

  ! Allocate vectors for storing HX values
  numBody = obs_numBody(obsSpaceData)
  allocate(HXmean(numBody))
  allocate(HXens(numBody, nEns))
  allocate(obsVal(numBody))

  ! extract observation value, Y
  call enkf_extractObsRealBodyColumn(obsVal, obsSpaceData, OBS_VAR)

  ! ------------------------------------------------------------
  ! Compute H(X) for all members either with TLM or nonlinear H()
  ! ------------------------------------------------------------
  if ( useTlmH ) then

    ! Compute ensemble mean column and use nonlinear H()
    call enkf_computeColumnsMean(column_mean, columns)

    write(*,*) ''
    write(*,*) 'midas-ensembleH: apply nonlinear H to ensemble mean'
    write(*,*) ''
    ! compute Y-H(X) in OBS_OMP
    call tmg_start(7,'OBSOPER')
    call inn_computeInnovation(column_mean, obsSpaceData )
    call tmg_stop(7)

    ! extract observation-minus-HXmean value, Y-HXmean
    call enkf_extractObsRealBodyColumn(HXmean, obsSpaceData, OBS_OMP)
    ! compute HXmean = Y - (Y-HXmean)
    HXmean(:) = obsVal(:) - HXmean(:)

    ! Compute ensemble perturbation columns and use TL H()
    call enkf_computeColumnsPerturbations(columns, column_mean)
    write(*,*) ''
    write(*,*) 'midas-ensembleH: apply tangent linear H to ensemble perturbations'
    write(*,*) ''
    do memberIndex = 1, nEns
      ! compute H(dX) in OBS_WORK
      call tmg_start(7,'OBSOPER')
      call oop_Htl(columns(memberIndex), column_mean, obsSpaceData, memberIndex)
      call tmg_stop(7)

      ! extract HXpert value
      call enkf_extractObsRealBodyColumn(HXens(:,memberIndex), obsSpaceData, OBS_WORK)
      ! Recombine mean and perturbation to get total HX values
      HXens(:,memberIndex) = HXmean(:) + HXens(:,memberIndex)
    end do

  else

    ! Compute HX for all members using the nonlinear H()
    do memberIndex = 1, nEns
      ! compute Y-H(X) in OBS_OMP
      write(*,*) ''
      write(*,*) 'midas-ensembleH: apply nonlinear H to ensemble member ', memberIndex
      write(*,*) ''
      beSilent = .true.
      if ( memberIndex == 1 ) beSilent = .false.

      call tmg_start(7,'OBSOPER')
      call inn_computeInnovation(columns(memberIndex), obsSpaceData, beSilent_opt=beSilent)
      call tmg_stop(7)

      ! extract observation-minus-HX value, Y-HX
      call enkf_extractObsRealBodyColumn(HXens(:,memberIndex), obsSpaceData, OBS_OMP)
      ! compute HX = Y - (Y-HX)
      HXens(:,memberIndex) = obsVal(:) - HXens(:,memberIndex)

    end do

  end if ! useTlmH

  ! Put y-mean(H(X)) in OBS_OMP
  HXmean(:) = 0.0d0
  do memberIndex = 1, nEns
    do bodyIndex = 1, numBody
      HXmean(bodyIndex) = HXmean(bodyIndex) + HXens(bodyIndex,memberIndex)
    end do 
  end do
  HXmean(:) = HXmean(:)/nEns
  do bodyIndex = 1, numBody
    call obs_bodySet_r(obsSpaceData, OBS_OMP, bodyIndex,  &
                       obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex)-HXmean(bodyIndex))
  end do

  ! Gather obsSpaceData and HXens onto task 0
  call tmg_start(8,'OBS&HX_MPICOMM')
  if ( .not. obsf_filesSplit() ) then 
    call obs_expandToMpiGlobal(obsSpaceData)
  end if
  call enkf_gatherHX(HXens,HXensT_mpiglobal)
  call tmg_stop(8)

  ! Output mpiglobal H(X) and obsSpaceData files
  call tmg_start(9,'WRITEHXOBS')
  call obsf_writeFiles( obsSpaceData, HXensT_mpiglobal_opt = HXensT_mpiglobal, &
                        asciDumpObs_opt = asciDumpObs )
  call tmg_stop(9)

  !
  !- MPI, tmg finalize
  !  
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  call tmg_stop(1)

  call tmg_terminate(mpi_myid, 'TMG_ENSEMBLEH' )
  call rpn_comm_finalize(ierr) 

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !- 7.  Ending
  !
  if ( mpi_myid == 0 ) write(*,*) ' --------------------------------'
  if ( mpi_myid == 0 ) write(*,*) ' MIDAS-ENSEMBLEH ENDS'
  if ( mpi_myid == 0 ) write(*,*) ' --------------------------------'

end program midas_ensembleH

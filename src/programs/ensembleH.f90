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
  use mpi_mod
  use mpivar_mod
  use mathPhysConstants_mod
  use fileNames_mod
  use gridStateVector_mod
  use columnData_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use timeCoord_mod
  use utilities_mod
  use ramDisk_mod
  use enkf_mod
  use burpFiles_mod
  use obsSpaceData_mod
  use obsErrors_mod
  use obsOperators_mod
  use innovation_mod
  implicit none

  type(struct_gsv)                     :: statevector_member
  type(struct_obs), target             :: obsSpaceData
  type(struct_columnData), allocatable :: columns(:)
  type(struct_columnData)              :: column_mean

  real(4), allocatable :: HX_mean_r4(:)
  real(4), allocatable :: HX_ens_r4(:,:)
  real(4), allocatable :: obsVal_r4(:)

  type(struct_vco), pointer :: vco_ens => null()
  type(struct_hco), pointer :: hco_ens => null()
  type(struct_hco), pointer :: hco_ens_core => null()

  integer              :: fclos, fnom, fstopc, newdate, ierr
  integer              :: memberIndex, stepIndex, numStep, procIndex, numBody
  integer              :: idate, itime, nulnam
  integer              :: dateStamp, dateStampObs
  integer              :: get_max_rss

  character(len=256)  :: ensFileName
  character(len=3)    :: obsColumnMode
  character(len=48)   :: obsMpiStrategy

  logical :: beSilent

  real(8), pointer    :: column_ptr(:)

  ! namelist variables
  character(len=2)   :: ctrlVarHumidity
  character(len=256) :: ensPathName, ensFileBaseName
  logical  :: useTlmH
  integer  :: nEns, date
  NAMELIST /NAMENSEMBLEH/nEns, date, ensPathName, ensFileBaseName, ctrlVarHumidity, useTlmH

  write(*,'(/,' //  &
        '3(" *****************"),/,' //                   &
        '14x,"-- START OF MIDAS-ENSEMBLEH             --",/,' //   &
        '14x,"-- Program for applying H to ensemble --",/, ' //  &
        '14x,"-- Revision number ",a," --",/,' //  &
        '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

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
  date            = 0
  ensPathName     = 'ensemble'
  ensFileBaseName = ''
  ctrlVarHumidity = 'LQ'
  useTlmH         = .true.

  ! Read the namelist
  nulnam = 0
  ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namensembleh, iostat=ierr)
  if ( ierr /= 0) call utl_abort('midas-ensembleH: Error reading namelist')
  if ( mpi_myid == 0 ) write(*,nml=namensembleh)
  ierr = fclos(nulnam)

  ! Check if there are enough mpi tasks to read all members
  if ( mpi_nprocs < nEns ) then 
    write(*,*) 'nEns, mpi_nprocs = ', nEns, mpi_nprocs 
    call utl_abort('midas-ensembleH: Not enough mpi processes available to read all ensemble members')
  end if

  ! Initialize the analysis date-time from the namelist (for now)
  idate   = date/100
  itime   = (date-idate*100)*1000000
  ierr    = newdate(dateStamp, idate, itime, 3)
  if ( mpi_myid == 0 ) write(*,*) 'dateStamp from namelist= ', dateStamp

  ! Setup timeCoord module
  call tim_setup
  call tim_setDatestamp(dateStamp)
  numStep = tim_nstepobs

  !- Initialize variables of the model states
  call gsv_setup

  !- Initialize the Ensemble grid
  if (mpi_myid == 0) write(*,*) ''
  if (mpi_myid == 0) write(*,*) 'midas-ensembleH: Set hco and vco parameters for ensemble grid'
  ! Use the first ensemble member to initialize the grid
  call fln_ensFileName( ensFileName, ensPathName, 1 )
  call hco_SetupFromFile( hco_ens, ensFileName, ' ', 'ENSFILEGRID')
  call vco_setupFromFile( vco_ens, ensFileName )

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !- Setup and read the background ensemble (keep distribution as members on native grid)
  call tmg_start(2,'READ_ENSEMBLE')
  write(*,*) ''
  write(*,*) 'midas-ensembleH: allocating stateVector_member'
  write(*,*) ''
  call gsv_allocate(statevector_member, numStep, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                    mpi_local_opt=.false., dataKind_in_opt=4, allocGZsfc_opt=.true.)
  if ( (mpi_myid+1) <= nEns ) then
    write(*,*) ''
    write(*,*) 'midas-ensembleH: reading member from file'
    write(*,*) ''
    memberIndex = mpi_myid + 1
    call enkf_readMember(stateVector_member, ensPathName, memberIndex, ctrlVarHumidity)
  end if
  call tmg_stop(2)

  write(*,*) ''
  write(*,*) 'midas-ensembleH: read in the observations, keeping them distributed as is (RR)'
  write(*,*) ''
  obsColumnMode = 'VAR'
  obsMpiStrategy = 'LIKESPLITFILES'
  call burp_setupFiles(dateStampObs,'OminusF')
  call inn_setupObs(obsSpaceData, obsColumnMode, obsMpiStrategy, 'analysis')
  call oop_setup('analysis')

  write(*,*) ''
  write(*,*) 'midas-ensembleH: allocate an ensemble of columnData objects'
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  write(*,*) ''
  call col_setup
  if ( useTlmH ) then
    write(*,*) 'midas-ensembleH: allocating column_mean'
    call col_setVco(column_mean, vco_ens)
    call col_allocate(column_mean, obs_numheader(obsSpaceData), mpi_local=.true.)
  end if
  allocate(columns(nEns))
  do memberIndex = 1, nEns
    beSilent = .true.
    if ( memberIndex == 1 ) beSilent = .false.
    write(*,*) 'midas-ensembleH: allocating member ', memberIndex
    call col_setVco(columns(memberIndex), vco_ens)
    call col_allocate(columns(memberIndex), obs_numheader(obsSpaceData), mpi_local=.true., beSilent_opt=beSilent, setToZero_opt=.false.)
    write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  end do

  write(*,*) ''
  write(*,*) 'midas-ensembleH: Interpolate ensemble members to obs locations/times and redistribute'
  write(*,*) ''
  call tim_sutimeinterp(obsSpaceData, numStep)
  call enkf_setupColumnsFromEnsemble(stateVector_member, obsSpaceData, columns)

  ! Initialize the observation error covariances
  call oer_setObsErrors(obsSpaceData, 'analysis') ! IN

  ! Allocate vectors for storing HX values
  numBody = obs_numBody(obsSpaceData)
  allocate(HX_mean_r4(numBody))
  allocate(HX_ens_r4(numBody, nEns))
  allocate(obsVal_r4(numBody))

  ! extract observation value, Y
  call enkf_extractObsBodyColumn(obsVal_r4, obsSpaceData, OBS_VAR)

  ! Compute H(X) for all members either with TLM or nonlinear H()
  if ( useTlmH ) then

    ! Compute ensemble mean column and use nonlinear H()
    call enkf_computeColumnsMean(column_mean, columns)

    write(*,*) ''
    write(*,*) 'midas-ensembleH: apply nonlinear H to ensemble mean'
    write(*,*) ''
    ! compute Y-H(X) in OBS_OMP
    call inn_computeInnovation(column_mean, obsSpaceData)

    ! extract observation-minus-HXmean value, Y-HXmean
    call enkf_extractObsBodyColumn(HX_mean_r4, obsSpaceData, OBS_OMP)
    ! compute HXmean = Y - (Y-HXmean)
    HX_mean_r4(:) = obsVal_r4(:) - HX_mean_r4(:)

    ! Compute ensemble perturbation columns and use TL H()
    call enkf_computeColumnsPerturbations(columns, column_mean)
    write(*,*) ''
    write(*,*) 'midas-ensembleH: apply tangent linear H to ensemble perturbations'
    write(*,*) ''
    do memberIndex = 1, nEns
      ! compute H(dX) in OBS_WORK
      call oop_Htl(columns(memberIndex), column_mean, obsSpaceData, memberIndex)

      ! extract HXpert value
      call enkf_extractObsBodyColumn(HX_ens_r4(:,memberIndex), obsSpaceData, OBS_WORK)
      ! Recombine mean and perturbation to get total HX values
      HX_ens_r4(:,memberIndex) = HX_mean_r4(:) + HX_ens_r4(:,memberIndex)
    end do

    ! Output H(X) files
    call enkf_writeHXensemble(HX_ens_r4, obsSpaceData)

  else

    ! Compute HX for all members using the nonlinear H()
    do memberIndex = 1, nEns
      ! compute Y-H(X) in OBS_OMP
      write(*,*) ''
      write(*,*) 'midas-ensembleH: apply nonlinear H to ensemble member ', memberIndex
      write(*,*) ''
      call inn_computeInnovation(columns(memberIndex), obsSpaceData)

      ! extract observation-minus-HX value, Y-HX
      call enkf_extractObsBodyColumn(HX_ens_r4(:,memberIndex), obsSpaceData, OBS_OMP)
      ! compute HX = Y - (Y-HX)
      HX_ens_r4(:,memberIndex) = obsVal_r4(:) - HX_ens_r4(:,memberIndex)
    end do

    ! Need to reconcile assimilation flag between all members
    ! ????????????????

    ! Output H(X) files
    call enkf_writeHXensemble(HX_ens_r4, obsSpaceData)

  end if





  ! Deallocate things
  if ( stateVector_member%allocated ) call gsv_deallocate(statevector_member)

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

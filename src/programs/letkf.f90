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

program midas_letkf
  ! :Purpose: Main program for the local ensemble transform Kalman filter (LETKF).
  !           Note that the actual calculation of the analyses is in the local
  !           subroutine letkf_computeAnalyses, *contained* in this program.
  !           Many aspects of this program are controlled throught the namelist
  !           block NAMLETKF.
  use mpi, only : mpi_statuses_ignore ! this is the mpi library module
  use mpi_mod
  use mathPhysConstants_mod
  use fileNames_mod
  use ensembleObservations_mod
  use ensembleStateVector_mod
  use gridStateVector_mod
  use columnData_mod
  use localizationFunction_mod
  use tovs_nl_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use analysisGrid_mod
  use timeCoord_mod
  use utilities_mod
  use ramDisk_mod
  use statetocolumn_mod
  use obsFiles_mod
  use obsSpaceData_mod
  use obsErrors_mod
  use obsOperators_mod
  use innovation_mod
  use gridVariableTransforms_mod
  use physicsFunctions_mod
  use varNameList_mod
  use enkf_mod
  use humidityLimits_mod
  implicit none

  type(struct_obs), target  :: obsSpaceData
  type(struct_ens), target  :: ensembleTrl4D
  type(struct_ens), pointer :: ensembleTrl
  type(struct_ens)          :: ensembleAnl
  type(struct_gsv)          :: stateVectorMeanTrl4D
  type(struct_gsv)          :: stateVectorMeanTrl
  type(struct_gsv)          :: stateVectorMeanAnl
  type(struct_gsv)          :: stateVectorMeanInc
  type(struct_gsv)          :: stateVectorWithZandP4D
  type(struct_gsv)          :: stateVectorStdDevTrl
  type(struct_gsv)          :: stateVectorStdDevAnl
  type(struct_gsv)          :: stateVectorMeanTrlPressure
  type(struct_gsv)          :: stateVectorMeanTrlPressure_1step
  type(struct_gsv)          :: stateVectorMeanAnlSfcPressure
  type(struct_gsv)          :: stateVectorDeterTrl4D
  type(struct_gsv)          :: stateVectorDeterTrl
  type(struct_gsv)          :: stateVectorDeterAnl
  type(struct_gsv)          :: stateVectorDeterInc
  type(struct_gsv)          :: stateVectorHeightSfc
  type(struct_columnData)   :: column

  type(struct_eob) :: ensObs, ensObs_mpiglobal

  type(struct_vco), pointer :: vco_ens => null()
  type(struct_hco), pointer :: hco_ens => null()
  type(struct_hco), pointer :: hco_ens_core => null()

  integer :: memberIndex, stepIndex, middleStepIndex, nLev_M, nLev_T
  integer :: myLonBegHalo, myLonEndHalo, myLatBegHalo, myLatEndHalo
  integer :: myLonBeg, myLonEnd, myLatBeg, myLatEnd, numK
  integer :: nulnam, nulFile, dateStamp, datePrint, timePrint, imode, ierr, nsize
  integer :: get_max_rss, fclos, fnom, fstopc, newdate
  integer, allocatable :: dateStampList(:), dateStampListInc(:)

  character(len=256) :: ensFileName, outFileName, deterFileName
  character(len=9)   :: obsColumnMode
  character(len=48)  :: obsMpiStrategy
  character(len=48)  :: midasMode
  character(len=4)   :: memberIndexStr

  logical :: deterExists

  real(4), pointer     :: logPres_M_r4(:,:,:)

  ! interpolation information for weights (in enkf_mod)
  type(struct_enkfInterpInfo) :: wInterpInfo

  ! namelist variables
  character(len=256) :: ensPathName ! absolute or relative path to ensemble directory
  integer  :: nEns             ! ensemble size
  integer  :: maxNumLocalObs   ! maximum number of obs in each local volume to assimilate
  integer  :: weightLatLonStep ! separation of lat-lon grid points for weight calculation
  integer  :: randomSeed       ! seed used for random perturbation additive inflation
  logical  :: updateMembers    ! true means compute and write out members analyses/increments
  logical  :: writeIncrements  ! write ens of increments and mean increment
  real(8)  :: hLocalize        ! horizontal localization radius (in km)
  real(8)  :: vLocalize        ! vertical localization radius (in units of ln(Pressure in Pa))
  real(8)  :: alphaRTPP        ! RTPP coefficient (between 0 and 1; 0 means no relaxation)
  real(8)  :: alphaRTPS        ! RTPS coefficient (between 0 and 1; 0 means no relaxation)
  real(8)  :: alphaRandomPert  ! Random perturbation additive inflation coeff (0->1)
  logical  :: imposeSaturationLimit ! switch for choosing to impose saturation limit of humidity
  logical  :: imposeRttovHuLimits   ! switch for choosing to impose the RTTOV limits on humidity
  character(len=20) :: obsTimeInterpType ! type of time interpolation to obs time
  character(len=20) :: mpiDistribution   ! type of mpiDistribution for weight calculation ('ROUNDROBIN' or 'TILES')
  NAMELIST /NAMLETKF/nEns, ensPathName, hLocalize, vLocalize,  &
                     maxNumLocalObs, weightLatLonStep,  &
                     updateMembers, writeIncrements,  &
                     alphaRTPP, alphaRTPS, randomSeed, alphaRandomPert,  &
                     imposeSaturationLimit, imposeRttovHuLimits, obsTimeInterpType, &
                     mpiDistribution


  write(*,'(/,' //  &
        '3(" *****************"),/,' //                   &
        '14x,"-- START OF MIDAS-LETKF               --",/,' //   &
        '14x,"-- Program for Local Ensemble Transform Kalman Filter --",/, ' //  &
        '14x,"-- Revision number ",a," --",/,' //  &
        '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  ! Some high-level configuration settings
  midasMode = 'analysis'
  obsColumnMode = 'ENKFMIDAS'
  obsMpiStrategy = 'LIKESPLITFILES'

  !
  !- 0. MPI, TMG and misc. initialization
  !
  call mpi_initialize
  call tmg_init(mpi_myid, 'TMG_LETKF' )

  call tmg_start(1,'MAIN')
  call tmg_start(2,'LETKF-preAnl')
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Avoid printing lots of stuff to listing for std file I/O
  ierr = fstopc('MSGLVL','ERRORS',0)

  ! Setup the ramdisk directory (if supplied)
  call ram_setup

  !
  !- 1. Set/Read values for the namelist NAMLETKF
  !

  !- 1.1 Setting default namelist variable values
  nEns                  = 10
  maxNumLocalObs        = 1000
  weightLatLonStep      = 1
  ensPathName           = 'ensemble'
  updateMembers         = .false.
  writeIncrements       = .false.
  hLocalize             = 500.0D0
  vLocalize             = -1.0D0
  alphaRTPP             =  0.0D0
  alphaRTPS             =  0.0D0
  randomSeed            =  -999
  alphaRandomPert       =  0.0D0
  imposeSaturationLimit = .false.
  imposeRttovHuLimits   = .false.
  obsTimeInterpType     = 'LINEAR'
  mpiDistribution       = 'ROUNDROBIN'

  !- 1.2 Read the namelist
  nulnam = 0
  ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namletkf, iostat=ierr)
  if ( ierr /= 0) call utl_abort('midas-letkf: Error reading namelist')
  if ( mpi_myid == 0 ) write(*,nml=namletkf)
  ierr = fclos(nulnam)

  !- 1.3 Some minor modifications of namelist values
  hLocalize = hLocalize * 1000.0D0 ! convert from km to m
  if (alphaRTPP < 0.0D0) alphaRTPP = 0.0D0
  if (alphaRTPS < 0.0D0) alphaRTPS = 0.0D0
  if (alphaRandomPert < 0.0D0) alphaRandomPert = 0.0D0

  !
  !- 2.  Initialization
  !

  !- 2.1 Read the observations
  call obsf_setup( dateStamp, midasMode )

  !- 2.2 Initialize date/time-related info

  ! Setup timeCoord module, set dateStamp with value from obs files
  call tim_setup()
  if (tim_nstepobsinc /= 1 .and. tim_nstepobsinc /= tim_nstepobs) then
    call utl_abort('midas-letkf: invalid value for namelist variable DSTEPOBSINC. ' // &
                   'Increments can be either 3D or have same number of time steps as trials')
  end if
  call tim_setDateStamp(dateStamp)
  allocate(dateStampList(tim_nstepobs))
  call tim_getstamplist(dateStampList,tim_nstepobs,tim_getDatestamp())
  allocate(dateStampListInc(tim_nstepobsinc))
  call tim_getstamplist(dateStampListInc,tim_nstepobsinc,tim_getDatestamp())

  write(*,*) 'midas-letkf: analysis dateStamp = ',tim_getDatestamp()

  !- 2.3 Initialize variables of the model states
  call gsv_setup

  !- 2.4 Initialize the Ensemble grid
  if (mpi_myid == 0) write(*,*) ''
  if (mpi_myid == 0) write(*,*) 'midas-letkf: Set hco and vco parameters for ensemble grid'
  call fln_ensFileName( ensFileName, ensPathName, memberIndex_opt=1 )
  call hco_SetupFromFile( hco_ens, ensFileName, ' ', 'ENSFILEGRID')
  call vco_setupFromFile( vco_ens, ensFileName )
  nLev_M = vco_getNumLev(vco_ens, 'MM')
  nLev_T = vco_getNumLev(vco_ens, 'TH')
  if (nLev_M /= nLev_T) call utl_abort('midas-letkf: nLev_M /= nLev_T - currently not supported')

  if ( hco_ens % global ) then
    call agd_SetupFromHCO( hco_ens ) ! IN
  else
    ! Setup the LAM analysis grid metrics
    hco_ens_core => hco_ens
    call agd_SetupFromHCO( hco_ens, hco_ens_core ) ! IN
  end if

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !- 2.5 Read in the observations and other obs-related set up

  ! Read the observations
  call inn_setupObs( obsSpaceData, obsColumnMode, obsMpiStrategy, midasMode,  &
                     obsClean_opt = .false. )

  ! Set up the obs operators
  call oop_setup(midasMode)

  ! Initialize the observation error covariances
  call oer_setObsErrors(obsSpaceData, midasMode) ! IN

  ! Allocate vectors for storing HX values
  call eob_allocate(ensObs, nEns, obs_numBody(obsSpaceData), obsSpaceData)

  !- 2.6 Initialize a single columnData object
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  call col_setup
  call col_setVco(column, vco_ens)
  call col_allocate(column, obs_numheader(obsSpaceData),  &
                    mpiLocal_opt=.true., setToZero_opt=.true.)
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !- 2.7 Read the sfc height from ensemble member 1
  call gsv_allocate( stateVectorHeightSfc, 1, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0','TT'/) )
  call gsv_readFromFile( stateVectorHeightSfc, ensFileName, ' ', ' ',  &
                         containsFullField_opt=.true., readHeightSfc_opt=.true. )

  !- 2.8 Allocate various statevectors related to ensemble mean
  call gsv_allocate( stateVectorMeanTrl4D, tim_nstepobs, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeightSfc_opt=.true., &
                     allocHeight_opt=.false., allocPressure_opt=.false. )
  call gsv_zero(stateVectorMeanTrl4D)
  call gsv_allocate( stateVectorMeanTrl, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeightSfc_opt=.true., &
                     allocHeight_opt=.false., allocPressure_opt=.false. )
  call gsv_zero(stateVectorMeanTrl)
  call gsv_allocate( stateVectorMeanInc, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeightSfc_opt=.true., &
                     allocHeight_opt=.false., allocPressure_opt=.false. )
  call gsv_zero(stateVectorMeanInc)
  call gsv_allocate( stateVectorMeanAnl, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeightSfc_opt=.true., &
                     allocHeight_opt=.false., allocPressure_opt=.false. )
  call gsv_zero(stateVectorMeanAnl)

  !- 2.9 Allocate statevector for storing state with heights and pressures allocated (for s2c_nl)
  call gsv_allocate( stateVectorWithZandP4D, tim_nstepobs, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeightSfc_opt=.true. )
  call gsv_zero(stateVectorWithZandP4D)

  !- 2.10 Allocate ensembles, read the Trl ensemble
  call ens_allocate(ensembleTrl4D, nEns, tim_nstepobs, hco_ens, vco_ens, dateStampList)
  call ens_readEnsemble(ensembleTrl4D, ensPathName, biPeriodic=.false.)

  !- 2.11 Compute ensemble mean and copy to meanTrl and meanAnl stateVectors
  call ens_computeMean(ensembleTrl4D)
  call ens_copyEnsMean(ensembleTrl4D, stateVectorMeanTrl4D)
  if (tim_nstepobsinc < tim_nstepobs) then
    call gsv_copy4Dto3D(stateVectorMeanTrl4D, stateVectorMeanTrl)
  else
    call gsv_copy(stateVectorMeanTrl4D, stateVectorMeanTrl)
  end if
  call gsv_copy(stateVectorMeanTrl, stateVectorMeanAnl)

  !- 2.12 If deterministic background exists, do allocation and then read it
  call fln_ensFileName( deterFileName, ensPathName, shouldExist_opt=.false. )
  inquire(file=deterFileName, exist=deterExists)
  if (deterExists) then
    write(*,*)
    write(*,*) 'midas-letkf: Deterministic background state found, will provide deterministic analysis.'
    write(*,*) 'filename = ', deterFileName
    call gsv_allocate( stateVectorDeterTrl4D, tim_nstepobs, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_zero(stateVectorDeterTrl4D)
    call gsv_allocate( stateVectorDeterTrl, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_zero(stateVectorDeterTrl)
    call gsv_allocate( stateVectorDeterInc, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_zero(stateVectorDeterInc)
    call gsv_allocate( stateVectorDeterAnl, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_zero(stateVectorDeterAnl)
    do stepIndex = 1, tim_nstepobs
      call gsv_readFromFile( stateVectorDeterTrl4D, deterFileName, ' ', ' ',  &
                             stepIndex_opt=stepIndex, containsFullField_opt=.true., &
                             readHeightSfc_opt=.true. )
    end do
    if (tim_nstepobsinc < tim_nstepobs) then
      call gsv_copy4Dto3D(stateVectorDeterTrl4D, stateVectorDeterTrl)
    else
      call gsv_copy(stateVectorDeterTrl4D, stateVectorDeterTrl)
    end if
    call gsv_copy(stateVectorDeterTrl, stateVectorDeterAnl)
  else
    write(*,*)
    write(*,*) 'midas-letkf: No deterministic background state present.'
    write(*,*)
  end if

  !
  !- 3. Compute HX values with results in ensObs
  !

  !- 3.1 If it exists, compute HX for deterministic background
  if (deterExists) then

    write(*,*) ''
    write(*,*) 'midas-letkf: apply nonlinear H to deterministic background'
    write(*,*) ''

    ! copy deterministic background to state with pressure and heights allocated
    call gsv_copy(stateVectorDeterTrl4D, stateVectorWithZandP4D, allowMismatch_opt=.true.)

    ! Compute Y-H(X) in OBS_OMP
    call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorWithZandP4D)
    call s2c_nl( stateVectorWithZandP4D, obsSpaceData, column, timeInterpType=obsTimeInterpType )

    call tmg_start(6,'LETKF-obsOperators')
    call inn_computeInnovation(column, obsSpaceData, beSilent_opt=.true.)
    call tmg_stop(6)

    ! Copy to ensObs: Y-HX and various other obs quantities
    call eob_setLatLonObs(ensObs)
    call eob_setDeterYb(ensObs)

  end if

  !- 3.2 Loop over all members and compute HX for each
  do memberIndex = 1, nEns

    write(*,*) ''
    write(*,*) 'midas-letkf: apply nonlinear H to ensemble member ', memberIndex
    write(*,*) ''

    ! copy 1 member to a stateVector
    call ens_copyMember(ensembleTrl4D, stateVectorWithZandP4D, memberIndex)

    call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorWithZandP4D)
    call s2c_nl( stateVectorWithZandP4D, obsSpaceData, column, timeInterpType=obsTimeInterpType, dealloc_opt=.false. )

    ! Compute Y-H(X) in OBS_OMP
    call tmg_start(6,'LETKF-obsOperators')
    call inn_computeInnovation(column, obsSpaceData, beSilent_opt=.true.)
    call tmg_stop(6)

    ! Copy to ensObs: Y-HX and various other obs quantities
    if (memberIndex == 1) call eob_setLatLonObs(ensObs)
    call eob_setYb(ensObs, memberIndex)

  end do

  !- 3.3 Set some additional information in ensObs and then communicate globally

  ! Set assimilation flag in ensObs using value after applying H to all members
  call eob_setAssFlag(ensObs)

  !  Compute and remove the mean of Yb
  call eob_calcRemoveMeanYb(ensObs)

  ! Put y-mean(H(X)) in OBS_OMP, for writing to obs files
  call eob_setMeanOMP(ensObs)

  ! Apply obs operators to ensemble mean background for several purposes
  write(*,*) ''
  write(*,*) 'midas-letkf: apply nonlinear H to ensemble mean background'
  write(*,*) ''
  call gsv_copy(stateVectorMeanTrl4D, stateVectorWithZandP4D, allowMismatch_opt=.true.)
  call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorWithZandP4D)
  call s2c_nl( stateVectorWithZandP4D, obsSpaceData, column, timeInterpType=obsTimeInterpType, dealloc_opt=.false. )
  call tvs_allocTransmission ! this will cause radiative transmission profiles to be stored for use in eob_setPres
  call tmg_start(6,'LETKF-obsOperators')
  call inn_computeInnovation(column, obsSpaceData, beSilent_opt=.false.)
  call tmg_stop(6)

  ! Set pressure for all obs for vertical localization, based on ensemble mean pressure and height
  call eob_setPres(ensObs, column)

  ! Clean and then globally communicate obs-related data to all mpi tasks
  call eob_clean(ensObs)
  call eob_allGather(ensObs,ensObs_mpiglobal)

  call tmg_stop(2)

  !
  !- 4. Final preparations for computing analyses
  !

  !- 4.1 Adjust tile limits for weight calculation to include interpolation halo
  numK = stateVectorMeanInc%nk
  myLonBeg = stateVectorMeanInc%myLonBeg
  myLonEnd = stateVectorMeanInc%myLonEnd
  myLatBeg = stateVectorMeanInc%myLatBeg
  myLatEnd = stateVectorMeanInc%myLatEnd
  myLonBegHalo = 1 + weightLatLonStep * floor(real(myLonBeg - 1)/real(weightLatLonStep))
  myLonEndHalo = min(stateVectorMeanInc%ni, 1 + weightLatLonStep * ceiling(real(myLonEnd - 1)/real(weightLatLonStep)))
  myLatBegHalo = 1 + weightLatLonStep * floor(real(myLatBeg - 1)/real(weightLatLonStep))
  myLatEndHalo = min(stateVectorMeanInc%nj, 1 + weightLatLonStep * ceiling(real(myLatEnd - 1)/real(weightLatLonStep)))
  write(*,*) 'midas-letkf: myLonBeg/End, myLatBeg/End (original)  = ',  &
             myLonBeg, myLonEnd, myLatBeg, myLatEnd
  write(*,*) 'midas-letkf: myLonBeg/End, myLatBeg/End (with Halo) = ',  &
             myLonBegHalo, myLonEndHalo, myLatBegHalo, myLatEndHalo
  write(*,*) 'midas-letkf: number of local gridpts where weights computed = ',  &
             ( 1 + ceiling(real(myLonEndHalo - myLonBegHalo) / real(weightLatLonStep)) ) *  &
             ( 1 + ceiling(real(myLatEndHalo - myLatBegHalo) / real(weightLatLonStep)) )

  !- 4.2 Copy trial ensemble to nstepobsinc time steps
  if (tim_nstepobsinc < tim_nstepobs) then
    allocate(ensembleTrl)
    call ens_allocate(ensembleTrl, nEns, tim_nstepobsinc, hco_ens, vco_ens, dateStampListInc)
    call ens_copy4Dto3D(ensembleTrl4D,ensembleTrl)
    call ens_deallocate(ensembleTrl4D)
  else
    ! number of timesteps is the same, so just point to it
    ensembleTrl => ensembleTrl4D
  end if

  !- 4.3 Copy trl ensemble to anl ensemble
  call ens_allocate(ensembleAnl, nEns, tim_nstepobsinc, hco_ens, vco_ens, dateStampListInc)
  call ens_copy(ensembleTrl,ensembleAnl)

  !- 4.4 Compute background ens mean 3D log pressure and make mpiglobal for vertical localization
  call gsv_allocate( stateVectorMeanTrlPressure, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0','P_M','P_T'/) )
  call gsv_zero(stateVectorMeanTrlPressure)
  call gsv_copy(stateVectorMeanTrl, stateVectorMeanTrlPressure, allowMismatch_opt=.true.)
  call gvt_transform(stateVectorMeanTrlPressure,'PsfcToP_nl')
  if (mpi_myid == 0) then
    call gsv_allocate( stateVectorMeanTrlPressure_1step, 1, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.false., &
                       dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0','P_M','P_T'/) )
  end if
  call gsv_transposeTilesToStep(stateVectorMeanTrlPressure_1step, stateVectorMeanTrlPressure, (tim_nstepobsinc+1)/2)
  call gsv_deallocate(stateVectorMeanTrlPressure)
  if (mpi_myid == 0) then
    logPres_M_r4 => gsv_getField3d_r4(stateVectorMeanTrlPressure_1step,'P_M')
    logPres_M_r4(:,:,:) = log(logPres_M_r4(:,:,:))
  else
    allocate(logPres_M_r4(stateVectorMeanTrlPressure%ni, stateVectorMeanTrlPressure%nj, nLev_M))
  end if
  nsize = stateVectorMeanTrlPressure%ni * stateVectorMeanTrlPressure%nj * nLev_M
  call rpn_comm_bcast(logPres_M_r4, nsize, 'mpi_real4', 0, 'GRID', ierr)

  !- 4.5 Setup for interpolating weights from coarse to full resolution
  call tmg_start(92,'LETKF-interpolateWeights')
  call enkf_setupInterpInfo(wInterpInfo, stateVectorMeanInc%hco, weightLatLonStep,  &
                            myLonBegHalo, myLonEndHalo,  &
                            myLatBegHalo, myLatEndHalo)
  call tmg_stop(92)

  !
  !- 5. Main calculation of ensemble and deterministic analyses
  !

  !- 5.1 Do actual analysis calculations in internal subroutine
  call tmg_start(3,'LETKF-doAnalysis')
  call letkf_computeAnalyses
  call tmg_stop(3)

  !- 5.2 Compute ensemble spread stddev Trl and Anl
  call gsv_allocate( stateVectorStdDevTrl, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeight_opt=.false., allocPressure_opt=.false. )
  call gsv_allocate( stateVectorStdDevAnl, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeight_opt=.false., allocPressure_opt=.false. )
  call ens_computeMean(ensembleTrl)
  call ens_computeMean(ensembleAnl)
  call ens_computeStdDev(ensembleTrl)
  call ens_computeStdDev(ensembleAnl)
  call ens_copyEnsStdDev(ensembleTrl, stateVectorStdDevTrl)
  call ens_copyEnsStdDev(ensembleAnl, stateVectorStdDevAnl)

  !
  !- 6. Post processing of analysis results
  !

  !- 6.1 Apply RTPS, if requested
  if (alphaRTPS > 0.0D0) then
    call enkf_RTPS(ensembleAnl, ensembleTrl, stateVectorStdDevAnl, stateVectorStdDevTrl, stateVectorMeanAnl, alphaRTPS)
  end if

  !- 6.2 Apply random additive inflation, if requested
  if (alphaRandomPert > 0.0D0) then
    ! If namelist value is -999, set random seed using the date (as in standard EnKF)
    if (randomSeed == -999) then
      imode = -3 ! stamp to printable date and time: YYYYMMDD, HHMMSShh
      ierr = newdate(dateStamp, datePrint, timePrint, imode)
      timePrint = timePrint/1000000
      datePrint =  datePrint*100 + timePrint
      ! remove the year and add 9
      randomSeed = 9 + datePrint - 1000000*(datePrint/1000000)
      write(*,*) 'midas-letkf: randomSeed for additive inflation set to ', randomSeed
    end if
    call tmg_start(101,'LETKF-randomPert')
    call enkf_addRandomPert(ensembleAnl, stateVectorMeanTrl, alphaRandomPert, randomSeed)
    call tmg_stop(101)
  end if

  !- 6.3 Impose limits on humidity, if requested
  if (imposeSaturationLimit .or. imposeRttovHuLimits) then 
    call tmg_start(102,'LETKF-imposeHulimits')
    ! impose limits on analysis ensemble
    call ens_clipHumidity(ensembleAnl, imposeSaturationLimit, imposeRttovHuLimits)
    ! and recompute analysis mean
    call ens_computeMean(ensembleAnl)
    call ens_copyEnsMean(ensembleAnl, stateVectorMeanAnl)
    ! and recompute mean increment
    call gsv_copy(stateVectorMeanAnl, stateVectorMeanInc)
    call gsv_add(stateVectorMeanTrl, stateVectorMeanInc, scaleFactor_opt=-1.0D0)

    ! impose limits on deterministic analysis
    if (deterExists) then
      if ( imposeSaturationLimit ) call qlim_gsvSaturationLimit(stateVectorDeterAnl)
      if ( imposeRttovHuLimits   ) call qlim_gsvRttovLimit(stateVectorDeterAnl)
      ! and recompute deterministic increment
      call gsv_copy(stateVectorDeterAnl, stateVectorDeterInc)
      call gsv_add(stateVectorDeterTrl, stateVectorDeterInc, scaleFactor_opt=-1.0D0)
    end if
    call tmg_stop(102)
  end if

  !- 6.4 Recompute the analysis spread stddev after inflation and humidity limits
  if (alphaRTPS > 0.0D0 .or. alphaRandomPert > 0.0D0) then
    call ens_computeStdDev(ensembleAnl)
    call ens_copyEnsStdDev(ensembleAnl, stateVectorStdDevAnl)
  end if

  !
  !- 7. Output everything
  !
  call tmg_start(4,'LETKF-writeOutput')

  !- 7.0 Prepare stateVector with only MeanAnl surface pressure and surface height
  call gsv_allocate( stateVectorMeanAnlSfcPressure, tim_nstepobsinc, hco_ens, vco_ens,   &
                     dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0'/) )
  call gsv_zero(stateVectorMeanAnlSfcPressure)
  call gsv_copy(stateVectorMeanAnl, stateVectorMeanAnlSfcPressure, allowMismatch_opt=.true.)
  call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorMeanAnlSfcPressure)

  !- 7.1 Output the ensemble spread stddev Trl and Anl
  call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp() )
  outFileName = trim(outFileName) // '_stddev'
  do stepIndex = 1, tim_nstepobsinc
    call gsv_writeToFile(stateVectorStdDevTrl, outFileName, 'ENSSTD_TRL',  &
                         typvar_opt='R', writeHeightSfc_opt=.false., &
                         stepIndex_opt=stepIndex, containsFullField_opt=.false.)
  end do
  do stepIndex = 1, tim_nstepobsinc
    call gsv_writeToFile(stateVectorStdDevAnl, outFileName, 'ENSSTD_ANL',  &
                         typvar_opt='R', writeHeightSfc_opt=.false., &
                         stepIndex_opt=stepIndex, containsFullField_opt=.false.)
  end do

  !- 7.2 Output the ensemble mean increment (include MeanAnl Psfc)
  if (writeIncrements) then
    call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp(), 0, ensFileNameSuffix_opt='inc' )
    do stepIndex = 1, tim_nstepobsinc
      call gsv_writeToFile(stateVectorMeanInc, outFileName, 'ENSMEAN_INC',  &
                           typvar_opt='R', writeHeightSfc_opt=.false., &
                           stepIndex_opt=stepIndex, containsFullField_opt=.false.)
      call gsv_writeToFile(stateVectorMeanAnlSfcPressure, outFileName, 'ENSMEAN_INC',  &
                           typvar_opt='A', writeHeightSfc_opt=.true., &
                           stepIndex_opt=stepIndex, containsFullField_opt=.true.)
    end do
  end if

  !- 7.3 Output the ensemble mean analysis state
  call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp(), 0 )
  !if (gsv_varExist(stateVectorMeanAnl,'LPR')) call gvt_transform(stateVectorMeanAnl, 'LPRtoPR')
  do stepIndex = 1, tim_nstepobsinc
    call gsv_writeToFile(stateVectorMeanAnl, outFileName, 'ENSMEAN_ANL',  &
                         typvar_opt='A', writeHeightSfc_opt=.false., &
                         stepIndex_opt=stepIndex, containsFullField_opt=.true.)
  end do

  if (deterExists) then
    !- 7.4 Output the deterministic increment (include MeanAnl Psfc)
    if (writeIncrements) then
      call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp(), ensFileNameSuffix_opt='inc' )
      do stepIndex = 1, tim_nstepobsinc
        call gsv_writeToFile(stateVectorDeterInc, outFileName, 'DETER_INC',  &
                             typvar_opt='R', writeHeightSfc_opt=.false., &
                             stepIndex_opt=stepIndex, containsFullField_opt=.false.)
        call gsv_writeToFile(stateVectorMeanAnlSfcPressure, outFileName, 'DETER_INC',  &
                             typvar_opt='A', writeHeightSfc_opt=.true., &
                             stepIndex_opt=stepIndex, containsFullField_opt=.true.)
      end do
    end if

    !- 7.5 Output the deterministic analysis state
    call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp() )
    !if (gsv_varExist(stateVectorDeterAnl,'LPR')) call gvt_transform(stateVectorDeterAnl, 'LPRtoPR')
    do stepIndex = 1, tim_nstepobsinc
      call gsv_writeToFile(stateVectorDeterAnl, outFileName, 'DETER_ANL',  &
                           typvar_opt='A', writeHeightSfc_opt=.false., &
                           stepIndex_opt=stepIndex, containsFullField_opt=.true.)
    end do
  end if

  !- 7.6 Output all ensemble member analyses
  if (updateMembers) then
    !if (gsv_varExist(varName='LPR')) call gvt_transform(ensembleAnl, 'LPRtoPR')
    call ens_writeEnsemble(ensembleAnl, '.', '', ' ', 'ENS_ANL', 'A',  &
                           numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                           containsFullField_opt=.true.)
  end if

  !- 7.7 Output all ensemble member increments (include MeanAnl Psfc)
  if (writeIncrements) then
    ! WARNING: Increment put in ensembleTrl for output
    call ens_add(ensembleAnl, ensembleTrl, scaleFactorInOut_opt=-1.0D0)
    call ens_writeEnsemble(ensembleTrl, '.', '', ' ', 'ENS_INC', 'R',  &
                           numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                           containsFullField_opt=.false.)
    ! Also write the reference (analysis) surface pressure to increment files
    do memberIndex = 1, nEns
      call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp(),  &
                               memberIndex, ensFileNameSuffix_opt='inc' )
      write(memberIndexStr,'(I4.4)') memberIndex
      do stepIndex = 1, tim_nstepobsinc
        call gsv_writeToFile(stateVectorMeanAnlSfcPressure, outFileName,  &
                             'ENS_INC' // memberIndexStr,  &
                             typvar_opt='A', writeHeightSfc_opt=.true., &
                             stepIndex_opt=stepIndex, containsFullField_opt=.true.)
      end do
    end do
  end if

  !- 7.8 Output obs files with mean OMP and OMA

  ! Compute Y-H(Xa_mean) in OBS_OMA
  write(*,*) ''
  write(*,*) 'midas-letkf: apply nonlinear H to ensemble mean analysis'
  write(*,*) ''
  if (tim_nstepobsinc < tim_nstepobs) then
    ! meanAnl is only 3D, so need to make 4D for s2c_nl
    middleStepIndex = (tim_nstepobs + 1) / 2
    call gsv_copy(stateVectorMeanAnl, stateVectorWithZandP4D, allowMismatch_opt=.true., stepIndexOut_opt=middleStepIndex)
    call gsv_3dto4d(stateVectorWithZandP4D)
  else
    call gsv_copy(stateVectorMeanAnl, stateVectorWithZandP4D, allowMismatch_opt=.true.)
  end if
  call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorWithZandP4D)
  call s2c_nl( stateVectorWithZandP4D, obsSpaceData, column, timeInterpType=obsTimeInterpType )
  call tmg_start(6,'LETKF-obsOperators')
  call inn_computeInnovation(column, obsSpaceData, destObsColumn_opt=OBS_OMA, beSilent_opt=.false.)
  call tmg_stop(6)

  ! Write (update) observation files
  call obsf_writeFiles( obsSpaceData )

  !
  !- 8. MPI, tmg finalize
  !  
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  call tmg_stop(4)
  call tmg_stop(1)

  call tmg_terminate(mpi_myid, 'TMG_LETKF' )
  call rpn_comm_finalize(ierr) 

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  if ( mpi_myid == 0 ) write(*,*) ' --------------------------------'
  if ( mpi_myid == 0 ) write(*,*) ' MIDAS-LETKF ENDS'
  if ( mpi_myid == 0 ) write(*,*) ' --------------------------------'

contains

subroutine letkf_computeAnalyses
  ! :Purpose: Local subroutine containing the code for computing
  !           the LETKF analyses for all ensemble members, ensemble
  !           mean and (if present) a deterministic state.
  implicit none

  integer :: memberIndex, memberIndex1, memberIndex2, procIndex, procIndexSend
  integer :: latIndex, lonIndex, stepIndex, kIndex, levIndex, levIndex2
  integer :: bodyIndex, localObsIndex, numLocalObs, numLocalObsFound
  integer :: countMaxExceeded, maxCountMaxExceeded, numGridPointWeights
  integer :: myNumLatLonRecv, myNumLatLonSend, numLatLonRecvMax
  integer :: numLatLonTotalUnique, latLonIndex
  integer :: sendTag, recvTag, nsize, numRecv, numSend
  real(8) :: anlLat, anlLon, anlLogPres, distance
  real(8) :: localization
  real(4) :: secondsBeg, secondsEnd, secondsTotal=0.0

  integer, allocatable :: localBodyIndices(:)
  integer, allocatable :: myLatIndexesRecv(:), myLonIndexesRecv(:)
  integer, allocatable :: myLatIndexesSend(:), myLonIndexesSend(:)
  integer, allocatable :: myNumProcIndexesSend(:)
  integer, allocatable :: myProcIndexesRecv(:), myProcIndexesSend(:,:)
  integer, allocatable :: requestIdRecv(:), requestIdSend(:)

  real(8), allocatable :: distances(:)
  real(8), allocatable :: PaInv(:,:), PaSqrt(:,:), Pa(:,:), YbTinvR(:,:)
  real(8), allocatable :: weightsTemp(:)
  real(8), allocatable :: weightsMembers(:,:,:,:), weightsMembersLatLon(:,:,:)
  real(8), allocatable :: weightsMean(:,:,:,:), weightsMeanLatLon(:,:,:)
  real(8), allocatable :: weightsDeter(:,:,:,:), weightsDeterLatLon(:,:,:)
  real(8), allocatable :: memberAnlPert(:)
  real(4), allocatable :: allSecondsTotal(:,:)

  real(4), pointer     :: meanTrl_ptr_r4(:,:,:,:), meanAnl_ptr_r4(:,:,:,:), meanInc_ptr_r4(:,:,:,:)
  real(4), pointer     :: deterTrl_ptr_r4(:,:,:,:), deterAnl_ptr_r4(:,:,:,:), deterInc_ptr_r4(:,:,:,:)
  real(4), pointer     :: memberTrl_ptr_r4(:,:,:,:), memberAnl_ptr_r4(:,:,:,:)

  logical :: firstTime = .true.

  !
  ! Set things up for the redistribution of work across mpi tasks
  !
  call setupMpiDistribution(myNumLatLonRecv, myNumLatLonSend,   &
                            myLatIndexesRecv, myLonIndexesRecv,   &
                            myLatIndexesSend, myLonIndexesSend,   &
                            myProcIndexesRecv, myProcIndexesSend, &
                            myNumProcIndexesSend)
  allocate(requestIdSend(3*myNumLatLonSend*maxval(myNumProcIndexesSend)))
  allocate(requestIdRecv(3*myNumLatLonRecv))

  !
  ! Compute gridded 3D ensemble weights
  !
  allocate(localBodyIndices(maxNumLocalObs))
  allocate(distances(maxNumLocalObs))
  allocate(YbTinvR(nEns,maxNumLocalObs))
  allocate(PaInv(nEns,nEns))
  allocate(PaSqrt(nEns,nEns))
  allocate(Pa(nEns,nEns))
  allocate(memberAnlPert(nEns))
  allocate(weightsTemp(nEns))
  weightsTemp(:) = 0.0d0
  ! Weights for mean analysis
  allocate(weightsMean(nEns,1,myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
  weightsMean(:,:,:,:) = 0.0d0
  allocate(weightsMeanLatLon(nEns,1,myNumLatLonSend))
  weightsMeanLatLon(:,:,:) = 0.0d0
  ! Weights for member analyses
  allocate(weightsMembers(nEns,nEns,myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
  weightsMembers(:,:,:,:) = 0.0d0
  allocate(weightsMembersLatLon(nEns,nEns,myNumLatLonSend))
  weightsMembersLatLon(:,:,:) = 0.0d0
  ! Weights for deterministic analysis
  if (deterExists) then
    allocate(weightsDeter(nEns,1,myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
    weightsDeter(:,:,:,:) = 0.0d0
    allocate(weightsDeterLatLon(nEns,1,myNumLatLonSend))
    weightsDeterLatLon(:,:,:) = 0.0d0
  end if
  call lfn_Setup(LocFunctionWanted='FifthOrder')

  ! Compute the weights for ensemble mean and members
  countMaxExceeded = 0
  maxCountMaxExceeded = 0
  numGridPointWeights = 0
  LEV_LOOP: do levIndex = 1, nLev_M
    write(*,*) 'computing ensemble updates for vertical level = ', levIndex

    !
    ! First post all recv instructions for communication of weights
    !
    call tmg_start(103,'LETKF-commWeights')
    numSend = 0
    numRecv = 0
    do latLonIndex = 1, myNumLatLonRecv
      latIndex = myLatIndexesRecv(latLonIndex)
      lonIndex = myLonIndexesRecv(latLonIndex)
      procIndex = myProcIndexesRecv(latLonIndex)
      recvTag = (latIndex-1)*stateVectorMeanInc%ni + lonIndex

      nsize = nEns
      numRecv = numRecv + 1
      call mpi_irecv( weightsMean(:,1,lonIndex,latIndex),  &
                      nsize, mpi_datyp_real8, procIndex-1, recvTag,  &
                      mpi_comm_grid, requestIdRecv(numRecv), ierr )
      if (deterExists) then
        numRecv = numRecv + 1
        recvTag = recvTag + stateVectorMeanInc%ni*stateVectorMeanInc%nj
        call mpi_irecv( weightsDeter(:,1,lonIndex,latIndex),  &
                        nsize, mpi_datyp_real8, procIndex-1, recvTag,  &
                        mpi_comm_grid, requestIdRecv(numRecv), ierr )
      end if
      if (updateMembers) then
        nsize = nEns*nEns
        numRecv = numRecv + 1
        recvTag = recvTag + stateVectorMeanInc%ni*stateVectorMeanInc%nj
        call mpi_irecv( weightsMembers(:,:,lonIndex,latIndex),  &
                        nsize, mpi_datyp_real8, procIndex-1, recvTag,  &
                        mpi_comm_grid, requestIdRecv(numRecv), ierr )
      end if
    end do
    call tmg_stop(103)

    LATLON_LOOP: do latLonIndex = 1, myNumLatLonSend
      latIndex = myLatIndexesSend(latLonIndex)
      lonIndex = myLonIndexesSend(latLonIndex)

      numGridPointWeights = numGridPointWeights + 1

      ! lat-lon of the grid point for which we are doing the analysis
      anlLat = hco_ens%lat2d_4(lonIndex,latIndex)
      anlLon = hco_ens%lon2d_4(lonIndex,latIndex)
      anlLogPres = logPres_M_r4(lonIndex,latIndex,levIndex)

      ! Get list of nearby observations and distances to gridpoint
      call tmg_start(9,'LETKF-getLocalBodyIndices')
      numLocalObs = eob_getLocalBodyIndices(ensObs_mpiglobal, localBodyIndices,     &
                                            distances, anlLat, anlLon, anlLogPres,  &
                                            hLocalize, vLocalize, numLocalObsFound)
      if (numLocalObsFound > maxNumLocalObs) then
        countMaxExceeded = countMaxExceeded + 1
        maxCountMaxExceeded = max(maxCountMaxExceeded, numLocalObsFound)
      end if
      call tmg_stop(9)

      call tmg_start(91,'LETKF-calcWeights')
      call cpu_time(secondsBeg)

      ! Extract initial quantities YbTinvR and first term of PaInv (YbTinvR*Yb)
      do localObsIndex = 1, numLocalObs
        bodyIndex = localBodyIndices(localObsIndex)

        ! Compute value of localization function
        call tmg_start(18,'LETKF-locFunction')
        ! Horizontal
        localization = lfn_Response(distances(localObsIndex),hLocalize)
        ! Vertical - use pressures at the grid point (not obs) location
        if (vLocalize > 0) then
          distance = abs( anlLogPres - ensObs_mpiglobal%logPres(bodyIndex) )
          localization = localization * lfn_Response(distance,vLocalize)
        end if
        call tmg_stop(18)

        do memberIndex = 1, nEns
          YbTinvR(memberIndex,localObsIndex) = ensObs_mpiglobal%Yb(memberIndex, bodyIndex) * &
                                               localization * ensObs_mpiglobal%varObsInv(bodyIndex)
        end do

        if (localObsIndex == 1) PaInv(:,:) = 0.0D0
        do memberIndex2 = 1, nEns
          do memberIndex1 = 1, nEns
            PaInv(memberIndex1,memberIndex2) = PaInv(memberIndex1,memberIndex2) +  &
                 YbTinvR(memberIndex1,localObsIndex) * ensObs_mpiglobal%Yb(memberIndex2, bodyIndex)
          end do
        end do

      end do ! localObsIndex

      ! Rest of the computation of local weights for this grid point
      if (numLocalObs > 0) then

        ! Add second term of PaInv
        do memberIndex = 1, nEns
          PaInv(memberIndex,memberIndex) = PaInv(memberIndex,memberIndex) + real(nEns - 1,8)
        end do

        ! Compute Pa and sqrt(Pa) matrices from PaInv
        Pa(:,:) = PaInv(:,:)
        call tmg_start(90,'LETKF-matInverse')
        if (updateMembers) then
          call utl_matInverse(Pa, nEns, inverseSqrt_opt=PaSqrt)
        else
          call utl_matInverse(Pa, nEns)
        end if
        call tmg_stop(90)

        ! Compute ensemble mean local weights as Pa * YbTinvR * (obs - meanYb)
        weightsTemp(:) = 0.0d0
        do localObsIndex = 1, numLocalObs
          bodyIndex = localBodyIndices(localObsIndex)
          do memberIndex = 1, nEns
            weightsTemp(memberIndex) = weightsTemp(memberIndex) +   &
                                       YbTinvR(memberIndex,localObsIndex) *  &
                                       ( ensObs_mpiglobal%obsValue(bodyIndex) - &
                                         ensObs_mpiglobal%meanYb(bodyIndex) )
          end do
        end do

        weightsMeanLatLon(:,1,latLonIndex) = 0.0d0
        do memberIndex2 = 1, nEns
          do memberIndex1 = 1, nEns
            weightsMeanLatLon(memberIndex1,1,latLonIndex) = weightsMeanLatLon(memberIndex1,1,latLonIndex) +  &
                 Pa(memberIndex1,memberIndex2)*weightsTemp(memberIndex2)
          end do
        end do

        if (deterExists) then
          ! Compute deterministic analysis local weights as Pa * YbTinvR * (obs - deterYb)
          weightsTemp(:) = 0.0d0
          do localObsIndex = 1, numLocalObs
            bodyIndex = localBodyIndices(localObsIndex)
            do memberIndex = 1, nEns
              weightsTemp(memberIndex) = weightsTemp(memberIndex) +   &
                                         YbTinvR(memberIndex,localObsIndex) *  &
                                         ( ensObs_mpiglobal%obsValue(bodyIndex) - &
                                           ensObs_mpiglobal%deterYb(bodyIndex) )
            end do
          end do

          weightsDeterLatLon(:,1,latLonIndex) = 0.0d0
          do memberIndex2 = 1, nEns
            do memberIndex1 = 1, nEns
              weightsDeterLatLon(memberIndex1,1,latLonIndex) = weightsDeterLatLon(memberIndex1,1,latLonIndex) +  &
                   Pa(memberIndex1,memberIndex2)*weightsTemp(memberIndex2)
            end do
          end do
        end if

        ! Compute ensemble perturbation weights: (1-alphaRTPP)*[(Nens-1)^1/2*PaSqrt]+alphaRTPP*I
        if (updateMembers) then
          weightsMembersLatLon(:,:,latLonIndex) = (1.0d0 - alphaRTPP) * sqrt(real(nEns - 1,8)) * PaSqrt(:,:)
          do memberIndex = 1, nEns
            weightsMembersLatLon(memberIndex,memberIndex,latLonIndex) = alphaRTPP +  &
                 weightsMembersLatLon(memberIndex,memberIndex,latLonIndex)
          end do
        end if

      else
        ! no observations near this grid point, set weights to zero
        weightsMeanLatLon(:,1,latLonIndex) = 0.0d0
        if (deterExists) weightsDeterLatLon(:,1,latLonIndex) = 0.0d0
      end if ! numLocalObs > 0

      call cpu_time(secondsEnd)
      secondsTotal = secondsTotal + (secondsEnd-secondsBeg)
      call tmg_stop(91)

      !
      ! Now post all send instructions (each lat-lon may be sent to multiple tasks)
      !
      call tmg_start(103,'LETKF-commWeights')
      latIndex = myLatIndexesSend(latLonIndex)
      lonIndex = myLonIndexesSend(latLonIndex)
      do procIndex = 1, myNumProcIndexesSend(latLonIndex)
        sendTag = (latIndex-1)*stateVectorMeanInc%ni + lonIndex
        procIndexSend = myProcIndexesSend(latLonIndex, procIndex)

        nsize = nEns
        numSend = numSend + 1
        call mpi_isend( weightsMeanLatLon(:,1,latLonIndex),  &
                        nsize, mpi_datyp_real8, procIndexSend-1, sendTag,  &
                        mpi_comm_grid, requestIdSend(numSend), ierr )
        if (deterExists) then
          numSend = numSend + 1
          sendTag = sendTag + stateVectorMeanInc%ni*stateVectorMeanInc%nj
          call mpi_isend( weightsDeterLatLon(:,1,latLonIndex),  &
                          nsize, mpi_datyp_real8, procIndexSend-1, sendTag,  &
                          mpi_comm_grid, requestIdSend(numSend), ierr )
        end if
        if (updateMembers) then
          nsize = nEns*nEns
          numSend = numSend + 1
          sendTag = sendTag + stateVectorMeanInc%ni*stateVectorMeanInc%nj
          call mpi_isend( weightsMembersLatLon(:,:,latLonIndex),  &
                          nsize, mpi_datyp_real8, procIndexSend-1, sendTag,  &
                          mpi_comm_grid, requestIdSend(numSend), ierr )
        end if
      end do
      call tmg_stop(103)

    end do LATLON_LOOP

    !
    ! Wait for communiations to finish before continuing
    !
    call tmg_start(103,'LETKF-commWeights')
    if (firstTime) write(*,*) 'numSend/Recv = ', numSend, numRecv
    firstTime = .false.

    if ( numRecv > 0 ) then
      call mpi_waitAll(numRecv, requestIdRecv(1:numRecv), MPI_STATUSES_IGNORE, ierr)
    end if

    if ( numSend > 0 ) then
      call mpi_waitAll(numSend, requestIdSend(1:numSend), MPI_STATUSES_IGNORE, ierr)
    end if

    call tmg_stop(103)

    !
    ! Interpolate weights from coarse to full resolution
    !
    call tmg_start(92,'LETKF-interpolateWeights')
    if (weightLatLonStep > 1) then

      call enkf_interpWeights(wInterpInfo, weightsMean)

      if (deterExists) call enkf_interpWeights(wInterpInfo, weightsDeter)

      if (updateMembers) call enkf_interpWeights(wInterpInfo, weightsMembers)

    end if
    call tmg_stop(92)

    call tmg_start(100,'LETKF-applyWeights')

    !
    ! Apply the weights to compute the ensemble mean and members
    !
    meanInc_ptr_r4 => gsv_getField_r4(stateVectorMeanInc)
    meanTrl_ptr_r4 => gsv_getField_r4(stateVectorMeanTrl)
    meanAnl_ptr_r4 => gsv_getField_r4(stateVectorMeanAnl)
    if (deterExists) then
      deterInc_ptr_r4 => gsv_getField_r4(stateVectorDeterInc)
      deterTrl_ptr_r4 => gsv_getField_r4(stateVectorDeterTrl)
      deterAnl_ptr_r4 => gsv_getField_r4(stateVectorDeterAnl)
    end if
    do latIndex = myLatBeg, myLatEnd
      LON_LOOP5: do lonIndex = myLonBeg, myLonEnd

        ! skip this grid point if all weights zero (no nearby obs)
        if (all(weightsMean(:,1,lonIndex,latIndex) == 0.0d0)) cycle LON_LOOP5

        ! Compute the ensemble mean increment and analysis
        do kIndex = 1, numK
          ! Only treat kIndex values that correspond with current levIndex
          if (vnl_varLevelFromVarname(gsv_getVarNameFromK(stateVectorMeanInc,kIndex)) == 'SF') then
            levIndex2 = nLev_M
          else
            levIndex2 = gsv_getLevFromK(stateVectorMeanInc,kIndex)
          end if
          if (levIndex2 /= levIndex) cycle
          memberTrl_ptr_r4 => ens_getOneLev_r4(ensembleTrl,kIndex)
          do stepIndex = 1, tim_nstepobsinc
            ! mean increment
            do memberIndex = 1, nEns
              meanInc_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) = meanInc_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) +  &
                   weightsMean(memberIndex,1,lonIndex,latIndex) *  &
                   (memberTrl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) - meanTrl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex))
            end do ! memberIndex
            ! mean analysis
            meanAnl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) = meanTrl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) +  &
                 meanInc_ptr_r4(lonIndex,latIndex,kIndex,stepIndex)
          end do ! stepIndex
        end do ! kIndex

        ! Compute the deterministic increment and analysis
        if (deterExists) then
          do kIndex = 1, numK
            ! Only treat kIndex values that correspond with current levIndex
            if (vnl_varLevelFromVarname(gsv_getVarNameFromK(stateVectorMeanInc,kIndex)) == 'SF') then
              levIndex2 = nLev_M
            else
              levIndex2 = gsv_getLevFromK(stateVectorMeanInc,kIndex)
            end if
            if (levIndex2 /= levIndex) cycle
            memberTrl_ptr_r4 => ens_getOneLev_r4(ensembleTrl,kIndex)
            do stepIndex = 1, tim_nstepobsinc
              ! deterministic increment
              do memberIndex = 1, nEns
                deterInc_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) = deterInc_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) +  &
                     weightsDeter(memberIndex,1,lonIndex,latIndex) *  &
                     (memberTrl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) - deterTrl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex))
              end do ! memberIndex
              ! deterministic analysis
              deterAnl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) = deterTrl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) +  &
                   deterInc_ptr_r4(lonIndex,latIndex,kIndex,stepIndex)
            end do ! stepIndex
          end do ! kIndex
        end if

        ! Compute the ensemble member analyses
        if (updateMembers) then
          do kIndex = 1, numK
            ! Only treat kIndex values that correspond with current levIndex
            if (vnl_varLevelFromVarname(gsv_getVarNameFromK(stateVectorMeanInc,kIndex)) == 'SF') then
              levIndex2 = nLev_M
            else
              levIndex2 = gsv_getLevFromK(stateVectorMeanInc,kIndex)
            end if
            if (levIndex2 /= levIndex) cycle
            memberTrl_ptr_r4 => ens_getOneLev_r4(ensembleTrl,kIndex)
            memberAnl_ptr_r4 => ens_getOneLev_r4(ensembleAnl,kIndex)
            do stepIndex = 1, tim_nstepobsinc
              ! Compute analysis member perturbation
              memberAnlPert(:) = 0.0d0
              do memberIndex2 = 1, nEns
                do memberIndex1 = 1, nEns
                  memberAnlPert(memberIndex2) = memberAnlPert(memberIndex2) + &
                       weightsMembers(memberIndex1,memberIndex2,lonIndex,latIndex) *  &
                       (memberTrl_ptr_r4(memberIndex1,stepIndex,lonIndex,latIndex) -  &
                        meanTrl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex))
                end do ! memberIndex1
              end do ! memberIndex2
              ! Add analysis member perturbation to mean analysis
              memberAnl_ptr_r4(:,stepIndex,lonIndex,latIndex) =  &
                   meanAnl_ptr_r4(lonIndex,latIndex,kIndex,stepIndex) + memberAnlPert(:)
            end do ! stepIndex
          end do ! kIndex
        end if ! updateMembers

      end do LON_LOOP5
    end do

    call tmg_stop(100)

  end do LEV_LOOP

  if (countMaxExceeded > 0) then
    write(*,*) 'midas-letkf: WARNING: Found more local obs than specified max number at ', &
               real(100*countMaxExceeded)/real(numGridPointWeights), '% of grid points.'
    write(*,*) '                      Maximum number found was ', maxCountMaxExceeded,  &
               ' which is greater than specified number ', maxNumLocalObs
    write(*,*) '                      Therefore will keep closest obs only.'
  end if

  call tmg_start(19,'LETKF-barr')
  call rpn_comm_barrier('GRID',ierr)
  call tmg_stop(19)

  ! Write out an ascii file with the order of timing of weight calculation
  allocate(allSecondsTotal(mpi_nprocs,2))
  call rpn_comm_gather( secondsTotal, 1, 'mpi_real4', allSecondsTotal(:,1), 1, 'mpi_real4', &
                        0, 'GRID', ierr )
  if (mpi_myid == 0) then
    do procIndex = 1, mpi_nprocs
      allSecondsTotal(procIndex,2) = real(procIndex)
    end do
    call utl_heapsort2d(allSecondsTotal)
    nulFile = 0
    ierr = fnom(nulFile, './procOrder.txt', 'FTN+SEQ+R/W', 0)
    do procIndex = 1, mpi_nprocs
      write(nulFile,*) allSecondsTotal(procIndex,:)
    end do
    ierr = fclos(nulFile)
  end if

end subroutine letkf_computeAnalyses


subroutine setupMpiDistribution(myNumLatLonRecv, myNumLatLonSend, &
                                myLatIndexesRecv, myLonIndexesRecv, &
                                myLatIndexesSend, myLonIndexesSend, &
                                myProcIndexesRecv, myProcIndexesSend, &
                                myNumProcIndexesSend)

  ! :Purpose: Setup for distribution of grid points over mpi tasks.

  implicit none

  ! Arguments
  integer              :: myNumLatLonRecv, myNumLatLonSend
  integer, allocatable :: myLatIndexesRecv(:), myLonIndexesRecv(:)
  integer, allocatable :: myLatIndexesSend(:), myLonIndexesSend(:)
  integer, allocatable :: myProcIndexesRecv(:), myProcIndexesSend(:,:)
  integer, allocatable :: myNumProcIndexesSend(:)

  ! Locals
  integer :: latIndex, lonIndex, procIndex, procIndexSend, latLonIndex
  integer :: numLatLonRecvMax, numLatLonTotalUnique
  integer, allocatable :: allLatIndexesRecv(:,:), allLonIndexesRecv(:,:)
  integer, allocatable :: allLatIndexesSend(:,:), allLonIndexesSend(:,:)
  integer, allocatable :: allNumLatLonRecv(:), allNumLatLonSend(:)

  if (trim(mpiDistribution) == 'TILES') then

    ! First, determine number of grid points needed locally (for recv-ing)
    myNumLatLonRecv = 0
    do latIndex = myLatBegHalo, myLatEndHalo
      LON_LOOP0: do lonIndex = myLonBegHalo, myLonEndHalo
        ! If this lat-lon is to be interpolated, then skip calculation
        if (wInterpInfo%numIndexes(lonIndex,latIndex) > 0) cycle LON_LOOP0
        myNumLatLonRecv = myNumLatLonRecv + 1
      end do LON_LOOP0
    end do

    write(*,*) 'midas-letkf: myNumLatLonRecv =', myNumLatLonRecv

    ! Determine list of grid point indexes where weights needed locally (for recv-ing)
    allocate(myLatIndexesRecv(myNumLatLonRecv))
    allocate(myLonIndexesRecv(myNumLatLonRecv))
    allocate(myProcIndexesRecv(myNumLatLonRecv))
    myNumLatLonRecv = 0
    do latIndex = myLatBegHalo, myLatEndHalo
      LON_LOOP1: do lonIndex = myLonBegHalo, myLonEndHalo
        ! If this lat-lon is to be interpolated, then skip calculation
        if (wInterpInfo%numIndexes(lonIndex,latIndex) > 0) cycle LON_LOOP1
        myNumLatLonRecv = myNumLatLonRecv + 1

        myLatIndexesRecv(myNumLatLonRecv) = latIndex
        myLonIndexesRecv(myNumLatLonRecv) = lonIndex
        myProcIndexesRecv(myNumLatLonRecv) = mpi_myid+1
      end do LON_LOOP1
    end do

    ! No communication, so send info equals recv info
    myNumLatLonSend = myNumLatLonRecv
    allocate(myLatIndexesSend(myNumLatLonSend))
    allocate(myLonIndexesSend(myNumLatLonSend))
    allocate(myProcIndexesSend(myNumLatLonSend,1))
    allocate(myNumProcIndexesSend(myNumLatLonSend))

    myLatIndexesSend(:) = myLatIndexesRecv(:)
    myLonIndexesSend(:) = myLonIndexesRecv(:)
    myProcIndexesSend(:,1) = myProcIndexesRecv(:)
    myNumProcIndexesSend(:) = 1

  else if (trim(mpiDistribution) == 'ROUNDROBIN') then

    ! First, determine number of grid points needed locally (for recv-ing)
    myNumLatLonRecv = 0
    do latIndex = myLatBegHalo, myLatEndHalo
      LON_LOOP2: do lonIndex = myLonBegHalo, myLonEndHalo
        ! If this lat-lon is to be interpolated, then skip calculation
        if (wInterpInfo%numIndexes(lonIndex,latIndex) > 0) cycle LON_LOOP2
        myNumLatLonRecv = myNumLatLonRecv + 1
      end do LON_LOOP2
    end do

    ! Communicate to all mpi tasks
    allocate(allNumLatLonRecv(mpi_nprocs))
    call rpn_comm_allgather(myNumLatLonRecv, 1, "mpi_integer",  &
                            allNumLatLonRecv, 1,"mpi_integer", "GRID", ierr)
    numLatLonRecvMax = maxval(allNumLatLonRecv)
    write(*,*) 'midas-letkf: allNumLatLonRecv =', allNumLatLonRecv(:)
    write(*,*) 'midas-letkf: numLatLonRecvSum =', sum(allNumLatLonRecv)
    write(*,*) 'midas-letkf: numLatLonRecvMax =', numLatLonRecvMax

    ! Determine list of grid point indexes where weights needed locally (for recv-ing)
    allocate(myLatIndexesRecv(numLatLonRecvMax))
    allocate(myLonIndexesRecv(numLatLonRecvMax))
    allocate(myProcIndexesRecv(numLatLonRecvMax))
    myLatIndexesRecv(:) = -1
    myLonIndexesRecv(:) = -1
    myProcIndexesRecv(:) = -1
    myNumLatLonRecv = 0
    do latIndex = myLatBegHalo, myLatEndHalo
      LON_LOOP3: do lonIndex = myLonBegHalo, myLonEndHalo
        ! If this lat-lon is to be interpolated, then skip calculation
        if (wInterpInfo%numIndexes(lonIndex,latIndex) > 0) cycle LON_LOOP3
        myNumLatLonRecv = myNumLatLonRecv + 1

        myLatIndexesRecv(myNumLatLonRecv) = latIndex
        myLonIndexesRecv(myNumLatLonRecv) = lonIndex
      end do LON_LOOP3
    end do

    ! Communicate to all mpi tasks this list of grid point lat-lon indexes
    allocate(allLatIndexesRecv(numLatLonRecvMax, mpi_nprocs))
    allocate(allLonIndexesRecv(numLatLonRecvMax, mpi_nprocs))
    call rpn_comm_allgather(myLatIndexesRecv, numLatLonRecvMax, "mpi_integer",  &
                            allLatIndexesRecv, numLatLonRecvMax, "mpi_integer",  &
                            "GRID", ierr)
    call rpn_comm_allgather(myLonIndexesRecv, numLatLonRecvMax, "mpi_integer",  &
                            allLonIndexesRecv, numLatLonRecvMax, "mpi_integer",  &
                            "GRID", ierr)

    ! From these lat-lon lists, create unique master list of all grid points where weights computed
    ! and assign to mpi tasks for doing the calculation and for send-ing
    allocate(myLatIndexesSend(numLatLonRecvMax))
    allocate(myLonIndexesSend(numLatLonRecvMax))
    myLatIndexesSend(:) = -1
    myLonIndexesSend(:) = -1
    numLatLonTotalUnique = 0
    myNumLatLonSend = 0
    do procIndex = 1, mpi_nprocs
      WEIGHTS1LEV_LOOP: do latLonIndex = 1, allNumLatLonRecv(procIndex)
        if (alreadyFound(allLatIndexesRecv, allLonIndexesRecv, latLonIndex, procIndex)) &
             cycle WEIGHTS1LEV_LOOP
        ! Count the total number of weights
        numLatLonTotalUnique = numLatLonTotalUnique + 1

        ! Round-robin distribution of master list across mpi tasks
        procIndexSend = 1 + mod(numLatLonTotalUnique-1, mpi_nprocs)

        ! Store the lat-lon indexes of the weights I am responsible for
        if (procIndexSend == (mpi_myid+1)) then
          myNumLatLonSend = myNumLatLonSend + 1
          myLatIndexesSend(myNumLatLonSend) =  &
               allLatIndexesRecv(latLonIndex, procIndex)
          myLonIndexesSend(myNumLatLonSend) =  &
               allLonIndexesRecv(latLonIndex, procIndex)
        end if
      end do WEIGHTS1LEV_LOOP
    end do
    write(*,*) 'midas-letkf: number of lat/lon points where weights to be computed =',  &
               numLatLonTotalUnique

    ! Communicate to all mpi tasks this list of grid point lat-lon indexes
    allocate(allNumLatLonSend(mpi_nprocs))
    call rpn_comm_allgather(myNumLatLonSend, 1, "mpi_integer",  &
                            allNumLatLonSend, 1,"mpi_integer", "GRID", ierr)
    allocate(allLatIndexesSend(numLatLonRecvMax, mpi_nprocs))
    allocate(allLonIndexesSend(numLatLonRecvMax, mpi_nprocs))
    call rpn_comm_allgather(myLatIndexesSend, numLatLonRecvMax, "mpi_integer",  &
                            allLatIndexesSend, numLatLonRecvMax, "mpi_integer",  &
                            "GRID", ierr)
    call rpn_comm_allgather(myLonIndexesSend, numLatLonRecvMax, "mpi_integer",  &
                            allLonIndexesSend, numLatLonRecvMax, "mpi_integer",  &
                            "GRID", ierr)

    ! Figure out which mpi tasks I will need to send my results to
    allocate(myProcIndexesSend(myNumLatLonSend,mpi_nprocs))
    allocate(myNumProcIndexesSend(myNumLatLonSend))
    myProcIndexesSend(:,:) = -1
    myNumProcIndexesSend(:) = 0
    do latLonIndex = 1, myNumLatLonSend
      do procIndex = 1, mpi_nprocs
        if ( any( (myLatIndexesSend(latLonIndex) == allLatIndexesRecv(1:allNumLatLonRecv(procIndex), procIndex)) .and.  &
                  (myLonIndexesSend(latLonIndex) == allLonIndexesRecv(1:allNumLatLonRecv(procIndex), procIndex)) ) ) then
          myNumProcIndexesSend(latLonIndex) = myNumProcIndexesSend(latLonIndex) + 1
          myProcIndexesSend(latLonIndex,myNumProcIndexesSend(latLonIndex)) = procIndex
        end if
      end do
    end do

    ! Figure out which mpi tasks I will receive the results from
    do latLonIndex = 1, myNumLatLonRecv
      do procIndex = 1, mpi_nprocs
        if ( any( (myLatIndexesRecv(latLonIndex) == allLatIndexesSend(1:allNumLatLonSend(procIndex), procIndex)) .and.  &
                  (myLonIndexesRecv(latLonIndex) == allLonIndexesSend(1:allNumLatLonSend(procIndex), procIndex)) ) ) then
          myProcIndexesRecv(latLonIndex) = procIndex
        end if
      end do
    end do

  else
    call utl_abort('midas-letkf: unknown MPI distribution selected')
  end if

  write(*,*) 'midas-letkf: lat/lon/proc indexes I need to receive:'
  do latLonIndex = 1, myNumLatLonRecv
    write(*,*) myLatIndexesRecv(latLonIndex), myLonIndexesRecv(latLonIndex),  &
               myProcIndexesRecv(latLonIndex)
  end do

  write(*,*) 'midas-letkf: number of lat/lon indexes I am responsible for =', myNumLatLonSend
  write(*,*) 'midas-letkf: the lat/lon/proc indexes I am responsible for:'
  do latLonIndex = 1, myNumLatLonSend
    write(*,*) myLatIndexesSend(latLonIndex), myLonIndexesSend(latLonIndex),  &
               myProcIndexesSend(latLonIndex,1:myNumProcIndexesSend(latLonIndex))
  end do

end subroutine setupMpiDistribution


function alreadyFound(allLatIndexesRecv, allLonIndexesRecv, latLonIndex, procIndex) result(found)
  implicit none
  ! Arguments:
  integer :: allLatIndexesRecv(:,:), allLonIndexesRecv(:,:)
  integer :: latLonIndex, procIndex
  logical :: found

  ! Locals:
  integer :: latLonIndex2, procIndex2, numLatLonRecvMax

  numLatLonRecvMax = size(allLatIndexesRecv, 1)

  ! check on all previous mpi tasks if this lat/lon has already been encountered
  found = .false.
  do procIndex2 = 1, procIndex-1
    WEIGHTS1LEV_LOOP2: do latLonIndex2 = 1, numLatLonRecvMax
      if (allLatIndexesRecv(latLonIndex2, procIndex2) < 0) cycle WEIGHTS1LEV_LOOP2
      if ( (allLatIndexesRecv(latLonIndex, procIndex) ==  &
            allLatIndexesRecv(latLonIndex2, procIndex2)) .and.  &
           (allLonIndexesRecv(latLonIndex, procIndex) ==  &
            allLonIndexesRecv(latLonIndex2, procIndex2)) ) then
        found = .true.
        exit WEIGHTS1LEV_LOOP2
      end if
    end do WEIGHTS1LEV_LOOP2
  end do

end function alreadyFound


end program midas_letkf

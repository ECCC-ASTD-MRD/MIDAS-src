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
  use mpi_mod
  use mathPhysConstants_mod
  use fileNames_mod
  use ensembleObservations_mod
  use ensembleStateVector_mod
  use gridStateVector_mod
  use columnData_mod
  use tovs_nl_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use analysisGrid_mod
  use timeCoord_mod
  use obsTimeInterp_mod
  use utilities_mod
  use ramDisk_mod
  use statetocolumn_mod
  use obsFiles_mod
  use obsSpaceData_mod
  use obsErrors_mod
  use obsFilter_mod
  use innovation_mod
  use gridVariableTransforms_mod
  use physicsFunctions_mod
  use enkf_mod
  use humidityLimits_mod
  use clib_interfaces_mod
  implicit none

  type(struct_obs), target  :: obsSpaceData
  type(struct_ens), target  :: ensembleTrl4D
  type(struct_ens), pointer :: ensembleTrl
  type(struct_ens)          :: ensembleAnl
  type(struct_ens)          :: ensembleTrlSubSample
  type(struct_ens)          :: ensembleAnlSubSample
  type(struct_gsv)          :: stateVectorMeanTrl4D
  type(struct_gsv)          :: stateVectorMeanTrl
  type(struct_gsv)          :: stateVectorMeanAnl
  type(struct_gsv)          :: stateVectorMeanAnlSubSample
  type(struct_gsv)          :: stateVectorMeanInc
  type(struct_gsv)          :: stateVectorMeanIncSubSample
  type(struct_gsv)          :: stateVectorWithZandP4D
  type(struct_gsv)          :: stateVectorStdDevTrl
  type(struct_gsv)          :: stateVectorStdDevAnl
  type(struct_gsv)          :: stateVectorStdDevAnlPert
  type(struct_gsv)          :: stateVectorMeanAnlSfcPres
  type(struct_gsv)          :: stateVectorMeanAnlSfcPresMpiGlb
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

  integer :: memberIndex, stepIndex, middleStepIndex, randomSeedRandomPert, randomSeedObs
  integer :: nulnam, dateStamp, datePrint, timePrint, imode, ierr
  integer :: get_max_rss, fclos, fnom, fstopc, newdate
  integer, allocatable :: dateStampList(:), dateStampListInc(:)

  character(len=256) :: ensFileName, outFileName, deterFileName
  character(len=9)   :: obsColumnMode
  character(len=48)  :: obsMpiStrategy
  character(len=48)  :: midasMode
  character(len=4)   :: memberIndexStr
  character(len=4), pointer :: varNames(:)
  character(len=12)  :: etiketMean='', etiketStd=''

  logical :: deterministicStateExists

  ! interpolation information for weights (in enkf_mod)
  type(struct_enkfInterpInfo) :: wInterpInfo

  ! namelist variables
  character(len=20)  :: algorithm ! name of the chosen LETKF algorithm: 'LETKF', 'CVLETKF'
  integer            :: numSubEns ! number of sub-ensembles to split the full ensemble
  character(len=256) :: ensPathName ! absolute or relative path to ensemble directory
  integer  :: nEns                ! ensemble size
  integer  :: maxNumLocalObs      ! maximum number of obs in each local volume to assimilate
  integer  :: weightLatLonStep    ! separation of lat-lon grid points for weight calculation
  integer  :: randomSeed          ! seed used for random perturbation additive inflation
  logical  :: modifyAmsubObsError ! reduce AMSU-B obs error stddev in tropics
  logical  :: backgroundCheck     ! apply additional background check using ensemble spread
  logical  :: huberize            ! apply huber norm quality control procedure
  logical  :: rejectHighLatIR     ! reject all IR observations at high latitudes
  logical  :: rejectRadNearSfc    ! reject radiance observations near the surface
  logical  :: writeSubSample      ! write sub-sample members for initializing medium-range fcsts
  real(8)  :: hLocalize           ! horizontal localization radius (in km)
  real(8)  :: vLocalize           ! vertical localization radius (in units of ln(Pressure in Pa))
  real(8)  :: alphaRTPP           ! RTPP coefficient (between 0 and 1; 0 means no relaxation)
  real(8)  :: alphaRTPS           ! RTPS coefficient (between 0 and 1; 0 means no relaxation)
  real(8)  :: alphaRandomPert     ! Random perturbation additive inflation coeff (0->1)
  real(8)  :: alphaRandomPertSubSample ! Random perturbation additive inflation coeff for medium-range fcsts
  logical  :: imposeSaturationLimit ! switch for choosing to impose saturation limit of humidity
  logical  :: imposeRttovHuLimits   ! switch for choosing to impose the RTTOV limits on humidity
  character(len=20) :: obsTimeInterpType ! type of time interpolation to obs time
  character(len=12) :: etiket0
  character(len=20) :: mpiDistribution   ! type of mpiDistribution for weight calculation ('ROUNDROBIN' or 'TILES')
  NAMELIST /NAMLETKF/algorithm, nEns, numSubEns, ensPathName, hLocalize, vLocalize,  &
                     maxNumLocalObs, weightLatLonStep,  &
                     modifyAmsubObsError, backgroundCheck, huberize, rejectHighLatIR, rejectRadNearSfc,  &
                     writeSubSample,  &
                     alphaRTPP, alphaRTPS, randomSeed, alphaRandomPert, alphaRandomPertSubSample,  &
                     imposeSaturationLimit, imposeRttovHuLimits, obsTimeInterpType, &
                     etiket0, mpiDistribution


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
  algorithm             = 'LETKF'
  ensPathName           = 'ensemble'
  nEns                  = 10
  numSubEns             = 2
  maxNumLocalObs        = 1000
  weightLatLonStep      = 1
  modifyAmsubObsError   = .false.
  backgroundCheck       = .false.
  huberize              = .false.
  rejectHighLatIR       = .false.
  rejectRadNearSfc      = .false.
  writeSubSample        = .false.
  hLocalize             = 500.0D0
  vLocalize             = -1.0D0
  alphaRTPP             =  0.0D0
  alphaRTPS             =  0.0D0
  randomSeed            =  -999
  alphaRandomPert       =  0.0D0
  alphaRandomPertSubSample =  -1.0D0
  imposeSaturationLimit = .false.
  imposeRttovHuLimits   = .false.
  obsTimeInterpType     = 'LINEAR'
  etiket0               = 'E26_0_0P'
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
  if (alphaRandomPertSubSample < 0.0D0) alphaRandomPertSubSample = 0.0D0
  if (trim(algorithm) /= 'LETKF' .and. trim(algorithm) /= 'CVLETKF' .and.  &
      trim(algorithm) /= 'CVLETKF-PERTOBS') then
    call utl_abort('midas-letkf: unknown LETKF algorithm: ' // trim(algorithm))
  end if

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
  if (vco_getNumLev(vco_ens, 'MM') /= vco_getNumLev(vco_ens, 'TH')) then
    call utl_abort('midas-letkf: nLev_M /= nLev_T - currently not supported')
  end if

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

  ! Initialize obs error covariances and set flag using 'util' column of stats_tovs
  call oer_setObsErrors(obsSpaceData, midasMode, useTovsUtil_opt=.true.) ! IN

  ! Call suprep again to filter out channels according to 'util' column of stats_tovs
  call filt_suprep(obsSpaceData)

  ! Allocate vectors for storing HX values
  call eob_allocate(ensObs, nEns, obs_numBody(obsSpaceData), obsSpaceData)
  call eob_zero(ensObs)

  ! Set lat, lon, obs values in ensObs
  call eob_setLatLonObs(ensObs)

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
  inquire(file=deterFileName, exist=deterministicStateExists)
  if (deterministicStateExists) then
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
  if (deterministicStateExists) then

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

    ! Copy to ensObs: Y-HX for deterministic background
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

    ! Copy to ensObs: Y-HX for this member
    call eob_setYb(ensObs, memberIndex)

  end do

  !- 3.3 Set some additional information in ensObs and additional quality 
  !      control before finally communicating ensObs globally

  ! Compute and remove the mean of Yb
  call eob_calcAndRemoveMeanYb(ensObs)

  ! Put HPHT in OBS_HPHT, for writing to obs files
  call eob_setHPHT(ensObs)

  ! Compute random observation perturbations
  if (trim(algorithm) == 'CVLETKF-PERTOBS') then
    randomSeedObs = 1 + mpi_myid
    call eob_calcRandPert(ensObs, randomSeedObs)
  end if

  ! Apply obs operators to ensemble mean background for several purposes
  write(*,*) ''
  write(*,*) 'midas-letkf: apply nonlinear H to ensemble mean background'
  write(*,*) ''
  call gsv_copy(stateVectorMeanTrl4D, stateVectorWithZandP4D, allowMismatch_opt=.true.)
  call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorWithZandP4D)
  call s2c_nl( stateVectorWithZandP4D, obsSpaceData, column, timeInterpType=obsTimeInterpType, dealloc_opt=.false. )
  call tvs_allocTransmission ! this will cause radiative transmission profiles to be stored for use in eob_setLogPres
  call tmg_start(6,'LETKF-obsOperators')
  call inn_computeInnovation(column, obsSpaceData, beSilent_opt=.false.)
  call tmg_stop(6)

  ! Put y-mean(H(X)) in OBS_OMP for writing to obs files (overwrites y-H(mean(X)))
  call eob_setMeanOMP(ensObs)

  ! Set pressure for all obs for vertical localization, based on ensemble mean pressure and height
  call eob_setLogPres(ensObs, column)

  ! Modify the obs error stddev for AMSUB in the tropics
  if (modifyAmsubObsError) call enkf_modifyAmsubObsError(obsSpaceData)

  ! Apply a background check (reject limit is set in the routine)
  if (backgroundCheck) call eob_backgroundCheck(ensObs)

  ! Set values of obs_sigi and obs_sigo before hubernorm modifies obs_oer
  call eob_setSigiSigo(ensObs)

  ! Apply huber norm quality control procedure (modifies obs_oer)
  if (huberize) call eob_huberNorm(ensObs)

  !- Reject all IR radiance observation in arctic and antarctic (.i.e |lat|>60. )
  if (rejectHighLatIR) call enkf_rejectHighLatIR(obsSpaceData)

  ! Reject radiance observations too close to the surface
  if (rejectRadNearSfc) call eob_rejectRadNearSfc(ensObs)

  ! Compute inverse of obs error variance (done here to use dynamic GPS-RO, GB-GPS based on mean O-P)
  call eob_setObsErrInv(ensObs)

  ! Clean and globally communicate obs-related data to all mpi tasks
  call eob_allGather(ensObs,ensObs_mpiglobal)

  ! Print number of assimilated obs per family to the listing
  write(*,*) 'oti_timeBinning: After extra filtering done in midas-letkf'
  call oti_timeBinning(obsSpaceData,tim_nstepobs)

  call tmg_stop(2)

  !
  !- 4. Final preparations for computing analyses
  !

  !- 4.1 Copy trial ensemble to nstepobsinc time steps
  if (tim_nstepobsinc < tim_nstepobs) then
    allocate(ensembleTrl)
    call ens_allocate(ensembleTrl, nEns, tim_nstepobsinc, hco_ens, vco_ens, dateStampListInc)
    call ens_copy4Dto3D(ensembleTrl4D,ensembleTrl)
    call ens_deallocate(ensembleTrl4D)
  else
    ! number of timesteps is the same, so just point to it
    ensembleTrl => ensembleTrl4D
  end if

  !- 4.2 Copy trl ensemble to anl ensemble
  call ens_allocate(ensembleAnl, nEns, tim_nstepobsinc, hco_ens, vco_ens, dateStampListInc)
  call ens_copy(ensembleTrl,ensembleAnl)

  !- 4.3 Setup for interpolating weights from coarse to full resolution
  call tmg_start(92,'LETKF-interpolateWeights')
  call enkf_setupInterpInfo(wInterpInfo, stateVectorMeanInc%hco, weightLatLonStep,  &
                            stateVectorMeanInc%myLonBeg, stateVectorMeanInc%myLonEnd,  &
                            stateVectorMeanInc%myLatBeg, stateVectorMeanInc%myLatEnd)
  call tmg_stop(92)

  !
  !- 5. Main calculation of ensemble and deterministic analyses
  !

  !- 5.1 Do actual analysis calculations for chosen flavour of LETKF
  call tmg_start(3,'LETKF-doAnalysis')
  call enkf_LETKFanalyses(algorithm, numSubEns,  &
                          ensembleAnl, ensembleTrl, ensObs_mpiglobal,  &
                          stateVectorMeanInc, stateVectorMeanTrl, stateVectorMeanAnl, &
                          stateVectorDeterInc, stateVectorDeterTrl, stateVectorDeterAnl, &
                          wInterpInfo, maxNumLocalObs,  &
                          hLocalize, vLocalize, alphaRTPP, mpiDistribution)
  call tmg_stop(3)

  !- 5.2 Compute ensemble spread stddev Trl and Anl
  call gsv_allocate( stateVectorStdDevTrl, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeight_opt=.false., allocPressure_opt=.false. )
  call gsv_allocate( stateVectorStdDevAnl, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeight_opt=.false., allocPressure_opt=.false. )
  call gsv_allocate( stateVectorStdDevAnlPert, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
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
    ! recompute the analysis spread stddev
    call ens_computeStdDev(ensembleAnl)
    call ens_copyEnsStdDev(ensembleAnl, stateVectorStdDevAnl)
  end if

  !- 6.2 If SubSample requested, copy sub-sample of analysis and trial members
  if (writeSubSample) then
    ! Copy sub-sampled analysis and trial ensemble members
    call enkf_selectSubSample(ensembleAnl, ensembleTrl,  &
                              ensembleAnlSubSample, ensembleTrlSubSample)

    ! Create subdirectory for outputting sub sample increments
    ierr = clib_mkdir_r('subspace')

    ! Allocate stateVectors to store and output sub-sampled ensemble mean analysis and increment
    call gsv_allocate( stateVectorMeanAnlSubSample, tim_nstepobsinc,  &
                       hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_zero(stateVectorMeanAnlSubSample)
    call gsv_allocate( stateVectorMeanIncSubSample, tim_nstepobsinc,  &
                       hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_zero(stateVectorMeanIncSubSample)

  end if

  !- 6.3 Apply random additive inflation, if requested
  if (alphaRandomPert > 0.0D0) then
    ! If namelist value is -999, set random seed using the date (as in standard EnKF)
    if (randomSeed == -999) then
      imode = -3 ! stamp to printable date and time: YYYYMMDD, HHMMSShh
      ierr = newdate(dateStamp, datePrint, timePrint, imode)
      timePrint = timePrint/1000000
      datePrint =  datePrint*100 + timePrint
      ! Remove the year and add 9
      randomSeedRandomPert = 9 + datePrint - 1000000*(datePrint/1000000)
      write(*,*) 'midas-letkf: randomSeed for additive inflation set to ', randomSeedRandomPert
    else
      randomSeedRandomPert = randomSeed
    end if
    call tmg_start(101,'LETKF-randomPert')
    call enkf_addRandomPert(ensembleAnl, stateVectorMeanTrl, alphaRandomPert, randomSeedRandomPert)
    call tmg_stop(101)
  end if

  !- 6.4 Impose limits on humidity, if requested
  if (imposeSaturationLimit .or. imposeRttovHuLimits) then 
    call tmg_start(102,'LETKF-imposeHulimits')
    ! Impose limits on analysis ensemble
    if (imposeSaturationLimit .or. imposeRttovHuLimits) then
      if (mpi_myid == 0) write(*,*) ''
      if (mpi_myid == 0) write(*,*) 'midas-letkf: limits will be imposed on the humidity of analysis ensemble'
      if (mpi_myid == 0 .and. imposeSaturationLimit ) write(*,*) '              -> Saturation Limit'
      if (mpi_myid == 0 .and. imposeRttovHuLimits   ) write(*,*) '              -> Rttov Limit'

      if ( imposeSaturationLimit ) call qlim_saturationLimit(ensembleAnl)
      if ( imposeRttovHuLimits   ) call qlim_rttovLimit     (ensembleAnl)
    end if

    ! And recompute analysis mean
    call ens_computeMean(ensembleAnl)
    call ens_copyEnsMean(ensembleAnl, stateVectorMeanAnl)
    ! And recompute mean increment
    call gsv_copy(stateVectorMeanAnl, stateVectorMeanInc)
    call gsv_add(stateVectorMeanTrl, stateVectorMeanInc, scaleFactor_opt=-1.0D0)

    ! Impose limits on deterministic analysis
    if (deterministicStateExists) then
      if ( imposeSaturationLimit ) call qlim_saturationLimit(stateVectorDeterAnl)
      if ( imposeRttovHuLimits   ) call qlim_rttovLimit(stateVectorDeterAnl)
      ! And recompute deterministic increment
      call gsv_copy(stateVectorDeterAnl, stateVectorDeterInc)
      call gsv_add(stateVectorDeterTrl, stateVectorDeterInc, scaleFactor_opt=-1.0D0)
    end if
    call tmg_stop(102)
  end if

  !- 6.5 Recompute the analysis spread stddev after inflation and humidity limits
  call ens_computeStdDev(ensembleAnl)
  call ens_copyEnsStdDev(ensembleAnl, stateVectorStdDevAnlPert)

  !- 6.6 If SubSample requested, do remaining processing and output of sub-sampled members
  if (writeSubSample) then

    ! Apply random additive inflation to sub-sampled ensemble, if requested
    if (alphaRandomPertSubSample > 0.0D0) then
      ! If namelist value is -999, set random seed using the date (as in standard EnKF)
      if (randomSeed == -999) then
        imode = -3 ! stamp to printable date and time: YYYYMMDD, HHMMSShh
        ierr = newdate(dateStamp, datePrint, timePrint, imode)
        timePrint = timePrint/1000000
        datePrint =  datePrint*100 + timePrint
        ! Remove the year and add 9
        randomSeedRandomPert = 9 + datePrint - 1000000*(datePrint/1000000)
        write(*,*) 'midas-letkf: randomSeed for additive inflation set to ', randomSeedRandomPert
      else
        randomSeedRandomPert = randomSeed
      end if
      call tmg_start(101,'LETKF-randomPert')
      call enkf_addRandomPert(ensembleAnlSubSample, stateVectorMeanTrl,  &
                              alphaRandomPertSubSample, randomSeedRandomPert)
      call tmg_stop(101)
    end if

    ! Compute analysis mean of sub-sampled ensemble
    call ens_computeMean(ensembleAnlSubSample)

    ! Shift members to have same mean as full ensemble and impose humidity limits, if requested
    call ens_recenter(ensembleAnlSubSample, stateVectorMeanAnl, recenteringCoeff=1.0D0)
    if (imposeSaturationLimit .or. imposeRttovHuLimits) then
      if (mpi_myid == 0) write(*,*) ''
      if (mpi_myid == 0) write(*,*) 'midas-letkf: limits will be imposed on the humidity of recentered ensemble'
      if (mpi_myid == 0 .and. imposeSaturationLimit ) write(*,*) '              -> Saturation Limit'
      if (mpi_myid == 0 .and. imposeRttovHuLimits   ) write(*,*) '              -> Rttov Limit'

      if ( imposeSaturationLimit ) call qlim_saturationLimit(ensembleAnlSubSample)
      if ( imposeRttovHuLimits   ) call qlim_rttovLimit     (ensembleAnlSubSample)
    end if

    ! Re-compute analysis mean of sub-sampled ensemble
    call ens_computeMean(ensembleAnlSubSample)
    call ens_copyEnsMean(ensembleAnlSubSample, stateVectorMeanAnlSubSample)

    ! And compute mean increment with respect to mean of full trial ensemble
    call gsv_copy(stateVectorMeanAnlSubSample, stateVectorMeanIncSubSample)
    call gsv_add(stateVectorMeanTrl, stateVectorMeanIncSubSample, scaleFactor_opt=-1.0D0)

  end if

  !
  !- 7. Output everything
  !
  call tmg_start(4,'LETKF-writeOutput')

  !- 7.0 Output obs files with mean OMP and OMA

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

  !- 7.1 Prepare stateVector with only MeanAnl surface pressure and surface height
  call gsv_allocate( stateVectorMeanAnlSfcPres, tim_nstepobsinc, hco_ens, vco_ens,   &
                     dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0'/) )
  call gsv_zero(stateVectorMeanAnlSfcPres)
  if (mpi_myid <= (nEns-1)) then
    call gsv_allocate( stateVectorMeanAnlSfcPresMpiGlb, tim_nstepobsinc, hco_ens, vco_ens,   &
                       dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.false., &
                       dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0'/) )
    call gsv_zero(stateVectorMeanAnlSfcPresMpiGlb)
  end if
  call gsv_copy(stateVectorMeanAnl, stateVectorMeanAnlSfcPres, allowMismatch_opt=.true.)
  call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorMeanAnlSfcPres)
  call gsv_transposeTilesToMpiGlobal(stateVectorMeanAnlSfcPresMpiGlb, stateVectorMeanAnlSfcPres)

  !- 7.2 Output ens stddev and mean in trialrms, analrms and analpertrms files

  ! determine middle timestep for output of these files
  middleStepIndex = (tim_nstepobsinc + 1) / 2

  ! output trialmean, trialrms
  call getRmsEtiket(etiketMean, etiketStd, 'F')
  call fln_ensFileName( outFileName, '.', shouldExist_opt=.false. )
  outFileName = trim(outFileName) // '_trialmean'
  call gsv_writeToFile(stateVectorMeanTrl, outFileName, trim(etiketMean),  &
                       typvar_opt='P', writeHeightSfc_opt=.false., numBits_opt=16,  &
                       stepIndex_opt=middleStepIndex, containsFullField_opt=.true.)
  call fln_ensFileName( outFileName, '.', shouldExist_opt=.false. )
  outFileName = trim(outFileName) // '_trialrms'
  call gsv_writeToFile(stateVectorStdDevTrl, outFileName, trim(etiketStd),  &
                       typvar_opt='P', writeHeightSfc_opt=.false., numBits_opt=16, &
                       stepIndex_opt=middleStepIndex, containsFullField_opt=.false.)

  ! output analmean, analrms
  call getRmsEtiket(etiketMean, etiketStd, 'A')
  call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp() )
  outFileName = trim(outFileName) // '_analmean'
  call gsv_writeToFile(stateVectorMeanAnl, outFileName, trim(etiketMean),  &
                       typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=16, &
                       stepIndex_opt=middleStepIndex, containsFullField_opt=.true.)
  call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp() )
  outFileName = trim(outFileName) // '_analrms'
  call gsv_writeToFile(stateVectorStdDevAnl, outFileName, trim(etiketStd),  &
                       typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=16, &
                       stepIndex_opt=middleStepIndex, containsFullField_opt=.false.)

  ! output analpertmean, analpertrms
  call getRmsEtiket(etiketMean, etiketStd, 'P')
  call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp() )
  outFileName = trim(outFileName) // '_analpertmean'
  call gsv_writeToFile(stateVectorMeanAnl, outFileName, trim(etiketMean),  &
                       typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=16, &
                       stepIndex_opt=middleStepIndex, containsFullField_opt=.true.)
  call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp() )
  outFileName = trim(outFileName) // '_analpertrms'
  call gsv_writeToFile(stateVectorStdDevAnlPert, outFileName, trim(etiketStd),  &
                       typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=16, &
                       stepIndex_opt=middleStepIndex, containsFullField_opt=.false.)

  !- 7.3 Output the ensemble mean increment (include MeanAnl Psfc) and analysis

  ! convert transformed to model variables for ensemble mean of analysis and trial
  call gvt_transform(stateVectorMeanAnl,'AllTransformedToModel',allowOverWrite_opt=.true.)
  call gvt_transform(stateVectorMeanTrl,'AllTransformedToModel',allowOverWrite_opt=.true.)
  ! and recompute mean increment for converted model variables (e.g. VIS and PR)
  nullify(varNames)
  call gsv_varNamesList(varNames, stateVectorMeanAnl)
  call gsv_deallocate( stateVectorMeanInc )
  call gsv_allocate( stateVectorMeanInc, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles',  &
                     dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=varNames )
  call gsv_copy(stateVectorMeanAnl, stateVectorMeanInc)
  call gsv_add(stateVectorMeanTrl, stateVectorMeanInc, scaleFactor_opt=-1.0D0)
  deallocate(varNames)

  ! output ensemble mean increment
  call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp(), 0, ensFileNameSuffix_opt='inc' )
  do stepIndex = 1, tim_nstepobsinc
    call gsv_writeToFile(stateVectorMeanInc, outFileName, 'ENSMEAN_INC',  &
                         typvar_opt='R', writeHeightSfc_opt=.false., numBits_opt=16, &
                         stepIndex_opt=stepIndex, containsFullField_opt=.false.)
    call gsv_writeToFile(stateVectorMeanAnlSfcPres, outFileName, 'ENSMEAN_INC',  &
                         typvar_opt='A', writeHeightSfc_opt=.true., &
                         stepIndex_opt=stepIndex, containsFullField_opt=.true.)
  end do

  ! output ensemble mean analysis state
  call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp(), 0 )
  do stepIndex = 1, tim_nstepobsinc
    call gsv_writeToFile(stateVectorMeanAnl, outFileName, 'ENSMEAN_ANL',  &
                         typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=16, &
                         stepIndex_opt=stepIndex, containsFullField_opt=.true.)
  end do

  !- 7.4 Output the deterministic increment (include MeanAnl Psfc) and analysis
  if (deterministicStateExists) then
    ! convert transformed to model variables for deterministic analysis and trial
    call gvt_transform(stateVectorDeterAnl,'AllTransformedToModel',allowOverWrite_opt=.true.)
    call gvt_transform(stateVectorDeterTrl,'AllTransformedToModel',allowOverWrite_opt=.true.)
    ! and recompute mean increment for converted model variables (e.g. VIS and PR)
    nullify(varNames)
    call gsv_varNamesList(varNames, stateVectorDeterAnl)
    call gsv_deallocate( stateVectorDeterInc )
    call gsv_allocate( stateVectorDeterInc, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles',  &
                       dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=varNames )
    call gsv_copy(stateVectorDeterAnl, stateVectorDeterInc)
    call gsv_add(stateVectorDeterTrl, stateVectorDeterInc, scaleFactor_opt=-1.0D0)
    deallocate(varNames)

    ! output the deterministic increment
    call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp(), ensFileNameSuffix_opt='inc' )
    do stepIndex = 1, tim_nstepobsinc
      call gsv_writeToFile(stateVectorDeterInc, outFileName, 'DETER_INC',  &
                           typvar_opt='R', writeHeightSfc_opt=.false., numBits_opt=16, &
                           stepIndex_opt=stepIndex, containsFullField_opt=.false.)
      call gsv_writeToFile(stateVectorMeanAnlSfcPres, outFileName, 'DETER_INC',  &
                           typvar_opt='A', writeHeightSfc_opt=.true., numBits_opt=16, &
                           stepIndex_opt=stepIndex, containsFullField_opt=.true.)
    end do

    ! output the deterministic analysis state
    call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp() )
    do stepIndex = 1, tim_nstepobsinc
      call gsv_writeToFile(stateVectorDeterAnl, outFileName, 'DETER_ANL',  &
                           typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=16, &
                           stepIndex_opt=stepIndex, containsFullField_opt=.true.)
    end do
  end if

  !- 7.5 Output all ensemble member analyses and increments
  ! convert transformed to model variables for analysis and trial ensembles
  call gvt_transform(ensembleAnl,'AllTransformedToModel',allowOverWrite_opt=.true.)
  call gvt_transform(ensembleTrl,'AllTransformedToModel',allowOverWrite_opt=.true.)
  call tmg_start(104,'LETKF-writeEns')
  call ens_writeEnsemble(ensembleAnl, '.', '', ' ', 'ENS_ANL', 'A',  &
                         numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                         containsFullField_opt=.true.)
  call tmg_stop(104)

  ! WARNING: Increment put in ensembleTrl for output
  call ens_add(ensembleAnl, ensembleTrl, scaleFactorInOut_opt=-1.0D0)
  call tmg_start(104,'LETKF-writeEns')
  call ens_writeEnsemble(ensembleTrl, '.', '', ' ', 'ENS_INC', 'R',  &
                         numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                         containsFullField_opt=.false.)
  ! Also write the reference (analysis) surface pressure to increment files
  call writeToAllMembers(stateVectorMeanAnlSfcPresMpiGlb, nEns,  &
                         etiket='ENS_INC', typvar='A', fileNameSuffix='inc',  &
                         ensPath='.')
  call tmg_stop(104)

  !- 7.6 Output the sub-sampled ensemble analyses and increments
  if (writeSubSample) then

    ! Output the ensemble mean increment (include MeanAnl Psfc)
    call fln_ensAnlFileName( outFileName, 'subspace', tim_getDateStamp(), 0, ensFileNameSuffix_opt='inc' )
    do stepIndex = 1, tim_nstepobsinc
      call gsv_writeToFile(stateVectorMeanIncSubSample, outFileName, 'ENSMEAN_INC',  &
                           typvar_opt='R', writeHeightSfc_opt=.false., numBits_opt=16, &
                           stepIndex_opt=stepIndex, containsFullField_opt=.false.)
      call gsv_writeToFile(stateVectorMeanAnlSfcPres, outFileName, 'ENSMEAN_INC',  &
                           typvar_opt='A', writeHeightSfc_opt=.true., &
                           stepIndex_opt=stepIndex, containsFullField_opt=.true.)
    end do

    ! Output the ensemble mean analysis state
    call fln_ensAnlFileName( outFileName, 'subspace', tim_getDateStamp(), 0 )
    do stepIndex = 1, tim_nstepobsinc
      call gsv_writeToFile(stateVectorMeanAnlSubSample, outFileName, 'ENSMEAN_ANL',  &
                           typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=16, &
                           stepIndex_opt=stepIndex, containsFullField_opt=.true.)
    end do

    ! Output the sub-sampled analysis ensemble members
    call tmg_start(104,'LETKF-writeEns')
    call ens_writeEnsemble(ensembleAnlSubSample, 'subspace', '', ' ', 'ENS_ANL', 'A',  &
                           numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                           containsFullField_opt=.true.)
    call tmg_stop(104)

    ! Output the sub-sampled ensemble increments (include MeanAnl Psfc)
    ! WARNING: Increment put in ensembleTrlSubSample for output
    call ens_add(ensembleAnlSubSample, ensembleTrlSubSample, scaleFactorInOut_opt=-1.0D0)
    call tmg_start(104,'LETKF-writeEns')
    call ens_writeEnsemble(ensembleTrlSubSample, 'subspace', '', ' ', 'ENS_INC', 'R',  &
                           numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                           containsFullField_opt=.false.)
    ! Also write the reference (analysis) surface pressure to increment files
    call writeToAllMembers(stateVectorMeanAnlSfcPresMpiGlb,  &
                           ens_getNumMembers(ensembleAnlSubSample),  &
                           etiket='ENS_INC', typvar='A', fileNameSuffix='inc',  &
                           ensPath='subspace')
    call tmg_stop(104)
  end if

  call tmg_stop(4)

  !
  !- 8. MPI, tmg finalize
  !  
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  call tmg_stop(1)

  call tmg_terminate(mpi_myid, 'TMG_LETKF' )
  call rpn_comm_finalize(ierr) 

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  if ( mpi_myid == 0 ) write(*,*) ' --------------------------------'
  if ( mpi_myid == 0 ) write(*,*) ' MIDAS-LETKF ENDS'
  if ( mpi_myid == 0 ) write(*,*) ' --------------------------------'


contains

  subroutine getRmsEtiket(etiketMean, etiketStd, etiketType)
    implicit none

    ! arguments:
    character(len=*) :: etiketMean
    character(len=*) :: etiketStd
    character(len=*) :: etiketType

    if (trim(etiketType) == 'F') then

      ! create trialrms etiket, e.g. E2AVGTRPALL E24_3GMP0256
      etiketStd(1:5) = etiket0(1:5)
      etiketStd(6:7) = 'GM'
      etiketStd(8:8) = etiket0(8:8)
      write(etiketStd(9:12),'(I4.4)') nEns
      etiketMean(1:2) = etiket0(1:2)
      etiketMean(3:7) = 'AVGTR'
      etiketMean(8:8) = etiket0(8:8)
      etiketMean(9:11) = 'ALL'

    else if (trim(etiketType) == 'A') then

      ! create analrms etiket, e.g. E2AVGANPALL, E24_3_0P0256
      etiketStd(1:8) = etiket0(1:8)
      write(etiketStd(9:12),'(I4.4)') nEns
      etiketMean(1:2) = etiket0(1:2)
      etiketMean(3:7)='AVGAN'
      etiketMean(8:8) = etiket0(8:8)
      etiketMean(9:11) = 'ALL'

    else if (trim(etiketType) == 'P') then

      ! create analpertrms etiket, e.g. E2AVGPTPALL, E24_3PTP0256
      etiketStd(1:5) = etiket0(1:5)
      etiketStd(6:7) = 'PT'
      etiketStd(8:8) = etiket0(8:8)
      write(etiketStd(9:12),'(I4.4)') nEns
      etiketMean(1:2) = etiket0(1:2)
      etiketMean(3:7) = 'AVGPT'
      etiketMean(8:8) = etiket0(8:8)
      etiketMean(9:11) = 'ALL'

    else
      call utl_abort('midas-letkf: unknown value of etiketType')
    end if

  end subroutine getRmsEtiket

  subroutine writeToAllMembers(stateVector, nEns, etiket, typvar,  &
                               fileNameSuffix, ensPath)
    implicit none

    ! arguments:
    type(struct_gsv) :: stateVector
    integer          :: nEns
    character(len=*) :: etiket
    character(len=*) :: typvar
    character(len=*) :: fileNameSuffix
    character(len=*) :: ensPath

    ! locals:
    integer :: memberIndex, stepIndex, writeFilePE(nEns)
    character(len=4) :: memberIndexStr

    do memberIndex = 1, nEns
      writeFilePE(memberIndex) = mod(memberIndex-1, mpi_nprocs)
    end do

    do memberIndex = 1, nEns

      if (mpi_myid == writeFilePE(memberIndex)) then

        call fln_ensAnlFileName( outFileName, ensPath, tim_getDateStamp(),  &
                                 memberIndex, ensFileNameSuffix_opt=fileNameSuffix )
        write(memberIndexStr,'(I4.4)') memberIndex

        do stepIndex = 1, tim_nstepobsinc
          call gsv_writeToFile(stateVector, outFileName,  &
                               trim(etiket) // memberIndexStr,  &
                               typvar_opt=trim(typvar), writeHeightSfc_opt=.true., &
                               stepIndex_opt=stepIndex, containsFullField_opt=.true., &
                               numBits_opt=16)
        end do

      end if

    end do

  end subroutine writeToAllMembers

end program midas_letkf

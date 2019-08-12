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
  use utilities_mod
  use ramDisk_mod
  use statetocolumn_mod
  use obsFiles_mod
  use obsSpaceData_mod
  use obsErrors_mod
  use obsFilter_mod
  use obsOperators_mod
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
  integer :: myLonBeg, myLonEnd, myLatBeg, myLatEnd, randomSeed2
  integer :: nulnam, nulFile, dateStamp, datePrint, timePrint, imode, ierr, nsize
  integer :: get_max_rss, fclos, fnom, fstopc, newdate
  integer, allocatable :: dateStampList(:), dateStampListInc(:)

  character(len=256) :: ensFileName, outFileName, deterFileName
  character(len=9)   :: obsColumnMode
  character(len=48)  :: obsMpiStrategy
  character(len=48)  :: midasMode
  character(len=4)   :: memberIndexStr

  logical :: deterExists

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
  logical  :: writeSubSample   ! write sub-sample members for initializing medium-range fcsts
  real(8)  :: hLocalize        ! horizontal localization radius (in km)
  real(8)  :: vLocalize        ! vertical localization radius (in units of ln(Pressure in Pa))
  real(8)  :: alphaRTPP        ! RTPP coefficient (between 0 and 1; 0 means no relaxation)
  real(8)  :: alphaRTPS        ! RTPS coefficient (between 0 and 1; 0 means no relaxation)
  real(8)  :: alphaRandomPert  ! Random perturbation additive inflation coeff (0->1)
  real(8)  :: alphaRandomPertSubSample ! Random perturbation additive inflation coeff for medium-range fcsts
  logical  :: imposeSaturationLimit ! switch for choosing to impose saturation limit of humidity
  logical  :: imposeRttovHuLimits   ! switch for choosing to impose the RTTOV limits on humidity
  character(len=20) :: obsTimeInterpType ! type of time interpolation to obs time
  character(len=20) :: mpiDistribution   ! type of mpiDistribution for weight calculation ('ROUNDROBIN' or 'TILES')
  NAMELIST /NAMLETKF/nEns, ensPathName, hLocalize, vLocalize,  &
                     maxNumLocalObs, weightLatLonStep,  &
                     updateMembers, writeIncrements, writeSubSample,  &
                     alphaRTPP, alphaRTPS, randomSeed, alphaRandomPert, alphaRandomPertSubSample,  &
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

  ! Initialize obs error covariances and set flag using 'util' column of stats_tovs
  call oer_setObsErrors(obsSpaceData, midasMode, useTovsUtil_opt=.true.) ! IN

  ! Call suprep again to filter out channels according to 'util' column of stats_tovs
  call filt_suprep(obsSpaceData)

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
  call enkf_LETKFanalyses(ensembleAnl, ensembleTrl, ensObs_mpiglobal,  &
                          stateVectorMeanInc, stateVectorMeanTrl, stateVectorMeanAnl, &
                          stateVectorDeterInc, stateVectorDeterTrl, stateVectorDeterAnl, &
                          wInterpInfo, maxNumLocalObs,  &
                          hLocalize, vLocalize, updateMembers, alphaRTPP, mpiDistribution)
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
      randomSeed2 = 9 + datePrint - 1000000*(datePrint/1000000)
      write(*,*) 'midas-letkf: randomSeed for additive inflation set to ', randomSeed2
    else
      randomSeed2 = randomSeed
    end if
    call tmg_start(101,'LETKF-randomPert')
    call enkf_addRandomPert(ensembleAnl, stateVectorMeanTrl, alphaRandomPert, randomSeed2)
    call tmg_stop(101)
  end if

  !- 6.4 Impose limits on humidity, if requested
  if (imposeSaturationLimit .or. imposeRttovHuLimits) then 
    call tmg_start(102,'LETKF-imposeHulimits')
    ! Impose limits on analysis ensemble
    call ens_clipHumidity(ensembleAnl, imposeSaturationLimit, imposeRttovHuLimits)
    ! And recompute analysis mean
    call ens_computeMean(ensembleAnl)
    call ens_copyEnsMean(ensembleAnl, stateVectorMeanAnl)
    ! And recompute mean increment
    call gsv_copy(stateVectorMeanAnl, stateVectorMeanInc)
    call gsv_add(stateVectorMeanTrl, stateVectorMeanInc, scaleFactor_opt=-1.0D0)

    ! Impose limits on deterministic analysis
    if (deterExists) then
      if ( imposeSaturationLimit ) call qlim_gsvSaturationLimit(stateVectorDeterAnl)
      if ( imposeRttovHuLimits   ) call qlim_gsvRttovLimit(stateVectorDeterAnl)
      ! And recompute deterministic increment
      call gsv_copy(stateVectorDeterAnl, stateVectorDeterInc)
      call gsv_add(stateVectorDeterTrl, stateVectorDeterInc, scaleFactor_opt=-1.0D0)
    end if
    call tmg_stop(102)
  end if

  !- 6.5 Recompute the analysis spread stddev after inflation and humidity limits
  if (alphaRTPS > 0.0D0 .or. alphaRandomPert > 0.0D0) then
    call ens_computeStdDev(ensembleAnl)
    call ens_copyEnsStdDev(ensembleAnl, stateVectorStdDevAnl)
  end if

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
        randomSeed2 = 9 + datePrint - 1000000*(datePrint/1000000)
        write(*,*) 'midas-letkf: randomSeed for additive inflation set to ', randomSeed2
      else
        randomSeed2 = randomSeed
      end if
      call tmg_start(101,'LETKF-randomPert')
      call enkf_addRandomPert(ensembleAnlSubSample, stateVectorMeanTrl,  &
                              alphaRandomPertSubSample, randomSeed2)
      call tmg_stop(101)
    end if

    ! Compute analysis mean of sub-sampled ensemble
    call ens_computeMean(ensembleAnlSubSample)

    ! Shift members to have same mean as full ensemble and impose humidity limits, if requested
    call ens_recenter(ensembleAnlSubSample, stateVectorMeanAnl, recenteringCoeff=1.0D0, &
                      imposeRttovHuLimits_opt=.true., &
                      imposeSaturationLimit_opt=.true.)

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
    do stepIndex = 1, tim_nstepobsinc
      call gsv_writeToFile(stateVectorDeterAnl, outFileName, 'DETER_ANL',  &
                           typvar_opt='A', writeHeightSfc_opt=.false., &
                           stepIndex_opt=stepIndex, containsFullField_opt=.true.)
    end do
  end if

  !- 7.6 Output all ensemble member analyses
  if (updateMembers) then
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

  !- 7.8 Output the sub-sampled ensemble analyses and increments
  if (writeSubSample) then

    ! Output the ensemble mean increment (include MeanAnl Psfc)
    call fln_ensAnlFileName( outFileName, 'subspace', tim_getDateStamp(), 0, ensFileNameSuffix_opt='inc' )
    do stepIndex = 1, tim_nstepobsinc
      call gsv_writeToFile(stateVectorMeanIncSubSample, outFileName, 'ENSMEAN_INC',  &
                           typvar_opt='R', writeHeightSfc_opt=.false., &
                           stepIndex_opt=stepIndex, containsFullField_opt=.false.)
      call gsv_writeToFile(stateVectorMeanAnlSfcPressure, outFileName, 'ENSMEAN_INC',  &
                           typvar_opt='A', writeHeightSfc_opt=.true., &
                           stepIndex_opt=stepIndex, containsFullField_opt=.true.)
    end do

    ! Output the ensemble mean analysis state
    call fln_ensAnlFileName( outFileName, 'subspace', tim_getDateStamp(), 0 )
    do stepIndex = 1, tim_nstepobsinc
      call gsv_writeToFile(stateVectorMeanAnlSubSample, outFileName, 'ENSMEAN_ANL',  &
                           typvar_opt='A', writeHeightSfc_opt=.false., &
                           stepIndex_opt=stepIndex, containsFullField_opt=.true.)
    end do

    ! Output the sub-sampled analysis ensemble members
    call ens_writeEnsemble(ensembleAnlSubSample, 'subspace', '', ' ', 'ENS_ANL', 'A',  &
                           numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                           containsFullField_opt=.true.)

    ! Output the sub-sampled ensemble increments (include MeanAnl Psfc)
    ! WARNING: Increment put in ensembleTrlSubSample for output
    call ens_add(ensembleAnlSubSample, ensembleTrlSubSample, scaleFactorInOut_opt=-1.0D0)
    call ens_writeEnsemble(ensembleTrlSubSample, 'subspace', '', ' ', 'ENS_INC', 'R',  &
                           numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                           containsFullField_opt=.false.)
    ! Also write the reference (analysis) surface pressure to increment files
    do memberIndex = 1, ens_getNumMembers(ensembleAnlSubSample)
      call fln_ensAnlFileName( outFileName, 'subspace', tim_getDateStamp(),  &
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

  !- 7.9 Output obs files with mean OMP and OMA

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


end program midas_letkf

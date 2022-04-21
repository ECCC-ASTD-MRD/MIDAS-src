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
  !           Note that the actual calculation of the analyses is in the
  !           subroutine enkf_LETKFanalyses.
  !           Many aspects of this program are controlled throught the namelist
  !           block NAMLETKF.
  use version_mod
  use mpi_mod
  use mathPhysConstants_mod
  use fileNames_mod
  use ensembleObservations_mod
  use ensembleStateVector_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use columnData_mod
  use tovs_nl_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use oceanMask_mod
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
  use enkf_mod
  use ensPostProcess_mod
  implicit none

  type(struct_obs), target  :: obsSpaceData
  type(struct_ens), target  :: ensembleTrl4D
  type(struct_ens), pointer :: ensembleTrl
  type(struct_ens)          :: ensembleAnl
  type(struct_gsv)          :: stateVectorMeanTrl4D
  type(struct_gsv)          :: stateVectorMeanAnl
  type(struct_gsv)          :: stateVectorWithZandP4D
  type(struct_gsv)          :: stateVectorHeightSfc
  type(struct_gsv)          :: stateVectorCtrlTrl
  type(struct_gsv)          :: stateVectorRecenter
  type(struct_columnData)   :: column

  type(struct_eob) :: ensObs, ensObs_mpiglobal

  type(struct_vco), pointer :: vco_ens => null()
  type(struct_hco), pointer :: hco_ens => null()
  type(struct_ocm)          :: oceanMask

  integer :: memberIndex, middleStepIndex, stepIndex, randomSeedObs
  integer :: nulnam, dateStamp, ierr
  integer :: get_max_rss, fclos, fnom, fstopc
  integer, allocatable :: dateStampList(:), dateStampListInc(:)

  character(len=256) :: ensFileName, ctrlFileName, recenterFileName
  character(len=9)   :: obsColumnMode
  character(len=48)  :: obsMpiStrategy
  character(len=48)  :: midasMode

  logical :: nwpFields   ! indicates if fields are on momentum and thermo levels
  logical :: oceanFields ! indicates if fields are on depth levels

  ! interpolation information for weights (in enkf_mod)
  type(struct_enkfInterpInfo) :: wInterpInfo

  ! namelist variables
  character(len=20)  :: algorithm  ! name of the chosen LETKF algorithm: 'LETKF', 'CVLETKF'
  logical            :: ensPostProcessing ! do all post-processing of analysis ensemble
  logical            :: recenterInputEns  ! read a deterministic state to recenter ensemble
  integer            :: numSubEns  ! number of sub-ensembles to split the full ensemble
  character(len=256) :: ensPathName ! absolute or relative path to ensemble directory
  integer  :: nEns                 ! ensemble size
  logical  :: randomShuffleSubEns  ! choose to randomly shuffle members into subensembles 
  integer  :: maxNumLocalObs       ! maximum number of obs in each local volume to assimilate
  integer  :: weightLatLonStep     ! separation of lat-lon grid points for weight calculation
  logical  :: modifyAmsubObsError  ! reduce AMSU-B obs error stddev in tropics
  logical  :: backgroundCheck      ! apply additional background check using ensemble spread
  logical  :: huberize             ! apply huber norm quality control procedure
  logical  :: rejectHighLatIR      ! reject all IR observations at high latitudes
  logical  :: rejectRadNearSfc     ! reject radiance observations near the surface
  logical  :: ignoreEnsDate        ! when reading ensemble, ignore the date
  logical  :: outputOnlyEnsMean    ! when writing ensemble, can choose to only write member zero
  real(8)  :: hLocalize(4)         ! horizontal localization radius (in km)
  real(8)  :: hLocalizePressure(3) ! pressures where horizontal localization changes (in hPa)
  real(8)  :: vLocalize            ! vertical localization radius (units: ln(Pressure in Pa) or meters)
  real(8)  :: minDistanceToLand    ! for ice/ocean DA: minimum distance to land for assimilating obs
  character(len=20) :: obsTimeInterpType ! type of time interpolation to obs time
  character(len=20) :: mpiDistribution   ! type of mpiDistribution for weight calculation ('ROUNDROBIN' or 'TILES')
  character(len=12) :: etiket_anl        ! etiket for output files
  NAMELIST /NAMLETKF/algorithm, ensPostProcessing, recenterInputEns, nEns, numSubEns, &
                     ensPathName, randomShuffleSubEns,  &
                     hLocalize, hLocalizePressure, vLocalize, minDistanceToLand,  &
                     maxNumLocalObs, weightLatLonStep,  &
                     modifyAmsubObsError, backgroundCheck, huberize, rejectHighLatIR, rejectRadNearSfc,  &
                     ignoreEnsDate, outputOnlyEnsMean, obsTimeInterpType, mpiDistribution, etiket_anl

  ! Some high-level configuration settings
  midasMode = 'analysis'
  obsColumnMode = 'ENKFMIDAS'
  obsMpiStrategy = 'LIKESPLITFILES'

  call ver_printNameAndVersion('letkf','Program for Local Ensemble Transform Kalman Filter')

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
  ensPostProcessing     = .false.
  recenterInputEns      = .false.
  ensPathName           = 'ensemble'
  nEns                  = 10
  numSubEns             = 2
  randomShuffleSubEns   = .false.
  maxNumLocalObs        = 1000
  weightLatLonStep      = 1
  modifyAmsubObsError   = .false.
  backgroundCheck       = .false.
  huberize              = .false.
  rejectHighLatIR       = .false.
  rejectRadNearSfc      = .false.
  ignoreEnsDate         = .false.
  outputOnlyEnsMean     = .false.
  hLocalize(:)          = -1.0D0
  hLocalizePressure     = (/14.0D0, 140.0D0, 400.0D0/)
  vLocalize             = -1.0D0
  minDistanceToLand     = -1.0D0
  obsTimeInterpType     = 'LINEAR'
  mpiDistribution       = 'ROUNDROBIN'
  etiket_anl            = 'ENS_ANL'

  !- 1.2 Read the namelist
  nulnam = 0
  ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namletkf, iostat=ierr)
  if ( ierr /= 0) call utl_abort('midas-letkf: Error reading namelist')
  if ( mpi_myid == 0 ) write(*,nml=namletkf)
  ierr = fclos(nulnam)

  !- 1.3 Some minor modifications of namelist values
  if (hLocalize(1) > 0.0D0 .and. hLocalize(2) < 0.0D0) then
    ! if only 1 value given for hLocalize, use it for entire column
    hLocalize(2:4) = hLocalize(1)
  end if
  hLocalize(:) = hLocalize(:) * 1000.0D0 ! convert from km to m
  hLocalizePressure(:) = log(hLocalizePressure(:) * MPC_PA_PER_MBAR_R8)

  if (minDistanceToLand > 0.0D0) then
    minDistanceToLand = minDistanceToLand * 1000.0D0 ! convert from km to m
  end if

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

  write(*,*) 'midas-letkf: analysis dateStamp = ',dateStamp

  !- 2.3 Initialize variables of the model states
  call gsv_setup

  !- 2.4 Initialize the Ensemble grid
  if (mpi_myid == 0) write(*,*) ''
  if (mpi_myid == 0) write(*,*) 'midas-letkf: Set hco and vco parameters for ensemble grid'
  call fln_ensFileName( ensFileName, ensPathName, memberIndex_opt=1, &
                        copyToRamDisk_opt=.false. )
  call hco_SetupFromFile( hco_ens, ensFileName, ' ', 'ENSFILEGRID')
  call vco_setupFromFile( vco_ens, ensFileName )
  if (vco_getNumLev(vco_ens, 'MM') /= vco_getNumLev(vco_ens, 'TH')) then
    call utl_abort('midas-letkf: nLev_M /= nLev_T - currently not supported')
  end if
  nwpFields   = (vco_getNumLev(vco_ens,'TH') > 0 .or. vco_getNumLev(vco_ens,'MM') > 0)
  oceanFields = (vco_getNumLev(vco_ens,'DP') > 0)
  if (.not.nwpFields .and. .not.oceanFields) then
    call utl_abort('midas-letkf: vertical coordinate does not contain nwp nor ocean fields')
  end if

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !- 2.5 Read in the observations and other obs-related set up

  ! Read the observations
  call inn_setupObs( obsSpaceData, hco_ens, obsColumnMode, obsMpiStrategy, midasMode,  &
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

  !- 2.7 Read the sfc height from ensemble member 1 - only if we are doing NWP
  if ( nwpFields ) then
    call gsv_allocate( stateVectorHeightSfc, 1, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0','TT'/) )
    call gio_readFromFile( stateVectorHeightSfc, ensFileName, ' ', ' ',  &
                           containsFullField_opt=.true., readHeightSfc_opt=.true. )
  end if

  !- 2.8 Allocate various statevectors related to ensemble mean
  call gsv_allocate( stateVectorMeanTrl4D, tim_nstepobs, hco_ens, vco_ens, &
                     dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeightSfc_opt=.true., &
                     allocHeight_opt=.false., allocPressure_opt=.false. )
  call gsv_zero(stateVectorMeanTrl4D)
  call gsv_allocate( stateVectorMeanAnl, tim_nstepobsinc, hco_ens, vco_ens, &
                     dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeightSfc_opt=.true., &
                     allocHeight_opt=.false., allocPressure_opt=.false. )
  call gsv_zero(stateVectorMeanAnl)

  !- 2.9 Allocate statevector for storing state with heights and pressures allocated (for s2c_nl)
  call gsv_allocate( stateVectorWithZandP4D, tim_nstepobs, hco_ens, vco_ens, &
                     dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeightSfc_opt=.true. )
  call gsv_zero(stateVectorWithZandP4D)

  !- 2.10 Allocate ensembles, read the Trl ensemble
  call ens_allocate(ensembleTrl4D, nEns, tim_nstepobs, hco_ens, vco_ens, dateStampList)
  call ens_readEnsemble(ensembleTrl4D, ensPathName, biPeriodic=.false., &
                        ignoreDate_opt=ignoreEnsDate)

  !- 2.11 If desired, read a deterministic state for recentering the ensemble
  if (recenterInputEns) then
    call gsv_allocate( stateVectorRecenter, tim_nstepobs, hco_ens, vco_ens, &
                       dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.false., &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_zero(stateVectorRecenter)
    call fln_ensTrlFileName( recenterFileName, './', tim_getDateStamp() )
    do stepIndex = 1, tim_nstepobs
      call gio_readFromFile( stateVectorRecenter, recenterFileName, ' ', ' ',  &
                             stepIndex_opt=stepIndex, containsFullField_opt=.true., &
                             readHeightSfc_opt=.false. )
    end do
    call ens_recenter( ensembleTrl4D, stateVectorRecenter, recenteringCoeff_opt=1.0d0 )
    call gsv_deallocate( stateVectorRecenter )
  end if
  
  !- 2.12 Compute ensemble mean and copy to meanTrl and meanAnl stateVectors
  call ens_computeMean(ensembleTrl4D)
  call ens_copyEnsMean(ensembleTrl4D, stateVectorMeanTrl4D)
  if (tim_nstepobsinc < tim_nstepobs) then
    call gsv_copy4Dto3D(stateVectorMeanTrl4D, stateVectorMeanAnl)
  else
    call gsv_copy(stateVectorMeanTrl4D, stateVectorMeanAnl)
  end if
  
  !
  !- 3. Compute HX values with results in ensObs
  !

  !- 3.1 Loop over all members and compute HX for each
  do memberIndex = 1, nEns

    write(*,*) ''
    write(*,*) 'midas-letkf: apply nonlinear H to ensemble member ', memberIndex
    write(*,*) ''

    ! copy 1 member to a stateVector
    call ens_copyMember(ensembleTrl4D, stateVectorWithZandP4D, memberIndex)

    ! copy the surface height field
    if (nwpFields) then
      call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorWithZandP4D)
    end if

    call s2c_nl( stateVectorWithZandP4D, obsSpaceData, column, hco_ens, &
                 timeInterpType=obsTimeInterpType, dealloc_opt=.false. )

    ! Compute Y-H(X) in OBS_OMP
    call tmg_start(6,'LETKF-obsOperators')
    call inn_computeInnovation(column, obsSpaceData, beSilent_opt=.true.)
    call tmg_stop(6)

    ! Copy to ensObs: Y-HX for this member
    call eob_setYb(ensObs, memberIndex)

  end do

  !- 3.2 Set some additional information in ensObs and additional quality 
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
  call gsv_copy(stateVectorMeanTrl4D, stateVectorWithZandP4D, allowVarMismatch_opt=.true.)
  if (nwpFields) then
    call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorWithZandP4D)
  end if
  call s2c_nl( stateVectorWithZandP4D, obsSpaceData, column, hco_ens, &
               timeInterpType=obsTimeInterpType, dealloc_opt=.false. )
  call tvs_allocTransmission(col_getNumLev(column,'TH')) ! this will cause radiative transmission profiles to be stored for use in eob_setVertLocation
  call tmg_start(6,'LETKF-obsOperators')
  call inn_computeInnovation(column, obsSpaceData, beSilent_opt=.false.)
  call tmg_stop(6)

  ! Put y-mean(H(X)) in OBS_OMP for writing to obs files (overwrites y-H(mean(X)))
  call eob_setMeanOMP(ensObs)

  ! Set vertical location for all obs for vertical localization (based on ensemble mean pressure and height)
  if (vLocalize > 0.0d0) then
    if (nwpFields) then
      call eob_setTypeVertCoord(ensObs,'logPressure')
    else if (oceanFields) then
      call eob_setTypeVertCoord(ensObs,'depth')
    end if
    call eob_setVertLocation(ensObs, column)
  end if

  ! Modify the obs error stddev for AMSUB in the tropics
  if (modifyAmsubObsError) call enkf_modifyAmsubObsError(obsSpaceData)

  ! Apply a background check (reject limit is set in the routine)
  if (backgroundCheck) call eob_backgroundCheck(ensObs)

  ! For ice/ocean DA: remove obs that are too close to land
  if (minDistanceToLand > 0.0D0) then
    call ens_copyMask(ensembleTrl4D,oceanMask)
    call eob_removeObsNearLand(ensObs, oceanMask, minDistanceToLand)
  end if

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
  call enkf_setupInterpInfo(wInterpInfo, stateVectorMeanAnl%hco, weightLatLonStep,  &
                            stateVectorMeanAnl%myLonBeg, stateVectorMeanAnl%myLonEnd,  &
                            stateVectorMeanAnl%myLatBeg, stateVectorMeanAnl%myLatEnd)
  call tmg_stop(92)

  !
  !- 5. Main calculation of ensemble analyses
  !
  call tmg_start(3,'LETKF-doAnalysis')
  call enkf_LETKFanalyses(algorithm, numSubEns, randomShuffleSubEns,  &
                          ensembleAnl, ensembleTrl, ensObs_mpiglobal,  &
                          stateVectorMeanAnl, &
                          wInterpInfo, maxNumLocalObs,  &
                          hLocalize, hLocalizePressure, vLocalize, mpiDistribution)
  call tmg_stop(3)

  !- 6. Output obs files with mean OMP and (unrecentered) OMA

  ! Compute Y-H(Xa_mean) in OBS_OMA
  write(*,*) ''
  write(*,*) 'midas-letkf: apply nonlinear H to ensemble mean analysis'
  write(*,*) ''
  if (tim_nstepobsinc < tim_nstepobs) then
    ! meanAnl is only 3D, so need to make 4D for s2c_nl
    middleStepIndex = (tim_nstepobs + 1) / 2
    call gsv_copy(stateVectorMeanAnl, stateVectorWithZandP4D, allowVarMismatch_opt=.true., stepIndexOut_opt=middleStepIndex)
    call gsv_3dto4d(stateVectorWithZandP4D)
  else
    call gsv_copy(stateVectorMeanAnl, stateVectorWithZandP4D, allowVarMismatch_opt=.true.)
  end if
  if (nwpFields) then
    call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorWithZandP4D)
  end if
  call s2c_nl( stateVectorWithZandP4D, obsSpaceData, column, hco_ens, &
               timeInterpType=obsTimeInterpType )
  call tmg_start(6,'LETKF-obsOperators')
  call inn_computeInnovation(column, obsSpaceData, destObsColumn_opt=OBS_OMA, beSilent_opt=.false.)
  call tmg_stop(6)

  ! Write (update) observation files
  call obsf_writeFiles( obsSpaceData )

  !- 7. Post processing of the analysis results (if desired) and write everything to files
  if (ensPostProcessing) then
    !- Allocate and read the Trl control member (used to compute control member increment for IAU)
    call gsv_allocate( stateVectorCtrlTrl, tim_nstepobsinc, hco_ens, vco_ens, &
                       dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    if (recenterInputEns) then
      !- Use the deterministic trial, if we are recentering the input ensemble
      ctrlFileName = trim(recenterFileName)
    else
      !- Otherwise, use member 0000
      call fln_ensFileName(ctrlFileName, ensPathName, memberIndex_opt=0, &
                           copyToRamDisk_opt=.false.)
    end if
    do stepIndex = 1, tim_nstepobsinc
      call gio_readFromFile( stateVectorCtrlTrl, ctrlFileName, ' ', ' ',  &
                             stepIndex_opt=stepIndex, containsFullField_opt=.true., &
                             readHeightSfc_opt=.false. )
    end do

    call tmg_start(8,'LETKF-postProcess')
    call epp_postProcess(ensembleTrl, ensembleAnl, stateVectorHeightSfc, stateVectorCtrlTrl, &
                         writeTrlEnsemble=.false., outputOnlyEnsMean_opt=outputOnlyEnsMean)
    call tmg_stop(8)
  else
    !- Just write the raw analysis ensemble to files
    if (mpi_myid == 0) then
      write(*,*) 'midas-letkf: No ensemble post-processing requested, so just write the raw analysis ensemble'
    end if
    call tmg_start(104,'LETKF-writeEns')
    if (.not. outputOnlyEnsMean) then
      call ens_writeEnsemble(ensembleAnl, '.', '', etiket_anl, 'A',  &
                             numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                             containsFullField_opt=.true.)
    end if
    call tmg_stop(104)

  end if

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

end program midas_letkf

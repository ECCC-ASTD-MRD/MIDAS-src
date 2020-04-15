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
  use enkf_mod
  implicit none

  type(struct_obs), target  :: obsSpaceData
  type(struct_ens), target  :: ensembleTrl4D
  type(struct_ens), pointer :: ensembleTrl
  type(struct_ens)          :: ensembleAnl
  type(struct_gsv)          :: stateVectorMeanTrl4D
  type(struct_gsv)          :: stateVectorMeanAnl
  type(struct_gsv)          :: stateVectorMeanInc
  type(struct_gsv)          :: stateVectorWithZandP4D
  type(struct_gsv)          :: stateVectorStdDevTrl
  type(struct_gsv)          :: stateVectorStdDevAnl
  type(struct_gsv)          :: stateVectorStdDevAnlPert
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

  character(len=256) :: ensFileName
  character(len=9)   :: obsColumnMode
  character(len=48)  :: obsMpiStrategy
  character(len=48)  :: midasMode
  character(len=4)   :: memberIndexStr

  ! interpolation information for weights (in enkf_mod)
  type(struct_enkfInterpInfo) :: wInterpInfo

  ! namelist variables
  character(len=20)  :: algorithm  ! name of the chosen LETKF algorithm: 'LETKF', 'CVLETKF'
  logical            :: ensPostProcessing ! do all post-processing of analysis ensemble
  integer            :: numSubEns  ! number of sub-ensembles to split the full ensemble
  character(len=256) :: ensPathName ! absolute or relative path to ensemble directory
  integer  :: nEns                 ! ensemble size
  integer  :: maxNumLocalObs       ! maximum number of obs in each local volume to assimilate
  integer  :: weightLatLonStep     ! separation of lat-lon grid points for weight calculation
  real(8)  :: alphaRTPP            ! RTPP coefficient (between 0 and 1; 0 means no relaxation)
  logical  :: modifyAmsubObsError  ! reduce AMSU-B obs error stddev in tropics
  logical  :: backgroundCheck      ! apply additional background check using ensemble spread
  logical  :: huberize             ! apply huber norm quality control procedure
  logical  :: rejectHighLatIR      ! reject all IR observations at high latitudes
  logical  :: rejectRadNearSfc     ! reject radiance observations near the surface
  real(8)  :: hLocalize(4)         ! horizontal localization radius (in km)
  real(8)  :: hLocalizePressure(3) ! pressures where horizontal localization changes (in hPa)
  real(8)  :: vLocalize            ! vertical localization radius (in units of ln(Pressure in Pa))
  character(len=20) :: obsTimeInterpType ! type of time interpolation to obs time
  character(len=20) :: mpiDistribution   ! type of mpiDistribution for weight calculation ('ROUNDROBIN' or 'TILES')
  NAMELIST /NAMLETKF/algorithm, ensPostProcessing, nEns, numSubEns, ensPathName,  &
                     hLocalize, hLocalizePressure, vLocalize,  &
                     maxNumLocalObs, weightLatLonStep,  &
                     modifyAmsubObsError, backgroundCheck, huberize, rejectHighLatIR, rejectRadNearSfc,  &
                     alphaRTPP, obsTimeInterpType, mpiDistribution


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
  ensPostProcessing     = .false.
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
  hLocalize(:)          = -1.0D0
  hLocalizePressure     = (/14.0D0, 140.0D0, 400.0D0/)
  vLocalize             = -1.0D0
  alphaRTPP             =  0.0D0
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
  if (hLocalize(1) > 0.0D0 .and. hLocalize(2) < 0.0D0) then
    ! if only 1 value given for hLocalize, use it for entire column
    hLocalize(2:4) = hLocalize(1)
  end if
  hLocalize(:) = hLocalize(:) * 1000.0D0 ! convert from km to m
  hLocalizePressure(:) = log(hLocalizePressure(:) * MPC_PA_PER_MBAR_R8)
  if (alphaRTPP < 0.0D0) alphaRTPP = 0.0D0
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

    call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorWithZandP4D)
    call s2c_nl( stateVectorWithZandP4D, obsSpaceData, column, timeInterpType=obsTimeInterpType, dealloc_opt=.false. )

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
  call enkf_setupInterpInfo(wInterpInfo, stateVectorMeanAnl%hco, weightLatLonStep,  &
                            stateVectorMeanAnl%myLonBeg, stateVectorMeanAnl%myLonEnd,  &
                            stateVectorMeanAnl%myLatBeg, stateVectorMeanAnl%myLatEnd)
  call tmg_stop(92)

  !
  !- 5. Main calculation of ensemble analyses
  !
  call tmg_start(3,'LETKF-doAnalysis')
  call enkf_LETKFanalyses(algorithm, numSubEns,  &
                          ensembleAnl, ensembleTrl, ensObs_mpiglobal,  &
                          stateVectorMeanAnl, &
                          wInterpInfo, maxNumLocalObs,  &
                          hLocalize, hLocalizePressure, vLocalize, alphaRTPP, mpiDistribution)
  call tmg_stop(3)

  !- 6. Output obs files with mean OMP and (unrecentered) OMA

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

  !- 7. Post processing of the analysis results (if desired) and write everything to files
  if (ensPostProcessing) then
    call tmg_start(8,'LETKF-postProcess')
    call enkf_postProcess(ensembleTrl, ensembleAnl, stateVectorHeightSfc)
    call tmg_stop(8)
  else
    ! just write the raw analysis ensemble to files
    if (mpi_myid == 0) then
      write(*,*) 'midas-letkf: No ensemble post-processing requested, so just write the raw analysis ensemble'
    end if
    call tmg_start(104,'LETKF-writeEns')
    call ens_writeEnsemble(ensembleAnl, '.', '', ' ', 'ENS_ANL', 'A',  &
                           numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                           containsFullField_opt=.true.)
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

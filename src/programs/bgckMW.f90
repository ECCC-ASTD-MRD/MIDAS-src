!--------------------------------------- LICENCE BEGIN -----------------------------------
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

program midas_bgckMW
  !
  !:Purpose: Main program for background check of passive microwave instruments. 
  !
  use bgckmicrowave_mod
  use MathPhysConstants_mod
  use utilities_mod
  use mpi_mod

  implicit none

  integer                       :: reportIndex                                        ! report index
  integer                       :: lastReportIndex                                    ! report index
  integer                       :: reportNumObs                                       ! number of obs in current report
  integer                       :: reportNumChannel                                   ! "      "   channels "      "
  integer                       :: nobs_tot                                           ! obs. total number
  integer                       :: n_bad_reps                                         ! bad reports number
  integer                       :: satIndexObserrFile                                 ! satellite index in obserror file
  logical                       :: resumeReport                                       ! logical to newInformationFlagify resume report
  logical                       :: ifLastReport                                       ! logical to newInformationFlagify last report
  character(len=9)              :: burpFileSatId                                      ! station id in burp file
  character(len=9), allocatable :: satelliteId(:)                                     ! satellite Id in stats error file
  real, allocatable             :: modelInterpTerrain(:)                              ! topo in standard file interpolated to obs point
  real, allocatable             :: modelInterpSeaIce(:)                               ! Glace de mer " "
  real, allocatable             :: modelInterpGroundIce(:)                            ! Glace de continent " "
  real,    allocatable          :: obsLatitude(:)                                     ! obs. point latitudes
  real,    allocatable          :: obsLongitude(:)                                    ! obs. point longitude
  integer, allocatable          :: satIdentifier(:)                                   ! Satellite identifier
  real,    allocatable          :: satZenithAngle(:)                                  ! sat. satZenithAngle angle
  integer, allocatable          :: landQualifierIndice(:)                             ! land qualifyer
  integer, allocatable          :: terrainTypeIndice(:)                               ! terrain type
  real,    allocatable          :: obsTemperatureBrillance(:)                         ! temperature de brillance
  real,    allocatable          :: ompTemperatureBrillance(:)                         ! o-p temperature de "
  real,    allocatable          :: biasCorr(:)                                        ! bias correction fo obsTemperatureBrillance
  integer, allocatable          :: satelliteScanPosition(:)                           ! scan position
  integer, allocatable          :: obsQcFlag1(:,:)                                    ! Obs Quality flag 1
  integer, allocatable          :: obsQcFlag2(:)                                      ! Obs Quality flag 2 
  integer, allocatable          :: obsChannels(:)                                     ! obsTemperatureBrillance channels
  integer, allocatable          :: ompChannels(:)                                     ! zomp channel
  integer, allocatable          :: obsFlags(:)                                        ! obs. flag
  integer, allocatable          :: satOrbit(:)                                        ! orbit
  integer, allocatable          :: obsGlobalMarker(:)                                 ! global marker
  integer                       :: IUTILST(mwbg_maxNumChan,mwbg_maxNumSat)            ! channel use option for each sat.
  real                          :: TOVERRST(mwbg_maxNumChan,mwbg_maxNumSat)           ! obs. error per channel for each sat
  integer                       :: rejectionCodArray(mwbg_maxNumTest,&
                                                     mwbg_maxNumChan,mwbg_maxNumSat)  ! number of rejection 
  !                                                                                                   per sat. per channl per test
  integer                       :: rejectionCodArray2(mwbg_maxNumTest,&
                                                      mwbg_maxNumChan,mwbg_maxNumSat) ! number of rejection per 
  !                                                                                                   sat. per channl per test
  !                                                                                                   for ATMS 2nd category of tests
  integer, allocatable          :: qcIndicator(:,:)                                   ! indicateur controle de qualite tovs par canal 
  !                                                                                    =0, ok,
  !                                                                                    >0, rejet,
  integer, allocatable          :: globalQcIndicator(:)                                ! indicateur global controle de qualite tovs. Code:
  !                                                                                    =0, ok,
  !                                                                                    >0, rejet d'au moins un canal
  integer, allocatable          :: newInformationFlag(:)                               ! ATMS Information flag (newInformationFlag) values 
  !                                                                                    (new BURP element  025174 in header). FOR AMSUA 
  !                                                                                    just fill with zeros 
  real,    allocatable          :: cloudLiquidWaterPath(:)                             ! cloud liquid water. NB: for AMSUA, 
  !                                                                                    cloudLiquidWaterPath=0.5(model_cloudLiquidWaterPath 
  !                                                                                    + obs_cloudLiquidWaterPath)
  real,    allocatable          :: atmScatteringIndex(:)                               ! scattering index
  integer, external             :: exdb, exfin, fnom, fclos
  integer                       :: ier, istat, nulnam
  ! Temporary arrays
  integer, allocatable          :: reportNumObsTmp(:)                                       ! number of obs in current report
  integer, allocatable          :: reportNumChannelTmp(:)                                   ! "      "   channels "      "
  real,    allocatable          :: obsLatitudeTmp(:,:)                                     ! obs. point latitudes
  real,    allocatable          :: obsLongitudeTmp(:,:)                                    ! obs. point longitude
  integer, allocatable          :: satIdentifierTmp(:,:)                                   ! Satellite identifier
  real,    allocatable          :: satZenithAngleTmp(:,:)                                  ! sat. satZenithAngle angle
  integer, allocatable          :: landQualifierIndiceTmp(:,:)                             ! land qualifyer
  integer, allocatable          :: terrainTypeIndiceTmp(:,:)                               ! terrain type
  real,    allocatable          :: obsTemperatureBrillanceTmp(:,:)                         ! temperature de brillance
  real,    allocatable          :: ompTemperatureBrillanceTmp(:,:)                         ! o-p temperature de "
  real,    allocatable          :: biasCorrTmp(:,:)                                        ! bias correction fo obsTemperatureBrillance
  integer, allocatable          :: satelliteScanPositionTmp(:,:)                           ! scan position
  integer, allocatable          :: obsQcFlag1Tmp(:,:,:)                                    ! Obs Quality flag 1
  integer, allocatable          :: obsQcFlag2Tmp(:,:)                                      ! Obs Quality flag 2 
  integer, allocatable          :: obsChannelsTmp(:,:)                                     ! obsTemperatureBrillance channels
  integer, allocatable          :: ompChannelsTmp(:,:)                                     ! zomp channel
  integer, allocatable          :: obsFlagsTmp(:,:)                                        ! obs. flag
  integer, allocatable          :: satOrbitTmp(:,:)                                        ! orbit
  integer, allocatable          :: obsGlobalMarkerTmp(:,:)                                 ! global marker
  logical, allocatable          :: listeOfGoodReport(:)
  logical, allocatable          :: listeOfResumeReport(:)
  ! namelist variables
  character(len=9)              :: instName                        ! instrument name
  character(len=90)             :: burpFileNameIn                  ! burp input file name
  character(len=90)             :: burpFileNameOut                 ! burp output file name
  character(len=9)              :: ETIKRESU                        ! resume etiket name
  character(len=128)            :: mglg_file                       ! glace de mer file
  character(len=128)            :: statsFile                       ! stats error file
  real                          :: clwQcThreshold                  ! 
  logical                       :: allowStateDepSigmaObs           !
  logical                       :: useUnbiasedObsForClw            !
  logical                       :: RESETQC                         ! reset Qc flags option
  logical                       :: debug                           ! debug mode
  logical                       :: writeModelLsqTT                 ! logical for writing lsq and tt in file
  logical                       :: writeEle25174                   ! logical for writing ele 25174 in file
  logical                       :: writeTbValuesToFile             ! logical for replacing missing tb value
  integer                       :: reportNumMax                    ! Max number of reports in file
  integer                       :: locationNumMax                  ! Max number of obs per report
  integer                       :: channelNumMax                   ! Max number of channel in report

  namelist /nambgck/instName, burpFileNameIn, burpFileNameOut, mglg_file, statsFile, &
                    writeModelLsqTT, writeEle25174, clwQcThreshold, allowStateDepSigmaObs, &
                    useUnbiasedObsForClw, debug, RESETQC, ETIKRESU, writeTbValuesToFile

  namelist/nammwobs/reportNumMax, locationNumMax, channelNumMax

  istat = exdb('midas-bgckMW','DEBUT','NON')

  write(*,'(/,' //                                                &
            '3(" *****************"),/,' //                       &
            '14x,"-- START OF MAIN PROGRAM midas-bgckMW: --",/,' //   &
            '14x,"-- BACKGROUND CHECK FOR MW OBSERVATIONS --",/, ' //&
            '14x,"-- Revision : ",a," --",/,' //       &
            '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  !- 1.0 mpi
  call mpi_initialize

  !- 1.1 timings
  call tmg_init(mpi_myid, 'TMG_BGCKMW' )
  call tmg_start(1,'MAIN')
 
  ! default nambgck namelist values
  instName              = 'AMSUA'
  burpFileNameIn        = './obsto_amsua'
  burpFileNameOut       = './obsto_amsua.out'
  mglg_file             = './fstmglg'
  statsFile             = './stats_amsua_assim'  
  writeModelLsqTT       = .false.
  writeEle25174         = .false.
  clwQcThreshold        = 0.3 
  allowStateDepSigmaObs = .false.
  useUnbiasedObsForClw  = .false.
  debug                 = .false.
  RESETQC               = .false.
  ETIKRESU              = '>>BGCKALT'
  writeTbValuesToFile   = .false.

  ! reading nambgck namelist
  nulnam = 0
  ier = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read(nulnam, nml=nambgck, iostat=ier)
  if ( ier /= 0 ) then 
    call utl_abort('midas_bgckmw: Error reading nambgck namelist')
  end if 
  write(*,nml=nambgck)
  ier = fclos(nulnam)

  mwbg_debug = debug
  mwbg_clwQcThreshold = clwQcThreshold 
  mwbg_allowStateDepSigmaObs = allowStateDepSigmaObs
  mwbg_useUnbiasedObsForClw = useUnbiasedObsForClw

  ! default nammwobs namelist
  reportNumMax = 1000
  locationNumMax = 3000
  channelNumMax = 22

  ! reading nambgck namelist
  nulnam = 0
  ier = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read(nulnam, nml=nammwobs, iostat=ier)
  if ( ier /= 0 ) then 
    call utl_abort('midas_bgckmw: Error reading nammwobs namelist')
  end if 
  write(*,nml=nammwobs)
  ier = fclos(nulnam)
  

  ! Initializations of counters (for total reports/locations in the file).
  n_bad_reps = 0
  nobs_tot = 0 
  
  !#################################################################################
  ! STEP 1 ) Lecture des statistiques d'erreur totale pour les  TOVS 
  !#################################################################################
  call mwbg_readStatTovs(statsFile, instName, satelliteId, IUTILST, TOVERRST)

  ! Allocation of the temporary big arrays to assign in all file data

  call mwbg_allocate1DLogicalArray(listeOfResumeReport,reportNumMax)
  call mwbg_allocate1DLogicalArray(listeOfGoodReport,reportNumMax)
  call mwbg_allocate1DIntegerArray(reportNumObsTmp,reportNumMax)
  call mwbg_allocate1DIntegerArray(reportNumChannelTmp,reportNumMax)
  call mwbg_allocate2DIntegerArray(satIdentifierTmp,reportNumMax,locationNumMax)
  call mwbg_allocate2DIntegerArray(landQualifierIndiceTmp,reportNumMax,locationNumMax)
  call mwbg_allocate2DIntegerArray(terrainTypeIndiceTmp,reportNumMax,locationNumMax)
  call mwbg_allocate2DRealArray(obsLatitudeTmp,reportNumMax,locationNumMax)
  call mwbg_allocate2DRealArray(obsLongitudeTmp,reportNumMax,locationNumMax)
  call mwbg_allocate2DRealArray(satZenithAngleTmp,reportNumMax,locationNumMax)
  
  call mwbg_allocate2DIntegerArray(satelliteScanPositionTmp,reportNumMax,locationNumMax)
  call mwbg_allocate3DIntegerArray(obsQcFlag1Tmp,reportNumMax,locationNumMax,3)
  call mwbg_allocate2DIntegerArray(obsQcFlag2Tmp,reportNumMax,locationNumMax)
  call mwbg_allocate2DIntegerArray(satOrbitTmp,reportNumMax,locationNumMax)
  call mwbg_allocate2DIntegerArray(obsGlobalMarkerTmp,reportNumMax,locationNumMax)

  call mwbg_allocate2DRealArray(obsTemperatureBrillanceTmp,reportNumMax,locationNumMax*channelNumMax)
  call mwbg_allocate2DRealArray(biasCorrTmp,reportNumMax,locationNumMax*channelNumMax)
  call mwbg_allocate2DRealArray(ompTemperatureBrillanceTmp,reportNumMax,locationNumMax*channelNumMax)
  call mwbg_allocate2DIntegerArray(obsChannelsTmp,reportNumMax,locationNumMax*channelNumMax)
  call mwbg_allocate2DIntegerArray(ompChannelsTmp,reportNumMax,locationNumMax*channelNumMax)
  call mwbg_allocate2DIntegerArray(obsFlagsTmp,reportNumMax,locationNumMax*channelNumMax)
  ! MAIN LOOP through all the reports in the BURP file
  listeOfGoodReport = .false.
  listeOfResumeReport = .false.
  reportIndex = 0
  REPORTS: do 

    reportIndex = reportIndex + 1
    !###############################################################################
    ! STEP 0 ) Read TOV Observations in big arrays
    !###############################################################################
    write(*,*) ' ==> mwbg_getData: '
    call mwbg_getData(burpFileNameIn, reportIndex, satIdentifier, satZenithAngle, landQualifierIndice, &
                      terrainTypeIndice, obsLatitude, obsLongitude, obsTemperatureBrillance, biasCorr, ompTemperatureBrillance, &
                      satelliteScanPosition, reportNumChannel, reportNumObs, obsQcFlag1, obsQcFlag2, obsChannels, ompChannels, &
                      obsFlags, satOrbit, obsGlobalMarker, resumeReport, ifLastReport, &
                      instName, burpFileSatId)

    if (resumeReport) listeOfResumeReport(reportIndex) = .true.
    if (ifLastReport) lastReportIndex = reportIndex
    if (.not. resumeReport) then
      if (ALL(ompTemperatureBrillance(:) == MPC_missingValue_R4)) then 
        n_bad_reps = n_bad_reps + 1  
        cycle REPORTS
      end if
      ! Increment total number of obs pts read
      nobs_tot = nobs_tot + reportNumObs
      ! copy obs content in temporary arrays
      listeOfGoodReport(reportIndex) = .true.
      reportNumObsTmp(reportIndex) = reportNumObs
      reportNumChannelTmp(reportIndex) = reportNumChannel
      satIdentifierTmp(reportIndex,1:size(satIdentifier))           = satIdentifier
      satZenithAngleTmp(reportIndex,1:size(satZenithAngle))          = satZenithAngle
      landQualifierIndiceTmp(reportIndex,1:size(landQualifierIndice))     = landQualifierIndice  
      terrainTypeIndiceTmp(reportIndex,1:size(terrainTypeIndice))       = terrainTypeIndice
      obsLatitudeTmp(reportIndex,1:size(obsLatitude))             = obsLatitude
      obsLongitudeTmp(reportIndex,1:size(obsLongitude))            = obsLongitude
      obsTemperatureBrillanceTmp(reportIndex,1:size(obsTemperatureBrillance)) = obsTemperatureBrillance 
      biasCorrTmp(reportIndex,1:size(biasCorr))                = biasCorr
      ompTemperatureBrillanceTmp(reportIndex,1:size(ompTemperatureBrillance)) = ompTemperatureBrillance
      satelliteScanPositionTmp(reportIndex,1:size(satelliteScanPosition))   = satelliteScanPosition 
      obsQcFlag1Tmp(reportIndex,1:size(obsQcFlag1,1),1:size(obsQcFlag1,2))            = obsQcFlag1
      obsQcFlag2Tmp(reportIndex,1:size(obsQcFlag2))              = obsQcFlag2 
      obsChannelsTmp(reportIndex,1:size(obsChannels))             = obsChannels
      ompChannelsTmp(reportIndex,1:size(ompChannels))             = ompChannels
      obsFlagsTmp(reportIndex,1:size(obsFlags))                = obsFlags
      satOrbitTmp(reportIndex,1:size(satOrbit))                = satOrbit
      obsGlobalMarkerTmp(reportIndex,1:size(obsGlobalMarker))         = obsGlobalMarker
    end if
    if (ifLastReport) exit REPORTS
  end do REPORTS
  write(*,*) ' Number of obs pts read from BURP file              = ', nobs_tot
  write(*,*) ' Number of BURP file reports                        = ', reportIndex
  write(*,*) ' Number of bad BURP file reports (all data flagged) = ', n_bad_reps 
  write(*,*) 'liste of resume report = ', listeOfResumeReport
  write (*,*) 'SAT IDENTIFIER = ', satIdentifierTmp
  ! QC LOOP 
  reportIndex = 0
  ifLastReport = .false.
  REPORTS2: do 
    reportIndex = reportIndex + 1
    ifLastReport = reportIndex == lastReportIndex
    if (.not. listeOfResumeReport(reportIndex)) then
      if (.not. listeOfGoodreport(reportIndex)) then
        write(*,*) 'Bad report index'
        cycle REPORTS2
      else
        !###############################################################################
        ! STEP 3) trouver l'indice du satellite                                        !
        !###############################################################################
        write(*,*) ' ==> mwbg_findSatelliteIndex: '
        call mwbg_findSatelliteIndex(burpFileSatId, satelliteId, satIndexObserrFile)
    
        !###############################################################################
        ! STEP 4) Interpolation de le champ MX(topogrpahy), MG et GL aux pts TOVS.
        !###############################################################################
        write(*,*) ' ==> mwbg_readGeophysicFieldsAndInterpolate: '
        call mwbg_readGeophysicFieldsAndInterpolate(instName, obsLatitudeTmp(reportIndex,:), &
                                                    obsLongitudeTmp(reportIndex,:), modelInterpTerrain, &
                                                    modelInterpGroundIce, modelInterpSeaIce)

        !###############################################################################
        ! STEP 5) Controle de qualite des TOVS. Data QC flags (obsFlags) are modified here!
        !###############################################################################
        write(*,*) ' ==> mwbg_tovCheck For: ', instName
        if (instName == 'AMSUA') then
          call mwbg_tovCheckAmsua(TOVERRST, IUTILST, satIdentifierTmp(reportIndex,:), landQualifierIndiceTmp(reportIndex,:), &
                                  satOrbitTmp(reportIndex,:), obsChannelsTmp(reportIndex,:), ompChannelsTmp(reportIndex,:), &
                                  obsTemperatureBrillanceTmp(reportIndex,:), biasCorrTmp(reportIndex,:), &
                                  ompTemperatureBrillanceTmp(reportIndex,:), qcIndicator, reportNumChannelTmp(reportIndex), &
                                  reportNumObsTmp(reportIndex), mwbg_realMisg, satIndexObserrFile, globalQcIndicator, &
                                  satelliteScanPositionTmp(reportIndex,:), modelInterpGroundIce, modelInterpTerrain, modelInterpSeaIce, &
                                  terrainTypeIndiceTmp(reportIndex,:), satZenithAngleTmp(reportIndex,:), obsFlagsTmp(reportIndex,:), &
                                  newInformationFlag, cloudLiquidWaterPath, atmScatteringIndex, rejectionCodArray, burpFileSatId, &
                                  RESETQC, obsLatitudeTmp(reportIndex,:))
        else if (instName == 'ATMS') then
          call mwbg_tovCheckAtms(TOVERRST, IUTILST, mglg_file, obsLatitudeTmp(reportIndex,:), obsLongitudeTmp(reportIndex,:), &
                                 landQualifierIndiceTmp(reportIndex,:), terrainTypeIndiceTmp(reportIndex,:), &
                                 satZenithAngleTmp(reportIndex,:), obsQcFlag2Tmp(reportIndex,:), obsQcFlag1Tmp(reportIndex,:,:), &
                                 satIdentifierTmp(reportIndex,:), satOrbitTmp(reportIndex,:), obsChannelsTmp(reportIndex,:), &
                                 ompChannelsTmp(reportIndex,:), obsTemperatureBrillanceTmp(reportIndex,:), biasCorrTmp(reportIndex,:), &
                                 ompTemperatureBrillanceTmp(reportIndex,:), qcIndicator, reportNumChannelTmp(reportIndex), &
                                 reportNumChannelTmp(reportIndex), reportNumObsTmp(reportIndex), mwbg_realMisg, satIndexObserrFile, &
                                 newInformationFlag, globalQcIndicator, satelliteScanPositionTmp(reportIndex,:), modelInterpTerrain, &
                                 obsGlobalMarkerTmp(reportIndex,:), obsFlagsTmp(reportIndex,:), cloudLiquidWaterPath, atmScatteringIndex, &
                                 rejectionCodArray, rejectionCodArray2, burpFileSatId, RESETQC, ifLastReport)
        else
          write(*,*) 'midas-bgckMW: instName = ', instName
          call utl_abort('midas-bgckMW: unknown instName')
        end if

        !###############################################################################
        ! STEP 6) Accumuler Les statistiques sur les rejets
        !###############################################################################
        write(*,*) ' ==> mwbg_qcStats For: ', instName
        call mwbg_qcStats(instName, qcIndicator, obsChannelsTmp(reportIndex,:), satIndexObserrFile, reportNumChannel, &
                          reportNumObsTmp(reportIndex), satelliteId, .FALSE., rejectionCodArray, rejectionCodArray2)

      end if 
    end if
    !###############################################################################
    ! STEP 7) Update the burpfile out burpFileNameIn
    !###############################################################################
    write(*,*) ' ==> mwbg_updateBurp For : ', instName
    call mwbg_updateBurp(burpFileNameIn, ReportIndex, ETIKRESU, obsTemperatureBrillanceTmp(reportIndex,:), cloudLiquidWaterPath, &
                         atmScatteringIndex, newInformationFlag, obsGlobalMarkerTmp(reportIndex,:), RESETQC, globalQcIndicator, &
                         landQualifierIndiceTmp(reportIndex,:), terrainTypeIndiceTmp(reportIndex,:), &
                         obsFlagsTmp(reportIndex,:), writeTbValuesToFile, writeModelLsqTT, writeEle25174, burpFileNameout)

    if (ifLastReport) exit REPORTS2

  end do REPORTS2

  write(*,*) ' Number of obs pts read from BURP file              = ', nobs_tot
  write(*,*) ' Number of BURP file reports                        = ', reportIndex
  write(*,*) ' Number of bad BURP file reports (all data flagged) = ', n_bad_reps

  !###############################################################################
  ! STEP 8) Print the statistics in listing file 
  !###############################################################################
  call mwbg_qcStats(instName, qcIndicator, obsChannelsTmp(reportIndex,:), satIndexObserrFile, reportNumChannel, &
                    reportNumObs, satelliteId, .TRUE., rejectionCodArray, rejectionCodArray2)

  call tmg_stop(1)
  call tmg_terminate(mpi_myid, 'TMG_BGCKMW' )

  call rpn_comm_finalize(ier)

  istat  = exfin('midas-bgckMW','FIN','NON')

end program midas_bgckMW

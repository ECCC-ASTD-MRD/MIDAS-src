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
  use obsSpaceData_mod

  implicit none
  
  type(struct_obs)              :: obsSpaceData                                       ! ObsSpace Data object
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
  real,    allocatable          :: obsTb(:)                         ! temperature de brillance
  real,    allocatable          :: ompTb(:)                         ! o-p temperature de "
  real,    allocatable          :: obsTbBiasCorr(:)                                        ! bias correction fo obsTb
  integer, allocatable          :: satScanPosition(:)                           ! scan position
  integer, allocatable          :: obsQcFlag1(:,:)                                    ! Obs Quality flag 1
  integer, allocatable          :: obsQcFlag2(:)                                      ! Obs Quality flag 2 
  integer, allocatable          :: obsChannels(:)                                     ! obsTb channels
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
  integer, allocatable          :: obsNumPerReport(:)                                 ! number of obs per report
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

  ! reading nammwobs namelist
  nulnam = 0
  ier = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read(nulnam, nml=nammwobs, iostat=ier)
  if ( ier /= 0 ) then 
    call utl_abort('midas_bgckmw: Error reading nammwobs namelist')
  end if 
  write(*,nml=nammwobs)
  ier = fclos(nulnam)
  
  ! Initialize obsSpaceData object
  call obs_class_initialize('ALL')
  call obs_initialize( obsSpaceData, mpi_local=.true. )

  ! Initializations of counters (for total reports/locations in the file).
  n_bad_reps = 0
  nobs_tot = 0 
  
  !#################################################################################
  ! Lecture des statistiques d'erreur totale pour les  TOVS 
  !#################################################################################
  call mwbg_readStatTovs(statsFile, instName, satelliteId, IUTILST, TOVERRST)

  ! Allocation of the arrays to use for ObsSpacedata 
  call utl_reAllocate(obsNumPerReport,reportNumMax)
  call utl_reAllocate(listeOfResumeReport,reportNumMax)
  call utl_reAllocate(listeOfGoodReport,reportNumMax)

  ! MAIN LOOP through all the reports in the BURP file
  listeOfGoodReport(:) = .true.
  listeOfResumeReport(:) = .false.
  reportIndex = 0

  REPORTS: do 
    reportIndex = reportIndex + 1
    !###############################################################################
    ! STEP 0 ) Read TOV Observations in big arrays
    !###############################################################################
    write(*,*) ' ==> mwbg_getData: '
    call mwbg_getData(burpFileNameIn, reportIndex, satIdentifier, satZenithAngle,   &
                      landQualifierIndice, terrainTypeIndice, obsLatitude,          &
                      obsLongitude, obsTb, obsTbBiasCorr, ompTb, satScanPosition,   &
                      reportNumChannel, reportNumObs, obsQcFlag1, obsQcFlag2,       &
                      obsChannels, ompChannels, obsFlags, satOrbit, obsGlobalMarker,&
                      resumeReport, ifLastReport, instName, burpFileSatId)

    if (resumeReport) listeOfResumeReport(reportIndex) = .true.
    if (ifLastReport) lastReportIndex = reportIndex
    if (.not. resumeReport) then
      if (ALL(ompTb(:) == MPC_missingValue_R4)) then 
        n_bad_reps = n_bad_reps + 1  
        listeOfGoodReport(reportIndex) = .false.
        cycle REPORTS
      end if
      ! fill array with number of location per report
      obsNumPerReport(reportIndex) = reportNumObs
      ! Increment total number of obs pts read
      nobs_tot = nobs_tot + reportNumObs
      ! Put vectors in ObsSpaceData
      write(*,*) ' ==> mwbg_copyObsToObsSpace: '
      call mwbg_copyObsToObsSpace(instName, reportNumObs, reportNumChannel,         &
                                  satIdentifier, satZenithAngle,landQualifierIndice,&
                                  terrainTypeIndice, obsLatitude, obsLongitude,     &
                                  satScanPosition, obsQcFlag1, satOrbit,            &
                                  obsGlobalMarker, burpFileSatId, obsTb,            &
                                  obsTbBiasCorr, ompTb, obsQcFlag2, obsChannels,    &
                                  ompChannels, obsFlags, obsSpaceData)

    end if
    if (ifLastReport) exit REPORTS
  end do REPORTS
  
  ! QC LOOP 
  reportIndex = 0
  ifLastReport = .false.

  QCLoop: do 
    reportIndex = reportIndex + 1
    ifLastReport = reportIndex == lastReportIndex
    if (listeOfGoodreport(reportIndex)) then
      if (.not. listeOfResumeReport(reportIndex)) then
        ! read burp arrays from ObspaceData Object
        write(*,*) ' ==> mwbg_readObsFromObsSpace: '
        call mwbg_readObsFromObsSpace(instName, obsNumPerReport(reportIndex),       &
                                 satIdentifier, satZenithAngle,landQualifierIndice, & 
                                 terrainTypeIndice, obsLatitude, obsLongitude,      &
                                 satScanPosition, obsQcFlag1, satOrbit,             &
                                 obsGlobalMarker, burpFileSatId, obsTb,             &
                                 obsTbBiasCorr, ompTb, obsQcFlag2, obsChannels,     &
                                 ompChannels, obsFlags, reportNumObs,               &
                                 reportNumChannel, obsSpaceData)

        !###############################################################################
        ! STEP 1) trouver l'indice du satellite                                        !
        !###############################################################################
        write(*,*) ' ==> mwbg_findSatelliteIndex: '
        call mwbg_findSatelliteIndex(burpFileSatId, satelliteId, satIndexObserrFile)
    
        !###############################################################################
        ! STEP 2) Interpolation de le champ MX(topogrpahy), MG et GL aux pts TOVS.
        !###############################################################################
        write(*,*) ' ==> mwbg_readGeophysicFieldsAndInterpolate: '
        call mwbg_readGeophysicFieldsAndInterpolate(instName, obsLatitude,            &
                                                obsLongitude, modelInterpTerrain,     &
                                                modelInterpGroundIce, modelInterpSeaIce)

        !###############################################################################
        ! STEP 3) Controle de qualite des TOVS. Data QC flags (obsFlags) are modified here!
        !###############################################################################
        write(*,*) ' ==> mwbg_tovCheck For: ', instName
        if (instName == 'AMSUA') then
          call mwbg_tovCheckAmsua(TOVERRST, IUTILST, satIdentifier, landQualifierIndice,&
                              satOrbit, obsChannels, ompChannels, obsTb, obsTbBiasCorr, & 
                              ompTb, qcIndicator, reportNumChannel, reportNumObs,       &
                              mwbg_realMisg, satIndexObserrFile, globalQcIndicator,     &
                              satScanPosition, modelInterpGroundIce, modelInterpTerrain,&
                              modelInterpSeaIce, terrainTypeIndice, satZenithAngle,     &
                              obsFlags, newInformationFlag, cloudLiquidWaterPath,       &
                              atmScatteringIndex, rejectionCodArray, burpFileSatId,     &
                              RESETQC, obsLatitude)
        else if (instName == 'ATMS') then
          call mwbg_tovCheckAtms(TOVERRST, IUTILST,mglg_file, obsLatitude, obsLongitude,&
                              landQualifierIndice, terrainTypeIndice, satZenithAngle,   &
                              obsQcFlag2, obsQcFlag1, satIdentifier, satOrbit,          &
                              obsChannels, ompChannels, obsTb, obsTbBiasCorr, ompTb,    &
                              qcIndicator, reportNumChannel, reportNumChannel,          &
                              reportNumObs, mwbg_realMisg, satIndexObserrFile,          &
                              newInformationFlag, globalQcIndicator, satScanPosition,   &
                              modelInterpTerrain, obsGlobalMarker, obsFlags,            &
                              cloudLiquidWaterPath,atmScatteringIndex,rejectionCodArray,&
                              rejectionCodArray2, burpFileSatId, RESETQC)
        else
          write(*,*) 'midas-bgckMW: instName = ', instName
          call utl_abort('midas-bgckMW: unknown instName')
        end if

        !###############################################################################
        ! STEP 4) Accumuler Les statistiques sur les rejets
        !###############################################################################
        write(*,*) ' ==> mwbg_qcStats For: ', instName
        call mwbg_qcStats(instName, qcIndicator, obsChannels, satIndexObserrFile,       &
                          reportNumChannel, reportNumObs, satelliteId, .FALSE.,         &
                          rejectionCodArray, rejectionCodArray2)
      end if

    !###############################################################################
    ! STEP 5) Update the burpfile out burpFileNameIn
    !###############################################################################
    write(*,*) ' ==> mwbg_updateBurp For : ', instName
    call mwbg_updateBurp(burpFileNameIn,reportIndex,ETIKRESU,obsTb,cloudLiquidWaterPath,&
                         atmScatteringIndex,newInformationFlag,obsGlobalMarker, RESETQC,& 
                         globalQcIndicator, landQualifierIndice, terrainTypeIndice,     &
                         obsFlags, writeTbValuesToFile, writeModelLsqTT, writeEle25174, &
                         burpFileNameout)
    end if
    if (ifLastReport) exit QCLoop
  end do QCLoop

  write(*,*) ' Number of obs pts read from BURP file              = ', nobs_tot
  write(*,*) ' Number of BURP file reports                        = ', reportIndex
  write(*,*) ' Number of bad BURP file reports (all data flagged) = ', n_bad_reps

  !###############################################################################
  ! Print the statistics in listing file 
  !###############################################################################
  call mwbg_qcStats(instName, qcIndicator, obsChannels, satIndexObserrFile,              &
                    reportNumChannel, reportNumObs, satelliteId,.TRUE.,rejectionCodArray,&
                    rejectionCodArray2)

  ! deasllocation
  deallocate(obsNumPerReport)
  deallocate(listeOfResumeReport)
  deallocate(listeOfGoodReport)
  ! deallocate obsSpaceData object
  call obs_finalize(obsSpaceData) ! deallocate obsSpaceData

  call tmg_stop(1)
  call tmg_terminate(mpi_myid, 'TMG_BGCKMW' )
  call rpn_comm_finalize(ier)
  istat  = exfin('midas-bgckMW','FIN','NON')

end program midas_bgckMW

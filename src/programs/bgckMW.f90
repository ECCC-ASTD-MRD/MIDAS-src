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
  integer                       :: reportNumObs                                       ! number of obs in current report
  integer                       :: reportNumChannel                                   ! "      "   channels "      "
  integer                       :: nobs_tot                                           ! obs. total number
  integer                       :: n_bad_reps                                         ! bad reports number
  integer                       :: satelliteIndexObserrorFile                         ! satellite index in obserror file
  logical                       :: resumeReport                                       ! logical to newInformationFlagify resume report
  logical                       :: ifLastReport                                       ! logical to newInformationFlagify last report
  character(len=9)              :: burpFileSatId                                      ! station id in burp file
  character(len=9), allocatable :: satelliteId(:)                                     ! satellite Id in stats error file
  real, allocatable             :: modelInterpTerrain(:)                              ! topo in standard file interpolated to obs point
  real, allocatable             :: modelInterpSeaIce(:)                               ! Glace de mer " "
  real, allocatable             :: modelInterpGroundIce(:)                            ! Glace de continent " "
  real,    allocatable          :: obsLatitude(:)                                     ! obs. point latitudes
  real,    allocatable          :: obsLongitude(:)                                    ! obs. point longitude
  integer, allocatable          :: satellitenewInformationFlagifier(:)                ! Satellite index in burp file
  real,    allocatable          :: satelliteZenithAngle(:)                            ! sat. satelliteZenithAngle angle
  integer, allocatable          :: landQualifierIndice(:)                             ! land qualifyer
  integer, allocatable          :: terrainTypeIndice(:)                               ! terrain type
  real,    allocatable          :: obsTemperatureBrillance(:)                         ! temperature de brillance
  real,    allocatable          :: btClear(:)                                         ! clear-sky observed radiance
  real,    allocatable          :: ompTemperatureBrillance(:)                         ! o-p temperature de "
  real,    allocatable          :: biasCorr(:)                                        ! bias correction fo obsTemperatureBrillance
  integer, allocatable          :: satelliteScanPosition(:)                           ! scan position
  integer, allocatable          :: obsQcFlag1(:,:)                                    ! Obs Quality flag 1
  integer, allocatable          :: obsQcFlag2(:)                                      ! Obs Quality flag 2 
  integer, allocatable          :: observationChannels(:)                             ! obsTemperatureBrillance channels
  integer, allocatable          :: ompChannels(:)                                     ! zomp channel
  integer, allocatable          :: observationFlags(:)                                ! obs. flag
  integer, allocatable          :: satelliteOrbitnewInformationFlagifier(:)           ! orbit
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
  integer, allocatable          :: qualityControlIndicator(:,:)                     ! indicateur controle de qualite tovs par canal 
  !                                                                =0, ok,
  !                                                                >0, rejet,
  integer, allocatable          :: globalQcIndicator(:)                      ! indicateur global controle de qualite tovs. Code:
  !                                                                =0, ok,
  !                                                                >0, rejet d'au moins un canal
  integer, allocatable          :: newInformationFlag(:)                        ! ATMS Information flag (newInformationFlag) values (new BURP element 
  !                                                                  025174 in header). FOR AMSUA just fill with zeros 
  real,    allocatable          :: cloudLiquidWaterPathObs(:)                       ! cloud liquid water path from observation.
  real,    allocatable          :: cloudLiquidWaterPathFG(:)                        ! cloud liquid water path from background.
  real,    allocatable          :: atmScatteringIndex(:)                        ! scattering index
  integer, external             :: exdb, exfin, fnom, fclos
  integer                       :: ier, istat, nulnam

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
  logical                       :: useAveragedClwForQC

  namelist /nambgck/instName, burpFileNameIn, burpFileNameOut, mglg_file, statsFile, &
                    writeModelLsqTT, writeEle25174, clwQcThreshold, allowStateDepSigmaObs, &
                    useUnbiasedObsForClw, debug, RESETQC, ETIKRESU, writeTbValuesToFile,
                    useAveragedClwForQC 

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
 
  ! default namelist values
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
  useAveragedClwForQC   = .false.

  ! reading namelist
  nulnam = 0
  ier = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read(nulnam, nml=nambgck, iostat=ier)
  if ( ier /= 0 ) then 
    call utl_abort('midas_bgckmw: Error reading namelist')
  end if 
  write(*,nml=nambgck)
  ier = fclos(nulnam)

  mwbg_debug = debug
  mwbg_clwQcThreshold = clwQcThreshold 
  mwbg_allowStateDepSigmaObs = allowStateDepSigmaObs
  mwbg_useUnbiasedObsForClw = useUnbiasedObsForClw
  mwbg_useAveragedClwForQC = useAveragedClwForQC

  if ( (      mwbg_allowStateDepSigmaObs .and. .not. mwbg_useAveragedClwForQC) .or. &
       (.not. mwbg_allowStateDepSigmaObs .and.       mwbg_useAveragedClwForQC) ) then
    write(*,*) 'Using state-dependant obs and averaged CLW for QC should be compatible and are not.'
    call abort()
  end if

  ! Initializations of counters (for total reports/locations in the file).
  n_bad_reps = 0
  nobs_tot = 0 
  
  !#################################################################################
  ! STEP 1 ) Lecture des statistiques d'erreur totale pour les  TOVS 
  !#################################################################################
  call mwbg_readStatTovs(statsFile, instName, satelliteId, IUTILST, TOVERRST)

  ! MAIN LOOP through all the reports in the BURP file
  reportIndex = 0
  REPORTS: do 

    reportIndex = reportIndex + 1
    !###############################################################################
    ! STEP 2 ) Lecture des observations TOVS dans le fichier burpFileNameIn
    !###############################################################################
    write(*,*) ' ==> mwbg_getData: '
    call mwbg_getData(burpFileNameIn, reportIndex, satellitenewInformationFlagifier, satelliteZenithAngle, landQualifierIndice, &
                      terrainTypeIndice, obsLatitude, obsLongitude, obsTemperatureBrillance, btClear, biasCorr, ompTemperatureBrillance, &
                      satelliteScanPosition, reportNumChannel, reportNumObs, obsQcFlag1, obsQcFlag2, observationChannels, ompChannels, &
                      observationFlags, satelliteOrbitnewInformationFlagifier, obsGlobalMarker, resumeReport, ifLastReport, &
                      instName, burpFileSatId)

    if (.not. resumeReport) then
      if (ALL(ompTemperatureBrillance(:) == MPC_missingValue_R4)) then 
        n_bad_reps = n_bad_reps + 1  
        cycle REPORTS
      end if
      ! Increment total number of obs pts read
      nobs_tot = nobs_tot + reportNumObs

      !###############################################################################
      ! STEP 3) trouver l'indice du satellite                                        !
      !###############################################################################
      write(*,*) ' ==> mwbg_findSatelliteIndex: '
      call mwbg_findSatelliteIndex(burpFileSatId, satelliteId, satelliteIndexObserrorFile)
    
      !###############################################################################
      ! STEP 4) Interpolation de le champ MX(topogrpahy), MG et GL aux pts TOVS.
      !###############################################################################
      write(*,*) ' ==> mwbg_readGeophysicFieldsAndInterpolate: '
      call mwbg_readGeophysicFieldsAndInterpolate(instName, obsLatitude, obsLongitude, modelInterpTerrain, modelInterpGroundIce, &
                                                  modelInterpSeaIce)

      !###############################################################################
      ! STEP 5) Controle de qualite des TOVS. Data QC flags (observationFlags) are modified here!
      !###############################################################################
      write(*,*) ' ==> mwbg_tovCheck For: ', instName
      if (instName == 'AMSUA') then
        call mwbg_tovCheckAmsua(TOVERRST, IUTILST, satellitenewInformationFlagifier, landQualifierIndice, &
                                satelliteOrbitnewInformationFlagifier, observationChannels, ompChannels, obsTemperatureBrillance, btClear, &
                                biasCorr, ompTemperatureBrillance, qualityControlIndicator, reportNumChannel, reportNumObs, &
                                mwbg_realMisg, satelliteIndexObserrorFile, globalQcIndicator, satelliteScanPosition, &
                                modelInterpGroundIce, modelInterpTerrain, modelInterpSeaIce, terrainTypeIndice, satelliteZenithAngle, &
                                observationFlags, newInformationFlag, cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, atmScatteringIndex, rejectionCodArray, &
                                burpFileSatId, RESETQC, obsLatitude)
      else if (instName == 'ATMS') then
        call mwbg_tovCheckAtms(TOVERRST, IUTILST, mglg_file, obsLatitude, obsLongitude, landQualifierIndice, terrainTypeIndice, &
                               satelliteZenithAngle, obsQcFlag2, obsQcFlag1, satellitenewInformationFlagifier, &
                               satelliteOrbitnewInformationFlagifier, observationChannels, ompChannels, obsTemperatureBrillance, &
                               biasCorr, ompTemperatureBrillance, qualityControlIndicator, reportNumChannel, reportNumChannel, &
                               reportNumObs, mwbg_realMisg, satelliteIndexObserrorFile, newInformationFlag, globalQcIndicator, &
                               satelliteScanPosition, modelInterpTerrain, obsGlobalMarker, observationFlags, cloudLiquidWaterPathObs, &
                               atmScatteringIndex, rejectionCodArray, rejectionCodArray2, burpFileSatId, RESETQC, ifLastReport)
      else
        write(*,*) 'midas-bgckMW: instName = ', instName
        call utl_abort('midas-bgckMW: unknown instName')
      end if

      !###############################################################################
      ! STEP 6) Accumuler Les statistiques sur les rejets
      !###############################################################################
      write(*,*) ' ==> mwbg_qcStats For: ', instName
      call mwbg_qcStats(instName, qualityControlIndicator, observationChannels, satelliteIndexObserrorFile, reportNumChannel, &
                        reportNumObs, satelliteId, .FALSE., rejectionCodArray, rejectionCodArray2)

    end if 

    !###############################################################################
    ! STEP 7) Update the burpfile out burpFileNameIn
    !###############################################################################
    write(*,*) ' ==> mwbg_updateBurp For : ', instName
    call mwbg_updateBurp(burpFileNameIn, ReportIndex, ETIKRESU, obsTemperatureBrillance, cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, atmScatteringIndex, &
                         newInformationFlag, obsGlobalMarker, RESETQC, globalQcIndicator, landQualifierIndice, terrainTypeIndice, &
                         observationFlags, writeTbValuesToFile, writeModelLsqTT, writeEle25174, burpFileNameout)

    if (ifLastReport) exit REPORTS

  end do REPORTS

  write(*,*) ' Number of obs pts read from BURP file              = ', nobs_tot
  write(*,*) ' Number of BURP file reports                        = ', reportIndex
  write(*,*) ' Number of bad BURP file reports (all data flagged) = ', n_bad_reps

  !###############################################################################
  ! STEP 8) Print the statistics in listing file 
  !###############################################################################
  call mwbg_qcStats(instName, qualityControlIndicator, observationChannels, satelliteIndexObserrorFile, reportNumChannel, &
                    reportNumObs, satelliteId, .TRUE., rejectionCodArray, rejectionCodArray2)

  call tmg_stop(1)
  call tmg_terminate(mpi_myid, 'TMG_BGCKMW' )

  call rpn_comm_finalize(ier)

  istat  = exfin('midas-bgckMW','FIN','NON')

end program midas_bgckMW

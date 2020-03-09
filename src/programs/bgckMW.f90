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
  !:Purpose: Main program for background check of microwave instruments. 
  !
  use bgckmicrowave_mod
  use MathPhysConstants_mod

  implicit none
  integer, parameter            :: nchanAtms=22                    ! ATMS channel number        
  integer, parameter            :: mxscan=96                       ! 
  integer, PARAMETER            :: MXSAT = 9                       ! max number of satellite
  integer, parameter            :: JPNSAT = 9                      ! 
  integer, parameter            :: JPCH = 50                       ! missing value
  real, parameter               :: zmisg=9.9e09                    !
  integer, parameter            :: MXCHN = 42                      !
  integer, parameter            :: JPMXREJ = 15                    ! test index
  integer, parameter            :: ipc=6                           ! 
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
  logical                       :: lutb                            ! logical for replacing missing tb value
  integer                       :: nobs_tot                        ! obs. total number
  integer                       :: n_bad_reps                      ! bad reports number
  integer                       :: reportIndex                     ! report index
  integer                       :: ntOut                           ! number of obs in current report
  integer                       :: nvalOut                         ! "      "   channels "      "
  integer                       :: exdb                            !
  integer                       :: exfin                           !
  integer                       :: ier                             !
  integer                       :: istat                           !
  integer                       :: nulnam                          !
  integer                       :: fnom                            !
  integer                       :: fclos                           !
  integer                       :: INOSAT                          !
  logical                       :: resumeReport                    ! logical to identify resume report
  logical                       :: ifLastReport                    ! logical to identify last report
  character(len=9)              :: STNID                           ! station id in burp file
  character(len=9), allocatable :: satelliteId(:)                  ! satellite Id in stats error file
  real, allocatable             :: MTINTRP(:)                      ! topo in standard file interpolated to obs point
  real, allocatable             :: GLINTRP(:)                      ! Glace de mer " "
  real, allocatable             :: MGINTRP(:)                      ! Glace de continent " "
  real,    allocatable          :: zlat(:)                         ! obs. point latitudes
  real,    allocatable          :: zlon(:)                         ! obs. point longitude
  integer, allocatable          :: ISAT(:)                         ! Satellite index in burp file
  real,    allocatable          :: zenith(:)                       ! sat. zenith angle
  integer, allocatable          :: ilq(:)                          ! land qualifyer
  integer, allocatable          :: itt(:)                          ! terrain type
  real,    allocatable          :: ztb(:)                          ! temperature de brillance
  real,    allocatable          :: ZOMP(:)                         ! o-p temperature de "
  real,    allocatable          :: biasCorr(:)                     ! bias correction fo ztb
  integer, allocatable          :: scanpos(:)                      ! scan position
  integer, allocatable          :: qcflag1(:,:)                    ! 
  integer, allocatable          :: qcflag2(:)                      ! 
  integer, allocatable          :: ican(:)                         ! ztb channels
  integer, allocatable          :: icanomp(:)                      ! zomp channel
  integer, allocatable          :: IMARQ(:)                        ! obs. flag
  integer, allocatable          :: IORBIT(:)                       ! orbit
  integer, allocatable          :: globMarq(:)                     ! global marker
  integer                       :: IUTILST(JPCH,JPNSAT)            ! channel use option for each sat.
  real                          :: TOVERRST(JPCH,JPNSAT)           ! obs. error per channel for each sat
  integer                       :: MREJCOD(JPMXREJ,MXCHN,MXSAT)    ! number of rejection per sat. per channl per test
  integer                       :: MREJCOD2(JPMXREJ,MXCHN,MXSAT)   ! number of rejection per sat. per channl per test
  !                                                                  for ATMS 2nd category of tests
  integer, allocatable          :: lsq(:)                          ! land sea qualifyer
  integer, allocatable          :: trn(:)                          ! terrain type
  integer, allocatable          :: icheck(:,:)                     ! indicateur controle de qualite tovs par canal 
  !                                                                =0, ok,
  !                                                                >0, rejet,
  integer, allocatable          :: ichkprf(:)                      ! indicateur global controle de qualite tovs. Code:
  !                                                                =0, ok,
  !                                                                >0, rejet d'au moins un canal
  integer, allocatable          :: ident(:)                        ! ATMS Information flag (ident) values (new BURP element 
  !                                                                  025174 in header). FOR AMSUA just fill with zeros 
  real,    allocatable          :: clw(:)                          ! cloud liquid water. NB: for AMSUA, clw=0.5(model_clw + obs_clw)
  real,    allocatable          :: scatw(:)                        ! scatering index
  integer                       :: seaIcePointNum                  ! Number of waterobs points converted to sea ice points
  integer                       :: clwMissingPointNum              ! Number of points where CLW/SI missing over water due bad data
  integer                       :: n_reps_tb2misg                  ! Number of BURP file reports where Tb set to zmisg
  integer                       :: cldcnt                          ! 
  integer                       :: flgcnt                          ! Total number of filtered obs pts
  integer                       :: rejcnt                          ! Number of problem obs pts (Tb err, QCfail)
  integer                       :: landcnt                         ! Number of obs pts found over land/ice
  integer                       :: iwvcnt                          ! Number of pts with Mean 183 Ghz Tb < 240K
  integer                       :: pcpcnt                          ! Number of scatter/precip obs
  integer                       :: drycnt                          ! Number of pts flagged for AMSU-B Dryness Index

  external EXDB,EXFIN
  
  namelist /nambgck/instName, burpFileNameIn, burpFileNameOut, mglg_file, statsFile, &
                    writeModelLsqTT, writeEle25174, clwQcThreshold, allowStateDepSigmaObs, &
                    useUnbiasedObsForClw, debug, RESETQC, ETIKRESU 

  istat = exdb('midas-bgckMW','DEBUT','NON')

  write(*,'(/,' //                                                &
            '3(" *****************"),/,' //                       &
            '14x,"-- START OF MAIN PROGRAM midas-bgckMW: --",/,' //   &
            '14x,"-- BACKGROUND CHECK FOR MW OBSERVATIONS --",/, ' //&
            '14x,"-- Revision : ",a," --",/,' //       &
            '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  ! Initializations of counters (for total reports/locations in the file).
  lutb = .False.
  flgcnt = 0
  landcnt = 0
  rejcnt = 0
  cldcnt = 0
  iwvcnt = 0
  pcpcnt = 0
  drycnt = 0
  seaIcePointNum = 0
  clwMissingPointNum = 0
  n_bad_reps = 0
  nobs_tot = 0 
  n_reps_tb2misg = 0 
 
  ! default namelist values
  instName = 'AMSUA'
  burpFileNameIn = './obsto_amsua'
  burpFileNameOut = './obsto_amsua.out'
  mglg_file  = './fstmglg'
  statsFile = './stats_amsua_assim'  
  writeModelLsqTT = .false.
  writeEle25174 = .false.
  clwQcThreshold = 0.3 
  allowStateDepSigmaObs = .false.
  useUnbiasedObsForClw = .false.
  debug = .false.
  RESETQC = .FALSE.
  ETIKRESU = '>>BGCKALT'

  ! reading namelist
  nulnam = 0
  ier = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read(nulnam, nml=nambgck, iostat=ier)
  if ( ier /= 0 ) then 
    write(*,*) 'midas_bgckmw: Error reading namelist'
    call abort()
  end if 
  write(*,nml=nambgck)
  ier = fclos(nulnam)

  mwbg_debug = debug
  mwbg_clwQcThreshold = clwQcThreshold 
  mwbg_allowStateDepSigmaObs = mwbg_allowStateDepSigmaObs
  mwbg_useUnbiasedObsForClw = useUnbiasedObsForClw
 
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
    call mwbg_getData(burpFileNameIn, reportIndex, ISAT, zenith, ilq, itt, &
                      zlat, zlon, ztb, biasCorr, ZOMP, scanpos, nvalOut, &
                      ntOut, qcflag1, qcflag2, ican, icanomp, IMARQ, IORBIT, &
                      globMarq, resumeReport, ifLastReport, instName, STNID)

    if (.not. resumeReport) then
      if (ALL(ZOMP(:) == MPC_missingValue_R4)) then 
        n_bad_reps = n_bad_reps + 1  
        write(*,*) 'Bad Report '
        cycle REPORTS
      end if
      ! Increment total number of obs pts read
      nobs_tot = nobs_tot + ntOut

      !###############################################################################
      ! STEP 3) trouver l'indice du satellite                                        !
      !###############################################################################
      write(*,*) ' ==> mwbg_findSatelliteIndex: '
      call mwbg_findSatelliteIndex(STNID, satelliteId, INOSAT)
    
      !###############################################################################
      ! STEP 4) Interpolation de le champ MX(topogrpahy), MG et GL aux pts TOVS.
      !###############################################################################
      write(*,*) ' ==> mwbg_readGeophysicFieldsAndInterpolate: '
      call mwbg_readGeophysicFieldsAndInterpolate(instName, zlat, zlon, MTINTRP, MGINTRP, &
                                                  GLINTRP)

      !###############################################################################
      ! STEP 5) Controle de qualite des TOVS. Data QC flags (IMARQ) are modified here!
      !###############################################################################
      write(*,*) ' ==> mwbg_tovCheck For: ', instName
      if (instName == 'AMSUA') then
        call mwbg_tovCheckAmsua(TOVERRST, IUTILST, ISAT, ilq, IORBIT, ican, ICANOMP, ztb, &
                                biasCorr, ZOMP, ICHECK, nvalOut, ntOut, ZMISG, INOSAT, &
                                ICHKPRF, scanpos, MGINTRP, MTINTRP, GLINTRP, itt, zenith, &
                                IMARQ, ident, clw, scatw, MREJCOD, STNID, RESETQC, &
                                ZLAT)
      else if (instName == 'ATMS') then
        call mwbg_tovCheckAtms(TOVERRST, IUTILST, mglg_file, zlat, zlon, ilq, itt, lsq, trn, &
                               zenith, qcflag2, qcflag1, ISAT, IORBIT, ICAN, ICANOMP, &
                               ztb, biasCorr, ZOMP, ICHECK, nvalOut, nvalOut, ntOut, zmisg,INOSAT, &
                               ident, ICHKPRF, scanpos, MTINTRP, globMarq, IMARQ, clw, scatw, MREJCOD, &
                               MREJCOD2, STNID, RESETQC,  ipc, n_reps_tb2misg, drycnt,&
                               landcnt, rejcnt, iwvcnt, pcpcnt, flgcnt, cldcnt, seaIcePointNum, &
                               clwMissingPointNum)
      end if

      !###############################################################################
      ! STEP 6) Accumuler Les statistiques sur les rejets
      !###############################################################################
      write(*,*) ' ==> mwbg_qcStats For: ', instName
      call mwbg_qcStats(instName, ICHECK, ican, INOSAT, nvalOut, &
                             ntOut, satelliteId, .FALSE., MREJCOD, MREJCOD2)

    end if 

    !###############################################################################
    ! STEP 7) Update the burpfile out burpFileNameIn
    !###############################################################################
    write(*,*) ' ==> mwbg_updateBurp For : ', instName
    call mwbg_updateBurp(burpFileNameIn, ReportIndex, ETIKRESU, ztb, clw, scatw, ident, &
                           globMarq, RESETQC, ICHKPRF, ilq, itt, IMARQ, lutb, &
                           writeModelLsqTT, writeEle25174, burpFileNameout)

    if (ifLastReport) exit REPORTS

  end do REPORTS

  write(*,*) ' Number of obs pts read from BURP file              = ', nobs_tot
  write(*,*) ' Number of BURP file reports                        = ', reportIndex
  write(*,*) ' Number of bad BURP file reports (all data flagged) = ', n_bad_reps

  !###############################################################################
  ! STEP 8) Print the statistics in listing file 
  !###############################################################################
  call mwbg_qcStats(instName, ICHECK, ican, INOSAT, nvalOut, &
                    ntOut, satelliteId, .TRUE., MREJCOD, MREJCOD2)
  
  istat  = exfin('midas-bgckMW','FIN','NON')

end program midas_bgckMW

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

program midas_bgckAllMW
  !
  ! :Purpose: Main program for background check of microwave instruments. 
  !
  use burp_module
  use bgckmicrowave_mod
  use MathPhysConstants_mod

  implicit none
  

  integer, parameter            :: nchanAtms=22
  integer, parameter            :: mxscan=96
  integer, PARAMETER            :: MXSAT = 9
  integer, parameter            :: JPNSAT = 9
  integer, parameter            :: JPCH = 50
  real, parameter               :: zmisg=9.9e09
  integer, parameter            :: MXCHN = 42   
  integer, parameter            :: JPMXREJ = 15  
  integer, parameter            :: ipc=6

  character(len=128)            :: mglg_file
  character(len=128)            :: statsFile
  logical                       :: lutb
  
  integer                       :: nobs_tot
  integer                       :: n_bad_reps
  integer                       :: n_reps
  integer                       :: reportIndex
  integer                       :: satNumber 
  integer                       :: ntOut 
  integer                       :: nvalOut
  integer                       :: exdb
  integer                       :: exfin
  integer                       :: ier
  integer                       :: istat
  integer                       :: nulnam
  integer                       :: INOSAT
  logical                       :: bad_report
  logical                       :: resumeReport
  logical                       :: ifLastReport
  logical                       :: RESETQC 
  logical                       :: debug
  logical                       :: writeModelLsqTT
  logical                       :: writeEle25174 
  integer, parameter            :: MISGINT = -1
  character(len=9)              :: STNID
  character(len=9), allocatable :: satelliteId(:)
  character(len=5)              :: instName
  character(len=90)             :: burpFileNameIn
  character(len=90)             :: burpFileNameOut
  character(len=9)              :: ETIKRESU
  real, allocatable             :: MTINTRP(:)
  real, allocatable             :: GLINTRP(:)
  real, allocatable             :: MGINTRP(:)
 
  real,    allocatable          :: zlat(:)
  real,    allocatable          :: zlon(:)
  integer, allocatable          :: ISAT(:)
  real,    allocatable          :: zenith(:)
  integer, allocatable          :: ilq(:)
  integer, allocatable          :: itt(:)
  real,    allocatable          :: ztb(:)
  real,    allocatable          :: ZOMP(:)
  real,    allocatable          :: biasCorr(:)
  integer, allocatable          :: scanpos(:)
  integer, allocatable          :: qcflag1(:,:)
  integer, allocatable          :: qcflag2(:)
  integer, allocatable          :: ican(:)
  integer, allocatable          :: icanomp(:)
  integer, allocatable          :: IMARQ(:)
  integer, allocatable          :: IORBIT(:)
  integer, allocatable          :: globMarq(:)
  integer                       :: IUTILST(JPCH,JPNSAT)
  real                          :: TOVERRST(JPCH,JPNSAT)
  integer                       :: MREJCOD(JPMXREJ,MXCHN,MXSAT)
  integer                       :: MREJCOD2(JPMXREJ,MXCHN,MXSAT)
  integer, allocatable          :: lsq(:)
  integer, allocatable          :: trn(:)
  integer, allocatable          :: icheck(:,:)
  integer, allocatable          :: ichkprf(:)
  logical, allocatable          :: waterobs(:)
  logical, allocatable          :: grossrej(:)
  logical, allocatable          :: cloudobs(:)
  logical, allocatable          :: iwvreject(:)
  logical, allocatable          :: precipobs(:)
  integer, allocatable          :: ident(:)
  real,    allocatable          :: zdi(:)
  real,    allocatable          :: clw(:)
  real,    allocatable          :: scatw(:)
  real,    allocatable          :: iwv(:)
  real,    allocatable          :: SeaIce(:)
  real,    allocatable          :: scatec(:)
  real,    allocatable          :: scatbg(:)
  logical, allocatable          :: lflagchn(:,:)
  logical, allocatable          :: lqc(:,:)
  real                          :: glbscanb(nchanAtms,mxscan,mxsat)
  real                          :: dglbscanb(nchanAtms,mxscan,mxsat)
  character(len=9)              :: csatid(mxsat)
  real, allocatable             :: ztbcor(:)
  integer                       :: seaIcePointNum
  integer                       :: clwMissingPointNum
  integer                       :: n_reps_tb2misg

  integer                       :: cldcnt


  integer                       :: ilnmx
  integer                       :: status
  integer                       :: nele
  integer                       :: nt
  integer                       :: blat
  integer                       :: blon
  integer                       :: idtyp
  integer                       :: flgcnt
  integer                       :: rejcnt
  integer                       :: landcnt
  integer                       :: iwvcnt
  integer                       :: pcpcnt
  integer                       :: drycnt
  integer                       :: indx1, indx2, ii, numbsat, iich, jj, kk

  external EXDB,EXFIN
  
  namelist /nambgck/instName, burpFileNameIn, burpFileNameOut, mglg_file, statsFile, &
                    writeModelLsqTT, writeEle25174, debug, RESETQC, ETIKRESU 

  istat = exdb('midas_bgck','DEBUT','NON')


  write(*,'(/,' //                                                &
            '3(" *****************"),/,' //                       &
            '14x,"-- START OF MAIN PROGRAM MIDAS-BGCKMW: --",/,' //   &
            '14x,"-- BACKGROUND CHECK FOR MW OBSERVATIONS --",/, ' //&
            '14x,"-- Revision : ",a," --",/,' //       &
            '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  ! default values
  instName = 'AMSUA'
  burpFileNameIn = './obsto_amsua'
  burpFileNameOut = './obsto_amsua.out'
  mglg_file  = './fstmglg'
  statsFile = './stats_amsua_assim'  
  writeModelLsqTT = .false.
  writeEle25174 = .false.
  debug = .false.
  RESETQC = .FALSE.
  ETIKRESU = '>>BGCKALT'


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
  reportIndex = 0
  nobs_tot = 0 
  n_reps_tb2misg = 0 
  
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

 
  !#################################################################################
  ! STEP 1 ) Lecture des statistiques d'erreur totale pour les  TOVS 
  !#################################################################################
  call mwbg_readStatTovs(statsFile, instName, satelliteId, IUTILST, TOVERRST)
  satNumber = size(satelliteId)

  ! MAIN LOOP through all the reports in the BURP file
  
  REPORTS: do 

    reportIndex = reportIndex + 1

    !###############################################################################
    ! STEP 2 ) Lecture des observations TOVS dans le fichier burpFileNameIn
    !###############################################################################
    call mwbg_getData(burpFileNameIn, reportIndex, ISAT, zenith, ilq, itt, &
                      zlat, zlon, ztb, biasCorr, ZOMP, scanpos, nvalOut, &
                      ntOut, qcflag1, qcflag2, ican, icanomp, IMARQ, IORBIT, &
                      globMarq, resumeReport, ifLastReport, instName, STNID)


    if ( .not. resumeReport ) then
      if ( ALL(ZOMP(:) == MPC_missingValue_R4 )) then 
        n_bad_reps = n_bad_reps + 1  
        write(*,*) 'Bad Report '
        cycle REPORTS
      end if
      ! Increment total number of obs pts read
      nobs_tot = nobs_tot + ntOut

      !###############################################################################
      ! STEP 3) trouver l'indice du satellite                                        !
      !###############################################################################
      call mwbg_findSatelliteIndex(STNID, satelliteId, SatNumber, INOSAT)
    
      !###############################################################################
      ! STEP 4) Interpolation de le champ MX(topogrpahy), MG et GL aux pts TOVS.
      !    N.B.: on examine ce champ sur une boite centree sur chaque obs.
 
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
      call mwbg_qcStats(instName, satNumber, ICHECK, ican, INOSAT, nvalOut, &
                             ntOut, satelliteId, .FALSE., MREJCOD, MREJCOD2)

    end if 

    !###############################################################################
    ! STEP 7) Update the burpfile out burpFileNameIn
    !###############################################################################

    write(*,*) ' ==> mwbg_updateBurp For : ', instName
    
    call mwbg_updateBurp(burpFileNameIn, ReportIndex, ETIKRESU, ztb, clw, scatw, ident, &
                           globMarq, RESETQC, ICHKPRF, ilq, itt, IMARQ, lutb, &
                           writeModelLsqTT, writeEle25174, burpFileNameout)

    if ( ifLastReport) exit REPORTS

  end do REPORTS

  write(*,*) ' Number of obs pts read from BURP file              = ', nobs_tot
  write(*,*) ' Number of BURP file reports                        = ', reportIndex
  write(*,*) ' Number of bad BURP file reports (all data flagged) = ', n_bad_reps

  !###############################################################################
  ! STEP 8) Print the statistics in listing file 
  !###############################################################################

  call mwbg_qcStats(instName, satNumber, ICHECK, ican, INOSAT, nvalOut, &
                         ntOut, satelliteId, .TRUE., MREJCOD, MREJCOD2)
  
  istat  = exfin('bgckMW','FIN','NON')

end program midas_bgckAllMW

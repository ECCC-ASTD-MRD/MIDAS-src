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

program midas_bgckmw
  !
  ! :Purpose: Main program for background check of microwave instruments. 
  !
  use burp_module
  use bgckmicrowave_mod
  use MathPhysConstants_mod

  implicit none

  integer, PARAMETER            :: MXSAT = 9
  integer, PARAMETER            :: MXVAL = 22
  integer, PARAMETER            :: MXNT = 3000
  integer                       :: nobs_tot
  integer                       :: n_bad_reps
  integer                       :: n_reps
  integer                       :: reportIndex
  integer                       :: satNumber 
  integer                       :: ntOut 
  integer                       :: nvalOut
  integer                       :: junk 
  integer                       :: exdb
  integer                       :: exfin
  integer                       :: ier
  integer                       :: I, ISTAT
  integer                       :: nulnam
  integer                       :: INOSAT
  logical                       :: bad_report
  logical                       :: resumeReport
  logical                       :: ifLastReport
  logical                       :: RESETQC 
  logical                       :: debug
  logical                       :: clwQcThreshold
  logical                       :: allowStateDepSigmaObs  
  real, parameter               :: zmisg=9.9e09
  integer, parameter            :: MISGINT = -1
  character(len=9)              :: STNID
  character(len=9), allocatable :: satelliteId(:)
  !character(len=9) :: satelliteId(MXSAT)

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
  real,    allocatable          :: clw_avg(:) 
  real,    allocatable          :: clw(:)
  real,    allocatable          :: scatw(:)
  integer, allocatable          :: icheck(:,:)
  integer, allocatable          :: ichkprf(:)

  EXTERNAL EXDB,EXFIN
  namelist /nambgck/ debug, RESETQC, ETIKRESU, clwQcThreshold

  JUNK = EXDB('BGCKMW','DEBUT','NON')

  write(*,'(/,' //                                                &
            '3(" *****************"),/,' //                       &
            '14x,"-- START OF MAIN PROGRAM MIDAS-BGCKMW: --",/,' //   &
            '14x,"-- BACKGROUND CHECK FOR MW OBSERVATIONS --",/, ' //&
            '14x,"-- Revision : ",a," --",/,' //       &
            '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  ! 1) Debut
  ! default values
  debug = .false.
  RESETQC = .FALSE.
  ETIKRESU = '>>BGCKALT'
  clwQcThreshold = 0.3
  allowStateDepSigmaObs = .false.
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
  mwbg_allowStateDepSigmaObs = allowStateDepSigmaObs

  burpFileNameIn = './obsto_amsua'
  burpFileNameOut = './obsto_amsua.out'


  ! 2) Lecture des statistiques d'erreur totale pour les  TOVS 
  call mwbg_readStatTovs('./stats_amsua_assim', 'AMSUA', satelliteId)
  satNumber = size(satelliteId)

  ! MAIN LOOP through all the reports in the BURP file
  n_bad_reps = 0
  reportIndex = 0
  nobs_tot = 0 
  
  REPORTS: do 

    reportIndex = reportIndex + 1
    !  Get all the required data from the blocks in the report (Rpt_in)
    call mwbg_getData(burpFileNameIn, reportIndex, ISAT, zenith, ilq, itt, &
                      zlat, zlon, ztb, biasCorr, ZOMP, scanpos, nvalOut, &
                      ntOut, qcflag1, qcflag2, ican, icanomp, IMARQ, IORBIT, &
                      globMarq, resumeReport, ifLastReport, 'AMSUA', STNID)


    if ( .not. resumeReport ) then
      if ( ALL(ZOMP(:) == MPC_missingValue_R4 )) then 
        n_bad_reps = n_bad_reps + 1  
        write(*,*) 'Bad Report '
        cycle REPORTS
      end if
      ! Increment total number of obs pts read
      nobs_tot = nobs_tot + ntOut
      ! trouver l'indice du satellite
      INOSAT = 0
      DO I = 1, satNumber
        if ( STNID .EQ. '^'//satelliteId(I) ) THEN
          INOSAT = I
          write(*,*)' SATELLITE = ', STNID
          write(*,*)'    INOSAT = ', INOSAT
        end if
      ENDDO
      if ( INOSAT .EQ. 0 ) THEN
        write(*,*) 'SATELLITE ',TRIM(STNID), &
                   ' NOT FOUND IN STATS FILE!'
        call ABORT()
      end if
    
      write(*,*) ' Interpolation'
      ! 5) Interpolation de la glace et le champ terre/mer du modele aux pts TOVS.
      ! N.B.: on examine ces champs sur une boite centree sur chaque obs.
      call mwbg_readGeophysicFieldsAndInterpolate('AMSUA', zlat, zlon, MTINTRP, &
                                                  MGINTRP, GLINTRP)

      ! 6) Controle de qualite des TOVS. Data QC flags (IMARQ) are modified here!
      call mwbg_tovCheckAmsua(ISAT, ilq, IORBIT, ican, ICANOMP, ztb, biasCorr, &
                              ZOMP, ICHECK, nvalOut, ntOut, ZMISG, INOSAT, ICHKPRF, &
                              scanpos, MGINTRP, MTINTRP, GLINTRP, itt, zenith, &
                              IMARQ, clw, clw_avg, scatw, STNID, RESETQC, ZLAT)

      write(*,*) ' Statistiques'
      ! Accumuler Les statistiques sur les rejets
      call mwbg_qcStats('AMSUA', satNumber, ICHECK, ican, INOSAT, nvalOut, &
                             ntOut, satelliteId, .FALSE.)
      ! set terrain type to sea ice given certain conditions 
      call mwbg_setTerrainTypeToSeaIce(GLINTRP, ilq, itt)
    
    end if 

    write(*,*) ' Mise a jour du fichier burp'
    ! 7) Mise a jour de rapport.
    call mwbg_updateBurpAmsua(burpFileNameIn, ReportIndex, ETIKRESU, clw_avg, scatw, &
                              globMarq, ICHKPRF, ilq, itt, & 
                              RESETQC, IMARQ, burpFileNameout)

    if ( ifLastReport) exit REPORTS

  end do REPORTS

  write(*,*) ' Number of obs pts read from BURP file              = ', nobs_tot
  write(*,*) ' Number of BURP file reports                        = ', reportIndex
  write(*,*) ' Number of bad BURP file reports (all data flagged) = ', n_bad_reps

  ! 10) Fin
  ! Imprimer les statistiques sur les rejets
  call mwbg_qcStats('AMSUA', satNumber, ICHECK, ican, INOSAT, nvalOut, &
                         ntOut, satelliteId, .TRUE.)
  
  junk  = exfin('bgckMW','FIN','NON')

end program midas_bgckmw

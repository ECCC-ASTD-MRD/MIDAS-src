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

program midas_inovQCAtms
  !
  ! :Purpose: Main program for background check of microwave instruments. 
  !

  !OBJET          Effectuer le controle de qualite des radiances level 1b 
  !               ATMS de NPP.
  !
  !CLES           IXENT   - FICHIER BURP CONTENANT LES DONNEES
  !                         TOVS REGROUPES                               (ENTREE)
  !               IRGEO   - FICHIER STANDARD CONTENANT LES CHAMPS 
  !                         GEOPHYSIQUES MF/MX                           (ENTREE)
  !               ISSTAT  - FICHIER SEQUENTIEL CONTENANT LES STATISTIQUES 
  !                         D'ERREUR TOTALE DES TOVS (NEW 2013 FORMAT)   (ENTREE)
  !               OXSRT   - FICHIER BURP CONTENANT LES DONNEES
  !                         TOVS REGROUPES AVEC MARQUEURS QC RAJOUTES    (SORTIE)
  !               DEBUG   - ACTIONNER LE MODE DEBUG
  !NOTES
  !
  !  Must be run after satqc_atms, satbcor and 3D-Var (O-P mode)
  !
  !        5 tests are done:                                                      QC flag bits set
  !          1) check for data rejected by SATQC_ATMS (QC flag bit 7 ON)          --> bit 9(,7)
  !          2) topography rejection for low-peaking channels (with MF/MX field), --> bits 9,18
  !          3) check for uncorrected radiance (QC flag bit 6 OFF),               --> bit 11           
  !          4) Innovation (O-P) based QC                                         --> bit 9,16
  !         5) channel blacklisting (from UTIL column in stats_atms_assim file)  --> bit 8
  !
  !       *** Array ITEST(5) in SUBROUTINE mwbg_tovCheckAtms used to select tests. ***

  !  (i)   Basic QC tests plus filtering for surface-sensitive channels, cloud water (CLW), 
  !         scattering index, Dryness Index are done in first QC program SATQC_ATMS.
  !         QC flag bit 7 is set for the rejected data. Test 1 of this program checks
  !         for such data and sets bit 9 ON (for 3D-Var rejection).
  !
  ! (ii)   This program sets data QC flag bit 9 ON in all "data reject" cases except
  !         for test 12 (blacklisting) where bit 8 is set ON.
  !
  !  To compile on Linux (pgi9xx) (on arxt10):
  !
  ! >  s.compile -src atms_inovqc.f -o atms_inovqc_Linux -librmn rmn_014_rc2 {-debug}
  !
  !  To RUN :
  ! 
  ! >  atms_inovqc_Linux -IXENT burpin -IRGEO MT_fst -ISSTAT stats_atms_assim -OXSRT burpout { -DEBUG oui }
  !
  !       burpin            =  BURP file containing ATMS Tb (Sat_QCd, bias corrected, and O-P) [bits 7 and 6 set]
  !       MT_fst            =  Standard file containing filtered model topo fields MF or MX (for TOPO check)
  !                            NOTE: GEM analysis (_000) and 3h trial (_180m) files contain these fields.
  !       stats_atms_assim  =  ATMS observation error file (new 2013 format)
  !                            NOTE: ERBGCK used for rogue O-P check, UTIL column for channel selection
  !       burpout           =  BURP file containing ATMS Tb [with bits 8, 9 and 11 set for 3dvar reject]
  !
  use bgckmicrowave_mod

  IMPLICIT NONE

  INTEGER MXELM
  INTEGER MXLAT, MXLON
  PARAMETER ( MXELM  =    30 )
  PARAMETER ( MXLAT  =     5 )
  PARAMETER ( MXLON  =     5 )

  integer ezsetopt, ezqkdef
  integer gdllsval, gdmg

  INTEGER FNOM,MRFOPN,MRFCLS,MRFPUT,MRBUPD
  INTEGER FCLOS,MRFOPR,FSTOUV
  INTEGER FSTINF,FSTPRM,FSTLIR,FSTFRM
  INTEGER MRFOPC,MRFMXL
  INTEGER HANDLE,ISTAT,NOMBRE,NVAL,NT
  INTEGER EXDB,EXFIN 
  INTEGER IER,IREC,IREC2,JUNK
  INTEGER I,ILNMX, NELE, J, JN, JL
  INTEGER IUNENT, IUNSRT, IUNGEO, IUNSTAT, INUMSAT, nulnam 
  INTEGER INO,INOMP,INOMPNA,INOSAT
  INTEGER IDUM,IDUM1,IDUM2,IDUM3,IDUM4,IDUM5,IDUM6,IDUM7
  INTEGER IDUM8,IDUM9,IDUM10,IDUM11,IDUM12,IDUM13
  INTEGER IDUM14,IDUM15,IDUM16,IDUM17,IDUM18
  INTEGER IG1,IG2,IG3,IG4
  INTEGER IG1R,IG2R,IG3R,IG4R
  INTEGER NI,NJ,NK,INDX,NLAT,NLON

  INTEGER TBLVAL    (MXELM*MXVAL*MXNT)
  INTEGER KTBLVALN  (MXELM*MXVAL*MXNT)
  INTEGER LSTELE    (MXELM)
  INTEGER KLISTEN   (MXELM)
  INTEGER ELDALT    (MXELM)
  INTEGER IDATA     (MXVAL*MXNT)
  INTEGER ITERMER   (MXNT)
  INTEGER ITERRAIN  (MXNT)
  INTEGER ISCNCNT   (MXNT)
  INTEGER ISCNPOS   (MXNT)
  INTEGER ISAT      (MXNT)
  INTEGER IORBIT    (MXNT)
  INTEGER IDENTF    (MXNT)
  INTEGER ICANO     (MXVAL*MXNT)
  INTEGER ICANOMP   (MXVAL*MXNT)
  INTEGER ICANOMPNA (MXVAL*MXNT)
  INTEGER ICHECK    (MXVAL*MXNT)
  INTEGER ICHKPRF   (MXNT)
  INTEGER IMARQ     (MXVAL*MXNT)

  CHARACTER *9   STNID
  CHARACTER *9   CSATID(MXSAT)

  CHARACTER *12  ETIKXX
  CHARACTER *9   ETIKRESU
  CHARACTER *4   NOMVXX,CLNOMVAR
  CHARACTER *2   TYPXX 
  CHARACTER *1   GRTYP

  INTEGER, ALLOCATABLE, DIMENSION(:) :: BUF1
  REAL, ALLOCATABLE, DIMENSION(:) :: MT

  REAL  DONIALT (MXELM*MXVAL*MXNT)
  REAL  PRVALN  (MXELM*MXVAL*MXNT)
  REAL  ZDATA   (MXVAL*MXNT)
  REAL  SATZEN  (MXNT)
  REAL  MTINTRP (MXNT)
  REAL  ZO      (MXVAL*MXNT)
  REAL  ZCOR    (MXVAL*MXNT)
  REAL  ZOMP    (MXVAL*MXNT)
  REAL  ZOMPNA  (MXVAL*MXNT)
  REAL  ZLATBOX (MXLAT*MXLON,MXNT)
  REAL  ZLONBOX (MXLAT*MXLON,MXNT)
  REAL  MTINTBOX(MXLAT*MXLON,MXNT)
  REAL  XLAT,XLON

  REAL  DLAT, DLON, TOPOFACT

  LOGICAL RESETQC, SKIPENR

  DATA IUNENT  / 10 /
  DATA IUNSRT  / 20 /
  DATA IUNGEO  / 50 /
  DATA IUNSTAT / 60 /
  DATA DLAT   / 0.4 /
  DATA DLON   / 0.6 /

  EXTERNAL CCARD,EXDB,EXFIN
  EXTERNAL FCLOS,FNOM

  LOGICAL DEBUG
  COMMON /DBGCOM/ DEBUG

  namelist /nambgck/ debug, RESETQC, ETIKRESU 

  ! 1) Debut
  IER = FNOM(IUNGEO,'./fstgzmx','STD+RND+R/O',0)

  JUNK = EXDB('ATMS_INOVQC','30NOV13','NON')

  DEBUG = .FALSE.
  RESETQC = .FALSE.
  ETIKRESU = '>>BGCKALT'

  ! read namelist
  nulnam = 0
  ier = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read(nulnam, nml=nambgck, iostat=ier)
  if ( ier /= 0 ) then
    write(*,*) 'midas_inovQCAtms: Error reading namelist'
    call abort()
  end if
  write(*,nml=nambgck)
  ier = fclos(nulnam)

  ! ouverture du fichier entree burp en mode lecture
  IER = FNOM(IUNENT,'obsatms','RND',0)
  IF(IER .NE. 0) THEN
     PRINT *,' ERREUR D ASSOCIATION DE FICHIER'
     STOP
  ENDIF
  NOMBRE =	MRFOPN(IUNENT,'READ')
  ISTAT =	MRFOPC('MSGLVL','ERROR')

  ! ouverture du fichier sortie burp 
  IER = FNOM(IUNSRT,'obsatms.out','RND',0)
  IF(IER .NE. 0) THEN
     PRINT *,' ERREUR D ASSOCIATION DE FICHIER'
     STOP
  ENDIF
  NOMBRE =	MRFOPN(IUNSRT,'CREATE')

  ! Lecture du fichier burp entree
  ILNMX = MRFMXL(IUNENT)
  IF (DEBUG) THEN
     write(*,*)'MRFMXL: ILNMX =',ILNMX
  ENDIF
  ALLOCATE ( buf1(ILNMX*2), STAT=ier)
  BUF1(1) = ILNMX*2

  ! Valeur manquante burp
  ISTAT = MRFOPR('MISSING',ZMISG)
  IF (DEBUG) THEN
    write(*,*)' MISSING VALUE =', ZMISG
  ENDIF

  ! 2) Lecture des statistiques d'erreur totale pour les  TOVS 
  IER = FNOM(IUNSTAT,'./stats_atms_assim','SEQ+FMT',0)
  IF(IER.LT.0)THEN
    write(*,*) '(" ATMS_INOVQC: Problem opening ", &
          "ATMS total error statistics file ", stats_atms_assim)'
    CALL ABORT ()
  END IF
  CALL mwbg_readStatTovsAtms(IUNSTAT,INUMSAT,CSATID)
  write(*,*) " SATID's = "
  DO I = 1, INUMSAT
    write(*,*) '  ', CSATID(I)
  ENDDO

  ! 3) Lecture des champs geophysiques (MF/MX) du modele
  IER = FSTOUV(IUNGEO,'RND')

  ! TOPOGRAPHIE (MF ou MX).
  !     MF est la topographie filtree avec unites en metres (filtered ME).
  !     MX est la topographie filtree avec unites en m2/s2  (geopotential topography).
  TOPOFACT = 1.0
  IREC = FSTINF(IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1,' ','MF')
  CLNOMVAR = 'MF'
  IF (IREC .LT. 0) THEN
    TOPOFACT = 9.80616
    IREC = FSTINF(IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1,' ','MX')
    CLNOMVAR = 'MX'
  ENDIF
  IF (IREC .LT. 0) THEN
    write(*,*) ' LA TOPOGRAPHIE (MF or MX) EST INEXISTANTE' 
    CALL ABORT()
  ELSE
    ALLOCATE ( MT(NI*NJ), STAT=ier)
    IER = FSTLIR(MT,IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1, &
         ' ',CLNOMVAR)
  ENDIF
      
  IER = FSTPRM ( IREC, IDUM1, IDUM2, IDUM3, IDUM4, & 
      IDUM5, IDUM6, IDUM7, IDUM8, IDUM9, IDUM10,  &
      IDUM11, TYPXX, NOMVXX, ETIKXX, GRTYP, IG1, &
      IG2, IG3, IG4, IDUM12, IDUM13, IDUM14,  &
      IDUM15, IDUM16, IDUM17, IDUM18 )
  
  write(*,*) ' GRTYP = ', grtyp 

  ! 4) Lire les donnees TOVS
  HANDLE = 0
1000 CALL mwbg_readTovsAtms(IUNENT,HANDLE,ISAT,ZMISG,BUF1,TBLVAL, &
           LSTELE,ELDALT,IDATA,ICANO,ICANOMP,ICANOMPNA,INO, &
           INOMP,INOMPNA,DONIALT,ZDATA,ZO,ZCOR,ZOMP,STNID, &
           ZOMPNA,IMARQ,MXELM,MXVAL,MXNT,NELE,NVAL,NT, &
           ISCNCNT,ISCNPOS,MXSCAN,ZLAT,ZLON,ITERMER, &
           IORBIT,SATZEN,ITERRAIN,SKIPENR,IDENTF)

  ! All data read?
  IF ( HANDLE .LE. 0 ) GO TO 2500

  ! Enregistrement resume?
  IF ( STNID(1:2) .EQ. '>>' ) THEN   
    IER=MRBUPD(IUNSRT,BUF1,-1,-1,ETIKRESU,-1,-1,-1,-1,-1, &
            -1,-1,-1,-1,-1,-1,0,-1,0)                                            
    ISTAT = MRFPUT(IUNSRT,0,BUF1)
  ENDIF  

  ! Sauter l'enregistrement?
  IF ( SKIPENR ) GO TO 1000
     
  ! trouver l'indice du satellite
  INOSAT = 0
  DO I = 1,MXSAT
    IF ( STNID .EQ. '^'//CSATID(I) ) THEN
      INOSAT = I
    ENDIF
  ENDDO
  IF ( INOSAT .EQ. 0 ) THEN
    write(*,*)'SATELLITE NON-VALIDE', STNID
    CALL ABORT()
  ENDIF

  ! 5) Interpolation de le champ MF/MX (topogrpahy) aux pts TOVS.
  !    N.B.: on examine ce champ sur une boite centree sur chaque obs.
  NLAT = (MXLAT-1)/2
  NLON = (MXLON-1)/2
  DO JN = 1, NT
    INDX = 0
    DO I = -NLAT, NLAT
      XLAT = ZLAT(JN) +I*DLAT
      XLAT = MAX(-90.0,MIN(90.0,XLAT))
      DO J = -NLON, NLON
        INDX = INDX + 1
        XLON = ZLON(JN) +J*DLON
        IF ( XLON .LT. -180. ) XLON = XLON + 360.
        IF ( XLON .GT.  180. ) XLON = XLON - 360.
        IF ( XLON .lt.    0. ) XLON = XLON + 360.
        ZLATBOX(INDX,JN) = XLAT
        ZLONBOX(INDX,JN) = XLON
      ENDDO
    ENDDO
  ENDDO

  ier  = ezsetopt('INTERP_DEGREE','LINEAR')
  gdmg = ezqkdef(ni,nj,grtyp,ig1,ig2,ig3,ig4,iungeo)
  ier=gdllsval (gdmg,mtintbox,mt,ZLATBOX,ZLONBOX,MXLAT*MXLON*NT)

  DO JN = 1, NT
    IF (DEBUG) THEN
      PRINT *, ' ------------------  '
      PRINT *, ' JN = ', JN
      PRINT *, '   '
      PRINT *, ' zlat,zlon = ', zlat(jn), zlon(jn)
      PRINT *, '   '
      PRINT *, ' ZLATBOX = '
      PRINT *,  (ZLATBOX(I,JN),I=1,MXLAT*MXLON)
      PRINT *, ' ZLONBOX = '
      PRINT *,  (ZLONBOX(I,JN),I=1,MXLAT*MXLON)
      PRINT *, ' MTINTBOX = '
      PRINT *,  (MTINTBOX(I,JN),I=1,MXLAT*MXLON)
    ENDIF
    MTINTRP(JN) = 0.0
    DO I=1,MXLAT*MXLON
      MTINTRP(JN) = MAX(MTINTRP(JN),MTINTBOX(I,JN)/TOPOFACT)
    ENDDO
    IF (DEBUG) THEN
      PRINT *, ' MTINTRP = ', MTINTRP(JN)
    ENDIF
  ENDDO

  ! 6) Controle de qualite des TOVS. Data QC flags (IMARQ) are modified here!
  CALL mwbg_tovCheckAtms(ISAT,IORBIT,ICANO,ICANOMP,ZO,ZCOR, &
               ZOMP,ICHECK,INO,INOMP,NT,ZMISG,INOSAT,IDENTF, &
               ICHKPRF,ISCNPOS,MTINTRP,IMARQ,STNID,RESETQC)

  ! Accumuler Les statistiques sur les rejets
  CALL mwbg_qcStatsAtms (INUMSAT,ICHECK,ICANO, &
                INOSAT,CSATID,INO, NT, .FALSE.)

  ! 7) Mise a jour des marqueurs.
  CALL mwbg_updatFlgAtms (BUF1,LSTELE,TBLVAL,ELDALT,DONIALT, &
                 IDATA,ZDATA,ICHKPRF,ICHECK, & 
                 RESETQC,IMARQ)

  ! Ecriture de l'enregistrement sur le fichier burp de sortie
  ISTAT = MRFPUT(IUNSRT,0,BUF1) 

  GO TO 1000  
     
  ! All data read!

2500  continue

  ! 9) Fin
  ! Imprimer les statistiques sur les rejets
  CALL mwbg_qcStatsAtms (INUMSAT,ICHECK,ICANO, &
               INOSAT,CSATID,INO, NT, .TRUE.)

  ! fermeture des fichiers 
9999  CONTINUE
  JUNK  = EXFIN('ATMS_INOVQC',' ','NON')
  ISTAT = MRFCLS(IUNENT)
  ISTAT = MRFCLS(IUNSRT)
  ISTAT = FSTFRM(IUNGEO)
  ISTAT = FCLOS (IUNGEO)
  ISTAT = FCLOS (IUNSTAT)

  STOP

end program midas_inovQCAtms

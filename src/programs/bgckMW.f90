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
  use bgckmicrowave_mod

  implicit none

  INTEGER MXELM, MXVAL, MXNT, NCLES, MXSAT, MXSCAN
  INTEGER MXLAT, MXLON
  PARAMETER ( NCLES  =     9 )
  PARAMETER ( MXELM  =    40 )
  PARAMETER ( MXVAL  =    50 )
  PARAMETER ( MXNT   =  2000 )
  PARAMETER ( MXSAT  =     9 )
  PARAMETER ( MXSCAN =    56 )
  PARAMETER ( MXLAT  =     5 )
  PARAMETER ( MXLON  =     5 )

  integer ezsetopt, ezsetval, ezqkdef
  integer gdllsval, gdxyfll, gdmg, gdmt, gdgl

  INTEGER FNOM,MRFOPN,MRFCLS,MRFPUT,MRBUPD
  INTEGER FCLOS,MRFOPR,FSTOUV
  INTEGER FSTINF,FSTPRM,FSTLIR,FSTFRM
  INTEGER MRFOPC,MRFMXL
  INTEGER HANDLE,ISTAT,NOMBRE,NVAL,NT
  INTEGER EXDB,EXFIN 
  INTEGER NPOSIT, IER,IREC,IREC2,JUNK
  INTEGER I,ILNMX, NELE, J, JN, JL
  INTEGER IUNENT, IUNSRT, IUNGEO, IUNSTAT, INUMSAT, nulnam 
  INTEGER INO,INOMP,INOMPNA,INOSAT
  INTEGER IDUM,IDUM1,IDUM2,IDUM3,IDUM4,IDUM5,IDUM6,IDUM7
  INTEGER IDUM8,IDUM9,IDUM10,IDUM11,IDUM12,IDUM13
  INTEGER IDUM14,IDUM15,IDUM16,IDUM17,IDUM18
  INTEGER IG1GL,IG2GL,IG3GL,IG4GL
  INTEGER IG1MG,IG2MG,IG3MG,IG4MG
  INTEGER IG1MT,IG2MT,IG3MT,IG4MT
  INTEGER NITIC,NJTAC,NIGL,NJGL,NK,INDX,NLAT,NLON
  INTEGER NIMG,NJMG,NIMT,NJMT

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
  INTEGER ICANO     (MXVAL*MXNT)
  INTEGER ICANOMP   (MXVAL*MXNT)
  INTEGER ICANOMPNA (MXVAL*MXNT)
  INTEGER ICHECK    (MXVAL*MXNT)
  INTEGER ICHKPRF   (MXNT)
  INTEGER IMARQ     (MXVAL*MXNT)

  CHARACTER(len=9)   STNID
  CHARACTER(len=9)   CSATID(MXSAT)

  CHARACTER(len=12)  ETIKXX
  CHARACTER(len=9)   ETIKRESU
  CHARACTER(len=4)   NOMVXX,CLNOMVAR
  CHARACTER(len=2)   TYPXX 
  CHARACTER(len=1)   GRTYPGL,GRTYPMG,GRTYPMT

  INTEGER, ALLOCATABLE, DIMENSION(:) :: BUF1
  REAL, ALLOCATABLE, DIMENSION(:) :: MG
  REAL, ALLOCATABLE, DIMENSION(:) :: MT
  REAL, ALLOCATABLE, DIMENSION(:) :: GL
  REAL, ALLOCATABLE, DIMENSION(:) :: TIC
  REAL, ALLOCATABLE, DIMENSION(:) :: TAC

  REAL  DONIALT (MXELM*MXVAL*MXNT)
  REAL  PRVALN  (MXELM*MXVAL*MXNT)
  REAL  ZDATA   (MXVAL*MXNT)
  REAL  ZLAT    (MXNT)
  REAL  ZLON    (MXNT)
  REAL  SATZEN  (MXNT)
  REAL  MGINTRP (MXNT)
  REAL  MTINTRP (MXNT)
  REAL  GLINTRP (MXNT)
  REAL  ZO      (MXVAL*MXNT)
  REAL  ZCOR    (MXVAL*MXNT)
  REAL  ZOMP    (MXVAL*MXNT)
  REAL  ZOMPNA  (MXVAL*MXNT)
  REAL  ZLATBOX (MXLAT*MXLON,MXNT)
  REAL  ZLONBOX (MXLAT*MXLON,MXNT)
  REAL  MGINTBOX(MXLAT*MXLON,MXNT)
  REAL  MTINTBOX(MXLAT*MXLON,MXNT)
  REAL  GLINTBOX(MXLAT*MXLON,MXNT)
  REAL  XLAT,XLON

  REAL  ZMISG, DLAT, DLON, TOPOFACT

  LOGICAL RESETQC, SKIPENR

  DATA IUNENT  / 10 /
  DATA IUNSRT  / 20 /
  DATA IUNGEO  / 50 /
  DATA IUNSTAT / 60 /
  DATA ZMISG  /9.9E09 /
  DATA DLAT   / 0.4 /
  DATA DLON   / 0.6 /

  EXTERNAL EXDB,EXFIN
  EXTERNAL FCLOS,FNOM

  LOGICAL DEBUG
  COMMON /DBGCOM/ DEBUG

  namelist /nambgck/ debug, RESETQC, ETIKRESU 

  ! 1) Debut
      IER = FNOM(6     ,'./out.txt' ,'SEQ'    ,0) 
      IER = FNOM(IUNGEO,'./fstglmg' ,'STD+RND+R/O',0)

      JUNK = EXDB('SATQC_AMSUA','07MAR14','NON')

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

  ! ouverture du fichier entree burp en mode lecture
  IER = FNOM(IUNENT,'./obsto_amsua','RND',0)
  IF(IER .NE. 0) THEN
     PRINT *,' ERREUR D ASSOCIATION DE FICHIER'
     STOP
  ENDIF
  NOMBRE =	MRFOPN(IUNENT,'READ')
  ISTAT =	MRFOPC('MSGLVL','ERROR')

  ! ouverture du fichier sortie burp 
  IER = FNOM(IUNSRT,'./obsto_amsua.out','RND',0)
  IF(IER .NE. 0) THEN
     PRINT *,' ERREUR D ASSOCIATION DE FICHIER'
     STOP
  ENDIF
  NOMBRE =	MRFOPN(IUNSRT,'CREATE')

  ! Lecture du fichier burp entree
  ILNMX = MRFMXL(IUNENT)
  IF (DEBUG) THEN
     WRITE(6,*)'MRFMXL: ILNMX =',ILNMX
  ENDIF
  ALLOCATE ( buf1(ILNMX*2), STAT=ier)
  BUF1(1) = ILNMX*2

  ! Valeur manquante burp
  ISTAT = MRFOPR('MISSING',ZMISG)
  IF (DEBUG) THEN
     WRITE(6,*)' MISSING VALUE =', ZMISG
  ENDIF

  ! 2) Lecture des statistiques d'erreur totale pour les  TOVS 
  IER = FNOM(IUNSTAT,'./stats_amsua_assim','SEQ+FMT',0)
  IF(IER.LT.0)THEN
     WRITE (6,*) '(" SATQC_AMSUA: Problem opening ", &
           "TOVS total error statistics file ", stats_amsua_assim)'               
     CALL ABORT ()
  END IF
  CALL mwbg_readStatTovs(IUNSTAT,INUMSAT,CSATID)
  WRITE(6,*) " SATID's = "
  DO I = 1, INUMSAT
     WRITE(6,*) '  ', CSATID(I)
  ENDDO

  ! 3) Lecture des champs geophysiques (GL,MG,MT) du modele
  IER = FSTOUV(IUNGEO,'RND')
  ! MG
  IREC = FSTINF(IUNGEO,NIMG,NJMG,NK,-1,' ',-1,-1,-1,' ','MG')
  IF (IREC .LT. 0) THEN
     WRITE (6,*) ' LE MASQUE TERRE-MER EST INEXISTANT' 
     CALL ABORT()
  ENDIF

  ALLOCATE ( mg(NIMG*NJMG), STAT=ier)
  IER = FSTLIR(MG,IUNGEO,NIMG,NJMG,NK,-1,' ',-1,-1,-1,&
               ' ','MG')

  IER = FSTPRM ( IREC, IDUM1, IDUM2, IDUM3, IDUM4, &
       IDUM5, IDUM6, IDUM7, IDUM8, IDUM9, IDUM10, &
       IDUM11, TYPXX, NOMVXX, ETIKXX, GRTYPMG, IG1MG,&
       IG2MG, IG3MG, IG4MG, IDUM12, IDUM13, IDUM14, &
       IDUM15, IDUM16, IDUM17, IDUM18 )

  ! TOPOGRAPHIE (ME ou MX).
  ! ME est la topographie avec unites en metres.
  ! MX est la topographie avec unites en m2/s2.
  TOPOFACT = 1.0
  IREC = FSTINF(IUNGEO,NIMT,NJMT,NK,-1,' ',-1,-1,-1,' ','ME')
  CLNOMVAR = 'ME'
  IF (IREC .LT. 0) THEN
     TOPOFACT = 9.80616
     IREC = FSTINF(IUNGEO,NIMT,NJMT,NK,-1,' ',-1,-1,-1,' ','MX')
     CLNOMVAR = 'MX'
  ENDIF
  IF (IREC .LT. 0) THEN
     WRITE (6,*) ' LA TOPOGRAPHIE EST INEXISTANTE' 
     CALL ABORT()
  ENDIF

  ALLOCATE ( MT(NIMT*NJMT), STAT=ier)
  IER = FSTLIR(MT,IUNGEO,NIMT,NJMT,NK,-1,' ',-1,-1,-1,&
               ' ',CLNOMVAR)

  IER = FSTPRM ( IREC, IDUM1, IDUM2, IDUM3, IDUM4, & 
       IDUM5, IDUM6, IDUM7, IDUM8, IDUM9, IDUM10, & 
       IDUM11, TYPXX, NOMVXX, ETIKXX, GRTYPMT, IG1MT, &
       IG2MT, IG3MT, IG4MT, IDUM12, IDUM13, IDUM14, &
       IDUM15, IDUM16, IDUM17, IDUM18 )

  ! GL
  IREC = FSTINF(IUNGEO,NIGL,NJGL,NK,-1,' ',-1,-1,-1,' ','GL')
  IF (IREC .LT. 0) THEN
     WRITE (6,*) 'LE CHAMP GL EST INEXISTANT' 
     CALL ABORT()
  ENDIF

  ALLOCATE ( gl(NIGL*NJGL), STAT=ier)
  IER = FSTLIR(GL,IUNGEO,NIGL,NJGL,NK,-1,' ',-1,-1,-1, &
               ' ','GL')

  IER = FSTPRM ( IREC, IDUM1, IDUM2, IDUM3, IDUM4, & 
       IDUM5, IDUM6, IDUM7, IDUM8, IDUM9, IDUM10, &
       IDUM11, TYPXX, NOMVXX, ETIKXX, GRTYPGL, IG1GL, &
       IG2GL, IG3GL, IG4GL, IDUM12, IDUM13, IDUM14, &
       IDUM15, IDUM16, IDUM17, IDUM18 )

  ! On regle les parametres pour l'interpolation avec EZSCINT
  ier = ezsetopt('INTERP_DEGREE','LINEAR')
  WRITE (*,*) ' GRILLE MG : ',grtypmg,nimg,njmg, &
              ig1mg,ig2mg,ig3mg,ig4mg
  gdmg= ezqkdef(nimg,njmg,grtypmg,ig1mg,ig2mg,ig3mg,ig4mg,iungeo)
  WRITE (*,*) ' GRILLE MT : ',grtypmt,nimt,njmt, &
              ig1mt,ig2mt,ig3mt,ig4mt
  gdmt= ezqkdef(nimt,njmt,grtypmt,ig1mt,ig2mt,ig3mt,ig4mt,iungeo)
  WRITE (*,*) ' GRILLE GL : ',grtypgl,nigl,njgl, &
              ig1gl,ig2gl,ig3gl,ig4gl
  gdgl= ezqkdef(nigl,njgl,grtypgl,ig1gl,ig2gl,ig3gl,ig4gl,iungeo)

  ! 4) Lire les donnees TOVS
  HANDLE = 0
1000 CALL mwbg_readTovs(IUNENT,HANDLE,ISAT,ZMISG,BUF1,TBLVAL, &
           LSTELE,ELDALT,IDATA,ICANO,ICANOMP,ICANOMPNA,INO, &
           INOMP,INOMPNA,DONIALT,ZDATA,ZO,ZCOR,ZOMP,STNID, &
           ZOMPNA,IMARQ,MXELM,MXVAL,MXNT,NELE,NVAL,NT, &
           ISCNCNT,ISCNPOS,MXSCAN,ZLAT,ZLON,ITERMER, &
           IORBIT,SATZEN,ITERRAIN,SKIPENR)

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
    WRITE(6,*) 'SATELLITE ',TRIM(STNID), &
               ' NOT FOUND IN STATS FILE!'
    CALL ABORT()
  ENDIF

  ! 5) Interpolation de la glace et le champ terre/mer du modele aux pts TOVS.
  ! N.B.: on examine ces champs sur une boite centree sur chaque obs.
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

  ier = gdllsval(gdmg,mgintbox,mg,ZLATBOX,ZLONBOX,MXLAT*MXLON*NT)
  IF (ier .lt. 0) THEN
     WRITE(*,*) 'ERROR in the interpolation of MG'
     CALL ABORT()
  ENDIF
  ier = gdllsval(gdmt,mtintbox,mt,ZLATBOX,ZLONBOX,MXLAT*MXLON*NT)
  IF (ier .lt. 0) THEN
     WRITE(*,*) 'ERROR in the interpolation of MT'
     CALL ABORT()
  ENDIF
  ier = gdllsval(gdgl,glintbox,gl,ZLATBOX,ZLONBOX,MXLAT*MXLON*NT)
  IF (ier .lt. 0) THEN
     WRITE(*,*) 'ERROR in the interpolation of GL'
     CALL ABORT()
  ENDIF

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
        PRINT *, ' MGINTBOX = '
        PRINT *,  (MGINTBOX(I,JN),I=1,MXLAT*MXLON)
        PRINT *, ' MTINTBOX = '
        PRINT *,  (MTINTBOX(I,JN),I=1,MXLAT*MXLON)
        PRINT *, ' GLINTBOX = '
        PRINT *,  (GLINTBOX(I,JN),I=1,MXLAT*MXLON)
     ENDIF
     MGINTRP(JN) = 0.0
     MTINTRP(JN) = 0.0
     GLINTRP(JN) = 0.0
     DO I=1,MXLAT*MXLON
        MGINTRP(JN) = MAX(MGINTRP(JN),MGINTBOX(I,JN))
        MTINTRP(JN) = MAX(MTINTRP(JN),MTINTBOX(I,JN)/TOPOFACT)
        GLINTRP(JN) = MAX(GLINTRP(JN),GLINTBOX(I,JN))
     ENDDO
     IF (DEBUG) THEN
        PRINT *, ' MGINTRP = ', MGINTRP(JN)
        PRINT *, ' MTINTRP = ', MTINTRP(JN)
        PRINT *, ' GLINTRP = ', GLINTRP(JN)
     ENDIF
  ENDDO

  ! 6) Controle de qualite des TOVS. Data QC flags (IMARQ) are modified here!
  CALL mwbg_tovCheckAmsua(ISAT,ITERMER,IORBIT,ICANO, &
                 ICANOMP,ZO,ZCOR,ZOMP,ICHECK, &
                 INO,INOMP,NT,ZMISG, &
                 INOSAT,ICHKPRF,ISCNPOS, &
                 MGINTRP,MTINTRP,GLINTRP,ITERRAIN,SATZEN, &
                 IMARQ,STNID, &
                 RESETQC,ZLAT)
  ! Accumuler Les statistiques sur les rejets
  CALL mwbg_qcStatsAmsua(INUMSAT,ICHECK,ICANO, &
                 INOSAT,CSATID,INO, NT,.FALSE.)

  ! 7) Mise a jour des marqueurs.
  CALL mwbg_UPDATFLG(BUF1,LSTELE,TBLVAL,ELDALT,DONIALT, &
                  IDATA,ZDATA,ICHKPRF,ITERMER,ICHECK, &
                  RESETQC,IMARQ)

  ! 8) Remplacer le terrain type en fichier burp.
  CALL mwbg_ADDTRRN(BUF1,LSTELE,TBLVAL,ELDALT,DONIALT, &
                IDATA,ZDATA,ITERMER,ITERRAIN, &
                GLINTRP,KTBLVALN,KLISTEN,PRVALN)

  ! 9) Option d'enlever les blocs SATEM.
  ! not done any more, j.h. june 2001
  ! Ecriture de l'enregistrement sur le fichier burp de sortie
  ISTAT = MRFPUT(IUNSRT,0,BUF1) 

  GO TO 1000  
     
  ! All data read!
2500  continue

  ! 10) Fin
  ! Imprimer les statistiques sur les rejets
  CALL mwbg_qcStatsAmsua(INUMSAT,ICHECK,ICANO, &
               INOSAT,CSATID,INO, NT,.TRUE.)

  ! fermeture des fichiers 
9999  CONTINUE
  JUNK  = EXFIN('SATQC_AMSUA',' ','NON')
  ISTAT = MRFCLS(IUNENT)
  ISTAT = MRFCLS(IUNSRT)
  ISTAT = FSTFRM(IUNGEO)
  ISTAT = FCLOS (IUNGEO)
  ISTAT = FCLOS (IUNSTAT)

  STOP

end program midas_bgckmw

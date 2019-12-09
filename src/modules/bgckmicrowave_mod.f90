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

module bgckmicrowave_mod
  ! MODULE bgckmicrowave_mod (prefix='mwbg' category='1. High-level functionality')
  !
  ! :Purpose: Variables for microwave background check and quality control.
  !
  !use 
  implicit none
  save
  private
  ! Public functions (methods)
  public :: mwbg_readStatTovs, mwbg_readTovs, mwbg_tovCheckAmsua, mwbg_qcStatsAmsua, mwbg_UPDATFLG, mwbg_ADDTRRN  

  integer, PARAMETER :: MXVAL = 50
  integer, PARAMETER :: MXNT = 2000

contains


  INTEGER FUNCTION ISRCHEQR (KLIST, KLEN, KENTRY)
! OBJET          Rechercher un element dans une liste (valeurs reelles).
! APPEL          INDX = ISRCHEQR (KLIST, KLEN, KENTRY)
!ARGUMENTS      indx    - output -  position de l'element recherche:
!                                   =0, element introuvable,
!                                   >0, position de l'element trouve,
!               klist   - input  -  la liste
!               klen    - input  -  longueur de la liste
!               kentry  - input  -  l'element recherche

  IMPLICIT NONE

  INTEGER  KLEN, JI

  REAL  KLIST(KLEN)
  REAL  KENTRY

  ISRCHEQR = 0
  DO JI=1,KLEN
     IF ( NINT(KLIST(JI)) .EQ. NINT(KENTRY) ) THEN
        ISRCHEQR = JI
        RETURN
     ENDIF
  ENDDO

  RETURN
  END


  FUNCTION ISRCHEQI (KLIST, KLEN, KENTRY) result(ISRCHEQI_out)

!OBJET          Rechercher un element dans une liste (valeurs entieres).
!ARGUMENTS      indx    - output -  position de l'element recherche:
!                                   =0, element introuvable,
!                                   >0, position de l'element trouve,
!               klist   - input  -  la liste
!               klen    - input  -  longueur de la liste
!               kentry  - input  -  l'element recherche

  IMPLICIT NONE

  integer :: ISRCHEQI_out
  INTEGER  KLEN, JI

  INTEGER  KLIST(KLEN)
  INTEGER  KENTRY

  ISRCHEQI_out = 0
  DO JI=1,KLEN
     IF ( KLIST(JI) .EQ. KENTRY ) THEN
        ISRCHEQI_out = JI
        RETURN
     ENDIF
  ENDDO

  RETURN
  END


  SUBROUTINE XTRBLK (KBTYP,KBFAM,KBUF,KLISTE,KTBLVAL,KDLISTE, &
                    PRVAL,KNELE,KNVAL,KNT,KBLKNO)
!OBJET          Extraire un bloc d'un fichier burp.
!*ARGUMENTS      kbtyp   - input  -  BTYP du bloc recherche
!               kbfam   - input  -  descripteur de la famille du bloc recherche
!               kbuf    - input  -  tableau contenant le rapport
!               kliste  - output -  liste des noms d'elements
!               ktblval - output -  champ de donnees (valeurs entieres)
!               kdliste - output -  liste des noms d'elements (decodes)
!               prval   - output -  champ de donnees (valeurs reelles)
!               knele   - output -  nombre d'elements
!               knval   - output -  nombre de donnees par element
!               knt     - output -  nombre de groupes KNELE X KNVAL
!               kblkno  - output -  numero d'ordre du bloc

  IMPLICIT NONE

  INTEGER KBUF    (:)
  INTEGER KLISTE  (:)
  INTEGER KDLISTE (:)
  INTEGER KTBLVAL (:)

  INTEGER KBLKNO,KNELE,KNVAL,KNT,KBTYP,ISTAT
  INTEGER MRBLOC,MRBXTR,MRBTYP,MRBCVT,MRBDCL,MRBPRM
  INTEGER IBFAM,IBDESC,IBTYP,INBIT,IBIT0,IDATYP
  INTEGER IBKNAT,IBKTYP,IBKSTP,JN,KBFAM

  REAL PRVAL(:)

  LOGICAL DEBUG
  COMMON /DBGCOM/ DEBUG

  KBLKNO = 0
  KBLKNO =  MRBLOC(KBUF,KBFAM,-1,KBTYP,KBLKNO)
  IF(KBLKNO .GT. 0) THEN
     IF ( DEBUG) THEN
        WRITE(6,*)'FOUND BLOCK WITH BTYP = ',KBTYP
     ENDIF

! extraction du bloc
     ISTAT = MRBXTR(KBUF,KBLKNO,KLISTE,KTBLVAL)

! extraction des parametres descripteurs du bloc
     ISTAT = MRBPRM (KBUF,KBLKNO,KNELE,KNVAL,KNT,IBFAM, &
                   IBDESC,IBTYP,INBIT,IBIT0,IDATYP)

! extrait bktyp et btyp
     ISTAT = MRBTYP(IBKNAT,IBKTYP,IBKSTP,IBTYP)
     IF ( DEBUG) THEN
        WRITE(6,*)'IBTYP=',IBTYP
        WRITE(6,*)'IBKNAT,IBKTYP,IBKSTP=',IBKNAT,IBKTYP,IBKSTP
        WRITE(6,*)'KNELE,KNVAL,KNT,INBIT,IBDESC,IBFAM,IDATYP=' &
                 ,KNELE,KNVAL,KNT,INBIT,IBDESC,IBFAM,IDATYP
     ENDIF

! extraction des informations du bloc et de la liste
     ISTAT = MRBXTR(KBUF,KBLKNO,KLISTE,KTBLVAL)
     ISTAT = MRBCVT(KLISTE,KTBLVAL,PRVAL,KNELE, &
                   KNVAL,KNT,0)
     ISTAT = MRBDCL(KLISTE,KDLISTE,KNELE)

     IF ( DEBUG) THEN
        WRITE(6,*) 'PRVAL:'
        WRITE(6,*) (PRVAL(JN),JN=1,KNELE*KNVAL*KNT)
        WRITE(6,*) 'KDLISTE:'
        WRITE(6,*) (KDLISTE(JN), JN=1,KNELE)
     ENDIF

  ENDIF

  RETURN
  END


  SUBROUTINE RMBLK (KBTYP,KBFAM,KBUF)
!OBJET          Enlever un bloc d'un fichier burp.
!ARGUMENTS      kbtyp   - input  -  BTYP du bloc
!               kbfam   - input  -  descripteur de la famille du bloc
!               kbuf    - in/out -  tableau contenant le rapport
  IMPLICIT NONE

  INTEGER KBUF    (:)

  INTEGER KBLKNO,KBTYP,ISTAT,KBFAM
  INTEGER MRBLOC,MRBDEL

  LOGICAL DEBUG
  COMMON /DBGCOM/ DEBUG

  KBLKNO = 0
  KBLKNO =  MRBLOC(KBUF,KBFAM,-1,KBTYP,KBLKNO)
  IF(KBLKNO .GT. 0) THEN
     IF ( DEBUG) THEN
        WRITE(6,*)'FOUND BLOCK WITH BTYP = ',KBTYP
     ENDIF

! enlever le bloc
     ISTAT = MRBDEL(KBUF,KBLKNO)
     IF(ISTAT .NE. 0) THEN
        WRITE(6,*)' PROBLEME AVEC RMBLK - BLOC ',KBLKNO,' ENLEVE'
     ENDIF
     IF ( DEBUG) THEN
        WRITE(6,*)'MRBDEL ISTAT = ',ISTAT
     ENDIF

  ENDIF

  RETURN
  END


  SUBROUTINE XTRDATA (KDLISTE,KTBLVAL,PRVAL,KNELE,KNVAL,KNT, &
                         KELEM,KDATA,PDATA,KPNTR)
!OBJET          Extraire les donnees de l'element specifie 
!               (uni ou multi niveaux).
!
!ARGUMENTS      kdliste - input  -  liste des noms d'elements (decodes)
!               ktblval - input  -  champ de donnees (valeurs entieres)
!               prval   - input  -  champ de donnees (valeurs reelles)
!               knele   - input  -  nombre d'elements
!               knval   - input  -  nombre de donnees par element
!               knt     - input  -  nombre de groupes KNELE X KNVAL
!               kelem   - input  -  element recherche 
!               kdata   - output -  donnees extraites (valeurs entieres)
!               pdata   - output -  donnees extraites (valeurs reelles)
!               kpntr   - output -  pointeur de l'element recherche:
!                                   =0, element introuvable,
!                                   >0, pointeur de l'element.
  IMPLICIT NONE

  INTEGER KDLISTE (:)
  INTEGER KTBLVAL (:)
  INTEGER KDATA(:)

  REAL PRVAL(:)
  REAL PDATA(:)

  INTEGER KNELE,KNVAL,KNT,JI,KPNTR
  INTEGER INDX,KELEM,JJ,IPOS

  LOGICAL DEBUG
  COMMON /DBGCOM/ DEBUG

  KPNTR = 0
  DO JI = 1, KNELE
     IF ( KDLISTE(JI) .EQ. KELEM ) THEN
        KPNTR = JI
        GO TO 100
     ENDIF
  ENDDO

! any missing elements?

100  CONTINUE
  IF ( KPNTR .EQ. 0  ) THEN
     IF (DEBUG) THEN
        WRITE(6,*) ' XTRDATA: No data for element ',KELEM
     ENDIF
     RETURN
  ELSE 
     IPOS = 0
     DO JI = 1, KNT 
        DO JJ = 1, KNVAL 
           IPOS = IPOS + 1
           INDX = (JI-1)*KNELE*KNVAL + (JJ-1)*KNELE + KPNTR 
           KDATA(IPOS) = KTBLVAL(INDX) 
           PDATA(IPOS) = PRVAL  (INDX)
        ENDDO
     ENDDO
     IF (DEBUG) THEN
        WRITE(6,*) ' XTRDATA: kdata =  ',(KDATA(JJ), &
                  JJ=1,KNVAL*KNT)
        WRITE(6,*) ' XTRDATA: pdata =  ',(PDATA(JJ), &
                  JJ=1,KNVAL*KNT)
     ENDIF
  ENDIF

  RETURN
  END


  SUBROUTINE REPDATA (KDLISTE,KTBLVAL,KNELE,KNVAL,KNT, &
                         KELEM,KDATA,KPNTR)
! OBJET          Remplacer les donnees de l'element specifie 
!               (uni ou multi niveaux). Valeurs entieres seulement.
!
!APPEL          CALL REPDATA (KDLISTE,KTBLVAL,KNELE,KNVAL,KNT,
!                             KELEM,KDATA,KPNTR)  
!
!ARGUMENTS      kdliste - input  -  liste des noms d'elements (decodes)
!               ktblval - in/out -  champ de donnees (valeurs entieres)
!               knele   - input  -  nombre d'elements
!               knval   - input  -  nombre de donnees par element
!               knt     - input  -  nombre de groupes KNELE X KNVAL
!               kelem   - input  -  element specifie 
!               kdata   - input  -  donnees extraites (valeurs entieres)
!               kpntr   - output -  pointeur de l'element a remplacer:
!                                   =0, element introuvable,
!                                   >0, pointeur de l'element.
    IMPLICIT NONE

    INTEGER KDLISTE (:)
    INTEGER KTBLVAL (:)
    INTEGER KDATA(:)

    INTEGER KNELE,KNVAL,KNT,JI,KPNTR
    INTEGER INDX,KELEM,JJ,IPOS

    LOGICAL DEBUG
    COMMON /DBGCOM/ DEBUG

    KPNTR = 0
    DO JI = 1, KNELE
       IF ( KDLISTE(JI) .EQ. KELEM ) KPNTR = JI
    ENDDO

! any missing elements?
    IF ( KPNTR .EQ. 0  ) THEN
       IF (DEBUG) THEN
          WRITE(6,*) ' XTRDATA: No data for element ',KELEM
       ENDIF
       RETURN
    ELSE 
       IPOS = 0
       DO JI = 1, KNT 
          DO JJ = 1, KNVAL 
             IPOS = IPOS + 1
             INDX = (JI-1)*KNELE*KNVAL + (JJ-1)*KNELE + KPNTR 
             KTBLVAL(INDX) =  KDATA(IPOS) 
          ENDDO
       ENDDO
       IF (DEBUG) THEN
          WRITE(6,*) ' REPDATA: ktblval =  ',(KTBLVAL(JJ), &
                    JJ=1,KNELE*KNVAL*KNT)
       ENDIF
    ENDIF

    RETURN
  END


  SUBROUTINE mwbg_tovCheckAmsua(KSAT,KTERMER,KORBIT,ICANO,ICANOMP, &
                          ZO,ZCOR,ZOMP,ICHECK,KNO,KNOMP, &
                          KNT,PMISG,KNOSAT,KCHKPRF,ISCNPOS, &
                          MGINTRP,MTINTRP,GLINTRP,ITERRAIN, &
                          SATZEN,IMARQ,STNID, &
                          RESETQC,ZLAT)
!OBJET          Effectuer le controle de qualite des radiances tovs.
!ARGUMENTS      ksat    - input  -  numero d'identificateur du satellite
!               ktermer - input  -  indicateur terre/mer
!               korbit  - input  -  numero d'orbite
!               icano   - input  -  canaux des observations
!               icanomp - input  -  canaux des residus (o-p)
!               zcor    - input  -  correction aux radiances
!               zo      - input  -  radiances
!               zomp    - input  -  residus (o-p)
!               icheck  - output -  indicateur controle de qualite tovs par canal 
!                                   =0, ok,
!                                   >0, rejet,
!               kno     - input  -  nombre de canaux des observations
!               knomp   - input  -  nombre de canaux des residus (o-p)
!               knt     - input  -  nombre de tovs
!               pmisg   - input  -  valeur manquante burp
!               knosat  - input  -  numero de satellite (i.e. indice)
!               kchkprf - output -  indicateur global controle de qualite tovs. Code:
!                                   =0, ok,
!                                   >0, rejet d'au moins un canal.
!               iscnpos - input  -  position sur le "scan"
!               mgintrp - input  -  masque terre/mer du modele
!               mtintrp - input  -  topographie du modele
!               glintrp - input  -  etendue de glace du modele
!               iterrain- input  -  indicateur du type de terrain
!               satzen  - input  -  angle zenith du satellite (deg.)
!               imarq   - in/out -  marqueurs des radiances
!               stnid   - input  -  identificateur du satellite
!               resetqc - input  -  reset du controle de qualite?
!               zlat    - input  -  latitude
!
!NOTES  
!               Quinze tests sont effectues menant aux erreurs suivantes:
!                  1) topography reject,
!                  2) invalid land/sea qualifier,
!                  3) invalid terrain type,
!                  4) invalid field of view number,
!                  5) satellite zenith angle out of range,
!                  6) inconsistent field of view and sat. zenith angle,
!                  7) inconsistent land/sea qualifier and model mask,
!                  8) inconsistent terrain type and model ice, (NOT USED)
!                  9) uncorrected radiance,
!                 10) rejected by RTTOV,
!                 11) radiance gross check failure,
!                 12) cloud liquid water reject,
!                 13) scattering index reject,
!                 14) radiance residual rogue check failure,
!                 15) channel reject (channel selection).
    IMPLICIT NONE

    INTEGER     MXNT
    PARAMETER ( MXNT   =  2000 )

    INTEGER MXCHN, MXSAT, MXSCAN, MXCLWREJ, MXCANPRED, MXSFCREJ2
    INTEGER MXSCANHIRS, MXSCANAMSU, MXSCATREJ, MXSFCREJ, NTESTS
    PARAMETER  ( MXCHN     = 42 )
    PARAMETER  ( MXSAT     =  9 )
    PARAMETER  ( MXSCAN    = 56 )
    PARAMETER  ( MXSCANHIRS= 56 )
    PARAMETER  ( MXSCANAMSU= 30 )
    PARAMETER  ( MXCLWREJ  =  6 )
    PARAMETER  ( MXSFCREJ  =  6 )
    PARAMETER  ( MXSFCREJ2 =  4 )
    PARAMETER  ( MXSCATREJ =  7 )
    PARAMETER  ( MXCANPRED =  9 )


    INTEGER JPNSAT,JPCH
 
    PARAMETER (JPNSAT =  9) 
    PARAMETER (JPCH = 50)

    INTEGER JPMXSFC, JPMXREJ
    PARAMETER (JPMXSFC =  2)
    PARAMETER (JPMXREJ = 15)
    
    INTEGER NCHNA   (     JPNSAT)
    INTEGER MLISCHNA(JPCH,JPNSAT)
    INTEGER IUTILST (JPCH,JPNSAT)
    REAL    TOVERRST(JPCH,JPNSAT)
    
    COMMON /COMTOVST/ NCHNA, MLISCHNA, TOVERRST, IUTILST

    INTEGER KNO,KNOMP,KNT,KNOSAT,MAXVAL
    INTEGER JI,JJ,INDX8,INDX12,INO,ICHN
    INTEGER JK,IBIT,JC,INDX,INDXCAN
    INTEGER ITRN

    INTEGER KSAT    (KNT)
    INTEGER KTERMER (:)
    INTEGER ISCNPOS (:)
    INTEGER KORBIT  (KNT)
    INTEGER ICANO   (MXVAL*MXNT)
    INTEGER KCANO   (KNO  ,KNT)
    INTEGER ICANOMP (MXVAL*MXNT)
    INTEGER KCANOMP (KNOMP,KNT)
    INTEGER ICHECK  (KNO  ,KNT)
    INTEGER KCHKPRF (KNT)
    INTEGER IMARQ    (MXVAL*MXNT)
    INTEGER KMARQ   (KNO  ,KNT)
    INTEGER ITERRAIN(KNT)
    INTEGER ICLWREJ (MXCLWREJ)
    INTEGER ISFCREJ (MXSFCREJ)
    INTEGER ISFCREJ2(MXSFCREJ2)
    INTEGER ISCATREJ(MXSCATREJ)

    REAL  PMISG,ZSEUILCLW,EPSILON,MISGINT,ZANGL,MISGRODY,ZSEUILSCAT
    REAL  APPROXIM, ANGDIF, XCHECKVAL
    REAL  ZO      (MXVAL*MXNT)
    REAL  PTBO    (KNO    ,KNT)
    REAL  ZCOR    (MXVAL*MXNT)
    REAL  PTBCOR  (KNO    ,KNT)
    REAL  ZOMP    (MXVAL*MXNT)
    REAL  PTBOMP  (KNOMP  ,KNT)
    REAL  MGINTRP (:)
    REAL  MTINTRP (:)
    REAL  GLINTRP (:)
    REAL  SATZEN  (:)
    REAL  ZLAT    (:)
    REAL  GROSSMIN(MXCHN)
    REAL  GROSSMAX(MXCHN) 
    REAL  ROGUEFAC(MXCHN)

    real tb23 (mxnt)
    real tb31 (mxnt)
    real tb50 (mxnt)
    real tb53 (mxnt)
    real tb89 (mxnt)
    real ice  (mxnt)
    real tpw  (mxnt)
    real clw  (mxnt)
    real scatl(mxnt)
    real scatw(mxnt)

    integer err (mxnt)
    integer rain(mxnt)
    integer snow(mxnt)

    CHARACTER *9   STNID

    LOGICAL LLFIRST,GROSSERROR,FULLREJCT,RESETQC,SFCREJCT

    INTEGER  MREJCOD,INTOT,INTOTRJF, INTOTRJP
    COMMON /STATS/  MREJCOD(JPMXREJ,MXCHN,MXSAT), &
                   INTOT(MXSAT), &
                   INTOTRJF(MXSAT),INTOTRJP(MXSAT)

    LOGICAL DEBUG
    COMMON /DBGCOM/ DEBUG

    SAVE LLFIRST

    DATA  LLFIRST / .TRUE. /
    DATA  EPSILON / 0.01   /
    DATA  MISGINT / -1     /
    DATA  MISGRODY / -99.     /
!      DATA  ROGUEFAC/ 3.0    / changed, jh, from 3 to 4, jan 2001
!      DATA  ROGUEFAC/ 4.0    /
    DATA  ROGUEFAC / 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
                     4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
                     4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 2.0, 2.0, 2.0, &
                     3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
                     4.0, 2.0/
    DATA  ICLWREJ  / 28, 29, 30, 31, 32, 42 /
    DATA  ISFCREJ  / 28, 29, 30, 31, 32, 42 /
    DATA  ISCATREJ / 28, 29, 30, 31, 32, 33, 42 /
    DATA  ISFCREJ2 / 28, 29, 30, 42 /
                   
    DATA GROSSMIN / 200., 190., 190., 180., 180., 180., 170., &
                    170., 180., 170., 170., 170., 180., 180., &
                    180., 180., 170., 180., 180., 000., 120., &
                    190., 180., 180., 180., 190., 200., 120., &
                    120., 160., 190., 190., 200., 190., 180., &
                    180., 180., 180., 190., 190., 200., 130./

    DATA GROSSMAX / 270., 250., 250., 250., 260., 280., 290., &
                    320., 300., 320., 300., 280., 320., 300., &
                    290., 280., 330., 350., 350., 000., 310., &
                    300., 250., 250., 270., 280., 290., 310., &
                    310., 310., 300., 300., 260., 250., 250., &
                    250., 260., 260., 270., 280., 290., 330./  

    ! copy the original input 1D array to 2D array. The 2D arrays are used in this s/r.
    DO JJ=1,KNT
      DO JI=1,KNO
        INDX = (JJ-1)*KNO + JI 
        KCANO(JI,JJ) = ICANO(INDX)
        PTBCOR(JI,JJ) = ZCOR(INDX)  
        PTBO(JI,JJ) = ZO(INDX)  
        KMARQ(JI,JJ) = IMARQ(INDX)  
      end do
    end do

    DO JJ=1,KNT
      DO JI=1,KNOMP
        INDX = (JJ-1)*KNOMP + JI 
        KCANOMP(JI,JJ) = ICANOMP(INDX)
        PTBOMP(JI,JJ) = ZOMP(INDX)  
      end do
    end do


! Initialisation, la premiere fois seulement!
    IF (LLFIRST) THEN
       DO JI = 1, JPMXREJ
          DO JJ = 1, MXCHN
             DO JK = 1, MXSAT
                MREJCOD(JI,JJ,JK) = 0
             ENDDO
          ENDDO
       ENDDO
       LLFIRST = .FALSE.
    ENDIF

! Verification de l'integrite des donnees, c'est-a-dire que:
!         i)  les dimensions des divers blocs de donnees concordent,
!         ii) les listes de canaux concordent.
    IF ( KNO .NE. KNOMP     ) THEN
       WRITE(6,*)'ERROR IN DIMENSIONS OF TOVS DATA'
       CALL ABORT()
    ENDIF

    DO JJ=1,KNT
       DO JI=1,KNO
          IF ( KCANO(JI,JJ) .NE. KCANOMP(JI,JJ) ) THEN
             WRITE(6,*)'INCONSISTENT CHANNEL LISTS FOR TOVS DATA'
             CALL ABORT()
          ENDIF
       ENDDO
    ENDDO

! Initialisations
    DO JJ=1,KNT
       DO JI=1,KNO
          ICHECK(JI,JJ) = 0
          IF ( RESETQC ) KMARQ(JI,JJ) = 0
       ENDDO
    ENDDO

!  Run Grody AMSU-A algorithms.
!     Grody parameters.
!     extract required channels:
!        23 Ghz = AMSU-A 1 = channel #28
!        31 Ghz = AMSU-A 2 = channel #29
!        50 Ghz = AMSU-A 3 = channel #30
!        53 Ghz = AMSU-A 5 = channel #32
!        89 Ghz = AMSU-A15 = channel #42
    DO JJ=1,KNT
       DO JI=1,KNO
          ichn = KCANO(JI,JJ)
          if ( ptbo(ji,jj) .ne. pmisg ) then
             if ( ptbcor(ji,jj) .ne. pmisg ) then
                if ( ichn .eq. 28 ) tb23(jj) = ptbo(ji,jj) &
                     - ptbcor(ji,jj)
                if ( ichn .eq. 29 ) tb31(jj) = ptbo(ji,jj) &
                     - ptbcor(ji,jj)
                if ( ichn .eq. 30 ) tb50(jj) = ptbo(ji,jj) &
                     - ptbcor(ji,jj)
                if ( ichn .eq. 32 ) tb53(jj) = ptbo(ji,jj) &
                     - ptbcor(ji,jj)
                if ( ichn .eq. 42 ) tb89(jj) = ptbo(ji,jj) &
                     - ptbcor(ji,jj)
             else
                if ( ichn .eq. 28 ) tb23(jj) = ptbo(ji,jj)
                if ( ichn .eq. 29 ) tb31(jj) = ptbo(ji,jj)
                if ( ichn .eq. 30 ) tb50(jj) = ptbo(ji,jj)
                if ( ichn .eq. 32 ) tb53(jj) = ptbo(ji,jj)
                if ( ichn .eq. 42 ) tb89(jj) = ptbo(ji,jj)
             endif
          else
             if ( ichn .eq. 28 ) tb23(jj) = 0.
             if ( ichn .eq. 29 ) tb31(jj) = 0.
             if ( ichn .eq. 30 ) tb50(jj) = 0.
             if ( ichn .eq. 32 ) tb53(jj) = 0.
             if ( ichn .eq. 42 ) tb89(jj) = 0.
          endif
       ENDDO
    ENDDO

    call grody (err, knt, tb23, tb31, tb50, tb53, tb89, &
                satzen, zlat, ktermer, ice, &
                tpw, clw, rain, snow, scatl, scatw)   

! 10) test 10: RTTOV reject check.                                 (single)
! Rejected datum flag has bit #9 on.
    IF (.NOT.RESETQC) THEN
       INO = 10
       DO JJ=1,KNT
          DO JI=1,KNO
            IF ( KCANO(JI,JJ) .NE. 20 ) THEN
               IBIT = AND(KMARQ(JI,JJ), 2**9)
               IF ( IBIT .NE. 0  ) THEN
                  ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
                  KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
                  MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                       MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
                  IF (DEBUG) THEN
                     WRITE(6,*)STNID(2:9),' RTTOV REJECT.', &
                              'CHANNEL=', KCANO(JI,JJ), &
                              ' IMARQ= ',KMARQ(JI,JJ)
                  ENDIF
               ENDIF
            ENDIF
          ENDDO
       ENDDO
    ENDIF

! 1) test 1: Topography check                                     (partial)
! Channel 6 is rejected for topography >  250m.
! Channel 7 is rejected for topography > 2000m.
    INO = 1
    DO JJ=1,KNT
       DO JI=1,KNO
         IF ( KCANO(JI,JJ) .EQ. 33 ) THEN
            IF ( MTINTRP(JJ) .GE. 250.  ) THEN
               ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
               KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
               KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**18)
               MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                    MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
               IF (DEBUG) THEN
                  WRITE(6,*)STNID(2:9),' TOPOGRAPHY REJECT.', &
                           'CHANNEL=', KCANO(JI,JJ), &
                           ' TOPO= ',MTINTRP(JJ)
               ENDIF
            ENDIF
         ELSEIF ( KCANO(JI,JJ) .EQ. 34 ) THEN
            IF ( MTINTRP(JJ) .GE. 2000.  ) THEN
               ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
               KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
               KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**18)
               MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                    MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
               IF (DEBUG) THEN
                  WRITE(6,*)STNID(2:9),' TOPOGRAPHY REJECT.', &
                           'CHANNEL=', KCANO(JI,JJ), &
                           ' TOPO= ',MTINTRP(JJ)
               ENDIF
            ENDIF
         ENDIF
       ENDDO
    ENDDO

! 2) test 2: "Land/sea qualifier" code check.                        (full)
! allowed values are: 0, land,
!                       1, sea,
!                       2, coast.
    INO = 2
    DO JJ=1,KNT
       IF ( KTERMER(JJ) .LT.  0  .OR. &
           KTERMER(JJ) .GT.  2        ) THEN
          DO JI=1,KNO
             ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
             KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
             KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
             MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
               MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
          ENDDO
          IF (DEBUG) THEN
             WRITE(6,*) STNID(2:9),'LAND/SEA QUALIFIER CODE', &
                     ' REJECT. KTERMER=', KTERMER(JJ)
          ENDIF
       ENDIF
    ENDDO

! 3) test 3: "Terrain type" code check.                              (full)
!   allowed values are: -1, missing,
!                        0, sea-ice,
!                        1, snow on land.
    INO = 3
    DO JJ=1,KNT
       IF ( ITERRAIN(JJ) .NE.  MISGINT ) THEN
          IF ( ITERRAIN(JJ) .LT.  0  .OR. &
              ITERRAIN(JJ) .GT.  1        ) THEN
             DO JI=1,KNO
                ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
                KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
                KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
                MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                  MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
             ENDDO
             IF (DEBUG) THEN
                WRITE(6,*)STNID(2:9),'TERRAIN TYPE CODE', &
                        ' REJECT. TERRAIN=', ITERRAIN(JJ)
             ENDIF
          ENDIF
       ENDIF
    ENDDO

! 4) test 4: Field of view number check.                             (full)
!
! Field of view acceptable range is [1,MXSCANAMSU]  for AMSU footprints.
    INO = 4
    DO JJ=1,KNT
       DO JI=1,KNO
          IF ( ISCNPOS(JJ) .LT. 1          .OR. &
              ISCNPOS(JJ) .GT. MXSCANAMSU       ) THEN
             ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
             KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
             KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
             MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
               MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
             IF (DEBUG) THEN
                WRITE(6,*)STNID(2:9),'FIELD OF VIEW NUMBER', &
                         ' REJECT. FIELD OF VIEW= ', &
                         ISCNPOS(JJ)
             ENDIF
          ENDIF
       ENDDO
    ENDDO

! 5) test 5: Satellite zenith angle check.                           (full)
! Satellite zenith angle acceptable range is [0.,60.].
    INO = 5
    DO JJ=1,KNT
       IF ( SATZEN(JJ) .NE.  PMISG ) THEN
          IF ( SATZEN(JJ) .LT.  0.  .OR. &
              SATZEN(JJ) .GT. 60.       ) THEN
             DO JI=1,KNO
                ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
                KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
                KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
                MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                  MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
             ENDDO
             IF (DEBUG) THEN
                WRITE(6,*)STNID(2:9),' SATELLITE ZENITH ANGLE', &
                         ' REJECT. SATZEN= ', &
                         SATZEN(JJ)
             ENDIF
          ENDIF
       ENDIF
    ENDDO

! 6) test 6: "Sat. zenith angle"/"field of view" consistency check.  (full)
! Acceptable difference between "Satellite zenith angle"  and
! "approximate angle computed from field of view number" is 1.8 degrees.
    ZANGL   = 3.92

    INO = 6
    DO JJ=1,KNT
       IF ( SATZEN (JJ) .NE.  PMISG   .AND. &
           ISCNPOS(JJ) .NE.  MISGINT       ) THEN
          APPROXIM = ABS((ISCNPOS(JJ)-MXSCANAMSU/2.-0.5)*ZANGL)
          ANGDIF = ABS(SATZEN (JJ)-APPROXIM)
          IF ( ANGDIF .GT. 1.8 ) THEN 
             DO JI=1,KNO
                ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
                KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
                KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
                MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                  MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
             ENDDO
             IF (DEBUG) THEN
                WRITE(6,*)STNID(2:9),' ANGLE/FIELD OF VIEW', &
                         ' INCONSISTENCY REJECT. SATZEN= ', &
                         SATZEN(JJ), &
                         ' FIELD OF VIEW= ',ISCNPOS(JJ), &
                         ' ANGDIF= ',ANGDIF  
             ENDIF
          ENDIF
       ENDIF
    ENDDO

! 7) test 7: "Land/sea qual."/"model land/sea" consistency check.    (full)
! Acceptable conditions are:
!       a) both over ocean (ktermer=1; mg<0.01), new threshold 0.20, jh dec 2000,
!       b) both over land  (ktermer=0; mg>0.80), new threshold 0.50, jh dec 2000.
! Other conditions are unacceptable.
    INO = 7
    DO JJ=1,KNT
       IF ( KTERMER (JJ) .NE.  MISGINT  ) THEN
          IF     ( KTERMER(JJ) .EQ. 1       .AND. &
                  MGINTRP(JJ) .LT. 0.20          ) THEN
          ELSEIF ( KTERMER(JJ) .EQ. 0       .AND. &
                  MGINTRP(JJ) .GT. 0.50          ) THEN
          ELSE
             DO JI=1,KNO
                ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
                KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
                KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
                MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                  MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
             ENDDO
             IF (DEBUG) THEN
                WRITE(6,*)STNID(2:9),' LAND/SEA QUALIFIER', &
                         ' INCONSISTENCY REJECT. KTERMER= ', &
                         KTERMER(JJ), &
                         ' MODEL MASK= ',MGINTRP(JJ)
             ENDIF
          ENDIF
       ENDIF
    ENDDO

! 8) test 8: "Terrain type"/"Land/sea qual."/"model ice" consistency check.           (full)
! Unacceptable conditions are:
!        a) terrain is sea-ice and model has no ice(iterrain=0; gl<0.01).
!        b) terrain is sea-ice and land/sea qualifier is land (iterrain=0; ktermer=0).
!        c) terrain is snow on land and land/sea qualifier is sea (iterrain=1; ktermer=1).
!        d) terrain is missing, land/sea qualifier is sea and model has ice(iterrain=-1; ktermer=1; gl>0.01). (enleve jh, jan 2001)
!
! skip this test, terrain type in burp file now comes from operational analysis. (jh, march 2003)
    go to 800 
    INO = 8
    DO JJ=1,KNT
       IF ( ITERRAIN (JJ) .NE.  MISGINT  ) THEN
          IF     ( ITERRAIN(JJ) .EQ. 0       .AND. &
                  GLINTRP (JJ) .LT. 0.01          ) THEN
             DO JI=1,KNO
                ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
                MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                  MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
             ENDDO
             IF (DEBUG) THEN
                WRITE(6,*)STNID(2:9),' TERRAIN TYPE/MODEL ICE', &
                         ' INCONSISTENCY REJECT. TERRAIN= ', &
                         ITERRAIN(JJ), &
                         ' MODEL ICE= ',GLINTRP(JJ)
             ENDIF
          ENDIF
          IF     ( ITERRAIN(JJ) .EQ. 0       .AND. &
                  KTERMER (JJ) .EQ. 0           ) THEN
             DO JI=1,KNO
                ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
                MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                  MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
             ENDDO
             IF (DEBUG) THEN
                WRITE(6,*)STNID(2:9),' TERRAIN TYPE/LAND?SEA QUAL.', &
                         ' INCONSISTENCY REJECT. TERRAIN= ', &
                         ITERRAIN(JJ), &
                         ' LAND/SEA= ',KTERMER(JJ)
             ENDIF
          ENDIF
          IF     ( ITERRAIN(JJ) .EQ. 1       .AND. &
                  KTERMER (JJ) .EQ. 1           ) THEN
             DO JI=1,KNO
                ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
                MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                  MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
             ENDDO
             IF (DEBUG) THEN
                WRITE(6,*)STNID(2:9),' TERRAIN TYPE/LAND?SEA QUAL.', &
                         ' INCONSISTENCY REJECT. TERRAIN= ', &
                         ITERRAIN(JJ), &
                         ' LAND/SEA= ',KTERMER(JJ)
             ENDIF
          ENDIF
       ELSE
          IF     ( KTERMER (JJ) .EQ. 1      .AND. &
                  GLINTRP (JJ) .GT. 0.01         ) THEN
!               DO JI=1,KNO
!                  ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
!                  MREJCOD(INO,KCANO(JI,JJ),KNOSAT) =
!                    MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
!               ENDDO
             IF (DEBUG) THEN
                WRITE(6,*)STNID(2:9),' TERRAIN TYPE MSG/MODEL ICE', &
                         ' INCONSISTENCY REJECT. TERRAIN= ', &
                         ITERRAIN(JJ), &
                         ' LAND/SEA= ',KTERMER(JJ)
             ENDIF
          ENDIF
       ENDIF
    ENDDO
800  continue

! 9) test 9: Uncorrected Tb check.                                 (single)
! Uncorrected datum (flag bit #6 off). In this case switch bit 11 ON.
    IF (.NOT.RESETQC) THEN
    INO = 9
    DO JJ=1,KNT
       DO JI=1,KNO
         IF ( KCANO(JI,JJ) .NE. 20 ) THEN
            IBIT = AND(KMARQ(JI,JJ), 2**6)
            IF ( IBIT .EQ. 0  ) THEN
               ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
               KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**11)
               MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                    MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
               IF (DEBUG) THEN
                  WRITE(6,*)STNID(2:9),' UNCORRECTED TB REJECT.', &
                           'CHANNEL=', KCANO(JI,JJ), &
                           ' IMARQ= ',KMARQ(JI,JJ)
               ENDIF
            ENDIF
         ENDIF
       ENDDO
    ENDDO
    ENDIF

! 11) test 11: Radiance observation "Gross" check.                 (single) 
!  Change this test from full to single. jh nov 2000.
    INO = 11
    DO JJ=1,KNT
       GROSSERROR = .FALSE.
       DO JI=1,KNO
         IF ( KCANO(JI,JJ) .NE. 20     .AND. &
             KCANO(JI,JJ) .GE.  1     .AND. &
             KCANO(JI,JJ) .LE.  MXCHN       ) THEN  
          IF ( PTBO(JI,JJ) .NE. PMISG .AND. &
           ( PTBO(JI,JJ).LT.GROSSMIN(KCANO(JI,JJ)).OR. &
             PTBO(JI,JJ).GT.GROSSMAX(KCANO(JI,JJ))     ) ) THEN
             GROSSERROR = .TRUE.
             ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
             KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
             KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
             MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                    MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
             IF (DEBUG) THEN
                WRITE(6,*)STNID(2:9),' GROSS CHECK REJECT.', &
                         'CHANNEL=', KCANO(JI,JJ), &
                         ' TB= ',PTBO(JI,JJ)
             ENDIF
          ENDIF
         ENDIF
       ENDDO
    ENDDO

! 12) test 12: Grody cloud liquid water check.                    (partial)
! For Cloud Liquid Water > 0.3, reject AMSUA-A channels 1-5 and 15.
    INO = 12
    ZSEUILCLW = 0.3
    DO JJ=1,KNT
       IF ( CLW(JJ) .NE.  MISGRODY  ) THEN
          IF ( CLW(JJ) .GT. ZSEUILCLW   ) THEN
             DO JI=1,KNO
                INDXCAN = ISRCHEQI (ICLWREJ,MXCLWREJ,KCANO(JI,JJ))
                IF ( INDXCAN.NE.0 )  THEN
                   ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
                   KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
                   KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
                   MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                           MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
                ENDIF
             ENDDO
             IF (DEBUG) THEN
!                 IF (.true.) THEN
                WRITE(6,*)STNID(2:9),'Grody cloud liquid water check', &
                         ' REJECT. CLW= ',CLW(JJ), &
                         ' SEUIL= ',ZSEUILCLW
             ENDIF
          ENDIF
       ENDIF
    ENDDO

! 13) test 13: Grody scattering index check.                      (partial)
! For Scattering Index > 9, reject AMSUA-A channels 1-6 and 15.
    INO = 13
    ZSEUILSCAT = 9.0
    DO JJ=1,KNT
       IF ( SCATW(JJ) .NE.  MISGRODY  ) THEN
          IF (  KTERMER (JJ) .EQ.  1 .AND. &
               ITERRAIN(JJ) .NE.  0 .AND. &   
               SCATW   (JJ) .GT. ZSEUILSCAT   ) THEN
             DO JI=1,KNO
                INDXCAN = ISRCHEQI (ISCATREJ,MXSCATREJ,KCANO(JI,JJ))
                IF ( INDXCAN.NE.0 )  THEN
                   ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
                   KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
                   KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
                   MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                           MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
                ENDIF
             ENDDO
             IF (DEBUG) THEN
!                 IF (.true.) THEN
                WRITE(6,*)STNID(2:9),'Grody scattering index check', &
                          ' REJECT. SCATW= ',SCATW(JJ), &
                          ' SEUIL= ',ZSEUILSCAT
             ENDIF
          ENDIF
       ENDIF
    ENDDO

! 14) test 14: "Rogue check" for (O-P) Tb residuals out of range. (single/full)
! Les observations, dont le residu (O-P) depasse par un facteur (roguefac) l'erreur totale des TOVS.
! N.B.: a reject by any of the 3 surface channels produces the rejection of AMSUA-A channels 1-5 and 15. 
    INO = 14
    DO JJ=1,KNT

       SFCREJCT = .FALSE.
       DO JI=1,KNO
          ICHN = KCANO(JI,JJ)
          IF ( ICHN .NE. 20 ) THEN
             XCHECKVAL = ROGUEFAC(ICHN)*TOVERRST(ICHN,KNOSAT) 
             IF ( PTBOMP(JI,JJ)      .NE. PMISG    .AND. &
                  ABS(PTBOMP(JI,JJ)) .GE. XCHECKVAL     ) THEN
                ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
                KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
                KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**16)
                MREJCOD(INO,ICHN,KNOSAT) = &
                    MREJCOD(INO,ICHN,KNOSAT) + 1 
                IF (DEBUG) THEN
                   WRITE(6,*)STNID(2:9),'ROGUE CHECK REJECT.NO.', &
                          ' OBS = ',JJ, &
                          ' CHANNEL= ',ICHN, &
                          ' CHECK VALUE= ',XCHECKVAL, &
                          ' TBOMP= ',PTBOMP(JI,JJ)
                ENDIF
                IF ( ICHN .EQ. 28 .OR. &
                     ICHN .EQ. 29 .OR. &
                     ICHN .EQ. 30      ) THEN
                   SFCREJCT = .TRUE.
                ENDIF
             ENDIF
          ENDIF
       ENDDO

       IF ( SFCREJCT ) THEN
          DO JI=1,KNO
             INDXCAN = ISRCHEQI (ISFCREJ,MXSFCREJ,KCANO(JI,JJ))
             IF ( INDXCAN .NE. 0 )  THEN
                IF ( ICHECK(JI,JJ) .NE. INO ) THEN
                   ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
                   KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
                   KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**16)
                   MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                            MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
                ENDIF
             ENDIF
          ENDDO
       ENDIF

    ENDDO

! 15) test 15: Channel Selection using array IUTILST(chan,sat)
!  IUTILST = 0 (blacklisted)
!            1 (assmilate)
!            2 (assimilate over open water only)
!
!  We also set QC flag bits 7 and 9 ON for channels with IUTILST=2 
!  over land or sea-ice
!    and 
!  we set QC flag bits 7 and 9 ON for channels 1-3,15 over land
!  or sea-ice REGARDLESS of IUTILST value (but IUTILST=0 always for
!  these unassimilated channels).
    INO = 15

    DO JJ=1,KNT
      ITRN = ITERRAIN(JJ)
      IF ( KTERMER (JJ) .EQ. 1    .AND. &
           ITERRAIN(JJ) .EQ. -1   .AND. &
           GLINTRP (JJ) .GE. 0.01       ) THEN
         ITRN = 0
      ENDIF        
      DO JI=1,KNO
          ICHN = KCANO(JI,JJ)
          INDXCAN = ISRCHEQI (ISFCREJ2,MXSFCREJ2,ICHN)
          IF ( INDXCAN .NE. 0 )  THEN
            IF ( KTERMER (JJ) .EQ. 0 .OR. ITRN .EQ. 0 )  THEN
              KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
              KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
            ENDIF
          ENDIF
          IF ( IUTILST(ICHN,KNOSAT) .NE. 1 ) THEN
            SFCREJCT = .FALSE.
            IF ( IUTILST(ICHN,KNOSAT) .EQ. 0 ) THEN
              SFCREJCT = .TRUE.
              KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**11)
            ELSE 
              IF ( KTERMER (JJ) .EQ. 0 .OR. ITRN .EQ. 0 )  THEN
                SFCREJCT = .TRUE.
                KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
                KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
              ENDIF
            ENDIF
            IF ( SFCREJCT ) THEN
              ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
              MREJCOD(INO,ICHN,KNOSAT) = & 
                 MREJCOD(INO,ICHN,KNOSAT) + 1 
              IF (DEBUG) THEN
                 WRITE(6,*)STNID(2:9),'CHANNEL REJECT: ', &
                        ' OBS = ',JJ, &
                        ' CHANNEL= ',ICHN
              ENDIF
            ENDIF
          ENDIF
        ENDDO
    ENDDO

    IF (DEBUG) THEN
       WRITE(6,*)'ICHECK = ',((ICHECK(JI,JJ),JI=1,KNO),JJ=1,KNT)
    ENDIF

!  Synthese de la controle de qualite au niveau de chaque point
!  d'observation. Code:
!            =0, aucun rejet,
!            >0, au moins un canal rejete.

    DO JJ=1,KNT
       KCHKPRF(JJ) = 0
       DO JI=1,KNO 
             KCHKPRF(JJ) = MAX(KCHKPRF(JJ),ICHECK(JI,JJ))
       ENDDO
    ENDDO

    IF (DEBUG) THEN
       WRITE(6,*)'KCHKPRF = ',(KCHKPRF(JJ),JJ=1,KNT)
    ENDIF 

    ! Copy the modified FLAG to the 1D array, used outside this s/r. 
    DO JJ=1,KNT
      DO JI=1,KNO
        INDX = (JJ-1)*KNO + JI 
        IMARQ(INDX) = KMARQ(JI,JJ)
      end do
    end do

    RETURN
  END


  SUBROUTINE mwbg_qcStatsAmsua(INUMSAT,ICHECK,KCANO,KNOSAT, &
                        CSATID,KNO,KNT,LDPRINT)
! OBJET          Cumuler ou imprimer des statistiques decriptives
!                des rejets tovs.
!ARGUMENTS      
!               inumsat - input  -  nombre actuel de satellites
!               icheck  - input  -  indicateur controle de qualite tovs par canal 
!                                   =0, ok,
!                                   >0, rejet,
!               kcano   - input  -  canaux des observations
!               knosat  - input  -  numero du satellite (i.e. index) 
!               csatid  - input  -  liste des noms de satellite 
!               kno     - input  -  nombre de canaux des observations
!               knt     - input  -  nombre de tovs
!               ldprint - input  -  mode: imprimer ou cumuler? 
    IMPLICIT NONE

    INTEGER MXCHN, MXSAT
    PARAMETER  ( MXCHN = 42 )
    PARAMETER  ( MXSAT =  9 )

    INTEGER     JPMXREJ
    PARAMETER ( JPMXREJ = 15)

    CHARACTER*9   CSATID (MXSAT)

    INTEGER  JI, JJ, JK, KNO, KNT, KNOSAT
    INTEGER  INTOTOBS, INTOTACC, INUMSAT

    INTEGER  ICHECK (KNO,KNT)
    INTEGER  KCANO  (KNO,KNT)

    INTEGER  MREJCOD,INTOT,INTOTRJF,INTOTRJP
    COMMON /STATS/  MREJCOD(JPMXREJ,MXCHN,MXSAT), &
                   INTOT(MXSAT), &
                   INTOTRJF(MXSAT),INTOTRJP(MXSAT)

    LOGICAL LLFIRST, LDPRINT, FULLREJCT, FULLACCPT

    SAVE  LLFIRST

    DATA  LLFIRST / .TRUE. /

! Initialize

    IF ( LLFIRST ) THEN
       DO JJ = 1, MXSAT
          INTOTRJF(JJ) = 0
          INTOTRJP(JJ) = 0
          INTOT(JJ)  = 0
       ENDDO
       LLFIRST = .FALSE.
    ENDIF

    IF (.NOT. LDPRINT ) THEN

! Accumulate statistics on rejects

       DO JJ = 1, KNT

          INTOT(KNOSAT) = INTOT(KNOSAT) + 1
! Fully accepted, fully rejected or partially rejected?
          FULLREJCT = .TRUE.
          FULLACCPT = .TRUE.
          DO JI = 1, KNO
             IF ( KCANO(JI,JJ) .NE. 20 ) THEN
                IF ( ICHECK(JI,JJ) .NE. 0 ) THEN
                   FULLACCPT = .FALSE.
                ELSE
                   FULLREJCT = .FALSE.
                ENDIF
             ENDIF
          ENDDO
          IF ( FULLREJCT ) THEN
             INTOTRJF(KNOSAT) = INTOTRJF(KNOSAT) + 1
          ENDIF
          IF ( .NOT.FULLREJCT .AND. .NOT.FULLACCPT ) THEN
             INTOTRJP(KNOSAT) = INTOTRJP(KNOSAT) + 1
          ENDIF
       ENDDO

    ELSE

! Print statistics
      DO JK = 1, INUMSAT

       INTOTOBS = INTOT(JK)
       INTOTACC = INTOTOBS - INTOTRJF(JK) - INTOTRJP(JK)

       WRITE(6,'(/////50("*"))')
       WRITE(6,'(     50("*")/)')
       WRITE(6,'(T5,"SUMMARY OF QUALITY CONTROL FOR ", &
        A8)') CSATID(JK) 
       WRITE(6,'(T5,"------------------------------------- ",/)')
       WRITE(6,'( &
        "   TOTAL NUMBER OF AMSU-A  = ",I10,/ &
        " - TOTAL FULL REJECTS      = ",I10,/ &
        " - TOTAL PARTIAL REJECTS   = ",I10,/ &
        "   ------------------------------------",/ &
        "   TOTAL FULLY ACCEPTED    = ",I10,/)') &
         INTOTOBS, INTOTRJF(JK), INTOTRJP(JK), INTOTACC

       WRITE(6,'(//,1x,114("-"))')
       WRITE(6,'(t10,"|",t47,"REJECTION CATEGORIES")')
       WRITE(6,'(" CHANNEL",t10,"|",105("-"))')
       WRITE(6,'(t10,"|",15i7)') (JI,JI=1,JPMXREJ)
       WRITE(6,'(1x,"--------|",105("-"))')
       DO JJ = 1, MXCHN
          WRITE(6,'(3X,I2,t10,"|",15I7)') JJ,(MREJCOD(JI,JJ,JK), &
                                     JI=1,JPMXREJ)
       ENDDO
       WRITE(6,'(1x,114("-"))')
     ENDDO

! Print legend
     PRINT *, ' '
     PRINT *, ' '
     PRINT *, ' -----------------------------------------------------'
     PRINT *, ' Definition of rejection categories:'
     PRINT *, ' -----------------------------------------------------'
     PRINT *, '  1 - topography reject'
     PRINT *, '  2 - invalid land/sea qualifier'
     PRINT *, '  3 - invalid terrain type'
     PRINT *, '  4 - invalid field of view number'
     PRINT *, '  5 - satellite zenith angle out of range '
     PRINT *, '  6 - inconsistent field of view and sat. zenith angle'
     PRINT *, '  7 - inconsistent land/sea qualifier and model mask'
     PRINT *, '  8 - inconsistent terrain type and land/sea', &
             ' qualifier/model ice (NOT DONE)'
     PRINT *, '  9 - uncorrected radiance'
     PRINT *, ' 10 - rejected by RTTOV'
     PRINT *, ' 11 - radiance gross check failure'
     PRINT *, ' 12 - cloud liquid water reject'
     PRINT *, ' 13 - scattering index reject'
     PRINT *, ' 14 - radiance residual rogue check failure'
     PRINT *, ' 15 - rejection by channel selection'
     PRINT *, ' -----------------------------------------------------'
     PRINT *, ' '

    ENDIF

    RETURN
  END


  SUBROUTINE mwbg_UPDATFLG(KBUF1,KLISTE,KTBLVAL,KDLISTE,PRVAL, &
                      KDATA,PDATA,KCHKPRF,KTERMER,ICHECK, &
                      RESETQC,IMARQ)
!OBJET          Allumer les bits des marqueurs pour les tovs rejetes.
!               Mettre a jour l'indicateur terre/mer qui a
!               possiblement ete corrige pour la glace marine.
!               Modifier le bktyp des donnees, marqueurs et (O-P) pourt
!               signifier "vu par AO". 
!
!APPEL          CALL   mwbg_UPDATFLG (KBUF1,KLISTE,KTBLVAL,KDLISTE,PRVAL,
!                                KDATA,PDATA,KCHKPRF,KTERMER,ICHECK,
!                                RESETQC,IMARQ)
!
!ARGUMENTS      kbuf1   - in/out -  tableau contenant le rapport
!               kliste  - input  -  liste des noms d'elements
!               ktblval - input  -  champ de donnees (valeurs entieres)
!               kdliste - input  -  liste des noms d'elements (decodes)
!               prval   - input  -  champ de donnees (valeurs reelles)
!               kdata   - input  -  donnees extraites (valeurs entieres)
!               pdata   - input  -  donnees extraites (valeurs reelles)
!               kchprf  - input  -  indicateur global controle de qualite tovs. Code:
!                                   =0, ok,
!                                   >0, rejet,
!               ktermer - input  -  indicateur terre/mer
!               icheck  - input  -  indicateur controle de qualite tovs au 
!                                   niveau de chaque canal
!               resetqc - input  -  reset the quality control flags before adding the new ones ? 
!               imarq   - input  -  modified flag values from mwbg_tovCheckAmsua
!!
    IMPLICIT NONE

    INTEGER KBUF1   (:)
    INTEGER KDLISTE (:)
    INTEGER KLISTE  (:)
    INTEGER KTBLVAL (:)
    INTEGER KDATA   (:)
    INTEGER KCHKPRF (:)
    INTEGER KTERMER (:)
    INTEGER ICHECK  (:)
    INTEGER IMARQ   (:)

    INTEGER IDUM1,IDUM2,IDUM3,IBFAM
    INTEGER IBDESC,IBTYP,INBIT,IBIT0,IDATYP
    INTEGER IBKNAT,IBKTYP,IBKSTP 
    INTEGER MRBPRM, MRBDEL, MRBADD, MRBTYP,MRBXTR
    INTEGER INELE,INVAL,INT,JI,IPNTR,MRBREP,MRBLOC
    INTEGER JJ,IBLKNO,ISTAT,IMELERAD

    REAL PRVAL (:)
    REAL PDATA (:)

    LOGICAL RESETQC

    LOGICAL DEBUG
    COMMON /DBGCOM/ DEBUG

! 1) Bloc info 3d: bloc 5120.
!    Modifier les marqueurs globaux de 24bits pour les donnees rejetees.

!  extraire le bloc
    CALL XTRBLK (5120,-1,KBUF1,KLISTE,KTBLVAL,KDLISTE,PRVAL, &
                INELE,INVAL,INT,IBLKNO)    
    IF(IBLKNO .LE. 0) THEN
       WRITE(6,*)'3D INFO BLOCK NOT FOUND'
       CALL ABORT()
    ENDIF

! extraire les marqueurs globaux de 24bits; element 55200
    CALL XTRDATA (KDLISTE,KTBLVAL,PRVAL,INELE,INVAL,INT, &
                 55200,KDATA,PDATA,IPNTR)    
    IF(IPNTR .EQ. 0) THEN
       WRITE(6,*)'GLOBAL FLAGS MISSING'
       CALL ABORT()
    ENDIF
    IF (DEBUG) THEN
        WRITE(6,*) ' OLD FLAGS = ', (KDATA(JJ),JJ=1,INVAL*INT)
    ENDIF

!  allumer la bit (6) indiquant que l'observation a un element
!  rejete par le controle de qualite de l'AO.
!  N.B.: si on est en mode resetqc, on remet le marqueur global a
!        sa valeur de defaut, soit 1024,  avant de faire la mise a jour.
    DO JI = 1, INT
       IF (RESETQC) THEN
          KDATA(JI) = 1024  
       ENDIF
       IF ( KCHKPRF(JI).NE.0  ) THEN
          KDATA(JI) = OR (KDATA(JI),2**6)
       ENDIF
    ENDDO
    IF (DEBUG) THEN
        WRITE(6,*) ' NEW FLAGS = ', (KDATA(JJ),JJ=1,INVAL*INT)
    ENDIF

! Remplacer les nouveaux marqueurs dans le tableau.
    CALL REPDATA (KDLISTE,KTBLVAL,INELE,INVAL,INT, &
                       55200,KDATA,IPNTR)

! Remplacer le bloc.
    ISTAT = MRBREP (KBUF1,IBLKNO,KTBLVAL)

! 2) Bloc info (general): bloc 3072
!    Modifier les indicateurs terre/mer possiblement corriges pour la glace
!    marine.

!  extraire le bloc
    CALL XTRBLK (3072,-1,KBUF1,KLISTE,KTBLVAL,KDLISTE,PRVAL, &
                INELE,INVAL,INT,IBLKNO)    
    IF(IBLKNO .LE. 0) THEN
       WRITE(6,*)'INFO BLOCK NOT FOUND'
       CALL ABORT()
    ENDIF

! extraire l'indicateur terre/mer; element 8012
    CALL XTRDATA (KDLISTE,KTBLVAL,PRVAL,INELE,INVAL,INT, &
                 8012,KDATA,PDATA,IPNTR)    
    IF(IPNTR .EQ. 0) THEN
          WRITE(6,*)'LAND/SEA INDICATOR MISSING'
       CALL ABORT()
    ENDIF
    IF (DEBUG) THEN
        WRITE(6,*) ' OLD LAND/SEA = ', (KDATA(JJ),JJ=1,INVAL*INT)
    ENDIF

! les indicateurs terre/mer corriges
    IF (DEBUG) THEN
        WRITE(6,*) ' NEW LAND/SEA = ', (KTERMER(JJ),JJ=1,INVAL*INT)
    ENDIF

! Remplacer les nouveaux indicateurs terre/mer dans le tableau.
    CALL REPDATA (KDLISTE,KTBLVAL,INELE,INVAL,INT, &
                       8012,KTERMER,IPNTR)

! Remplacer le bloc.
    ISTAT = MRBREP (KBUF1,IBLKNO,KTBLVAL)

! 3) Bloc multi niveaux de radiances: bloc 9218, 9248, 9264.
!    Modifier le bktyp pour signifier "vu par AO".

!  localiser le bloc
    IBLKNO =  MRBLOC(KBUF1,-1,-1,9218,0)   
    IF(IBLKNO .LE. 0) THEN
       IBLKNO =  MRBLOC(KBUF1,-1,-1,9248,0)  
    ENDIF          
    IF(IBLKNO .LE. 0) THEN
       IBLKNO =  MRBLOC(KBUF1,-1,-1,9264,0)  
    ENDIF        
    IF(IBLKNO .LE. 0) THEN
       WRITE(6,*)'RADIANCE DATA BLOCK NOT FOUND'
       CALL ABORT()
    ENDIF

! extraction du bloc

! Remplacer le bloc et modifier le bktyp pour signifier "vu par AO".
    ISTAT = MRBXTR(KBUF1,IBLKNO,KLISTE,KTBLVAL)
    ISTAT = MRBPRM (KBUF1,IBLKNO,INELE,INVAL,INT,IBFAM, &
                   IBDESC,IBTYP,INBIT,IBIT0,IDATYP)
    ISTAT = MRBDEL(KBUF1,IBLKNO)

    ISTAT = MRBTYP(IBKNAT,IBKTYP,IBKSTP,IBTYP)
    IBKTYP = IBKTYP + 4
    IBTYP = MRBTYP(IBKNAT,IBKTYP,IBKSTP,-1)
    ISTAT = MRBADD (KBUF1,IBLKNO,INELE,INVAL,INT,IBFAM, &
                   IBDESC,IBTYP,INBIT,IBIT0,IDATYP, &
                   KLISTE,KTBLVAL)

! 4) Bloc marqueurs multi niveaux de radiances: bloc 15362, 15392, 15408.
!    Modifier les marqueurs de 13bits associes a chaque radiance.
!    Modifier le bktyp pour signifier "vu par AO".

! extraire le bloc
    CALL XTRBLK (15362,-1,KBUF1,KLISTE,KTBLVAL,KDLISTE,PRVAL, &
                INELE,INVAL,INT,IBLKNO)    
    IF(IBLKNO .LE. 0) THEN
       CALL XTRBLK (15392,-1,KBUF1,KLISTE,KTBLVAL,KDLISTE,PRVAL, &
                   INELE,INVAL,INT,IBLKNO)    
    ENDIF        
    IF(IBLKNO .LE. 0) THEN
       CALL XTRBLK (15408,-1,KBUF1,KLISTE,KTBLVAL,KDLISTE,PRVAL, &
                   INELE,INVAL,INT,IBLKNO)    
    ENDIF    
    IF(IBLKNO .LE. 0) THEN
       WRITE(6,*)'RADIANCE DATA FLAG BLOCK NOT FOUND'
       CALL ABORT()
    ENDIF

! extraire les marqueurs de 13bits des radiances; element 212163 (LEVEL 1B)
    IMELERAD =  212163 
    CALL XTRDATA (KDLISTE,KTBLVAL,PRVAL,INELE,INVAL,INT, &
                 IMELERAD,KDATA,PDATA,IPNTR) 
    IF(IPNTR .EQ. 0) THEN
       WRITE(6,*)'RADIANCE DATA FLAGS MISSING'
       CALL ABORT()
    ENDIF
    IF (DEBUG) THEN
        WRITE(6,*) ' OLD FLAGS = ', (KDATA(JJ),JJ=1,INVAL*INT)
    ENDIF

! update data flags
    DO JI = 1, INVAL*INT
        KDATA(JI) = IMARQ(JI)
    ENDDO
    IF (DEBUG) THEN
        WRITE(6,*) ' ICHECK = ', (ICHECK(JJ),JJ=1,INVAL*INT)
        WRITE(6,*) ' NEW FLAGS = ', (KDATA(JJ),JJ=1,INVAL*INT)
    ENDIF

! Remplacer les nouveaux marqueurs dans le tableau.
    CALL REPDATA (KDLISTE,KTBLVAL,INELE,INVAL,INT, &
                 IMELERAD,KDATA,IPNTR)

! Remplacer le bloc et modifier le bktyp pour signifier "vu par AO".
    ISTAT = MRBPRM (KBUF1,IBLKNO,IDUM1,IDUM2,IDUM3,IBFAM, &
                   IBDESC,IBTYP,INBIT,IBIT0,IDATYP)
    ISTAT = MRBDEL(KBUF1,IBLKNO)

    ISTAT = MRBTYP(IBKNAT,IBKTYP,IBKSTP,IBTYP)
    IBKTYP = IBKTYP + 4
    IBTYP = MRBTYP(IBKNAT,IBKTYP,IBKSTP,-1)
    ISTAT = MRBADD (KBUF1,IBLKNO,INELE,INVAL,INT,IBFAM, &
                   IBDESC,IBTYP,INBIT,IBIT0,IDATYP, &
                   KLISTE,KTBLVAL)

! 5) Bloc multi niveaux de residus de radiances (O-P): bloc 9322, 9226, 9258, 9274, bfam 14
!    Modifier le bktyp pour signifier "vu par AO".

!  localiser le bloc
    IBLKNO =  MRBLOC(KBUF1,-1,-1,9322,0)   
    IF(IBLKNO .LE. 0) THEN
       IBLKNO =  MRBLOC(KBUF1,-1,-1,9226,0)  
    ENDIF          
    IF(IBLKNO .LE. 0) THEN
       IBLKNO =  MRBLOC(KBUF1,-1,-1,9258,0)  
    ENDIF        
    IF(IBLKNO .LE. 0) THEN
       IBLKNO =  MRBLOC(KBUF1,-1,-1,9274,0)  
    ENDIF   
    IF(IBLKNO .LE. 0) THEN
       WRITE(6,*)'WARNING: (O-P) RADIANCE BLOCK NOT FOUND'
       CALL ABORT()
    ENDIF

! Remplacer le bloc et modifier le bktyp pour signifier "vu par AO".
    ISTAT = MRBXTR(KBUF1,IBLKNO,KLISTE,KTBLVAL)
    ISTAT = MRBPRM (KBUF1,IBLKNO,INELE,INVAL,INT,IBFAM, &
                   IBDESC,IBTYP,INBIT,IBIT0,IDATYP)
    ISTAT = MRBDEL(KBUF1,IBLKNO)

    ISTAT = MRBTYP(IBKNAT,IBKTYP,IBKSTP,IBTYP)
    IBKTYP = IBKTYP + 4
    IBTYP = MRBTYP(IBKNAT,IBKTYP,IBKSTP,-1)
    ISTAT = MRBADD (KBUF1,IBLKNO,INELE,INVAL,INT,IBFAM, &
                   IBDESC,IBTYP,INBIT,IBIT0,IDATYP, &
                   KLISTE,KTBLVAL)

    RETURN
  END


  SUBROUTINE mwbg_readTovs(IUNENT,HANDLE,ISAT,ZMISG,BUF1,TBLVAL, &
       LSTELE,ELDALT,IDATA,ICANO,ICANOMP,ICANOMPNA,INO, &
       INOMP,INOMPNA,DONIALT,ZDATA,ZO,ZCOR,ZOMP,STNID, &
       ZOMPNA,IMARQ,MXELM,MXVAL,MXNT,NELE,NVAL,NT, &
       ISCNCNT,ISCNPOS,MXSCAN,ZLAT,ZLON,ITERMER, &
       IORBIT,SATZEN,ITERRAIN,SKIPENR)
!OBJET          Lire les donnees ATOVS.
!
!ARGUMENTS      iunent    - input  -  unite logique du fichier burp
!               handle    - in/out -  "handle" du fichier burp
!               isat      - output -  numero d'identificateur du satellite 
!               zmisg     - input  -  valeur manquante burp 
!               buf1      - input  -  tableau contenant le rapport burp
!               tblval    - input  -  champ de donnees (valeurs entieres)
!               lstele    - input  -  liste des noms d'elements
!               eldalt    - input  -  liste des noms d'elements (decodes)
!               idata     - input  -  vecteur de travail
!               icano     - output -  canaux des radiances
!               icanomp   - output -  canaux des residus
!               icanompna - output -  canaux des residus au nadir
!               ino       - output -  nombre de canaux (radiances)
!               inomp     - output -  nombre de canaux (residus)
!               inompna   - output -  nombre de canaux (residus au nadir)
!               donialt   - input  -  champ de donnees (valeurs reelles)
!               zdata     - input  -  vecteur de travail
!               zo        - output -  radiances
!               zcor      - output -  correction aux radiances
!               zomp      - output -  residus des radiances
!               stnid     - output -  etiquette du satellite
!               zompna    - output -  residus des radiances au nadir
!               imarq     - output -  marqueurs des radiances
!               mxelm     - input  -  nombre maximum d'elements
!               mxval     - input  -  nombre maximum de donnees par element
!               mxnt      - input  -  nombre maximum de groupes mxelm X mxval
!               nele      - output -  nombre d'elements
!               nval      - output -  nombre de donnees par element
!               nt        - output -  nombre de groupes NELE X NVAL
!               iscncnt   - output -  satellite location counter
!               iscnpos   - output -  position sur le "scan"
!               mxscan    - input  -  limite superieure
!               zlat      - output -  latitude
!               zlon      - output -  longitude
!               itermer   - output -  indicateur terre/mer
!               iorbit    - output -  numero d'orbite
!               satzen    - output -  angle zenith du satellite (deg.)
!               iterrain  - output -  indicateur du type de terrain
!               skipenr   - output -  sauter l'enregistrement?
    IMPLICIT NONE

    INTEGER MXELM, MXVAL, MXNT, MXSCAN

    INTEGER  IUNENT, HANDLE

    INTEGER MRFGET
    INTEGER MRBHDR,MRFLOC
    INTEGER TEMPS,DATE,OARS,RUNN
    INTEGER IDTYP,LATI,LONG,DX,DY,ELEV,DRND,NBLOCS,BLKNO
    INTEGER ISTAT,NVAL,NT,MISGINT
    INTEGER FLGS,SUP,XAUX
    INTEGER I,JJ,IPNTR,NELE,INO,INOMP,INOMPNA

    INTEGER BUF1     (:)
    INTEGER TBLVAL   (MXELM*MXVAL*MXNT)
    INTEGER LSTELE   (MXELM)
    INTEGER ELDALT   (MXELM)
    INTEGER IDATA    (MXVAL*MXNT)
    INTEGER ISCNCNT  (MXNT)
    INTEGER ISCNPOS  (MXNT)
    INTEGER ICANO    (MXVAL*MXNT)
    INTEGER ICANOMP  (MXVAL*MXNT)
    INTEGER ICANOMPNA(MXVAL*MXNT)
    INTEGER IMARQ    (MXVAL*MXNT)
    INTEGER ISAT     (MXNT)
    INTEGER ITERMER  (MXNT)
    INTEGER IORBIT   (MXNT)
    INTEGER ITERRAIN (MXNT)

    REAL ZMISG

    REAL  DONIALT (MXELM*MXVAL*MXNT)
    REAL  ZDATA   (MXVAL*MXNT)
    REAL  ZO      (MXVAL*MXNT)
    REAL  ZCOR    (MXVAL*MXNT)
    REAL  ZOMP    (MXVAL*MXNT)
    REAL  ZOMPNA  (MXVAL*MXNT)
    REAL  ZLAT    (MXNT)
    REAL  ZLON    (MXNT)
    REAL  SATZEN  (MXNT)

    CHARACTER*9 STNID 

    LOGICAL  SKIPENR

    DATA MISGINT  /   -1    /

    LOGICAL DEBUG
    COMMON /DBGCOM/ DEBUG

    SKIPENR = .FALSE.

    HANDLE = MRFLOC(IUNENT,HANDLE,'*********',-1,-1,-1,-1, &
                   -1,SUP,0)

    IF(HANDLE .GT. 0) THEN
       IF (DEBUG) THEN
          WRITE(6,*)'PROCESSING BRP FILE:HANDLE=',HANDLE
       ENDIF
       ISTAT = MRFGET(HANDLE,BUF1)

! obtenir les parametres descripteurs de l'enregistrement lu
       ISTAT = MRBHDR(BUF1,TEMPS,FLGS,STNID,IDTYP,LATI,LONG,DX,DY, &
                      ELEV,DRND,DATE,OARS,RUNN,NBLOCS,SUP,0,XAUX,0)

! sauter les enregistrements resume
       IF ( STNID(1:2) .EQ. '>>' ) THEN
          SKIPENR = .TRUE.
          RETURN
       ENDIF

! initialisation
       DO I = 1,MXELM*MXVAL*MXNT
          DONIALT(I) = ZMISG
       ENDDO
   
       DO  I = 1,MXNT
          SATZEN  (I) = ZMISG
          ITERRAIN(I) = MISGINT
       ENDDO

! 1) Bloc info 3d: bloc 5120

! extraire le bloc
       CALL XTRBLK (5120,-1,BUF1,LSTELE,TBLVAL,ELDALT,DONIALT, &
                   NELE,NVAL,NT,BLKNO)    
       IF(BLKNO .LE. 0) THEN
          WRITE(6,*)'3D INFO BLOCK NOT FOUND'
          CALL ABORT()
       ENDIF

! extraire la latitude; element 5002
       CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                    5002,IDATA,ZLAT,IPNTR)    
       IF(IPNTR .EQ. 0) THEN
          WRITE(6,*)'LATITUDE MISSING'
          CALL ABORT()
       ENDIF
       IF (DEBUG) THEN
           WRITE(6,*) ' LATITUDE = ', (ZLAT(JJ),JJ=1,NVAL*NT)
       ENDIF

! extraire la longitude; element 6002
       CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                    6002,IDATA,ZLON,IPNTR)    
       IF(IPNTR .EQ. 0) THEN
          WRITE(6,*)'LONGITUDE MISSING'
          CALL ABORT()
       ENDIF
       IF (DEBUG) THEN
           WRITE(6,*) ' LONGITUDE = ', (ZLON(JJ),JJ=1,NVAL*NT)
       ENDIF

! 2) Bloc info (general): bloc 3072

! extraire le bloc
       CALL XTRBLK (3072,-1,BUF1,LSTELE,TBLVAL,ELDALT,DONIALT, &
                   NELE,NVAL,NT,BLKNO)    
       IF(BLKNO .LE. 0) THEN
          WRITE(6,*)'GENERAL INFO BLOCK NOT FOUND'
          CALL ABORT()
       ENDIF

! extraire l'indicateur d'identification du satellite; element 0 01 007.
       CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                    1007,ISAT,ZDATA,IPNTR) 
       IF(IPNTR .EQ. 0) THEN
          WRITE(6,*)' SATELLITE IDENTIFIER MISSING'
          CALL ABORT()
       ENDIF
       IF (DEBUG) THEN
           WRITE(6,*) ' ISAT = ', (ISAT(JJ),JJ=1,NVAL*NT)
       ENDIF

! extraire l'indicateur terre/mer; element 8012.
       CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                    8012,ITERMER,ZDATA,IPNTR)    
       IF(IPNTR .EQ. 0) THEN
          WRITE(6,*)'LAND/SEA INDICATOR MISSING'
          CALL ABORT()
       ENDIF
       IF (DEBUG) THEN
           WRITE(6,*) ' ITERMER = ', (ITERMER(JJ),JJ=1,NVAL*NT)
       ENDIF

! extraire technique de traitement.

!       not done anymore, j.h. june 2001

! extraire le "beam position (field of view no.)"; element 5043.
       CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                    5043,ISCNPOS,ZDATA,IPNTR)
       IF(IPNTR .EQ. 0) THEN
          WRITE(6,*)'INFORMATION ON SCAN POSITION MISSING'
          CALL ABORT()
       ENDIF
       IF (DEBUG) THEN
           WRITE(6,*) ' ISCNPOS = ', (ISCNPOS(JJ),JJ=1,NVAL*NT)
       ENDIF

! extraire le "satellite zenith angle"; element 7024.
       CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                    7024,IDATA,SATZEN,IPNTR)
       IF(IPNTR .EQ. 0) THEN
          WRITE(6,*)'SATELLITE ZENITH ANGLE MISSING'
          CALL ABORT()
       ENDIF
       IF (DEBUG) THEN
          WRITE(6,*) ' SATZEN = ', (SATZEN(JJ),JJ=1,NVAL*NT)
       ENDIF

! le "terrain type": pour les donnees 1b; element 13039.
       CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                    13039,ITERRAIN,ZDATA,IPNTR)
       IF(IPNTR .EQ. 0 .AND. DEBUG) THEN
          WRITE(6,*)'WARNING: TERRAIN TYPE MISSING'
       ENDIF
       IF (DEBUG) THEN
           WRITE(6,*) ' ITERRAIN = ', (ITERRAIN(JJ),JJ=1,NVAL*NT)
       ENDIF

! extraire le numero d'orbite; element 05040.
       CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                    05040,IORBIT,ZDATA,IPNTR)
       IF(IPNTR .EQ. 0) THEN
          WRITE(6,*)'ORBIT NUMBER MISSING'
          CALL ABORT()
       ENDIF
       IF (DEBUG) THEN
           WRITE(6,*) ' IORBIT = ', (IORBIT(JJ),JJ=1,NVAL*NT)
       ENDIF

! 3) Bloc multi niveaux de radiances: bloc 9218, 9248, 9264.

! extraire le bloc
       CALL XTRBLK (9218,-1,BUF1,LSTELE,TBLVAL,ELDALT,DONIALT, &
                   NELE,NVAL,NT,BLKNO)     
       IF(BLKNO .LE. 0) THEN
          CALL XTRBLK (9248,-1,BUF1,LSTELE,TBLVAL,ELDALT,DONIALT, &
                      NELE,NVAL,NT,BLKNO)     
       ENDIF       
       IF(BLKNO .LE. 0) THEN
          CALL XTRBLK (9264,-1,BUF1,LSTELE,TBLVAL,ELDALT,DONIALT, &
                      NELE,NVAL,NT,BLKNO)     
       ENDIF  
       IF(BLKNO .LE. 0) THEN
          WRITE(6,*)'RADIANCE DATA BLOCK NOT FOUND'
          CALL ABORT()
       ENDIF

! extraire les canaux; element 2150 
       CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                    2150,ICANO,ZDATA,IPNTR)
       IF(IPNTR .EQ. 0) THEN
          WRITE(6,*)'CHANNELS FOR RADIANCE OBS. MISSING'
          CALL ABORT()
       ENDIF
       INO = NVAL
       IF (DEBUG) THEN
           WRITE(6,*) ' INO  = ',  INO
           WRITE(6,*) ' ICANO= ', (ICANO(JJ),JJ=1,NVAL*NT)
       ENDIF

! extraire les radiances.
       CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                    12163,IDATA,ZO,IPNTR) 
       IF(IPNTR .EQ. 0) THEN
          WRITE(6,*)'RADIANCE OBS. MISSING'
          CALL ABORT()
       ENDIF
       IF (DEBUG) THEN
           WRITE(6,*) ' ZO = ', (ZO(JJ),JJ=1,NVAL*NT)
       ENDIF

! extraire la correction aux radiances, element 012233.
       CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                    012233,IDATA,ZCOR,IPNTR)
       IF(IPNTR .EQ. 0) THEN
          IF (DEBUG) THEN
             WRITE(6,*)'RADIANCE CORRECTION MISSING'
          ENDIF
          DO I = 1, NVAL*NT
             ZCOR(I) = ZMISG
          ENDDO
       ENDIF
       IF (DEBUG) THEN
           WRITE(6,*) ' ZCOR = ', (ZCOR(JJ),JJ=1,NVAL*NT)
       ENDIF                                

! 4) Bloc marqueurs multi niveaux de radiances: bloc 15362, 15392, 15408.

!  extraire le bloc
      CALL XTRBLK (15362,-1,BUF1,LSTELE,TBLVAL,ELDALT,DONIALT, &
                  NELE,NVAL,NT,BLKNO)    
      IF(BLKNO .LE. 0) THEN
         CALL XTRBLK (15392,-1,BUF1,LSTELE,TBLVAL,ELDALT,DONIALT, &
                     NELE,NVAL,NT,BLKNO)    
      ENDIF        
      IF(BLKNO .LE. 0) THEN
         CALL XTRBLK (15408,-1,BUF1,LSTELE,TBLVAL,ELDALT,DONIALT, &
                     NELE,NVAL,NT,BLKNO)    
      ENDIF    
      IF(BLKNO .LE. 0) THEN
         WRITE(6,*)'RADIANCE DATA FLAG BLOCK NOT FOUND'
         CALL ABORT()
      ENDIF

! extraire les marqueurs de 13bits des radiances; element 212163 
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   212163,IMARQ,ZDATA,IPNTR)
      IF(IPNTR .EQ. 0) THEN
         WRITE(6,*)'RADIANCE DATA FLAGS MISSING'
         CALL ABORT()
      ENDIF
      IF (DEBUG) THEN
          WRITE(6,*) ' RADIANCE FLAGS = ', (IMARQ(JJ),JJ=1,NVAL*NT)
      ENDIF

! 5) Bloc multi niveaux de residus de radiances (O-P): bloc 9322, 9226, 9258, 9274, bfam 14

! extraire le bloc
       CALL XTRBLK (9322,14,BUF1,LSTELE,TBLVAL,ELDALT,DONIALT, &
                   NELE,NVAL,NT,BLKNO)  
       IF(BLKNO .LE. 0) THEN
          CALL XTRBLK (9226,14,BUF1,LSTELE,TBLVAL,ELDALT,DONIALT, &
                      NELE,NVAL,NT,BLKNO)  
       ENDIF 
       IF(BLKNO .LE. 0) THEN
          CALL XTRBLK (9258,14,BUF1,LSTELE,TBLVAL,ELDALT,DONIALT, &
                      NELE,NVAL,NT,BLKNO)  
       ENDIF 
       IF(BLKNO .LE. 0) THEN
          CALL XTRBLK (9274,14,BUF1,LSTELE,TBLVAL,ELDALT,DONIALT, &
                      NELE,NVAL,NT,BLKNO)  
       ENDIF
       IF(BLKNO .LE. 0) THEN
          WRITE(6,*)'WARNING: (O-P) RADIANCE BLOCK NOT FOUND'
          SKIPENR = .TRUE.
          RETURN
       ENDIF

! extraire les canaux; element 2150 
       CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                    2150,ICANOMP,ZDATA,IPNTR) 
       IF(IPNTR .EQ. 0) THEN
          WRITE(6,*)'CHANNELS FOR (O-P) RADIANCE MISSING'
          CALL ABORT()
       ENDIF
       INOMP = NVAL
       IF (DEBUG) THEN
           WRITE(6,*) ' INOMP   = ',  INOMP
           WRITE(6,*) ' ICANOMP = ', (ICANOMP(JJ),JJ=1,NVAL*NT)
       ENDIF

! extraire les residus (O-P) en radiances.
       CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                    12163,IDATA,ZOMP,IPNTR)
       IF(IPNTR .EQ. 0) THEN
          WRITE(6,*)'(O-P) RADIANCES MISSING'
          CALL ABORT()
       ENDIF
       IF (DEBUG) THEN
           WRITE(6,*) ' ZOMP = ', (ZOMP(JJ),JJ=1,NVAL*NT)
       ENDIF

! 6) Bloc multi niveaux de residus de radiances (O-P) au nadir: bloc 9322 ou 9226, bfam 32.

! not done any more, j.h. june 2001
    ENDIF

    RETURN
  END


  SUBROUTINE mwbg_ADDTRRN(KBUF1,KLISTE,KTBLVAL,KDLISTE,PRVAL, &
                     KDATA,PDATA,KTERMER,ITERRAIN, &
                     GLINTRP,KTBLVALN,KLISTEN,PRVALN)
!OBJET          Ajouter l'element "terrain type" au bloc info d'un
!               fichier burp.
!
!ARGUMENTS      kbuf1   - in/out -  tableau contenant le rapport
!               kliste  - input  -  liste des noms d'elements
!               ktblval - input  -  champ de donnees (valeurs entieres)
!               kdliste - input  -  liste des noms d'elements (decodes)
!               prval   - input  -  champ de donnees (valeurs reelles)
!               kdata   - input  -  donnees extraites (valeurs entieres)
!               pdata   - input  -  donnees extraites (valeurs reelles)
!               ktermer - input  -  indicateur terre/mer
!               iterrain- input  -  indicateur du type de terrain
!               glintrp - input  -  etendue de glace du modele
!               ktblvaln- in/out -  nouveau champ de donnees (valeurs entieres)
!               kdlisten- in/out -  nouvelle liste des noms d'elements (decodes)
!               prvaln  - in/out -  nouveau champ de donnees (valeurs reelles)
    IMPLICIT NONE

    INTEGER KBUF1    (:)
    INTEGER KDLISTE  (:)
    INTEGER KLISTE   (:)
    INTEGER KLISTEN  (:)
    INTEGER KTBLVAL  (:)
    INTEGER KTBLVALN (:)
    INTEGER KDATA    (:)
    INTEGER KTERMER  (:)
    INTEGER ITERRAIN (:)

    INTEGER INELE,INVAL,INT,JI,IPNTR,MRBREP
    INTEGER JJ,IBLKNO,ISTAT,INELEN
    INTEGER IDUM1,IDUM2,IDUM3,IBFAM
    INTEGER IBDESC,IBTYP,INBIT,IBIT0,IDATYP   
    INTEGER MRBPRM, MRBDEL, MRBADD

    REAL PRVAL   (:)
    REAL PRVALN  (:)
    REAL PDATA   (:)
    REAL GLINTRP (:)

    LOGICAL DEBUG
    COMMON /DBGCOM/ DEBUG

! 1) Bloc info (general): bloc 3072
!    Ajouter le "terrain type".

! extraire le bloc
    CALL XTRBLK (3072,-1,KBUF1,KLISTE,KTBLVAL,KDLISTE,PRVAL, &
                INELE,INVAL,INT,IBLKNO)    
    IF(IBLKNO .LE. 0) THEN
       WRITE(6,*)'INFO BLOCK NOT FOUND'
       CALL ABORT()
    ENDIF

!  Dans les conditions suivantes:
!      1) l'indicateur terre/mer indique l'ocean (ktermer=1),
!      2) le "terrain type" est manquant (iterrain=-1),
!      3) le modele indique de la glace (gl >= 0.01),
!  on specifie "sea ice" pour le "terrain type" (iterrain=0).
    IF (DEBUG) THEN
        WRITE(6,*) ' OLD TERRAIN TYPE = ', &
                  (ITERRAIN(JJ),JJ=1,INT)
        WRITE(6,*) ' KTERMER = ', &
                  (KTERMER(JJ),JJ=1,INT)
        WRITE(6,*) ' GLINTRP = ', &
                  (GLINTRP(JJ),JJ=1,INT)
    ENDIF
    DO JI = 1, INT
       IF ( KTERMER (JI) .EQ. 1    .AND. &
           ITERRAIN(JI) .EQ. -1   .AND. &
           GLINTRP (JI) .GE. 0.01       ) THEN
          ITERRAIN(JI) = 0
       ENDIF
    ENDDO
    IF (DEBUG) THEN
        WRITE(6,*) ' NEW TERRAIN TYPE = ', &
                  (ITERRAIN(JJ),JJ=1,INT)
    ENDIF

!  Le "terrain type" (element 13039) existe-t-il dans le bloc? 
!       si non, ajouter cet element, 
!       si oui, remplacer cet element.

    CALL XTRDATA (KDLISTE,KTBLVAL,PRVAL,INELE,INVAL,INT, &
                 13039,KDATA,PDATA,IPNTR)    
    IF(IPNTR .EQ. 0) THEN
       IF (DEBUG) THEN
           WRITE(6,*)'TERRAIN TYPE MISSING IN BLOCK. ADD IT'
           WRITE(6,*) ' KTBLVAL  = ', (KTBLVAL (JJ),JJ=1, &
                                       INELE*INVAL*INT)
       ENDIF
       CALL ADDDAT (KLISTE,KTBLVAL,PRVAL,INELE,INVAL,INT, &
                   13039,ITERRAIN,PDATA, &
                   KLISTEN,KTBLVALN,PRVALN,INELEN,0)
       IF (DEBUG) THEN
           WRITE(6,*) ' KTBLVALN = ', (KTBLVALN(JJ),JJ=1, &
                                        INELEN*INVAL*INT)
       ENDIF
       ISTAT = MRBPRM (KBUF1,IBLKNO,IDUM1,IDUM2,IDUM3,IBFAM, &
                   IBDESC,IBTYP,INBIT,IBIT0,IDATYP)
       ISTAT = MRBDEL(KBUF1,IBLKNO)
       ISTAT = MRBADD (KBUF1,IBLKNO,INELEN,INVAL,INT,IBFAM, &
                   IBDESC,IBTYP,INBIT,IBIT0,IDATYP, &
                   KLISTEN,KTBLVALN)
    ELSE
       CALL REPDATA (KDLISTE,KTBLVAL,INELE,INVAL,INT, &
                    13039,ITERRAIN,IPNTR)
       ISTAT = MRBREP (KBUF1,IBLKNO,KTBLVAL)
    ENDIF

! les indicateurs terre/mer corriges
    IF (DEBUG) THEN
        WRITE(6,*) ' NEW LAND/SEA = ', (KTERMER(JJ),JJ=1,INVAL*INT)
    ENDIF

! Remplacer le bloc.
!    ISTAT = MRBREP (KBUF1,IBLKNO,KTBLVAL)


    RETURN
  END


  SUBROUTINE ADDDAT (KLISTE,KTBLVAL,PRVAL,KNELE,KNVAL,KNT, &
                    KELEM,KDATA,PDATA,KLISTEN,KTBLVALN, &
                    PRVALN,KNELEN,KMODE)
!OBJET          Remplacer les donnees de l'element specifie 
!               (uni ou multi niveaux). Valeurs entieres seulement.
!
!APPEL          CALL ADDDAT   (KLISTE,KTBLVAL,PRVAL,KNELE,KNVAL,KNT,
!                              KELEM,KDATA,PDATA,KLISTEN,KTBLVALN,
!                              PRVALN,KNELEN,KMODE)
!
!ARGUMENTS      kliste  - input  -  liste des noms d'elements (non decodes)
!               ktblval - in/out -  champ de donnees (valeurs entieres)
!               prval   - in/out -  champ de donnees (valeurs reelles)
!               knele   - input  -  nombre d'elements
!               knval   - input  -  nombre de donnees par element
!               knt     - input  -  nombre de groupes KNELE X KNVAL
!               kelem   - input  -  element specifie 
!               kdata   - input  -  donnees extraites (valeurs entieres)
!               pdata   - input  -  donnees extraites (valeurs relles)
!               klisten - output -  nouvelle liste des noms d'elements (non decodes)
!               ktblvaln- output -  nouveau champ de donnees (valeurs entieres)
!               prvaln  - output -  nouveau champ de donnees (valeurs reelles)
!               knelen  - output -  nouveau nombre d'elements
!               kmode   - input  -  =0, mrbcvt pas appele (ajoute kdata),
!                                   =1, mrbcvt appele (ajoute pdata).
    IMPLICIT NONE

    INTEGER KLISTE   (:)
    INTEGER KTBLVAL  (:)
    INTEGER KLISTEN  (:)
    INTEGER KTBLVALN (:)
    INTEGER KDATA    (:)

    INTEGER KNELE,KNVAL,KNT,JI,MRBCVT,KNELEN
    INTEGER INDX,KELEM,JJ,IPOS,ISTAT,JK,INDXP,MRBCOL
    INTEGER KMODE

    REAL PDATA (:)
    REAL PRVAL (:)
    REAL PRVALN(:)

    LOGICAL DEBUG

    DATA  DEBUG /.FALSE./

    KNELEN = KNELE + 1
    DO JK = 1, KNELE 
       KLISTEN(JK) =  KLISTE(JK)
    ENDDO
    ISTAT = MRBCOL(KELEM,KLISTEN(KNELEN),1)

    IPOS = 0
    INDX = 0
    INDXP= 0
    DO JI = 1, KNT 
       DO JJ = 1, KNVAL 
          INDXP = INDXP + 1
          DO JK = 1, KNELE 
             INDX = INDX + 1
             IPOS = IPOS + 1
             IF (KMODE.EQ.0) THEN
                KTBLVALN(IPOS) =  KTBLVAL(INDX)
             ELSE
                PRVALN(IPOS) =  PRVAL(INDX)
             ENDIF 
          ENDDO
          IPOS = IPOS + 1
          IF (KMODE.EQ.0) THEN
             KTBLVALN(IPOS) =  KDATA(INDXP)
          ELSE
             PRVALN(IPOS) =  PDATA(INDXP) 
          ENDIF 
       ENDDO
    ENDDO
    IF (KMODE.NE.0) THEN
       ISTAT = MRBCVT(KLISTEN,KTBLVALN,PRVALN,KNELEN, &
                     KNVAL,KNT,1)
    ENDIF
    IF (DEBUG) THEN
       WRITE(6,*) ' ADDDAT: KTBLVALN =  ',(KTBLVALN(JJ), &
                 JJ=1,KNELEN*KNVAL*KNT)
    ENDIF

    RETURN
  END


  SUBROUTINE GRODY (ier, ni, tb23, tb31, tb50, tb53, tb89, pangl, &
                   plat, ilansea, ice, tpw, clw, rain, snow, scatl, &
                   scatw )
!OBJET          Compute the following parameters using 5 AMSU-A
!               channels:
!                  - sea ice, 
!                  - total precipitable water, 
!                  - cloud liquid water, 
!                  - ocean/land rain, 
!                  - snow cover/glacial ice,
!                  - scattering index (sur la terre et sur l'eau).
!               The four channels used are: 23Ghz, 31Ghz, 50Ghz and 89Ghz.
!
!REGERENCES     N. Grody, NOAA/NESDIS, ....
!
!APPEL          CALL   GRODY (ier, ni, tb23, tb31, tb50, tb53, tb89, pangl, plat,
!                             ilansea, ice, tpw, clw, rain, snow, scatl, scatw) 
!
!ARGUMENTS      ier     - output - error return code:
!                                  0, ok,  
!                                  1, input parameter out of range. 
!               ni      - input  -  number of points to process
!               tb23    - input  -  23Ghz brightness temperature (K)
!               tb31    - input  -  31Ghz brightness temperature (K)
!               tb50    - input  -  50Ghz brightness temperature (K)
!               tb53    - input  -  53Ghz brightness temperature (K)
!               tb89    - input  -  89Ghz brightness temperature (K)
!               pangl   - input  -  satellite zenith angle (deg.)
!               plat    - input  -  lalitude (deg.)
!               ilansea - input  -  land/sea indicator (0=land;1=ocean)
!               ice     - output -  sea ice concentration (0-100%)
!               tpw     - output -  total precipitable water (0-70mm)
!               clw     - output -  cloud liquid water (0-3mm)
!               rain    - output -  rain identification (0=no rain; 1=rain)
!               snow    - output -  snow cover and glacial ice identification: 
!                                   (0=no snow; 1=snow; 2=glacial ice)
!               scatl   - output -  scattering index over land
!               scatw   - output -  scattering index over water
!
!! Notes: In the case where an output parameter cannot be calculated, the
!!        value of this parameter is to to the missing value, i.e. -99.

    IMPLICIT NONE

    integer ni, i

    integer ier    (:)
    integer ilansea(:)
    integer rain   (:)
    integer snow   (:)

    real zmisg, siw, sil, df1, df2, df3, a, b, c, d, e23
    real ei, cosz, tt, scat, sc31, abslat, t23, t31, t50, t89
    real sc50, par, t53
    real dif285t23, dif285t31, epsilon

    real tb23  (:)
    real tb31  (:)
    real tb50  (:)
    real tb53  (:)
    real tb89  (:)
    real pangl (:)
    real plat  (:)
    real ice   (:)
    real tpw   (:)
    real clw   (:)
    real scatl (:)
    real scatw (:)

    data zmisg   / -99.     /
    data epsilon /   1.E-30 /

! 1) Initialise output parameters:
    do i = 1, ni
       ice  (i) = zmisg
       tpw  (i) = zmisg
       clw  (i) = zmisg
       scatl(i) = zmisg
       scatw(i) = zmisg
       rain (i) = nint(zmisg)
       snow (i) = nint(zmisg)
    enddo

! 2) Validate input parameters:
    do i = 1, ni

       if ( tb23(i)    .lt. 120.  .or. &
           tb23(i)    .gt. 350.  .or. &
           tb31(i)    .lt. 120.  .or. &
           tb31(i)    .gt. 350.  .or. &
           tb50(i)    .lt. 120.  .or. &
           tb50(i)    .gt. 350.  .or. &
           tb53(i)    .lt. 120.  .or. &
           tb53(i)    .gt. 350.  .or. &
           tb89(i)    .lt. 120.  .or. &
           tb89(i)    .gt. 350.  .or. &
           pangl(i)   .lt. -90.  .or. &
           pangl(i)   .gt.  90.  .or. &
           plat(i)    .lt. -90.  .or. &
           plat(i)    .gt.  90.  .or. &
           ilansea(i) .lt.   0   .or. &
           ilansea(i) .gt.   1        ) then
          ier(i) = 1
       else
          ier(i) = 0
       endif

    enddo

!_3) Compute parameters:
    do i = 1, ni

      if ( ier(i) .eq. 0 ) then

         abslat = abs(plat(i))
         cosz   = cosd(pangl(i))
         t23 = tb23(i)
         t31 = tb31(i)
         t50 = tb50(i)
         t53 = tb53(i)
         t89 = tb89(i)
         dif285t23=max(285.-t23,epsilon)
         dif285t31=max(285.-t31,epsilon)

! scattering indices:
         siw = -113.2 + (2.41 - 0.0049*t23)*t23 &
                     + 0.454*t31 &
                     - t89 
         sil = t23 - t89 

         scatl (i) = sil
         scatw (i) = siw

! discriminate functions:
         df1 =  2.85 + 0.020*t23 - 0.028*t50 ! used to identify (also remove) sea ice
         df2 =  5.10 + 0.078*t23 - 0.096*t50 ! used to identify (also remove) warm deserts
         df3 = 10.20 + 0.036*t23 - 0.074*t50 ! used to identify (also remove) cold deserts

         if ( ilansea(i) .eq. 1 ) then

! Ocean Parameters

!_3.1) Sea ice:
            if ( abslat .lt. 50. ) then
               ice(i) = 0.0
            else
               if ( df1 .lt. 0.45 ) then
                  ice(i) = 0.0
               else
                  a =  1.7340  - 0.6236*cosz
                  b =  0.0070  + 0.0025*cosz
                  c = -0.00106 
                  d = -0.00909
                  e23 = a + b*t31 + c*t23 + d*t50 ! theoretical 23Ghz sfc emissivity (0.3-1.)
                  if ( (t23-t31) .ge. 5. ) then   ! fov contains multiyear or new ice/water
                     ei = 0.88
                  else
                     ei = 0.95
                  endif
                  ice(i) = 100*(e23-0.45)/(ei-0.45) ! sea-ice concentration within fov (0-100%) 
                  ice(i) = min(100.,max(0.,ice(i)))/100.   !jh (0.-1.)
               endif
            endif

! 3.2) Total precipitable water:
            ! identify and remove sea ice
            if ( abslat .gt. 50.  .and. &
                df1    .gt.  0.2        ) then  
               tpw(i) = zmisg
            else
               a =  247.920  - (69.235 - 44.177*cosz)*cosz
               b = -116.270
               c =   73.409
               tpw(i) = a + b*log(dif285t23) & 
                         + c*log(dif285t31)
               tpw(i) = tpw(i)*cosz           ! theoretical total precipitable water (0-70mm)
               tpw(i) = 0.942*tpw(i) - 2.17   ! corrected   total precipitable water 
               tpw(i) = min(70.,max(0.,tpw(i)))   ! jh     
            endif

!_3.3) Cloud liquid water:
            ! identify and remove sea ice
            if ( abslat .gt. 50.  .and. &
                df1    .gt.  0.0        ) then  
               clw(i) = zmisg
            else
               a =  8.240 - (2.622 - 1.846*cosz)*cosz
               b =  0.754
               c = -2.265
               clw(i) = a + b*log(dif285t23) & 
                         + c*log(dif285t31)
               clw(i) = clw(i)*cosz           ! theoretical cloud liquid water (0-3mm)
               clw(i) = clw(i) - 0.03         ! corrected   cloud liquid water 
               clw(i) = min(3.,max(0.,clw(i)))   ! jh       
            endif

!_3.4) Ocean rain: 0=no rain; 1=rain.
            ! identify and remove sea ice
            if ( abslat .gt. 50.  .and. &
                df1    .gt.  0.0        ) then  
               rain(i) = nint(zmisg)
            else                                   ! remove non-precipitating clouds
               if ( clw(i) .gt. 0.3 .or. &
                   siw    .gt. 9.0      ) then 
                  rain(i) = 1
               else
                  rain(i) = 0
               endif
            endif

         else

! Land Parameters

! 3.5) Rain  over land: 0=no rain; 1=rain.
            tt = 168. + 0.49*t89
            if ( sil .ge. 3. ) then
               rain(i) = 1
            else 
               rain(i) = 0
            endif
            
            ! remove snow cover
            if ( t23 .le. 261. .and. &
                t23 .lt. tt         ) then
               rain(i) = 0
            endif

            ! remove warm deserts
            if ( t89 .gt. 273.  .or. &
                df2 .lt.   0.6      ) then
               rain(i) = 0
            endif

! 3.6) Snow cover and glacial ice: 0=no snow; 1=snow; 2=glacial ice.
            tt = 168. + 0.49*t89
            scat = t23 - t89
            sc31 = t23 - t31
            sc50 = t31 - t50
            par  = t50 - t53

            ! re-frozen snow
            if ( t89  .lt. 255.  .and. &
                scat .lt. sc31        ) then
               scat = sc31
            endif

            ! identify glacial ice
            if ( scat .lt.   3.  .and. &
                t23  .lt. 215.        ) then
               snow(i) = 2
            endif
            if ( scat .ge. 3.       ) then
               snow(i) = 1
            else
               snow(i) = 0
            endif

            ! remove precipitation
            if ( t23 .ge. 262.  .or. &
                t23 .ge. tt           ) then
               snow(i) = 0
            endif
            if ( df3 .le. 0.35         ) then    ! remove deserts
               snow(i) = 0
            endif

            ! high elevation deserts
            if ( scat .lt.  15.  .and. &
                sc31 .lt.   3.  .and. &
                par  .gt.   2.        ) then
               snow(i) = 0
            endif

            ! remove frozen ground
            if ( scat .lt.   9.  .and. &
                sc31 .lt.   3.  .and. &
                sc50 .lt.   0.        ) then
               snow(i) = 0
            endif

         endif

      endif
      if ( .false. ) then
      print *, ' '
      print *, ' i,tb23(i),tb31(i),tb50(i),tb89(i),pangl(i),plat(i), &
                ilansea(i) = ', &
                i,tb23(i),tb31(i),tb50(i),tb89(i),pangl(i),plat(i), &
                ilansea(i)
      print *, ' ier(i),ice(i),tpw(i),clw(i),rain(i),snow(i)=', &
                ier(i),ice(i),tpw(i),clw(i),rain(i),snow(i)
      if ( i .eq. 100 )stop
      endif

    enddo

    return
  end


  SUBROUTINE mwbg_readStatTovs(ILUTOV,INUMSAT,CSATID)
!OBJET          Lire les statistiques de l'erreur totale pour les TOVS.
!
!ARGUMENTS      ilutov  - input  -  unite logique du fichier stats des TOVS
!               inumsat - output -  nombre de satellites
!               csatid  - output -  identificateur de satellite 
    IMPLICIT NONE

    INTEGER JPNSAT,JPCH
 
    PARAMETER (JPNSAT =  9) 
    PARAMETER (JPCH = 50)

    INTEGER JPMXSFC
    PARAMETER (JPMXSFC =  2)

    INTEGER NCHNA    (JPNSAT)
    INTEGER MLISCHNA (JPCH,JPNSAT)
    REAL    TOVERRST (JPCH,JPNSAT)
    INTEGER IUTILST  (JPCH,JPNSAT)
    COMMON /COMTOVST/ NCHNA, MLISCHNA, TOVERRST, IUTILST

    INTEGER ILUTOV, JI, JJ, JK, JL, JM, I, ICHN, NULOUT
    INTEGER INUMSAT, INDX, IPOS

    INTEGER NUMCHNIN(JPNSAT), ISATID(JPNSAT)

    REAL*8  TOVERRIN(JPCH,2,JPNSAT)
    REAL*8  ZDUM

    CHARACTER*132  CLDUM
    CHARACTER*17   CSATSTR

    CHARACTER*9   CSATID(JPNSAT)
    CHARACTER*12  CTYPSTAT(2)

    DATA NULOUT /  6 /

    DATA CTYPSTAT     / 'Monitoring',  'Assimilation'  /  

    WRITE(NULOUT,FMT=9000)
 9000 FORMAT(//,10x,"-mwbg_readStatTovs: reading total error statistics" &
          ," required for TOVS processing")
!
! 1. Initialize
 100  CONTINUE
    DO JL = 1, JPNSAT
        NCHNA(JL) = 0
        NUMCHNIN(JL) = 0
        ISATID(JL) = 0
        DO JI = 1, JPCH
          TOVERRIN(JI,1,JL) = 0.0
          TOVERRIN(JI,2,JL) = 0.0
          IUTILST (JI,JL) = 0
          TOVERRST(JI,JL) = 0.0
        ENDDO
    ENDDO

! 2. Open the file

 200  CONTINUE
! .... not done here anymore, jh, august 2000

! 3. Print the file contents

 300  CONTINUE

    WRITE(NULOUT,'(20X,"ASCII dump of stats_tovs file: "//)')
    DO JI = 1, 9999999
       READ (ILUTOV,'(A)',ERR=900,END=400) CLDUM
       WRITE(NULOUT,'(A)')   CLDUM
    ENDDO

! 4. Read number of satellites

 400  CONTINUE

    REWIND(ILUTOV)
    READ (ILUTOV,*,ERR=900)
    READ (ILUTOV,*,ERR=900) INUMSAT
    READ (ILUTOV,*,ERR=900)

! 5. Read the satellite identification, the number of channels,
! the observation errors and the utilization flags

 500  CONTINUE

    DO JL = 1, INUMSAT
       READ (ILUTOV,*,ERR=900)
       READ (ILUTOV,'(A)',ERR=900) CLDUM
       CSATSTR = TRIM(ADJUSTL(CLDUM))

! Get satellite (e.g. NOAA18) from satellite/instrument (e.g. NOAA18 AMSUA)
       INDX = INDEX(CSATSTR,'AMSUA')
       IF ( INDX .GT. 3 ) THEN
         IF ( INDEX(CSATSTR,'EOS-2') .GT. 0 ) THEN
           CSATID(JL) = 'AQUA'
         ELSE
           IPOS = INDX-2
           CSATID(JL) = CSATSTR(1:IPOS)
         ENDIF
       ELSE
         WRITE ( NULOUT, '(" mwbg_readStatTovs: Non-AMSUA ", &
                    "instrument found in stats file!")' )
         WRITE ( NULOUT,'(A)' ) CLDUM
         CALL ABORT ()
       ENDIF
       READ (ILUTOV,*,ERR=900)
       READ (ILUTOV,*,ERR=900) ISATID(JL), NUMCHNIN(JL)
       DO JI = 1, 3
          READ (ILUTOV,*,ERR=900)
       ENDDO
! Set errors to ERRBGCK column values
       DO JI = 1, NUMCHNIN(JL)
          READ (ILUTOV,*,ERR=900) ICHN, &
                   TOVERRIN(ICHN,1,JL), &
                   TOVERRIN(ICHN,2,JL), &
                   IUTILST (ICHN,JL), ZDUM
          TOVERRST(ICHN,JL) = TOVERRIN(ICHN,1,JL)
       ENDDO
       READ (ILUTOV,*,ERR=900)
    ENDDO

 510  CONTINUE

! 6. Print error stats for assimilated channels

 600  CONTINUE

    WRITE(NULOUT,'(//5X,"Total errors for TOVS data"/)') 
    DO JL = 1, INUMSAT
       INDX = 0
       DO JI = 1, JPCH
          IF ( IUTILST(JI,JL) .NE. 0 ) THEN
            NCHNA(JL) = NCHNA(JL) + 1
            INDX = INDX + 1
            MLISCHNA(INDX,JL) = JI
          ENDIF
       ENDDO
       DO JK = 1, 2
          WRITE(NULOUT,'(/7X,"Satellite: ",A,5X,A)') &
            CSATID(JL), CTYPSTAT(JK)
          WRITE(NULOUT,'(7X,"Channels   : ",30(T22,27I4/))') & 
           (MLISCHNA(JI,JL),JI=1,NCHNA(JL))
          WRITE(NULOUT,'(7X,"Total errors: ",30(T22,27f4.1/))') &
           (TOVERRIN(MLISCHNA(JI,JL),JK,JL), &
            JI=1,NCHNA(JL))
       ENDDO
    ENDDO
     
! 7. Close the file
     
 700  CONTINUE

! .... not done here anymore, jh, august 2000

    RETURN

! Read error

900   WRITE ( NULOUT, '(" mwbg_readStatTovs: Problem reading ", &
                     "TOVS total error stats file ")' ) 
    CALL ABORT ()

    RETURN
  END
      

end module bgckmicrowave_mod

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
  use burp_module

  implicit none
  save
  private

  ! public variables
  public :: mwbg_debug, mwbg_clwThreshold

  ! Public functions
  public :: mwbg_readStatTovs, mwbg_readStatTovsAtms, mwbg_readTovs, mwbg_readTovsAtms, mwbg_tovCheckAmsua, mwbg_tovCheckAtms, mwbg_qcStatsAmsua, mwbg_qcStatsAtms, mwbg_UPDATFLG, mwbg_updateBurpAmsua, mwbg_updatFlgAtms, mwbg_ADDTRRN, mwbg_ADDTRRNF90, mwbg_getDataAmsua, mwbg_writeBlocksAmsua

  logical :: mwbg_debug, mwbg_clwThreshold

  integer, parameter :: JPNSAT = 9
  integer, parameter :: JPCH = 50
  integer, parameter :: MXCHN = 42 
  integer, parameter :: JPMXREJ = 15
  integer, parameter :: MXSAT = 9
  integer, PARAMETER :: MXVAL = 22
  integer, PARAMETER :: MXNT = 3000
  integer, parameter :: nchanAtms=22
  integer, parameter :: mxscan=96
  real, parameter    :: zmisg=9.9e09
  integer, parameter :: MISGINT = -1

  INTEGER :: NCHNA(JPNSAT), MLISCHNA(JPCH,JPNSAT), IUTILST(JPCH,JPNSAT)
  REAL    :: TOVERRST(JPCH,JPNSAT)

  INTEGER :: MREJCOD(JPMXREJ,MXCHN,MXSAT), INTOT(MXSAT), INTOTRJF(MXSAT), INTOTRJP(MXSAT)

  ! ATMS QC programs
  ! public variables
  public :: mwbg_modlsqtt, mwbg_useUnbiasedObsForClw 

  ! Public functions
  public :: mwbg_getData, mwbg_landIceMaskAtms, mwbg_grossValueCheck
  public :: mwbg_firstQcCheckAtms, mwbg_nrlFilterAtms, mwbg_writeBlocks

  integer  :: error

  logical :: mwbg_modlsqtt, mwbg_useUnbiasedObsForClw 

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

    KBLKNO = 0
    KBLKNO =  MRBLOC(KBUF,KBFAM,-1,KBTYP,KBLKNO)
    IF(KBLKNO .GT. 0) THEN
      IF ( mwbg_debug ) THEN
         WRITE(*,*)'FOUND BLOCK WITH BTYP = ',KBTYP
      ENDIF

      ! extraction du bloc
      ISTAT = MRBXTR(KBUF,KBLKNO,KLISTE,KTBLVAL)

      ! extraction des parametres descripteurs du bloc
      ISTAT = MRBPRM (KBUF,KBLKNO,KNELE,KNVAL,KNT,IBFAM, &
                     IBDESC,IBTYP,INBIT,IBIT0,IDATYP)

      ! extrait bktyp et btyp
      ISTAT = MRBTYP(IBKNAT,IBKTYP,IBKSTP,IBTYP)
      IF ( mwbg_debug ) THEN
         WRITE(*,*)'IBTYP=',IBTYP
         WRITE(*,*)'IBKNAT,IBKTYP,IBKSTP=',IBKNAT,IBKTYP,IBKSTP
         WRITE(*,*)'KNELE,KNVAL,KNT,INBIT,IBDESC,IBFAM,IDATYP=' &
                  ,KNELE,KNVAL,KNT,INBIT,IBDESC,IBFAM,IDATYP
      ENDIF

      ! extraction des informations du bloc et de la liste
      ISTAT = MRBXTR(KBUF,KBLKNO,KLISTE,KTBLVAL)
      ISTAT = MRBCVT(KLISTE,KTBLVAL,PRVAL,KNELE, &
                    KNVAL,KNT,0)
      ISTAT = MRBDCL(KLISTE,KDLISTE,KNELE)

      IF ( mwbg_debug ) THEN
         WRITE(*,*) 'PRVAL:'
         WRITE(*,*) (PRVAL(JN),JN=1,KNELE*KNVAL*KNT)
         WRITE(*,*) 'KDLISTE:'
         WRITE(*,*) (KDLISTE(JN), JN=1,KNELE)
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

    KBLKNO = 0
    KBLKNO =  MRBLOC(KBUF,KBFAM,-1,KBTYP,KBLKNO)
    IF(KBLKNO .GT. 0) THEN
      IF ( mwbg_debug ) THEN
        WRITE(*,*)'FOUND BLOCK WITH BTYP = ',KBTYP
      ENDIF

      ! enlever le bloc
      ISTAT = MRBDEL(KBUF,KBLKNO)
      IF(ISTAT .NE. 0) THEN
        WRITE(*,*)' PROBLEME AVEC RMBLK - BLOC ',KBLKNO,' ENLEVE'
      ENDIF
      IF ( mwbg_debug ) THEN
        WRITE(*,*)'MRBDEL ISTAT = ',ISTAT
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
      IF ( mwbg_debug ) THEN
        WRITE(*,*) ' XTRDATA: No data for element ',KELEM
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
      IF ( mwbg_debug ) THEN
        WRITE(*,*) ' XTRDATA: kdata =  ',(KDATA(JJ), &
                  JJ=1,KNVAL*KNT)
        WRITE(*,*) ' XTRDATA: pdata =  ',(PDATA(JJ), &
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

    KPNTR = 0
    DO JI = 1, KNELE
      IF ( KDLISTE(JI) .EQ. KELEM ) KPNTR = JI
    ENDDO

    ! any missing elements?
    IF ( KPNTR .EQ. 0  ) THEN
      IF ( mwbg_debug ) THEN
        WRITE(*,*) ' XTRDATA: No data for element ',KELEM
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
      IF ( mwbg_debug ) THEN
        WRITE(*,*) ' REPDATA: ktblval =  ',(KTBLVAL(JJ), &
                   JJ=1,KNELE*KNVAL*KNT)
      ENDIF
    ENDIF

    RETURN
  END


  SUBROUTINE mwbg_tovCheckAmsua(KSAT, KTERMER, KORBIT, ICANO, ICANOMP, ZO, ZCOR, &
                                ZOMP, ICHECK, KNO, KNT, PMISG, KNOSAT, KCHKPRF, &
                                ISCNPOS, MGINTRP, MTINTRP, GLINTRP, ITERRAIN, SATZEN, &
                                IMARQ, clw, scatw, STNID, RESETQC, ZLAT)
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
    !               clw     - output -  retrieved cloud liquid water
    !               scatw   - output -  scattering index over water
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

    INTEGER MXCLWREJ, MXCANPRED, MXSFCREJ2
    INTEGER MXSCANHIRS, MXSCANAMSU, MXSCATREJ, MXSFCREJ, NTESTS
    PARAMETER  ( MXSCANHIRS= 56 )
    PARAMETER  ( MXSCANAMSU= 30 )
    PARAMETER  ( MXCLWREJ  =  6 )
    PARAMETER  ( MXSFCREJ  =  6 )
    PARAMETER  ( MXSFCREJ2 =  4 )
    PARAMETER  ( MXSCATREJ =  7 )
    PARAMETER  ( MXCANPRED =  9 )

    INTEGER JPMXSFC
    PARAMETER (JPMXSFC =  2)
    
    INTEGER KNO,KNT,KNOSAT,MAXVAL
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
    INTEGER KCANOMP (KNO,KNT)
    INTEGER ICHECK  (KNO  ,KNT)
    INTEGER KCHKPRF (KNT)
    INTEGER IMARQ    (MXVAL*MXNT)
    INTEGER KMARQ   (KNO  ,KNT)
    INTEGER ITERRAIN(KNT)
    INTEGER ICLWREJ (MXCLWREJ)
    INTEGER ISFCREJ (MXSFCREJ)
    INTEGER ISFCREJ2(MXSFCREJ2)
    INTEGER ISCATREJ(MXSCATREJ)

    REAL  PMISG,EPSILON,ZANGL,MISGRODY,ZSEUILSCAT
    REAL  APPROXIM, ANGDIF, XCHECKVAL
    REAL  ZO      (MXVAL*MXNT)
    REAL  PTBO    (KNO    ,KNT)
    REAL  ZCOR    (MXVAL*MXNT)
    REAL  PTBCOR  (KNO    ,KNT)
    REAL  ZOMP    (MXVAL*MXNT)
    REAL  PTBOMP  (KNO  ,KNT)
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

    SAVE LLFIRST

    DATA  LLFIRST / .TRUE. /
    DATA  EPSILON / 0.01   /
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

        if ( ICANO(INDX) /= ICANOMP(INDX) ) then
          WRITE(*,*)'ERROR IN DIMENSIONS OF TOVS DATA'
          CALL ABORT()
        end if

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

    DO JJ=1,KNT
      DO JI=1,KNO
        IF ( KCANO(JI,JJ) .NE. KCANOMP(JI,JJ) ) THEN
          WRITE(*,*)'INCONSISTENT CHANNEL LISTS FOR TOVS DATA'
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

    ! 10) test 10: RTTOV reject check (single)
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
              IF ( mwbg_DEBUG ) THEN
                WRITE(*,*)STNID(2:9),' RTTOV REJECT.', &
                          'CHANNEL=', KCANO(JI,JJ), &
                          ' IMARQ= ',KMARQ(JI,JJ)
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    ! 1) test 1: Topography check (partial)
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
            IF ( mwbg_DEBUG ) THEN
              WRITE(*,*)STNID(2:9),' TOPOGRAPHY REJECT.', &
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
            IF ( mwbg_DEBUG ) THEN
              WRITE(*,*)STNID(2:9),' TOPOGRAPHY REJECT.', &
                        'CHANNEL=', KCANO(JI,JJ), &
                        ' TOPO= ',MTINTRP(JJ)
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    ! 2) test 2: "Land/sea qualifier" code check (full)
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
        IF ( mwbg_DEBUG ) THEN
          WRITE(*,*) STNID(2:9),'LAND/SEA QUALIFIER CODE', &
                   ' REJECT. KTERMER=', KTERMER(JJ)
        ENDIF
      ENDIF
    ENDDO

    ! 3) test 3: "Terrain type" code check (full)
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
          IF ( mwbg_debug ) THEN
            WRITE(*,*)STNID(2:9),'TERRAIN TYPE CODE', &
                     ' REJECT. TERRAIN=', ITERRAIN(JJ)
          ENDIF
        ENDIF
      ENDIF
    ENDDO

    ! 4) test 4: Field of view number check (full)
    !
    ! Field of view acceptable range is [1,MXSCANAMSU]  for AMSU footprints.
    INO = 4
    DO JJ=1,KNT
      DO JI=1,KNO
        IF ( ISCNPOS(JJ) .LT. 1 .OR. &
            ISCNPOS(JJ) .GT. MXSCANAMSU ) THEN
          ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
          KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
          KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
          MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
            MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
          IF ( mwbg_debug ) THEN
            WRITE(*,*)STNID(2:9),'FIELD OF VIEW NUMBER', &
                      ' REJECT. FIELD OF VIEW= ', ISCNPOS(JJ)
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    ! 5) test 5: Satellite zenith angle check (full)
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
          IF ( mwbg_debug ) THEN
            WRITE(*,*)STNID(2:9),' SATELLITE ZENITH ANGLE', &
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
          IF ( mwbg_debug ) THEN
            WRITE(*,*)STNID(2:9),' ANGLE/FIELD OF VIEW', &
                      ' INCONSISTENCY REJECT. SATZEN= ', &
                      SATZEN(JJ), ' FIELD OF VIEW= ',ISCNPOS(JJ), &
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
          IF ( mwbg_debug ) THEN
            WRITE(*,*)STNID(2:9),' LAND/SEA QUALIFIER', &
                      ' INCONSISTENCY REJECT. KTERMER= ', &
                      KTERMER(JJ), ' MODEL MASK= ',MGINTRP(JJ)
          ENDIF
        ENDIF
      ENDIF
    ENDDO

    ! 8) test 8: "Terrain type"/"Land/sea qual."/"model ice" consistency check. (full)
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
          IF ( mwbg_debug ) THEN
            WRITE(*,*)STNID(2:9),' TERRAIN TYPE/MODEL ICE', &
                      ' INCONSISTENCY REJECT. TERRAIN= ', &
                      ITERRAIN(JJ), ' MODEL ICE= ',GLINTRP(JJ)
          ENDIF
        ENDIF
        IF ( ITERRAIN(JJ) .EQ. 0       .AND. &
               KTERMER (JJ) .EQ. 0           ) THEN
          DO JI=1,KNO
            ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
            MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
               MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
          ENDDO
          IF ( mwbg_debug ) THEN
            WRITE(*,*)STNID(2:9),' TERRAIN TYPE/LAND?SEA QUAL.', &
                      ' INCONSISTENCY REJECT. TERRAIN= ', &
                      ITERRAIN(JJ), ' LAND/SEA= ',KTERMER(JJ)
          ENDIF
        ENDIF
        IF ( ITERRAIN(JJ) .EQ. 1       .AND. &
                KTERMER (JJ) .EQ. 1           ) THEN
          DO JI=1,KNO
            ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
            MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
               MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
          ENDDO
          IF ( mwbg_debug ) THEN
            WRITE(*,*)STNID(2:9),' TERRAIN TYPE/LAND?SEA QUAL.', &
                      ' INCONSISTENCY REJECT. TERRAIN= ', &
                      ITERRAIN(JJ), ' LAND/SEA= ',KTERMER(JJ)
          ENDIF
        ENDIF
      ELSE
        IF     ( KTERMER (JJ) .EQ. 1      .AND. &
                GLINTRP (JJ) .GT. 0.01         ) THEN
  !             DO JI=1,KNO
  !                ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
  !                MREJCOD(INO,KCANO(JI,JJ),KNOSAT) =
  !                  MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
  !             ENDDO
           IF ( mwbg_debug ) THEN
              WRITE(*,*)STNID(2:9),' TERRAIN TYPE MSG/MODEL ICE', &
                       ' INCONSISTENCY REJECT. TERRAIN= ', &
                       ITERRAIN(JJ), &
                       ' LAND/SEA= ',KTERMER(JJ)
           ENDIF
        ENDIF
      ENDIF
    ENDDO
800  continue

    ! 9) test 9: Uncorrected Tb check (single)
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
              IF ( mwbg_debug ) THEN
                WRITE(*,*)STNID(2:9),' UNCORRECTED TB REJECT.', &
                           'CHANNEL=', KCANO(JI,JJ), ' IMARQ= ',KMARQ(JI,JJ)
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    ! 11) test 11: Radiance observation "Gross" check (single) 
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
            IF ( mwbg_debug ) THEN
              WRITE(*,*)STNID(2:9),' GROSS CHECK REJECT.', &
                        'CHANNEL=', KCANO(JI,JJ), ' TB= ',PTBO(JI,JJ)
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    ! 12) test 12: Grody cloud liquid water check (partial)
    ! For Cloud Liquid Water > clwThreshold, reject AMSUA-A channels 1-5 and 15.
    INO = 12
    DO JJ=1,KNT
      IF ( CLW(JJ) .NE.  MISGRODY  ) THEN
        IF ( CLW(JJ) .GT. mwbg_clwThreshold   ) THEN
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
          IF ( mwbg_debug ) THEN
  !          IF (.true.) THEN
            WRITE(*,*)STNID(2:9),'Grody cloud liquid water check', &
                      ' REJECT. CLW= ',CLW(JJ), ' SEUIL= ',mwbg_clwThreshold
          ENDIF
        ENDIF
      ENDIF
    ENDDO

    ! 13) test 13: Grody scattering index check (partial)
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
          IF ( mwbg_debug ) THEN
  !          IF (.true.) THEN
            WRITE(*,*)STNID(2:9),'Grody scattering index check', &
                       ' REJECT. SCATW= ',SCATW(JJ), ' SEUIL= ',ZSEUILSCAT
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
            IF ( mwbg_debug ) THEN
              WRITE(*,*)STNID(2:9),'ROGUE CHECK REJECT.NO.', &
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
              IF ( mwbg_debug ) THEN
                 WRITE(*,*)STNID(2:9),'CHANNEL REJECT: ', &
                        ' OBS = ',JJ, &
                        ' CHANNEL= ',ICHN
              ENDIF
            ENDIF
          ENDIF
        ENDDO
    ENDDO

    IF ( mwbg_debug ) THEN
       WRITE(*,*)'ICHECK = ',((ICHECK(JI,JJ),JI=1,KNO),JJ=1,KNT)
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

    IF ( mwbg_debug ) THEN
      WRITE(*,*)'KCHKPRF = ',(KCHKPRF(JJ),JJ=1,KNT)
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

    CHARACTER*9   CSATID (MXSAT)

    INTEGER  JI, JJ, JK, KNO, KNT, KNOSAT
    INTEGER  INTOTOBS, INTOTACC, INUMSAT

    INTEGER  ICHECK (KNO,KNT)
    INTEGER  KCANO  (KNO,KNT)

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

        WRITE(*,'(/////50("*"))')
        WRITE(*,'(     50("*")/)')
        WRITE(*,'(T5,"SUMMARY OF QUALITY CONTROL FOR ", &
         A8)') CSATID(JK) 
        WRITE(*,'(T5,"------------------------------------- ",/)')
        WRITE(*,'( &
         "   TOTAL NUMBER OF AMSU-A  = ",I10,/ &
         " - TOTAL FULL REJECTS      = ",I10,/ &
         " - TOTAL PARTIAL REJECTS   = ",I10,/ &
         "   ------------------------------------",/ &
         "   TOTAL FULLY ACCEPTED    = ",I10,/)') &
          INTOTOBS, INTOTRJF(JK), INTOTRJP(JK), INTOTACC

        WRITE(*,'(//,1x,114("-"))')
        WRITE(*,'(t10,"|",t47,"REJECTION CATEGORIES")')
        WRITE(*,'(" CHANNEL",t10,"|",105("-"))')
        WRITE(*,'(t10,"|",15i7)') (JI,JI=1,JPMXREJ)
        WRITE(*,'(1x,"--------|",105("-"))')
        DO JJ = 1, MXCHN
           WRITE(*,'(3X,I2,t10,"|",15I7)') JJ,(MREJCOD(JI,JJ,JK), &
                                      JI=1,JPMXREJ)
        ENDDO
        WRITE(*,'(1x,114("-"))')
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

    ! 1) Bloc info 3d: bloc 5120.
    !    Modifier les marqueurs globaux de 24bits pour les donnees rejetees.

    !  extraire le bloc
    CALL XTRBLK (5120,-1,KBUF1,KLISTE,KTBLVAL,KDLISTE,PRVAL, &
                INELE,INVAL,INT,IBLKNO)    
    IF(IBLKNO .LE. 0) THEN
      WRITE(*,*)'3D INFO BLOCK NOT FOUND'
      CALL ABORT()
    ENDIF

    ! extraire les marqueurs globaux de 24bits; element 55200
    CALL XTRDATA (KDLISTE,KTBLVAL,PRVAL,INELE,INVAL,INT, &
                 55200,KDATA,PDATA,IPNTR)    
    IF(IPNTR .EQ. 0) THEN
      WRITE(*,*)'GLOBAL FLAGS MISSING'
      CALL ABORT()
    ENDIF
    IF ( mwbg_debug ) THEN
      WRITE(*,*) ' OLD FLAGS = ', (KDATA(JJ),JJ=1,INVAL*INT)
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
    IF ( mwbg_debug ) THEN
       WRITE(*,*) ' NEW FLAGS = ', (KDATA(JJ),JJ=1,INVAL*INT)
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
      WRITE(*,*)'INFO BLOCK NOT FOUND'
      CALL ABORT()
    ENDIF

    ! extraire l'indicateur terre/mer; element 8012
    CALL XTRDATA (KDLISTE,KTBLVAL,PRVAL,INELE,INVAL,INT, &
                 8012,KDATA,PDATA,IPNTR)    
    IF(IPNTR .EQ. 0) THEN
      WRITE(*,*)'LAND/SEA INDICATOR MISSING'
      CALL ABORT()
    ENDIF
    IF ( mwbg_debug ) THEN
      WRITE(*,*) ' OLD LAND/SEA = ', (KDATA(JJ),JJ=1,INVAL*INT)
    ENDIF

    ! les indicateurs terre/mer corriges
    IF ( mwbg_debug ) THEN
      WRITE(*,*) ' NEW LAND/SEA = ', (KTERMER(JJ),JJ=1,INVAL*INT)
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
      WRITE(*,*)'RADIANCE DATA BLOCK NOT FOUND'
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
      WRITE(*,*)'RADIANCE DATA FLAG BLOCK NOT FOUND'
      CALL ABORT()
    ENDIF

    ! extraire les marqueurs de 13bits des radiances; element 212163 (LEVEL 1B)
    IMELERAD =  212163 
    CALL XTRDATA (KDLISTE,KTBLVAL,PRVAL,INELE,INVAL,INT, &
                 IMELERAD,KDATA,PDATA,IPNTR) 
    IF(IPNTR .EQ. 0) THEN
      WRITE(*,*)'RADIANCE DATA FLAGS MISSING'
      CALL ABORT()
    ENDIF
    IF ( mwbg_debug ) THEN
       WRITE(*,*) ' OLD FLAGS = ', (KDATA(JJ),JJ=1,INVAL*INT)
    ENDIF

    ! update data flags
    DO JI = 1, INVAL*INT
      KDATA(JI) = IMARQ(JI)
    ENDDO
    IF ( mwbg_debug ) THEN
      WRITE(*,*) ' ICHECK = ', (ICHECK(JJ),JJ=1,INVAL*INT)
      WRITE(*,*) ' NEW FLAGS = ', (KDATA(JJ),JJ=1,INVAL*INT)
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
      WRITE(*,*)'WARNING: (O-P) RADIANCE BLOCK NOT FOUND'
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


  SUBROUTINE mwbg_updateBurpAmsua(clw, scatw, KCHKPRF, KTERMER, ITERRAIN, &
                                  GLINTRP, ICHECK, RESETQC, IMARQ, &
                                  rpt, rpt_out)
    !OBJET          Allumer les bits des marqueurs pour les tovs rejetes.
    !               Mettre a jour l'indicateur terre/mer qui a
    !               possiblement ete corrige pour la glace marine.
    !               Modifier le bktyp des donnees, marqueurs et (O-P) pourt
    !               signifier "vu par AO". 
    !
    !ARGUMENTS      kchprf  - input  -  indicateur global controle de qualite tovs. Code:
    !                                   =0, ok,
    !                                   >0, rejet,
    !               ktermer - input  -  indicateur terre/mer
    !               iterrain- input  -  indicateur du type de terrain
    !               glintrp - input  -  etendue de glace du modele
    !               icheck  - input  -  indicateur controle de qualite tovs au 
    !                                   niveau de chaque canal
    !               resetqc - input  -  reset the quality control flags before adding the new ones ? 
    !               imarq   - input  -  modified flag values from mwbg_tovCheckAmsua
    !               rpt     - input  -  tableau contenant le rapport
    !               rpt_out - output -  report to write
    IMPLICIT NONE

    real    :: clw(:)
    real    :: scatw(:)
    INTEGER :: KCHKPRF (:)
    INTEGER :: ITERRAIN(:)
    INTEGER :: KTERMER (:)
    REAL    :: GLINTRP (:)
    INTEGER :: ICHECK  (:)
    LOGICAL :: RESETQC
    INTEGER :: IMARQ   (:)
    type(BURP_RPT) :: rpt
    type(BURP_RPT) :: rpt_out

    type(BURP_BLOCK) :: blk, blk_copy 

    INTEGER :: KDATA(MXVAL*MXNT)

    integer :: error, ref_blk, my_nt,  my_nval, my_nele, my_idtyp
    integer :: my_btyp, my_bktyp, my_bfam, new_bktyp 
    integer :: indice, indice1, indice2, kk, jj, JI, j, ipos, idata  

    Call BURP_Init(blk, B2=blk_copy, IOSTAT=error)

    ! Read and modify the blocks in rpt and add them to rpt_out
    ref_blk = 0
    BLOCKS: do

      ref_blk = BURP_Find_Block(rpt,BLOCK= blk,SEARCH_FROM= ref_blk,IOSTAT= error)
      if (error /= burp_noerr) call abort()
      
      if (ref_blk < 0) Exit

      Call BURP_Get_Property(blk, &
                  NELE   = my_nele, &
                  NVAL   = my_nval, &       ! 1 or number of channels (obs per location) if Tb data/flag block
                  NT     = my_nt, &         ! 1 or number of locations in block
                  BTYP   = my_btyp, &
                  BKTYP  = my_bktyp, & 
                  BFAM   = my_bfam, &
                  IOSTAT = error)
      if (error /= burp_noerr) call abort()

      ! 1) Bloc info 3d: bloc 5120.
      !    Modifier les marqueurs globaux de 24bits pour les donnees rejetees.
      if (my_btyp == 5120) then     

        ! extraire les marqueurs globaux de 24bits; element 55200
        indice = BURP_Find_Element(blk,55200,IOSTAT=error)
        if ( indice > 0 ) then
          j = 1
          do kk = 1, my_nt
            kdata(kk) = BURP_Get_Tblval(blk,indice,j,kk,error)
          end do
        else
          write(*,*) 'GLOBAL FLAGS missing in 3D block (btyp=5120).'
          call abort()
        endif
        IF (mwbg_debug) THEN
          write(*,*) ' OLD FLAGS = ', (KDATA(JJ),JJ=1,my_nt)
        ENDIF

        ! allumer la bit (6) indiquant que l'observation a un element
        ! rejete par le controle de qualite de l'AO.
        !  N.B.: si on est en mode resetqc, on remet le marqueur global a
        !        sa valeur de defaut, soit 1024,  avant de faire la mise a jour.
        DO JI = 1, my_nt
          IF (RESETQC) THEN
            KDATA(JI) = 1024  
          ENDIF
          IF ( KCHKPRF(JI).NE.0  ) THEN
            KDATA(JI) = OR (KDATA(JI),2**6)
          ENDIF
        ENDDO
        IF (mwbg_debug) THEN
          write(*,*) ' KCHKPRF   = ', (KCHKPRF(JJ),JJ=1,my_nt)
          write(*,*) ' NEW FLAGS = ', (KDATA(JJ),JJ=1,my_nt)
        ENDIF

        ! Remplacer les nouveaux marqueurs dans le tableau.
        j = 1
        do kk = 1, my_nt
          idata = kdata(kk)
          Call BURP_Set_Tblval(blk,indice,j,kk,idata)
        end do

        blk_copy = blk
        Call BURP_Write_Block(rpt_out,BLOCK=blk_copy,CONVERT_BLOCK=.FALSE.,IOSTAT=error)
        if (error /= burp_noerr)  call abort()

      ! 2) Bloc info (general): bloc 3072
      !    Modifier les indicateurs terre/mer possiblement corriges pour la glace
      !    marine.
      !    Add new elements CLW and scatw
      else if (my_btyp == 3072) then

        ! Add new elements CLW and scatw
        Call BURP_Resize_Block(blk, ADD_NELE = 2, IOSTAT = error)
        if (error /= burp_noerr)  call abort()
        
        Call BURP_Set_Element(blk,NELE_IND=my_nele+1,ELEMENT=13209,IOSTAT=error)
        Call BURP_Set_Element(blk,NELE_IND=my_nele+2,ELEMENT=13208,IOSTAT=error)
        Call BURP_Encode_Block(blk)   ! encode the element numbers in the block

        j = 1
        do kk =1, my_nt
          Call BURP_Set_Rval(blk,NELE_IND=my_nele+1,NVAL_IND=j,NT_IND=kk,RVAL=clw(kk),IOSTAT=error)
          Call BURP_Set_Rval(blk,NELE_IND=my_nele+2,NVAL_IND=j,NT_IND=kk,RVAL=scatw(kk),IOSTAT=error)
        end do
        Call BURP_Convert_Block(blk)

        ! extraire l'indicateur terre/mer; element 8012
        indice = BURP_Find_Element(blk,8012,IOSTAT=error)
        if ( indice > 0 ) then
          j = 1
          do kk = 1, my_nt
            kdata(kk) = BURP_Get_Tblval(blk,indice,j,kk,error)
          end do
        else
          write(*,*) 'LAND/SEA INDICATOR MISSING in INFO block (btyp=3072).'
          call abort()
        endif
        IF (mwbg_debug) THEN
          write(*,*) ' OLD LAND/SEA = ', (KDATA(JJ),JJ=1,my_nt)
          WRITE(*,*) ' NEW LAND/SEA = ', (KTERMER(JJ),JJ=1,my_nt)
        ENDIF

        ! Remplacer les nouveaux indicateurs terre/mer dans le tableau.
        j = 1
        do kk = 1, my_nt
          idata = KTERMER(kk)
          Call BURP_Set_Tblval(blk,indice,j,kk,idata)
        end do

        !  Dans les conditions suivantes:
        !      1) l'indicateur terre/mer indique l'ocean (ktermer=1),
        !      2) le "terrain type" est manquant (iterrain=-1),
        !      3) le modele indique de la glace (gl >= 0.01),
        !  on specifie "sea ice" pour le "terrain type" (iterrain=0).
        IF ( mwbg_debug ) THEN
          WRITE(*,*) ' OLD TERRAIN TYPE = ', (ITERRAIN(JJ),JJ=1,my_nt)
          WRITE(*,*) ' KTERMER = ', (KTERMER(JJ),JJ=1,my_nt)
          WRITE(*,*) ' GLINTRP = ', (GLINTRP(JJ),JJ=1,my_nt)
        ENDIF
        DO JI = 1, my_nt
          IF ( KTERMER (JI) == 1 .and. ITERRAIN(JI) == -1 .and. GLINTRP (JI) >= 0.01 ) &
                ITERRAIN(JI) = 0
        ENDDO
        IF ( mwbg_debug ) THEN
          WRITE(*,*) ' NEW TERRAIN TYPE = ', (ITERRAIN(JJ),JJ=1,my_nt)
        ENDIF

        indice = BURP_Find_Element(blk,13039,IOSTAT=error)
        if ( indice <= 0 ) then
          write(*,*) 'TERRAIN TYPE MISSING in INFO block (btyp=3072).'
          call abort()
        endif

        ! Update terrain type dans le tableau.
        j = 1
        do kk = 1, my_nt
          idata = ITERRAIN(kk)
          Call BURP_Set_Tblval(blk,indice,j,kk,idata)
        end do

        blk_copy = blk
        Call BURP_Write_Block(rpt_out,BLOCK=blk_copy,CONVERT_BLOCK=.FALSE.,IOSTAT=error)
        if (error /= burp_noerr)  call abort()

      ! 3) Bloc multi niveaux de radiances: bloc 9218, 9248, 9264.
      !    Modifier le bktyp pour signifier "vu par AO".
      else if ( (my_btyp == 9218 .or. my_btyp == 9248 .or. my_btyp ==9264) .and. &
                my_bfam == 0 ) then 

        ! Remplacer le bloc et modifier le bktyp pour signifier "vu par AO".
        new_bktyp = my_bktyp + 4
        Call BURP_Set_Property(blk,BKTYP=new_bktyp)

        blk_copy = blk
        Call BURP_Write_Block(rpt_out,BLOCK=blk_copy,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
        if (error /= burp_noerr)  call abort()

      ! 4) Bloc marqueurs multi niveaux de radiances: bloc 15362, 15392, 15408.
      !    Modifier les marqueurs de 13bits associes a chaque radiance.
      !    Modifier le bktyp pour signifier "vu par AO".
      else if ( (my_btyp == 15362 .or. my_btyp == 15392 .or. my_btyp == 15408) .and. &
                my_bfam == 0 ) then 

        ! extraire les marqueurs de 13bits des radiances; element 212163 (LEVEL 1B)
        indice = BURP_Find_Element(blk,212163,IOSTAT=error)
        if ( indice > 0 ) then
          ipos = 0
          do kk = 1, my_nt
            do j = 1, my_nval
              ipos = ipos + 1
              KDATA(ipos) = BURP_Get_Tblval(blk,indice,j,kk,error)
            end do
          end do
        else
          write(*,*) 'ERREUR - Element 212163 missing in flag block.'
          call abort()
        endif 
        IF (mwbg_debug) THEN
          write(*,*) ' OLD FLAGS = ', (KDATA(JJ),JJ=1,my_nval*my_nt)
        ENDIF

        ! update data flags
        DO JI = 1, my_nval * my_nt
          KDATA(JI) = IMARQ(JI)
        ENDDO
        IF (mwbg_debug) THEN
          write(*,*) ' ICHECK = ', (ICHECK(JJ),JJ=1,my_nval*my_nt)
          write(*,*) ' NEW FLAGS = ', (KDATA(JJ),JJ=1,my_nval*my_nt)
        ENDIF

        ! Remplacer les nouveaux marqueurs dans le tableau.
        ipos = 0
        do kk =1, my_nt
          do j = 1, my_nval
            ipos = ipos + 1
            Call BURP_Set_Tblval(blk,indice,j,kk,KDATA(ipos))
          enddo
        enddo

        ! Remplacer le bloc et modifier le bktyp pour signifier "vu par AO".
        new_bktyp = my_bktyp + 4
        Call BURP_Set_Property(blk,BKTYP=new_bktyp,IOSTAT=error)
        if (error /= burp_noerr)  call abort()

        blk_copy = blk
        Call BURP_Write_Block(rpt_out,BLOCK=blk_copy,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
        if (error /= burp_noerr)  call abort()

      ! 5) Bloc multi niveaux de residus de radiances (O-P): bloc 9322, 9226, 9258, 9274, bfam 14
      !    Modifier le bktyp pour signifier "vu par AO".
      else if ( (my_btyp == 9322 .or. my_btyp == 9226 .or. my_btyp == 9258 .or. &
                my_btyp == 9274) .and. my_bfam == 14 ) then 

        ! Remplacer le bloc et modifier le bktyp pour signifier "vu par AO".
        new_bktyp = my_bktyp + 4
        Call BURP_Set_Property(blk,BKTYP=new_bktyp,IOSTAT=error)
        if (error /= burp_noerr)  call abort()

        blk_copy = blk
        Call BURP_Write_Block(rpt_out,BLOCK=blk_copy,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
        if (error /= burp_noerr)  call abort()

      ! OTHER BLOCK 
      else

        blk_copy = blk
        Call BURP_Write_Block(rpt_out,BLOCK=blk_copy,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
        if (error /= burp_noerr)  call abort()

      end if

    enddo BLOCKS

  end subroutine mwbg_updateBurpAmsua


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
    INTEGER ISTAT,NVAL,NT
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

    SKIPENR = .FALSE.

    HANDLE = MRFLOC(IUNENT,HANDLE,'*********',-1,-1,-1,-1, &
                   -1,SUP,0)

    IF(HANDLE .GT. 0) THEN
      IF ( mwbg_debug ) THEN
         WRITE(*,*)'PROCESSING BRP FILE:HANDLE=',HANDLE
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
        WRITE(*,*)'3D INFO BLOCK NOT FOUND'
        CALL ABORT()
      ENDIF

      ! extraire la latitude; element 5002
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   5002,IDATA,ZLAT,IPNTR)    
      IF(IPNTR .EQ. 0) THEN
        WRITE(*,*)'LATITUDE MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_debug ) THEN
        WRITE(*,*) ' LATITUDE = ', (ZLAT(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! extraire la longitude; element 6002
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   6002,IDATA,ZLON,IPNTR)    
      IF(IPNTR .EQ. 0) THEN
        WRITE(*,*)'LONGITUDE MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_debug ) THEN
        WRITE(*,*) ' LONGITUDE = ', (ZLON(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! 2) Bloc info (general): bloc 3072

      ! extraire le bloc
      CALL XTRBLK (3072,-1,BUF1,LSTELE,TBLVAL,ELDALT,DONIALT, &
                  NELE,NVAL,NT,BLKNO)    
      IF(BLKNO .LE. 0) THEN
        WRITE(*,*)'GENERAL INFO BLOCK NOT FOUND'
        CALL ABORT()
      ENDIF

      ! extraire l'indicateur d'identification du satellite; element 0 01 007.
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   1007,ISAT,ZDATA,IPNTR) 
      IF(IPNTR .EQ. 0) THEN
        WRITE(*,*)' SATELLITE IDENTIFIER MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_debug ) THEN
        WRITE(*,*) ' ISAT = ', (ISAT(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! extraire l'indicateur terre/mer; element 8012.
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   8012,ITERMER,ZDATA,IPNTR)    
      IF(IPNTR .EQ. 0) THEN
        WRITE(*,*)'LAND/SEA INDICATOR MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_debug ) THEN
        WRITE(*,*) ' ITERMER = ', (ITERMER(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! extraire technique de traitement.

      ! not done anymore, j.h. june 2001

      ! extraire le "beam position (field of view no.)"; element 5043.
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   5043,ISCNPOS,ZDATA,IPNTR)
      IF(IPNTR .EQ. 0) THEN
        WRITE(*,*)'INFORMATION ON SCAN POSITION MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_debug ) THEN
        WRITE(*,*) ' ISCNPOS = ', (ISCNPOS(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! extraire le "satellite zenith angle"; element 7024.
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   7024,IDATA,SATZEN,IPNTR)
      IF(IPNTR .EQ. 0) THEN
        WRITE(*,*)'SATELLITE ZENITH ANGLE MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_debug ) THEN
        WRITE(*,*) ' SATZEN = ', (SATZEN(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! le "terrain type": pour les donnees 1b; element 13039.
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   13039,ITERRAIN,ZDATA,IPNTR)
      IF(IPNTR .EQ. 0 .AND. mwbg_debug ) THEN
        WRITE(*,*)'WARNING: TERRAIN TYPE MISSING'
      ENDIF
      IF ( mwbg_debug ) THEN
        WRITE(*,*) ' ITERRAIN = ', (ITERRAIN(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! extraire le numero d'orbite; element 05040.
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   05040,IORBIT,ZDATA,IPNTR)
      IF(IPNTR .EQ. 0) THEN
        WRITE(*,*)'ORBIT NUMBER MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_debug ) THEN
        WRITE(*,*) ' IORBIT = ', (IORBIT(JJ),JJ=1,NVAL*NT)
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
        WRITE(*,*)'RADIANCE DATA BLOCK NOT FOUND'
        CALL ABORT()
      ENDIF

      ! extraire les canaux; element 2150 
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   2150,ICANO,ZDATA,IPNTR)
      IF(IPNTR .EQ. 0) THEN
        WRITE(*,*)'CHANNELS FOR RADIANCE OBS. MISSING'
        CALL ABORT()
      ENDIF
      INO = NVAL
      IF ( mwbg_debug ) THEN
        WRITE(*,*) ' INO  = ',  INO
        WRITE(*,*) ' ICANO= ', (ICANO(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! extraire les radiances.
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   12163,IDATA,ZO,IPNTR) 
      IF(IPNTR .EQ. 0) THEN
        WRITE(*,*)'RADIANCE OBS. MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_debug ) THEN
        WRITE(*,*) ' ZO = ', (ZO(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! extraire la correction aux radiances, element 012233.
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   012233,IDATA,ZCOR,IPNTR)
      IF(IPNTR .EQ. 0) THEN
        IF ( mwbg_debug ) THEN
          WRITE(*,*)'RADIANCE CORRECTION MISSING'
        ENDIF
        DO I = 1, NVAL*NT
          ZCOR(I) = ZMISG
        ENDDO
      ENDIF
      IF ( mwbg_debug ) THEN
        WRITE(*,*) ' ZCOR = ', (ZCOR(JJ),JJ=1,NVAL*NT)
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
        WRITE(*,*)'RADIANCE DATA FLAG BLOCK NOT FOUND'
        CALL ABORT()
      ENDIF

      ! extraire les marqueurs de 13bits des radiances; element 212163 
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   212163,IMARQ,ZDATA,IPNTR)
      IF(IPNTR .EQ. 0) THEN
        WRITE(*,*)'RADIANCE DATA FLAGS MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_debug ) THEN
        WRITE(*,*) ' RADIANCE FLAGS = ', (IMARQ(JJ),JJ=1,NVAL*NT)
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
        WRITE(*,*)'WARNING: (O-P) RADIANCE BLOCK NOT FOUND'
        SKIPENR = .TRUE.
        RETURN
      ENDIF

      ! extraire les canaux; element 2150 
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   2150,ICANOMP,ZDATA,IPNTR) 
      IF(IPNTR .EQ. 0) THEN
        WRITE(*,*)'CHANNELS FOR (O-P) RADIANCE MISSING'
        CALL ABORT()
      ENDIF
      INOMP = NVAL
      IF ( mwbg_debug ) THEN
        WRITE(*,*) ' INOMP   = ',  INOMP
        WRITE(*,*) ' ICANOMP = ', (ICANOMP(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! extraire les residus (O-P) en radiances.
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   12163,IDATA,ZOMP,IPNTR)
      IF(IPNTR .EQ. 0) THEN
        WRITE(*,*)'(O-P) RADIANCES MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_debug ) THEN
        WRITE(*,*) ' ZOMP = ', (ZOMP(JJ),JJ=1,NVAL*NT)
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

    ! 1) Bloc info (general): bloc 3072
    !    Ajouter le "terrain type".

    ! extraire le bloc
    CALL XTRBLK (3072,-1,KBUF1,KLISTE,KTBLVAL,KDLISTE,PRVAL, &
                INELE,INVAL,INT,IBLKNO)    
    IF(IBLKNO .LE. 0) THEN
      WRITE(*,*)'INFO BLOCK NOT FOUND'
      CALL ABORT()
    ENDIF

    !  Dans les conditions suivantes:
    !      1) l'indicateur terre/mer indique l'ocean (ktermer=1),
    !      2) le "terrain type" est manquant (iterrain=-1),
    !      3) le modele indique de la glace (gl >= 0.01),
    !  on specifie "sea ice" pour le "terrain type" (iterrain=0).
    IF ( mwbg_debug ) THEN
      WRITE(*,*) ' OLD TERRAIN TYPE = ', &
                (ITERRAIN(JJ),JJ=1,INT)
      WRITE(*,*) ' KTERMER = ', &
                (KTERMER(JJ),JJ=1,INT)
      WRITE(*,*) ' GLINTRP = ', &
                (GLINTRP(JJ),JJ=1,INT)
    ENDIF
    DO JI = 1, INT
      IF ( KTERMER (JI) .EQ. 1    .AND. &
          ITERRAIN(JI) .EQ. -1   .AND. &
          GLINTRP (JI) .GE. 0.01       ) THEN
        ITERRAIN(JI) = 0
      ENDIF
    ENDDO
    IF ( mwbg_debug ) THEN
      WRITE(*,*) ' NEW TERRAIN TYPE = ', &
                  (ITERRAIN(JJ),JJ=1,INT)
    ENDIF

    !  Le "terrain type" (element 13039) existe-t-il dans le bloc? 
    !       si non, ajouter cet element, 
    !       si oui, remplacer cet element.
    CALL XTRDATA (KDLISTE,KTBLVAL,PRVAL,INELE,INVAL,INT, &
                 13039,KDATA,PDATA,IPNTR)    
    IF(IPNTR .EQ. 0) THEN
      IF ( mwbg_debug ) THEN
        WRITE(*,*)'TERRAIN TYPE MISSING IN BLOCK. ADD IT'
        WRITE(*,*) ' KTBLVAL  = ', (KTBLVAL (JJ),JJ=1, &
                                      INELE*INVAL*INT)
      ENDIF
      CALL ADDDAT (KLISTE,KTBLVAL,PRVAL,INELE,INVAL,INT, &
                  13039,ITERRAIN,PDATA, &
                  KLISTEN,KTBLVALN,PRVALN,INELEN,0)
      IF ( mwbg_debug ) THEN
        WRITE(*,*) ' KTBLVALN = ', (KTBLVALN(JJ),JJ=1, &
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
    IF ( mwbg_debug ) THEN
      WRITE(*,*) ' NEW LAND/SEA = ', (KTERMER(JJ),JJ=1,INVAL*INT)
    ENDIF

    ! Remplacer le bloc.
    !    ISTAT = MRBREP (KBUF1,IBLKNO,KTBLVAL)

    RETURN
  END


  SUBROUTINE mwbg_ADDTRRNF90(KTERMER, ITERRAIN, GLINTRP, rpt)
    !OBJET          Ajouter l'element "terrain type" au bloc info d'un
    !               fichier burp.
    !
    !ARGUMENTS      ktermer - input  -  indicateur terre/mer
    !               iterrain- input  -  indicateur du type de terrain
    !               glintrp - input  -  etendue de glace du modele
    !               rpt     - in/out -  tableau contenant le rapport
    IMPLICIT NONE

    INTEGER :: KTERMER(:)
    INTEGER :: ITERRAIN(:)
    REAL :: GLINTRP(:)
    type(BURP_RPT) :: rpt

    type(BURP_BLOCK) :: blk, blk_copy 

    INTEGER :: KDATA(MXVAL*MXNT)
    REAL :: PDATA(MXVAL*MXNT)

    integer :: error, ref_blk, my_nt,  my_nval, my_nele, my_idtyp
    integer :: indice, indice1, indice2, kk, jj, JI, j, ipos, idata  
    integer :: bktyp, new_bktyp 

    Call BURP_Init(blk, IOSTAT=error)

    ! 1) Bloc info (general): bloc 3072
    !    Ajouter le "terrain type".

    ! extraire le bloc
    ref_blk = 0
    ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 3072, &
                    IOSTAT      = error)
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      write(*,*) 'ERREUR -  INFO block (btyp=3072) not found.'
      call abort()
    endif

    call BURP_Get_Property(blk, &
                     NELE = my_nele, &
                     NT   = my_nt, &    ! number of locations in the box (report)
                     NVAL = my_nval, &
                     IOSTAT=error)
    if (error /= burp_noerr)  call abort()

    !  Dans les conditions suivantes:
    !      1) l'indicateur terre/mer indique l'ocean (ktermer=1),
    !      2) le "terrain type" est manquant (iterrain=-1),
    !      3) le modele indique de la glace (gl >= 0.01),
    !  on specifie "sea ice" pour le "terrain type" (iterrain=0).
    IF ( mwbg_debug ) THEN
      WRITE(*,*) ' OLD TERRAIN TYPE = ', (ITERRAIN(JJ),JJ=1,my_nt)
      WRITE(*,*) ' KTERMER = ', (KTERMER(JJ),JJ=1,my_nt)
      WRITE(*,*) ' GLINTRP = ', (GLINTRP(JJ),JJ=1,my_nt)
    ENDIF
    DO JI = 1, my_nt
      IF ( KTERMER (JI) == 1 .and. ITERRAIN(JI) == -1 .and. GLINTRP (JI) >= 0.01 ) &
            ITERRAIN(JI) = 0
    ENDDO
    IF ( mwbg_debug ) THEN
      WRITE(*,*) ' NEW TERRAIN TYPE = ', (ITERRAIN(JJ),JJ=1,my_nt)
    ENDIF

    indice = BURP_Find_Element(blk,13039,IOSTAT=error)
    if ( indice <= 0 ) then
      write(*,*) 'TERRAIN TYPE MISSING in INFO block (btyp=3072).'
      call abort()
    endif

    ! Update terrain type dans le tableau.
    j = 1
    do kk = 1, my_nt
      idata = ITERRAIN(kk)
      Call BURP_Set_Tblval(blk,indice,j,kk,idata)
    end do
    blk_copy = blk

    Call BURP_Delete_Block(rpt,BLOCK=blk,IOSTAT = error)
    if (error /= burp_noerr)  call abort()

    Call BURP_Write_Block(rpt,BLOCK=blk_copy,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
    if (error /= burp_noerr)  call abort()

    ! les indicateurs terre/mer corriges
    IF ( mwbg_debug ) THEN
      WRITE(*,*) ' NEW LAND/SEA = ', (KTERMER(JJ),JJ=1,my_nt)
    ENDIF

  END subroutine mwbg_ADDTRRNF90


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
    IF ( mwbg_debug ) THEN
      WRITE(*,*) ' ADDDAT: KTBLVALN =  ',(KTBLVALN(JJ), &
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
    ! Notes: In the case where an output parameter cannot be calculated, the
    !        value of this parameter is to to the missing value, i.e. -99.

    IMPLICIT NONE

    integer ni, i

    integer ier    (:)
    integer ilansea(:)
    integer rain   (:)
    integer snow   (:)

    real zmisgLocal, siw, sil, df1, df2, df3, a, b, c, d, e23
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

    data zmisgLocal   / -99.     /
    data epsilon /   1.E-30 /

    ! 1) Initialise output parameters:
    do i = 1, ni
      ice  (i) = zmisgLocal
      tpw  (i) = zmisgLocal
      clw  (i) = zmisgLocal
      scatl(i) = zmisgLocal
      scatw(i) = zmisgLocal
      rain (i) = nint(zmisgLocal)
      snow (i) = nint(zmisgLocal)
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

    !3) Compute parameters:
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

          !3.1) Sea ice:
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
            tpw(i) = zmisgLocal
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

          !3.3) Cloud liquid water:
          ! identify and remove sea ice
          if ( abslat .gt. 50.  .and. &
              df1    .gt.  0.0        ) then  
            clw(i) = zmisgLocal
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

          !3.4) Ocean rain: 0=no rain; 1=rain.
          ! identify and remove sea ice
          if ( abslat .gt. 50.  .and. &
              df1    .gt.  0.0        ) then  
            rain(i) = nint(zmisgLocal)
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

    INTEGER JPMXSFC
    PARAMETER (JPMXSFC =  2)

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


  SUBROUTINE mwbg_tovCheckAtms(KSAT, KORBIT, KCANO, KCANOMP, PTBO, PTBCOR, &
                               PTBOMP, ICHECK, KNO, KNOMP, KNT, PMISG, KNOSAT, IDENTF, &
                               KCHKPRF, ISCNPOS, MTINTRP, IMARQ, STNID, RESETQC)

    !OBJET          Effectuer le controle de qualite des radiances tovs.
    !
    !
    !ARGUMENTS      ksat    - input  -  numero d'identificateur du satellite
    !               korbit  - input  -  numero d'orbite
    !               kcano   - input  -  canaux des observations
    !               kcanomp - input  -  canaux des residus (o-p)
    !               ptbcor  - input  -  correction aux radiances
    !               ptbo    - input  -  radiances
    !               ptbomp  - input  -  residus (o-p)
    !               icheck  - output -  indicateur controle de qualite tovs par canal 
    !                                   =0, ok,
    !                                   >0, rejet,
    !               kno     - input  -  nombre de canaux des observations
    !               knomp   - input  -  nombre de canaux des residus (o-p)
    !               knt     - input  -  nombre de tovs
    !               pmisg   - input  -  valeur manquante burp
    !               knosat  - input  -  numero de satellite (i.e. indice)
    !               identf  - input  -  special ATMS QC integer flag
    !               kchkprf - output -  indicateur global controle de qualite tovs. Code:
    !                                   =0, ok,
    !                                   >0, rejet d'au moins un canal.
    !               iscnpos - input  -  position sur le "scan"
    !               mtintrp - input  -  topographie du modele
    !               imarq   - in/out -  marqueurs des radiances
    !               stnid   - input  -  identificateur du satellite
    !               resetqc - input  -  reset du controle de qualite?
    !
    !NOTES  
    !        5 tests are done:  
    !          1) check for data rejected by first bgckAtms program (SATQC_ATMS standalone program) (QC flag bit 7 ON)          --> bit 9(,7)
    !          2) topography rejection for low-peaking channels (with MF/MX field), --> bits 9,18
    !          3) check for uncorrected radiance (QC flag bit 6 OFF),               --> bit 11           
    !          4) Innovation (O-P) based QC                                         --> bit 9
    !          5) channel blacklisting (from UTIL column in stats_atms_assim file)  --> bit 8
    IMPLICIT NONE

    INTEGER MXCHN
    INTEGER MXSCANAMSU, MXSFCREJ, MXCH2OMPREJ,  MXTOPO
    PARAMETER  ( MXCHN      = 22 )
    PARAMETER  ( MXSCANAMSU = 96 )
    PARAMETER  ( MXSFCREJ   =  8 )
    PARAMETER  ( MXCH2OMPREJ=  6 )
    PARAMETER  ( MXTOPO     =  5 )

    !   Number of tests
    INTEGER     JPMXREJ
    PARAMETER ( JPMXREJ = 5)

    INTEGER JPNSAT,JPCH
 
    PARAMETER (JPNSAT =  3) 
    PARAMETER (JPCH   = 22)
    
    INTEGER NCHNA    (JPNSAT)
    INTEGER MLISCHNA (JPCH,JPNSAT)
    REAL    TOVERRST (JPCH,JPNSAT)
    INTEGER IUTILST  (JPCH,JPNSAT)
    COMMON /COMTOVST/ NCHNA, MLISCHNA, TOVERRST, IUTILST

    INTEGER KNO,KNOMP,KNT,KNOSAT,MAXVAL
    INTEGER JI,JJ,INDX8,INDX12,INO,ICHN
    INTEGER JK,IBIT,JC,INDX,INDXCAN,INDXTOPO

    INTEGER KSAT    (KNT)
    INTEGER ISCNPOS (KNT)
    INTEGER KORBIT  (KNT)
    INTEGER KCANO   (KNO  ,KNT)
    INTEGER KCANOMP (KNOMP,KNT)
    INTEGER ICHECK  (KNO  ,KNT)
    INTEGER KCHKPRF (KNT)
    INTEGER IMARQ   (KNO  ,KNT)
    INTEGER ISFCREJ (MXSFCREJ)
    INTEGER ICH2OMPREJ(MXCH2OMPREJ)
    INTEGER IDENTF  (KNT)
    INTEGER B7CHCK  (KNO  ,KNT)

    REAL  PMISG
    REAL  APPROXIM, XCHECKVAL
    REAL  PTBO    (KNO    ,KNT)
    REAL  PTBCOR  (KNO    ,KNT)
    REAL  PTBOMP  (KNOMP  ,KNT)
    REAL  MTINTRP (KNT)
    REAL  ROGUEFAC(MXCHN)

    INTEGER ITEST  (JPMXREJ)

    INTEGER ICHTOPO (MXTOPO)
    REAL    ZCRIT   (MXTOPO)

    CHARACTER *9   STNID

    LOGICAL LLFIRST,GROSSERROR,FULLREJCT,RESETQC,SFCREJCT,CH2OMPREJCT

    INTEGER  MREJCOD,INTOT,INTOTRJF,INTOTRJP,MREJCOD2
    COMMON /STATS/  MREJCOD(JPMXREJ,MXCHN,MXSAT), &
                   INTOT(MXSAT), &
                   INTOTRJF(MXSAT),INTOTRJP(MXSAT), &
                   MREJCOD2(JPMXREJ,MXCHN,MXSAT)

    LOGICAL DEBUG
    COMMON /DBGCOM/ DEBUG

    SAVE LLFIRST

    DATA  LLFIRST / .TRUE. /

    DATA  ROGUEFAC / 2.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 4.0, &
                    4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 2.0, &
                    2.0, 4.0, 4.0, 4.0, 4.0, 4.0 /

    ! Channel sets for rejection in test 9 
    ! These LT channels are rejected if O-P fails rogue check for window ch. 1, 2, or 3
    DATA  ISFCREJ    / 1, 2, 3, 4, 5, 6, 16, 17 /
    !   These AMSU-B channels are rejected if ch. 17 O-P fails rogue check over OPEN WATER only    
    DATA  ICH2OMPREJ / 17, 18, 19, 20, 21, 22 /

    !  Data for TOPOGRAPHY CHECK
    !   Channel AMSUA-6 (atms ch 7) is rejected for topography  >  250m.
    !   Channel AMSUA-7 (atms ch 8) is rejected for topography  > 2000m.
    !   Channel AMSUB-3 (atms ch 22) is rejected for topography > 2500m.
    !                    atms ch 21  is rejected for topography > 2250m.
    !   Channel AMSUB-4 (atms ch 20) is rejected for topography > 2000m.
    DATA ICHTOPO  /     7,     8,    20,    21,    22  /
    DATA ZCRIT    /   250., 2000., 2000., 2250., 2500. /

    !  Test selection (0=skip test, 1=do test)
    !             1  2  3  4  5 
    DATA ITEST  / 1, 1, 1, 1, 1 /
       
    ! Initialisation, la premiere fois seulement!
    IF (LLFIRST) THEN
      DO JI = 1, JPMXREJ
        DO JJ = 1, MXCHN
          DO JK = 1, MXSAT
            MREJCOD(JI,JJ,JK)  = 0
            MREJCOD2(JI,JJ,JK) = 0
          ENDDO
        ENDDO
      ENDDO
      LLFIRST = .FALSE.
    ENDIF

    !  Verification de l'integrite des donnees, c'est-a-dire que:
    !         i)  les dimensions des divers blocs de donnees concordent,
    !         ii) les listes de canaux concordent.
    IF ( KNO .NE. KNOMP     ) THEN
      write(*,*)'ERROR IN DIMENSIONS OF TOVS DATA'
      CALL ABORT()
    ENDIF

    DO JJ=1,KNT
      DO JI=1,KNO
        IF ( KCANO(JI,JJ) .NE. KCANOMP(JI,JJ) ) THEN
          write(*,*)'INCONSISTENT CHANNEL LISTS FOR TOVS DATA'
          CALL ABORT()
        ENDIF
      ENDDO
    ENDDO

    !  Initialisations
    DO JJ=1,KNT
      DO JI=1,KNO
        ICHECK(JI,JJ) = 0
        B7CHCK(JI,JJ) = 0
        IF ( RESETQC ) IMARQ(JI,JJ) = 0
      ENDDO
    ENDDO

    ! 1) test 1: Check flag bit 7 on from the first bgckAtms program
    !  Includes observations flagged for cloud liquid water, scattering index,
    !  dryness index plus failure of several QC checks.
    INO = 1
    if ( itest(ino) .eq. 1 ) then
      DO JJ=1,KNT
        DO JI=1,KNO
          IBIT = AND(IMARQ(JI,JJ), 2**7)
          IF ( IBIT .NE. 0  ) THEN
            ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
            B7CHCK(JI,JJ) = 1
            IMARQ(JI,JJ) = OR(IMARQ(JI,JJ),2**9)
            MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                 MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
            IF (DEBUG) THEN
              write(*,*)STNID(2:9),' first bgckAtms program REJECT.', &
                        'CHANNEL=', KCANO(JI,JJ), &
                        ' IMARQ= ',IMARQ(JI,JJ)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    endif

    ! 2) test 2: Topography check (partial)
    INO = 2
    if ( itest(ino) .eq. 1 ) then
      DO JJ=1,KNT
        DO JI=1,KNO
          INDXTOPO = ISRCHEQI(ICHTOPO, MXTOPO, KCANO(JI,JJ))
          IF ( INDXTOPO .GT. 0 ) THEN
            IF ( MTINTRP(JJ) .GE. ZCRIT(INDXTOPO) ) THEN
              ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
              IMARQ(JI,JJ) = OR(IMARQ(JI,JJ),2**9)
              IMARQ(JI,JJ) = OR(IMARQ(JI,JJ),2**18)
              MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                   MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
              IF ( B7CHCK(JI,JJ) .EQ. 0 ) THEN
                MREJCOD2(INO,KCANO(JI,JJ),KNOSAT) = &
                   MREJCOD2(INO,KCANO(JI,JJ),KNOSAT)+ 1                 
              ENDIF
              IF (DEBUG) THEN
                write(*,*)STNID(2:9),' TOPOGRAPHY REJECT.', &
                          'CHANNEL=', KCANO(JI,JJ), &
                          ' TOPO= ',MTINTRP(JJ)
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    endif

    ! 3) test 3: Uncorrected Tb check (single)
    !  Uncorrected datum (flag bit #6 off). In this case switch bit 11 ON.
    IF (.NOT.RESETQC) THEN
      INO = 3
      if ( itest(ino) .eq. 1 ) then
        DO JJ=1,KNT
          DO JI=1,KNO
            IBIT = AND(IMARQ(JI,JJ), 2**6)
            IF ( IBIT .EQ. 0  ) THEN
              ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
              IMARQ(JI,JJ) = OR(IMARQ(JI,JJ),2**11)
              MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                 MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
              IF ( B7CHCK(JI,JJ) .EQ. 0 ) THEN
                MREJCOD2(INO,KCANO(JI,JJ),KNOSAT) = &
                    MREJCOD2(INO,KCANO(JI,JJ),KNOSAT)+ 1                 
              ENDIF
              IF (DEBUG) THEN
                write(*,*)STNID(2:9),' UNCORRECTED TB REJECT.', &
                          'CHANNEL=', KCANO(JI,JJ), &
                          ' IMARQ= ',IMARQ(JI,JJ)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      endif
    ENDIF

    ! 4) test 4: "Rogue check" for (O-P) Tb residuals out of range. (single/full)
    !             Also, over WATER remove CH.17-22 if CH.17 |O-P|>5K (partial)
    !  Les observations, dont le residu (O-P) depasse par un facteur (roguefac) 
    !   l'erreur totale des TOVS.
    !  N.B.: a reject by any of the 3 amsua surface channels 1-3 produces the 
    !           rejection of ATMS sfc/tropospheric channels 1-6 and 16-17.
    !  OVER OPEN WATER
    !    ch. 17 Abs(O-P) > 5K produces rejection of all ATMS amsub channels 17-22.
    INO = 4
    if ( itest(ino) .eq. 1 ) then
      DO JJ=1,KNT
        SFCREJCT = .FALSE.
        CH2OMPREJCT = .FALSE.
        DO JI=1,KNO
          ICHN = KCANO(JI,JJ)
          XCHECKVAL = ROGUEFAC(ICHN) * &
                     TOVERRST(ICHN,KNOSAT) 
          IF ( PTBOMP(JI,JJ)      .NE. PMISG    .AND. &
              ABS(PTBOMP(JI,JJ)) .GE. XCHECKVAL     ) THEN
            ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
            IMARQ(JI,JJ) = OR(IMARQ(JI,JJ),2**9)
            IMARQ(JI,JJ) = OR(IMARQ(JI,JJ),2**16)
            MREJCOD(INO,KCANO(JI,JJ),KNOSAT) =  &
               MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
            IF ( B7CHCK(JI,JJ) .EQ. 0 ) THEN
              MREJCOD2(INO,KCANO(JI,JJ),KNOSAT) = &
                 MREJCOD2(INO,KCANO(JI,JJ),KNOSAT)+ 1                 
            ENDIF
            IF (DEBUG) THEN
              write(*,*)STNID(2:9),'ROGUE CHECK REJECT.NO.', &
                     ' OBS = ',JJ, &
                     ' CHANNEL= ',KCANO(JI,JJ), &
                     ' CHECK VALUE= ',XCHECKVAL, &
                     ' TBOMP= ',PTBOMP(JI,JJ)
            ENDIF
            IF ( ICHN .EQ. 1 .OR. &
                ICHN .EQ. 2 .OR. &
                ICHN .EQ. 3    ) THEN
              SFCREJCT = .TRUE.
            ENDIF
          ENDIF
          IF ( ICHN .EQ. 17 .AND. PTBOMP(JI,JJ) .NE. PMISG .AND. &
              ABS(PTBOMP(JI,JJ)) .GT. 5.0 ) THEN
            CH2OMPREJCT = .TRUE.
          ENDIF
        ENDDO

        IF ( SFCREJCT ) THEN
          DO JI=1,KNO
            INDXCAN = ISRCHEQI (ISFCREJ,MXSFCREJ,KCANO(JI,JJ))
            IF ( INDXCAN .NE. 0 )  THEN
              IF ( ICHECK(JI,JJ) .NE. INO ) THEN
                ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
                IMARQ(JI,JJ) = OR(IMARQ(JI,JJ),2**9)
                IMARQ(JI,JJ) = OR(IMARQ(JI,JJ),2**16)
                MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                        MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
                IF ( B7CHCK(JI,JJ) .EQ. 0 ) THEN
                  MREJCOD2(INO,KCANO(JI,JJ),KNOSAT) = &
                     MREJCOD2(INO,KCANO(JI,JJ),KNOSAT)+ 1                 
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDIF

        !  amsub channels 17-22 obs are rejected if, for ch17 ABS(O-P) > 5K
        !    Apply over open water only (bit 0 ON in QC integer identf).
        !    Only apply if obs not rejected in this test already.
        IBIT = AND(IDENTF(JJ), 2**0)
        IF ( CH2OMPREJCT .AND. (IBIT .NE. 0) ) THEN
          DO JI=1,KNO
            INDXCAN = ISRCHEQI (ICH2OMPREJ,MXCH2OMPREJ,KCANO(JI,JJ))
            IF ( INDXCAN .NE. 0 )  THEN
              IF ( ICHECK(JI,JJ) .NE. INO ) THEN
                ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
                IMARQ(JI,JJ) = OR(IMARQ(JI,JJ),2**9)
                IMARQ(JI,JJ) = OR(IMARQ(JI,JJ),2**16)
                MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                        MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
                IF ( B7CHCK(JI,JJ) .EQ. 0 ) THEN
                  MREJCOD2(INO,KCANO(JI,JJ),KNOSAT) = &
                     MREJCOD2(INO,KCANO(JI,JJ),KNOSAT)+ 1                 
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDIF

      ENDDO
    endif

    ! 5) test 5: Channel selection using array IUTILST(chan,sat)
    !  IUTILST = 0 (blacklisted)
    !            1 (assmilate)
    INO = 5
    if ( itest(ino) .eq. 1 ) then
      DO JJ=1,KNT
        DO JI=1,KNO
           ICHN = KCANO(JI,JJ)
           IF ( IUTILST(ICHN,KNOSAT) .EQ. 0 ) THEN
             IMARQ(JI,JJ) = OR(IMARQ(JI,JJ),2**8)
             MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                   MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
             IF (DEBUG) THEN
               write(*,*)STNID(2:9),'CHANNEL REJECT: ', &
                      ' OBS = ',JJ, &
                      ' CHANNEL= ',ICHN                  
             ENDIF
           ENDIF
        ENDDO
      ENDDO      
    endif

    IF (DEBUG) THEN
      write(*,*) 'ICHECK = ',((ICHECK(JI,JJ),JI=1,KNO),JJ=1,KNT)
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
      write(*,*)'KCHKPRF = ',(KCHKPRF(JJ),JJ=1,KNT)
    ENDIF 

    RETURN
  END


  SUBROUTINE mwbg_qcStatsAtms(INUMSAT, ICHECK, KCANO, KNOSAT, CSATID, KNO, &
                              KNT, LDPRINT)
    !OBJET          Cumuler ou imprimer des statistiques decriptives
    !               des rejets tovs.
    !
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

    INTEGER MXCHN
    PARAMETER  ( MXCHN = 22 )

    INTEGER     JPMXREJ
    PARAMETER ( JPMXREJ = 5)

    CHARACTER*9   CSATID (MXSAT)

    INTEGER  JI, JJ, JK, KNO, KNT, KNOSAT
    INTEGER  INTOTOBS, INTOTACC, INUMSAT

    INTEGER  ICHECK (KNO,KNT)
    INTEGER  KCANO  (KNO,KNT)

    INTEGER  MREJCOD,INTOT,INTOTRJF,INTOTRJP,MREJCOD2
    COMMON /STATS/  MREJCOD(JPMXREJ,MXCHN,MXSAT), &
                   INTOT(MXSAT), &
                   INTOTRJF(MXSAT),INTOTRJP(MXSAT), &
                   MREJCOD2(JPMXREJ,MXCHN,MXSAT)

    LOGICAL LLFIRST, LDPRINT, FULLREJCT, FULLACCPT

    SAVE  LLFIRST

    DATA  LLFIRST / .TRUE. /

    ! Initialize
    IF ( LLFIRST ) THEN
      DO JJ = 1, MXSAT
        INTOTRJF(JJ) = 0
        INTOTRJP(JJ) = 0
        INTOT(JJ)  = 0
        LLFIRST = .FALSE.
      ENDDO
    ENDIF

    IF (.NOT. LDPRINT ) THEN

      ! Accumulate statistics on rejects
      DO JJ = 1, KNT
        INTOT(KNOSAT) = INTOT(KNOSAT) + 1

        ! Full accepted, fully rejected or partially rejected?
        FULLREJCT = .TRUE.
        FULLACCPT = .TRUE.
        DO JI = 1, KNO
           IF ( ICHECK(JI,JJ) .NE. 0 ) THEN
              FULLACCPT = .FALSE.
           ELSE
              FULLREJCT = .FALSE.
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

        write(*,'(/////50("*"))')
        write(*,'(     50("*")/)')
        write(*,'(T5,"SUMMARY OF QUALITY CONTROL FOR ", &
         A8)') CSATID(JK) 
        write(*,'(T5,"------------------------------------- ",/)')
        write(*,'( &
        "   TOTAL NUMBER OF ATMS    = ",I10,/ &
        " - TOTAL FULL REJECTS      = ",I10,/ &
        " - TOTAL PARTIAL REJECTS   = ",I10,/ &
        "   ------------------------------------",/ &
        "   TOTAL FULLY ACCEPTED    = ",I10,/)') &
          INTOTOBS, INTOTRJF(JK), INTOTRJP(JK), INTOTACC

        write(*,'(//,1x,59("-"))')
        write(*,'(t10,"|",t19,"1. REJECTION CATEGORIES")')
        write(*,'(" CHANNEL",t10,"|",50("-"))')
        write(*,'(t10,"|",5i7)') (JI,JI=1,JPMXREJ)
        write(*,'(1x,"--------|",50("-"))')
        DO JJ = 1, MXCHN
          write(*,'(3X,I2,t10,"|",5I7)') JJ,(MREJCOD(JI,JJ,JK), &
                                      JI=1,JPMXREJ)
        ENDDO
        write(*,'(1x,59("-"))')

        write(*,'(//,1x,59("-"))')
        write(*,'(t10,"|",t19,"2. QC2 REJECT CATEGORIES")')
        write(*,'(" CHANNEL",t10,"|",50("-"))')
        write(*,'(t10,"|",5i7)') (JI,JI=1,JPMXREJ)
        write(*,'(1x,"--------|",50("-"))')
        DO JJ = 1, MXCHN
          write(*,'(3X,I2,t10,"|",5I7)') JJ,(MREJCOD2(JI,JJ,JK), &
                                      JI=1,JPMXREJ)
        ENDDO
        write(*,'(1x,59("-"))')
      ENDDO

      ! Print legend
      PRINT *, ' '
      PRINT *, ' '
      PRINT *, ' -----------------------------------------------------'
      PRINT *, ' Definition of rejection categories: '
      PRINT *, ' -----------------------------------------------------'
      PRINT *, '  1 - first bgckAtms program reject [bit 7]'
      PRINT *, '  2 - topography reject'
      PRINT *, '  3 - uncorrected radiance'
      PRINT *, '  4 - innovation (O-P) based reject'
      PRINT *, '  5 - rejection by channel selection'
      PRINT *, ' -----------------------------------------------------'
      PRINT *, ' '
      PRINT *, ' QC2 REJECT numbers in Table 2 are for data that '
      PRINT *, ' passed test 1 (data with QC flag bit 7 OFF)'
      PRINT *, ' '

    ENDIF

    RETURN
  END


  SUBROUTINE mwbg_updatFlgAtms(KCHKPRF, ICHECK, RESETQC, IMARQ, rpt)
    !OBJET          Allumer les bits des marqueurs pour les tovs rejetes.
    !               Modifier le bktyp des donnees, marqueurs et (O-P) pourt
    !               signifier "vu par AO". 
    !
    !ARGUMENTS      kchprf  - input  -  indicateur global controle de qualite tovs. Code:
    !                                   =0, ok,
    !                                   >0, rejet,
    !               icheck  - input  -  indicateur controle de qualite tovs au 
    !                                   niveau de chaque canal
    !               resetqc - input  -  reset the quality control flags before adding the new ones ? 
    !               imarq   - input  -  modified flag values from mwbg_tovCheckAtms
    !               rpt     - in/out -  tableau contenant le rapport
    IMPLICIT NONE

    INTEGER KCHKPRF (:)
    INTEGER ICHECK  (:)
    LOGICAL RESETQC
    INTEGER IMARQ   (:)
    type(BURP_RPT) :: rpt

    type(BURP_BLOCK) :: blk, blk_copy 

    integer :: KDATA (MXVAL*MXNT)

    integer :: error, ref_blk, my_nt,  my_nval, my_nele, my_idtyp
    integer :: indice, indice1, indice2, kk, jj, JI, j, ipos, idata  
    integer :: bktyp, new_bktyp 

    Call BURP_Init(blk, IOSTAT=error)

    ! 1) Bloc info 3d: bloc 5120.
    !    Modifier les marqueurs globaux de 24bits pour les donnees rejetees.

    ! extraire le bloc
    ref_blk = 0
    ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 5120, &
                    IOSTAT      = error)
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      write(*,*) 'ERREUR -  Location/time (3D) block (btyp=5120) not found.'
      call abort()
    endif

    call BURP_Get_Property(blk, &
                     NELE = my_nele, &
                     NT   = my_nt, &    ! number of locations in the box (report)
                     NVAL = my_nval, &
                     IOSTAT=error)
    if (error /= burp_noerr)  call abort()

    ! extraire les marqueurs globaux de 24bits; element 55200
    indice = BURP_Find_Element(blk,55200,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk = 1, my_nt
        kdata(kk) = BURP_Get_Tblval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'GLOBAL FLAGS missing in 3D block (btyp=5120).'
      call abort()
    endif
    IF (mwbg_debug) THEN
      write(*,*) ' OLD FLAGS = ', (KDATA(JJ),JJ=1,my_nt)
    ENDIF

    ! allumer la bit (6) indiquant que l'observation a un element
    ! rejete par le controle de qualite de l'AO.
    !  N.B.: si on est en mode resetqc, on remet le marqueur global a
    !        sa valeur de defaut, soit 1024,  avant de faire la mise a jour.
    DO JI = 1, my_nt
      IF (RESETQC) THEN
        KDATA(JI) = 1024  
      ENDIF
      IF ( KCHKPRF(JI).NE.0  ) THEN
        KDATA(JI) = OR (KDATA(JI),2**6)
      ENDIF
    ENDDO
    IF (mwbg_debug) THEN
      write(*,*) ' KCHKPRF   = ', (KCHKPRF(JJ),JJ=1,my_nt)
      write(*,*) ' NEW FLAGS = ', (KDATA(JJ),JJ=1,my_nt)
    ENDIF

    ! Remplacer les nouveaux marqueurs dans le tableau.
    j = 1
    do kk = 1, my_nt
      idata = kdata(kk)
      Call BURP_Set_Tblval(blk,indice,j,kk,idata)
    end do
    blk_copy = blk

    Call BURP_Delete_Block(rpt,BLOCK=blk,IOSTAT = error)
    if (error /= burp_noerr)  call abort()

    Call BURP_Write_Block(rpt,BLOCK=blk_copy,CONVERT_BLOCK=.FALSE.,IOSTAT=error)
    if (error /= burp_noerr)  call abort()

    ! 3) Bloc multi niveaux de radiances: bloc 9218, 9248, 9264.
    !    Modifier le bktyp pour signifier "vu par AO".

    ! localiser le bloc
    ref_blk = 0
    ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 9218, &
                    IOSTAT      = error)
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      ref_blk = 0
      ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 9248, &
                    IOSTAT      = error)
    end if
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      ref_blk = 0
      ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0 , &
                    BTYP        = 9264, &
                    IOSTAT      = error)
      if (ref_blk < 0) then
        write(*,*) 'ERREUR - RADIANCE DATA BLOCK NOT FOUND'
        call abort()
      endif
    endif

    call BURP_Get_Property(blk, &
                     NELE = my_nele, &
                     NT   = my_nt, &    ! number of locations in the box (report)
                     NVAL = my_nval, &
                     BKTYP = bktyp, & 
                     IOSTAT=error)
    if (error /= burp_noerr)  call abort()

    ! Remplacer le bloc et modifier le bktyp pour signifier "vu par AO".
    new_bktyp = BKTYP + 4
    Call BURP_Set_Property(blk,BKTYP=new_bktyp)
    blk_copy = blk

    Call BURP_Delete_Block(rpt,BLOCK=blk,IOSTAT = error)
    if (error /= burp_noerr)  call abort()

    Call BURP_Write_Block(rpt,BLOCK=blk_copy,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
    if (error /= burp_noerr)  call abort()

    ! 4) Bloc marqueurs multi niveaux de radiances: bloc 15362, 15392, 15408.
    !    Modifier les marqueurs de 13bits associes a chaque radiance.
    !    Modifier le bktyp pour signifier "vu par AO".

    ! extraire le bloc
    ref_blk = 0
    ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 15362, &
                    IOSTAT      = error)
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      ref_blk = 0
      ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 15392, &
                    IOSTAT      = error)
    end if
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      ref_blk = 0
      ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0 , &
                    BTYP        = 15408, &
                    IOSTAT      = error)
      if (ref_blk < 0) then
        write(*,*) 'ERREUR - RADIANCE DATA BLOCK NOT FOUND'
        call abort()
      endif
    endif

    ! extraire les marqueurs de 13bits des radiances; element 212163 (LEVEL 1B)
    call BURP_Get_Property(blk, &
                      NELE = my_nele, &
                      NT   = my_nt, &    ! number of locations in the box (report)
                      NVAL = my_nval, &
                      BKTYP = bktyp, & 
                      IOSTAT=error)
    if (error /= burp_noerr)  call abort()

    indice = BURP_Find_Element(blk,212163,IOSTAT=error)
    if ( indice > 0 ) then
      ipos = 0
      do kk = 1, my_nt
        do j = 1, my_nval
          ipos = ipos + 1
          KDATA(ipos) = BURP_Get_Tblval(blk,indice,j,kk,error)
        end do
      end do
    else
      write(*,*) 'ERREUR - Element 212163 missing in flag block.'
      call abort()
    endif 
    IF (mwbg_debug) THEN
      write(*,*) ' OLD FLAGS = ', (KDATA(JJ),JJ=1,my_nval*my_nt)
    ENDIF

    ! update data flags
    DO JI = 1, my_nval * my_nt
      KDATA(JI) = IMARQ(JI)
    ENDDO
    IF (mwbg_debug) THEN
      write(*,*) ' ICHECK = ', (ICHECK(JJ),JJ=1,my_nval*my_nt)
      write(*,*) ' NEW FLAGS = ', (KDATA(JJ),JJ=1,my_nval*my_nt)
    ENDIF

    ! Remplacer les nouveaux marqueurs dans le tableau.
    ipos = 0
    do kk =1, my_nt
      do j = 1, my_nval
        ipos = ipos + 1
        Call BURP_Set_Tblval(blk,indice,j,kk,KDATA(ipos))
      enddo
    enddo

    ! Remplacer le bloc et modifier le bktyp pour signifier "vu par AO".
    new_bktyp = BKTYP + 4
    Call BURP_Set_Property(blk,BKTYP=new_bktyp,IOSTAT=error)
    if (error /= burp_noerr)  call abort()
    blk_copy = blk

    Call BURP_Delete_Block(rpt,BLOCK=blk,IOSTAT = error)
    if (error /= burp_noerr)  call abort()

    Call BURP_Write_Block(rpt,BLOCK=blk_copy,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
    if (error /= burp_noerr)  call abort()

    ! 5) Bloc multi niveaux de residus de radiances (O-P): bloc 9322, 9226, 9258, 9274, bfam 14
    !    Modifier le bktyp pour signifier "vu par AO".

    ! localiser le bloc
    ref_blk = 0
    ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 14, &
                    BTYP        = 9322, &
                    IOSTAT      = error)
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      ref_blk = 0
      ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 14, &
                    BTYP        = 9226, &
                    IOSTAT      = error)
    end if
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      ref_blk = 0
      ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 14, &
                    BTYP        = 9258, &
                    IOSTAT      = error)
    end if
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      ref_blk = 0
      ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 14, &
                    BTYP        = 9274, &
                    IOSTAT      = error)
      if (ref_blk < 0) then
        write(*,*) 'ERREUR - OMP DATA block (btyp 9322 or 9226 or 9258 or 9274) not found.'
        call abort()
      endif
    endif

    call BURP_Get_Property(blk, &
                      NELE = my_nele, &
                      NT   = my_nt, &    ! number of locations in the box (report)
                      NVAL = my_nval, &
                      BKTYP = bktyp, & 
                      IOSTAT=error)
    if (error /= burp_noerr)  call abort()

    ! Remplacer le bloc et modifier le bktyp pour signifier "vu par AO".
    new_bktyp = BKTYP + 4
    Call BURP_Set_Property(blk,BKTYP=new_bktyp)
    blk_copy = blk

    Call BURP_Delete_Block(rpt,BLOCK=blk,IOSTAT = error)
    if (error /= burp_noerr)  call abort()

    Call BURP_Write_Block(rpt,BLOCK=blk_copy,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
    if (error /= burp_noerr)  call abort()

  END subroutine mwbg_updatFlgAtms


  SUBROUTINE mwbg_readTovsAtms(IUNENT,HANDLE,ISAT,ZMISG,BUF1,TBLVAL, &
       LSTELE,ELDALT,IDATA,ICANO,ICANOMP,ICANOMPNA,INO, &
       INOMP,INOMPNA,DONIALT,ZDATA,ZO,ZCOR,ZOMP,STNID, &
       ZOMPNA,IMARQ,MXELM,MXVAL,MXNT,NELE,NVAL,NT, &
       ISCNCNT,ISCNPOS,MXSCAN,ZLAT,ZLON,ITERMER, &
       IORBIT,SATZEN,ITERRAIN,SKIPENR,IDENTF)
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
    !               identf    - output -  ATMS special QC flag integer (ele 025174)
    IMPLICIT NONE

    INTEGER MXELM, MXVAL, MXNT, MXSCAN

    INTEGER  IUNENT, HANDLE

    INTEGER MRFGET
    INTEGER MRBHDR,MRFLOC
    INTEGER TEMPS,DATE,OARS,RUNN
    INTEGER IDTYP,LATI,LONG,DX,DY,ELEV,DRND,NBLOCS,BLKNO
    INTEGER ISTAT,NVAL,NT
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
    INTEGER IDENTF   (MXNT)

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

    SKIPENR = .FALSE.

    HANDLE = MRFLOC(IUNENT,HANDLE,'*********',-1,-1,-1,-1, &
                   -1,SUP,0)

    IF(HANDLE .GT. 0) THEN
      IF ( mwbg_DEBUG ) THEN
        write(*,*)'PROCESSING BRP FILE:HANDLE=',HANDLE
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
        ISAT(I)     = MISGINT
      ENDDO

      ! 1) Bloc info 3d: bloc 5120

      ! extraire le bloc
      CALL XTRBLK (5120,-1,BUF1,LSTELE,TBLVAL,ELDALT,DONIALT, &
                  NELE,NVAL,NT,BLKNO)    
      IF(BLKNO .LE. 0) THEN
        write(*,*)'3D INFO BLOCK NOT FOUND'
        CALL ABORT()
      ENDIF

      ! extraire la latitude; element 5002
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   5002,IDATA,ZLAT,IPNTR)    
      IF(IPNTR .EQ. 0) THEN
        write(*,*)'LATITUDE MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_DEBUG ) THEN
        write(*,*) ' LATITUDE = ', (ZLAT(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! extraire la longitude; element 6002
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   6002,IDATA,ZLON,IPNTR)    
      IF(IPNTR .EQ. 0) THEN
        write(*,*)'LONGITUDE MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_DEBUG ) THEN
        write(*,*) ' LONGITUDE = ', (ZLON(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! 2) Bloc info (general): bloc 3072

      ! extraire le bloc
      CALL XTRBLK (3072,-1,BUF1,LSTELE,TBLVAL,ELDALT,DONIALT, &
                  NELE,NVAL,NT,BLKNO)    
      IF(BLKNO .LE. 0) THEN
        write(*,*)'GENERAL INFO BLOCK NOT FOUND'
        CALL ABORT()
      ENDIF

      ! extraire l'indicateur d'identification du satellite; element 0 01 007.
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   1007,ISAT,ZDATA,IPNTR) 
      IF(IPNTR .EQ. 0) THEN
        write(*,*)' SATELLITE IDENTIFIER MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_DEBUG ) THEN
        write(*,*) ' ISAT = ', (ISAT(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! extraire l'indicateur terre/mer; element 8012.
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   8012,ITERMER,ZDATA,IPNTR)    
      IF(IPNTR .EQ. 0) THEN
        write(*,*)'LAND/SEA INDICATOR MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_debug ) THEN
        write(*,*) ' ITERMER = ', (ITERMER(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! extraire technique de traitement.
      !       not done anymore, j.h. june 2001
      ! extraire le "beam position (field of view no.)"; element 5043.
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   5043,ISCNPOS,ZDATA,IPNTR)
      IF(IPNTR .EQ. 0) THEN
        write(*,*)'INFORMATION ON SCAN POSITION MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_debug ) THEN
        write(*,*) ' ISCNPOS = ', (ISCNPOS(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! extraire le "satellite zenith angle"; element 7024.
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   7024,IDATA,SATZEN,IPNTR)
      IF(IPNTR .EQ. 0) THEN
        write(*,*)'SATELLITE ZENITH ANGLE MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_debug ) THEN
        write(*,*) ' SATZEN = ', (SATZEN(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! le "terrain type": pour les donnees 1b; element 13039.
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   13039,ITERRAIN,ZDATA,IPNTR)
      IF(IPNTR .EQ. 0 .AND. mwbg_debug) THEN
        write(*,*)'WARNING: TERRAIN TYPE MISSING'
      ENDIF
      IF ( mwbg_debug ) THEN
        write(*,*) ' ITERRAIN = ', (ITERRAIN(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! extraire le numero d'orbite; element 05040.
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   05040,IORBIT,ZDATA,IPNTR)
      IF(IPNTR .EQ. 0) THEN
        write(*,*)'ORBIT NUMBER MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_debug ) THEN
        write(*,*) ' IORBIT = ', (IORBIT(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! extraire le special qc flag integer; element 25174.
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   25174,IDENTF,ZDATA,IPNTR)
      IF(IPNTR .EQ. 0) THEN
        write(*,*)'SPECIAL QC FLAG INTEGER MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_debug ) THEN
        write(*,*) ' IDENTF = ', (IDENTF(JJ),JJ=1,NVAL*NT)
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
        write(*,*)'RADIANCE DATA BLOCK NOT FOUND'
        CALL ABORT()
      ENDIF

      ! extraire les canaux; element 2150 or 5042 
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   2150,ICANO,ZDATA,IPNTR)
      IF(IPNTR .EQ. 0) THEN
        CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   5042,ICANO,ZDATA,IPNTR)
      ENDIF
      IF(IPNTR .EQ. 0) THEN
        write(*,*)'CHANNELS FOR RADIANCE OBS. MISSING'
        CALL ABORT()
      ENDIF
      INO = NVAL
      IF ( mwbg_debug ) THEN
        write(*,*) ' INO  = ',  INO
        write(*,*) ' ICANO= ', (ICANO(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! extraire les radiances.
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   12163,IDATA,ZO,IPNTR) 
      IF(IPNTR .EQ. 0) THEN
        write(*,*)'RADIANCE OBS. MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_debug ) THEN
        write(*,*) ' ZO = ', (ZO(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! extraire la correction aux radiances, element 012233.
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   012233,IDATA,ZCOR,IPNTR)
      IF(IPNTR .EQ. 0) THEN
        IF ( mwbg_debug ) THEN
          write(*,*)'RADIANCE CORRECTION MISSING'
        ENDIF
        DO I = 1, NVAL*NT
          ZCOR(I) = ZMISG
        ENDDO
      ENDIF
      IF ( mwbg_debug ) THEN
        write(*,*) ' ZCOR = ', (ZCOR(JJ),JJ=1,NVAL*NT)
      ENDIF                                

      ! 4) Bloc marqueurs multi niveaux de radiances: bloc 15362, 15392, 15408.

      ! extraire le bloc
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
        write(*,*)'RADIANCE DATA FLAG BLOCK NOT FOUND'
        CALL ABORT()
      ENDIF

      ! extraire les marqueurs de 13bits des radiances; element 212163 
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   212163,IMARQ,ZDATA,IPNTR)
      IF(IPNTR .EQ. 0) THEN
        write(*,*)'RADIANCE DATA FLAGS MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_debug ) THEN
        write(*,*) ' RADIANCE FLAGS = ', (IMARQ(JJ),JJ=1,NVAL*NT)
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
        write(*,*)'WARNING: (O-P) RADIANCE BLOCK NOT FOUND'
        SKIPENR = .TRUE.
        RETURN
      ENDIF

      ! extraire les canaux; element 2150 or 5042
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   2150,ICANOMP,ZDATA,IPNTR) 
      IF(IPNTR .EQ. 0) THEN
        CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   5042,ICANOMP,ZDATA,IPNTR)
      ENDIF
      IF(IPNTR .EQ. 0) THEN
        write(*,*)'CHANNELS FOR (O-P) RADIANCE MISSING'
        CALL ABORT()
      ENDIF
      INOMP = NVAL
      IF ( mwbg_debug ) THEN
        write(*,*) ' INOMP   = ',  INOMP
        write(*,*) ' ICANOMP = ', (ICANOMP(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! extraire les residus (O-P) en radiances.
      CALL XTRDATA (ELDALT,TBLVAL,DONIALT,NELE,NVAL,NT, &
                   12163,IDATA,ZOMP,IPNTR)
      IF(IPNTR .EQ. 0) THEN
        write(*,*)'(O-P) RADIANCES MISSING'
        CALL ABORT()
      ENDIF
      IF ( mwbg_debug ) THEN
        write(*,*) ' ZOMP = ', (ZOMP(JJ),JJ=1,NVAL*NT)
      ENDIF

      ! 6) Bloc multi niveaux de residus de radiances (O-P) au nadir: bloc 9322 ou 9226, bfam 32.

      ! not done any more, j.h. june 2001
    ENDIF

    RETURN
  END


  SUBROUTINE mwbg_readStatTovsAtms(ILUTOV,INUMSAT,CSATID)
    !OBJET          Read stats_atms_assim file
    !
    !ARGUMENTS      ilutov  - input  -  unite logique du fichier stats des TOVS
    !               inumsat - output -  nombre de satellites
    !               csatid  - output -  identificateur de satellite 
    IMPLICIT NONE

    INTEGER JPNSAT,JPCH
 
    PARAMETER (JPNSAT =  3) 
    PARAMETER (JPCH = 22)

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
9000 FORMAT(//,10x,"-mwbg_readStatTovsAtms: reading total error statistics" &
          ," required for TOVS processing")

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

200  CONTINUE

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
    !    .  the observation errors and the utilization flags
500  CONTINUE

    DO JL = 1, INUMSAT
      READ (ILUTOV,*,ERR=900)
      READ (ILUTOV,'(A)',ERR=900) CLDUM
      CSATSTR = TRIM(ADJUSTL(CLDUM))

      !        Get satellite (e.g. NPP) from satellite/instrument (e.g. NPP ATMS)
      INDX = INDEX(CSATSTR,'ATMS')
      IF ( INDX .GT. 3 ) THEN
        IPOS = INDX-2
        CSATID(JL) = CSATSTR(1:IPOS)
      ELSE
        WRITE ( NULOUT, '(" mwbg_readStatTovsAtms: Non-ATMS ", &
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
        WRITE(NULOUT,'(/7X,"Satellite: ",A,5X,A)')  &
          CSATID(JL), CTYPSTAT(JK)
        WRITE(NULOUT,'(7X,"Channels   : ",30(T22,27I4/))') &
         (MLISCHNA(JI,JL),JI=1,NCHNA(JL))
        WRITE(NULOUT,'(7X,"Total errors: ",30(T22,27f4.1/))') &
         (TOVERRIN(MLISCHNA(JI,JL),JK,JL), &
          JI=1,NCHNA(JL))
      ENDDO
    ENDDO
     
700  CONTINUE

    RETURN

    ! Read error
900   WRITE ( NULOUT, '(" mwbg_readStatTovsAtms: Problem reading ", &
                       "TOVS total error stats file ")' ) 
    CALL ABORT ()

    RETURN
  END


  subroutine mwbg_getData(reportIndex, rpt, ISAT, zenith, ilq, itt, zlat, zlon, ztb, &
                          biasCorr, ZOMP, scanpos, nvalOut, ntOut, qcflag1, qcflag2, &
                          ican, icanomp, IMARQ, IORBIT)
    !--------------------------------------------------------------------------------------
    ! Object:   This routine extracts the needed data from the blocks in the report:
    !             reportIndex      = report index
    !             rpt              = report
    !             ISAT(nt)         = satellite identifier
    !             zenith(nt)       = satellite zenith angle (btyp=3072,ele=7024)
    !             ilq(nt)          = land/sea qualifier     (btyp=3072,ele=8012)
    !             itt(nt)          = terrain-type (ice)     (btyp=3072,ele=13039)
    !             zlat(nt)                                  (btyp=5120,ele=5002)
    !             zlon(nt)                                  (btyp=5120,ele=6002)
    !             ztb(nt*nval),    = brightness temperature (btyp=9248/9264,ele=12163) 
    !             biasCorr(nt*nval) = bias correction 
    !             ZOMP(nt*nval)    = OMP values
    !             scanpos(nt)      = scan position (fov)    (btyp=3072,ele=5043) 
    !             nvalOut          = number of channels     (btyp=9248/9264)
    !             ntOut            = number of locations    (btyp=5120,etc.)
    !             qcflag1(nt,3)    = flag values for btyp=3072 block ele 033078, 033079, 033080
    !             qcflag2(nt*nval) = flag values for btyp=9248 block ele 033081
    !             ican(nt*nval)    = channel numbers btyp=9248 block ele 5042 (= 1-22)
    !             icanomp(nt*nval) = omp channel numbers btyp= block ele  (= 1-22)
    !             IMARQ(nt*nval)   = data flags
    !             IORBIT(nt*nval)  = orbit number

    ! NOTE:  reportIndex = report number (from MAIN program) **** DO NOT MODIFY ****
    !       kk = variable for loops over locations (nt)
    !        j = variable for loops over nval (nval = 1 or nchanAtms)

    integer                     :: reportIndex
    type(BURP_RPT)              :: rpt
    integer, dimension(:)       :: scanpos
    real, dimension(:)          :: ztb, biasCorr, ZOMP
    real, dimension(:)          :: zlat,zlon,zenith
    integer, dimension(:)       :: ilq,itt,ISAT, IMARQ
    integer, dimension(:)       :: ican,icanomp,qcflag2,IORBIT
    integer, dimension(:,:)     :: qcflag1
    integer :: nvalOut, ntOut

    type(BURP_BLOCK)       :: blk

    integer :: error, ref_blk, my_nt,  my_nval, my_nele, my_idtyp
    integer :: indice, indice1, indice2, indice3, indice4, kk, j
    integer :: ind_lat,ind_lon

    integer :: ipos

    ! 1) Get the lat,lon from time/location block    BTYP = 5120  (also get nt)
    ref_blk = 0
    ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 5120, &
                    IOSTAT      = error)
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      write(*,*) 'ERREUR -  Location/time (3D) block (btyp=5120) not found in report number ', reportIndex
      call abort()
    endif

    call BURP_Get_Property(blk, &
                     NELE = my_nele, &
                     NT   = my_nt, &    ! number of locations in the box (report)
                     NVAL = my_nval, IOSTAT=error)
    if (error /= burp_noerr)  call abort()
   
    ind_lat = BURP_Find_Element(blk,5002,IOSTAT=error)
    ind_lon = BURP_Find_Element(blk,6002,IOSTAT=error)
    
    if ( (ind_lat > 0) .AND. (ind_lon > 0) ) then
      j = 1
      do kk =1, my_nt
        zlat(kk) = BURP_Get_Rval(blk,ind_lat,j,kk,error)
        zlon(kk) = BURP_Get_Rval(blk,ind_lon,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - lat, lon elements (5002,6002) missing in 3D block (btyp=5120). Report = ', reportIndex
      call abort()
    endif

    ! 2) Get info elements from the INFO block   BTYP = 3072
    ref_blk = 0
    ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 3072, &
                    IOSTAT      = error)
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      write(*,*) 'ERREUR - INFO block (btyp=3072) not found in report number ', reportIndex
      call abort()
    endif

    call BURP_Get_Property(blk, &
                      NELE = my_nele, &
                      NT   = my_nt, &    ! number of locations in the box (report)
                      NVAL = my_nval, IOSTAT=error)
    if (error /= burp_noerr)  call abort()

    indice = BURP_Find_Element(blk,1007,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        ISAT(kk) = BURP_Get_Tblval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - SATELLITE IDENTIFIER MISSING (ele=1007). Report = ', reportIndex
      call abort()
    endif

    indice = BURP_Find_Element(blk,5040,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        IORBIT(kk) = BURP_Get_Tblval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - ORBIT NUMBER MISSING (ele=5040). Report = ', reportIndex
      call abort()
    endif

    indice = BURP_Find_Element(blk,7024,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        zenith(kk) = BURP_Get_Rval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - satellite zenith angle missing in INFO block (ele=7024). Report = ', reportIndex
      call abort()
    endif

    indice = BURP_Find_Element(blk,8012,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        ilq(kk) = BURP_Get_Tblval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - land/sea qualifier missing in INFO block (ele=8012). Report = ', reportIndex
      call abort()
    endif   

    indice = BURP_Find_Element(blk,13039,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        itt(kk) = BURP_Get_Tblval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - terrain-type missing in INFO block (ele=13039). Report = ', reportIndex
      call abort()
    endif   

    indice = BURP_Find_Element(blk,5043,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        scanpos(kk) = BURP_Get_Tblval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - scan position missing in INFO block (ele=5043). Report = ', reportIndex
      call abort()
    endif   

    indice = BURP_Find_Element(blk,33078,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        qcflag1(kk,1) = BURP_Get_Tblval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - Geolocation quality code missing in INFO block (ele=33078). Report = ', reportIndex
      call abort()
    endif  

    indice = BURP_Find_Element(blk,33079,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        qcflag1(kk,2) = BURP_Get_Tblval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - Granule level QC flag missing in INFO block (ele=33079). Report = ', reportIndex
      call abort()
    endif  

    indice = BURP_Find_Element(blk,33080,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        qcflag1(kk,3) = BURP_Get_Tblval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - Scan level QC flag missing in INFO block (ele=33080). Report = ', reportIndex
      call abort()
    endif  

    ! 3) Get data from the DATA block     BTYP = 9248 or 9264    (also get nval = nchanAtms)
    ref_blk = 0
    ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 9248, &
                    IOSTAT      = error)
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      ref_blk = 0
      ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 9264, &
                    IOSTAT      = error)
      if (ref_blk < 0) then
        write(*,*) 'ERREUR - DATA block (btyp 9248 or 9264) not found in report number ', reportIndex
        call abort()
      endif
    endif

    call BURP_Get_Property(blk, &
                      NELE = my_nele, &
                      NT   = my_nt, &    ! number of locations in the box (report)
                      NVAL = my_nval, IOSTAT=error)
    if (error /= burp_noerr)  call abort()

    nvalOut = my_nval    ! set nvalOut (#channels) for MAIN program
    ntOut = my_nt

    indice1 = BURP_Find_Element(blk,12163,IOSTAT=error)
    indice2 = BURP_Find_Element(blk,33081,IOSTAT=error)
    indice3 = BURP_Find_Element(blk,5042,IOSTAT=error)
    indice4 = BURP_Find_Element(blk,12233,IOSTAT=error)
   
    if ( indice1 > 0 .and. indice2 > 0 .and. indice3 > 0 .and. indice4 > 0 ) then
      ipos = 0
      do kk = 1, my_nt
        do j = 1, my_nval
          ipos = ipos + 1
          ztb(ipos)     = BURP_Get_Rval(blk,indice1,j,kk,error)
          qcflag2(ipos) = BURP_Get_Tblval(blk,indice2,j,kk,error)
          ican(ipos)    = BURP_Get_Tblval(blk,indice3,j,kk,error)
          biasCorr(ipos)= BURP_Get_Rval(blk,indice4,j,kk,error)
        end do
      end do
    else
      write(*,*) 'ERREUR - Elements are missing in DATA block. Report = ', reportIndex
      if ( indice1 <= 0 )  write(*,*) '       Tb data             (012163) are missing!'
      if ( indice2 <= 0 )  write(*,*) '       Data level QC flags (033081) are missing!'
      if ( indice3 <= 0 )  write(*,*) '       Channel numbers     (005042) are missing!'
      if ( indice4 <= 0 )  write(*,*) '       Bias correction     (012233) are missing!'
      call abort()
    endif 

    ! 4) Get OMP data from the DATA block     BTYP = 9248 or 9264 and bfma = 14
    ref_blk = 0
    ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 14, &
                    BTYP        = 9322, &
                    IOSTAT      = error)
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      ref_blk = 0
      ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 14, &
                    BTYP        = 9226, &
                    IOSTAT      = error)
    end if
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      ref_blk = 0
      ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 14, &
                    BTYP        = 9258, &
                    IOSTAT      = error)
    end if
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      ref_blk = 0
      ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 14, &
                    BTYP        = 9274, &
                    IOSTAT      = error)
      if (ref_blk < 0) then
        write(*,*) 'ERREUR - OMP DATA block (btyp 9322 or 9226 or 9258 or 9274) not found in report number ', reportIndex
        call abort()
      endif
    endif

    call BURP_Get_Property(blk, &
                      NELE = my_nele, &
                      NT   = my_nt, &    ! number of locations in the box (report)
                      NVAL = my_nval, IOSTAT=error)
    if (error /= burp_noerr)  call abort()

    indice1 = BURP_Find_Element(blk,5042,IOSTAT=error)
    indice2 = BURP_Find_Element(blk,12163,IOSTAT=error)
    if ( indice1 > 0 .and. indice2 > 0 ) then
      ipos = 0
      do kk = 1, my_nt
        do j = 1, my_nval
          ipos = ipos + 1
          icanomp(ipos) = BURP_Get_Tblval(blk,indice1,j,kk,error)
          ZOMP(ipos)    = BURP_Get_Rval(blk,indice2,j,kk,error)
        end do
      end do
    else
      write(*,*) 'ERREUR - Elements are missing in OMP DATA block. Report = ', reportIndex
      if ( indice1 <= 0 )  write(*,*) '       Channel numbers     (005042) are missing!'
      if ( indice2 <= 0 )  write(*,*) '       Tb data             (012163) are missing!'
      call abort()
    endif 

    ! 5) Bloc marqueurs multi niveaux de radiances: bloc 15362, 15392, 15408.
    ref_blk = 0
    ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 15362, &
                    IOSTAT      = error)
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      ref_blk = 0
      ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 15392, &
                    IOSTAT      = error)
    end if
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      ref_blk = 0
      ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 15408, &
                    IOSTAT      = error)
      if (ref_blk < 0) then
        write(*,*) 'ERREUR - flag block (btyp 15362 or 15392 or 15408) not found in report number ', reportIndex
        call abort()
      endif
    end if

    call BURP_Get_Property(blk, &
                      NELE = my_nele, &
                      NT   = my_nt, &    ! number of locations in the box (report)
                      NVAL = my_nval, IOSTAT=error)
    if (error /= burp_noerr)  call abort()

    indice = BURP_Find_Element(blk,212163,IOSTAT=error)
   
    if ( indice > 0 ) then
      ipos = 0
      do kk = 1, my_nt
        do j = 1, my_nval
          ipos = ipos + 1
          IMARQ(ipos) = BURP_Get_Tblval(blk,indice,j,kk,error)
        end do
      end do
    else
      write(*,*) 'ERREUR - Element 212163 missing in flag block. Report = ', reportIndex
      call abort()
    endif 

    return

  end subroutine mwbg_getData


  subroutine mwbg_getDataAmsua(reportIndex, rpt, ISAT, zenith, ilq, itt, zlat, zlon, ztb, &
                          biasCorr, ZOMP, scanpos, nvalOut, ntOut, qcflag2, &
                          ican, icanomp, IMARQ, IORBIT, bad_report)
    !--------------------------------------------------------------------------------------
    ! Object:   This routine extracts the needed data from the blocks in the report:
    !             reportIndex      = report index
    !             rpt              = report
    !             ISAT(nt)         = satellite identifier
    !             zenith(nt)       = satellite zenith angle (btyp=3072,ele=7024)
    !             ilq(nt)          = land/sea qualifier     (btyp=3072,ele=8012)
    !             itt(nt)          = terrain-type (ice)     (btyp=3072,ele=13039)
    !             zlat(nt)                                  (btyp=5120,ele=5002)
    !             zlon(nt)                                  (btyp=5120,ele=6002)
    !             ztb(nt*nval),    = brightness temperature (btyp=9248/9264,ele=12163) 
    !             biasCorr(nt*nval) = bias correction 
    !             ZOMP(nt*nval)    = OMP values
    !             scanpos(nt)      = scan position (fov)    (btyp=3072,ele=5043) 
    !             nvalOut          = number of channels     (btyp=9248/9264)
    !             ntOut            = number of locations    (btyp=5120,etc.)
    !             qcflag2(nt*nval) = flag values for btyp=9248 block ele 033081
    !             ican(nt*nval)    = channel numbers btyp=9248 block ele 5042 (= 1-22)
    !             icanomp(nt*nval) = omp channel numbers btyp= block ele  (= 1-22)
    !             IMARQ(nt*nval)   = data flags
    !             IORBIT(nt*nval)  = orbit number

    ! NOTE:  reportIndex = report number (from MAIN program) **** DO NOT MODIFY ****
    !       kk = variable for loops over locations (nt)
    !        j = variable for loops over nval (nval = 1 or nchanAtms)

    integer                     :: reportIndex
    type(BURP_RPT)              :: rpt
    integer, dimension(:)       :: scanpos
    real, dimension(:)          :: ztb, biasCorr, ZOMP
    real, dimension(:)          :: zlat,zlon,zenith
    integer, dimension(:)       :: ilq,itt,ISAT, IMARQ
    integer, dimension(:)       :: ican,icanomp,qcflag2,IORBIT
    integer :: nvalOut, ntOut
    logical, intent(out) :: bad_report

    type(BURP_BLOCK)       :: blk

    integer :: error, ref_blk, my_nt,  my_nval, my_nele, my_idtyp
    integer :: indice, indice1, indice2, indice3, indice4, kk, j
    integer :: ind_lat,ind_lon

    integer :: ipos


    zenith(:) = ZMISG
    itt(:) = MISGINT
    bad_report = .false.

    ! 1) Get the lat,lon from time/location block    BTYP = 5120  (also get nt)
    ref_blk = 0
    ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 5120, &
                    IOSTAT      = error)
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      write(*,*) 'ERREUR -  Location/time (3D) block (btyp=5120) not found in report number ', reportIndex
      call abort()
    endif

    call BURP_Get_Property(blk, &
                     NELE = my_nele, &
                     NT   = my_nt, &    ! number of locations in the box (report)
                     NVAL = my_nval, IOSTAT=error)
    if (error /= burp_noerr)  call abort()

    do kk =1, my_nt
    end do
   
    ind_lat = BURP_Find_Element(blk,5002,IOSTAT=error)
    ind_lon = BURP_Find_Element(blk,6002,IOSTAT=error)
    
    if ( (ind_lat > 0) .AND. (ind_lon > 0) ) then
      j = 1
      do kk =1, my_nt
        zlat(kk) = BURP_Get_Rval(blk,ind_lat,j,kk,error)
        zlon(kk) = BURP_Get_Rval(blk,ind_lon,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - lat, lon elements (5002,6002) missing in 3D block (btyp=5120). Report = ', reportIndex
      call abort()
    endif

    ! 2) Get info elements from the INFO block   BTYP = 3072
    ref_blk = 0
    ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 3072, &
                    IOSTAT      = error)
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      write(*,*) 'ERREUR - INFO block (btyp=3072) not found in report number ', reportIndex
      call abort()
    endif

    call BURP_Get_Property(blk, &
                      NELE = my_nele, &
                      NT   = my_nt, &    ! number of locations in the box (report)
                      NVAL = my_nval, IOSTAT=error)
    if (error /= burp_noerr)  call abort()

    indice = BURP_Find_Element(blk,1007,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        ISAT(kk) = BURP_Get_Tblval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - SATELLITE IDENTIFIER MISSING (ele=1007). Report = ', reportIndex
      call abort()
    endif

    indice = BURP_Find_Element(blk,5040,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        IORBIT(kk) = BURP_Get_Tblval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - ORBIT NUMBER MISSING (ele=5040). Report = ', reportIndex
      call abort()
    endif

    indice = BURP_Find_Element(blk,7024,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        zenith(kk) = BURP_Get_Rval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - satellite zenith angle missing in INFO block (ele=7024). Report = ', reportIndex
      call abort()
    endif

    indice = BURP_Find_Element(blk,8012,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        ilq(kk) = BURP_Get_Tblval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - land/sea qualifier missing in INFO block (ele=8012). Report = ', reportIndex
      call abort()
    endif   

    indice = BURP_Find_Element(blk,13039,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        itt(kk) = BURP_Get_Tblval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - terrain-type missing in INFO block (ele=13039). Report = ', reportIndex
      call abort()
    endif   

    indice = BURP_Find_Element(blk,5043,IOSTAT=error)
    if ( indice > 0 ) then
      j = 1
      do kk =1, my_nt
        scanpos(kk) = BURP_Get_Tblval(blk,indice,j,kk,error)
      end do
    else
      write(*,*) 'ERREUR - scan position missing in INFO block (ele=5043). Report = ', reportIndex
      call abort()
    endif   

    ! 3) Get data from the DATA block     BTYP = 9248 or 9264    (also get nval = nchanAtms)
    ref_blk = 0
    ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 9248, &
                    IOSTAT      = error)
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      ref_blk = 0
      ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 9264, &
                    IOSTAT      = error)
      if (ref_blk < 0) then
        write(*,*) 'ERREUR - DATA block (btyp 9248 or 9264) not found in report number ', reportIndex
        call abort()
      endif
    endif

    call BURP_Get_Property(blk, &
                      NELE = my_nele, &
                      NT   = my_nt, &    ! number of locations in the box (report)
                      NVAL = my_nval, IOSTAT=error)
    if (error /= burp_noerr)  call abort()

    nvalOut = my_nval    ! set nvalOut (#channels) for MAIN program
    ntOut = my_nt

    indice = BURP_Find_Element(blk,2150,IOSTAT=error)
    if ( indice > 0 ) then
      ipos = 0
      do kk = 1, my_nt
        do j = 1, my_nval
          ipos = ipos + 1
          ican(ipos) = BURP_Get_Tblval(blk,indice,j,kk,error)
        end do
      end do
    else
      write(*,*) ' ERREUR - Elements Channel numbers (002150) are missing in DATA block! Report = ', reportIndex
      call abort()
    end if

    indice = BURP_Find_Element(blk,12163,IOSTAT=error)
    if ( indice > 0 ) then
      ipos = 0
      do kk = 1, my_nt
        do j = 1, my_nval
          ipos = ipos + 1
          ztb(ipos) = BURP_Get_Rval(blk,indice,j,kk,error)
        end do
      end do
    else
      write(*,*) ' ERREUR - Elements Tb data (012163) are missing in DATA block! Report = ', reportIndex
      call abort()
    end if

    indice = BURP_Find_Element(blk,33032,IOSTAT=error)
    if ( indice > 0 ) then
      ipos = 0
      do kk = 1, my_nt
        do j = 1, my_nval
          ipos = ipos + 1
          qcflag2(ipos) = BURP_Get_Tblval(blk,indice,j,kk,error)
        end do
      end do
    else
      write(*,*) ' ERREUR - Elements Data level QC flags (033032) are missing in DATA block! Report = ', reportIndex
      call abort()
    end if

    indice = BURP_Find_Element(blk,12233,IOSTAT=error)
    if ( indice > 0 ) then
      ipos = 0
      do kk = 1, my_nt
        do j = 1, my_nval
          ipos = ipos + 1
          biasCorr(ipos)= BURP_Get_Rval(blk,indice,j,kk,error)
        end do
      end do
    else
      if ( mwbg_debug ) write(*,*) ' ERREUR - Elements Bias correction (012233) are missing in DATA block! Report = ', reportIndex
      biasCorr(:) = zmisg
    end if

    ! 5) Bloc marqueurs multi niveaux de radiances: bloc 15362, 15392, 15408.
    ref_blk = 0
    ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 15362, &
                    IOSTAT      = error)
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      ref_blk = 0
      ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 15392, &
                    IOSTAT      = error)
    end if
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      ref_blk = 0
      ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 0, &
                    BTYP        = 15408, &
                    IOSTAT      = error)
      if (ref_blk < 0) then
        write(*,*) 'ERREUR - flag block (btyp 15362 or 15392 or 15408) not found in report number ', reportIndex
        call abort()
      endif
    end if

    call BURP_Get_Property(blk, &
                      NELE = my_nele, &
                      NT   = my_nt, &    ! number of locations in the box (report)
                      NVAL = my_nval, IOSTAT=error)
    if (error /= burp_noerr)  call abort()

    indice = BURP_Find_Element(blk,212163,IOSTAT=error)
   
    if ( indice > 0 ) then
      ipos = 0
      do kk = 1, my_nt
        do j = 1, my_nval
          ipos = ipos + 1
          IMARQ(ipos) = BURP_Get_Tblval(blk,indice,j,kk,error)
        end do
      end do
    else
      write(*,*) 'ERREUR - Element 212163 missing in flag block. Report = ', reportIndex
      call abort()
    endif 

    ! 4) Get OMP data from the DATA block     BTYP = 9248 or 9264 and bfma = 14
    ref_blk = 0
    ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 14, &
                    BTYP        = 9322, &
                    IOSTAT      = error)
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      ref_blk = 0
      ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 14, &
                    BTYP        = 9226, &
                    IOSTAT      = error)
    end if
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      ref_blk = 0
      ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 14, &
                    BTYP        = 9258, &
                    IOSTAT      = error)
    end if
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      ref_blk = 0
      ref_blk = BURP_Find_Block(rpt, &
                    BLOCK       = blk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = 14, &
                    BTYP        = 9274, &
                    IOSTAT      = error)
    end if
    if (error /= burp_noerr) call abort()
    if (ref_blk < 0) then
      write(*,*) 'ERREUR - OMP DATA block (btyp 9322 or 9226 or 9258 or 9274) not found in report number ', reportIndex
      bad_report = .true.
      return
    endif

    call BURP_Get_Property(blk, &
                      NELE = my_nele, &
                      NT   = my_nt, &    ! number of locations in the box (report)
                      NVAL = my_nval, IOSTAT=error)
    if (error /= burp_noerr)  call abort()

    indice1 = BURP_Find_Element(blk,2150,IOSTAT=error)

    indice2 = BURP_Find_Element(blk,12163,IOSTAT=error)

    if ( indice1 > 0 .and. indice2 > 0 ) then
      ipos = 0
      do kk = 1, my_nt
        do j = 1, my_nval
          ipos = ipos + 1
          icanomp(ipos) = BURP_Get_Tblval(blk,indice1,j,kk,error)
          ZOMP(ipos)    = BURP_Get_Rval(blk,indice2,j,kk,error)
        end do
      end do
    else
      write(*,*) 'ERREUR - Elements are missing in OMP DATA block. Report = ', reportIndex
      if ( indice1 <= 0 )  write(*,*) '       Channel numbers     (005042) are missing!'
      if ( indice2 <= 0 )  write(*,*) '       Tb data             (012163) are missing!'
      call abort()
    endif 

    return

  end subroutine mwbg_getDataAmsua


  subroutine mwbg_writeBlocks(reportIndex, ztb, lsq, trn, riwv, rclw, ident, &
                              logicalFlags, IMARQ, lutb, rpt, rpt_out)
    ! Object:   This routine modifies the blocks in the input Report (rpt) 
    integer :: reportIndex
    real,    intent(in), dimension(:)   :: ztb
    integer, intent(in), dimension(:)   :: lsq
    integer, intent(in), dimension(:)   :: trn
    real,    intent(in), dimension(:)   :: riwv
    real,    intent(in), dimension(:)   :: rclw
    integer, intent(in), dimension(:)   :: ident
    logical, intent(in), dimension(:,:) :: logicalFlags
    integer, intent(inout), dimension(:):: IMARQ
    logical :: lutb
    type(BURP_RPT)         :: rpt
    type(BURP_RPT)         :: rpt_out

    type(BURP_BLOCK)       :: blk, blk_copy

    integer :: error, ref_blk, my_nt,  my_nval, my_nele, my_btyp, my_bfam, iidata
    integer :: indice, indice1, indice2, j, kk, ipos

    Call BURP_Init(blk, B2=blk_copy, IOSTAT=error)

    ! 1) Read and modify the blocks in rpt and add them to rpt_out
    ref_blk = 0
    BLOCKS: do

      ref_blk = BURP_Find_Block(rpt,BLOCK= blk,SEARCH_FROM= ref_blk,IOSTAT= error)
      if (error /= burp_noerr) call abort()
      
      if (ref_blk < 0) Exit

      Call BURP_Get_Property(blk, &
                  NELE   = my_nele, &
                  NVAL   = my_nval, &       ! 1 or number of channels (obs per location) if Tb data/flag block
                  NT     = my_nt, &         ! 1 or number of locations in block
                  BTYP   = my_btyp, &
                  BFAM   = my_bfam, &
                  IOSTAT = error)
      if (error /= burp_noerr) call abort()

      ! 3D Block
      if (my_btyp == 5120) then     
        ! Set bit 6 in 24-bit global flags if any data rejected

        ! Extract the global flags, element 55200
        indice = BURP_Find_Element(blk,55200,IOSTAT=error)
        if ( indice > 0 ) then
          j = 1
          do kk =1, my_nt
            iidata = BURP_Get_Tblval(blk,indice,j,kk,error)
            if ( ANY(logicalFlags(kk,:)) ) iidata = IBSET(iidata,6)
            Call BURP_Set_Tblval(blk,indice,j,kk,iidata)
          end do
        else
          write(*,*) 'ERREUR - Global flag missing in 3D block (ele=55200). Report = ', reportIndex
          call abort()
        endif 
          
        blk_copy = blk
        Call BURP_Write_Block(rpt_out,BLOCK=blk_copy,CONVERT_BLOCK=.FALSE.,IOSTAT=error)
        if (error /= burp_noerr)  call abort()

      ! INFO block (mix of integer and real data)
      elseif (my_btyp == 3072) then


        ! Add new elements QC indent flag (ident), CLW (rclw) and ECMWF_SI (riwv)
        !   OPTION: replace land-sea qualifier (lsq) and terrain type (trn) with internal values
        !   NOTE: after setting real values (Rval), call BURP_Convert_Block()!!

        Call BURP_Resize_Block(blk, ADD_NELE = 3, IOSTAT = error)
        if (error /= burp_noerr)  call abort()
        
        Call BURP_Set_Element(blk, NELE_IND = my_nele+1, ELEMENT = 25174, IOSTAT = error)
        Call BURP_Set_Element(blk, NELE_IND = my_nele+2, ELEMENT = 13209, IOSTAT = error)
        Call BURP_Set_Element(blk, NELE_IND = my_nele+3, ELEMENT = 13208, IOSTAT = error)
        Call BURP_Encode_Block(blk)   ! encode the element numbers in the block

        j = 1
        do kk =1, my_nt
          iidata = ident(kk)
          Call BURP_Set_Rval  (blk, NELE_IND=my_nele+2,NVAL_IND=j,NT_IND=kk, RVAL=rclw(kk),IOSTAT=error)
          Call BURP_Set_Rval  (blk, NELE_IND=my_nele+3,NVAL_IND=j,NT_IND=kk, RVAL=riwv(kk),IOSTAT=error)
        end do
      
        Call BURP_Convert_Block(blk)

        j = 1
        do kk =1, my_nt
          iidata = ident(kk)
          Call BURP_Set_Tblval(blk, NELE_IND=my_nele+1,NVAL_IND=j,NT_IND=kk, TBLVAL=iidata,IOSTAT= error)
        end do

        if (mwbg_modlsqtt) then
          indice1 = BURP_Find_Element(blk,  8012, IOSTAT=error)
          if (error /= burp_noerr)  call abort()
          indice2 = BURP_Find_Element(blk, 13039, IOSTAT=error)
          if (error /= burp_noerr)  call abort()
          if ( indice1 > 0 .and. indice2 > 0 ) then
            j = 1
            do kk =1, my_nt
              iidata = lsq(kk)
              Call BURP_Set_Tblval(blk,indice1,j,kk,iidata,error)
              if (error /= burp_noerr)  call abort()
              iidata = trn(kk)
              Call BURP_Set_Tblval(blk,indice2,j,kk,iidata,error)
              if (error /= burp_noerr)  call abort()
            end do
          else
            write(*,*) 'ERREUR - land/sea qualifier (ele=8012) and/or terrain type (ele=13039) not found in INFO block. Report = ', reportIndex
            call abort()
          endif                     
        endif

        blk_copy = blk
        Call BURP_Write_Block(rpt_out,BLOCK=blk_copy,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
        if (error /= burp_noerr)  call abort()      

      !  DATA block
      elseif (my_btyp == 9248 .or. my_btyp ==9264) then 
        ! Modify Tb data if any data (ztb) were set to zmisg (lutb=.true.)

        if (lutb) then
          indice = BURP_Find_Element(blk, 12163, IOSTAT=error)
          if (error /= burp_noerr)  call abort()
          if ( indice > 0 ) then
            ipos = 0
            do kk =1, my_nt
              do j = 1, my_nval
                ipos = ipos + 1
                Call BURP_Set_Rval(blk,NELE_IND=indice,NVAL_IND=j,NT_IND=kk,RVAL=ztb(ipos),IOSTAT=error)  
                if (error /= burp_noerr)  call abort()
              enddo
            enddo
          else
            write(*,*) 'ERREUR - Cannot find Tb (ele=12163) in DATA block!. Report = ', reportIndex
            call abort()
          endif
!          if ( mwbg_debug ) Call BURP_TO_STDOUT(blk, CONVERT =.FALSE.)
          blk_copy = blk
          Call BURP_Write_Block(rpt_out,BLOCK=blk_copy,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
          if (error /= burp_noerr)  call abort() 
        else
          Call BURP_Write_Block(rpt_out,BLOCK=blk,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
          if (error /= burp_noerr)  call abort() 
        endif

      ! FLAG block
      elseif (my_btyp == 15392 .or. my_btyp == 15408) then 
        ! Modify data flag values (set bit 7) for rejected data    

        indice = BURP_Find_Element(blk, 212163, IOSTAT=error)
        if (error /= burp_noerr)  call abort()
        if ( indice > 0 ) then 
          ipos = 0
          do kk =1, my_nt
            do j = 1, my_nval
              ipos = ipos + 1
              iidata = BURP_Get_Tblval(blk,indice,j,kk,error)

              if (logicalFlags(kk,j)) then
                iidata = IBSET(iidata,7)
                IMARQ(ipos) = iidata
              end if

              Call BURP_Set_Tblval(blk,indice,j,kk,iidata)
            enddo
          enddo
        else
          write(*,*) 'ERREUR - Data QC flags (ele=212163) not found in FLAG block. Report = ', reportIndex
          call abort()      
        endif

        blk_copy = blk
        Call BURP_Write_Block(rpt_out,BLOCK=blk_copy,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
        if (error /= burp_noerr)  call abort()

      ! OTHER BLOCK 
      else 
        Call BURP_Write_Block(rpt_out,BLOCK=blk,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
        if (error /= burp_noerr)  call abort() 
        
      endif      

    enddo BLOCKS

    return

  end subroutine mwbg_writeBlocks


  subroutine mwbg_writeBlocksAmsua(clw, scatw, reportIndex, rpt, rpt_out)
    ! Object:   This routine modifies the blocks in the input Report (rpt) 
    real, dimension(:)     :: clw, scatw
    integer                :: reportIndex
    type(BURP_RPT)         :: rpt
    type(BURP_RPT)         :: rpt_out

    type(BURP_BLOCK)       :: blk, blk_copy

    integer :: error, ref_blk, my_nt,  my_nval, my_nele, my_btyp, my_bfam, iidata
    integer :: indice, indice1, indice2, j, kk, ipos

    Call BURP_Init(blk, B2=blk_copy, IOSTAT=error)

    ! 1) Read and modify the blocks in rpt and add them to rpt_out
    ref_blk = 0
    BLOCKS: do

      ref_blk = BURP_Find_Block(rpt,BLOCK= blk,SEARCH_FROM= ref_blk,IOSTAT= error)
      if (error /= burp_noerr) call abort()
      
      if (ref_blk < 0) Exit

      Call BURP_Get_Property(blk, &
                  NELE   = my_nele, &
                  NVAL   = my_nval, &       ! 1 or number of channels (obs per location) if Tb data/flag block
                  NT     = my_nt, &         ! 1 or number of locations in block
                  BTYP   = my_btyp, &
                  BFAM   = my_bfam, &
                  IOSTAT = error)
      if (error /= burp_noerr) call abort()

      ! 3D Block
      if (my_btyp == 5120) then     
        blk_copy = blk
        Call BURP_Write_Block(rpt_out,BLOCK=blk_copy,CONVERT_BLOCK=.FALSE.,IOSTAT=error)
        if (error /= burp_noerr)  call abort()

      ! INFO block (mix of integer and real data)
      elseif (my_btyp == 3072) then
        ! Add new elements CLW and scatw
        Call BURP_Resize_Block(blk, ADD_NELE = 2, IOSTAT = error)
        if (error /= burp_noerr)  call abort()
        
        Call BURP_Set_Element(blk,NELE_IND=my_nele+1,ELEMENT=13209,IOSTAT=error)
        Call BURP_Set_Element(blk,NELE_IND=my_nele+2,ELEMENT=13208,IOSTAT=error)
        Call BURP_Encode_Block(blk)   ! encode the element numbers in the block

        j = 1
        do kk =1, my_nt
          Call BURP_Set_Rval(blk,NELE_IND=my_nele+1,NVAL_IND=j,NT_IND=kk,RVAL=clw(kk),IOSTAT=error)
          Call BURP_Set_Rval(blk,NELE_IND=my_nele+2,NVAL_IND=j,NT_IND=kk,RVAL=scatw(kk),IOSTAT=error)
        end do
        Call BURP_Convert_Block(blk)

        blk_copy = blk
        Call BURP_Write_Block(rpt_out,BLOCK=blk_copy,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
        if (error /= burp_noerr)  call abort()      

      !  DATA block
      elseif (my_btyp == 9248 .or. my_btyp ==9264) then 

        blk_copy = blk
        Call BURP_Write_Block(rpt_out,BLOCK=blk_copy,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
        if (error /= burp_noerr)  call abort() 

      ! FLAG block
      elseif (my_btyp == 15392 .or. my_btyp == 15408) then 

        blk_copy = blk
        Call BURP_Write_Block(rpt_out,BLOCK=blk_copy,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
        if (error /= burp_noerr)  call abort()

      ! OTHER BLOCK 
      else 
        blk_copy = blk
        Call BURP_Write_Block(rpt_out,BLOCK=blk_copy,CONVERT_BLOCK=.TRUE.,IOSTAT=error)
        if (error /= burp_noerr)  call abort() 
        
      endif      

    enddo BLOCKS

    return

  end subroutine mwbg_writeBlocksAmsua


  subroutine mwbg_landIceMaskAtms(mglg_file,npts,zlat,zlon,zlq,ztt,waterobs)
    ! Adapted from: land_ice_mask_ssmis.ftn90 of mwbg_ssmis (D. Anselmo, S. Macpherson)
    !
    ! Object:   This routine sets waterobs array by performing a land/ice proximity check using
    !           using analysis MG and LG (or GL) fields used by the model which produces the trial field.
    !           The purpose of this check is to remove obs that reside close to coasts or ice,
    !           and so whose TBs may be contaminated.
    !           The GEM Global (glbhyb2) analysis contains MG and LG fields (on different grids).
    !
    !           NOTE: The 0.1 deg binary ice field check from land_ice_mask_ssmis.ftn90
    !           was removed. The land/sea qualifier (zlq) and terrain type (ztt) are modified
    !           to indicate proximity to land and sea-ice but are NOT changed in output BURP file.
    !
    !           In the application of this check, a 5x5 mesh, with spacing defined by rlat_km and
    !           rlon_km, is positioned with its center over an obs pt (2 grid pts on either side
    !           of the obs pt; size of mesh is equal to 4*rlat_km x 4*rlon_km). The values of MG
    !           and LG are evaluated at the grid points of this mesh. The maximum value of each
    !           determines whether the obs pt is too close to ice or land to be retained.
    !           **NOTE: the threshold value for MG has a very strong effect on the distance
    !                   from land that is permitted for an obs to be retained
    !
    !
    !      Maximum FOV             x---x---x---x---x     ^
    !         = 75km x 75km        |   |   |   |   |     |
    !         for Meso-sphere CHs  x---x---x---x---x     |
    !         = 74km x 47km        |   |   |   |   |     |
    !         for 19 GHz           x---x---o---x---x     | = 4*rlat_km
    !                              |   |   |   |   |     | = 4*40 km
    !                           ^  x---x---x---x---x     | = 160 km = 80 km north & south
    !                   rlat_km |  |   |   |   |   |     |
    !                           v  x---x---x---x---x     v
    !                                          <--->
    !                                         rlon_km
    !
    !                              <--------------->
    !                                 = 4*rlon_km
    !                                 = 4*40 km
    !                                 = 160 km = 80 km east & west
    !
    !
    !               MG value = 1.0  ==>  LAND       MG value = 0.0  ==>  OCEAN
    !               LG value = 1.0  ==>  ICE        LG value = 0.0  ==>  NO ICE
    !
    !
    ! Version:      Date:      Comment:
    ! --------      -----      --------
    !   0.1       16/08/12     Original adapted code.      S. Macpherson  
    !   0.2       01/03/14     Open mglg_file in R/O mode  S. Macpherson
    !
    !--------------------------------------------------------------------
    !  Variable Definitions
    !  --------------------
    ! mglg_file  - input  -  name of file holding model MG and LG (or GL) fields
    ! npts       - input  -  number of input obs pts in report
    ! zlat       - input  -  array holding lat values for all obs pts in report
    ! zlon       - input  -  array holding lon values for all obs pts in report
    ! zlq        - in/out -  array holding land/sea qualifier values for all obs
    !                        pts of report (0 = land, 1 = sea)
    ! ztt        - in/out -  array holding terrain-type values for all obs pts
    !                        of current report (-1 land/open water, 0 = ice)
    ! waterobs   - output -  logical array identifying for each obs in current report
    !                        whether it is over open water, far from coast/ice
    ! mxlat      -internal-  number of grid pts in lat. direction for mesh
    ! mxlon      -internal-  number of grid pts in lon. direction for mesh
    ! rlat_km    -internal-  spacing desired between mesh grid points in km
    !                        along lat. direction
    ! rlon_km    -internal-  spacing desired between mesh grid points in km
    !                        along lon. direction
    ! dlat       -internal-  spacing between mesh grid points along lon. direction
    !                        in degrees computed from rlat_km
    ! dlon       -internal-  spacing between mesh grid points along lon. direction
    !                        in degrees computed from rlon_km
    ! rkm_per_deg -internal- distance in km per degree
    !                           = Earth radius * PI/180.0
    !                           = 6371.01 km * PI/180.0
    !                           = 111.195 km
    ! nlat,nlon  -internal-  used to define the lat/lon of the grid pts of mesh
    ! zlatbox    -internal-  lat values at all grid pts of mesh for all obs pts
    ! zlonbox    -internal-  lon values at all grid pts of mesh for all obs pts
    ! latmesh    -internal-  lat values at all grid pts of mesh for 1 obs pt
    ! lonmesh    -internal-  lon values at all grid pts of mesh for 1 obs pt
    ! mgintob    -internal-  interpolated MG values at all grid pts of mesh for 1 obs pt
    ! lgintob    -internal-  interpolated LG values at all grid pts of mesh for 1 obs pt
    ! mgintrp    -internal-  max. interpolated MG value on mesh for all obs pts
    ! lgintrp    -internal-  max. interpolated LG value on mesh for all obs pts
    ! MGthresh   -internal-  maximum allowable land fraction for obs to be kept
    ! LGthresh   -internal-  maximum allowable ice  fraction for obs to be kept
    !--------------------------------------------------------------------
    implicit none

    ! Arguments:
    character(len=128), intent(in) :: mglg_file

    integer, intent(in)                   :: npts
    real,    intent(in),     dimension(:) :: zlat,zlon
    integer, intent(inout),  dimension(:) :: zlq, ztt

    logical, intent(out), dimension(:) :: waterobs

    ! Locals:
    integer, parameter :: mxlat=5,mxlon=5
    integer, parameter :: iungeo=50

    integer :: ier,key,istat
    integer :: ni,nj,nk,nilg,njlg
    integer :: ig1,ig2,ig3,ig4,ig1lg,ig2lg,ig3lg,ig4lg
    integer :: idum4,idum5,idum6,idum7,idum8,idum9,idum10,idum11
    integer :: idum12,idum13,idum14,idum15,idum16,idum17,idum18

    integer :: indx,ii,jj,kk
    integer :: nlat,nlon

    integer, dimension(10) :: alloc_status
  
    real, parameter :: pi=3.141592654
    real, parameter :: MGthresh=0.01,LGthresh=0.01
    real, parameter :: rlat_km=40.0,rlon_km=40.0
    real, parameter :: rkm_per_deg=111.195

    real :: xlat,xlatrad,xlon,rii,rjj
    real :: dlat,dlon

    character(len=12) :: etikxx
    character(len=4)  :: nomvxx
    character(len=2)  :: typxx
    character(len=1)  :: grtyp,grtyplg
  
    logical  :: llg

    ! F90 allocatable arrays:
    real, allocatable, dimension(:)   :: mg,lg
    real, allocatable, dimension(:)   :: latmesh,lonmesh
    real, allocatable, dimension(:)   :: mgintob,lgintob
    real, allocatable, dimension(:,:) :: zlatbox,zlonbox
    real, allocatable, dimension(:)   :: mgintrp,lgintrp
  
    ! RMNLIB interpolating functions:
    integer :: ezsetopt,ezqkdef
    integer :: gdllsval,gdid,gdidlg

    ! Define FORTRAN FST functions:
    integer, external :: fstinf,fstprm,fstlir
    integer, external :: fstouv,fstfrm,fstinl,fstvoi

    integer :: idum1,idum2,idum3

    ! Allocate space for arrays holding values on mesh grid pts.
    alloc_status(:) = 0
    allocate ( latmesh(mxlat*mxlon), stat=alloc_status(1) )
    allocate ( lonmesh(mxlat*mxlon), stat=alloc_status(2) )
    allocate ( mgintob(mxlat*mxlon), stat=alloc_status(3) )
    allocate ( lgintob(mxlat*mxlon), stat=alloc_status(4) )
    allocate ( zlatbox(mxlat*mxlon,npts), stat=alloc_status(5) )
    allocate ( zlonbox(mxlat*mxlon,npts), stat=alloc_status(6) )
    if( any(alloc_status /= 0) ) then
      write(*,*) 'mwbg_landIceMaskAtms: Memory allocation error '
      call abort()
    endif

    ! Open FST file.
    ier = fnom( iungeo,mglg_file,'STD+RND+R/O',0 )
    ier = fstouv( iungeo,'RND' )

    ! Read MG field.
    key = fstinf(iungeo,ni,nj,nk,-1,' ',-1,-1,-1,' ' ,'MG')
    if ( key <  0 ) then
      write(*,*) 'mwbg_landIceMaskAtms: The MG field is MISSING '
      call abort()
    end if

    allocate ( mg(ni*nj), stat=alloc_status(7) )
    if( any(alloc_status /= 0) ) then
      write(*,*) 'mwbg_landIceMaskAtms: Memory allocation error '
      call abort()
    endif
    ier = fstlir(mg,iungeo,ni,nj,nk,-1,' ',-1,-1,-1,' ','MG')

    ier = fstprm(key,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,    &
                idum9,idum10,idum11,typxx,nomvxx,etikxx,grtyp,ig1,ig2,  &
                ig3,ig4,idum12,idum13,idum14,idum15,idum16,idum17,      &
                idum18)


    ! Read LG field. Use GL field as backup.
    ! **CAUTION**: Discontinuities in GL field may cause interpolation problems! LG field is preferable.
    llg=.false.
    key = fstinf(iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ' ,'LG')
    if ( key <  0 ) then
      key = fstinf(iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ' ,'GL')
      if ( key <  0 ) then
        write(*,*) 'mwbg_landIceMaskAtms: No ice (LG or GL) fields found. Aborting! '
        call abort()
      else
        !write(*,*) 'mwbg_landIceMaskAtms: The GL field was found and will be used.'
      endif
    else
      llg=.true.
    endif

    allocate ( lg(nilg*njlg), stat=alloc_status(8) )
    if( any(alloc_status /= 0) ) then
      write(*,*) 'mwbg_landIceMaskAtms: Memory allocation error '
      call abort()
    endif
    
    if ( llg ) then
      ier = fstlir(lg,iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ','LG')
    else
      ier = fstlir(lg,iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ','GL')
    endif

    ier = fstprm(key,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,          &
                idum9,idum10,idum11,typxx,nomvxx,etikxx,grtyplg,ig1lg,ig2lg,  &
                ig3lg,ig4lg,idum12,idum13,idum14,idum15,idum16,idum17,        &
                idum18)

    ! For each obs pt, define a grid of artificial pts surrounding it.
    nlat = ( mxlat - 1 ) / 2
    nlon = ( mxlon - 1 ) / 2

    dlat = rlat_km / rkm_per_deg
    do kk = 1, npts
      indx = 0

      do ii = -nlat, nlat
        rii = float(ii)
        xlat = zlat(kk) + rii*dlat
        xlat = max( -90.0, min(90.0,xlat) )
        xlatrad = xlat*pi/180.0

        do jj = -nlon, nlon
          dlon = rlon_km / ( rkm_per_deg*cos(xlatrad) )
          rjj = float(jj)
          indx = indx + 1
          xlon = zlon(kk) + rjj*dlon
          if ( xlon < -180. ) xlon = xlon + 360.
          if ( xlon >  180. ) xlon = xlon - 360.
          if ( xlon <    0. ) xlon = xlon + 360.
          zlatbox(indx,kk) = xlat
          zlonbox(indx,kk) = xlon
        end do

      end do
    end do


    ! Interpolate values from MG and LG field to grid pts of mesh centred over each obs pt.
    ! Determine for each obs pt, the max interpolated MG and LG value within the box
    ! surrounding it.
    ier    = ezsetopt('INTERP_DEGREE','LINEAR')
    gdid   = ezqkdef(ni,nj,grtyp,ig1,ig2,ig3,ig4,iungeo)
    gdidlg = ezqkdef(nilg,njlg,grtyplg,ig1lg,ig2lg,ig3lg,ig4lg,iungeo)

    allocate ( mgintrp(npts), stat=alloc_status(9) )
    allocate ( lgintrp(npts), stat=alloc_status(10) )
    if( any(alloc_status /= 0) ) then
      write(*,*) 'mwbg_landIceMaskAtms: Memory allocation error '
      call abort()
    endif

    mgintrp(:) = 0.0
    lgintrp(:) = 0.0
    do kk = 1, npts

      latmesh = zlatbox(:,kk)
      lonmesh = zlonbox(:,kk)

      ier  = gdllsval(gdid,mgintob,mg,latmesh,lonmesh,mxlat*mxlon)
      ier  = gdllsval(gdidlg,lgintob,lg,latmesh,lonmesh,mxlat*mxlon)

      mgintrp(kk) = maxval(mgintob(:))
      lgintrp(kk) = maxval(lgintob(:))

    end do

    !  Initialize all obs as being over land and free of ice or snow.
    !  Determine which obs are over open water.
    waterobs(:) = .false.   ! not over open water
    ztt(:) = -1             ! no ice (reset terain type)
    zlq(:) = 0              ! land   (reset land/sea qualifier)

    do kk = 1, npts
      if ( mgintrp(kk) < MGthresh ) zlq(kk) = 1  ! ocean point away from coast
      if ( lgintrp(kk) >= LGthresh .and. zlq(kk) == 1 ) ztt(kk) = 0  ! sea-ice affected point
      if ( lgintrp(kk)  < LGthresh .and. zlq(kk) == 1 ) then
        waterobs(kk) = .true.  ! water point not in close proximity to land or sea-ice
      end if
    end do

    ! Deallocate arrays and close FST file.
    alloc_status(:) = 0
    deallocate ( mgintrp, stat=alloc_status(1) )
    deallocate ( lgintrp, stat=alloc_status(2) )
    deallocate ( mg,      stat=alloc_status(3) )
    deallocate ( lg,      stat=alloc_status(4) )
    deallocate ( latmesh, stat=alloc_status(5) )
    deallocate ( lonmesh, stat=alloc_status(6) )
    deallocate ( mgintob, stat=alloc_status(7) )
    deallocate ( lgintob, stat=alloc_status(8) )
    deallocate ( zlatbox, stat=alloc_status(9) )
    deallocate ( zlonbox, stat=alloc_status(10) )
    if( any(alloc_status /= 0) ) then
      write(*,*) 'mwbg_landIceMaskAtms: Memory deallocation error '
      call abort()
    endif
    ier = fstfrm(iungeo)
    ier = fclos(iungeo)

    return
  end subroutine mwbg_landIceMaskAtms


  subroutine mwbg_grossValueCheck(npts,ztb,grossrej)
    !  Object: Check Tbs for values that are missing or outside physical limits.
    !          **NOTE: REJECT ALL CHANNELS OF ONE IS FOUND TO BE BAD.
    !
    ! Variable Definitions:
    ! npts            - input  -  number of obs pts to process
    ! ztb             - input  -  Tbs from input BURP file
    ! grossrej        - output -  logical array defining which obs are to be rejected
    implicit none

    ! Arguments
    integer, intent(in) :: npts

    real,    intent(in),  dimension(:) :: ztb
    logical, intent(out), dimension(:) :: grossrej

    ! Locals
    integer :: ii, indx1, indx2

    grossrej(1:npts) = .true.
    indx1 = 1
    do ii = 1, npts

      indx2 = ii*nchanAtms
      if ( all( ztb(indx1:indx2) > 50.0 ) .and. all( ztb(indx1:indx2) < 380.0 ) ) then
        grossrej(ii) = .false.
      end if
      indx1 = indx2 + 1

    end do

    return
  end subroutine mwbg_grossValueCheck


  subroutine mwbg_firstQcCheckAtms(zenith, ilq, itt, zlat, zlon, ztb, scanpos, stnid,&
                                   nval, nt, lqc, grossrej, lsq, trn, qcflag1, qcflag2, &
                                   ican, blat, blon, lutb)
    !  This routine performs basic quality control checks on the data. It sets array
    !  lqc(nt,nchanAtms) elements to .true. to flag data with failed checks.
    !
    !  The 7 QC checks are:
    !                 1) Invalid land/sea qualifier or terrain type,
    !                 2) Invalid field of view number,
    !                 3) Satellite zenith angle missing or out of range, (> 75 deg),
    !                 4) lat,lon check (lat,lon = O(-90.), 0(-180.))
    !                 5) Change in (computed) lsq,trn from (input) ilq,itt (from MG,LG fields)
    !                      ilq= 0,1 (from hi-res land/sea mask interpolated to obs point [CMDA])
    !                      itt=-1,0 (from hi-res ice analysis  interpolated to obs point [CMDA])
    !                      lsq= 0,1 (from max interp MG (0.0 to 1.0) in box surrounding obs point)
    !                      trn=-1,0 (from max interp LG (0.0 to 1.0) in box surrounding obs point)
    !                 6) ATMS quality flag check (qual. flag elements 33078,33079,33080,33081)

    !
    !  In most cases, lqc(ii,nchanAtms) is set to .true. for all channels at point ii
    !  if the check detects a problem. In addition, Tb (ztb) is set to missing_value 
    !  for checks 3 and 4 fails.
    implicit none

    integer, intent(in), dimension(:)    :: ilq, itt, scanpos
    integer, intent(in), dimension(:)    :: ican, qcflag2
    integer, intent(in), dimension(:,:)  :: qcflag1
    integer, intent(in)                  :: nval, nt
    integer, intent(in), dimension(:)    :: lsq, trn

    logical, intent(in), dimension(:)    :: grossrej     ! dim(nt), true if 1 or more Tb fail gross error check
    logical, intent(out)                 :: lutb         ! true if Tb(ztb) are set to missing_value

    real, intent(in), dimension(:)       :: zlat, zlon
    real, intent(inout), dimension(:)    :: ztb, zenith
    
    integer, intent(in)                  :: blat, blon   ! NT box lat,lon (header)
     
    logical, intent(inout), dimension(:,:) :: lqc        ! dim(nt,nchanAtms), lqc = .false. on input
     
    character(len=9), intent(in)         :: stnid
     
    !  Locals
    integer :: ii, jj, indx1, icount
    logical :: fail, fail1, fail2

    write(*,*) '============================================================================'
    write(*,*) 'mwbg_firstQcCheckAtms: Processing data box: Stnid, lat, lon = ', stnid, blat, blon
    write(*,*) ' '

    lutb = .false.

    ! Global rejection checks

    ! Check if number of channels is correct
    if ( nval /= nchanAtms ) then
      write(*,*) 'WARNING: Number of channels (',nval, ') is not equal to nchanAtms (', nchanAtms,')'
      write(*,*) '         All data flagged as bad and returning to calling routine!'
      lqc(:,:) = .true.  ! flag all data in report as bad
      return
    endif

    ! Check for errors in channel numbers (should be 1-22 for each location ii)
    indx1 = 1
    fail = .false.
    do ii = 1,nt
      do jj = 1,nchanAtms
        if ( ican(indx1+jj-1) /= jj ) fail = .true.
      enddo
      indx1 = indx1 + nchanAtms
    enddo
    if ( fail ) then
      write(*,*) 'WARNING: Bad channel number(s) detected!'
      write(*,*) '         All data flagged as bad and returning to calling routine!'
      write(*,*) '  ican(nt*nchanAtms) array = ', ican(:)
      lqc(:,:) = .true.  ! flag all data in report as bad
      return
    endif

    ! 1) invalid land/sea qualifier or terrain type
    !  ilq = 0 (land),     1 (sea)
    !  itt = 0 (sea-ice), -1 otherwise
    !  lsq = 1 (sea, away from land/coast [MG]),      0 otherwise
    !  trn = 0 (over or near analyzed sea-ice [LG]), -1 otherwise
    do ii = 1,nt
      fail = .false.
      if ( ilq(ii) < 0  .or. ilq(ii) > 2 ) fail = .true.
      if ( itt(ii) < -1 .or. itt(ii) > 1 ) fail = .true.
      if ( fail ) then
        write(*,*) 'WARNING: Invalid land/sea qualifier or terrain type!'
        write(*,*) '  ilq, itt, (lat, lon) = ', ilq(ii), itt(ii), '(',zlat(ii), zlon(ii),')'
      endif
      if ( ilq(ii) == 0 .and. itt(ii) == 0 ) then
        fail = .true.
        write(*,*) 'WARNING: Sea ice point (itt=0) at land point (ilq=0)!'
        write(*,*) ' lat, lon =  ', zlat(ii), zlon(ii)
      endif
      if ( fail ) lqc(ii,:) = .true.
    enddo

    do ii = 1,nt
      fail = .false.
      if ( lsq(ii) < 0  .or. lsq(ii) > 2 ) fail = .true.
      if ( trn(ii) < -1 .or. trn(ii) > 1 ) fail = .true.
      if ( fail ) then
        write(*,*) 'WARNING: Invalid model-based (MG/LG) land/sea qualifier or terrain type!'
        write(*,*) '  lsq, trn, (lat, lon) = ', lsq(ii), trn(ii), '(',zlat(ii), zlon(ii),')'
      endif
      if ( fail ) lqc(ii,:) = .true.
    enddo
 
    ! 2) invalid field of view number
    do ii = 1,nt
      fail = .false.
      if ( scanpos(ii) < 1  .or. scanpos(ii) > mxscan ) then
        fail = .true.
        write(*,*) 'WARNING: Invalid field of view! scanpos, lat, lon = ', scanpos(ii), zlat(ii), zlon(ii)
      endif
      if ( fail ) lqc(ii,:) = .true.
    enddo

    ! 3) satellite zenith angle missing or out of range (> 75 deg)
    !  If bad zenith, then set Tb (and zenith) = missing value
    indx1 = 1
    do ii = 1,nt
      fail = .false.
      if ( zenith(ii) > 75.0 .or. zenith(ii) < 0. ) then
        fail = .true.
        write(*,*) 'WARNING: Bad or missing zenith angle! zenith, lat, lon = ', zenith(ii), zlat(ii), zlon(ii)
        zenith(ii) = zmisg
        lutb = .true.
      endif
      do jj = 1,nchanAtms
        if ( fail ) then
          lqc(ii,jj) = .true.
          ztb(indx1+jj-1) = zmisg
        endif
      enddo
      indx1 = indx1 + nchanAtms
    enddo

    ! 4) Lat,lon check
    ! Check for undecoded BURP file integer values of lat,lon = 0,0
    ! (usually associated with missing zenith angle and erroneous Tb=330K)

    icount = 0
    indx1 = 1
    do ii = 1,nt
      fail = .false.
      if ( zlat(ii) == -90.0  .and. zlon(ii) == -180.0 ) then
        fail = .true.
        icount =  icount + 1
        lutb = .true.
      endif
      do jj = 1,nchanAtms
        if ( fail ) then
          lqc(ii,jj) = .true.
          ztb(indx1+jj-1) = zmisg
        endif
      enddo
      indx1 = indx1 + nchanAtms
    enddo
    if ( icount > 0 ) write(*,*) 'WARNING: Bad lat,lon pair(s) detected. Number of locations = ', icount

    icount = 0
    indx1 = 1
    do ii = 1,nt
      fail = .false.
      if ( abs(zlat(ii)) > 90.0  .or. abs(zlon(ii)) > 180.0 ) then
        fail = .true.
        icount =  icount + 1
        lutb = .true.
      endif
      do jj = 1,nchanAtms
        if ( fail ) then
          lqc(ii,jj) = .true.
          ztb(indx1+jj-1) = zmisg
        endif
      enddo
      indx1 = indx1 + nchanAtms
    enddo
    if ( icount > 0 ) write(*,*) 'WARNING: Lat or lon out of range! Number of locations = ', icount

    !  5) Change in land/sea qualifier or terrain-type based on MG,LG fields
    icount = 0
    do ii = 1,nt
      fail = .false.
      if ( (ilq(ii) /= lsq(ii)) .or. (itt(ii) /= trn(ii)) ) fail = .true.
      if ( fail ) then
        icount =  icount + 1
      endif
    enddo
    if ( icount > 0 ) write(*,*) 'INFO: Num. pts with land/sea qualifier or terrain type changed (MG,LG) = ', icount

    ! 6) ATMS quality flag check (qual. flag elements 33078,33079,33080,33081)

    !  33078 Geolocation quality code     qcflag1(ii,1)  code value = 0-15 (0= OK, 15=misg)
    !  33079 Granule level quality flags  qcflag1(ii,2)  16 bit flag  (start bit 6(2^5)=32) (misg=2^16-1 = 65535)
    !  33080 Scan level quality flags     qcflag1(ii,3)  20 bit flag  (start bit 7(2^6)=64) (misg=2^20-1) 
    !  33081 Channel data quality flags   qcflag2        12 bit flag  (start bit 3(2^2)=4)  (misg=2^12-1)
    !
    !  See http://www.wmo.int/pages/prog/www/WMOCodes/WMO306_vI2/2010edition/BUFRver16/BUFR_16_0_0_TableD.pdf

    indx1 = 1
    do ii = 1,nt 
      fail1 = .false.
      fail = .false.
      if ( (qcflag1(ii,1) > 0) .or. (qcflag1(ii,2) >= 32) .or. (qcflag1(ii,3) >= 64) ) then
        write(*,*) 'WARNING: INFO BLOCK QC flag(s) indicate problem with data'
        write(*,*) ' ele33078 = ',qcflag1(ii,1),' ele33079 = ',qcflag1(ii,2),' ele33080 = ', qcflag1(ii,3)
        write(*,*) ' lat, lon = ', zlat(ii), zlon(ii)
        fail1 = .true.
        if ( grossrej(ii) ) write(*,*) ' NOTE: grossrej is also true for this point!'
      endif
      do jj = 1,nchanAtms
        fail2 = .false.
        if ( qcflag2(indx1+jj-1) >= 4 ) then
          !write(*,*) 'WARNING: DATA BLOCK QC flag ele33081 = ', qcflag2(indx1+jj-1)
          !write(*,*) '    Lat, lon, channel = ', zlat(ii), zlon(ii), ican(indx1+jj-1)
          fail2 = .true.
          fail = .true.
          !if ( (.not. fail1) .and. grossrej(ii) ) write(*,*) ' NOTE: grossrej is also true for this point!'
        endif
        if ( fail2 .or. fail1 ) lqc(ii,jj) = .true.
      enddo
      if ( fail ) write(*,*) 'WARNING: DATA BLOCK QC flag ele33081 >= 4 for one or more channels! lat, lon = ', zlat(ii), zlon(ii)
      indx1 = indx1 + nchanAtms
    enddo
     
    write(*,*) 'mwbg_firstQcCheckAtms: Total number of data processed in this box = ', nt*nchanAtms
    write(*,*) '         Total number of data flagged in this box   = ', COUNT(lqc)
    write(*,*) ' '

    return
  end subroutine mwbg_firstQcCheckAtms


  subroutine mwbg_nrlFilterAtms(ier, ni, tb23, bcor23, tb31, bcor31, tb50, bcor50, &
                                tb89, bcor89, tb165, bcor165, pangl, plat, ilansea, iglace, &
                                waterobs, grossrej, clw, si_ecmwf, si_bg, iNumSeaIce, iRej, &
                                SeaIce)
    !OBJET          Compute the following parameters using 5 ATMS channels:
    !                  - sea ice, 
    !                  - cloud liquid water (clw), 
    !                  - 2 scattering indices (si) (ECMWF, Bennartz-Grody)
    !               The five channels used are: 23Ghz, 31Ghz, 50Ghz, 89Ghz, and 165Ghz.
    !
    !NOTES*
    !                o  open water points are converted to sea-ice points if sea ice concentration >= 0.55
    !                   and iglace (itt or terrain type) is changed accordingly
    !                o  clw are missing when out-of-range parameters/Tb detected or grossrej = .true.
    !                o  clw and si only computed over open water away from coasts and sea-ice
    !                o  clw and si = -99.0 where value cannot be computed.
    !
    !REFERENCES     Ben Ruston, NRL Monterey
    !                  JCSDA Seminar 12/12/12: Impact of NPP Satellite Assimilation in the U.S. Navy Global Modeling System
    !
    !
    !ARGUMENTS      ier         - output - error return code for each location:
    !                                        0, ok,  
    !                                        1, input parameter out of range or grossrej=.true. 
    !               ni          - input  -  number of points to process (= NT)
    !               tb23        - input  -  23Ghz brightness temperature (K) -- ch. 1
    !               tb31        - input  -  31Ghz brightness temperature (K) -- ch. 2
    !               tb50        - input  -  50Ghz brightness temperature (K) -- ch. 3
    !               tb89        - input  -  89Ghz brightness temperature (K) -- ch. 16
    !               tb165       - input  -  165Ghz brightness temperature (K) -- ch. 17
    !               pangl       - input  -  satellite zenith angle (deg.)
    !               plat        - input  -  latitude (deg.)
    !               ilansea     - input  -  land/sea indicator (0=land, 1=ocean)
    !               iglace      - in/out -  terrain type (0=ice, -1 otherwise)
    !               waterobs    - in/out -  .true. if open water point (away from coasts and sea-ice)
    !               grossrej    - input  -  .true. if any channel had a gross error from mwbg_grossValueCheck
    !               clw         - output -  cloud liquid water (kg/m**2) from tb23 & tb31
    !               si_ecmwf    - output -  ECMWF scattering index from tb89 & tb165
    !               si_bg       - output -  Bennartz-Grody scattering index from tb89 & tb165
    !               iNumSeaIce  - in/out -  running counter for number of open water points
    !                                       with sea-ice detected (from algorithm)
    !               iRej        - in/out -  running counter for number of locations with bad
    !                                       pangl, plat, ilansea, or with grossrej=true
    !               SeaIce      - output -  computed sea-ice fraction from tb23 & tb50 
    !
    !               ice         - internal -  sea ice
    !             
    !
    ! Notes: In the case where an output parameter cannot be calculated, the
    !        value of this parameter is set to the missing value, i.e. -99.
    !
    implicit none

    integer    ::  i

    integer, intent(in)                   ::  ni
    integer, intent(inout)                ::  iNumSeaIce
    integer, intent(in), dimension(:)     ::  ilansea
    integer, intent(out), dimension(:)    ::  ier
    integer, intent(inout), dimension(:)  ::  iglace
    integer, intent(inout)                ::  iRej
    
    
    logical, intent(in), dimension(:)     ::  grossrej
    logical, intent(inout), dimension(:)  ::  waterobs

    real, intent(in),  dimension(:)  ::  tb23, tb31, tb50, tb89, tb165, pangl, plat
    real, intent(in),  dimension(:)  ::  bcor23, bcor31, bcor50, bcor89, bcor165
    real, intent(out), dimension(:)  ::  clw, si_ecmwf, si_bg, SeaIce 

    real, dimension(ni)              ::  ice

    real       ::  aa, deltb, abslat, cosz
    real       ::  t23, t31, t50, t89, t165

    real, parameter  ::  rmisg = -99.

    ier = 0

    ! 1) Initialise parameters:
    do i = 1, ni
      ice(i)      = rmisg
      clw(i)      = rmisg
      si_ecmwf(i) = rmisg
      si_bg(i)    = rmisg
      SeaIce(i)   = 0.0
    enddo

    ! 2) Validate input parameters:
    do i = 1, ni
      if ( pangl(i)   .lt.   0.  .or. &
           pangl(i)   .gt.  70.  .or. &
           plat(i)    .lt. -90.  .or. & 
           plat(i)    .gt.  90.  .or. &  
           ilansea(i) .lt.   0   .or. & 
           ilansea(i) .gt.   1        ) then
         ier(i) = 1
      endif

      ! Skip computations for points where all data are rejected  (bad Tb ANY channel)       
      if ( grossrej(i) ) ier(i) = 1 

    enddo

    ! 3) Compute parameters:
    do i = 1, ni

      if ( ier(i) .eq. 0 ) then

        abslat = abs(plat(i))
        cosz   = cosd(pangl(i))

        if ( mwbg_useUnbiasedObsForClw ) then
          t23 = tb23(i)
          t31 = tb31(i)
          t50 = tb50(i)
          t89 = tb89(i)
          t165 = tb165(i)
        else
          t23 = tb23(i) - bcor23(i)
          t31 = tb31(i) - bcor31(i)
          t50 = tb50(i) - bcor50(i)
          t89 = tb89(i) - bcor89(i)
          t165 = tb165(i) - bcor165(i)
        end if
        deltb = t89 - t165

        ! Check for sea-ice over water points. Set terrain type to 0 if ice>=0.55 detected.
        if ( ilansea(i) .eq. 1 ) then  ! water point

          if ( abslat .lt. 50. ) then
            ice(i) = 0.0
          else
            ice(i) = 2.85 + 0.020*t23 - 0.028*t50
          endif
          
          SeaIce(i) = ice(i)
          
          if ( ice(i) .ge. 0.55 .and. waterobs(i) ) then
            iNumSeaIce = iNumSeaIce + 1
            waterobs(i) = .false.
            iglace(i) = 0
          endif
          
        endif

        ! Compute CLW and Scattering Indices (over open water only)
        if ( waterobs(i) ) then
          if ( t23 .lt. 284. .and. t31 .lt. 284. ) then
            aa = 8.24 - (2.622 - 1.846*cosz)*cosz
            clw(i) = aa + 0.754*alog(285.0-t23) - 2.265*alog(285.0-t31)
            clw(i) = clw(i)*cosz
            if ( clw(i) .lt. 0.0 ) clw(i) = 0.0
          endif
          si_ecmwf(i) = deltb - (-46.94 + 0.248*pangl(i))
          si_bg(i)    = deltb - (-39.201 + 0.1104*pangl(i))
        endif

      else  ! ier(i) .eq. 1 case
         iRej = iRej + 1

      endif ! if ( ier(i) .eq. 0 )

      if ( mwbg_debug .and. (i .le. 100) ) then
        write(*,*) ' '
        write(*,*) ' i,tb23(i),tb31(i),tb50(i),tb89(i),tb165(i),pangl(i),plat(i), ilansea(i) = ', &
     &             i,tb23(i),tb31(i),tb50(i),tb89(i),tb165(i),pangl(i),plat(i), ilansea(i)
        write(*,*) ' ier(i),ice(i),clw(i),si_ecmwf(i),si_bg(i),iglace(i),waterobs(i) =',ier(i),ice(i),&
     &             clw(i),si_ecmwf(i),si_bg(i),iglace(i),waterobs(i)
      endif

    enddo   ! i loop over ni points

    return

  end subroutine mwbg_nrlFilterAtms


  subroutine mwbg_readCoeff(sats,chans,fovbias,coeff,nsat,nchan,nfov,npred,cinstrum,maxpred,iun,coeff_file,ptypes)
    ! max # of satellites, # of channels (24), # of FOV (60)
    ! RETURNS:

    !   sats(nsat)            = satellite names
    !   chans(nsat,nchan(i))  = channel numbers of each channel of each satellite i
    !   npred(nsat,nchan(i))  = number of predictors for each channel of each satellite i
    !   fovbias(i,j,k)        = bias for satellite i, channel j, FOV k   k=1,nfov
    !     if FOV not considered for instrument, nfov = 1 and fovbias is global bias for channel
    !   coeff(i,j,1)          = regression constant
    !   coeff(i,j,2), ..., coeff(i,j,npred(i,j)) = predictor coefficients

    !   nsat, nchan, nfov, cinstrum (output) are determined from file
    !   if returned nsat = 0, coeff_file was empty

    !   maxpred (input) is max number of predictors
    implicit none

    ! IN
    integer                        :: maxpred, iun
    character(len=90)              :: coeff_file

    ! OUT
    character(len=9), dimension(mxsat) :: sats        ! dim(maxsat), satellite names
    integer*4, dimension(mxsat,nchanAtms)   :: chans       ! dim(maxsat, maxchan), channel numbers
    real, dimension(mxsat,nchanAtms,mxscan) :: fovbias     ! dim(maxsat,maxchan,maxfov), bias as F(fov)
    real, dimension(mxsat,nchanAtms,maxpred+1) :: coeff    ! dim(maxsat,maxchan,maxpred+1)
    integer                            :: nsat, nfov, nbscan
    integer, dimension(mxsat)          :: nchan       ! dim(maxsat), number of channels
    integer, dimension(mxsat,nchanAtms)     :: npred       ! dim(maxsat, maxchan), number of predictors
    character(len=5)                   :: cinstrum    ! string: instrument (e.g. SSMIS)
    character(len=2), dimension(mxsat,nchanAtms,maxpred)  :: ptypes ! dim(maxsat,maxchan,maxpred)

    ! LOCAL
    character(len=8)               :: sat
    character(len=120)             :: line
    integer*4                      :: chan
    integer                        :: ndata, nbfov, nbpred, i, j, k, ier, istat, ii
    logical                        :: newsat
    real                           :: dummy

    coeff    = 0.0
    fovbias  = 0.0
    sats     = 'XXXXXXXXX'
    cinstrum = 'XXXXX'
    chans    = 0
    npred    = 0
    nsat     = 0
    nchan    = 0
    nfov     = 0
    ptypes   = 'XX'

    nbscan = mxscan

    ier = FNOM(iun,coeff_file,'FMT',0)

    IF (ier == 0) THEN

      WRITE(*,*)
      WRITE(*,*) 'Bias correction coefficient file open = ', coeff_file

      READ(iun,*,IOSTAT=istat)
      IF ( istat < 0 ) THEN
        WRITE(*,*) '  ERROR- File appears empty.'
        RETURN
      END IF
      REWIND(iun)

      ii = 0

      ! Loop over the satellites/channels in the file
      do
        read(iun,'(A)',IOSTAT=istat) line
        if ( istat < 0 ) EXIT
        if ( line(1:3) == 'SAT' ) then
          newsat = .true.
          read(line,'(T53,A8,1X,A5,1X,I6,1X,I8,1X,I2,1X,I3)',IOSTAT=istat) sat, cinstrum, chan, ndata, nbpred, nbfov
          if ( istat /= 0 ) then
            write(*,*) ' ERROR - reading data from SATELLITE line in coeff file!'
            return
          endif
          do i = 1, mxsat
            if ( trim(sats(i)) == trim(sat) ) then
              newsat = .false.
              ii = i
            endif
          end do
          if ( newsat ) then
            ii = ii + 1
            if ( ii > mxsat ) then
              write(*,*) ' ERROR - max number of satellites exceeded in coeff file!'
              return
            endif
            sats(ii) = sat
            if (ii > 1) nchan(ii-1) = j
            j = 1
          else
            j = j + 1
          endif
          chans(ii, j) = chan
          npred(ii, j) = nbpred
          if ( nbpred > maxpred ) then
            write(*,*) ' ERROR - max number of predictors exceeded in coeff file!'
            return
          endif
          read(iun,'(A)',IOSTAT=istat) line
          if ( line(1:3) /= 'PTY' ) then
            write(*,*) ' ERROR - list of predictors is missing in coeff file!'
            return
          endif
          if ( nbpred > 0 ) then
            read(line,'(T8,6(1X,A2))',IOSTAT=istat) (ptypes(ii,j,k),k=1,nbpred)
            if ( istat /= 0 ) then
              write(*,*) ' ERROR - reading predictor types from PTYPES line in coeff file!'
              return
            endif
          endif
          read(iun,*,IOSTAT=istat) (fovbias(ii,j,k),k=1,nbfov)
          if ( istat /= 0 ) then
            write(*,*) ' ERROR - reading fovbias in coeff file!'
            return
          endif
          if ( nbpred > 0 ) then
            read(iun,*,IOSTAT=istat) (coeff(ii,j,k),k=1,nbpred+1)
          else
            read(iun,*,IOSTAT=istat) dummy
          endif
          if ( istat /= 0 ) then
            write(*,*) ' ERROR - reading coeff in coeff file!'
            return
          endif

        else
          EXIT
        endif

      end do

      if ( ii == 0 ) then
        write(*,*) ' ERROR - No data read from coeff file!'
        return
      endif

      nsat      = ii
      nfov      = nbfov
      nchan(ii) = j
      if ( nbscan /= 0 ) then
        if ( nfov /= mxscan ) then
          write(*,*) ' INFO - Number of FOV in coeff file (nfov) does not equal default value (mxscan).'
          write(*,*) '         nfov = ', nfov
          write(*,*) '       mxscan = ', mxscan
          !call abort()
        endif
      else ! nbscan = 0 case
        if ( nfov /= 1 ) then
          write(*,*) ' INFO - Number of FOV in coeff file (nfov) does not equal default value (1).'
          write(*,*) '         nfov = ', nfov
          !call abort()
        endif
      endif

      write(*,*) ' '
      write(*,*) ' ------------- BIAS CORRECTION COEFFICIENT FILE ------------------ '
      write(*,*) ' '
      write(*,*) ' Number of satellites =     ', nsat
      write(*,*) ' Number of FOV =            ', nfov
      write(*,*) ' Max number of predictors = ', maxval(npred)
      write(*,*) ' '
      do i = 1, nsat
        write(*,*) '  Satellite = ' // sats(i)
        write(*,*) '     Number of channels = ', nchan(i)
        write(*,*) '     predictors, fovbias, coeff for each channel: '
        do j = 1, nchan(i)
          write(*,*) i, chans(i,j)
          if ( npred(i,j) > 0 ) then 
            write(*,'(6(1X,A2))') (ptypes(i,j,k),k=1,npred(i,j))
          else
            write(*,'(A)') 'No predictors'
          endif
          write(*,*) (fovbias(i,j,k),k=1,nfov)
          write(*,*) (coeff(i,j,k),k=1,npred(i,j)+1)
        end do
      end do
      write(*,*) ' '

      ier = FCLOS(iun)

    ELSE

      write(*,*) 'mwbg_readCoeff: ERROR - Problem opening the coeff file!'

    ENDIF

    RETURN

  end subroutine mwbg_readCoeff


end module bgckmicrowave_mod

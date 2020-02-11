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
  use MathPhysConstants_mod

  implicit none
  save
  private

  ! public variables
  public :: mwbg_debug, mwbg_clwQcThreshold, mwbg_allowStateDepSigmaObs
  public :: mwbg_modlsqtt, mwbg_useUnbiasedObsForClw 

  ! Public functions/subroutines
  public :: mwbg_readStatTovs, mwbg_tovCheckAmsua, mwbg_qcStats
  public :: mwbg_updateBurpAmsua
  public :: mwbg_readGeophysicFieldsAndInterpolate
  public :: mwbg_burpErrorHistory
  public :: mwbg_setTerrainTypeToSeaIce
  public :: mwbg_readStatTovsAtms, mwbg_tovCheckAtms, mwbg_qcStatsAtms
  public :: mwbg_updatFlgAtms, mwbg_getData, mwbg_getBurpReportAdresses, mwbg_landIceMaskAtms
  public :: mwbg_grossValueCheck, mwbg_firstQcCheckAtms, mwbg_nrlFilterAtms
  public :: mwbg_writeBlocks

  logical :: mwbg_debug, mwbg_allowStateDepSigmaObs
  real    :: mwbg_clwQcThreshold
  logical :: mwbg_modlsqtt, mwbg_useUnbiasedObsForClw 

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
  real    :: clwThreshArr(JPCH,JPNSAT,2)
  real    :: sigmaObsErr(JPCH,JPNSAT,2)
  integer :: useStateDepSigmaObs(JPCH,JPNSAT)

  INTEGER :: MREJCOD(JPMXREJ,MXCHN,MXSAT), INTOT(MXSAT), INTOTRJF(MXSAT), INTOTRJP(MXSAT)
  INTEGER :: MREJCOD2(JPMXREJ,MXCHN,MXSAT)

  integer  :: error

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
  END FUNCTION ISRCHEQR


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
  END FUNCTION ISRCHEQI


  SUBROUTINE mwbg_tovCheckAmsua(KSAT, KTERMER, KORBIT, ICANO, ICANOMP, ZO, ZCOR, &
                                ZOMP, ICHECK, KNO, KNT, PMISG, KNOSAT, KCHKPRF, &
                                ISCNPOS, MGINTRP, MTINTRP, GLINTRP, ITERRAIN, SATZEN, &
                                IMARQ, clw, clw_avg, scatw, STNID, RESETQC, ZLAT)
    !:Purpose:          Effectuer le controle de qualite des radiances tovs.
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
    implicit none 
    !Arguments:
    integer, intent(in)                    :: KSAT(KNT)          ! numero d'identificateur du satellite
    integer, intent(in)                    :: KTERMER(KNT)      ! indicateur terre/mer
    integer, intent(in)                    :: ISCNPOS(KNT)      ! position sur le "scan"
    integer, intent(in)                    :: KORBIT(KNT)      ! numero d'orbite
    integer, intent(in)                    :: ICANO(KNO*KNT)     ! canaux des observations
    integer, intent(in)                    :: ITERRAIN(KNT)      ! indicateur du type de terrain
    integer, intent(in)                    :: ICANOMP(KNO*KNT)  ! canaux des residus (o-p)
    integer, intent(in)                    :: KNO                ! nombre de canaux des observations 
    integer, intent(in)                    :: KNT                ! nombre de tovs
    integer, intent(in)                    :: KNOSAT             ! numero de satellite (i.e. indice)
    integer, intent(inout)                 :: IMARQ(KNO*KNT)     ! marqueurs des radiances
    real, intent(in)                       :: ZO(KNO*KNT)        ! radiances
    real, intent(in)                       :: ZCOR(KNO*KNT)      ! correction aux radiances
    real, intent(in)                       :: ZOMP(KNO*KNT)  ! residus (o-p)
    real, intent(in)                       :: MGINTRP(KNT)      ! masque terre/mer du modele
    real, intent(in)                       :: MTINTRP(KNT)      ! topographie du modele
    real, intent(in)                       :: GLINTRP(KNT)      ! etendue de glace du modele
    real, intent(in)                       :: SATZEN(KNT)      ! angle zenith du satellite (deg.)
    real, intent(in)                       :: ZLAT(KNT)      ! latitude
    real, intent(in)                       :: PMISG              ! missing value
    character *9, intent(in)               :: STNID              ! identificateur du satellite
    logical, intent(in)                    :: RESETQC            ! reset du controle de qualite?
    integer,allocatable, intent(out)       :: ICHECK(:,:)    ! indicateur controle de qualite tovs par canal 
    !                                                              =0, ok,
    !                                                              >0, rejet,
    real,allocatable,  intent(out)        :: clw(:)         ! retrieved cloud liquid water
    real,allocatable,  intent(out)        :: clw_avg(:)     ! Averaged retrieved cloud liquid water, 
    !                                                              from observation and background
    real,allocatable, intent(out)         :: scatw(:)       ! scattering index over water

    integer,allocatable, intent(out)      :: KCHKPRF(:)       ! indicateur global controle de qualite tovs. Code:
    !                                                              =0, ok,
    !                                                              >0, rejet d'au moins un canal

    !locals
    integer, parameter                     :: MXSCANHIRS= 56 
    integer, parameter                     :: MXSCANAMSU= 30 
    integer, parameter                     :: MXCLWREJ  =  6 
    integer, parameter                     :: MXSFCREJ  =  6 
    integer, parameter                     :: MXSFCREJ2 =  4 
    integer, parameter                     :: MXSCATREJ =  7 
    integer, parameter                     :: MXCANPRED =  9 
    integer, parameter                     :: JPMXSFC = 2
    real, parameter                        :: cloudyClwThreshold = 0.3
    integer                                :: KMARQ   (KNO,KNT)
    integer                                :: KCANO   (KNO,KNT)
    integer                                :: KCANOMP (KNO,KNT)
    real                                   :: PTBO    (KNO,KNT)
    real                                   :: PTBCOR  (KNO,KNT)
    real                                   :: PTBOMP  (KNO,KNT)

    integer                                :: MAXVAL
    integer                                :: JI
    integer                                :: JJ
    integer                                :: INDX8
    integer                                :: INDX12
    integer                                :: INO
    integer                                :: ICHN
    integer                                :: JK
    integer                                :: IBIT
    integer                                :: JC
    integer                                :: INDX
    integer                                :: INDXCAN
    integer                                :: ITRN
    integer                                :: alloc_status 
    integer                                :: ICLWREJ (MXCLWREJ)
    integer                                :: ISFCREJ (MXSFCREJ)
    integer                                :: ISFCREJ2(MXSFCREJ2)
    integer                                :: ISCATREJ(MXSCATREJ)
    real                                   :: EPSILON
    real                                   :: ZANGL
    real                                   :: MISGRODY
    real                                   :: ZSEUILSCAT
    real                                   :: APPROXIM
    real                                   :: ANGDif
    real                                   :: XCHECKVAL
    real                                   :: GROSSMIN(MXCHN)
    real                                   :: GROSSMAX(MXCHN) 
    real                                   :: ROGUEFAC(MXCHN)
    real                                   :: tb23 (KNT)
    real                                   :: tb31 (KNT)
    real                                   :: tb50 (KNT)
    real                                   :: tb53 (KNT)
    real                                   :: tb89 (KNT)
    real                                   :: tb23_P (KNT)
    real                                   :: tb31_P (KNT)
    real                                   :: tb50_P (KNT)
    real                                   :: tb53_P (KNT)
    real                                   :: tb89_P (KNT)
    real                                   :: ice  (KNT)
    real                                   :: tpw  (KNT)
    real                                   :: scatl(KNT)
    integer                                :: err (KNT)
    integer                                :: rain(KNT)
    integer                                :: snow(KNT)
    logical                                :: GROSSERROR
    logical                                :: FULLREJCT
    logical                                :: SFCREJCT
    logical, save                          :: LLFIRST


    data  LLFIRST / .TRUE. /
    data  EPSILON / 0.01   /
    data  MISGRODY / -99.     /
    !      data  ROGUEFAC/ 3.0    / changed, jh, from 3 to 4, jan 2001
    !      data  ROGUEFAC/ 4.0    /
    data  ROGUEFAC / 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
>>>>>>> Issue #308: modifications to microwave module and bgckMW prog.
                     4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
                     4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 2.0, 2.0, 2.0, &
                     3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
                     4.0, 2.0/
    data  ICLWREJ  / 28, 29, 30, 31, 32, 42 /
    data  ISFCREJ  / 28, 29, 30, 31, 32, 42 /
    data  ISCATREJ / 28, 29, 30, 31, 32, 33, 42 /
    data  ISFCREJ2 / 28, 29, 30, 42 /
                   
    data GROSSMIN / 200., 190., 190., 180., 180., 180., 170., &
                    170., 180., 170., 170., 170., 180., 180., &
                    180., 180., 170., 180., 180., 000., 120., &
                    190., 180., 180., 180., 190., 200., 120., &
                    120., 160., 190., 190., 200., 190., 180., &
                    180., 180., 180., 190., 190., 200., 130./

    data GROSSMAX / 270., 250., 250., 250., 260., 280., 290., &
                    320., 300., 320., 300., 280., 320., 300., &
                    290., 280., 330., 350., 350., 000., 310., &
                    300., 250., 250., 270., 280., 290., 310., &
                    310., 310., 300., 300., 260., 250., 250., &
                    250., 260., 260., 270., 280., 290., 330./  

    ! Allocation
    alloc_status = 0
    if(allocated(clw_avg)) deallocate(clw_avg)
    allocate(clw_avg(knt), stat = alloc_status)
    if (alloc_status /= 0) then
      write(*,*) ' Allocation Error in sub. mwbg_tovCheckAmsua '
      call abort()
    end if 
    alloc_status = 0
    if(allocated(clw)) deallocate(clw)
    allocate(clw(knt), stat = alloc_status)
    if (alloc_status /= 0) then
      write(*,*) ' Allocation Error in sub. mwbg_tovCheckAmsua '
      call abort()
    end if 
    alloc_status = 0
    if(allocated(scatw)) deallocate(scatw)
    allocate(scatw(knt), stat = alloc_status)
    if (alloc_status /= 0) then
      write(*,*) ' Allocation Error in sub. mwbg_tovCheckAmsua '
      call abort()
    end if 
    alloc_status = 0
    if(allocated(kchkprf)) deallocate(kchkprf)
    allocate(kchkprf(knt), stat = alloc_status)
    if (alloc_status /= 0) then
      write(*,*) ' Allocation Error in sub. mwbg_tovCheckAmsua '
      call abort()
    end if 
    alloc_status = 0
    if(allocated(icheck)) deallocate(icheck)
    allocate(icheck(kno,knt), stat = alloc_status)
    if (alloc_status /= 0) then
      write(*,*) ' Allocation Error in sub. mwbg_tovCheckAmsua '
      call abort()
    end if 

    ! copy the original input 1D array to 2D array. The 2D arrays are used in this s/r.
    do JJ=1,KNT
      do JI=1,KNO
        INDX = (JJ-1)*KNO + JI 
        KCANO(JI,JJ) = ICANO(INDX)
        PTBCOR(JI,JJ) = ZCOR(INDX)  
        PTBO(JI,JJ) = ZO(INDX)  
        KMARQ(JI,JJ) = IMARQ(INDX)  

        if ( ICANO(INDX) /= ICANOMP(INDX) ) then
          write(*,*)'ERROR IN DIMENSIONS OF TOVS data'
          CALL ABORT()
        end if

        KCANOMP(JI,JJ) = ICANOMP(INDX)
        PTBOMP(JI,JJ) = ZOMP(INDX)  
      end do
    end do

    ! Initialisation, la premiere fois seulement!
    if (LLFIRST) then
       do JI = 1, JPMXREJ
          do JJ = 1, MXCHN
             do JK = 1, MXSAT
                MREJCOD(JI,JJ,JK) = 0
             end do
          end do
       end do
       LLFIRST = .FALSE.
    end if

    do JJ=1,KNT
      do JI=1,KNO
        if ( KCANO(JI,JJ) .NE. KCANOMP(JI,JJ) ) then
          write(*,*)'INCONSISTENT CHANNEL LISTS FOR TOVS data'
          CALL ABORT()
        end if
      end do
    end do

    ! Initialisations
    do JJ=1,KNT
      do JI=1,KNO
        ICHECK(JI,JJ) = 0
        if ( RESETQC ) KMARQ(JI,JJ) = 0
      end do
    end do

    !  Run Grody AMSU-A algorithms.
    !     Grody parameters.
    !     extract required channels:
    !        23 Ghz = AMSU-A 1 = channel #28
    !        31 Ghz = AMSU-A 2 = channel #29
    !        50 Ghz = AMSU-A 3 = channel #30
    !        53 Ghz = AMSU-A 5 = channel #32
    !        89 Ghz = AMSU-A15 = channel #42
    do JJ=1,KNT
      do JI=1,KNO
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

          if ( ichn .eq. 28 ) tb23_P(jj) = ptbo(ji,jj) - ptbomp(ji,jj)
          if ( ichn .eq. 29 ) tb31_P(jj) = ptbo(ji,jj) - ptbomp(ji,jj)
          if ( ichn .eq. 30 ) tb50_P(jj) = ptbo(ji,jj) - ptbomp(ji,jj)
          if ( ichn .eq. 32 ) tb53_P(jj) = ptbo(ji,jj) - ptbomp(ji,jj)
          if ( ichn .eq. 42 ) tb89_P(jj) = ptbo(ji,jj) - ptbomp(ji,jj)
        else
          if ( ichn .eq. 28 ) tb23(jj) = 0.
          if ( ichn .eq. 29 ) tb31(jj) = 0.
          if ( ichn .eq. 30 ) tb50(jj) = 0.
          if ( ichn .eq. 32 ) tb53(jj) = 0.
          if ( ichn .eq. 42 ) tb89(jj) = 0.

          if ( ichn .eq. 28 ) tb23_P(jj) = 0.  
          if ( ichn .eq. 29 ) tb31_P(jj) = 0. 
          if ( ichn .eq. 30 ) tb50_P(jj) = 0. 
          if ( ichn .eq. 32 ) tb53_P(jj) = 0. 
          if ( ichn .eq. 42 ) tb89_P(jj) = 0. 
        endif
      end do
    end do

    call grody (err, knt, tb23, tb31, tb50, tb53, tb89, &
                tb23_P, tb31_P, tb50_P, tb53_P, tb89_P, &
                satzen, zlat, ktermer, ice, tpw, clw, clw_avg, &
                rain, snow, scatl, scatw)   

    ! 10) test 10: RTTOV reject check (single)
    ! Rejected datum flag has bit #9 on.
    if (.NOT.RESETQC) then
      INO = 10
      do JJ=1,KNT
        do JI=1,KNO
          if ( KCANO(JI,JJ) .NE. 20 ) then
            IBIT = AND(KMARQ(JI,JJ), 2**9)
            if ( IBIT .NE. 0  ) then
              ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
              KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
              MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                   MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
              if ( mwbg_DEBUG ) then
                write(*,*)STNID(2:9),' RTTOV REJECT.', &
                          'CHANNEL=', KCANO(JI,JJ), &
                          ' IMARQ= ',KMARQ(JI,JJ)
              end if
            end if
          end if
        end do
      end do
    end if

    ! 1) test 1: Topography check (partial)
    ! Channel 6 is rejected for topography >  250m.
    ! Channel 7 is rejected for topography > 2000m.
    INO = 1
    do JJ=1,KNT
      do JI=1,KNO
        if ( KCANO(JI,JJ) .EQ. 33 ) then
          if ( MTINTRP(JJ) .GE. 250.  ) then
            ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
            KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
            KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**18)
            MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                 MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
            if ( mwbg_DEBUG ) then
              write(*,*)STNID(2:9),' TOPOGRAPHY REJECT.', &
                        'CHANNEL=', KCANO(JI,JJ), &
                        ' TOPO= ',MTINTRP(JJ)
            end if
          end if
        ELSEif ( KCANO(JI,JJ) .EQ. 34 ) then
          if ( MTINTRP(JJ) .GE. 2000.  ) then
            ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
            KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
            KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**18)
            MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                 MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
            if ( mwbg_DEBUG ) then
              write(*,*)STNID(2:9),' TOPOGRAPHY REJECT.', &
                        'CHANNEL=', KCANO(JI,JJ), &
                        ' TOPO= ',MTINTRP(JJ)
            end if
          end if
        end if
      end do
    end do

    ! 2) test 2: "Land/sea qualifier" code check (full)
    ! allowed values are: 0, land,
    !                       1, sea,
    !                       2, coast.
    INO = 2
    do JJ=1,KNT
      if ( KTERMER(JJ) .LT.  0  .OR. &
          KTERMER(JJ) .GT.  2        ) then
        do JI=1,KNO
          ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
          KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
          KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
          MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
            MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
        end do
        if ( mwbg_DEBUG ) then
          write(*,*) STNID(2:9),'LAND/SEA QUALifIER CODE', &
                   ' REJECT. KTERMER=', KTERMER(JJ)
        end if
      end if
    end do

    ! 3) test 3: "Terrain type" code check (full)
    !   allowed values are: -1, missing,
    !                        0, sea-ice,
    !                        1, snow on land.
    INO = 3
    do JJ=1,KNT
      if ( ITERRAIN(JJ) .NE.  MISGINT ) then
        if ( ITERRAIN(JJ) .LT.  0  .OR. &
           ITERRAIN(JJ) .GT.  1        ) then
          do JI=1,KNO
            ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
            KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
            KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
            MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
              MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
          end do
          if ( mwbg_debug ) then
            write(*,*)STNID(2:9),'TERRAIN TYPE CODE', &
                     ' REJECT. TERRAIN=', ITERRAIN(JJ)
          end if
        end if
      end if
    end do

    ! 4) test 4: Field of view number check (full)
    !
    ! Field of view acceptable range is [1,MXSCANAMSU]  for AMSU footprints.
    INO = 4
    do JJ=1,KNT
      do JI=1,KNO
        if ( ISCNPOS(JJ) .LT. 1 .OR. &
            ISCNPOS(JJ) .GT. MXSCANAMSU ) then
          ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
          KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
          KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
          MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
            MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
          if ( mwbg_debug ) then
            write(*,*)STNID(2:9),'FIELD OF VIEW NUMBER', &
                      ' REJECT. FIELD OF VIEW= ', ISCNPOS(JJ)
          end if
        end if
      end do
    end do

    ! 5) test 5: Satellite zenith angle check (full)
    ! Satellite zenith angle acceptable range is [0.,60.].
    INO = 5
    do JJ=1,KNT
      if ( SATZEN(JJ) .NE.  PMISG ) then
        if ( SATZEN(JJ) .LT.  0.  .OR. &
           SATZEN(JJ) .GT. 60.       ) then
          do JI=1,KNO
            ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
            KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
            KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
            MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
               MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
          end do
          if ( mwbg_debug ) then
            write(*,*)STNID(2:9),' SATELLITE ZENITH ANGLE', &
                      ' REJECT. SATZEN= ', &
                      SATZEN(JJ)
          end if
        end if
      end if
    end do

    ! 6) test 6: "Sat. zenith angle"/"field of view" consistency check.  (full)
    ! Acceptable difference between "Satellite zenith angle"  and
    ! "approximate angle computed from field of view number" is 1.8 degrees.
    ZANGL   = 3.92

    INO = 6
    do JJ=1,KNT
      if ( SATZEN (JJ) .NE.  PMISG   .AND. &
         ISCNPOS(JJ) .NE.  MISGINT       ) then
        APPROXIM = ABS((ISCNPOS(JJ)-MXSCANAMSU/2.-0.5)*ZANGL)
        ANGDif = ABS(SATZEN (JJ)-APPROXIM)
        if ( ANGDif .GT. 1.8 ) then 
          do JI=1,KNO
            ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
            KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
            KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
            MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
               MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
          end do
          if ( mwbg_debug ) then
            write(*,*)STNID(2:9),' ANGLE/FIELD OF VIEW', &
                      ' INCONSISTENCY REJECT. SATZEN= ', &
                      SATZEN(JJ), ' FIELD OF VIEW= ',ISCNPOS(JJ), &
                      ' ANGDif= ',ANGDif  
          end if
        end if
      end if
    end do

    ! 7) test 7: "Land/sea qual."/"model land/sea" consistency check.    (full)
    ! Acceptable conditions are:
    !       a) both over ocean (ktermer=1; mg<0.01), new threshold 0.20, jh dec 2000,
    !       b) both over land  (ktermer=0; mg>0.80), new threshold 0.50, jh dec 2000.
    ! Other conditions are unacceptable.
    INO = 7
    do JJ=1,KNT
      if ( KTERMER (JJ) .NE.  MISGINT  ) then
        if     ( KTERMER(JJ) .EQ. 1       .AND. &
                MGINTRP(JJ) .LT. 0.20          ) then
        ELSEif ( KTERMER(JJ) .EQ. 0       .AND. &
                MGINTRP(JJ) .GT. 0.50          ) then
        ELSE
          do JI=1,KNO
            ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
            KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
            KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
            MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
               MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
          end do
          if ( mwbg_debug ) then
            write(*,*)STNID(2:9),' LAND/SEA QUALifIER', &
                      ' INCONSISTENCY REJECT. KTERMER= ', &
                      KTERMER(JJ), ' MODEL MASK= ',MGINTRP(JJ)
          end if
        end if
      end if
    end do

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
    do JJ=1,KNT
      if ( ITERRAIN (JJ) .NE.  MISGINT  ) then
        if     ( ITERRAIN(JJ) .EQ. 0       .AND. &
                GLINTRP (JJ) .LT. 0.01          ) then
          do JI=1,KNO
            ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
            MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
               MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
          end do
          if ( mwbg_debug ) then
            write(*,*)STNID(2:9),' TERRAIN TYPE/MODEL ICE', &
                      ' INCONSISTENCY REJECT. TERRAIN= ', &
                      ITERRAIN(JJ), ' MODEL ICE= ',GLINTRP(JJ)
          end if
        end if
        if ( ITERRAIN(JJ) .EQ. 0       .AND. &
               KTERMER (JJ) .EQ. 0           ) then
          do JI=1,KNO
            ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
            MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
               MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
          end do
          if ( mwbg_debug ) then
            write(*,*)STNID(2:9),' TERRAIN TYPE/LAND?SEA QUAL.', &
                      ' INCONSISTENCY REJECT. TERRAIN= ', &
                      ITERRAIN(JJ), ' LAND/SEA= ',KTERMER(JJ)
          end if
        end if
        if ( ITERRAIN(JJ) .EQ. 1       .AND. &
                KTERMER (JJ) .EQ. 1           ) then
          do JI=1,KNO
            ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
            MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
               MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
          end do
          if ( mwbg_debug ) then
            write(*,*)STNID(2:9),' TERRAIN TYPE/LAND?SEA QUAL.', &
                      ' INCONSISTENCY REJECT. TERRAIN= ', &
                      ITERRAIN(JJ), ' LAND/SEA= ',KTERMER(JJ)
          end if
        end if
      ELSE
        if     ( KTERMER (JJ) .EQ. 1      .AND. &
                GLINTRP (JJ) .GT. 0.01         ) then
  !             do JI=1,KNO
  !                ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
  !                MREJCOD(INO,KCANO(JI,JJ),KNOSAT) =
  !                  MREJCOD(INO,KCANO(JI,JJ),KNOSAT) + 1
  !             end do
           if ( mwbg_debug ) then
              write(*,*)STNID(2:9),' TERRAIN TYPE MSG/MODEL ICE', &
                       ' INCONSISTENCY REJECT. TERRAIN= ', &
                       ITERRAIN(JJ), &
                       ' LAND/SEA= ',KTERMER(JJ)
           end if
        end if
      end if
    end do
800  continue

    ! 9) test 9: Uncorrected Tb check (single)
    ! Uncorrected datum (flag bit #6 off). In this case switch bit 11 ON.
    if (.NOT.RESETQC) then
      INO = 9
      do JJ=1,KNT
        do JI=1,KNO
          if ( KCANO(JI,JJ) .NE. 20 ) then
            IBIT = AND(KMARQ(JI,JJ), 2**6)
            if ( IBIT .EQ. 0  ) then
              ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
              KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**11)
              MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                    MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
              if ( mwbg_debug ) then
                write(*,*)STNID(2:9),' UNCORRECTED TB REJECT.', &
                           'CHANNEL=', KCANO(JI,JJ), ' IMARQ= ',KMARQ(JI,JJ)
              end if
            end if
          end if
        end do
      end do
    end if

    ! 11) test 11: Radiance observation "Gross" check (single) 
    !  Change this test from full to single. jh nov 2000.
    INO = 11
    do JJ=1,KNT
      GROSSERROR = .FALSE.
      do JI=1,KNO
        if ( KCANO(JI,JJ) .NE. 20     .AND. &
            KCANO(JI,JJ) .GE.  1     .AND. &
            KCANO(JI,JJ) .LE.  MXCHN       ) then  
          if ( PTBO(JI,JJ) .NE. PMISG .AND. &
             ( PTBO(JI,JJ).LT.GROSSMIN(KCANO(JI,JJ)).OR. &
               PTBO(JI,JJ).GT.GROSSMAX(KCANO(JI,JJ))     ) ) then
            GROSSERROR = .TRUE.
            ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
            KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
            KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
            MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                   MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
            if ( mwbg_debug ) then
              write(*,*)STNID(2:9),' GROSS CHECK REJECT.', &
                        'CHANNEL=', KCANO(JI,JJ), ' TB= ',PTBO(JI,JJ)
            end if
          end if
        end if
      end do
    end do

    ! 12) test 12: Grody cloud liquid water check (partial)
    ! For Cloud Liquid Water > clwQcThreshold, reject AMSUA-A channels 1-5 and 15.
    INO = 12
    do JJ=1,KNT
      if ( CLW(JJ) .NE.  MISGRODY  ) then
        if ( CLW(JJ) .GT. mwbg_clwQcThreshold   ) then
          do JI=1,KNO
            INDXCAN = ISRCHEQI (ICLWREJ,MXCLWREJ,KCANO(JI,JJ))
            if ( INDXCAN.NE.0 )  then
              ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
              KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
              KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
              MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                       MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
            end if
          end do
          if ( mwbg_debug ) then
            write(*,*)STNID(2:9),'Grody cloud liquid water check', &
                      ' REJECT. CLW= ',CLW(JJ), ' SEUIL= ',mwbg_clwQcThreshold
          end if
        end if

        ! trun on bit=23 for cloud-affected radiances (to be used in gen_bias_corr)
        if ( CLW(JJ) > cloudyClwThreshold ) then
          do JI = 1,KNO
            INDXCAN = ISRCHEQI(ICLWREJ,MXCLWREJ,KCANO(JI,JJ))
            if ( INDXCAN /= 0 ) KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**23)
          end do
          if ( mwbg_debug ) then
            write(*,*) STNID(2:9),' Grody cloud liquid water check', &
                      ' cloud-affected obs. CLW= ',CLW(JJ), ', threshold= ',cloudyClwThreshold 
          end if
        end if

      end if
    end do

    ! 13) test 13: Grody scattering index check (partial)
    ! For Scattering Index > 9, reject AMSUA-A channels 1-6 and 15.
    INO = 13
    ZSEUILSCAT = 9.0
    do JJ=1,KNT
      if ( SCATW(JJ) .NE.  MISGRODY  ) then
        if (  KTERMER (JJ) .EQ.  1 .AND. &
             ITERRAIN(JJ) .NE.  0 .AND. &   
             SCATW   (JJ) .GT. ZSEUILSCAT   ) then
          do JI=1,KNO
            INDXCAN = ISRCHEQI (ISCATREJ,MXSCATREJ,KCANO(JI,JJ))
            if ( INDXCAN.NE.0 )  then
              ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
              KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
              KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
              MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                       MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
            end if
          end do
          if ( mwbg_debug ) then
            write(*,*)STNID(2:9),'Grody scattering index check', &
                       ' REJECT. SCATW= ',SCATW(JJ), ' SEUIL= ',ZSEUILSCAT
          end if
        end if
      end if
    end do

    ! 14) test 14: "Rogue check" for (O-P) Tb residuals out of range. (single/full)
    ! Les observations, dont le residu (O-P) depasse par un facteur (roguefac) l'erreur totale des TOVS.
    ! N.B.: a reject by any of the 3 surface channels produces the rejection of AMSUA-A channels 1-5 and 15. 
    INO = 14
    do JJ=1,KNT
      SFCREJCT = .FALSE.
      do JI=1,KNO
        ICHN = KCANO(JI,JJ)
        if ( ICHN .NE. 20 ) then
          XCHECKVAL = ROGUEFAC(ICHN)*TOVERRST(ICHN,KNOSAT) 
          if ( PTBOMP(JI,JJ)      .NE. PMISG    .AND. &
              ABS(PTBOMP(JI,JJ)) .GE. XCHECKVAL     ) then
            ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
            KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
            KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**16)
            MREJCOD(INO,ICHN,KNOSAT) = &
                MREJCOD(INO,ICHN,KNOSAT) + 1 
            if ( mwbg_debug ) then
              write(*,*)STNID(2:9),'ROGUE CHECK REJECT.NO.', &
                      ' OBS = ',JJ, &
                      ' CHANNEL= ',ICHN, &
                      ' CHECK VALUE= ',XCHECKVAL, &
                      ' TBOMP= ',PTBOMP(JI,JJ)
            end if
            if ( ICHN .EQ. 28 .OR. &
                 ICHN .EQ. 29 .OR. &
                 ICHN .EQ. 30      ) then
              SFCREJCT = .TRUE.
            end if
          end if
        end if
      end do

      if ( SFCREJCT ) then
        do JI=1,KNO
          INDXCAN = ISRCHEQI (ISFCREJ,MXSFCREJ,KCANO(JI,JJ))
          if ( INDXCAN .NE. 0 )  then
            if ( ICHECK(JI,JJ) .NE. INO ) then
               ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
               KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
               KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**16)
               MREJCOD(INO,KCANO(JI,JJ),KNOSAT) = &
                         MREJCOD(INO,KCANO(JI,JJ),KNOSAT)+ 1
            end if
          end if
        end do
      end if

    end do

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

    do JJ=1,KNT
      ITRN = ITERRAIN(JJ)
      if ( KTERMER (JJ) .EQ. 1    .AND. &
           ITERRAIN(JJ) .EQ. -1   .AND. &
           GLINTRP (JJ) .GE. 0.01       ) then
         ITRN = 0
      end if        
      do JI=1,KNO
          ICHN = KCANO(JI,JJ)
          INDXCAN = ISRCHEQI (ISFCREJ2,MXSFCREJ2,ICHN)
          if ( INDXCAN .NE. 0 )  then
            if ( KTERMER (JJ) .EQ. 0 .OR. ITRN .EQ. 0 )  then
              KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
              KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
            end if
          end if
          if ( IUTILST(ICHN,KNOSAT) .NE. 1 ) then
            SFCREJCT = .FALSE.
            if ( IUTILST(ICHN,KNOSAT) .EQ. 0 ) then
              SFCREJCT = .TRUE.
              KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**11)
            ELSE 
              if ( KTERMER (JJ) .EQ. 0 .OR. ITRN .EQ. 0 )  then
                SFCREJCT = .TRUE.
                KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**9)
                KMARQ(JI,JJ) = OR(KMARQ(JI,JJ),2**7)
              end if
            end if
            if ( SFCREJCT ) then
              ICHECK(JI,JJ) = MAX(ICHECK(JI,JJ),INO)
              MREJCOD(INO,ICHN,KNOSAT) = & 
                 MREJCOD(INO,ICHN,KNOSAT) + 1 
              if ( mwbg_debug ) then
                 write(*,*)STNID(2:9),'CHANNEL REJECT: ', &
                        ' OBS = ',JJ, &
                        ' CHANNEL= ',ICHN
              end if
            end if
          end if
        end do
    end do

    if ( mwbg_debug ) then
       write(*,*)'ICHECK = ',((ICHECK(JI,JJ),JI=1,KNO),JJ=1,KNT)
    end if

    !  Synthese de la controle de qualite au niveau de chaque point
    !  d'observation. Code:
    !            =0, aucun rejet,
    !            >0, au moins un canal rejete.

    do JJ=1,KNT
      KCHKPRF(JJ) = 0
      do JI=1,KNO 
        KCHKPRF(JJ) = MAX(KCHKPRF(JJ),ICHECK(JI,JJ))
      end do
    end do

    if ( mwbg_debug ) then
      write(*,*)'KCHKPRF = ',(KCHKPRF(JJ),JJ=1,KNT)
    end if 

    ! Copy the modified FLAG to the 1D array, used outside this s/r. 
    do JJ=1,KNT
      do JI=1,KNO
        INDX = (JJ-1)*KNO + JI 
        IMARQ(INDX) = KMARQ(JI,JJ)
      end do
    end do

    return
  end subroutine mwbg_tovCheckAmsua


  subroutine mwbg_qcStats(instName, satNumber, ICHECK, KCANO, KNOSAT, &
                              KNO, KNT, satelliteId, LDprint)
    !:Purpose:          Cumuler ou imprimer des statistiques decriptives
    !                   des rejets tovs.
    implicit none 
    !Arguments:
    character(*), intent(in)               :: instName               ! Instrument Name
    integer, intent(in)                    :: satNumber              ! Nombre de satellite portant instrument
    integer, intent(in)                    :: ICHECK(KNO,KNT)        ! indicateur controle de qualite tovs par canal 
    !                                                                  =0, ok,
    !                                                                  >0, rejet,
    integer, intent(in)                    :: KCANO(KNO,KNT)         ! canaux des observations
    integer, intent(in)                    :: KNOSAT                 ! numero d'identificateur du satellite
    integer, intent(in)                    :: KNO                    ! nombre de canaux des observations 
    integer, intent(in)                    :: KNT                    ! nombre de tovs
    character *9, intent(in)               :: satelliteId(satNumber) ! identificateur du satellite
    logical, intent(in)                    :: LDprint                ! mode: imprimer ou cumuler?
    !Locals
    integer                                :: JI
    integer                                :: JJ
    integer                                :: JK
    integer                                :: INTOTOBS
    integer                                :: INTOTACC


    logical, save                          :: LLFIRST = .True.
    logical                                :: FULLREJCT
    logical                                :: FULLACCPT
    logical                                :: legendRejetAmsua
    logical                                :: legendRejetAtms

    legendRejetAmsua = .false.
    legendRejetAtms = .false.

    if (instName == 'AMSUA') then
      legendRejetAmsua = .true.
    else if (instName == 'ATMS') then
      legendRejetAtms = .true.
    end if 
    ! Initialize
    if ( LLFIRST ) then
      do JJ = 1, MXSAT
        INTOTRJF(JJ) = 0
        INTOTRJP(JJ) = 0
        INTOT(JJ)  = 0
      end do
      LLFIRST = .false.
    end if

    if (.NOT. LDprint ) then

      ! Accumulate statistics on rejects
      do JJ = 1, KNT

        INTOT(KNOSAT) = INTOT(KNOSAT) + 1

        ! Fully accepted, fully rejected or partially rejected?
        FULLREJCT = .true.
        FULLACCPT = .true.
        do JI = 1, KNO
          if ( KCANO(JI,JJ) .NE. 20 ) then
            if ( ICHECK(JI,JJ) .NE. 0 ) then
              FULLACCPT = .false.
            else
              FULLREJCT = .false.
            end if
          end if
        end do
        if ( FULLREJCT ) then
          INTOTRJF(KNOSAT) = INTOTRJF(KNOSAT) + 1
        end if
        if ( .NOT.FULLREJCT .AND. .NOT.FULLACCPT ) then
          INTOTRJP(KNOSAT) = INTOTRJP(KNOSAT) + 1
        end if
      end do

    else

      ! Print statistics
      do JK = 1, satNumber

        INTOTOBS = INTOT(JK)
        INTOTACC = INTOTOBS - INTOTRJF(JK) - INTOTRJP(JK)

        write(*,'(/////50("*"))')
        write(*,'(     50("*")/)')
        write(*,'(T5,"SUMMARY OF QUALITY CONTROL FOR ", &
         A8)') satelliteId(JK) 
        write(*,'(T5,"------------------------------------- ",/)')
        write(*,'( &
         "   TOTAL NUMBER OF INST  = ",I10,/ &
         " - TOTAL FULL REJECTS      = ",I10,/ &
         " - TOTAL PARTIAL REJECTS   = ",I10,/ &
         "   ------------------------------------",/ &
         "   TOTAL FULLY ACCEPTED    = ",I10,/)') &
          INTOTOBS, INTOTRJF(JK), INTOTRJP(JK), INTOTACC

        write(*,'(//,1x,114("-"))')
        write(*,'(t10,"|",t47,"REJECTION CATEGORIES")')
        write(*,'(" CHANNEL",t10,"|",105("-"))')
        write(*,'(t10,"|",15i7)') (JI,JI=1,JPMXREJ)
        write(*,'(1x,"--------|",105("-"))')
        do JJ = 1, MXCHN
           write(*,'(3X,I2,t10,"|",15I7)') JJ,(MREJCOD(JI,JJ,JK), &
                                      JI=1,JPMXREJ)
        end do
        write(*,'(1x,114("-"))')
      end do

      ! Print legend
      if (legendRejetAmsua) then
        print *, ' '
        print *, ' '
        print *, ' -----------------------------------------------------'
        print *, ' Definition of rejection categories:'
        print *, ' -----------------------------------------------------'
        print *, '  1 - topography reject'
        print *, '  2 - invalid land/sea qualifier'
        print *, '  3 - invalid terrain type'
        print *, '  4 - invalid field of view number'
        print *, '  5 - satellite zenith angle out of range '
        print *, '  6 - inconsistent field of view and sat. zenith angle'
        print *, '  7 - inconsistent land/sea qualifier and model mask'
        print *, '  8 - inconsistent terrain type and land/sea', &
                 ' qualifier/model ice (NOT doNE)'
        print *, '  9 - uncorrected radiance'
        print *, ' 10 - rejected by RTTOV'
        print *, ' 11 - radiance gross check failure'
        print *, ' 12 - cloud liquid water reject'
        print *, ' 13 - scattering index reject'
        print *, ' 14 - radiance residual rogue check failure'
        print *, ' 15 - rejection by channel selection'
        print *, ' -----------------------------------------------------'
        print *, ' ' 
      else if (legendRejetAtms) then 
        print *, ' '
        print *, ' '
        print *, ' -----------------------------------------------------'
        print *, ' Definition of rejection categories: '
        print *, ' -----------------------------------------------------'
        print *, '  1 - first bgckAtms program reject [bit 7]'
        print *, '  2 - topography reject'
        print *, '  3 - uncorrected radiance'
        print *, '  4 - innovation (O-P) based reject'
        print *, '  5 - rejection by channel selection'
        print *, ' -----------------------------------------------------'
        print *, ' '
        print *, ' QC2 REJECT numbers in Table 2 are for data that '
        print *, ' passed test 1 (data with QC flag bit 7 OFF)'
        print *, ' '
      end if

    end if

    return
  end subroutine mwbg_qcStats


  subroutine resetQcCases(RESETQC, KCHKPRF, globMarq)
    !:Purpose:        allumer la bit (6) indiquant que l'observation a un element
    !                 rejete par le controle de qualite de l'AO.
    !                 N.B.: si on est en mode resetqc, on remet le marqueur global a
    !                 sa valeur de defaut, soit 1024,  avant de faire la mise a jour.
    implicit none
    !Arguments:
    logical,              intent(in)     :: RESETQC       !reset the quality control flags before adding the new ones ?
    integer,              intent(in)     :: KCHKPRF(:)    !indicateur global controle de qualite tovs. Code:
    integer,              intent(inout)  :: globMarq(:)   !Marqueurs globaux  
    !Locals
    
    integer                              :: dataNum 
    integer                              :: dataIndex
    logical                              :: debug

    debug = mwbg_debug
    dataNum = size(globMarq)
    DO dataIndex = 1, dataNum
      IF (RESETQC) THEN
        globMarq(dataIndex) = 1024  
      ENDIF
      IF ( KCHKPRF(dataIndex).NE.0  ) THEN
        globMarq(dataIndex) = OR (globMarq(dataIndex),2**6)
      ENDIF
    ENDDO
    IF (debug) THEN
      write(*,*) ' KCHKPRF   = ', (KCHKPRF(dataIndex),dataIndex=1,dataNum)
      write(*,*) ' NEW FLAGS = ', (globMarq(dataIndex),dataIndex=1,dataNum)
    ENDIF

  end  subroutine resetQcCases

  subroutine mwbg_setTerrainTypeToSeaIce(GLINTRP, KTERMER, ITERRAIN)

    !:Purpose:       Dans les conditions suivantes:
    !                1) l'indicateur terre/mer indique l'ocean (ktermer=1),
    !                2) le "terrain type" est manquant (iterrain=-1),
    !                3) le modele indique de la glace (gl >= 0.01),
    !                on specifie "sea ice" pour le "terrain type" (iterrain=0).
    implicit none 
    !Arguments:
    real,                 intent(in)     :: GLINTRP(:)    !sea ice
    integer,              intent(in)     :: KTERMER(:)    !land sea qualifier
    integer,              intent(inout)  :: ITERRAIN(:)   !terrain type
    !Locals
    
    integer                              :: dataNum 
    integer                              :: dataIndex
    logical                              :: debug

    debug = mwbg_debug
    dataNum = size(ITERRAIN)

    IF ( mwbg_debug ) THEN
      WRITE(*,*) ' OLD TERRAIN TYPE = ', (ITERRAIN(dataIndex),dataIndex=1,dataNum)
      WRITE(*,*) ' KTERMER = ', (KTERMER(dataIndex),dataIndex=1,dataNum)
      WRITE(*,*) ' GLINTRP = ', (GLINTRP(dataIndex),dataIndex=1,dataNum)
    ENDIF
    DO dataIndex = 1, dataNum
      IF ( KTERMER (dataIndex) == 1 .and. ITERRAIN(dataIndex) == -1 .and. GLINTRP (dataIndex) >= 0.01 ) &
           ITERRAIN(dataIndex) = 0
    ENDDO
    IF ( mwbg_debug ) THEN
      WRITE(*,*) ' NEW TERRAIN TYPE = ', (ITERRAIN(dataIndex),dataIndex=1,dataNum)
    ENDIF
    
  end  subroutine mwbg_setTerrainTypeToSeaIce

  subroutine mwbg_addIntegerElementBurpBlock(burpBlock, burpElement, burpArr, burpChannelNum, burpLocationNum, burpNele)

    !:Purpose: This subroutine takes a block, the btyp and the Intger element 
    !          and that we want to be added to the given block
    !          The element can be an array of element

    implicit none 
    ! Arguments:
    type(BURP_BLOCK),           intent(inout)       :: burpBlock               ! burp report
    integer,                    intent(in)          :: burpElement             ! burp element 
    integer,                    intent(in)          :: burpArr(:)              ! burp integer array to add
    integer,                    intent(in)          :: burpLocationNum         ! nlocs
    integer,                    intent(in)          :: burpChannelNum          ! nchannels
    integer,                    intent(in)          :: burpNele 
    
    ! Locals
    integer                                         :: ref_blk 
    logical                                         :: debug
    integer                                         :: burpLocationIndex
    integer                                         :: dataIndex
    integer                                         :: burpChannelIndex
    integer                                         :: error                   ! error status

    debug = mwbg_debug
    debug = .true. 
    Call BURP_Resize_Block(burpBlock, ADD_NELE = 1, IOSTAT = error)
    if (error /= burp_noerr)  call abort()

    ! Add new elements CLW and scatw
    Call BURP_Set_Element(burpBlock,NELE_IND=burpNele,ELEMENT=burpElement,IOSTAT=error)
    call BURP_Encode_Block(burpBlock)   ! encode the element numbers in the block
        
    dataIndex = 1 
    do burpChannelIndex = 1,burpChannelNum  
      do burpLocationIndex  =1, burpLocationNum
        Call BURP_Set_Tblval(burpBlock,NELE_IND=burpNele,NVAL_IND=burpChannelIndex,NT_IND=burpLocationIndex,TBLVAL=burpArr(dataIndex),IOSTAT=error)
        dataIndex = dataIndex + 1        
      end do
    end do 
    Call BURP_Convert_Block(burpBlock)

  end subroutine mwbg_addIntegerElementBurpBlock

  subroutine mwbg_addRealElementBurpBlock(burpBlock, burpElement, burpArr, burpChannelNum, burpLocationNum, burpNele)

    !:Purpose: This subroutine takes a block, the btyp and the Real element 
    !          and that we want to be added to the given block
    !          The element can be an array of element

    implicit none 
    ! Arguments:
    type(BURP_BLOCK),           intent(inout)       :: burpBlock               ! burp report
    integer,                    intent(in)          :: burpElement             ! burp element 
    real,                       intent(in)          :: burpArr(:)              ! burp Real array to add
    integer,                    intent(in)          :: burpLocationNum         ! nlocs
    integer,                    intent(in)          :: burpChannelNum          ! nchannels
    integer,                    intent(in)          :: burpNele 
    
    ! Locals
    integer                                         :: ref_blk 
    logical                                         :: debug
    integer                                         :: burpLocationIndex
    integer                                         :: dataIndex
    integer                                         :: burpChannelIndex
    integer                                         :: error                   ! error status

    debug = mwbg_debug
    debug = .true. 
    Call BURP_Resize_Block(burpBlock, ADD_NELE = 1, IOSTAT = error)
    if (error /= burp_noerr)  call abort()

    ! Add new elements CLW and scatw
    Call BURP_Set_Element(burpBlock,NELE_IND=burpNele,ELEMENT=burpElement,IOSTAT=error)
    call BURP_Encode_Block(burpBlock)   ! encode the element numbers in the block
    
    dataIndex = 1 
    do burpChannelIndex = 1,burpChannelNum  
      do burpLocationIndex  =1, burpLocationNum
        Call BURP_Set_Rval(burpBlock,NELE_IND=burpNele,NVAL_IND=burpChannelIndex,NT_IND=burpLocationIndex,RVAL=burpArr(dataIndex),IOSTAT=error)
        dataIndex = dataIndex + 1        
      end do
    end do 
    Call BURP_Convert_Block(burpBlock)

  end subroutine mwbg_addRealElementBurpBlock


  subroutine mwbg_copyIntegerElementToBurpBlock(burpBlock, burpElement, burpArr, burpArrName, burpChannelNum, burpLocationNum)

    !:Purpose: This subroutine takes a block, the btyp and the Real element 
    !          and that we want to be added to the given block
    !          The element can be an array of element

    implicit none 
    ! Arguments:
    type(BURP_BLOCK),           intent(inout)       :: burpBlock               ! burp report
    integer,                    intent(in)          :: burpElement             ! burp element 
    integer,                    intent(in)          :: burpArr(:)              ! burp Real array to add
    character (*),              intent(in)          :: burpArrName             ! burp array name
    integer,                    intent(in)          :: burpLocationNum         ! nlocs
    integer,                    intent(in)          :: burpChannelNum          ! nchannels
    
    ! Locals
    integer                                         :: ref_blk 
    integer                                         :: indice 
    logical                                         :: debug
    integer                                         :: burpLocationIndex
    integer                                         :: burpChannelIndex
    integer                                         :: error                  
    integer                                         :: idata 

    debug = mwbg_debug
    debug = .true. 

    indice = BURP_Find_Element(burpBlock,burpElement,IOSTAT=error)
    if (error /= burp_noerr)  call abort()
    if ( indice > 0 ) then
      do burpChannelIndex = 1, burpChannelNum
        do burpLocationIndex = 1, burpLocationNum
            idata = burpArr(burpLocationIndex)
            Call BURP_Set_Tblval(burpBlock,indice,burpChannelIndex,burpLocationIndex,idata)
        end do
      end do
    else
      write(*,*)'Element ', burpArrName, ' MISSING in 3D block (btyp=5120).'
      call abort()
    end if

  end subroutine mwbg_copyIntegerElementToBurpBlock

  
  subroutine modifyBurpBktypAndWriteReport(burpReport, burpBlock, newBktyp, convertBlock)
    !
    !:Purpose:This subroutine takes a block, and its current bktyp 
    !          and then modify the bktyp then write the block to the 
    !          report

    implicit none 
    ! Arguments:
    type(BURP_BLOCK),           intent(inout)       :: burpBlock               ! burp block
    type(BURP_RPT),             intent(inout)       :: burpReport              ! burp report
    integer,                    intent(in)          :: newBktyp                ! burp element 
    logical,                    intent(in)          :: convertBlock            ! convert block yes or no
    
    ! Locals
    integer                                         :: error                  
    type(BURP_BLOCK)                                :: burpBlock_copy          ! burp block

    ! Remplacer le bloc et modifier le bktyp pour signifier "vu par AO".
    Call BURP_Set_Property(burpBlock,BKTYP=newBktyp,IOSTAT=error)
    burpBlock_copy = burpBlock
    Call BURP_Write_Block(burpReport,BLOCK = burpBlock_copy ,CONVERT_BLOCK=convertBlock,IOSTAT=error)
    if (error /= burp_noerr)  call abort()
    

  end subroutine modifyBurpBktypAndWriteReport


  SUBROUTINE mwbg_updateBurpAmsua(burpFileNameIn, ReportIndex, ETIKRESU, clw, scatw, &
                                  globMarq, KCHKPRF, KTERMER, ITERRAIN, & 
                                  RESETQC, IMARQ, burpFileNameout)
    !--------------------------------------------------------------------------------------
    !:Purpose:      Allumer les bits des marqueurs pour les tovs rejetes.
    !               Mettre a jour l'indicateur terre/mer qui a
    !               possiblement ete corrige pour la glace marine.
    !               Modifier le bktyp des donnees, marqueurs et (O-P) pourt
    !               signifier "vu par AO". 
    implicit none
    !Arguments:
    !
    character(len=90),    intent(in)     :: burpFileNameIn     !burp input Obs File
    integer,              intent(in)     :: ReportIndex        !report eportIndex
    character(len=9),     intent(in)     :: ETIKRESU           !resume report etiket
    real,                 intent(in)     :: clw(:)             !cloud liquid water path 
    real,                 intent(in)     :: scatw(:)           !scattering index 
    integer,              intent(inout)  :: globMarq(:)        !Marqueurs globaux  
    integer,              intent(in)     :: KCHKPRF(:)         !indicateur global controle de qualite tovs. Code:
                                                               !=0, ok,
                                                               !>0, rejet, 
    integer,              intent(in)     :: KTERMER(:)         !indicateur terre/mer 
    integer,              intent(inout)  :: ITERRAIN(:)        !indicateur du type de terrain 
                                                               !niveau de chaque canal  
    logical,              intent(in)     :: RESETQC            !reset the quality control flags before adding the new ones ?
    integer,              intent(in)     :: IMARQ(:)           !modified flag values from mwbg_tovCheckAmsua  
    character(len=90),    intent(in)     :: burpFileNameout    ! burp output Obs File

    !Locals:
    type(BURP_FILE), save                :: File_in
    type(BURP_FILE), save                :: File_out
    type(BURP_RPT)                       :: reportIn
    type(BURP_RPT)                       :: reportOut

    type(BURP_BLOCK)                     :: blk 
    type(BURP_BLOCK)                     :: blk_copy 

    logical, save                        :: ifFirstCall=.True.
    logical                              :: resumeReport 
    character(len=9)                     :: idStn
    integer, allocatable, save           :: adresses(:)        ! First Loop over all reports to get their adresses and save them
    integer, save                        :: nb_rpts
    integer, save                        :: nsize

    integer                              :: error
    integer                              :: ref_blk
    integer                              :: my_nt
    integer                              :: my_nval
    integer                              :: my_nele
    integer                              :: my_idtyp
    integer                              :: idtyp
    integer                              :: nlocs
    integer                              :: blat
    integer                              :: blon
    integer                              :: nblocs
    integer                              :: handle
    integer                              :: iun_burpin
    integer                              :: my_btyp
    integer                              :: my_bktyp 
    integer                              :: my_bfam
    integer                              :: new_bktyp
    integer                              :: indice 
    integer                              :: kk
    integer                              :: jj
    integer                              :: JI
    integer                              :: ipos
    integer                              :: idata 
    integer                              :: j

    if (ifFirstCall) then
      write(*,*) 'mwbg_updateBurpAmsua First Call : Initialisation'
      ! LOOP OVER ALL REPORTS OF THE INPUT FILE, APPLY PROCESSING, AND WRITE TO OUTPUT FILE.
      call mwbg_getBurpReportAdresses (burpFileNameIn, adresses)
      ! initialisation
      Call BURP_Init(File_in, F2= File_out, IOSTAT=error)
      ! ouverture du fichier burp d'entree et de sortie
      Call BURP_Init(reportIn, R2=reportOut, IOSTAT=error)
      ! Number of reports and maximum report size from input BURP file
      Call BURP_New(File_in,  FILENAME = burpFileNameIn,  MODE= FILE_ACC_READ,   IOSTAT= error)
      Call BURP_Get_Property(File_in, NRPTS=nb_rpts, IO_UNIT= iun_burpin)
      nsize = MRFMXL(iun_burpin)
      ! Add nsize to report size to accomodate modified (larger) data blocks
      nsize = nsize*3
      Call BURP_New(File_out, FILENAME= burpFileNameOut, MODE= FILE_ACC_CREATE, IOSTAT= error)
      ifFirstCall = .False.
    end if


    write(*,*) 'Get reportIn' 
    Call BURP_Get_Report(File_in, REPORT= reportIn, REF= adresses(ReportIndex), IOSTAT= error) 
    write(*,*) 'BURP_Get_Property reportIn'
    if (error /= burp_noerr) call mwbg_burpErrorHistory(file_in, reportIn) 
    Call BURP_Get_Property(reportIn,STNID=idStn,IDTYP=idtyp,ELEV=nlocs,LATI=blat,LONG=blon,NBLK=nblocs,HANDLE=handle)

    if ( idStn(1:2) .eq. ">>" ) then
      write(*,*) 'This is a resume Report: Copy ETIK and RETURN' 
      resumeReport = .True.
      ! change the header
      Call BURP_Set_Property(reportIn,STNID=ETIKRESU)  
      Call BURP_Write_Report(File_out,reportIn,IOSTAT=error)
      if (error /= burp_noerr) call mwbg_burpErrorHistory(File_in, reportIn, File_out, reportOut)

      if (reportIndex == nb_rpts) then
        write(*,*) ' reportIndex =  nb_rpts : END OF FILE'
        Call BURP_Free(File_in, F2 = File_out,IOSTAT=error)
        Call BURP_Free(reportIn, R2 = reportOut,IOSTAT=error)
        write(*,*) 'FREE FILE_IN and FILE_OUT'
      end if
      
      return 
    end if  

    write(*,*) ' Observation Report: Proceed to writing actions' 
    ! Create new report (reportOut) to contain modified blocks from reportIn
    Call BURP_New(reportOut, Alloc_Space = nsize,  IOSTAT=error)
    if (error /= burp_noerr) call mwbg_burpErrorHistory(file_in, reportIn, file_out, reportOut)
    
    write(*,*) 'BURP_INIT_REPORT_WRITE reportOut to File_out'
    ! initiliser pour ecriture a File_out
    Call BURP_INIT_Report_Write(File_out,reportOut,IOSTAT=error)
    if (error /= burp_noerr) call mwbg_burpErrorHistory(File_in, reportIn, File_out, reportOut)
    write(*,*) 'BURP_Copy_Header TO reportOut'
    !  copier le header du rapport 
    Call BURP_Copy_Header(TO= reportOut, FROM= reportIn)

    Call BURP_Init(blk, B2=blk_copy, IOSTAT=error)

    ! Read and modify the blocks in rpt and add them to reportOut
    ref_blk = 0
    BLOCKS: do

      ref_blk = BURP_Find_Block(reportIn,BLOCK= blk,SEARCH_FROM= ref_blk,IOSTAT= error)
      if (error /= burp_noerr) call abort()
      
      if (ref_blk < 0) Exit BLOCKS

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

        IF (mwbg_debug) THEN
          write(*,*) ' OLD FLAGS = ', (globMarq(JJ),JJ=1,my_nt)
        ENDIF

        call resetQcCases(RESETQC, KCHKPRF, globMarq)
        ! Remplacer les nouveaux marqueurs dans le tableau.
        call mwbg_copyIntegerElementToBurpBlock(blk, 55200, globMarq, "MARQUEURS_GLOBAUX", my_nval, my_nt)
        call modifyBurpBktypAndWriteReport(reportOut, blk, my_bktyp,.False.)        
      ! 2) Bloc info (general): bloc 3072
      !    Modifier les indicateurs terre/mer possiblement corriges pour la glace
      !    marine.
      !    Add new elements CLW and scatw
      else if (my_btyp == 3072) then
        ! Add new elements CLW

        call mwbg_addRealElementBurpBlock(blk, 13209, clw, 1, my_nt, my_nele+1)
        ! Add new elements scatw
        call mwbg_addRealElementBurpBlock(blk, 13208, scatw, 1, my_nt, my_nele+2)

        ! Remplacer les nouveaux indicateurs terre/mer dans le tableau.
        call mwbg_copyIntegerElementToBurpBlock(blk, 8012, KTERMER, "Land_Sea_INDICATOR", my_nval, my_nt)
        ! Update terrain type dans le tableau.
        call mwbg_copyIntegerElementToBurpBlock(blk, 13039, ITERRAIN, "TERRAIN_TYPE", my_nval, my_nt)
        ! finally copy blk
        call modifyBurpBktypAndWriteReport(reportOut, blk, my_bktyp,.False.)        

      ! 3) Bloc multi niveaux de radiances: bloc 9218, 9248, 9264.
      !    Modifier le bktyp pour signifier "vu par AO".
      else if ( (my_btyp == 9218 .or. my_btyp == 9248 .or. my_btyp ==9264) .and. &
                my_bfam == 0 ) then 

        ! Remplacer le bloc et modifier le bktyp pour signifier "vu par AO".
        call modifyBurpBktypAndWriteReport(reportOut, blk, my_bktyp+4, .True.)        
      ! 4) Bloc marqueurs multi niveaux de radiances: bloc 15362, 15392, 15408.
      !    Modifier les marqueurs de 13bits associes a chaque radiance.
      !    Modifier le bktyp pour signifier "vu par AO".
      else if ( (my_btyp == 15362 .or. my_btyp == 15392 .or. my_btyp == 15408) .and. &
                my_bfam == 0 ) then 

        ! extraire les marqueurs de 13bits des radiances; element 212163 (LEVEL 1B)
        indice = BURP_Find_Element(blk,212163,IOSTAT=error)
        ! Remplacer les nouveaux marqueurs dans le tableau.
        ipos = 0
        do kk =1, my_nt
          do j = 1, my_nval
            ipos = ipos + 1
            Call BURP_Set_Tblval(blk,indice,j,kk,IMARQ(ipos))
          enddo
        enddo
        ! Remplacer le bloc et modifier le bktyp pour signifier "vu par AO".
        call modifyBurpBktypAndWriteReport(reportOut, blk, my_bktyp+4, .True.)        
      ! 5) Bloc multi niveaux de residus de radiances (O-P): bloc 9322, 9226, 9258, 9274, bfam 14
      !    Modifier le bktyp pour signifier "vu par AO".
      else if ( (my_btyp == 9322 .or. my_btyp == 9226 .or. my_btyp == 9258 .or. &
                my_btyp == 9274) .and. my_bfam == 14 ) then 
        ! Remplacer le bloc et modifier le bktyp pour signifier "vu par AO".
        call modifyBurpBktypAndWriteReport(reportOut, blk, my_bktyp+4, .True.)
      ! OTHER BLOCK 
      else
        call modifyBurpBktypAndWriteReport(reportOut, blk, my_bktyp, .True.)

      end if

    enddo BLOCKS
    write(*,*) ' Update file'

    ! Write the modified report to the output file
    Call BURP_Write_Report(File_out,reportOut,IOSTAT=error)
    if (error /= burp_noerr) call mwbg_burpErrorHistory(File_in, reportIn, File_out, reportOut)
    write(*,*) ' New Report'

    if (reportIndex == nb_rpts) then
      write(*,*) ' reportIndex =  nb_rpts : END OF FILE'
      Call BURP_Free(File_in, F2 = File_out,IOSTAT=error)
      Call BURP_Free(reportIn, R2 = reportOut,IOSTAT=error)
      write(*,*) 'FREE FILE_IN and FILE_OUT'
    end if 

  end subroutine mwbg_updateBurpAmsua



  SUBROUTINE GRODY (ier, ni, tb23, tb31, tb50, tb53, tb89, &
                   tb23_P, tb31_P, tb50_P, tb53_P, tb89_P, &
                   pangl, plat, ilansea, ice, tpw, clw, clw_avg, &
                   rain, snow, scatl, scatw)
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
    !               tb23_P  - input  -  23Ghz brightness temperature from background (K)
    !               tb31_P  - input  -  31Ghz brightness temperature from background (K)
    !               tb50_P  - input  -  50Ghz brightness temperature from background (K)
    !               tb53_P  - input  -  53Ghz brightness temperature from background (K)
    !               tb89_P  - input  -  89Ghz brightness temperature from background (K)
    !               pangl   - input  -  satellite zenith angle (deg.)
    !               plat    - input  -  lalitude (deg.)
    !               ilansea - input  -  land/sea indicator (0=land;1=ocean)
    !               ice     - output -  sea ice concentration (0-100%)
    !               tpw     - output -  total precipitable water (0-70mm)
    !               clw     - output -  cloud liquid water (0-3mm)
    !               clw_avg - output -  averaged cloud liquid water from obs and 
    !                                   background (0-3mm)
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
    real dif285t23_P, dif285t31_P

    real tb23  (:)
    real tb31  (:)
    real tb50  (:)
    real tb53  (:)
    real tb89  (:)
    real tb23_P(:)
    real tb31_P(:)
    real tb50_P(:)
    real tb53_P(:)
    real tb89_P(:)
    real pangl (:)
    real plat  (:)
    real ice   (:)
    real tpw   (:)
    real clw   (:)
    real clw_P
    real clw_avg(:)
    real scatl (:)
    real scatw (:)

    data zmisgLocal   / -99.     /
    data epsilon /   1.E-30 /

    ! 1) Initialise output parameters:
    do i = 1, ni
      ice  (i) = zmisgLocal
      tpw  (i) = zmisgLocal
      clw  (i) = zmisgLocal
      clw_avg(i) = zmisgLocal
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
        dif285t23  =max(285.-t23,epsilon)
        dif285t23_P=max(285.-tb23_P(i),epsilon)
        dif285t31  =max(285.-t31,epsilon)
        dif285t31_P=max(285.-tb31_P(i),epsilon)

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

          !3.3) Cloud liquid water from obs (clw) and background state (clw_P):
          ! identify and remove sea ice
          if ( abslat .gt. 50.  .and. &
              df1    .gt.  0.0        ) then  
            clw(i) = zmisgLocal
            clw_avg(i) = zmisgLocal
          else
            a =  8.240 - (2.622 - 1.846*cosz)*cosz
            b =  0.754
            c = -2.265
            clw(i) = a + b*log(dif285t23) & 
                      + c*log(dif285t31)
            clw(i) = clw(i)*cosz           ! theoretical cloud liquid water (0-3mm)
            clw(i) = clw(i) - 0.03         ! corrected   cloud liquid water 
            clw(i) = min(3.,max(0.,clw(i)))   ! jh       

            clw_P = a + b*log(dif285t23_P) & 
                      + c*log(dif285t31_P)
            clw_P = clw_P*cosz           ! theoretical cloud liquid water (0-3mm)
            clw_P = clw_P - 0.03         ! corrected   cloud liquid water 
            clw_P = min(3.,max(0.,clw_P))   ! jh       

            ! averaged CLW from observation and background
            clw_avg(i) = 0.5 * (clw(i) + clw_P)
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

      if ( mwbg_DEBUG .and. i <= 100 ) then
        print *, 'GRODY: i,tb23(i),tb31(i),tb50(i),tb89(i),pangl(i),plat(i), &
                  ilansea(i) = ', &
                  i,tb23(i),tb31(i),tb50(i),tb89(i),pangl(i),plat(i), &
                  ilansea(i)
        print *, 'GRODY: ier(i),ice(i),tpw(i),clw(i),clw_avg(i),rain(i),snow(i)=', &
                  ier(i),ice(i),tpw(i),clw(i),clw_avg(i),rain(i),snow(i)
      endif

    enddo

    return
  end SUBROUTINE GRODY


  subroutine mwbg_readStatTovs(statFileName, instName, satelliteId)
    !:Purpose:       Lire les statistiques de l'erreur totale pour les TOVS.
    !
    implicit none
    !Arguments:
    character(*), intent(in) :: statFileName  ! fichier stats des TOVS
    character(*), intent(in) :: instName      ! Instrument Name
    character(len = 9), allocatable, intent(out) :: satelliteId(:)    ! Identificateur de satellite 

    !Locals

    integer :: iunStat
    integer :: ier 
    integer :: ISTAT 
    integer :: fnom 
    integer :: ICHN 
    integer :: satNumber 
    integer :: jpnsatIndex
    integer :: jpnchanIndex
    integer, parameter :: casesNum = 2
    integer :: satIndex 
    integer :: casesIndex 
    integer :: INDX, IPOS
    integer :: alloc_status 
    real*8  :: tovErrorIn(JPCH,2,JPNSAT)
    real    :: tovErrorStatus(JPCH,JPNSAT)
    integer :: channelInNum(JPNSAT)
    integer :: identifSatId(JPNSAT)
    real*8 :: ZDUM
    
    character(132) :: CLDUM
    character(17) :: CSATSTR
    character(12)  :: satelliteStatsType(casesNum)

    DATA satelliteStatsType     / 'Monitoring',  'Assimilation'  /  
    
    iunStat = 0
    ier = fnom(iunStat,statFileName,'SEQ+FMT',0)
    if (ier .LT. 0) then
      write (*,*) '" bgckMW: Problem opening ", &
                   "TOVS total error statistics file "', statFileName                
      call ABORT ()
    end if



    write(*,FMT=9000)
9000 FORMAT(//,10x,"-mwbg_readStatTovs: reading total error statistics" &
          ," required for TOVS processing")

    ! 1. Initialize
100  CONTINUE
    do jpnsatIndex = 1, JPNSAT
      NCHNA(jpnsatIndex) = 0
      channelInNum(jpnsatIndex) = 0
      identifSatId(jpnsatIndex) = 0
      do jpnchanIndex = 1, JPCH
        tovErrorIn(jpnchanIndex,1,jpnsatIndex) = 0.0
        tovErrorIn(jpnchanIndex,2,jpnsatIndex) = 0.0
        IUTILST (jpnchanIndex,jpnsatIndex) = 0
        tovErrorStatus(jpnchanIndex,jpnsatIndex) = 0.0
      end do
    end do

    ! 2. Open the file
200  CONTINUE
    ! .... not done here anymore, jh, august 2000
    write(*,*) 'mwbg_readStatTovs: reading total error statistics required for ', &
              'TOVS processing'

    write(*,'(20X,"ASCII dump of stats_tovs file: "//)')
    do jpnchanIndex = 1, 9999999
      read (iunStat,'(A)',ERR=900,END=400) CLDUM
      write(*,'(A)')   CLDUM
    end do

    ! 4. Read number of satellites
400  CONTINUE

    rewind(iunStat)
    read (iunStat,*,ERR=900)
    read (iunStat,*,ERR=900) satNumber
    read (iunStat,*,ERR=900)
    alloc_status = 0
    if (allocated(satelliteId)) deallocate(satelliteId)
    allocate(satelliteId(satNumber), stat=alloc_status)
    if (alloc_status /= 0) then
      write(*,*) 'ERROR - allocate(satelliteId) in mwbg_readStatTovs'
      call abort()
    endif

    ! Read the satellite identification, the number of channels,
    ! the observation errors and the utilization flags
500  CONTINUE

    do jpnsatIndex = 1, satNumber
      read (iunStat,*,ERR=900)
      read (iunStat,'(A)',ERR=900) CLDUM
      CSATSTR = TRIM(ADJUSTL(CLDUM))

      ! Get satellite (e.g. NOAA18) from satellite/instrument (e.g. NOAA18 AMSUA)
      INDX = INDEX(CSATSTR, instName)
      if ( INDX .GT. 3 ) THEN
        if ( INDEX(CSATSTR,'EOS-2') .GT. 0 ) THEN
          satelliteId(jpnsatIndex) = 'AQUA'
        else
          IPOS = INDX-2
          satelliteId(jpnsatIndex) = CSATSTR(1:IPOS)
        end if
      else
        write ( *,*) '(" mwbg_readStatTovs: Non- "', instName, &
                   ' instrument found in stats file!")'
        write ( *,'(A)' ) CLDUM
        call ABORT ()
      end if
      read (iunStat,*,ERR=900)
      read (iunStat,*,ERR=900) identifSatId(jpnsatIndex), channelInNum(jpnsatIndex)
      do jpnchanIndex = 1, 3
        read (iunStat,*,ERR=900)
      end do

      ! Set errors to ERRBGCK column values
      do jpnchanIndex = 1, channelInNum(jpnsatIndex)
        read (iunStat,*,ERR=900) ICHN, &
                  tovErrorIn(ICHN,1,jpnsatIndex), &
                  tovErrorIn(ICHN,2,jpnsatIndex), &
                  IUTILST (ICHN,jpnsatIndex), ZDUM
        TOVERRST(ICHN,jpnsatIndex) = tovErrorIn(ICHN,1,jpnsatIndex)
      end do
      read (iunStat,*,ERR=900)
    end do

510  CONTINUE
    ! 6. Print error stats for assimilated channels

600  CONTINUE

    write(*,'(//5X,"Total errors for TOVS data"/)') 
    do jpnsatIndex = 1, satNumber
      INDX = 0
      do jpnchanIndex = 1, JPCH
        if ( IUTILST(jpnchanIndex,jpnsatIndex) .NE. 0 ) THEN
          NCHNA(jpnsatIndex) = NCHNA(jpnsatIndex) + 1
          INDX = INDX + 1
          MLISCHNA(INDX,jpnsatIndex) = jpnchanIndex
        end if
      end do
      do casesIndex = 1, CasesNum
        write(*,'(/7X,"Satellite: ",A,5X,A)') &
          satelliteId(jpnsatIndex), satelliteStatsType(casesIndex)
        write(*,'(7X,"Channels   : ",30(T22,27I4/))') & 
         (MLISCHNA(jpnchanIndex,jpnsatIndex),jpnchanIndex=1,NCHNA(jpnsatIndex))
        write(*,'(7X,"Total errors: ",30(T22,27f4.1/))') &
         (tovErrorIn(MLISCHNA(jpnchanIndex,jpnsatIndex),casesIndex,jpnsatIndex), &
          jpnchanIndex=1,NCHNA(jpnsatIndex))
      end do
    end do
     
    ! 7. Close the file
700  CONTINUE
    ISTAT = FCLOS (iunStat)
    ! .... not done here anymore, jh, august 2000
    write(*,*) " SATID's = "
    do jpnsatIndex = 1, satNumber
      write(*,*) '  ', satelliteId(jpnsatIndex)
    end do

    return

    ! Read error
900   write ( *, '(" mwbg_readStatTovs: Problem reading ", &
                     "TOVS total error stats file ")' ) 
    ISTAT = FCLOS (iunStat)
    call ABORT ()

    return
  end subroutine mwbg_readStatTovs

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

    INTEGER MXSCANAMSU, MXSFCREJ, MXCH2OMPREJ,  MXTOPO
    PARAMETER  ( MXSCANAMSU = 96 )
    PARAMETER  ( MXSFCREJ   =  8 )
    PARAMETER  ( MXCH2OMPREJ=  6 )
    PARAMETER  ( MXTOPO     =  5 )

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
            IF ( mwbg_debug ) THEN
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
              IF ( mwbg_debug ) THEN
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
              IF ( mwbg_debug ) THEN
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
            IF ( mwbg_debug ) THEN
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
             IF ( mwbg_debug ) THEN
               write(*,*)STNID(2:9),'CHANNEL REJECT: ', &
                      ' OBS = ',JJ, &
                      ' CHANNEL= ',ICHN                  
             ENDIF
           ENDIF
        ENDDO
      ENDDO      
    endif

    IF ( mwbg_debug ) THEN
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

    IF ( mwbg_debug ) THEN
      write(*,*)'KCHKPRF = ',(KCHKPRF(JJ),JJ=1,KNT)
    ENDIF 

    RETURN
  END SUBROUTINE mwbg_tovCheckAtms


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
  END SUBROUTINE mwbg_qcStatsAtms


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


  SUBROUTINE mwbg_readStatTovsAtms(ILUTOV,INUMSAT,CSATID)
    !OBJET          Read stats_atms_assim file
    !
    !ARGUMENTS      ilutov  - input  -  unite logique du fichier stats des TOVS
    !               inumsat - output -  nombre de satellites
    !               csatid  - output -  identificateur de satellite 
    IMPLICIT NONE

    INTEGER ILUTOV, JI, JJ, JK, JL, JM, I, ICHN
    INTEGER INUMSAT, INDX, IPOS

    INTEGER NUMCHNIN(JPNSAT), ISATID(JPNSAT)

    REAL*8  TOVERRIN(JPCH,2,JPNSAT)
    REAL*8  ZDUM

    CHARACTER*132  CLDUM
    CHARACTER*17   CSATSTR

    CHARACTER*9   CSATID(JPNSAT)
    CHARACTER*12  CTYPSTAT(2)

    DATA CTYPSTAT     / 'Monitoring',  'Assimilation'  /  

    WRITE(*,*) 'mwbg_readStatTovsAtms: reading total error statistics required for TOVS processing'

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

    WRITE(*,*) 'ASCII dump of stats_tovs file:'
    DO JI = 1, 9999999
      READ (ILUTOV,'(A)',ERR=900,END=400) CLDUM
      WRITE(*,'(A)') CLDUM
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
        WRITE (*,*) 'mwbg_readStatTovsAtms: Non-ATMS instrument found in stats file!'
        WRITE (*,'(A)') CLDUM
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

    WRITE(*,*) 'Total errors for TOVS data' 
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
        WRITE(*,'(2(A),5X,A)') 'Satellite: ', CSATID(JL), CTYPSTAT(JK)
        WRITE(*,'(A,30(T22,27I4/))') 'Channels   : ', (MLISCHNA(JI,JL),JI=1,NCHNA(JL))
        WRITE(*,'(A,30(T22,27f4.1/))') 'Total errors: ', &
          (TOVERRIN(MLISCHNA(JI,JL),JK,JL), JI=1,NCHNA(JL))
      ENDDO
    ENDDO
     
700  CONTINUE

    RETURN

    ! Read error
900 WRITE (*,*) 'mwbg_readStatTovsAtms: Problem reading TOVS total error stats file'
    CALL ABORT ()

    RETURN
  END SUBROUTINE mwbg_readStatTovsAtms
  
  
  SUBROUTINE mwbg_getBurpReportAdresses (brp_in, adresses)

    !:Purpose: Initial scan of file to get number of reports and number of data locations.
    !          Store address of each report in array adresses(nb_rpts) for main REPORTS loop
 
    implicit none 
    !Arguments:
    character(len=90),              intent(in)      :: brp_in        ! burp file
    integer,          allocatable,  intent(out)     :: adresses(:)   ! Location of the reports

    !integer,           intent(out)     :: nsize      ! 3* Number of the largest report in file

    !Locals:
    type(BURP_FILE)                    :: File_in  ! observation burp file
    type(BURP_RPT)                     :: rpt_in    ! observation burp file
    integer                            :: nobs_tot      ! Number of obs. in burp file
    integer                            :: nb_rpts      ! Number of report in burp file
    integer                            :: ref_rpt
    integer                            :: ibrptime
    integer                            :: irun
    integer                            :: iun_burpin
    integer                            :: nlocs 
    integer                            :: nsize 
    integer                            :: error 
    integer                            :: alloc_status

    integer, parameter                 :: mxnt = 3000
    character(len=9)                   :: DataIdType
    character(len=20)                  :: opt_missing

    !nsize = 0
   
    ! initialisation
    Call BURP_Init(File_in, IOSTAT=error)
    Call BURP_Init(Rpt_in,  IOSTAT=error)
    ! ouverture du fichier burp d'entree et de sortie
    Call BURP_New(File_in, FILENAME= brp_in, MODE= FILE_ACC_READ, IOSTAT= error)
    ! Number of reports and maximum report size from input BURP file
    Call BURP_Get_Property(File_in, NRPTS=nb_rpts, IO_UNIT= iun_burpin)
    if (nb_rpts.le.1) then
      write(*,*) 'The input BURP file ''', trim(brp_in), ''' is empty!'
      stop
    end if
    nsize = MRFMXL(iun_burpin)
    write(*,*)
    write(*,*) 'Number of reports containing observations = ', nb_rpts-1
    write(*,*) 'Size of largest report = ', nsize
    write(*,*)

    ! Add nsize to report size to accomodate modified (larger) data blocks
    !nsize = nsize*3
    if(allocated(adresses)) deallocate(adresses)
    allocate(adresses(nb_rpts), stat=alloc_status)
  
    if (alloc_status /= 0) then
      write(*,*) 'ERROR - allocate(adresses(nb_rpts)). alloc_status =' , alloc_status
      call abort()
    endif
  
    adresses(:) = 0
    ref_rpt = 0
    nb_rpts = 0
    nobs_tot = 0
    do
      ref_rpt = BURP_Find_Report(File_in, REPORT= Rpt_in, SEARCH_FROM= ref_rpt, IOSTAT= error)
      if (error /= burp_noerr) call mwbg_burpErrorHistory(file_in, rpt_in)
      if (ref_rpt < 0) Exit
      Call BURP_Get_Property(Rpt_in,TEMPS=ibrptime,ELEV=nlocs,STNID=DataIdType,RUNN=irun)  
      ! ELEV= the number of locations in the data box (for grouped data) ==> nt in each block
      if ( DataIdType(1:2) .eq. ">>" ) then
        write(*,*) 'Type de fichier a l_entree = ',DataIdType 
        if (DataIdType .ne. ">>DERIALT") then
          write(*,*) 'WARNING - le type de fichier devrait etre >>DERIALT'
        endif
      else if (DataIdType(1:1) .eq. "^" ) then
        if ( nlocs > mxnt ) then
          write(*,*) 'ERROR: Number of locations (nlocs) in report ',nb_rpts+1, ' exceeds limit (mxnt)!'
          write(*,*) '       nlocs = ', nlocs
          write(*,*) '       mxnt  = ', mxnt
          call  mwbg_burpErrorHistory(File_in, Rpt_in)
        endif
        nobs_tot = nobs_tot + nlocs
      endif
      nb_rpts = nb_rpts+1
      adresses(nb_rpts) = ref_rpt
    
    end do
    write(*,*) ' Scan 1: Number of reports in input BURP file (nb_rpts) = ', nb_rpts
    write(*,*) '         Number of data locations (nobs_tot)             = ', nobs_tot
    ! if no reports ABORT
    if ( nb_rpts == 0 ) call mwbg_burpErrorHistory(File_in, Rpt_in)
    ! if no observations STOP
    if ( nobs_tot == 0 ) call mwbg_burpErrorHistory(File_in, Rpt_in)

  END SUBROUTINE  mwbg_getBurpReportAdresses


  subroutine mwbg_burpErrorHistory (burpFile, burpRpt, burpFile_opt, burpRpt_opt)

    ! :Purpose: Send and abort signal when there is a fail in reading burp file
    !          The history of all the actions in this attempt is listed
    !          Free the memory of the FIle and report instance

    implicit none
    ! Arguments 
    type(BURP_FILE),           intent(inout)                    :: burpFile  ! observation burp file
    type(BURP_RPT),            intent(inout)                    :: burpRpt   ! observation burp file
    type(BURP_FILE), optional, intent(inout)                    :: burpFile_opt  ! observation burp file
    type(BURP_RPT),  optional, intent(inout)                    :: burpRpt_opt   ! observation burp file

    write(*,*) BURP_STR_ERROR()
    write(*,*) "history"
    Call BURP_STR_ERROR_HISTORY()
    Call BURP_Free(burpFile)
    Call BURP_Free(burpRpt)
    if (present(burpFile_opt) .and. present(burpRpt_opt)) then
      Call BURP_Free(burpFile_opt)
      Call BURP_Free(burpRpt_opt)
    end if
    call abort()

  end subroutine mwbg_burpErrorHistory


  subroutine mwbg_readGeophysicFieldsAndInterpolate(instName, zlat, zlon, MTINTRP, MGINTRP, GLINTRP)

    implicit none

    !:Purpose: Reads Modele Geophysical variables and save for the first time
    !         TOPOGRAPHIE (MF ou MX):
    !             MF est la topographie filtree avec unites en metres (filtered ME).
    !             MX est la topographie filtree avec unites en m2/s2  (geopotential topography).
    !         Glace de Mer (GL)
    !         Masque Terre-Mer (MG)
    !         Then Interpolate Those variables to observation location
    !Arguments: 
    character(*),       intent(in)   :: instName       ! Instrument Name
    real,               intent(in)   :: zlat(:)        ! Obseravtion Lats
    real,               intent(in)   :: zlon(:)        ! Observation Lons
    real, allocatable,  intent(out)  :: MGINTRP(:)     ! Glace de mer interpolees au pt d'obs.
    real, allocatable,  intent(out)  :: MTINTRP(:)     ! topographie filtree (en metres) et interpolees
    real ,allocatable,  intent(out)  :: GLINTRP(:)     ! Glace de mer interpolees au pt d'obs.
  
    ! Locals:
    real, allocatable, save  :: GL(:)                  ! Modele Glace de Mer (GL)
    real, allocatable, save  :: MG(:)                  ! Modele Masque Terre-Mer (MG)
    real, allocatable, save  :: MT(:)                  ! Modele Topographie (MT)
    real,              save  :: TOPOFACT               ! Facteur x topo pour avoir des unites en metre
    logical,           save  :: ifFirstCall = .True.   ! If .True. we read GL, MT and MG
    integer,           save  ::  gdmt                  ! topo interpolation param
    integer,           save  ::  gdmg                  ! mask terre-mer interpolation param
    integer,           save  ::  gdgl                  ! glace interpolation param
    integer                  ::  gdllsval          
    integer                  :: IUNGEO
    logical                  :: readGlaceMask 
    logical                  :: debug 
    integer                  :: ier, irec 
    integer                  :: ezqkdef, ezsetopt
    integer                  :: FSTINF,FSTPRM,FCLOS
    integer                  :: FSTLIR,FSTFRM, FNOM, FSTOUV
    integer                  :: NI, NJ, NK, IG1, IG2, IG3, IG4
    integer                  :: IDUM,IDUM1,IDUM2,IDUM3,IDUM4
    integer                  :: IDUM5,IDUM6,IDUM7,IDUM8
    integer                  :: IDUM9,IDUM10,IDUM11,IDUM12,IDUM13
    integer                  :: IDUM14,IDUM15,IDUM16,IDUM17,IDUM18
    character(len=12)        :: ETIKXX
    character(len=4)         :: CLNOMVAR
    character(len=4)         :: NOMVXX
    character(len=2)         :: TYPXX 
    character(len=1)         :: GRTYP
    character(len=128)       :: fileGeoName 
    INTEGER                  :: NLAT
    INTEGER                  :: NLON
    INTEGER, PARAMETER       :: MXLON = 5
    INTEGER, PARAMETER       :: MXLAT = 5
    INTEGER, PARAMETER       :: MXELM = 40
    REAL,    PARAMETER       :: DLAT = 0.4
    REAL,    PARAMETER       :: DLON = 0.6
    REAL                     :: XLAT
    REAL                     :: XLON
    real, allocatable        :: ZLATBOX (:,:)
    real, allocatable        :: ZLONBOX (:,:)
    real, allocatable        :: MGINTBOX(:,:)
    real, allocatable        :: MTINTBOX(:,:)
    real, allocatable        :: GLINTBOX(:,:)
    integer                  :: dataIndex
    integer                  :: boxPointIndex
    integer                  :: latIndex
    integer                  :: lonIndex
    integer                  :: zlatNum
    integer                  :: zlonNum
    integer                  :: dataNum
    integer                  :: boxPointNum

    ! STEP 0: CHECK IF ZLAT AND ZLON ARE SAME DIMENSION
    zlatNum = size(zlat)
    zlonNum = size(zlon)
    if (zlatNum .ne. zlonNum) then
      write(*,*) 'ERREUR: OBSERVATION ZLAT and ZLON should have SAME LENGTH'
      call abort()
    else 
      dataNum = zlatNum
    end if

    ! STEP 1: READ MT, GL and MG from the FST FILE 
    debug = mwbg_debug 
    if (instName == 'ATMS') then
      fileGeoName = './fstgzmx'
      readGlaceMask = .False.
    else if (instName == 'AMSUA') then
      fileGeoName = './fstglmg'
      readGlaceMask = .True.
    end if 
    if(ifFirstCall) then
      IUNGEO = 0 
      IER = FNOM(IUNGEO,fileGeoName,'STD+RND+R/O',0)

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
        IF(ALLOCATED(MT)) DEALLOCATE(MT)
        ALLOCATE ( MT(NI*NJ), STAT=ier)
        IER = FSTLIR(MT,IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1, &
           ' ',CLNOMVAR)
      ENDIF
      
      IER = FSTPRM ( IREC, IDUM1, IDUM2, IDUM3, IDUM4, & 
          IDUM5, IDUM6, IDUM7, IDUM8, IDUM9, IDUM10,  &
          IDUM11, TYPXX, NOMVXX, ETIKXX, GRTYP, IG1, &
          IG2, IG3, IG4, IDUM12, IDUM13, IDUM14,  &
          IDUM15, IDUM16, IDUM17, IDUM18 )
       WRITE (*,*) ' GRILLE MT : ',grtyp,ni,nj, &
                ig1,ig2,ig3,ig4
      ier  = ezsetopt('INTERP_DEGREE','LINEAR')  
      gdmt = ezqkdef(ni,nj,grtyp,ig1,ig2,ig3,ig4,iungeo)

      if (readGlaceMask) then 
        ! MG
        IREC = FSTINF(IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1,' ','MG')
        IF (IREC .LT. 0) THEN
          WRITE (*,*) ' LE MASQUE TERRE-MER EST INEXISTANT' 
          CALL ABORT()
        ENDIF

        IF(ALLOCATED(MG)) DEALLOCATE(MG)
        ALLOCATE ( MG(NI*NJ), STAT=ier)
        IER = FSTLIR(MG,IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1,&
                 ' ','MG')

        IER = FSTPRM ( IREC, IDUM1, IDUM2, IDUM3, IDUM4, &
             IDUM5, IDUM6, IDUM7, IDUM8, IDUM9, IDUM10, &
             IDUM11, TYPXX, NOMVXX, ETIKXX, GRTYP, IG1,&
             IG2, IG3, IG4, IDUM12, IDUM13, IDUM14, &
             IDUM15, IDUM16, IDUM17, IDUM18 )
        WRITE (*,*) ' GRILLE MG : ',grtyp,ni,nj, &
                ig1,ig2,ig3,ig4
        ier  = ezsetopt('INTERP_DEGREE','LINEAR')  
        gdmg = ezqkdef(ni,nj,grtyp,ig1,ig2,ig3,ig4,iungeo)
        ! GL
        IREC = FSTINF(IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1,' ','GL')
        IF (IREC .LT. 0) THEN
          WRITE (*,*) 'LE CHAMP GLACE DE MER EST INEXISTANT' 
          CALL ABORT()
        ENDIF

        IF(ALLOCATED(GL)) DEALLOCATE(GL)
        ALLOCATE ( GL(NI*NJ), STAT=ier)
        IER = FSTLIR(GL,IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1, &
                 ' ','GL')

        IER = FSTPRM ( IREC, IDUM1, IDUM2, IDUM3, IDUM4, & 
             IDUM5, IDUM6, IDUM7, IDUM8, IDUM9, IDUM10, &
             IDUM11, TYPXX, NOMVXX, ETIKXX, GRTYP, IG1, &
             IG2, IG3, IG4, IDUM12, IDUM13, IDUM14, &
             IDUM15, IDUM16, IDUM17, IDUM18 )
        WRITE (*,*) ' GRILLE GL : ',grtyp,ni,nj, &
                ig1,ig2,ig3,ig4
        ier  = ezsetopt('INTERP_DEGREE','LINEAR')  
        gdgl = ezqkdef(ni,nj,grtyp,ig1,ig2,ig3,ig4,iungeo)
      else 
        gdgl = -1
        gdmg = -1
      end if 
      IER = FSTFRM(IUNGEO)
      IER = FCLOS(IUNGEO)
      ifFirstCall = .False. 
    end if

    ! STEP 3:  Interpolation de la glace et le champ terre/mer du modele aux pts TOVS.
    ! N.B.: on examine ces champs sur une boite centree sur chaque obs.
    boxPointNum = MXLAT*MXLON
    IF(ALLOCATED(ZLATBOX)) DEALLOCATE(ZLATBOX)
    ALLOCATE (ZLATBOX(boxPointNum, dataNum) , STAT=ier) 
    IF(ALLOCATED(ZLONBOX)) DEALLOCATE(ZLONBOX)
    ALLOCATE (ZLONBOX(boxPointNum, dataNum) , STAT=ier) 
    IF(ALLOCATED(MTINTBOX)) DEALLOCATE(MTINTBOX)
    ALLOCATE (MTINTBOX(boxPointNum, dataNum) , STAT=ier) 
    IF(ALLOCATED(GLINTBOX)) DEALLOCATE(GLINTBOX)
    ALLOCATE (GLINTBOX(boxPointNum, dataNum) , STAT=ier) 
    IF(ALLOCATED(MGINTBOX)) DEALLOCATE(MGINTBOX)
    ALLOCATE (MGINTBOX(boxPointNum, dataNum) , STAT=ier) 
    NLAT = (MXLAT-1)/2
    NLON = (MXLON-1)/2
    DO dataIndex = 1, dataNum
      boxPointIndex = 0
      DO latIndex = -NLAT, NLAT
        XLAT = ZLAT(dataIndex) +latIndex*DLAT
        XLAT = MAX(-90.0,MIN(90.0,XLAT))
        DO lonIndex = -NLON, NLON
          boxPointIndex = boxPointIndex + 1
          XLON = ZLON(dataIndex) +lonIndex*DLON
          IF ( XLON .LT. -180. ) XLON = XLON + 360.
          IF ( XLON .GT.  180. ) XLON = XLON - 360.
          IF ( XLON .lt.    0. ) XLON = XLON + 360.
           ZLATBOX(boxPointIndex,dataIndex) = XLAT
           ZLONBOX(boxPointIndex,dataIndex) = XLON
         ENDDO
      ENDDO
    ENDDO
    ier = ezsetopt('INTERP_DEGREE','LINEAR')
    ier = gdllsval(gdmt,mtintbox,mt,ZLATBOX,ZLONBOX,boxPointNum*dataNum)
    IF (ier .lt. 0) THEN
      WRITE(*,*) 'ERROR in the interpolation of MT'
      CALL ABORT()
    END IF
    if(readGlaceMask) then   
      ier = gdllsval(gdmg,mgintbox,mg,ZLATBOX,ZLONBOX,boxPointNum*dataNum)
      IF (ier .lt. 0) THEN
        WRITE(*,*) 'ERROR in the interpolation of MG'
        CALL ABORT()
      ENDIF
      ier = gdllsval(gdgl,glintbox,gl,ZLATBOX,ZLONBOX,boxPointNum*dataNum)
      IF (ier .lt. 0) THEN
        WRITE(*,*) 'ERROR in the interpolation of GL'
        CALL ABORT()
      END IF
    end if 

    IF(ALLOCATED(MTINTRP)) DEALLOCATE(MTINTRP)
    ALLOCATE (MTINTRP(dataNum) , STAT=ier) 
    IF(ALLOCATED(MGINTRP)) DEALLOCATE(MGINTRP)
    ALLOCATE (MGINTRP(dataNum) , STAT=ier) 
    IF(ALLOCATED(GLINTRP)) DEALLOCATE(GLINTRP)
    ALLOCATE (GLINTRP(dataNum) , STAT=ier) 
    DO dataIndex = 1, dataNum
      IF (DEBUG) THEN
        PRINT *, ' ------------------  '
        PRINT *, ' dataIndex = ', dataIndex
        PRINT *, '   '
        PRINT *, ' zlat,zlon = ', zlat(dataIndex), zlon(dataIndex)
        PRINT *, '   '
        PRINT *, ' ZLATBOX = '
        PRINT *,  (ZLATBOX(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
        PRINT *, ' ZLONBOX = '
        PRINT *,  (ZLONBOX(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
        PRINT *, ' MGINTBOX = '
        PRINT *,  (MGINTBOX(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
        PRINT *, ' MTINTBOX = '
        PRINT *,  (MTINTBOX(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
        PRINT *, ' GLINTBOX = '
        PRINT *,  (GLINTBOX(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
      END IF
      MGINTRP(dataIndex) = 0.0
      MTINTRP(dataIndex) = 0.0
      GLINTRP(dataIndex) = 0.0
      DO boxPointIndex=1,MXLAT*MXLON
        MTINTRP(dataIndex) = MAX(MTINTRP(dataIndex),MTINTBOX(boxPointIndex,dataIndex)/TOPOFACT)
        if(readGlaceMask) then      
          MGINTRP(dataIndex) = MAX(MGINTRP(dataIndex),MGINTBOX(boxPointIndex,dataIndex))
          GLINTRP(dataIndex) = MAX(GLINTRP(dataIndex),GLINTBOX(boxPointIndex,dataIndex))
        end if
      END DO
      IF (DEBUG) THEN
        PRINT *, ' MGINTRP = ', MGINTRP(dataIndex)
        PRINT *, ' MTINTRP = ', MTINTRP(dataIndex)
        PRINT *, ' GLINTRP = ', GLINTRP(dataIndex)
      END IF
    END DO
  end subroutine mwbg_readGeophysicFieldsAndInterpolate


  SUBROUTINE readBurpInteger (repIndex, burpRpt, burpBlkTypList, burpFam, burpEle, error, burpArr, &
                              burpArrName, burpLocationNum, burpChannelNum, abortIfMissing)
    !:Purpose: This subroutine takes the report, the family and the element 
    !          and read (out) the corresponding INTEGER Array from the current block
    !          In some cases, if the array does not exist, it will be filled with MSING

    implicit none 
    ! Arguments:
    integer,                    intent(in)          :: repIndex                ! report index
    type(BURP_RPT),             intent(in)          :: burpRpt                 ! burp report
    integer,                    intent(in)          :: burpBlkTypList(:)       ! burp block TYPE
    integer,                    intent(in)          :: burpFam                 ! burp family
    integer,                    intent(in)          :: burpEle                 ! burp Element num
    character (*),              intent(in)          :: burpArrName             ! burp array name
    logical      ,              intent(in)          :: abortIfMissing          ! abort if the array is missing
    integer, allocatable,       intent(out)         :: burpArr(:)              ! burp INTEGER array read
    integer,                    intent(out)         :: error                   ! error status
    integer,                    intent(out)         :: burpLocationNum          ! my_nt value to be returned
    integer,                    intent(out)         :: burpChannelNum         ! my_nval value to be returned
    
    ! Locals
    type(BURP_BLOCK)                                :: burpBlk
    integer                                         :: positionIndex 
    integer                                         :: burpLocationIndex
    integer                                         :: burpChannelIndex
    integer                                         :: burpBlkTypListIndex
    integer                                         :: burpReadIndice
    integer                                         :: burpNele 
    integer                                         :: ref_blk 
    integer                                         :: burpBlkTypListNum
    integer                                         :: alloc_status 
    logical                                         :: debug

    debug = mwbg_debug
    debug = .true.

    burpBlkTypListNum = size(burpBlkTypList)
    burpBlkTypListIndex = 1
    ref_blk = 0
    do while ((ref_blk <= 0 ) .and. (burpBlkTypListIndex <= burpBlkTypListNum))
      ref_blk = 0
      write(*,*) 'ref_blk = ', ref_blk
      write(*,*) 'burpFam = ', burpFam
      write(*,*) 'burpBlkTypList(burpBlkTypListIndex) = ', burpBlkTypList(burpBlkTypListIndex)
      write(*,*) ' error = ', error
      ref_blk = BURP_Find_Block(burpRpt, &
                    BLOCK       = burpBlk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = burpFam, &
                    BTYP        = burpBlkTypList(burpBlkTypListIndex), &
                    IOSTAT      = error)
      if (error /= burp_noerr) call abort()
      burpBlkTypListIndex = burpBlkTypListIndex + 1
    end do
    if (ref_blk < 0) then
      if ( abortIfMissing ) then
        write(*,*) ' ERREUR - Elements ', burpArrName, ' BLOCK (', burpEle ,') are missing in Report = ', repIndex
        call abort()
      else
        burpArr(:) = MPC_missingValue_R4
        return
      end if
    end if
    call BURP_Get_Property(burpBlk, &
                     NELE = burpNele, &
                     NT   = burpLocationNum, & 
                     NVAL = burpChannelNum, IOSTAT=error)
    if (error /= burp_noerr)  call abort()
    alloc_status = 0
    if(allocated(burpArr)) deallocate(burpArr)
    if(debug) then
      write(*,*)' VAR = ', burpArrName
      write(*,*)' DIMS: NLOCs = ', burpLocationNum, '  NVALs = ',  burpChannelNum
    end if 
    allocate(burpArr(burpLocationNum*burpChannelNum), stat=alloc_status)
    if( alloc_status /= 0)  then
      write(*,*) ' midas_bgckAtms: Memory allocation error  in readBurpReal Subroutine'
      call abort()
    endif
    burpReadIndice = BURP_Find_Element(burpBlk,burpEle,IOSTAT=error)
    if ( burpReadIndice > 0 ) then
      positionIndex = 0
      do burpLocationIndex = 1, burpLocationNum
        do burpChannelIndex = 1, burpChannelNum
          positionIndex = positionIndex + 1
          burpArr(positionIndex) = BURP_Get_Tblval(burpBlk,burpReadIndice,burpChannelIndex,burpLocationIndex,error)
        end do
      end do
    else if ( burpReadIndice < 0 ) then  
      if ( abortIfMissing ) then
        write(*,*) ' ERREUR - Elements ', burpArrName, ' data (', burpEle ,') are missing in DATA block! Report = ', repIndex
        call abort()
      else
        burpArr(:) = MPC_missingValue_INT
      end if
    end if
  END SUBROUTINE readBurpInteger


  SUBROUTINE readBurpReal (repIndex, burpRpt, burpBlkTypList, burpFam, burpEle, error, burpArr, &
                           burpArrName, burpLocationNum, burpChannelNum, abortIfMissing)
    !:Purpose: This subroutine takes the report, the family and the element 
    !          and read (out) the corresponding REAL Array from the current block
    !          In some cases, if the array does not exist, it will be filled with MSING

    implicit none 
    ! Arguments:
    integer,                    intent(in)          :: repIndex                ! report index
    type(BURP_RPT),             intent(in)          :: burpRpt                 ! burp report
    integer,                    intent(in)          :: burpBlkTypList(:)       ! burp block TYPE
    integer,                    intent(in)          :: burpFam                 ! burp family
    integer,                    intent(in)          :: burpEle                 ! burp Element num
    character (*),              intent(in)          :: burpArrName             ! burp array name
    logical      ,              intent(in)          :: abortIfMissing          ! abort if the array is missing
    real, allocatable,          intent(out)         :: burpArr(:)              ! burp REAL array read
    integer,                    intent(out)         :: error                   ! error status
    integer,                    intent(out)         :: burpLocationNum          ! my_nt value to be returned
    integer,                    intent(out)         :: burpChannelNum         ! my_nval value to be returned
    
    ! Locals
    type(BURP_BLOCK)                                :: burpBlk
    integer                                         :: positionIndex 
    integer                                         :: burpLocationIndex
    integer                                         :: burpChannelIndex
    integer                                         :: burpBlkTypListIndex
    integer                                         :: burpReadIndice
    integer                                         :: burpNele 
    integer                                         :: ref_blk 
    integer                                         :: burpBlkTypListNum
    integer                                         :: alloc_status 
    logical                                         :: debug

    debug = mwbg_debug
    debug = .true. 

    burpBlkTypListNum = size(burpBlkTypList)
    burpBlkTypListIndex = 1
    ref_blk = 0
    do while ((ref_blk <= 0 ) .and. (burpBlkTypListIndex <= burpBlkTypListNum))
      ref_blk = 0
      write(*,*) 'ref_blk = ', ref_blk
      write(*,*) 'burpFam = ', burpFam
      write(*,*) 'burpBlkTypList(burpBlkTypListIndex) = ', burpBlkTypList(burpBlkTypListIndex)
      write(*,*) ' error = ', error
      ref_blk = BURP_Find_Block(burpRpt, &
                    BLOCK       = burpBlk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = burpFam, &
                    BTYP        = burpBlkTypList(burpBlkTypListIndex), &
                    IOSTAT      = error)
      if (error /= burp_noerr) call abort()
      burpBlkTypListIndex = burpBlkTypListIndex + 1
    end do
    if (ref_blk < 0) then
      if ( abortIfMissing ) then
        write(*,*) ' ERREUR - Elements ', burpArrName, ' BLOCK (', burpEle ,') are missing in Report = ', repIndex
        call abort()
      else
        burpArr(:) = MPC_missingValue_R4
        return
      end if
    end if
    call BURP_Get_Property(burpBlk, &
                     NELE = burpNele, &
                     NT   = burpLocationNum, & 
                     NVAL = burpChannelNum, IOSTAT=error)
    if (error /= burp_noerr)  call abort()
    alloc_status = 0
    if(allocated(burpArr)) deallocate(burpArr)
    if(debug) then
      write(*,*)' VAR = ', burpArrName
      write(*,*)' DIMS: NLOCs = ', burpLocationNum, '  NVALs = ',  burpChannelNum
    end if 
    allocate(burpArr(burpLocationNum*burpChannelNum), stat=alloc_status)
    if( alloc_status /= 0)  then
      write(*,*) ' midas_bgckAtms: Memory allocation error  in readBurpReal Subroutine'
      call abort()
    endif
    burpReadIndice = BURP_Find_Element(burpBlk,burpEle,IOSTAT=error)
    if ( burpReadIndice > 0 ) then
      positionIndex = 0
      do burpLocationIndex = 1, burpLocationNum
        do burpChannelIndex = 1, burpChannelNum
          positionIndex = positionIndex + 1
          burpArr(positionIndex) = BURP_Get_Rval(burpBlk,burpReadIndice,burpChannelIndex,burpLocationIndex,error)
        end do
      end do
    else if ( burpReadIndice < 0 ) then  
      if ( abortIfMissing ) then
        write(*,*) ' ERREUR - Elements ', burpArrName, ' data (', burpEle ,') are missing in DATA block! Report = ', repIndex
        call abort()
      else
        burpArr(:) = MPC_missingValue_R4
      end if
    end if
  END SUBROUTINE readBurpReal


  subroutine mwbg_getData(burpFileNameIn, reportIndex, ISAT, zenith, ilq, itt, &
                          zlat, zlon, ztb, biasCorr, ZOMP, scanpos, nvalOut, &
                          ntOut, qcflag1, qcflag2, ican, icanomp, IMARQ, IORBIT, &
                          globMarq, resumeReport, iFlastReport, InstName, STNID)
    !--------------------------------------------------------------------------------------
    !:Purpose:   This routine extracts the needed data from the blocks in the report:
    !             rpt              = report

    ! NOTE:  reportIndex = report number (from MAIN program) **** DO NOT MODIFY ****
    !       kk = variable for loops over locations (nt)
    !        j = variable for loops over nval (nval = 1 or nchanAtms)
    !Arguments:
    !
    character(len=90),    intent(in)     :: burpFileNameIn
    integer,              intent(in)     :: reportIndex    ! report index
    !type(BURP_RPT),       intent(in)     :: rpt           ! report
    character(*),         intent(in)     :: InstName       ! Instrument Name
    integer, allocatable, intent(out)    :: ISAT(:)        ! satellite identifier
    real   , allocatable, intent(out)    :: zenith(:)      ! satellite zenith angle (btyp=3072,ele=7024) 
    integer, allocatable, intent(out)    :: ilq(:)         ! land/sea qualifier     (btyp=3072,ele=8012)
    integer, allocatable, intent(out)    :: itt(:)         ! terrain-type (ice)     (btyp=3072,ele=13039)
    real   , allocatable, intent(out)    :: zlat(:)        ! latitude values (btyp=5120,ele=5002)
    real   , allocatable, intent(out)    :: zlon(:)        ! longitude values (btyp=5120,ele=6002)
    real   , allocatable, intent(out)    :: ztb(:)         ! brightness temperature (btyp=9248/9264,ele=12163) 
    real   , allocatable, intent(out)    :: biasCorr(:)    ! bias correction 
    real   , allocatable, intent(out)    :: ZOMP(:)        ! OMP values
    integer, allocatable, intent(out)    :: scanpos(:)     ! scan position (fov)    (btyp=3072,ele=5043)
    integer,              intent(out)    :: nvalOut        ! number of channels     (btyp=9248/9264)
    integer,              intent(out)    :: ntOut          ! number of locations    (btyp=5120,etc.)
    integer, allocatable, intent(out)    :: qcflag1(:,:)   ! flag values for btyp=3072 block ele 033078, 033079, 033080
    integer, allocatable, intent(out)    :: qcflag2(:)     ! flag values for btyp=9248 block ele 033081      
    integer, allocatable, intent(out)    :: ican(:)        ! channel numbers btyp=9248 block ele 5042 (= 1-22)
    integer, allocatable, intent(out)    :: icanomp(:)     ! omp channel numbers btyp= block ele  (= 1-22)
    integer, allocatable, intent(out)    :: IMARQ(:)       ! data flags
    integer, allocatable, intent(out)    :: IORBIT(:)      ! orbit number
    integer, allocatable, intent(out)    :: globMarq(:)    ! global Marqueur Data
    logical,              intent(out)    :: resumeReport   ! True if resume Report is read
    logical,              intent(out)    :: ifLastReport   ! True if last Report is read
    character(*),         intent(out)    :: STNID          ! Platform Name
    ! Locals
    type(BURP_FILE), save                :: File_in
    type(BURP_RPT)                       :: reportIn
    integer, allocatable                 :: qcflag1FirstColomn(:)
    integer, allocatable                 :: qcflag1SecondColomn(:)
    integer, allocatable                 :: qcflag1ThirdColomn(:)
    integer                              :: allocNtOut
    integer                              :: allocNvalOut
    integer                              :: burpLocationNum
    integer                              :: burpChannelNum
    integer                              :: eleChannel     
    integer                              :: eleDataQcFlag   
    integer                              :: eleQcFlag1(3) 
    integer                              :: idtyp
    integer                              :: nlocs 
    integer                              :: blat 
    integer                              :: blon 
    integer                              :: iun_burpin
    integer                              :: nblocs 
    integer                              :: handle 
    character(len=9)                     :: idStn0
    character(len=20)                    :: opt_missing
    logical, save                        :: ifFirstCall = .True.  ! If .True., will open file_in and do some work then set to .False.
    integer, allocatable, save           :: adresses(:)              ! First Loop over all reports to get their adresses and save them
    integer, save                        :: nb_rpts
    
    ! 0) Logical block to assign some attributs depending on InstName
    if (InstName == 'ATMS') then
      eleChannel       = 5042
      eleDataQcFlag    = 33081
      eleQcFlag1(1)    = 33078 
      eleQcFlag1(2)    = 33079 
      eleQcFlag1(3)    = 33080 
    else if (InstName == 'AMSUA') then
      eleChannel       = 2150
      eleDataQcFlag    = 33032
      eleQcFlag1(:)    = -1 
    else
      write(*,*) 'ERREUR - Instrument Name not Recognized '
      call abort()
    end if


    if (ifFirstCall) then
      write(*,*) 'mwbg_getData First Call : Initialisation'
      ! Set BURP "missing value" for reals
  ! Set BURP "missing value" for reals
      opt_missing = 'MISSING'
      Call BURP_Set_Options(REAL_OPTNAME=opt_missing,REAL_OPTNAME_VALUE=zmisg)

      !opt_missing = 'MISSING'
      !Call BURP_Set_Options(REAL_OPTNAME=opt_missing,REAL_OPTNAME_VALUE=MPC_missingValue_R4)
      ! LOOP OVER ALL REPORTS OF THE INPUT FILE, APPLY PROCESSING, AND WRITE TO OUTPUT FILE.
      call mwbg_getBurpReportAdresses (burpFileNameIn, adresses)
      ! initialisation
      Call BURP_Init(File_in, IOSTAT=error)
      ! ouverture du fichier burp d'entree et de sortie
      Call BURP_Init(reportIn,  IOSTAT=error)
      ! Number of reports and maximum report size from input BURP file
      Call BURP_New(File_in,  FILENAME = burpFileNameIn,  MODE= FILE_ACC_READ,   IOSTAT= error)
      Call BURP_Get_Property(File_in, NRPTS=nb_rpts, IO_UNIT= iun_burpin)
      write(*,*)
      write(*,*) 'Number of reports containing observations = ', nb_rpts-1
      write(*,*)
      ifFirstCall = .False.
    end if
    if (nb_rpts.lt.1) then
      write(*,*) 'The input BURP file ''', trim(burpFileNameIn), ''' is empty!'
      call abort()
    end if    
    ! last report check
    ifLastReport = reportIndex == nb_rpts
    
    Call BURP_Get_Report(File_in, REPORT= reportIn, REF= adresses(reportIndex), IOSTAT= error) 
    if (error /= burp_noerr) call mwbg_burpErrorHistory(file_in, reportIn) 
    Call BURP_Get_Property(reportIn,STNID=idStn0,IDTYP=idtyp,ELEV=nlocs,LATI=blat,LONG=blon,NBLK=nblocs,HANDLE=handle)
    if ( idStn0(1:2) .eq. ">>" ) then
      resumeReport = .True.
      write(*,*) 'Resume Report: Will return to next Repport'
      return 
    else 
      resumeReport = .False.
      write(*,*) 'Resume Report = ', resumeReport
      STNID = idStn0 
    end if
    
    write(*,*) 'NLOCs = ', nlocs
    write(*,*) 'report Index = ', reportIndex

    !  Get OMP data from the DATA block     BTYP =  9322 or 9226 or 9258 or 9274 and bfma = 14
    call readBurpReal (reportIndex, reportIn, (/9322,9226,9258,9274/), 14, 12163, error, ZOMP, 'Omp_Data', &
                   burpLocationNum, burpChannelNum, abortIfMissing = .FALSE.) 
    
    if ( ALL(ZOMP(:) == MPC_missingValue_R4 )) then
      return
    end if

    call readBurpInteger (reportIndex, reportIn, (/9322,9226,9258,9274/), 14, eleChannel, error, ICANOMP, 'OMP_Channels', &
                          burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 

    !  Get the lat,lon from time/location block    BTYP = 5120  (also get nt)
    call readBurpReal(reportIndex, reportIn, (/5120/), 0, 5002, error, zlat, 'LAT', &
                      burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 

    call readBurpReal(reportIndex, reportIn, (/5120/), 0, 6002, error, zlon, 'LON', &
                      burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 

    !  Get info elements from the INFO block   BTYP = 3072
    call readBurpInteger(reportIndex, reportIn, (/3072/), 0, 1007, error, ISAT, 'Sat_Identifier', &
                         burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
    call readBurpInteger(reportIndex, reportIn, (/3072/), 0, 5040, error, IORBIT, 'Orbit_Number', &
                         burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
    call readBurpReal(reportIndex, reportIn, (/3072/), 0, 7024, error, ZENITH, 'Zenith_Angle', &
                      burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
    call readBurpInteger(reportIndex, reportIn, (/3072/), 0, 8012, error, ILQ, 'LandSea_Qualifier', &
                         burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
    call readBurpInteger(reportIndex, reportIn, (/3072/), 0, 13039, error, ITT, 'Terrain_Type', &
                         burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
    call readBurpInteger(reportIndex, reportIn, (/3072/), 0, 5043, error, SCANPOS, 'Scan_Position', &
                         burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 

    !  Get info elements 33078 33079 33080 if needed
    if (ALL(eleQcFlag1 /= -1)) then
      call readBurpInteger(reportIndex, reportIn, (/3072/), 0, eleQcFlag1(1), error, qcflag1FirstColomn, 'Geoloc_Quality_QcFlag', &
                           burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
      call readBurpInteger(reportIndex, reportIn, (/3072/), 0, eleQcFlag1(2), error, qcflag1SecondColomn, 'Granule_Level_QcFlag', &
                           burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
      call readBurpInteger(reportIndex, reportIn, (/3072/), 0, eleQcFlag1(3), error, qcflag1ThirdColomn, 'Scan_Level_QcFlag', &
                           burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
      if(allocated(qcflag1)) deallocate(qcflag1)
      allocate(qcflag1(size(qcflag1FirstColomn),3))
      qcflag1(:,1) = qcflag1FirstColomn
      qcflag1(:,2) = qcflag1SecondColomn
      qcflag1(:,3) = qcflag1ThirdColomn
    end if

    !  Get data from the DATA block     BTYP = 9248 or 9264    (also get nval = nchanAtms)
    call readBurpReal (reportIndex, reportIn, (/9248,9264/), 0, 12163, error, ztb, 'Tb_data', &
                       burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
    nvalOut = burpChannelNum    ! set nvalOut (#channels) for MAIN program
    ntOut = burpLocationNum     ! set ntOut (#locations) for MAIN program
    call readBurpReal (reportIndex, reportIn, (/9248,9264/), 0, 12233, error, biasCorr, 'Bias_Corr_data', &
                       burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
    call readBurpInteger (reportIndex, reportIn, (/9248,9264/), 0, eleChannel, error, ICAN, 'Channel_Numbers', &
                          burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
    call readBurpInteger (reportIndex, reportIn, (/9248,9264/), 0, eleDataQcFlag, error, qcflag2, 'Data_level_Qc_Flag', &
                          burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 

    !  Bloc marqueurs multi niveaux de radiances: bloc 15362, 15392, 15408.
    call readBurpInteger (reportIndex, reportIn, (/15362,15392,15408/), 0, 212163, error, IMARQ, 'IMARQ', &
                          burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.)

    !  Bloc info 3d: bloc 5120.
    call readBurpInteger (reportIndex, reportIn, (/5120/), 0, 55200, error, globMarq, 'Marqueurs_Gobaux', &
                          burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.)

  end subroutine mwbg_getData



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


  function calcStateDepObsErr_r4(clwThresh1,clwThresh2,sigmaThresh1,sigmaThresh2,clw_avg) result(sigmaObsErrUsed)
    implicit none
    real :: clwThresh1
    real :: clwThresh2
    real :: sigmaThresh1
    real :: sigmaThresh2
    real :: clw_avg
    real :: sigmaObsErrUsed

    if ( clw_avg <= clwThresh1 ) then
      sigmaObsErrUsed = sigmaThresh1
    else if ( clw_avg >  clwThresh1 .and. & 
                  clw_avg <= clwThresh2 ) then
      sigmaObsErrUsed = sigmaThresh1 + &
                      (sigmaThresh2 - sigmaThresh1) / &
                      (clwThresh2 - clwThresh1) * &
                      (clw_avg - clwThresh1) 
    else
      sigmaObsErrUsed = sigmaThresh2
    end if

  end function calcStateDepObsErr_r4


end module bgckmicrowave_mod

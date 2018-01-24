!-------------------------------------- LICENCE BEGIN ------------------------------------
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

!--------------------------------------------------------------------------
!! MODULE varqc (prefix="")
!!
!! *Purpose*: Procedures related to variational quality control including
!!            hard-coded values that determine how quickly the observation
!!            weight begins to be reduced
!!
!--------------------------------------------------------------------------
module varqc_mod
  use MathPhysConstants_mod
  use EarthConstants_mod
  use codtyp_mod
  use bufr
  use obsSpaceData_mod
  use columnData_mod
  use rmatrix_mod ,only : rmat_lnondiagr
  use varNameList_mod
  implicit none
  save
  private

  ! public procedures
  public :: vqc_setup, vqc_tl, vqc_ad, vqc_listrej


  contains

  subroutine vqc_setup(obsSpaceData)
    !
    !s/r vqc_setup - SET CERTAIN PARAMETERS FOR THE ASYMMETRIC
    !                CHECK AND FOR VARIATIONAL QUALITY CONTROL
    !
    !Author  : B. BRASNETT CMDA   JUNE 1999
    !Revision:
    !          S. Pellerin ARMA/SMC Nov. 2002
    !             . Elemination of nincrem variable
    !          R. Sarrazin/B. Brasnett CMC March 2004
    !             . tighten rejection on satwinds
    !          J. Hale CMC Sept. 2005
    !             . added MHS (codtyp=182).
    !          S. Macpherson ARMA/CMC Sept. 2007
    !             . add parameters for GB-GPS ZTD
    !          S. Macpherson ARMA/CMC March 2013
    !             . modified GPS ZTD parameters to increase QC-Var rejections
    !          Y.J. Rochon ARQI/AQRD, Feb 2015
    !             . Added ZACH and ZDCH for chemical constituents
    !             . Added call to bufr_IsAtmosConstituent(ITYP)
    !               with addition of 'use varNamelist_mod'
    !
    implicit none
    type(struct_obs) :: obsSpaceData

    integer jdata,kindic,iter,jjo,idata,idatend,idburp
    integer ityp,iass,ifld,iother,jj,istyp,ilev
    real(8) zagz,zahu,zatt,zduv,zdgz,zdtt,zdhu,zauv,zslev,zapn,zdpn
    real(8) zabt,zdbt,zabtb,zdbtb,zach,zdch
    real(8) zlev,zjo,zval,zspdo,zspdf,zofcst,zoval,zdiff,zaasym,zoer
    real(8) zfcst,zlat,zlon,zprior,zaps,zdps,zauvra,zattra,zattsym
    logical llok
    real(8) zazd, zdzd
    !
    !_____prior probabilities for winds:ZAUV
    !     prior probabilities for scalar variables: zagz and zahu
    !     standard deviation multiple for background check for winds: zduv
    !     standard deviation multiple for background check for heights: zdgz
    !     standard deviation multiple for background check for humidity: zdhu
    !
    zagz   = 1.d-12
    zatt   = 5.d-2
!    zattra = 1.d-2
    zattra = 0.005d0
!    zauv   = 0.01d0
    zauv   = 0.02d0
!    zauvra = 1.d-3
    zauvra = 1.d-5
!    zahu   = 0.01d0
    zahu   = 0.05d0
    zapn = 1.d-4
    zaps = 1.d-4
    zabt   = 1.0d-1
    zabtb  = 1.0d-1
    zattsym = 1.d-1
    zazd = 2.0d-2
    zach = 1.0d-3
    zduv = 5.d0
    zdgz = 5.d0
    zdtt = 5.d0
    zdhu = 5.d0
    zdpn = 5.d0
    zdps = 5.d0
    zdbt = 3.d0
    zdbtb = 3.d0
    zdzd = 3.d0
    zdch = 10.d0

    DO JJO = 1, obs_numheader(obsSpaceData)
       IDATA = obs_headElem_i(obsSpaceData,OBS_RLN,JJO)
       IDATEND = obs_headElem_i(obsSpaceData,OBS_NLV,JJO) + IDATA - 1
       IDBURP = obs_headElem_i(obsSpaceData,OBS_ITY,JJO)
       ZLAT = obs_headElem_r(obsSpaceData,OBS_LAT,JJO)*MPC_DEGREES_PER_RADIAN_R8
       ZLON = obs_headElem_r(obsSpaceData,OBS_LON,JJO)*MPC_DEGREES_PER_RADIAN_R8
       DO JDATA = IDATA, IDATEND
          ITYP = obs_bodyElem_i(obsSpaceData,OBS_VNM,JDATA)
          IASS = obs_bodyElem_i(obsSpaceData,OBS_ASS,JDATA)
          ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,JDATA)*MPC_MBAR_PER_PA_R8
          ZOER = obs_bodyElem_r(obsSpaceData,OBS_OER,JDATA)
          ZVAL = obs_bodyElem_r(obsSpaceData,OBS_VAR,JDATA)
          ZFCST= ZVAL - obs_bodyElem_r(obsSpaceData,OBS_OMP,JDATA)

          IF (ITYP .EQ. BUFR_NETS .OR. ITYP .EQ. BUFR_NEPS .OR.  &
              ITYP .EQ. BUFR_NEPN .OR. ITYP .EQ. BUFR_NESS .OR.  &
              ITYP .EQ. BUFR_NEUS .OR. ITYP .EQ. BUFR_NEVS .OR.  &
              ITYP .EQ. BUFR_NEZD) THEN
             LLOK = (obs_bodyElem_i(obsSpaceData,OBS_ASS,JDATA) .EQ. 1)
          ELSE
             LLOK = (IASS .EQ. 1) .AND. ((obs_bodyElem_i(obsSpaceData,OBS_XTR,JDATA) .EQ.0) &
                .OR. ((obs_bodyElem_i(obsSpaceData,OBS_XTR,JDATA) .EQ. 2) .AND. &
                      (obs_bodyElem_i(obsSpaceData,OBS_VNM,JDATA) .EQ. BUFR_NEGZ)))
          ENDIF

          IF (LLOK) THEN
             IF (ITYP .EQ. BUFR_NEUU .OR. ITYP .EQ. BUFR_NEVV .OR.  &
                 ITYP .EQ. BUFR_NEUS .OR. ITYP .EQ. BUFR_NEVS) THEN
                ZAASYM = 1.0d0
                IOTHER = -1
                IF (ITYP .EQ. BUFR_NEUU .OR. ITYP .EQ. BUFR_NEUS) THEN
                  DO JJ=IDATA,JDATA
                    ISTYP = obs_bodyElem_i(obsSpaceData,OBS_VNM,JJ)
                    ZSLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,JJ)*MPC_MBAR_PER_PA_R8
                    IF ((ISTYP .EQ. BUFR_NEVV .OR. ISTYP .EQ. BUFR_NEVS)  &
                         .AND. ZLEV .EQ. ZSLEV) THEN
                      IOTHER = JJ
                    ENDIF
                  ENDDO
                ELSE
                  DO JJ=IDATA,JDATA
                    ISTYP = obs_bodyElem_i(obsSpaceData,OBS_VNM,JJ)
                    ZSLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,JJ)*MPC_MBAR_PER_PA_R8
                    IF ((ISTYP .EQ. BUFR_NEUU .OR. ISTYP .EQ. BUFR_NEUS)  &
                         .AND. ZLEV .EQ. ZSLEV) THEN
                      IOTHER = JJ
                    ENDIF
                  ENDDO
                ENDIF
                IF (IOTHER .NE. -1) THEN
                  ZOER = obs_bodyElem_r(obsSpaceData,OBS_OER,JDATA)
                  ZOVAL = obs_bodyElem_r(obsSpaceData,OBS_VAR,IOTHER)
                  ZOFCST = ZOVAL-obs_bodyElem_r(obsSpaceData,OBS_OMP,IOTHER)
                  ZSPDO = SQRT(ZOVAL*ZOVAL + ZVAL*ZVAL)
                  ZSPDF = SQRT(ZOFCST*ZOFCST + ZFCST*ZFCST)
                  ZDIFF = ZSPDO - ZSPDF
                  ILEV = NINT(ZLEV)
                  !
                  !___tighten rejection criterion for satob winds
                  !
                  IF (IDBURP .EQ. 88 .OR. IDBURP .EQ. 188) THEN
                    IF (ZDIFF .LT. 0.0d0 .AND. ABS(ZLAT) .GT. 20.d0 .AND.  &
                         ILEV .LT. 550) THEN
                       ZAASYM = 0.7d0*MAX(ZSPDF,1.0d0)
                       ZPRIOR = ZAASYM*ZAUV
                       IF (ZPRIOR .GT. 0.99d0) ZAASYM = 0.99d0/ZAUV
                    ELSE
                       !
                       !___zaasym used to specify new criterion for satwinds
                       !   than are not included in the asymmetric test
                       !
                       ZAASYM = 10.d0
                    ENDIF
                  ENDIF

                ENDIF
                !
                !__INITIAL VALUE OF GAMMA FOR QCVAR (WINDS)
                !
                call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(1.d0 -  &
                   (1.d0-(ZAUV*ZAASYM))*(1.d0-(ZAUV*ZAASYM)))*(2.d0*MPC_PI_r8)/  &
                  ((1.d0-(ZAUV*ZAASYM))*(1.d0-(ZAUV*ZAASYM))*  &
                                    (2.d0*ZDUV)*(2.d0*ZDUV)))
                IF ((IDBURP .GE. 32  .AND. IDBURP .LE. 38) .OR.  &
                   (IDBURP .GE. 135 .AND. IDBURP .LE. 142) .OR.  &
                   (IDBURP .EQ. 132) )  &
                call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,  &
                       (1.d0 - (1.d0-ZAUVRA)*(1.d0-ZAUVRA))  &
                       *(2.d0*MPC_PI_R8)/((1.d0-ZAUVRA)*(1.d0-ZAUVRA)*  &
                                   (2.d0*ZDUV)*(2.d0*ZDUV)))
             ELSEIF (ITYP .EQ. BUFR_NEGZ) THEN
                !
                ! INITIAL VALUE OF GAMMA FOR QCVAR (HEIGHTS)
                !
                call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZAGZ*  &
                    SQRT(2.d0*MPC_PI_R8))/((1.d0-ZAGZ)*(2.d0*ZDGZ)))
             ELSEIF (ITYP .EQ. BUFR_NETT) THEN
                !
                ! INITIAL VALUE OF GAMMA FOR QCVAR (TEMPERATURES)
                !
                call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZATT*  &
                    SQRT(2.d0*MPC_PI_R8))/((1.d0-ZATT)*(2.d0*ZDTT)))
                IF ((IDBURP .GE. 32  .AND. IDBURP .LE. 38) .OR.  &
                    (IDBURP .GE. 135 .AND. IDBURP .LE. 142) .OR.  &
                    (IDBURP .EQ. 132) )  &
                  call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZATTRA*  &
                    SQRT(2.d0*MPC_PI_R8))/((1.d0-ZATTRA)*(2.d0*ZDTT)))
             ELSEIF (ITYP .EQ. BUFR_NETS) THEN
                !
                ! INITIAL VALUE OF GAMMA FOR QCVAR (SCREEN-LEVEL TEMPERATURES)
                !
                call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZATT*  &
                    SQRT(2.d0*MPC_PI_R8))/((1.d0-ZATT)*(2.d0*ZDTT)))
                !
                ! ASYMMETRIC TEST FOR SHIP TEMPERATURES
                !
                IF ((IDBURP .EQ. 13) .AND. (ZVAL .GT. ZFCST)) THEN
                  call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZATTSYM*  &
                      SQRT(2.d0*MPC_PI_R8))/((1.d0-ZATTSYM)*(2.d0*ZDTT)))
                ENDIF
             ELSEIF (ITYP .EQ. BUFR_NEPN) THEN
                !
                ! INITIAL VALUE OF GAMMA FOR QCVAR (MSL PRESSURE)
                !
                call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZAPN*  &
                    SQRT(2.d0*MPC_PI_R8))/((1.d0-ZAPN)*(2.d0*ZDPN)))
             ELSEIF (ITYP .EQ. BUFR_NEPS) THEN
                !
                ! INITIAL VALUE OF GAMMA FOR QCVAR (STATION PRESSURE)
                !
                call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZAPS*  &
                    SQRT(2.d0*MPC_PI_R8))/((1.d0-ZAPS)*(2.d0*ZDPS)))
             ELSEIF ( ITYP .EQ. 12062 .OR. ITYP .EQ. 12063 .OR.  &
                      ITYP .EQ. 12163) THEN
                !
                ! INITIAL VALUE OF GAMMA FOR BRIGHTNESS TEMPERATURES
                ! TOVS AMSU-A + AIRS + IASI + GEORAD !!!!! 
                !
                call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZABT*  &
                    SQRT(2.d0*MPC_PI_R8))/((1.d0-ZABT)*(2.d0*ZDBT)))
                IF (IDBURP .EQ. 181 .OR. IDBURP .EQ. 182) THEN
                  !
                  ! INITIAL VALUE OF GAMMA FOR TOVS AMSU-B (181) AND MHS (182)
                  !
                  call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZABTB*  &
                        SQRT(2.d0*MPC_PI_R8))/((1.d0-ZABTB)*(2.d0*ZDBTB)))
                ENDIF
             ELSEIF (ITYP .EQ. BUFR_NEZD) THEN
                !
                ! INITIAL VALUE OF GAMMA FOR QCVAR (GPS ZENITH DELAY)
                !
                call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZAZD*  &
                    SQRT(2.d0*MPC_PI_R8))/((1.d0-ZAZD)*(2.d0*ZDZD)))
             ELSEIF (bufr_IsAtmosConstituent(ITYP)) then
                !
                ! INITIAL VALUE OF GAMMA FOR THE CH (constituents) family
                !
                  call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZACH* &
                       SQRT(2.d0*MPC_PI_R8))/((1.d0-ZACH)*(2.d0*ZDCH)))
             ELSE
                !
                ! INITIAL VALUE OF GAMMA FOR QCVAR (OTHERS)
                !
                call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZAHU*  &
                    SQRT(2.d0*MPC_PI_R8))/((1.d0-ZAHU)*(2.d0*ZDHU)))
             ENDIF
          ENDIF
       END DO

    END DO

  end subroutine vqc_setup


  subroutine vqc_tl(obsSpaceData)
    implicit none
    !
    !Purpose : 1) Modify Jo [OBS_JOBS] according to
    !             Andersson and Jarvinen 1999, Variational quality control,
    !             Q.J.R., 125, pp. 697-722.
    !          2) Save the values of (1-Wqc) in OBS_QCV
    !             for gradient factorization and postalt flag criterion.
    !
    !Author  : S. Pellerin, ARMA, January 2009
    !          Generalisation of QCVAR originally embeded in observation
    !          operators from P. Koclas, J. Halle and J. St-James
    !
    type(struct_obs) :: obsSpaceData
    integer :: index_body,istyp,jj,index_header,ityp,index_body_start,ierr,index_family
    real*8 :: zgami,zjon,zqcarg,zppost,zlev,zslev
    logical :: lluv
    logical :: includeFlag

    BODY: do index_body = 1, obs_numbody(obsSpaceData)
      includeFlag = (obs_bodyElem_i(obsSpaceData,OBS_ASS,index_body).eq.1).and.  &
                    (obs_getFamily(obsSpaceData,bodyIndex=index_body).ne.'RO')
      ! pas de qcvar pour  les radiances en mode matrice R non diagonale
      if (rmat_lnondiagr) includeFlag = includeFlag .and.  &
        (obs_getFamily(obsSpaceData,bodyIndex=index_body).ne.'TO') 

      if (includeFlag) then
        index_header = obs_bodyElem_i(obsSpaceData,OBS_HIND,index_body)
        index_body_start = obs_headElem_i(obsSpaceData,OBS_RLN,INDEX_HEADER)
        ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY)
        zgami = obs_bodyElem_r(obsSpaceData,OBS_POB,index_body)
        ITYP = obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY)
        LLUV = ((ITYP .EQ. BUFR_NEUU .OR. ITYP .EQ. BUFR_NEUS) .AND.  &
               col_varExist('UU')) .OR. ((ITYP .EQ. BUFR_NEVV .OR.  &
               ITYP .EQ. BUFR_NEVS).AND.col_varExist('VV'))
        IF (LLUV) THEN
          IF (ITYP .EQ. BUFR_NEUU .OR. ITYP .EQ. BUFR_NEUS)THEN
            !
            ! In order to calculate the contribution to Jo from a wind, the o-a
            ! must be available for both u and v components. Hence, loop over only
            ! data for which o-a has already been calculated
            !
            DO JJ=INDEX_BODY_START, INDEX_BODY
              ISTYP = obs_bodyElem_i(obsSpaceData,OBS_VNM,JJ)
              ZSLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,JJ)
              IF ((ISTYP .EQ. BUFR_NEVV .OR.  &
                   ISTYP .EQ. BUFR_NEVS) .AND.  &
                   ZSLEV .EQ. ZLEV) THEN
                ZJON=obs_bodyElem_r(obsSpaceData,OBS_JOBS,INDEX_BODY)+  &
                     obs_bodyElem_r(obsSpaceData,OBS_JOBS,JJ)
                ZQCARG = ZGAMI + EXP(-1.0D0*ZJON)
                ZPPOST = ZGAMI/ZQCARG
                !
                ! Store the value of o-a multiplied by one minus the posterior
                ! probability of gross error (needed for the adjoint calculations)
                !
                call obs_bodySet_r(obsSpaceData,OBS_QCV,INDEX_BODY, ZPPOST)
                call obs_bodySet_r(obsSpaceData,OBS_QCV,JJ, ZPPOST)

                call obs_bodySet_r(obsSpaceData,OBS_JOBS,INDEX_BODY,-LOG(ZQCARG/(ZGAMI+1.D0))/2.D0)
                call obs_bodySet_r(obsSpaceData,OBS_JOBS,JJ, -LOG(ZQCARG/(ZGAMI+1.D0))/2.D0)
              ENDIF
            ENDDO
          ELSE ! ITYP
            DO JJ=INDEX_BODY_START, INDEX_BODY
              ISTYP = obs_bodyElem_i(obsSpaceData,OBS_VNM,JJ)
              ZSLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,JJ)
              IF ((ISTYP .EQ. BUFR_NEUU .OR.  &
                   ISTYP .EQ. BUFR_NEUS) .AND.  &
                   ZSLEV .EQ. ZLEV) THEN
                ZJON=obs_bodyElem_r(obsSpaceData,OBS_JOBS,INDEX_BODY)+  &
                     obs_bodyElem_r(obsSpaceData,OBS_JOBS,JJ)
                ZQCARG = ZGAMI + EXP(-1.0D0*ZJON)
                ZPPOST = ZGAMI/ZQCARG
                call obs_bodySet_r(obsSpaceData,OBS_QCV,INDEX_BODY, ZPPOST)
                call obs_bodySet_r(obsSpaceData,OBS_QCV,JJ, ZPPOST)
                call obs_bodySet_r(obsSpaceData,OBS_JOBS,INDEX_BODY,-LOG(ZQCARG/(ZGAMI+1.D0))/2.D0)
                call obs_bodySet_r(obsSpaceData,OBS_JOBS,JJ, -LOG(ZQCARG/(ZGAMI+1.D0))/2.D0)
              ENDIF
            enddo
          endif !ITYP
        else ! LLUV
          zjon = obs_bodyElem_r(obsSpaceData,OBS_JOBS,index_body)
          zqcarg = zgami + exp(-1.0D0*zjon)
          zppost = zgami/zqcarg
          call obs_bodySet_r(obsSpaceData,OBS_QCV,index_body, zppost)
          call obs_bodySet_r(obsSpaceData,OBS_JOBS,INDEX_BODY, - log(zqcarg/(zgami+1.D0)))
        endif ! LLUV

      endif ! includeFlag

    enddo BODY

  end subroutine vqc_tl


  subroutine vqc_ad(obsSpaceData)
    !
    !Purpose : Factorizes Grad(Jo) according to Andersson and Jarvinen
    !          1999, Variational quality control, Q.J.R., 125,
    !          pp. 697-722.
    !          It uses the value of (1-Wqc) saved in OBS_QCV
    !          in vqc_tl
    !
    !Author  : S. Pellerin, ARMA, January 2009
    !          Generalisation of QCVAR originally embeded in adjoint of
    !          observation operators from P. Koclas, J. Halle and
    !          J. St-James
    !
    implicit none

    type(struct_obs) :: obsSpaceData
    integer :: index_body
    logical :: includeFlag

    do index_body=1,obs_numbody(obsSpaceData)
      includeFlag = (obs_bodyElem_i(obsSpaceData,OBS_ASS,index_body).eq.1) .and.  &
                    (obs_getFamily(obsSpaceData,bodyIndex=index_body).ne.'RO')
      ! pas de qcvar pour les radiances en mode matrice R non diagonale
      if (rmat_lnondiagr) includeFlag = includeFlag .and.  &
         (obs_getFamily(obsSpaceData,bodyIndex=index_body).ne.'TO')

      if (includeFlag) then
        call obs_bodySet_r(obsSpaceData,OBS_WORK,index_body,  &
               obs_bodyElem_r(obsSpaceData,OBS_WORK,index_body)  &
               *(1.d0 - obs_bodyElem_r(obsSpaceData,OBS_QCV,index_body)))
      endif
    enddo

  end subroutine vqc_ad


  subroutine vqc_listrej(lobsSpaceData)
    !
    !PURPOSE: LIST ALL OBSERVATIONS REJECTED BY THE VARIATIONAL QC
    !         SET QC FLAGS CONSISTENT WITH VARQC DECISIONS
    !         SET GLOBAL FLAG INDICATING REPORT CONTAINS REJECTED OBSERVATION
    !           AS REQUIRED
    !
    !
    !AUTHOR: B. BRASNETT (CMDA/MSC) MARCH 2000
    !
    !*REVISION: Y.J. Rochon ARQI Jan 2015
    !          - Additions of CH families with CODTYP=195 (remote sounding)
    !            and 196 (in-situ). 
    !          - Changed meters to hectometers (hm) for CH output
    !          - Output ZVAR, ZFCST and ZANA as percent difference relative to 
    !            ZFCST for CH constituents 
    !          - Extended output of obs_headElem_i(lobsSpaceData,OBS_ONM,INDEX_HEADER)
    !            from I5 to I7.
    !          Y. Rochon and M. Sitwell, June 2016
    !          - Introduction of codtypname to prevent recursive write error and 
    !            change of A16 to A21,1X, for consistency with max codtyp_get_name. 
    implicit none
    type(struct_obs) :: lobsSpaceData
    INTEGER, parameter :: numFamily = 10
    CHARACTER(len=2), parameter :: listFamily(numFamily) = (/'UA','AI','SF','SW','PR','RO','GP','SC','TO','CH'/)
    INTEGER, parameter :: NUMITEM=14
    INTEGER ICOUNT(NUMITEM,numFamily)

    INTEGER JFAM,JITEM,INDEX_HEADER
    INTEGER INDEX_BODY,INDEX_BODY2,ISTART,ITYP,ISTYP
    INTEGER ISPDO,IDIRO,ISPDF,IDIRF,ISPDA,IDIRA,ILEV
    INTEGER IDBURP,IOTHER
    REAL*8    ZVAR,ZFCST,ZANA,ZPOST,ZLAT,ZLON,ZUU,ZVV
    REAL*8    SPD,DEG,ZLEV,ZSLEV,ZCUT
    CHARACTER *4 CLITM(NUMITEM),CLDESC
    CHARACTER *2 CLUNITS
    CHARACTER *21 CODTYPNAME
    LOGICAL LLOK,LLELREJ
    !
    !  ------NOTE----------
    ! CURRENTLY SUPPORTED FAMILIES OF DATA 'UA' 'AI' 'SF' 'HU' 'TO' 'GO' 'GP'
    !
    DATA ZCUT /0.75D0/
    DATA CLITM(1), CLITM(2), CLITM(3), CLITM(4), CLITM(5), CLITM(6)  &
           /'WND',    'WND',    'HGT',    'TMP',    'DPD',   'STNP'/
    DATA CLITM(7), CLITM(8), CLITM(9), CLITM(10), CLITM(11), CLITM(12)  &
          /'MSLP',   'TSFC',   'SDPD',   'SWND',   'SWND',   'BTMP'/
    DATA CLITM(13),CLITM(14) /'ZTD', 'CHM'/


    DO JFAM=1,numFamily
      DO JITEM=1,NUMITEM
        ICOUNT(JITEM,JFAM) = 0
      ENDDO
    ENDDO
    WRITE(*,*) 'LIST OF DATA REJECTED BY VARIATIONAL QUALITY CONTROL'
    !
    !_____LOOP OVER ALL REPORTS, ONE FAMILY AT A TIME
    !
    FAMILY: DO JFAM = 1,numFamily

      WRITE(*,*) ' '
      IF (listFamily(JFAM) .NE. 'TO') THEN
        WRITE(*,*) 'IDENT     TYPE DESCRIPTION           LAT   LONG    LEVEL  ITEM   OBSVD     BKGND     ANAL  POST PROB REPORT'
      ELSE
        WRITE(*,*) 'IDENT     TYPE DESCRIPTION           LAT   LONG    CHNL   ITEM   OBSVD     BKGND     ANAL  POST PROB REPORT'
      ENDIF
      ! loop over all header indices of the family
      call obs_set_current_header_list(lobsSpaceData,listFamily(jfam))
      HEADER: do
        index_header = obs_getHeaderIndex(lobsSpaceData)
        if (index_header < 0) exit HEADER

        LLELREJ = .FALSE.
        IDBURP  = obs_headElem_i(lobsSpaceData,OBS_ITY,INDEX_HEADER)
        ZLAT = obs_headElem_r(lobsSpaceData,OBS_LAT,INDEX_HEADER)*MPC_DEGREES_PER_RADIAN_R8
        ZLON = obs_headElem_r(lobsSpaceData,OBS_LON,INDEX_HEADER)*MPC_DEGREES_PER_RADIAN_R8

        ! loop over all body indices for this index_header
        call obs_set_current_body_list(lobsSpaceData, index_header)
        BODY: do 
           index_body = obs_getBodyIndex(lobsSpaceData)
           if (index_body < 0) exit BODY
           ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,INDEX_BODY)
           IF (ITYP .EQ. BUFR_NETS .OR. ITYP .EQ. BUFR_NEPS .OR.  &
               ITYP .EQ. BUFR_NEPN .OR. ITYP .EQ. BUFR_NESS .OR.  &
               ITYP .EQ. BUFR_NEUS .OR. ITYP .EQ. BUFR_NEVS .OR.  &
               ITYP .EQ. BUFR_NEZD) THEN
              LLOK = (obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY) .EQ. 1)
           ELSE
              LLOK = (obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY) .EQ. 1 .AND.  &
                      obs_bodyElem_i(lobsSpaceData,OBS_XTR,INDEX_BODY) .EQ. 0) .OR.  &
                     (obs_bodyElem_i(lobsSpaceData,OBS_XTR,INDEX_BODY) .EQ. 2 .AND.  &
                      obs_bodyElem_i(lobsSpaceData,OBS_VNM,INDEX_BODY) .EQ. BUFR_NEGZ)
           ENDIF
           IF (LLOK) THEN
             zpost = obs_bodyElem_r(lobsSpaceData,OBS_QCV,INDEX_BODY)
             ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)*MPC_MBAR_PER_PA_R8
             IF (obs_bodyElem_i(lobsSpaceData,OBS_VCO,INDEX_BODY) .EQ. 2) THEN
               ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)*MPC_MBAR_PER_PA_R8
               CLUNITS = 'MB'
             ELSE IF (obs_bodyElem_i(lobsSpaceData,OBS_VCO,INDEX_BODY) .EQ. 1) THEN
               !
               ! VERTICAL COORDINATE IS HEIGHT
               !
               ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)
               CLUNITS = ' M'
               IF (listFamily(JFAM).EQ.'CH') THEN
                  ZLEV = ZLEV*0.01
                  CLUNITS = 'HM'
               END IF
             ELSE IF (obs_bodyElem_i(lobsSpaceData,OBS_VCO,INDEX_BODY) .EQ. -1) THEN
               !
               ! TOVS CHANNEL NUMBER
               !
               ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)
               CLUNITS = '  '
             ELSE 
               ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)
               CLUNITS = '  '
             ENDIF
             ZVAR = obs_bodyElem_r(lobsSpaceData,OBS_VAR,INDEX_BODY)
             ZFCST= ZVAR - obs_bodyElem_r(lobsSpaceData,OBS_OMP,INDEX_BODY)
             ZANA = ZVAR - obs_bodyElem_r(lobsSpaceData,OBS_OMA,INDEX_BODY)
             !
             !_____________TREAT WINDS AS SPECIAL CASE
             !              BUFR_NEUU       = 011003 (U COMPONENT)           (m/s)
             !              BUFR_NEVV       = 011004 (V COMPONENT)           (m/s)
             !              BUFR_NEUS       = 011215 (U COMPONENT AT 10 M)   (m/s)
             !              BUFR_NEVS       = 011216 (V COMPONENT AT 10 M)   (m/s)
             !
             IF (((ITYP.EQ.BUFR_NEUU .OR. ITYP .EQ. BUFR_NEUS) .AND.  &
                  col_varExist('UU')).OR.  &
                 ((ITYP.EQ.BUFR_NEVV .OR. ITYP .EQ. BUFR_NEVS) .AND.  &
                  col_varExist('VV'))) THEN
               IOTHER = -1
               IF (ITYP .EQ. BUFR_NEUU .OR. ITYP .EQ. BUFR_NEUS) THEN
                 ISTART=obs_headElem_i(lobsSpaceData,OBS_RLN,INDEX_HEADER)
                 DO INDEX_BODY2=ISTART,INDEX_BODY
                   ISTYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,INDEX_BODY2)
                   IF (obs_bodyElem_i(lobsSpaceData,OBS_VCO,INDEX_BODY2) .EQ. 2) THEN
                     ZSLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY2)*MPC_MBAR_PER_PA_R8
                   ENDIF
                   IF (obs_bodyElem_i(lobsSpaceData,OBS_VCO,INDEX_BODY2) .EQ. 1) THEN
                     ZSLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY2)
                   ENDIF
                   IF ((ISTYP .EQ. BUFR_NEVV .OR. ISTYP .EQ. BUFR_NEVS) .AND. &
                        ZLEV .EQ. ZSLEV) THEN
                     IOTHER = INDEX_BODY2
                   ENDIF
                 ENDDO
               ELSE
                 ISTART=obs_headElem_i(lobsSpaceData,OBS_RLN,INDEX_HEADER)
                 DO INDEX_BODY2=ISTART,INDEX_BODY
                   ISTYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,INDEX_BODY2)
                   IF (obs_bodyElem_i(lobsSpaceData,OBS_VCO,INDEX_BODY2) .EQ. 2) THEN
                     ZSLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY2)*MPC_MBAR_PER_PA_R8
                   ENDIF
                   IF (obs_bodyElem_i(lobsSpaceData,OBS_VCO,INDEX_BODY2) .EQ. 1) THEN
                     ZSLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY2)
                   ENDIF
                   IF ((ISTYP .EQ. BUFR_NEUU .OR. ISTYP .EQ. BUFR_NEUS) .AND.  &
                        ZLEV .EQ. ZSLEV) THEN
                     IOTHER = INDEX_BODY2
                   ENDIF
                 ENDDO
               ENDIF
               IF (IOTHER .NE. -1) THEN
                 IF (ITYP .EQ. BUFR_NEVV .OR. ITYP .EQ. BUFR_NEUS) THEN
                   ZUU = ZVAR
                   ZVV = obs_bodyElem_r(lobsSpaceData,OBS_VAR,IOTHER)
                 ELSE
                   ZVV = ZVAR
                   ZUU = obs_bodyElem_r(lobsSpaceData,OBS_VAR,IOTHER)
                 ENDIF
                 SPD = SQRT(ZUU*ZUU + ZVV*ZVV)
                 ISPDO = NINT(SPD)
                 IF (ZUU.EQ.0.D0 .AND. ZVV.EQ.0.D0)THEN
                   IDIRO = 999
                 ELSE
                   DEG = 270. - ATAN2(ZVV,ZUU)*MPC_DEGREES_PER_RADIAN_R8
                   IF (DEG .GT. 360.D0)DEG = DEG - 360.D0
                   IF (DEG .LT. 0.D0)  DEG = DEG + 360.D0
                   IDIRO = NINT(DEG)
                 ENDIF
                 IF (ITYP .EQ. BUFR_NEUU .OR. ITYP .EQ. BUFR_NEUS) THEN
                   ZUU = ZFCST
                   ZVV = obs_bodyElem_r(lobsSpaceData,OBS_VAR,IOTHER) - &
                         obs_bodyElem_r(lobsSpaceData,OBS_OMP,IOTHER)
                 ELSE
                   ZVV = ZFCST
                   ZUU = obs_bodyElem_r(lobsSpaceData,OBS_VAR,IOTHER) - &
                         obs_bodyElem_r(lobsSpaceData,OBS_OMP,IOTHER)
                 ENDIF
                 SPD=SQRT(ZUU*ZUU + ZVV*ZVV)
                 ISPDF = NINT(SPD)
                 IF (ZUU .EQ. 0.D0 .AND. ZVV .EQ. 0.D0) THEN
                   IDIRF = 999
                 ELSE
                   DEG = 270.D0 - ATAN2(ZVV,ZUU)*MPC_DEGREES_PER_RADIAN_R8
                   IF (DEG .GT. 360.D0)DEG = DEG - 360.D0
                   IF (DEG .LT. 0.D0)  DEG = DEG + 360.D0
                   IDIRF = NINT(DEG)
                 ENDIF
                 IF (ITYP .EQ. BUFR_NEUU .OR. ITYP .EQ. BUFR_NEUS) THEN
                   ZUU = ZANA
                   ZVV = obs_bodyElem_r(lobsSpaceData,OBS_VAR,IOTHER) -  &
                         obs_bodyElem_r(lobsSpaceData,OBS_OMA,IOTHER)
                 ELSE
                   ZVV = ZANA
                   ZUU = obs_bodyElem_r(lobsSpaceData,OBS_VAR,IOTHER) -  &
                         obs_bodyElem_r(lobsSpaceData,OBS_OMA,IOTHER)
                 ENDIF
                 SPD=SQRT(ZUU*ZUU + ZVV*ZVV)
                 ISPDA = NINT(SPD)
                 IF (ZUU.EQ.0.D0 .AND. ZVV.EQ.0.D0) THEN
                   IDIRA = 999
                 ELSE
                   DEG = 270.D0 - ATAN2(ZVV,ZUU)*MPC_DEGREES_PER_RADIAN_R8
                   IF (DEG .GT. 360.D0)DEG = DEG - 360.D0
                   IF (DEG .LT. 0.D0)  DEG = DEG + 360.D0
                   IDIRA = NINT(DEG)
                 ENDIF
                 ILEV = NINT(ZLEV)
                 IF (ZPOST .GT. ZCUT) THEN
                   LLELREJ = .TRUE.
                   call obs_bodySet_i(lobsSpaceData,OBS_FLG,INDEX_BODY,  &
                     IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,INDEX_BODY),9))
                   call obs_bodySet_i(lobsSpaceData,OBS_FLG,IOTHER,  &
                     IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,IOTHER),9))
                   call obs_bodySet_i(lobsSpaceData,OBS_FLG,INDEX_BODY,  &
                     IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,INDEX_BODY),17))
                   call obs_bodySet_i(lobsSpaceData,OBS_FLG,IOTHER,  &
                     IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,IOTHER),17))
                   IF (ITYP .EQ. BUFR_NEUU .OR.  &
                       ITYP .EQ. BUFR_NEVV) THEN
                     ICOUNT(1,JFAM) = ICOUNT(1,JFAM) + 1
                   ENDIF
                   IF (ITYP .EQ. BUFR_NEUS .OR. ITYP .EQ. BUFR_NEVS) THEN
                     ICOUNT(10,JFAM) = ICOUNT(10,JFAM) + 1
                   ENDIF
                   codtypname=codtyp_get_name(IDBURP)
                   WRITE(*,620) obs_elem_c(lobsSpaceData,'STID',INDEX_HEADER),  &
                     IDBURP,codtypname,ZLAT,ZLON,  &
                     ILEV,CLUNITS,IDIRO,ISPDO,IDIRF,ISPDF,IDIRA,  &
                     ISPDA,ZPOST,obs_headElem_i(lobsSpaceData,OBS_ONM,  &
                     INDEX_HEADER)
 620               FORMAT(A9,1X,i3,1x,A21,1X,F5.1,2X,F5.1,2X,I4,A2,2X,  &
                     'WND',3X,I3,'/',I3,3X,I3,'/',I3,3X,I3,  &
                     '/',I3,2X,F7.4,1X,I8)

                 ENDIF
               ENDIF
             ELSE
               ILEV = NINT(ZLEV)
               IF (ZPOST .GT. ZCUT) THEN
                 LLELREJ = .TRUE.
                 IF (listFamily(JFAM).EQ.'CH') THEN
                    ICOUNT(14,JFAM)=ICOUNT(14,JFAM)+1
                    CLDESC=CLITM(14)
                    if (obs_headElem_i(lobsSpaceData,OBS_CHM,INDEX_HEADER).ge.0) then

                       ! CONVERT TO PERCENTAGE DIFFERENCE FROM THE FORECAST

                       if (abs(ZFCST).lt.1.D-10*max(1.D-30,abs(ZVAR-ZFCST))) then
                          write(*,*) 'vqc_listrej: WARNING for CH obs. BKGND and or OBSVD value out of bounds for output.'
                          write(*,*) 'OBSVD and ANAL set to 0.0 below in addition to BKGND.'
                          ZVAR=0.0D0
                          ZANA=0.0D0
                       else
                          ZVAR=(ZVAR-ZFCST)/ZFCST*100.
                          ZANA=(ZANA-ZFCST)/ZFCST*100.
                       end if
                       ZFCST=0.0
                       CLDESC=vnl_varnameFromVarnum(ITYP,obs_headElem_i(lobsSpaceData,OBS_CHM,INDEX_HEADER))
                         
                    else if (ITYP .EQ. BUFR_NETT) THEN

                       ! CONVERT FROM KELVIN TO CELCIUS

                       CLDESC = CLITM(4)
                       ZVAR = ZVAR - MPC_K_C_DEGREE_OFFSET_R8
                       ZFCST = ZFCST - MPC_K_C_DEGREE_OFFSET_R8
                       ZANA  = ZANA - MPC_K_C_DEGREE_OFFSET_R8
                    end if
                 ELSE IF (ITYP .EQ. BUFR_NEGZ) THEN
                   CLDESC = CLITM(3)
                   ICOUNT(3,JFAM) = ICOUNT(3,JFAM) + 1
                   !
                   ! CONVERT FROM GEOPOTENTIAL TO GEOPOTENTIAL HEIGHT
                   !
                   ZVAR = ZVAR/RG
                   ZFCST = ZFCST/RG
                   ZANA  = ZANA/RG
                 ElSE IF (ITYP .EQ. BUFR_NETT) THEN
                   CLDESC = CLITM(4)
                   ICOUNT(4,JFAM) = ICOUNT(4,JFAM) + 1
                   !
                   ! CONVERT FROM KELVIN TO CELCIUS
                   !
                   ZVAR = ZVAR - MPC_K_C_DEGREE_OFFSET_R8
                   ZFCST = ZFCST - MPC_K_C_DEGREE_OFFSET_R8
                   ZANA  = ZANA - MPC_K_C_DEGREE_OFFSET_R8
                 ELSE IF (ITYP .EQ. BUFR_NEES) THEN
                   CLDESC = CLITM(5)
                   ICOUNT(5,JFAM) = ICOUNT(5,JFAM) + 1
                 ELSE IF (ITYP .EQ. BUFR_NEPS) THEN
                   CLDESC = CLITM(6)
                   ICOUNT(6,JFAM) = ICOUNT(6,JFAM) + 1
                   !
                   ! CONVERT FROM PASCALS TO MILLIBARS
                   !
                   ZVAR = ZVAR*MPC_MBAR_PER_PA_R8
                   ZANA = ZANA*MPC_MBAR_PER_PA_R8
                   ZFCST = ZFCST*MPC_MBAR_PER_PA_R8
                 ELSE IF (ITYP .EQ. BUFR_NEPN) THEN
                   CLDESC = CLITM(7)
                   ICOUNT(7,JFAM) = ICOUNT(7,JFAM) + 1
                   !
                   ! CONVERT FROM PASCALS TO MILLIBARS
                   !
                   ZVAR = ZVAR*MPC_MBAR_PER_PA_R8
                   ZANA = ZANA*MPC_MBAR_PER_PA_R8
                   ZFCST = ZFCST*MPC_MBAR_PER_PA_R8
                 ELSE IF (ITYP .EQ. BUFR_NETS) THEN
                   CLDESC = CLITM(8)
                   ICOUNT(8,JFAM) = ICOUNT(8,JFAM) + 1
                   !
                   ! CONVERT FROM KELVIN TO CELCIUS
                   !
                   ZVAR = ZVAR - MPC_K_C_DEGREE_OFFSET_R8
                   ZFCST = ZFCST - MPC_K_C_DEGREE_OFFSET_R8
                   ZANA  = ZANA - MPC_K_C_DEGREE_OFFSET_R8
                 ELSE IF (ITYP .EQ. BUFR_NESS) THEN
                   CLDESC = CLITM(9)
                   ICOUNT(9,JFAM) = ICOUNT(9,JFAM) + 1
                 ELSE IF ( ITYP.EQ.BUFR_NBT1 .OR. ITYP.EQ.BUFR_NBT2 .OR.  &
                      ITYP.EQ.BUFR_NBT3      ) THEN
                   CLDESC = CLITM(12)
                   ICOUNT(12,JFAM) = ICOUNT(12,JFAM) + 1
                 ELSE IF (ITYP .EQ. BUFR_NEZD) THEN
                   CLDESC = CLITM(13)
                   ICOUNT(13,JFAM) = ICOUNT(13,JFAM) + 1
                   !
                   ! CONVERT ZTD FROM METRES TO MILLIMETRES
                   !
                   ZVAR = ZVAR * 1000.D0
                   ZFCST = ZFCST * 1000.D0
                   ZANA  = ZANA * 1000.D0
                 ENDIF

                 call obs_bodySet_i(lobsSpaceData,OBS_FLG,INDEX_BODY,  &
                   IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,INDEX_BODY),9))
                 call obs_bodySet_i(lobsSpaceData,OBS_FLG,INDEX_BODY,  &
                   IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,INDEX_BODY),17))
                 codtypname=codtyp_get_name(IDBURP)
                 WRITE(*,630) obs_elem_c(lobsSpaceData,'STID',INDEX_HEADER),IDBURP,  &
                   codtypname,ZLAT,ZLON,ILEV,CLUNITS,CLDESC, &
                   ZVAR,ZFCST,ZANA,ZPOST,obs_headElem_i(lobsSpaceData,OBS_ONM,INDEX_HEADER)
 630               FORMAT(A9,1X,I3,1X,A21,1X,F5.1,2X,F5.1,2X,I4,A2,2X,  &
                    A4,2X,F7.1,3X,F7.1,3X,F7.1,2X,F7.4,1X,I8,1X)
               ENDIF
             ENDIF
           ENDIF
        ENDDO BODY
        !
        ! NOW SET THE GLOBAL FLAGS IN THE BURP RECORD HEADER
        !
        IF (LLELREJ) THEN
          call obs_headSet_i(lobsSpaceData,OBS_ST1,INDEX_HEADER, IBSET( obs_headElem_i(lobsSpaceData,OBS_ST1,INDEX_HEADER), 06))
        ENDIF
      ENDDO HEADER
    ENDDO FAMILY


    WRITE(*,640)
 640  FORMAT(//)
    WRITE(*,*) ' REJECTED DATA ACCORDING TO FAMILY OF REPORT.'
    WRITE(*,670)(listFamily(JFAM),JFAM=1,numFamily)
 670  FORMAT(5X,11(7X,A2))
    DO JITEM=1,NUMITEM
      IF ( .NOT. (JITEM .EQ. 2 .OR. JITEM .EQ. 11)) THEN
        WRITE(*,680)CLITM(JITEM),(ICOUNT(JITEM,JFAM),JFAM=1,numFamily)
 680      FORMAT(1X,A4,11(4X,I5))
      ENDIF
    ENDDO
    WRITE(*,640)

  end subroutine vqc_listrej


end module varqc_mod
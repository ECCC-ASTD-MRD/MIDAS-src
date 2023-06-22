
module varQC_mod
  ! MODULE varQC_mod (prefix='vqc' category='1. High-level functionality')
  !
  ! :Purpose: Procedures related to variational quality control including
  !           hard-coded values that determine how quickly the observation
  !           weight begins to be reduced
  !
  use MathPhysConstants_mod
  use earthConstants_mod
  use codtyp_mod
  use bufr_mod
  use obsSpaceData_mod
  use columnData_mod
  use rmatrix_mod ,only : rmat_lnondiagr
  use varNameList_mod
  use obsFamilyList_mod
  
  implicit none
  save
  private

  ! public procedures
  public :: vqc_setup, vqc_NlTl, vqc_ad, vqc_listrej


  contains

  subroutine vqc_setup(obsSpaceData)
    !
    ! :Purpose: To set certain parameters for the asymmetric check
    !           and for variational quality control
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData ! obsSpaceData object

    ! Locals:
    integer jdata, jjo, idata, idatend, idburp
    integer ityp, iass, iother, jj, istyp, ilev
    real(8) zagz, zahu, zatt, zauv, zabt, zabtb, zapn, zaps, zazd, zach, zatm, zaice, zapr
    real(8) zdgz, zdhu, zdtt, zduv, zdbt, zdbtb, zdpn, zdps, zdzd, zdch, zdtm, zdice, zdpr
    real(8) zattra, zauvra, zattsym, zdvis, zavis
    real(8) zslev, zlev, zval, zspdo, zspdf, zofcst, zoval, zdiff, zaasym, zoer
    real(8) zfcst, zlat, zlon, zprior
    logical llok

    !
    !_____prior probabilities for winds:ZAUV
    !     prior probabilities for scalar variables: zagz and zahu
    !     standard deviation multiple for background check for winds: zduv
    !     standard deviation multiple for background check for heights: zdgz
    !     standard deviation multiple for background check for humidity: zdhu
    !
    zagz  = 1.0d-12
!    zahu  = 1.0d-2
    zahu  = 5.0d-2
    zatt  = 5.0d-2
!    zauv  = 1.0d-2
    zauv  = 2.0d-2
    zabt  = 1.0d-1
    zabtb = 1.0d-1
    zapn  = 1.0d-4
    zaps  = 1.0d-4
    zazd  = 2.0d-2
    zach  = 1.0d-3
    zatm  = 1.0d-12
    zaice = 0.0d0
    zapr  = 0.0d0

    zattra = 0.005d0
    zattsym = 1.d-1
    zauvra = 1.d-5
    zavis= 1.0d-3

    zdgz  = 5.d0
    zdhu  = 5.d0
    zdtt  = 5.d0
    zduv  = 5.d0
    zdbt  = 3.d0
    zdbtb = 3.d0
    zdpn  = 5.d0
    zdps  = 5.d0
    zdzd  = 3.d0
    zdch  = 1.d1
    zdtm  = 5.d0
    zdice = 5.d2
    zdpr  = 5.d2
    zdvis= 5.d0

    do jjo = 1, obs_numheader(obsSpaceData)

       IDATA   = obs_headElem_i( obsSpaceData, OBS_RLN, jjo )
       IDATEND = obs_headElem_i( obsSpaceData, OBS_NLV, jjo ) + IDATA - 1
       IDBURP  = obs_headElem_i( obsSpaceData, OBS_ITY, jjo )
       ZLAT    = obs_headElem_r( obsSpaceData, OBS_LAT, jjo ) * MPC_DEGREES_PER_RADIAN_R8
       ZLON    = obs_headElem_r( obsSpaceData, OBS_LON, jjo ) * MPC_DEGREES_PER_RADIAN_R8

       do JDATA = IDATA, IDATEND

          ityp  = obs_bodyElem_i( obsSpaceData, OBS_VNM, JDATA )
          IASS  = obs_bodyElem_i( obsSpaceData, OBS_ASS, JDATA )
          ZLEV  = obs_bodyElem_r( obsSpaceData, OBS_PPP, JDATA ) * MPC_MBAR_PER_PA_R8
          ZOER  = obs_bodyElem_r( obsSpaceData, OBS_OER, JDATA )
          ZVAL  = obs_bodyElem_r( obsSpaceData, OBS_VAR, JDATA )

          ZFCST = ZVAL - obs_bodyElem_r( obsSpaceData, OBS_OMP,JDATA)

          if (ityp == BUFR_NETS .or. ityp == BUFR_NEPS .or.  &
              ityp == BUFR_NEPN .or. ityp == BUFR_NESS .or.  &
              ityp == BUFR_NEUS .or. ityp == BUFR_NEVS .or.  &
              ityp == BUFR_NEZD .or. ityp == bufr_vis  .or.  &
              ityp == bufr_logVis .or. ityp == bufr_gust .or. &
              ityp == bufr_radarPrecip .or. ityp == bufr_logRadarPrecip ) then
             LLOK = (obs_bodyElem_i(obsSpaceData,OBS_ASS,JDATA) == obs_assimilated)
          else
             LLOK = (IASS == obs_assimilated) .and. ((obs_bodyElem_i(obsSpaceData,OBS_XTR,JDATA) ==0) &
                .or. ((obs_bodyElem_i(obsSpaceData,OBS_XTR,JDATA) == 2) .and. &
                      (obs_bodyElem_i(obsSpaceData,OBS_VNM,JDATA) == BUFR_NEGZ)))
          end if

          if (LLOK) then
             if (ityp == BUFR_NEUU .or. ityp == BUFR_NEVV .or.  &
                 ityp == BUFR_NEUS .or. ityp == BUFR_NEVS) then
                ZAASYM = 1.0d0
                IOTHER = -1
                if (ityp == BUFR_NEUU .or. ityp == BUFR_NEUS) then
                  do JJ = IDATA, JDATA
                    ISTYP = obs_bodyElem_i( obsSpaceData, OBS_VNM, JJ )
                    ZSLEV = obs_bodyElem_r( obsSpaceData, OBS_PPP, JJ ) * MPC_MBAR_PER_PA_R8
                    if ((ISTYP == BUFR_NEVV .or. ISTYP == BUFR_NEVS)  &
                         .and. ZLEV == ZSLEV) then
                      IOTHER = JJ
                    end if
                  end do
                else
                  do JJ=IDATA,JDATA
                    ISTYP = obs_bodyElem_i(obsSpaceData,OBS_VNM,JJ)
                    ZSLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,JJ)*MPC_MBAR_PER_PA_R8
                    if ((ISTYP == BUFR_NEUU .or. ISTYP == BUFR_NEUS)  &
                         .and. ZLEV == ZSLEV) then
                      IOTHER = JJ
                    end if
                  end do
                end if
                if (IOTHER /= -1) then
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
                  if (IDBURP == 88 .or. IDBURP == 188) then
                    if (ZDIFF < 0.0d0 .and. ABS(ZLAT) > 20.d0 .and.  &
                         ILEV < 550) then
                       ZAASYM = 0.7d0*MAX(ZSPDF,1.0d0)
                       ZPRIOR = ZAASYM*ZAUV
                       if (ZPRIOR > 0.99d0) ZAASYM = 0.99d0/ZAUV
                    else
                       !
                       !___zaasym used to specify new criterion for satwinds
                       !   than are not included in the asymmetric test
                       !
                       ZAASYM = 10.d0
                    end if
                  end if

                end if
                !
                !__INITIAL VALUE OF GAMMA FOR QCVAR (WINDS)
                !
                call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(1.d0 -  &
                   (1.d0-(ZAUV*ZAASYM))*(1.d0-(ZAUV*ZAASYM)))*(2.d0*MPC_PI_r8)/  &
                  ((1.d0-(ZAUV*ZAASYM))*(1.d0-(ZAUV*ZAASYM))*  &
                                    (2.d0*ZDUV)*(2.d0*ZDUV)))
                if ((IDBURP >= 32  .and. IDBURP <= 38) .or.  &
                   (IDBURP >= 135 .and. IDBURP <= 142) .or.  &
                   (IDBURP == 132) )  &
                call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,  &
                       (1.d0 - (1.d0-ZAUVRA)*(1.d0-ZAUVRA))  &
                       *(2.d0*MPC_PI_R8)/((1.d0-ZAUVRA)*(1.d0-ZAUVRA)*  &
                                   (2.d0*ZDUV)*(2.d0*ZDUV)))
             else if (ityp == BUFR_NEGZ) then
                !
                ! INITIAL VALUE OF GAMMA FOR QCVAR (HEIGHTS)
                !
                call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZAGZ*  &
                    SQRT(2.d0*MPC_PI_R8))/((1.d0-ZAGZ)*(2.d0*ZDGZ)))
             else if (ityp == BUFR_NETT) then
                !
                ! INITIAL VALUE OF GAMMA FOR QCVAR (TEMPERATURES)
                !
                call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZATT*  &
                    SQRT(2.d0*MPC_PI_R8))/((1.d0-ZATT)*(2.d0*ZDTT)))
                if ((IDBURP >= 32  .and. IDBURP <= 38) .or.  &
                    (IDBURP >= 135 .and. IDBURP <= 142) .or.  &
                    (IDBURP == 132) )  &
                  call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZATTRA*  &
                    SQRT(2.d0*MPC_PI_R8))/((1.d0-ZATTRA)*(2.d0*ZDTT)))
             else if (ityp == BUFR_NETS) then
                !
                ! INITIAL VALUE OF GAMMA FOR QCVAR (SCREEN-LEVEL TEMPERATURES)
                !
                call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZATT*  &
                    SQRT(2.d0*MPC_PI_R8))/((1.d0-ZATT)*(2.d0*ZDTT)))
                !
                ! ASYMMETRIC TEST FOR SHIP TEMPERATURES
                !
                if ((IDBURP == 13) .and. (ZVAL > ZFCST)) then
                  call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZATTSYM*  &
                      SQRT(2.d0*MPC_PI_R8))/((1.d0-ZATTSYM)*(2.d0*ZDTT)))
                end if
             else if (ityp == BUFR_NEPN) then
                !
                ! INITIAL VALUE OF GAMMA FOR QCVAR (MSL PRESSURE)
                !
                call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZAPN*  &
                    SQRT(2.d0*MPC_PI_R8))/((1.d0-ZAPN)*(2.d0*ZDPN)))
             else if (ityp == BUFR_NEPS) then
                !
                ! INITIAL VALUE OF GAMMA FOR QCVAR (STATION PRESSURE)
                !
                call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZAPS*  &
                    SQRT(2.d0*MPC_PI_R8))/((1.d0-ZAPS)*(2.d0*ZDPS)))
             else if ( ityp == 12062 .or. ityp == 12063 .or.  &
                      ityp == 12163) then
                !
                ! INITIAL VALUE OF GAMMA FOR BRIGHTNESS TEMPERATURES
                ! TOVS AMSU-A + AIRS + IASI + GEORAD !!!!! 
                !
                call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZABT*  &
                    SQRT(2.d0*MPC_PI_R8))/((1.d0-ZABT)*(2.d0*ZDBT)))
                if (IDBURP == 181 .or. IDBURP == 182 .or. &
                    IDBURP == 200) then
                  !
                  ! INITIAL VALUE OF GAMMA FOR TOVS AMSU-B (181) AND MHS (182)
                  ! AND MWHS2(200)
                  !
                  call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZABTB*  &
                        SQRT(2.d0*MPC_PI_R8))/((1.d0-ZABTB)*(2.d0*ZDBTB)))
                end if
             else if (ityp == BUFR_NEZD) then
                !
                ! INITIAL VALUE OF GAMMA FOR QCVAR (GPS ZENITH DELAY)
                !
                call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZAZD*  &
                    SQRT(2.d0*MPC_PI_R8))/((1.d0-ZAZD)*(2.d0*ZDZD)))
             else if (bufr_IsAtmosConstituent(ityp)) then
                !
                ! INITIAL VALUE OF GAMMA FOR THE CH (constituents) family
                !
                  call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZACH* &
                       SQRT(2.d0*MPC_PI_R8))/((1.d0-ZACH)*(2.d0*ZDCH)))
             else if (ityp == bufr_sst ) then
                !
                ! INITIAL VALUE OF GAMMA FOR SST
                !
                call obs_bodySet_r( obsSpaceData, OBS_POB, JDATA, ( zatm *  &
                    sqrt( 2.d0 * MPC_PI_R8 )) / (( 1.d0 - zatm ) * ( 2.d0 * zdtm )))
             else if (ityp == BUFR_ICEC ) then
                !
                ! INITIAL VALUE OF GAMMA FOR SEA ICE
                !
                call obs_bodySet_r( obsSpaceData, OBS_POB, JDATA, ( zaice *  &
                    sqrt( 2.d0 * MPC_PI_R8 )) / (( 1.d0 - zaice ) * ( 2.d0 * zdice )))
             else if (ityp == BUFR_radarPrecip .or. ityp == BUFR_logRadarPrecip ) then
                !
                ! INITIAL VALUE OF GAMMA FOR RADAR PRECIPITATION
                !
                call obs_bodySet_r( obsSpaceData, OBS_POB, JDATA, ( zapr *  &
                    sqrt( 2.d0 * MPC_PI_R8 )) / (( 1.d0 - zapr ) * ( 2.d0 * zdpr )))
             else if (ityp == bufr_gust ) then
                !
                ! INITIAL VALUE OF GAMMA FOR WIND GUST
                !
                call obs_bodySet_r( obsSpaceData, OBS_POB, JDATA, ( zauv *  &
                    sqrt( 2.d0 * MPC_PI_R8 )) / (( 1.d0 - zauv ) * ( 2.d0 * zduv )))
              else if (ityp == bufr_vis .or. ityp == bufr_logVis) then
                !
                ! INITIAL VALUE OF GAMMA FOR VISIBILITY
                !
                call obs_bodySet_r( obsSpaceData, OBS_POB, JDATA, ( zavis *  &
                    sqrt( 2.d0 * MPC_PI_R8 )) / (( 1.d0 - zavis ) * ( 2.d0 * zdvis )))
             else
                !
                ! INITIAL VALUE OF GAMMA FOR QCVAR (OTHERS)
                !
                call obs_bodySet_r(obsSpaceData,OBS_POB,JDATA,(ZAHU*  &
                    SQRT(2.d0*MPC_PI_R8))/((1.d0-ZAHU)*(2.d0*ZDHU)))
             end if
          end if
       end do

    end do

  end subroutine vqc_setup


  subroutine vqc_NlTl(obsSpaceData)
    !
    ! :Purpose: 1) Modify Jo [OBS_JOBS] according to
    !              Andersson and Jarvinen 1999, Variational quality control,
    !              Q.J.R., 125, pp. 697-722.
    !           2) Save the values of (1-Wqc) in OBS_QCV
    !              for gradient factorization and postalt flag criterion.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData ! obsSpaceData object

    ! Locals:
    integer :: index_body,istyp,jj,index_header,ityp,index_body_start
    real*8 :: zgami,zjon,zqcarg,zppost,zlev,zslev
    logical :: lluv
    logical :: includeFlag

    BODY: do index_body = 1, obs_numbody(obsSpaceData)
      includeFlag = (obs_bodyElem_i(obsSpaceData,OBS_ASS,index_body) == obs_assimilated).and.  &
                    (obs_getFamily(obsSpaceData,bodyIndex_opt=index_body).ne.'RO')
      ! pas de qcvar pour  les radiances en mode matrice R non diagonale
      if (rmat_lnondiagr) includeFlag = includeFlag .and.  &
        (obs_getFamily(obsSpaceData,bodyIndex_opt=index_body) /= 'TO') 

      if (includeFlag) then
        index_header = obs_bodyElem_i(obsSpaceData,OBS_HIND,index_body)
        index_body_start = obs_headElem_i(obsSpaceData,OBS_RLN,INDEX_HEADER)
        ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY)
        zgami = obs_bodyElem_r(obsSpaceData,OBS_POB,index_body)
        ityp = obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY)
        LLUV = ((ityp == BUFR_NEUU .or. ityp == BUFR_NEUS) .and.  &
               col_varExist(varName='UU')) .or. ((ityp == BUFR_NEVV .or.  &
               ityp == BUFR_NEVS).and.col_varExist(varName='VV'))
        if (LLUV) then
          if (ityp == BUFR_NEUU .or. ityp == BUFR_NEUS)then
            !
            ! In order to calculate the contribution to Jo from a wind, the o-a
            ! must be available for both u and v components. Hence, loop over only
            ! data for which o-a has already been calculated
            !
            do JJ=INDEX_BODY_START, INDEX_BODY
              ISTYP = obs_bodyElem_i(obsSpaceData,OBS_VNM,JJ)
              ZSLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,JJ)
              if ((ISTYP == BUFR_NEVV .or.  &
                   ISTYP == BUFR_NEVS) .and.  &
                   ZSLEV == ZLEV) then
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
              end if
            end do
          else ! ityp
            do JJ=INDEX_BODY_START, INDEX_BODY
              ISTYP = obs_bodyElem_i(obsSpaceData,OBS_VNM,JJ)
              ZSLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,JJ)
              if ((ISTYP == BUFR_NEUU .or.  &
                   ISTYP == BUFR_NEUS) .and.  &
                   ZSLEV == ZLEV) then
                ZJON=obs_bodyElem_r(obsSpaceData,OBS_JOBS,INDEX_BODY)+  &
                     obs_bodyElem_r(obsSpaceData,OBS_JOBS,JJ)
                ZQCARG = ZGAMI + EXP(-1.0D0*ZJON)
                ZPPOST = ZGAMI/ZQCARG
                call obs_bodySet_r(obsSpaceData,OBS_QCV,INDEX_BODY, ZPPOST)
                call obs_bodySet_r(obsSpaceData,OBS_QCV,JJ, ZPPOST)
                call obs_bodySet_r(obsSpaceData,OBS_JOBS,INDEX_BODY,-LOG(ZQCARG/(ZGAMI+1.D0))/2.D0)
                call obs_bodySet_r(obsSpaceData,OBS_JOBS,JJ, -LOG(ZQCARG/(ZGAMI+1.D0))/2.D0)
              end if
            enddo
          endif !ityp
        else ! LLUV
          zjon = obs_bodyElem_r(obsSpaceData,OBS_JOBS,index_body)
          zqcarg = zgami + exp(-1.0D0*zjon)
          zppost = zgami/zqcarg
          call obs_bodySet_r(obsSpaceData,OBS_QCV,index_body, zppost)
          call obs_bodySet_r(obsSpaceData,OBS_JOBS,INDEX_BODY, - log(zqcarg/(zgami+1.D0)))
        endif ! LLUV

      endif ! includeFlag

    enddo BODY

  end subroutine vqc_NlTl


  subroutine vqc_ad(obsSpaceData)
    !
    ! :Purpose: Factorizes Grad(Jo) according to Andersson and Jarvinen
    !           1999, Variational quality control, Q.J.R., 125, pp. 697-722.
    !           It uses the value of (1-Wqc) saved in OBS_QCV in vqc_NlTl
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData ! obsSpaceData object

    ! Locals:
    integer :: index_body
    logical :: includeFlag

    do index_body=1,obs_numbody(obsSpaceData)
      includeFlag = (obs_bodyElem_i(obsSpaceData,OBS_ASS,index_body) == obs_assimilated) .and.  &
                    (obs_getFamily(obsSpaceData,bodyIndex_opt=index_body) /= 'RO')
      ! pas de qcvar pour les radiances en mode matrice R non diagonale
      if (rmat_lnondiagr) includeFlag = includeFlag .and.  &
         (obs_getFamily(obsSpaceData,bodyIndex_opt=index_body) /= 'TO')

      if (includeFlag) then
        call obs_bodySet_r(obsSpaceData,OBS_WORK,index_body,  &
               obs_bodyElem_r(obsSpaceData,OBS_WORK,index_body)  &
               *(1.d0 - obs_bodyElem_r(obsSpaceData,OBS_QCV,index_body)))
      endif
    enddo

  end subroutine vqc_ad


  subroutine vqc_listrej(lobsSpaceData)
    !
    ! :Purpose: List all observations rejected by the variational QC
    !           Set QC flags consistent with VARQC decisions
    !           Set global flag indicating report contains rejected observations
    !           as required
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: lobsSpaceData

    ! Locals:
    integer, parameter :: numitem = 16
    integer :: icount( numitem, ofl_numFamily )
    integer :: jfam, jitem, headerIndex
    integer :: bodyIndex, bodyIndex2, ISTART,ityp,ISTYP
    integer :: ISPDO,IDIRO,ISPDF,IDIRF,ISPDA,IDIRA,ILEV
    integer :: IDBURP,IOTHER,obsONM
    real(8) :: ZVAR,ZFCST,ZANA,ZPOST,ZLAT,ZLON,ZUU,ZVV
    real(8) :: SPD,DEG,ZLEV,ZSLEV,ZCUT
    character(len=4)  :: CLITM(NUMITEM),CLDESC
    character(len=2)  :: CLUNITS
    character(len=21) :: CODTYPNAME
    character(len=12) :: stnId
    logical :: LLOK,LLELREJ

    !  ------NOTE----------
    ! CURRENTLY SUPPORTED FAMILIES OF DATA 'UA' 'AI' 'SF' 'HU' 'TO' 'GO' 'GP'

    DATA ZCUT /0.75D0/
    DATA CLITM(1), CLITM(2), CLITM(3), CLITM(4), CLITM(5), CLITM(6)  &
           /'WND',    'WND',    'HGT',    'TMP',    'DPD',   'STNP'/
    DATA CLITM(7), CLITM(8), CLITM(9), CLITM(10), CLITM(11), CLITM(12)  &
          /'MSLP',   'TSFC',   'SDPD',   'SWND',   'SWND',   'BTMP'/
    DATA CLITM(13), CLITM(14), CLITM(15) , CLITM(16)  &
           /'ZTD',   'CHM',     'SST',    'ICEC'/


    do jfam = 1, ofl_numFamily
      do jitem = 1, numitem
        icount( jitem, jfam ) = 0
      end do
    end do
    write(*,*) 'LIST OF DATA REJECTED BY VARIATIONAL QUALITY CONTROL'
    !
    !_____LOOP OVER ALL REPORTS, ONE FAMILY AT A TIME
    !
    FAMILY: do jfam = 1, ofl_numFamily

      write(*,*) ' '
      if (ofl_familyList(jfam) /= 'TO' ) then
        write(*,*) 'IDENT     TYPE DESCRIPTION           LAT   LONG    LEVEL  ITEM   OBSVD     BKGND     ANAL  POST PROB REPORT'
      else
        write(*,*) 'IDENT     TYPE DESCRIPTION           LAT   LONG    CHNL   ITEM   OBSVD     BKGND     ANAL  POST PROB REPORT'
      end if

      ! loop over all header indices of the family
      call obs_set_current_header_list( lobsSpaceData, ofl_familyList( jfam ) )

      HEADER: do
        headerIndex = obs_getHeaderIndex( lobsSpaceData )
        if (headerIndex < 0) exit HEADER

        llelrej = .false.
        idburp  = obs_headElem_i( lobsSpaceData, OBS_ITY, headerIndex )
        zlat    = obs_headElem_r( lobsSpaceData, OBS_LAT, headerIndex ) * MPC_DEGREES_PER_RADIAN_R8
        zlon    = obs_headElem_r( lobsSpaceData, OBS_LON, headerIndex ) * MPC_DEGREES_PER_RADIAN_R8

        ! loop over all body indices for this headerIndex
        call obs_set_current_body_list( lobsSpaceData, headerIndex )

        BODY: do 
           bodyIndex = obs_getBodyIndex( lobsSpaceData )
           if ( bodyIndex < 0) exit BODY
           ityp = obs_bodyElem_i( lobsSpaceData, OBS_VNM, bodyIndex )
           if (ityp == BUFR_NETS .or. ityp == BUFR_NEPS .or.  &
               ityp == BUFR_NEPN .or. ityp == BUFR_NESS .or.  &
               ityp == BUFR_NEUS .or. ityp == BUFR_NEVS .or.  &
               ityp == BUFR_NEZD .or. ityp == bufr_sst  .or. ityp == BUFR_ICEC .or.  &
               ityp == bufr_vis  .or. ityp == bufr_logVis  .or. ityp == bufr_gust .or.  &
               ityp == bufr_radarPrecip .or. ityp == bufr_logRadarPrecip ) then
              llok = (obs_bodyElem_i( lobsSpaceData, OBS_ASS, bodyIndex ) == obs_assimilated )
           else
              llok = (obs_bodyElem_i( lobsSpaceData, OBS_ASS, bodyIndex ) == obs_assimilated .and.  &
                      obs_bodyElem_i( lobsSpaceData, OBS_XTR, bodyIndex ) == 0) .or.  &
                     (obs_bodyElem_i( lobsSpaceData, OBS_XTR, bodyIndex ) == 2 .and.  &
                      obs_bodyElem_i( lobsSpaceData, OBS_VNM, bodyIndex ) == BUFR_NEGZ)
           end if
           if ( llok ) then
             zpost = obs_bodyElem_r( lobsSpaceData, OBS_QCV, bodyIndex )
             zlev  = obs_bodyElem_r( lobsSpaceData, OBS_PPP, bodyIndex ) * MPC_MBAR_PER_PA_R8
             if ( obs_bodyElem_i( lobsSpaceData, OBS_VCO, bodyIndex ) == 2) then
               zlev = obs_bodyElem_r( lobsSpaceData, OBS_PPP, bodyIndex ) * MPC_MBAR_PER_PA_R8
               CLUNITS = 'MB'
             else if (obs_bodyElem_i(lobsSpaceData,OBS_VCO,bodyIndex) == 1) then
               !
               ! VERTICAL COORDINATE IS HEIGHT
               !
               ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,bodyIndex)
               CLUNITS = ' M'
               if (ofl_familyList(JFAM)=='CH') then
                  ZLEV = ZLEV*0.01
                  CLUNITS = 'HM'
               end if
             else if (obs_bodyElem_i(lobsSpaceData,OBS_VCO,bodyIndex) == -1) then
               !
               ! TOVS CHANNEL NUMBER
               !
               ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,bodyIndex)
               CLUNITS = '  '
             else 
               ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,bodyIndex)
               CLUNITS = '  '
             end if

             ZVAR = obs_bodyElem_r(lobsSpaceData,OBS_VAR,bodyIndex)

             ZFCST= ZVAR - obs_bodyElem_r(lobsSpaceData,OBS_OMP,bodyIndex)
             ZANA = ZVAR - obs_bodyElem_r(lobsSpaceData,OBS_OMA,bodyIndex)
             !
             !_____________TREAT WINDS AS SPECIAL CASE
             !              BUFR_NEUU       = 011003 (U COMPONENT)           (m/s)
             !              BUFR_NEVV       = 011004 (V COMPONENT)           (m/s)
             !              BUFR_NEUS       = 011215 (U COMPONENT AT 10 M)   (m/s)
             !              BUFR_NEVS       = 011216 (V COMPONENT AT 10 M)   (m/s)
             !
             if (((ityp==BUFR_NEUU .or. ityp == BUFR_NEUS) .and.  &
                  col_varExist(varName='UU')).or.  &
                 ((ityp==BUFR_NEVV .or. ityp == BUFR_NEVS) .and.  &
                  col_varExist(varName='VV'))) then
               IOTHER = -1
               if (ityp == BUFR_NEUU .or. ityp == BUFR_NEUS) then
                 ISTART=obs_headElem_i(lobsSpaceData,OBS_RLN,headerIndex)
                 do bodyIndex2=ISTART,bodyIndex
                   ISTYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,bodyIndex2)
                   if (obs_bodyElem_i(lobsSpaceData,OBS_VCO,bodyIndex2) == 2) then
                     ZSLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,bodyIndex2)*MPC_MBAR_PER_PA_R8
                   end if
                   if (obs_bodyElem_i(lobsSpaceData,OBS_VCO,bodyIndex2) == 1) then
                     ZSLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,bodyIndex2)
                   end if
                   if ((ISTYP == BUFR_NEVV .or. ISTYP == BUFR_NEVS) .and. &
                        ZLEV == ZSLEV) then
                     IOTHER = bodyIndex2
                   end if
                 end do
               else
                 ISTART=obs_headElem_i(lobsSpaceData,OBS_RLN,headerIndex)
                 do bodyIndex2=ISTART,bodyIndex
                   ISTYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,bodyIndex2)
                   if (obs_bodyElem_i(lobsSpaceData,OBS_VCO,bodyIndex2) == 2) then
                     ZSLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,bodyIndex2)*MPC_MBAR_PER_PA_R8
                   end if
                   if (obs_bodyElem_i(lobsSpaceData,OBS_VCO,bodyIndex2) == 1) then
                     ZSLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,bodyIndex2)
                   end if
                   if ((ISTYP == BUFR_NEUU .or. ISTYP == BUFR_NEUS) .and.  &
                        ZLEV == ZSLEV) then
                     IOTHER = bodyIndex2
                   end if
                 end do
               end if
               if (IOTHER /= -1) then
                 if (ityp == BUFR_NEVV .or. ityp == BUFR_NEUS) then
                   ZUU = ZVAR
                   ZVV = obs_bodyElem_r(lobsSpaceData,OBS_VAR,IOTHER)
                 else
                   ZVV = ZVAR
                   ZUU = obs_bodyElem_r(lobsSpaceData,OBS_VAR,IOTHER)
                 end if
                 SPD = SQRT(ZUU*ZUU + ZVV*ZVV)
                 ISPDO = NINT(SPD)
                 if (ZUU==0.D0 .and. ZVV==0.D0)then
                   IDIRO = 999
                 else
                   DEG = 270. - ATAN2(ZVV,ZUU)*MPC_DEGREES_PER_RADIAN_R8
                   if (DEG > 360.D0)DEG = DEG - 360.D0
                   if (DEG < 0.D0)  DEG = DEG + 360.D0
                   IDIRO = NINT(DEG)
                 end if
                 if (ityp == BUFR_NEUU .or. ityp == BUFR_NEUS) then
                   ZUU = ZFCST
                   ZVV = obs_bodyElem_r(lobsSpaceData,OBS_VAR,IOTHER) - &
                         obs_bodyElem_r(lobsSpaceData,OBS_OMP,IOTHER)
                 else
                   ZVV = ZFCST
                   ZUU = obs_bodyElem_r(lobsSpaceData,OBS_VAR,IOTHER) - &
                         obs_bodyElem_r(lobsSpaceData,OBS_OMP,IOTHER)
                 end if
                 SPD=SQRT(ZUU*ZUU + ZVV*ZVV)
                 ISPDF = NINT(SPD)
                 if (ZUU == 0.D0 .and. ZVV == 0.D0) then
                   IDIRF = 999
                 else
                   DEG = 270.D0 - ATAN2(ZVV,ZUU)*MPC_DEGREES_PER_RADIAN_R8
                   if (DEG > 360.D0)DEG = DEG - 360.D0
                   if (DEG < 0.D0)  DEG = DEG + 360.D0
                   IDIRF = NINT(DEG)
                 end if
                 if (ityp == BUFR_NEUU .or. ityp == BUFR_NEUS) then
                   ZUU = ZANA
                   ZVV = obs_bodyElem_r(lobsSpaceData,OBS_VAR,IOTHER) -  &
                         obs_bodyElem_r(lobsSpaceData,OBS_OMA,IOTHER)
                 else
                   ZVV = ZANA
                   ZUU = obs_bodyElem_r(lobsSpaceData,OBS_VAR,IOTHER) -  &
                         obs_bodyElem_r(lobsSpaceData,OBS_OMA,IOTHER)
                 end if
                 SPD=SQRT(ZUU*ZUU + ZVV*ZVV)
                 ISPDA = NINT(SPD)
                 if (ZUU==0.D0 .and. ZVV==0.D0) then
                   IDIRA = 999
                 else
                   DEG = 270.D0 - ATAN2(ZVV,ZUU)*MPC_DEGREES_PER_RADIAN_R8
                   if (DEG > 360.D0)DEG = DEG - 360.D0
                   if (DEG < 0.D0)  DEG = DEG + 360.D0
                   IDIRA = NINT(DEG)
                 end if
                 ILEV = NINT(ZLEV)
                 if (ZPOST > ZCUT) then
                   LLELREJ = .TRUE.
                   call obs_bodySet_i(lobsSpaceData,OBS_FLG,bodyIndex,  &
                     IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,bodyIndex),9))
                   call obs_bodySet_i(lobsSpaceData,OBS_FLG,IOTHER,  &
                     IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,IOTHER),9))
                   call obs_bodySet_i(lobsSpaceData,OBS_FLG,bodyIndex,  &
                     IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,bodyIndex),17))
                   call obs_bodySet_i(lobsSpaceData,OBS_FLG,IOTHER,  &
                     IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,IOTHER),17))
                   if (ityp == BUFR_NEUU .or.  &
                       ityp == BUFR_NEVV) then
                     ICOUNT(1,JFAM) = ICOUNT(1,JFAM) + 1
                   end if
                   if (ityp == BUFR_NEUS .or. ityp == BUFR_NEVS) then
                     ICOUNT(10,JFAM) = ICOUNT(10,JFAM) + 1
                   end if
                   codtypname=codtyp_get_name(IDBURP)
                   stnId = obs_elem_c(lobsSpaceData,'STID',headerIndex)
                   obsONM = obs_headElem_i(lobsSpaceData,OBS_ONM,headerIndex)
                   write(*,620) stnId,IDBURP,codtypname,ZLAT,ZLON,  &
                     ILEV,CLUNITS,IDIRO,ISPDO,IDIRF,ISPDF,IDIRA,  &
                     ISPDA,ZPOST,obsONM
 620               FORMAT(A9,1X,i3,1x,A21,1X,F5.1,2X,F5.1,2X,I4,A2,2X,  &
                     'WND',3X,I3,'/',I3,3X,I3,'/',I3,3X,I3,  &
                     '/',I3,2X,F7.4,1X,I8)

                 end if
               end if
             else
               ILEV = NINT(ZLEV)
               if (ZPOST > ZCUT) then
                 LLELREJ = .TRUE.
                 if (ofl_familyList(JFAM)=='CH') then
                    ICOUNT(14,JFAM)=ICOUNT(14,JFAM)+1
                    CLDESC=CLITM(14)
                    if (obs_headElem_i(lobsSpaceData,OBS_CHM,headerIndex).ge.0) then

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
                       CLDESC=vnl_varnameFromVarnum(ityp,obs_headElem_i(lobsSpaceData,OBS_CHM,headerIndex))
                         
                    else if (ityp == BUFR_NETT) then

                       ! CONVERT FROM KELVIN TO CELCIUS

                       CLDESC = CLITM(4)
                       ZVAR = ZVAR - MPC_K_C_DEGREE_OFFSET_R8
                       ZFCST = ZFCST - MPC_K_C_DEGREE_OFFSET_R8
                       ZANA  = ZANA - MPC_K_C_DEGREE_OFFSET_R8
                    end if
                 else if (ityp == BUFR_NEGZ) then
                   CLDESC = CLITM(3)
                   ICOUNT(3,JFAM) = ICOUNT(3,JFAM) + 1
                   !
                   ! CONVERT FROM GEOPOTENTIAL TO GEOPOTENTIAL HEIGHT
                   !
                   ZVAR = ZVAR/ec_rg
                   ZFCST = ZFCST/ec_rg
                   ZANA  = ZANA/ec_rg
                 else if (ityp == BUFR_NETT) then
                   CLDESC = CLITM(4)
                   ICOUNT(4,JFAM) = ICOUNT(4,JFAM) + 1
                   !
                   ! CONVERT FROM KELVIN TO CELCIUS
                   !
                   ZVAR = ZVAR - MPC_K_C_DEGREE_OFFSET_R8
                   ZFCST = ZFCST - MPC_K_C_DEGREE_OFFSET_R8
                   ZANA  = ZANA - MPC_K_C_DEGREE_OFFSET_R8
                 else if (ityp == BUFR_NEES) then
                   CLDESC = CLITM(5)
                   ICOUNT(5,JFAM) = ICOUNT(5,JFAM) + 1
                 else if (ityp == BUFR_NEPS) then
                   CLDESC = CLITM(6)
                   ICOUNT(6,JFAM) = ICOUNT(6,JFAM) + 1
                   !
                   ! CONVERT FROM PASCALS TO MILLIBARS
                   !
                   ZVAR = ZVAR*MPC_MBAR_PER_PA_R8
                   ZANA = ZANA*MPC_MBAR_PER_PA_R8
                   ZFCST = ZFCST*MPC_MBAR_PER_PA_R8
                 else if (ityp == BUFR_NEPN) then
                   CLDESC = CLITM(7)
                   ICOUNT(7,JFAM) = ICOUNT(7,JFAM) + 1
                   !
                   ! CONVERT FROM PASCALS TO MILLIBARS
                   !
                   ZVAR = ZVAR*MPC_MBAR_PER_PA_R8
                   ZANA = ZANA*MPC_MBAR_PER_PA_R8
                   ZFCST = ZFCST*MPC_MBAR_PER_PA_R8
                 else if (ityp == BUFR_NETS) then
                   CLDESC = CLITM(8)
                   ICOUNT(8,JFAM) = ICOUNT(8,JFAM) + 1
                   !
                   ! CONVERT FROM KELVIN TO CELCIUS
                   !
                   ZVAR = ZVAR - MPC_K_C_DEGREE_OFFSET_R8
                   ZFCST = ZFCST - MPC_K_C_DEGREE_OFFSET_R8
                   ZANA  = ZANA - MPC_K_C_DEGREE_OFFSET_R8
                 else if (ityp == BUFR_NESS) then
                   CLDESC = CLITM(9)
                   ICOUNT(9,JFAM) = ICOUNT(9,JFAM) + 1
                 else if ( ityp==BUFR_NBT1 .or. ityp==BUFR_NBT2 .or.  &
                      ityp==BUFR_NBT3      ) then
                   CLDESC = CLITM(12)
                   ICOUNT(12,JFAM) = ICOUNT(12,JFAM) + 1
                 else if (ityp == BUFR_NEZD) then
                   CLDESC = CLITM(13)
                   ICOUNT(13,JFAM) = ICOUNT(13,JFAM) + 1
                   !
                   ! CONVERT ZTD FROM METRES TO MILLIMETRES
                   !
                   ZVAR = ZVAR * 1000.D0
                   ZFCST = ZFCST * 1000.D0
                   ZANA  = ZANA * 1000.D0
                 else if (ityp == bufr_sst) then
                   CLDESC = CLITM(15)
                   ICOUNT(15,JFAM) = ICOUNT(15,JFAM) + 1
                 else if (ityp == BUFR_ICEC) then
                   CLDESC = CLITM(16)
                   ICOUNT(16,JFAM) = ICOUNT(16,JFAM) + 1
                 end if

                 call obs_bodySet_i(lobsSpaceData,OBS_FLG,bodyIndex,  &
                   IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,bodyIndex),9))
                 call obs_bodySet_i(lobsSpaceData,OBS_FLG,bodyIndex,  &
                   IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,bodyIndex),17))
                 codtypname=codtyp_get_name(IDBURP)
                 stnId = obs_elem_c(lobsSpaceData,'STID',headerIndex)
                 obsONM = obs_headElem_i(lobsSpaceData,OBS_ONM,headerIndex)
                 write(*,630) stnId,IDBURP,  &
                   codtypname,ZLAT,ZLON,ILEV,CLUNITS,CLDESC, &
                   ZVAR,ZFCST,ZANA,ZPOST,obsONM
 630               FORMAT(A9,1X,I3,1X,A21,1X,F5.1,2X,F5.1,2X,I4,A2,2X,  &
                    A4,2X,F7.1,3X,F7.1,3X,F7.1,2X,F7.4,1X,I8,1X)
               end if
             end if
           end if
        end do BODY
        !
        ! NOW SET THE GLOBAL FLAGS IN THE BURP RECORD HEADER
        !
        if (LLELREJ) then
          call obs_headSet_i(lobsSpaceData,OBS_ST1,headerIndex, IBSET( obs_headElem_i(lobsSpaceData,OBS_ST1,headerIndex), 06))
        end if
      end do HEADER
    end do FAMILY


    write(*,640)
 640  FORMAT(//)
    write(*,*) ' REJECTED DATA ACCORDING TO FAMILY OF REPORT.'
    write(*,670)(ofl_familyList(JFAM),JFAM=1,ofl_numFamily)
 670  FORMAT(5X,11(7X,A2))
    do JITEM=1,NUMITEM
      if ( .NOT. (JITEM == 2 .or. JITEM == 11)) then
        write(*,680)CLITM(JITEM),(ICOUNT(JITEM,JFAM),JFAM=1,ofl_numFamily)
 680      FORMAT(1X,A4,11(4X,I5))
      end if
    end do
    write(*,640)

  end subroutine vqc_listrej


end module varQC_mod

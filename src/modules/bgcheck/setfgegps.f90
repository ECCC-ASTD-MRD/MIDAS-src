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

!--------------------------------------------------------------------------
!! *Purpose*: Set FGE for all GPS ZTD observations using
!!            Jacobians from ZTD observation operator
!!
!! OPTION: Test ZTD operators (compares H(x+dx)-H(x) with (dH/dx)*dx
!!         when LTESTOP = .true.)
!!
!! @author S. Macpherson *ARMA/MSC  December 2004
!!
!! Revisions:
!!
!!        -S. Macpherson *ARMA/MSC  18 March 2010
!!           - add optional NL, TL and AD operator tests
!!        -S. Macpherson *ARMA/MSC   August 2010
!!           - use new GPS ZTD observation operator (from GPS-RO modules)
!!        -S. Macpherson *ARMA/MSC   December 2012
!!           - update from Rev189 to Rev213
!!           - use new ZTD-specific GPS modules modgps04profilezd, modgps08ztdop
!!           - LTESTOP option now set in 3dvar namelist
!!           - if numGPSZTD=0, does nothing and returns
!!        -S. Macpherson *ARMA/MSC   November 2014
!!           - add surface pressure (P0) argument to call gps_structztd()
!!        -M. Bani Shahabadi Dec 2018
!!           - use the calculated height in tt2phi in the gps_structztd_v2()
!!
!!v     *********************************************************************
!!v     ****                   9 October 2015                            ****
!!v     ****                                                             ****
!!v     **** NOTE: Effective Rev644M, this routine is no longer used!    ****
!!v     ****       FGE for ZTD is no longer needed for background check. ****
!!v     ****       Routine is only called when LTESTOP=.true., in which  ****
!!v     ****       case the operator test only is done.                  ****
!!v     ****                                                             ****
!!v     *********************************************************************
!!
!--------------------------------------------------------------------------
      SUBROUTINE SETFGEGPS(lcolumn,lcolumng,lobsSpaceData)
      use EarthConstants_mod
      use MathPhysConstants_mod
      use bufr_mod
      use columnData_mod
      use obsSpaceData_mod
      use verticalCoord_mod
      use gps_mod
      IMPLICIT NONE
!   lcolumn  contains background errors for control variables on model levels
!   lcolumng contains lo-res first guess profiles at obs locations
      type(struct_columnData) :: lcolumn, lcolumng
      type(struct_obs) :: lobsSpaceData
      type(struct_vco), pointer :: vco_anl
      REAL*8 ZLAT, Lat
      REAL*8 ZLON, Lon
      REAL*8, allocatable :: ZPP(:)
      REAL*8, allocatable :: ZDP(:)
      REAL*8, allocatable :: ZTT(:)
      REAL*8, allocatable :: ZHU(:)
      REAL*8, allocatable :: ZGZ(:)
      REAL*8, allocatable :: ZGZ2(:)
      REAL*8, allocatable :: ZTTB(:)
      REAL*8, allocatable :: ZHUB(:)
      REAL*8, allocatable :: ZQQB(:)
      REAL*8, allocatable :: ZQQ(:)
      REAL*8, allocatable :: ZTTB_P(:)
      REAL*8, allocatable :: ZQQB_P(:)
      REAL*8, allocatable :: ZGZ_P(:)
      REAL*8, allocatable :: RZHUB_P(:)
      REAL*8, allocatable :: ZPP_P(:)
      
      REAL*8 ZP0
      REAL*8 ZP0B, ZP0B_P
      REAL*8 ZMT, ZTOP, ZBOT
 
      REAL*8 JAC(ngpscvmx)
      REAL*8 DX (ngpscvmx)

      REAL*8 ZOER, ZLEV, ZTDOBS, ZVAR, ZPSMOD
      REAL*8 ZJP0, ZLSUM
      REAL*8 DELTAH_NL, DELTAH_TL
      REAL*8 PERTFAC, ZTDM
      REAL*8 ZDZMIN, ZSUMTEST

      INTEGER INDEX_HEADER, FIRST_HEADER
      INTEGER IDATYP, ITYP
      INTEGER IDATA, IDATEND, INDEX_BODY
      INTEGER JL, JK, NFLEV_T, ILYR, IOBS
      INTEGER INOBS_OPT, INOBS_JAC, icount, status, iversion

      LOGICAL  ASSIM, LLOK, LSTAG
      CHARACTER*9  STN_JAC
      
      CHARACTER(len=2) :: varLevel
      
      TYPE(GPS_PROFILEZD)    :: PRF, PRFP
      TYPE(GPS_DIFF)         :: ZTDopv, ZTDopvP

      IF (numGPSZTD .EQ. 0) RETURN

!C
!C     * 1.  Initializations
!C     *     ---------------
!C
      NFLEV_T = col_getNumLev(lcolumng,'TH')
      allocate(ZPP(NFLEV_T))
      allocate(ZDP(NFLEV_T))
      allocate(ZTT(NFLEV_T))
      allocate(ZHU(NFLEV_T))
      allocate(ZGZ(NFLEV_T))
      allocate(ZTTB(NFLEV_T))
      allocate(ZHUB(NFLEV_T))
      allocate(ZQQB(NFLEV_T))
      allocate(ZQQ(NFLEV_T))

!c     Number of locations/sites for observation operator test
      INOBS_OPT = 50
!c     Number of locations/sites for Jacobian printout
      INOBS_JAC  = 5
!c     Factor to multiply background errors for perturbation vector
      PERTFAC = 0.75d0
!C
      STN_JAC = 'FSL_BRFT '
!c      
      ZDZMIN = DZMIN
!c
      vco_anl => col_getVco(lcolumng)
      status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=iversion)
      if (iversion .eq. 5002) then
         LSTAG = .TRUE. 
         WRITE(*,*)'VERTICAL COORD OF ANALYSIS FIELDS IS STAGGERED'
         WRITE(*,*)'VCODE= ',iversion,' LSTAG= ',LSTAG
      else
         LSTAG = .FALSE.
         WRITE(*,*)'VERTICAL COORD OF ANALYSIS FIELDS IS NOT STAGGERED'
         WRITE(*,*)'VCODE= ',iversion,' LSTAG= ',LSTAG
      endif

!C
      IF ( .NOT.LTESTOP ) THEN

      first_header=-1
      icount = 0
!C
      ! loop over all header indices of the 'GP' family
      call obs_set_current_header_list(lobsSpaceData,'GP')
      HEADER: do
         index_header = obs_getHeaderIndex(lobsSpaceData)
         if (index_header < 0) exit HEADER
         if (first_header .eq. -1) first_header = index_header
!C     
!C     *    .   Process only zenith delay data (codtyp 189 and BUFR_NEZD)
!C
               IDATYP = obs_headElem_i(lobsSpaceData,OBS_ITY,INDEX_HEADER)
               IF ( IDATYP .EQ. 189 ) THEN
!C
!C                 Loop over data in the observations
!C
                  IDATA   = obs_headElem_i(lobsSpaceData,OBS_RLN,INDEX_HEADER)
                  IDATEND = obs_headElem_i(lobsSpaceData,OBS_NLV,INDEX_HEADER) + IDATA - 1
                  ASSIM = .FALSE.
!C
!C                 Scan for requested assimilations, and count them.
!C
                  DO INDEX_BODY= IDATA, IDATEND
                     ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,INDEX_BODY)
                     LLOK = ( (ITYP .EQ. BUFR_NEZD) .AND. (obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY) .EQ. 1) )
                     IF ( LLOK ) THEN
                        ASSIM = .TRUE.
                        ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)
                        icount = icount + 1
                     ENDIF
                  ENDDO
!C
!C     *           If assimilations are requested, apply the AD observation operator
!C
                  IF (ASSIM) THEN
!C     
!C     *        LR background profile and background errors at the observation location x :
!C
                     Lat  = obs_headElem_r(lobsSpaceData,OBS_LAT,INDEX_HEADER)
                     Lon  = obs_headElem_r(lobsSpaceData,OBS_LON,INDEX_HEADER)
                     ZLAT = Lat * MPC_DEGREES_PER_RADIAN_R8
                     ZLON = Lon * MPC_DEGREES_PER_RADIAN_R8
                     ZP0B = col_getElem(lcolumng,1,INDEX_HEADER,'P0')
                     DO JL = 1, NFLEV_T
                       ZPP(JL)  = col_getPressure(lcolumng,JL,INDEX_HEADER,'TH')
!C                     Get ZDP = dP/dP0
                       ZDP(JL)  = col_getPressureDeriv(lcolumng,JL,INDEX_HEADER,'TH')
                       ZTTB(JL) = col_getElem(lcolumng,JL,INDEX_HEADER,'TT')- 273.15d0
                       ZTT(JL)  = col_getElem(lcolumn,JL,INDEX_HEADER,'TT')
                       DX(JL)   = ZTT(JL)
                       ZHUB(JL) = col_getElem(lcolumng,JL,INDEX_HEADER,'HU')
                       ZQQB(JL) = ZHUB(JL)
                       ZHU(JL)  = col_getElem(lcolumn,JL,INDEX_HEADER,'HU')
                       DX(NFLEV_T+JL) = ZHU(JL)
                       ZGZ(JL)  = col_getHeight(lcolumng,JL,INDEX_HEADER,'TH')
                       DX(2*NFLEV_T+JL) = col_getHeight(lcolumn,JL,INDEX_HEADER,'TH')
                     ENDDO
                     ZP0  = col_getElem(lcolumn,1,INDEX_HEADER,'P0')
                     DX(3*NFLEV_T+1) = ZP0
                     ZMT  = ZGZ(NFLEV_T)
                     CALL gps_structztd_v2(NFLEV_T,Lat,Lon,ZMT,ZP0B,ZPP,ZDP,ZTTB,ZHUB,ZGZ,LBEVIS,IREFOPT,PRF)
                     CALL gps_ztdopv(ZLEV,PRF,LBEVIS,ZDZMIN,ZTDopv,ZPSMOD,IZTDOP)
                     JAC = ZTDopv%DVar
!c
                     DO INDEX_BODY= IDATA, IDATEND
                        ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,INDEX_BODY)
                        IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY).EQ.1 .AND. ITYP.EQ.BUFR_NEZD ) THEN
!C
!C     *                    Observation error    SDERR
!c                           ZOER = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)

!C     *                    Observation height (m)
                           ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)

                           ZLSUM  = 0.0d0
!C
                           DO JL = 1, 3*NFLEV_T+1
                             ZLSUM = ZLSUM + (JAC(JL)*DX(JL))**2
                           ENDDO
                           call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,SQRT(ZLSUM))

                           IF (icount .LE. INOBS_JAC) THEN
!                           IF ( obs_elem_c(lobsSpaceData,'STID',INDEX_HEADER) .EQ. STN_JAC ) THEN
                             WRITE(*,'(A11,A9)') 'SETFGEGPS: ',obs_elem_c(lobsSpaceData,'STID',INDEX_HEADER)
                             WRITE(*,*) '  ZTD, ZTD FGE = ', ZTDopv%Var, SQRT(ZLSUM)
                             WRITE(*,'(A11,A9,3(1x,f7.2))')   &
                               'SETFGEGPS: ',obs_elem_c(lobsSpaceData,'STID',INDEX_HEADER),ZLAT,ZLON,ZLEV
                             WRITE(*,*) 'JL JACT JACQ FGE_T FGE_LQ QQ'
                             DO JL = 1, NFLEV_T
                               WRITE(*,'(1X,I2,5(1x,E13.6))') JL,JAC(JL),JAC(JL+NFLEV_T)/ZQQB(JL),ZTT(JL),ZHU(JL),ZQQB(JL)
                             ENDDO                         
                             WRITE(*,*) 'JACPS FGE_PS'
                             WRITE(*,'(2(1x,E13.6))') JAC(3*NFLEV_T+1), ZP0
                           ENDIF

                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF

      ENDDO HEADER
      
      ENDIF

!c---------------------------------------------------------------------------------------------------------------

      IF ( LTESTOP ) THEN
      
      allocate(ZTTB_P(NFLEV_T))
      allocate(ZQQB_P(NFLEV_T))
      allocate(ZGZ2(NFLEV_T))
      allocate(ZGZ_P(NFLEV_T))
      allocate(ZPP_P(NFLEV_T))

      icount = 0
      ZSUMTEST = 0
!C
      ! loop over all header indices of the 'GP' family
      call obs_set_current_header_list(lobsSpaceData,'GP')
      HEADER2: DO
         index_header = obs_getHeaderIndex(lobsSpaceData)
         if (index_header < 0) exit HEADER2
         if (icount > INOBS_OPT ) exit HEADER2

!C       Loop over data in the observations
!C
         IDATA   = obs_headElem_i(lobsSpaceData,OBS_RLN,INDEX_HEADER)
         IDATEND = obs_headElem_i(lobsSpaceData,OBS_NLV,INDEX_HEADER) + IDATA - 1
!C     
!C       LR background profile and background errors at the observation location x :
!C
         Lat  = obs_headElem_r(lobsSpaceData,OBS_LAT,INDEX_HEADER)
         Lon  = obs_headElem_r(lobsSpaceData,OBS_LON,INDEX_HEADER)
         ZLAT = Lat * MPC_DEGREES_PER_RADIAN_R8
         ZLON = Lon * MPC_DEGREES_PER_RADIAN_R8
         ZP0B = col_getElem(lcolumng,1,INDEX_HEADER,'P0')
         DO JL = 1, NFLEV_T
            ZPP(JL)  = col_getPressure(lcolumng,JL,INDEX_HEADER,'TH')
!C          Get ZDP = dP/dP0
            ZDP(JL)  = col_getPressureDeriv(lcolumng,JL,INDEX_HEADER,'TH')
            ZTTB(JL) = col_getElem(lcolumng,JL,INDEX_HEADER,'TT')- 273.15d0
            ZTT(JL)  = col_getElem(lcolumn,JL,INDEX_HEADER,'TT') * PERTFAC
            ZQQB(JL) = col_getElem(lcolumng,JL,INDEX_HEADER,'HU')
            ZQQ(JL)  = col_getElem(lcolumn,JL,INDEX_HEADER,'HU') * PERTFAC
            ZGZ(JL)  = col_getHeight(lcolumng,JL,INDEX_HEADER,'TH')
            ZGZ2(JL)  = col_getHeight(lcolumn,JL,INDEX_HEADER,'TH') * PERTFAC
         ENDDO
         ZP0  = col_getElem(lcolumn,1,INDEX_HEADER,'P0') * PERTFAC
         ZMT  = ZGZ(NFLEV_T)

         DO JL = 1, NFLEV_T
             DX (      JL) = ZTT(JL)
             DX (NFLEV_T+JL) = ZQQ(JL)
             DX (2*NFLEV_T+JL) = ZGZ2(JL)
         ENDDO
         DX (3*NFLEV_T+1) = ZP0

         ZTDOBS = -1.0d0
         DO INDEX_BODY = IDATA, IDATEND
           ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,INDEX_BODY)
           IOBS = obs_bodyElem_i(lobsSpaceData,OBS_HIND,INDEX_BODY)
           IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY).EQ.1 .AND. ITYP .EQ. BUFR_NEZD ) THEN
             varLevel = vnl_varLevelFromVarnum(ITYP)
             ZTDOBS  = obs_bodyElem_r(lobsSpaceData,OBS_VAR,INDEX_BODY)
             ZLEV    = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)
             ILYR    = obs_bodyElem_i(lobsSpaceData,OBS_LYR,INDEX_BODY)
             ZTOP    = col_getHeight(lcolumng,ILYR,IOBS,varLevel)
             if ( ILYR .LT. NFLEV_T ) then
               ZBOT    = col_getHeight(lcolumng,ILYR+1,IOBS,varLevel)
             else
               ZBOT    = ZTOP
             endif
             icount  = icount + 1
           ENDIF
         ENDDO

         IF ( ZTDOBS .GT. 0.d0 ) THEN
!c         Create the pertubation control vector
           DO JL = 1, NFLEV_T
             ZPP_P(JL)  = ZPP(JL)  + ZDP(JL)*ZP0
             ZTTB_P(JL) = ZTTB(JL) + ZTT(JL)
             ZQQB_P(JL) = ZQQB(JL) + ZQQ(JL)
             ZGZ_P(JL) = ZGZ(JL) + ZGZ2(JL)
           ENDDO
           ZP0B_P = ZP0B + ZP0
!C
!C         Non-linear observation operator --> delta_H = H(x+delta_x) - H(x)
!c
           CALL gps_structztd_v2(NFLEV_T,Lat,Lon,ZMT,ZP0B,ZPP,ZDP,ZTTB,ZQQB,ZGZ,LBEVIS,IREFOPT,PRF)
           CALL gps_structztd_v2(NFLEV_T,Lat,Lon,ZMT,ZP0B_P,ZPP_P,ZDP,ZTTB_P,ZQQB_P,ZGZ_P,LBEVIS,IREFOPT,PRFP)
           CALL gps_ztdopv(ZLEV,PRF,LBEVIS,ZDZMIN,ZTDopv,ZPSMOD,IZTDOP)
           JAC  = ZTDopv%DVar
           ZTDM = ZTDopv%Var
           CALL gps_ztdopv(ZLEV,PRFP,LBEVIS,ZDZMIN,ZTDopvP,ZPSMOD,IZTDOP)
           DELTAH_NL = ZTDopvP%Var - ZTDopv%Var
!c
!c         Linear  --> delta_H = dH/dx * delta_x
!c
           DELTAH_TL = 0.0d0
           DO JL = 1, 3*NFLEV_T+1
             DELTAH_TL = DELTAH_TL + JAC(JL)*DX(JL)
           ENDDO
!c
           WRITE(*,*) 'SETFGEGPS: GPS ZTD OBSOP TEST FOR SITE ', obs_elem_c(lobsSpaceData,'STID',INDEX_HEADER)
           WRITE(*,*) ' '
           WRITE(*,*) '  DZ (M), MODEL LEVEL ABOVE = ', ZLEV-ZMT, ILYR
           WRITE(*,*) '  ZLEV (M), ZTOP (M), ZBOT (M) = ', ZLEV, ZTOP, ZBOT
           WRITE(*,*) '  ZTD OBS (MM)            = ', ZTDOBS*1000.d0
           WRITE(*,*) '  ZTD_MOD                 = ', ZTDM*1000.d0
           WRITE(*,*) '  DELTAH_NL, DELTAH_TL = ', DELTAH_NL*1000.d0, DELTAH_TL*1000.d0
           WRITE(*,*) ' '
           WRITE(*,*) '  DELTAH_TL/DELTAH_NL = ', DELTAH_TL/DELTAH_NL
           WRITE(*,*) ' '  
           
           ZSUMTEST = ZSUMTEST + (DELTAH_TL/DELTAH_NL)
           
         ENDIF

      ENDDO HEADER2
      
      WRITE(*,*) ' '
      WRITE(*,*) 'SETFGEGPS: ----- GPS ZTD OBSOP TEST SUMMARY -----'
      WRITE(*,*) '           NUMBER OF TESTS (sites) = ', icount
      WRITE(*,*) '           AVG DELTAH_TL/DELTAH_NL = ', ZSUMTEST/FLOAT(icount)
      WRITE(*,*) ' '  

      deallocate(ZTTB_P)
      deallocate(ZQQB_P)
      deallocate(ZGZ2)
      deallocate(ZGZ_P)
      deallocate(ZPP_P)

      ENDIF
!-----------------------------------------------------------------------------------------------------------

      deallocate(ZPP)
      deallocate(ZDP)
      deallocate(ZTT)
      deallocate(ZHU)
      deallocate(ZGZ)
      deallocate(ZTTB)
      deallocate(ZHUB)
      deallocate(ZQQB)
      deallocate(ZQQ)

      RETURN
      END SUBROUTINE SETFGEGPS

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
!! *Purpose*: Construct the FIRST GUESS ERROR VARIANCES from the
!!            diff-calculated dependencies and the primary errors.
!!
!! @author J.M. Aparicio *MSC/ARMA Nov 2004
!!
!!          Adapted Nov 2012 for both refractivity and bending angle data
!!
!--------------------------------------------------------------------------
      SUBROUTINE SETFGEDIF(CDFAM,lcolumng,lobsSpaceData)
      use EarthConstants_mod
      use MathPhysConstants_mod
      use gps_mod
      use obsSpaceData_mod
      use columnData_mod
      IMPLICIT NONE
!C
      type(struct_columnData) :: lcolumng
      type(struct_obs)        :: lobsSpaceData
!C
      INTEGER INDEX_HEADER, IDATYP, INDEX_BODY, iProfile
      CHARACTER*2 CDFAM
      REAL*8 zLat, Lat, sLat
      REAL*8 zLon, Lon
      REAL*8 zAzm, Azm
      INTEGER IAZM, ISAT
      REAL*8 Rad, Geo, WFGPS
      REAL*8, allocatable :: zPP(:)
      REAL*8, allocatable :: zDP(:)
      REAL*8, allocatable :: zTT(:)
      REAL*8, allocatable :: zHU(:)
      REAL*8, allocatable :: zUU(:)
      REAL*8, allocatable :: zVV(:)
      INTEGER JF, stat
      INTEGER JL, JJ
      REAL*8 ZP0, ZMT
      REAL*8 HNH1, ZFGE, ZERR
      INTEGER JV, NGPSLEV, NWNDLEV
      LOGICAL  ASSIM, LFIRST, FIRSTHEADER

      INTEGER NH, NH1

!      REAL*8 JAC(ngpscvmx)
      REAL*8 DV (ngpscvmx)
      TYPE(GPS_PROFILE)           :: PRF
      REAL*8       , allocatable :: H   (:),AZMV(:)
      TYPE(GPS_DIFF), allocatable :: RSTV(:),RSTVP(:),RSTVM(:)
      type(struct_vco), pointer  :: vco_anl

      WRITE(*,*)'ENTER SETFGEDIFF'
!C
!C     * 1.  Initializations
!C     *     ---------------
!C
      NGPSLEV=col_getNumLev(lcolumng,'TH')
      NWNDLEV=col_getNumLev(lcolumng,'MM')
      LFIRST=.FALSE.
      if ( .NOT.allocated(gps_vRO_Jacobian) ) then
         LFIRST = .TRUE.
         allocate(zPP (NGPSLEV))
         allocate(zDP (NGPSLEV))
         allocate(zTT (NGPSLEV))
         allocate(zHU (NGPSLEV))
         allocate(zUU (NGPSLEV))
         allocate(zVV (NGPSLEV))

         allocate(gps_vRO_Jacobian(gps_numROProfiles,GPSRO_MAXPRFSIZE,2*NGPSLEV+1))
         allocate(gps_vRO_lJac    (gps_numROProfiles))
         gps_vRO_lJac=.false.

         allocate( H    (GPSRO_MAXPRFSIZE) )
         allocate( AZMV (GPSRO_MAXPRFSIZE) )
         allocate( RSTV (GPSRO_MAXPRFSIZE) )
!C         IF (LEVELGPSRO.EQ.1) THEN
!C            allocate( RSTVP(GPSRO_MAXPRFSIZE) )
!C            allocate( RSTVM(GPSRO_MAXPRFSIZE) )
!C         ENDIF
      endif

      vco_anl => col_getVco(lcolumng)
!C
!C    Loop over all header indices of the 'RO' family:
!C
      call obs_set_current_header_list(lobsSpaceData,CDFAM)
      FIRSTHEADER=.TRUE.
      HEADER: do
         index_header = obs_getHeaderIndex(lobsSpaceData)
         if (index_header < 0) exit HEADER
!C
!C     * Process only refractivity data (codtyp 169)
!C
         IDATYP = obs_headElem_i(lobsSpaceData,OBS_ITY,INDEX_HEADER)
         IF ( IDATYP .EQ. 169 ) THEN
!C
!C     *    Scan for requested data values of the profile, and count them
!C
            ASSIM = .FALSE.
            NH = 0
            call obs_set_current_body_list(lobsSpaceData, INDEX_HEADER)
            BODY: do
               INDEX_BODY = obs_getBodyIndex(lobsSpaceData)
               if (INDEX_BODY < 0) exit BODY
               IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY).EQ.1 ) THEN
                  ASSIM = .TRUE.
                  NH = NH + 1
               ENDIF
            ENDDO BODY
!C
!C     *    If assimilations are requested, prepare and apply the observation operator
!C
            IF (ASSIM) THEN
               iProfile=gps_iprofile_from_index(INDEX_HEADER)
!C
!C     *       Profile at the observation location:
!C
               if (.not.gps_vRO_lJac(iProfile)) then
!C
!C     *          Basic geometric variables of the profile:
!C
                  zLat = obs_headElem_r(lobsSpaceData,OBS_LAT,INDEX_HEADER)
                  zLon = obs_headElem_r(lobsSpaceData,OBS_LON,INDEX_HEADER)
                  IAZM = obs_headElem_i(lobsSpaceData,OBS_AZA,INDEX_HEADER)
                  ISAT = obs_headElem_i(lobsSpaceData,OBS_SAT,INDEX_HEADER)
                  Rad  = obs_headElem_r(lobsSpaceData,OBS_TRAD,INDEX_HEADER)
                  Geo  = obs_headElem_r(lobsSpaceData,OBS_GEOI,INDEX_HEADER)
                  zAzm = 0.01d0*IAZM / MPC_DEGREES_PER_RADIAN_R8
                  zMT  = col_getHeight(lcolumng,NGPSLEV,INDEX_HEADER,'TH')
                  WFGPS= 0.d0
                  DO JJ=1,NUMGPSSATS
                     IF (ISAT.EQ.IGPSSAT(JJ)) WFGPS=WGPS(JJ)
                  ENDDO
                  Lat  = zLat * MPC_DEGREES_PER_RADIAN_R8
                  Lon  = zLon * MPC_DEGREES_PER_RADIAN_R8
                  Azm  = zAzm * MPC_DEGREES_PER_RADIAN_R8
                  sLat = sin(zLat)
                  zMT  = zMT * RG / gpsgravitysrf(sLat)
                  zP0  = col_getElem(lcolumng,1,INDEX_HEADER,'P0')
                  DO JL = 1, NGPSLEV
!C
!C     *             Profile x
!C
                     zPP(JL) = col_getPressure(lcolumng,JL,INDEX_HEADER,'TH')
!C     *             True implementation of zDP (dP/dP0)
                     zDP(JL) = col_getPressureDeriv(lcolumng,JL,INDEX_HEADER,'TH')
                     zTT(JL) = col_getElem(lcolumng,JL,INDEX_HEADER,'TT') - p_TC
                     zHU(JL) = col_getElem(lcolumng,JL,INDEX_HEADER,'HU')
                     zUU(JL) = 0.d0
                     zVV(JL) = 0.d0
                  ENDDO
                  DO JL = 1, NWNDLEV
                     zUU(JL) = col_getElem(lcolumng,JL,INDEX_HEADER,'UU')
                     zVV(JL) = col_getElem(lcolumng,JL,INDEX_HEADER,'VV')
                  ENDDO
                  zUU(NGPSLEV) = zUU(NWNDLEV)
                  zVV(NGPSLEV) = zUU(NWNDLEV)
!C     
!C     *          GPS profile structure:
!C
                  call gps_struct1sw(ngpslev,zLat,zLon,zAzm,zMT,Rad,geo,zP0,zPP,zDP,zTT,zHU,zUU,zVV,prf)
!C
!C     *          Prepare the vector of all the observations:
!C
                  NH1 = 0
                  call obs_set_current_body_list(lobsSpaceData, INDEX_HEADER)
                  BODY_2: do
                     INDEX_BODY = obs_getBodyIndex(lobsSpaceData)
                     if (INDEX_BODY < 0) exit BODY_2
                     IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY).EQ.1 ) THEN
                        NH1      = NH1 + 1
                        H(NH1)   = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)
                        AZMV(NH1)= zAzm
                     ENDIF
                  ENDDO BODY_2
!C
!C     *          Apply the observation operator:
!C
                  IF (LEVELGPSRO.EQ.1) THEN
                     CALL GPS_BNDOPV1(H      , AZMV, NH, PRF, RSTV)
!C                     CALL GPS_BNDOPV1(H+WFGPS, AZMV, NH, PRF, RSTVP)
!C                     CALL GPS_BNDOPV1(H-WFGPS, AZMV, NH, PRF, RSTVM)
!C                     do nh1 = 1, nh
!C                        RSTV(nh1)=(RSTVP(nh1)+RSTV(nh1)+RSTVM(nh1))/3.d0
!C                     enddo
                  ELSE
                     CALL GPS_REFOPV (H,       NH, PRF, RSTV)
                  ENDIF
                  DO NH1=1,NH
                     gps_vRO_Jacobian(iProfile,NH1,:)= RSTV(NH1)%DVAR(1:2*NGPSLEV+1)
                  ENDDO
                  gps_vRO_lJac(iProfile)=.true.
               endif
!C
!C     *       Local error
!C
               DO JL = 1, NGPSLEV
                  DV (        JL) = 1.d0
                  DV (NGPSLEV+JL) = 1.d0
               ENDDO
               DV (2*NGPSLEV+1)   = 2.d0
!C
!C     *       Perform the H(xb)DV operation:
!C
               NH1 = 0
               call obs_set_current_body_list(lobsSpaceData, INDEX_HEADER)
               BODY_3: do
                  INDEX_BODY = obs_getBodyIndex(lobsSpaceData)
                  if (INDEX_BODY < 0) exit BODY_3
                  IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY).EQ.1 ) THEN
                     NH1 = NH1 + 1
!C
!C     *             Observation jacobian
!C
!                     JAC = RSTV(NH1)%DVAR
!C
!C     *             Evaluate sqrt( H(xb)DV **2 )
!C
                     ZFGE = 0.d0
                     DO JV = 1, 2*PRF%NGPSLEV+1
                        ZFGE = ZFGE + (gps_vRO_Jacobian(iProfile,NH1,JV) * DV(JV))**2
                     ENDDO
                     ZFGE = SQRT(ZFGE)
                     ZERR = obs_bodyElem_r(lobsSpaceData,OBS_OER,INDEX_BODY)
!C     
!C     *             FIRST GUESS ERROR VARIANCE
!C
                     call obs_bodySet_r(lobsSpaceData,OBS_HPHT,INDEX_BODY,ZFGE)
                     IF (FIRSTHEADER) THEN
11                      FORMAT(A12,2I5,F12.2,3F16.8)
                        WRITE(*,11)'SETFGEDIFFGE',NH1,NH,H(NH1),RSTV(NH1)%VAR,ZFGE,ZERR
                     ENDIF
                  ENDIF
               ENDDO BODY_3
            ENDIF
         ENDIF
         FIRSTHEADER = .FALSE.
      ENDDO HEADER

      IF (LFIRST) THEN
!C         IF (LEVELGPSRO.EQ.1) THEN
!C            deallocate( RSTVM )
!C            deallocate( RSTVP )
!C         ENDIF
         deallocate( RSTV )
         deallocate( AZMV )
         deallocate( H    )

         deallocate(zVV)
         deallocate(zUU)
         deallocate(zHU)
         deallocate(zTT)
         deallocate(zDP)
         deallocate(zPP)
         deallocate(gps_vRO_Jacobian)
      ENDIF

      WRITE(*,*)'EXIT SETFGEDIFF'
      RETURN
      END SUBROUTINE SETFGEDIF

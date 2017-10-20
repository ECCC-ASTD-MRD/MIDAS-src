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
!!  *Purpose*: Set background check flag on GPSRO data if ABS(O-P)/P is too large
!!
!!  @author P. KOCLAS. Mar 2008
!!  Modified: 
!!            -J.M. Aparicio, Dec 2012.
!!                -Simplified and adapted to both refractivity and bending angle data
!!
!--------------------------------------------------------------------------
      SUBROUTINE BGCGPSRO(lcolumnhr,lobsSpaceData)
      use EarthConstants_mod
      use MathPhysConstants_mod
      use obsSpaceData_mod
      use columnData_mod
      use verticalCoord_mod
      use gps_mod
      IMPLICIT NONE

      type(struct_columnData) :: lcolumnhr
      type(struct_obs) :: lobsSpaceData
      type(struct_vco), pointer :: vco_trl
      REAL*8 HNH1, ZOBS, ZMHX, ZOMF, ZREF, ZOER, Rad

      INTEGER INDEX_HEADER, JD
      INTEGER IDATYP
      INTEGER IDATA   , IDATEND, INDEX_BODY
      INTEGER JL, JH,  NGPSLEV
      INTEGER NH, NH1, stat, iversion

      LOGICAL  LSTAG


      LSTAG = .FALSE.
      vco_trl => col_getVco(lcolumnhr)
      stat = vgd_get(vco_trl%vgrid,key='ig_1 - vertical coord code',value=iversion)
      if (iversion .eq. 5002) LSTAG = .TRUE. 
      
      WRITE(*,*)'ENTER BGCSGPSRO'
!C
!C     * 1.  Initializations
!C     *     ---------------
!C
      NGPSLEV=col_getNumLev(lcolumnhr,'TH')
!C
!C     Loop over all files
!C
      ! loop over all header indices of the 'RO' family
      call obs_set_current_header_list(lobsSpaceData,'RO')
      HEADER: do
         index_header = obs_getHeaderIndex(lobsSpaceData)
         if (index_header < 0) exit HEADER
!C     
!C     * Process only refractivity data (codtyp 169)
!C
         IDATYP = obs_headElem_i(lobsSpaceData,OBS_ITY,INDEX_HEADER)
         IF ( IDATYP .EQ. 169 ) THEN
!C
!C     *    Basic geometric variables of the profile:
!C
            Rad  = obs_headElem_r(lobsSpaceData,OBS_TRAD,INDEX_HEADER)
!C
!C     *    Loops over data in the observation
!C
            IDATA   = obs_headElem_i(lobsSpaceData,OBS_RLN,INDEX_HEADER)
            IDATEND = obs_headElem_i(lobsSpaceData,OBS_NLV,INDEX_HEADER) + IDATA - 1
!C
!C     *    Scan for requested assimilations, and count them
!C
            DO INDEX_BODY= IDATA, IDATEND
               IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY).EQ.1 ) THEN
                  HNH1 = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)
                  IF (LEVELGPSRO.EQ.1) HNH1 = HNH1-Rad
!C
!C     *          Increment OMF = Y - H(x)
!C
                  ZOMF = obs_bodyElem_r(lobsSpaceData,OBS_OMP,INDEX_BODY)
!C
!C     *          Observation value    Y
!C
                  ZOBS = obs_bodyElem_r(lobsSpaceData,OBS_VAR,INDEX_BODY)
                  ZOER = obs_bodyElem_r(lobsSpaceData,OBS_OER,INDEX_BODY)
                  ZMHX = ZOBS-ZOMF
!C
!C     *          Reference order of magnitude value:
!C
                  IF (LEVELGPSRO.EQ.1) THEN
                     ZREF = 0.025d0*exp(-HNH1/6500.d0)
                  ELSE
                     IF (NUMGPSSATS .GE. 1) THEN
                       ZREF = 300.d0*exp(-HNH1/6500.d0)
                     ELSE
                       ZREF = ZMHX
                     ENDIF
                  ENDIF
!C                           
!C     *          OMF Tested criteria:
!C
                  IF (NUMGPSSATS .GE. 1) THEN
                    IF (DABS(ZOMF)/ZREF.GT.BGCKBAND .OR. DABS(ZOMF)/ZOER.GT.3.d0) THEN
                      call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body),16))
                      call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body),9))
                    ENDIF
                  ELSE
                    IF (DABS(ZOMF)/ZREF.GT.BGCKBAND) THEN
                      call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body),16))
                      call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body),9))
                      WRITE(*,'(A40,F10.0,3F12.4)') ' REJECT BGCSGPSRO H  O  P (O-P/ZREF) =',HNH1,ZOBS,ZMHX,(ZOMF)/ZREF
                    ENDIF                  
                  ENDIF
                  
               ENDIF
            ENDDO
         ENDIF

      ENDDO HEADER

      WRITE(*,*)'EXIT BGCSGPSRO'
      RETURN
      
      END SUBROUTINE BGCGPSRO
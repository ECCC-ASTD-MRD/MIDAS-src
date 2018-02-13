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
!! *Purpose*:  Calculate a background check for a data family
!!             AND SET the appropriate quality control flags in
!!             obsSpaceData
!!
!! @author P. Koclas *CMC/CMDA  January 1999
!!
!! Revisions:
!!
!!        -S. Macpherson *ARMA/MRD  9 October 2015
!!           - Changed background check for GB-GPS ZTD observations:
!!v               Changed
!!v                  ZBGCHK=(ZOMP)**2/(ZFGE**2 + ZOER**2)
!!v               to
!!v                  ZBGCHK=(ZOMP**2)/(ZFGE**2) 
!!v               where ZFGE = StdDev(O-P) estimated from ZWD and set in
!!v               s/r SETERRGPSGB (stored in obsSpaceData col. OBS_HPHT).
!!
!!        -Y.J. Rochon *ARQD/ARQI March 2016
!!          - Allow ZFGE**2 + ZOER**2 .LT. 1.D-5 for CDFAM='CH'
!!          - Set output format 124 for 'CH'
!!
!!
!!        -S. Laroche  *ARMA/MRD October 2017
!!           -New background check for AMVs that checks u and v
!!            simultaneously.
!!
!--------------------------------------------------------------------------
      SUBROUTINE BGCDATA(PJO,CDFAM,lobsSpaceData,new_bgck_sw)

      use MathPhysConstants_mod
      use bufr_mod
      use obsSpaceData_mod
      use gps_mod
      use utilities_mod
      IMPLICIT NONE
!*
!C NOTE 1: YSFERRWGT IN MODGPSZTD_MOD (FROM NML FILE) IS USED HERE FOR ERROR WEIGHTING
!C         OF TIME SERIES (FGAT) GPS MET OBSERVATIONS PS, TS, DPDS. IT IS APPLIED
!C         (ALONG WITH YZDERRWGT FOR ZTD) IN S/R SETERR AS A MULT. FACTOR TO ERRORS.
!C         YZDERRWGT AND YSFERRWGT = 1 FOR NORMAL 3D-VAR (3D-THINNING).
!C

!C
      type(struct_obs) :: lobsSpaceData
      REAL*8 PJO,zsum,zsumo
      CHARACTER*2 CDFAM
      LOGICAL :: new_bgck_sw
      INTEGER IFLAG,INAM,ISETFLAG,INDEX_HEADER,IDBURP
      INTEGER IBEGIN,ILAST,ierr,fclos,fnom,ITYP
      INTEGER J,J2,JJ,JD,INDEX_BODY,icoun
      INTEGER INOBS, INREJ, INZOBS, INZREJ
      INTEGER INPOBS, INTOBS, INDOBS, INPREJ, INTREJ, INDREJ
      REAL*8 ZOER,ZOMP,ZFGE,ZBGCHK,ZVAR,ZLEV,ZLAT,ZLON,ZSOP
      LOGICAL LLOK, LLZD, LMODIF1020

      INTEGER i_ass,i_vco,i_vnm,i_oth,istyp,index_body_u,index_body_v,index_body_start
      REAL*8 uu_d,uu_r,uu_f,vv_d,vv_r,vv_f,duv2,duv2_lim,zslev
      LOGICAL found_u,found_v

      lmodif1020=.false.

      WRITE(*,*)' '
      WRITE(*,*)' ------------------------------'
      WRITE(*,*)'  BACKGROUND CHECK FOR'
      WRITE(*,*)'       ',CDFAM, ' DATA'
      WRITE(*,*)' ------------------------------'
      WRITE(*,*)' '
      WRITE(*,'(a55,a74)')'  STNID     LATITU LONGITU  ID Elem        Level        ',  &
      ' Value        Sigmao       Sigmap         O-P       SigmaOP         qcflag  '
      WRITE(*,'(a55,a74)')'  -----     ------ -------  -- ----        -----        ',  &
      ' -----        ------       ------         ---       -------         ------  '
      icoun=0
      zsum=0.D0
      zsumo=0.D0
      PJO=-99.99D0

      INOBS=0
      INREJ=0
!C
!C*    Initialize counters for GP family observations
!C
      IF (CDFAM .EQ. 'GP') THEN
        INZOBS=0
        INPOBS=0
        INTOBS=0
        INDOBS=0
        INZREJ=0
        INPREJ=0
        INTREJ=0
        INDREJ=0
      ENDIF

      IF (.not. new_bgck_sw .or. CDFAM .ne. 'SW') THEN

        ! loop over all header indices of the CDFAM family
        call obs_set_current_header_list(lobsSpaceData,CDFAM)
        HEADER: do
          index_header = obs_getHeaderIndex(lobsSpaceData)
          if (index_header < 0) exit HEADER
!C
!C*    1. Computation of (HX - Z)**2/(SIGMAo**2 +SIGMAp**2)
!C     .  ----------------------------------------------------
!C
          ! loop over all body indices for this index_header
          call obs_set_current_body_list(lobsSpaceData, index_header)
          BODY: do 
            index_body = obs_getBodyIndex(lobsSpaceData)
            if (index_body < 0) exit BODY

            ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body)
            LLOK=( obs_bodyElem_i(lobsSpaceData,OBS_ASS,index_body) .EQ. 1)
            IF ( LLOK ) THEN
              INOBS = INOBS + 1
              IF (CDFAM .EQ. 'GP') THEN
                IF (ITYP .EQ. BUFR_NEZD) INZOBS = INZOBS+1
                IF (ITYP .EQ. BUFR_NEPS) INPOBS = INPOBS+1
                IF (ITYP .EQ. BUFR_NETS) INTOBS = INTOBS+1
                IF (ITYP .EQ. BUFR_NESS) INDOBS = INDOBS+1
              ENDIF
              ZVAR = obs_bodyElem_r(lobsSpaceData,OBS_VAR,index_body)
              ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,index_body)
              ZOER = obs_bodyElem_r(lobsSpaceData,OBS_OER,index_body)
              ZLAT  = obs_headElem_r(lobsSpaceData,OBS_LAT,index_header) &
                                                  * MPC_DEGREES_PER_RADIAN_R8
              ZLON  = obs_headElem_r(lobsSpaceData,OBS_LON,index_header) &
                                                  * MPC_DEGREES_PER_RADIAN_R8
!C
!C            BACK GROUND CHECK
!C
              ZOMP  = obs_bodyElem_r(lobsSpaceData,OBS_OMP,index_body)
!C            NOTE: For GB-GPS ZTD observations, ZFGE is not the FGE but rather Std(O-P)
!C                  set in s/r SETERRGPSGB.
              ZFGE  = obs_bodyElem_r(lobsSpaceData,OBS_HPHT,index_body)
              IF ( CDFAM .EQ. 'GP' ) THEN
                IF (ITYP .EQ. BUFR_NEZD) THEN
                  ZOER = ZOER/YZDERRWGT
                ELSE
                  ZOER = ZOER/YSFERRWGT
                ENDIF
                IF (ZOER .LT. 1.D-3 .AND. ITYP .NE. BUFR_NEZD) THEN
                  WRITE(*,*)' Problem for GP STNID ZOER= ' ,  &
                    obs_elem_c(lobsSpaceData,'STID',index_header),ZOER
                  CALL utl_abort('BGCDATA: PROBLEM WITH OER.')
                ENDIF
                IF (ZFGE .LT. 1.D-3 ) THEN
                  WRITE(*,*)' Problem for GP STNID FGE= ',  &
                     obs_elem_c(lobsSpaceData,'STID',index_header),ZFGE
                  ZFGE=1.D-3
                ENDIF
              ELSE IF (CDFAM.EQ.'CH') THEN
                IF ( ZFGE**2 + ZOER**2 .LT. 1.D-60)THEN
                  WRITE(*,*) ' Problem for STNID FGE ZOER=',  &
                             obs_elem_c(lobsSpaceData,'STID',index_header),ZFGE,ZOER
                  ZFGE=1.D-30
                  ZOER=1.D-30
                ENDIF
              ELSE
                IF ( ZFGE**2 + ZOER**2 .LT. 1.D-5)THEN
                  WRITE(*,*) ' Problem for STNID FGE ZOER=',  &
                             obs_elem_c(lobsSpaceData,'STID',index_header),ZFGE,ZOER
                  ZFGE=1.D-5
                  ZOER=1.D-5
                ENDIF
              ENDIF
              IF ( CDFAM .EQ. 'GP' .AND. ITYP .EQ. BUFR_NEZD ) THEN
                ZBGCHK=(ZOMP**2)/(ZFGE**2)
                ZSOP = ZFGE
              ELSE
                ZBGCHK=(ZOMP)**2/(ZFGE**2 + ZOER**2)
                ZSOP = SQRT(ZFGE**2 + ZOER**2)
              ENDIF
!C
!C            UPDATE QUALITY CONTROL FLAGS
!C            ( ELEMENT FLAGS + GLOBAL HEADER FLAGS)
!C
              INAM =obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body)
              if( inam.eq.12192 .and. obs_bodyElem_i(lobsSpaceData,OBS_XTR,index_body).eq.0)then
                zsum=zsum+zfge*zfge
                zsumo=zsumo+zoer*zoer
                icoun=icoun+1
              endif
              IDBURP=obs_headElem_i(lobsSpaceData,OBS_ITY,INDEX_HEADER)
              IFLAG=ISETFLAG(CDFAM,IDBURP,INAM,ZLEV,ZBGCHK,LMODIF1020)
!C
!C            CONVERT ZTD VALUES FROM M TO MM FOR PRINTOUT
!C
              LLZD = .FALSE.
              IF ( CDFAM .EQ. 'GP' .AND. ITYP .EQ. BUFR_NEZD) THEN
                ZVAR = ZVAR * 1000.0D0
                ZOER = ZOER * 1000.0D0
                ZFGE = ZFGE * 1000.0D0
                ZOMP = ZOMP * 1000.0D0
                ZSOP = ZSOP * 1000.0D0
                LLZD = .TRUE.
              ENDIF

              IF (IFLAG .GE. 2 .OR. (LLZD .AND. LTESTOP)) THEN
                IF (CDFAM.NE.'CH') THEN 
                  WRITE(*,122)  &
                     obs_elem_c(lobsSpaceData,'STID',index_header),  &
                     zlat,zlon,IDBURP,INAM,ZLEV,ZVAR,ZOER  &
                     ,ZFGE,ZOMP,ZSOP,ZBGCHK,IFLAG
                ELSE 
                  WRITE(*,124)  &
                     obs_elem_c(lobsSpaceData,'STID',index_header),  &
                     zlat,zlon,IDBURP,INAM,ZLEV,ZVAR,ZOER  &
                     ,ZFGE,ZOMP,ZSOP,ZBGCHK,IFLAG
                END IF
                IF (IFLAG .GE. 2) INREJ = INREJ + 1
                IF (IFLAG .GE. 2 .AND. CDFAM .EQ. 'GP') THEN
                  IF (ITYP .EQ. BUFR_NEZD) INZREJ = INZREJ+1
                  IF (ITYP .EQ. BUFR_NEPS) INPREJ = INPREJ+1
                  IF (ITYP .EQ. BUFR_NETS) INTREJ = INTREJ+1
                  IF (ITYP .EQ. BUFR_NESS) INDREJ = INDREJ+1
                ENDIF
              ENDIF

              IF ( IFLAG .EQ. 1 ) THEN
                call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body),13))
              ELSEIF ( IFLAG .EQ. 2 ) THEN
                call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body),14))
                call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body),16))
                call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body),09))
                call obs_headSet_i(lobsSpaceData,OBS_ST1,index_header,ibset(obs_headElem_i(lobsSpaceData,OBS_ST1,index_header),06))

              ELSEIF ( IFLAG .EQ. 3 ) THEN
                call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body),15))
                call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body),16))
                call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body),09))
                call obs_headSet_i(lobsSpaceData,OBS_ST1,index_header,ibset(obs_headElem_i(lobsSpaceData,OBS_ST1,index_header),06))
              ENDIF
            ENDIF
          ENDDO BODY
124   FORMAT(2x,a9,1x,f6.2,1x,f7.2,1x,I3,1x,I5,7(2x,G11.2),I3 )
122   FORMAT(2x,a9,1x,f6.2,1x,f7.2,1x,I3,1x,I5,7(2x,F11.2),I3 )

        ENDDO HEADER

      ELSE !IF (.not. new_bgck_sw .or. CDFAM .ne. 'SW') THEN

        call obs_set_current_body_list(lobsSpaceData, 'SW')
        bodyuv: do
          index_body = obs_getBodyIndex(lobsSpaceData)
          if (index_body < 0) exit bodyuv

          ! Only process pressure level observations flagged to be assimilated
          i_ass = obs_bodyElem_i (lobsSpaceData,OBS_ASS,index_body)
          i_vco = obs_bodyElem_i (lobsSpaceData,OBS_VCO,index_body)

          if(i_ass.ne.1 .or. i_vco.ne.2) cycle bodyuv
 
          index_header     = obs_bodyElem_i(lobsSpaceData,OBS_HIND,index_body)
          index_body_start = obs_headElem_i(lobsSpaceData,OBS_RLN,index_header)

          i_vnm = obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body)
          zlev  = obs_bodyElem_r(lobsSpaceData,OBS_PPP,index_body)

          found_u = .false.
          found_v = .false.

          if (i_vnm .EQ. BUFR_NEUU) then
      
            index_body_u = index_body
            found_u = .true.

            do i_oth =index_body_start,index_body
              istyp = obs_bodyElem_i(lobsSpaceData,OBS_VNM,i_oth)
              zslev = obs_bodyElem_r(lobsSpaceData,OBS_PPP,i_oth)
              if ( istyp .eq. BUFR_NEVV .and. zslev .eq. zlev ) then
                index_body_v = i_oth
                found_v = .true.
              end if
            end do

          else ! i_vnm

            index_body_v = index_body
            found_v = .true.
 
            do i_oth =index_body_start,index_body
              istyp = obs_bodyElem_i(lobsSpaceData,OBS_VNM,i_oth)
              zslev = obs_bodyElem_r(lobsSpaceData,OBS_PPP,i_oth)
              IF ( istyp .eq. BUFR_NEUU .and. zslev .eq. zlev ) then
                index_body_u = i_oth
                found_u = .true.
              end if
            end do

          end if !i_vnm

          if(found_u .and. found_v) then

            uu_d = obs_bodyElem_r(lobsSpaceData,OBS_OMP,index_body_u)
            vv_d = obs_bodyElem_r(lobsSpaceData,OBS_OMP,index_body_v)
            uu_r = obs_bodyElem_r(lobsSpaceData,OBS_OER,index_body_u)
            vv_r = obs_bodyElem_r(lobsSpaceData,OBS_OER,index_body_v)
            uu_f = obs_bodyElem_r(lobsSpaceData,OBS_HPHT,index_body_u)
            vv_f = obs_bodyElem_r(lobsSpaceData,OBS_HPHT,index_body_v)

            duv2 = uu_d**2 + vv_d**2

            iflag = 0
            duv2_lim = (uu_r**2 + uu_f**2 + vv_r**2 + vv_f**2)*1
            if(duv2 > duv2_lim) then
              iflag = 1
            end if
            duv2_lim = (uu_r**2 + uu_f**2 + vv_r**2 + vv_f**2)*2
            if(duv2 > duv2_lim) then
              iflag = 2
            end if
            duv2_lim = (uu_r**2 + uu_f**2 + vv_r**2 + vv_f**2)*4
            if(duv2 > duv2_lim) then
              iflag = 3
            end if

            inobs = inobs + 2
            if (iflag .ge. 2) inrej = inrej + 2

            if ( iflag .eq. 1 ) then
              call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body_u,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body_u),13))
              call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body_v,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body_v),13))

            else if ( iflag .eq. 2 ) then
              call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body_u,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body_u),14))
              call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body_u,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body_u),16))
              call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body_u,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body_u),09))
              call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body_v,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body_v),14))
              call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body_v,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body_v),16))
              call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body_v,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body_v),09))
              call obs_headSet_i(lobsSpaceData,OBS_ST1,index_header,ibset(obs_headElem_i(lobsSpaceData,OBS_ST1,index_header),06))

            else if ( iflag .eq. 3 ) then
              call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body_u,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body_u),15))
              call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body_u,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body_u),16))
              call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body_u,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body_u),09))
              call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body_v,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body_v),15))
              call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body_v,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body_v),16))
              call obs_bodySet_i(lobsSpaceData,OBS_FLG,index_body_v,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,index_body_v),09))
              call obs_headSet_i(lobsSpaceData,OBS_ST1,index_header,ibset(obs_headElem_i(lobsSpaceData,OBS_ST1,index_header),06))
            end if

          end if !if(found_u .and. found_v) then

        end do bodyuv

      end if !IF (.not. new_bgck_sw .or. CDFAM .ne. 'SW') THEN

      IF ( INOBS .GT. 0 ) THEN
        WRITE(*,*)' '
        WRITE(*,*) '  BGCDATA: FINISHED BGCHECK OF ',CDFAM, ' DATA'
        WRITE(*,123) 'BGCDATA:   ',INREJ, ' OBSERVATIONS REJECTED OUT OF ', INOBS
        WRITE(*,*)' '
      END IF

      IF ( (INOBS .GT. 0) .AND. (CDFAM .EQ. 'GP') ) THEN
        WRITE(*,*)' '
        WRITE(*,*) '  BGCDATA:    REPORT FOR GP FAMILY OF OBSERVATIONS'
        WRITE(*,123) 'BGCDATA:   ',INZREJ, ' ZTD  OBSERVATIONS REJECTED OUT OF ', INZOBS
        WRITE(*,123) 'BGCDATA:   ',INPREJ, ' PSFC OBSERVATIONS REJECTED OUT OF ', INPOBS
        WRITE(*,123) 'BGCDATA:   ',INTREJ, ' TSFC OBSERVATIONS REJECTED OUT OF ', INTOBS
        WRITE(*,123) 'BGCDATA:   ',INDREJ, ' DPDS OBSERVATIONS REJECTED OUT OF ', INDOBS
        WRITE(*,*)' '
      END IF

123   FORMAT(2X,A,I0,A,I0)

      WRITE(*,*)' '
      WRITE(*,*)' ---------------------------'
      WRITE(*,*)'           DONE BGCDATA     '
      WRITE(*,*)' ---------------------------'
      WRITE(*,*)' '
      if ( icoun .gt. 0) then
         write(*,*) ' icoun meanfge=',icoun,zsum/icoun
         write(*,*) ' icoun meanzoer',icoun,zsumo/icoun
      endif

      END SUBROUTINE BGCDATA

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
!! MODULE backgroundCheck_mod (prefix='bgck' category='1. High-level functionality')
!!
!! *Purpose*: Contains subroutines related to the background check 
!!
!--------------------------------------------------------------------------
module backgroundCheck_mod
  use mathPhysConstants_mod
  use bufr_mod
  use obsSpaceData_mod
  use gps_mod
  use utilities_mod
  use mpi_mod
  use mpivar_mod
  use columnData_mod
  use obsSpaceDiag_mod
  use earthConstants_mod
  use verticalCoord_mod
  use computeHBHT_mod
  implicit none
  private

  ! public procedures
  public :: bgck_bgcheck_conv

  contains

!--------------------------------------------------------------------------
!! *Purpose*: Do background check on all conventionnal observations
!!
!! @author P. Koclas *CMC/CMDA  Nov 1998
!!
!! Revision:  M. Sitwell (ARQI/AQRD) May 2015, March 2016
!!           - Added call to BGCDATA for chemical constituents
!!           - Added loop over bgfam list
!!
!!           Y. Rochon (ARQI/AQRD) June 2016
!!           - Added call to osd_ObsSpaceDiag
!!
!!           S. Laroche (ARMA/MRD) October 2017
!!           - New option NEW_BGCK_SW for AMVs
!!
!--------------------------------------------------------------------------
subroutine bgck_bgcheck_conv( columng, columnhr, obsSpaceData )

  IMPLICIT NONE

  type(struct_obs)        :: obsSpaceData  ! Observation-related data
  type(struct_columnData) :: columng       !
  type(struct_columnData) :: columnhr      ! 

  integer :: j, jdata
  real(8) :: zjo

  integer            :: nulNam, ier, fnom, fclos
  character(len=256) :: namFile
  logical            :: NEW_BGCK_SW

  character(len=2), dimension(12) :: bgfam = (/ 'UA', 'AI', 'HU', 'SF', 'ST', 'SW', 'SC', 'PR', 'GP', 'CH', 'TM', 'AL' /)
      
  call tmg_start(3,'BGCHECK_CONV')

  write(*,FMT=9000)
9000 FORMAT(//,3(" **********"),/," BEGIN CONVENTIONNAL BACKGROUND CHECK",/,3(" **********"),/)

  NEW_BGCK_SW = .false.

  NAMELIST /NAMBGCKCONV/NEW_BGCK_SW
  namFile=trim("flnml")
  nulNam=0
  ier = FNOM( NULNAM, NAMFILE, 'R/O', 0 )

  read( nulNam, nml = NAMBGCKCONV, IOSTAT = ier )
  if ( ier /= 0 ) then
    write(*,*) 'bgcheck_conv: No valid namelist NAMBGCKCONV found'
  end if

  ier = fclos(nulNam)

  write(*,*) 'new_bgck_sw = ',new_bgck_sw


!     CALCULATE HBHT (sigma_B in observation space)
!     ----------------------------------------------
!
  call hbht_compute( columng,   &   ! IN
                     columnhr,  &   ! IN
                     obsSpaceData ) ! INOUT

!
!     DO A BACKGROUND CHECK ON ALL THE OBSERVATIONS
!     ----------------------------------------------

  do j = 1, size( bgfam )
    ! For SW only, old and new background check schemes controlled by "new_bgck_sw"
    if ( obs_famExist( obsSpaceData, bgfam(j) )) CALL bgck_data( ZJO, bgfam(j), obsSpaceData, new_bgck_sw )
  end do

  if (obs_famExist(obsSpaceData,'RO')) CALL bgck_gpsro( columnhr , obsSpaceData )

! Conduct obs-space post-processing diagnostic tasks (some diagnostic 
! computations controlled by NAMOSD namelist in flnml)

  call osd_ObsSpaceDiag( obsSpaceData, columng, analysisMode_opt = .false. )


  call tmg_stop(3)

end subroutine bgck_bgcheck_conv


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
      SUBROUTINE bgck_data(PJO,CDFAM,lobsSpaceData,new_bgck_sw)

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
      INTEGER IFLAG,INAM,INDEX_HEADER,IDBURP
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
            LLOK=( obs_bodyElem_i(lobsSpaceData,OBS_ASS,index_body) == obs_assimilated)
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

              ! Protect against zfge values that are too small
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

              ! Calculate zbgchk from zfge
              IF ( CDFAM .EQ. 'GP' .AND. ITYP .EQ. BUFR_NEZD ) THEN
                ZBGCHK=(ZOMP**2)/(ZFGE**2)
                ZSOP = ZFGE
              ELSE
                ZBGCHK=(ZOMP)**2/(ZFGE**2 + ZOER**2)
                ZSOP = SQRT(ZFGE**2 + ZOER**2)
              ENDIF

              ! UPDATE QUALITY CONTROL FLAGS, based on zbgchk
              ! ( ELEMENT FLAGS + GLOBAL HEADER FLAGS)
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

          if(i_ass /= obs_assimilated .or. i_vco.ne.2) cycle bodyuv
 
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

      END SUBROUTINE bgck_data


!--------------------------------------------------------------------------
!!  *Purpose*: Set background check flag on GPSRO data if ABS(O-P)/P is too large
!!
!!  @author P. KOCLAS. Mar 2008
!!  Modified: 
!!            -J.M. Aparicio, Dec 2012.
!!                -Simplified and adapted to both refractivity and bending angle data
!!
!--------------------------------------------------------------------------
      SUBROUTINE bgck_gpsro(lcolumnhr,lobsSpaceData)
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
               IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY) == obs_assimilated ) THEN
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
      
      END SUBROUTINE bgck_gpsro


!--------------------------------------------------------------------------
!! *Purpose*:  Set BACKGROUND CHECK FLAGS According to values set in a table.
!!             Original values in table come from ecmwf.
!!
!! @author P. Koclas *CMC/CMSV  September 1998
!!
!! Arguments:
!!     -KVNAM= VARIABLE NAME ( BURP )
!!     -KODTYP=BURP CODE TYPE
!!     -CDFAM= FAMILY  NAME ( 'UA' , 'AI'   ...etc.. )
!!     -ZLEV = LEVEL
!!     -zbgchk=NORMALIZED BACKGROUND DEPARTURE
!!     -lmodif1020=switch to activate special criteria for backound check (*ua 10-20 mb)
!!
!--------------------------------------------------------------------------
      function isetflag(cdfam,kodtyp,kvnam,zlev,zbgchk,lmodif1020)

      IMPLICIT NONE
      integer isetflag
      integer kvnam,kodtyp
      real*8 zlev
      real*8 zbgchk
      character*2 cdfam
      logical lmodif1020

      real*8 zgzcrit(3),zttcrit(3),zuvcrit(3),zescrit(3),zdzcrit(3),zalcrit(3)
      real*8 zpscrit(3),zpncrit(3),ztscrit(3),zswcrit(3),zzdcrit(3),zviscrit(3)
      real*8 zchcrit(3)

      isetflag=0

! ASSUME CVCORD = GEMHYB
         zttcrit(1) = 9.00D0
         zttcrit(2) = 16.00D0
         zttcrit(3) = 25.00D0

         zuvcrit(1) = 10.00D0
         zuvcrit(2) = 20.00D0
         zuvcrit(3) = 30.00D0
 
         zalcrit(1) = 10.00D0
         zalcrit(2) = 20.00D0
         zalcrit(3) = 30.00D0

         zescrit(1) = 10.00D0
         zescrit(2) = 20.00D0
         zescrit(3) = 30.00D0

         zpscrit(1) = 9.00D0
         zpscrit(2) = 16.00D0
         zpscrit(3) = 25.00D0

         zpncrit(1) = 10.00D0
         zpncrit(2) = 20.00D0
         zpncrit(3) = 30.00D0

         zswcrit(1) = 10.00D0
         zswcrit(2) = 20.00D0
         zswcrit(3) = 30.00D0

         ztscrit(1) = 5.00D0
         ztscrit(2) = 25.00D0
         ztscrit(3) = 30.00D0

         zdzcrit(1) = 2.25D0
         zdzcrit(2) = 5.06D0
         zdzcrit(3) = 7.56D0

         zgzcrit(1) = 12.25D0
         zgzcrit(2) = 25.00D0
         zgzcrit(3) = 36.00D0

         zzdcrit(1) = 9.00D0
         zzdcrit(2) = 16.00D0
         zzdcrit(3) = 25.00D0

         zchcrit(1) = 9.00D0
         zchcrit(2) = 16.00D0
         zchcrit(3) = 25.00D0

         zviscrit(1) = 10.00D0
         zviscrit(2) = 20.00D0
         zviscrit(3) = 30.00D0

         if ( kodtyp .eq. 37 ) then
           zuvcrit(2)=25.D0
         else
           zuvcrit(2)=20.D0
         endif
!C
!C     SET FLAG FOR HEIGHTS
!C
      if ( kvnam .eq. BUFR_NEGZ ) then
         if (      zbgchk .gt. zgzcrit(1) .and. zbgchk .lt. zgzcrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zgzcrit(2) .and. zbgchk .lt. zgzcrit(3) ) then
            isetflag=2
         else if ( zbgchk .ge. zgzcrit(3) )then
           isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR TEMPERATURE
!C
      if ( kvnam .eq. BUFR_NETT ) then
         if (      zbgchk .gt. zttcrit(1) .and. zbgchk .lt. zttcrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zttcrit(2) .and. zbgchk .lt. zttcrit(3) ) then
            isetflag=2
         else if ( zbgchk .ge. zttcrit(3) )then
           isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR SATEMS
!C
      if ( kvnam .eq. BUFR_NEDZ ) then
         if (      zbgchk .gt. zdzcrit(1) .and. zbgchk .lt. zdzcrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zdzcrit(2) .and. zbgchk .lt. zdzcrit(3) ) then
            isetflag=2
         else if ( zbgchk .ge. zdzcrit(3) )then
           isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR WIND COMPONENTS
!C
      if ( kvnam .eq. BUFR_NEUU .or. kvnam .eq. BUFR_NEVV ) then
         if (      zbgchk .gt. zuvcrit(1) .and. zbgchk .lt. zuvcrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zuvcrit(2) .and. zbgchk .lt. zuvcrit(3) ) then
            isetflag=2
         else if ( zbgchk .ge. zuvcrit(3) )then
           isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR ALADIN HLOS WIND OBSERVATIONS
!C
      if ( kvnam .eq. BUFR_NEAL ) then
         if      ( zbgchk >= zalcrit(1) .and. zbgchk < zalcrit(2) ) then
            isetflag=1
         else if ( zbgchk >= zalcrit(2) .and. zbgchk < zalcrit(3) ) then
            isetflag=2
         else if ( zbgchk >= zalcrit(3) )then
            isetflag=3
         end if
      end if
!C
!C     SET FLAG FOR SURFACE WIND COMPONENTS
!C
      if ( kvnam .eq. BUFR_NEUS .or. kvnam .eq. BUFR_NEVS ) then
         if (      zbgchk .gt. zswcrit(1) .and. zbgchk .lt. zswcrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zswcrit(2) .and. zbgchk .lt. zswcrit(3) ) then
            isetflag=2
         else if ( zbgchk .ge. zswcrit(3) )then
           isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR SURFACE WIND GUST
!C
      if ( kvnam == bufr_gust ) then
         if (      zbgchk >= zswcrit(1) .and. zbgchk < zswcrit(2) ) then
           isetflag=1
         else if ( zbgchk >= zswcrit(2) .and. zbgchk < zswcrit(3) ) then
            isetflag=2
         else if ( zbgchk >= zswcrit(3) )then
           isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR DEW POINT DEPRESSION
!C
      if ( kvnam .eq. BUFR_NEES ) then
         if (      zbgchk .gt. zescrit(1) .and. zbgchk .lt. zescrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zescrit(2) .and. zbgchk .lt. zescrit(3) ) then
           isetflag=2
         else if ( zbgchk .ge. zescrit(3) )then
           isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR SURFACE PRESSURE
!C
      if ( kvnam .eq. BUFR_NEPS ) then
         if (      zbgchk .gt. zpscrit(1) .and. zbgchk .lt. zpscrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zpscrit(2) .and. zbgchk .lt. zpscrit(3) ) then
           isetflag=2
         else if ( zbgchk .ge. zpscrit(3) )then
           isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR MEAN SEA LEVEL PRESSURE
!C
      if ( kvnam .eq. BUFR_NEPN ) then
         if (      zbgchk .gt. zpncrit(1) .and. zbgchk .lt. zpncrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zpncrit(2) .and. zbgchk .lt. zpncrit(3) ) then
           isetflag=2
         else if ( zbgchk .ge. zpncrit(3) )then
           isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR SURFACE TEMPERATURE
!C
      if ( kvnam .eq. BUFR_NETS ) then
         if (      zbgchk .gt. ztscrit(1) .and. zbgchk .lt. ztscrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. ztscrit(2) .and. zbgchk .lt. ztscrit(3) ) then
           isetflag=2
         else if ( zbgchk .ge. ztscrit(3) )then
           isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR VISIBILITY
!C
      if ( kvnam .eq. bufr_vis ) then
         if (      zbgchk .gt. zviscrit(1) .and. zbgchk .lt. zviscrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zviscrit(2) .and. zbgchk .lt. zviscrit(3) ) then
           isetflag=2
         else if ( zbgchk .ge. zviscrit(3) )then
           isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR GB-GPS ZENITH DELAY
!C
      if ( kvnam .eq. BUFR_NEZD ) then
         if (      zbgchk .gt. zzdcrit(1) .and. zbgchk .lt. zzdcrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zzdcrit(2) .and. zbgchk .lt. zzdcrit(3) ) then
           isetflag=2
         else if ( zbgchk .ge. zzdcrit(3) )then
           isetflag =3
         endif
      endif
!C
!C     SET FLAG FOR CHEMICAL CONSTITUENTS
!C
      if ( cdfam .eq. 'CH' ) then
         if (      zbgchk .gt. zchcrit(1) .and. zbgchk .lt. zchcrit(2) ) then
           isetflag=1
         else if ( zbgchk .gt. zchcrit(2) .and. zbgchk .lt. zchcrit(3) ) then
           isetflag=2
         else if ( zbgchk .ge. zchcrit(3) )then
           isetflag =3
         endif
      endif

      return
      end function isetflag


end module backgroundCheck_mod

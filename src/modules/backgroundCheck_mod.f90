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

module backgroundCheck_mod
  ! MODULE backgroundCheck_mod (prefix='bgck' category='1. High-level functionality')
  !
  ! :Purpose: Performs the background check on all conventional observations
  !
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
  use horizontalCoord_mod
  use obsSpaceErrorStdDev_mod
  use obsFamilyList_mod
  
  implicit none
  private

  ! public procedures
  public :: bgck_bgcheck_conv

  integer,          parameter :: numFamilyToProcess = 11
  character(len=2), parameter :: familyListToProcess(numFamilyToProcess)= (/ &
                                 'UA','AI','SF','SC','SW','PR','RO','GP','RA', &
                                 'CH','AL'/)

  contains

  !--------------------------------------------------------------------------
  ! bgck_bgcheck_conv
  !--------------------------------------------------------------------------
  subroutine bgck_bgcheck_conv( columnTrlOnAnlIncLev, columnTrlOnTrlLev, hco_anl, obsSpaceData )
     !
     !:Purpose: Do background check on all conventional observations
     !
     implicit none

     ! Arguments:
     type(struct_obs)          :: obsSpaceData  ! Observation-related data
     type(struct_columnData)   :: columnTrlOnAnlIncLev       !
     type(struct_columnData)   :: columnTrlOnTrlLev      ! 
     type(struct_hco), pointer :: hco_anl

     ! Locals:
     integer :: familyIndex
     real(8) :: zjo
     integer            :: nulNam, ier, fnom, fclos
     character(len=256) :: namFile
     logical            :: NEW_BGCK_SW, noObsToProcess

     write(*,FMT=9000)
9000 FORMAT(//,3(" **********"),/," BEGIN CONVENTIONNAL BACKGROUND CHECK",/,3(" **********"),/)

     ! Check if any observations are present for conventional background check
     noObsToProcess = .true.
     do familyIndex = 1, numFamilyToProcess
       if (obs_famExist(obsSpaceData,familyListToProcess(familyIndex))) then
         noObsToProcess = .false.
       end if
     end do
     if (noObsToProcess) then
       write(*,*) 'bgcheck_conv: No observations to process'
       return
     end if
     
     call tmg_start(3,'BGCHECK_CONV')

     NEW_BGCK_SW = .false.

     NAMELIST /NAMBGCKCONV/NEW_BGCK_SW
     namFile=trim("flnml")
     nulNam=0
     ier = FNOM( NULNAM, NAMFILE, 'R/O', 0 )
     read( nulNam, nml = NAMBGCKCONV, IOSTAT = ier )
     if ( ier /= 0 ) write(*,*) 'bgcheck_conv: No valid namelist NAMBGCKCONV found'
     ier = fclos(nulNam)
     
     write(*,*) 'new_bgck_sw = ',new_bgck_sw
 
     ! Obtain or calc OmP-error std dev when requested and possible.
     ! Otherwise calc HBHT contribution (sigma_B in observation space)  
     ! -------------------------------------------------------------------- 

     call ose_computeStddev( columnTrlOnAnlIncLev, hco_anl, & ! IN
                             obsSpaceData )      ! INOUT

     ! DO A BACKGROUND CHECK ON ALL THE OBSERVATIONS
     ! ----------------------------------------------

     do familyIndex = 1, ofl_numFamily
       ! For SW only, old and new background check schemes controlled by "new_bgck_sw"
       if ( obs_famExist(obsSpaceData,ofl_familyList(familyIndex)) ) then
         call bgck_data(zjo, ofl_familyList(familyIndex), obsSpaceData, new_bgck_sw)
       end if
     end do

     if (obs_famExist(obsSpaceData,'RO')) CALL bgck_gpsro( columnTrlOnTrlLev , obsSpaceData )

     ! Conduct obs-space post-processing diagnostic tasks (some diagnostic 
     ! computations controlled by NAMOSD namelist in flnml)

     call osd_ObsSpaceDiag(obsSpaceData, columnTrlOnAnlIncLev, hco_anl, analysisMode_opt = .false.)

     call tmg_stop(3)

  end subroutine bgck_bgcheck_conv

  !--------------------------------------------------------------------------
  ! bgck_data
  !--------------------------------------------------------------------------
  subroutine bgck_data(PJO,CDFAM,lobsSpaceData,new_bgck_sw)
      !
      !:Purpose:  Calculate a background check for a data family and set the
      !           appropriate quality-control flags in obsSpaceData
      !

      IMPLICIT NONE
     
      ! NOTE 1: YSFERRWGT IN MODGPSZTD_MOD (FROM NML FILE) IS USED HERE FOR ERROR WEIGHTING
      !         OF TIME SERIES (FGAT) GPS MET OBSERVATIONS PS, TS, DPDS. IT IS APPLIED
      !         (ALONG WITH YZDERRWGT FOR ZTD) IN S/R SETERR AS A MULT. FACTOR TO ERRORS.
      !         YZDERRWGT AND YSFERRWGT = 1 FOR NORMAL 3D-VAR (3D-THINNING).

      type(struct_obs) :: lobsSpaceData
      real(8) :: PJO,zsum,zsumo
      CHARACTER(len=2) :: CDFAM
      LOGICAL :: new_bgck_sw
      INTEGER :: IFLAG,INAM,INDEX_HEADER,IDBURP
      INTEGER :: ITYP
      INTEGER :: INDEX_BODY,icoun
      INTEGER :: INOBS, INREJ, INZOBS, INZREJ
      INTEGER :: INPOBS, INTOBS, INDOBS, INPREJ, INTREJ, INDREJ
      real(8) :: ZOER,ZOMP,ZFGE,ZOMPER,ZBGCHK,ZVAR,ZLEV,ZLAT,ZLON,ZSOP
      LOGICAL :: LLOK, LLZD, LMODIF1020
      character(len=12) :: stnid

      INTEGER :: i_ass,i_vco,i_vnm,i_oth,istyp,index_body_u,index_body_v,index_body_start
      real(8) :: uu_d,uu_r,uu_f,vv_d,vv_r,vv_f,duv2,duv2_lim,zslev
      LOGICAL :: found_u,found_v

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

      ! Initialize counters for GP family observations

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

          stnid = obs_elem_c(lobsSpaceData,'STID',index_header)

          ! 1. Computation of (HX - Z)**2/(SIGMAo**2 +SIGMAp**2)
          ! ----------------------------------------------------
          
          ! loop over all body indices for this index_header
          call obs_set_current_body_list(lobsSpaceData, index_header)
          BODY: do 
            index_body = obs_getBodyIndex(lobsSpaceData)
            if (index_body < 0) exit BODY

            ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body)
            LLOK=( obs_bodyElem_i(lobsSpaceData,OBS_ASS,index_body) == obs_assimilated)
            IF ( LLOK ) THEN
              INOBS = INOBS + 1
              IF (CDFAM == 'GP') THEN
                IF (ITYP == BUFR_NEZD) INZOBS = INZOBS+1
                IF (ITYP == BUFR_NEPS) INPOBS = INPOBS+1
                IF (ITYP == BUFR_NETS) INTOBS = INTOBS+1
                IF (ITYP == BUFR_NESS) INDOBS = INDOBS+1
              ENDIF
              ZVAR = obs_bodyElem_r(lobsSpaceData,OBS_VAR,index_body)
              ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,index_body)              
              ZLAT  = obs_headElem_r(lobsSpaceData,OBS_LAT,index_header) &
                                                  * MPC_DEGREES_PER_RADIAN_R8
              ZLON  = obs_headElem_r(lobsSpaceData,OBS_LON,index_header) &
                                                  * MPC_DEGREES_PER_RADIAN_R8

              ! BACKGROUND CHECK

              ZOMP  = obs_bodyElem_r(lobsSpaceData,OBS_OMP,index_body)

              ! Get error std dev
              
              ! std(O-P) (valid/available if > 0.0)
              zomper = obs_bodyElem_r(lobsSpaceData,OBS_OMPE,index_body)
                                         
              ! std(O)
              zoer = obs_bodyElem_r(lobsSpaceData,OBS_OER,index_body)
              ! std(P)
              zfge  = obs_bodyElem_r(lobsSpaceData,OBS_HPHT,index_body)
              ! NOTE: For GB-GPS ZTD observations (GP family), ZFGE is not the FGE but
              !       rather Std(O-P) set in s/r SETERRGPSGB.
              if ( cdfam == 'GP' .and. ityp == BUFR_NEZD ) then
                if ( zomper <= 0.0d0 .and. zfge > 0.0d0 ) zomper = zfge
              end if
              
              ! Protect against error std dev values that are too small
              IF ( CDFAM == 'GP' ) THEN
                IF (ITYP == BUFR_NEZD) THEN
                  ZOER = ZOER/YZDERRWGT
                ELSE
                  ZOER = ZOER/YSFERRWGT
                ENDIF
                IF ( ZOER < 1.0d-3 .AND. ITYP /= BUFR_NEZD ) THEN
                  WRITE(*,*)' Problem for GP STNID ZOER= ' ,  &
                    stnid, ZOER
                  CALL utl_abort('BGCDATA: PROBLEM WITH OER.')
                ENDIF
                IF ( ZFGE < 1.0d-3 ) THEN
                  WRITE(*,*)' Problem for GP STNID FGE= ',  &
                     stnid, ZFGE
                  ZFGE=1.0d-3
                ENDIF
                if ( zomper > 0.0d0 .and. zomper < 1.0d-3 ) zomper = 1.0d-3
              ELSE IF ( CDFAM == 'CH' ) THEN
                IF ( ZFGE**2 + ZOER**2 < 1.0d-60)THEN
                  WRITE(*,*) ' Problem for STNID FGE ZOER=',  &
                           stnid, ZFGE, ZOER
                  ZFGE=1.0d30
                  ZOER=1.0d-30
                ENDIF
                if ( zomper > 0.0 .and. zomper < 1.0d-30 ) zomper = 1.0d-30
              ELSE
                if ( zfge < 0.0d0 ) zfge = 0.0d0
                IF ( ZFGE**2 + ZOER**2 < 1.0d-5)THEN
                  WRITE(*,*) ' Problem for STNID FGE ZOER=',  &
                             stnid, ZFGE, ZOER
                  ZFGE=1.0d-5
                  ZOER=1.0d-5
                ENDIF
                if ( zomper > 0.0 .and. zomper < 1.0D-5 ) zomper = 1.0d-5
              ENDIF
              
              ! Calculate zbgchk
              
              if ( zomper > 0.0d0 ) then
                zbgchk = (zomp**2)/(zomper**2) 
                zsop = zomper
              else
                ZBGCHK=(ZOMP)**2/(ZFGE**2 + ZOER**2)
                ZSOP = SQRT(ZFGE**2 + ZOER**2)
              end if
              
              ! UPDATE QUALITY CONTROL FLAGS, based on zbgchk
              ! ( ELEMENT FLAGS + GLOBAL HEADER FLAGS)
              INAM =obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body)
              if( inam.eq.12192 .and. obs_bodyElem_i(lobsSpaceData,OBS_XTR,index_body).eq.0)then
                zsum=zsum+zfge*zfge
                zsumo=zsumo+zoer*zoer
                icoun=icoun+1
              endif
              IDBURP=obs_headElem_i(lobsSpaceData,OBS_ITY,INDEX_HEADER)
              IFLAG=ISETFLAG(CDFAM,IDBURP,INAM,ZBGCHK)

              ! CONVERT ZTD VALUES FROM M TO MM FOR PRINTOUT

              LLZD = .FALSE.
              IF ( CDFAM .EQ. 'GP' .AND. ITYP .EQ. BUFR_NEZD) THEN
                ZVAR = ZVAR * 1000.0D0
                ZOER = ZOER * 1000.0D0
                ZFGE = ZFGE * 1000.0D0
                ZOMP = ZOMP * 1000.0D0
                ZSOP = ZSOP * 1000.0D0
                LLZD = .TRUE.
              ENDIF

              stnid = obs_elem_c(lobsSpaceData,'STID',index_header)
              IF (IFLAG .GE. 2 .OR. (LLZD .AND. LTESTOP)) THEN
                IF (CDFAM.NE.'CH') THEN 
                  WRITE(*,122)  &
                     stnid,zlat,zlon,IDBURP,INAM,ZLEV,ZVAR,ZOER  &
                     ,ZFGE,ZOMP,ZSOP,ZBGCHK,IFLAG
                ELSE 
                  WRITE(*,124)  &
                     stnid,zlat,zlon,IDBURP,INAM,ZLEV,ZVAR,ZOER  &
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

124       FORMAT(2x,a9,1x,f6.2,1x,f7.2,1x,I3,1x,I5,7(2x,G11.2),I3 )
122       FORMAT(2x,a9,1x,f6.2,1x,f7.2,1x,I3,1x,I5,7(2x,F11.2),I3 )

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

  end subroutine bgck_data

  !--------------------------------------------------------------------------
  ! bgck_gpsro
  !--------------------------------------------------------------------------
  subroutine bgck_gpsro(columnTrlOnTrlLev,lobsSpaceData)
      !
      !:Purpose: Set background-check flag on GPSRO data if ABS(O-P)/P is too
      !          large
      !
      IMPLICIT NONE

      type(struct_columnData) :: columnTrlOnTrlLev
      type(struct_obs) :: lobsSpaceData
      type(struct_vco), pointer :: vco_trl
      real(8) :: HNH1, ZOBS, ZMHX, ZOMF, ZREF, ZOER, Rad

      INTEGER :: INDEX_HEADER
      INTEGER :: IDATYP
      INTEGER :: IDATA   , IDATEND, INDEX_BODY
      INTEGER :: NGPSLEV
      INTEGER :: stat, iversion

      LOGICAL :: LSTAG

      LSTAG = .FALSE.
      vco_trl => col_getVco(columnTrlOnTrlLev)
      stat = vgd_get(vco_trl%vgrid,key='ig_1 - vertical coord code',value=iversion)
      if (iversion .eq. 5002) LSTAG = .TRUE. 
      
      WRITE(*,*)'ENTER BGCSGPSRO'
      
      ! 1.  Initializations
      !     ---------------

      NGPSLEV=col_getNumLev(columnTrlOnTrlLev,'TH')

      ! Loop over all files

      ! loop over all header indices of the 'RO' family
      call obs_set_current_header_list(lobsSpaceData,'RO')
      HEADER: do
         index_header = obs_getHeaderIndex(lobsSpaceData)
         if (index_header < 0) exit HEADER
     
         ! Process only refractivity data (codtyp 169)

         IDATYP = obs_headElem_i(lobsSpaceData,OBS_ITY,INDEX_HEADER)
         IF ( IDATYP .EQ. 169 ) THEN

            ! Basic geometric variables of the profile
            
            Rad  = obs_headElem_r(lobsSpaceData,OBS_TRAD,INDEX_HEADER)

            ! Loops over data in the observation
            
            IDATA   = obs_headElem_i(lobsSpaceData,OBS_RLN,INDEX_HEADER)
            IDATEND = obs_headElem_i(lobsSpaceData,OBS_NLV,INDEX_HEADER) + IDATA - 1

            ! Scan for requested assimilations, and count them

            DO INDEX_BODY= IDATA, IDATEND
               IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY) == obs_assimilated ) THEN
                  HNH1 = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)
                  IF (LEVELGPSRO.EQ.1) HNH1 = HNH1-Rad

                  ! Increment OMF = Y - H(x)

                  ZOMF = obs_bodyElem_r(lobsSpaceData,OBS_OMP,INDEX_BODY)

                  ! Observation value    Y

                  ZOBS = obs_bodyElem_r(lobsSpaceData,OBS_VAR,INDEX_BODY)
                  ZOER = obs_bodyElem_r(lobsSpaceData,OBS_OER,INDEX_BODY)
                  ZMHX = ZOBS-ZOMF

                  ! Reference order of magnitude value:

                  IF (LEVELGPSRO.EQ.1) THEN
                     ZREF = 0.025d0*exp(-HNH1/6500.d0)
                  ELSE
                     IF (.NOT.gpsroBNorm) THEN
                       ZREF = 300.d0*exp(-HNH1/6500.d0)
                     ELSE
                       ZREF = ZMHX
                     ENDIF
                  ENDIF
                           
                  ! OMF Tested criteria:

                  IF (.NOT.gpsroBNorm) THEN
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
      
  end subroutine bgck_gpsro

  !--------------------------------------------------------------------------
  ! isetflag
  !--------------------------------------------------------------------------
  function isetflag(cdfam,kodtyp,kvnam,zbgchk)
      !
      !:Purpose: Set BACKGROUND-CHECK FLAGS According to values set in a table.
      !          Original values in table come from ecmwf.
      !

      IMPLICIT NONE
      integer isetflag

      ! Arguments:
      character(len=2) :: cdfam ! FAMILY  NAME ( 'UA' , 'AI'   ...etc.. )
      integer :: kodtyp ! BURP CODE TYPE
      integer :: kvnam  ! VARIABLE NAME ( BURP )
      real(8) :: zbgchk ! NORMALIZED BACKGROUND DEPARTURE

      ! Locals:      
      real(8), parameter :: zsacrit(3) = (/ 10.00D0, 20.00D0, 30.00D0 /)
      real(8), parameter :: zttcrit(3) = (/  9.00D0, 16.00D0, 25.00D0 /)
      real(8), parameter :: zalcrit(3) = (/ 10.00D0, 20.00D0, 30.00D0 /)
      real(8), parameter :: zescrit(3) = (/ 10.00D0, 20.00D0, 30.00D0 /)
      real(8), parameter :: zpscrit(3) = (/  9.00D0, 16.00D0, 25.00D0 /)
      real(8), parameter :: zpncrit(3) = (/ 10.00D0, 20.00D0, 30.00D0 /)
      real(8), parameter :: zswcrit(3) = (/ 10.00D0, 20.00D0, 30.00D0 /)
      real(8), parameter :: ztscrit(3) = (/  5.00D0, 25.00D0, 30.00D0 /)
      real(8), parameter :: zdzcrit(3) = (/  2.25D0,  5.06D0,  7.56D0 /)
      real(8), parameter :: zgzcrit(3) = (/ 12.25D0, 25.00D0, 36.00D0 /)
      real(8), parameter :: zzdcrit(3) = (/  9.00D0, 16.00D0, 25.00D0 /)
      real(8), parameter :: zchcrit(3) = (/  9.00D0, 16.00D0, 25.00D0 /)
      real(8), parameter :: zLogViscrit(3) = (/ 10.00D0, 20.00D0, 30.00D0 /)

      !real(8) :: zuvcrit(3) = (/ 10.00D0, 20.00D0, 30.00D0 /)
      real(8) :: zuvcrit(3)
      
      zuvcrit(1) = 10.00D0
      zuvcrit(3) = 30.00D0
      if ( kodtyp == 37 ) then
        zuvcrit(2)=25.00D0
      else
        zuvcrit(2)=20.00D0
      end if
      
      isetflag=0

      if ( kvnam == BUFR_NEGZ ) then
      
         ! SET FLAG FOR HEIGHTS
         
         if ( zbgchk >= zgzcrit(1) .and. zbgchk < zgzcrit(2) ) then
           isetflag=1
         else if ( zbgchk >= zgzcrit(2) .and. zbgchk < zgzcrit(3) ) then
           isetflag=2
         else if ( zbgchk >= zgzcrit(3) )then
           isetflag =3
         endif

      else if ( kvnam == BUFR_NETT ) then
      
         ! SET FLAG FOR TEMPERATURE

         if ( zbgchk >= zttcrit(1) .and. zbgchk < zttcrit(2) ) then
           isetflag=1
         else if ( zbgchk >= zttcrit(2) .and. zbgchk < zttcrit(3) ) then
           isetflag=2
         else if ( zbgchk >= zttcrit(3) )then
           isetflag =3
         endif

      else if ( kvnam == BUFR_NEDZ ) then
     
         ! SET FLAG FOR SATEMS
         
         if ( zbgchk >= zdzcrit(1) .and. zbgchk < zdzcrit(2) ) then
           isetflag=1
         else if ( zbgchk >= zdzcrit(2) .and. zbgchk < zdzcrit(3) ) then
           isetflag=2
         else if ( zbgchk >= zdzcrit(3) )then
           isetflag =3
         endif

      else if ( kvnam == BUFR_NEFS ) then

         ! SET FLAG FOR WIND SPEED

         if ( zbgchk >= zsacrit(1) .and. zbgchk < zsacrit(2) ) then
           isetflag=1
         else if ( zbgchk >= zsacrit(2) .and. zbgchk < zsacrit(3) ) then
           isetflag=2
         else if ( zbgchk >= zsacrit(3) )then
           isetflag =3
         endif

      else if ( kvnam == BUFR_NEUU .or. kvnam == BUFR_NEVV ) then

         ! SET FLAG FOR WIND COMPONENTS

         if ( zbgchk >= zuvcrit(1) .and. zbgchk < zuvcrit(2) ) then
           isetflag=1
         else if ( zbgchk >= zuvcrit(2) .and. zbgchk < zuvcrit(3) ) then
           isetflag=2
         else if ( zbgchk >= zuvcrit(3) )then
           isetflag =3
         endif
         
      else if ( kvnam == BUFR_NEAL ) then

         ! SET FLAG FOR ALADIN HLOS WIND OBSERVATIONS
 
         if ( zbgchk >= zalcrit(1) .and. zbgchk < zalcrit(2) ) then
           isetflag=1
         else if ( zbgchk >= zalcrit(2) .and. zbgchk < zalcrit(3) ) then
           isetflag=2
         else if ( zbgchk >= zalcrit(3) )then
           isetflag=3
         end if

      else if ( kvnam == BUFR_NEUS .or. kvnam == BUFR_NEVS ) then

         ! SET FLAG FOR SURFACE WIND COMPONENTS

         if ( zbgchk >= zswcrit(1) .and. zbgchk < zswcrit(2) ) then
           isetflag=1
         else if ( zbgchk >= zswcrit(2) .and. zbgchk < zswcrit(3) ) then
           isetflag=2
         else if ( zbgchk >= zswcrit(3) )then
           isetflag =3
         endif

      else if ( kvnam == bufr_gust ) then
      
         ! SET FLAG FOR SURFACE WIND GUST

         if ( zbgchk >= zswcrit(1) .and. zbgchk < zswcrit(2) ) then
           isetflag=1
         else if ( zbgchk >= zswcrit(2) .and. zbgchk < zswcrit(3) ) then
           isetflag=2
         else if ( zbgchk >= zswcrit(3) )then
           isetflag =3
         endif

      else if ( kvnam == bufr_nees ) then

         ! SET FLAG FOR DEW POINT DEPRESSION

         if ( zbgchk >= zescrit(1) .and. zbgchk < zescrit(2) ) then
           isetflag=1
         else if ( zbgchk >= zescrit(2) .and. zbgchk < zescrit(3) ) then
           isetflag=2
         else if ( zbgchk >= zescrit(3) )then
           isetflag =3
         endif

      else if ( kvnam == bufr_neps ) then

         ! SET FLAG FOR SURFACE PRESSURE
         
         if ( zbgchk >= zpscrit(1) .and. zbgchk < zpscrit(2) ) then
           isetflag=1
         else if ( zbgchk >= zpscrit(2) .and. zbgchk < zpscrit(3) ) then
           isetflag=2
         else if ( zbgchk >= zpscrit(3) )then
           isetflag =3
         endif

      else if ( kvnam == bufr_nepn ) then

         ! SET FLAG FOR MEAN SEA LEVEL PRESSURE

         if ( zbgchk >= zpncrit(1) .and. zbgchk < zpncrit(2) ) then
           isetflag=1
         else if ( zbgchk >= zpncrit(2) .and. zbgchk < zpncrit(3) ) then
           isetflag=2
         else if ( zbgchk >= zpncrit(3) )then
           isetflag =3
         endif

      else if ( kvnam == bufr_nets ) then

         ! SET FLAG FOR SURFACE TEMPERATURE

         if ( zbgchk >= ztscrit(1) .and. zbgchk < ztscrit(2) ) then
           isetflag=1
         else if ( zbgchk >= ztscrit(2) .and. zbgchk < ztscrit(3) ) then
           isetflag=2
         else if ( zbgchk >= ztscrit(3) )then
           isetflag =3
         endif

      else if ( kvnam .eq. BUFR_NESS ) then
      
         ! SET FLAG FOR SURFACE DEW POINT DEPRESSION

         if ( zbgchk >= zescrit(1) .and. zbgchk < zescrit(2) ) then
           isetflag=1
         else if ( zbgchk >= zescrit(2) .and. zbgchk < zescrit(3) ) then
           isetflag=2
         else if ( zbgchk >= zescrit(3) )then
           isetflag =3
         endif

      else if ( kvnam == bufr_vis ) then

         ! SET FLAG FOR VISIBILITY

         if ( zbgchk >= zLogViscrit(1) .and. zbgchk < zLogViscrit(2) ) then
           isetflag=1
         else if ( zbgchk >= zLogViscrit(2) .and. zbgchk < zLogViscrit(3) ) then
           isetflag=2
         else if ( zbgchk >= zLogViscrit(3) )then
           isetflag =3
         endif

      else if ( kvnam == bufr_nezd ) then

         ! SET FLAG FOR GB-GPS ZENITH DELAY

         if ( zbgchk >= zzdcrit(1) .and. zbgchk < zzdcrit(2) ) then
           isetflag=1
         else if ( zbgchk >= zzdcrit(2) .and. zbgchk < zzdcrit(3) ) then
           isetflag=2
         else if ( zbgchk >= zzdcrit(3) )then
           isetflag =3
         endif

       else if ( cdfam .eq. 'CH' ) then
      
         ! SET FLAG FOR CHEMICAL CONSTITUENTS

         if ( zbgchk >= zchcrit(1) .and. zbgchk < zchcrit(2) ) then
           isetflag=1
         else if ( zbgchk >= zchcrit(2) .and. zbgchk < zchcrit(3) ) then
           isetflag=2
         else if ( zbgchk >= zchcrit(3) )then
           isetflag =3
         endif
         
      endif

      return
  end function isetflag
      
end module backgroundCheck_mod

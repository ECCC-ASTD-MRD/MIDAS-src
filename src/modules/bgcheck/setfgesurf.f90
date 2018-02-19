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
!! *Purpose*: Interpolate vertically the contents of "column" to
!!            the pressure levels of the observations.
!!            A linear interpolation in ln(p) is performed.
!!
!! @author P. Koclas *CMC/AES  September 2000
!!
!--------------------------------------------------------------------------
      SUBROUTINE SETFGESURF(lcolumn,lcolumng,lobsSpaceData)
      use EarthConstants_mod
      use MathPhysConstants_mod
      use bufr_mod
      use columnData_mod
      use obsSpaceData_mod
      IMPLICIT NONE

      type(struct_columnData) :: lcolumn,lcolumng
      type(struct_obs) :: lobsSpaceData
      INTEGER IPB,IPT,IDIM
      INTEGER INDEX_HEADER,IK
      INTEGER INDEX_BODY,ITYP
      REAL*8 ZWB,ZWT
      REAL*8 ZLEV,ZPT,ZPB,ZHHH
      CHARACTER(len=2) :: cfam,varLevel
      LOGICAL LLOK

      ! loop over all body rows
      BODY: do index_body=1,obs_numbody(lobsSpaceData)
        cfam=obs_getFamily(lobsSpaceData,bodyIndex=index_body)
        if(cfam.eq.'SF' .or. &
           cfam.eq.'UA' .or. &
           cfam.eq.'SC' .or. &
           cfam.eq.'GP') then
!
!         Process all data within the domain of the model (excluding GB-GPS
!         ZTD data)
!
          LLOK=.FALSE.
          IF ( obs_bodyElem_i(lobsSpaceData,OBS_VCO,index_body) == 1 ) THEN
            ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body)
            IF (ITYP == BUFR_NETS .OR. ITYP == BUFR_NEPS .OR.  &
                ITYP == BUFR_NEPN .OR. ITYP == BUFR_NESS .OR.  &
                ITYP == BUFR_NEUS .OR. ITYP == BUFR_NEVS .OR.  &
                ITYP == BUFR_NEFS .OR. ITYP == BUFR_NEDS ) THEN
              LLOK=(obs_bodyElem_i(lobsSpaceData,OBS_ASS,index_body) .EQ. 1 )
            ELSEIF ( ITYP == BUFR_NEZD ) THEN
              ! make sure total zenith delay (from ground-based GPS) not treated
              LLOK=.FALSE.
            ELSE
              LLOK=(obs_bodyElem_i(lobsSpaceData,OBS_ASS,index_body) .EQ. 1 .AND.  &
                    obs_bodyElem_i(lobsSpaceData,OBS_XTR,index_body) .ge. 0)
              if(llok) write(*,*) 'setfgesurf: WARNING!!! unknown obs seen'
              if(llok) write(*,*) 'setfgesurf: ityp=',ityp,', cfam=',cfam
            ENDIF

            IF ( LLOK ) THEN
              INDEX_HEADER = obs_bodyElem_i(lobsSpaceData,OBS_HIND,index_body)
              ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body)
              varLevel = vnl_varLevelFromVarnum(ityp)
              idim=1
              if (varLevel.eq.'SF') idim=0
              IK   = obs_bodyElem_i(lobsSpaceData,OBS_LYR,index_body)
              ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,index_body)
              ZHHH = ZLEV * GRAV

              IF (ITYP == BUFR_NETS .OR. ITYP == BUFR_NEPS .OR.   &
                  ITYP == BUFR_NEPN .OR. ITYP == BUFR_NESS .OR.   &
                  ITYP == BUFR_NEUS .OR. ITYP == BUFR_NEVS    ) THEN

                IPT  = IK + col_getOffsetFromVarno(lcolumng,ityp)
                IPB  = IPT+1
                call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,col_getElem(lcolumn,IPB,INDEX_HEADER))
              ELSE
                IPT  = IK + col_getOffsetFromVarno(lcolumng,ityp)
                IPB  = IPT+1
                ZPT  = col_getHeight(lcolumng,IK,INDEX_HEADER,varLevel)
                ZPB  = col_getHeight(lcolumng,IK+1,INDEX_HEADER,varLevel)
                ZWB  = idim*(ZPT-ZHHH)/(ZPT-ZPB)
                ZWT  = 1.d0 - ZWB
                IF ( obs_bodyElem_i(lobsSpaceData,OBS_XTR,index_body) .eq. 0) then
                  call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,   &
                    zwb*col_getElem(lcolumn,IPB,INDEX_HEADER) + ZWT*col_getElem(lcolumn,IPT,INDEX_HEADER))
                ELSE
                  call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,   &
                    col_getElem(lcolumn,IK + col_getOffsetFromVarno(lcolumng,ityp),INDEX_HEADER))
                ENDIF
                if(obs_elem_c(lobsSpaceData,'STID',index_header) .eq. '99999999') then
                  write(*,*) 'setfgesurf:stn,ityp,xtr,ipt,ipb,zwt,zwb'   &
                       ,obs_elem_c(lobsSpaceData,'STID',index_header),ityp,    &
                        obs_bodyElem_i(lobsSpaceData,OBS_XTR,index_body),ipt,ipb,zwt,zwb
                  write(*,*) 'setfgesurf:gobs(ipb),gobs(ipt),fge'   &
                       ,col_getElem(lcolumn,IPB,INDEX_HEADER),col_getElem(lcolumn,IPT,INDEX_HEADER)      &
                       ,obs_bodyElem_r(lobsSpaceData,OBS_HPHT,index_body)
                endif
              ENDIF
            ENDIF
          ENDIF

        ENDIF

      ENDDO BODY

      RETURN
      END SUBROUTINE SETFGESURF

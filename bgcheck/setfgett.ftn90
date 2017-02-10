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
!!            the pressure levels of the observations. Then
!!            compute THE FIRST GUESS ERROR VARIANCES
!!            A linear interpolation in ln(p) is performed.
!!
!! @author P. Koclas *CMC/CMSV November 1998
!!
!--------------------------------------------------------------------------
      SUBROUTINE SETFGETT(lcolumn,lcolumng,lobsSpaceData)
      use bufr
      use columnData_mod
      use obsSpaceData_mod
      IMPLICIT NONE
      type(struct_columnData) :: lcolumn,lcolumng
      type(struct_obs) :: lobsSpaceData
      INTEGER IPB,IPT
      INTEGER INDEX_HEADER,ITYP,IK,IBEGIN,ILAST
      INTEGER J,INDEX_BODY
      REAL*8 ZWB,ZWT
      REAL*8 ZLEV,ZPB,ZPT
      character(len=2) :: varLevel

      ! loop over all body rows
      BODY: do index_body=1,obs_numbody(lobsSpaceData)
!
!*    1. Computation of sigmap
!     .  -----------------------------
!
         IF ( (obs_bodyElem_i(lobsSpaceData,OBS_ASS,index_body) .EQ. 1) .and.    &
              (obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body).EQ. BUFR_NETT) ) THEN

            IF ( (obs_bodyElem_i(lobsSpaceData,OBS_XTR,index_body) .NE. 0) .and.    &
                 (obs_bodyElem_i(lobsSpaceData,OBS_VCO,index_body) .EQ. 2) ) THEN
               ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body)
               varLevel = vnl_varLevelFromVarnum(ityp)
               IK=col_getNumLev(lcolumng,varLevel)-1
               INDEX_HEADER = obs_bodyElem_i(lobsSpaceData,OBS_HIND,index_body)
               IPT  = IK + col_getOffsetFromVarno(lcolumng,ityp)
               IPB  = IPT +1
               call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,col_getElem(lcolumn,IPB,INDEX_HEADER))
            ELSE
               INDEX_HEADER = obs_bodyElem_i(lobsSpaceData,OBS_HIND,index_body)
               ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,index_body)
               IK   = obs_bodyElem_i(lobsSpaceData,OBS_LYR,index_body)
               IPT  = IK
               IPB  = IPT+1
               ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body)
               varLevel = vnl_varLevelFromVarnum(ityp)
               ZPT  = col_getPressure(lcolumng,IK,INDEX_HEADER,varLevel)
               ZPB  = col_getPressure(lcolumng,IK+1,INDEX_HEADER,varLevel)
               ZWB  = LOG(ZLEV/ZPT)/LOG(ZPB/ZPT)
               ZWT  = 1.0D0 - ZWB
!
!              FIRST GUESS ERROR VARIANCE
!
               call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,   &
                  (ZWB*col_getElem(lcolumn,IPB,INDEX_HEADER,'TT') + ZWT*col_getElem(lcolumn,IPT,INDEX_HEADER,'TT')))
            ENDIF

         ENDIF

      END DO BODY

      RETURN
      END SUBROUTINE SETFGETT

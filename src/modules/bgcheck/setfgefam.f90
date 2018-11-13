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
!--------------------------------------------------------------------------
      SUBROUTINE SETFGEFAM(CDFAM,lcolumn,lcolumng,lobsSpaceData)
      use bufr_mod
      use columnData_mod
      use obsSpaceData_mod
      use utilities_mod
      use EarthConstants_mod 
      IMPLICIT NONE
      type(struct_columnData) :: lcolumn,lcolumng
      type(struct_obs) :: lobsSpaceData
      CHARACTER*2 CDFAM
      INTEGER IPB,IPT
      INTEGER INDEX_HEADER,ITYP,IK
      INTEGER INDEX_BODY
      REAL*8 ZWB,ZWT
      REAL*8 ZLEV,ZPB,ZPT
      character(len=2) :: varLevel

      ! loop over all header indices of the CDFAM family
      call obs_set_current_header_list(lobsSpaceData,CDFAM)
      HEADER: do
         index_header = obs_getHeaderIndex(lobsSpaceData)
         if (index_header < 0) exit HEADER

         ! loop over all body indices for this index_header
         call obs_set_current_body_list(lobsSpaceData, index_header)
         BODY: do 
            index_body = obs_getBodyIndex(lobsSpaceData)
            if (index_body < 0) exit BODY
!
!*    1. Computation of sigmap
!     .  -----------------------------
!
            IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,index_body) .EQ. 1 .AND.   &
                 obs_bodyElem_i(lobsSpaceData,OBS_VCO,index_body) .EQ. 2      ) then
            IF  (obs_bodyElem_i(lobsSpaceData,OBS_XTR,index_body) .NE. 0) THEN
               ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body)
               varLevel = vnl_varLevelFromVarnum(ityp)
               IK   = col_getNumLev(LCOLUMNG,varLevel)
               IPB  = IK + col_getOffsetFromVarno(lcolumng,ityp)
               if(ITYP .ne. BUFR_NEGZ) then
                 call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,col_getElem(lcolumn,IPB,INDEX_HEADER))
               else
                 call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,col_getHeight(lcolumn,IK,INDEX_HEADER,'TH')*RG)
               endif
            ELSE
               ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body)
               varLevel = vnl_varLevelFromVarnum(ityp)
               ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,index_body)
               IK   = obs_bodyElem_i(lobsSpaceData,OBS_LYR,index_body)
               IPT  = IK + col_getOffsetFromVarno(lcolumng,ityp)
               IPB  = IPT+1
               ZPT  = col_getPressure(lcolumng,IK,INDEX_HEADER,varLevel)
               ZPB  = col_getPressure(lcolumng,IK+1,INDEX_HEADER,varLevel)
               ZWB  = LOG(ZLEV/ZPT)/LOG(ZPB/ZPT)
               ZWT  = 1.0D0 - ZWB
!
!              FIRST GUESS ERROR VARIANCE
!
               if(ITYP .ne. BUFR_NEGZ) then
                 call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,   &
                      (ZWB*col_getElem(lcolumn,IPB,INDEX_HEADER) + ZWT*col_getElem(lcolumn,IPT,INDEX_HEADER)))
               else
                 call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,   &
                      (ZWB*col_getHeight(lcolumn,IK+1,INDEX_HEADER,'TH')*RG + ZWT*col_getHeight(lcolumn,IK,INDEX_HEADER,'TH')*RG))
               endif
               if(obs_bodyElem_r(lobsSpaceData,OBS_HPHT,index_body).le.0.d0) then
                 write(*,*) 'SETFGEFAM: CDFAM = ',CDFAM
                 write(*,*) 'SETFGEFAM: IPB,IPT,ZWB,ZWT,ITYP,ZLEV=',IPB,IPT,ZWB,ZWT,ITYP,ZLEV
                 write(*,*) 'SETFGEFAM: lcolumn_all(IPB,INDEX_HEADER)=',col_getElem(lcolumn,IPB,INDEX_HEADER)
                 write(*,*) 'SETFGEFAM: lcolumn_all(IPT,INDEX_HEADER)=',col_getElem(lcolumn,IPT,INDEX_HEADER)
                 write(*,*) 'SETFGEFAM: get_height(IK+1,INDEX_HEADER)=',col_getHeight(lcolumn,IK+1,INDEX_HEADER,'TH')*RG
                 write(*,*) 'SETFGEFAM: get_height(IK  ,INDEX_HEADER)=',col_getHeight(lcolumn,IK  ,INDEX_HEADER,'TH')*RG
                 CALL utl_abort('SETFGEFAM: First-guess stdev bad value')
               endif
            ENDIF
            ENDIF

         END DO BODY

      END DO HEADER

      RETURN
      END SUBROUTINE SETFGEFAM

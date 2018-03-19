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
!! MODULE obsUtil_mod (prefix="obs")
!!
!! *Purpose*: Common routines used by burpfiles_mod and sqlitefiles_mod
!!            
!!
!--------------------------------------------------------------------------

module obsUtil_mod
  
  use obsSpaceData_mod
  use bufr_mod
  use MathPhysConstants_mod
  use codePrecision_mod
  use EarthConstants_mod, only:  GRAV

  implicit none
  save
  private
  
  ! public procedures
  public :: VINT3DFD, SETASSFLG, FLAGUVTOFD_OBSDAT, FDTOUV_OBSDAT, ADJUST_HUM_GZ, ADJUST_SFVCOORD, set_err_gbgps, cvt_obs_instrum

  contains

  SUBROUTINE FLAGUVTOFD_OBSDAT(obsSpaceData)
!
!**s/r FLAGUVTOFD_OBSDAT  - Update WIND DIRECTION AND SPEED FLAGS
!
!
!Author  : P. Koclas *CMC/CMDA  April 2013
!
!
!Arguments
!
      IMPLICIT NONE
!
      type(struct_obs) :: obsSpaceData
      INTEGER :: IUU,IVV,IFF,IDD
      INTEGER :: FLAGU,FLAGV,NEWFLAG
      INTEGER :: INDEX_HEADER,ISTART,IEND,jwintyp
      INTEGER :: INDEX_BODY,INDEX_BODY2
      REAL*8  :: ZLEVU
      LOGICAL ::  LLOK
      CHARACTER*9 :: STID
!-----------------------------------------------------------------------
!
      WIND_TYPE: do jwintyp=1,2

         if (jwintyp == 1) then
            IUU=BUFR_NEUU
            IVV=BUFR_NEVV
            IDD=BUFR_NEDD
            IFF=BUFR_NEFF
         else
            IUU=BUFR_NEUS
            IVV=BUFR_NEVS
            IDD=BUFR_NEDS
            IFF=BUFR_NEFS
         end if
!
!
!
         BODY: DO INDEX_BODY=1,obs_numBody(obsSpaceData)

            LLOK= ( obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY) == IUU)

             FLAGU=-1
    !----------------
            IF ( LLOK ) THEN
    !----------------
               INDEX_HEADER = obs_bodyElem_i(obsSpaceData,OBS_HIND,INDEX_BODY)
               ISTART       = obs_headElem_i(obsSpaceData,OBS_RLN,INDEX_HEADER)
               IEND=obs_headElem_i(obsSpaceData,OBS_NLV,INDEX_HEADER) +ISTART-1
               STID=obs_elem_c(obsSpaceData,'STID',INDEX_HEADER)


               ZLEVU = obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY)
!
!****************************************************************************
!  GET FLAG OF U COMPONENT
!***********************************************************************
!
               FLAGU=obs_bodyElem_i(obsSpaceData,OBS_FLG,INDEX_BODY)

               BODY_2: DO INDEX_BODY2=ISTART,IEND
                  IF ( ( obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY2) == IVV) &
                 .AND. ( obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY2) == ZLEVU) ) THEN
!
!****************************************************************************
!  GET FLAG OF V COMPONENT
!***********************************************************************
!
                     FLAGV= obs_bodyElem_i(obsSpaceData,OBS_FLG,INDEX_BODY2)
                     NEWFLAG =IOR(FLAGU,FLAGV)
!   
                  END IF
               END DO BODY_2
!
!***********************************************************************
!                UPDATE FLAGS OF DIRECTION AN SPEED
!***********************************************************************
!
               BODY_2_2: DO INDEX_BODY2=ISTART,IEND
       !===============================================
                  IF ((obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY2) == IDD) &
                 .AND. obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY2) == ZLEVU ) THEN

                     NEWFLAG =IOR(FLAGU,FLAGV)
                     call obs_bodySet_i(obsSpaceData, OBS_FLG, INDEX_BODY2, NEWFLAG) 

                  END IF
	       !===============================================

       !===============================================
                  IF ((obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY2) == IFF) &
                 .AND. obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY2) == ZLEVU ) THEN

                     NEWFLAG =IOR(FLAGU,FLAGV)
                     call obs_bodySet_i(obsSpaceData,OBS_FLG,INDEX_BODY2, NEWFLAG)
                  END IF
	       !===============================================
               END DO BODY_2_2

	    !----------------
            END IF
	    !----------------

         END DO BODY

      END DO WIND_TYPE

  END SUBROUTINE FLAGUVTOFD_OBSDAT


  SUBROUTINE VINT3DFD(elem_i,obsSpaceData)
      !
      ! s/r VINT3DFD  - Computation of DIRECTION AND SPEED RESIDUALS
      !
      ! Author  : P. Koclas *CMC/AES  September 1999
      ! Revision:
      !     1.0  P. Koclas CMC :  September 2000
      !                 -remove quality control flag and (ff dd) component initializtions
      !          JM Belanger CMDA/SMC  Jan 2001
      !                   . 32 bits conversion
      !
      !     Purpose:  -Compute direction and speed residuals from u and
      !                v residuals.
      !
      implicit none

      type(struct_obs) :: obsSpaceData
      integer, intent(in) :: elem_i
      INTEGER IUU,IVV,IFF,IDD
      INTEGER INDEX_HEADER,ISTART,IEND,jwintyp
      INTEGER INDEX_BODY,INDEX_BODY2
      REAL*8 ZLEVU
      REAL*8 MODUL,ANG,UU,VV
      LOGICAL LLOK

      WIND_TYPE: do jwintyp=1,2

         if (jwintyp == 1) then
            IUU=BUFR_NEUU
            IVV=BUFR_NEVV
            IDD=BUFR_NEDD
            IFF=BUFR_NEFF
         else
            IUU=BUFR_NEUS
            IVV=BUFR_NEVS
            IDD=BUFR_NEDS
            IFF=BUFR_NEFS
         end if

         ! Process all data within the domain of the model

         BODY: DO INDEX_BODY=1,obs_numBody(obsSpaceData)
            LLOK= (obs_bodyElem_i(obsSpaceData,OBS_ASS,INDEX_BODY) == 1)  &
            .AND. (obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY) == IUU)
            IF ( LLOK ) THEN
               INDEX_HEADER = obs_bodyElem_i(obsSpaceData,OBS_HIND,INDEX_BODY)
               ISTART=obs_headElem_i(obsSpaceData,OBS_RLN,INDEX_HEADER)
               IEND=obs_headElem_i(obsSpaceData,OBS_NLV,INDEX_HEADER) +ISTART-1
               ZLEVU = obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY)
               UU=-obs_bodyElem_r(obsSpaceData,elem_i,INDEX_BODY) +  &
                   obs_bodyElem_r(obsSpaceData,OBS_VAR,INDEX_BODY)
               BODY_2: DO INDEX_BODY2=ISTART,IEND
                  IF ((obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY2) == IVV)  &
                 .AND.(obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY2) == ZLEVU)) THEN
                   VV=-obs_bodyElem_r(obsSpaceData,elem_i,INDEX_BODY2) +  &
                       obs_bodyElem_r(obsSpaceData,OBS_VAR,INDEX_BODY2)

                     ! 1-calculate angle

                     MODUL=SQRT((UU**2)+(VV**2))
                     IF (MODUL == 0.) THEN
                        ANG=0.0D0
                     ELSE
                        ANG=ATAN2(VV,UU)
                        ANG= (270.0D0 - ANG  * MPC_DEGREES_PER_RADIAN_R8 )

                        ! 2-Change to meteorological definition of wind direction.

                        IF (ANG > 360.0D0) ANG=ANG-360.0D0
                        IF (ANG <= 0.0D0)   ANG=ANG+360.0D0
                     END IF
   
                  END IF
               END DO BODY_2

               ! insert resduals into obsSpaceData

               BODY_2_2: DO INDEX_BODY2=ISTART,IEND
                  IF ((obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY2) == IDD)  &
                 .AND. obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY2) == ZLEVU ) THEN

                     call obs_bodySet_r(obsSpaceData, elem_i, INDEX_BODY2,    &
                          obs_bodyElem_r(obsSpaceData,OBS_VAR,INDEX_BODY2) - ANG )

                     IF ( obs_bodyElem_r(obsSpaceData,elem_i,INDEX_BODY2) >  180.0d0)  &
                        call obs_bodySet_r(obsSpaceData, elem_i, INDEX_BODY2,   &
                                       obs_bodyElem_r(obsSpaceData,elem_i,INDEX_BODY2)-360.0d0)
                     IF ( obs_bodyElem_r(obsSpaceData,elem_i,INDEX_BODY2) <= -180.0d0)  &
                        call obs_bodySet_r(obsSpaceData, elem_i, INDEX_BODY2,  &
                                       obs_bodyElem_r(obsSpaceData,elem_i,INDEX_BODY2)+360.0d0)

                      call obs_bodySet_r(obsSpaceData, elem_i, INDEX_BODY2, -1.0d0*  &
                                         obs_bodyElem_r(obsSpaceData,elem_i,INDEX_BODY2))

                      call obs_bodySet_r(obsSpaceData,OBS_OER,INDEX_BODY2,1.0d0)
                      call obs_bodySet_i(obsSpaceData,OBS_ASS,INDEX_BODY2, 1)
                      call obs_bodySet_i(obsSpaceData,OBS_FLG,INDEX_BODY2, 0)
                  END IF
                  IF ((obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY2) == IFF)  &
                 .AND. obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY2) == ZLEVU ) THEN
                     call obs_bodySet_r(obsSpaceData,elem_i, INDEX_BODY2,   &
                          obs_bodyElem_r(obsSpaceData,OBS_VAR,INDEX_BODY2) - MODUL)
                     call obs_bodySet_r(obsSpaceData,OBS_OER,INDEX_BODY2,1.0d0)
                     call obs_bodySet_i(obsSpaceData,OBS_ASS,INDEX_BODY2, 1)
                     call obs_bodySet_i(obsSpaceData,OBS_FLG,INDEX_BODY2, 0)
                  END IF
               END DO BODY_2_2
            END IF

         END DO BODY

      END DO WIND_TYPE

  END SUBROUTINE VINT3DFD

  subroutine setassflg(obsSpaceData)
    ! Purpose:  Set banco quality control bit #12 for all data assimilated
    !           by current analysis.
    implicit none

    type(struct_obs) :: obsSpaceData
    integer :: index_body

    ! Process all data
    do index_body=1,obs_numBody(obsSpaceData)
      if (obs_bodyElem_i(obsSpaceData,OBS_ASS,index_body) == 1)  then
        call obs_bodySet_i(obsSpaceData,OBS_FLG,index_body,ibset( obs_bodyElem_i(obsSpaceData,OBS_FLG,index_body), 12 ))
      end if
    end do

  end subroutine setassflg

  subroutine fdtouv_obsdat(obsdat,start,end,ppmis)
!
!---------------------------------------------------------------
!
! Author  : P. Koclas, CMC/CMDA December  2012
!           CONVERT DD , FF  WINDS TO
!            UU (est-west),  VV (north-south) COMPONENTS
!
!    ARGUMENTS:
!                 INPUT:
!
!                       -obsdat     : CMA_table INSTANCE 
!                       -START     : FIRST OBERVATION
!                       -END       : LAST  OBERVATION
!                       -PPMIS     : MISSING VALUE  
!
!        **************************************************
!         IT IS ASSUMED THAT CMA CONTAINS ENTRIES   FOR 
!          UU AND VV  with observed values = missing value
!        **************************************************
!
!---------------------------------------------------------------
!
    implicit none
    type (struct_obs), intent(inout) :: obsdat
    real(4)        :: ppmis
    integer(4)     :: start,end

    integer(4)     :: varno,varno2,varno4
    real(4)        :: obsuv
    integer(4)     :: jo,rln,nlv,j,j2,j4,jpos,ilem
    integer(4)     :: ddflag,ffflag,newflag,uuflag,vvflag
    integer(4)     :: ilemf,ilemu,ilemv,indu_misg,indv_misg,indum,indvm
    logical        :: llmisdd,llmisff,llmis,lluv_misg,llu_misg,llv_misg
    logical        :: lluv_present,llu_present,llv_present
    real(obs_real) :: uu,vv,dd,ff
    real(obs_real) :: level_dd,level4,level,level_uu

    ffflag = 0  ! bhe 

    HEADER1: do jo = start, end
        
      rln = obs_headElem_i(obsdat,OBS_RLN,jo)
      nlv = obs_headElem_i(obsdat,OBS_NLV,jo)

      BODY1: do j = rln, nlv + rln -1
        dd = ppmis
        ff = ppmis

        varno = obs_bodyElem_i(obsdat,OBS_VNM,j)
        llmisdd = .true.

        if ( varno /= bufr_nedd .and. varno /= bufr_neds ) cycle BODY1

        if( varno == bufr_neds) then
          ilemf = bufr_nefs
          ilemu = bufr_neus
          ilemv = bufr_nevs
        else
          ilemf = bufr_neff
          ilemu = bufr_neuu
          ilemv = bufr_nevv
        end if

        dd       = obs_bodyElem_r(obsdat,OBS_VAR,j)
        ddflag   = obs_bodyElem_i(obsdat,OBS_FLG,j)
        level_dd = obs_bodyElem_r(obsdat,OBS_PPP,j)
          
        llu_misg = .false.
        llv_misg = .false.
        llu_present = .false.
        llv_present = .false.
        indum = -1
        indvm = -1

        ! FIND IF  U AND V ARE ALREADY IN CMA
        UVINOBSDAT: do j4 = j, nlv + rln -1

          level4 = obs_bodyElem_r(obsdat, OBS_PPP,j4)
          if (level4 == level_dd) then
            varno4 = obs_bodyElem_i(obsdat, OBS_VNM,j4)
            select case (varno4)
              case (11003,11004,11002,11001,11215,11216,11011,11012)

                obsuv = obs_bodyElem_r(obsdat, OBS_VAR,j4)
                if (  (varno4 == ilemu)     .and.  (obsuv /= ppmis) ) then
                  llu_present = .true.
                  indum = j4
                else if ( (varno4 == ilemv) .and. (obsuv /= ppmis) ) then
                  llv_present = .true.
                  indvm = j4
                end if
                
                if (  (varno4 == ilemu)     .and. (obsuv == ppmis) ) then
                  llu_misg = .true.
                  indu_misg = j4
                else if ( (varno4 == ilemv)  .and. (obsuv == ppmis) ) then
                  llv_misg = .true.
                  indv_misg = j4
                end if

            end select
          end if

        end do UVINOBSDAT

        lluv_misg = (llu_misg .and. llv_misg)
        lluv_present = (llu_present .and. llv_present)

        if ( lluv_misg) then

          CALCUV: do j2 = j, nlv + rln -1

            llmisff = .true.
            llmisdd = .true.
            llmis   = .true.
            level = obs_bodyElem_r(obsdat,OBS_PPP,j2)
            if ( level /= level_dd) cycle
            varno2 = obs_bodyElem_i(obsdat,OBS_VNM,j2)
            if ( varno2 == ilemf ) then

              ff = obs_bodyElem_r(obsdat,OBS_VAR,j2)
              ffflag = obs_bodyElem_i(obsdat,OBS_FLG,j2)
              IF ( (dd == 0.  .and. ff > 0.) .or. ( dd > 360. .or. dd  < 0.) ) then
                llmisdd = .true.
                llmisff = .true.
              else if ( dd == ppmis .or. ff == ppmis) then
                llmisdd = .true.
                llmisff = .true.
              else
                llmisdd = .false.
                llmisff = .false.
              end if

              ! IF SPEED = 0 CALM WIND IS ASSUMED.
              if (ff == 0.0) then
                dd = 0.
              end if
                   
              dd = dd + 180.
              if ( dd > 360.) dd = dd - 360.
              dd = dd * mpc_radians_per_degree_r8
                
              ! U,V COMPONENTS ARE
              uu = ff * sin(dd)
              vv = ff * cos(dd)
              if  ( ( llmisdd == .true.) .or. ( llmisff == .true. ) ) then
                llmis = .true.
                if ( indu_misg > 0 .or. indv_misg > 0 ) then
                  call obs_bodySet_i(obsdat,OBS_VNM,INDU_MISG,-1)
                  call obs_bodySet_i(obsdat,OBS_VNM,INDV_MISG,-1)
                end if
              else
                llmis = .false.
              end if

            end if
            newflag = ior(ddflag,ffflag)

            if ( indum > 0 .or. indvm > 0 ) then
              call obs_bodySet_i(obsdat,OBS_VNM,indu_misg,-1)
              call obs_bodySet_i(obsdat,OBS_VNM,indv_misg,-1)
            end if
            if (llmis) then
              if ( indum > 0 .or. indvm > 0 ) then
                call obs_bodySet_i(obsdat,OBS_FLG,induM,newflag)
                call obs_bodySet_i(obsdat,OBS_FLG,indvM,newflag)
              end if
            else if (.not.llmis) then
              call obs_bodySet_r(obsdat,OBS_VAR,indu_misg,uu)
              call obs_bodySet_i(obsdat,OBS_FLG,indu_misg,newflag)

              call obs_bodySet_r(obsdat,OBS_VAR,indv_misg,vv)
              call obs_bodySet_i(obsdat,OBS_FLG,indv_misg,newflag)
            end if

          end do CALCUV

        else                       

          if ( lluv_present ) then
            call obs_bodySet_i(obsdat,OBS_VNM,indu_misg,-1)
            call obs_bodySet_i(obsdat,OBS_VNM,indv_misg,-1)
          ELSE
            if (indum > 0) then
              call obs_bodySet_i(obsdat,OBS_VNM,indum,-1)
            end if
            if (indvm > 0) then
              call obs_bodySet_i(obsdat,OBS_VNM,indvm,-1)
            end if
          end if

        end if

      end do BODY1

    end do HEADER1


    do jo = start, end

      rln = obs_headElem_i(obsdat,OBS_RLN,jo)
      nlv = obs_headElem_i(obsdat,OBS_NLV,jo)

      do j = rln, nlv + rln -1

        llmisdd = .true.
        varno = obs_bodyElem_i(obsdat,OBS_VNM,j)
        level = obs_bodyElem_r(obsdat,OBS_PPP,j)

        select case (varno)

          case (bufr_neuu)
            ilem = bufr_nevv
          case (bufr_neus)
            ilem = bufr_nevs        
          case default
            cycle

        end select

        Jpos = -1

        ! TRANSFER THE FLAG BITS  FROM ONE WIND COMPONENT TO THE OTHER
        do j4 = rln, nlv + rln -1
          uu       = obs_bodyElem_r(obsdat,OBS_VAR,j4)
          level_uu = obs_bodyElem_r(obsdat,OBS_PPP,j4)
          Jpos = -1
        
          if ( level_uu == level .and. uu == ppmis ) then
            call obs_bodySet_i(obsdat,OBS_VNM,j4,-1)
          end if

          if ( level_uu == level .and. uu /= ppmis ) then
            uuflag = obs_bodyElem_i(obsdat,OBS_FLG,j4)
            varno2 = obs_bodyElem_i(obsdat,OBS_VNM,j4)

            if ( ilem == varno2 ) then
              vvflag  = obs_bodyElem_i(obsdat,OBS_FLG,j)
              newflag = ior(uuflag,vvflag)
              call obs_bodySet_i(obsdat,OBS_FLG,j, newflag)
              call obs_bodySet_i(obsdat,OBS_FLG,j4,newflag)
              jpos = j4
              exit
            end if

          end if

        end do !j4

        ! ELIMINATE ENTRIES WHERE ONE COMPONENT OF WIND (UU OR VV) IS MISSING
        if (jpos < 0) then
          write(*,*) ' eliminate winds for station : ', obs_elem_c(obsdat,'STID',JO),  &
            obs_bodyElem_i(obsdat,OBS_VNM,J), obs_bodyElem_r(obsdat,OBS_PPP,J)
          call obs_bodySet_i(obsdat,OBS_VNM,j,-1)
        end if

      end do !j

    end do !jo

  end subroutine fdtouv_obsdat
  real function surfvcord(ilem,idtyp)
      implicit none
      integer :: ilem,idtyp,type
      real :: vcordsf2
!***********************************************************************
!
!      PURPOSE: SEt vertical coordinate for surface data.
!
!       AUTHOR:   P. KOCLAS (CMC/CMDA) December 2011
!
!       Revision : 
!
!    ARGUMENTS:
!               INPUT:
!                      -ILEMP   : BURP ELEMENT NUMBER
!                      -IDTYP   : BURP CODETYPE
!
!               OUTPUT:
!                      -SURFVCORD
!
!
!***********************************************************************
!
!
!     GENERATE TABLES TO ADJUST VERTICAL COORDINATE OF SURFACE DATA
!
!     DEFAULT VALUE 
      vcordsf2=0.

      select case(idtyp)
        case(135,136,137,138,32,34,35,37,38,159,160,161,162)
          ! UPPER AIR LAND
          type=3

        case(139,140,141,142,33,36)
          ! UPPER AIR SHIP
          type=4

        case(12,14,146)
          ! SYNOPS
          type=1

        case(13,18,145,147)
          ! SHIPS
          type=2

        case(254)
          ! SCATTEROMETER WINDS
          type=5

        case default
          type=-99

      end select

      select case(type)
         case (1)
           select case(ilem)
             case (bufr_neds,bufr_nefs,bufr_neus,bufr_nevs)
               ! us,vs,ffs,dds
               vcordsf2=10.0

             case (bufr_nepn)
               vcordsf2=0.0

             case (bufr_neps)
               vcordsf2=0.0

             case (bufr_nets)
               vcordsf2=1.5

             case (bufr_nees,bufr_ness)
               vcordsf2=1.5

           end select

         case (2)
           select case(ilem)
             ! us,vs,ffs,dds
             case (bufr_neds,bufr_nefs,bufr_neus,bufr_nevs)
               vcordsf2=20.0

             case (bufr_nepn)
               vcordsf2=0.0

             case (bufr_neps)
               vcordsf2=0.0

             case (bufr_nets)
               vcordsf2=11.5

             case (bufr_nees,bufr_ness)
               vcordsf2=11.5
           end select

         case (3)
           select case(ilem)
             case (bufr_neds,bufr_nefs,bufr_neus,bufr_nevs)
               vcordsf2=10.0

             case (bufr_nepn)
               vcordsf2=0.0

             case (bufr_neps)
               vcordsf2=0.0

             case (bufr_nets)
               vcordsf2=1.5

             case (bufr_nees)
               vcordsf2=0.0

             case (bufr_ness)
               vcordsf2=1.5

           end select

         case (4)
           select case(ilem)
             case (bufr_neds,bufr_nefs,bufr_neus,bufr_nevs)
               vcordsf2=20.0

             case (bufr_nepn)
               vcordsf2=0.0

             case (bufr_neps)
               vcordsf2=0.0

             case (bufr_nets)
               vcordsf2=1.5

             case (bufr_nees)
               vcordsf2=0.0

             case (bufr_ness)
               vcordsf2=1.5

           end select

         case (5)
           select case(ilem)
             case (bufr_neds,bufr_nefs,bufr_neus,bufr_nevs)
               vcordsf2=10.0

           end select

      end select

      surfvcord = vcordsf2

  end function  surfvcord



  subroutine  adjust_sfvcoord(obsdat,start,end)
!
!**s/r ADJUST_SFVCOORD  - Computation of HEIGHT ASSIGNED TO SURFACE OBSERVATIONS
!
!
!Author  : P. Koclas *CMC/CMDA  April 2013
!Revision:
!          S. Macpherson *ARMA  Oct 2013
!              -- add GB-GPS (GP family) element BUFR_NEZD (ele 15031)
!              -- NOTE that for GP data, ELEV = GPS Antenna Height so
!                 no adjustment is needed (SFC_VCO=0).
!
!*    Purpose:  -Compute  HEIGHT ASSIGNED TO SURFACE OBSERVATIONS
!                and INSERT INTO CMA.
!
!
!Arguments
!               INPUT:
!                  -OBSDAT    : instance of obsspace_data module object
!                  -START     : FIRST OBERVATION
!                  -END       : LAST  OBERVATION
!
      implicit none
      integer  :: start,end
      integer  :: j,jo,rln,nlv
      integer  :: varno,codtyp,ity
      real     :: sfc_vco,elev
      real(obs_real) :: ppp
      type (struct_obs), intent(inout):: obsdat

      write(*,*)'   ADJUST_SFVCOORD '

      do jo = start, end
        rln = obs_headElem_i(obsdat,OBS_RLN,jo)
        nlv = obs_headElem_i(obsdat,OBS_NLV,jo)
        ity = obs_headElem_i(obsdat,OBS_ITY,jo)
        codtyp = ity
        elev = obs_headElem_r(obsdat,OBS_ALT,jo)
        do j = rln, nlv + rln -1

          varno = obs_bodyElem_i(obsdat,OBS_VNM,j)
          select case(varno)
            case(bufr_neds,bufr_nefs,bufr_neus,bufr_nevs,bufr_nets,bufr_ness,bufr_nepn,bufr_neps,bufr_nehs,bufr_nezd)

             sfc_vco= surfvcord(varno,codtyp)
             if ( varno /= bufr_nepn) then
                ppp = elev + sfc_vco
                call obs_bodySet_r(obsdat,OBS_PPP,j,ppp)
                call obs_bodySet_i(obsdat,OBS_VCO,j,1)
             else
                 ppp = 0.
                call obs_bodySet_r(obsdat,OBS_PPP,j,ppp)
                call obs_bodySet_i(obsdat,OBS_VCO,j,1)
             end if
          end select
        end do

      end do

      write(*,*)' DONE   ADJUST_SFVCOORD '

  end subroutine adjust_sfvcoord

  subroutine  adjust_hum_gz(obsdat,start,end)
!**s/r ADJUST_HUM_GZ  - Adjust  t-td and GZ in obsdat
!
!
!Author  : P. Koclas *CMC/CMDA  April 2013
!Revision:
!
!*    Purpose:  - Adjust  t-td values to zesmax=30. in obsdat
!                 set Z to GZ                       in obsdat
!
!
!Arguments
! 
!               INPUT:
!                  -OBSDAT    : instance of obsspace_data module object
!                  -START     : FIRST OBERVATION
!                  -END       : LAST  OBERVATION
!
!
      implicit none
      integer  :: start,end

      integer  :: j,jo,rln,nlv
      integer  :: varno
      type (struct_obs), intent(inout):: obsdat

      real(obs_real)    :: zesmax,gz,obsv
      real              :: rmin

      zesmax = 30.0

      write(*,*)'   ADJUST_HUM_GZ '

      do jo = start, end
        rln = obs_headElem_i(obsdat,OBS_RLN,JO)
        nlv = obs_headElem_i(obsdat,OBS_NLV,JO)

        do j = rln, nlv + rln -1

          varno = obs_bodyElem_i(obsdat,OBS_VNM,j)
          select case(varno)
            case(bufr_nees,bufr_ness)
             obsv = obs_bodyElem_r(obsdat,OBS_VAR,j)
             if ( obsv > zesmax) then
                obsv = zesmax
             end if
             call obs_bodySet_r(obsdat,OBS_VAR,j, obsv )
            case(bufr_negz)
             obsv = obs_bodyElem_r(obsdat,OBS_VAR,j)
             gz = obsv*grav
             call obs_bodySet_r(obsdat,OBS_VAR,j,gz )
          end select

        end do

      end do

      write(*,*)' DONE   ADJUST_HUM_GZ '

  end subroutine adjust_hum_gz

  subroutine  set_err_gbgps(obsdat,start,end)
!**s/r SET_ERR_GBGPS  - SET INITIAL ERROR FRO GROUND BASED GPS
!
!
!Author  : P. Koclas *CMC/CMDA  July 2013
!Revision:
!
!*    Purpose:  - PUT 15032 observation element as error of 15031 element  in obsdat
!
!
!Arguments
! 
!               INPUT:
!                  -OBSDAT    : instance of obsspace_data module object
!                  -START     : FIRST OBERVATION
!                  -END       : LAST  OBERVATION
!
      implicit none
      real(obs_real)    :: obsv
      integer  :: start,end

      integer  :: j,jo,rln,nlv
      integer  :: varno
      type (struct_obs), intent(inout):: obsdat

      write(*,*)'   SET_ERR_GBGPS '

      do jo = start, end
        rln = obs_headElem_i(obsdat,OBS_RLN,jo)
        nlv = obs_headElem_i(obsdat,OBS_NLV,jo)

        obsv = real(MPC_missingValue_R8,obs_real)
        do j = rln, nlv + rln -1

          varno = obs_bodyElem_i(obsdat,OBS_VNM,j)
          if ( varno == 15032 ) then
             obsv = obs_bodyElem_r(obsdat,OBS_VAR,j)
             call obs_bodySet_i(obsdat,OBS_VNM,j,999 )
             exit
          end if

        end do
        do j = rln, nlv + rln -1

          varno = obs_bodyElem_i(obsdat,OBS_VNM,j)
          if ( varno == 15031 .and. obsv /= real(MPC_missingValue_R8,obs_real)) then
             call obs_bodySet_r(obsdat,OBS_OER,j,obsv)
             exit
          end if

        end do

      end do

      write(*,*)' DONE   SET_ERR_GBGPS '

  end subroutine set_err_gbgps

  integer function cvt_obs_instrum(sensor)
    !
    ! func CVT_OBS_INSTRUM : Map burp satellite sensor indicator (element #2048) to
    !                         the corresponding burp satellite instrument (element
    !                         #2019).
    !
    ! Author  : J. Halle CMDA/SMC May 2002
    ! Revision:
    !           J.W. Blezius ARMA Feb 2011 - converted from subroutine to a function
    !
    !     Purpose:  Map burp satellite sensor indicator (element #2048) to
    !               burp satellite instrument (element #2019). This is a more
    !               complete common element, allowing for future expansion.
    !
    !               Table of BURP satellite sensor indicator element #002048
    !               --------------------------------------------------------
    !               Satellite sensor     BURP satellite sensor indicator
    !               ----------------     -------------------------------
    !               HIRS                 0
    !               MSU                  1
    !               SSU                  2
    !               AMSUA                3
    !               AMSUB                4
    !               AVHRR                5
    !               SSMI                 6
    !               NSCAT                7
    !               SEAWINDS             8
    !               Reserved             9-14
    !               Missing value        15
    !

    implicit none
    integer :: sensor      ! BURP satellite sensor indicator (element #2048)
    integer :: instrument  ! BURP satellite instrument       (element #2019)

    select case (sensor)
      case (000);   instrument=606  ! HIRS
      case (001);   instrument=623  ! MSU
      case (002);   instrument=627  ! SSU
      case (003);   instrument=570  ! AMSUA
      case (004);   instrument=574  ! AMSUB
      case (005);   instrument=591  ! AVHRR
      case (006);   instrument=905  ! SSMI
      case (007);   instrument=312  ! NSCAT
      case (008);   instrument=313  ! SEAWINDS
      case (015);   instrument=2047 ! Missing value
      case default; instrument=2047 ! Unrecognized value
    end select

    cvt_obs_instrum = instrument
    return

  end function cvt_obs_instrum


end module obsUtil_mod

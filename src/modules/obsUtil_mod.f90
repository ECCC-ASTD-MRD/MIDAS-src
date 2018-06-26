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
!! MODULE obsUtil_mod (prefix="obsu")
!!
!! *Purpose*: Common routines used by burpfiles_mod and sqlitefiles_mod
!!            
!!
!--------------------------------------------------------------------------

module obsUtil_mod
  
  use obsSpaceData_mod
  use bufr_mod
  use mathPhysConstants_mod
  use codePrecision_mod
  use earthConstants_mod, only:  grav
  use codtyp_mod

  implicit none
  save
  private
  public :: obsu_computeDirectionSpeedResiduals, obsu_setassflg, obsu_updateFlagWindDirectionSpeed
  public :: obsu_windDirectionToUV, obsu_adjustHumGZ, obsu_computeVertCoordSurfObs, obsu_setGbgpsError, obsu_cvt_obs_instrum
  public :: obsu_scaleFSO

  contains

  subroutine obsu_updateFlagWindDirectionSpeed(obsSpaceData)
    implicit none
    ! argument
    type(struct_obs) :: obsSpaceData
    ! locals
    integer          :: iuu, ivv, iff, idd, flagu, flagv, newflag
    integer          :: headerIndex, headerIndexStart, headerIndexEnd, jWindType, bodyIndex, bodyIndex2
    real(obs_real)   :: zlevu
    logical          :: llok
    character(len=9) :: stid

    WIND_TYPE: do jWindType = 1, 2
      if (jWindType == 1 ) then
        iuu=bufr_neuu
        ivv=bufr_nevv
        idd=bufr_nedd
        iff=bufr_neff
      else
        iuu=bufr_neus
        ivv=bufr_nevs
        idd=bufr_neds
        iff=bufr_nefs
      end if

      BODY: do bodyIndex = 1, obs_numBody(obsSpaceData)

        llok = ( obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) == iuu )
        flagu=-1
        if ( llok ) then
          headerIndex      = obs_bodyElem_i(obsSpaceData, OBS_HIND, bodyIndex   )
          headerIndexStart = obs_headElem_i(obsSpaceData, OBS_RLN , headerIndex )
          headerIndexEnd   = obs_headElem_i(obsSpaceData, OBS_NLV , headerIndex ) + headerIndexStart - 1
          stid             = obs_elem_c    (obsSpaceData, 'STID'  , headerIndex )
          zlevu            = obs_bodyElem_r(obsSpaceData, OBS_PPP , bodyIndex   )
          flagu            = obs_bodyElem_i(obsSpaceData, OBS_FLG , bodyIndex   ) !  GET FLAG OF U COMPONENT

          BODY2: do bodyIndex2 = headerIndexStart, headerIndexEnd

            if  ( obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2 ) == ivv .and. &
                  obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2 ) == zlevu ) then
              flagv = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex2 ) !  GET FLAG OF V COMPONENT
              newflag = ior(flagu, flagv)
            end if

          end do BODY2

          ! UPDATE FLAGS OF DIRECTION AN SPEED
          BODY2_2: do bodyIndex2 = headerIndexStart, headerIndexEnd

            if ( obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2 ) == idd .and. &
                 obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2 ) == zlevu ) then
              newflag =IOR(flagu,flagv)
              call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex2, newflag ) 
            end if

            if  ( obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex2) == iff   .and. &
                  obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex2) == zlevu ) then
              newflag =IOR(flagu,flagv)
              call obs_bodySet_i(obsSpaceData,OBS_FLG,bodyIndex2, newflag)
            end if

          end do BODY2_2

        end if ! llok

      end do BODY

    end do WIND_TYPE

  end subroutine obsu_updateFlagWindDirectionSpeed


  subroutine obsu_computeDirectionSpeedResiduals(elem_i,obsSpaceData)
    implicit none
    ! arguments
    type(struct_obs)    :: obsSpaceData
    integer, intent(in) :: elem_i
    ! locals
    integer :: iuu, ivv, iff, idd
    integer :: headerIndex, headerIndexStart, headerIndexEnd, jWindType, bodyIndex, bodyIndex2
    real(obs_real) :: zlevu, module, ang, uu, vv
    logical :: llok

    WIND_TYPE: do jWindType = 1, 2
      if (jWindType == 1) then
        iuu = bufr_neuu
        ivv = bufr_nevv
        idd = bufr_nedd
        iff = bufr_neff
      else
        iuu = bufr_neus
        ivv = bufr_nevs
        idd = bufr_neds
        iff = bufr_nefs
      end if

      ! Process all data within the domain of the model
      BODY: do bodyIndex = 1, obs_numBody(obsSpaceData)

        llok = ( obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex ) == 1 .and. &
                 obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex) == iuu )

        if ( llok ) then
          headerIndex      = obs_bodyElem_i(obsSpaceData, OBS_HIND, bodyIndex  )
          headerIndexStart = obs_headElem_i(obsSpaceData, OBS_RLN , headerIndex )
          headerIndexEnd   = obs_headElem_i(obsSpaceData, OBS_NLV , headerIndex ) + headerIndexStart - 1
          zlevu            = obs_bodyElem_r(obsSpaceData, OBS_PPP , bodyIndex )
          uu = - obs_bodyElem_r(obsSpaceData, elem_i, bodyIndex ) + obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
         
           BODY2: do bodyIndex2 = headerIndexStart, headerIndexEnd
     
             if  ( obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2 ) == ivv .and. &
                  obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2 ) == zlevu ) then
              vv = -obs_bodyElem_r(obsSpaceData, elem_i, bodyIndex2 ) + obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex2 )
              ! 1-calculate angle
              module = sqrt((uu**2)+(vv**2))

              if (module == 0.) then
                ang = 0.0d0
              else
                ang = atan2(vv,uu)
                ang = (270.0d0 - ang  * MPC_DEGREES_PER_RADIAN_R8 )
                ! 2-Change to meteorological definition of wind direction.
                if (ang > 360.0d0 ) ang = ang - 360.0d0
                if (ang <= 0.0d0  ) ang = ang + 360.0d0
              end if

            end if

          end do BODY2
          
          ! insert resduals into obsSpaceData
          BODY2_2: do bodyIndex2 = headerIndexStart, headerIndexEnd

            if  ( obs_bodyElem_i( obsSpaceData, OBS_VNM, bodyIndex2 ) == idd .and. &
                  obs_bodyElem_r( obsSpaceData, OBS_PPP, bodyIndex2 ) == zlevu ) then
              call obs_bodySet_r( obsSpaceData, elem_i, bodyIndex2, &
                                  obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex2 ) - ang )
              
              if ( obs_bodyElem_r( obsSpaceData,elem_i,bodyIndex2) >  180.0d0)  &
                call obs_bodySet_r( obsSpaceData, elem_i, bodyIndex2, &
                                    obs_bodyElem_r(obsSpaceData, elem_i, bodyIndex2 ) - 360.0d0 )
              if ( obs_bodyElem_r(obsSpaceData,elem_i,bodyIndex2) <= -180.0d0)  &
                call obs_bodySet_r( obsSpaceData, elem_i, bodyIndex2, &
                                    obs_bodyElem_r(obsSpaceData, elem_i, bodyIndex2 ) + 360.0d0 )
              call obs_bodySet_r( obsSpaceData, elem_i, bodyIndex2, &
                                  - 1.0d0 * obs_bodyElem_r(obsSpaceData,elem_i,bodyIndex2 ) )
              call obs_bodySet_r( obsSpaceData, OBS_OER, bodyIndex2, 1.0d0 )
              call obs_bodySet_i( obsSpaceData, OBS_ASS, bodyIndex2, 1 )
              call obs_bodySet_i( obsSpaceData, OBS_FLG, bodyIndex2, 0 )
            end if
            if  ( obs_bodyElem_i( obsSpaceData, OBS_VNM, bodyIndex2 ) == iff .and. &
                  obs_bodyElem_r( obsSpaceData, OBS_PPP, bodyIndex2 ) == zlevu ) then
              call obs_bodySet_r( obsSpaceData, elem_i,  bodyIndex2, &
                                  obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex2 ) - module )
              call obs_bodySet_r( obsSpaceData, OBS_OER, bodyIndex2, 1.0d0 )
              call obs_bodySet_i( obsSpaceData, OBS_ASS, bodyIndex2, 1)
              call obs_bodySet_i( obsSpaceData, OBS_FLG, bodyIndex2, 0)
            end if

          end do BODY2_2

        end if

      end do BODY

    end do WIND_TYPE

  end subroutine obsu_computeDirectionSpeedResiduals


  subroutine obsu_setassflg(obsSpaceData)
    ! Purpose:  Set banco quality control bit #12 for all data assimilated
    !           by current analysis.
    implicit none
    ! argument
    type(struct_obs) :: obsSpaceData
    ! local
    integer :: bodyIndex

    ! Process all data
    BODY: do bodyIndex = 1, obs_numBody(obsSpaceData)

      if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex ) == 1 ) &
        call obs_bodySet_i( obsSpaceData, OBS_FLG, bodyIndex, &
                            ibset( obs_bodyElem_i(obsSpaceData, OBS_FLG,bodyIndex), 12 ))
    end do BODY

  end subroutine obsu_setassflg


  subroutine obsu_windDirectionToUV(obsdat, headerIndexStart, headerIndexEnd, missingValue )
    implicit none
    ! arguments
    type (struct_obs), intent(inout) :: obsdat
    integer          , intent(in)    :: headerIndexStart, headerIndexEnd
    real             , intent(in)    :: missingValue
    ! locals
    integer        :: varno, varno2, varno4
    real           :: obsuv
    integer        :: headerIndex, bodyIndex, bodyIndexStart, bodyIndexEnd, bodyIndex2,jpos, ilem
    integer        :: ddflag, ffflag, newflag, uuflag, vvflag
    integer        :: ilemf, ilemu, ilemv, indu_misg, indv_misg, indum, indvm
    logical        :: llmisdd, llmisff, llmis, lluv_misg, llu_misg, llv_misg
    logical        :: lluv_present, llu_present, llv_present
    real(obs_real) :: uu, vv, dd, ff
    real(obs_real) :: level_dd, level4, level, level_uu
    character(len=*), parameter :: my_name = 'obsu_windDirectionToUV'
    character(len=*), parameter :: my_warning = '****** '// my_name //' WARNING: '
    character(len=*), parameter :: my_error   = '******** '// my_name //' ERROR: '

    ffflag = 0  ! bhe 
    write(*,'(a,2i8)')   my_name//': first and last observations: ', headerIndexStart, headerIndexEnd
    write(*,'(a,f12.4)') my_name//': missing value: ', missingValue
    write(*,*)'****************************************************' 
   
    HEADER1: do headerIndex = headerIndexStart, headerIndexEnd
      
      bodyIndexStart = obs_headElem_i(obsdat, OBS_RLN, headerIndex )
      bodyIndexEnd   = obs_headElem_i(obsdat, OBS_NLV, headerIndex ) + bodyIndexStart - 1
      
      BODY1: do bodyIndex = bodyIndexStart, bodyIndexEnd 

        dd = missingValue
        ff = missingValue
        varno = obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex)
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
      
        dd       = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex )
        ddflag   = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex )
        level_dd = obs_bodyElem_r(obsdat, OBS_PPP, bodyIndex )
        llu_misg = .false.
        llv_misg = .false.
        llu_present = .false.
        llv_present = .false.
        indum = -1
        indvm = -1
      
        ! FIND IF  U AND V ARE ALREADY IN CMA
        UVINOBSDAT: do bodyIndex2 = bodyIndex, bodyIndexEnd

          level4 = obs_bodyElem_r(obsdat, OBS_PPP, bodyIndex2 )

          if (level4 == level_dd) then
            varno4 = obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex2 )

            select case (varno4)
              case (bufr_neuu, bufr_nevv, bufr_neff, bufr_nedd, bufr_neus, bufr_nevs, bufr_neds, bufr_nefs )
                obsuv = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex2 )
                if ( varno4 == ilemu .and. obsuv /= missingValue ) then
                  llu_present = .true.
                  indum = bodyIndex2
                else if ( varno4 == ilemv .and. obsuv /= missingValue ) then
                  llv_present = .true.
                  indvm = bodyIndex2
                end if
                if ( varno4 == ilemu .and. obsuv == missingValue ) then
                  llu_misg = .true.
                  indu_misg = bodyIndex2
                else if ( varno4 == ilemv .and. obsuv == missingValue ) then
                  llv_misg = .true.
                  indv_misg = bodyIndex2
                end if
            end select

          end if

        end do UVINOBSDAT

        lluv_misg = (llu_misg .and. llv_misg)
        lluv_present = (llu_present .and. llv_present)

        if ( lluv_misg ) then

          CALCUV: do bodyIndex2 = bodyIndex, bodyIndexEnd
            llmisff = .true.
            llmisdd = .true.
            llmis   = .true.
            level = obs_bodyElem_r(obsdat, OBS_PPP, bodyIndex2 )

            if ( level /= level_dd) cycle CALCUV

            varno2 = obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex2 )

            if ( varno2 == ilemf ) then

              ff = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex2 )
              ffflag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex2 )

              if ( dd == 0.d0 .and. ff > 0. .or. dd > 360. .or. dd < 0. ) then
                llmisdd = .true.
                llmisff = .true.
              else if ( dd == missingValue .or. ff == missingValue ) then
                llmisdd = .true.
                llmisff = .true.
              else
                llmisdd = .false.
                llmisff = .false.
              end if

              ! IF SPEED = 0 CALM WIND IS ASSUMED.
              if (ff == 0.d0) dd = 0.d0
              dd = dd + 180.
              if ( dd > 360.) dd = dd - 360.
              dd = dd * mpc_radians_per_degree_r8
              ! U,V COMPONENTS ARE
              uu = ff * sin(dd)
              vv = ff * cos(dd)

              if ( llmisdd == .true. .or. llmisff == .true. ) then

                llmis = .true.
                if ( indu_misg > 0 .or. indv_misg > 0 ) then
                  call obs_bodySet_i(obsdat, OBS_VNM, indu_misg, -1 )
                  call obs_bodySet_i(obsdat, OBS_VNM, indv_misg, -1 )
                end if

              else

                llmis = .false.

              end if

            end if

            newflag = ior(ddflag, ffflag )

            if ( indum > 0 .or. indvm > 0 ) then
              call obs_bodySet_i(obsdat, OBS_VNM, indu_misg, -1 )
              call obs_bodySet_i(obsdat, OBS_VNM, indv_misg, -1 )
            end if

            if ( llmis ) then
              if ( indum > 0 .or. indvm > 0 ) then
                call obs_bodySet_i(obsdat, OBS_FLG, induM, newflag )
                call obs_bodySet_i(obsdat, OBS_FLG, indvM, newflag )
              end if
            else if ( .not. llmis ) then
              call obs_bodySet_r(obsdat, OBS_VAR, indu_misg, uu )
              call obs_bodySet_i(obsdat, OBS_FLG, indu_misg, newflag )
              call obs_bodySet_r(obsdat, OBS_VAR, indv_misg, vv )
              call obs_bodySet_i(obsdat, OBS_FLG, indv_misg, newflag )
            end if

          end do CALCUV

        else      
                 
          if ( lluv_present ) then
            call obs_bodySet_i(obsdat, OBS_VNM, indu_misg, -1 )
            call obs_bodySet_i(obsdat, OBS_VNM, indv_misg, -1 )
          else
            if ( indum > 0 ) then
              call obs_bodySet_i(obsdat, OBS_VNM, indum, -1 )
            end if
            if ( indvm > 0 ) then
              call obs_bodySet_i(obsdat, OBS_VNM, indvm, -1 )
            end if

          end if

        end if

      end do BODY1

    end do HEADER1

    HEADER: do headerIndex = headerIndexStart, headerIndexEnd

      bodyIndexStart = obs_headElem_i(obsdat, OBS_RLN, headerIndex )
      bodyIndexEnd   = obs_headElem_i(obsdat, OBS_NLV, headerIndex ) + bodyIndexStart - 1
      
      BODY: do bodyIndex = bodyIndexStart, bodyIndexEnd
        llmisdd = .true.
        varno = obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex )
        level = obs_bodyElem_r(obsdat, OBS_PPP, bodyIndex )

        select case (varno)
          case (bufr_neuu)
            ilem = bufr_nevv
          case (bufr_neus)
            ilem = bufr_nevs        
          case default
            cycle BODY
        end select

        Jpos = -1

        ! TRANSFER THE FLAG BITS  FROM ONE WIND COMPONENT TO THE OTHER
        BODY2: do bodyIndex2 = bodyIndexStart, bodyIndexEnd

          uu       = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex2 )
          level_uu = obs_bodyElem_r(obsdat, OBS_PPP, bodyIndex2 )
          Jpos = -1

          if ( level_uu == level .and. uu == missingValue ) call obs_bodySet_i(obsdat, OBS_VNM, bodyIndex2, -1 )

          if ( level_uu == level .and. uu /= missingValue ) then

            uuflag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex2 )
            varno2 = obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex2 )

            if ( ilem == varno2 ) then
              vvflag  = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex )
              newflag = ior( uuflag, vvflag )
              call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, newflag)
              call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex2, newflag)
              jpos = bodyIndex2
              exit BODY2
            end if

          end if

        end do BODY2

        ! ELIMINATE ENTRIES WHERE ONE COMPONENT OF WIND (UU OR VV) IS MISSING
        if (jpos < 0) then

          write(*,*) ' eliminate winds for station : ', obs_elem_c(obsdat,'STID',headerIndex ),  &
          obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex ), obs_bodyElem_r(obsdat, OBS_PPP, bodyIndex )
          call obs_bodySet_i(obsdat, OBS_VNM, bodyIndex, -1 )

        end if

      end do BODY

    end do HEADER

  end subroutine obsu_windDirectionToUV

  real function surfvcord(ilem, codtyp)
    implicit none
    ! arguments
    integer :: ilem,codtyp
    ! locals
    integer :: type
    real    :: vcordsf2

    vcordsf2 = 0.
    if ( codtyp == codtyp_get_codtyp('temppilot')      .or. codtyp == codtyp_get_codtyp('tempsynop')       .or. &
         codtyp == codtyp_get_codtyp('pilotsynop')     .or. codtyp == codtyp_get_codtyp('temppilotsynop')  .or. &
         codtyp == codtyp_get_codtyp('pilot')          .or. codtyp == codtyp_get_codtyp('pilotmobil')      .or. &
         codtyp == codtyp_get_codtyp('temp')           .or. codtyp == codtyp_get_codtyp('tempdrop')        .or. &
         codtyp == codtyp_get_codtyp('tempmobil')      .or. codtyp == codtyp_get_codtyp('temppilotmobil')  .or. &
         codtyp == codtyp_get_codtyp('tempsynopmobil') .or. codtyp == codtyp_get_codtyp('pilotsynopmobil') .or. &
         codtyp == codtyp_get_codtyp('temppilotsynopmobil') ) type = 3  ! UPPER AIR LAND
    
    if ( codtyp == codtyp_get_codtyp('temppilotship')  .or. codtyp == codtyp_get_codtyp('tempshipship')    .or. &
         codtyp == codtyp_get_codtyp('tempsshipship')  .or. codtyp == codtyp_get_codtyp('pilotshipship')   .or. &
         codtyp == codtyp_get_codtyp('pilotship')      .or. codtyp == codtyp_get_codtyp('tempship') ) type = 4 ! UPPER AIR SHIP
 
    if ( codtyp == codtyp_get_codtyp('synopnonauto')   .or. codtyp == codtyp_get_codtyp('synopmobil')      .or. &
         codtyp == codtyp_get_codtyp('asynopauto') ) type = 1 ! SYNOPS

    if ( codtyp == codtyp_get_codtyp('shipnonauto')    .or. codtyp == codtyp_get_codtyp('drifter')         .or. &
         codtyp == codtyp_get_codtyp('synoppatrol')    .or. codtyp == codtyp_get_codtyp('ashipauto') ) type = 2 ! SHIPS

    if ( codtyp == codtyp_get_codtyp('ascat') ) type = 5 ! SCATTEROMETER WINDS

    select case(type)
      case (1)
        select case(ilem)
          case (bufr_neds, bufr_nefs, bufr_neus, bufr_nevs)
            ! us,vs,ffs,dds
            vcordsf2 = 10.0
          case (bufr_nepn)
            vcordsf2 = 0.0
          case (bufr_neps)
            vcordsf2 = 0.0
          case (bufr_nets)
            vcordsf2 = 1.5
          case (bufr_nees, bufr_ness)
            vcordsf2 = 1.5
        end select
      case (2)
        select case(ilem)
          ! us,vs,ffs,dds
          case (bufr_neds, bufr_nefs, bufr_neus, bufr_nevs)
            vcordsf2 = 20.0
          case (bufr_nepn)
            vcordsf2 = 0.0
          case (bufr_neps)
            vcordsf2 = 0.0
          case (bufr_nets)
            vcordsf2 = 11.5
          case (bufr_nees, bufr_ness)
            vcordsf2 = 11.5
       end select
     case (3)
       select case(ilem)
         case (bufr_neds, bufr_nefs, bufr_neus, bufr_nevs)
           vcordsf2 = 10.0
         case (bufr_nepn)
           vcordsf2 = 0.0
         case (bufr_neps)
           vcordsf2 = 0.0
         case (bufr_nets)
           vcordsf2 = 1.5
         case (bufr_nees)
           vcordsf2 = 0.0
         case (bufr_ness)
           vcordsf2 = 1.5
       end select
     case (4)
       select case(ilem)
         case (bufr_neds, bufr_nefs, bufr_neus, bufr_nevs)
           vcordsf2 = 20.0
         case (bufr_nepn)
           vcordsf2 = 0.0
         case (bufr_neps)
           vcordsf2 = 0.0
         case (bufr_nets)
           vcordsf2 = 1.5
         case (bufr_nees)
           vcordsf2 = 0.0
         case (bufr_ness)
           vcordsf2 = 1.5
        end select
      case (5)

      select case(ilem)
        case (bufr_neds, bufr_nefs, bufr_neus, bufr_nevs)
          vcordsf2 = 10.0
      end select
    end select

    surfvcord = vcordsf2

  end function  surfvcord


  subroutine  obsu_computeVertCoordSurfObs(obsdat, headerIndexStart, headerIndexEnd )
    implicit none
    ! arguments
    type (struct_obs), intent(inout) :: obsdat
    integer          , intent(in)    :: headerIndexStart, headerIndexEnd
    ! locals
    integer  :: bodyIndex, headerIndex, bodyIndexStart, bodyIndexEnd
    integer  :: varno, codtyp
    real           :: sfc_vco, elev
    real(obs_real) :: ppp

    HEADER: do headerIndex = headerIndexStart, headerIndexEnd
        
      bodyIndexStart = obs_headElem_i(obsdat, OBS_RLN, headerIndex )
      bodyIndexEnd   = obs_headElem_i(obsdat, OBS_NLV, headerIndex ) + bodyIndexStart - 1
      codtyp = obs_headElem_i(obsdat, OBS_ITY, headerIndex )
      elev = obs_headElem_r(obsdat, OBS_ALT, headerIndex )

      BODY: do bodyIndex = bodyIndexStart, bodyIndexEnd
        varno = obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex )

        select case(varno)
          case(bufr_neds, bufr_nefs, bufr_neus, bufr_nevs, bufr_nets, bufr_ness, bufr_nepn, bufr_neps, bufr_nehs, bufr_nezd )
            sfc_vco = surfvcord(varno, codtyp )
            if ( varno /= bufr_nepn) then
              ppp = elev + sfc_vco
              call obs_bodySet_r(obsdat, OBS_PPP, bodyIndex, ppp )
              call obs_bodySet_i(obsdat, OBS_VCO, bodyIndex, 1 )
            else
              ppp = 0.
              call obs_bodySet_r(obsdat,OBS_PPP, bodyIndex, ppp )
              call obs_bodySet_i(obsdat,OBS_VCO, bodyIndex, 1 )
            end if
        end select

      end do BODY

    end do HEADER

  end subroutine obsu_computeVertCoordSurfObs


  subroutine  obsu_adjustHumGZ(obsdat, headerIndexStart, headerIndexEnd )
    implicit none
    ! arguments
    type (struct_obs), intent(inout) :: obsdat
    integer          , intent(in)    :: headerIndexStart, headerIndexEnd
    ! locals
    integer  :: bodyIndex, headerIndex, bodyIndexStart, bodyIndexEnd
    integer  :: varno
    real(obs_real)    :: zesmax,gz,obsv
    real              :: rmin

    zesmax = 30.0

    HEADER: do headerIndex = headerIndexStart, headerIndexEnd 

      bodyIndexStart = obs_headElem_i(obsdat, OBS_RLN, headerIndex )
      bodyIndexEnd   = obs_headElem_i(obsdat, OBS_NLV, headerIndex ) + bodyIndexStart - 1
      
      BODY: do bodyIndex = bodyIndexStart, bodyIndexEnd

        varno = obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex )

        select case(varno )
          case(bufr_nees, bufr_ness )
            obsv = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex )
            if ( obsv > zesmax ) obsv = zesmax
            call obs_bodySet_r(obsdat, OBS_VAR, bodyIndex, obsv )
          case(bufr_negz )
            obsv = obs_bodyElem_r(obsdat,OBS_VAR, bodyIndex )
            gz = obsv * grav
            call obs_bodySet_r(obsdat, OBS_VAR, bodyIndex, gz )
        end select

      end do BODY

    end do HEADER

  end subroutine obsu_adjustHumGZ


  subroutine  obsu_setGbgpsError(obsdat, headerIndexStart, headerIndexEnd )
    implicit none
    ! arguments
    type (struct_obs), intent(inout) :: obsdat
    integer          , intent(in)    :: headerIndexStart, headerIndexEnd
    ! locals
    real(obs_real)    :: obsv
    integer  :: bodyIndex, headerIndex,bodyIndexStart, bodyIndexEnd
    integer  :: varno

    HEADER: do headerIndex = headerIndexStart, headerIndexEnd
 
      bodyIndexStart = obs_headElem_i(obsdat, OBS_RLN, headerIndex )
      bodyIndexEnd   = obs_headElem_i(obsdat, OBS_NLV, headerIndex ) + bodyIndexStart - 1
      obsv = real(MPC_missingValue_R8, obs_real)
      
      BODY: do bodyIndex = bodyIndexStart, bodyIndexEnd 

        varno = obs_bodyElem_i(obsdat,OBS_VNM,bodyIndex )

        if ( varno == bufr_nefe ) then
          obsv = obs_bodyElem_r(obsdat,OBS_VAR, bodyIndex )
          call obs_bodySet_i(obsdat, OBS_VNM, bodyIndex, 999 )
          exit BODY
        end if

      end do BODY

      BODY2: do bodyIndex = bodyIndexStart, bodyIndexEnd

        varno = obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex )

        if ( varno == bufr_nezd .and. obsv /= real(MPC_missingValue_R8, obs_real)) then
          call obs_bodySet_r(obsdat, OBS_OER, bodyIndex, obsv )
          exit BODY2
        end if

      end do BODY2

    end do HEADER

  end subroutine obsu_setGbgpsError

  
  integer function obsu_cvt_obs_instrum(sensor)
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
    implicit none
    integer :: sensor      ! BURP satellite sensor indicator (element #2048)
    integer :: instrument  ! BURP satellite instrument       (element #2019)

    select case (sensor)
      case (000);   instrument = 606  ! HIRS
      case (001);   instrument = 623  ! MSU
      case (002);   instrument = 627  ! SSU
      case (003);   instrument = 570  ! AMSUA
      case (004);   instrument = 574  ! AMSUB
      case (005);   instrument = 591  ! AVHRR
      case (006);   instrument = 905  ! SSMI
      case (007);   instrument = 312  ! NSCAT
      case (008);   instrument = 313  ! SEAWINDS
      case (015);   instrument = 2047 ! Missing value
      case default; instrument = 2047 ! Unrecognized value
    end select

    obsu_cvt_obs_instrum = instrument

  end function obsu_cvt_obs_instrum

  subroutine obsu_scaleFSO(obsdat)
    implicit none

    type (struct_obs), intent(inout):: obsdat
    integer  :: headerIndex, bodyIndex, bodyIndexBeg, bodyIndexEnd
    real(8)  :: FSOVal

    do headerIndex = 1, obs_numHeader(obsdat)

      bodyIndexBeg = obs_headElem_i(obsdat,OBS_RLN,headerIndex)
      bodyIndexEnd = obs_headElem_i(obsdat,OBS_NLV,headerIndex) + bodyIndexBeg - 1

      do bodyIndex = bodyIndexBeg, bodyIndexEnd
        if ( obs_bodyElem_i(obsdat,OBS_ASS,bodyIndex) == 1 ) then
          FSOVal = obs_bodyElem_r(obsdat,OBS_FSO,bodyIndex)
          ! FSO value is quite small so giving scale factor 1.0e6 for storage in burp
          ! but treating GBGPS variable 15031 differently from the others
          if ( obs_bodyElem_i(obsdat,OBS_VNM,bodyIndex)  /= 15031 ) then
            call obs_bodySet_r(obsdat,OBS_FSO,bodyIndex, FSOVal*1.0e6)
          end if
         end if
      end do

    end do

  end subroutine obsu_scaleFSO

end module obsUtil_mod

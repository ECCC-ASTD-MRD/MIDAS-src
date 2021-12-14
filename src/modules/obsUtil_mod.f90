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

module obsUtil_mod
  ! MODULE obsUtil_mod (prefix='obsu' category='6. Observation input/output')
  !
  ! :Purpose: Common routines used by burpfiles_mod and sqlitefiles_mod
  !  
  
  use obsSpaceData_mod
  use bufr_mod
  use mathPhysConstants_mod
  use codePrecision_mod
  use codtyp_mod
  use obsVariableTransforms_mod
  use utilities_mod

  implicit none
  save
  private
  public :: obsu_setassflg, obsu_updateSourceVariablesFlag
  public :: obsu_computeVertCoordSurfObs, obsu_setGbgpsError, obsu_cvt_obs_instrum
  
contains

  !--------------------------------------------------------------------------
  ! obsu_updateSourceVariablesFlag
  !--------------------------------------------------------------------------
  subroutine obsu_updateSourceVariablesFlag(obsSpaceData)
    implicit none
    ! argument
    type(struct_obs) :: obsSpaceData

    ! locals
    integer          :: transformBufrCode, transformBufrCodeExtra
    integer          :: flag, flagExtra, mergedFlag
    integer          :: sourceBufrCode, sourceBufrCodeExtra
    integer          :: headerIndex, bodyIndexStart, bodyIndexEnd, bodyIndex, bodyIndex2
    real(pre_obsReal)   :: level

    body : do bodyIndex = 1, obs_numBody(obsSpaceData) 

      if (ovt_isTransformedVariable(obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex))) then

        ! We have found a transformed variable
        transformBufrCode = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)

        ! We will assign his assimilation flag to its source variable
        sourceBufrCode = ovt_getSourceBufrCode(transformBufrCode)

        headerIndex      = obs_bodyElem_i(obsSpaceData, OBS_HIND, bodyIndex   )
        bodyIndexStart   = obs_headElem_i(obsSpaceData, OBS_RLN , headerIndex )
        bodyIndexEnd     = obs_headElem_i(obsSpaceData, OBS_NLV , headerIndex ) + bodyIndexStart - 1
        level            = obs_bodyElem_r(obsSpaceData, OBS_PPP , bodyIndex   )
        flag             = obs_bodyElem_i(obsSpaceData, OBS_FLG , bodyIndex   )

        if (ovt_isWindObs(sourceBufrCode)) then

           ! Get the flag of the companion transformed variable
          transformBufrCodeExtra = ovt_getDestinationBufrCode(sourceBufrCode, &
                                                                   extra_opt=.true.)

          flagExtra = -999
          do bodyIndex2 = bodyIndexStart, bodyIndexEnd
            if ( obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2 ) == transformBufrCodeExtra .and. &
                 obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2 ) == level ) then
              flagExtra = obs_bodyElem_i(obsSpaceData, OBS_FLG , bodyIndex2)
            end if
          end do

          ! Combine flags
          if (flagExtra /= -999) then
            mergedFlag = ior(flag,flagExtra)
          else
            call utl_abort('obsu_updateSourceVariablesFlag: could not find the wind companion variable')
          end if
          
          ! Find sourceBufrCode and sourceBufrCodeExtra and update their flag
          sourceBufrCodeExtra = ovt_getSourceBufrCode(transformBufrCode,&
                                                                      extra_opt=.true.)
          do bodyIndex2 = bodyIndexStart, bodyIndexEnd
            if ( obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2 ) == sourceBufrCode .and. &
                 obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2 ) == level ) then
              call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex2, mergedFlag ) 
            end if
            if ( obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2 ) == sourceBufrCodeExtra .and. &
                 obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2 ) == level ) then
              call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex2, mergedFlag ) 
            end if
          end do
          
        else
          
          ! Find sourceBufrCode and update its flag
          do bodyIndex2 = bodyIndexStart, bodyIndexEnd
            if ( obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2 ) == sourceBufrCode .and. &
                 obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2 ) == level ) then
              call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex2, flag ) 
            end if
          end do

        end if
        
      end if
      
    end do body

  end subroutine obsu_updateSourceVariablesFlag

  !--------------------------------------------------------------------------
  ! obsu_setassflg
  !--------------------------------------------------------------------------
  subroutine obsu_setassflg(obsSpaceData)
    !
    ! :Purpose: Set banco quality control bit #12 for all data assimilated
    !           by current analysis.
    !
    implicit none
    ! argument
    type(struct_obs) :: obsSpaceData

    ! local
    integer :: bodyIndex

    ! Process all data
    BODY: do bodyIndex = 1, obs_numBody(obsSpaceData)

      if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex ) == obs_assimilated ) &
        call obs_bodySet_i( obsSpaceData, OBS_FLG, bodyIndex, &
                            ibset( obs_bodyElem_i(obsSpaceData, OBS_FLG,bodyIndex), 12 ))
    end do BODY

  end subroutine obsu_setassflg

  !--------------------------------------------------------------------------
  ! surfvcord
  !--------------------------------------------------------------------------
  real function surfvcord( varno, codtyp )

    implicit none
    ! Arguments
    integer    :: varno
    integer    :: codtyp
    ! locals
    character(len=9)  :: family

    family = codtypfam(codtyp)
    surfvcord = 0.0

    select case(family)
      case ('synop')
        select case(varno)
          case (bufr_neds, bufr_nefs, bufr_neus, bufr_nevs, bufr_gust)
            surfvcord = 10.0
          case (bufr_nets, bufr_nees, bufr_ness, bufr_vis, bufr_logVis)
            surfvcord = 1.5
        end select
      case ('ship')
        select case(varno)
          case (bufr_neds, bufr_nefs, bufr_neus, bufr_nevs, bufr_gust)
            surfvcord = 20.0
          case (bufr_nets, bufr_nees, bufr_ness, bufr_vis, bufr_logVis)
            surfvcord = 11.5
        end select
      case ('upairland')
        select case(varno)
          case (bufr_neds, bufr_nefs, bufr_neus, bufr_nevs, bufr_gust)
            surfvcord = 10.0
          case (bufr_nets, bufr_ness, bufr_vis, bufr_logVis)
            surfvcord = 1.5
        end select
      case ('upairship')
        select case(varno)
          case (bufr_neds, bufr_nefs, bufr_neus, bufr_nevs, bufr_gust)
            surfvcord = 20.0
          case (bufr_nets, bufr_ness, bufr_vis, bufr_logVis)
            surfvcord = 1.5
        end select
      case ('scatwinds')
        select case(varno)
          case (bufr_neds, bufr_nefs, bufr_neus, bufr_nevs)
            surfvcord = 10.0
        end select
    end select

  end function  surfvcord

  !--------------------------------------------------------------------------
  ! codtypfam
  !--------------------------------------------------------------------------
  function codtypfam(codtyp) result(family)
    implicit none
    integer          :: codtyp
    character(len=9) :: family

    if (       codtyp == codtyp_get_codtyp('synopnonauto')   .or. codtyp == codtyp_get_codtyp('synopmobil')      &
       .or.    codtyp == codtyp_get_codtyp('asynopauto') ) then
       family = 'synop'
    else if (  codtyp == codtyp_get_codtyp('shipnonauto')    .or. codtyp == codtyp_get_codtyp('drifter')         &
       .or.    codtyp == codtyp_get_codtyp('synoppatrol')    .or. codtyp == codtyp_get_codtyp('ashipauto') ) then
       family = 'ship'
    else if (  codtyp == codtyp_get_codtyp('temppilot')      .or. codtyp == codtyp_get_codtyp('tempsynop')       &
       .or.    codtyp == codtyp_get_codtyp('pilotsynop')     .or. codtyp == codtyp_get_codtyp('temppilotsynop')  &
       .or.    codtyp == codtyp_get_codtyp('pilot')          .or. codtyp == codtyp_get_codtyp('pilotmobil')      &
       .or.    codtyp == codtyp_get_codtyp('temp')           .or. codtyp == codtyp_get_codtyp('tempdrop')        &
       .or.    codtyp == codtyp_get_codtyp('tempmobil')      .or. codtyp == codtyp_get_codtyp('temppilotmobil')  &
       .or.    codtyp == codtyp_get_codtyp('tempsynopmobil') .or. codtyp == codtyp_get_codtyp('pilotsynopmobil') &
       .or.    codtyp == codtyp_get_codtyp('temppilotsynopmobil') ) then
       family = 'upairland'
    else if (  codtyp == codtyp_get_codtyp('temppilotship')  .or. codtyp ==codtyp_get_codtyp('tempshipship')     &
       .or.    codtyp == codtyp_get_codtyp('tempsshipship')  .or. codtyp ==codtyp_get_codtyp('pilotshipship')    &
       .or.    codtyp == codtyp_get_codtyp('pilotship')      .or. codtyp ==codtyp_get_codtyp('tempship') ) then
       family = 'upairship'
    else if (  codtyp == codtyp_get_codtyp('ascat') ) then
       family = 'scatwinds'
    else
       family = 'other'
    end if

  end function codtypfam

  !--------------------------------------------------------------------------
  ! obsu_computeVertCoordSurfObs
  !--------------------------------------------------------------------------
  subroutine  obsu_computeVertCoordSurfObs(obsdat, headerIndexStart, headerIndexEnd )
    !
    implicit none

    ! Arguments
    type (struct_obs), intent(inout) :: obsdat
    integer          , intent(in)    :: headerIndexStart
    integer          , intent(in)    :: headerIndexEnd

    ! locals
    integer  :: bodyIndex, headerIndex, bodyIndexStart, bodyIndexEnd
    integer  :: varno, codtyp
    real           :: sfc_vco, elev
    real(pre_obsReal) :: ppp

    HEADER: do headerIndex = headerIndexStart, headerIndexEnd
        
      bodyIndexStart = obs_headElem_i(obsdat, OBS_RLN, headerIndex )
      bodyIndexEnd   = obs_headElem_i(obsdat, OBS_NLV, headerIndex ) + bodyIndexStart - 1
      codtyp = obs_headElem_i(obsdat, OBS_ITY, headerIndex )
      elev = obs_headElem_r(obsdat, OBS_ALT, headerIndex )

      BODY: do bodyIndex = bodyIndexStart, bodyIndexEnd
        varno = obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex )

        select case(varno)
          case(bufr_neds, bufr_nefs, bufr_neus, bufr_nevs, bufr_nets, bufr_ness, bufr_nepn, bufr_neps, bufr_nehs, &
               bufr_nezd, bufr_vis, bufr_logVis, bufr_gust )
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

  !--------------------------------------------------------------------------
  ! obsu_setGbgpsError
  !--------------------------------------------------------------------------
  subroutine  obsu_setGbgpsError(obsdat, headerIndexStart, headerIndexEnd )
    !
    implicit none

    ! Arguments
    type (struct_obs), intent(inout) :: obsdat
    integer          , intent(in)    :: headerIndexStart
    integer          , intent(in)    :: headerIndexEnd

    ! locals
    real(pre_obsReal)    :: obsv
    integer  :: bodyIndex, headerIndex,bodyIndexStart, bodyIndexEnd
    integer  :: varno

    HEADER: do headerIndex = headerIndexStart, headerIndexEnd
 
      bodyIndexStart = obs_headElem_i(obsdat, OBS_RLN, headerIndex )
      bodyIndexEnd   = obs_headElem_i(obsdat, OBS_NLV, headerIndex ) + bodyIndexStart - 1
      obsv = real(MPC_missingValue_R8, pre_obsReal)
      
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

        if ( varno == bufr_nezd .and. obsv /= real(MPC_missingValue_R8, pre_obsReal)) then
          call obs_bodySet_r(obsdat, OBS_OER, bodyIndex, obsv )
          exit BODY2
        end if

      end do BODY2

    end do HEADER

  end subroutine obsu_setGbgpsError

  !--------------------------------------------------------------------------
  ! obsu_cvt_obs_instrum
  !--------------------------------------------------------------------------
  integer function obsu_cvt_obs_instrum(sensor)
    !
    ! :Purpose: Map burp satellite sensor indicator (element #2048) to
    !           burp satellite instrument (element #2019). This is a more
    !           complete common element, allowing for future expansion.
    !
    ! :Table of BURP satellite sensor indicator element #002048:
    ! ==================  =============================== 
    ! Satellite sensor    BURP satellite sensor indicator
    ! ==================  =============================== 
    !               HIRS                 0
    !                MSU                 1
    !                SSU                 2
    !              AMSUA                 3
    !              AMSUB                 4
    !              AVHRR                 5
    !               SSMI                 6
    !              NSCAT                 7
    !           SEAWINDS                 8
    !           Reserved                9-14
    !      Missing value                15
    ! ==================  =============================== 
    !
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

end module obsUtil_mod

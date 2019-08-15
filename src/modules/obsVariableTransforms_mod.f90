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

module obsVariableTransforms_mod
  ! MODULE obsVariableTransforms_mod (prefix='ovt' category='3. High-level transformations')
  !
  ! :Purpose: To store various functions for variable transforms using inputs
  !           from obsSpaceData. Outputs are also placed ObsSpaceData.
  !  
  
  use obsSpaceData_mod
  use bufr_mod
  use codePrecision_mod
  use mathPhysConstants_mod
  use earthConstants_mod, only: grav
  use codtyp_mod
  use utilities_mod
  use obsFilter_mod

  implicit none
  save
  private

  public :: ovt_setup, ovt_transformObsValues, ovt_transformResiduals
  public :: ovt_getTransformVariableBufrCode, ovt_getSourceVariableBufrCode, ovt_variableBufrCodeSkipped
  public :: ovt_isWindObs, ovt_isTransformedVariable, ovt_adjustHumGZ

  integer  , parameter   :: nTransformSupported = 2
  integer  , parameter   :: nVariableBufrCodeMax= 2

  type :: struct_ovt
    character(len=48) :: transformName                 (nTransformSupported)
    integer           :: nVariableBufrCode             (nTransformSupported)
    integer           :: sourceVariableBufrCode        (nTransformSupported,nVariableBufrCodeMax)
    integer           :: sourceVariableBufrCodeExtra   (nTransformSupported,nVariableBufrCodeMax)
    integer           :: transformVariableBufrCode     (nTransformSupported,nVariableBufrCodeMax)
    integer           :: transformVariableBufrCodeExtra(nTransformSupported,nVariableBufrCodeMax)
    logical           :: wind                          (nTransformSupported) = .false.
    logical           :: activeTransform               (nTransformSupported) = .false.
    integer           :: nSkippedVariableBufrCodes    = 0
    integer           :: skippedVariableBufrCodes(50) = -999
    logical           :: setup = .false.
  end type struct_ovt

  ! Module's internal data
  type(struct_ovt) :: ovt

contains

  !--------------------------------------------------------------------------
  ! ovt_initStructure
  !--------------------------------------------------------------------------
  subroutine ovt_initStructure
    implicit none

    integer :: transformIndex

    ! Upper and surface winds
    transformIndex = 1
    ovt%transformName                 (transformIndex)   = 'windSpeedDirectionToUV'
    ovt%nVariableBufrCode             (transformIndex)   = 2
    ovt%sourceVariableBufrCode        (transformIndex,:) = (/bufr_nedd, bufr_neds/) ! direction
    ovt%sourceVariableBufrCodeExtra   (transformIndex,:) = (/bufr_neff, bufr_nefs/) ! speed
    ovt%transformVariableBufrCode     (transformIndex,:) = (/bufr_neuu, bufr_neus/) ! u-wind
    ovt%transformVariableBufrCodeExtra(transformIndex,:) = (/bufr_nevv, bufr_nevs/) ! v-wind
    ovt%wind                          (transformIndex)   = .true.

    ! log of visibility
    transformIndex = 2
    ovt%transformName            (transformIndex)   = 'visToLogVis'
    ovt%nVariableBufrCode        (transformIndex)   = 1
    ovt%sourceVariableBufrCode   (transformIndex,:) = (/bufr_vis   , -1/) ! visibility
    ovt%transformVariableBufrCode(transformIndex,:) = (/bufr_logVis, -1/) ! log(visibility)
    ovt%wind                     (transformIndex)   = .false.

    ! Common
    ovt%nSkippedVariableBufrCodes   = 2
    ovt%skippedVariableBufrCodes(1) = bufr_neff  ! because we still want to be able ...
    ovt%skippedVariableBufrCodes(2) = bufr_nefs  ! ... to assimilate wind speed obs alone

  end subroutine ovt_initStructure

  !--------------------------------------------------------------------------
  ! ovt_setup
  !--------------------------------------------------------------------------
  subroutine ovt_setup(variableBufrCodeReaded)
    implicit none

    integer, intent(in) :: variableBufrCodeReaded(:)

    integer, allocatable :: variableBufrCodeAssimilated(:)

    integer :: nVariableBufrCodeReaded, nVariableBufrCodeAssimilated
    integer :: readVariableBufrCodeIndex, transformIndex, assimVariableBufrCodeIndex

    logical :: variableTransformNeeded, foundTransformation
    logical, save :: firstTime = .true.

    if (firstTime) then
      call ovt_initStructure
      firstTime = .false.
    end if

    nVariableBufrCodeAssimilated = filt_nVariableBufrCodeAssimilated()
    nVariableBufrCodeReaded = size(variableBufrCodeReaded)
    variableTransformNeeded = .false.

    allocate(variableBufrCodeAssimilated(nVariableBufrCodeAssimilated))
    call filt_getVariableBufrCodeAssimilated(variableBufrCodeAssimilated)

    do readVariableBufrCodeIndex = 1, nVariableBufrCodeReaded
    
      ! Check if a transform is neeeded
      if (filt_variableBufrCodeAssimilated(variableBufrCodeReaded(readVariableBufrCodeIndex)) .or. &
          variableBufrCodeReaded(readVariableBufrCodeIndex) == bufr_neff             .or. &
          variableBufrCodeReaded(readVariableBufrCodeIndex) == bufr_nefs ) then
        cycle ! No transformation needed. Move on.
      end if

      ! Find and activate the appropriate transform
      foundTransformation = .false.

      transformLoop : do transformIndex = 1, nTransformSupported
        assimVariableBufrCodeLoop : do assimVariableBufrCodeIndex = 1, nVariableBufrCodeAssimilated
          if (any(ovt%sourceVariableBufrCode     (transformIndex,:) == variableBufrCodeReaded(readVariableBufrCodeIndex))        .and. &
              any(ovt%transformVariableBufrCode(transformIndex,:) == variableBufrCodeAssimilated(assimVariableBufrCodeIndex)) )  then
            if (.not. ovt%activeTransform(transformIndex)) then
              write(*,*) 'ovt_setup: transform activated : ', trim(ovt%transformName(transformIndex))
              ovt%activeTransform(transformIndex) = .true.
            end if
            foundTransformation = .true.
            exit transformLoop
          end if
        end do assimVariableBufrCodeLoop
      end do transformLoop

      if (.not. foundTransformation) then
        write(*,*)
        write(*,*) 'ovt_setup: !WARNING! No transform found for the readed bufr code = ', variableBufrCodeReaded(readVariableBufrCodeIndex)
        write(*,*) '           We are assuming that this observation is readed but not assimilated.'
        write(*,*) '           Please consider removing this bufr code from the readed observation list.'
        ovt%nSkippedVariableBufrCodes = ovt%nSkippedVariableBufrCodes + 1
        ovt%skippedVariableBufrCodes(ovt%nSkippedVariableBufrCodes) = variableBufrCodeReaded(readVariableBufrCodeIndex)
      end if

    end do

    ovt%setup = .true.

    deallocate(variableBufrCodeAssimilated)

  end subroutine ovt_setup

  !--------------------------------------------------------------------------
  ! ovt_variableBufrCodeSkipped
  !--------------------------------------------------------------------------
  function ovt_variableBufrCodeSkipped(sourceVariableBufrCode) result(skip)
    implicit none

    integer, intent(in) :: sourceVariableBufrCode
    logical :: skip
    
    integer :: transformIndex, variableBufrCodeIndex

    skip = .false.
    do variableBufrCodeIndex = 1, ovt%nSkippedVariableBufrCodes
      if (ovt%skippedVariableBufrCodes(variableBufrCodeIndex) == sourceVariableBufrCode) then
        skip = .true.
        exit
      end if
    end do

  end function ovt_variableBufrCodeSkipped

  !--------------------------------------------------------------------------
  ! ovt_getTransformVariableBufrCode
  !--------------------------------------------------------------------------
  function ovt_getTransformVariableBufrCode(sourceVariableBufrCode,extra_opt) result(transformVariableBufrCode)
    implicit none

    integer, intent(in) :: sourceVariableBufrCode
    integer :: transformVariableBufrCode
    logical, optional :: extra_opt
    
    logical :: extra
    integer :: transformIndex, variableBufrCodeIndex

    if (present(extra_opt)) then
      extra = extra_opt
    else
      extra = .false. ! default
    end if

    transformVariableBufrCode = -1

    transformLoop : do transformIndex = 1, nTransformSupported
      variableBufrCodeLoop : do variableBufrCodeIndex = 1, ovt%nVariableBufrCode(transformIndex)
        if (ovt%sourceVariableBufrCode(transformIndex,variableBufrCodeIndex) == sourceVariableBufrCode) then
          if (extra) then
            transformVariableBufrCode = ovt%transformVariableBufrCodeExtra(transformIndex,variableBufrCodeIndex)
          else
            transformVariableBufrCode = ovt%transformVariableBufrCode(transformIndex,variableBufrCodeIndex)
          end if
          exit transformLoop
        end if
      end do variableBufrCodeLoop
    end do transformLoop

    if (transformVariableBufrCode == -1) then
      write(*,*)
      write(*,*) 'ovt_getTransformVariableBufrCode: source variableBufrCode = ', sourceVariableBufrCode
      call utl_abort('ovt_getTransformVariableBufrCode: found no associated variableBufrCode for the above source variable bufr code')
    end if

  end function ovt_getTransformVariableBufrCode

  !--------------------------------------------------------------------------
  ! ovt_getSourceVariableBufrCode
  !--------------------------------------------------------------------------
  function ovt_getSourceVariableBufrCode(transformVariableBufrCode,extra_opt) result(sourceVariableBufrCode)
    implicit none

    integer, intent(in) :: transformVariableBufrCode
    integer :: sourceVariableBufrCode
    logical, optional :: extra_opt
    
    logical :: extra
    integer :: transformIndex, transformVariableBufrCodeIndex

    if (present(extra_opt)) then
      extra = extra_opt
    else
      extra = .false. ! default
    end if

    sourceVariableBufrCode = -1

    transformLoop : do transformIndex = 1, nTransformSupported
      transformVariableBufrCodeLoop : do transformVariableBufrCodeIndex = 1, ovt%nVariableBufrCode(transformIndex)
        if (ovt%transformVariableBufrCode(transformIndex,transformVariableBufrCodeIndex) == transformVariableBufrCode) then
          if (extra) then
            sourceVariableBufrCode = ovt%sourceVariableBufrCodeExtra(transformIndex,transformVariableBufrCodeIndex)
          else
            sourceVariableBufrCode = ovt%sourceVariableBufrCode(transformIndex,transformVariableBufrCodeIndex)
          end if
          exit transformLoop
        end if
      end do transformVariableBufrCodeLoop
    end do transformLoop

    if (sourceVariableBufrCode == -1) then
      write(*,*)
      write(*,*) 'ovt_getSourceVariableBufrCode: tranform variable bufr code = ', transformVariableBufrCode
      call utl_abort('ovt_getSourceVariableBufrCode: found no associated variable bufr code for the above transform variable bufr code')
    end if

  end function ovt_getSourceVariableBufrCode

  !--------------------------------------------------------------------------
  ! ovt_isWindObs
  !--------------------------------------------------------------------------
  function ovt_isWindObs(sourceVariableBufrCode) result(wind)
    implicit none

    integer, intent(in) :: sourceVariableBufrCode
    logical :: wind

    integer :: transformIndex, variableBufrCodeIndex

    transformLoop : do transformIndex = 1, nTransformSupported
      variableBufrCodeLoop : do variableBufrCodeIndex = 1, ovt%nVariableBufrCode(transformIndex)
        if (ovt%sourceVariableBufrCode(transformIndex,variableBufrCodeIndex) == sourceVariableBufrCode) then
          wind = ovt%wind(transformIndex)
          exit transformLoop
        end if
      end do variableBufrCodeLoop
    end do transformLoop

  end function ovt_isWindObs

  !--------------------------------------------------------------------------
  ! ovt_isTransformedVariable
  !--------------------------------------------------------------------------
  function ovt_isTransformedVariable(variableBufrCode) result(transformed)
    implicit none

    integer, intent(in) :: variableBufrCode
    logical :: transformed

    integer :: transformIndex, variableBufrCodeIndex

    transformed = .false.

    transformLoop : do transformIndex = 1, nTransformSupported
      if (ovt%activeTransform(transformIndex)) then
        variableBufrCodeLoop : do variableBufrCodeIndex = 1, ovt%nVariableBufrCode(transformIndex)
          if (ovt%transformVariableBufrCode(transformIndex,variableBufrCodeIndex) == variableBufrCode) then
            transformed = .true.
            exit transformLoop
          end if
        end do variableBufrCodeLoop
      end if
    end do transformLoop

  end function ovt_isTransformedVariable

  !--------------------------------------------------------------------------
  ! ovt_transformObsValues
  !--------------------------------------------------------------------------
  subroutine ovt_transformObsValues(obsSpaceData, headerIndexStart, headerIndexEnd)
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsSpaceData
    integer          , intent(in)    :: headerIndexStart 
    integer          , intent(in)    :: headerIndexEnd

    integer :: transformIndex

    if (.not. ovt%setup) then
      call utl_abort('ovt_transformObsValues: this module has not been setup')
    end if

    do transformIndex = 1, nTransformSupported

      if (ovt%activeTransform(transformIndex)) then
        select case(trim(ovt%transformName(transformIndex)))
        case ('windSpeedDirectionToUV')
          call ovt_windSpeedDirectionToUV(obsSpaceData, headerIndexStart, headerIndexEnd)
        case ('visToLogVis')
          call ovt_visToLogVis           (obsSpaceData, headerIndexStart, headerIndexEnd)
        case default
          call utl_abort('ovt_transformObsValues: Unsupported function ' // trim(ovt%transformName(transformIndex)))
        end select
      end if

    end do

  end subroutine ovt_transformObsValues

  !--------------------------------------------------------------------------
  ! ovt_TransformResiduals
  !--------------------------------------------------------------------------
  subroutine ovt_transformResiduals(obsSpaceData, residualTypeID)
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsSpaceData
    integer          , intent(in)    :: residualTypeID

    integer :: transformIndex

    if (.not. ovt%setup) then
      call utl_abort('ovt_transformResiduals: this module has not been setup')
    end if

    do transformIndex = 1, nTransformSupported

      if (ovt%activeTransform(transformIndex)) then
        select case(trim(ovt%transformName(transformIndex)))
        case ('windSpeedDirectionToUV')
            call ovt_UVtoWindSpeedDirection_residual(obsSpaceData, residualTypeID)
          case ('visToLogVis') 
            call ovt_visToLogVis_residual           (obsSpaceData, residualTypeID)
        case default
          call utl_abort('ovt_transformResiduals: Unsupported function ' // trim(ovt%transformName(transformIndex)))
        end select
      end if

    end do
    
  end subroutine ovt_transformResiduals

  !--------------------------------------------------------------------------
  ! ovt_windSpeedDirectionToUV
  !--------------------------------------------------------------------------
  subroutine ovt_windSpeedDirectionToUV(obsSpaceData, headerIndexStart, headerIndexEnd)
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsSpaceData
    integer          , intent(in)    :: headerIndexStart 
    integer          , intent(in)    :: headerIndexEnd

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

    ffflag = 0

    HEADER1: do headerIndex = headerIndexStart, headerIndexEnd
      
      bodyIndexStart = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex )
      bodyIndexEnd   = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex ) + bodyIndexStart - 1
      
      BODY1: do bodyIndex = bodyIndexStart, bodyIndexEnd 

        dd = obs_missingValue
        ff = obs_missingValue
        varno = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex)
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
      
        dd       = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
        ddflag   = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex )
        level_dd = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex )
        llu_misg = .false.
        llv_misg = .false.
        llu_present = .false.
        llv_present = .false.
        indum = -1
        indvm = -1

        ! FIND IF  U AND V ARE ALREADY IN CMA
        UVINOBSSPACEDATA: do bodyIndex2 = bodyIndex, bodyIndexEnd

          level4 = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2 )

          if (level4 == level_dd) then
            varno4 = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2 )
            
            select case (varno4)
            case (bufr_neuu, bufr_nevv, bufr_neff, bufr_nedd, bufr_neus, bufr_nevs, bufr_neds, bufr_nefs )
              obsuv = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex2 )
              if ( varno4 == ilemu .and. obsuv /= obs_missingValue ) then
                llu_present = .true.
                indum = bodyIndex2
              else if ( varno4 == ilemv .and. obsuv /= obs_missingValue ) then
                llv_present = .true.
                indvm = bodyIndex2
              end if
              if ( varno4 == ilemu .and. obsuv == obs_missingValue ) then
                llu_misg = .true.
                indu_misg = bodyIndex2
              else if ( varno4 == ilemv .and. obsuv == obs_missingValue ) then
                llv_misg = .true.
                indv_misg = bodyIndex2
              end if
            end select

          end if

        end do UVINOBSSPACEDATA

        lluv_misg = (llu_misg .and. llv_misg)
        lluv_present = (llu_present .and. llv_present)

        if ( lluv_misg ) then

          CALCUV: do bodyIndex2 = bodyIndex, bodyIndexEnd
            llmisff = .true.
            llmisdd = .true.
            llmis   = .true.
            level = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2 )

            if ( level /= level_dd) cycle CALCUV

            varno2 = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2 )

            if ( varno2 == ilemf ) then

              ff = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex2 )
              ffflag = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex2 )

              if ( dd == 0.d0 .and. ff > 0. .or. dd > 360. .or. dd < 0. ) then
                llmisdd = .true.
                llmisff = .true.
              else if ( dd == obs_missingValue .or. ff == obs_missingValue ) then
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
                  call obs_bodySet_i(obsSpaceData, OBS_VNM, indu_misg, -1 )
                  call obs_bodySet_i(obsSpaceData, OBS_VNM, indv_misg, -1 )
                end if
              else
                llmis = .false.
              end if

            end if

            newflag = ior(ddflag, ffflag )

            if ( indum > 0 .or. indvm > 0 ) then
              call obs_bodySet_i(obsSpaceData, OBS_VNM, indu_misg, -1 )
              call obs_bodySet_i(obsSpaceData, OBS_VNM, indv_misg, -1 )
            end if

            if ( llmis ) then
              if ( indum > 0 .or. indvm > 0 ) then
                call obs_bodySet_i(obsSpaceData, OBS_FLG, induM, newflag )
                call obs_bodySet_i(obsSpaceData, OBS_FLG, indvM, newflag )
              end if
            else if ( .not. llmis ) then
              call obs_bodySet_r(obsSpaceData, OBS_VAR, indu_misg, uu )
              call obs_bodySet_i(obsSpaceData, OBS_FLG, indu_misg, newflag )
              call obs_bodySet_r(obsSpaceData, OBS_VAR, indv_misg, vv )
              call obs_bodySet_i(obsSpaceData, OBS_FLG, indv_misg, newflag )
            end if

          end do CALCUV

        else      
                 
          if ( lluv_present ) then
            call obs_bodySet_i(obsSpaceData, OBS_VNM, indu_misg, -1 )
            call obs_bodySet_i(obsSpaceData, OBS_VNM, indv_misg, -1 )
          else
            if ( indum > 0 ) then
              call obs_bodySet_i(obsSpaceData, OBS_VNM, indum, -1 )
            end if
            if ( indvm > 0 ) then
              call obs_bodySet_i(obsSpaceData, OBS_VNM, indvm, -1 )
            end if

          end if

        end if

      end do BODY1

    end do HEADER1

    HEADER: do headerIndex = headerIndexStart, headerIndexEnd

      bodyIndexStart = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex )
      bodyIndexEnd   = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex ) + bodyIndexStart - 1
      
      BODY: do bodyIndex = bodyIndexStart, bodyIndexEnd
        llmisdd = .true.
        varno = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex )
        level = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex )

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

          uu       = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex2 )
          level_uu = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2 )
          Jpos = -1

          if ( level_uu == level .and. uu == obs_missingValue ) call obs_bodySet_i(obsSpaceData, OBS_VNM, bodyIndex2, -1 )

          if ( level_uu == level .and. uu /= obs_missingValue ) then

            uuflag = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex2 )
            varno2 = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2 )

            if ( ilem == varno2 ) then
              vvflag  = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex )
              newflag = ior( uuflag, vvflag )
              call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex, newflag)
              call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex2, newflag)
              jpos = bodyIndex2
              exit BODY2
            end if

          end if

        end do BODY2

        ! ELIMINATE ENTRIES WHERE ONE COMPONENT OF WIND (UU OR VV) IS MISSING
        if (jpos < 0) then

          write(*,*) ' eliminate winds for station : ', obs_elem_c(obsSpaceData,'STID',headerIndex ),  &
          obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex ), obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex )
          call obs_bodySet_i(obsSpaceData, OBS_VNM, bodyIndex, -1 )

        end if

      end do BODY

    end do HEADER

  end subroutine ovt_windSpeedDirectionToUV

  !--------------------------------------------------------------------------
  ! ovt_UVtoWindSpeedDirection_residual
  !--------------------------------------------------------------------------
  subroutine ovt_UVtoWindSpeedDirection_residual(obsSpaceData, elem_i)
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

        llok = ( obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex ) == obs_assimilated .and. &
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
                                    obs_bodyElem_r(obsSpaceData, elem_i, bodyIndex2 ) - real(360.0d0,OBS_REAL) )
              if ( obs_bodyElem_r(obsSpaceData,elem_i,bodyIndex2) <= -180.0d0)  &
                call obs_bodySet_r( obsSpaceData, elem_i, bodyIndex2, &
                                    obs_bodyElem_r(obsSpaceData, elem_i, bodyIndex2 ) + real(360.0d0,OBS_REAL) )
              call obs_bodySet_r( obsSpaceData, elem_i, bodyIndex2, &
                                  - real(1.0d0,OBS_REAL) * obs_bodyElem_r(obsSpaceData,elem_i,bodyIndex2 ) )
              call obs_bodySet_r( obsSpaceData, OBS_OER, bodyIndex2, real(1.0d0,OBS_REAL) )
              call obs_bodySet_i( obsSpaceData, OBS_ASS, bodyIndex2, obs_assimilated )
              call obs_bodySet_i( obsSpaceData, OBS_FLG, bodyIndex2, 0 )
            end if
            if  ( obs_bodyElem_i( obsSpaceData, OBS_VNM, bodyIndex2 ) == iff .and. &
                  obs_bodyElem_r( obsSpaceData, OBS_PPP, bodyIndex2 ) == zlevu ) then
              call obs_bodySet_r( obsSpaceData, elem_i,  bodyIndex2, &
                                  obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex2 ) - module )
              call obs_bodySet_r( obsSpaceData, OBS_OER, bodyIndex2, real(1.0d0,OBS_REAL) )
              call obs_bodySet_i( obsSpaceData, OBS_ASS, bodyIndex2, obs_assimilated)
              call obs_bodySet_i( obsSpaceData, OBS_FLG, bodyIndex2, 0)
            end if

          end do BODY2_2

        end if

      end do BODY

    end do WIND_TYPE

  end subroutine ovt_UVtoWindSpeedDirection_residual

  !--------------------------------------------------------------------------
  ! ovt_visToLogVis
  !--------------------------------------------------------------------------
  subroutine ovt_visToLogVis(obsSpaceData, headerIndexStart, headerIndexEnd)
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsSpaceData
    integer          , intent(in)    :: headerIndexStart 
    integer          , intent(in)    :: headerIndexEnd

    ! locals
    integer        :: headerIndex, bodyIndex, bodyIndexStart, bodyIndexEnd, bodyIndex2
    integer        :: visFlag, logVisFlag
    real(obs_real) :: visObs, visLevel, logVisObs, level

    ! Loop through headers
    header: do headerIndex = headerIndexStart, headerIndexEnd
      
      bodyIndexStart = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex )
      bodyIndexEnd   = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex ) + bodyIndexStart - 1
      
      ! Find each visibily report
      body: do bodyIndex = bodyIndexStart, bodyIndexEnd 

        if (obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex) /= bufr_vis) cycle body

        visObs   = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
        visFlag  = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex )
        visLevel = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex )

        ! Find the associated logVis body created earlier in the proper reading routine
        body2: do bodyIndex2 = bodyIndex, bodyIndexEnd

          level = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2 )

          if ( level /= visLevel   ) cycle body2

          if (obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2 ) /= bufr_logVis) cycle body2

          if (visObs == obs_missingValue) then
            logVisObs = visObs
          else
            ! vis -> log(vis)
            logVisObs = log(max(min(visObs,MPC_MAXIMUM_VIS_R8),MPC_MINIMUM_VIS_R8))
          end if
          call obs_bodySet_r(obsSpaceData, OBS_VAR, bodyIndex2, logVisObs)

          logVisFlag  = visFlag
          call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex2, logVisFlag )

          exit body2

        end do body2

      end do body

    end do header

  end subroutine ovt_visToLogVis

  !--------------------------------------------------------------------------
  ! ovt_visToLogVis_redidual
  !--------------------------------------------------------------------------
  subroutine ovt_visToLogVis_residual(obsSpaceData, residualTypeID)
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsSpaceData
    integer          , intent(in)    :: residualTypeID

    ! locals
    integer        :: headerIndex, bodyIndex, bodyIndexStart, bodyIndexEnd, bodyIndex2
    real(obs_real) :: visObs, logVisLevel, logVisObs, level
    real(obs_real) :: visResidual, logVisResidual

    ! Find each log of visibily assimilated observations
    body: do bodyIndex = 1, obs_numBody(obsSpaceData)
      
      if ( obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) /= obs_assimilated .or. &
           obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex) /= bufr_logVis ) cycle

      headerIndex    = obs_bodyElem_i(obsSpaceData, OBS_HIND, bodyIndex  )
      bodyIndexStart = obs_headElem_i(obsSpaceData, OBS_RLN , headerIndex )
      bodyIndexEnd   = obs_headElem_i(obsSpaceData, OBS_NLV , headerIndex ) + bodyIndexStart - 1
      logVisLevel    = obs_bodyElem_r(obsSpaceData, OBS_PPP , bodyIndex )

      logVisObs      = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
      logVisResidual = obs_bodyElem_r(obsSpaceData, residualTypeID, bodyIndex )

      ! Find the associated vis body
      body2: do bodyIndex2 = bodyIndexStart, bodyIndexEnd

        level = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2 )

        if (level /= logVisLevel) cycle body2
        if (obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2) /= bufr_vis) cycle body2

        ! log(vis) -> vis
        visObs      = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex2 )
        visResidual = visObs - exp(logVisObs-logVisResidual) ! o-p or o-a
        call obs_bodySet_r(obsSpaceData, residualTypeID, bodyIndex2, visResidual)

        ! Set the obs error to missing
        call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex2, obs_missingValue)

        ! Set flags
        call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex2, obs_assimilated)
        call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex2, 0)

        exit body2

      end do body2
      
    end do body

  end subroutine ovt_visToLogVis_residual

  !--------------------------------------------------------------------------
  ! ovt_adjustHumGZ
  !--------------------------------------------------------------------------
  subroutine  ovt_adjustHumGZ(obsSpaceData, headerIndexStart, headerIndexEnd )
    implicit none

    ! arguments
    type (struct_obs), intent(inout) :: obsSpaceData
    integer          , intent(in)    :: headerIndexStart
    integer          , intent(in)    :: headerIndexEnd

    ! locals
    integer  :: bodyIndex, headerIndex, bodyIndexStart, bodyIndexEnd
    integer  :: varno
    real(obs_real)    :: ESmax,gz,obsValue

    ESmax = 30.0

    HEADER: do headerIndex = headerIndexStart, headerIndexEnd 

      bodyIndexStart = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex )
      bodyIndexEnd   = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex ) + bodyIndexStart - 1
      
      BODY: do bodyIndex = bodyIndexStart, bodyIndexEnd

        varno = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex )

        select case(varno)
        case(bufr_nees, bufr_ness)
          obsValue = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
          if ( obsValue > ESmax ) obsValue = ESmax
          call obs_bodySet_r(obsSpaceData, OBS_VAR, bodyIndex, obsValue )
        case(bufr_negz)
          obsValue = obs_bodyElem_r(obsSpaceData,OBS_VAR, bodyIndex )
          gz = obsValue * grav
          call obs_bodySet_r(obsSpaceData, OBS_VAR, bodyIndex, gz )
        end select

      end do BODY

    end do HEADER

  end subroutine ovt_adjustHumGZ
  
end module obsVariableTransforms_mod

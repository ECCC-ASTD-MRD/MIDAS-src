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
  public :: ovt_getDestinationVariableBufrCode, ovt_getSourceVariableBufrCode, ovt_variableBufrCodeSkipped
  public :: ovt_isWindObs, ovt_isTransformedVariable, ovt_adjustHumGZ

  integer, parameter :: nTransformSupported  = 2

  type :: struct_transform
    character(len=48)    :: name
    integer              :: nBufrCode
    integer, allocatable :: sourceBufrCode        (:)
    integer, allocatable :: sourceBufrCodeExtra   (:)
    integer, allocatable :: destinationBufrCode     (:)
    integer, allocatable :: destinationBufrCodeExtra(:)
    logical              :: wind   = .false.
    logical              :: active = .false.
  end type struct_transform

  type(struct_transform) :: transform(nTransformSupported)

  integer                :: nSkippedBufrCodes    = 0
  integer                :: skippedBufrCodes(50) = -999

  logical                :: initialized                = .false.

contains

  !--------------------------------------------------------------------------
  ! ovt_initStructure
  !--------------------------------------------------------------------------
  subroutine ovt_initStructure
    !
    ! :Purpose: To set the transforms handled by this module
    ! 
    implicit none

    integer :: transformIndex

    ! Upper and surface winds
    transformIndex = 1
    transform(transformIndex)%name = 'windSpeedDirectionToUV'
    transform(transformIndex)%nBufrCode = 2
    allocate(transform(transformIndex)%sourceBufrCode          (transform(transformIndex)%nBufrCode))
    allocate(transform(transformIndex)%sourceBufrCodeExtra     (transform(transformIndex)%nBufrCode))
    allocate(transform(transformIndex)%destinationBufrCode     (transform(transformIndex)%nBufrCode))
    allocate(transform(transformIndex)%destinationBufrCodeExtra(transform(transformIndex)%nBufrCode))
    transform(transformIndex)%sourceBufrCode          (:) = (/bufr_nedd, bufr_neds/) ! direction
    transform(transformIndex)%sourceBufrCodeExtra     (:) = (/bufr_neff, bufr_nefs/) ! speed
    transform(transformIndex)%destinationBufrCode     (:) = (/bufr_neuu, bufr_neus/) ! u-wind
    transform(transformIndex)%destinationBufrCodeExtra(:) = (/bufr_nevv, bufr_nevs/) ! v-wind
    transform(transformIndex)%wind = .true.

    ! log of visibility
    transformIndex = 2
    transform(transformIndex)%name = 'visToLogVis'
    transform(transformIndex)%nBufrCode = 1
    allocate(transform(transformIndex)%sourceBufrCode     (transform(transformIndex)%nBufrCode))
    allocate(transform(transformIndex)%destinationBufrCode(transform(transformIndex)%nBufrCode))
    transform(transformIndex)%sourceBufrCode     (:) = (/bufr_vis   /) ! visibility
    transform(transformIndex)%destinationBufrCode(:) = (/bufr_logVis/) ! log(visibility)
    transform(transformIndex)%wind = .false.

    ! Skipped variables
    nSkippedBufrCodes   = 2
    skippedBufrCodes(1) = bufr_neff  ! because we still want to be able ...
    skippedBufrCodes(2) = bufr_nefs  ! ... to assimilate wind speed obs alone

  end subroutine ovt_initStructure

  !--------------------------------------------------------------------------
  ! ovt_setup
  !--------------------------------------------------------------------------
  subroutine ovt_setup(bufrCodeReaded)
    !
    ! :Purpose: To determine which transform must be actived
    !
    implicit none

    integer, intent(in) :: bufrCodeReaded(:)

    integer, allocatable :: bufrCodeAssimilated(:)

    integer :: nBufrCodeReaded, nBufrCodeAssimilated
    integer :: readBufrCodeIndex, transformIndex, assimBufrCodeIndex

    logical :: variableTransformNeeded, foundTransformation
    logical, save :: firstTime = .true.

    if (firstTime) then
      call ovt_initStructure
      firstTime = .false.
    end if

    nBufrCodeAssimilated = filt_nVariableBufrCodeAssimilated()
    nBufrCodeReaded = size(bufrCodeReaded)
    variableTransformNeeded = .false.

    allocate(bufrCodeAssimilated(nBufrCodeAssimilated))
    call filt_getVariableBufrCodeAssimilated(bufrCodeAssimilated)

    do readBufrCodeIndex = 1, nBufrCodeReaded
    
      ! Check if a transform is neeeded
      if (filt_variableBufrCodeAssimilated(bufrCodeReaded(readBufrCodeIndex)) .or. &
          bufrCodeReaded(readBufrCodeIndex) == bufr_neff             .or. &
          bufrCodeReaded(readBufrCodeIndex) == bufr_nefs ) then
        cycle ! No transformation needed. Move on.
      end if

      ! Find and activate the appropriate transform
      foundTransformation = .false.

      transformLoop : do transformIndex = 1, nTransformSupported
        assimBufrCodeLoop : do assimBufrCodeIndex = 1, nBufrCodeAssimilated
          if (any(transform(transformIndex)%sourceBufrCode(:) == bufrCodeReaded(readBufrCodeIndex))        .and. &
              any(transform(transformIndex)%destinationBufrCode(:) == bufrCodeAssimilated(assimBufrCodeIndex)) )  then
            if (.not. transform(transformIndex)%active) then
              write(*,*) 'ovt_setup: transform activated : ', trim(transform(transformIndex)%name)
              transform(transformIndex)%active = .true.
            end if
            foundTransformation = .true.
            exit transformLoop
          end if
        end do assimBufrCodeLoop
      end do transformLoop

      if (.not. foundTransformation) then
        write(*,*)
        write(*,*) 'ovt_setup: !WARNING! No transform found for the readed bufr code = ', bufrCodeReaded(readBufrCodeIndex)
        write(*,*) '           We are assuming that this observation is readed but not assimilated.'
        write(*,*) '           Please consider removing this bufr code from the readed observation list.'
        nSkippedBufrCodes = nSkippedBufrCodes + 1
        skippedBufrCodes(nSkippedBufrCodes) = bufrCodeReaded(readBufrCodeIndex)
      end if

    end do

    initialized = .true.

    deallocate(bufrCodeAssimilated)

  end subroutine ovt_setup

  !--------------------------------------------------------------------------
  ! ovt_variableBufrCodeSkipped
  !--------------------------------------------------------------------------
  function ovt_variableBufrCodeSkipped(sourceBufrCode) result(skip)
    !
    ! :Purpose: To determine if an observation must be ignored or not
    !
    implicit none

    integer, intent(in) :: sourceBufrCode
    logical :: skip
    
    integer :: transformIndex, bufrCodeIndex

    if (.not. initialized) then
      call utl_abort(' ovt_variableBufrCodeSkipped: this module has not been setup')
    end if

    skip = .false.
    do bufrCodeIndex = 1, nSkippedBufrCodes
      if (skippedBufrCodes(bufrCodeIndex) == sourceBufrCode) then
        skip = .true.
        exit
      end if
    end do

  end function ovt_variableBufrCodeSkipped

  !--------------------------------------------------------------------------
  ! ovt_getDestinationVariableBufrCode
  !--------------------------------------------------------------------------
  function ovt_getDestinationVariableBufrCode(sourceBufrCode,extra_opt) result(destinationBufrCode)
    !
    ! :Purpose: To get the bufr code of the transformed/destination variable based on the 
    !           bufr code of the source variable
    !
    implicit none

    integer, intent(in) :: sourceBufrCode
    integer :: destinationBufrCode
    logical, optional :: extra_opt
    
    logical :: extra
    integer :: transformIndex, bufrCodeIndex

    if (.not. initialized) then
      call utl_abort(' ovt_getDestinationVariableBufrCode: this module has not been setup')
    end if

    if (present(extra_opt)) then
      extra = extra_opt
    else
      extra = .false. ! default
    end if

    destinationBufrCode = -1

    transformLoop : do transformIndex = 1, nTransformSupported
      bufrCodeLoop : do bufrCodeIndex = 1, transform(transformIndex)%nBufrCode
        if (transform(transformIndex)%sourceBufrCode(bufrCodeIndex) == sourceBufrCode) then
          if (extra) then
            destinationBufrCode = transform(transformIndex)%destinationBufrCodeExtra(bufrCodeIndex)
          else
            destinationBufrCode = transform(transformIndex)%destinationBufrCode(bufrCodeIndex)
          end if
          exit transformLoop
        end if
      end do bufrCodeLoop
    end do transformLoop

    if (destinationBufrCode == -1) then
      write(*,*)
      write(*,*) 'ovt_getDestinationVariableBufrCode: source bufrCode = ', sourceBufrCode
      call utl_abort('ovt_getDestinationVariableBufrCode: found no associated bufrCode for the above source variable bufr code')
    end if

  end function ovt_getDestinationVariableBufrCode

  !--------------------------------------------------------------------------
  ! ovt_getSourceVariableBufrCode
  !--------------------------------------------------------------------------
  function ovt_getSourceVariableBufrCode(destinationBufrCode,extra_opt) result(sourceBufrCode)
    !
    ! :Purpose: To get the bufr code of the source variable based on the 
    !           bufr code of the destination/transformed variable
    !
    implicit none

    integer, intent(in) :: destinationBufrCode
    integer :: sourceBufrCode
    logical, optional :: extra_opt
    
    logical :: extra
    integer :: transformIndex, destinationBufrCodeIndex

    if (.not. initialized) then
      call utl_abort(' ovt_getSourceVariableBufrCode: this module has not been setup')
    end if

    if (present(extra_opt)) then
      extra = extra_opt
    else
      extra = .false. ! default
    end if

    sourceBufrCode = -1

    transformLoop : do transformIndex = 1, nTransformSupported
      destinationBufrCodeLoop : do destinationBufrCodeIndex = 1, transform(transformIndex)%nBufrCode
        if (transform(transformIndex)%destinationBufrCode(destinationBufrCodeIndex) == destinationBufrCode) then
          if (extra) then
            sourceBufrCode = transform(transformIndex)%sourceBufrCodeExtra(destinationBufrCodeIndex)
          else
            sourceBufrCode = transform(transformIndex)%sourceBufrCode(destinationBufrCodeIndex)
          end if
          exit transformLoop
        end if
      end do destinationBufrCodeLoop
    end do transformLoop

    if (sourceBufrCode == -1) then
      write(*,*)
      write(*,*) 'ovt_getSourceVariableBufrCode: tranform variable bufr code = ', destinationBufrCode
      call utl_abort('ovt_getSourceVariableBufrCode: found no associated variable bufr code for the above transform variable bufr code')
    end if

  end function ovt_getSourceVariableBufrCode

  !--------------------------------------------------------------------------
  ! ovt_isWindObs
  !--------------------------------------------------------------------------
  function ovt_isWindObs(sourceBufrCode) result(wind)
    !
    ! :Purpose: To determine if a bufr code is wind related
    !
    implicit none

    integer, intent(in) :: sourceBufrCode
    logical :: wind

    integer :: transformIndex, bufrCodeIndex

    if (.not. initialized) then
      call utl_abort(' ovt_isWindObs: this module has not been setup')
    end if

    transformLoop : do transformIndex = 1, nTransformSupported
      bufrCodeLoop : do bufrCodeIndex = 1, transform(transformIndex)%nBufrCode
        if (transform(transformIndex)%sourceBufrCode(bufrCodeIndex) == sourceBufrCode) then
          wind = transform(transformIndex)%wind
          exit transformLoop
        end if
      end do bufrCodeLoop
    end do transformLoop

  end function ovt_isWindObs

  !--------------------------------------------------------------------------
  ! ovt_isTransformedVariable
  !--------------------------------------------------------------------------
  function ovt_isTransformedVariable(bufrCode) result(transformed)
    !
    ! :Purpose: To determine if a bufr code is a transfomed variabled
    !
    implicit none

    integer, intent(in) :: bufrCode
    logical :: transformed

    integer :: transformIndex, bufrCodeIndex

    if (.not. initialized) then
      call utl_abort(' ovt_isTransformedVariable: this module has not been setup')
    end if

    transformed = .false.

    transformLoop : do transformIndex = 1, nTransformSupported
      if (transform(transformIndex)%active) then
        bufrCodeLoop : do bufrCodeIndex = 1, transform(transformIndex)%nBufrCode
          if (transform(transformIndex)%destinationBufrCode(bufrCodeIndex) == bufrCode) then
            transformed = .true.
            exit transformLoop
          end if
        end do bufrCodeLoop
      end if
    end do transformLoop

  end function ovt_isTransformedVariable

  !--------------------------------------------------------------------------
  ! ovt_transformObsValues
  !--------------------------------------------------------------------------
  subroutine ovt_transformObsValues(obsSpaceData, headerIndexStart, headerIndexEnd)
    !
    ! :Purpose: To perform observation variable transforms
    !
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsSpaceData
    integer          , intent(in)    :: headerIndexStart 
    integer          , intent(in)    :: headerIndexEnd

    integer :: transformIndex

    if (.not. initialized) then
      call utl_abort('ovt_transformObsValues: this module has not been setup')
    end if

    do transformIndex = 1, nTransformSupported

      if (transform(transformIndex)%active) then
        select case(trim(transform(transformIndex)%name))
        case ('windSpeedDirectionToUV')
          call ovt_windSpeedDirectionToUV(obsSpaceData, headerIndexStart, headerIndexEnd)
        case ('visToLogVis')
          call ovt_visToLogVis           (obsSpaceData, headerIndexStart, headerIndexEnd)
        case default
          call utl_abort('ovt_transformObsValues: Unsupported function ' // trim(transform(transformIndex)%name))
        end select
      end if

    end do

  end subroutine ovt_transformObsValues

  !--------------------------------------------------------------------------
  ! ovt_TransformResiduals
  !--------------------------------------------------------------------------
  subroutine ovt_transformResiduals(obsSpaceData, residualTypeID)
    !
    ! :Purpose: To compute the o-p or o-a of the source variable(s)
    !
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsSpaceData
    integer          , intent(in)    :: residualTypeID

    integer :: transformIndex

    if (.not. initialized) then
      call utl_abort('ovt_transformResiduals: this module has not been setup')
    end if

    do transformIndex = 1, nTransformSupported

      if (transform(transformIndex)%active) then
        select case(trim(transform(transformIndex)%name))
        case ('windSpeedDirectionToUV')
            call ovt_UVtoWindSpeedDirection_residual(obsSpaceData, residualTypeID)
          case ('visToLogVis') 
            call ovt_visToLogVis_residual           (obsSpaceData, residualTypeID)
        case default
          call utl_abort('ovt_transformResiduals: Unsupported function ' // trim(transform(transformIndex)%name))
        end select
      end if

    end do
    
  end subroutine ovt_transformResiduals

  !--------------------------------------------------------------------------
  ! ovt_windSpeedDirectionToUV
  !--------------------------------------------------------------------------
  subroutine ovt_windSpeedDirectionToUV(obsSpaceData, headerIndexStart, headerIndexEnd)
    !
    ! :Purpose: To transform wind observation in terms of speed and direction to u and v
    !
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsSpaceData
    integer          , intent(in)    :: headerIndexStart 
    integer          , intent(in)    :: headerIndexEnd

    ! locals
    integer        :: bufrCode, bufrCode2, bufrCode4
    real(obs_real) :: obsValue
    integer        :: headerIndex, bodyIndex, bodyIndexStart, bodyIndexEnd, bodyIndex2,bodyIndexFound, bufrCodeAssociated
    integer        :: directionFlag, speedFlag, combinedFlag, uWindFlag, vWindFlag
    integer        :: bufrCodeSpeed, bufrCodeUWind, bufrCodeVWind, indu_missing, indv_missing, indum, indvm
    logical        :: direction_missing, speed_missing, missing, uv_missing, uWind_missing, vWind_missing
    !logical        :: uv_present, uWind_present, vWind_present
    real(obs_real) :: uWind, vWind, direction, speed
    real(obs_real) :: level_direction, level4, level, level_uWind

    speedFlag = 0

    header1: do headerIndex = headerIndexStart, headerIndexEnd
      
      bodyIndexStart = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex )
      bodyIndexEnd   = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex ) + bodyIndexStart - 1
      
      body1: do bodyIndex = bodyIndexStart, bodyIndexEnd 

        direction = obs_missingValue
        speed     = obs_missingValue
        bufrCode  = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex)
        direction_missing = .true.
      
        if ( bufrCode /= bufr_nedd .and. bufrCode /= bufr_neds ) cycle body1

        if( bufrCode == bufr_neds) then
          ! Surface obs
          bufrCodeSpeed = bufr_nefs
          bufrCodeUWind = bufr_neus
          bufrCodeVWind = bufr_nevs
        else
          ! Upper air obs
          bufrCodeSpeed = bufr_neff
          bufrCodeUWind = bufr_neuu
          bufrCodeVWind = bufr_nevv
        end if
      
        direction       = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
        directionFlag   = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex )
        level_direction = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex )
        uWind_missing = .true.
        vWind_missing = .true.
        !uWind_present = .false.
        !vWind_present = .false.
        indum = -1
        indvm = -1

        ! check if u and v are present in obsSpaceData
        do bodyIndex2 = bodyIndex, bodyIndexEnd

          level4 = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2 )

          if (level4 == level_direction) then
            bufrCode4 = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2 )
            
            select case (bufrCode4)
            case (bufr_neuu, bufr_nevv, bufr_neff, bufr_nedd, bufr_neus, bufr_nevs, bufr_neds, bufr_nefs )
              obsValue = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex2 )
              if ( bufrCode4 == bufrCodeUWind .and. obsValue /= obs_missingValue ) then
                call utl_abort('ovt_windSpeedDirectionToUV: u-wind is aleady present. This should not happen!')
                !uWind_present = .true.
                uWind_missing = .false.
                indum = bodyIndex2
              else if ( bufrCode4 == bufrCodeVWind .and. obsValue /= obs_missingValue ) then
                call utl_abort('ovt_windSpeedDirectionToUV: v-wind is aleady present. This should not happen!')
                !vWind_present = .true.
                vWind_missing = .false.
                indvm = bodyIndex2
              end if
              if ( bufrCode4 == bufrCodeUWind .and. obsValue == obs_missingValue ) then
                !uWind_missing = .true.
                indu_missing = bodyIndex2
              else if ( bufrCode4 == bufrCodeVWind .and. obsValue == obs_missingValue ) then
                !vWind_missing = .true.
                indv_missing = bodyIndex2
              end if
            end select

          end if

        end do

        uv_missing = (uWind_missing .and. vWind_missing)
        !uv_present = (uWind_present .and. vWind_present)

        if ( uv_missing ) then

          calcuv: do bodyIndex2 = bodyIndex, bodyIndexEnd
            speed_missing = .true.
            direction_missing = .true.
            missing   = .true.
            level = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2 )

            if ( level /= level_direction) cycle calcuv

            bufrCode2 = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2 )

            if ( bufrCode2 == bufrCodeSpeed ) then

              speed = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex2 )
              speedFlag = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex2 )

              if ( direction == 0.d0 .and. speed > 0. .or. direction > 360. .or. direction < 0. ) then
                direction_missing = .true.
                speed_missing = .true.
              else if ( direction == obs_missingValue .or. speed == obs_missingValue ) then
                direction_missing = .true.
                speed_missing = .true.
              else
                direction_missing = .false.
                speed_missing = .false.
              end if

              if (speed == 0.d0) direction = 0.d0
              direction = direction + 180.
              if ( direction > 360.) direction = direction - 360.
              direction = direction * mpc_radians_per_degree_r8

              ! (speed,direction) -> (u,v)
              uWind = speed * sin(direction)
              vWind = speed * cos(direction)

              if ( direction_missing == .true. .or. speed_missing == .true. ) then
                missing = .true.
                if ( indu_missing > 0 .or. indv_missing > 0 ) then
                  call obs_bodySet_i(obsSpaceData, OBS_VNM, indu_missing, -1 )
                  call obs_bodySet_i(obsSpaceData, OBS_VNM, indv_missing, -1 )
                end if
              else
                missing = .false.
              end if

            end if

            combinedFlag = ior(directionFlag, speedFlag )

            if ( indum > 0 .or. indvm > 0 ) then
              call obs_bodySet_i(obsSpaceData, OBS_VNM, indu_missing, -1 )
              call obs_bodySet_i(obsSpaceData, OBS_VNM, indv_missing, -1 )
            end if

            if ( missing ) then
              if ( indum > 0 .or. indvm > 0 ) then
                call obs_bodySet_i(obsSpaceData, OBS_FLG, induM, combinedFlag )
                call obs_bodySet_i(obsSpaceData, OBS_FLG, indvM, combinedFlag )
              end if
            else if ( .not. missing ) then
              call obs_bodySet_r(obsSpaceData, OBS_VAR, indu_missing, uWind )
              call obs_bodySet_i(obsSpaceData, OBS_FLG, indu_missing, combinedFlag )
              call obs_bodySet_r(obsSpaceData, OBS_VAR, indv_missing, vWind )
              call obs_bodySet_i(obsSpaceData, OBS_FLG, indv_missing, combinedFlag )
            end if

          end do calcuv

        else      
 
          call utl_abort('ovt_windSpeedDirectionToUV: u-wind and/or v-wind are aleady present. This should not happen!')

        end if

      end do body1

    end do header1

    header: do headerIndex = headerIndexStart, headerIndexEnd

      bodyIndexStart = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex )
      bodyIndexEnd   = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex ) + bodyIndexStart - 1
      
      body: do bodyIndex = bodyIndexStart, bodyIndexEnd
        direction_missing = .true.
        bufrCode = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex )
        level = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex )

        select case (bufrCode)
          case (bufr_neuu)
            bufrCodeAssociated = bufr_nevv
          case (bufr_neus)
            bufrCodeAssociated = bufr_nevs        
          case default
            cycle body
        end select

        bodyIndexFound = -1

        ! TRANSFER THE FLAG BITS  FROM ONE WIND COMPONENT TO THE OTHER
        body2: do bodyIndex2 = bodyIndexStart, bodyIndexEnd

          uWind       = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex2 )
          level_uWind = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2 )
          bodyIndexFound = -1

          if ( level_uWind == level .and. uWind == obs_missingValue ) call obs_bodySet_i(obsSpaceData, OBS_VNM, bodyIndex2, -1 )

          if ( level_uWind == level .and. uWind /= obs_missingValue ) then

            uWindFlag = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex2 )
            bufrCode2 = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2 )

            if ( bufrCodeAssociated == bufrCode2 ) then
              vWindFlag  = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex )
              combinedFlag = ior( uWindFlag, vWindFlag )
              call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex, combinedFlag)
              call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex2, combinedFlag)
              bodyIndexFound = bodyIndex2
              exit body2
            end if

          end if

        end do body2

        ! ELIMINATE ENTRIES WHERE ONE COMPONENT OF WIND (UWIND OR VWIND) IS MISSING
        if (bodyIndexFound < 0) then

          write(*,*) ' ovt_windSpeedDirectionToUV: eliminate winds for station ', obs_elem_c(obsSpaceData,'STID',headerIndex ),  &
          obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex ), obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex )
          call obs_bodySet_i(obsSpaceData, OBS_VNM, bodyIndex, -1 )

        end if

      end do body

    end do header

  end subroutine ovt_windSpeedDirectionToUV

  !--------------------------------------------------------------------------
  ! ovt_UVtoWindSpeedDirection_residual
  !--------------------------------------------------------------------------
  subroutine ovt_UVtoWindSpeedDirection_residual(obsSpaceData, elem_i)
    !
    ! :Purpose: To transform wind residuals in terms of u and v to speed and direction
    !
    implicit none

    ! arguments
    type(struct_obs)    :: obsSpaceData
    integer, intent(in) :: elem_i

    ! locals
    integer :: iuu, ivv, iff, idd
    integer :: headerIndex, headerIndexStart, headerIndexEnd, jWindType, bodyIndex, bodyIndex2
    real(obs_real) :: zlevu, module, ang, uu, vv
    logical :: ok

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
      body: do bodyIndex = 1, obs_numBody(obsSpaceData)

        ok = ( obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex ) == obs_assimilated .and. &
                 obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex) == iuu )

        if ( ok ) then
          headerIndex      = obs_bodyElem_i(obsSpaceData, OBS_HIND, bodyIndex  )
          headerIndexStart = obs_headElem_i(obsSpaceData, OBS_RLN , headerIndex )
          headerIndexEnd   = obs_headElem_i(obsSpaceData, OBS_NLV , headerIndex ) + headerIndexStart - 1
          zlevu            = obs_bodyElem_r(obsSpaceData, OBS_PPP , bodyIndex )
          uu = - obs_bodyElem_r(obsSpaceData, elem_i, bodyIndex ) + obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
         
           body2: do bodyIndex2 = headerIndexStart, headerIndexEnd
     
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

          end do body2
          
          ! insert resduals into obsSpaceData
          body2_2: do bodyIndex2 = headerIndexStart, headerIndexEnd

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

          end do body2_2

        end if

      end do body

    end do WIND_TYPE

  end subroutine ovt_UVtoWindSpeedDirection_residual

  !--------------------------------------------------------------------------
  ! ovt_visToLogVis
  !--------------------------------------------------------------------------
  subroutine ovt_visToLogVis(obsSpaceData, headerIndexStart, headerIndexEnd)
    !
    ! :Purpose: To transform visibily observation to log(visibility)
    !
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsSpaceData
    integer          , intent(in)    :: headerIndexStart 
    integer          , intent(in)    :: headerIndexEnd

    ! locals
    integer        :: headerIndex, bodyIndex, bodyIndexStart, bodyIndexEnd, bodyIndex2
    integer        :: visFlag, logVisFlag
    real(obs_real) :: visObs, visLevel, logVisObs, level
    logical        :: logVisFound

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
        logVisFound = .false.
        body2: do bodyIndex2 = bodyIndex, bodyIndexEnd

          level = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2 )

          if ( level /= visLevel ) cycle body2

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

          logVisFound = .true.
          exit body2

        end do body2

        if (.not. logVisFound) then
          call utl_abort('ovt_visToLogVis: logVis bodyIndex not found!')
        end if

      end do body

    end do header

  end subroutine ovt_visToLogVis

  !--------------------------------------------------------------------------
  ! ovt_visToLogVis_redidual
  !--------------------------------------------------------------------------
  subroutine ovt_visToLogVis_residual(obsSpaceData, residualTypeID)
    !
    ! :Purpose: To transform log(visibily) residuals to visibility
    !
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsSpaceData
    integer          , intent(in)    :: residualTypeID

    ! locals
    integer        :: headerIndex, bodyIndex, bodyIndexStart, bodyIndexEnd, bodyIndex2
    real(obs_real) :: visObs, logVisLevel, logVisObs, level
    real(obs_real) :: visResidual, logVisResidual
    logical        :: visFound

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
      visFound = .false.
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

        visFound = .true.
        exit body2

      end do body2
      
      if (.not. visFound) then
        call utl_abort('ovt_visToLogVis_residual: vis bodyIndex not found!')
      end if

    end do body

  end subroutine ovt_visToLogVis_residual

  !--------------------------------------------------------------------------
  ! ovt_adjustHumGZ
  !--------------------------------------------------------------------------
  subroutine  ovt_adjustHumGZ(obsSpaceData, headerIndexStart, headerIndexEnd )
    !
    ! :Purpose: To apply a threshold on dew-point departure values and to 
    !           transform geopotential height values to geopotential
    !
    implicit none

    ! arguments
    type (struct_obs), intent(inout) :: obsSpaceData
    integer          , intent(in)    :: headerIndexStart
    integer          , intent(in)    :: headerIndexEnd

    ! locals
    integer  :: bodyIndex, headerIndex, bodyIndexStart, bodyIndexEnd
    integer  :: bufrCode

    real(obs_real), parameter :: ESmax = 30.0
    real(obs_real) :: gz, obsValue

    do headerIndex = headerIndexStart, headerIndexEnd 

      bodyIndexStart = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex )
      bodyIndexEnd   = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex ) + bodyIndexStart - 1
      
      do bodyIndex = bodyIndexStart, bodyIndexEnd

        bufrCode = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex )

        select case(bufrCode)
        case(bufr_nees, bufr_ness)
          obsValue = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
          if ( obsValue > ESmax ) obsValue = ESmax
          call obs_bodySet_r(obsSpaceData, OBS_VAR, bodyIndex, obsValue )
        case(bufr_negz)
          obsValue = obs_bodyElem_r(obsSpaceData,OBS_VAR, bodyIndex )
          gz = obsValue * grav
          call obs_bodySet_r(obsSpaceData, OBS_VAR, bodyIndex, gz )
        end select

      end do

    end do

  end subroutine ovt_adjustHumGZ
  
end module obsVariableTransforms_mod

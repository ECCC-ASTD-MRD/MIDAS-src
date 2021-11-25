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
  public :: ovt_getDestinationBufrCode, ovt_getSourceBufrCode, ovt_bufrCodeSkipped
  public :: ovt_isWindObs, ovt_isTransformedVariable, ovt_adjustHumGZ

  integer, parameter :: nTransformSupported  = 3

  type :: struct_transform
    character(len=48)    :: name
    integer              :: nBufrCode
    integer, allocatable :: sourceBufrCode          (:)
    integer, allocatable :: sourceBufrCodeExtra     (:)
    integer, allocatable :: destinationBufrCode     (:)
    integer, allocatable :: destinationBufrCodeExtra(:)
    logical              :: wind   = .false.
    logical              :: active = .false.
  end type struct_transform

  type(struct_transform) :: transform(nTransformSupported)

  integer                :: nSkippedBufrCodes    = 0
  integer                :: skippedBufrCodes(50) = -999

  logical                :: initialized = .false.

contains

  !--------------------------------------------------------------------------
  ! ovt_initStructure
  !--------------------------------------------------------------------------
  subroutine ovt_initStructure
    !
    ! :Purpose: To set the transforms handled by this module
    ! 
    implicit none

    ! Locals:
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

    ! log of precipitation rate
    transformIndex = 3
    transform(transformIndex)%name = 'precipToLogPrecip'
    transform(transformIndex)%nBufrCode = 1
    allocate(transform(transformIndex)%sourceBufrCode     (transform(transformIndex)%nBufrCode))
    allocate(transform(transformIndex)%destinationBufrCode(transform(transformIndex)%nBufrCode))
    transform(transformIndex)%sourceBufrCode     (:) = (/bufr_radarPrecip   /) ! precipitation
    transform(transformIndex)%destinationBufrCode(:) = (/bufr_logRadarPrecip/) ! log(precipitation)
    transform(transformIndex)%wind = .false.

    ! Skipped variables
    nSkippedBufrCodes   = 2
    skippedBufrCodes(1) = bufr_neff  ! because we still want to be able ...
    skippedBufrCodes(2) = bufr_nefs  ! ... to assimilate wind speed obs alone

  end subroutine ovt_initStructure

  !--------------------------------------------------------------------------
  ! ovt_setup
  !--------------------------------------------------------------------------
  subroutine ovt_setup(bufrCodeRead)
    !
    ! :Purpose: To determine which transform must be actived
    !
    implicit none

    ! Argument:
    integer, intent(in) :: bufrCodeRead(:)          ! The list of bufr code read

    ! Locals:
    integer, allocatable :: bufrCodeAssimilated(:)  ! The list of bufr code assimilated

    integer :: nBufrCodeRead, nBufrCodeAssimilated
    integer :: readBufrCodeIndex, transformIndex, assimBufrCodeIndex

    logical :: variableTransformNeeded, foundTransformation
    logical, save :: firstTime = .true.

    if (firstTime) then
      call ovt_initStructure
      firstTime = .false.
    end if

    nBufrCodeAssimilated = filt_nBufrCodeAssimilated()
    nBufrCodeRead = size(bufrCodeRead)
    variableTransformNeeded = .false.

    allocate(bufrCodeAssimilated(nBufrCodeAssimilated))
    call filt_getBufrCodeAssimilated(bufrCodeAssimilated)

    do readBufrCodeIndex = 1, nBufrCodeRead
    
      ! Check if a transform is neeeded
      if (filt_bufrCodeAssimilated(bufrCodeRead(readBufrCodeIndex)) .or. &
          bufrCodeRead(readBufrCodeIndex) == bufr_neff              .or. &
          bufrCodeRead(readBufrCodeIndex) == bufr_nefs ) then
        cycle ! No transformation needed. Move on.
              ! Note that this is where we decide that wind speed will be ignored, 
              ! because we do all the appropriate wind manipulations when we encounter direction.
              ! This strategy allow to assimilate lonely speed observations, like SAR winds.
      end if

      ! Find and activate the appropriate transform
      foundTransformation = .false.

      transformLoop : do transformIndex = 1, nTransformSupported
        assimBufrCodeLoop : do assimBufrCodeIndex = 1, nBufrCodeAssimilated
          if (any(transform(transformIndex)%sourceBufrCode(:) == bufrCodeRead(readBufrCodeIndex))        .and. &
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
        if ( .not. any( skippedBufrCodes(:) == bufrCodeRead(readBufrCodeIndex) ) ) then
          write(*,*)
          write(*,*) 'ovt_setup: !WARNING! No transform found for the read bufr code = ', bufrCodeRead(readBufrCodeIndex)
          write(*,*) '           We are assuming that this observation is read but not assimilated.'
          write(*,*) '           Please consider removing this bufr code from the read observation list'
          write(*,*) '           unless the bufr element is read for another purpose.'
          nSkippedBufrCodes = nSkippedBufrCodes + 1
          skippedBufrCodes(nSkippedBufrCodes) = bufrCodeRead(readBufrCodeIndex)
        end if
      end if

    end do

    initialized = .true.

    deallocate(bufrCodeAssimilated)

  end subroutine ovt_setup

  !--------------------------------------------------------------------------
  ! ovt_bufrCodeSkipped
  !--------------------------------------------------------------------------
  function ovt_bufrCodeSkipped(sourceBufrCode) result(skip)
    !
    ! :Purpose: To NEVER activate a variable transform for this bufr code,
    !           even when this bufr_code is read but not found in the assimilated list.
    !           So far, this function is only used to "skipped" wind speed reports
    !           because we do all the appropriate wind manipulations when we encounter direction.
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: sourceBufrCode ! The input bufr code
    logical :: skip                       ! The decision
    
    ! Locals:
    integer :: bufrCodeIndex

    if (.not. initialized) then
      call utl_abort(' ovt_bufrCodeSkipped: this module has not been setup')
    end if

    skip = .false.
    do bufrCodeIndex = 1, nSkippedBufrCodes
      if (skippedBufrCodes(bufrCodeIndex) == sourceBufrCode) then
        skip = .true.
        exit
      end if
    end do

  end function ovt_bufrCodeSkipped

  !--------------------------------------------------------------------------
  ! ovt_getDestinationBufrCode
  !--------------------------------------------------------------------------
  function ovt_getDestinationBufrCode(sourceBufrCode,extra_opt) result(destinationBufrCode)
    !
    ! :Purpose: To get the bufr code of the transformed/destination variable based on the 
    !           bufr code of the source variable
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: sourceBufrCode    ! The input source bufr code
    integer :: destinationBufrCode           ! The returned destination/transform bufr code
    logical, optional :: extra_opt           ! Should we look in the "extra" bufr code list or not
    
    ! Locals:
    logical :: extra
    integer :: transformIndex, bufrCodeIndex

    if (.not. initialized) then
      call utl_abort(' ovt_getDestinationBufrCode: this module has not been setup')
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
      write(*,*) 'ovt_getDestinationBufrCode: source bufrCode = ', sourceBufrCode
      call utl_abort('ovt_getDestinationBufrCode: found no associated bufrCode for the above source variable bufr code')
    end if

  end function ovt_getDestinationBufrCode

  !--------------------------------------------------------------------------
  ! ovt_getSourceBufrCode
  !--------------------------------------------------------------------------
  function ovt_getSourceBufrCode(destinationBufrCode,extra_opt) result(sourceBufrCode)
    !
    ! :Purpose: To get the bufr code of the source variable based on the 
    !           bufr code of the destination/transformed variable
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: destinationBufrCode  ! The input destination/transform bufr code
    integer :: sourceBufrCode                   ! The returned source bufr code
    logical, optional :: extra_opt              ! Should we look in the "extra" bufr code list or not
    
    ! Locals:
    logical :: extra
    integer :: transformIndex, destinationBufrCodeIndex

    if (.not. initialized) then
      call utl_abort(' ovt_getSourceBufrCode: this module has not been setup')
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
      write(*,*) 'ovt_getSourceBufrCode: tranform variable bufr code = ', destinationBufrCode
      call utl_abort('ovt_getSourceBufrCode: found no associated variable bufr code for the above transform variable bufr code')
    end if

  end function ovt_getSourceBufrCode

  !--------------------------------------------------------------------------
  ! ovt_isWindObs
  !--------------------------------------------------------------------------
  function ovt_isWindObs(sourceBufrCode) result(wind)
    !
    ! :Purpose: To determine if a bufr code is wind related
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: sourceBufrCode ! The input source bufr code
    logical :: wind                       ! Is this bufr code linked to wind or not

    ! Locals:
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

    ! Arguments:
    integer, intent(in) :: bufrCode    ! The input bufr code
    logical :: transformed             ! Is this a bufr code associated to a transform variable or not

    ! Locals:
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
    type (struct_obs), intent(inout) :: obsSpaceData      ! The observation database
    integer          , intent(in)    :: headerIndexStart  ! The initial header index to analyse
    integer          , intent(in)    :: headerIndexEnd    ! The final header index to analyse

    ! Locals:
    integer :: transformIndex

    if (obs_numHeader(obsSpaceData) == 0) return

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
        case ('precipToLogPrecip')
          call ovt_precipToLogPrecip     (obsSpaceData, headerIndexStart, headerIndexEnd)
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
    type (struct_obs), intent(inout) :: obsSpaceData    ! The observation database 
    integer          , intent(in)    :: residualTypeID  ! The residual type ID (o-p or o-a)

    ! Local:
    integer :: transformIndex

    if (obs_numHeader(obsSpaceData) == 0) return

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
        case ('precipToLogPrecip') 
          call ovt_precipToLogPrecip_residual     (obsSpaceData, residualTypeID)
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
    type (struct_obs), intent(inout) :: obsSpaceData      ! The observation database
    integer          , intent(in)    :: headerIndexStart  ! The initial header index to analyse
    integer          , intent(in)    :: headerIndexEnd    ! The final header index to analyse

    ! Locals:
    integer        :: bufrCode, bufrCode2, bufrCode3
    integer        :: headerIndex, bodyIndex, bodyIndexStart, bodyIndexEnd, bodyIndex2, bufrCodeAssociated
    integer        :: directionFlag, speedFlag, combinedFlag, uWindFlag, vWindFlag
    integer        :: speedBufrCode, uWindBufrCode, vWindBufrCode, uWindbodyIndex, vWindBodyIndex

    logical        :: direction_missing, speed_missing
    logical        :: uWind_present, vWind_present

    real(pre_obsReal) :: uWind, vWind, direction, speed
    real(pre_obsReal) :: level_direction, level, level2, level3

    speedFlag = 0

    ! Loop through headers
    header: do headerIndex = headerIndexStart, headerIndexEnd
      
      bodyIndexStart = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex )
      bodyIndexEnd   = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex ) + bodyIndexStart - 1
      
      ! Find the wind direction report
      body: do bodyIndex = bodyIndexStart, bodyIndexEnd 

        direction = obs_missingValue_R
        speed     = obs_missingValue_R
        bufrCode  = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex)
        direction_missing = .true.
      
        if ( bufrCode /= bufr_nedd .and. bufrCode /= bufr_neds ) cycle body

        if( bufrCode == bufr_neds) then
          ! Surface obs
          speedBufrCode = bufr_nefs
          uWindBufrCode = bufr_neus
          vWindBufrCode = bufr_nevs
        else
          ! Upper air obs
          speedBufrCode = bufr_neff
          uWindBufrCode = bufr_neuu
          vWindBufrCode = bufr_nevv
        end if
      
        direction       = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
        directionFlag   = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex )
        level_direction = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex )
        uWind_present = .false.
        vWind_present = .false.

        ! Check if u and v are present in obsSpaceData
        do bodyIndex2 = bodyIndexStart, bodyIndexEnd

          level3 = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2 )

          if (level3 == level_direction) then
            bufrCode3 = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2 )
            
            if ( bufrCode3 == uWindBufrCode ) then
              uWind_present = .true.
              uWindbodyIndex = bodyIndex2
            else if ( bufrCode3 == vWindBufrCode) then
              vWind_present = .true.
              vWindBodyIndex = bodyIndex2
            end if

          end if

        end do

        if (.not. uWind_present .or. .not. vWind_present) then
          call utl_abort('ovt_windSpeedDirectionToUV: uWind and/or vWind bodyIndex not found!')
        end if

        ! Find the speed report and compute uWind and vWind
        calcuv: do bodyIndex2 = bodyIndexStart, bodyIndexEnd
          speed_missing     = .true.
          direction_missing = .true.
          level = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2 )

          if ( level /= level_direction) cycle calcuv

          bufrCode2 = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2 )

          if ( bufrCode2 == speedBufrCode ) then

            speed     = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex2 )
            speedFlag = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex2 )

            if ( direction == 0.d0 .and. speed > 0. .or. direction > 360. .or. direction < 0. ) then
              direction_missing = .true.
              speed_missing     = .true.
            else if ( direction == obs_missingValue_R .or. speed == obs_missingValue_R ) then
              direction_missing = .true.
              speed_missing     = .true.
            else
              direction_missing = .false.
              speed_missing     = .false.
            end if

            if ( direction_missing .or. speed_missing ) then

              call obs_bodySet_i(obsSpaceData, OBS_VNM, uWindbodyIndex, -1 )
              call obs_bodySet_i(obsSpaceData, OBS_VNM, vWindBodyIndex, -1 )

            else

              if (speed == 0.d0) direction = 0.d0
              direction = direction + 180.
              if ( direction > 360.) direction = direction - 360.
              direction = direction * mpc_radians_per_degree_r8
            
              ! (speed,direction) -> (u,v)
              uWind = speed * sin(direction)
              vWind = speed * cos(direction)
          
              combinedFlag = ior(directionFlag, speedFlag )

              call obs_bodySet_r(obsSpaceData, OBS_VAR, uWindbodyIndex, uWind )
              call obs_bodySet_i(obsSpaceData, OBS_FLG, uWindbodyIndex, combinedFlag )
              call obs_bodySet_r(obsSpaceData, OBS_VAR, vWindBodyIndex, vWind )
              call obs_bodySet_i(obsSpaceData, OBS_FLG, vWindBodyIndex, combinedFlag )

            end if

          end if

        end do calcuv

      end do body

    end do header

    !
    !- Merge uWind and vWind flags (JFC: not sure why this is needed because this was already done above)
    !
    header2: do headerIndex = headerIndexStart, headerIndexEnd

      bodyIndexStart = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
      bodyIndexEnd   = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) + bodyIndexStart - 1
    
      ! Search uWind component
      body2: do bodyIndex = bodyIndexStart, bodyIndexEnd
        bufrCode = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex)
        level    = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex)

        select case (bufrCode)
          case (bufr_neuu)
            bufrCodeAssociated = bufr_nevv
          case (bufr_neus)
            bufrCodeAssociated = bufr_nevs        
          case default
            cycle body2
        end select

        uWindbodyIndex = bodyIndex

        ! Eleminate entries where uWind is missing
        uWind = obs_bodyElem_r(obsSpaceData, OBS_VAR, uWindbodyIndex)
        if ( uWind == obs_missingValue_R ) then
          call obs_bodySet_i(obsSpaceData, OBS_VNM, uWindbodyIndex, -1)
        end if

        ! Search the associated vWind component
        vWindBodyIndex = -1
        body3: do bodyIndex2 = bodyIndexStart, bodyIndexEnd
          bufrCode2 = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2)
          level2    = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2)

          if ( bufrCode2 /= bufrCodeAssociated .or. level2 /= level ) cycle

          vWindBodyIndex = bodyIndex2

          if ( uWind == obs_missingValue_R ) then
            call obs_bodySet_i(obsSpaceData, OBS_VNM, vWindBodyIndex, -1)
          else
            uWindFlag = obs_bodyElem_i(obsSpaceData, OBS_FLG, uWindbodyIndex)
            vWindFlag = obs_bodyElem_i(obsSpaceData, OBS_FLG, vWindBodyIndex)
            combinedFlag = ior( uWindFlag, vWindFlag )
            call obs_bodySet_i(obsSpaceData, OBS_FLG, uWindbodyIndex, combinedFlag)
            call obs_bodySet_i(obsSpaceData, OBS_FLG, vWindBodyIndex, combinedFlag)
          end if

          exit body3

        end do body3

        ! Eleminate entries where vWind is missing
        if (vWindBodyIndex < 0 .and. uWind /= obs_missingValue_R) then
          call obs_bodySet_i(obsSpaceData, OBS_VNM, uWindbodyIndex, -1)
        end if

      end do body2

    end do header2

  end subroutine ovt_windSpeedDirectionToUV

  !--------------------------------------------------------------------------
  ! ovt_UVtoWindSpeedDirection_residual
  !--------------------------------------------------------------------------
  subroutine ovt_UVtoWindSpeedDirection_residual(obsSpaceData, residualTypeID)
    !
    ! :Purpose: To transform wind residuals in terms of u and v to speed and direction
    !
    implicit none

    ! Arguments:
    type(struct_obs)    :: obsSpaceData    ! The observation database
    integer, intent(in) :: residualTypeID  ! The residual type ID (o-p or o-a)

    ! Locals:
    integer :: uWindBufrCode, vWindBufrCode, speedBufrCode, directionBufrCode
    integer :: headerIndex, headerIndexStart, headerIndexEnd, windTypeIndex, bodyIndex, bodyIndex2

    real(pre_obsReal) :: uWindLevel, speed, direction, uWind, vWind

    windType: do windTypeIndex = 1, 2
      if (windTypeIndex == 1) then
        uWindBufrCode = bufr_neuu
        vWindBufrCode = bufr_nevv
        directionBufrCode = bufr_nedd
        speedBufrCode = bufr_neff
      else
        uWindBufrCode = bufr_neus
        vWindBufrCode = bufr_nevs
        directionBufrCode = bufr_neds
        speedBufrCode = bufr_nefs
      end if

      ! Process all data within the domain of the model
      body: do bodyIndex = 1, obs_numBody(obsSpaceData)

        if ( obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) == obs_assimilated .and. &
             obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex) == uWindBufrCode ) then

          headerIndex      = obs_bodyElem_i(obsSpaceData, OBS_HIND, bodyIndex  )
          headerIndexStart = obs_headElem_i(obsSpaceData, OBS_RLN , headerIndex )
          headerIndexEnd   = obs_headElem_i(obsSpaceData, OBS_NLV , headerIndex ) + headerIndexStart - 1
          uWindLevel            = obs_bodyElem_r(obsSpaceData, OBS_PPP , bodyIndex )
          uWind = - obs_bodyElem_r(obsSpaceData, residualTypeID, bodyIndex ) + obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
         
           body2: do bodyIndex2 = headerIndexStart, headerIndexEnd
     
             if  ( obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2) == vWindBufrCode .and. &
                   obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2) == uWindLevel ) then
               
              vWind = -obs_bodyElem_r(obsSpaceData, residualTypeID, bodyIndex2 ) + obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex2 )

              ! Calculate angle
              speed = sqrt((uWind**2)+(vWind**2))

              if (speed == 0.) then
                direction = 0.0d0
              else
                direction = atan2(vWind,uWind)
                direction = (270.0d0 - direction  * MPC_DEGREES_PER_RADIAN_R8 )
                ! Change to meteorological definition of wind direction.
                if (direction > 360.0d0 ) direction = direction - 360.0d0
                if (direction <= 0.0d0  ) direction = direction + 360.0d0
              end if

            end if

          end do body2
          
          ! insert resduals into obsSpaceData
          body2_2: do bodyIndex2 = headerIndexStart, headerIndexEnd

            if  ( obs_bodyElem_i( obsSpaceData, OBS_VNM, bodyIndex2) == directionBufrCode .and. &
                  obs_bodyElem_r( obsSpaceData, OBS_PPP, bodyIndex2) == uWindLevel ) then
              call obs_bodySet_r( obsSpaceData, residualTypeID, bodyIndex2, &
                                  obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex2 ) - direction )
              
              if ( obs_bodyElem_r( obsSpaceData,residualTypeID,bodyIndex2) >  180.0d0)  &
                call obs_bodySet_r( obsSpaceData, residualTypeID, bodyIndex2, &
                                    obs_bodyElem_r(obsSpaceData, residualTypeID, bodyIndex2 ) - real(360.0d0,pre_obsReal) )
              if ( obs_bodyElem_r(obsSpaceData,residualTypeID,bodyIndex2) <= -180.0d0)  &
                call obs_bodySet_r( obsSpaceData, residualTypeID, bodyIndex2, &
                                    obs_bodyElem_r(obsSpaceData, residualTypeID, bodyIndex2 ) + real(360.0d0,pre_obsReal) )
              call obs_bodySet_r( obsSpaceData, residualTypeID, bodyIndex2, &
                                  - real(1.0d0,pre_obsReal) * obs_bodyElem_r(obsSpaceData,residualTypeID,bodyIndex2 ) )
              call obs_bodySet_r( obsSpaceData, OBS_OER, bodyIndex2, real(1.0d0,pre_obsReal) )
              call obs_bodySet_i( obsSpaceData, OBS_ASS, bodyIndex2, obs_assimilated )
              call obs_bodySet_i( obsSpaceData, OBS_FLG, bodyIndex2, 0 )
            end if
            if  ( obs_bodyElem_i( obsSpaceData, OBS_VNM, bodyIndex2 ) == speedBufrCode .and. &
                  obs_bodyElem_r( obsSpaceData, OBS_PPP, bodyIndex2 ) == uWindLevel ) then
              call obs_bodySet_r( obsSpaceData, residualTypeID,  bodyIndex2, &
                                  obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex2 ) - speed )
              call obs_bodySet_r( obsSpaceData, OBS_OER, bodyIndex2, real(1.0d0,pre_obsReal) )
              call obs_bodySet_i( obsSpaceData, OBS_ASS, bodyIndex2, obs_assimilated)
              call obs_bodySet_i( obsSpaceData, OBS_FLG, bodyIndex2, 0)
            end if

          end do body2_2

        end if

      end do body

    end do windType

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
    type (struct_obs), intent(inout) :: obsSpaceData      ! The observation database
    integer          , intent(in)    :: headerIndexStart  ! The initial header index to analyse
    integer          , intent(in)    :: headerIndexEnd    ! The final header index to analyse

    ! Locals:
    integer        :: headerIndex, bodyIndex, bodyIndexStart, bodyIndexEnd, bodyIndex2
    integer        :: visFlag, logVisFlag
    real(pre_obsReal) :: visObs, visLevel, logVisObs, level
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

          if (visObs == obs_missingValue_R) then
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
  ! ovt_visToLogVis_residual
  !--------------------------------------------------------------------------
  subroutine ovt_visToLogVis_residual(obsSpaceData, residualTypeID)
    !
    ! :Purpose: To transform log(visibily) residuals to visibility
    !
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsSpaceData    ! The observation database
    integer          , intent(in)    :: residualTypeID  ! The residual type ID (o-p or o-a)

    ! Locals:
    integer        :: headerIndex, bodyIndex, bodyIndexStart, bodyIndexEnd, bodyIndex2
    real(pre_obsReal) :: visObs, logVisLevel, logVisObs, level
    real(pre_obsReal) :: visResidual, logVisResidual
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
        call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex2, obs_missingValue_R)

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
  ! ovt_precipToLogPrecip
  !--------------------------------------------------------------------------
  subroutine ovt_precipToLogPrecip(obsSpaceData, headerIndexStart, headerIndexEnd)
    !
    ! :Purpose: To transform precipitation observation to log(precipitation)
    !
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsSpaceData      ! The observation database
    integer          , intent(in)    :: headerIndexStart  ! The initial header index to analyse
    integer          , intent(in)    :: headerIndexEnd    ! The final header index to analyse

    ! Locals:
    integer        :: headerIndex, bodyIndex, bodyIndexStart, bodyIndexEnd, bodyIndex2
    integer        :: precipFlag, logPrecipFlag
    real(pre_obsReal) :: precipObs, precipLevel, logPrecipObs, level
    logical        :: logPrecipFound

    ! Loop through headers
    header: do headerIndex = headerIndexStart, headerIndexEnd
      
      bodyIndexStart = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex )
      bodyIndexEnd   = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex ) + bodyIndexStart - 1
      
      ! Find each precipitation report
      body: do bodyIndex = bodyIndexStart, bodyIndexEnd 

        if (obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex) /= bufr_radarPrecip) cycle body

        precipObs   = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
        precipFlag  = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex )
        precipLevel = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex )

        ! Find the associated logPrecip body created earlier in the proper reading routine
        logPrecipFound = .false.
        body2: do bodyIndex2 = bodyIndex, bodyIndexEnd

          level = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2 )

          if ( level /= precipLevel ) cycle body2

          if (obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2 ) /= bufr_logRadarPrecip) cycle body2

          if (precipObs == obs_missingValue_R) then
            logPrecipObs = precipObs
          else
            ! precip -> log(precip)
            logPrecipObs = log(MPC_MINIMUM_PR_R8 + max(0.0d0,precipObs))
          end if
          call obs_bodySet_r(obsSpaceData, OBS_VAR, bodyIndex2, logPrecipObs)

          logPrecipFlag  = precipFlag
          call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex2, logPrecipFlag )

          logPrecipFound = .true.
          exit body2

        end do body2

        if (.not. logPrecipFound) then
          call utl_abort('ovt_precipToLogPrecip: logPrecip bodyIndex not found!')
        end if

      end do body

    end do header

  end subroutine ovt_precipToLogPrecip

  !--------------------------------------------------------------------------
  ! ovt_precipToLogPrecip_residual
  !--------------------------------------------------------------------------
  subroutine ovt_precipToLogPrecip_residual(obsSpaceData, residualTypeID)
    !
    ! :Purpose: To transform log(precip) residuals to precip
    !
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsSpaceData    ! The observation database
    integer          , intent(in)    :: residualTypeID  ! The residual type ID (o-p or o-a)

    ! Locals:
    integer        :: headerIndex, bodyIndex, bodyIndexStart, bodyIndexEnd, bodyIndex2
    real(pre_obsReal) :: precipObs, logPrecipLevel, logPrecipObs, level
    real(pre_obsReal) :: precipResidual, logPrecipResidual
    logical        :: precipFound

    ! Find each log of precipitation assimilated observations
    body: do bodyIndex = 1, obs_numBody(obsSpaceData)
      
      if ( obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) /= obs_assimilated .or. &
           obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex) /= bufr_logRadarPrecip ) cycle

      headerIndex    = obs_bodyElem_i(obsSpaceData, OBS_HIND, bodyIndex  )
      bodyIndexStart = obs_headElem_i(obsSpaceData, OBS_RLN , headerIndex )
      bodyIndexEnd   = obs_headElem_i(obsSpaceData, OBS_NLV , headerIndex ) + bodyIndexStart - 1
      logPrecipLevel    = obs_bodyElem_r(obsSpaceData, OBS_PPP , bodyIndex )

      logPrecipObs      = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
      logPrecipResidual = obs_bodyElem_r(obsSpaceData, residualTypeID, bodyIndex )

      ! Find the associated precip body
      precipFound = .false.
      body2: do bodyIndex2 = bodyIndexStart, bodyIndexEnd

        level = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex2 )

        if (level /= logPrecipLevel) cycle body2
        if (obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2) /= bufr_radarPrecip) cycle body2

        ! log(precip) -> precip
        precipObs      = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex2 )
        precipResidual = precipObs - max(0.0D0,  &               ! o-p or o-a
                         exp(logPrecipObs - logPrecipResidual - MPC_MINIMUM_PR_R8))
        call obs_bodySet_r(obsSpaceData, residualTypeID, bodyIndex2, precipResidual)

        ! Set the obs error to missing
        call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex2, obs_missingValue_R)

        ! Set flags
        call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex2, obs_assimilated)
        call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex2, 0)

        precipFound = .true.
        exit body2

      end do body2
      
      if (.not. precipFound) then
        call utl_abort('ovt_precipToLogPrecip_residual: precip bodyIndex not found!')
      end if

    end do body

  end subroutine ovt_precipToLogPrecip_residual

  !--------------------------------------------------------------------------
  ! ovt_adjustHumGZ
  !--------------------------------------------------------------------------
  subroutine  ovt_adjustHumGZ(obsSpaceData, headerIndexStart, headerIndexEnd )
    !
    ! :Purpose: To apply a threshold on dew-point departure values and to 
    !           transform geopotential height values to geopotential
    !
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsSpaceData      ! The observation database
    integer          , intent(in)    :: headerIndexStart  ! The initial header index to analyse
    integer          , intent(in)    :: headerIndexEnd    ! The final header index to analyse

    ! Locals:
    integer  :: bodyIndex, headerIndex, bodyIndexStart, bodyIndexEnd
    integer  :: bufrCode

    real(pre_obsReal), parameter :: ESmax = 30.0
    real(pre_obsReal) :: gz, obsValue

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

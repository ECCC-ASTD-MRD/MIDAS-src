
module obsFlags_mod
  ! MODULE obsFlags_mod (prefix='flg' category='8. Low-level utilities and constants')
  !
  !:Purpose:  Contains a list of all possible bit numbers used when setting or
  !           interpreting the bits for indicating the status of observations.
  !
  use obsSpaceData_mod
  use midasMpi_mod
  use utilities_mod

  implicit none
  save
  private

  ! Public variables (parameters)
  !public :: flg_

  ! Public procedures
  public :: flg_bitNumtoRefValue, flg_refValueTobitNum
  public :: flg_flagIsOn

  ! interface for testing bit value
  interface flg_flagIsOn
    module procedure flg_flagIsOnFromObsData
    module procedure flg_flagIsOnFromFlagInt
  end interface flg_flagIsOn

contains

  !--------------------------------------------------------------------------
  ! flg_flagIsOn
  !--------------------------------------------------------------------------
  function flg_flagIsOnFromObsData(obsSpaceData, bodyIndex,  &
                                   bitNumber_opt, refValue_opt) result(flagIsOn)
    !
    ! :Purpose: Check if a specific flag bit is on.
    !
    implicit none

    ! Arguments:
    integer,                    intent(in) :: bodyIndex
    type(struct_obs),           intent(in) :: obsSpaceData
    integer,          optional, intent(in) :: bitNumber_opt
    integer,          optional, intent(in) :: refValue_opt
    ! Result:
    logical                       :: flagIsOn

    ! Locals:
    integer :: bitNumber, flagInteger

    if (present(bitNumber_opt) .and. present(refValue_opt)) then
      call utl_abort('flg_flagIsOn: only one of bitNumber_opt and refValue_opt can be set')
    end if

    if (.not.present(bitNumber_opt) .and. .not.present(refValue_opt)) then
      call utl_abort('flg_flagIsOn: either bitNumber_opt or refValue_opt must be set')
    end if

    if (present(bitNumber_opt)) then
      bitNumber = bitNumber_opt
    else
      bitNumber = flg_refValueToBitNum(refValue_opt)
    end if

    flagInteger = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex)
    flagIsOn = btest(flagInteger, bitNumber)

  end function flg_flagIsOnFromObsData

  !--------------------------------------------------------------------------
  ! flg_flagIsOn
  !--------------------------------------------------------------------------
  function flg_flagIsOnFromFlagInt(bitNumber_opt, refValue_opt, &
                                   flagInteger) result(flagIsOn)
    !
    ! :Purpose: Check if a specific flag bit is on.
    !
    implicit none

    ! Arguments:
    integer,                    intent(in) :: flagInteger
    integer,          optional, intent(in) :: bitNumber_opt
    integer,          optional, intent(in) :: refValue_opt
    ! Result:
    logical                       :: flagIsOn

    ! Locals:
    integer :: bitNumber

    if (present(bitNumber_opt) .and. present(refValue_opt)) then
      call utl_abort('flg_flagIsOn: only one of bitNumber_opt and refValue_opt can be set')
    end if

    if (.not.present(bitNumber_opt) .and. .not.present(refValue_opt)) then
      call utl_abort('flg_flagIsOn: either bitNumber_opt or refValue_opt must be set')
    end if

    if (present(bitNumber_opt)) then
      bitNumber = bitNumber_opt
    else
      bitNumber = flg_refValueToBitNum(refValue_opt)
    end if

    flagIsOn = btest(flagInteger, bitNumber)

  end function flg_flagIsOnFromFlagInt

  !--------------------------------------------------------------------------
  ! flg_bitNumtoRefValue
  !--------------------------------------------------------------------------
  function flg_bitNumtoRefValue(bitNumber) result(refValue)
    !
    ! :Purpose: Convert bit numbers to reference values using the magic formula.
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: bitNumber
    ! Result:
    integer             :: refValue

    refValue = 13 - bitNumber

  end function flg_bitNumtoRefValue

  !--------------------------------------------------------------------------
  ! flg_refValueTobitNum
  !--------------------------------------------------------------------------
  function flg_refValueTobitNum(refValue) result(bitNumber)
    !
    ! :Purpose: Convert bit numbers to reference values using the magic formula.
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: refValue
    ! Result:
    integer             :: bitNumber

    bitNumber = 13 - refValue

  end function flg_refValueTobitNum

end module obsFlags_mod

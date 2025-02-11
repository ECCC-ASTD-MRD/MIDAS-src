
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
  public :: flg_flagIsOn, flg_flagIsOnFromRefValue
  public :: flg_resetFlag, flg_setFlag, flg_clearFlag

  ! interface for testing bit value
  interface flg_flagIsOn
    module procedure flg_flagIsOnFromObsData
    module procedure flg_flagIsOnFromFlagInt
    module procedure flg_flagsAreOnFromObsData
    module procedure flg_flagsAreOnFromFlagInt
  end interface flg_flagIsOn

  interface flg_setFlag
    module procedure flg_setFlagFromObsData
    module procedure flg_setFlagsFromObsData
    module procedure flg_setFlagFromFlagInt
    module procedure flg_setFlagsFromFlagInt
  end interface flg_setFlag

  interface flg_clearFlag
    module procedure flg_clearFlagFromObsData
    module procedure flg_clearFlagFromFlagInt
  end interface flg_clearFlag

contains

  !--------------------------------------------------------------------------
  ! flg_flagIsOnFromRefValue
  !--------------------------------------------------------------------------
  function flg_flagIsOnFromRefValue(obsSpaceData, bodyIndex,  &
                                    refValue) result(flagIsOn)
    !
    ! :Purpose: Check if a specific flag bit is on.
    !
    implicit none

    ! Arguments:
    integer,          intent(in) :: bodyIndex
    type(struct_obs), intent(in) :: obsSpaceData
    integer,          intent(in) :: refValue
    ! Result:
    logical                      :: flagIsOn

    ! Locals:
    integer :: bitNumber, flagInteger

    bitNumber = flg_refValueToBitNum(refValue)
    flagInteger = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex)
    flagIsOn = btest(flagInteger, bitNumber)

  end function flg_flagIsOnFromRefValue

  !--------------------------------------------------------------------------
  ! flg_flagIsOnFromObsData
  !--------------------------------------------------------------------------
  function flg_flagIsOnFromObsData(obsSpaceData, bodyIndex,  &
                                   bitNumber) result(flagIsOn)
    !
    ! :Purpose: Check if a specific flag bit is on.
    !
    implicit none

    ! Arguments:
    integer,          intent(in) :: bodyIndex
    type(struct_obs), intent(in) :: obsSpaceData
    integer,          intent(in) :: bitNumber
    ! Result:
    logical                      :: flagIsOn

    ! Locals:
    integer :: flagInteger

    flagInteger = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex)
    flagIsOn = btest(flagInteger, bitNumber)

  end function flg_flagIsOnFromObsData

  !--------------------------------------------------------------------------
  ! flg_flagIsOnFromFlagInt
  !--------------------------------------------------------------------------
  function flg_flagIsOnFromFlagInt(flagInteger, bitNumber) result(flagIsOn)
    !
    ! :Purpose: Check if a specific flag bit is on.
    !
    implicit none

    ! Arguments:
    integer,          intent(in) :: flagInteger
    integer,          intent(in) :: bitNumber
    ! Result:
    logical                      :: flagIsOn

    flagIsOn = btest(flagInteger, bitNumber)

  end function flg_flagIsOnFromFlagInt

  !--------------------------------------------------------------------------
  ! flg_flagsAreOnFromObsData
  !--------------------------------------------------------------------------
  function flg_flagsAreOnFromObsData(booleanOper, obsSpaceData, bodyIndex,  &
                                     bitNumbers) result(flagsAreOn)
    !
    ! :Purpose: Check if a specific flag bit is on.
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: booleanOper
    type(struct_obs), intent(in) :: obsSpaceData
    integer,          intent(in) :: bodyIndex
    integer,          intent(in) :: bitNumbers(:)
    ! Result:
    logical                      :: flagsAreOn

    ! Locals:
    integer :: flagInteger, numBitNumbers, bitIndex

    numBitNumbers = size(bitNumbers)

    flagInteger = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex)

    if (trim(booleanOper) == 'OR') then
      flagsAreOn = .false.
      do bitIndex = 1, numBitNumbers
        flagsAreOn = flagsAreOn .or. btest(flagInteger, bitNumbers(bitIndex))
      end do
    else
      write(*,*) 'booleanOper = ***', trim(booleanOper), '***'
      call utl_abort('flg_flagsAreOnFromObsData: Only OR is implemented so far')
    end if

  end function flg_flagsAreOnFromObsData

  !--------------------------------------------------------------------------
  ! flg_flagsAreOnFromFlagInt
  !--------------------------------------------------------------------------
  function flg_flagsAreOnFromFlagInt(booleanOper, flagInteger, &
                                     bitNumbers) result(flagsAreOn)
    !
    ! :Purpose: Check if a specific flag bit is on.
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: booleanOper
    integer,          intent(in) :: flagInteger
    integer,          intent(in) :: bitNumbers(:)
    ! Result:
    logical                      :: flagsAreOn

    ! Locals:
    integer :: numBitNumbers, bitIndex

    numBitNumbers = size(bitNumbers)

    if (trim(booleanOper) == 'OR') then
      flagsAreOn = .false.
      do bitIndex = 1, numBitNumbers
        flagsAreOn = flagsAreOn .or. btest(flagInteger, bitNumbers(bitIndex))
      end do
    else
      write(*,*) 'booleanOper = ***', trim(booleanOper), '***'
      call utl_abort('flg_flagsAreOnFromFlagInt: Only OR is implemented so far')
    end if

  end function flg_flagsAreOnFromFlagInt

  !--------------------------------------------------------------------------
  ! flg_resetFlag
  !--------------------------------------------------------------------------
  subroutine flg_resetFlag(obsSpaceData, bodyIndex)
    !
    ! :Purpose: Reset the value of the flag integer to zero.
    !
    implicit none

    ! Arguments:
    integer,                    intent(in)    :: bodyIndex
    type(struct_obs),           intent(inout) :: obsSpaceData

    call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex, 0)
 
  end subroutine flg_resetFlag

  !--------------------------------------------------------------------------
  ! flg_setFlagFromObsData
  !--------------------------------------------------------------------------
  subroutine flg_setFlagFromObsData(obsSpaceData, bodyIndex,  &
                                    bitNumber)
    !
    ! :Purpose: Set the value of a specific flag bit.
    !
    implicit none

    ! Arguments:
    integer,          intent(in)    :: bodyIndex
    type(struct_obs), intent(inout) :: obsSpaceData
    integer,          intent(in)    :: bitNumber

    ! Locals:
    integer :: flagIntegerOrig, flagIntegerUpdated

    flagIntegerOrig = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex)
    flagIntegerUpdated = ibset(flagIntegerOrig, bitNumber)
    call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex, flagIntegerUpdated)
 
  end subroutine flg_setFlagFromObsData

  !--------------------------------------------------------------------------
  ! flg_setFlagsFromObsData
  !--------------------------------------------------------------------------
  subroutine flg_setFlagsFromObsData(obsSpaceData, bodyIndex, bitNumbers)
    !
    ! :Purpose: Set the value of multiple specified flag bits.
    !
    implicit none

    ! Arguments:
    integer,                    intent(in)    :: bodyIndex
    type(struct_obs),           intent(inout) :: obsSpaceData
    integer,                    intent(in)    :: bitNumbers(:)

    ! Locals:
    integer :: flagIntegerOrig, flagIntegerUpdated
    integer :: numBitsToSet, bitIndex

    numBitsToSet = size(bitNumbers)

    do bitIndex = 1, numBitsToSet
      flagIntegerOrig = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex)
      flagIntegerUpdated = ibset(flagIntegerOrig, bitNumbers(bitIndex))
      call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex, flagIntegerUpdated)
    end do
 
  end subroutine flg_setFlagsFromObsData

  !--------------------------------------------------------------------------
  ! flg_setFlagFromFlagInt
  !--------------------------------------------------------------------------
  subroutine flg_setFlagFromFlagInt(flagInteger, bitNumber)
    !
    ! :Purpose: Set the value of a specific flag bit.
    !
    implicit none

    ! Arguments:
    integer,          intent(inout) :: flagInteger
    integer,          intent(in)    :: bitNumber

    flagInteger = ibset(flagInteger, bitNumber)
 
  end subroutine flg_setFlagFromFlagInt

  !--------------------------------------------------------------------------
  ! flg_setFlagsFromFlagInt
  !--------------------------------------------------------------------------
  subroutine flg_setFlagsFromFlagInt(flagInteger, bitNumbers)
    !
    ! :Purpose: Set the value of multiple specified flag bits.
    !
    implicit none

    ! Arguments:
    integer,                    intent(inout) :: flagInteger
    integer,                    intent(in)    :: bitNumbers(:)

    ! Locals:
    integer :: numBitsToSet, bitIndex

    numBitsToSet = size(bitNumbers)

    do bitIndex = 1, numBitsToSet
      flagInteger = ibset(flagInteger, bitNumbers(bitIndex))
    end do
 
  end subroutine flg_setFlagsFromFlagInt

  !--------------------------------------------------------------------------
  ! flg_clearFlagFromObsData
  !--------------------------------------------------------------------------
  subroutine flg_clearFlagFromObsData(obsSpaceData, bodyIndex,  &
                                      bitNumber)
    !
    ! :Purpose: Clear the value of a specific flag bit.
    !
    implicit none

    ! Arguments:
    integer,          intent(in)    :: bodyIndex
    type(struct_obs), intent(inout) :: obsSpaceData
    integer,          intent(in)    :: bitNumber

    ! Locals:
    integer :: flagIntegerOrig, flagIntegerUpdated

    flagIntegerOrig = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex)
    flagIntegerUpdated = ibclr(flagIntegerOrig, bitNumber)
    call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex, flagIntegerUpdated)
 
  end subroutine flg_clearFlagFromObsData

  !--------------------------------------------------------------------------
  ! flg_clearFlagFromFlagInt
  !--------------------------------------------------------------------------
  subroutine flg_clearFlagFromFlagInt(flagInteger, bitNumber)
    !
    ! :Purpose: Clear the value of a specific flag bit.
    !
    implicit none

    ! Arguments:
    integer,          intent(inout) :: flagInteger
    integer,          intent(in)    :: bitNumber

    flagInteger = ibclr(flagInteger, bitNumber)
 
  end subroutine flg_clearFlagFromFlagInt

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

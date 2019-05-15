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

module thinning_mod
  ! MODULE thinning (prefix='thn' category='1. High-level functionality')
  !
  ! :Purpose: Using observation-type-specific algorithms, set bit 11 of 'flag'
  !           on data that are not to be assimilated.
  !
  ! :Note:    This module is intended to group all of the thinning methods in a
  !           single fortran module.
  !
  !           So far, only aladin wind data are treated.
  !
  use codePrecision_mod
  use bufr_mod
  use obsSpaceData_mod
  use utilities_mod
  implicit none
  private
  public :: thn_thinAladin

contains

  subroutine thn_thinAladin(obsdat)
    implicit none

    ! ARGUMENTS
    type(struct_obs), intent(inout) :: obsdat


    ! NAMELIST VARIABLES
    integer :: keepNthVertical ! keep every nth vertical datum

    namelist /thin_aladin/keepNthVertical


    ! LOCAL VARIABLES
    integer :: nulnam
    integer :: fnom, fclos, ierr



    keepNthVertical=-1

    ! Read the namelist for Aladin observations
    nulnam = 0
    ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    if(ierr.ne.0) call utl_abort('thn_thinAladin: Error opening file flnml')
    read(nulnam,nml=thin_aladin,iostat=ierr)
    if(ierr.ne.0) call utl_abort('thn_thinAladin: Error reading namelist')
    write(*,nml=thin_aladin)
    ierr=fclos(nulnam)

    if(keepNthVertical > 0) then
      call keepNthObs(obsdat, 'AL', keepNthVertical)
    end if

  end subroutine thn_thinAladin


!_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
!_/
!_/ The following methods are intended to be general algorithms that may be
!_/ called by any of the observation-type-specific thinning methods.
!_/
!_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/


!--------------------------------------------------------------------------
!!
!! *Purpose*: Of the observations in a column that have not already been
!!            rejected, keep every nth observation and throw out the rest.
!!
!!            Set bit 11 of OBS_FLG on observations that are to be rejected.
!!
!--------------------------------------------------------------------------
  subroutine keepNthObs(obsdat, familyType, keepNthVertical)
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in) :: familyType
    integer,          intent(in) :: keepNthVertical

    integer, parameter :: BIT9 =int(Z'200')
    integer, parameter :: BIT11=int(Z'800')
    integer, parameter :: PROFILE_NOT_FOUND=-1
    integer :: headerIndex, bodyIndex, bodyIndex2, bodyIndexStart, bodyIndexEnd
    integer :: flag
    integer :: countKeepN ! count to keep every Nth observation in the column
    integer :: newProfileId
    real(OBS_REAL) :: level

    countKeepN=0

    ! Loop over all body indices (columns) of the family of interest and
    ! thin each column independently of the others
    call obs_set_current_body_list(obsdat, familyType)
    BODY: do 
      bodyIndex = obs_getBodyIndex(obsdat)
      if (bodyIndex < 0) exit BODY

      ! If the datum is not being assimilated, ignore it
      if(            obs_bodyElem_i(obsdat, OBS_ASS,  bodyIndex) == obs_notAssimilated) cycle BODY

      ! If datum already rejected, ignore it
      flag         = obs_bodyElem_i(obsdat, OBS_FLG,  bodyIndex  )
      if(btest(flag,BIT9))cycle BODY

      headerIndex  = obs_bodyElem_i(obsdat, OBS_HIND, bodyIndex  )
      newProfileId = obs_headElem_i(obsdat, OBS_PRFL, headerIndex)

      countKeepN=countKeepN + 1
      if(     countKeepN == keepNthVertical &
         .OR. new_column() &
        )then
        ! Reset the counter and keep this observation
        countKeepN=0
      else
        ! Reject this observation
        call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(flag, BIT11))
      end if

    end do BODY


  contains
    function new_column()
      ! Determine whether the current observation begins a new vertical column
      ! (Assume that observations are in chronological order)
      !
      ! Note:  This method has been written with aladin observations in mind.
      !        It might be necessary to generalize the method.
      implicit none
      logical :: new_column

      integer, save :: previousProfileId=huge(previousProfileId)

      if(newProfileId == PROFILE_NOT_FOUND)then
        ! The profile ID for this element is missing.
        ! Assume that it is the same as the previous element
        newProfileId = previousProfileId
      end if

      if(newProfileId /= previousProfileId)then
        previousProfileId = newProfileId
        new_column=.true.
      else
        new_column=.false.
      end if
    end function new_column
  end subroutine keepNthObs

end module thinning_mod

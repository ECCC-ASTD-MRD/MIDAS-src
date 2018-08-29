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
!! MODULE thinning (prefix="thn")
!!
!! *Purpose*: Using observation-type specific algorithms, set bit 11 of 'flag'
!!            on data that are not to be assimilated.
!!
!! *Note*:    This module is intended to group all of the thinning methods in a
!!            single fortran module.
!!
!!            So far, only aladin wind data are treated.
!!
!--------------------------------------------------------------------------
module thinning_mod
  use bufr_mod
  use obsSpaceData_mod
  use utilities_mod
  implicit none
  private
  public thn_thinAladin

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



    ! Read the namelist for Aladin observations
    nulnam = 0
    ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    if(ierr.ne.0) call utl_abort('midas-obsSelect: Error opening file flnml')
    read(nulnam,nml=thin_aladin,iostat=ierr)
    if(ierr.ne.0) call utl_abort('midas-obsSelect: Error reading namelist')
    write(*,nml=thin_aladin)
    ierr=fclos(nulnam)

    call keepNthObs(obsdat, 'AL', keepNthVertical)

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
    integer :: headerIndex, bodyIndex, bodyIndex2, bodyIndexStart, bodyIndexEnd
    integer :: flag
    integer :: countKeepN ! count to keep every Nth observation in the column
    integer :: newProfileId
    real :: level

    countKeepN=0

    ! Loop over all header indices of the family of interest
    call obs_set_current_header_list(obsdat, familyType)
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER

      ! Loop over all body indices (still in the same family)
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY

        ! if datum already rejected, ignore it
        flag           = obs_bodyElem_i(obsdat, OBS_FLG , bodyIndex  )
        if(btest(flag,BIT9))cycle BODY

        headerIndex    = obs_bodyElem_i(obsdat, OBS_HIND, bodyIndex  )
        bodyIndexStart = obs_headElem_i(obsdat, OBS_RLN , headerIndex)
        bodyIndexEnd =   obs_headElem_i(obsdat, OBS_NLV , headerIndex) &
                       + bodyIndexStart - 1
        level          = obs_bodyElem_r(obsdat, OBS_PPP , bodyIndex  )

        ! Obtain supplementary parameters that are stored as observations
        BODY_SUPP: do bodyIndex2 = bodyIndexStart, bodyIndexEnd
          if (       obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex2 ) == BUFR_NEPR &
              .and.  obs_bodyElem_r(obsdat, OBS_PPP, bodyIndex2 ) == level &
             ) then
            newProfileId = obs_bodyElem_i(obsdat, OBS_VAR, bodyIndex2 )
          end if
        end do BODY_SUPP

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
    end do HEADER


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

      if(newProfileId /= previousProfileId)then
        previousProfileId = newProfileId
        new_column=.true.
      else
        new_column=.false.
      end if
    end function new_column
  end subroutine keepNthObs

END module thinning_mod

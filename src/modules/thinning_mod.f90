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
  use mpi_mod
  use bufr_mod
  use mathPhysConstants_mod
  use obsSpaceData_mod
  use codtyp_mod
  use utilities_mod
  implicit none
  private

  public :: thn_thinAladin, thn_thinIASI

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

    ! Read the namelist for Aladin observations (if it exists)
    if (utl_isNamelistPresent('thin_aladin','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_thinAladin: Error opening file flnml')
      read(nulnam,nml=thin_aladin,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinAladin: Error reading namelist')
      if (mpi_myid == 0) write(*,nml=thin_aladin)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'thn_thinAladin: Namelist block thin_aladin is missing in the namelist.'
      write(*,*) '                The default value will be taken.'
    end if

    if (keepNthVertical > 0) then
      call thn_keepNthObs(obsdat, 'AL', keepNthVertical)
    end if

  end subroutine thn_thinAladin

  
  subroutine thn_thinIASI(obsdat)
    implicit none

    ! ARGUMENTS
    type(struct_obs), intent(inout) :: obsdat

    ! NAMELIST VARIABLES
    integer :: deltmax ! time window by bin (from bin center to bin edge) (in minutes)
    integer :: deltax  ! thinning (dimension of box sides) (in km)
    integer :: deltrad ! radius around box center for chosen obs (in km)
    namelist /thin_iasi/deltmax, deltax, deltrad

    ! LOCAL VARIABLES
    integer :: nulnam
    integer :: fnom, fclos, ierr

    ! Default namelist values
    deltmax = 22
    deltax  = 150
    deltrad = 45

    ! Read the namelist for Aladin observations (if it exists)
    if (utl_isNamelistPresent('thin_iasi','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_thinIASI: Error opening file flnml')
      read(nulnam,nml=thin_iasi,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinIASI: Error reading namelist')
      write(*,nml=thin_iasi)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'thn_thinIASI: Namelist block thin_iasi is missing in the namelist.'
      write(*,*) '                The default value will be taken.'
    end if

    call thn_thinByLatLonBoxes(obsdat, deltmax, deltax, deltrad, 'TO',  &
                               codtyp_opt=codtyp_get_codtyp('iasi'))

  end subroutine thn_thinIASI

!_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
!_/
!_/ The following methods are intended to be general algorithms that may be
!_/ called by any of the observation-type-specific thinning methods.
!_/
!_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/


  subroutine thn_keepNthObs(obsdat, familyType, keepNthVertical)
    !
    ! :Purpose: Of the observations in a column that have not already been
    !           rejected, keep every nth observation and throw out the rest.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.
    !
    implicit none

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

    countKeepN=0

    ! Loop over all body indices (columns) of the family of interest and
    ! thin each column independently of the others
    call obs_set_current_body_list(obsdat, familyType)
    BODY: do 
      bodyIndex = obs_getBodyIndex(obsdat)
      if (bodyIndex < 0) exit BODY

      ! If the datum is not being assimilated, ignore it
      if (obs_bodyElem_i(obsdat, OBS_ASS, bodyIndex) == obs_notAssimilated) cycle BODY

      ! If datum already rejected, ignore it
      flag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
      if (btest(flag,BIT9)) cycle BODY

      headerIndex  = obs_bodyElem_i(obsdat, OBS_HIND, bodyIndex  )
      newProfileId = obs_headElem_i(obsdat, OBS_PRFL, headerIndex)

      countKeepN=countKeepN + 1
      if ( countKeepN == keepNthVertical .or. &
           new_column() ) then
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

      if (newProfileId == PROFILE_NOT_FOUND) then
        ! The profile ID for this element is missing.
        ! Assume that it is the same as the previous element
        newProfileId = previousProfileId
      end if

      if (newProfileId /= previousProfileId) then
        previousProfileId = newProfileId
        new_column=.true.
      else
        new_column=.false.
      end if
    end function new_column

  end subroutine thn_keepNthObs


  subroutine thn_thinByLatLonBoxes(obsdat, deltmax, deltax, deltrad, familyType,  &
                                   codtyp_opt)
    !
    ! :Purpose: Only keep the observation closest to the center of each
    !           lat-lon (and time) box.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    integer, intent(in)             :: deltmax
    integer, intent(in)             :: deltax
    integer, intent(in)             :: deltrad
    character(len=*), intent(in)    :: familyType
    integer, optional, intent(in)   :: codtyp_opt

    ! Locals:
    integer :: headerIndex, bodyIndex, obsDate, obsTime, obsFlag
    real(8) :: obsLatInDegrees, obsLonInDegrees 

    write(*,*) 'thn_thinByLatLonBoxes: Starting'

    ! loop over all header indices of the specified family
    call obs_set_current_header_list(obsdat,trim(familyType))
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER

      if (present(codtyp_opt)) then
        if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) /= codtyp_opt) cycle HEADER
      end if

      obsLonInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obs_headElem_r(obsdat, OBS_LON, headerIndex)
      obsLatInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obs_headElem_r(obsdat, OBS_LAT, headerIndex)
      obsDate = obs_headElem_i(obsdat, OBS_DAT, headerIndex)
      obsTime = obs_headElem_i(obsdat, OBS_ETM, headerIndex)

      ! loop over all body indices for this headerIndex
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY

        obsFlag = obs_headElem_i(obsdat, OBS_FLG, headerIndex)

      end do BODY

    end do HEADER

    write(*,*) 'thn_thinByLatLonBoxes: Finished'

  end subroutine thn_thinByLatLonBoxes

end module thinning_mod

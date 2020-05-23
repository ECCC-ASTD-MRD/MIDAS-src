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
  use timeCoord_mod
  use codtyp_mod
  use physicsFunctions_mod
  use utilities_mod
  implicit none
  private

  public :: thn_thinAladin, thn_thinIASI

contains

  !--------------------------------------------------------------------------
  ! thn_thinAladin
  !--------------------------------------------------------------------------
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

  !--------------------------------------------------------------------------
  ! thn_thinIASI
  !--------------------------------------------------------------------------
  subroutine thn_thinIASI(obsdat)
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat

    ! Namelist variables
    integer :: deltmax ! time window by bin (from bin center to bin edge) (in minutes)
    integer :: deltax  ! thinning (dimension of box sides) (in km)
    integer :: deltrad ! radius around box center for chosen obs (in km)
    namelist /thin_iasi/deltmax, deltax, deltrad

    ! Locals:
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
      write(*,*) '              The default value will be taken.'
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


  !--------------------------------------------------------------------------
  ! thn_keepNthObs
  !--------------------------------------------------------------------------
  subroutine thn_keepNthObs(obsdat, familyType, keepNthVertical)
    !
    ! :Purpose: Of the observations in a column that have not already been
    !           rejected, keep every nth observation and throw out the rest.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in) :: familyType
    integer,          intent(in) :: keepNthVertical

    ! Locals:
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
      if (btest(flag,9)) cycle BODY

      headerIndex  = obs_bodyElem_i(obsdat, OBS_HIND, bodyIndex  )
      newProfileId = obs_headElem_i(obsdat, OBS_PRFL, headerIndex)

      countKeepN=countKeepN + 1
      if ( countKeepN == keepNthVertical .or. &
           new_column() ) then
        ! Reset the counter and keep this observation
        countKeepN=0
      else
        ! Reject this observation
        call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(flag,11))
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

  !--------------------------------------------------------------------------
  ! thn_thinByLatLonBoxes
  !--------------------------------------------------------------------------
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
    integer :: nblat, nblon, latIndex, numChannels, delMinutes
    integer :: lonBinIndex, latBinIndex, timeBinIndex
    integer :: ierr, nsize, procIndex
    real(4) :: latr, length, distance
    real(8) :: lonBoxCenterInRad, latBoxCenterInRad
    real(8) :: obsLatInRad, obsLonInRad
    real(8) :: obsLatInDegrees, obsLonInDegrees, stepObsIndex
    integer, parameter :: LAT_LENGTH = 10000 
    integer, parameter :: LON_LENGTH = 40000 
    real(4), allocatable :: latdeg(:)
    integer, allocatable :: ngrd(:)
    integer, allocatable :: headerIndexKeep(:,:,:), numChannelsKeep(:,:,:)
    integer, allocatable :: allHeaderIndex(:,:,:,:), allNumChannels(:,:,:,:)
    integer, allocatable :: delMinutesKeep(:,:,:)
    integer, allocatable :: allDelMinutes(:,:,:,:)
    integer, allocatable :: procIndexKeep(:,:,:)
    real(4), allocatable :: distanceKeep(:,:,:)
    real(4), allocatable :: allDistance(:,:,:,:)
    logical, allocatable :: rejectThisHeader(:)
    logical :: keepThisObs

    write(*,*) 'thn_thinByLatLonBoxes: Starting'

    ! Initial setup
    nblat = nint( 2. * real(LAT_LENGTH) / real(deltax) )
    nblon = nint(      real(LON_LENGTH) / real(deltax) )
    allocate(headerIndexKeep(nblat,nblon,tim_nstepobs))
    allocate(numChannelsKeep(nblat,nblon,tim_nstepobs))
    allocate(distanceKeep(nblat,nblon,tim_nstepobs))
    allocate(delMinutesKeep(nblat,nblon,tim_nstepobs))
    allocate(latdeg(nblat))
    allocate(ngrd(nblat))
    headerIndexKeep(:,:,:) = -1
    numChannelsKeep(:,:,:) = 0
    distanceKeep(:,:,:)    = 0.0
    delMinutesKeep(:,:,:)  = deltmax
    latdeg(:)              = 0.0
    ngrd(:)                = 0

    allocate(allHeaderIndex(nblat,nblon,tim_nstepobs,mpi_nprocs))
    allocate(allNumChannels(nblat,nblon,tim_nstepobs,mpi_nprocs))
    allocate(allDistance(nblat,nblon,tim_nstepobs,mpi_nprocs))
    allocate(allDelMinutes(nblat,nblon,tim_nstepobs,mpi_nprocs))
    allocate(procIndexKeep(nblat,nblon,tim_nstepobs))
    procIndexKeep(:,:,:) = -1

    ! set spatial boxes properties
    ! latdeg(:) : latitude (deg) of northern side of the box
    ! ngrd(:)   : number of longitudinal boxes at this latitude
    do latIndex = 1, nblat
      latdeg(latIndex) = (latIndex*180./nblat) - 90.
      if ( latdeg(latIndex) <= 0.0 ) then
        latr = latdeg(latIndex) * MPC_PI_R8 / 180.
      else
        latr = latdeg(latIndex-1) * MPC_PI_R8 / 180.
      end if
      length = LON_LENGTH * cos(latr)
      ngrd(latIndex)   = nint(length/deltax)
    end do

    ! loop over all header indices of the specified family
    call obs_set_current_header_list(obsdat,trim(familyType))
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER

      if (present(codtyp_opt)) then
        if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) /= codtyp_opt) cycle HEADER
      end if

      obsLonInRad = obs_headElem_r(obsdat, OBS_LON, headerIndex)
      obsLatInRad = obs_headElem_r(obsdat, OBS_LAT, headerIndex)
      obsLonInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obsLonInRad
      obsLatInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obsLatInRad
      obsDate = obs_headElem_i(obsdat, OBS_DAT, headerIndex)
      obsTime = obs_headElem_i(obsdat, OBS_ETM, headerIndex)

      numChannels = 0

      ! loop over all body indices for this headerIndex
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY

        if (obs_bodyElem_i(obsdat, OBS_ASS, bodyIndex) == obs_assimilated) then
          numChannels = numChannels + 1
        end if

      end do BODY

      ! Determine the lat and lon bin indexes
      do latIndex = 1, nblat
        if ( obsLatInDegrees <= (latdeg(latIndex)+0.000001) ) then
          latBinIndex = latIndex
          exit
        end if
      end do
      lonBinIndex = int( obsLonInDegrees/(360.0/ngrd(latBinIndex)) ) + 1
      if ( lonBinIndex > ngrd(latBinIndex) ) lonBinIndex = ngrd(latBinIndex)

      ! Determine the time bin index
      call tim_getStepObsIndex(stepObsIndex, tim_getDatestamp(), &
                               obsDate, obsTime, tim_nstepobs)
      timeBinIndex = nint(stepObsIndex)
      delMinutes = nint(60.0 * tim_dstepobs * abs(real(timeBinIndex) - stepObsIndex))

      ! Determine distance from box center
      latBoxCenterInRad = MPC_RADIANS_PER_DEGREE_R8 * &
                          (latdeg(latBinIndex) - 0.5 * (180./nblat))
      lonBoxCenterInRad = MPC_RADIANS_PER_DEGREE_R8 * &
                          ((360. / ngrd(latBinIndex)) * (lonBinIndex - 0.5))
      distance = 1.0d-3 * phf_calcDistance(latBoxCenterInRad, lonBoxCenterInRad, &
                                       obsLatInRad, obsLonInRad )

      ! Apply thinning criteria
      keepThisObs = .false.

      ! keep if distance to box center smaller than limit and 
      ! time under limit fixed at input and maximise number of channels 
      if ( (distance < deltrad)                                              .and. &
           (numChannels >= numChannelsKeep(latBinIndex,lonBinIndex,timeBinIndex)) .and. &
           (delMinutes <= deltmax) ) keepThisObs = .true.
     
      ! keep the closest to bin central time
      if ( numChannels == numChannelsKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
        if ( headerIndexKeep(latBinIndex,lonBinIndex,timeBinIndex) == -1 ) then
          keepThisObs = .false.
        else
          if ( delMinutes > delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex) ) keepThisObs = .false.
         
          if ( delMinutes == delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
            if ( distance > distanceKeep(latBinIndex,lonBinIndex,timeBinIndex) ) keepThisObs = .false.
          end if

        end if
      end if

      ! save the observation if thinning criteria is satisfied
      if ( keepThisObs ) then
        headerIndexKeep(latBinIndex,lonBinIndex,timeBinIndex) = headerIndex
        distanceKeep(latBinIndex,lonBinIndex,timeBinIndex)    = distance
        delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex)  = delMinutes
        numChannelsKeep(latBinIndex,lonBinIndex,timeBinIndex) = numChannels
      end if

    end do HEADER

    ! communicate results to all other mpi tasks
    nsize = nblat * nblon * tim_nstepobs
    call rpn_comm_allgather(distanceKeep, nsize, 'mpi_real4',  &
                            allDistance,  nsize, 'mpi_real4', 'grid', ierr)
    call rpn_comm_allgather(delMinutesKeep, nsize, 'mpi_integer',  &
                            allDelMinutes,  nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(numChannelsKeep, nsize, 'mpi_integer',  &
                            allNumChannels,  nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(headerIndexKeep, nsize, 'mpi_integer',  &
                            allheaderIndex,  nsize, 'mpi_integer', 'grid', ierr)

    ! reset arrays that store info about kept obs
    headerIndexKeep(:,:,:) = -1
    numChannelsKeep(:,:,:) = 0
    distanceKeep(:,:,:)    = 0.0
    delMinutesKeep(:,:,:)  = deltmax

    do timeBinIndex = 1, tim_nstepobs
      do lonBinIndex = 1, nblon
        do latBinIndex = 1, nblat

          ! Apply thinning criteria to results from all mpi tasks
          do procIndex = 1, mpi_nprocs

            headerIndex = allHeaderIndex(latBinIndex,lonBinIndex,timeBinIndex,procIndex)
            distance    = allDistance(latBinIndex,lonBinIndex,timeBinIndex,procIndex)
            delMinutes  = allDelMinutes(latBinIndex,lonBinIndex,timeBinIndex,procIndex)
            numChannels = allNumChannels(latBinIndex,lonBinIndex,timeBinIndex,procIndex)
            
            keepThisObs = .false.

            ! keep if distance to box center smaller than limit and 
            ! time under limit fixed at input and maximise number of channels 
            if ( (distance < deltrad)                                              .and. &
                 (numChannels >= numChannelsKeep(latBinIndex,lonBinIndex,timeBinIndex)) .and. &
                 (delMinutes <= deltmax) ) keepThisObs = .true.
     
            ! keep the closest to bin central time
            if ( numChannels == numChannelsKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
              if ( headerIndexKeep(latBinIndex,lonBinIndex,timeBinIndex) == -1 ) then
                keepThisObs = .false.
              else
                if ( delMinutes > delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex) ) keepThisObs = .false.
         
                if ( delMinutes == delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
                  if ( distance > distanceKeep(latBinIndex,lonBinIndex,timeBinIndex) ) keepThisObs = .false.
                end if

              end if
            end if

            ! save the observation if thinning criteria is satisfied
            if ( keepThisObs ) then
              procIndexKeep(latBinIndex,lonBinIndex,timeBinIndex)   = procIndex
              headerIndexKeep(latBinIndex,lonBinIndex,timeBinIndex) = headerIndex
              distanceKeep(latBinIndex,lonBinIndex,timeBinIndex)    = distance
              delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex)  = delMinutes
              numChannelsKeep(latBinIndex,lonBinIndex,timeBinIndex) = numChannels
            end if

          end do ! procIndex

        end do
      end do
    end do

    ! determine which headerIndex values are rejected
    allocate(rejectThisHeader(obs_numheader(obsdat)))
    rejectThisHeader(:) = .true.
    do timeBinIndex = 1, tim_nstepobs
      do lonBinIndex = 1, nblon
        do latBinIndex = 1, nblat
          if (procIndexKeep(latBinIndex,lonBinIndex,timeBinIndex) == mpi_myid+1) then
            headerIndex = headerIndexKeep(latBinIndex,lonBinIndex,timeBinIndex)
            rejectThisHeader(headerIndex) = .false.
          end if
        end do
      end do
    end do

    ! modify flags, by setting bit 11 for those rejected
    call obs_set_current_header_list(obsdat,trim(familyType))
    HEADER2: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER2

      if (present(codtyp_opt)) then
        if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) /= codtyp_opt) cycle HEADER2
      end if

      if (.not. rejectThisHeader(headerIndex)) cycle HEADER2

      ! loop over all body indices for this headerIndex
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY2: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY2

        obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
        call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,11))
      end do BODY2

    end do HEADER2

    write(*,*) 'thn_thinByLatLonBoxes: Finished'

  end subroutine thn_thinByLatLonBoxes

end module thinning_mod

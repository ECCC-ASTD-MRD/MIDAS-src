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
  use mpi_mod
  use bufr_mod
  use mathPhysConstants_mod
  use earthConstants_mod
  use obsSpaceData_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use gridStateVector_mod
  use codtyp_mod
  use physicsFunctions_mod
  use utilities_mod
  use kdtree2_mod
  implicit none
  private

  public :: thn_thinHyper, thn_thinTovs, thn_thinCSR
  public :: thn_thinRaobs, thn_thinAircraft, thn_thinScat, thn_thinSatWinds
  public :: thn_thinGbGps, thn_thinGpsRo, thn_thinAladin

  integer, external :: get_max_rss

contains

  !--------------------------------------------------------------------------
  ! thn_thinRaobs
  !--------------------------------------------------------------------------
  subroutine thn_thinRaobs(obsdat)
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat

    ! Locals:
    integer :: nulnam
    integer :: fnom, fclos, ierr

    ! Namelist variables
    logical :: verticalThinningES !
    logical :: ecmwfRejetsES      !
    logical :: rejectTdZeroC      !

    namelist /thin_raobs/ verticalThinningES, ecmwfRejetsES, rejectTdZeroC

    ! return if no aircraft obs
    if (.not. obs_famExist(obsdat,'UA')) return

    ! Default values for namelist variables
    verticalThinningES = .true.
    ecmwfRejetsES = .true.
    rejectTdZeroC = .true.

    ! Read the namelist for Radiosonde observations (if it exists)
    if (utl_isNamelistPresent('thin_raobs','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_thinRaobs: Error opening file flnml')
      read(nulnam,nml=thin_raobs,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinRaobs: Error reading thin_raobs namelist')
      if (mpi_myid == 0) write(*,nml=thin_raobs)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'thn_thinRaobs: Namelist block thin_raobs is missing in the namelist.'
      write(*,*) '               The default value will be taken.'
      if (mpi_myid == 0) write(*,nml=thin_raobs)
    end if

    call thn_radiosonde(obsdat, verticalThinningES, ecmwfRejetsES, rejectTdZeroC)

  end subroutine thn_thinRaobs

  !--------------------------------------------------------------------------
  ! thn_thinAircraft
  !--------------------------------------------------------------------------
  subroutine thn_thinAircraft(obsdat)
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat

    ! Locals:
    integer :: nulnam
    integer :: fnom, fclos, ierr

    ! Namelist variables
    integer :: deltmax ! maximum time difference (in minutes)

    namelist /thin_aircraft/deltmax

    ! return if no aircraft obs
    if (.not. obs_famExist(obsdat,'AI')) return

    ! Default values for namelist variables
    deltmax = 90

    ! Read the namelist for Aircraft observations (if it exists)
    if (utl_isNamelistPresent('thin_aircraft','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_thinAircraft: Error opening file flnml')
      read(nulnam,nml=thin_aircraft,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinAircraft: Error reading thin_aircraft namelist')
      if (mpi_myid == 0) write(*,nml=thin_aircraft)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'thn_thinAircraft: Namelist block thin_aircraft is missing in the namelist.'
      write(*,*) '                  The default value will be taken.'
      if (mpi_myid == 0) write(*,nml=thin_aircraft)
    end if

    call thn_aircraftByBoxes(obsdat, 'AI', deltmax)

  end subroutine thn_thinAircraft

  !--------------------------------------------------------------------------
  ! thn_thinSatWinds
  !--------------------------------------------------------------------------
  subroutine thn_thinSatWinds(obsdat)
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat

    ! Locals:
    integer :: nulnam
    integer :: fnom, fclos, ierr

    ! Namelist variables
    integer :: deltemps ! number of time bins between adjacent observations
    integer :: deldist  ! minimal distance in km between adjacent observations

    namelist /thin_satwind/deltemps, deldist

    ! return if no satwind obs
    if (.not. obs_famExist(obsdat,'SW')) return

    ! Default values for namelist variables
    deltemps = 6
    deldist  = 200

    ! Read the namelist for SatWinds observations (if it exists)
    if (utl_isNamelistPresent('thin_satwind','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_thinSatWinds: Error opening file flnml')
      read(nulnam,nml=thin_satwind,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinSatWinds: Error reading thin_satwind namelist')
      if (mpi_myid == 0) write(*,nml=thin_satwind)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'thn_thinSatWinds: Namelist block thin_satwind is missing in the namelist.'
      write(*,*) '                  The default value will be taken.'
      if (mpi_myid == 0) write(*,nml=thin_satwind)
    end if

    call thn_satWindsByDistance(obsdat, 'SW', deltemps, deldist)

  end subroutine thn_thinSatWinds

  !--------------------------------------------------------------------------
  ! thn_thinGpsRo
  !--------------------------------------------------------------------------
  subroutine thn_thinGpsRo(obsdat)
    ! :Purpose: Main thinning subroutine GPS radio-occultation obs.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat

    ! Locals:
    integer :: nulnam
    integer :: fnom, fclos, ierr

    ! Namelist variables
    real(8) :: heightMin     ! niveau a partir du quel on accepte les donnees
    real(8) :: heightMax     ! niveau a partir du quel on rejette les donnees
    real(8) :: heightSpacing ! epaisseur minimale entre deux niveaux

    namelist /thin_gpsro/heightMin, heightMax, heightSpacing

    ! return if no gb-gps obs
    if (.not. obs_famExist(obsdat,'RO')) return

    ! Default values for namelist variables
    heightMin     = 1000.0d0
    heightMax     = 40000.0d0
    heightSpacing = 750.0d0

    ! Read the namelist for GpsRo observations (if it exists)
    if (utl_isNamelistPresent('thin_gpsro','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_thinGpsRo: Error opening file flnml')
      read(nulnam,nml=thin_gpsro,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinGpsRo: Error reading thin_gpsro namelist')
      if (mpi_myid == 0) write(*,nml=thin_gpsro)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'thn_thinGpsRo: Namelist block thin_gpsro is missing in the namelist.'
      write(*,*) '               The default value will be taken.'
      if (mpi_myid == 0) write(*,nml=thin_gpsro)
    end if

    call thn_gpsroVertical(obsdat, heightMin, heightMax, heightSpacing)

  end subroutine thn_thinGpsRo

  !--------------------------------------------------------------------------
  ! thn_thinGbGps
  !--------------------------------------------------------------------------
  subroutine thn_thinGbGps(obsdat)
    ! :Purpose: Main thinning subroutine ground-based GPS obs.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat

    ! Locals:
    integer :: nulnam
    integer :: fnom, fclos, ierr

    ! Namelist variables
    integer :: deltemps     ! number of time bins between adjacent observations
    integer :: deldist      ! minimal distance in km between adjacent observations
    logical :: removeUncorrected ! remove obs that are not bias corrected (bit 6)

    namelist /thin_gbgps/deltemps, deldist, removeUncorrected

    ! return if no gb-gps obs
    if (.not. obs_famExist(obsdat,'GP')) return

    ! Default values for namelist variables
    deltemps     = 8
    deldist      = 50
    removeUncorrected = .false.

    ! Read the namelist for GbGps observations (if it exists)
    if (utl_isNamelistPresent('thin_gbgps','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_thinGbGps: Error opening file flnml')
      read(nulnam,nml=thin_gbgps,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinGbGps: Error reading thin_gbgps namelist')
      if (mpi_myid == 0) write(*,nml=thin_gbgps)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'thn_thinGbGps: Namelist block thin_gbgps is missing in the namelist.'
      write(*,*) '               The default value will be taken.'
      if (mpi_myid == 0) write(*,nml=thin_gbgps)
    end if

    call thn_gbgpsByDistance(obsdat, deltemps, deldist, removeUncorrected)

  end subroutine thn_thinGbGps

  !--------------------------------------------------------------------------
  ! thn_thinAladin
  !--------------------------------------------------------------------------
  subroutine thn_thinAladin(obsdat)
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat

    ! Locals:
    integer :: nulnam
    integer :: fnom, fclos, ierr

    ! Namelist variables
    integer :: keepNthVertical ! keep every nth vertical datum

    namelist /thin_aladin/keepNthVertical

    ! return if no Aladin obs
    if (.not. obs_famExist(obsdat,'AL')) return

    ! Default values for namelist variables
    keepNthVertical=-1

    ! Read the namelist for Aladin observations (if it exists)
    if (utl_isNamelistPresent('thin_aladin','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_thinAladin: Error opening file flnml')
      read(nulnam,nml=thin_aladin,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinAladin: Error reading thin_aladin namelist')
      if (mpi_myid == 0) write(*,nml=thin_aladin)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'thn_thinAladin: Namelist block thin_aladin is missing in the namelist.'
      write(*,*) '                The default value will be taken.'
      if (mpi_myid == 0) write(*,nml=thin_aladin)
    end if

    if (keepNthVertical > 0) then
      call thn_keepNthObs(obsdat, 'AL', keepNthVertical)
    end if

  end subroutine thn_thinAladin

  !--------------------------------------------------------------------------
  ! thn_thinCSR
  !--------------------------------------------------------------------------
  subroutine thn_thinCSR(obsdat)
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat

    ! Locals:
    integer :: nulnam
    integer :: fnom, fclos, ierr

    ! Namelist variables
    integer :: deltax     ! thinning (dimension of box sides) (in km)
    integer :: deltrad    ! radius around box center for chosen obs (in km)

    namelist /thin_csr/deltax, deltrad

    ! return if no TOVS obs
    if (.not. obs_famExist(obsdat,'TO')) return

    ! Default namelist values
    deltax  = 150
    deltrad = 45

    ! Read the namelist for CSR observations (if it exists)
    if (utl_isNamelistPresent('thin_csr','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_thinCSR: Error opening file flnml')
      read(nulnam,nml=thin_csr,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinCSR: Error reading thin_csr namelist')
      if (mpi_myid == 0) write(*,nml=thin_csr)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'thn_thinCSR: Namelist block thin_csr is missing in the namelist.'
      write(*,*) '             The default value will be taken.'
      if (mpi_myid == 0) write(*,nml=thin_csr)
    end if

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_csrByLatLonBoxes(obsdat, deltax, deltrad)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine thn_thinCSR

  !--------------------------------------------------------------------------
  ! thn_thinScat
  !--------------------------------------------------------------------------
  subroutine thn_thinScat(obsdat)
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat

    ! Locals:
    integer :: nulnam
    integer :: fnom, fclos, ierr

    ! Namelist variables
    integer :: deltax     ! thinning (dimension of box sides) (in km)
    integer :: deltmax    ! temporal thinning resolution (in minutes)

    namelist /thin_scat/deltax, deltmax

    ! return if no scat obs
    if (.not. obs_famExist(obsdat,'SC')) return

    ! Default namelist values
    deltax = 100
    deltmax = 90

    ! Read the namelist for Scat observations (if it exists)
    if (utl_isNamelistPresent('thin_scat','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_thinScat: Error opening file flnml')
      read(nulnam,nml=thin_scat,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinScat: Error reading thin_scat namelist')
      if (mpi_myid == 0) write(*,nml=thin_scat)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'thn_thinScat: Namelist block thin_scat is missing in the namelist.'
      write(*,*) '              The default value will be taken.'
      if (mpi_myid == 0) write(*,nml=thin_scat)
    end if

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_scatByLatLonBoxes(obsdat, deltax, deltmax)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine thn_thinScat

  !--------------------------------------------------------------------------
  ! thn_thinTovs
  !--------------------------------------------------------------------------
  subroutine thn_thinTovs(obsdat)
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat

    ! Locals:
    integer :: nulnam
    integer :: fnom, fclos, ierr

    ! Namelist variables
    integer :: delta    ! 

    namelist /thin_tovs/delta

    ! return if no TOVS obs
    if (.not. obs_famExist(obsdat,'TO')) return

    ! Default namelist values
    delta = 100

    ! Read the namelist for TOVS observations (if it exists)
    if (utl_isNamelistPresent('thin_tovs','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_thinTovs: Error opening file flnml')
      read(nulnam,nml=thin_tovs,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinTovs: Error reading thin_tovs namelist')
      if (mpi_myid == 0) write(*,nml=thin_tovs)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'thn_thinTovs: Namelist block thin_tovs is missing in the namelist.'
      write(*,*) '              The default value will be taken.'
      if (mpi_myid == 0) write(*,nml=thin_tovs)
    end if

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_tovsFilt(obsdat, delta, codtyp_get_codtyp('amsua'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_tovsFilt(obsdat, delta, codtyp_get_codtyp('amsub'), &
                      codtyp2_opt=codtyp_get_codtyp('mhs'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_tovsFilt(obsdat, delta, codtyp_get_codtyp('atms'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_tovsFilt(obsdat, delta, codtyp_get_codtyp('mwhs2'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine thn_thinTovs

  !--------------------------------------------------------------------------
  ! thn_thinHyper
  !--------------------------------------------------------------------------
  subroutine thn_thinHyper(obsdat)
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat

    ! Locals:
    integer :: nulnam
    integer :: fnom, fclos, ierr

    ! Namelist variables
    logical :: removeUnCorrected ! indicate if obs without bias correction should be removed
    integer :: deltmax           ! time window by bin (from bin center to bin edge) (in minutes)
    integer :: deltax            ! thinning (dimension of box sides) (in km)
    integer :: deltrad           ! radius around box center for chosen obs (in km)
    namelist /thin_hyper/removeUnCorrected, deltmax, deltax, deltrad

    ! return if no TOVS obs
    if (.not. obs_famExist(obsdat,'TO')) return

    ! Default namelist values
    removeUnCorrected = .true.
    deltmax = 22
    deltax  = 150
    deltrad = 45

    ! Read the namelist for Aladin observations (if it exists)
    if (utl_isNamelistPresent('thin_hyper','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_thinHyper: Error opening file flnml')
      read(nulnam,nml=thin_hyper,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinHyper: Error reading thin_hyper namelist')
      if (mpi_myid == 0) write(*,nml=thin_hyper)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'thn_thinHyper: Namelist block thin_hyper is missing in the namelist.'
      write(*,*) '               The default value will be taken.'
      if (mpi_myid == 0) write(*,nml=thin_hyper)
    end if

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_hyperByLatLonBoxes(obsdat, removeUnCorrected, deltmax, deltax, deltrad, &
                                'TO', codtyp_get_codtyp('airs'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_hyperByLatLonBoxes(obsdat, removeUnCorrected, deltmax, deltax, deltrad, &
                                'TO', codtyp_get_codtyp('iasi'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_hyperByLatLonBoxes(obsdat, removeUnCorrected, deltmax, deltax, deltrad, &
                                'TO', codtyp_get_codtyp('cris'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_hyperByLatLonBoxes(obsdat, removeUnCorrected, deltmax, deltax, deltrad, &
                                'TO', codtyp_get_codtyp('crisfsr'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine thn_thinHyper

!_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
!_/
!_/ The following methods are intended to be general algorithms that may be
!_/ called by any of the observation-type-specific thinning methods.
!_/
!_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  !--------------------------------------------------------------------------
  ! thn_gpsroVertical
  !--------------------------------------------------------------------------
  subroutine thn_gpsroVertical(obsdat, heightMin, heightMax, heightSpacing)
    !
    ! :Purpose: Original method for thinning GPSRO data by vertical distance.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    real(8),          intent(in)    :: heightMin
    real(8),          intent(in)    :: heightMax
    real(8),          intent(in)    :: heightSpacing

    ! Local parameters:
    integer, parameter :: gpsroVarNo = BUFR_NERF

    ! Locals:
    integer :: countObs, countObsInMpi, headerIndex, bodyIndex, ierr
    integer :: numLev, levIndex, obsVarNo, obsFlag
    real(8) :: nextHeightMin
    logical :: rejectObs
    integer, allocatable :: bodyIndexList(:)
    real(8), allocatable :: obsHeights(:)

    write(*,*)
    write(*,*) 'thn_gpsroVertical: Starting'
    write(*,*)

    ! Check if any observations to be treated
    countObs = 0
    call obs_set_current_header_list(obsdat,'RO')
    HEADER0: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER0
      countObs = countObs + 1
    end do HEADER0

    call rpn_comm_allReduce(countObs, countObsInMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    if (countObsInMpi == 0) then
      write(*,*) 'thn_gpsroVertical: no gpsro observations present'
      return
    end if

    write(*,*) 'thn_gpsroVertical: number of obs initial = ', &
               countObs, countObsInMpi

    call obs_set_current_header_list(obsdat,'RO')
    HEADER1: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER1

      ! count number of levels for this profile
      numLev = 0
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY1: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY1

        obsVarNo = obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex)
        if (obsVarNo /= gpsroVarNo) cycle BODY1

        numLev = numLev + 1

      end do BODY1

      allocate(obsHeights(numLev))
      allocate(bodyIndexList(numLev))

      ! extract altitudes for this profile
      levIndex = 0
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY2: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY2
        
        obsVarNo = obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex)
        if (obsVarNo /= gpsroVarNo) cycle BODY2

        levIndex = levIndex + 1
        obsHeights(levIndex) = obs_bodyElem_r(obsdat, OBS_PPP, bodyIndex)
        bodyIndexList(levIndex) = bodyIndex

      end do BODY2

      ! ensure altitudes are in ascending order
      call thn_QsortReal8(obsHeights,bodyIndexList)

      ! apply vertical thinning
      nextHeightMin = heightMin
      LEVELS: do levIndex = 1, numLev
        
        if ( obsHeights(levIndex) >= nextHeightMin .and. &
             obsHeights(levIndex) < heightMax ) then
          nextHeightMin = obsHeights(levIndex) + heightSpacing
          rejectObs = .false.
        else
          rejectObs = .true.
        end if

        if (rejectObs) then
          bodyIndex = bodyIndexList(levIndex)
          obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
          call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,11))
        end if
        
      end do LEVELS

      deallocate(obsHeights)
      deallocate(bodyIndexList)

    end do HEADER1

    write(*,*)
    write(*,*) 'thn_gpsroVertical: Finished'
    write(*,*)

  end subroutine thn_gpsroVertical

  !--------------------------------------------------------------------------
  ! thn_radiosonde
  !--------------------------------------------------------------------------
  subroutine thn_radiosonde(obsdat, verticalThinningES, ecmwfRejetsES, rejectTdZeroC)
    !
    ! :Purpose: Original method for thinning radiosonde data vertically.
    !           We assume that each vertical level is stored in obsSpaceData
    !           with a separate headerIndex. That is, the 4D representation.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    logical,          intent(in)    :: verticalThinningES
    logical,          intent(in)    :: ecmwfRejetsES
    logical,          intent(in)    :: rejectTdZeroC

    ! Locals:
    type(struct_hco), pointer :: hco_sfc
    type(struct_vco), pointer :: vco_sfc
    type(struct_gsv)          :: stateVectorPsfc
    integer :: nulnam, fnom, fclos, ezgdef, ezsint, ezdefset, ezsetopt
    integer :: ierr, nltot, nltotMpi, levStnIndex, countLevel,nlev_max
    integer :: nbstn, nbstnMpi, countProfile, lastProfileIndex
    integer :: profileIndex, headerIndex, bodyIndex, levIndex, stepIndex
    integer :: ig1obs, ig2obs, ig3obs, ig4obs, lalo
    integer :: obsFlag
    real(4) :: obsValue, obsOmp, obsStepIndex
    real(4) :: zig1, zig2, zig3, zig4, zpresa, zpresb
    real(8) :: obsStepIndex_r8
    integer, allocatable :: niv_stn(:), stn_type(:), date(:), temps(:)
    integer, allocatable :: date_ini(:), temps_ini(:), temps_lch(:), flgs_h(:)
    integer, allocatable :: flgs_t(:,:), flgs_o(:,:)
    real(4), allocatable :: lat_stn(:), lon_stn(:), int_surf(:,:)
    real(4), allocatable :: obsv(:,:), ombv(:,:), mod_pp(:,:)
    real(4), pointer     :: surfPressure(:,:,:,:)
    character(len=9), allocatable :: ids_stn(:)
    character(len=20) :: trlmFileName
    character(len=2)  :: fileNumber
    logical :: upperAirObs

    integer :: iniv,istn
    integer :: countAcc_dd, countRej_dd, countAccMpi_dd, countRejMpi_dd
    integer :: countAcc_ff, countRej_ff, countAccMpi_ff, countRejMpi_ff
    integer :: countAcc_tt, countRej_tt, countAccMpi_tt, countRejMpi_tt
    integer :: countAcc_es, countRej_es, countAccMpi_es, countRejMpi_es

    ! Local parameters:
    integer, parameter :: nbvar=5, nbtrj=2, nbmod=300

    ! Namelist variables:
    real(8) :: rprefinc
    real(8) :: rptopinc
    real(8) :: rcoefinc
    real(4) :: vlev(nbmod)
    integer :: numlev
    namelist /namgem/rprefinc, rptopinc, rcoefinc, numlev, vlev

    ! Check if any observations to be treated and count number of "profiles"
    nltot = 0
    nbstn = 0
    lastProfileIndex = -1
    call obs_set_current_header_list(obsdat,'UA')
    HEADER0: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER0

      ! skip if this headerIndex doesn't contain upper air obs
      upperAirObs = .false.
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY0: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY0

        select case (obs_bodyElem_i(obsdat,OBS_VNM,bodyIndex))
        case (bufr_neuu, bufr_nevv, bufr_nett, bufr_nees)
          upperAirObs = .true.
        end select
      end do BODY0
      profileIndex = obs_headElem_i(obsdat, obs_prfl, headerIndex)
      if (.not. upperAirObs) cycle HEADER0

      nltot = nltot + 1

      profileIndex = obs_headElem_i(obsdat, obs_prfl, headerIndex)
      if (profileIndex /= lastProfileIndex) then
        lastProfileIndex = profileIndex
        nbstn = nbstn + 1
      end if

    end do HEADER0

    call rpn_comm_allReduce(nltot, nltotMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    if (nltotMpi == 0) then
      write(*,*) 'thn_radiosonde: no UA obs observations present'
      return
    end if

    call rpn_comm_allReduce(nbstn, nbstnMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)

    write(*,*) 'thn_radiosonde: number of obs initial = ', &
               nltot, nltotMpi
    write(*,*) 'thn_radiosonde: number of profiles    = ', &
               nbstn, nbstnMpi

    ! Allocate some quanitities needed for each profile
    allocate(niv_stn(nbstn+1))
    allocate(stn_type(nbstn))
    allocate(date(nbstn))
    allocate(temps(nbstn))
    allocate(date_ini(nbstn))
    allocate(temps_ini(nbstn))
    allocate(temps_lch(nbstn))
    allocate(flgs_h(nbstn))
    allocate(lat_stn(nbstn))
    allocate(lon_stn(nbstn))
    allocate(ids_stn(nbstn))
    allocate(mod_pp(nbstn,nbmod))
    allocate(flgs_t(nbtrj,nltot))
    allocate(flgs_o(nbvar,nltot))
    allocate(obsv(nbvar,nltot))
    allocate(ombv(nbvar,nltot))
    niv_stn(:) = 0
    date_ini(:) = -1
    flgs_t(:,:) = -1
    flgs_o(:,:) = 0
    obsv(:,:) = -999.0
    ombv(:,:) = -999.0

    ! Fill in the arrays for each profile
    levStnIndex = 0
    countProfile = 0
    lastProfileIndex = -1
    nlev_max = -1
    call obs_set_current_header_list(obsdat,'UA')
    HEADER1: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER1

      ! skip if this headerIndex doesn't contain upper air obs
      upperAirObs = .false.
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY1: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY1

        select case (obs_bodyElem_i(obsdat,OBS_VNM,bodyIndex))
        case (bufr_neuu, bufr_nevv, bufr_nett, bufr_nees)
          upperAirObs = .true.
        end select
      end do BODY1
      if (.not. upperAirObs) cycle HEADER1

      levStnIndex = levStnIndex + 1

      profileIndex = obs_headElem_i(obsdat, obs_prfl, headerIndex)
      if (profileIndex /= lastProfileIndex) then
        lastProfileIndex = profileIndex
        countProfile = countProfile + 1
        countLevel = 0
      end if

      countLevel = countLevel + 1
      niv_stn(countProfile+1) = niv_stn(countProfile) + countLevel

      if (countLevel > nlev_max) nlev_max = countLevel

      ! Get some information from the first header in this profile
      if (countLevel == 1) then
        date_ini(countProfile)  = obs_headElem_i(obsdat,obs_dat,headerIndex)
        temps_ini(countProfile) = obs_headElem_i(obsdat,obs_etm,headerIndex)
        lat_stn(countProfile)   = obs_headElem_r(obsdat,obs_lat,headerIndex) * &
                                  MPC_DEGREES_PER_RADIAN_R8
        lon_stn(countProfile)   = obs_headElem_r(obsdat,obs_lon,headerIndex) * &
                                  MPC_DEGREES_PER_RADIAN_R8
        lat_stn(countProfile)   = 0.01*nint(100.0*lat_stn(countProfile))
        lon_stn(countProfile)   = 0.01*nint(100.0*lon_stn(countProfile))
      end if
      date(countProfile)      = obs_headElem_i(obsdat,obs_hdd,headerIndex)
      temps(countProfile)     = obs_headElem_i(obsdat,obs_hdt,headerIndex)
      temps_lch(countProfile) = obs_headElem_i(obsdat,obs_lch,headerIndex)
      stn_type(countProfile)  = obs_headElem_i(obsdat,obs_rtp,headerIndex)
      flgs_h(countProfile)    = obs_headElem_i(obsdat,obs_st1,headerIndex)
      ids_stn(countProfile)   = obs_elem_c(obsdat,'STID',headerIndex)

      call obs_set_current_body_list(obsdat, headerIndex)
      BODY2: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY2

        obsFlag  = obs_bodyElem_i(obsdat,obs_flg,bodyIndex)
        obsValue = obs_bodyElem_r(obsdat,obs_var,bodyIndex)
        obsOmp   = obs_bodyElem_r(obsdat,obs_omp,bodyIndex)
        select case (obs_bodyElem_i(obsdat,OBS_VNM,bodyIndex))
        case (bufr_nedd)
          flgs_o(1,levStnIndex) = obsFlag
          obsv(1,levStnIndex)   = obsValue
        case (bufr_neuu)
          ombv(1,levStnIndex)   = obsOmp
        case (bufr_neff)
          flgs_o(2,levStnIndex) = obsFlag
          obsv(2,levStnIndex)   = obsValue
        case (bufr_nevv)
          ombv(2,levStnIndex)   = obsOmp
        case (bufr_nett)
          flgs_o(3,levStnIndex) = obsFlag
          obsv(3,levStnIndex)   = obsValue
          ombv(3,levStnIndex)   = obsOmp
        case (bufr_nees)
          flgs_o(4,levStnIndex) = obsFlag
          obsv(4,levStnIndex)   = obsValue
          ombv(4,levStnIndex)   = obsOmp
        case (4015)
          flgs_t(1,levStnIndex) = obsFlag
        case (5001)
          flgs_t(2,levStnIndex) = obsFlag
        end select
        obsv(5,levStnIndex) = obs_bodyElem_r(obsdat,obs_ppp,bodyIndex)*0.01
      end do BODY2

    end do HEADER1

    ! Default namelist values
    numlev = 80
    vlev(:) = -1
    rprefinc = 0.0d0
    rptopinc = 0.0d0
    rcoefinc = 0.0d0
    ! Read the namelist defining the vertical levels
    if (utl_isNamelistPresent('namgem','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_radiosonde: Error opening file flnml')
      read(nulnam,nml=namgem,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_radiosonde: Error reading namgem namelist')
      if (mpi_myid == 0) write(*,nml=namgem)
      ierr = fclos(nulnam)
    else
      call utl_abort('thn_radiosonde: Namelist block namgem is missing in the namelist.')
    end if

    ! Read the trial surface pressure
    nullify(vco_sfc)
    nullify(hco_sfc)
    trlmFileName = './trlm_01'
    call vco_setupFromFile(vco_sfc, trlmFileName)
    call hco_setupFromFile(hco_sfc, trlmFileName, ' ')
    call gsv_allocate( stateVectorPsfc, tim_nstepobs, hco_sfc, vco_sfc, &
                       datestamp_opt=tim_getDatestamp(), mpi_local_opt=.false., &
                       dataKind_opt=4, varNames_opt=(/'P0'/), &
                       hInterpolateDegree_opt='LINEAR' )
    do stepIndex = 1, tim_nstepobs
      write(fileNumber,'(I2.2)') stepIndex
      trlmFileName = './trlm_' // trim(fileNumber)
      call gsv_readFromFile( stateVectorPsfc, trlmFileName, ' ', ' ',  &
                             stepIndex_opt=stepIndex, containsFullField_opt=.true. )
    end do
    call gsv_getField(stateVectorPsfc,surfPressure)

    ! Setup the interpolation to obs locations of surface pressure
    zig1 = 0.0
    zig2 = 0.0
    zig3 = 1.0
    zig4 = 1.0
    ierr = ezsetopt('INTERP_DEGREE', 'LINEAR')
    call cxgaig('L',ig1obs,ig2obs,ig3obs,ig4obs,zig1,zig2,zig3,zig4)
    lalo = ezgdef(nbstn,1,'Y','L',ig1obs,ig2obs,ig3obs,ig4obs,lon_stn,lat_stn)
    ierr = ezdefset(lalo,hco_sfc%ezScintId)

    ! Do the interpolation of surface pressure to obs locations
    allocate(int_surf(nbstn,tim_nstepobs))
    do stepIndex = 1, tim_nstepobs
      ierr = ezsint(int_surf(:,stepIndex),surfPressure(:,:,1,stepIndex))
    end do
    call gsv_deallocate( stateVectorPsfc )

    ! Compute pressure profile at each obs location
    do countProfile = 1, nbstn

      ! Calculate the stepIndex corresponding to the launch time
      call tim_getStepObsIndex(obsStepIndex_r8, tim_getDatestamp(), &
                               date_ini(countProfile), temps_ini(countProfile), tim_nstepobs)
      obsStepIndex = nint(obsStepIndex_r8)

      ! Calculate pressure levels for each station based on GEM3 vertical coordinate
      do levIndex  = 1, numLev
        zpresb = ((vlev(levIndex) - rptopinc/rprefinc) / (1.0-rptopinc/rprefinc))**rcoefinc
        zpresa = rprefinc * (vlev(levIndex)-zpresb)
        mod_pp(countProfile,levIndex) = 0.01*(zpresa + zpresb*int_surf(countProfile,obsStepIndex))
      end do

    end do

    call check_duplicated_stations( ids_stn, niv_stn, lat_stn, lon_stn, date, &
                                    flgs_t, flgs_o, obsv, ombv, temps_lch, flgs_h, &
                                    nbvar, nltot, nbstn )
    
    call thinning_model( flgs_o, obsv, mod_pp, nbvar, numlev, nltot, &
                         nbstn, nlev_max, niv_stn )

    if( verticalThinningES ) then
      call thinning_es( flgs_o, obsv, mod_pp, numlev, nltot, nbstn, nlev_max, niv_stn )
    end if

    if( ecmwfRejetsES ) then
      call backlisting_ecmwf( flgs_o, obsv, ids_stn, stn_type, mod_pp, numlev, nltot, nbstn, niv_stn )
    end if

    countAcc_dd=0;  countRej_dd=0
    countAcc_ff=0;  countRej_ff=0
    countAcc_tt=0;  countRej_tt=0
    countAcc_es=0;  countRej_es=0
    do istn=1,nbstn
      do iniv=niv_stn(istn)+1,niv_stn(istn+1)
        if(btest(flgs_o(1,iniv),11)) then
          countRej_dd = countRej_dd + 1
        else
          countAcc_dd = countAcc_dd + 1
        end if
        if(btest(flgs_o(2,iniv),11)) then
          countRej_ff = countRej_ff + 1
        else
          countAcc_ff = countAcc_ff + 1
        end if
        if(btest(flgs_o(3,iniv),11)) then
          countRej_tt = countRej_tt + 1
        else
          countAcc_tt = countAcc_tt + 1
        end if
        if(btest(flgs_o(4,iniv),11) .or. btest(flgs_o(4,iniv),8)) then
          countRej_es = countRej_es + 1
        else
          countAcc_es = countAcc_es + 1
        end if
      end do
    end do

    call rpn_comm_allReduce(countRej_dd, countRejMpi_dd, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(countRej_ff, countRejMpi_ff, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(countRej_tt, countRejMpi_tt, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(countRej_es, countRejMpi_es, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(countAcc_dd, countAccMpi_dd, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(countAcc_ff, countAccMpi_ff, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(countAcc_tt, countAccMpi_tt, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(countAcc_es, countAccMpi_es, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)

    write(*,*) 'DD Rej/Acc = ',countRejMpi_dd, countAccMpi_dd
    write(*,*) 'FF Rej/Acc = ',countRejMpi_ff, countAccMpi_ff
    write(*,*) 'TT Rej/Acc = ',countRejMpi_tt, countAccMpi_tt
    write(*,*) 'ES Rej/Acc = ',countRejMpi_es, countAccMpi_es

    ! Deallocate arrays

  end subroutine thn_radiosonde

  !---------------------------------------------------------------------
  ! Following are several subroutines needed by thn_radiosonde
  !---------------------------------------------------------------------

  subroutine check_duplicated_stations ( ids_stn, niv_stn, lat_stn, lon_stn, date, &
                                         flgs_t, flgs_o, obsv, ombv, temps_lch, flgs_h, &
                                         nbvar, nltot, ntot_stn )
    ! Purpose : Check duplicated stations and select the best TAC/BUFR profiles
    !           
    ! Version 2.0
    !
    ! Author: S. Laroche (ARMA)  June 2013
    !
    ! Revision: 2.0 S. Laroche (ARMA) Decembre 2016
    !
    implicit none

    integer,           intent(in)    :: nbvar, nltot, ntot_stn
    integer,           intent(in)    :: date(:), temps_lch(:), flgs_h(:)
    real,              intent(in)    :: lat_stn(:), lon_stn(:)
    real,              intent(in)    :: obsv(:,:), ombv(:,:)
    integer,           intent(in)    :: flgs_t(:,:)
    integer,           intent(in)    :: flgs_o(:,:)
    integer,           intent(inout) :: niv_stn(:)
    character (len=9), intent(inout) :: ids_stn(:)

    logical :: condition, same_profile, stnid_not_found
    integer :: ilev, jlev, klev, istn, jstn, kstn, iele, iniv, nniv, grt_nval, pos_stn, duplicate, duplicate_total, select_ij 
    integer :: ij_bufr, ij_tac, n_checked, nb_stnid, nb_same
    integer, parameter :: max_nb_stnid  = 2000
    character (len=9)  :: stnid_found(max_nb_stnid)
    integer            :: stnid_numbr(max_nb_stnid)
    integer :: n_choix(5), cloche(30), select_critere

    n_choix(:) = 0
    cloche(:) = 0

    ! Check for duplication of TAC or BUFR profiles

    duplicate_total = 0

    do istn=1,ntot_stn

      if( ids_stn(istn) /= 'NOT_VALID' ) then

        grt_nval = niv_stn(istn+1) - niv_stn(istn) + 1
        pos_stn = istn
        duplicate = 0

        ! Verify if there exists two records with the same stnid, date, temps, flgs in header.
        ! Keep the one with the greatest number of vertical levels

        do jstn=istn+1,ntot_stn

          if( ids_stn(jstn) /= 'NOT_VALID' ) then

            condition = ids_stn(jstn) == ids_stn(istn) .and. &
                        date(jstn) == date(istn) .and. &
                        temps_lch(jstn) == temps_lch(istn) .and. &
                        lat_stn(jstn) == lat_stn(istn) .and. &
                        lon_stn(jstn) == lon_stn(istn) .and. &
                        flgs_h(jstn) == flgs_h(istn)

            if( condition ) then
              duplicate = duplicate + 1
              if( (niv_stn(jstn+1) - niv_stn(jstn) + 1) > grt_nval ) then
                grt_nval = niv_stn(jstn+1) - niv_stn(jstn) + 1
                pos_stn  = jstn
              end if
            end if

          end if

        end do ! jstn

        ! Invalid all duplicated station except the one with the greatest number of vertical levels

        if(duplicate > 0) then
          do jstn=istn,ntot_stn
            if( ids_stn(jstn) /= 'NOT_VALID' .and. ids_stn(jstn) == ids_stn(istn) .and. pos_stn /= jstn ) then
              write(*,'(2A20,I10,2F9.2,I10)') 'Station duplique ',ids_stn(jstn),niv_stn(jstn+1)-niv_stn(jstn), lat_stn(jstn), lon_stn(jstn), temps_lch(jstn)
              ids_stn(jstn) = 'NOT_VALID'
            end if
          end do
          duplicate_total = duplicate_total + duplicate
        end if

      end if

    end do ! istn

    ! Select best profile for collocated TAC or BUFR reports

    if (mpi_myid == 0) open(unit=11, file='./selected_stations_tac_burf.txt', status='UNKNOWN')

    n_checked = 0

    do istn=1,ntot_stn

      if( ids_stn(istn) /= 'NOT_VALID' ) then

        do jstn=istn+1,ntot_stn

          if( ids_stn(jstn) /= 'NOT_VALID' ) then

            condition = (ids_stn(jstn) == ids_stn(istn)) .and. (btest(flgs_h(istn),23) .neqv. btest(flgs_h(jstn),23))

            if( condition ) then

              if(temps_lch(jstn) == temps_lch(istn)) then
                same_profile = .true.
              else
                call check_if_same_profile(istn,jstn,flgs_h,flgs_o,obsv,niv_stn,nbvar,nltot,ids_stn,same_profile)
              end if

              if( same_profile ) then

                n_checked  = n_checked + 1

                call compare_profiles(istn,jstn,ids_stn,flgs_h,flgs_t,flgs_o,obsv,ombv,niv_stn,nbvar,nltot,cloche,select_critere,select_ij)

                n_choix(select_critere) = n_choix(select_critere) + 1

                ij_bufr = jstn
                if(      btest(flgs_h(istn),23) ) ij_bufr = istn
                ij_tac  = jstn
                if( .not.btest(flgs_h(istn),23) ) ij_tac = istn

                if (mpi_myid == 0) write(11,'(A9,2F10.3,3I10)') ids_stn(istn),lat_stn(istn),lon_stn(istn) &
                     ,niv_stn(ij_bufr+1)-niv_stn(ij_bufr)+1, niv_stn(ij_tac+1)-niv_stn(ij_tac)+1,select_critere

                if( select_ij == istn) ids_stn(jstn) = 'NOT_VALID'
                if( select_ij == jstn) ids_stn(istn) = 'NOT_VALID'

              end if ! same_profile

            end if ! condition 

          end if ! ids_stn(jstn) /= 'NOT_VALID'

        end do ! jstn

      end if ! ids_stn(istn) /= 'NOT_VALID' 

    end do ! istn

    write(*,*)
    write(*,*)
    write(*,'(a30,I10)') 'nb of total duplicates '  ,duplicate_total
    write(*,'(a30,I10)') 'nb TAC vs BURF checked'   ,n_checked
    write(*,'(a30,I10)') 'nb not selected'          ,n_choix(1)
    write(*,'(a30,I10)') 'nb suspicious traj BUFR'  ,n_choix(2)
    write(*,'(a30,I10)') 'nb Energy higher   BUFR'  ,n_choix(3)
    write(*,'(a30,I10)') 'nb variables lower BUFR'  ,n_choix(4)
    write(*,'(a30,I10)') 'nb TAC  selected'         ,n_choix(2) + n_choix(3) + n_choix(4)
    write(*,'(a30,I10)') 'nb BUFR selected'         ,n_choix(5)
    write(*,*)
    write(*,*)

    do iele=1,29
      write(*,'(a15,f4.2,a3,f4.2,I5)') 'nb e ratio ',(iele/10.)-.05,' - ',(iele/10.)+.05,cloche(iele)
    end do
    write(*,*)
    write(*,'(a30,I5)') 'nb of e_tot(2)) very small  ',cloche(30)
    write(*,*)

    if (mpi_myid == 0) close(unit=11)

    ! Check whether there is still duplications

    nb_stnid = 0

    do istn=1,ntot_stn

      if( ids_stn(istn) /= 'NOT_VALID' ) then

        if (nb_stnid < max_nb_stnid ) then
          stnid_not_found=.true.
          if ( nb_stnid == 0) then
            nb_stnid = nb_stnid + 1
            stnid_found(nb_stnid) = ids_stn(istn)
            stnid_numbr(nb_stnid) = istn
          else
            nb_same = 0
            do jstn=1,nb_stnid
              if ( stnid_found(jstn) == ids_stn(istn) ) then
                stnid_not_found=.false.
                nb_same = nb_same + 1
                kstn = stnid_numbr(jstn)
              end if
            end do
            if ( stnid_not_found ) then
              nb_stnid = nb_stnid + 1
              stnid_found(nb_stnid) = ids_stn(istn)
              stnid_numbr(nb_stnid) = istn
            else
              write(*,*) 'Multi profiles found : ',ids_stn(istn),ids_stn(kstn),nb_same,istn,kstn,btest(flgs_h(istn),23),btest(flgs_h(kstn),23)
              write(*,'(a30,2i10,2f10.2,i10)') 'date, lch, lat lon ',date(istn),temps_lch(istn),lat_stn(istn),lon_stn(istn),niv_stn(istn+1)-niv_stn(istn)+1
              write(*,'(a30,2i10,2f10.2,i10)') 'date, lch, lat lon ',date(kstn),temps_lch(kstn),lat_stn(kstn),lon_stn(kstn),niv_stn(kstn+1)-niv_stn(kstn)+1
              nniv = MIN(niv_stn(istn+1)-niv_stn(istn)+1,niv_stn(kstn+1)-niv_stn(kstn)+1)
              do iniv=1,nniv
                ilev = niv_stn(istn)+iniv
                klev = niv_stn(kstn)+iniv
              end do
              write(*,*)

            end if
          end if
        else
          write(*,*) 'nb_stnid >= ',max_nb_stnid
          call qqexit(1)
        end if

      end if

    end do

  end subroutine check_duplicated_stations

  subroutine check_if_same_profile(istn,jstn,flgs_h,flgs_o,obsv,niv_stn,nbvar,nltot,ids_stn,same_profile)
    ! Purpose : Check duplicated stations and select the best TAC/BUFR profiles
    !           
    ! Version 2.0
    !
    ! Author: S. Laroche (ARMA)  June 2013
    !
    ! Revision: 2.0 S. Laroche (ARMA) Decembre 2016
    ! Purpose:  
    !        
    implicit none

    integer,           intent(in)    :: istn,jstn,nbvar,nltot
    integer,           intent(in)    :: flgs_o(:,:)
    real,              intent(in)    :: obsv(:,:)
    integer,           intent(in)    :: niv_stn(:),flgs_h(:)
    logical,           intent(out)   :: same_profile
    character (len=9), intent(inout) :: ids_stn(:)

    logical :: condition,debug_show_thinning_profiles
    integer, parameter :: nlev_man = 16
    integer :: ivar, iniv, ilev, lev_i, lev_j, n_sum
    real    :: sum_val, min_del_p, min_del_p1, min_del_p2

    ! Standard levels
    real, dimension(nlev_man) :: man_levels
    man_levels = (/ 1000.,925.,850.,700.,500.,400.,300.,250.,200.,150.,100.,70.,50.,30.,20.,10./) 

    sum_val = 0.0
    n_sum   = 0

    same_profile = .false.

    do iniv=1,nlev_man

      min_del_p1 = 1000.
      do ilev=niv_stn(istn)+1,niv_stn(istn+1)
        if( ABS(man_levels(iniv) - obsv(5,ilev)) < min_del_p1 ) then
          lev_i = ilev
          min_del_p1 = ABS(man_levels(iniv) - obsv(5,ilev))
        end if
      end do
      min_del_p2 = 1000.
      do ilev=niv_stn(jstn)+1,niv_stn(jstn+1)
        if( ABS(man_levels(iniv) - obsv(5,ilev)) < min_del_p2 ) then
          lev_j = ilev
          min_del_p2 = ABS(man_levels(iniv) - obsv(5,ilev))
        end if
      end do

      if( min_del_p1 < 1.0 .and. min_del_p2 < 1.0 ) then
        do ivar=2,4,2

          if ( obsv(ivar,lev_i) /= -999.0 .and. obsv(ivar,lev_j) /= -999.0 ) then
            sum_val = sum_val + ABS(obsv(ivar,lev_i) - obsv(ivar,lev_j))
            n_sum   = n_sum + 1
          end if

        end do

      end if

    end do !iniv

    if( n_sum > 0 ) then
      if(sum_val/n_sum < 0.5) same_profile = .true.
    end if

  end subroutine check_if_same_profile

  subroutine compare_profiles(istn,jstn,ids_stn,flgs_h,flgs_t,flgs_o,obsv,ombv,niv_stn,nbvar,nltot,cloche,select_critere,select_ij)
    ! Purpose : 
    !           
    ! Version 2.0
    !
    ! Author: S. Laroche (ARMA)  Decembre 2016
    !
    ! Revision: 2.0 S. Laroche (ARMA) Decembre 2016
    !
    implicit none

    integer,           intent(in)    :: istn,jstn,nbvar,nltot
    integer,           intent(in)    :: flgs_t(:,:)
    integer,           intent(in)    :: flgs_o(:,:)
    real,              intent(in)    :: obsv(:,:),ombv(:,:)
    integer,           intent(in)    :: niv_stn(:),flgs_h(:)
    integer,           intent(inout) :: cloche(:)
    integer,           intent(out)   :: select_critere,select_ij
    character (len=9), intent(inout) :: ids_stn(:)

    logical :: condition,tac_bufr,trajectoire_OK
    integer :: iexp,ivar,iniv,ilev,iele,ijtt,ij_bufr,ij_tac,nb_time_f,nb_posi_f,nb_traj
    integer :: nb_values(2)
    real ::  p_bot,v_bot,d_pres,t_pres,d_ombv,e_norm,p_lo,p_hi
    real ::  pourcent_traj_incorrect,factor_tolerance_NE
    real ::  e_var(nbvar,2),e_tot(2),p_bas(nbvar,2),p_haut(nbvar,2),max_ombv(nbvar,2)

    pourcent_traj_incorrect = 10.0
    factor_tolerance_NE     =  1.4
    iele = 1
    max_ombv(:,:)  = 0.0
    select_critere = 1

    tac_bufr       = .true.
    trajectoire_OK = .true.

    if( btest(flgs_h(istn),23)      .and.      btest(flgs_h(jstn),23) ) tac_bufr = .false.
    if( .not.btest(flgs_h(istn),23) .and. .not.btest(flgs_h(jstn),23) ) tac_bufr = .false.

    if( tac_bufr ) then

      ij_bufr = jstn
      if(      btest(flgs_h(istn),23) ) ij_bufr = istn
      ij_tac  = jstn
      if( .not.btest(flgs_h(istn),23) ) ij_tac  = istn
      select_ij = ij_tac

      ! 1. Evalue si la trajectoire native est correct

      if( btest(flgs_h(ij_bufr),14) ) then

        nb_time_f = 0
        nb_posi_f = 0
        nb_traj   = niv_stn(ij_bufr+1) - niv_stn(ij_bufr) + 1

        do iniv=niv_stn(ij_bufr)+1,niv_stn(ij_bufr+1)
          if( btest(flgs_t(1,iniv),4) ) nb_time_f = nb_time_f + 1
          if( btest(flgs_t(2,iniv),4) ) nb_posi_f = nb_posi_f + 1
        end do

        if(100.*nb_time_f/nb_traj > pourcent_traj_incorrect .or. 100.*nb_posi_f/nb_traj > pourcent_traj_incorrect) then
          trajectoire_OK = .false.
          select_ij = ij_tac
          select_critere = 2
          !    write(*,'(2a15,2f10.3,i10)') 'time posi ',ids_stn(ij_bufr),100.*nb_time_f/nb_traj,100.*nb_posi_f/nb_traj,nb_traj
        end if

      end if

      if(trajectoire_OK) then

        ! Cherche les pressions (bas et haut) des extrimites des profils (TAC et BUFR)

        do iexp=1,2

          if(iexp == 1) ijtt = ij_bufr
          if(iexp == 2) ijtt = ij_tac

          do ivar=1,3

            p_bas(ivar,iexp) = 1000.0
            do iniv=niv_stn(ijtt)+1,niv_stn(ijtt+1)
              condition = ombv(ivar,iniv) /= -999.0 .and. &
                   .not.btest(flgs_o(ivar,iniv),18) .and. .not.btest(flgs_o(ivar,iniv),16) .and. &
                   .not.btest(flgs_o(ivar,iniv),9)  .and. .not.btest(flgs_o(ivar,iniv),8)  .and. &
                   .not.btest(flgs_o(ivar,iniv),2)  .and. .not.btest(flgs_o(ivar,iniv),11)
              if( condition ) then
                p_bas(ivar,iexp) = obsv(5,iniv)
                exit
              end if
 
            end do
            p_haut(ivar,iexp) = p_bas(ivar,iexp)
 
            if(iniv < niv_stn(ijtt+1) ) then

              do ilev=iniv+1,niv_stn(ijtt+1)

                condition = ombv(ivar,ilev) /= -999.0 .and. &
                     .not.btest(flgs_o(ivar,ilev),18) .and. .not.btest(flgs_o(ivar,ilev),16) .and. &
                     .not.btest(flgs_o(ivar,ilev),9)  .and. .not.btest(flgs_o(ivar,ilev),8)  .and. &
                     .not.btest(flgs_o(ivar,ilev),2)  .and. .not.btest(flgs_o(ivar,ilev),11)
                if( condition ) then
                  p_haut(ivar,iexp) = obsv(5,ilev)
                end if

              end do

            end if

          end do !ivar

        end do !iexp

        ! Calcul de la norme energie equivalente (NE)

        do iexp=1,2

          if(iexp == 1) ijtt = ij_bufr
          if(iexp == 2) ijtt = ij_tac

          nb_values(iexp) = 0

          do ivar=1,3

            p_lo =  p_bas(ivar,1)
            if( p_lo > p_bas(ivar,2)  ) p_lo = p_bas(ivar,2)
            p_hi =  p_haut(ivar,1)
            if( p_hi < p_haut(ivar,2) ) p_hi = p_haut(ivar,2)

            t_pres           = 0.0
            e_norm           = 0.0
            e_var(ivar,iexp) = 0.0

            do iniv=niv_stn(ijtt)+1,niv_stn(ijtt+1)

              condition = ombv(ivar,iniv) /= -999.0 .and. &
                   .not.btest(flgs_o(ivar,iniv),18) .and. .not.btest(flgs_o(ivar,iniv),16) .and. &
                   .not.btest(flgs_o(ivar,iniv),9)  .and. .not.btest(flgs_o(ivar,iniv),8)  .and. &
                   .not.btest(flgs_o(ivar,iniv),2)  .and. .not.btest(flgs_o(ivar,iniv),11)
              if( condition ) then
                nb_values(iexp) = nb_values(iexp) + 1
                if( obsv(5,iniv) <= p_lo ) then
                  p_bot = obsv(   5,iniv)
                  v_bot = ombv(ivar,iniv)
                  exit
                end if
              end if

            end do

            if(iniv < niv_stn(ijtt+1) ) then

              do ilev=iniv+1,niv_stn(ijtt+1)

                condition = ombv(ivar,ilev) /= -999.0 .and. &
                     .not.btest(flgs_o(ivar,ilev),18) .and. .not.btest(flgs_o(ivar,ilev),16) .and. &
                     .not.btest(flgs_o(ivar,ilev),9)  .and. .not.btest(flgs_o(ivar,ilev),8)  .and. &
                     .not.btest(flgs_o(ivar,ilev),2)  .and. .not.btest(flgs_o(ivar,ilev),11)
                if( condition ) then
                  nb_values(iexp) = nb_values(iexp) + 1
                  if( obsv(5,ilev) >= p_hi ) then
                    d_pres = log(p_bot) - log(obsv(5,ilev))
                    d_ombv = ( v_bot + ombv(ivar,ilev) ) / 2 
                    e_norm = e_norm + d_pres * d_ombv**2
                    t_pres = t_pres + d_pres
                    p_bot  = obsv(   5,ilev)
                    v_bot  = ombv(ivar,ilev)
                    if( ABS(d_ombv) > max_ombv(ivar,iexp) ) max_ombv(ivar,iexp) =  ABS(d_ombv)
                  end if
                end if

              end do

            end if
            if(t_pres > 0.) e_var(ivar,iexp) = e_norm/t_pres

          end do !ivar

          e_tot(iexp) = e_var(1,iexp) + e_var(2,iexp) + (1005./300.)*e_var(3,iexp)

        end do !iexp

        ! 2. Choisi le profil TAC (iexp == 2) si le nombre de niveaux ou il y a des donnees
        !    est plus grand ou egal au nombre dans le profil BUFR (iexp == 1)

        if ( nb_values(2) >= nb_values(1) ) then

          select_ij      = ij_tac
          select_critere = 4

        else

          ! 3. Choisi le profil BUFR si la norme energie equivalente du profil ne depasse pas 'factor_tolerance_NE'

          if(ABS(e_tot(2)) > 1.e-6) iele = NINT(10.*e_tot(1)/e_tot(2)) + 1
          if(iele > 30) iele = 30 

          cloche(iele) = cloche(iele) + 1

          if( e_tot(1) <= factor_tolerance_NE*e_tot(2) ) then

            select_ij      = ij_bufr
            select_critere = 5

          else

            select_ij      = ij_tac
            select_critere = 3

          end if

        end if

      end if ! trajectoire_OK

    else

      write(*,*) 'MEME TYPE DE FICHIER' 

    end if !tac_bufr

  end subroutine compare_profiles

  subroutine thinning_model( flgs_o, obsv, mod_pp, nbvar, nblev, nltot, ntot_stn, nlev_max, niv_stn )
    ! Purpose : 
    !           
    ! Version 2.0
    !
    ! Author: S. Laroche (ARMA)  June 2013
    !
    ! Revision: 2.0 S. Laroche (ARMA) Decembre 2016
    !
    implicit none

    integer,           intent(in)    :: nbvar, nblev, nltot, ntot_stn, nlev_max
    integer,           intent(in)    :: niv_stn(:)
    integer,           intent(inout) :: flgs_o(:,:)
    real,              intent(in)    :: obsv(:,:)
    real,              intent(in)    :: mod_pp(:,:)

    logical :: condition,debug_show_thinning_profiles
    integer :: istn,ilev,iniv,iman,ivar,ivar_dr,ivar_uv,ivar_pp,status
    integer :: ielm,jelm,nbelm,niv_valide,grt_nbelm,srt_nbelm
    real, dimension(16) :: man_levels
    real                :: p_top,p_bot,del_p,del_pmin
    integer, dimension(:), allocatable :: niv_elm, val_elm, elm_elm

    ! Standard levels including 925 hPa
    man_levels(1:16) = (/ 1000.,925.,850.,700.,500.,400.,300.,250.,200.,150.,100.,70.,50.,30.,20.,10./) 

    ivar_dr = 1
    ivar_uv = 2
    ivar_pp = 5
    debug_show_thinning_profiles = .false.

    allocate (  niv_elm(nlev_max), STAT=status )
    allocate (  val_elm(nlev_max), STAT=status )
    allocate (  elm_elm(nlev_max), STAT=status )

    do istn=1,ntot_stn

      ! If one flgs_o is set to 0 but not the other then make the flgs_os consistent
      ! and flgs_o wind direction and module if one of these wind variables is missing

      do iniv=niv_stn(istn)+1,niv_stn(istn+1)

        if( flgs_o(ivar_dr,iniv) == 0  .and. flgs_o(ivar_uv,iniv) /= 0 ) flgs_o(ivar_dr,iniv) = flgs_o(ivar_uv,iniv)
        if( flgs_o(ivar_dr,iniv) /= 0  .and. flgs_o(ivar_uv,iniv) == 0 ) flgs_o(ivar_uv,iniv) = flgs_o(ivar_dr,iniv)

        if( (obsv(ivar_dr,iniv)  < 0. .and. obsv(ivar_uv,iniv) >= 0.) .or. &
             (obsv(ivar_dr,iniv) >= 0. .and. obsv(ivar_uv,iniv)  < 0.) ) then
          flgs_o(ivar_dr,iniv) = ibset(flgs_o(ivar_dr,iniv),11)
          flgs_o(ivar_uv,iniv) = ibset(flgs_o(ivar_uv,iniv),11)
        end if

      end do

      do ilev=1,nblev

        ! Pressures at the bottom and top of the model layer

        if(ilev == 1) then
          p_top = mod_pp(istn,ilev)
          p_bot = exp(0.5*alog(mod_pp(istn,ilev+1)*mod_pp(istn,ilev)))
        else if(ilev == nblev) then
          p_top = exp(0.5*alog(mod_pp(istn,ilev-1)*mod_pp(istn,ilev)))
          p_bot = mod_pp(istn,ilev)
        else
          p_top = exp(0.5*alog(mod_pp(istn,ilev-1)*mod_pp(istn,ilev)))
          p_bot = exp(0.5*alog(mod_pp(istn,ilev+1)*mod_pp(istn,ilev)))
        end if

        ! Select the observation levels between p_top and p_bot

        nbelm = 0
        do iniv=niv_stn(istn)+1,niv_stn(istn+1)

          if(obsv(ivar_pp,iniv) > p_top  .and. obsv(ivar_pp,iniv) <= p_bot) then

            nbelm = nbelm + 1
            niv_elm(nbelm) = iniv
            val_elm(nbelm) = 0

            ! val_elm(nbelm) contains the number of valid observations at each level selected

            do ivar=1,nbvar
              condition = obsv(ivar,iniv) >= 0. .and. &
                   .not.btest(flgs_o(ivar,iniv),18) .and. .not.btest(flgs_o(ivar,iniv),16) .and. &
                   .not.btest(flgs_o(ivar,iniv),9)  .and. .not.btest(flgs_o(ivar,iniv),8)  .and. &
                   .not.btest(flgs_o(ivar,iniv),2)  .and. .not.btest(flgs_o(ivar,iniv),11)
              if ( condition ) val_elm(nbelm) = val_elm(nbelm) + 1
            end do

          end if

        end do ! iniv

        if( nbelm > 0 ) then

          ! Sort the observation levels by greatest numbers of valid observations

          do ielm=1,nbelm
            elm_elm(ielm) = ielm
          end do
          do ielm=1,nbelm
            do jelm=ielm,nbelm
              if( val_elm(jelm) > val_elm(ielm) ) then
                srt_nbelm     = val_elm(ielm)
                val_elm(ielm) = val_elm(jelm)
                val_elm(jelm) = srt_nbelm
                srt_nbelm     = elm_elm(ielm)
                elm_elm(ielm) = elm_elm(jelm)
                elm_elm(jelm) = srt_nbelm 
              end if
            end do
          end do
          grt_nbelm = 1
          if( nbelm > 1 ) then
            do ielm=2,nbelm
              if(val_elm(ielm) == val_elm(1)) grt_nbelm = grt_nbelm + 1
            end do
          end if

          do ivar=1,nbvar

            niv_valide = 0

            ! Rule #1 : get valid the observation at the standard level if one is found in the layer

            do ielm=1,nbelm

              iniv = niv_elm(ielm)

              condition = obsv(ivar,iniv) >= 0. .and. &
                   .not.btest(flgs_o(ivar,iniv),18) .and. .not.btest(flgs_o(ivar,iniv),16) .and. &
                   .not.btest(flgs_o(ivar,iniv),9)  .and. .not.btest(flgs_o(ivar,iniv),8)  .and. &
                   .not.btest(flgs_o(ivar,iniv),2)  .and. .not.btest(flgs_o(ivar,iniv),11)

              do iman=1,16
                if(obsv(ivar_pp,iniv) == man_levels(iman) .and. condition ) niv_valide = iniv
              end do

            end do ! ielm

            if( niv_valide == 0 ) then

              ! Rule #2 : get the closest valid observation to the model level and most complete

              del_pmin = p_bot - p_top
              do ielm=1,grt_nbelm

                iniv = niv_elm(elm_elm(ielm))

                del_p = abs(obsv(ivar_pp,iniv) - mod_pp(istn,ilev))
                condition = obsv(ivar,iniv) >= 0. .and. &
                     .not.btest(flgs_o(ivar,iniv),18) .and. .not.btest(flgs_o(ivar,iniv),16) .and. &
                     .not.btest(flgs_o(ivar,iniv),9)  .and. .not.btest(flgs_o(ivar,iniv),8)  .and. &
                     .not.btest(flgs_o(ivar,iniv),2)  .and. .not.btest(flgs_o(ivar,iniv),11)

                if ( del_p < del_pmin .and. condition) then
                  del_pmin = del_p
                  niv_valide = iniv
                end if

              end do ! ielm

            end if  ! niv_valide == 0

            if( niv_valide == 0 ) then

              ! Rule #3 : get the closest valid observation if not previously selected

              del_pmin = p_bot - p_top
              do ielm=1,nbelm

                iniv = niv_elm(ielm)

                del_p = abs(obsv(ivar_pp,iniv) - mod_pp(istn,ilev))
                condition = obsv(ivar,iniv) >= 0. .and. &
                     .not.btest(flgs_o(ivar,iniv),18) .and. .not.btest(flgs_o(ivar,iniv),16) .and. &
                     .not.btest(flgs_o(ivar,iniv),9)  .and. .not.btest(flgs_o(ivar,iniv),8)  .and. &
                     .not.btest(flgs_o(ivar,iniv),2)  .and. .not.btest(flgs_o(ivar,iniv),11)

                if ( del_p < del_pmin .and. condition) then
                  del_pmin = del_p
                  niv_valide = iniv
                end if

              end do ! ielm

            end if  ! niv_valide == 0

            ! Apply thinning flgs_o to all observations except the one on level niv_valide

            do ielm=1,nbelm

              iniv = niv_elm(ielm)

              if( iniv /= niv_valide .and. obsv(ivar,iniv) >= 0. ) flgs_o(ivar,iniv) = ibset(flgs_o(ivar,iniv),11)

            end do

          end do !ivar

        end if !nbelm > 0

      end do !ilev

      ! Following lines for debugging purposes

      if( debug_show_thinning_profiles ) then

        ilev = nblev
        p_bot = mod_pp(istn,ilev)
        write (*,*) ' '
        write (*,'(a40,I8,a40)') '      ==================== Station no. ',istn,' =============================='
        write (*,*) ' '
        write (*,'(a80,f10.2,i10)') 'lowest model layer-----------------------------------------------------------> ',p_bot,ilev
        do while(obsv(ivar_pp,niv_stn(istn)+1) < p_bot)
          ilev = ilev - 1
          p_bot = exp(0.5*alog(mod_pp(istn,ilev+1)*mod_pp(istn,ilev)))
        end do
        if(ilev == nblev) ilev = nblev - 1

        do iniv=niv_stn(istn)+1,niv_stn(istn+1)

          p_bot = exp(0.5*alog(mod_pp(istn,ilev+1)*mod_pp(istn,ilev)))
          do while(obsv(ivar_pp,iniv) < p_bot)
            write (*,'(a80,f10.2,i10)') 'model layer------------------------------------------------------------------> ',p_bot,ilev
            ilev = ilev - 1
            p_bot = exp(0.5*alog(mod_pp(istn,ilev+1)*mod_pp(istn,ilev)))
          end do
          write (*,'(5f8.2,4i8)') obsv(ivar_pp,iniv),  obsv(1,iniv),  obsv(2,iniv),  obsv(3,iniv),  obsv(4,iniv) &
               ,flgs_o(1,iniv),flgs_o(2,iniv),flgs_o(3,iniv),flgs_o(4,iniv)

        end do

      end if

    end do !istn

    deallocate ( niv_elm, val_elm, elm_elm, STAT=status )

  end subroutine thinning_model

  subroutine thinning_es( flgs_o, obsv, mod_pp, nblev, nltot, ntot_stn, nlev_max, niv_stn )
    ! Purpose : Check duplicated stations and select the best TAC/BUFR profiles
    !           
    ! Version 2.0
    !
    ! Author: S. Laroche (ARMA)  June 2013
    !
    ! Revision: 1.1 E. Lapalme (ARMA) August 2014
    !           Bugfix: replace nblev by nltot_es in second do loop
    !
    ! Revision: 2.0 S. Laroche (ARMA) Decembre 2016
    !
    implicit none

    integer,           intent(in)    :: nblev, nltot, ntot_stn, nlev_max
    integer,           intent(in)    :: niv_stn(:)
    integer,           intent(inout) :: flgs_o(:,:)
    real,              intent(in)    :: obsv(:,:)
    real,              intent(in)    :: mod_pp(:,:)

    integer, parameter  :: nbniv_es = 42
    logical :: condition
    integer :: istn, ilev, iniv, ielm, ivar_es, ivar_pp, nbelm, niv_valide, status
    real, dimension(nbniv_es) :: es_levels
    real :: del_p, del_pmin, p_bot, p_top
    integer, dimension(:), allocatable :: niv_elm

    ! Selected pressure levels for additional thinning to humitidy observations

    es_levels(1:nbniv_es) = (/ 1025.,1000.,975.,950.,925.,900.,875.,850.,825., &
                               800.,775.,750.,725.,700.,675.,650.,625.,600.,575., &
                               550.,525.,500.,475.,450.,425.,400.,375.,350.,325., &
                               300.,275.,250.,225.,200.,175.,150.,100., 70., 50., &
                               30., 20., 10./)

    ivar_es = 4 ; ivar_pp = 5

    allocate (  niv_elm(nlev_max), STAT=status )

    do istn=1,ntot_stn

      ! For each station, retain the nearest observations to selected pressure levels in es_levels

      do ilev=1,nbniv_es

        ! Pressures at the bottom and top of the ES layer

        if(ilev == 1) then
          p_bot = es_levels(ilev)
          p_top = exp(0.5*alog(es_levels(ilev+1)*es_levels(ilev)))
        else if(ilev == nbniv_es) then
          p_bot = exp(0.5*alog(es_levels(ilev-1)*es_levels(ilev)))
          p_top = es_levels(ilev)
        else
          p_bot = exp(0.5*alog(es_levels(ilev-1)*es_levels(ilev)))
          p_top = exp(0.5*alog(es_levels(ilev+1)*es_levels(ilev)))
        end if
   
        ! Select the levels between p_top and p_bot

        nbelm = 0

        do iniv=niv_stn(istn)+1,niv_stn(istn+1)

          if(obsv(ivar_pp,iniv) > p_top  .and. obsv(ivar_pp,iniv) <= p_bot) then
            nbelm  = nbelm + 1
            niv_elm(nbelm) = iniv
          end if

        end do ! iniv

        ! Get the closest valid observation to the selected pressure level

        niv_valide = 0
        del_pmin   = p_bot - p_top

        do ielm=1,nbelm

          iniv = niv_elm(ielm)

          del_p = abs(obsv(ivar_pp,iniv) - es_levels(ilev))
          condition = .not.btest(flgs_o(ivar_es,iniv),18) .and. .not.btest(flgs_o(ivar_es,iniv),16) .and. &
               .not.btest(flgs_o(ivar_es,iniv),9)  .and. .not.btest(flgs_o(ivar_es,iniv),8)  .and. &
               .not.btest(flgs_o(ivar_es,iniv),2)  .and. .not.btest(flgs_o(ivar_es,iniv),11) .and. &
               obsv(ivar_es,iniv) >= 0.          .and. del_p < del_pmin

          if ( condition ) then
            del_pmin = del_p
            niv_valide = iniv
          end if

        end do !ielm

        ! Apply thinning flgs_o to all observations except the one selected with niv_valide

        do ielm=1,nbelm

          iniv = niv_elm(ielm)

          if( iniv /= niv_valide .and. obsv(ivar_es,iniv) >= 0. ) flgs_o(ivar_es,iniv) = ibset(flgs_o(ivar_es,iniv),11)

        end do

      end do !ilev

    end do !istn

    deallocate ( niv_elm, STAT=status )

  end subroutine thinning_es

  subroutine backlisting_ecmwf( flgs_o, obsv, ids_stn, stn_type, mod_pp, nblev, nltot, ntot_stn, niv_stn )
    ! Purpose :
    !           
    ! Version 2.0
    !
    ! Author: S. Laroche (ARMA)  June 2013
    !
    ! Revision: 2.0 S. Laroche (ARMA) February 2017
    !              -Refonte du code pour la premiere mise en operation
    !              -Mise a jour des types sondes RS92 et RS41 (type = 1)
    !              -Retrait de la serie 80 et 90
    !
    implicit none

    integer,           intent(in)    :: nblev, nltot, ntot_stn
    integer,           intent(in)    :: niv_stn(:), stn_type(:)
    integer,           intent(inout) :: flgs_o(:,:)
    real,              intent(in)    :: obsv(:,:)
    real,              intent(in)    :: mod_pp(:,:)
    character (len=9), intent(in)    :: ids_stn(:)

    logical :: condition_tt, condition_es, rejet_es
    integer :: ierr, istn, ilev, iniv, ivar_es, ivar_tt, ivar_pp
    integer :: type, n_type_0_p, n_type_0_t, n_type_1_p, n_type_1_t, n_es_total
    integer :: n_type_0_pMpi, n_type_0_tMpi, n_type_1_pMpi, n_type_1_tMpi, n_es_totalMpi

    n_type_0_p = 0 ;n_type_0_t = 0 ; n_type_1_p = 0 ; n_type_1_t = 0 ; n_es_total = 0

    ivar_tt = 3 ; ivar_es = 4 ; ivar_pp = 5

    do istn=1, ntot_stn

      type   =  0

      if ( stn_type(istn) /= -1 ) then

        if ( stn_type(istn) == 14  .or.  &
             stn_type(istn) == 24  .or.  &
             stn_type(istn) == 41  .or.  &
             stn_type(istn) == 42  .or.  &
             stn_type(istn) == 52  .or.  &
             stn_type(istn) == 79  .or.  &
             stn_type(istn) == 80  .or.  &
             stn_type(istn) == 81  .or.  &
             stn_type(istn) == 83  .or.  &
             stn_type(istn) == 141 .or.  &
             stn_type(istn) == 142 ) then

          type = 1

        end if

      end if

      do iniv=niv_stn(istn)+1, niv_stn(istn+1)

        condition_tt =  &
             .not.btest(flgs_o(ivar_tt,iniv),18) .and. .not.btest(flgs_o(ivar_tt,iniv),16) .and. &
             .not.btest(flgs_o(ivar_tt,iniv),9)  .and. .not.btest(flgs_o(ivar_tt,iniv),8)  .and. &
             .not.btest(flgs_o(ivar_tt,iniv),2)  .and. .not.btest(flgs_o(ivar_tt,iniv),11) .and. &
             obsv(ivar_tt,iniv) >= 0.0
        condition_es =  &
             .not.btest(flgs_o(ivar_es,iniv),18) .and. .not.btest(flgs_o(ivar_es,iniv),16) .and. &
             .not.btest(flgs_o(ivar_es,iniv),9)  .and. .not.btest(flgs_o(ivar_es,iniv),8)  .and. &
             .not.btest(flgs_o(ivar_es,iniv),2)  .and. .not.btest(flgs_o(ivar_es,iniv),11) .and. &
             obsv(ivar_es,iniv) >= 0.0

        if( condition_tt ) then

          rejet_es = .false. 
          if(condition_es) n_es_total = n_es_total + 1

          if ( type == 0 ) then

            if ( .not.rejet_es .and. obsv(ivar_pp,iniv) < 300.0 ) then
              if(condition_es) n_type_0_p = n_type_0_p + 1
              rejet_es = .true. 
            end if
            if ( .not.rejet_es .and. obsv(ivar_tt,iniv) < 233.15 ) then
              if(condition_es) n_type_0_t = n_type_0_t + 1
              rejet_es = .true. 
            end if
    
          else if  ( type == 1 ) then

            if ( .not.rejet_es .and. obsv(ivar_pp,iniv) < 100.0 ) then
              if(condition_es) n_type_1_p = n_type_1_p + 1
              rejet_es = .true. 
            end if
            if ( .not.rejet_es .and. obsv(ivar_tt,iniv) < 213.15 ) then
              if(condition_es)  n_type_1_t = n_type_1_t + 1
              rejet_es = .true. 
            end if

          end if

          if( rejet_es ) flgs_o(ivar_es,iniv) = ibset(flgs_o(ivar_es,iniv),8)

        else

          flgs_o(ivar_es,iniv) = flgs_o(ivar_tt,iniv)

        end if

      end do !iniv

    end do !istn

    
    call rpn_comm_allReduce(n_es_total, n_es_totalMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(n_type_0_p, n_type_0_pMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(n_type_0_t, n_type_0_tMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(n_type_1_p, n_type_1_pMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(n_type_1_t, n_type_1_tMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    write(*,*)
    write(*,*) ' Rejet des donnees ES inspire de ECMWF'
    write(*,*)
    write(*,'(a30,I10)') 'nb total donnees es        ',n_es_totalMpi
    write(*,'(a30,I10)') 'nb rejet type 0 pression   ',n_type_0_pMpi
    write(*,'(a30,I10)') 'nb rejet type 0 temperature',n_type_0_tMpi
    write(*,'(a30,I10)') 'nb rejet type 1 pression   ',n_type_1_pMpi
    write(*,'(a30,I10)') 'nb rejet type 1 temperature',n_type_1_tMpi
    write(*,*)
    if (n_es_totalMpi > 0) then
      write(*,'(a30,I10,f10.2)') 'nb rejet total et %', &
           n_type_0_pMpi+n_type_0_tMpi+n_type_1_pMpi+n_type_1_tMpi, &
           100.0*(n_type_0_pMpi+n_type_0_tMpi+n_type_1_pMpi+n_type_1_tMpi)/n_es_totalMpi
    end if
    write(*,*)

  end subroutine backlisting_ecmwf

  
  !--------------------------------------------------------------------------
  ! thn_gbGpsByDistance
  !--------------------------------------------------------------------------
  subroutine thn_gbGpsByDistance(obsdat, deltemps, deldist, removeUncorrected)
    !
    ! :Purpose: Original method for thinning GB-GPS data by the distance method.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    integer,          intent(in)    :: deltemps
    integer,          intent(in)    :: deldist
    logical,          intent(in)    :: removeUncorrected

    ! Local parameters:
    real(4), parameter :: normZtdScore = 50.0 ! normalization factor for zdscores
    character(len=3), parameter :: winpos='mid' ! Preference to obs close to middle of window

    ! Locals:
    integer :: ierr, numHeader, numHeaderMpi, numHeaderMaxMpi, bodyIndex, headerIndex
    integer :: countObs, countObsOutMpi, countObsInMpi
    integer :: obsDate, obsTime, obsVarno, ztdObsFlag, obsFlag, nsize
    integer :: bgckCount, bgckCountMpi, blackListCount, blackListCountMpi
    integer :: unCorrectCount, unCorrectCountMpi, badTimeCount, badTimeCountMpi
    integer :: numSelected, middleStep
    integer :: obsIndex1, obsIndex2, headerIndex1, headerIndex2
    integer :: headerIndexBeg, headerIndexEnd
    real(4) :: thinDistance, deltaLat, deltaLon, obsLat1, obsLat2
    real(4) :: normFormalErr, missingFormalErr, formalError, finalZtdScore, ztdScore
    real(8) :: obsLonInDegrees, obsLatInDegrees, obsStepIndex_r8
    character(len=12)  :: stnId
    logical :: skipThisObs, thisStnIdNoaa
    integer, allocatable :: quality(:), qualityMpi(:)
    integer, allocatable :: obsLonBurpFile(:), obsLatBurpFile(:)
    integer, allocatable :: obsLonBurpFileMpi(:), obsLatBurpFileMpi(:)
    integer, allocatable :: obsStepIndex(:), obsStepIndexMpi(:)
    integer, allocatable :: headerIndexSelected(:), headerIndexSorted(:)
    logical, allocatable :: valid(:), validMpi(:)

    write(*,*)
    write(*,*) 'thn_gbGpsByDistance: Starting'
    write(*,*)

    numHeader = obs_numHeader(obsdat)
    call rpn_comm_allReduce(numHeader, numHeaderMaxMpi, 1, 'mpi_integer', &
                            'mpi_max', 'grid', ierr)
    numHeaderMpi = numHeaderMaxMpi * mpi_nprocs

    ! Check if any observations to be treated
    countObs = 0
    call obs_set_current_header_list(obsdat,'GP')
    HEADER0: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER0
      countObs = countObs + 1
    end do HEADER0

    call rpn_comm_allReduce(countObs, countObsInMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    if (countObsInMpi == 0) then
      write(*,*) 'thn_gbGpsByDistance: no gb-gps observations present'
      return
    end if

    write(*,*) 'thn_gbGpsByDistance: number of obs initial = ', &
               countObs, countObsInMpi

    thinDistance = real(deldist)
    write(*,*)
    write(*,*) 'Minimun thinning distance ',thinDistance

    middleStep   = nint( ((tim_windowSize/2.0d0) - tim_dstepobs/2.d0) / &
                   tim_dstepobs) + 1

    write(*,*)
    write(*,*) 'Number of time bins                     = ', tim_nstepobs
    write(*,*) 'Minimum number of time bins between obs = ', deltemps
    write(*,*) 'Central time bin                        = ', middleStep
    write(*,*)

    ! Allocations: 
    allocate(obsLatBurpFile(numHeaderMaxMpi))
    allocate(obsLonBurpFile(numHeaderMaxMpi))
    allocate(obsStepIndex(numHeaderMaxMpi))
    allocate(quality(numHeaderMaxMpi))
    allocate(valid(numHeaderMaxMpi))

    allocate(obsLatBurpFileMpi(numHeaderMpi))
    allocate(obsLonBurpFileMpi(numHeaderMpi))
    allocate(obsStepIndexMpi(numHeaderMpi))
    allocate(qualityMpi(numHeaderMpi))

    allocate(headerIndexSorted(numHeaderMpi))
    allocate(headerIndexSelected(numHeaderMpi))
    allocate(validMpi(numHeaderMpi))

    validMpi(:) = .false.

    quality(:)  = 9999
    obsLatBurpFile(:)    = 0
    obsLonBurpFile(:)    = 0
    obsStepIndex(:)    = 0

    qualityMpi(:)        = 9999
    obsLatBurpFileMpi(:) = 0
    obsLonBurpFileMpi(:) = 0
    obsStepIndexMpi(:)   = 0

    badTimeCount = 0
    bgckCount = 0
    blackListCount = 0
    unCorrectCount = 0

    ! First pass through observations
    call obs_set_current_header_list(obsdat,'GP')
    HEADER1: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER1

      ! convert and store stnId as integer array
      stnId = obs_elem_c(obsdat,'STID',headerIndex)

      ! get latitude and longitude
      obsLonInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obs_headElem_r(obsdat, OBS_LON, headerIndex)
      obsLatInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obs_headElem_r(obsdat, OBS_LAT, headerIndex)
      obsLonBurpFile(headerIndex) = nint(100.0*obsLonInDegrees)
      obsLatBurpFile(headerIndex) = 9000 + nint(100.0*obsLatInDegrees)
      if (obsLonBurpFile(headerIndex) >= 18000) then
        obsLonBurpFile(headerIndex) = obsLonBurpFile(headerIndex) - 18000
      else
        obsLonBurpFile(headerIndex) = obsLonBurpFile(headerIndex) + 18000
      end if

      ! get step bin
      obsDate = obs_headElem_i(obsdat, OBS_DAT, headerIndex)
      obsTime = obs_headElem_i(obsdat, OBS_ETM, headerIndex)
      call tim_getStepObsIndex(obsStepIndex_r8, tim_getDatestamp(), &
                               obsDate, obsTime, tim_nstepobs)
      obsStepIndex(headerIndex) = nint(obsStepIndex_r8)

      ! thisStnIdNoaa==TRUE means obs has collocated GPS met (Psfc) observations
      thisStnIdNoaa = ((index(stnId,'FSL_') + index(stnId,'-NOAA') + index(stnId,'-UCAR')) > 0)
      ! normalization factor (Units = mm) for ZTD formal error (formalError) 
      if (thisStnIdNoaa) then
        normFormalErr  = 15.0
        missingFormalErr =  7.0
      else
        normFormalErr  = 5.0
        missingFormalErr = 3.0
      end if

      ! get ztd flag and formal error value, element 15032
      formalError = -1.0
      ztdScore = -1.0
      ztdObsFlag = -1
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY1: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY1
        obsVarno = obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex)
        if (obsVarno == bufr_nezd) then
          ! convert units from m to mm
          formalError = 1000.0*obs_bodyElem_r(obsdat, OBS_OER, bodyIndex)
        end if
        if (obsVarno == bufr_ztdScore) then
          ztdScore = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex)
        end if
        if (obsVarno == bufr_nezd) then
          ztdObsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
        end if
      end do BODY1
      if (formalError == -1.0 .or. formalError == 1000.0*MPC_missingValue_R4) then
        formalError = missingFormalErr
      end if
      if (ztdScore == -1.0) then
        ztdScore = 999.0
      end if
      if (ztdObsFlag == -1) then
        call utl_abort('thn_gbGpsByDistance: ztd not found')
      end if

      ! ZTD quality estimate using monitoring zdscore and the formal error
      finalZtdScore = 80*(ztdScore/normZtdScore) + 20*(formalError/normFormalErr)
    
      ! Give preference to FSL/UCAR ZTD observations (usually include collocated GPS met Psfc)
      if (.not. thisStnIdNoaa) finalZtdScore = finalZtdScore + 5.0
    
      ! Give preference to obs near middle or end of the assimilation window
      if (winpos == 'mid') then
        finalZtdScore = finalZtdScore + 25.0*float(abs(middleStep-obsStepIndex(headerIndex))) / &
                        float(middleStep-1)
      else if (winpos == 'end') then
        finalZtdScore = finalZtdScore + 25.0*float(tim_nstepobs-obsStepIndex(headerIndex)) / &
                        float(tim_nstepobs-1)
      end if

      ! Quality (lower is better), typical values 20->100, if no zdscore then > 1600
      quality(headerIndex) = nint(finalZtdScore)

      ! obs is outside time window
      if(obsStepIndex(headerIndex) == -1.0d0) then
        badTimeCount = badTimeCount + 1
        quality(headerIndex) = 9999
      end if

      ! ZTD O-P failed background/topography checks, ZTD is blacklisted, ZTD is not bias corrected
      if (       btest(ztdObsFlag,16) ) bgckCount = bgckCount + 1
      if (       btest(ztdObsFlag,8) )  blackListCount = blackListCount + 1
      if ( btest(ztdObsFlag,16) .or. btest(ztdObsFlag,18) .or. &
           btest(ztdObsFlag,8) ) then
        quality(headerIndex) = 9999
      end if
      if (removeUncorrected) then
        if ( .not. btest(ztdObsflag,6) ) then
          unCorrectCount = unCorrectCount + 1
          quality(headerIndex) = 9999
        end if
      end if

    end do HEADER1

    ! Gather needed information from all MPI tasks
    nsize = numHeaderMaxMpi
    call rpn_comm_allgather(quality,    nsize, 'mpi_integer',  &
                            qualityMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsLatBurpFile,    nsize, 'mpi_integer',  &
                            obsLatBurpFileMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsLonBurpFile,    nsize, 'mpi_integer',  &
                            obsLonBurpFileMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsStepIndex,    nsize, 'mpi_integer',  &
                            obsStepIndexMpi, nsize, 'mpi_integer', 'grid', ierr)

    do obsIndex1 = 1, numHeaderMpi
      headerIndexSorted(obsIndex1)  = obsIndex1
    end do

    call thn_QsortInt(qualityMpi,headerIndexSorted)

    numSelected       = 0   ! number of obs selected so far
    OBS_LOOP: do obsIndex1 = 1, numHeaderMpi

      if ( qualityMpi(obsIndex1) /= 9999 ) then

        headerIndex1 = headerIndexSorted(obsIndex1)

        ! Check if any of the obs already selected are close in space/time to this obs
        ! If no, then keep (select) this obs        
        if( numSelected >= 1 ) then
          skipThisObs = .false.
          LOOP2: do obsIndex2 = 1, numSelected
            headerIndex2 = headerIndexSelected(obsIndex2)
            if ( abs(obsStepIndexMpi(headerIndex1) -  &
                     obsStepIndexMpi(headerIndex2)) < deltemps  ) then
              deltaLat = abs(obsLatBurpFileMpi(headerIndex1) -  &
                             obsLatBurpFileMpi(headerIndex2))/100.
              deltaLon = abs(obsLonBurpFileMpi(headerIndex1) -  &
                             obsLonBurpFileMpi(headerIndex2))/100.
              if (deltaLon > 180.) deltaLon = 360. - deltaLon
              obsLat1 = ((obsLatBurpFileMpi(headerIndex1) - 9000)/100.)
              obsLat2  = ((obsLatBurpFileMpi(headerIndex2) - 9000)/100.)
              if ( thn_distanceArc(deltaLat,deltaLon,obsLat1,obsLat2) < thinDistance ) then
                skipThisObs = .true.
                exit LOOP2
              end if
            end if
          end do LOOP2
        else
          skipThisObs = .false.
        end if

        if (.not. skipThisObs) then
          numSelected = numSelected + 1
          headerIndexSelected(numSelected) = headerIndex1
        end if

      end if

    end do OBS_LOOP

    do obsIndex1 = 1, numSelected
      obsIndex2 = headerIndexSelected(obsIndex1)
      validMpi(obsIndex2) = .true.
    end do

    ! Update local copy of valid from global mpi version
    headerIndexBeg = 1 + mpi_myid * numHeaderMaxMpi
    headerIndexEnd = headerIndexBeg + numHeaderMaxMpi - 1
    valid(:) = validMpi(headerIndexBeg:headerIndexEnd)

    countObs = count(valid)
    call rpn_comm_allReduce(countObs, countObsOutMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    write(*,*) 'thn_gbGpsByDistance: number of obs after thinning = ', &
               countObs, countObsOutMpi

    ! modify the obs flags and count number of obs kept for each stnId
    call obs_set_current_header_list(obsdat,'GP')
    HEADER3: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER3

      ! do not keep this obs: set bit 11 and jump to the next obs
      if (.not. valid(headerIndex)) then
        call obs_set_current_body_list(obsdat, headerIndex)
        BODY3: do 
          bodyIndex = obs_getBodyIndex(obsdat)
          if (bodyIndex < 0) exit BODY3
          obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
          call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,11))
        end do BODY3
        cycle HEADER3
      end if

    end do HEADER3

    call rpn_comm_allReduce(badTimeCount, badTimeCountMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(unCorrectCount, unCorrectCountMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(blackListCount, blackListCountMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(bgckCount, bgckCountMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)

    write(*,*)
    write(*,'(a50,i10)') 'Number of input obs                  = ', countObsInMpi
    write(*,'(a50,i10)') 'Total number of rejected/thinned obs = ', countObsInMpi - countObsOutMpi
    write(*,'(a50,i10)') 'Number of output obs                 = ', countObsOutMpi
    write(*,*)
    write(*,*)
    write(*,'(a50,i10)') 'Number of blacklisted obs            = ', blackListCountMpi
    write(*,'(a50,i10)') 'Number of obs without bias correction= ', unCorrectCountMpi
    write(*,'(a60,4i10)') 'Number of rejects outside time window, BGCK, thinning ', &
         badTimeCountMpi, bgckCountMpi, &
         countObsInMpi - countObsOutMpi - &
         bgckCountMpi - badTimeCountMpi - blackListCountMpi - unCorrectCountMpi

    write(*,*)
    write(*,*) 'thn_gbGpsByDistance: Finished'
    write(*,*)

  end subroutine thn_gbGpsByDistance

  !--------------------------------------------------------------------------
  ! thn_satWindsByDistance
  !--------------------------------------------------------------------------
  subroutine thn_satWindsByDistance(obsdat, familyType, deltemps, deldist)
    !
    ! :Purpose: Original method for thinning SatWinds data by the distance method.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in)    :: familyType
    integer,          intent(in)    :: deltemps
    integer,          intent(in)    :: deldist

    ! Local parameters:
    integer, parameter :: numStnIdMax = 100
    integer, parameter :: numLayers = 11
    real(4), parameter :: layer(numLayers) = (/ 100000., 92500., 85000., 70000., &
                                                50000., 40000., 30000., 25000., &
                                                20000., 15000., 10000. /)

    ! Locals:
    integer :: ierr, numHeader, numHeaderMaxMpi, bodyIndex, headerIndex, stnIdIndex
    integer :: numStnId, stnIdIndexFound, lenStnId, charIndex
    integer :: obsDate, obsTime, layerIndex, obsVarno, obsFlag, uObsFlag, vObsFlag
    integer :: bgckCount, bgckCountMpi, missingCount, missingCountMpi, nsize
    integer :: countObs, countObsOutMpi, countObsInMpi, numSelected, numHeaderMpi
    integer :: obsIndex1, obsIndex2, headerIndex1, headerIndex2
    integer :: headerIndexBeg, headerIndexEnd, mpiTaskId
    real(4) :: thinDistance, deltaLat, deltaLon, obsLat1, obsLat2
    real(4) :: obsPressure
    real(8) :: obsLonInDegrees, obsLatInDegrees
    real(8) :: obsStepIndex_r8, deltaPress, deltaPressMin
    character(len=12)  :: stnId, stnidList(numStnIdMax)
    logical :: obsAlreadySameStep, skipThisObs
    integer :: numObsStnIdOut(numStnIdMax)
    integer :: numObsStnIdInMpi(numStnIdMax), numObsStnIdOutMpi(numStnIdMax)
    integer, allocatable :: stnIdInt(:,:), stnIdIntMpi(:,:), obsMethod(:), obsMethodMpi(:)
    integer, allocatable :: quality(:), qualityMpi(:)
    integer, allocatable :: obsLonBurpFile(:), obsLatBurpFile(:)
    integer, allocatable :: obsLonBurpFileMpi(:), obsLatBurpFileMpi(:)
    integer, allocatable :: obsStepIndex(:), obsStepIndexMpi(:)
    integer, allocatable :: obsLayerIndex(:), obsLayerIndexMpi(:)
    integer, allocatable :: headerIndexSorted(:), headerIndexSelected(:)
    logical, allocatable :: valid(:), validMpi(:), validMpi2(:)

    write(*,*)
    write(*,*) 'thn_satWindsByDistance: Starting'
    write(*,*)

    numHeader = obs_numHeader(obsdat)
    call rpn_comm_allReduce(numHeader, numHeaderMaxMpi, 1, 'mpi_integer', &
                            'mpi_max', 'grid', ierr)

    ! Check if any observations to be treated
    countObs = 0
    call obs_set_current_header_list(obsdat,trim(familyType))
    HEADER0: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER0
      countObs = countObs + 1
    end do HEADER0

    call rpn_comm_allReduce(countObs, countObsInMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    if (countObsInMpi == 0) then
      write(*,*) 'thn_satWindsByDistance: no satwind observations present'
      return
    end if

    write(*,*) 'thn_satWindsByDistance: number of obs initial = ', &
               countObs, countObsInMpi

    thinDistance = real(deldist)
    write(*,*) 'Minimun thinning distance ',thinDistance

    ! Allocations:
    allocate(valid(numHeaderMaxMpi))
    allocate(quality(numHeaderMaxMpi))
    allocate(obsLatBurpFile(numHeaderMaxMpi))
    allocate(obsLonBurpFile(numHeaderMaxMpi))
    allocate(obsStepIndex(numHeaderMaxMpi))
    allocate(obsLayerIndex(numHeaderMaxMpi))
    allocate(obsMethod(numHeaderMaxMpi))
    lenStnId = len(stnId)
    allocate(stnIdInt(lenStnId,numHeaderMaxMpi))

    ! Initializations:
    valid(:) = .false.
    quality(:) = 0
    obsLatBurpFile(:) = 0
    obsLonBurpFile(:) = 0
    obsLayerIndex(:) = 0
    obsStepIndex(:) = 0
    obsMethod(:) = 0
    stnIdInt(:,:) = 0
    numObsStnIdInMpi(:) = 0
    numObsStnIdOut(:) = 0
    numObsStnIdOutMpi(:) = 0
    bgckCount = 0
    missingCount = 0

    ! First pass through observations
    numStnId = 0
    call obs_set_current_header_list(obsdat,trim(familyType))
    HEADER1: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER1

      ! convert and store stnId as integer array
      stnId = obs_elem_c(obsdat,'STID',headerIndex)
      do charIndex = 1, lenStnId
        stnIdInt(charIndex,headerIndex) = iachar(stnId(charIndex:charIndex))
      end do

      ! get latitude and longitude
      obsLonInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obs_headElem_r(obsdat, OBS_LON, headerIndex)
      obsLatInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obs_headElem_r(obsdat, OBS_LAT, headerIndex)
      obsLonBurpFile(headerIndex) = nint(100.0*obsLonInDegrees)
      obsLatBurpFile(headerIndex) = 9000 + nint(100.0*obsLatInDegrees)

      ! get step bin
      obsDate = obs_headElem_i(obsdat, OBS_DAT, headerIndex)
      obsTime = obs_headElem_i(obsdat, OBS_ETM, headerIndex)
      call tim_getStepObsIndex(obsStepIndex_r8, tim_getDatestamp(), &
                               obsDate, obsTime, tim_nstepobs)
      obsStepIndex(headerIndex) = nint(obsStepIndex_r8)

      ! find layer (assumes 1 level only per headerIndex)
      obsPressure = -1.0
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY1: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY1
        
        if (obsPressure <= 0.0) then
          obsPressure = obs_bodyElem_r(obsdat, OBS_PPP, bodyIndex)
          exit BODY1
        end if
      end do BODY1
      ! modify obsPressure to be consistent with operational pgm
      obsPressure = 100.0*nint(obsPressure/100.0)
      deltaPressMin = abs( log(obsPressure) - log(layer(1)) )
      obsLayerIndex(headerIndex) = 1
      do layerIndex = 2, numLayers
        deltaPress = abs( log(obsPressure) - log(layer(layerIndex)) )
        if ( deltaPress < deltaPressMin ) then
          deltaPressMin = deltaPress
          obsLayerIndex(headerIndex) = layerIndex
        end if
      end do

      ! extract additional information
      obsMethod(headerIndex) = obs_headElem_i(obsdat, OBS_SWMT, headerIndex)

      ! set the observation quality based on QI1
      quality(headerIndex) = obs_headElem_i(obsdat, OBS_SWQ1, headerIndex)

      ! find observation flags (assumes 1 level only per headerIndex)
      uObsFlag = -1
      vObsFlag = -1
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY2: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY2
        obsVarno = obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex)
        if (obsVarno == bufr_neuu) then
          uObsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
        else if (obsVarno == bufr_nevv) then
          vObsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
        end if
      end do BODY2

      ! modify quality based on flags
      if (uObsFlag /= -1 .and. vObsFlag /= -1) then
        if ( btest(uObsFlag,16) .or. btest(vObsFlag,16) ) bgckCount = bgckCount + 1
        if ( btest(uObsFlag,16) .or. btest(vObsFlag,16) .or. &
             btest(uObsFlag,18) .or. btest(vObsFlag,18) ) then
          quality(headerIndex) = 0
        end if
      else
        quality(headerIndex) = 0
        MissingCount = MissingCount + 1
      end if

    end do HEADER1

    ! Gather needed information from all MPI tasks
    allocate(validMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(validMpi2(numHeaderMaxMpi*mpi_nprocs))
    allocate(qualityMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsLatBurpFileMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsLonBurpFileMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsStepIndexMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsLayerIndexMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsMethodMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(stnIdIntMpi(lenStnId,numHeaderMaxMpi*mpi_nprocs))

    nsize = numHeaderMaxMpi
    call rpn_comm_allgather(quality,    nsize, 'mpi_integer',  &
                            qualityMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsLatBurpFile,    nsize, 'mpi_integer',  &
                            obsLatBurpFileMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsLonBurpFile,    nsize, 'mpi_integer',  &
                            obsLonBurpFileMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsStepIndex,    nsize, 'mpi_integer',  &
                            obsStepIndexMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsLayerIndex,    nsize, 'mpi_integer',  &
                            obsLayerIndexMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsMethod,    nsize, 'mpi_integer',  &
                            obsMethodMpi, nsize, 'mpi_integer', 'grid', ierr)
    nsize = lenStnId * numHeaderMaxMpi
    call rpn_comm_allgather(stnIdInt,    nsize, 'mpi_integer',  &
                            stnIdIntMpi, nsize, 'mpi_integer', 'grid', ierr)

    ! build a global list of stnId over all mpi tasks
    numHeaderMpi = numHeaderMaxMpi * mpi_nprocs
    HEADER: do headerIndex = 1, numHeaderMpi
      if (all(stnIdIntMpi(:,headerIndex) == 0)) cycle HEADER

      ! Station ID converted back to character string
      do charIndex = 1, lenStnId
        stnId(charIndex:charIndex) = achar(stnIdIntMpi(charIndex,headerIndex))
      end do

      if (numStnId < numStnIdMax ) then
        stnIdIndexFound = -1
        do stnIdIndex = 1, numStnId
          if ( stnidList(stnIdIndex) == stnid ) stnIdIndexFound = stnIdIndex
        end do
        if ( stnIdIndexFound == -1 ) then
          numStnId = numStnId + 1
          stnidList(numStnId) = stnid
          stnIdIndexFound = numStnId
        end if
        numObsStnIdInMpi(stnIdIndexFound) = numObsStnIdInMpi(stnIdIndexFound) + 1
      else
        call utl_abort('thn_satWindsByDistance: numStnId too large')
      end if
    end do HEADER

    ! Thinning procedure
    allocate(headerIndexSorted(numHeaderMaxMpi*mpi_nprocs))
    do headerIndex = 1, numHeaderMaxMpi*mpi_nprocs
      headerIndexSorted(headerIndex) = headerIndex
    end do
    allocate(headerIndexSelected(numHeaderMaxMpi*mpi_nprocs))
    headerIndexSelected(:) = 0

    call thn_QsortInt(qualityMpi,headerIndexSorted)

    validMpi(:) = .false.
    call tmg_start(144,'bruteThinning')
    STNIDLOOP: do stnIdIndex = 1, numStnId
      write(*,*) 'thn_satWindsByDistance: applying thinning for: ', &
                 trim(stnidList(stnIdIndex))

      LAYERLOOP: do layerIndex = 1, numLayers

        ! do selection of obs for 1 satellite and layer at a time, on separate mpi tasks
        mpiTaskId = mod((stnIdIndex-1)*numLayers + layerIndex - 1, mpi_nprocs)
        if (mpi_myid /= mpiTaskId) cycle LAYERLOOP

        numSelected       = 0
        OBSLOOP1: do obsIndex1 = 1, numHeaderMpi
          headerIndex1 = headerIndexSorted(numHeaderMpi-obsIndex1+1)

          ! only consider obs with current layer being considered
          if (obsLayerIndexMpi(headerIndex1) /= layerIndex) cycle OBSLOOP1

          ! only consider obs with high quality
          if (qualityMpi(numHeaderMpi-obsIndex1+1) <= 10) cycle OBSLOOP1

          ! only consider obs from current satellite
          do charIndex = 1, lenStnId
            stnId(charIndex:charIndex) = achar(stnIdIntMpi(charIndex,headerIndex1))
          end do
          if (stnidList(stnIdIndex) /= stnId) cycle OBSLOOP1

          ! On compte le nombre d'observations qui sont deja
          ! selectionnees avec les memes parametres 'obsStepIndex' et 'obsLayerIndex'
          ! que l'observation consideree ici.
          call tmg_start(145,'countLoop')
          obsAlreadySameStep = .false.
          OBSLOOP2: do obsIndex2 = 1, numSelected
            headerIndex2 = headerIndexSelected(obsIndex2)
            if ( obsStepIndexMpi(headerIndex1) == obsStepIndexMpi(headerIndex2) ) then
              if ( (obsLatBurpFileMpi(headerIndex1) == obsLatBurpFileMpi(headerIndex2)) .and. &
                   (obsLonBurpFileMpi(headerIndex1) == obsLonBurpFileMpi(headerIndex2)) ) then
                ! Si une observation selectionnee porte deja le meme lat, lon, layer et step.
                call tmg_stop(145)
                cycle OBSLOOP1
              end if
              obsAlreadySameStep = .true.
              exit OBSLOOP2
            end if
          end do OBSLOOP2
          call tmg_stop(145)

          if ( obsAlreadySameStep ) then
            call tmg_start(146,'distanceLoop')
            ! Calcule les distances entre la donnee courante et toutes celles choisies 
            ! precedemment.
            skipThisObs = .false.
            OBSLOOP3: do obsIndex2 = 1, numSelected
              headerIndex2 = headerIndexSelected(obsIndex2)
              if ( abs( obsStepIndexMpi(headerIndex1) - &
                        obsStepIndexMpi(headerIndex2) ) < deltemps ) then
                deltaLat = abs( obsLatBurpFileMpi(headerIndex1) - &
                                obsLatBurpFileMpi(headerIndex2) ) / 100.
                deltaLon = abs( obsLonBurpFileMpi(headerIndex1) - &
                                obsLonBurpFileMpi(headerIndex2) ) / 100.
                if(deltaLon > 180.) deltaLon = 360. - deltaLon
                obsLat1 = ((obsLatBurpFileMpi(headerIndex1) - 9000)/100.)
                obsLat2 = ((obsLatBurpFileMpi(headerIndex2) - 9000)/100.)
                if ( thn_distanceArc(deltaLat,deltaLon,obsLat1,obsLat2) < thinDistance ) then
                  skipThisObs = .true.
                  exit OBSLOOP3
                end if
              end if
            end do OBSLOOP3
            call tmg_stop(146)

            if ( .not. skipThisObs ) then

              ! On selectionne la donnee si toutes celles choisies sont au-dela
              ! de thinDistance. Cet evaluation est faite dan la boucle
              ! 'check_list' precedante.
              numSelected = numSelected + 1
              headerIndexSelected(numSelected) = headerIndex1

            end if

          else

            ! On selectionne la donnee s'il y en a aucune choisie dans les intervalles 
            ! layer et step.
            numSelected = numSelected + 1
            headerIndexSelected(numSelected) = headerIndex1

          end if

        end do OBSLOOP1

        do obsIndex1 = 1, numSelected
          validMpi(headerIndexSelected(obsIndex1)) = .true.
        end do
      
      end do LAYERLOOP
    end do STNIDLOOP
    call tmg_stop(144)

    ! communicate values of validMpi computed on each mpi task
    call tmg_start(147,'reduceKeepObs')
    nsize = numHeaderMaxMpi * mpi_nprocs
    call rpn_comm_allReduce(validMpi, validMpi2, nsize, 'mpi_logical', &
                            'mpi_lor','grid',ierr)
    call tmg_stop(147)

    ! Update local copy of valid from global mpi version
    headerIndexBeg = 1 + mpi_myid * numHeaderMaxMpi
    headerIndexEnd = headerIndexBeg + numHeaderMaxMpi - 1
    valid(:) = validMpi2(headerIndexBeg:headerIndexEnd)
    
    countObs = count(valid)
    call rpn_comm_allReduce(countObs, countObsOutMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    write(*,*) 'thn_satWindsByDistance: number of obs after thinning = ', &
               countObs, countObsOutMpi

    ! modify the obs flags and count number of obs kept for each stnId
    call obs_set_current_header_list(obsdat,trim(familyType))
    HEADER3: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER3

      ! do not keep this obs: set bit 11 and jump to the next obs
      if (.not. valid(headerIndex)) then
        call obs_set_current_body_list(obsdat, headerIndex)
        BODY3: do 
          bodyIndex = obs_getBodyIndex(obsdat)
          if (bodyIndex < 0) exit BODY3
          obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
          call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,11))
        end do BODY3
        cycle HEADER3
      end if

      stnId = obs_elem_c(obsdat,'STID',headerIndex)
      stnIdIndexFound = -1
      do stnIdIndex = 1, numStnId
        if (stnidList(stnIdIndex) == stnId) stnIdIndexFound = stnIdIndex
      end do
      if (stnIdIndexFound == -1) call utl_abort('stnid not found in list')
      numObsStnIdOut(stnIdIndexFound) = numObsStnIdOut(stnIdIndexFound) + 1
    end do HEADER3

    call rpn_comm_allReduce(numObsStnIdOut, numObsStnIdOutMpi, &
                            numStnIdMax, 'mpi_integer', 'mpi_sum', 'grid', ierr)
    call rpn_comm_allReduce(bgckCount, bgckCountMpi, 1, &
                            'mpi_integer', 'mpi_sum', 'grid', ierr)

    ! Print counts
    write(*,*)
    write(*,'(a30, i10)') ' Number of obs in  = ', countObsInMpi
    write(*,'(a30, i10)') ' Total number of reject = ', countObsInMpi - countObsOutMpi
    write(*,'(a30, i10)') ' Number of obs out = ', countObsOutMpi
    write(*,*)
    write(*,'(a30,4i10)') ' Number of rejects from bgck, thinned, missing obs ', &
         bgckCountMpi, &
         countObsInMpi - countObsOutMpi - bgckCountMpi - missingCountMpi, &
         missingCountMpi
    write(*,*)
    write(*,'(a30,i10)') 'Number of satellites found = ',numStnId
    write(*,*)

    write(*,'(a30,2a15)') 'Satellite', 'nb AMVs in'
    write(*,*)
    do stnIdIndex = 1, numStnId
      write(*,'(a30,2i15)') stnidList(stnIdIndex), numObsStnIdInMpi(stnIdIndex)
    end do
    write(*,*)
    write(*,'(a30,2i10,f10.4)') 'Total number of obs in : ',sum(numObsStnIdInMpi)

    write(*,*)
    write(*,'(a30,2a15)') 'Satellite', 'nb AMVs out'
    write(*,*)
    do stnIdIndex = 1, numStnId
      write(*,'(a30,2i15)') stnidList(stnIdIndex), numObsStnIdOutMpi(stnIdIndex)
    end do
    write(*,*)
    write(*,'(a30,2i10,f10.4)') 'Total number of obs out : ',sum(numObsStnIdOutMpi)

    ! Deallocations:
    deallocate(valid)
    deallocate(quality)
    deallocate(obsLatBurpFile)
    deallocate(obsLonBurpFile)
    deallocate(obsStepIndex)
    deallocate(obsLayerIndex)
    deallocate(obsMethod)
    deallocate(stnIdInt)
    deallocate(validMpi)
    deallocate(validMpi2)
    deallocate(qualityMpi)
    deallocate(obsLatBurpFileMpi)
    deallocate(obsLonBurpFileMpi)
    deallocate(obsStepIndexMpi)
    deallocate(obsLayerIndexMpi)
    deallocate(obsMethodMpi)
    deallocate(stnIdIntMpi)
    deallocate(headerIndexSorted)
    deallocate(headerIndexSelected)

    write(*,*)
    write(*,*) 'thn_satWindsByDistance: Finished'
    write(*,*)

  end subroutine thn_satWindsByDistance

  !--------------------------------------------------------------------------
  ! thn_QsortInt
  !--------------------------------------------------------------------------
  recursive subroutine thn_QsortInt(A,B)
    implicit none

    integer, intent(inout) :: A(:)
    integer, intent(inout) :: B(:)
    integer :: iq

    if (size(A) > 1) then
      call thn_QsortIntpartition(A,B,iq)
      call thn_QsortInt(A(:iq-1),B(:iq-1))
      call thn_QsortInt(A(iq:),B(iq:))
    end if

  end subroutine thn_QsortInt

  !--------------------------------------------------------------------------
  ! thn_QsortIntpartition
  !--------------------------------------------------------------------------
  subroutine thn_QsortIntpartition(A,B,marker)
    implicit none

    integer, intent(inout) :: A(:)
    integer, intent(inout) :: B(:)
    integer, intent(out) :: marker
    integer :: i, j, tmpi
    integer :: temp
    integer :: x      ! pivot point

    x = A(1)
    i= 0
    j= size(A) + 1

    do
      j = j-1
      do
        if (A(j) <= x) exit
        j = j-1
      end do
      i = i+1
      do
        if (A(i) >= x) exit
        i = i+1
      end do
      if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
        tmpi = B(i)
        B(i) = B(j)
        B(j) = tmpi
      else if (i == j) then
        marker = i+1
        return
      else
        marker = i
        return
      end if
    end do

  end subroutine thn_QsortIntpartition

  !--------------------------------------------------------------------------
  ! thn_QsortReal8
  !--------------------------------------------------------------------------
  recursive subroutine thn_QsortReal8(A,B)
    implicit none

    real(8), intent(inout) :: A(:)
    integer, intent(inout) :: B(:)
    integer :: iq

    if (size(A) > 1) then
      call thn_QsortReal8partition(A,B,iq)
      call thn_QsortReal8(A(:iq-1),B(:iq-1))
      call thn_QsortReal8(A(iq:),B(iq:))
    end if

  end subroutine thn_QsortReal8

  !--------------------------------------------------------------------------
  ! thn_QsortReal8partition
  !--------------------------------------------------------------------------
  subroutine thn_QsortReal8partition(A,B,marker)
    implicit none

    real(8), intent(inout) :: A(:)
    integer, intent(inout) :: B(:)
    integer, intent(out) :: marker
    integer :: i, j, tmpi
    real(8) :: temp
    real(8) :: x      ! pivot point

    x = A(1)
    i= 0
    j= size(A) + 1

    do
      j = j-1
      do
        if (A(j) <= x) exit
        j = j-1
      end do
      i = i+1
      do
        if (A(i) >= x) exit
        i = i+1
      end do
      if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
        tmpi = B(i)
        B(i) = B(j)
        B(j) = tmpi
      else if (i == j) then
        marker = i+1
        return
      else
        marker = i
        return
      end if
    end do

  end subroutine thn_QsortReal8partition

  !--------------------------------------------------------------------------
  ! thn_distanceArc
  !--------------------------------------------------------------------------
  real function thn_distanceArc( deltaLat, deltaLon, lat1, lat2 )
    implicit none

    real, intent(in) :: deltaLat, deltaLon, lat1, lat2

    real, parameter :: PI = 3.141592654
    real, parameter :: RT = 6374.893
    real :: lat1_r, lat2_r, deltaLat_r, deltaLon_r, term_a

    lat1_r = lat1*PI/180.
    lat2_r = lat2*PI/180.
    deltaLat_r = deltaLat*PI/180.
    deltaLon_r = deltaLon*PI/180.

    term_a  = sin(deltaLat_r/2)*sin(deltaLat_r/2) +  &
              cos(lat1_r)*cos(lat2_r)*sin(deltaLon_r/2)*sin(deltaLon_r/2)
    if(term_a < 0.0) term_a = 0.0
    if(term_a > 1.0) term_a = 1.0

    thn_distanceArc = 2*RT*asin(sqrt(term_a))

  end function thn_distanceArc

  !--------------------------------------------------------------------------
  ! thn_aircraftByBoxes
  !--------------------------------------------------------------------------
  subroutine thn_aircraftByBoxes(obsdat, familyType, deltmax)
    !
    ! :Purpose: Original method for thinning aircraft data by lat-lon boxes.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in)    :: familyType
    integer,          intent(in)    :: deltmax

    ! Locals:
    character(len=20) :: trlmFileName
    character(len=2)  :: fileNumber
    type(struct_hco), pointer :: hco_thinning
    type(struct_vco), pointer :: vco_sfc
    type(struct_gsv)          :: stateVectorPsfc
    integer :: numLon, numLat, nsize, headerIndexBeg, headerIndexEnd
    integer :: nulnam, ierr, lonIndex, latIndex, levIndex, stepIndex, codtyp
    integer :: obsLonIndex, obsLatIndex, obsLevIndex, obsStepIndex
    integer :: numHeader, numHeaderMaxMpi, headerIndex, bodyIndex
    integer :: aiTypeCount(4), aiTypeCountMpi(4)
    integer :: obsFlag, obsVarno, obsDate, obsTime, countObs, countObsMpi
    logical :: ttMissing, huMissing, uuMissing, vvMissing, ddMissing, ffMissing
    real(8) :: zpresa, zpresb, obsPressure, delMinutes, obsStepIndex_r8
    real(8) :: obsLonInRad, obsLatInRad, obsLonInDegrees, obsLatInDegrees
    real(8) :: deltaLon, deltaLat, deltaPress, midDistance, score
    real(4), allocatable :: gridLat(:), gridLon(:)
    integer, allocatable :: rejectCount(:), rejectCountMpi(:)
    real(8), allocatable :: gridPressure(:,:,:,:)
    real(8), pointer     :: surfPressure(:,:,:,:)
    logical, allocatable :: valid(:), validMpi(:), isAircraft(:)
    integer, allocatable :: obsLatIndexVec(:), obsLonIndexVec(:)
    integer, allocatable :: obsTimeIndexVec(:), obsLevIndexVec(:)
    integer, allocatable :: obsLatIndexMpi(:), obsLonIndexMpi(:)
    integer, allocatable :: obsTimeIndexMpi(:), obsLevIndexMpi(:)
    logical, allocatable :: obsUVPresent(:), obsTTPresent(:)
    logical, allocatable :: obsUVPresentMpi(:), obsTTPresentMpi(:)
    real(4), allocatable :: obsDistance(:), obsUU(:), obsVV(:), obsTT(:)
    real(4), allocatable :: obsDistanceMpi(:), obsUUMpi(:), obsVVMpi(:), obsTTMpi(:)
    integer, allocatable :: handlesGrid(:,:,:), numObsGrid(:,:,:)
    real(4), allocatable :: minScoreGrid(:,:,:), minDistGrid(:,:,:), maxDistGrid(:,:,:)
    real(4), allocatable :: uuSumGrid(:,:,:), vvSumGrid(:,:,:), ttSumGrid(:,:,:)

    integer, external :: fnom, fclos
    integer, parameter :: maxLev = 500

    ! Namelist variables:
    real(8) :: rprefinc
    real(8) :: rptopinc
    real(8) :: rcoefinc
    real(4) :: vlev(maxLev)
    integer :: numlev
    namelist /namgem/rprefinc, rptopinc, rcoefinc, numlev, vlev

    write(*,*)
    write(*,*) 'thn_aircraftByBoxes: Starting'
    write(*,*)

    numHeader = obs_numHeader(obsdat)
    call rpn_comm_allReduce(numHeader, numHeaderMaxMpi, 1, 'mpi_integer', &
                            'mpi_max','grid',ierr)

    allocate(valid(numHeaderMaxMpi))
    allocate(isAircraft(numHeaderMaxMpi))
    valid(:) = .false.
    aiTypeCount(:) = 0

    call obs_set_current_header_list(obsdat,trim(familyType))
    HEADER0: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER0

      ! check observation type
      codtyp = obs_headElem_i(obsdat, OBS_ITY, headerIndex)
      if ( codtyp == codtyp_get_codtyp('airep') ) then
        aiTypeCount(1) = aiTypeCount(1) + 1
        valid(headerIndex) = .true.
      else if ( codtyp == codtyp_get_codtyp('amdar')  ) then
        aiTypeCount(2) = aiTypeCount(2) + 1
        valid(headerIndex) = .true.
      else if ( codtyp == codtyp_get_codtyp('acars') ) then
        aiTypeCount(3) = aiTypeCount(3) + 1
        valid(headerIndex) = .true.
      else if ( codtyp == codtyp_get_codtyp('ads') ) then
        aiTypeCount(4) = aiTypeCount(4) + 1
        valid(headerIndex) = .true.
      end if
    end do HEADER0
    isAircraft(:) = valid(:)

    ! Return if no aircraft obs to thin
    allocate(validMpi(numHeaderMaxMpi*mpi_nprocs))
    nsize = numHeaderMaxMpi
    call rpn_comm_allgather(valid,    nsize, 'mpi_logical',  &
                            validMpi, nsize, 'mpi_logical', 'grid', ierr)
    if (count(validMpi(:)) == 0) then
      write(*,*) 'thn_aircraftByBoxes: no aircraft observations present'
      return
    end if

    write(*,*) 'thn_aircraftByBoxes: numHeader, numHeaderMaxMpi = ', &
         numHeader, numHeaderMaxMpi

    countObs = count(valid(:))
    call rpn_comm_allReduce(countObs, countObsMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    write(*,*) 'thn_aircraftByBoxes: number of obs initial = ', countObs, countObsMpi

    ! Setup horizontal thinning grid
    nullify(hco_thinning)
    call hco_SetupFromFile(hco_thinning, './analysisgrid_thinning_ai', &
                           'ANALYSIS', 'Analysis')

    ! Default namelist values
    numlev = 80
    vlev(:) = -1
    rprefinc = 0.0d0
    rptopinc = 0.0d0
    rcoefinc = 0.0d0
    ! Read the namelist defining the vertical levels
    if (utl_isNamelistPresent('namgem','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_aircraftByBoxes: Error opening file flnml')
      read(nulnam,nml=namgem,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_aircraftByBoxes: Error reading namgem namelist')
      if (mpi_myid == 0) write(*,nml=namgem)
      ierr = fclos(nulnam)
    else
      call utl_abort('thn_aircraftByBoxes: Namelist block namgem is missing in the namelist.')
    end if

    ! Setup thinning grid parameters
    numLon = hco_thinning%ni
    numLat = hco_thinning%nj
    allocate(gridLat(numLat))
    allocate(gridLon(numLon))
    gridLon(:) = hco_thinning%lon(:) * MPC_DEGREES_PER_RADIAN_R8
    gridLat(:) = hco_thinning%lat(:) * MPC_DEGREES_PER_RADIAN_R8
    write(*,*) 'thinning grid lats = '
    write(*,*) gridLat(:)
    write(*,*) 'thinning grid lons = '
    write(*,*) gridLon(:)
    write(*,*) 'thinning grid vlev = '
    write(*,*) vlev(1:numLev)

    ! Allocate vectors
    allocate(rejectCount(tim_nstepObs))
    allocate(obsLatIndexVec(numHeaderMaxMpi))
    allocate(obsLonIndexVec(numHeaderMaxMpi))
    allocate(obsLevIndexVec(numHeaderMaxMpi))
    allocate(obsTimeIndexVec(numHeaderMaxMpi))
    allocate(obsDistance(numHeaderMaxMpi))
    allocate(obsUU(numHeaderMaxMpi))
    allocate(obsVV(numHeaderMaxMpi))
    allocate(obsTT(numHeaderMaxMpi))
    allocate(obsUVPresent(numHeaderMaxMpi))
    allocate(obsTTPresent(numHeaderMaxMpi))

    ! Allocate mpi global vectors
    allocate(rejectCountMpi(tim_nstepObs))
    allocate(obsLatIndexMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsLonIndexMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsLevIndexMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsTimeIndexMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsDistanceMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsUUMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsVVMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsTTMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsUVPresentMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsTTPresentMpi(numHeaderMaxMpi*mpi_nprocs))

    ! Initialize vectors
    rejectCount(:) = 0
    obsLatIndexVec(:) = 0
    obsLonIndexVec(:) = 0
    obsLevIndexVec(:) = 0
    obsTimeIndexVec(:) = 0
    obsDistance(:) = 0.
    obsUU(:) = 0.
    obsVV(:) = 0.
    obsTT(:) = 0.
    obsUVPresent(:) = .true.
    obsTTPresent(:) = .true.

    ! Read and interpolate the trial surface pressure
    nullify(vco_sfc)
    trlmFileName = './trlm_01'
    call vco_setupFromFile(vco_sfc, trlmFileName)
    call gsv_allocate( stateVectorPsfc, tim_nstepobs, hco_thinning, vco_sfc, &
                       datestamp_opt=tim_getDatestamp(), mpi_local_opt=.false., &
                       varNames_opt=(/'P0'/), hInterpolateDegree_opt='LINEAR' )
    do stepIndex = 1, tim_nstepobs
      write(fileNumber,'(I2.2)') stepIndex
      trlmFileName = './trlm_' // trim(fileNumber)
      call gsv_readFromFile( stateVectorPsfc, trlmFileName, ' ', ' ',  &
                             stepIndex_opt=stepIndex, containsFullField_opt=.true. )
    end do

    ! compute pressure of each model level using p0
    allocate(gridPressure(numLon,numLat,numLev,tim_nstepobs))
    gridPressure(:,:,:,:) = -1.0
    call gsv_getField(stateVectorPsfc,surfPressure)
    do stepIndex = 1, tim_nstepobs
      do levIndex  = 1, numLev
        zpresb = ( (vlev(levIndex) - rptopinc/rprefinc) /  &
                   (1.0D0-rptopinc/rprefinc) )**rcoefinc
        zpresa = rprefinc * (vlev(levIndex)-zpresb)
        do latIndex = 1, numLat
          do lonIndex = 1, numLon
            gridPressure(lonIndex,latIndex,levIndex,stepIndex) =  &
                 zpresa + zpresb*surfPressure(lonIndex,latIndex,1,stepIndex)
          end do
        end do
      end do
    end do
    call gsv_deallocate(stateVectorPsfc)

    call obs_set_current_header_list(obsdat,trim(familyType))
    HEADER1: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER1

      ! find time difference
      obsDate = obs_headElem_i(obsdat, OBS_DAT, headerIndex)
      obsTime = obs_headElem_i(obsdat, OBS_ETM, headerIndex)
      call tim_getStepObsIndex(obsStepIndex_r8, tim_getDatestamp(), &
                               obsDate, obsTime, tim_nstepobs)
      obsStepIndex = nint(obsStepIndex_r8)
      delMinutes = abs(nint(60.0 * tim_dstepobs * abs(real(obsStepIndex) - obsStepIndex_r8)))

      ! check time window
      if ( delMinutes > deltmax ) then
        valid(headerIndex) = .false.
        rejectCount(obsStepIndex) = rejectCount(obsStepIndex) + 1
      end if

      ! Only accept obs below 175hPa
      obsPressure = -1.0d0
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY
        
        if (obsPressure < 0.0d0) then
          obsPressure = obs_bodyElem_r(obsdat, OBS_PPP, bodyIndex)
          if ( obsPressure < 17500.0 .or. obsPressure > 110000.0 ) valid(headerIndex) = .false.
        end if
      end do BODY

      ttMissing = .true.
      huMissing = .true.
      uuMissing = .true.
      vvMissing = .true.
      ddMissing = .true.
      ffMissing = .true.      
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY2: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY2
        
        obsFlag  = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
        obsVarno = obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex)

        ! find number of elements availables
        if (obsVarno == BUFR_NETT) then
          if ( .not. (btest(obsFlag,18) .or. btest(obsFlag,16) .or. &
                      btest(obsFlag,9)  .or. btest(obsFlag,8)  .or. &
                      btest(obsFlag,2)) ) then
            ttMissing = .false.
            obsTT(headerIndex) = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex)
          end if
        else if (obsVarno == BUFR_NEES) then
          if ( .not. (btest(obsFlag,18) .or. btest(obsFlag,16) .or. &
                      btest(obsFlag,9)  .or. btest(obsFlag,8)  .or. &
                      btest(obsFlag,2)) ) huMissing = .false.
        else if (obsVarno == BUFR_NEUU) then
          if ( .not. (btest(obsFlag,18) .or. btest(obsFlag,16) .or. &
                      btest(obsFlag,9)  .or. btest(obsFlag,8)  .or. &
                      btest(obsFlag,2)) ) then
            uuMissing = .false.
            obsUU(headerIndex) = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex)
          end if
        else if (obsVarno == BUFR_NEVV) then
          if ( .not. (btest(obsFlag,18) .or. btest(obsFlag,16) .or. &
                      btest(obsFlag,9)  .or. btest(obsFlag,8)  .or. &
                      btest(obsFlag,2)) ) then
            vvMissing = .false.
            obsVV(headerIndex) = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex)
          end if
        else if (obsVarno == BUFR_NEDD) then
          if ( .not. (btest(obsFlag,18) .or. btest(obsFlag,16) .or. &
                      btest(obsFlag,9)  .or. btest(obsFlag,8)  .or. &
                      btest(obsFlag,2)) ) ddMissing = .false.
        else if (obsVarno == BUFR_NEFF) then
          if ( .not. (btest(obsFlag,18) .or. btest(obsFlag,16) .or. &
                      btest(obsFlag,9)  .or. btest(obsFlag,8)  .or. &
                      btest(obsFlag,2)) ) ffMissing = .false.
        end if

      end do BODY2

      ! wind components are rejected if speed or direction 
      if (ddMissing .or. ffMissing) then
        uuMissing = .true.
        vvMissing = .true.
      end if

      ! eliminate records with nothing to assimilate
      if (ttMissing .and. huMissing .and. uuMissing .and. vvMissing) then
        rejectCount(obsStepIndex) = rejectCount(obsStepIndex) + 1
        valid(headerIndex) = .false.
      end if

      if( uuMissing .or. vvMissing ) obsUVPresent(headerIndex) = .false.
      if( ttMissing )                obsTTPresent(headerIndex) = .false.

      ! Lat and Lon for each observation
      obsLonInRad = obs_headElem_r(obsdat, OBS_LON, headerIndex)
      obsLatInRad = obs_headElem_r(obsdat, OBS_LAT, headerIndex)
      obsLonInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obsLonInRad
      obsLatInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obsLatInRad

      ! latitude index
      deltaLat = abs(gridLat(1) - obsLatInDegrees)
      obsLatIndex = 1
      do latIndex = 2, numLat
        if (abs(gridLat(latIndex) - obsLatInDegrees) < deltaLat) then
          deltaLat = abs(gridLat(latIndex) - obsLatInDegrees)
          obsLatIndex = latIndex
        end if
      end do

      ! longitude index
      deltaLon = abs(gridLon(1) - obsLonInDegrees)
      obsLonIndex = 1
      do lonIndex = 2, numLon
        if ( abs(gridLon(lonIndex) - obsLonInDegrees) < deltaLon ) then
          deltaLon = abs(gridLon(lonIndex) - obsLonInDegrees)
          obsLonIndex = lonIndex
        end if
      end do

      ! layer index
      deltaPress = abs(gridPressure(obsLonIndex,obsLatIndex,1,obsStepIndex) - obsPressure)
      obsLevIndex = 1
      do levIndex = 2, numLev
        if ( abs(gridPressure(obsLonIndex,obsLatIndex,levIndex,obsStepIndex) - obsPressure) < &
             deltaPress ) then
          deltaPress = abs( gridPressure(obsLonIndex,obsLatIndex,levIndex,obsStepIndex) - &
                            obsPressure )
          obsLevIndex = levIndex
        end if
      end do

      obsLatIndexVec(headerIndex) = obsLatIndex
      obsLonIndexVec(headerIndex) = obsLonIndex
      obsLevIndexVec(headerIndex) = obsLevIndex
      obsTimeIndexVec(headerIndex) = obsStepIndex
      obsDistance(headerIndex) = sqrt((deltaLon*3)**2 + (deltaLat*3)**2 + (deltaPress/100.0)**2)

    end do HEADER1

    call rpn_comm_allReduce(aiTypeCount, aiTypeCountMpi, 4, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(rejectCount, rejectCountMpi, tim_nstepObs, 'mpi_integer', &
                            'mpi_sum','grid',ierr)

    write(*,*)
    write(*,'(a50,i10)') ' Total number of obs = ', sum(aiTypeCountMpi(:))
    write(*,*)
    do stepIndex = 1, tim_nstepobs
      write(*,'(a50,2i10)')' Number of rejects for bin = ', stepIndex, rejectCountMpi(stepIndex)
    end do
    write(*,'(a50,i10)')' Total number of rejects = ', sum(rejectCountMpi)
    write(*,*)
    write(*,'(a50,i10)') '====nb AIREP = ', aiTypeCountMpi(1)
    write(*,'(a50,i10)') '====nb AMDAR = ', aiTypeCountMpi(2)
    write(*,'(a50,i10)') '====nb ACARS = ', aiTypeCountMpi(3)
    write(*,'(a50,i10)') '====nb ADS = ',   aiTypeCountMpi(4)
    write(*,*)

    allocate(handlesGrid(numLat,numLon,numLev))
    allocate(minScoreGrid(numLat,numLon,numLev))
    allocate(minDistGrid(numLat,numLon,numLev))
    allocate(maxDistGrid(numLat,numLon,numLev))
    allocate(numObsGrid(numLat,numLon,numLev))
    allocate(uuSumGrid(numLat,numLon,numLev))
    allocate(vvSumGrid(numLat,numLon,numLev))
    allocate(ttSumGrid(numLat,numLon,numLev))

    ! Make all inputs to the following tests mpiglobal
    nsize = numHeaderMaxMpi
    call rpn_comm_allgather(valid,    nsize, 'mpi_logical',  &
                            validMpi, nsize, 'mpi_logical', 'grid', ierr)
    call rpn_comm_allgather(obsLatIndexVec, nsize, 'mpi_integer',  &
                            obsLatIndexMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsLonIndexVec, nsize, 'mpi_integer',  &
                            obsLonIndexMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsLevIndexVec, nsize, 'mpi_integer',  &
                            obsLevIndexMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsTimeIndexVec,nsize, 'mpi_integer',  &
                            obsTimeIndexMpi,nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsDistance,    nsize, 'mpi_real4',  &
                            obsDistanceMpi, nsize, 'mpi_real4', 'grid', ierr)
    call rpn_comm_allgather(obsUU,    nsize, 'mpi_real4',  &
                            obsUUMpi, nsize, 'mpi_real4', 'grid', ierr)
    call rpn_comm_allgather(obsVV,    nsize, 'mpi_real4',  &
                            obsVVMpi, nsize, 'mpi_real4', 'grid', ierr)
    call rpn_comm_allgather(obsTT,    nsize, 'mpi_real4',  &
                            obsTTMpi, nsize, 'mpi_real4', 'grid', ierr)
    call rpn_comm_allgather(obsUVPresent,    nsize, 'mpi_logical',  &
                            obsUVPresentMpi, nsize, 'mpi_logical', 'grid', ierr)
    call rpn_comm_allgather(obsTTPresent,    nsize, 'mpi_logical',  &
                            obsTTPresentMpi, nsize, 'mpi_logical', 'grid', ierr)

    STEP: do stepIndex = 1, tim_nstepobs
      write(*,'(a50,i10)' ) ' Process bin number = ', stepIndex

      handlesGrid(:,:,:) = -1
      minScoreGrid(:,:,:) = 1000000.
      minDistGrid(:,:,:) = 1000000.
      maxDistGrid(:,:,:) = 0.
      numObsGrid(:,:,:) = 0 
      uuSumGrid(:,:,:) = 0. 
      vvSumGrid(:,:,:) = 0. 
      ttSumGrid(:,:,:) = 0. 

      ! Calcul des distances min et max du centre la boite des rapports 
      ! contenus dans les boites
      do headerIndex = 1, numHeaderMaxMpi*mpi_nprocs
        if( .not. validMpi(headerIndex) ) cycle
        if( obsTimeIndexMpi(headerIndex) /= stepIndex ) cycle
        latIndex = obsLatIndexMpi(headerIndex)
        lonIndex = obsLonIndexMpi(headerIndex)
        levIndex = obsLevIndexMpi(headerIndex)
        if ( obsDistanceMpi(headerIndex) < minDistGrid(latIndex,lonIndex,levIndex) ) then
          minDistGrid(latIndex,lonIndex,levIndex) = obsDistanceMpi(headerIndex)
        end if
        if ( obsDistanceMpi(headerIndex) > maxDistGrid(latIndex,lonIndex,levIndex) ) then
          maxDistGrid(latIndex,lonIndex,levIndex) = obsDistanceMpi(headerIndex)
        end if
      end do

      ! Calcul des sommes de u, v et t des observations situees a une distance midDistance
      ! du centre de la boite
      do headerIndex = 1, numHeaderMaxMpi*mpi_nprocs
        if( .not. validMpi(headerIndex) ) cycle
        if( obsTimeIndexMpi(headerIndex) /= stepIndex ) cycle
        latIndex = obsLatIndexMpi(headerIndex)
        lonIndex = obsLonIndexMpi(headerIndex)
        levIndex = obsLevIndexMpi(headerIndex)
        midDistance = ( minDistGrid(latIndex,lonIndex,levIndex) + &
                        maxDistGrid(latIndex,lonIndex,levIndex) )/2.
        if ( (obsDistanceMpi(headerIndex) < midDistance) .and. &
             obsTTPresentMpi(headerIndex) .and. &
             obsUVPresentMpi(headerIndex) ) then
          numObsGrid(latIndex,lonIndex,levIndex) =  &
               numObsGrid(latIndex,lonIndex,levIndex) + 1
          uuSumGrid(latIndex,lonIndex,levIndex) =  &
               uuSumGrid(latIndex,lonIndex,levIndex) + obsUUMpi(headerIndex)
          vvSumGrid(latIndex,lonIndex,levIndex) =  &
               vvSumGrid(latIndex,lonIndex,levIndex) + obsVVMpi(headerIndex)
          ttSumGrid(latIndex,lonIndex,levIndex) =  &
               ttSumGrid(latIndex,lonIndex,levIndex) + obsTTMpi(headerIndex)
        end if
      end do

      ! Calcul la moyenne de u, v et t s'il y a plus de 3 rapports dans la boite
      do latIndex = 1, numLat
        do lonIndex = 1, numLon
          do levIndex = 1, numLev
            if(numObsGrid(latIndex,lonIndex,levIndex) >= 3) then
              uuSumGrid(latIndex,lonIndex,levIndex) =  &
                   uuSumGrid(latIndex,lonIndex,levIndex)/numObsGrid(latIndex,lonIndex,levIndex)
              vvSumGrid(latIndex,lonIndex,levIndex) =  &
                   vvSumGrid(latIndex,lonIndex,levIndex)/numObsGrid(latIndex,lonIndex,levIndex)
              ttSumGrid(latIndex,lonIndex,levIndex) =  &
                   ttSumGrid(latIndex,lonIndex,levIndex)/numObsGrid(latIndex,lonIndex,levIndex)
            end if
          end do
        end do
      end do

      ! S'il y a plus de 3 rapports dans la boite, le rapport dont le score est le plus
      ! petit est retenu. S'il y a 2 rapports ou moins, le rapport le plus pres du centre
      ! de la boite est retenu.
      do headerIndex = 1, numHeaderMaxMpi*mpi_nprocs
        if( .not. validMpi(headerIndex) ) cycle
        if( obsTimeIndexMpi(headerIndex) /= stepIndex ) cycle
        latIndex = obsLatIndexMpi(headerIndex)
        lonIndex = obsLonIndexMpi(headerIndex)
        levIndex = obsLevIndexMpi(headerIndex)

        if(numObsGrid(latIndex,lonIndex,levIndex) >= 3) then

          midDistance = ( minDistGrid(latIndex,lonIndex,levIndex) + &
                          maxDistGrid(latIndex,lonIndex,levIndex) )/2.
          if ((obsDistanceMpi(headerIndex) < midDistance) .and. &
               obsTTPresentMpi(headerIndex) .and. &
               obsUVPresentMpi(headerIndex) ) then
            score = sqrt( (uuSumGrid(latIndex,lonIndex,levIndex) - &
                           obsUUMpi(headerIndex))**2/(1.4**2) +   &
                          (vvSumGrid(latIndex,lonIndex,levIndex) - &
                           obsVVMpi(headerIndex))**2/(1.4**2) ) + &
                    (ttSumGrid(latIndex,lonIndex,levIndex) - &
                     obsTTMpi(headerIndex))**2/(0.9**2)

            if ( handlesGrid(latIndex,lonIndex,levIndex) /= -1 ) then
              if ( score >= minScoreGrid(latIndex,lonIndex,levIndex) ) then
                validMpi(headerIndex) = .false.
              end if
            end if
          
            if ( validMpi(headerIndex) ) then
              if ( handlesGrid(latIndex,lonIndex,levIndex) /= -1 ) then
                validMpi(handlesGrid(latIndex,lonIndex,levIndex)) = .false.
              end if
              minScoreGrid(latIndex,lonIndex,levIndex) = score
              validMpi(headerIndex) = .true.
              handlesGrid(latIndex,lonIndex,levIndex) = headerIndex
            end if

          else
            validMpi(headerIndex) = .false.
          end if

        else ! if(numObsGrid(latIndex,lonIndex,levIndex) < 3)

          if ( handlesGrid(latIndex,lonIndex,levIndex) /= -1 ) then
            if ( obsDistanceMpi(headerIndex) > minDistGrid(latIndex,lonIndex,levIndex) ) then
              validMpi(headerIndex) = .false.
            end if
          end if

          if ( validMpi(headerIndex) ) then
            if ( handlesGrid(latIndex,lonIndex,levIndex) /= -1 ) then
              validMpi(handlesGrid(latIndex,lonIndex,levIndex)) = .false.
            end if
            validMpi(headerIndex) = .true.
            handlesGrid(latIndex,lonIndex,levIndex) = headerIndex
          end if

        end if

      end do ! headerIndex

    end do  STEP

    ! Update local copy of valid from global mpi version
    headerIndexBeg = 1 + mpi_myid * numHeaderMaxMpi
    headerIndexEnd = headerIndexBeg + numHeaderMaxMpi - 1
    valid(:) = validMpi(headerIndexBeg:headerIndexEnd)

    deallocate(handlesGrid)
    deallocate(minScoreGrid)
    deallocate(minDistGrid)
    deallocate(maxDistGrid)
    deallocate(numObsGrid)
    deallocate(uuSumGrid)
    deallocate(vvSumGrid)
    deallocate(ttSumGrid)

    write(*,*)
    write(*,'(a50,i10)') ' Number of obs in  = ', sum(aiTypeCountMpi(:))
    write(*,'(a50,i10)') ' Number of obs out = ', count(validMpi(:))
    write(*,'(a50,i10)') ' Number of obs not out = ', &
         sum(aiTypeCountMpi(:)) - count(validMpi(:))
    write(*,*)

    ! Modify the flags for rejected observations
    do headerIndex = 1, numHeader
      ! skip observation if we're not supposed to consider it
      if (.not. isAirCraft(headerIndex)) cycle
     
      if (.not. valid(headerIndex)) then
        call obs_set_current_body_list(obsdat, headerIndex)
        BODY3: do 
          bodyIndex = obs_getBodyIndex(obsdat)
          if (bodyIndex < 0) exit BODY3
        
          obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
          call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,11))

        end do BODY3
      end if
    end do

    ! Deallocation
    deallocate(valid)
    deallocate(isAircraft)
    deallocate(validMpi)
    deallocate(gridLat)
    deallocate(gridLon)
    deallocate(rejectCount)
    deallocate(obsLatIndexVec)
    deallocate(obsLonIndexVec)
    deallocate(obsLevIndexVec)
    deallocate(obsTimeIndexVec)
    deallocate(obsDistance)
    deallocate(obsUU)
    deallocate(obsVV)
    deallocate(obsTT)
    deallocate(obsUVPresent)
    deallocate(obsTTPresent)
    deallocate(rejectCountMpi)
    deallocate(obsLatIndexMpi)
    deallocate(obsLonIndexMpi)
    deallocate(obsLevIndexMpi)
    deallocate(obsTimeIndexMpi)
    deallocate(obsDistanceMpi)
    deallocate(obsUUMpi)
    deallocate(obsVVMpi)
    deallocate(obsTTMpi)
    deallocate(obsUVPresentMpi)
    deallocate(obsTTPresentMpi)
    deallocate(gridPressure)

    write(*,*)
    write(*,*) 'thn_aircraftByBoxes: Finished'
    write(*,*)

  end subroutine thn_aircraftByBoxes

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

    write(*,*)
    write(*,*) 'thn_keepNthObs: Starting'
    write(*,*)

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

    write(*,*)
    write(*,*) 'thn_keepNthObs: Finished'
    write(*,*)

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
  ! thn_tovsFilt
  !--------------------------------------------------------------------------
  subroutine thn_tovsFilt(obsdat, delta, codtyp, codtyp2_opt)
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    integer, intent(in)             :: delta
    integer, intent(in)             :: codtyp
    integer, optional, intent(in)   :: codtyp2_opt

    ! Locals:
    integer :: numLat, numLon, headerIndex, headerIndexKeep, latIndex, lonIndex, latIndex2
    integer :: gridIndex, numGridLonsTotal, obsTime, obsDate, numHeader, numHeaderMaxMpi, ierr
    integer :: bodyIndex, stepIndex, obsIndex, obsFov
    integer :: loscan, hiscan, obsFlag, numObs, nsize, minLonBurpFileMpi(mpi_nprocs)
    integer :: procIndex, procIndexKeep, minLonBurpFile, countObs, countObsMpi
    integer :: countQc, countKept, countOther, countKeptMpi, countQcMpi, countGridPoints
    real(4) :: obsLatInRad, obsLonInRad, obsLat, obsLon, distance
    real(4) :: obsLatInDegrees, obsLonInDegrees, minDistance, minDistanceMpi(mpi_nprocs)
    real(4) :: rejectRate, gridLat, gridLon
    real(4) :: percentTotal, percentQc, percentOther, percentKept
    real(8), allocatable :: stepObsIndex(:)
    real(4), allocatable :: gridLats(:), gridLatsAll(:), gridLonsAll(:), obsDistance(:)
    logical, allocatable :: valid(:)
    integer, allocatable :: numGridLons(:), numObsGrid(:), obsGridIndex(:)
    integer, allocatable :: obsLonBurpFile(:), obsLatBurpFile(:), numObsAssim(:)
    integer, allocatable :: headerIndexList(:), headerIndexList2(:)
    integer, allocatable :: obsIndexGrid(:), obsIndexLink(:)

    ! Local parameters:
    integer, parameter :: latLength=10000
    integer, parameter :: lonLength=40000
    integer, parameter :: mxscanamsua=30
    integer, parameter :: mxscanamsub=90
    integer, parameter :: mxscanatms =96
    integer, parameter :: mxscanmwhs2=98

    write(*,*)
    write(*,*) 'thn_tovsFilt: Starting, ', trim(codtyp_get_name(codtyp))
    write(*,*)

    numHeader = obs_numHeader(obsdat)
    call rpn_comm_allReduce(numHeader, numHeaderMaxMpi, 1, 'mpi_integer', &
                            'mpi_max','grid',ierr)
    write(*,*) 'thn_tovsFilt: numHeader, numHeaderMaxMpi = ', numHeader, numHeaderMaxMpi

    ! Check if we have any observations to process
    allocate(valid(numHeaderMaxMpi))
    valid(:) = .false.
    do headerIndex = 1, numHeader
      if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) == codtyp) then
        valid(headerIndex) = .true.
      else if (present(codtyp2_opt)) then
        if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) == codtyp2_opt) then
          valid(headerIndex) = .true.
        end if
      end if
    end do
    if (count(valid(:)) == 0) then
      write(*,*) 'thn_tovsFilt: no observations for this instrument'
      deallocate(valid)
      return
    end if

    countObs = count(valid(:))
    call rpn_comm_allReduce(countObs, countObsMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    write(*,*) 'thn_tovsFilt: countObs initial                        = ', &
               countObs, countObsMpi

    ! Remove RARS obs that are also present from a global originating centre
    call thn_removeRarsDuplicates(obsdat, valid)

    countObs = count(valid(:))
    call rpn_comm_allReduce(countObs, countObsMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    write(*,*) 'thn_tovsFilt: countObs after thn_removeRarsDuplicates = ', &
               countObs, countObsMpi

    numLat = 2*latLength/delta
    numLon = lonLength/delta

    ! Allocations
    allocate(gridLats(numLat))
    allocate(gridLatsAll(numLat*numLon))
    allocate(gridLonsAll(numLat*numLon))
    allocate(numObsGrid(numLat*numLon))
    allocate(obsIndexGrid(numLat*numLon))
    allocate(obsIndexLink(numHeader))
    allocate(headerIndexList(numHeader))
    allocate(headerIndexList2(numHeader))
    allocate(numGridLons(numLat))
    allocate(obsGridIndex(numHeader))
    allocate(obsLonBurpFile(numHeader))
    allocate(obsLatBurpFile(numHeader))
    allocate(numObsAssim(numHeader))
    allocate(obsDistance(numHeader))
    allocate(stepObsIndex(numHeader))

    ! Initialize arrays
    gridLats(:) = 0.0
    gridLatsAll(:) = 0.0
    gridLonsAll(:) = 0.0
    obsDistance(:) = 0.0
    numObsGrid(:) = 0
    obsIndexGrid(:) = 0
    obsIndexLink(:) = 0
    headerIndexList(:) = 0
    headerIndexList2(:) = 0
    numGridLons(:) = 0
    obsGridIndex(:) = 0
    obsLonBurpFile(:) = 0
    obsLatBurpFile(:) = 0
    numObsAssim(:) = 0
    stepObsIndex(:) = 0.0d0

    ! Set up the grid used for thinning
    numGridLonsTotal = 0
    do latIndex = 1, numLat
      gridLats(latIndex)    = (latIndex*180./numLat) - 90.
      distance              = lonLength * cos(gridLats(latIndex) * mpc_pi_r4 / 180.)
      numGridLons(latIndex) = nint(distance/delta)
      numGridLons(latIndex) = max(numGridLons(latIndex),1)
      numGridLonsTotal      = numGridLonsTotal + numGridLons(latIndex)
    end do

    gridIndex = 0
    do latIndex = 1, numLat
      do lonIndex = 1, numGridLons(latIndex)
        gridLatsAll(gridIndex+lonIndex) = gridLats(latIndex)
        gridLonsAll(gridIndex+lonIndex) = (lonIndex-1) * 360. / numGridLons(latIndex)
      end do
      gridIndex = gridIndex + numGridLons(latIndex)
    end do

    ! Loop over all observation locations
    do headerIndex = 1, numHeader
      if ( .not. valid(headerIndex) ) cycle

      ! Lat and Lon for each observation
      obsLonInRad = obs_headElem_r(obsdat, OBS_LON, headerIndex)
      obsLatInRad = obs_headElem_r(obsdat, OBS_LAT, headerIndex)

      obsLonInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obsLonInRad
      obsLatInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obsLatInRad
      obsLonBurpFile(headerIndex) = nint(100.0*(obsLonInDegrees - 180.0))
      if(obsLonBurpFile(headerIndex) < 0) then
        obsLonBurpFile(headerIndex) = obsLonBurpFile(headerIndex) + 36000
      end if
      obsLatBurpFile(headerIndex) = 9000+nint(100.0*obsLatInDegrees)

      obsDate = obs_headElem_i(obsdat, OBS_DAT, headerIndex)
      obsTime = obs_headElem_i(obsdat, OBS_ETM, headerIndex)
      call tim_getStepObsIndex(stepObsIndex(headerIndex), tim_getDatestamp(), &
                               obsDate, obsTime, tim_nstepobs)

      ! Associate each observation to a grid point
      obsLat = (obsLatBurpFile(headerIndex) - 9000.) / 100.
      obsLon = obsLonBurpFile(headerIndex) / 100.
      do latIndex = 1, numLat-1
        if (obsLat <  gridLats(1))      obsLat = gridLats(1)
        if (obsLat >= gridLats(numLat)) obsLat = gridLats(numLat) - 0.5
        if (obsLat >= gridLats(latIndex) .and. obsLat < gridLats(latIndex+1)) then
          gridIndex = 1
          do latIndex2 = 1, latIndex-1
            gridIndex = gridIndex + numGridLons(latIndex2)
          end do
          gridIndex = gridIndex + ifix(obsLon/(360./numGridLons(latIndex)))
          exit
        end if
      end do
      obsGridIndex(headerIndex) = gridIndex
      numObsGrid(gridIndex) = numObsGrid(gridIndex) + 1

    end do ! headerIndex

    if      ( codtyp == codtyp_get_codtyp('amsua') ) then
      loscan   = 1
      hiscan   = mxscanamsua
    else if ( codtyp == codtyp_get_codtyp('amsub') ) then
      loscan   = 1
      hiscan   = mxscanamsub
    else if ( codtyp == codtyp_get_codtyp('atms') ) then
      loscan   = 2
      hiscan   = mxscanatms - 1
    else if ( codtyp == codtyp_get_codtyp('mwhs2') ) then
      loscan   = 1
      hiscan   = mxscanmwhs2
    end if

    countQc = 0
    numObsAssim(:) = 0
    do headerIndex = 1, numHeader

      if ( .not. valid(headerIndex) ) cycle

      ! Look at the obs flags
      rejectRate = 0.

      call obs_set_current_body_list(obsdat, headerIndex)
      BODY: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY
        
        ! If not a blacklisted channel (note that bit 11 is set in 
        ! satqc_amsu*.f for blacklisted channels)
        obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
        if ( .not. btest(obsFlag,11) ) then
          numObsAssim(headerIndex) = numObsAssim(headerIndex) + 1
          if ( btest(obsFlag,9) ) then
            rejectRate = rejectRate + 1.0
          end if
        end if
      end do BODY

      ! fixer le % de rejets a 100% si aucun canal n'est assimilable         
      if ( rejectRate == 0. .and. numObsAssim(headerIndex) == 0 ) then
        rejectRate = 1.
      else
        rejectRate = rejectRate / max(numObsAssim(headerIndex),1)  
      end if

      obsFov = obs_headElem_i(obsdat, OBS_FOV, headerIndex)
      if ( rejectRate >= 0.80 ) then
        countQc = countQc + 1
        valid(headerIndex) = .false.
      else if (obsFov < loscan .or.  &
               obsFov > hiscan ) then
        countQc = countQc + 1
        valid(headerIndex) = .false.
      end if

    end do

    countObs = count(valid(:))
    call rpn_comm_allReduce(countObs, countObsMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    write(*,*) 'thn_tovsFilt: countObs after QC                       = ', &
               countObs, countObsMpi

    ! Calculate distance of obs. from center of its grid box center
    do headerIndex = 1, numHeader
      if ( .not. valid(headerIndex) ) cycle

      gridIndex = obsGridIndex(headerIndex)
      if (numObsGrid(gridIndex) /= 0) then
        latIndex = (gridLatsAll(gridIndex)+90.)/(180./numLat)
        gridLat = gridLatsAll(gridIndex) + 0.5*(180./numLat)
        gridLon = gridLonsAll(gridIndex) + 0.5*360./numGridLons(latIndex)
        obsLat = (obsLatBurpFile(headerIndex) - 9000.) / 100.
        obsLon = obsLonBurpFile(headerIndex) / 100.
        obsDistance(headerIndex) = thn_separation(obsLon,obsLat,gridLon,gridLat) * &
             float(latLength) / 90.
      end if
    end do

    ! Create a linked list of observations (link to grid point)
    obsIndexGrid(:) = 0
    obsIndexLink(:) = 0
    countObs = 0
    do headerIndex = 1, numHeader
      if ( .not. valid(headerIndex) ) cycle

      gridIndex = obsGridIndex(headerIndex)
      if (numObsGrid(gridIndex) /= 0) then
        countObs = countObs + 1
        headerIndexList(countObs) = headerIndex
        obsIndexLink(countObs) = obsIndexGrid(gridIndex)
        obsIndexGrid(gridIndex) = countObs
      end if
    end do

    ! Loop over stepObs
    do stepIndex = 1, tim_nstepobs

      ! Loop over all grid points
      countGridPoints = 0
      do gridIndex = 1, numGridLonsTotal
        if (numObsGrid(gridIndex) /= 0) then
          countGridPoints = countGridPoints + 1
        end if
        numObs = 0
        obsIndex = obsIndexGrid(gridIndex)
        do
          if (obsIndex == 0) exit
          headerIndex = headerIndexList(obsIndex)
          if ( obsGridIndex(headerIndex) == gridIndex  .and. &
               valid(headerIndex)               .and. &
               nint(stepObsIndex(headerIndex)) == stepIndex ) then
            numObs = numObs + 1
            headerIndexList2(numObs) = headerIndex
          end if

          obsIndex = obsIndexLink(obsIndex)
        end do

        minDistance = 1000000.             

        ! Choose the obs closest to the grid point
        do obsIndex = 1, numObs
          if (obsDistance(headerIndexList2(obsIndex)) < minDistance) then
            minDistance = obsDistance(headerIndexList2(obsIndex))
            minLonBurpFile = obsLonBurpFile(headerIndexList2(obsIndex))
            headerIndexKeep = headerIndexList2(obsIndex)
          end if
        end do

        ! Check for multiple obs with same distance to grid point
        if (numObs > 0) then
          if ( count(obsDistance(headerIndexList2(1:numObs)) == minDistance) > 1 ) then
            ! resolve ambiguity by choosing obs with min value of lon
            minLonBurpFile = 10000000
            do obsIndex = 1, numObs
              if (obsDistance(headerIndexList2(obsIndex)) == minDistance) then
                if (obsLonBurpFile(headerIndexList2(obsIndex)) < minLonBurpFile) then
                  minLonBurpFile = obsLonBurpFile(headerIndexList2(obsIndex))
                  headerIndexKeep = headerIndexList2(obsIndex)
                end if
              end if
            end do
          end if
        end if

        do obsIndex = 1, numObs
          valid(headerIndexList2(obsIndex)) = .false.
        end do
        if (numObs > 0 .and. minDistance <= 75. ) then
          valid(headerIndexKeep) = .true.
        end if

        ! Communicate the distance of chosen observation among all mpi tasks
        call rpn_comm_allgather(minDistance,    1, 'mpi_real4',  &
                                minDistanceMpi, 1, 'mpi_real4', 'grid', ierr)

        ! Choose the closest to the center of the box among all mpi tasks
        minDistance = 1000000.
        do procIndex = 1, mpi_nprocs
          if (minDistanceMpi(procIndex) < minDistance) then
            minDistance = minDistanceMpi(procIndex)
            procIndexKeep = procIndex
          end if
        end do

        ! Adjust flags to only keep 1 observation among all mpi tasks
        if (minDistance < 1000000.) then
          if ( count(minDistanceMpi(:) == minDistance) > 1 ) then
            ! resolve ambiguity by choosing obs with min value of lon
            call rpn_comm_allgather(minLonBurpFile,    1, 'mpi_integer',  &
                                    minLonBurpFileMpi, 1, 'mpi_integer', 'grid', ierr)
            minLonBurpFile = 10000000
            do procIndex = 1, mpi_nprocs
              if (minDistanceMpi(procIndex) == minDistance) then
                if (minLonBurpFileMpi(procIndex) < minLonBurpFile) then
                  minLonBurpFile = minLonBurpFileMpi(procIndex)
                  procIndexKeep = procIndex
                end if
              end if
            end do            
          end if
          if (numObs > 0) then
            if (mpi_myid /= (procIndexKeep-1)) then
              valid(headerIndexKeep) = .false.
            end if
          end if
        end if

      end do ! gridIndex
    end do ! stepIndex

    countObs = count(valid(:))
    call rpn_comm_allReduce(countObs, countObsMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    write(*,*) 'thn_tovsFilt: obsCount after thinning                 = ', &
               countObs, countObsMpi

    ! modify the observation flags in obsSpaceData
    countObs = 0
    do headerIndex = 1, numHeader
      ! skip observation if we're not supposed to consider it
      if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) /= codtyp) then
        if (present(codtyp2_opt)) then
          if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) /= codtyp2_opt) then
            cycle
          end if
        else
          cycle
        end if
      end if
     
      countObs = countObs + 1

      if (.not. valid(headerIndex)) then
        call obs_set_current_body_list(obsdat, headerIndex)
        BODY2: do 
          bodyIndex = obs_getBodyIndex(obsdat)
          if (bodyIndex < 0) exit BODY2
        
          obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
          call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,11))

        end do BODY2
      end if
    end do

    ! print a summary to the listing
    countKept = count(valid)
    call rpn_comm_allReduce(countKept, countKeptMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(countQc,   countQcMpi,   1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(countObs, countObsMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)

    countOther = countObsMpi - countKeptMpi - countQcMpi

    percentTotal = 100.0
    percentQc    = (float(countQcMpi)   / float(countObsMpi)) * 100.0
    percentOther = (float(countOther)   / float(countObsMpi)) * 100.0
    percentKept  = (float(countKeptMpi) / float(countObsMpi)) * 100.0
         
    write(*,100)
100 format(/,' SOMMAIRE DES RESULTATS',/)
    write(*,200) countObsMpi, percentTotal, countQcMpi, percentQc, &
                 countOther, percentOther, countKeptMpi, percentKept, &
                 delta, numGridLonsTotal, countGridPoints
200 format(' NB.STNS TOTAL AU DEBUT EN ENTREE        =    ',I7, &
           '  =  ',F6.1,'%',/, &
           ' NB.STNS REJETEES AU CONTROLE QUALITATIF =  - ',I7, &
           '  =  ',F6.1,'%',/, &
           ' NB.STNS REJETEES POUR LES AUTRES RAISONS=  - ',I7, &
           '  =  ',F6.1,'%',/, &
           '                                              ', &
           '-------- = -------',/, &
           ' NB.STNS GARDEES A LA FIN                =    ',I7, &
           '  =  ',F6.1,'%',/////, &
           ' NB.PTS GRILLE AVEC LA RESOLUTION ',I5,' KM   =    ', &
           I7,/, &
           ' NB.PTS GRILLE TRAITES                       =    ',I7)

    ! end of sommair

    deallocate(valid)
    deallocate(gridLats)
    deallocate(gridLatsAll)
    deallocate(gridLonsAll)
    deallocate(numObsGrid)
    deallocate(obsIndexGrid)
    deallocate(obsIndexLink)
    deallocate(headerIndexList)
    deallocate(headerIndexList2)
    deallocate(numGridLons)
    deallocate(obsGridIndex)
    deallocate(obsLonBurpFile)
    deallocate(obsLatBurpFile)
    deallocate(numObsAssim)
    deallocate(obsDistance)
    deallocate(stepObsIndex)

    write(*,*)
    write(*,*) 'thn_tovsFilt: Finished'
    write(*,*)

  end subroutine thn_tovsFilt

  !--------------------------------------------------------------------------
  ! thn_removeRarsDuplicates
  !--------------------------------------------------------------------------
  subroutine thn_removeRarsDuplicates(obsdat, valid)
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    logical,          intent(inout) :: valid(:)

    ! Locals:
    integer :: nsize, ierr, lenStnId, headerIndex, headerIndex1, headerIndex2
    integer :: numHeader, numHeaderMaxMpi, charIndex, headerIndexBeg, headerIndexEnd
    integer :: obsDate, obsTime
    real(4) :: obsLatInRad, obsLonInRad
    real(8) :: dlhours
    logical :: global1, global2
    integer, allocatable :: centreOrig(:), centreOrigMpi(:)
    integer, allocatable :: obsFov(:), obsFovMpi(:)
    integer, allocatable :: obsDateStamp(:), obsDateStampMpi(:)
    integer, allocatable :: stnIdInt(:,:), stnIdIntMpi(:,:)
    logical, allocatable :: validMpi(:)
    character(len=12)    :: stnId
    
    ! Locals related to kdtree2:
    type(kdtree2), pointer            :: tree
    integer, parameter                :: maxNumSearch = 100
    integer                           :: numFoundSearch, resultIndex
    type(kdtree2_result)              :: searchResults(maxNumSearch)
    real(kdkind)                      :: maxRadius = 100.d6
    real(kdkind)                      :: refPosition(3)
    real(kdkind), allocatable         :: obsPosition3d(:,:)
    real(kdkind), allocatable         :: obsPosition3dMpi(:,:)

    ! Local parameters:
    integer, parameter :: centreOrigGlobal(3)=(/53, 74, 160/)

    ! Externals:
    integer, external    :: newdate

    numHeader = obs_numHeader(obsdat)
    call rpn_comm_allReduce(numHeader, numHeaderMaxMpi, 1, 'mpi_integer', &
                            'mpi_max','grid',ierr)

    ! Allocations
    allocate(obsPosition3d(3,numHeaderMaxMpi))
    allocate(obsPosition3dMpi(3,numHeaderMaxMpi*mpi_nprocs))
    allocate(centreOrig(numHeaderMaxMpi))
    allocate(centreOrigMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsFov(numHeaderMaxMpi))
    allocate(obsFovMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsDateStamp(numHeaderMaxMpi))
    allocate(obsDateStampMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(validMpi(numHeaderMaxMpi*mpi_nprocs))
    lenStnId = len(stnId)
    allocate(stnIdInt(lenStnId,numHeaderMaxMpi))
    allocate(stnIdIntMpi(lenStnId,numHeaderMaxMpi*mpi_nprocs))

    ! Some initializations
    centreOrig(:) = 0
    obsPosition3d(:,:) = 0.0

    ! Loop over all observation locations
    do headerIndex = 1, numHeader
      if ( .not. valid(headerIndex) ) cycle

      ! Originating centre of data
      centreOrig(headerIndex) = obs_headElem_i(obsdat, OBS_ORI, headerIndex)

      ! Station ID converted to integer array
      stnId = obs_elem_c(obsdat,'STID',headerIndex)
      do charIndex = 1, lenStnId
        stnIdInt(charIndex,headerIndex) = iachar(stnId(charIndex:charIndex))
      end do

      ! Date stamp for each observation
      obsDate = obs_headElem_i(obsdat, OBS_DAT, headerIndex)
      obsTime = obs_headElem_i(obsdat, OBS_ETM, headerIndex)
      ierr = newdate(obsDateStamp(headerIndex), obsDate, obsTime*10000+2900, 3)

      ! Field of View for each observation
      obsFov(headerIndex) = obs_headElem_i(obsdat, OBS_FOV, headerIndex)

      ! Lat and Lon for each observation
      obsLonInRad = obs_headElem_r(obsdat, OBS_LON, headerIndex)
      obsLatInRad = obs_headElem_r(obsdat, OBS_LAT, headerIndex)

      ! 3D location array for kdtree
      obsPosition3d(1,headerIndex) = RA * sin(obsLonInRad) * cos(obsLatInRad)
      obsPosition3d(2,headerIndex) = RA * cos(obsLonInRad) * cos(obsLatInRad)
      obsPosition3d(3,headerIndex) = RA *                    sin(obsLatInRad)
    end do

    nsize = 3 * numHeaderMaxMpi
    call rpn_comm_allgather(obsPosition3d,    nsize, 'mpi_real8',  &
                            obsPosition3dMpi, nsize, 'mpi_real8', 'grid', ierr)
    nsize = numHeaderMaxMpi
    call rpn_comm_allgather(valid,    nsize, 'mpi_logical',  &
                            validMpi, nsize, 'mpi_logical', 'grid', ierr)
    call rpn_comm_allgather(centreOrig,    nsize, 'mpi_integer',  &
                            centreOrigMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsFov,    nsize, 'mpi_integer',  &
                            obsFovMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsDateStamp,    nsize, 'mpi_integer',  &
                            obsDateStampMpi, nsize, 'mpi_integer', 'grid', ierr)
    nsize = lenStnId * numHeaderMaxMpi
    call rpn_comm_allgather(stnIdInt,    nsize, 'mpi_integer',  &
                            stnIdIntMpi, nsize, 'mpi_integer', 'grid', ierr)
    nullify(tree)

    tree => kdtree2_create(obsPosition3dMpi, sort=.true., rearrange=.true.)
    HEADER1: do headerIndex1 = 1, mpi_nprocs*numHeaderMaxMpi
        
      if ( .not. validMpi(headerIndex1) ) cycle HEADER1

      ! Find all obs within 10km
      refPosition(:) = obsPosition3dMpi(:,headerIndex1)
      call kdtree2_r_nearest(tp=tree, qv=refPosition, r2=maxRadius, nfound=numFoundSearch, &
                             nalloc=maxNumSearch, results=searchResults)
      if (numFoundSearch >= maxNumSearch) then
        call utl_abort('thn_tovsFilt: the parameter maxNumSearch must be increased')
      end if
      if (numFoundSearch == 0) then
        call utl_abort('thn_tovsFilt: no match found. This should not happen!!!')
      end if

      ! Loop over all of these nearby locations
      HEADER2: do resultIndex = 1, numFoundSearch
        headerIndex2 = searchResults(resultIndex)%idx

        if ( .not. validMpi(headerIndex2) ) cycle HEADER2

        ! Certaines stations locales nous envoient 
        ! le mauvais numero d'orbite. On ne peut donc
        ! pas s'y fier.
        ! Il faut comparer le temps de la reception des
        ! donnees
        if ( centreOrigMpi(headerIndex1) /= centreOrigMpi(headerIndex2) ) then
          if ( obsFovMpi(headerIndex1) == obsFovMpi(headerIndex2) ) then
            if ( all(stnIdIntMpi(:,headerIndex1) ==  stnIdIntMpi(:,headerIndex2)) ) then
            
              ! Difference (in hours) between obs time
              call difdatr(obsDateStampMpi(headerIndex1),obsDateStampMpi(headerIndex2),dlhours)

              ! Si la difference est moins de 6 minutes,
              ! on peut avoir affaire a un rars

              if ( abs(dlhours) <= 0.1 ) then

                ! si l'element_i est global, on doit le garder et rejeter l'element_j
                global1 = any(centreOrigGlobal(:) == centreOrigMpi(headerIndex1))
                if (global1) then 
                  validMpi(headerIndex2) = .false.
                else
                  ! toutefois, ca ne signifie pas que l'element_j est un rars
                  ! VERIFIER SI LA STATION 2 EST RARS
                  global2 = any(centreOrigGlobal(:) == centreOrigMpi(headerIndex2))

                  ! Si l'element_j est global, rejeter l'element_i
                  ! Si les 2 elements sont rars, garder le 1er
                  if (global2) then 
                    validMpi(headerIndex1) = .false.
                    cycle HEADER1
                  else
                    validMpi(headerIndex2) = .false.
                  end if
                end if

              end if ! abs(dlhours) <= 0.1
              
            end if ! STID1 == STID2
          end if ! FOV1 == FOV2
        end if ! centreOrig1 /= centreOrig2

      end do HEADER2
    end do HEADER1
    call kdtree2_destroy(tree)

    ! update local copy of 'valid' array
    headerIndexBeg = 1 + mpi_myid * numHeaderMaxMpi
    headerIndexEnd = headerIndexBeg + numHeaderMaxMpi - 1
    valid(:) = validMpi(headerIndexBeg:headerIndexEnd)

    deallocate(obsPosition3d)
    deallocate(obsPosition3dMpi)
    deallocate(centreOrig)
    deallocate(centreOrigMpi)
    deallocate(obsFov)
    deallocate(obsFovMpi)
    deallocate(obsDateStamp)
    deallocate(obsDateStampMpi)
    deallocate(validMpi)
    deallocate(stnIdInt)
    deallocate(stnIdIntMpi)

  end subroutine thn_removeRarsDuplicates

  !--------------------------------------------------------------------------
  ! thn_scatByLatLonBoxes
  !--------------------------------------------------------------------------
  subroutine thn_scatByLatLonBoxes(obsdat, deltax, deltmax)
    !
    ! :Purpose: Only keep the observation closest to the center of each
    !           lat-lon (and time) box for SCAT observations.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    integer, intent(in)             :: deltax
    integer, intent(in)             :: deltmax

    ! Locals parameters:
    integer, parameter :: latLength = 10000 ! Earth dimension parameters
    integer, parameter :: lonLength = 40000 ! Earth dimension parameters
    integer, parameter :: numStnIdMax = 100

    ! Locals:
    integer :: bodyIndex, charIndex, nsize, lenStnId
    integer :: timeRejectCount, flagRejectCount, timeRejectCountMpi, flagRejectCountMpi
    integer :: uObsFlag, vObsFlag, obsVarno, stnIdIndex, numStnId, stnIdIndexFound
    integer :: numLat, numLon, latIndex, lonIndex, stepIndex, obsFlag
    integer :: ierr, headerIndex, numHeader, numHeaderMaxMpi
    integer :: headerIndexBeg, headerIndexEnd
    integer :: countObs, countObsInMpi, countObsOutMpi
    integer :: obsLonBurpFile, obsLatBurpFile, obsDate, obsTime
    integer :: numObsStnIdOut(numStnIdMax)
    integer :: numObsStnIdInMpi(numStnIdMax), numObsStnIdOutMpi(numStnIdMax)
    real(4) :: latInRadians, distance, obsLat, obsLon, gridLat, gridLon
    real(8) :: obsLatInDegrees, obsLonInDegrees, obsStepIndex_r8
    logical :: change
    real(4), allocatable :: gridLats(:), gridLatsMid(:), gridLonsMid(:,:)
    integer, allocatable :: headerIndexGrid(:,:,:), delMinutesGrid(:,:,:)
    real(4), allocatable :: distanceGrid(:,:,:)
    integer, allocatable :: obsLatIndex(:), obsLonIndex(:), obsStepIndex(:), numGridLons(:)
    real(4), allocatable :: obsDistance(:)
    integer, allocatable :: obsLatIndexMpi(:), obsLonIndexMpi(:), obsStepIndexMpi(:)
    integer, allocatable :: obsDelMinutes(:), obsDelMinutesMpi(:)
    integer, allocatable :: stnIdInt(:,:), stnIdIntMpi(:,:)
    real(4), allocatable :: obsDistanceMpi(:)
    logical, allocatable :: valid(:), validMpi(:)
    character(len=5)     :: stnIdTrim
    character(len=12)    :: stnId, stnidList(numStnIdMax)
    character(len=12), allocatable :: stnIdGrid(:,:,:)

    write(*,*)
    write(*,*) 'thn_scatByLatLonBoxes: Starting'
    write(*,*)

    numHeader = obs_numHeader(obsdat)
    call rpn_comm_allReduce(numHeader, numHeaderMaxMpi, 1, 'mpi_integer', &
                            'mpi_max','grid',ierr)
    write(*,*) 'thn_scatByLatLonBoxes: numHeader, numHeaderMaxMpi = ', &
               numHeader, numHeaderMaxMpi

    ! Check if we have any observations to process
    allocate(valid(numHeaderMaxMpi))
    valid(:) = .false.
    call obs_set_current_header_list(obsdat,'SC')
    HEADER0: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER0
      valid(headerIndex) = .true.
    end do HEADER0
    countObs = count(valid(:))
    call rpn_comm_allReduce(countObs, countObsInMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    if (countObsInMpi == 0) then
      write(*,*) 'thn_scatByLatLonBoxes: no observations for this instrument'
      deallocate(valid)
      return
    end if

    numLat = nint(2.*real(latLength)/real(deltax))
    numLon = nint(real(lonLength)/real(deltax))

    write(*,*)
    write(*,*) 'Number of horizontal boxes : ', numLon
    write(*,*) 'Number of vertical boxes   : ', numLat
    write(*,*) 'Number of temporal bins    : ', tim_nstepobs
    write(*,*)

    write(*,*) 'thn_scatByLatLonBoxes: countObs initial                   = ', &
               countObs, countObsInMpi

    ! Allocate arrays
    allocate(gridLats(numLat))
    allocate(gridLatsMid(numLat))
    allocate(gridLonsMid(numLat,numLon))
    allocate(numGridLons(numLat))
    allocate(stnIdGrid(numLat,numLon,tim_nstepobs))
    allocate(distanceGrid(numLat,numLon,tim_nstepobs))
    allocate(headerIndexGrid(numLat,numLon,tim_nstepobs))
    allocate(delMinutesGrid(numLat,numLon,tim_nstepobs))

    allocate(obsLatIndex(numHeaderMaxMpi))
    allocate(obsLonIndex(numHeaderMaxMpi))
    allocate(obsStepIndex(numHeaderMaxMpi))
    allocate(obsDistance(numHeaderMaxMpi))
    allocate(obsDelMinutes(numHeaderMaxMpi))
    lenStnId = len(stnId)
    allocate(stnIdInt(lenStnId,numHeaderMaxMpi))

    ! Allocation for MPI gather
    allocate(validMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsLatIndexMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsLonIndexMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsStepIndexMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsDistanceMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsDelMinutesMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(stnIdIntMpi(lenStnId,numHeaderMaxMpi*mpi_nprocs))

    gridLats(:)            = 0.
    gridLatsMid(:)         = 0.
    gridLonsMid(:,:)       = 0.
    numGridLons(:)         = 0
    stnIdGrid(:,:,:)       = ''
    distanceGrid(:,:,:)    = -1.0
    headerIndexGrid(:,:,:) = -1

    timeRejectCount = 0
    flagRejectCount = 0

    ! set spatial boxes properties
    do latIndex = 1, numLat
      gridLats(latIndex) = (latIndex*180./numLat) - 90.
      gridLatsMid(latIndex) = gridLats(latIndex) - (90./numLat)
      if (gridLats(latIndex) <= 0.0) then
        latInRadians = gridLats(latIndex) * MPC_PI_R8 / 180.
      else
        latInRadians = gridLats(latIndex-1) * MPC_PI_R8 / 180.
      end if
      distance = lonLength * cos(latInRadians)
      numGridLons(latIndex) = nint(distance/deltax)
      do lonIndex = 1, numGridLons(latIndex)
        gridLonsMid(latIndex,lonIndex) =  &
             (lonIndex * 36000 /  numGridLons(latIndex)) -  &
             (18000 / numGridLons(latIndex))
        gridLonsMid(latIndex,lonIndex) = 0.01 * gridLonsMid(latIndex,lonIndex)
      end do
    end do

    ! Station ID converted to integer array
    stnIdInt(:,:) = 0
    HEADER1: do headerIndex = 1, numHeader
      if (.not. valid(headerIndex)) cycle HEADER1

      stnId = obs_elem_c(obsdat,'STID',headerIndex)
      do charIndex = 1, lenStnId
        stnIdInt(charIndex,headerIndex) = iachar(stnId(charIndex:charIndex))
      end do
    end do HEADER1

    nsize = numHeaderMaxMpi
    call rpn_comm_allgather(valid,    nsize, 'mpi_logical',  &
                            validMpi, nsize, 'mpi_logical', 'grid', ierr)
    nsize = lenStnId * numHeaderMaxMpi
    call rpn_comm_allgather(stnIdInt,    nsize, 'mpi_integer',  &
                            stnIdIntMpi, nsize, 'mpi_integer', 'grid', ierr)

    ! build a global list of stnId over all mpi tasks
    numStnId = 0
    numObsStnIdInMpi(:) = 0
    HEADER2: do headerIndex = 1, numHeaderMaxMpi * mpi_nprocs
      if (all(stnIdIntMpi(:,headerIndex) == 0)) cycle HEADER2
      if (.not.validMpi(headerIndex)) cycle HEADER2

      ! Station ID converted back to character string
      do charIndex = 1, lenStnId
        stnId(charIndex:charIndex) = achar(stnIdIntMpi(charIndex,headerIndex))
      end do

      if (numStnId < numStnIdMax ) then
        stnIdIndexFound = -1
        do stnIdIndex = 1, numStnId
          if ( stnidList(stnIdIndex) == stnid ) stnIdIndexFound = stnIdIndex
        end do
        if ( stnIdIndexFound == -1 ) then
          numStnId = numStnId + 1
          stnidList(numStnId) = stnid
          stnIdIndexFound = numStnId
        end if
        numObsStnIdInMpi(stnIdIndexFound) = numObsStnIdInMpi(stnIdIndexFound) + 1
      else
        call utl_abort('thn_scatByLatLonBoxes: numStnId too large')
      end if
    end do HEADER2

    ! Initial pass through all observations
    HEADER3: do headerIndex = 1, numHeader
      if (.not. valid(headerIndex)) cycle HEADER3

      ! Station ID converted to integer array
      stnId = obs_elem_c(obsdat,'STID',headerIndex)
      do charIndex = 1, lenStnId
        stnIdInt(charIndex,headerIndex) = iachar(stnId(charIndex:charIndex))
      end do

      ! Obs lat-lon
      obsLonInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obs_headElem_r(obsdat, OBS_LON, headerIndex)
      obsLatInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obs_headElem_r(obsdat, OBS_LAT, headerIndex)
      obsLonBurpFile = nint(100.0*(obsLonInDegrees - 180.0))
      if(obsLonBurpFile < 0) obsLonBurpFile = obsLonBurpFile + 36000
      obsLatBurpFile = 9000+nint(100.0*obsLatInDegrees)

      ! compute box indices
      obsLat = (obsLatBurpFile - 9000.)/100.
      do latIndex = 1, numLat
        if ( obsLat <= (gridLats(latIndex) + 0.000001) ) then
          obsLatIndex(headerIndex) = latIndex
          exit
        end if
      end do

      if (obsLonBurpFile >= 18000) then
        obsLonBurpFile = obsLonBurpFile - 18000
      else
        obsLonBurpFile = obsLonBurpFile + 18000
      end if
      obsLonIndex(headerIndex) =  &
           int( obsLonBurpFile /  &
                (36000. / real(numGridLons(obsLatIndex(headerIndex)))) ) + 1
      if ( obsLonIndex(headerIndex) > numGridLons(obsLatIndex(headerIndex)) ) then
        obsLonIndex(headerIndex) = numGridLons(obsLatIndex(headerIndex))
      end if

      ! compute spatial distances
      obsLat = (obsLatBurpFile - 9000.)/100.
      obsLon = obsLonBurpFile/100.
      obsDistance(headerIndex) =  100.0 * &
           ( ((gridLatsMid(obsLatIndex(headerIndex))-obsLat))**2 +  &
             ((gridLonsMid(obsLatIndex(headerIndex),obsLonIndex(headerIndex))-obsLon))**2)**0.5

      ! calcul de la bin temporelle dans laquelle se trouve l'observation
      obsDate = obs_headElem_i(obsdat, OBS_DAT, headerIndex)
      obsTime = obs_headElem_i(obsdat, OBS_ETM, headerIndex)
      call tim_getStepObsIndex(obsStepIndex_r8, tim_getDatestamp(), &
                               obsDate, obsTime, tim_nstepobs)
      obsStepIndex(headerIndex) = nint(obsStepIndex_r8)
      obsDelMinutes(headerIndex) = nint( 60.0 * tim_dstepobs *  &
           abs(real(obsStepIndex(headerIndex)) - obsStepIndex_r8) )

      ! check time window
      if ( obsDelMinutes(headerIndex) > deltmax ) then
        timeRejectCount = timeRejectCount + 1
        valid(headerIndex) = .false.
      end if

      ! find observation flags (assumes 1 level only per headerIndex)
      uObsFlag = -1
      vObsFlag = -1
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY3: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY3
        obsVarno = obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex)
        if (obsVarno == bufr_neus) then
          uObsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
        else if (obsVarno == bufr_nevs) then
          vObsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
        end if
      end do BODY3

      ! modify valid based on flags
      if (uObsFlag /= -1 .and. vObsFlag /= -1) then
        if ( btest(uObsFlag,16) .or. btest(vObsFlag,16) .or. &
             btest(uObsFlag,18) .or. btest(vObsFlag,18) ) then
          flagRejectCount = flagRejectCount + 1
          valid(headerIndex) = .false.
        end if
      else
        valid(headerIndex) = .false.
      end if

    end do HEADER3

    countObs = count(valid(:))
    call rpn_comm_allReduce(countObs, countObsOutMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    write(*,*) 'thn_scatByLatLonBoxes: countObs after QC and time tests = ', &
               countObs, countObsOutMpi

    ! Gather data from all MPI tasks
    nsize = numHeaderMaxMpi
    call rpn_comm_allgather(valid,    nsize, 'mpi_logical',  &
                            validMpi, nsize, 'mpi_logical', 'grid', ierr)
    call rpn_comm_allgather(obsLatIndex,    nsize, 'mpi_integer',  &
                            obsLatIndexMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsLonIndex,    nsize, 'mpi_integer',  &
                            obsLonIndexMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsStepIndex,    nsize, 'mpi_integer',  &
                            obsStepIndexMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsDelMinutes,    nsize, 'mpi_integer',  &
                            obsDelMinutesMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsDistance,    nsize, 'mpi_real4',  &
                            obsDistanceMpi, nsize, 'mpi_real4', 'grid', ierr)
    
    ! Apply thinning algorithm
    HEADER4: do headerIndex = 1, numHeaderMaxMpi*mpi_nprocs
      if (.not. validMpi(headerIndex)) cycle HEADER4

      change = .true.

      latIndex  = obsLatIndexMpi(headerIndex)
      lonIndex  = obsLonIndexMpi(headerIndex)
      stepIndex = obsStepIndexMpi(headerIndex)

      ! Station ID converted back to character string
      do charIndex = 1, lenStnId
        stnId(charIndex:charIndex) = achar(stnIdIntMpi(charIndex,headerIndex))
      end do
      stnIdTrim = stnId(2:6)

      if ( stnIdGrid(latIndex,lonIndex,stepIndex) /= '' ) then

        ! This is an ASCAT observation
        if (stnIdTrim == 'METOP') then

          ! si l'obs retenue precedemment etait un ASCAT, on poursuit l'investigation
          if ( stnIdTrim == stnIdGrid(latIndex,lonIndex,stepIndex) ) then
          
            ! si la difference temporelle est plus grande que celle deja retenue
            if ( obsDelMinutesMpi(headerIndex) >  &
                 delMinutesGrid(latIndex,lonIndex,stepIndex) ) then
              change = .false.
            else
              ! si la distance au centre de la boite est plus grande que celle retenue 
              if ( (obsDelMinutesMpi(headerIndex) ==  &
                    delMinutesGrid(latIndex,lonIndex,stepIndex)) .and. &
                   (obsDistanceMpi(headerIndex) >=  &
                    distanceGrid(latIndex,lonIndex,stepIndex)) ) then
                change = .false.
              end if                    
            end if

          else

            ! si l'obs retenue precedemment etait autre que METOP
            change = .true.

          end if

        else ! satellites autre que METOP

          ! si l'obs retenue precedemment etait autre qu'un METOP, on poursuit l'investigation
          if ( stnIdTrim == stnIdGrid(latIndex,lonIndex,stepIndex) ) then
                 
            ! si la difference temporelle est plus grande que celle deja retenue, skip it
            if ( obsDelMinutesMpi(headerIndex) >  &
                 delMinutesGrid(latIndex,lonIndex,stepIndex) ) then
              change = .false.
            else
              ! si la distance au centre de la boite est plus grande que celle retenue, skip it 
              if ( (obsDelMinutesMpi(headerIndex) ==  &
                    delMinutesGrid(latIndex,lonIndex,stepIndex)) .and. &
                   (obsDistanceMpi(headerIndex) >=  &
                    distanceGrid(latIndex,lonIndex,stepIndex)) ) then
                change = .false.
              end if
            end if

          else

            ! si l'obs retenue precedemment etait un METOP
            change = .false.

          end if

        end if ! METOP

      end if

      ! update list of data to save
      if ( .not. change ) then
        ! keep previously accepted obs, so reject current obs
        validMpi(headerIndex) = .false.
      else
        ! reject previously accepted obs
        if ( headerIndexGrid(latIndex,lonIndex,stepIndex) /= -1 ) then
          validMpi(headerIndexGrid(latIndex,lonIndex,stepIndex)) = .false.
        end if

        ! keep current obs
        validMpi(headerIndex) = .true.
        headerIndexGrid(latIndex,lonIndex,stepIndex) = headerIndex
        stnIdGrid(latIndex,lonIndex,stepIndex) = stnIdTrim
        delMinutesGrid(latIndex,lonIndex,stepIndex) = obsDelMinutesMpi(headerIndex)
        distanceGrid(latIndex,lonIndex,stepIndex) = obsDistanceMpi(headerIndex)
      end if

    end do HEADER4

    ! update local copy of 'valid' array
    headerIndexBeg = 1 + mpi_myid * numHeaderMaxMpi
    headerIndexEnd = headerIndexBeg + numHeaderMaxMpi - 1
    valid(:) = validMpi(headerIndexBeg:headerIndexEnd)

    countObs = count(valid(:))
    call rpn_comm_allReduce(countObs, countObsOutMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    write(*,*) 'thn_scatByLatLonBoxes: countObs after choosing 1 per box  = ', &
               countObs, countObsOutMpi

    ! modify the observation flags in obsSpaceData and count obs for each stnId
    numObsStnIdOut(:) = 0
    call obs_set_current_header_list(obsdat,'SC')
    HEADER5: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER5
     
      if (.not. valid(headerIndex)) then
        call obs_set_current_body_list(obsdat, headerIndex)
        BODY5: do 
          bodyIndex = obs_getBodyIndex(obsdat)
          if (bodyIndex < 0) exit BODY5
        
          obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
          call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,11))

        end do BODY5
        cycle HEADER5
      end if

      ! count number of obs kept for each stnId
      stnId = obs_elem_c(obsdat,'STID',headerIndex)
      stnIdIndexFound = -1
      do stnIdIndex = 1, numStnId
        if (stnidList(stnIdIndex) == stnId) stnIdIndexFound = stnIdIndex
      end do
      if (stnIdIndexFound == -1) call utl_abort('stnid not found in list')
      numObsStnIdOut(stnIdIndexFound) = numObsStnIdOut(stnIdIndexFound) + 1
    end do HEADER5

    call rpn_comm_allReduce(numObsStnIdOut, numObsStnIdOutMpi, &
                            numStnIdMax, 'mpi_integer', 'mpi_sum', 'grid', ierr)
    call rpn_comm_allReduce(timeRejectCount, timeRejectCountMpi, 1, &
                            'mpi_integer', 'mpi_sum', 'grid', ierr)
    call rpn_comm_allReduce(flagRejectCount, flagRejectCountMpi, 1, &
                            'mpi_integer', 'mpi_sum', 'grid', ierr)

    write(*,*)
    write(*,'(a,i6)') ' Number of obs in input ', countObsInMpi
    write(*,'(a,i6)') ' Number of obs in output ', countObsOutMpi
    write(*,'(a,i6)') ' Number of obs not selected due to time ', timeRejectCountMpi
    write(*,'(a,i6)') ' Number of obs not selected due to topo ', flagRejectCountMpi
    write(*,*)
    
    write(*,'(a40,i10)' ) 'Number of satellites found = ', numStnId
    write(*,*)
  
    write(*,'(a40,2a15)' ) 'Satellite', 'nb SCAT in'
    write(*,*)
    do stnIdIndex = 1, numStnId
      write(*,'(a40,2i15)') stnidList(stnIdIndex), numObsStnIdInMpi(stnIdIndex)
    end do
    write(*,*)
    write(*,'(a40,2i10,f10.4)' ) 'Total number of obs in : ', sum(numObsStnIdInMpi(:))
  
    write(*,*)
    write(*,'(a40,2a15)' ) 'Satellite', 'nb SCAT out'
    write(*,*)
    do stnIdIndex = 1, numStnId
      write(*,'(a40,2i15)') stnidList(stnIdIndex), numObsStnIdOutMpi(stnIdIndex)
    end do
    write(*,*)
    write(*,'(a40,2i10,f10.4)' ) 'Total number of obs out : ', sum(numObsStnIdOutMpi(:))

    ! Deallocations:
    deallocate(valid)
    deallocate(gridLats)
    deallocate(gridLatsMid)
    deallocate(gridLonsMid)
    deallocate(numGridLons)
    deallocate(stnIdGrid)
    deallocate(distanceGrid)
    deallocate(headerIndexGrid)
    deallocate(delMinutesGrid)

    deallocate(obsLatIndex)
    deallocate(obsLonIndex)
    deallocate(obsStepIndex)
    deallocate(obsDistance)
    deallocate(obsDelMinutes)
    deallocate(stnIdInt)

    deallocate(validMpi)
    deallocate(obsLatIndexMpi)
    deallocate(obsLonIndexMpi)
    deallocate(obsStepIndexMpi)
    deallocate(obsDistanceMpi)
    deallocate(obsDelMinutesMpi)
    deallocate(stnIdIntMpi)

    write(*,*)
    write(*,*) 'thn_scatByLatLonBoxes: Finished'
    write(*,*)

  end subroutine thn_scatByLatLonBoxes

  !--------------------------------------------------------------------------
  ! thn_csrByLatLonBoxes
  !--------------------------------------------------------------------------
  subroutine thn_csrByLatLonBoxes(obsdat, deltax, deltrad)
    !
    ! :Purpose: Only keep the observation closest to the center of each
    !           lat-lon (and time) box for CSR observations.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    integer, intent(in)             :: deltax
    integer, intent(in)             :: deltrad

    ! Locals parameters:
    integer, parameter :: latLength = 10000 ! Earth dimension parameters
    integer, parameter :: lonLength = 40000 ! Earth dimension parameters
    integer, parameter :: maxNumChan = 15    ! nb max de canaux

    ! Locals:
    integer :: bodyIndex, channelIndex, charIndex, nsize, lenStnId
    integer :: numLat, numLon, latIndex, lonIndex, stepIndex, obsFlag
    integer :: ierr, headerIndex, numHeader, numHeaderMaxMpi, channelList(maxNumChan)
    integer :: headerIndexBeg, headerIndexEnd, countObs, countObsMpi
    integer :: obsLonBurpFile, obsLatBurpFile, obsDate, obsTime
    real(4) :: latInRadians, distance, obsLat, obsLon, gridLat, gridLon
    real(8) :: obsLatInDegrees, obsLonInDegrees, obsStepIndex_r8
    logical :: change
    real(4), allocatable :: gridLats(:)
    integer, allocatable :: numChanAssimGrid(:,:,:), headerIndexGrid(:,:,:)
    real(4), allocatable :: angleGrid(:,:,:), distanceGrid(:,:,:), cloudGrid(:,:,:,:)
    integer, allocatable :: obsLatIndex(:), obsLonIndex(:), obsStepIndex(:), numGridLons(:)
    real(4), allocatable :: obsCloud(:,:), obsAngle(:), obsDistance(:)
    integer, allocatable :: obsLatIndexMpi(:), obsLonIndexMpi(:), obsStepIndexMpi(:)
    integer, allocatable :: stnIdInt(:,:), stnIdIntMpi(:,:), numChannel(:), numChannelMpi(:)
    real(4), allocatable :: obsCloudMpi(:,:), obsAngleMpi(:), obsDistanceMpi(:)
    logical, allocatable :: valid(:), validMpi(:), channelAssim(:,:), channelAssimMpi(:,:)
    character(len=12) :: stnId
    character(len=12), allocatable :: stnIdGrid(:,:,:)

    write(*,*)
    write(*,*) 'thn_csrByLatLonBoxes: Starting'
    write(*,*)

    numHeader = obs_numHeader(obsdat)
    call rpn_comm_allReduce(numHeader, numHeaderMaxMpi, 1, 'mpi_integer', &
                            'mpi_max','grid',ierr)
    write(*,*) 'thn_csrByLatLonBoxes: numHeader, numHeaderMaxMpi = ', &
               numHeader, numHeaderMaxMpi

    ! Check if we have any observations to process
    allocate(valid(numHeaderMaxMpi))
    valid(:) = .false.
    do headerIndex = 1, numHeader
      if ( obs_headElem_i(obsdat, OBS_ITY, headerIndex) == &
           codtyp_get_codtyp('radianceclear') ) then
        valid(headerIndex) = .true.
      end if
    end do
    countObs = count(valid(:))
    call rpn_comm_allReduce(countObs, countObsMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    if (countObsMpi == 0) then
      write(*,*) 'thn_csrByLatLonBoxes: no observations for this instrument'
      deallocate(valid)
      return
    end if

    numLat = nint(2.*real(latLength)/real(deltax))
    numLon = nint(real(lonLength)/real(deltax))

    write(*,*)
    write(*,*) 'Number of horizontal boxes : ', numLon
    write(*,*) 'Number of vertical boxes   : ', numLat
    write(*,*) 'Number of temporal bins    : ', tim_nstepobs
    write(*,*)

    write(*,*) 'thn_csrByLatLonBoxes: countObs initial                   = ', &
               countObs, countObsMpi

    ! Allocate arrays
    allocate(gridLats(numLat))
    allocate(numGridLons(numLat))
    allocate(stnIdGrid(numLat,numLon,tim_nstepobs))
    allocate(numChanAssimGrid(numLat,numLon,tim_nstepobs))
    allocate(angleGrid(numLat,numLon,tim_nstepobs))
    allocate(distanceGrid(numLat,numLon,tim_nstepobs))
    allocate(headerIndexGrid(numLat,numLon,tim_nstepobs))
    allocate(cloudGrid(maxNumChan,numLat,numLon,tim_nstepobs))
    allocate(obsLatIndex(numHeaderMaxMpi))
    allocate(obsLonIndex(numHeaderMaxMpi))
    allocate(obsStepIndex(numHeaderMaxMpi))
    allocate(numChannel(numHeaderMaxMpi))
    allocate(channelAssim(maxNumChan,numHeaderMaxMpi))
    allocate(obsAngle(numHeaderMaxMpi))
    allocate(obsCloud(maxNumChan,numHeaderMaxMpi))
    allocate(obsDistance(numHeaderMaxMpi))

    gridLats(:)             = 0.
    numGridLons(:)          = 0
    stnIdGrid(:,:,:)        = ''
    numChanAssimGrid(:,:,:) = -1
    angleGrid(:,:,:)        = -1.0
    distanceGrid(:,:,:)     = -1.0
    headerIndexGrid(:,:,:)  = -1
    cloudGrid(:,:,:,:)      = -1.0
    channelAssim(:,:)       = .false.
    numChannel(:)           = 0

    ! set spatial boxes properties

    do latIndex = 1, numLat
      gridLats(latIndex) = (latIndex*180./numLat) - 90.
      if (gridLats(latIndex) <= 0.0) then
        latInRadians = gridLats(latIndex) * MPC_PI_R8 / 180.
      else
        latInRadians = gridLats(latIndex-1) * MPC_PI_R8 / 180.
      end if
      distance = lonLength * cos(latInRadians)
      numGridLons(latIndex) = nint(distance/deltax)
    end do

    ! Initial pass through all observations
    HEADER1: do headerIndex = 1, numHeader
      if (.not. valid(headerIndex)) cycle HEADER1

      obsLonInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obs_headElem_r(obsdat, OBS_LON, headerIndex)
      obsLatInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obs_headElem_r(obsdat, OBS_LAT, headerIndex)
      obsLonBurpFile = nint(100.0*(obsLonInDegrees - 180.0))
      if(obsLonBurpFile < 0) obsLonBurpFile = obsLonBurpFile + 36000
      obsLatBurpFile = 9000+nint(100.0*obsLatInDegrees)

      ! compute box indices
      do latIndex = 1, numLat
        if ( (obsLatBurpFile - 9000.)/100. <= (gridLats(latIndex) + 0.000001) ) then
          obsLatIndex(headerIndex) = latIndex
          exit
        end if
      end do

      obsLonIndex(headerIndex) = int(obsLonBurpFile /  &
           (36000. / numGridLons(obsLatIndex(headerIndex)))) + 1
      if ( obsLonIndex(headerIndex) > numGridLons(obsLatIndex(headerIndex)) ) then
        obsLonIndex(headerIndex) = numGridLons(obsLatIndex(headerIndex))
      end if

      ! compute spatial distances
      ! position of the observation
      obsLat = (obsLatBurpFile - 9000.) / 100.
      obsLon = obsLonBurpFile / 100.

      ! position of the box center
      gridLat = gridLats(obsLatIndex(headerIndex)) - 0.5 * (180./numLat)
      gridLon = (360. / numGridLons(obsLatIndex(headerIndex))) *  &
                (obsLonIndex(headerIndex) - 0.5)

      ! spatial separation
      obsDistance(headerIndex) = thn_separation(obsLon,obsLat,gridLon,gridLat) * &
                                 latLength / 90.

      ! calcul de la bin temporelle dans laquelle se trouve l'observation
      obsDate = obs_headElem_i(obsdat, OBS_DAT, headerIndex)
      obsTime = obs_headElem_i(obsdat, OBS_ETM, headerIndex)
      call tim_getStepObsIndex(obsStepIndex_r8, tim_getDatestamp(), &
                               obsDate, obsTime, tim_nstepobs)
      obsStepIndex(headerIndex) = nint(obsStepIndex_r8)

      ! check if distance too far from box center
      if (obsDistance(headerIndex) > real(deltrad)) then
        valid(headerIndex) = .false.
        cycle HEADER1
      end if

    end do HEADER1

    countObs = count(valid(:))
    call rpn_comm_allReduce(countObs, countObsMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    write(*,*) 'thn_csrByLatLonBoxes: countObs after deltrad test        = ', &
               countObs, countObsMpi

    ! Second pass through all observations
    HEADER2: do headerIndex = 1, numHeader
      if (.not. valid(headerIndex)) cycle HEADER2
      
      ! get the zenith angle
      obsAngle(headerIndex) = obs_headElem_r(obsdat, OBS_SZA, headerIndex)

      ! Keep obs only if at least one channel not rejected based on tests in suprep
      valid(headerIndex) = .false.
      channelIndex = 0
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY1: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY1

        if (obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex) /= bufr_nbt3) then
          cycle BODY1
        end if

        numChannel(headerIndex) = numChannel(headerIndex) + 1
        channelIndex = channelIndex + 1
        channelList(channelIndex) = nint(obs_bodyElem_r(obsdat, OBS_PPP, bodyIndex))

        if (obs_bodyElem_i(obsdat, OBS_ASS, bodyIndex) == obs_assimilated) then
          valid(headerIndex) = .true.
          channelAssim(channelIndex,headerIndex) = .true.
        end if
      end do BODY1

      ! extract cloud fraction
      CHANNELS: do channelIndex = 1, numChannel(headerIndex)
        obsCloud(channelIndex, headerIndex) = -1.0

        ! search for the cloud information for this channel
        call obs_set_current_body_list(obsdat, headerIndex)
        BODY2: do
          bodyIndex = obs_getBodyIndex(obsdat)
          if (bodyIndex < 0) exit BODY2

          ! skip if not this is not cloud
          if (obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex) /= bufr_cloudInSeg) then
            cycle BODY2
          end if
          ! check if channel number matches
          if (nint(obs_bodyElem_r(obsdat, OBS_PPP, bodyIndex)) == channelList(channelIndex)) then
            obsCloud(channelIndex, headerIndex) = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex)
            cycle CHANNELS
          end if
        end do BODY2
        if (obsCloud(channelIndex, headerIndex) == -1.0) then
          call utl_abort('thn_csrByLatLonBoxes: could not find cloud fraction in obsSpaceData')
        end if
      end do CHANNELS

    end do HEADER2

    countObs = count(valid(:))
    call rpn_comm_allReduce(countObs, countObsMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    write(*,*) 'thn_csrByLatLonBoxes: countObs after rejection flag test = ', &
               countObs, countObsMpi

    ! Allocation for MPI gather
    allocate(validMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsLatIndexMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsLonIndexMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsStepIndexMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(numChannelMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsAngleMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(obsDistanceMpi(numHeaderMaxMpi*mpi_nprocs))
    allocate(channelAssimMpi(maxNumChan,numHeaderMaxMpi*mpi_nprocs))
    allocate(obsCloudMpi(maxNumChan,numHeaderMaxMpi*mpi_nprocs))
    lenStnId = len(stnId)
    allocate(stnIdInt(lenStnId,numHeaderMaxMpi))
    allocate(stnIdIntMpi(lenStnId,numHeaderMaxMpi*mpi_nprocs))

    ! Station ID converted to integer array
    HEADER3: do headerIndex = 1, numHeader
      if (.not. valid(headerIndex)) cycle HEADER3

      stnId = obs_elem_c(obsdat,'STID',headerIndex)
      do charIndex = 1, lenStnId
        stnIdInt(charIndex,headerIndex) = iachar(stnId(charIndex:charIndex))
      end do
    end do HEADER3

    ! Gather data from all MPI tasks
    nsize = numHeaderMaxMpi
    call rpn_comm_allgather(valid,    nsize, 'mpi_logical',  &
                            validMpi, nsize, 'mpi_logical', 'grid', ierr)
    call rpn_comm_allgather(obsLatIndex,    nsize, 'mpi_integer',  &
                            obsLatIndexMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsLonIndex,    nsize, 'mpi_integer',  &
                            obsLonIndexMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsStepIndex,    nsize, 'mpi_integer',  &
                            obsStepIndexMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(numChannel,    nsize, 'mpi_integer',  &
                            numChannelMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsAngle,    nsize, 'mpi_real4',  &
                            obsAngleMpi, nsize, 'mpi_real4', 'grid', ierr)
    call rpn_comm_allgather(obsDistance,    nsize, 'mpi_real4',  &
                            obsDistanceMpi, nsize, 'mpi_real4', 'grid', ierr)

    nsize = maxNumChan * numHeaderMaxMpi
    call rpn_comm_allgather(channelAssim,    nsize, 'mpi_integer',  &
                            channelAssimMpi, nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(obsCloud,    nsize, 'mpi_real4',  &
                            obsCloudMpi, nsize, 'mpi_real4', 'grid', ierr)

    nsize = lenStnId * numHeaderMaxMpi
    call rpn_comm_allgather(stnIdInt,    nsize, 'mpi_integer',  &
                            stnIdIntMpi, nsize, 'mpi_integer', 'grid', ierr)
    
    ! Apply thinning algorithm
    HEADER4: do headerIndex = 1, numHeaderMaxMpi*mpi_nprocs
      if (.not. validMpi(headerIndex)) cycle HEADER4

      change = .true.

      latIndex  = obsLatIndexMpi(headerIndex)
      lonIndex  = obsLonIndexMpi(headerIndex)
      stepIndex = obsStepIndexMpi(headerIndex)

      ! Station ID converted back to character string
      do charIndex = 1, lenStnId
        stnId(charIndex:charIndex) = achar(stnIdIntMpi(charIndex,headerIndex))
      end do

      if ( stnIdGrid(latIndex,lonIndex,stepIndex) /= '' ) then

        ! on veut previlegier les profils avec le plus de canaux assimiles
        if ( count(channelAssimMpi(:,headerIndex)) <  &
             numChanAssimGrid(latIndex,lonIndex,stepIndex) ) change = .false.

        ! en cas d'egalite, on doit regarder d'autres conditions pour faire un choix
        if ( count(channelAssimMpi(:,headerIndex)) ==  &
             numChanAssimGrid(latIndex,lonIndex,stepIndex) ) then

          ! si le profil actuel est d'un autre instrument que celui deja considere
          ! choisir celui qui a le plus petit angle satellite
          if ( stnid /= stnIdGrid(latIndex,lonIndex,stepIndex) ) then
            if ( obsAngleMpi(headerIndex)  >  &
                 angleGrid(latIndex,lonIndex,stepIndex) ) change = .false.

            ! en cas d'egalite de l'angle,
            ! choisir le profil le plus pres du centre de la boite
            if ( ( obsAngleMpi(headerIndex) ==  &
                   angleGrid(latIndex,lonIndex,stepIndex) ) .and. &
                 ( obsDistanceMpi(headerIndex) >  &
                   distanceGrid(latIndex,lonIndex,stepIndex) ) ) change = .false.

          ! si le profil actuel est du meme instrument que celui deja considere
          ! choisir celui dont tous les canaux assimiles ont respectivement
          ! moins de fraction nuageuse que celui deja considere
          else
            do channelIndex = 1, numChannelMpi(headerIndex)
              if ( channelAssimMpi(channelIndex,headerIndex) .and. &
                   ( obsCloudMpi(channelIndex,headerIndex) >  &
                     cloudGrid(channelIndex,latIndex,lonIndex,stepIndex) ) ) change = .false.
            end do

            ! en cas d'egalite de la fraction nuageuse pour chaque canal present,
            ! choisir le profil le plus pres du centre de la boite
            do channelIndex = 1, numChannelMpi(headerIndex)
              if ( channelAssimMpi(channelIndex,headerIndex) .and. &
                   ( obsCloudMpi(channelIndex,headerIndex) <  &
                     cloudGrid(channelIndex,latIndex,lonIndex,stepIndex) ) ) exit
              if ( ( channelIndex == numChannelMpi(headerIndex) ) .and. &
                   ( obsDistanceMpi(headerIndex) >  &
                     distanceGrid(latIndex,lonIndex,stepIndex) ) ) change = .false.
            end do
          end if
        end if
      end if

      ! update list of data to save
      if ( .not. change ) then
        ! keep previously accepted obs, so reject current obs
        validMpi(headerIndex) = .false.
      else
        ! reject previously accepted obs
        if ( headerIndexGrid(latIndex,lonIndex,stepIndex) /= -1 ) then
          validMpi(headerIndexGrid(latIndex,lonIndex,stepIndex)) = .false.
        end if

        ! keep current obs
        validMpi(headerIndex) = .true.
        headerIndexGrid(latIndex,lonIndex,stepIndex) = headerIndex
        numChanAssimGrid(latIndex,lonIndex,stepIndex) =  &
             count(channelAssimMpi(:,headerIndex))
        stnIdGrid(latIndex,lonIndex,stepIndex) = stnid
        angleGrid(latIndex,lonIndex,stepIndex) = obsAngleMpi(headerIndex)
        cloudGrid(:,latIndex,lonIndex,stepIndex) = obsCloudMpi(:,headerIndex)
        distanceGrid(latIndex,lonIndex,stepIndex) = obsDistanceMpi(headerIndex)
      end if

    end do HEADER4

    ! update local copy of 'valid' array
    headerIndexBeg = 1 + mpi_myid * numHeaderMaxMpi
    headerIndexEnd = headerIndexBeg + numHeaderMaxMpi - 1
    valid(:) = validMpi(headerIndexBeg:headerIndexEnd)

    countObs = count(valid(:))
    call rpn_comm_allReduce(countObs, countObsMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    write(*,*) 'thn_csrByLatLonBoxes: countObs after choosing 1 per box  = ', &
               countObs, countObsMpi

    ! modify the observation flags in obsSpaceData
    HEADER5: do headerIndex = 1, numHeader
      ! skip observation if we're not supposed to consider it
      if ( obs_headElem_i(obsdat, OBS_ITY, headerIndex) /= &
           codtyp_get_codtyp('radianceclear') ) then
        cycle HEADER5
      end if
     
      if (.not. valid(headerIndex)) then
        call obs_set_current_body_list(obsdat, headerIndex)
        BODY3: do 
          bodyIndex = obs_getBodyIndex(obsdat)
          if (bodyIndex < 0) exit BODY3
        
          obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
          call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,11))

        end do BODY3
      end if
    end do HEADER5

    deallocate(valid)
    deallocate(gridLats)
    deallocate(numGridLons)
    deallocate(stnIdGrid)
    deallocate(numChanAssimGrid)
    deallocate(angleGrid)
    deallocate(distanceGrid)
    deallocate(headerIndexGrid)
    deallocate(cloudGrid)
    deallocate(obsLatIndex)
    deallocate(obsLonIndex)
    deallocate(obsStepIndex)
    deallocate(numChannel)
    deallocate(channelAssim)
    deallocate(obsAngle)
    deallocate(obsCloud)
    deallocate(obsDistance)
    deallocate(validMpi)
    deallocate(obsLatIndexMpi)
    deallocate(obsLonIndexMpi)
    deallocate(obsStepIndexMpi)
    deallocate(numChannelMpi)
    deallocate(obsAngleMpi)
    deallocate(obsDistanceMpi)
    deallocate(channelAssimMpi)
    deallocate(obsCloudMpi)
    deallocate(stnIdInt)
    deallocate(stnIdIntMpi)

    write(*,*)
    write(*,*) 'thn_csrByLatLonBoxes: Finished'
    write(*,*)

  end subroutine thn_csrByLatLonBoxes

  !--------------------------------------------------------------------------
  ! thn_hyperByLatLonBoxes
  !--------------------------------------------------------------------------
  subroutine thn_hyperByLatLonBoxes(obsdat, removeUnCorrected, &
                                    deltmax, deltax, deltrad,  &
                                    familyType, codtyp)
    !
    ! :Purpose: Only keep the observation closest to the center of each
    !           lat-lon (and time) box.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    logical, intent(in)             :: removeUnCorrected
    integer, intent(in)             :: deltmax
    integer, intent(in)             :: deltax
    integer, intent(in)             :: deltrad
    character(len=*), intent(in)    :: familyType
    integer, intent(in)             :: codtyp

    ! Locals:
    integer :: headerIndex, bodyIndex, obsDate, obsTime, obsFlag
    integer :: numLat, numLon, latIndex, numChannels, delMinutes
    integer :: lonBinIndex, latBinIndex, timeBinIndex
    integer :: ierr, nsize, procIndex, countHeader
    real(4) :: latInRadians, length, distance
    real(8) :: lonBoxCenterInDegrees, latBoxCenterInDegrees
    real(8) :: obsLatInRad, obsLonInRad, obsLat, obsLon
    real(8) :: obsLatInDegrees, obsLonInDegrees, obsStepIndex_r8
    integer, parameter :: latLength = 10000
    integer, parameter :: lonLength = 40000
    real(4), allocatable :: gridLats(:)
    integer, allocatable :: numGridLons(:)
    integer, allocatable :: headerIndexKeep(:,:,:), numChannelsKeep(:,:,:)
    integer, allocatable :: headerIndexKeepMpi(:,:,:,:), numChannelsKeepMpi(:,:,:,:)
    integer, allocatable :: delMinutesKeep(:,:,:), delMinutesKeepMpi(:,:,:,:)
    integer, allocatable :: procIndexKeep(:,:,:)
    real(4), allocatable :: distanceKeep(:,:,:), distanceKeepMpi(:,:,:,:)
    logical, allocatable :: rejectThisHeader(:)
    logical :: keepThisObs
    integer :: obsLonBurpFile, obsLatBurpFile
    character(len=12) :: stnid

    write(*,*)
    write(*,*) 'thn_hyperByLatLonBoxes: Starting, ', trim(codtyp_get_name(codtyp))
    write(*,*)

    ! Initial setup
    numLat = nint( 2. * real(latLength) / real(deltax) )
    numLon = nint(      real(lonLength) / real(deltax) )
    allocate(headerIndexKeep(numLat,numLon,tim_nstepobs))
    allocate(numChannelsKeep(numLat,numLon,tim_nstepobs))
    allocate(distanceKeep(numLat,numLon,tim_nstepobs))
    allocate(delMinutesKeep(numLat,numLon,tim_nstepobs))
    allocate(gridLats(numLat))
    allocate(numGridLons(numLat))
    headerIndexKeep(:,:,:) = -1
    numChannelsKeep(:,:,:) = 0
    distanceKeep(:,:,:)    = 0.0
    delMinutesKeep(:,:,:)  = deltmax
    gridLats(:)              = 0.0
    numGridLons(:)                = 0

    allocate(headerIndexKeepMpi(numLat,numLon,tim_nstepobs,mpi_nprocs))
    allocate(numChannelsKeepMpi(numLat,numLon,tim_nstepobs,mpi_nprocs))
    allocate(distanceKeepMpi(numLat,numLon,tim_nstepobs,mpi_nprocs))
    allocate(delMinutesKeepMpi(numLat,numLon,tim_nstepobs,mpi_nprocs))
    allocate(procIndexKeep(numLat,numLon,tim_nstepobs))
    procIndexKeep(:,:,:) = -1

    ! set spatial boxes properties
    ! gridLats(:) : latitude (deg) of northern side of the box
    ! numGridLons(:)   : number of longitudinal boxes at this latitude
    do latIndex = 1, numLat
      gridLats(latIndex) = (latIndex*180./numLat) - 90.
      if ( gridLats(latIndex) <= 0.0 ) then
        latInRadians = gridLats(latIndex) * MPC_PI_R8 / 180.
      else
        latInRadians = gridLats(latIndex-1) * MPC_PI_R8 / 180.
      end if
      length = lonLength * cos(latInRadians)
      numGridLons(latIndex)   = nint(length/deltax)
    end do

    countHeader = 0

    ! loop over all header indices of the specified family
    call obs_set_current_header_list(obsdat,trim(familyType))
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER

      if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) /= codtyp) cycle HEADER

      countHeader = countHeader + 1

      obsLonInRad = obs_headElem_r(obsdat, OBS_LON, headerIndex)
      obsLatInRad = obs_headElem_r(obsdat, OBS_LAT, headerIndex)
      obsDate = obs_headElem_i(obsdat, OBS_DAT, headerIndex)
      obsTime = obs_headElem_i(obsdat, OBS_ETM, headerIndex)

      obsLonInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obsLonInRad
      obsLatInDegrees = MPC_DEGREES_PER_RADIAN_R8 * obsLatInRad
      obsLonBurpFile = nint(100.0*(obsLonInDegrees - 180.0))
      if(obsLonBurpFile < 0) obsLonBurpFile = obsLonBurpFile + 36000
      obsLatBurpFile = 9000+nint(100.0*obsLatInDegrees)

      numChannels = 0

      ! loop over all body indices for this headerIndex
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY

        ! mark for rejection if not bias corrected (bit 6 not set)
        if (removeUnCorrected) then
          obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
          if (.not. btest(obsFlag,6)) then
            call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,11))
            call obs_bodySet_i(obsdat, OBS_ASS, bodyIndex, obs_notAssimilated)
          end if
        end if

        ! count the number of accepted channels
        if (obs_bodyElem_i(obsdat, OBS_ASS, bodyIndex) == obs_assimilated) then
          numChannels = numChannels + 1
        end if

      end do BODY

      ! Determine the lat and lon bin indexes
      do latIndex = 1, numLat
        if ( obsLatInDegrees <= (gridLats(latIndex)+0.000001) ) then
          latBinIndex = latIndex
          exit
        end if
      end do
      lonBinIndex = int( obsLonBurpFile/(36000.0/numGridLons(latBinIndex)) ) + 1
      if ( lonBinIndex > numGridLons(latBinIndex) ) lonBinIndex = numGridLons(latBinIndex)

      ! Determine the time bin index
      call tim_getStepObsIndex(obsStepIndex_r8, tim_getDatestamp(), &
                               obsDate, obsTime, tim_nstepobs)
      timeBinIndex = nint(obsStepIndex_r8)
      delMinutes = nint(60.0 * tim_dstepobs * abs(real(timeBinIndex) - obsStepIndex_r8))

      ! Determine distance from box center
      latBoxCenterInDegrees = gridLats(latBinIndex) - 0.5 * (180./numLat)
      lonBoxCenterInDegrees = (360. / numGridLons(latBinIndex)) * (lonBinIndex - 0.5)
      obsLat = (obsLatBurpFile - 9000.) / 100.
      obsLon = obsLonBurpFile / 100.
      distance = 1.0d-3 * phf_calcDistance(MPC_RADIANS_PER_DEGREE_R8 * latBoxCenterInDegrees, &
                                           MPC_RADIANS_PER_DEGREE_R8 * lonBoxCenterInDegrees, &
                                           MPC_RADIANS_PER_DEGREE_R8 * obsLat, &
                                           MPC_RADIANS_PER_DEGREE_R8 * obsLon )

      stnid = obs_elem_c(obsdat,'STID',headerIndex)

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
          if ( delMinutes > delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
            keepThisObs = .false.
          end if
         
          if ( delMinutes == delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
            if ( distance > distanceKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
              keepThisObs = .false.
            end if
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

    ! return if no observations for this instrument
    if (countHeader == 0) then
      call deallocLocals()
      write(*,*) 'thn_hyperByLatLonBoxes: no observations for this instrument'
      return
    end if

    ! communicate results to all other mpi tasks
    nsize = numLat * numLon * tim_nstepobs
    call rpn_comm_allgather(distanceKeep,     nsize, 'mpi_real4',  &
                            distanceKeepMpi,  nsize, 'mpi_real4', 'grid', ierr)
    call rpn_comm_allgather(delMinutesKeep,     nsize, 'mpi_integer',  &
                            delMinutesKeepMpi,  nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(numChannelsKeep,     nsize, 'mpi_integer',  &
                            numChannelsKeepMpi,  nsize, 'mpi_integer', 'grid', ierr)
    call rpn_comm_allgather(headerIndexKeep,     nsize, 'mpi_integer',  &
                            headerIndexKeepMpi,  nsize, 'mpi_integer', 'grid', ierr)

    ! reset arrays that store info about kept obs
    headerIndexKeep(:,:,:) = -1
    numChannelsKeep(:,:,:) = 0
    distanceKeep(:,:,:)    = 0.0
    delMinutesKeep(:,:,:)  = deltmax

    do timeBinIndex = 1, tim_nstepobs
      do lonBinIndex = 1, numLon
        do latBinIndex = 1, numLat

          ! Apply thinning criteria to results from all mpi tasks
          do procIndex = 1, mpi_nprocs

            headerIndex = headerIndexKeepMpi(latBinIndex,lonBinIndex,timeBinIndex,procIndex)
            distance    = distanceKeepMpi(latBinIndex,lonBinIndex,timeBinIndex,procIndex)
            delMinutes  = delMinutesKeepMpi(latBinIndex,lonBinIndex,timeBinIndex,procIndex)
            numChannels = numChannelsKeepMpi(latBinIndex,lonBinIndex,timeBinIndex,procIndex)
            
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
                if ( delMinutes > delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
                  keepThisObs = .false.
                end if
         
                if ( delMinutes == delMinutesKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
                  if ( distance > distanceKeep(latBinIndex,lonBinIndex,timeBinIndex) ) then
                    keepThisObs = .false.
                  end if
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
      do lonBinIndex = 1, numLon
        do latBinIndex = 1, numLat
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

      if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) /= codtyp) cycle HEADER2

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
    deallocate(rejectThisHeader)

    call deallocLocals()

    write(*,*)
    write(*,*) 'thn_hyperByLatLonBoxes: Finished'
    write(*,*)

  contains

    subroutine deallocLocals()
      implicit none

      deallocate(headerIndexKeep)
      deallocate(numChannelsKeep)
      deallocate(distanceKeep)
      deallocate(delMinutesKeep)
      deallocate(gridLats)
      deallocate(numGridLons)
      deallocate(headerIndexKeepMpi)
      deallocate(numChannelsKeepMpi)
      deallocate(distanceKeepMpi)
      deallocate(delMinutesKeepMpi)
      deallocate(procIndexKeep)
      
    end subroutine deallocLocals

  end subroutine thn_hyperByLatLonBoxes

  !--------------------------------------------------------------------------
  ! thn_separation
  !--------------------------------------------------------------------------
  function thn_separation(xlon1,xlat1,xlon2,xlat2)
    implicit none

    ! Arguments:
    real(4) :: thn_separation
    real(4) :: xlat1, xlat2, xlon1, xlon2

    ! Locals:
    real(4) :: cosval, degrad, raddeg

    raddeg = 180.0/3.14159265358979
    degrad = 1.0/raddeg
    cosval = sin(xlat1*degrad) * sin(xlat2*degrad) + &
             cos(xlat1*degrad) * cos(xlat2*degrad) * &
             cos((xlon1-xlon2) * degrad)

    if (cosval < -1.0d0) then
      cosval = -1.0d0
    else if (cosval > 1.0d0) then
      cosval = 1.0d0
    end if
    thn_separation = acos(cosval) * raddeg

  end function thn_separation

end module thinning_mod

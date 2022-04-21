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
  use gridStateVectorFileIO_mod
  use codtyp_mod
  use physicsFunctions_mod
  use utilities_mod
  use kdtree2_mod
  implicit none
  private

  public :: thn_thinHyper, thn_thinTovs, thn_thinCSR
  public :: thn_thinRaobs, thn_thinAircraft, thn_thinScat, thn_thinSatWinds
  public :: thn_thinSurface, thn_thinGbGps, thn_thinGpsRo, thn_thinAladin

  integer, external :: get_max_rss

contains

  !--------------------------------------------------------------------------
  ! thn_thinSurface
  !--------------------------------------------------------------------------
  subroutine thn_thinSurface(obsdat)
    ! :Purpose: Main subroutine for thinning of surface obs.

    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat

    ! Locals:
    integer :: nulnam
    integer :: fnom, fclos, ierr

    ! Namelist variables
    logical :: doThinning   ! if false, we return immediately
    real(8) :: step         ! time resolution (in hours)
    integer :: deltmax      ! maximum time difference (in minutes)
    logical :: useBlackList ! signal if blacklist file should be read and used

    namelist /thin_surface/doThinning, step, deltmax, useBlackList
    
    ! return if no surface obs
    if (.not. obs_famExist(obsdat,'SF')) return

    ! Default values for namelist variables
    doThinning = .false.
    step    = 6.0d0
    deltmax = 90
    useBlackList = .true.

    ! Read the namelist for Surface observations (if it exists)
    if (utl_isNamelistPresent('thin_surface','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('thn_thinSurface: Error opening file flnml')
      read(nulnam,nml=thin_surface,iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinSurface: Error reading namelist')
      if (mpi_myid == 0) write(*,nml=thin_surface)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'thn_thinSurface: Namelist block thin_surface is missing in the namelist.'
      write(*,*) '                 The default value will be taken.'
      if (mpi_myid == 0) write(*,nml=thin_surface)
    end if

    if (.not. doThinning) return

    call thn_surfaceInTime(obsdat, step, deltmax, useBlackList)

  end subroutine thn_thinSurface

  !--------------------------------------------------------------------------
  ! thn_thinRaobs
  !--------------------------------------------------------------------------
  subroutine thn_thinRaobs(obsdat)
    ! :Purpose: Main thinning subroutine Radiosonde obs.

    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat

    ! Locals:
    integer :: nulnam
    integer :: fnom, fclos, ierr

    ! Namelist variables
    logical :: verticalThinningES !
    logical :: ecmwfRejetsES      !

    namelist /thin_raobs/ verticalThinningES, ecmwfRejetsES

    ! return if no aircraft obs
    if (.not. obs_famExist(obsdat,'UA')) return

    ! Default values for namelist variables
    verticalThinningES = .true.
    ecmwfRejetsES = .true.

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

    call thn_radiosonde(obsdat, verticalThinningES, ecmwfRejetsES)

  end subroutine thn_thinRaobs

  !--------------------------------------------------------------------------
  ! thn_thinAircraft
  !--------------------------------------------------------------------------
  subroutine thn_thinAircraft(obsdat)
    ! :Purpose: Main thinning subroutine for aircraft obs.

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
    ! :Purpose: Main thinning subroutine for satellite winds (AMVs).

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
    integer :: gpsroVarNo    ! bufr element id to be used

    namelist /thin_gpsro/heightMin, heightMax, heightSpacing, gpsroVarNo

    ! return if no gb-gps obs
    if (.not. obs_famExist(obsdat,'RO')) return

    ! Default values for namelist variables
    heightMin     = 1000.0d0
    heightMax     = 40000.0d0
    heightSpacing = 750.0d0
    gpsroVarNo    = BUFR_NERF ! default is refractivity

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

    call thn_gpsroVertical(obsdat, heightMin, heightMax, heightSpacing, gpsroVarNo)

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
    ! :Purpose: Main thinning subroutine for Aladin winds obs.

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
    ! :Purpose: Main thinning subroutine for geostationary radiances (CSR).

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
    ! :Purpose: Main thinning subroutine for scatterometer winds.

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
    ! :Purpose: Main thinning subroutine for AMSU and ATMS obs.

    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat

    ! Locals:
    integer :: nulnam
    integer :: fnom, fclos, ierr

    ! Namelist variables
    integer :: delta    !
    integer :: deltrad

    namelist /thin_tovs/delta, deltrad

    ! return if no TOVS obs
    if (.not. obs_famExist(obsdat,'TO')) return

    ! Default namelist values
    delta   = 100
    deltrad = 75

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
    call thn_tovsFilt(obsdat, delta, deltrad, codtyp_get_codtyp('amsua'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_tovsFilt(obsdat, delta, deltrad, codtyp_get_codtyp('amsub'), &
                      codtyp2_opt=codtyp_get_codtyp('mhs'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_tovsFilt(obsdat, delta, deltrad, codtyp_get_codtyp('atms'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call thn_tovsFilt(obsdat, delta, deltrad, codtyp_get_codtyp('mwhs2'))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine thn_thinTovs

  !--------------------------------------------------------------------------
  ! thn_thinHyper
  !--------------------------------------------------------------------------
  subroutine thn_thinHyper(obsdat)
    ! :Purpose: Main thinning subroutine for hyperspectral infrared radiances.

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
  ! thn_surfaceInTime
  !--------------------------------------------------------------------------
  subroutine thn_surfaceInTime(obsdat, step, deltmax, useBlackList)
    !
    ! :Purpose: Original method for thinning surface data in time.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    real(8),          intent(in)    :: step
    integer,          intent(in)    :: deltmax
    logical,          intent(in)    :: useBlackList

    ! Drifter removal parameters:
    ! Remove incomplete DRIFTER reports (using listEleBadDrifter)?
    logical, parameter :: removeBadDrifters = .true.
    ! Minimum required number elements (1 is consistent with bextrep in ops)
    integer, parameter :: numEleMinBadDrifter = 1
    ! List of required elements for codtyp 18 (see ops derivate program)
    integer, parameter :: listEleBadDrifter(4) = (/ 10051, 11011, 11012, 12004 /)

    ! Selection parameters:
    integer, parameter :: numListCodtyp = 8 ! number of elements in listCodtyp
    ! List of codtyps to keep (what about SYNOP mobil? SA+SYNOP?)
    integer, parameter :: listCodtyp(numListCodtyp) = &
         (/ 12, 146, 13, 147, 18, 143, 144, 15 /)
    character(len=13), parameter :: listCodtypName(numListCodtyp) = &
         (/ 'SYNOP', 'ASYNOP', 'SHIP', 'ASHIP', &
            'DRIFTER', 'SWOB_regular', 'ASWOB_regular', 'METAR' /)
    ! Codtyps to which list_ele_select will be applied
    integer, parameter :: listCodtypSelect(3) = (/ 15, 143, 144 /)
    ! Elements to select (flags for all other elements will have bit 11 set)
    integer, parameter :: listEleSelect(10) = &    ! P, T, Td, U, V
         (/ 8194, 10004, 10051, 11011, 11012, 11215, 11216, 12004, 12006, 12203 /) 

    ! BlackList parameters:
    character(len=*), parameter :: blackListFileName = 'blacklist_sf'
    character(len=6), parameter :: blacklistMode = 'normal' 
    integer, parameter :: numColBlacklist = 5  ! number of columns in blacklist file
    integer, parameter :: numEleBlacklist = 11 ! number of elements in listEleBlacklist
    integer, parameter :: listEleBlacklist(numEleBlacklist) = &
         (/ 10004, 10004, 10051, 10051, 11011, 11012, 11215, 11216, 12004, 12006, 12203 /)
    integer, parameter :: listColBlackList(numEleBlacklist) = &
         (/     1,     2,     1,     2,     5,     5,     5,     5,     3,     4,     4 /)

    ! Locals:
    integer :: countObsIn, countObsInMpi, countObsOut
    integer :: countObsInAllMpi(mpi_nprocs), countObsInMyOffset
    integer :: numElements, codtyp, obsDateStamp, numStep
    integer :: listIndex, obsIndex, obsIndex2, headerIndex, bodyIndex, procIndex
    integer :: ierr, istat, nulfile, numRowBlacklist
    integer :: elemIndex, rowIndex, colIndex, obsFlag, obsVarNo
    integer :: rowBlackList(numColBlacklist) ! blacklist work array
    integer :: countObsInPerCodtyp(numListCodtyp), countObsOutPerCodtyp(numListCodtyp)
    integer :: numEleInPerCodtyp(numListCodtyp), numEleOutPerCodtyp(numListCodtyp)
    integer :: numBit8InPerCodtyp(numListCodtyp), numBit8OutPerCodtyp(numListCodtyp)
    integer :: numBit11InPerCodtyp(numListCodtyp), numBit11OutPerCodtyp(numListCodtyp)
    integer :: numBit8In, numBit8Out, numBit11In, numBit11Out
    integer :: numRemovedCodtyp, numRemovedCodTypPriority
    integer :: numRemovedTime, numRemovedDelt
    integer :: numEleIn, numEleOut, numRepeat1, numRepeat2, numRemovedDrifter
    integer, allocatable :: dataBlacklist(:,:) ! blacklist data matrix
    integer, allocatable :: obsLon(:), obsLat(:), obsDate(:)
    integer, allocatable :: obsLonMpi(:), obsLatMpi(:)
    integer, allocatable :: obsTime(:), obsStepIndex(:), obsStepIndexMpi(:)
    integer, allocatable :: obsCodtypIndex(:), obsDelT(:)
    integer, allocatable :: obsCodtypIndexMpi(:), obsDelTMpi(:)
    real(8) :: obsStepIndex_r8, deltaHours
    logical, allocatable :: valid(:), validMpi(:)
    character(len=9), allocatable :: obsStnid(:), obsStnidMpi(:), stnidBlacklist(:)

    integer, external :: newdate

    ! Check if any observations to be treated
    countObsIn = 0
    call obs_set_current_header_list(obsdat,'SF')
    HEADER0: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER0
      countObsIn = countObsIn + 1
    end do HEADER0

    call rpn_comm_allReduce(countObsIn, countObsInMpi, 1, 'mpi_integer', &
                            'mpi_sum', 'grid', ierr)
    write(*,*)
    if (countObsInMpi == 0) then
      write(*,*) 'thn_surfaceInTime: no surface observations present'
      return
    else
      write(*,*) 'thn_surfaceInTime: countObs initial = ', countObsIn, countObsInMpi
    end if
    write(*,*)

    ! Compute number of time steps in the window for thinning
    numStep = 2*nint(((tim_windowsize - step)/2.d0)/step) + 1
    write(*,*) 'thn_surfaceInTime: step, numStep = ', real(step), numStep
    write(*,*)
 
    ! Print some values to the listing
    write(*,*) 'Codtyps to which selection will be applied:'
    do listIndex = 1, size(listCodtypSelect)
      write(*,*) listCodtypSelect(listIndex)
    end do
    write(*,*) 'Elements to select from above codtyps:'
    do listIndex = 1, size(listEleSelect)
      write(*,*) listEleSelect(listIndex)
    end do
    write(*,*) 'Remove incomplete DRIFTER reports: ', removeBadDrifters
    write(*,*)

    ! Allocate arrays
    allocate(valid(countObsIn))
    allocate(obsCodtypIndex(countObsIn))
    allocate(obsLon(countObsIn))
    allocate(obsLat(countObsIn))
    allocate(obsDate(countObsIn))
    allocate(obsTime(countObsIn))
    allocate(obsDelT(countObsIn))
    allocate(obsStepIndex(countObsIn))
    allocate(obsStnid(countObsIn))

    ! Initialize dynamic arrays
    valid(:) =    .true.    ! array which keeps track of which reports to keep
    obsCodtypIndex(:) = 999 ! codtyp index array (999 means non-existent)
    obsLon(:)           = 0 ! longitude corresponding to each report
    obsLat(:)           = 0 ! latitude corresponding to each report
    obsDate(:)          = 0 ! DATE yyyymmdd corresponding to each report
    obsTime(:)          = 0 ! time hhmm corresponding to each report
    obsDelT(:)       = 9999 ! delta t array (departure from nearest bin time)
    obsStepIndex(:)     = 0 ! temporal bin corresponding to each report
    obsStnid(:)        = '' ! stnid corresponding to each report

    ! Initialize counters and counter arrays
    numEleIn                 = 0 ! number of elements in input file
    numBit8In                = 0 ! number of input elements flagged (blacklist)
    numBit11In               = 0 ! number of input elements flagged (selection)
    numEleOut                = 0 ! number of elements in output file
    numBit8Out               = 0 ! number of output elements flagged (blacklist)
    numBit11Out              = 0 ! number of output elements flagged (selection)
    numRemovedCodtyp         = 0 ! number of obs removed due to network
    numRemovedTime           = 0 ! number of obs removed (desired time window)
    numRemovedCodTypPriority = 0 ! number of obs removed due to codtyp
    numRemovedDelt           = 0 ! number of obs removed due to delta t
    numRemovedDrifter        = 0 ! number of incomplete drifter reports removed
    numRepeat1               = 0 ! number of obs same lon/lat/date/time
    numRepeat2               = 0 ! number of obs same lon/lat/date/time/codtyp

    countObsInPerCodtyp(:)  = 0 ! number of input reports for each codtyp
    numEleInPerCodtyp(:)    = 0 ! number of input elements for each codtyp
    numBit8InPerCodtyp(:)   = 0 ! number of input flags bit 8 for each codtyp
    numBit11InPerCodtyp(:)  = 0 ! number of input flags bit 11 for each codtyp
    countObsOutPerCodtyp(:) = 0 ! number of output reports for each codtyp
    numEleOutPerCodtyp(:)   = 0 ! number of output elements for each codtyp
    numBit8OutPerCodtyp(:)  = 0 ! number of output flags bit 8 for each codtyp
    numBit11OutPerCodtyp(:) = 0 ! number of output flags bit 11 for each codtyp

    ! Extract needed information from obsSpaceData
    obsIndex = 0
    call obs_set_current_header_list(obsdat,'SF')
    HEADER1: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER1
      obsIndex = obsIndex + 1

      ! Get index in codtyp list
      codtyp = obs_headElem_i(obsdat, OBS_ITY, headerIndex)
      do listIndex = 1, numListCodtyp
        if (codtyp == listCodtyp(listIndex)) then
          obsCodtypIndex(obsIndex) = listIndex
        end if
      end do

      ! Remove observation if not in codtyp list
      if (obsCodtypIndex(obsIndex) == 999) then
        valid(obsIndex) = .false.
        numRemovedCodtyp = numRemovedCodtyp + 1
        cycle HEADER1
      end if

      ! Check if DRIFTER reports contain required element(s), otherwise cycle
      if ((codtyp == 18) .and. removeBadDrifters) then
        numElements = 0
        call obs_set_current_body_list(obsdat, headerIndex)
        BODY1: do 
          bodyIndex = obs_getBodyIndex(obsdat)
          if (bodyIndex < 0) exit BODY1
          obsVarNo = obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex)
          if (any(listEleBadDrifter(:) == obsVarNo)) numElements = numElements + 1
        enddo BODY1
        if (numElements < numEleMinBadDrifter) then
          valid(obsIndex) = .false.
          numRemovedDrifter = numRemovedDrifter + 1
          cycle HEADER1
        end if
      end if

      ! Count the number of elements
      numElements = 0
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY2: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY2

        obsVarNo = obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex)
        if (obsVarNo == -1) cycle BODY2

        numElements = numElements + 1

        obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)

        ! Count input flags with bit 8 set
        if (btest(obsFlag, 8)) then
          numBit8InPerCodtyp(obsCodtypIndex(obsIndex)) = &
               numBit8InPerCodtyp(obsCodtypIndex(obsIndex)) + 1
          numBit8In = numBit8In + 1
        end if

        ! Count input flags with bit 11 set
        if (btest(obsFlag, 11)) then
          numBit11InPerCodtyp(obsCodtypIndex(obsIndex)) = &
               numBit11InPerCodtyp(obsCodtypIndex(obsIndex)) + 1
          numBit11In = numBit11In + 1
        end if

      end do BODY2

      numEleIn = numEleIn + numElements

      ! Counts per codtyp
      listIndex = obsCodtypIndex(obsIndex)
      countObsInPerCodtyp(listIndex) = countObsInPerCodtyp(listIndex) + 1
      numEleInPerCodtyp(listIndex) = numEleInPerCodtyp(listIndex) + numElements

      obsStnid(obsIndex) = obs_elem_c(obsdat,'STID',headerIndex)

      obsLon(obsIndex) = nint(100.0 * MPC_DEGREES_PER_RADIAN_R8 * &
                              obs_headElem_r(obsdat, OBS_LON, headerIndex))
      obsLat(obsIndex) = nint(100.0 * MPC_DEGREES_PER_RADIAN_R8 * &
                              obs_headElem_r(obsdat, OBS_LAT, headerIndex)) + 9000

      obsDate(obsIndex) = obs_headElem_i(obsdat, OBS_DAT, headerIndex)
      obsTime(obsIndex) = obs_headElem_i(obsdat, OBS_ETM, headerIndex)

      call tim_getStepObsIndex(obsStepIndex_r8, tim_getDatestamp(), &
                               obsDate(obsIndex), obsTime(obsIndex), numStep)
      obsStepIndex(obsIndex) = nint(obsStepIndex_r8)
      if (numStep > 1) then
        obsDelT(obsIndex) = &
             nint( 60.0 * step * &
                   (obsStepIndex_r8 - real(obsStepIndex(obsIndex))) )
      else
        ierr = newdate(obsDateStamp, obsDate(obsIndex), obsTime(obsIndex)*10000, 3)
        call difdatr(obsDateStamp,tim_getDateStamp(),deltaHours)
        obsDelT(obsIndex) = nint(60.0*deltaHours)
      end if

      ! Reject if time difference larger than deltmax
      if (abs(obsDelT(obsIndex)) > deltmax) then
        valid(obsIndex) = .false.
        numRemovedTime = numRemovedTime + 1
      endif

      write(*,*) 'DEBUG: ', obsLon(obsIndex), obsLat(obsIndex), obsDate(obsIndex), obsTime(obsIndex), obsStnid(obsIndex)
    end do HEADER1

    call utl_allReduce(numRemovedDrifter)
    call utl_allReduce(numRemovedTime)
    call utl_allReduce(numRemovedCodtyp)

    ! Read blacklist file
    if (useBlackList) then
      write(*,*) 'Opening blacklist file'
      numRowBlacklist = 0
      nulfile = 0
      open (unit=nulfile, file=blackListFileName, status='OLD', iostat=ierr)
      if (ierr /= 0) then
        write(*,*) 'Cannot open blacklist file ', trim(blackListFileName)
        call utl_abort('thn_surfaceInTime')
      end if
      read(nulfile, iostat=istat, fmt='(i6)') numRowBlacklist
      write(*,*) 'thn_surfaceInTime: Number of stations in blacklist: ', numRowBlacklist

      allocate(stnidBlacklist(numRowBlacklist))
      allocate(dataBlacklist(numRowBlacklist, numEleBlacklist))
      stnidBlacklist(:) = '' ! array of stnid values in blacklist file
      dataBlacklist(:,:)    = 0 ! blacklist matrix for stnids and elements

      do rowIndex = 1, numRowBlacklist
        read(nulfile, iostat=istat, fmt='(x,a8,x,5(x,i1))') &
             stnidBlacklist(rowIndex), &
             (rowBlackList(colIndex), colIndex = 1, numColBlacklist)
        do elemIndex = 1, numEleBlacklist
          if ( (blacklistMode == 'severe') .or. &
               (rowBlackList(listColBlackList(elemIndex))==1)) then
            dataBlacklist(rowIndex,elemIndex) = 1
          end if
        end do
      end do
      write(*,*) 'thn_surfaceInTime: Closing blacklist file'
      write(*,*)
      close (unit=nulfile)
    end if

    ! Gather array information over all mpi tasks
    allocate(validMpi(countObsInMpi))
    allocate(obsCodtypIndexMpi(countObsInMpi))
    allocate(obsLonMpi(countObsInMpi))
    allocate(obsLatMpi(countObsInMpi))
    allocate(obsDelTMpi(countObsInMpi))
    allocate(obsStepIndexMpi(countObsInMpi))
    allocate(obsStnidMpi(countObsInMpi))

    call intArrayToMpi(obsLon, obsLonMpi)
    call intArrayToMpi(obsLat, obsLatMpi)
    call intArrayToMpi(obsCodtypIndex, obsCodtypIndexMpi)
    call intArrayToMpi(obsStepIndex, obsStepIndexMpi)
    call intArrayToMpi(obsDelT, obsDelTMpi)
    call logicalArrayToMpi(valid, validMpi)
    call stringArrayToMpi(obsStnid, obsStnidMpi)
    
    ! Apply the thinning algorithm
    do obsIndex = 1, countObsInMpi
      
      ! If current report OK so far
      if (validMpi(obsIndex)) then
        ! Loop over previously-read reports
        do obsIndex2 = 1, obsIndex - 1
          ! If other report OK so far and both reports in same bin
          if ( validMpi(obsIndex2) .and. &
               (obsStepIndexMpi(obsIndex) == obsStepIndexMpi(obsIndex2)) ) then
            ! If reports are spatially colocated or have same stnid
            if ( ( (obsLonMpi(obsIndex) == obsLonMpi(obsIndex2)) .and. &
                   (obsLatMpi(obsIndex) == obsLatMpi(obsIndex2)) ) .or. &
                 (obsStnidMpi(obsIndex) == obsStnidMpi(obsIndex2)) ) then
              ! If both reports have same codtyp
              if (obsCodtypIndexMpi(obsIndex) == obsCodtypIndexMpi(obsIndex2)) then
                ! If current report closer to bin time
                if (abs(obsDelTMpi(obsIndex)) < abs(obsDelTMpi(obsIndex2))) then
                  ! Reject other report
                  validMpi(obsIndex2) = .false.
                  numRemovedDelt = numRemovedDelt + 1
                else if (abs(obsDelTMpi(obsIndex)) > abs(obsDelTMpi(obsIndex2))) then
                  ! other report closer to bin time, so reject current report
                  validMpi(obsIndex) = .false.
                  numRemovedDelt = numRemovedDelt + 1
                else if (obsDelTMpi(obsIndex) >= obsDelTMpi(obsIndex2)) then
                  ! both reports equally far from bin time
                  ! If other report not more recent then reject it
                  validMpi(obsIndex2) = .false.
                  if (obsDelTMpi(obsIndex) == obsDelTMpi(obsIndex2)) then
                    numRepeat2 = numRepeat2 + 1
                  else
                    numRemovedDelt = numRemovedDelt + 1
                  end if
                else
                  ! current report less recent, reject it
                  validMpi(obsIndex) = .false.
                  numRemovedDelt = numRemovedDelt + 1
                end if ! delta t
              else
                ! Reports do not have same codtyp
                ! If current report has higher codtyp precedence
                if (obsCodtypIndexMpi(obsIndex) < obsCodtypIndexMpi(obsIndex2)) then
                  ! Reject other report
                  validMpi(obsIndex2) = .false.
                else
                  ! Other report has higher codtyp precedence, so reject current report
                  validMpi(obsIndex) = .false.
                end if
                numRemovedCodTypPriority = numRemovedCodTypPriority + 1
                if (obsDelTMpi(obsIndex) == obsDelTMpi(obsIndex2)) then
                  numRepeat1 = numRepeat1 + 1
                end if
              end if ! codtyp
            end if ! lon, lat, ibin
          end if ! other report OK
        end do ! obsIndex2
      end if ! current report OK

    end do ! obsIndex

    ! Transfer mpi global array 'valid' to local array
    call rpn_comm_allGather(countObsIn,       1, 'mpi_integer',  &
                            countObsInAllMpi, 1, 'mpi_integer', &
                            'GRID', ierr )
    countObsInMyOffset = 0
    do procIndex = 1, mpi_myid
      countObsInMyOffset = countObsInMyOffset + countObsInAllMpi(procIndex)
    end do
    do obsIndex = 1, countObsIn
      valid(obsIndex) = validMpi(obsIndex + countObsInMyOffset)
    end do

    ! Do counts of kepts observations
    obsIndex = 0
    countObsOut = 0
    call obs_set_current_header_list(obsdat,'SF')
    HEADER2: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER2
      obsIndex = obsIndex + 1

      ! If rejected, set bit 11 for all data flags
      if (.not. valid(obsIndex)) then
        call obs_set_current_body_list(obsdat, headerIndex)
        BODY3: do
          bodyIndex = obs_getBodyIndex(obsdat)
          if (bodyIndex < 0) exit BODY3
          obsFlag  = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
          call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,11))
        end do BODY3
        
        cycle HEADER2
      end if

      codtyp   = obs_headElem_i(obsdat, OBS_ITY, headerIndex)

      ! Count the number of elements
      numElements = 0
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY4: do
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY4

        obsVarNo = obs_bodyElem_i(obsdat, OBS_VNM, bodyIndex)
        if (obsVarNo == -1) cycle BODY4

        numElements = numElements + 1

        ! Set bit 8 according to blacklist, if blacklist present
        if (useBlackList) then
          ! Traverse rows of blacklist matrix
          do rowIndex = 1, numRowBlacklist
            ! If stnid found in blacklist, flag elements as appropriate
            if (obsStnid(obsIndex) == stnidBlacklist(rowIndex)) then
              ! Traverse columns of blacklist matrix
              do elemIndex = 1, numEleBlacklist
                ! If element is to be blacklisted, set bit 8
                if ( obsVarNo == listEleBlacklist(elemIndex) .and. &
                     dataBlacklist(rowIndex,elemIndex) == 1 ) then
                  obsFlag  = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
                  call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,8))
                end if
              end do ! elemIndex
            end if ! stnid
          end do ! rowIndex
        end if ! useBlackList

        ! Set bit 11 according to requested codtyps and elements
        if (any(listCodtypSelect(:) == codtyp)) then
          ! If current element not in select list, set bit 11
          if (.not. any(listEleSelect(:) == obsVarNo)) then
            obsFlag  = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
write(*,*) 'Setting bit 11 for codtyp, elem = ', codtyp, obsVarNo
            call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,11))
          end if
        end if

        obsFlag  = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)

        ! Count output flags with bit 8 set
        if (btest(obsFlag, 8)) then
          numBit8OutPerCodtyp(obsCodtypIndex(obsIndex)) = &
               numBit8OutPerCodtyp(obsCodtypIndex(obsIndex)) + 1
          numBit8Out = numBit8Out + 1
        end if

        ! Count output flags with bit 11 set
        if (btest(obsFlag, 11)) then
          numBit11OutPerCodtyp(obsCodtypIndex(obsIndex)) = &
               numBit11OutPerCodtyp(obsCodtypIndex(obsIndex)) + 1
          numBit11Out = numBit11Out + 1
        end if

      end do BODY4

      numEleOut = numEleOut + numElements
      countObsOut = countObsOut + 1

      ! Counts per codtyp
      listIndex = obsCodtypIndex(obsIndex)
      countObsOutPerCodtyp(listIndex) = countObsOutPerCodtyp(listIndex) + 1
      numEleOutPerCodtyp(listIndex) = numEleOutPerCodtyp(listIndex) + numElements

    end do HEADER2

    ! Output statistics to screen

    ! numRepeat1 should include the case where codtyps are the same
    numRepeat1 = numRepeat1 + numRepeat2

    write(*,'(a)') ' Number of reports in input file'
    do listIndex = 1, numListCodtyp
      call utl_allReduce(countObsInPerCodtyp(listIndex))
      write(*,'(i4,3a,i7)') listCodtyp(listIndex), ' (', &
           listCodtypName(listIndex), '): ', countObsInPerCodtyp(listIndex)
    end do

    write(*,*)
    write(*,'(a)') ' Number of elements in input file'
    do listIndex = 1, numListCodtyp
      call utl_allReduce(numEleInPerCodtyp(listIndex))
      call utl_allReduce(numBit8InPerCodtyp(listIndex))
      call utl_allReduce(numBit11InPerCodtyp(listIndex))
      write(*,'(i4,3a,i7,a,i7,a,i7,a)') listCodtyp(listIndex), ' (', &
           listCodtypName(listIndex), '): ', numEleInPerCodtyp(listIndex), ' (', &
           numBit8InPerCodtyp(listIndex), ' bit 8, ', &
           numBit11InPerCodtyp(listIndex), ' bit 11)'
    end do

    write(*,*)
    write(*,'(a)') ' Number of reports in output file'
    do listIndex = 1, numListCodtyp
      call utl_allReduce(countObsOutPerCodtyp(listIndex))
      write(*,'(i4,3a,i7)') listCodtyp(listIndex), ' (', &
           listCodtypName(listIndex), '): ', countObsOutPerCodtyp(listIndex)
    end do

    write(*,*)
    write(*,'(a)') ' Number of elements in output file'
    do listIndex = 1, numListCodtyp
      call utl_allReduce(numEleOutPerCodtyp(listIndex))
      call utl_allReduce(numBit8OutPerCodtyp(listIndex))
      call utl_allReduce(numBit11OutPerCodtyp(listIndex))
      write(*,'(i4,3a,i7,a,i7,a,i7,a)') listCodtyp(listIndex), ' (', &
           listCodtypName(listIndex), '): ', numEleOutPerCodtyp(listIndex), ' (', &
           numBit8OutPerCodtyp(listIndex), ' bit 8, ', &
           numBit11OutPerCodtyp(listIndex), ' bit 11)'
    end do

    write(*,*)
    call utl_allReduce(countObsIn)
    call utl_allReduce(countObsOut)
    write(*,'(a,i7)') 'Total number of reports in input file:   ', countObsIn
    write(*,'(a,i7)') 'Total number of reports in output file:  ', countObsOut
    call utl_allReduce(numEleIn)
    call utl_allReduce(numBit8In)
    call utl_allReduce(numBit11In)
    call utl_allReduce(numEleOut)
    call utl_allReduce(numBit8Out)
    call utl_allReduce(numBit11Out)
    write(*,'(a,i7,a,i7,a,i7,a)') 'Total number of elements in input file:  ', numEleIn, &
          ' (', numBit8In, ' bit 8, ', numBit11In, ' bit 11)'
    write(*,'(a,i7,a,i7,a,i7,a)') 'Total number of elements in output file: ', numEleOut, &
          ' (', numBit8Out, ' bit 8, ', numBit11Out, ' bit 11)'
    write(*,*)
    write(*,'(a,i7)') 'Number of repeated reports (lon/lat/date/time):        ', numRepeat1
    write(*,'(a,f6.2)') 'Above count as a percentage of total reports in:        ', 100.0 * numRepeat1 / countObsIn
    write(*,'(a,i7)') 'Number of repeated reports (lon/lat/date/time/codtyp): ', numRepeat2
    write(*,'(a,f6.2)') 'Above count as a percentage of total reports in:        ', 100.00 * numRepeat2 / countObsIn
    write(*,*)
    write(*,'(a,i7)') 'Number of reports removed due to codtyp:               ', numRemovedCodtyp
    write(*,'(a,i7)') 'Number of reports removed due to time:                 ', numRemovedTime
    write(*,'(a,i7)') 'Number of reports removed using codtyp precedence:     ', numRemovedCodTypPriority
    write(*,'(a,i7)') 'Number of reports removed using delta t:               ', numRemovedDelt
    write(*,'(a,i7)') 'Number of incomplete drifter reports removed:          ', numRemovedDrifter
    write(*,*)
 
  end subroutine thn_surfaceInTime
  
  !--------------------------------------------------------------------------
  ! thn_gpsroVertical
  !--------------------------------------------------------------------------
  subroutine thn_gpsroVertical(obsdat, heightMin, heightMax, heightSpacing, &
                               gpsroVarNo)
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
    integer,          intent(in)    :: gpsroVarNo

    ! Locals:
    integer :: countObs, countObsMpi
    integer :: countObsReject, countObsRejectMpi, countObsTotal, countObsTotalMpi
    integer :: headerIndex, bodyIndex, ierr
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

    call rpn_comm_allReduce(countObs, countObsMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    if (countObsMpi == 0) then
      write(*,*) 'thn_gpsroVertical: no gpsro observations present'
      return
    end if

    countObsTotal  = 0
    countObsReject = 0
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
        countObsTotal = countObsTotal + 1

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
          countObsReject = countObsReject + 1
        end if
        
      end do LEVELS

      deallocate(obsHeights)
      deallocate(bodyIndexList)

    end do HEADER1

    call rpn_comm_allReduce(countObsTotal, countObsTotalMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(countObsReject, countObsRejectMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    write(*,*)' Number of GPS-RO elements, total     --->', countObsTotalMpi
    write(*,*)' Number of GPS-RO elements, rejected  --->', countObsRejectMpi
    write(*,*)' Number of GPS-RO elements, kept      --->', countObsTotalMpi - &
                                                            countObsRejectMpi

    write(*,*)
    write(*,*) 'thn_gpsroVertical: Finished'
    write(*,*)

  end subroutine thn_gpsroVertical

  !--------------------------------------------------------------------------
  ! thn_radiosonde
  !--------------------------------------------------------------------------
  subroutine thn_radiosonde(obsdat, verticalThinningES, ecmwfRejetsES)
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

    ! Locals:
    type(struct_hco), pointer :: hco_sfc
    type(struct_vco), pointer :: vco_sfc
    type(struct_gsv)          :: stateVectorPsfc
    integer :: fnom, fclos, ezgdef, ezsint, ezdefset, ezsetopt
    integer :: ierr, nulnam, numLevStn, numLevStnMpi, countLevel, numLevStnMax
    integer :: numStation, numStationMpi, stationIndex, stationIndexMpi, lastProfileIndex
    integer :: profileIndex, headerIndex, bodyIndex, levIndex, stepIndex, varIndex
    integer :: levStnIndex, levStnIndexMpi, obsFlag
    integer :: ig1obs, ig2obs, ig3obs, ig4obs, obsGridID
    real(4) :: obsValue, obsOmp, obsStepIndex
    real(4) :: zig1, zig2, zig3, zig4, zpresa, zpresb
    real(8) :: obsStepIndex_r8
    integer, allocatable :: obsLevOffset(:), obsType(:), obsHeadDate(:), obsHeadTime(:)
    integer, allocatable :: obsDate(:), obsTime(:), obsLaunchTime(:), stationFlags(:)
    integer, allocatable :: trajFlags(:,:), obsFlags(:,:)
    integer, allocatable :: obsLevOffsetMpi(:), obsHeadDateMpi(:)
    integer, allocatable :: obsLaunchTimeMpi(:), stationFlagsMpi(:)
    integer, allocatable :: trajFlagsMpi(:,:), obsFlagsMpi(:,:)
    integer, allocatable :: procIndexes(:), procIndexesMpi(:)
    real(4), allocatable :: obsLat(:), obsLon(:), surfPresInterp(:,:)
    real(4), allocatable :: obsValues(:,:), oMinusB(:,:), presInterp(:,:)
    real(4), allocatable :: obsLatMpi(:), obsLonMpi(:)
    real(4), allocatable :: obsValuesMpi(:,:), oMinusBMpi(:,:)
    real(4), pointer     :: surfPressure(:,:,:,:)
    character(len=9), allocatable :: stnId(:), stnIdMpi(:)
    character(len=20) :: trlmFileName
    character(len=2)  :: fileNumber
    logical :: upperAirObs
    integer :: countAcc_dd, countRej_dd, countAccMpi_dd, countRejMpi_dd
    integer :: countAcc_ff, countRej_ff, countAccMpi_ff, countRejMpi_ff
    integer :: countAcc_tt, countRej_tt, countAccMpi_tt, countRejMpi_tt
    integer :: countAcc_es, countRej_es, countAccMpi_es, countRejMpi_es

    ! Local parameters:
    integer, parameter :: numVars=5, numTraj=2, maxNumLev=300

    ! Namelist variables:
    real(8) :: rprefinc
    real(8) :: rptopinc
    real(8) :: rcoefinc
    real(4) :: vlev(maxNumLev)
    integer :: numlev
    namelist /namgem/rprefinc, rptopinc, rcoefinc, numlev, vlev

    ! Check if any observations to be treated and count number of "profiles"
    numLevStn = 0
    numStation = 0
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

      numLevStn = numLevStn + 1

      profileIndex = obs_headElem_i(obsdat, obs_prfl, headerIndex)
      if (profileIndex /= lastProfileIndex) then
        lastProfileIndex = profileIndex
        numStation = numStation + 1
      end if

    end do HEADER0

    call rpn_comm_allReduce(numLevStn, numLevStnMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    if (numLevStnMpi == 0) then
      write(*,*) 'thn_radiosonde: no UA obs observations present'
      return
    end if

    call rpn_comm_allReduce(numStation, numStationMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)

    write(*,*) 'thn_radiosonde: number of obs initial = ', &
               numLevStn, numLevStnMpi
    write(*,*) 'thn_radiosonde: number of profiles    = ', &
               numStation, numStationMpi

    ! Allocate some quanitities needed for each profile
    allocate(obsLevOffset(numStation+1))
    allocate(obsType(numStation))
    allocate(obsHeadDate(numStation))
    allocate(obsHeadTime(numStation))
    allocate(obsDate(numStation))
    allocate(obsTime(numStation))
    allocate(obsLaunchTime(numStation))
    allocate(stationFlags(numStation))
    allocate(obsLat(numStation))
    allocate(obsLon(numStation))
    allocate(stnId(numStation))
    allocate(presInterp(numStation,maxNumLev))
    allocate(trajFlags(numTraj,numLevStn))
    allocate(obsFlags(numVars,numLevStn))
    allocate(obsValues(numVars,numLevStn))
    allocate(oMinusB(numVars,numLevStn))
    obsLevOffset(:) = 0
    obsDate(:) = -1
    trajFlags(:,:) = -1
    obsFlags(:,:) = 0
    obsValues(:,:) = -999.0
    oMinusB(:,:) = -999.0

    ! Fill in the arrays for each profile
    levStnIndex = 0
    stationIndex = 0
    lastProfileIndex = -1
    numLevStnMax = -1
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
        stationIndex = stationIndex + 1
        countLevel = 0
      end if

      countLevel = countLevel + 1
      obsLevOffset(stationIndex+1) = obsLevOffset(stationIndex) + countLevel

      if (countLevel > numLevStnMax) numLevStnMax = countLevel

      ! Get some information from the first header in this profile
      if (countLevel == 1) then
        obsDate(stationIndex)  = obs_headElem_i(obsdat,obs_dat,headerIndex)
        obsTime(stationIndex)  = obs_headElem_i(obsdat,obs_etm,headerIndex)
        obsLat(stationIndex)   = obs_headElem_r(obsdat,obs_lat,headerIndex) * &
                                 MPC_DEGREES_PER_RADIAN_R8
        obsLon(stationIndex)   = obs_headElem_r(obsdat,obs_lon,headerIndex) * &
                                 MPC_DEGREES_PER_RADIAN_R8
        obsLat(stationIndex)   = 0.01*nint(100.0*obsLat(stationIndex))
        obsLon(stationIndex)   = 0.01*nint(100.0*obsLon(stationIndex))
      end if
      obsHeadDate(stationIndex)   = obs_headElem_i(obsdat,obs_hdd,headerIndex)
      obsHeadTime(stationIndex)   = obs_headElem_i(obsdat,obs_hdt,headerIndex)
      obsLaunchTime(stationIndex) = obs_headElem_i(obsdat,obs_lch,headerIndex)
      if (obsLaunchTime(stationIndex) == mpc_missingValue_int) then
        obsLaunchTime(stationIndex) = obsHeadTime(stationIndex)
      end if
      obsType(stationIndex)      = obs_headElem_i(obsdat,obs_rtp,headerIndex)
      stationFlags(stationIndex) = obs_headElem_i(obsdat,obs_st1,headerIndex)
      stnId(stationIndex)        = obs_elem_c(obsdat,'STID',headerIndex)

      trajFlags(1,levStnIndex) = obs_headElem_i(obsdat,obs_tflg,headerIndex)
      trajFlags(2,levStnIndex) = obs_headElem_i(obsdat,obs_lflg,headerIndex)

      call obs_set_current_body_list(obsdat, headerIndex)
      BODY2: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY2

        obsFlag  = obs_bodyElem_i(obsdat,obs_flg,bodyIndex)
        obsValue = obs_bodyElem_r(obsdat,obs_var,bodyIndex)
        obsOmp   = obs_bodyElem_r(obsdat,obs_omp,bodyIndex)
        select case (obs_bodyElem_i(obsdat,OBS_VNM,bodyIndex))
        case (bufr_nedd)
          obsFlags(1,levStnIndex)  = obsFlag
          obsValues(1,levStnIndex) = obsValue
        case (bufr_neuu)
          oMinusB(1,levStnIndex)   = obsOmp
        case (bufr_neff)
          obsFlags(2,levStnIndex)  = obsFlag
          obsValues(2,levStnIndex) = obsValue
        case (bufr_nevv)
          oMinusB(2,levStnIndex)   = obsOmp
        case (bufr_nett)
          obsFlags(3,levStnIndex)  = obsFlag
          obsValues(3,levStnIndex) = obsValue
          oMinusB(3,levStnIndex)   = obsOmp
        case (bufr_nees)
          obsFlags(4,levStnIndex)  = obsFlag
          obsValues(4,levStnIndex) = obsValue
          oMinusB(4,levStnIndex)   = obsOmp
        end select
        obsValues(5,levStnIndex) = 0.01 * obs_bodyElem_r(obsdat,obs_ppp,bodyIndex)
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
      call gio_readFromFile( stateVectorPsfc, trlmFileName, ' ', ' ',  &
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
    obsGridID = ezgdef(numStation,1,'Y','L',ig1obs,ig2obs,ig3obs,ig4obs,obsLon,obsLat)
    ierr = ezdefset(obsGridID,hco_sfc%ezScintId)

    ! Do the interpolation of surface pressure to obs locations
    allocate(surfPresInterp(numStation,tim_nstepobs))
    do stepIndex = 1, tim_nstepobs
      ierr = ezsint(surfPresInterp(:,stepIndex),surfPressure(:,:,1,stepIndex))
    end do
    call gsv_deallocate(stateVectorPsfc)

    ! Compute pressure profile at each obs location
    do stationIndex = 1, numStation

      ! Calculate the stepIndex corresponding to the launch time
      call tim_getStepObsIndex(obsStepIndex_r8, tim_getDatestamp(), &
                               obsDate(stationIndex), obsTime(stationIndex), tim_nstepobs)
      obsStepIndex = nint(obsStepIndex_r8)
      if (obsStepIndex < 0) then
        obsStepIndex = (tim_nstepobs+1)/2
        write(*,*) 'thn_radiosonde: Obs outside the assimilation window, set to middle of window'
      end if

      ! Calculate pressure levels for each station based on GEM3 vertical coordinate
      do levIndex  = 1, numLev
        zpresb = ((vlev(levIndex) - rptopinc/rprefinc) /  &
                 (1.0-rptopinc/rprefinc))**rcoefinc
        zpresa = rprefinc * (vlev(levIndex)-zpresb)
        presInterp(stationIndex,levIndex) =  &
             0.01 * (zpresa + zpresb*surfPresInterp(stationIndex,obsStepIndex))
      end do

    end do

    ! communicate several arrays to all MPI tasks for raobs_check_duplicated_stations
    allocate(obsLevOffsetMpi(numStationMpi+1))
    allocate(obsHeadDateMpi(numStationMpi))
    allocate(obsLaunchTimeMpi(numStationMpi))
    allocate(stationFlagsMpi(numStationMpi))
    allocate(obsLatMpi(numStationMpi))
    allocate(obsLonMpi(numStationMpi))
    allocate(trajFlagsMpi(numTraj,numLevStnMpi))
    allocate(obsFlagsMpi(numVars,numLevStnMpi))
    allocate(obsValuesMpi(numVars,numLevStnMpi))
    allocate(oMinusBMpi(numVars,numLevStnMpi))
    allocate(stnIdMpi(numStationMpi))

    call intArrayToMpi(obsLevOffset, obsLevOffsetMpi, is_obsLevOffset_opt=.true.)
    call intArrayToMpi(obsHeadDate, obsHeadDateMpi)
    call intArrayToMpi(obsLaunchTime, obsLaunchTimeMpi)
    call intArrayToMpi(stationFlags, stationFlagsMpi)
    call realArrayToMpi(obsLat, obsLatMpi)
    call realArrayToMpi(obsLon, obsLonMpi)
    do varIndex = 1, numTraj
      call intArrayToMpi(trajFlags(varIndex,:), trajFlagsMpi(varIndex,:))
    end do
    do varIndex = 1, numVars
      call intArrayToMpi(obsFlags(varIndex,:), obsFlagsMpi(varIndex,:))
      call realArrayToMpi(obsValues(varIndex,:), obsValuesMpi(varIndex,:))
      call realArrayToMpi(oMinusB(varIndex,:), oMinusBMpi(varIndex,:))
    end do
    call stringArrayToMpi(stnId,stnIdMpi)

    ! set stnIdMpi to 'NOT_VALID' for rejected duplicate stations
    call raobs_check_duplicated_stations( stnIdMpi, obsLevOffsetMpi, obsLatMpi, obsLonMpi, &
                                          obsHeadDateMpi, trajFlagsMpi, &
                                          obsFlagsMpi, obsValuesMpi, oMinusBMpi, &
                                          obsLaunchTimeMpi, stationFlagsMpi, &
                                          numVars, numStationMpi )

    ! modify obs flags based on stnIdMpi
    do stationIndexMpi = 1, numStationMpi
      if (stnIdMpi(stationIndexMpi) == 'NOT_VALID') then
        do levStnIndex = obsLevOffsetMpi(stationIndexMpi)+1, &
                         obsLevOffsetMpi(stationIndexMpi+1)
          do varIndex = 1, numVars
            obsFlagsMpi(varIndex,levStnIndex) =  &
                 ibset(obsFlagsMpi(varIndex,levStnIndex),11)
          end do
        end do
      end if
    end do

    ! copy global mpi flags to local copy
    allocate(procIndexes(numStation))
    allocate(procIndexesMpi(numStationMpi))
    procIndexes(:) = mpi_myid
    call intArrayToMpi(procIndexes, procIndexesMpi)
    levStnIndex = 0
    do stationIndexMpi = 1, numStationMpi
      if (procIndexesMpi(stationIndexMpi) == mpi_myid) then
        do levStnIndexMpi = obsLevOffsetMpi(stationIndexMpi)+1, &
                            obsLevOffsetMpi(stationIndexMpi+1)
          levStnIndex = levStnIndex + 1
          obsFlags(:,levStnIndex) = obsFlagsMpi(:,levStnIndexMpi)
        end do
      end if
    end do
    deallocate(procIndexes)
    deallocate(procIndexesMpi)

    call raobs_thinning_model( obsFlags, obsValues, presInterp, numVars, numlev, &
                               numStation, numLevStnMax, obsLevOffset )

    if ( verticalThinningES ) then
      call raobs_thinning_es( obsFlags, obsValues, numStation, numLevStnMax, obsLevOffset )
    end if

    if ( ecmwfRejetsES ) then
      call raobs_blacklisting_ecmwf( obsFlags, obsValues, obsType, numStation, obsLevOffset )
    end if

    countAcc_dd=0;  countRej_dd=0
    countAcc_ff=0;  countRej_ff=0
    countAcc_tt=0;  countRej_tt=0
    countAcc_es=0;  countRej_es=0
    do stationIndex = 1, numStation
      do levStnIndex = obsLevOffset(stationIndex)+1, obsLevOffset(stationIndex+1)
        if (btest(obsFlags(1,levStnIndex),11)) then
          countRej_dd = countRej_dd + 1
        else
          countAcc_dd = countAcc_dd + 1
        end if
        if (btest(obsFlags(2,levStnIndex),11)) then
          countRej_ff = countRej_ff + 1
        else
          countAcc_ff = countAcc_ff + 1
        end if
        if (btest(obsFlags(3,levStnIndex),11)) then
          countRej_tt = countRej_tt + 1
        else
          countAcc_tt = countAcc_tt + 1
        end if
        if (btest(obsFlags(4,levStnIndex),11) .or. btest(obsFlags(4,levStnIndex),8)) then
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

    write(*,*)
    write(*,*) 'DD Rej/Acc = ',countRejMpi_dd, countAccMpi_dd
    write(*,*) 'FF Rej/Acc = ',countRejMpi_ff, countAccMpi_ff
    write(*,*) 'TT Rej/Acc = ',countRejMpi_tt, countAccMpi_tt
    write(*,*) 'ES Rej/Acc = ',countRejMpi_es, countAccMpi_es
    write(*,*)

    ! Replace obs flags in obsSpaceData with obsFlags
    levStnIndex = 0
    call obs_set_current_header_list(obsdat,'UA')
    HEADER2: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER2

      ! skip if this headerIndex doesn't contain upper air obs
      upperAirObs = .false.
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY3: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY3

        select case (obs_bodyElem_i(obsdat,OBS_VNM,bodyIndex))
        case (bufr_neuu, bufr_nevv, bufr_nett, bufr_nees)
          upperAirObs = .true.
        end select
      end do BODY3
      if (.not. upperAirObs) cycle HEADER2

      levStnIndex = levStnIndex + 1

      call obs_set_current_body_list(obsdat, headerIndex)
      BODY4: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY4

        select case (obs_bodyElem_i(obsdat,OBS_VNM,bodyIndex))
        case (bufr_neuu)
          call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, obsFlags(1,levStnIndex))
        case (bufr_nevv)
          call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, obsFlags(2,levStnIndex))
        case (bufr_nett)
          call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, obsFlags(3,levStnIndex))
        case (bufr_nees)
          call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, obsFlags(4,levStnIndex))
        end select
      end do BODY4

    end do HEADER2

    ! Deallocate arrays
    deallocate(surfPresInterp)
    deallocate(obsLevOffset)
    deallocate(obsType)
    deallocate(obsHeadDate)
    deallocate(obsHeadTime)
    deallocate(obsDate)
    deallocate(obsTime)
    deallocate(obsLaunchTime)
    deallocate(stationFlags)
    deallocate(obsLat)
    deallocate(obsLon)
    deallocate(stnId)
    deallocate(presInterp)
    deallocate(trajFlags)
    deallocate(obsFlags)
    deallocate(obsValues)
    deallocate(oMinusB)
    deallocate(obsLevOffsetMpi)
    deallocate(obsHeadDateMpi)
    deallocate(obsLaunchTimeMpi)
    deallocate(stationFlagsMpi)
    deallocate(obsLatMpi)
    deallocate(obsLonMpi)
    deallocate(trajFlagsMpi)
    deallocate(obsFlagsMpi)
    deallocate(obsValuesMpi)
    deallocate(oMinusBMpi)
    deallocate(stnIdMpi)

  end subroutine thn_radiosonde

  !---------------------------------------------------------------------
  ! Following are several subroutines needed by thn_radiosonde
  !---------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! stringArrayToMpi
  !--------------------------------------------------------------------------
  subroutine stringArrayToMpi(array, arrayMpi)
    ! :Purpose: Do the equivalent of mpi_allgatherv for a string array
    
    implicit none

    ! Arguments:
    character(len=*) :: array(:)
    character(len=*) :: arrayMpi(:)

    ! Locals:
    integer :: ierr, arrayIndex, charIndex, lenString
    integer :: nsize, nsizeMpi, allnsize(mpi_nprocs)
    integer, allocatable :: stringInt(:), stringIntMpi(:)

    nsize = size(array)
    call rpn_comm_allgather( nsize,    1, 'mpi_integer',  &
                             allnsize, 1, 'mpi_integer', &
                             'GRID', ierr )
    nsizeMpi = sum(allnsize(:))

    allocate(stringInt(nsize))
    allocate(stringIntMpi(nsizeMpi))

    lenString = len(array(1))
    do charIndex = 1, lenString
      do arrayIndex = 1, nsize
        stringInt(arrayIndex) = iachar(array(arrayIndex)(charIndex:charIndex))
      end do

      call intArrayToMpi(stringInt, stringIntMpi)

      do arrayIndex = 1, nsizeMpi
        arrayMpi(arrayIndex)(charIndex:charIndex) = achar(stringIntMpi(arrayIndex))
      end do
    end do

    deallocate(stringInt)
    deallocate(stringIntMpi)
    
  end subroutine stringArrayToMpi

  !--------------------------------------------------------------------------
  ! intArrayToMpi
  !--------------------------------------------------------------------------
  subroutine intArrayToMpi(array, arrayMpi, is_obsLevOffset_opt)
    ! :Purpose: Do the equivalent of mpi_allgatherv for an integer array,
    !           but with special treatment if array of obsLevOffset.

    implicit none

    ! Arguments:
    integer :: array(:)
    integer :: arrayMpi(:)
    logical, optional :: is_obsLevOffset_opt

    ! Locals:
    integer :: ierr, procIndex, arrayIndex
    integer :: nsize, nsizeMpi, allnsize(mpi_nprocs), displs(mpi_nprocs)
    logical :: is_obsLevOffset
    integer, allocatable :: numLevels(:), numLevelsMpi(:)

    if (present(is_obsLevOffset_opt)) then
      is_obsLevOffset = is_obsLevOffset_opt
    else
      is_obsLevOffset = .false.
    end if

    if (is_obsLevOffset) then

      ! special treatment is requirement for the variable "obsLevOffset"

      nsize = size(array) - 1
      call rpn_comm_allgather( nsize,    1, 'mpi_integer',  &
                               allnsize, 1, 'mpi_integer', &
                               'GRID', ierr )
      nsizeMpi = sum(allnsize(:))

      if ( mpi_myid == 0 ) then
        displs(1) = 0
        do procIndex = 2, mpi_nprocs
          displs(procIndex) = displs(procIndex-1) + allnsize(procIndex-1)
        end do
      else
        displs(:) = 0
      end if

      allocate(numLevels(nsize))
      allocate(numLevelsMpi(nsizeMpi))
      do arrayIndex = 1, nsize
        numLevels(arrayIndex) = array(arrayIndex+1) - array(arrayIndex)
      end do
      call rpn_comm_gatherv( numLevels   , nsize, 'mpi_integer', &
                             numLevelsMpi, allnsize, displs, 'mpi_integer',  &
                             0, 'GRID', ierr )

      call rpn_comm_bcast(numLevelsMpi, nsizeMpi, 'mpi_integer',  &
                          0, 'GRID', ierr)
      arrayMpi(1) = 0
      do arrayIndex = 1, nsizeMpi
        arrayMpi(arrayIndex+1) = arrayMpi(arrayIndex) + numLevelsMpi(arrayIndex)
      end do
      deallocate(numLevels)
      deallocate(numLevelsMpi)

    else

      nsize = size(array)
      call rpn_comm_allgather( nsize,    1, 'mpi_integer',  &
                               allnsize, 1, 'mpi_integer', &
                               'GRID', ierr )
      nsizeMpi = sum(allnsize(:))

      if ( mpi_myid == 0 ) then
        displs(1) = 0
        do procIndex = 2, mpi_nprocs
          displs(procIndex) = displs(procIndex-1) + allnsize(procIndex-1)
        end do
      else
        displs(:) = 0
      end if

      call rpn_comm_gatherv( array   , nsize, 'mpi_integer', &
                             arrayMpi, allnsize, displs, 'mpi_integer',  &
                             0, 'GRID', ierr )

      call rpn_comm_bcast(arrayMpi, nsizeMpi, 'mpi_integer',  &
                          0, 'GRID', ierr)      
    end if

  end subroutine intArrayToMpi

  !--------------------------------------------------------------------------
  ! realArrayToMpi
  !--------------------------------------------------------------------------
  subroutine realArrayToMpi(array, arrayMpi)
    ! :Purpose: Do the equivalent of mpi_allgatherv for a real array,

    implicit none

    ! Arguments:
    real(4) :: array(:)
    real(4) :: arrayMpi(:)

    ! Locals:
    integer :: ierr, procIndex
    integer :: nsize, nsizeMpi, allnsize(mpi_nprocs), displs(mpi_nprocs)

    nsize = size(array)
    call rpn_comm_allgather( nsize,    1, 'mpi_integer',  &
                             allnsize, 1, 'mpi_integer', &
                             'GRID', ierr )
    nsizeMpi = sum(allnsize(:))

    if ( mpi_myid == 0 ) then
      displs(1) = 0
      do procIndex = 2, mpi_nprocs
        displs(procIndex) = displs(procIndex-1) + allnsize(procIndex-1)
      end do
    else
      displs(:) = 0
    end if
    
    call rpn_comm_gatherv( array   , nsize, 'mpi_real4', &
                           arrayMpi, allnsize, displs, 'mpi_real4',  &
                           0, 'GRID', ierr )

    call rpn_comm_bcast(arrayMpi, nsizeMpi, 'mpi_real4',  &
                        0, 'GRID', ierr)
    
  end subroutine realArrayToMpi

  !--------------------------------------------------------------------------
  ! logicalArrayToMpi
  !--------------------------------------------------------------------------
  subroutine logicalArrayToMpi(array, arrayMpi)
    ! :Purpose: Do the equivalent of mpi_allgatherv for a logical array,

    implicit none

    ! Arguments:
    logical :: array(:)
    logical :: arrayMpi(:)

    ! Locals:
    integer :: ierr, procIndex
    integer :: nsize, nsizeMpi, allnsize(mpi_nprocs), displs(mpi_nprocs)

    nsize = size(array)
    call rpn_comm_allgather( nsize,    1, 'mpi_integer',  &
                             allnsize, 1, 'mpi_integer', &
                             'GRID', ierr )
    nsizeMpi = sum(allnsize(:))

    if ( mpi_myid == 0 ) then
      displs(1) = 0
      do procIndex = 2, mpi_nprocs
        displs(procIndex) = displs(procIndex-1) + allnsize(procIndex-1)
      end do
    else
      displs(:) = 0
    end if
    
    call rpn_comm_gatherv( array   , nsize, 'mpi_logical', &
                           arrayMpi, allnsize, displs, 'mpi_logical',  &
                           0, 'GRID', ierr )

    call rpn_comm_bcast(arrayMpi, nsizeMpi, 'mpi_logical',  &
                        0, 'GRID', ierr)
    
  end subroutine logicalArrayToMpi

  !--------------------------------------------------------------------------
  ! raobs_check_duplicated_stations
  !--------------------------------------------------------------------------
  subroutine raobs_check_duplicated_stations ( stnId, obsLevOffset, obsLat, obsLon,  &
                                               obsHeadDate, trajFlags, obsFlags, obsValues, &
                                               oMinusB, obsLaunchTime, stationFlags, &
                                               numVars, numStation )
    ! :Purpose: Check duplicated stations and select the best TAC/BUFR profiles

    implicit none

    integer,           intent(in)    :: numVars, numStation
    integer,           intent(in)    :: obsHeadDate(:), obsLaunchTime(:), stationFlags(:)
    real(4),           intent(in)    :: obsLat(:), obsLon(:)
    real(4),           intent(in)    :: obsValues(:,:), oMinusB(:,:)
    integer,           intent(in)    :: trajFlags(:,:)
    integer,           intent(in)    :: obsFlags(:,:)
    integer,           intent(in)    :: obsLevOffset(:)
    character (len=9), intent(inout) :: stnId(:)

    logical :: condition, sameProfile, stnidNotFound
    integer :: stationIndex, stationIndex2, stationIndex3, catIndex
    integer :: greaterNumVal, numDuplicate, numDuplicateTotal, selectStationIndex 
    integer :: bufrStationIndex, tacStationIndex, numChecked, numStnid, numSame
    integer, parameter :: maxNumStnid  = 5000
    character (len=9)  :: stnidList(maxNumStnid)
    integer            :: stationIndexList(maxNumStnid)
    integer :: numCriteria(5), cloche(30), selectCriteria

    numCriteria(:) = 0
    cloche(:) = 0

    ! Check for duplication of TAC or BUFR profiles

    numDuplicateTotal = 0
    do stationIndex = 1, numStation

      if ( stnId(stationIndex) /= 'NOT_VALID' ) then

        greaterNumVal = obsLevOffset(stationIndex+1) - obsLevOffset(stationIndex) + 1
        selectStationIndex = stationIndex
        numDuplicate = 0

        ! Verify if there exists two records with the same stnid, date, time, flgs.
        ! Keep the one with the greatest number of vertical levels
        do stationIndex2 = stationIndex+1, numStation

          if ( stnId(stationIndex2) /= 'NOT_VALID' ) then
            condition = stnId(stationIndex2) == stnId(stationIndex) .and. &
                        obsHeadDate(stationIndex2) == obsHeadDate(stationIndex) .and. &
                        obsLaunchTime(stationIndex2) == obsLaunchTime(stationIndex) .and. &
                        obsLat(stationIndex2) == obsLat(stationIndex) .and. &
                        obsLon(stationIndex2) == obsLon(stationIndex) .and. &
                        stationFlags(stationIndex2) == stationFlags(stationIndex)

            if ( condition ) then
              numDuplicate = numDuplicate + 1
              if ( (obsLevOffset(stationIndex2+1) - obsLevOffset(stationIndex2) + 1) > &
                  greaterNumVal ) then
                greaterNumVal = obsLevOffset(stationIndex2+1) - &
                                obsLevOffset(stationIndex2) + 1
                selectStationIndex  = stationIndex2
              end if
            end if
          end if

        end do ! stationIndex2

        ! Invalid all duplicated station except one with greatest number of levels
        if (numDuplicate > 0) then
          do stationIndex2 = stationIndex, numStation
            if ( stnId(stationIndex2) /= 'NOT_VALID' .and. &
                stnId(stationIndex2) == stnId(stationIndex) .and. &
                selectStationIndex /= stationIndex2 ) then
              write(*,'(2A20,I10,2F9.2,I10)') 'Station duplique ', stnId(stationIndex2), &
                   obsLevOffset(stationIndex2+1)-obsLevOffset(stationIndex2), &
                   obsLat(stationIndex2), obsLon(stationIndex2), &
                   obsLaunchTime(stationIndex2)
              stnId(stationIndex2) = 'NOT_VALID'
            end if
          end do
          numDuplicateTotal = numDuplicateTotal + numDuplicate
        end if

      end if

    end do ! stationIndex

    ! Select best profile for collocated TAC or BUFR reports
    if (mpi_myid == 0) then
      open(unit=11, file='./selected_stations_tac_bufr.txt', status='UNKNOWN')
    end if

    numChecked = 0
    do stationIndex = 1, numStation

      if ( stnId(stationIndex) /= 'NOT_VALID' ) then

        do stationIndex2 = stationIndex+1, numStation

          if ( stnId(stationIndex2) /= 'NOT_VALID' ) then
            condition = (stnId(stationIndex2) == stnId(stationIndex)) .and. &
                        ( btest(stationFlags(stationIndex),23) .neqv. &
                          btest(stationFlags(stationIndex2),23) )

            if ( condition ) then

              if (obsLaunchTime(stationIndex2) == obsLaunchTime(stationIndex)) then
                sameProfile = .true.
              else
                call raobs_check_if_same_profile(stationIndex,stationIndex2, &
                                                 obsValues,obsLevOffset,sameProfile)
              end if

              if ( sameProfile ) then

                numChecked  = numChecked + 1

                call raobs_compare_profiles(stationIndex, stationIndex2, &
                                            stationFlags, trajFlags, obsFlags, obsValues, &
                                            oMinusB, obsLevOffset, &
                                            numVars, cloche, selectCriteria, &
                                            selectStationIndex)

                numCriteria(selectCriteria) = numCriteria(selectCriteria) + 1

                bufrStationIndex = stationIndex2
                if (      btest(stationFlags(stationIndex),23) ) bufrStationIndex = stationIndex
                tacStationIndex  = stationIndex2
                if ( .not.btest(stationFlags(stationIndex),23) ) tacStationIndex = stationIndex

                if (mpi_myid == 0) write(11,'(A9,2F10.3,3I10)') stnId(stationIndex), &
                     obsLat(stationIndex), obsLon(stationIndex), &
                     obsLevOffset(bufrStationIndex+1)-obsLevOffset(bufrStationIndex)+1, &
                     obsLevOffset(tacStationIndex+1)-obsLevOffset(tacStationIndex)+1, selectCriteria

                if ( selectStationIndex == stationIndex) stnId(stationIndex2) = 'NOT_VALID'
                if ( selectStationIndex == stationIndex2) stnId(stationIndex) = 'NOT_VALID'

              end if ! sameProfile

            end if ! condition 

          end if ! stnId(stationIndex2) /= 'NOT_VALID'

        end do ! stationIndex2

      end if ! stnId(stationIndex) /= 'NOT_VALID' 

    end do ! stationIndex

    write(*,*)
    write(*,*)
    write(*,'(a30,I10)') 'nb of total duplicates '  ,numDuplicateTotal
    write(*,'(a30,I10)') 'nb TAC vs BUFR checked'   ,numChecked
    write(*,'(a30,I10)') 'nb not selected'          ,numCriteria(1)
    write(*,'(a30,I10)') 'nb suspicious traj BUFR'  ,numCriteria(2)
    write(*,'(a30,I10)') 'nb Energy higher   BUFR'  ,numCriteria(3)
    write(*,'(a30,I10)') 'nb variables lower BUFR'  ,numCriteria(4)
    write(*,'(a30,I10)') 'nb TAC  selected'         ,numCriteria(2) + &
                                    numCriteria(3) + numCriteria(4)
    write(*,'(a30,I10)') 'nb BUFR selected'         ,numCriteria(5)
    write(*,*)
    write(*,*)

    do catIndex = 1, 29
      write(*,'(a15,f4.2,a3,f4.2,I5)') 'nb e ratio ',(catIndex/10.)-.05,' - ', &
                                       (catIndex/10.)+.05,cloche(catIndex)
    end do
    write(*,*)
    write(*,'(a30,I5)') 'nb of energyTot(2)) very small  ',cloche(30)
    write(*,*)

    if (mpi_myid == 0) close(unit=11)

    ! Check whether there is still duplications
    numStnid = 0
    do stationIndex = 1, numStation

      if ( stnId(stationIndex) /= 'NOT_VALID' ) then

        if (numStnid < maxNumStnid ) then
          stnidNotFound=.true.
          if ( numStnid == 0) then
            numStnid = numStnid + 1
            stnidList(numStnid) = stnId(stationIndex)
            stationIndexList(numStnid) = stationIndex
          else
            numSame = 0
            do stationIndex2 = 1, numStnid
              if ( stnidList(stationIndex2) == stnId(stationIndex) ) then
                stnidNotFound=.false.
                numSame = numSame + 1
                stationIndex3 = stationIndexList(stationIndex2)
              end if
            end do
            if ( stnidNotFound ) then
              numStnid = numStnid + 1
              stnidList(numStnid) = stnId(stationIndex)
              stationIndexList(numStnid) = stationIndex
            else
              write(*,*) 'Multi profiles found : ',stnId(stationIndex), &
                   stnId(stationIndex3),numSame,stationIndex,stationIndex3, &
                   btest(stationFlags(stationIndex),23),btest(stationFlags(stationIndex3),23)
              write(*,'(a30,2i10,2f10.2,i10)') 'date, lch, lat lon ', &
                   obsHeadDate(stationIndex),obsLaunchTime(stationIndex), &
                   obsLat(stationIndex), obsLon(stationIndex), &
                   obsLevOffset(stationIndex+1)-obsLevOffset(stationIndex)+1
              write(*,'(a30,2i10,2f10.2,i10)') 'date, lch, lat lon ', &
                   obsHeadDate(stationIndex3),obsLaunchTime(stationIndex3), &
                   obsLat(stationIndex3), obsLon(stationIndex3), &
                   obsLevOffset(stationIndex3+1)-obsLevOffset(stationIndex3)+1
              write(*,*)

            end if
          end if
        else
          write(*,*) 'numStnid >= ',maxNumStnid
          call utl_abort('raobs_check_duplicated_stations')
        end if

      end if

    end do

  end subroutine raobs_check_duplicated_stations

  !--------------------------------------------------------------------------
  ! raobs_check_if_same_profile
  !--------------------------------------------------------------------------
  subroutine raobs_check_if_same_profile( stationIndex, stationIndex2, &
                                          obsValues, obsLevOffset, sameProfile )
    ! :Purpose: Check if two raobs profiles are the same.

    implicit none

    integer, intent(in)    :: stationIndex, stationIndex2
    real(4), intent(in)    :: obsValues(:,:)
    integer, intent(in)    :: obsLevOffset(:)
    logical, intent(out)   :: sameProfile

    integer :: varIndex, levIndex, levStnIndex, levStnIndex1, levStnIndex2, numSum
    real(4) :: valSum, minDeltaP1, minDeltaP2

    ! Standard levels
    integer, parameter :: numStdLevels = 16
    real(4) :: standardLevels(numStdLevels)
    standardLevels = (/ 1000.,925.,850.,700.,500.,400.,300.,250.,200.,150., &
                        100.,70.,50.,30.,20.,10./) 

    valSum = 0.0
    numSum   = 0

    sameProfile = .false.

    do levIndex = 1, numStdLevels

      minDeltaP1 = 1000.
      do levStnIndex = obsLevOffset(stationIndex)+1, obsLevOffset(stationIndex+1)
        if ( abs(standardLevels(levIndex) - obsValues(5,levStnIndex)) < minDeltaP1 ) then
          levStnIndex1 = levStnIndex
          minDeltaP1 = abs(standardLevels(levIndex) - obsValues(5,levStnIndex))
        end if
      end do
      minDeltaP2 = 1000.
      do levStnIndex = obsLevOffset(stationIndex2)+1, obsLevOffset(stationIndex2+1)
        if ( abs(standardLevels(levIndex) - obsValues(5,levStnIndex)) < minDeltaP2 ) then
          levStnIndex2 = levStnIndex
          minDeltaP2 = abs(standardLevels(levIndex) - obsValues(5,levStnIndex))
        end if
      end do

      if ( minDeltaP1 < 1.0 .and. minDeltaP2 < 1.0 ) then
        do varIndex = 2, 4, 2

          if ( obsValues(varIndex,levStnIndex1) /= -999.0 .and. &
               obsValues(varIndex,levStnIndex2) /= -999.0 ) then
            valSum = valSum + abs(obsValues(varIndex,levStnIndex1) - &
                                  obsValues(varIndex,levStnIndex2))
            numSum = numSum + 1
          end if

        end do

      end if

    end do !levIndex

    if (numSum > 0) then
      if (valSum/numSum < 0.5) sameProfile = .true.
    end if

  end subroutine raobs_check_if_same_profile

  !--------------------------------------------------------------------------
  ! raobs_compare_profiles
  !--------------------------------------------------------------------------
  subroutine raobs_compare_profiles( stationIndex, stationIndex2, stationFlags, trajFlags, &
                                     obsFlags, obsValues, oMinusB, obsLevOffset, numVars, &
                                     cloche, selectCriteria, selectStationIndex )
    ! :Purpose: Perform a comparison between two raobs profiles.

    implicit none

    integer,           intent(in)    :: stationIndex, stationIndex2, numVars
    integer,           intent(in)    :: trajFlags(:,:)
    integer,           intent(in)    :: obsFlags(:,:)
    real(4),           intent(in)    :: obsValues(:,:), oMinusB(:,:)
    integer,           intent(in)    :: obsLevOffset(:), stationFlags(:)
    integer,           intent(inout) :: cloche(:)
    integer,           intent(out)   :: selectCriteria, selectStationIndex

    logical :: condition, tacAndBufr, trajInfoOk
    integer :: raobFormatIndex, varIndex, levStnIndex, levStnIndex2, catIndex
    integer :: bufrStationIndex, tacStationIndex, countTimeFlag, countLatFlag, countTraj
    integer :: countValues(2), thisStationIndex
    real(4) :: presBottom, ombBottom, deltaPres, sumDeltaPres, oMinusBavg
    real(4) :: sumEnergy, presLower, presUpper
    real(4) :: percentTrajBad, toleranceFactor
    real(4) :: sumEnergyVar(numVars,2), energyTot(2)
    real(4) :: presLowerTacBufr(numVars,2), presUpperTacBufr(numVars,2)

    percentTrajBad = 10.0
    toleranceFactor     =  1.4
    catIndex = 1
    selectCriteria = 1

    tacAndBufr       = .true.
    trajInfoOk = .true.

    if ( btest(stationFlags(stationIndex),23) .and. &
         btest(stationFlags(stationIndex2),23) ) tacAndBufr = .false.
    if ( .not.btest(stationFlags(stationIndex),23) .and. &
         .not.btest(stationFlags(stationIndex2),23) ) tacAndBufr = .false.

    if ( tacAndBufr ) then

      bufrStationIndex = stationIndex2
      if (      btest(stationFlags(stationIndex),23) ) bufrStationIndex = stationIndex
      tacStationIndex  = stationIndex2
      if ( .not.btest(stationFlags(stationIndex),23) ) tacStationIndex  = stationIndex
      selectStationIndex = tacStationIndex

      ! 1. Evalue si la trajectoire native est correct

      if ( btest(stationFlags(bufrStationIndex),14) ) then

        countTimeFlag = 0
        countLatFlag = 0
        countTraj = obsLevOffset(bufrStationIndex+1) - obsLevOffset(bufrStationIndex) + 1

        do levStnIndex = obsLevOffset(bufrStationIndex)+1, obsLevOffset(bufrStationIndex+1)
          if ( btest(trajFlags(1,levStnIndex),4) ) countTimeFlag = countTimeFlag + 1
          if ( btest(trajFlags(2,levStnIndex),4) ) countLatFlag  = countLatFlag + 1
        end do

        if (100.*countTimeFlag/countTraj > percentTrajBad .or.  &
            100.*countLatFlag/countTraj  > percentTrajBad) then
          trajInfoOk = .false.
          selectStationIndex = tacStationIndex
          selectCriteria = 2
        end if

      end if

      if (trajInfoOk) then

        ! Cherche les pressions (bas et haut) des extremites des profils (TAC et BUFR)

        do raobFormatIndex = 1, 2

          if (raobFormatIndex == 1) thisStationIndex = bufrStationIndex
          if (raobFormatIndex == 2) thisStationIndex = tacStationIndex

          do varIndex = 1, 3

            presLowerTacBufr(varIndex,raobFormatIndex) = 1000.0
            do levStnIndex = obsLevOffset(thisStationIndex)+1, &
                             obsLevOffset(thisStationIndex+1)
              condition = oMinusB(varIndex,levStnIndex) /= -999.0 .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),18) .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),16) .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),9)  .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),8)  .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),2)  .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),11)
              if ( condition ) then
                presLowerTacBufr(varIndex,raobFormatIndex) = obsValues(5,levStnIndex)
                exit
              end if
 
            end do
            presUpperTacBufr(varIndex,raobFormatIndex) = &
                 presLowerTacBufr(varIndex,raobFormatIndex)
 
            if (levStnIndex < obsLevOffset(thisStationIndex+1) ) then

              do levStnIndex2 = levStnIndex+1, obsLevOffset(thisStationIndex+1)

                condition = oMinusB(varIndex,levStnIndex2) /= -999.0 .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex2),18) .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex2),16) .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex2),9)  .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex2),8)  .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex2),2)  .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex2),11)
                if ( condition ) then
                  presUpperTacBufr(varIndex,raobFormatIndex) = obsValues(5,levStnIndex2)
                end if

              end do

            end if

          end do !varIndex

        end do !raobFormatIndex

        ! Calcul de la norme energie equivalente (NE)

        do raobFormatIndex = 1, 2

          if (raobFormatIndex == 1) thisStationIndex = bufrStationIndex
          if (raobFormatIndex == 2) thisStationIndex = tacStationIndex

          countValues(raobFormatIndex) = 0

          do varIndex = 1, 3

            presLower =  presLowerTacBufr(varIndex,1)
            if ( presLower > presLowerTacBufr(varIndex,2) ) then
              presLower = presLowerTacBufr(varIndex,2)
            end if
            presUpper =  presUpperTacBufr(varIndex,1)
            if ( presUpper < presUpperTacBufr(varIndex,2) ) then
              presUpper = presUpperTacBufr(varIndex,2)
            end if

            sumDeltaPres = 0.0
            sumEnergy = 0.0
            sumEnergyVar(varIndex,raobFormatIndex) = 0.0

            do levStnIndex = obsLevOffset(thisStationIndex)+1, &
                             obsLevOffset(thisStationIndex+1)

              condition = oMinusB(varIndex,levStnIndex) /= -999.0 .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),18) .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),16) .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),9)  .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),8)  .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),2)  .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),11)

              if ( condition ) then
                countValues(raobFormatIndex) = countValues(raobFormatIndex) + 1
                if ( obsValues(5,levStnIndex) <= presLower ) then
                  presBottom = obsValues(   5,levStnIndex)
                  ombBottom = oMinusB(varIndex,levStnIndex)
                  exit
                end if
              end if

            end do

            if (levStnIndex < obsLevOffset(thisStationIndex+1) ) then

              do levStnIndex2 = levStnIndex+1, obsLevOffset(thisStationIndex+1)

                condition = oMinusB(varIndex,levStnIndex2) /= -999.0 .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex2),18) .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex2),16) .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex2),9)  .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex2),8)  .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex2),2)  .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex2),11)
                if ( condition ) then
                  countValues(raobFormatIndex) = countValues(raobFormatIndex) + 1
                  if ( obsValues(5,levStnIndex2) >= presUpper ) then
                    deltaPres = log(presBottom) - log(obsValues(5,levStnIndex2))
                    oMinusBavg = ( ombBottom + oMinusB(varIndex,levStnIndex2) ) / 2 
                    sumEnergy = sumEnergy + deltaPres * oMinusBavg**2
                    sumDeltaPres = sumDeltaPres + deltaPres
                    presBottom = obsValues(   5,levStnIndex2)
                    ombBottom = oMinusB(varIndex,levStnIndex2)
                  end if
                end if

              end do

            end if
            if (sumDeltaPres > 0.) then
              sumEnergyVar(varIndex,raobFormatIndex) = sumEnergy/sumDeltaPres
            end if

          end do ! varIndex

          energyTot(raobFormatIndex) = sumEnergyVar(1,raobFormatIndex) + &
                                       sumEnergyVar(2,raobFormatIndex) + &
                                       (1005./300.)*sumEnergyVar(3,raobFormatIndex)

        end do ! raobFormatIndex

        ! 2. Choisi le profil TAC (raobFormatIndex==2) si le nombre de niveaux ou
        !    il y a des donnees est plus grand ou egal au nombre dans le profil
        !    BUFR (raobFormatIndex==1)

        if ( countValues(2) >= countValues(1) ) then

          selectStationIndex = tacStationIndex
          selectCriteria = 4

        else

          ! 3. Choisi le profil BUFR si la norme energie equivalente du profil
          !    ne depasse pas 'toleranceFactor'

          if (abs(energyTot(2)) > 1.e-6) catIndex = nint(10.*energyTot(1)/energyTot(2)) + 1
          if (catIndex > 30) catIndex = 30 

          cloche(catIndex) = cloche(catIndex) + 1

          if ( energyTot(1) <= toleranceFactor*energyTot(2) ) then

            selectStationIndex = bufrStationIndex
            selectCriteria = 5

          else

            selectStationIndex = tacStationIndex
            selectCriteria = 3

          end if

        end if

      end if ! trajInfoOk

    else

      write(*,*) 'MEME TYPE DE FICHIER' 

    end if !tacAndBufr

  end subroutine raobs_compare_profiles

  !--------------------------------------------------------------------------
  ! raobs_thinning_model
  !--------------------------------------------------------------------------
  subroutine raobs_thinning_model( obsFlags, obsValues, presInterp, numVars, numLev, &
                                   numStation, numLevStnMax, obsLevOffset )
    ! :Purpose: Perform raobs thinning by comparing with a set of model levels.

    implicit none

    integer,           intent(in)    :: numVars, numLev, numStation, numLevStnMax
    integer,           intent(in)    :: obsLevOffset(:)
    integer,           intent(inout) :: obsFlags(:,:)
    real(4),           intent(in)    :: obsValues(:,:)
    real(4),           intent(in)    :: presInterp(:,:)

    logical :: condition, debug
    integer :: stationIndex, levIndex, levStnIndex, stdLevelIndex, varIndex
    integer :: varIndexDD, varIndexFF, varIndexPres
    integer :: levSelectIndex, levSelectIndex2, numLevSelect, levStnIndexValid
    integer :: numLevSelectBest, tempInt
    real(4) :: presTop, presBottom, deltaPres, deltaPresMin
    integer, allocatable :: levStnIndexList(:), numValidObs(:), listIndex(:)

    ! Standard levels including 925 hPa
    integer, parameter :: numStdLevels = 16
    real(4) :: standardLevels(numStdLevels)
    standardLevels(1:numStdLevels) = (/ 1000.,925.,850.,700.,500.,400.,300.,250.,200.,150., &
                                        100.,70.,50.,30.,20.,10./) 

    varIndexDD = 1
    varIndexFF = 2
    varIndexPres = 5
    debug = .false.

    allocate(levStnIndexList(numLevStnMax))
    allocate(numValidObs(numLevStnMax))
    allocate(listIndex(numLevStnMax))

    do stationIndex = 1, numStation

      ! If one obsFlags is set to 0 but not the other then make the obsFlagss consistent
      ! and obsFlags wind direction and module if one of these wind variables is missing

      do levStnIndex = obsLevOffset(stationIndex)+1, obsLevOffset(stationIndex+1)

        if ( obsFlags(varIndexDD,levStnIndex) == 0  .and. &
             obsFlags(varIndexFF,levStnIndex) /= 0 ) then
          obsFlags(varIndexDD,levStnIndex) = obsFlags(varIndexFF,levStnIndex)
        end if
        if ( obsFlags(varIndexDD,levStnIndex) /= 0  .and. &
             obsFlags(varIndexFF,levStnIndex) == 0 ) then
          obsFlags(varIndexFF,levStnIndex) = obsFlags(varIndexDD,levStnIndex)
        end if

        if ( (obsValues(varIndexDD,levStnIndex)  < 0. .and. &
              obsValues(varIndexFF,levStnIndex) >= 0.) .or. &
             (obsValues(varIndexDD,levStnIndex) >= 0. .and. &
              obsValues(varIndexFF,levStnIndex)  < 0.) ) then
          obsFlags(varIndexDD,levStnIndex) = ibset(obsFlags(varIndexDD,levStnIndex),11)
          obsFlags(varIndexFF,levStnIndex) = ibset(obsFlags(varIndexFF,levStnIndex),11)
        end if

      end do

      do levIndex = 1, numLev

        ! Pressures at the bottom and top of the model layer

        if (levIndex == 1) then
          presTop = presInterp(stationIndex,levIndex)
          presBottom = exp( 0.5*alog(presInterp(stationIndex,levIndex+1) * &
                       presInterp(stationIndex,levIndex)) )
        else if (levIndex == numLev) then
          presTop = exp( 0.5*alog(presInterp(stationIndex,levIndex-1) * &
                       presInterp(stationIndex,levIndex)) )
          presBottom = presInterp(stationIndex,levIndex)
        else
          presTop = exp( 0.5*alog(presInterp(stationIndex,levIndex-1) * &
                       presInterp(stationIndex,levIndex)) )
          presBottom = exp( 0.5*alog(presInterp(stationIndex,levIndex+1) * &
                       presInterp(stationIndex,levIndex)) )
        end if

        ! Select the observation levels between presTop and presBottom

        numLevSelect = 0
        do levStnIndex = obsLevOffset(stationIndex)+1, obsLevOffset(stationIndex+1)

          if (obsValues(varIndexPres,levStnIndex) > presTop  .and. &
              obsValues(varIndexPres,levStnIndex) <= presBottom) then

            numLevSelect = numLevSelect + 1
            levStnIndexList(numLevSelect) = levStnIndex
            numValidObs(numLevSelect) = 0

            ! numValidObs(numLevSelect) contains number of valid obs at each level selected

            do varIndex = 1, numVars
              condition = obsValues(varIndex,levStnIndex) >= 0. .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),18) .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),16) .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),9)  .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),8)  .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),2)  .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),11)
              if ( condition ) numValidObs(numLevSelect) = numValidObs(numLevSelect) + 1
            end do

          end if

        end do ! levStnIndex

        if ( numLevSelect > 0 ) then

          ! Sort the observation levels by greatest numbers of valid observations

          do levSelectIndex = 1, numLevSelect
            listIndex(levSelectIndex) = levSelectIndex
          end do
          do levSelectIndex = 1, numLevSelect
            do levSelectIndex2 = levSelectIndex, numLevSelect
              if ( numValidObs(levSelectIndex2) > numValidObs(levSelectIndex) ) then
                tempInt     = numValidObs(levSelectIndex)
                numValidObs(levSelectIndex) = numValidObs(levSelectIndex2)
                numValidObs(levSelectIndex2) = tempInt
                tempInt     = listIndex(levSelectIndex)
                listIndex(levSelectIndex) = listIndex(levSelectIndex2)
                listIndex(levSelectIndex2) = tempInt 
              end if
            end do
          end do
          numLevSelectBest = 1
          if ( numLevSelect > 1 ) then
            do levSelectIndex = 2, numLevSelect
              if (numValidObs(levSelectIndex) == numValidObs(1)) then
                numLevSelectBest = numLevSelectBest + 1
              end if
            end do
          end if

          do varIndex = 1, numVars

            levStnIndexValid = 0

            ! Rule #1 : get valid the observation at the standard level if one
            !           is found in the layer

            do levSelectIndex = 1, numLevSelect

              levStnIndex = levStnIndexList(levSelectIndex)

              condition = obsValues(varIndex,levStnIndex) >= 0. .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),18) .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),16) .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),9)  .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),8)  .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),2)  .and. &
                   .not.btest(obsFlags(varIndex,levStnIndex),11)

              do stdLevelIndex = 1, numStdLevels
                if ( (obsValues(varIndexPres,levStnIndex) == &
                      standardLevels(stdLevelIndex)) .and. &
                     condition ) then
                  levStnIndexValid = levStnIndex
                end if
              end do

            end do ! levSelectIndex

            if ( levStnIndexValid == 0 ) then

              ! Rule #2 : get the closest valid observation to the model level
              !           and most complete

              deltaPresMin = presBottom - presTop
              do levSelectIndex = 1, numLevSelectBest

                levStnIndex = levStnIndexList(listIndex(levSelectIndex))

                deltaPres = abs(obsValues(varIndexPres,levStnIndex) - &
                                presInterp(stationIndex,levIndex))
                condition = obsValues(varIndex,levStnIndex) >= 0. .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex),18) .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex),16) .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex),9)  .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex),8)  .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex),2)  .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex),11)

                if ( deltaPres < deltaPresMin .and. condition) then
                  deltaPresMin = deltaPres
                  levStnIndexValid = levStnIndex
                end if

              end do ! levSelectIndex

            end if  ! levStnIndexValid == 0

            if ( levStnIndexValid == 0 ) then

              ! Rule #3 : get the closest valid observation if not previously selected

              deltaPresMin = presBottom - presTop
              do levSelectIndex = 1, numLevSelect

                levStnIndex = levStnIndexList(levSelectIndex)

                deltaPres = abs(obsValues(varIndexPres,levStnIndex) - &
                                presInterp(stationIndex,levIndex))
                condition = obsValues(varIndex,levStnIndex) >= 0. .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex),18) .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex),16) .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex),9)  .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex),8)  .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex),2)  .and. &
                     .not.btest(obsFlags(varIndex,levStnIndex),11)

                if ( deltaPres < deltaPresMin .and. condition) then
                  deltaPresMin = deltaPres
                  levStnIndexValid = levStnIndex
                end if

              end do ! levSelectIndex

            end if  ! levStnIndexValid == 0

            ! Apply thinning obsFlags to all obs except the one on level levStnIndexValid

            do levSelectIndex = 1, numLevSelect

              levStnIndex = levStnIndexList(levSelectIndex)

              if ( levStnIndex /= levStnIndexValid .and. &
                   obsValues(varIndex,levStnIndex) >= 0. ) then
                obsFlags(varIndex,levStnIndex) = ibset(obsFlags(varIndex,levStnIndex),11)
              end if

            end do

          end do !varIndex

        end if !numLevSelect > 0

      end do !levIndex

      ! Following lines for debugging purposes

      if ( debug ) then

        levIndex = numLev
        presBottom = presInterp(stationIndex,levIndex)
        write (*,*) ' '
        write (*,'(a40,I8,a40)') '      ==================== Station no. ', &
             stationIndex,' =============================='
        write (*,*) ' '
        write (*,'(a80,f10.2,i10)') 'lowest model layer-------------> ',presBottom,levIndex
        do while(obsValues(varIndexPres,obsLevOffset(stationIndex)+1) < presBottom)
          levIndex = levIndex - 1
          presBottom = exp(0.5 * alog(presInterp(stationIndex,levIndex+1) * &
                           presInterp(stationIndex,levIndex)))
        end do
        if (levIndex == numLev) levIndex = numLev - 1

        do levStnIndex = obsLevOffset(stationIndex)+1, obsLevOffset(stationIndex+1)

          presBottom = exp( 0.5*alog(presInterp(stationIndex,levIndex+1) * &
                       presInterp(stationIndex,levIndex)) )
          do while(obsValues(varIndexPres,levStnIndex) < presBottom)
            write (*,'(a80,f10.2,i10)') 'model layer---------------> ',presBottom,levIndex
            levIndex = levIndex - 1
            presBottom = exp( 0.5*alog(presInterp(stationIndex,levIndex+1) * &
                         presInterp(stationIndex,levIndex)) )
          end do
          write (*,'(5f8.2,4i8)') obsValues(varIndexPres,levStnIndex), &
               obsValues(1,levStnIndex), obsValues(2,levStnIndex), &
               obsValues(3,levStnIndex), obsValues(4,levStnIndex), &
               obsFlags(1,levStnIndex), obsFlags(2,levStnIndex), &
               obsFlags(3,levStnIndex), obsFlags(4,levStnIndex)

        end do

      end if

    end do !stationIndex

    deallocate(levStnIndexList, numValidObs, listIndex)

  end subroutine raobs_thinning_model

  !--------------------------------------------------------------------------
  ! raobs_thinning_es
  !--------------------------------------------------------------------------
  subroutine raobs_thinning_es( obsFlags, obsValues, numStation, &
                                numLevStnMax, obsLevOffset )
    ! :Purpose: Perform thinning of T-Td raobs observations.

    implicit none

    integer,           intent(in)    :: numStation, numLevStnMax
    integer,           intent(in)    :: obsLevOffset(:)
    integer,           intent(inout) :: obsFlags(:,:)
    real(4),           intent(in)    :: obsValues(:,:)

    integer, parameter  :: numLevES = 42
    logical :: condition
    integer :: stationIndex, levIndex, levStnIndex, levSelectIndex, varIndexES, varIndexPres
    integer :: numLevSelect, levStnIndexValid
    real(4) :: levelsES(numLevES), deltaPres, deltaPresMin, presBottom, presTop
    integer, allocatable :: levStnIndexList(:)

    ! Selected pressure levels for additional thinning to humitidy observations

    levelsES(1:numLevES) = (/ 1025.,1000.,975.,950.,925.,900.,875.,850.,825., &
                              800.,775.,750.,725.,700.,675.,650.,625.,600.,575., &
                              550.,525.,500.,475.,450.,425.,400.,375.,350.,325., &
                              300.,275.,250.,225.,200.,175.,150.,100., 70., 50., &
                              30., 20., 10./)

    varIndexES = 4
    varIndexPres = 5

    allocate(levStnIndexList(numLevStnMax))

    do stationIndex = 1, numStation

      ! Retain nearest observations to selected pressure levels in levelsES

      do levIndex = 1, numLevES

        ! Pressures at the bottom and top of the ES layer

        if (levIndex == 1) then
          presBottom = levelsES(levIndex)
          presTop = exp(0.5*alog(levelsES(levIndex+1)*levelsES(levIndex)))
        else if (levIndex == numLevES) then
          presBottom = exp(0.5*alog(levelsES(levIndex-1)*levelsES(levIndex)))
          presTop = levelsES(levIndex)
        else
          presBottom = exp(0.5*alog(levelsES(levIndex-1)*levelsES(levIndex)))
          presTop = exp(0.5*alog(levelsES(levIndex+1)*levelsES(levIndex)))
        end if
   
        ! Select the levels between presTop and presBottom

        numLevSelect = 0

        do levStnIndex = obsLevOffset(stationIndex)+1, obsLevOffset(stationIndex+1)

          if ( obsValues(varIndexPres,levStnIndex) > presTop .and. &
               obsValues(varIndexPres,levStnIndex) <= presBottom ) then
            numLevSelect  = numLevSelect + 1
            levStnIndexList(numLevSelect) = levStnIndex
          end if

        end do ! levStnIndex

        ! Get the closest valid observation to the selected pressure level

        levStnIndexValid = 0
        deltaPresMin = presBottom - presTop

        do levSelectIndex = 1, numLevSelect

          levStnIndex = levStnIndexList(levSelectIndex)

          deltaPres = abs(obsValues(varIndexPres,levStnIndex) - levelsES(levIndex))
          condition = .not.btest(obsFlags(varIndexES,levStnIndex),18) .and. &
               .not.btest(obsFlags(varIndexES,levStnIndex),16) .and. &
               .not.btest(obsFlags(varIndexES,levStnIndex),9)  .and. &
               .not.btest(obsFlags(varIndexES,levStnIndex),8)  .and. &
               .not.btest(obsFlags(varIndexES,levStnIndex),2)  .and. &
               .not.btest(obsFlags(varIndexES,levStnIndex),11) .and. &
               obsValues(varIndexES,levStnIndex) >= 0.         .and. &
               deltaPres < deltaPresMin

          if ( condition ) then
            deltaPresMin = deltaPres
            levStnIndexValid = levStnIndex
          end if

        end do !levSelectIndex

        ! Apply thinning obsFlags to all obs except the one selected with levStnIndexValid

        do levSelectIndex = 1, numLevSelect

          levStnIndex = levStnIndexList(levSelectIndex)

          if ( (levStnIndex /= levStnIndexValid) .and. &
               (obsValues(varIndexES,levStnIndex) >= 0.) ) then
            obsFlags(varIndexES,levStnIndex) = ibset(obsFlags(varIndexES,levStnIndex),11)
          end if

        end do

      end do !levIndex

    end do !stationIndex

    deallocate(levStnIndexList)

  end subroutine raobs_thinning_es

  !--------------------------------------------------------------------------
  ! raobs_blacklisting_ecmwf
  !--------------------------------------------------------------------------
  subroutine raobs_blacklisting_ecmwf( obsFlags, obsValues, obsType, numStation, obsLevOffset )
    ! :Purpose: Perform filtering of T-Td raobs observations based on
    !           approach inspired by ECMWF approach

    implicit none

    integer, intent(in)    :: numStation
    integer, intent(in)    :: obsLevOffset(:), obsType(:)
    integer, intent(inout) :: obsFlags(:,:)
    real(4), intent(in)    :: obsValues(:,:)

    logical :: conditionTT, conditionES, rejectES
    integer :: ierr, stationIndex, levStnIndex, varIndexES, varIndexTT, varIndexPres
    integer :: type, countReject0p, countReject0t, countReject1p, countReject1t
    integer :: countReject0pMpi, countReject0tMpi, countReject1pMpi, countReject1tMpi
    integer :: countTotalES, countTotalESMpi

    countReject0p = 0
    countReject0t = 0
    countReject1p = 0
    countReject1t = 0
    countTotalES = 0

    varIndexTT = 3
    varIndexES = 4
    varIndexPres = 5

    do stationIndex = 1, numStation
      type = 0

      if ( obsType(stationIndex) /= -1 ) then

        if ( obsType(stationIndex) == 14  .or.  &
             obsType(stationIndex) == 24  .or.  &
             obsType(stationIndex) == 41  .or.  &
             obsType(stationIndex) == 42  .or.  &
             obsType(stationIndex) == 52  .or.  &
             obsType(stationIndex) == 79  .or.  &
             obsType(stationIndex) == 80  .or.  &
             obsType(stationIndex) == 81  .or.  &
             obsType(stationIndex) == 83  .or.  &
             obsType(stationIndex) == 141 .or.  &
             obsType(stationIndex) == 142 ) then

          type = 1

        end if

      end if

      do levStnIndex = obsLevOffset(stationIndex)+1, obsLevOffset(stationIndex+1)

        conditionTT =  &
             .not.btest(obsFlags(varIndexTT,levStnIndex),18) .and. &
             .not.btest(obsFlags(varIndexTT,levStnIndex),16) .and. &
             .not.btest(obsFlags(varIndexTT,levStnIndex),9)  .and. &
             .not.btest(obsFlags(varIndexTT,levStnIndex),8)  .and. &
             .not.btest(obsFlags(varIndexTT,levStnIndex),2)  .and. &
             .not.btest(obsFlags(varIndexTT,levStnIndex),11) .and. &
             obsValues(varIndexTT,levStnIndex) >= 0.0
        conditionES =  &
             .not.btest(obsFlags(varIndexES,levStnIndex),18) .and. &
             .not.btest(obsFlags(varIndexES,levStnIndex),16) .and. &
             .not.btest(obsFlags(varIndexES,levStnIndex),9)  .and. &
             .not.btest(obsFlags(varIndexES,levStnIndex),8)  .and. &
             .not.btest(obsFlags(varIndexES,levStnIndex),2)  .and. &
             .not.btest(obsFlags(varIndexES,levStnIndex),11) .and. &
             obsValues(varIndexES,levStnIndex) >= 0.0

        if ( conditionTT ) then

          rejectES = .false. 
          if (conditionES) countTotalES = countTotalES + 1

          if ( type == 0 ) then

            if ( .not.rejectES .and. obsValues(varIndexPres,levStnIndex) < 300.0 ) then
              if (conditionES) countReject0p = countReject0p + 1
              rejectES = .true. 
            end if
            if ( .not.rejectES .and. obsValues(varIndexTT,levStnIndex) < 233.15 ) then
              if (conditionES) countReject0t = countReject0t + 1
              rejectES = .true. 
            end if
    
          else if  ( type == 1 ) then

            if ( .not.rejectES .and. obsValues(varIndexPres,levStnIndex) < 100.0 ) then
              if (conditionES) countReject1p = countReject1p + 1
              rejectES = .true. 
            end if
            if ( .not.rejectES .and. obsValues(varIndexTT,levStnIndex) < 213.15 ) then
              if (conditionES)  countReject1t = countReject1t + 1
              rejectES = .true. 
            end if

          end if

          if ( rejectES ) then
            obsFlags(varIndexES,levStnIndex) = ibset(obsFlags(varIndexES,levStnIndex),8)
          end if

        else

          obsFlags(varIndexES,levStnIndex) = obsFlags(varIndexTT,levStnIndex)

        end if

      end do !levStnIndex

    end do !stationIndex

    
    call rpn_comm_allReduce(countTotalES, countTotalESMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(countReject0p, countReject0pMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(countReject0t, countReject0tMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(countReject1p, countReject1pMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    call rpn_comm_allReduce(countReject1t, countReject1tMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    write(*,*)
    write(*,*) ' Rejet des donnees ES inspire de ECMWF'
    write(*,*)
    write(*,'(a30,I10)') 'nb total donnees es        ',countTotalESMpi
    write(*,'(a30,I10)') 'nb rejet type 0 pression   ',countReject0pMpi
    write(*,'(a30,I10)') 'nb rejet type 0 temperature',countReject0tMpi
    write(*,'(a30,I10)') 'nb rejet type 1 pression   ',countReject1pMpi
    write(*,'(a30,I10)') 'nb rejet type 1 temperature',countReject1tMpi
    write(*,*)
    if (countTotalESMpi > 0) then
      write(*,'(a30,I10,f10.2)') 'nb rejet total et %', &
           countReject0pMpi + countReject0tMpi + countReject1pMpi + countReject1tMpi, &
           100.0 * (countReject0pMpi + countReject0tMpi + countReject1pMpi + &
                    countReject1tMpi) / countTotalESMpi
    end if
    write(*,*)

  end subroutine raobs_blacklisting_ecmwf
  
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
    integer :: ztdScoreCount, ztdScoreCountMpi
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
    ztdScoreCount = 0

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
      ! ZTD quality is unknown (missing ztd score)
      if (ztdScore == 999.0) then
        ztdScoreCount = ztdScoreCount + 1
        quality(headerIndex) = 9999
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
    call rpn_comm_allReduce(ztdScoreCount, ztdScoreCountMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)

    write(*,*)
    write(*,'(a50,i10)') 'Number of input obs                  = ', countObsInMpi
    write(*,'(a50,i10)') 'Total number of rejected/thinned obs = ', countObsInMpi - countObsOutMpi
    write(*,'(a50,i10)') 'Number of output obs                 = ', countObsOutMpi
    write(*,*)
    write(*,*)
    write(*,'(a50,i10)') 'Number of blacklisted obs            = ', blackListCountMpi
    write(*,'(a50,i10)') 'Number of obs without bias correction= ', unCorrectCountMpi
    write(*,'(a50,i10)') 'Number of obs with no ztdScore       = ', ztdScoreCountMpi
    write(*,'(a60,4i10)') 'Number of rejects outside time window, BGCK, thinning ', &
         badTimeCountMpi, bgckCountMpi, &
         countObsInMpi - countObsOutMpi - &
         bgckCountMpi - badTimeCountMpi - blackListCountMpi - unCorrectCountMpi - ztdScoreCountMpi

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
    call rpn_comm_allReduce(missingCount, missingCountMpi, 1, &
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
    ! :Purpose: Quick sort algorithm for integer data.

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
    ! :Purpose: Subroutine called in quick sort for integers.

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
    ! :Purpose: Quick sort algorithm for real8 data.

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
    ! :Purpose: Subroutine called for quick sort of real8 data.

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
  real(4) function thn_distanceArc( deltaLat, deltaLon, lat1, lat2 )
    ! :Purpose: Compute arc distance.

    implicit none

    real(4), intent(in) :: deltaLat, deltaLon, lat1, lat2

    real(4), parameter :: PI = 3.141592654
    real(4), parameter :: RT = 6374.893
    real(4) :: lat1_r, lat2_r, deltaLat_r, deltaLon_r, term_a

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
      call gio_readFromFile( stateVectorPsfc, trlmFileName, ' ', ' ',  &
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
    ! :Purpose: Of the observations in a column that have not already been
    !           rejected, keep every nth observation and throw out the rest.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.

    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in) :: familyType
    integer,          intent(in) :: keepNthVertical

    ! Locals:
    integer, parameter :: PROFILE_NOT_FOUND=-1
    integer :: headerIndex, bodyIndex
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

      ! If datum already rejected, ignore it
      flag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
      if ( btest(flag,9) .or. &
           btest(flag,11) ) cycle BODY

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
      ! :Purpose: Determine whether the current observation begins a new vertical column
      !           (Assume that observations are in chronological order)
      !
      ! :Note:  This method has been written with aladin observations in mind.
      !         It might be necessary to generalize the method.
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
  subroutine thn_tovsFilt(obsdat, delta, deltrad, codtyp, codtyp2_opt)
    ! :Purpose: Thinning algorithm used for AMSU and ATMS radiance obs.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.

    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    integer, intent(in)             :: delta
    integer, intent(in)             :: deltrad
    integer, intent(in)             :: codtyp
    integer, optional, intent(in)   :: codtyp2_opt

    ! Locals:
    integer :: numLat, numLon, headerIndex, headerIndexKeep, latIndex, lonIndex, latIndex2
    integer :: gridIndex, numGridLonsTotal, obsTime, obsDate, numHeader, numHeaderMaxMpi, ierr
    integer :: bodyIndex, stepIndex, obsIndex, obsFov
    integer :: loscan, hiscan, obsFlag, numObs, minLonBurpFileMpi(mpi_nprocs)
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
    character(len=codtyp_name_length) :: instrumName

    ! Local parameters:
    integer, parameter :: latLength=10000
    integer, parameter :: lonLength=40000
    integer, parameter :: mxscanamsua=30
    integer, parameter :: mxscanamsub=90
    integer, parameter :: mxscanatms =96
    integer, parameter :: mxscanmwhs2=98

    instrumName = codtyp_get_name(codtyp)
    write(*,*)
    write(*,*) 'thn_tovsFilt: Starting, ', trim(instrumName)
    write(*,*)

    numHeader = obs_numHeader(obsdat)
    call rpn_comm_allReduce(numHeader, numHeaderMaxMpi, 1, 'mpi_integer', &
                            'mpi_max','grid',ierr)

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

    countObs = count(valid(:))
    call rpn_comm_allReduce(countObs, countObsMpi, 1, 'mpi_integer', &
                            'mpi_sum','grid',ierr)
    if (countObsMpi == 0) then
      write(*,*) 'thn_tovsFilt: no observations for this instrument'
      deallocate(valid)
      return
    end if

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
        if (numObs > 0 .and. minDistance <= real(deltRad) ) then
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
    write(*,*) 'thn_tovsFilt: countObs after thinning                 = ', &
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
    ! :Purpose: Remove duplicate TOVS observations due to RARS.

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
      obsPosition3d(1,headerIndex) = ec_ra * sin(obsLonInRad) * cos(obsLatInRad)
      obsPosition3d(2,headerIndex) = ec_ra * cos(obsLonInRad) * cos(obsLatInRad)
      obsPosition3d(3,headerIndex) = ec_ra *                    sin(obsLatInRad)
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
    ! :Purpose: Only keep the observation closest to the center of each
    !           lat-lon (and time) box for SCAT observations.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.

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
    real(4) :: latInRadians, distance, obsLat, obsLon
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
    write(*,*) 'thn_scatByLatLonBoxes: countObs after QC and time tests   = ', &
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
    write(*,'(a,i6)') 'scatByLatLonBoxes: Number of obs in input  = ', countObsInMpi
    write(*,'(a,i6)') 'scatByLatLonBoxes: Number of obs in output = ', countObsOutMpi
    write(*,'(a,i6)') 'scatByLatLonBoxes: Number of obs not selected due to time = ', &
         timeRejectCountMpi
    write(*,'(a,i6)') 'scatByLatLonBoxes: Number of obs not selected due to topo = ', &
         flagRejectCountMpi
    write(*,*)
    
    write(*,'(a40,i10)' ) 'Number of satellites found = ', numStnId
    write(*,*)
  
    write(*,'(a40,a15)' ) 'Satellite', 'nb SCAT in'
    write(*,*)
    do stnIdIndex = 1, numStnId
      write(*,'(a40,i15)') stnidList(stnIdIndex), numObsStnIdInMpi(stnIdIndex)
    end do
    write(*,*)
    write(*,'(a40,i15)' ) 'Total number of obs in : ', sum(numObsStnIdInMpi(:))
  
    write(*,*)
    write(*,'(a40,a15)' ) 'Satellite', 'nb SCAT out'
    write(*,*)
    do stnIdIndex = 1, numStnId
      write(*,'(a40,i15)') stnidList(stnIdIndex), numObsStnIdOutMpi(stnIdIndex)
    end do
    write(*,*)
    write(*,'(a40,i15)' ) 'Total number of obs out : ', sum(numObsStnIdOutMpi(:))

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
    ! :Purpose: Only keep the observation closest to the center of each
    !           lat-lon (and time) box for CSR observations.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.

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
    obsAngle(:)             = 0.0
    obsCloud(:,:)           = 0.0
    obsDistance(:)          = 0.0

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

        obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
        if ( .not.btest(obsFlag,8) .and. &
             .not.btest(obsFlag,9) .and. &
             .not.btest(obsFlag,11) ) then
          valid(headerIndex) = .true.
          channelAssim(channelIndex,headerIndex) = .true.
        end if
      end do BODY1

      ! read the new element 20081 in obsspacedata OBS_CLA
      CHANNELS: do channelIndex = 1, numChannel(headerIndex)
        obsCloud(channelIndex, headerIndex) = -1.0

        ! search for the cloud information for this channel
        call obs_set_current_body_list(obsdat, headerIndex)
        BODY2: do
          bodyIndex = obs_getBodyIndex(obsdat)
          if (bodyIndex < 0) exit BODY2

          ! check if channel number matches
          if (nint(obs_bodyElem_r(obsdat, OBS_PPP, bodyIndex)) == channelList(channelIndex)) then
            obsCloud(channelIndex, headerIndex) = real(obs_bodyElem_i(obsdat, OBS_CLA, bodyIndex))
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

    ! Initialize arrays
    obsLatIndexMpi(:)    = 0
    obsLonIndexMpi(:)    = 0
    obsStepIndexMpi(:)   = 0
    numChannelMpi(:)     = 0
    obsAngleMpi(:)       = 0.0
    obsDistanceMpi(:)    = 0.0
    channelAssimMpi(:,:) = 0
    obsCloudMpi(:,:)     = 0.0
    stnIdInt(:,:)        = 0
    stnIdIntMpi(:,:)     = 0

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
    ! :Purpose: Only keep the observation closest to the center of each
    !           lat-lon (and time) box.
    !           Set bit 11 of OBS_FLG on observations that are to be rejected.

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
    integer :: numHeader, numHeaderMpi, numHeaderMaxMpi, lenStnId 
    integer :: numLat, numLon, latIndex, numChannels, delMinutes
    integer :: lonBinIndex, latBinIndex, timeBinIndex, charIndex
    integer :: ierr, nsize, procIndex, countHeader, countHeaderMpi
    real(4) :: latInRadians, length, distance
    real(8) :: lonBoxCenterInDegrees, latBoxCenterInDegrees
    real(8) :: obsLatInRad, obsLonInRad, obsLat, obsLon
    real(8) :: obsLatInDegrees, obsLonInDegrees, obsStepIndex_r8
    integer, parameter :: latLength = 10000
    integer, parameter :: lonLength = 40000
    real(4), allocatable :: gridLats(:)
    integer, allocatable :: numGridLons(:)
    integer, allocatable :: stnIdInt(:,:), stnIdIntMpi(:,:)
    integer, allocatable :: headerIndexKeep(:,:,:), numChannelsKeep(:,:,:)
    integer, allocatable :: headerIndexKeepMpi(:,:,:,:), numChannelsKeepMpi(:,:,:,:)
    integer, allocatable :: delMinutesKeep(:,:,:), delMinutesKeepMpi(:,:,:,:)
    integer, allocatable :: procIndexKeep(:,:,:)
    real(4), allocatable :: distanceKeep(:,:,:), distanceKeepMpi(:,:,:,:)
    logical, allocatable :: rejectThisHeader(:)
    logical :: keepThisObs, stnIdFoundInList
    integer :: obsLonBurpFile, obsLatBurpFile
    integer, parameter :: numStnIdMax = 10
    integer :: numStnId, stnIdIndexFound, stnIdIndex, countMpi, numObsStnId(numStnIdMax)
    character(len=12) :: stnid, stnidList(numStnIdMax)
    character(len=codtyp_name_length) :: instrumName

    instrumName = codtyp_get_name(codtyp)
    write(*,*)
    write(*,*) 'thn_hyperByLatLonBoxes: Starting, ', trim(instrumName)
    write(*,*)

    numHeader = obs_numHeader(obsdat)
    call rpn_comm_allReduce(numHeader, numHeaderMaxMpi, 1, 'mpi_integer', &
                            'mpi_max', 'grid', ierr)

    lenStnId = len(stnId)
    allocate(stnIdInt(lenStnId,numHeaderMaxMpi))
    allocate(stnIdIntMpi(lenStnId,numHeaderMaxMpi*mpi_nprocs))
    stnIdInt(:,:) = 0
    stnIdIntMpi(:,:) = 0

    ! loop over all header indices of the specified family and get integer stnId
    countHeader = 0
    call obs_set_current_header_list(obsdat,trim(familyType))
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER
      if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) /= codtyp) cycle HEADER

      countHeader = countHeader + 1

      ! convert and store stnId as integer array
      stnId = obs_elem_c(obsdat,'STID',headerIndex)
      do charIndex = 1, lenStnId
        stnIdInt(charIndex,headerIndex) = iachar(stnId(charIndex:charIndex))
      end do
    end do HEADER

    ! return if no observations for this instrument
    call rpn_comm_allReduce(countHeader, countHeaderMpi, 1, 'mpi_integer', &
                            'mpi_sum', 'grid', ierr)
    if (countHeaderMpi == 0) then
      write(*,*) 'thn_hyperByLatLonBoxes: no observations for this instrument'
      return
    end if

    ! Gather stnIdInt from all MPI tasks
    nsize = lenStnId * numHeaderMaxMpi
    call rpn_comm_allgather(stnIdInt,    nsize, 'mpi_integer',  &
                            stnIdIntMpi, nsize, 'mpi_integer', 'grid', ierr)

    ! build a global stnIdList
    numStnId = 0
    numHeaderMpi = numHeaderMaxMpi * mpi_nprocs
    HEADER4: do headerIndex = 1, numHeaderMpi
      if (all(stnIdIntMpi(:,headerIndex) == 0)) cycle HEADER4

      ! Station ID converted back to character string
      do charIndex = 1, lenStnId
        stnId(charIndex:charIndex) = achar(stnIdIntMpi(charIndex,headerIndex))
      end do

      stnIdFoundInList = .false.
      if ( numStnId == 0 ) then
        stnIdFoundInList = .true.
        numStnId = numStnId + 1
        stnidList(numStnId) = stnId
      else
        do stnIdIndex = 1, numStnId
          if ( stnidList(stnIdIndex) == stnId ) cycle HEADER4
        end do
      end if

      if ( .not. stnIdFoundInList ) then
        numStnId = numStnId + 1
        stnidList(numStnId) = stnId
      end if

      if ( numStnId >= numStnIdMax ) then
        call utl_abort('thn_hyperByLatLonBoxes: numStnId too large')
      end if

    end do HEADER4

    ! loop over local headers to find numObs for each stnId
    call obs_set_current_header_list(obsdat,trim(familyType))
    HEADER0: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER0
      if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) /= codtyp) cycle HEADER0

      stnId = obs_elem_c(obsdat,'STID',headerIndex)

      do stnIdIndex = 1, numStnId
        if ( stnidList(stnIdIndex) == stnId ) then
          numObsStnId(stnIdIndex) = numObsStnId(stnIdIndex) + 1
        end if
      end do

    end do HEADER0

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

    write(*,'(a)') ' '
    write(*,'(a)') ' == Gridbox properties == '
    write(*,'(a)') ' '
    write(*,'(a,i8)') ' Number of horizontal boxes : ', numLon
    write(*,'(a,i8)') ' Number of vertical boxes   : ', numLat
    write(*,'(a,i8)') ' Number of temporal bins    : ', tim_nstepobs

    ! print some statistics before thinning

    instrumName = codtyp_get_name(codtyp)
    write(*,*)
    write(*,'(a)')        ' == Input file == '
    write(*,*)
    write(*,'(a,a4)')     ' Instrument = ', trim(instrumName)
    write(*,'(a,i4)')     ' Codtyp     = ', codtyp

    write(*,*)
    write(*,'(a,2i8)')     ' Number of valid satellite profiles  = ', countHeader, countHeaderMpi
    write(*,*)

    do stnIdIndex = 1, numStnId
      call rpn_comm_allReduce(numObsStnId(stnIdIndex), countMpi, 1, 'mpi_integer', &
                              'mpi_sum', 'grid', ierr)
      write(*,'(a9,a,2i8)')  stnidList(stnIdIndex), ' :  ', numObsStnId(stnIdIndex), countMpi
    end do
    
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

    ! loop over all header indices of the specified family
    call obs_set_current_header_list(obsdat,trim(familyType))
    HEADER1: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER1

      if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) /= codtyp) cycle HEADER1

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
          end if
        end if

        ! count the number of accepted channels
        obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
        if ( .not. btest(obsFlag,8) .and. &
             .not. btest(obsFlag,9) .and. &
             .not. btest(obsFlag,11) ) then
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

    end do HEADER1

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

    ! modify flags, by setting bit 11 for those rejected, and count obs
    countHeader = 0
    numObsStnId(:) = 0
    call obs_set_current_header_list(obsdat,trim(familyType))
    HEADER2: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER2

      if (obs_headElem_i(obsdat, OBS_ITY, headerIndex) /= codtyp) cycle HEADER2

      if (.not. rejectThisHeader(headerIndex)) then

        stnid = obs_elem_c(obsdat,'STID',headerIndex)
        stnIdIndexFound = -1
        do stnIdIndex = 1, numStnId
          if ( stnidList(stnIdIndex) == stnid ) stnIdIndexFound = stnIdIndex
        end do
        if ( stnIdIndexFound == -1 ) then
          call utl_abort('thn_hyperByLatLonBoxes: Problem with stnId')
        end if
        numObsStnId(stnIdIndexFound) = numObsStnId(stnIdIndexFound) + 1
        countHeader = countHeader + 1

        cycle HEADER2
      end if

      ! loop over all body indices for this headerIndex
      call obs_set_current_body_list(obsdat, headerIndex)
      BODY2: do 
        bodyIndex = obs_getBodyIndex(obsdat)
        if (bodyIndex < 0) exit BODY2

        obsFlag = obs_bodyElem_i(obsdat, OBS_FLG, bodyIndex)
        call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, ibset(obsFlag,11))
      end do BODY2

    end do HEADER2

    call rpn_comm_allReduce(countHeader, countHeaderMpi, 1, 'mpi_integer', &
                            'mpi_sum', 'grid', ierr)

    write(*,*)
    write(*,'(a)')        ' == Output file == '

    write(*,*)
    write(*,'(a,2i8)')     ' Number of valid satellite profiles  = ', countHeader, countHeaderMpi
    write(*,*)

    do stnIdIndex = 1, numStnId
      call rpn_comm_allReduce(numObsStnId(stnIdIndex), countMpi, 1, 'mpi_integer', &
                              'mpi_sum', 'grid', ierr)
      write(*,'(a9,a,2i8)')  stnidList(stnIdIndex), ' :  ', numObsStnId(stnIdIndex), countMpi
    end do
    
    deallocate(rejectThisHeader)
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
    deallocate(stnIdIntMpi)
    deallocate(stnIdInt)

    write(*,*)
    write(*,*) 'thn_hyperByLatLonBoxes: Finished'
    write(*,*)

  end subroutine thn_hyperByLatLonBoxes

  !--------------------------------------------------------------------------
  ! thn_separation
  !--------------------------------------------------------------------------
  function thn_separation(xlon1,xlat1,xlon2,xlat2)
    ! :Purpose: Compute the separation distance for some thinning algorithms.

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


module satWind_mod
  ! MODULE satWind_mod (prefix='swd' category'=8. Low-level utilities and constants')
  !
  !:Purpose: utilities realted SW (AMV) data
  !
  use midasMpi_mod
  use utilities_mod

  implicit none
  save
  private
  
  ! Public procedures
  public :: swd_readSwqi

contains

  !--------------------------------------------------------------------------
  ! swd_read_swqi
  !--------------------------------------------------------------------------
  subroutine swd_readSwqi(SWname,QIvalue)
    !
    !:Purpose: read NAMSW block in the namelist
    !
    implicit none

    ! Arguments:
    character(len=*), allocatable, intent(out) :: SWname(:)
    character(len=*), allocatable, intent(out) :: QIvalue(:)

    ! Locals:
    integer :: ierr
    integer :: nsats, isat
    integer, parameter             :: maxSat = 99
    character(len=20), allocatable :: SWQIArray(:)

    ! Namelist variables
    character(len=20) :: SWQI(maxSat)

    namelist /NAMSW/ SWQI

    ! Defeault values for namelist variables
    SWQI(:)  = ''
    SWQI(1)  = 'METSAT7:qi1'
    SWQI(2)  = 'METSAT8:qi1'
    SWQI(3)  = 'METSAT9:qi1'
    SWQI(4)  = 'METSAT10:qi1'
    SWQI(5)  = 'METSAT11:qi1'
    SWQI(6)  = 'HMWARI-8:qi1'
    SWQI(7)  = 'HMWARI-9:qi1'
    SWQI(8)  = 'GOES13:qi1'
    SWQI(9)  = 'GOES15:qi1'
    SWQI(10) = 'GOES16:qi1'
    SWQI(11) = 'GOES17:qi1'
    SWQI(12) = 'GOES18:qi1'
    SWQI(13) = 'NOAA15:qi1'
    SWQI(14) = 'NOAA16:qi1'
    SWQI(15) = 'NOAA18:qi1'
    SWQI(16) = 'NOAA19:qi1'
    SWQI(17) = 'NOAA20:qi1'
    SWQI(18) = 'NOAA21:qi1'
    SWQI(19) = 'NPP:qi1'
    SWQI(20) = 'AQUA:qi1'
    SWQI(21) = 'TERRA:qi1'
    SWQI(22) = 'METOP-1:qi1'
    SWQI(23) = 'METOP-2:qi1'
    SWQI(24) = 'METOP-3:qi1'
    SWQI(25) = 'METOP1-3:qi1'
    SWQI(26) = 'GEO-POL:qi1'

    ! Read the namelist for SatWinds observations
    if (utl_isNamelistPresent('NAMSW','./flnml')) then
      call utl_tmg_start(181,'low-level--readNML')
      read(utl_flnml, nml=NAMSW, iostat=ierr)
      if (ierr /= 0) call utl_abort('thn_thinSatWinds: Error reading NAMSW namelist')
      call utl_tmg_stop(181)
    else
      write(*,*)
      write(*,*) 'thn_thinSatWinds: Namelist block NAMSW is missing in the namelist.'
      write(*,*) '                  The default value will be taken.'
    end if
    if (mmpi_myid == 0) write(*,nml=NAMSW)

    nsats = getNumSats(maxSat,SWQI)
    allocate(SWname(nsats))
    allocate(QIvalue(nsats))
    !call SplitString(nsats,SWQI,SWname,QIvalue)
    do isat = 1, nsats
      call utl_splitString(SWQI(isat),':',SWQIArray)
      SWname(isat) = SWQIArray(1)
      QIvalue(isat) = SWQIArray(2)
      deallocate(SWQIArray)
    end do

  end subroutine swd_readSwqi

  !--------------------------------------------------------------------------
  ! getNumSats
  !--------------------------------------------------------------------------
  integer function getNumSats(maxSat,vars)
    !
    !:Purpose: count the number of satellites, i.e. count the number of non ''
    !
    implicit none

    ! Arguments:
    integer,          intent(in) :: maxSat
    character(len=*), intent(in) :: vars(maxSat)

    ! Locals:
    integer                       :: varIndex

    getNumSats = 0

    do varIndex = 1, maxSat
      if (trim(vars(varIndex)) /= '') getNumSats = getNumSats + 1
    end do

  end function getNumSats

end module satWind_mod

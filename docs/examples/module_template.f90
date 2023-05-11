
module verticalCoord_mod
  ! MODULE verticalcoord (prefix='vco' category='7. Low-level data objects and utilities')
  !
  !:Purpose: Derived type and procedures related to the vertical levels including
  !          a pointer to the associated VGRID descriptor
  use mpi_mod
  use mathPhysConstants_mod
  use vGrid_Descriptors
  use utilities_mod
  implicit none
  private

  ! public derived type
  public :: struct_vco
  ! public procedures
  public :: vco_setupManual, vco_equal, vco_allocate
  ! public entities accessed through inheritance
  public :: vgd_get, vgd_levels, vgd_ok, vgd_dpidpis, vgd_write

  type struct_vco
    logical :: initialized=.false.
    integer :: nlev_T= 0
    integer :: nlev_M= 0
    integer :: ip1_sfc  ! ip1 value for the surface (Vcode=5005)
    integer, pointer :: ip1_T(:) => null()
    integer, pointer :: ip1_M(:) => null()  ! encoded IP1 levels (Thermo/Moment)
    type(vgrid_descriptor) :: vgrid
    character(len=8) :: setuptype
  end type struct_vco

  contains

  subroutine vco_allocate(vco)
    implicit none
    ! arguments
    type(struct_vco), pointer :: vco
    ! locals
    integer :: numLev, statusSum, status

    statusSum = 0

    numLev = vco_getNumLev(vco,'MM')
    allocate (vco%ip1_M(numLev),stat = status)
    statusSum = statusSum + status

    numLev = vco_getNumLev(vco,'TH')
    allocate (vco%ip1_T(numLev),stat = status)
    statusSum = statusSum + status

    if(statusSum /= 0 ) then
      call utl_abort('vco_allocate: problem with allocate in vco')
    end if

  end subroutine vco_allocate


  subroutine vco_setupManual(vco,ip1,numLev)
    implicit none
    ! arguments
    type(struct_vco), pointer :: vco
    integer, intent(in) :: numLev
    integer, intent(in) :: ip1(numlev)
    ! locals
    integer :: ip1Sfc
    character(len=10) :: blk_S

    write(*,*) 
    write(*,*) 'vco_setupManual: Creating an adhoc verticalgrid using'
    write(*,*) '                   number of level = ', numLev
    write(*,*) '                   ip1             = ', ip1

    if ( associated(vco) ) then
      call utl_abort('vco_setupManual: the supplied vco pointer is not null!')
    end if

    allocate(vco)

    vco%setupType = 'Manual'
 
    vco%nlev_T = numLev
    vco%nlev_M = numLev

    call vco_allocate(vco)

    vco%ip1_T(:) = ip1(:)
    vco%ip1_M(:) = ip1(:)

    ! determine IP1 of sfc (hyb=1.0)
    call convip(ip1Sfc, 1.0, 5, 2, blk_s, .false.)

    vco%ip1_sfc = ip1Sfc  ! ip1 value for the surface (Vcode=5005)

    vco%initialized = .true.

  end subroutine vco_setupManual


  function vco_equal(vco1,vco2) result(equal)
    implicit none
    ! arguments
    type(struct_vco), pointer :: vco1, vco2
    logical                   :: equal

    equal = .true.

    if ( trim(vco1%setupType) == 'fromFile' .and. trim(vco2%setupType) == 'fromFile' ) then
      equal = equal .and. (vco1%vgrid == vco2%vgrid)
      if (.not. equal) then
        write(*,*) 'vco_equal: vgrid not equal'
        return
      end if
    end if

    ! Even if vgrid defined, not enough just to compare vgrid, must compare everything
    equal = equal .and. (vco1%nlev_T == vco2%nlev_T)
    if (.not. equal) then
      write(*,*) 'vco_equal: nlev_T not equal', vco1%nlev_T, vco2%nlev_T
      return
    end if
    equal = equal .and. (vco1%nlev_M == vco2%nlev_M)
    if (.not. equal) then
      write(*,*) 'vco_equal: nlev_M not equal', vco1%nlev_M, vco2%nlev_M
      return
    end if
    equal = equal .and. all(vco1%ip1_T(:) == vco2%ip1_T(:))
    if (.not. equal) then
      write(*,*) 'vco_equal: ip1_T not equal'
      return
    end if
    equal = equal .and. all(vco1%ip1_M(:) == vco2%ip1_M(:))
    if (.not. equal) then
      write(*,*) 'vco_equal: ip1_M not equal'
      return
    end if
    equal = equal .and. (vco1%ip1_sfc == vco2%ip1_sfc)
    if (.not. equal) then
      write(*,*) 'vco_equal: ip1_sfc not equal'
      return
    end if

  end function vco_equal

end module verticalCoord_mod

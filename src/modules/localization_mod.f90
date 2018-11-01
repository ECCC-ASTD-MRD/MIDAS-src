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
!! module localization (prefix="loc")
!!
!! *Purpose*: Master module for the computation of localized 3D gridpoint amplitude 
!!            fields for each ensemble member from a given (1D) control vector
!!
!--------------------------------------------------------------------------
module localization_mod
  use mpi_mod
  use mpivar_mod
  use utilities_mod
  use localizationSpectral_mod
  use localizationFunction_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use ensembleStatevector_mod
  implicit none
  save
  private

  ! public procedures
  public :: loc_setup, loc_Lsqrt, loc_LsqrtAd, loc_finalize
  public :: loc_getLocInfo, loc_getNumLocActive
  public :: loc_reducetompilocal, loc_reducetompilocal_r4
  public :: loc_expandtompiglobal, loc_expandtompiglobal_r4

  public :: struct_loc

  type :: struct_loc
     logical             :: initialized = .false.
     integer             :: id
     integer             :: nEns
     integer             :: nEnsOverDimension
     integer             :: cvDim
     character(len=64)   :: locType
     type(struct_vco), pointer :: vco => null()
     type(struct_hco), pointer :: hco => null()
  end type struct_loc

  integer, parameter :: nMaxLoc = 10
  integer            :: nLocAlreadyAllocated = 0
  type(struct_loc), target :: loc(nMaxLoc)

  logical, parameter :: verbose = .false.

CONTAINS

!--------------------------------------------------------------------------
! loc_setup
!--------------------------------------------------------------------------
  subroutine loc_setup(hco_loc, vco_loc, nEns, pressureProfile, ntrunc, locType, &
                       locMode, horizLengthScale1, horizLengthScale2, vertLengthScale, cvDim_out, &
                       id_out)
    implicit none
  
    type(struct_hco), pointer, intent(in) :: hco_loc
    type(struct_vco), pointer, intent(in) :: vco_loc

    integer, intent(in) :: nEns
    integer, intent(in) :: nTrunc

    real(8), intent(in) :: pressureProfile(vco_loc%nLev_M)
    real(8), intent(in) :: horizLengthScale1
    real(8), intent(in) :: horizLengthScale2
    real(8), intent(in) :: vertLengthScale

    character(len=*), intent(in) :: locType, locMode

    integer, intent(out) :: cvDim_out
    integer, intent(out) :: id_out

    integer :: id, nEnsOverDimension

    call tmg_start(130,'LOC_SETUP')

    if (verbose) write(*,*) 'Entering loc_Setup'

    !
    !- 1.  ID allocation
    !
    nLocAlreadyAllocated = nLocAlreadyAllocated + 1
    if (nLocAlreadyAllocated <= nMaxLoc) then
       id = nLocAlreadyAllocated
       id_out = id
       write(*,*)
       write(*,*) "loc_setup: Setting localization id = ",id
    else
       call utl_abort('loc_setup: Too many localizations!!!')
    end if

    !
    !- 2.  Setup
    !
    loc(id)%locType = trim(locType)
    loc(id)%hco => hco_loc
    loc(id)%vco => vco_loc

    select case (trim(loc(id)%locType))
    case('spectral')
       if (mpi_myid == 0) write(*,*)
       if (mpi_myid == 0) write(*,*) 'loc_setup: LocType = ', trim(locType)
       call lsp_setup(hco_loc, nEns, vco_loc%nLev_M, pressureProfile, ntrunc, locType,& ! IN
                      locMode, horizLengthScale1, horizLengthScale2, vertLengthScale, & ! IN
                      cvDim_out, id, nEnsOverDimension)                                 ! OUT
    case default
       write(*,*)
       write(*,*) 'locType = ', trim(locType)
       call utl_abort('loc_setup: unknown locType')
    end select

    loc(id)%id    = id
    loc(id)%cvDim = cvDim_out
    loc(id)%nEnsOverDimension = nEnsOverDimension 

    !
    !- 3.  Ending
    !
    loc(id)%initialized = .true.

    call tmg_stop(130)

  end subroutine loc_setup

!--------------------------------------------------------------------------
! loc_Lsqrt
!--------------------------------------------------------------------------
  subroutine loc_Lsqrt(id, controlVector, ensAmplitude, stepIndex)
    implicit none

    integer, intent(in)  :: id, stepIndex
    real(8), intent(in)  :: controlVector(:)
    type(struct_ens)     :: ensAmplitude

    if (verbose) write(*,*) 'Entering loc_Lsqrt'
    call idcheck(id)

    select case (trim(loc(id)%locType))
    case('spectral')
       call lsp_Lsqrt(loc(id)%id, controlVector, & ! IN
                      ensAmplitude,              & ! OUT
                      stepIndex)                   ! IN
    case default
       call utl_abort('loc_Lsqrt: unknown locType')
    end select

  end subroutine loc_Lsqrt

!--------------------------------------------------------------------------
! loc_LsqrtAd
!--------------------------------------------------------------------------
  subroutine loc_LsqrtAd(id, ensAmplitude, controlVector, stepIndex)
    implicit none

    integer, intent(in)   :: id, stepIndex
    real(8), intent(out)  :: controlVector(:)
    type(struct_ens)      :: ensAmplitude

    if (verbose) write(*,*) 'Entering loc_LsqrtAd'
    call idcheck(id)

    select case (trim(loc(id)%locType))
    case('spectral')
       call lsp_LsqrtAd(loc(id)%id,   & ! IN
                        ensAmplitude, & ! INOUT
                        controlVector,& ! OUT
                        stepIndex )     ! IN
    case default
       call utl_abort('loc_LsqrtAd: unknown locType')
    end select

  end subroutine loc_LsqrtAd

!--------------------------------------------------------------------------
! loc_finalize
!--------------------------------------------------------------------------
  subroutine loc_finalize(id)
    implicit none

    integer, intent(in)  :: id

    if (verbose) write(*,*) 'Entering loc_finalize'
    call idcheck(id)

    select case (trim(loc(id)%locType))
    case('spectral')
       call lsp_finalize(loc(id)%id)
    case default
       call utl_abort('loc_finalize: unknown locType')
    end select

  end subroutine loc_finalize

!--------------------------------------------------------------------------
! loc_reduceToMPILocal
!--------------------------------------------------------------------------
  subroutine loc_reduceToMPILocal(id,cv_mpilocal,cv_mpiglobal)
    implicit none

    integer, intent(in) :: id

    real(8), intent(out) :: cv_mpilocal(:)
    real(8), intent(in)  :: cv_mpiglobal(:)

    if (verbose) write(*,*) 'Entering loc_reduceToMPILocal'
    call idcheck(id)

    select case (trim(loc(id)%locType))
    case('spectral')
       call lsp_reduceToMPILocal(loc(id)%id,     & ! IN
                                 cv_mpilocal,    & ! OUT
                                 cv_mpiglobal)     ! IN
    case default
       call utl_abort('loc_reduceToMPILocal: unknown locType')
    end select

 end subroutine loc_reduceToMPILocal

!--------------------------------------------------------------------------
! loc_reduceToMPILocal_r4
!--------------------------------------------------------------------------
  subroutine loc_reduceToMPILocal_r4(id,cv_mpilocal,cv_mpiglobal)
    implicit none
    integer, intent(in)  :: id
    real(4), intent(out) :: cv_mpilocal(:)
    real(4), intent(in)  :: cv_mpiglobal(:)

    if (verbose) write(*,*) 'Entering loc_reduceToMPILocal_r4'
    call idcheck(id)

    select case (trim(loc(id)%locType))
    case('spectral')
       call lsp_reduceToMPILocal_r4(loc(id)%id,    & ! IN
                                    cv_mpilocal,   & ! OUT
                                    cv_mpiglobal)    ! IN
    case default
       call utl_abort('loc_reduceToMPILocal_r4: unknown locType')
    end select

 end subroutine loc_reduceToMPILocal_r4

!--------------------------------------------------------------------------
! loc_expandToMPIGlobal
!--------------------------------------------------------------------------
  subroutine loc_expandToMPIGlobal(id,cv_mpilocal,cv_mpiglobal)
    implicit none

    integer, intent(in)  :: id
    real(8), intent(in)  :: cv_mpilocal(:)
    real(8), intent(out) :: cv_mpiglobal(:)

    if (verbose) write(*,*) 'Entering loc_expandToMPIGlobal'
    call idcheck(id)
    
    select case (trim(loc(id)%locType))
    case('spectral')
       call lsp_expandToMPIGlobal(loc(id)%id, cv_mpilocal,  & ! IN
                                  cv_mpiglobal)               ! OUT
    case default
       call utl_abort('loc_expandToMPIGlobal: unknown locType')
    end select

  end subroutine loc_expandToMPIGlobal

!--------------------------------------------------------------------------
! loc_expandToMPIGlobal_r4
!--------------------------------------------------------------------------
  subroutine loc_expandToMPIGlobal_r4(id,cv_mpilocal,cv_mpiglobal)
    implicit none

    integer, intent(in)  :: id
    real(4), intent(in)  :: cv_mpilocal(:)
    real(4), intent(out) :: cv_mpiglobal(:)

    if (verbose) write(*,*) 'Entering loc_expandToMPIGlobal_r4'
    call idcheck(id)
    
    select case (trim(loc(id)%locType))
    case('spectral')
       call lsp_expandToMPIGlobal_r4(loc(id)%id, cv_mpilocal, & ! IN
                                     cv_mpiglobal)              ! OUT
    case default
       call utl_abort('loc_expandToMPIGlobal_r4: unknown locType')
    end select

  end subroutine  loc_expandToMPIGlobal_r4

!--------------------------------------------------------------------------
!   IDCHECK
!--------------------------------------------------------------------------
  subroutine idcheck(id)
    implicit none

    integer, intent(in) :: id

    if ( .not. loc(id)%initialized) then
       write(*,*)
       write(*,*) "transform ID ", id
       call utl_abort('loc_IDCHECK: Unknown transform ID')
    end if

  end subroutine idcheck

!--------------------------------------------------------------------------
!   loc_getNumLocActive
!--------------------------------------------------------------------------
  function loc_getNumLocActive() result(numLocActive)
    implicit none

    integer :: numLocActive

    if (verbose) write(*,*) 'Entering loc_getNumLocActive'

    numLocActive = nLocAlreadyAllocated

  end function loc_getNumLocActive

!--------------------------------------------------------------------------
!   loc_getLocInfo
!--------------------------------------------------------------------------
  function loc_getLocInfo(id) result(locInfo_ptr)
    implicit none

    integer, intent(in) :: id
    type(struct_loc), pointer :: locInfo_ptr

    if (verbose) write(*,*) 'Entering loc_getLocInfo'
    call idcheck(id)

    locInfo_ptr => loc(id)

  end function loc_getLocInfo

end module localization_mod

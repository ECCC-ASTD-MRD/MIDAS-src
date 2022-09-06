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

module localization_mod
  ! MODULE localization_mod (prefix='loc' category='2. B and R matrices')
  !
  ! :Purpose: Master module for the computation of localized 3D gridpoint
  !           amplitude fields for each ensemble member from a given (1D)
  !           control vector
  !
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

  ! public structure
  public :: struct_loc
  ! public procedures
  public :: loc_setup, loc_Lsqrt, loc_LsqrtAd, loc_finalize
  public :: loc_reducetompilocal, loc_reducetompilocal_r4
  public :: loc_expandtompiglobal, loc_expandtompiglobal_r4

  type :: struct_loc
     logical             :: initialized = .false.
     type(struct_lsp), pointer :: lsp => null()
     integer             :: nEns
     integer             :: nEnsOverDimension
     integer             :: cvDim
     character(len=64)   :: locType
     type(struct_vco), pointer :: vco => null()
     type(struct_hco), pointer :: hco => null()
  end type struct_loc

  logical, parameter :: verbose = .false.

CONTAINS

  !--------------------------------------------------------------------------
  ! loc_setup
  !--------------------------------------------------------------------------
  subroutine loc_setup(loc, cvDim_out, hco_loc, vco_loc, nEns, pressureProfile, ntrunc, locType, &
                       locMode, horizLengthScale1, horizLengthScale2, vertLengthScale)
    implicit none

    type(struct_loc) :: loc
    integer, intent(out) :: cvDim_out

    type(struct_hco), pointer, intent(in) :: hco_loc
    type(struct_vco), pointer, intent(in) :: vco_loc

    integer, intent(in) :: nEns
    integer, intent(in) :: nTrunc

    real(8), intent(in) :: pressureProfile(:)
    real(8), intent(in) :: horizLengthScale1
    real(8), intent(in) :: horizLengthScale2
    real(8), intent(in) :: vertLengthScale

    character(len=*), intent(in) :: locType, locMode

    integer :: nEnsOverDimension, nLev

    if (verbose) write(*,*) 'Entering loc_Setup'

    !
    !- 2.  Setup
    !
    loc%locType = trim(locType)
    loc%hco => hco_loc
    loc%vco => vco_loc

    if ( loc%vco%Vcode == 5002 .or. loc%vco%Vcode == 5005 ) then
      if (loc%vco%nLev_M > 0) then
        nLev = loc%vco%nLev_M
      else
        nLev = loc%vco%nLev_T
      end if
    else !  vco_anl%Vcode == 0
      nLev = 1
    end if

    select case (trim(loc%locType))
    case('spectral')
       if (mpi_myid == 0) write(*,*)
       if (mpi_myid == 0) write(*,*) 'loc_setup: LocType = ', trim(locType)
       call lsp_setup(hco_loc, nEns, nLev, pressureProfile, ntrunc, locType,          & ! IN
                      locMode, horizLengthScale1, horizLengthScale2, vertLengthScale, & ! IN
                      cvDim_out, loc%lsp, nEnsOverDimension)                            ! OUT
    case default
       write(*,*)
       write(*,*) 'locType = ', trim(locType)
       call utl_abort('loc_setup: unknown locType')
    end select

    loc%cvDim = cvDim_out
    loc%nEnsOverDimension = nEnsOverDimension 

    !
    !- 3.  Ending
    !
    loc%initialized = .true.

  end subroutine loc_setup

  !--------------------------------------------------------------------------
  ! loc_Lsqrt
  !--------------------------------------------------------------------------
  subroutine loc_Lsqrt(loc, controlVector, ensAmplitude, stepIndex)
    implicit none

    type(struct_loc)     :: loc

    integer, intent(in)  :: stepIndex
    real(8), intent(in)  :: controlVector(:)
    type(struct_ens)     :: ensAmplitude

    if (verbose) write(*,*) 'Entering loc_Lsqrt'

    select case (trim(loc%locType))
    case('spectral')
      call lsp_Lsqrt(loc%lsp, controlVector, & ! IN
                     ensAmplitude,          & ! OUT
                     stepIndex)               ! IN
    case default
       call utl_abort('loc_Lsqrt: unknown locType')
    end select

  end subroutine loc_Lsqrt

  !--------------------------------------------------------------------------
  ! loc_LsqrtAd
  !--------------------------------------------------------------------------
  subroutine loc_LsqrtAd(loc, ensAmplitude, controlVector, stepIndex)
    implicit none

    type(struct_loc)      :: loc

    integer, intent(in)   :: stepIndex
    real(8), intent(out)  :: controlVector(:)
    type(struct_ens)      :: ensAmplitude

    if (verbose) write(*,*) 'Entering loc_LsqrtAd'

    select case (trim(loc%locType))
    case('spectral')
       call lsp_LsqrtAd(loc%lsp,       & ! IN
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
  subroutine loc_finalize(loc)
    implicit none

    type(struct_loc) :: loc

    if (verbose) write(*,*) 'Entering loc_finalize'

    select case (trim(loc%locType))
    case('spectral')
      call lsp_finalize(loc%lsp)
    case default
      call utl_abort('loc_finalize: unknown locType')
    end select

  end subroutine loc_finalize

  !--------------------------------------------------------------------------
  ! loc_reduceToMPILocal
  !--------------------------------------------------------------------------
  subroutine loc_reduceToMPILocal(loc,cv_mpilocal,cv_mpiglobal)
    implicit none

    type(struct_loc)     :: loc

    real(8), intent(out) :: cv_mpilocal(:)
    real(8), intent(in)  :: cv_mpiglobal(:)

    if (verbose) write(*,*) 'Entering loc_reduceToMPILocal'

    select case (trim(loc%locType))
    case('spectral')
      call lsp_reduceToMPILocal(loc%lsp,      & ! IN
                                cv_mpilocal, & ! OUT
                                cv_mpiglobal)  ! IN
    case default
      call utl_abort('loc_reduceToMPILocal: unknown locType')
    end select

  end subroutine loc_reduceToMPILocal
 
  !--------------------------------------------------------------------------
  ! loc_reduceToMPILocal_r4
  !--------------------------------------------------------------------------
  subroutine loc_reduceToMPILocal_r4(loc,cv_mpilocal,cv_mpiglobal)
    implicit none

    type(struct_loc)     :: loc

    real(4), intent(out) :: cv_mpilocal(:)
    real(4), intent(in)  :: cv_mpiglobal(:)

    if (verbose) write(*,*) 'Entering loc_reduceToMPILocal_r4'

    select case (trim(loc%locType))
    case('spectral')
      call lsp_reduceToMPILocal_r4(loc%lsp,       & ! IN
                                   cv_mpilocal,  & ! OUT
                                   cv_mpiglobal)   ! IN
    case default
      call utl_abort('loc_reduceToMPILocal_r4: unknown locType')
    end select

  end subroutine loc_reduceToMPILocal_r4
 
  !--------------------------------------------------------------------------
  ! loc_expandToMPIGlobal
  !--------------------------------------------------------------------------
  subroutine loc_expandToMPIGlobal(loc,cv_mpilocal,cv_mpiglobal)
    implicit none

    type(struct_loc)     :: loc

    real(8), intent(in)  :: cv_mpilocal(:)
    real(8), intent(out) :: cv_mpiglobal(:)

    if (verbose) write(*,*) 'Entering loc_expandToMPIGlobal'
    
    select case (trim(loc%locType))
    case('spectral')
      call lsp_expandToMPIGlobal(loc%lsp, cv_mpilocal,  & ! IN
                                 cv_mpiglobal)           ! OUT
    case default
      call utl_abort('loc_expandToMPIGlobal: unknown locType')
    end select

  end subroutine loc_expandToMPIGlobal

  !--------------------------------------------------------------------------
  ! loc_expandToMPIGlobal_r4
  !--------------------------------------------------------------------------
  subroutine loc_expandToMPIGlobal_r4(loc,cv_mpilocal,cv_mpiglobal)
    implicit none

    type(struct_loc)     :: loc

    real(4), intent(in)  :: cv_mpilocal(:)
    real(4), intent(out) :: cv_mpiglobal(:)

    if (verbose) write(*,*) 'Entering loc_expandToMPIGlobal_r4'
    
    select case (trim(loc%locType))
    case('spectral')
      call lsp_expandToMPIGlobal_r4(loc%lsp, cv_mpilocal, & ! IN
                                    cv_mpiglobal)          ! OUT
    case default
      call utl_abort('loc_expandToMPIGlobal_r4: unknown locType')
    end select

  end subroutine  loc_expandToMPIGlobal_r4

end module localization_mod

!-------------------------------------- LICENCE BEGIN ------------------------------------
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

module calcHeightAndPressure_mod
  ! MODULE czp_calcHeightAndPressure (prefix='czp' category='4. Data Object transformations')
  !
  ! :Purpose: Subroutines for computing height and/or pressure on statevectors
  !           and columns depending on the vgrid kind.
  !           Nonlinear, tangent-linear and adjoint versions of these
  !           transformations are included in separate subroutines.
  !           Depending on the vertical representation of the state or column,
  !           pressure or height values are either computed or retrieved using
  !           the [vgrid](https://gitlab.science.gc.ca/RPN-SI/vgrid) library.
  !           When computation is required (for instance to compute height on a
  !           GEM-P, represented on pressure coordinates), thermodynamical 
  !           variables are required, typically `P0`, `TT` and `HU`.
  !           Height and pressure values are obtained for both thermodynamical
  !           and momentum levels and labeled `Z_T` (`P_T`) and `Z_M` (`P_M`).
  !
  use codePrecision_mod
  use midasMpi_mod
  use mathPhysConstants_mod
  use earthConstants_mod
  use physicsFunctions_mod
  use verticalCoord_mod
  use gridstatevector_mod
  use columnData_mod
  use utilities_mod
  use message_mod
  use varnamelist_mod
  use gps_mod
  use HorizontalCoord_mod
  implicit none
  save
  private

  ! public procedures
  public :: czp_calcZandP_nl, czp_calcZandP_tl, czp_calcZandP_ad
  public :: czp_calcHeight_nl, czp_calcHeight_tl, czp_calcHeight_ad
  public :: czp_calcPressure_nl, czp_calcPressure_tl, czp_calcPressure_ad
  public :: czp_calcReturnHeight_gsv_nl, czp_calcReturnPressure_gsv_nl
  public :: czp_calcReturnHeight_col_nl, czp_calcReturnPressure_col_nl
  public :: czp_fetchProfile

  interface czp_calcZandP_nl
    module procedure calcZandP_gsv_nl
    module procedure calcZandP_col_nl
  end interface czp_calcZandP_nl
  interface czp_calcZandP_tl
    module procedure calcZandP_gsv_tl
    module procedure calcZandP_col_tl
  end interface czp_calcZandP_tl
  interface czp_calcZandP_ad
    module procedure calcZandP_gsv_ad
    module procedure calcZandP_col_ad
  end interface czp_calcZandP_ad

  interface czp_calcHeight_nl
    module procedure calcHeight_gsv_nl
    module procedure calcHeight_col_nl
  end interface czp_calcHeight_nl
  interface czp_calcHeight_tl
    module procedure calcHeight_gsv_tl
    module procedure calcHeight_col_tl
  end interface czp_calcHeight_tl
  interface czp_calcHeight_ad
    module procedure calcHeight_gsv_ad
    module procedure calcHeight_col_ad
  end interface czp_calcHeight_ad

  interface czp_calcPressure_nl
    module procedure calcPressure_gsv_nl
    module procedure calcPressure_col_nl
  end interface czp_calcPressure_nl
  interface czp_calcPressure_tl
    module procedure calcPressure_gsv_tl
    module procedure calcPressure_col_tl
  end interface czp_calcPressure_tl
  interface czp_calcPressure_ad
    module procedure calcPressure_gsv_ad
    module procedure calcPressure_col_ad
  end interface czp_calcPressure_ad

  ! private module variables
  real(8), allocatable :: coeff_M_TT_gsv(:,:,:,:), coeff_M_HU_gsv(:,:,:,:)
  real(8), allocatable :: coeff_T_TT_gsv(:,:,:),   coeff_T_HU_gsv(:,:,:)
  real(8), allocatable :: coeff_M_P0_delPM_gsv(:,:,:,:)
  real(8), allocatable :: coeff_M_P0_dP_delPT_gsv(:,:,:,:)
  real(8), allocatable :: coeff_M_P0_dP_delP0_gsv(:,:,:,:)
  real(8), allocatable :: coeff_T_P0_delP1_gsv(:,:,:)
  real(8), allocatable :: coeff_T_P0_dP_delPT_gsv(:,:,:)
  real(8), allocatable :: coeff_T_P0_dP_delP0_gsv(:,:,:)

  real(8), allocatable :: coeff_M_TT_col(:,:), coeff_M_HU_col(:,:)
  real(8), allocatable :: coeff_T_TT_col(:),   coeff_T_HU_col(:)
  real(8), allocatable :: coeff_M_P0_delPM_col(:,:)
  real(8), allocatable :: coeff_M_P0_dP_delPT_col(:,:)
  real(8), allocatable :: coeff_M_P0_dP_delP0_col(:,:)
  real(8), allocatable :: coeff_T_P0_delP1_col(:),   coeff_T_P0_dP_delPT_col(:)
  real(8), allocatable :: coeff_T_P0_dP_delP0_col(:)

contains
  !---------------------------------------------------------------------
  ! subroutines operating on struct_gsv
  !---------------------------------------------------------------------

  !---------------------------------------------------------
  ! calcZandP_gsv_nl
  !---------------------------------------------------------
  subroutine calcZandP_gsv_nl(statevector)
    !
    ! :Purpose: Pressure and height computation on the grid in proper order
    !           depending on the vgrid kind.
    !           Depending on the vcode, the routine will check the existence of
    !           P_* (vcode=500x) or Z_* (vcode=2100x) first and proceed with
    !           pressure (height) computation.  Then, if the other variables are
    !           also present, it will secondly compute height (pressure).
    !           Hence if only P_* (Z_*) is present, only these are computed.
    !           If the first variable P_* (Z_*) is not present, nothing is done.
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout) :: statevector      ! statevector that will contain the Z_*/P_* fields

    ! Locals
    integer                   :: Vcode

    call msg('calcZandP_gsv_nl (czp)', 'START', verb_opt=2)

    Vcode = gsv_getVco(statevector)%vcode

    if (Vcode == 5002 .or. Vcode == 5005) then
      ! if P_T, P_M not allocated : do nothing
      if (gsv_varExist(statevector, 'P_*')) then
        call calcPressure_gsv_nl(statevector)
        if (gsv_varExist(statevector, 'Z_*')) then
          call calcHeight_gsv_nl(statevector)
        end if
      end if
    else if (Vcode == 21001) then
      ! if Z_T, Z_M not allocated : do nothing
      if (gsv_varExist(statevector, 'Z_*')) then
          call calcHeight_gsv_nl(statevector)
        if (gsv_varExist(statevector, 'P_*')) then
          call calcPressure_gsv_nl(statevector)
        end if
      end if
    end if

    call msg('calcZandP_gsv_nl (czp)', 'END', verb_opt=2)
  end subroutine calcZandP_gsv_nl

  !---------------------------------------------------------
  ! calcZandP_gsv_tl
  !---------------------------------------------------------
  subroutine calcZandP_gsv_tl(statevector, statevectorRef)
    !
    ! :Purpose: Pressure and height increment computation on the grid in proper
    !           order depending on the vgrid kind.
    !           Depending on the vcode, the routine will check the existence of
    !           P_* (vcode=500x) or Z_* (vcode=2100x) first and proceed with
    !           pressure (height) computation.  Then, if the other variables are
    !           also present, it will secondly compute height (pressure).
    !           Hence if only P_* (Z_*) is present, only these are computed.
    !           If the first variable P_* (Z_*) is not present, nothing is done.
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout) :: statevector      ! statevector that will contain the Z_*/P_* increments
    type(struct_gsv), intent(in)    :: statevectorRef   ! statevector containing needed reference fields

    ! Locals
    type(struct_vco), pointer :: vco
    integer                   :: Vcode


    call msg('calcZandP_gsv_tl (czp)', 'START', verb_opt=2)

    vco => gsv_getVco(statevector)
    Vcode = vco%vcode

    if (Vcode == 5002 .or. Vcode == 5005) then
      ! if P_T, P_M not allocated : do nothing
      if (gsv_varExist(statevector, 'P_*')) then

        if ( .not. gsv_containsNonZeroValues(stateVectorRef) ) then
          call utl_abort('calcZandP_gsv_tl: stateVectorRef not initialized')
        end if
        call calcPressure_gsv_tl(statevector, statevectorRef)

        if (gsv_varExist(statevector, 'Z_*')) then
          call calcHeight_gsv_tl(statevector, statevectorRef)
        end if

      end if
    else if (Vcode == 21001) then
      ! if Z_T, Z_M not allocated : do nothing
      if (gsv_varExist(statevector, 'Z_*')) then

        if ( .not. gsv_containsNonZeroValues(stateVectorRef) ) then
          call utl_abort('calcZandP_gsv_tl: stateVectorRef not initialized')
        end if
        call calcHeight_gsv_tl(statevector, statevectorRef)

        if (gsv_varExist(statevector, 'P_*')) then
          call calcPressure_gsv_tl(statevector, statevectorRef)
        end if

      end if
    end if

    call msg('calcZandP_gsv_tl (czp)', 'END', verb_opt=2)
  end subroutine calcZandP_gsv_tl

  !---------------------------------------------------------
  ! calcZandP_gsv_ad
  !---------------------------------------------------------
  subroutine calcZandP_gsv_ad(statevector, statevectorRef)
    !
    ! :Purpose: Pressure and height increment adjoint computation on the grid
    !           in proper order depending on the vgrid kind
    !           Depending on the vcode, the routine will check the existence of
    !           Z_* (vcode=500x) or P_* (vcode=2100x) first and proceed with
    !           height (pressure) adjoint computation.  Then, if the other
    !           variables are also present, it will secondly proceed with
    !           adjoint computation of pressure (height).
    !           Hence if only Z_* (P_*) is present, only these are computed.
    !           If the first variable Z_* (P_*) is not present, nothing is done.
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout) :: statevector      ! statevector that will contain the Z_*/P_* increments
    type(struct_gsv), intent(in)    :: statevectorRef   ! statevector containing needed reference fields

    ! Locals
    type(struct_vco), pointer :: vco
    integer                   :: Vcode


    call msg('calcZandP_gsv_ad (czp)', 'START', verb_opt=2)

    vco => gsv_getVco(statevector)
    Vcode = vco%vcode

    if (Vcode == 5002 .or. Vcode == 5005) then
      ! if Z_T, Z_M not allocated : do nothing
      if (gsv_varExist(statevector, 'Z_*')) then

        if ( .not. gsv_containsNonZeroValues(stateVectorRef) ) then
          call utl_abort('calcZandP_gsv_ad: stateVectorRef not initialized')
        end if
        call calcHeight_gsv_ad(statevector, statevectorRef)

        if (gsv_varExist(statevector, 'P_*')) then
          call calcPressure_gsv_ad(statevector, statevectorRef)
        end if

      end if
    else if (Vcode == 21001) then
      ! if P_T, P_M not allocated : do nothing
      if (gsv_varExist(statevector, 'P_*')) then

        if ( .not. gsv_containsNonZeroValues(stateVectorRef) ) then
          call utl_abort('calcZandP_gsv_ad: stateVectorRef not initialized')
        end if
        call calcPressure_gsv_ad(statevector, statevectorRef)

        if (gsv_varExist(statevector, 'Z_*')) then
          call calcHeight_gsv_ad(statevector, statevectorRef)
        end if

      end if
    end if

    call msg('calcZandP_gsv_ad (czp)', 'END', verb_opt=2)
  end subroutine calcZandP_gsv_ad

  !---------------------------------------------------------
  ! calcHeight_gsv_nl
  !---------------------------------------------------------
  subroutine calcHeight_gsv_nl(statevector)
    !
    ! :Purpose: Compute or retrieve heights and store values in statevector.
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout) :: statevector

    ! Locals
    integer :: Vcode

    real(4), pointer :: ptr_PT_r4(:,:,:,:), ptr_PM_r4(:,:,:,:)
    real(8), pointer :: ptr_PT_r8(:,:,:,:), ptr_PM_r8(:,:,:,:)
    real(4), pointer :: ptr_ZT_r4(:,:,:,:), ptr_ZM_r4(:,:,:,:)
    real(8), pointer :: ptr_ZT_r8(:,:,:,:), ptr_ZM_r8(:,:,:,:)

    call utl_tmg_start(172,'low-level--czp_calcHeight_nl')
    call msg('calcHeight_gsv_nl (czp)', 'START', verb_opt=2)

    Vcode = gsv_getVco(statevector)%vcode
    if (Vcode == 5005 .or. Vcode == 5002) then
      if ( gsv_getDataKind(statevector) == 4 ) then
        call gsv_getField(statevector, ptr_PT_r4, 'P_T')
        call gsv_getField(statevector, ptr_PM_r4, 'P_M')
        call gsv_getField(statevector, ptr_ZT_r4, 'Z_T')
        call gsv_getField(statevector, ptr_ZM_r4, 'Z_M')
        call calcHeight_gsv_nl_vcode500x( statevector, &
                                          PTin_r4_opt=ptr_PT_r4, &
                                          PMin_r4_opt=ptr_PM_r4, &
                                          ZTout_r4_opt=ptr_ZT_r4, &
                                          ZMout_r4_opt=ptr_ZM_r4)
      else
        call gsv_getField(statevector, ptr_PT_r8, 'P_T')
        call gsv_getField(statevector, ptr_PM_r8, 'P_M')
        call gsv_getField(statevector, ptr_ZT_r8, 'Z_T')
        call gsv_getField(statevector, ptr_ZM_r8, 'Z_M')
        call calcHeight_gsv_nl_vcode500x( statevector, &
                                          PTin_r8_opt=ptr_PT_r8, &
                                          PMin_r8_opt=ptr_PM_r8, &
                                          ZTout_r8_opt=ptr_ZT_r8, &
                                          ZMout_r8_opt=ptr_ZM_r8)
      end if

    else if (Vcode == 21001) then
      ! Development notes (@mad001) 
      !   probably some some gsv_varExist(statevector,.) needed for GEM-H
      if ( gsv_getDataKind(statevector) == 4 ) then
        call gsv_getField(statevector, ptr_ZT_r4, 'Z_T')
        call gsv_getField(statevector, ptr_ZM_r4, 'Z_M')
        call calcHeight_gsv_nl_vcode2100x_r4(statevector, ptr_ZT_r4, ptr_ZM_r4)
      else
        call gsv_getField(statevector, ptr_ZT_r8, 'Z_T')
        call gsv_getField(statevector, ptr_ZM_r8, 'Z_M')
        call calcHeight_gsv_nl_vcode2100x_r8(statevector, ptr_ZT_r8, ptr_ZM_r8)
      end if
    end if

    if ( gsv_getDataKind(statevector) == 4 ) then
      call msg('calcHeight_gsv_nl (czp)', &
             new_line('')//'Z_M = '&
           //str(ptr_ZM_r4(statevector%myLonBeg,statevector%myLatBeg,:,1), vertical_opt=.false.) &
           //new_line('')//'Z_T = '&
           //str(ptr_ZT_r4( statevector%myLonBeg,statevector%myLatBeg,:,1), vertical_opt=.false.), &
           verb_opt=2)
    else
      call msg('calcHeight_gsv_nl (czp)', &
             new_line('')//'Z_M = '&
           //str(ptr_ZM_r8(statevector%myLonBeg,statevector%myLatBeg,:,1), vertical_opt=.false.) &
           //new_line('')//'Z_T = '&
           //str(ptr_ZT_r8( statevector%myLonBeg,statevector%myLatBeg,:,1), vertical_opt=.false.), &
           verb_opt=2)
    end if

    call msg('calcHeight_gsv_nl (czp)', 'END', verb_opt=2)
    call utl_tmg_stop(172)
  end subroutine calcHeight_gsv_nl
  
  !---------------------------------------------------------
  ! czp_calcReturnHeight_gsv_nl
  !---------------------------------------------------------
  subroutine czp_calcReturnHeight_gsv_nl( statevector, &
                                          PTin_r4_opt, PMin_r4_opt, &
                                          PTin_r8_opt, PMin_r8_opt, &
                                          ZTout_r4_opt, ZMout_r4_opt, &
                                          ZTout_r8_opt, ZMout_r8_opt)
    !
    ! :Purpose: Compute or retrieve heights and return values in pointer arguments.
    !           Proceeds to vcode dispatching.
    !
    implicit none

    ! Arguments
    type(struct_gsv),           intent(in)    :: statevector
    real(4), optional, pointer, intent(in)    :: PTin_r4_opt(:,:,:,:)
    real(4), optional, pointer, intent(in)    :: PMin_r4_opt(:,:,:,:)
    real(8), optional, pointer, intent(in)    :: PTin_r8_opt(:,:,:,:)
    real(8), optional, pointer, intent(in)    :: PMin_r8_opt(:,:,:,:)
    real(4), optional, pointer, intent(inout) :: ZTout_r4_opt(:,:,:,:)
    real(4), optional, pointer, intent(inout) :: ZMout_r4_opt(:,:,:,:)
    real(8), optional, pointer, intent(inout) :: ZTout_r8_opt(:,:,:,:)
    real(8), optional, pointer, intent(inout) :: ZMout_r8_opt(:,:,:,:)

    ! Locals
    integer :: Vcode

    call utl_tmg_start(172,'low-level--czp_calcHeight_nl')
    call msg('czp_calcReturnHeight_gsv_nl', 'START', verb_opt=2)

    Vcode = gsv_getVco(statevector)%vcode
    if (Vcode == 5005 .or. Vcode == 5002) then
      if ( gsv_getDataKind(statevector) == 4 ) then
        if ( .not. (present(PTin_r4_opt) .and. present(PMin_r4_opt))) then
          call utl_abort('czp_calcReturnHeight_gsv_nl: dataKind=4: P{T,M}out_r4_opt expected')
        end if
        if ( .not. (present(ZTout_r4_opt) .and. present(ZMout_r4_opt))) then
          call utl_abort('czp_calcReturnHeight_gsv_nl: dataKind=4: Z{T,M}out_r4_opt expected')
        end if
        call calcHeight_gsv_nl_vcode500x( statevector, &
                                          PTin_r4_opt=PTin_r4_opt, & 
                                          PMin_r4_opt=PMin_r4_opt, & 
                                          ZTout_r4_opt=ZTout_r4_opt, &
                                          ZMout_r4_opt=ZMout_r4_opt)
      else ! datakind = 8
        if ( .not. (present(PTin_r8_opt) .and. present(PMin_r8_opt))) then
          call utl_abort('czp_calcReturnHeight_gsv_nl: dataKind=8: P{T,M}out_r8_opt expected')
        end if
        if ( .not. (present(ZTout_r8_opt) .and. present(ZMout_r8_opt))) then
          call utl_abort('czp_calcReturnHeight_gsv_nl: dataKind=8: Z{T,M}out_r8_opt expected')
        end if
        call calcHeight_gsv_nl_vcode500x( statevector, &
                                          PTin_r8_opt=PTin_r8_opt, &
                                          PMin_r8_opt=PMin_r8_opt, &
                                          ZTout_r8_opt=ZTout_r8_opt, &
                                          ZMout_r8_opt=ZMout_r8_opt)
      end if

    else if (Vcode == 21001) then
      if ( gsv_getDataKind(statevector) == 4 ) then
        if ( .not. (present(ZTout_r4_opt) .and. present(ZMout_r4_opt))) then
          call utl_abort('czp_calcReturnHeight_gsv_nl: dataKind=4: Z{T,M}_r4 expected')
        end if
        call calcHeight_gsv_nl_vcode2100x_r4(statevector, ZTout_r4_opt, ZMout_r4_opt)
      else
        if ( .not. (present(ZTout_r8_opt) .and. present(ZMout_r8_opt))) then
          call utl_abort('czp_calcReturnHeight_gsv_nl: dataKind=4: Z{T,M}_r4 expected')
        end if
        call calcHeight_gsv_nl_vcode2100x_r8(statevector, ZTout_r8_opt, ZMout_r8_opt)
      end if
    end if

    call msg('czp_calcReturnHeight_gsv_nl', 'END', verb_opt=2)
    call utl_tmg_stop(172) 
  end subroutine czp_calcReturnHeight_gsv_nl

  !---------------------------------------------------------
  ! calcHeight_gsv_nl_vcode2100x_r4
  !---------------------------------------------------------
  subroutine calcHeight_gsv_nl_vcode2100x_r4(statevector, Z_T, Z_M)
    !
    ! :Purpose: Retrieve heights for GEM-H statevector, return height values 
    !           in pointer arguments.
    !           real(4) version
    !
    implicit none

    ! Arguments
    type(struct_gsv),  intent(in)    :: statevector
    real(4), pointer,  intent(inout) :: Z_T(:,:,:,:)
    real(4), pointer,  intent(inout) :: Z_M(:,:,:,:)

    ! Locals
    integer ::  numStep, stepIndex, status
    real(kind=8), pointer       :: Hsfc(:,:)
    real(kind=4), allocatable   :: Hsfc4(:,:)
    real(kind=4), pointer       :: GZHeight_out(:,:,:)

    call msg('calcHeight_gsv_nl_vcode2100x_r4 (czp)', 'START', verb_opt=4)

    if ( .not. gsv_varExist(statevector,'Z_*')) then
      call utl_abort('calcHeight_gsv_nl_vcode2100x_r4 (czp): Z_T/Z_M do not exist in statevector!')
    end if

    allocate(Hsfc4( statevector%myLonBeg:statevector%myLonEnd, &
                    statevector%myLatBeg:statevector%myLatEnd))
    Hsfc => gsv_getHeightSfc(statevector)
    Hsfc4 = real(Hsfc,4)

    numStep = statevector%numStep

    do stepIndex = 1, numStep

      ! Z_M
      nullify(GZHeight_out)
      status = vgd_levels(statevector%vco%vgrid, &
                          ip1_list=statevector%vco%ip1_M, &
                          levels=GZHeight_out, &
                          sfc_field=Hsfc4, &
                          in_log=.false.)
      if( status .ne. VGD_OK ) then
          call utl_abort('calcHeight_gsv_nl_vcode2100x_r4 (czp): ERROR with vgd_levels')
      end if
      Z_M(:,:,:,stepIndex) = gz2alt_r4(statevector, GZHeight_out)
      deallocate(GZHeight_out)

      ! Z_T
      nullify(GZHeight_out)
      status = vgd_levels(statevector%vco%vgrid, &
                          ip1_list=statevector%vco%ip1_T, &
                          levels=GZHeight_out, &
                          sfc_field=Hsfc4, &
                          in_log=.false.)
      if( status .ne. VGD_OK ) then
          call utl_abort('calcHeight_gsv_nl_vcode2100x_r4 (czp): ERROR with vgd_levels')
      end if
      Z_T(:,:,:,stepIndex) = gz2alt_r4(statevector, GZHeight_out)
      deallocate(GZHeight_out)

    end do
    deallocate(Hsfc4)

    call msg('calcHeight_gsv_nl_vcode2100x_r4 (czp)', 'END', verb_opt=4)
  end subroutine calcHeight_gsv_nl_vcode2100x_r4

  !---------------------------------------------------------
  ! gz2alt_r4
  !---------------------------------------------------------
  function gz2alt_r4(statevector, gzHeight) result(alt)
    !
    ! :Purpose: Iterative conversion of geopotential height to geometric
    !           altitude.  (solution proposed by J. Aparicio)
    !           real(4) version.
    !
    implicit none

    ! Arguments
    real(kind=4), allocatable           :: alt(:,:,:)
    type(struct_gsv),      intent(in)   :: statevector
    real(kind=4), pointer, intent(in)   :: gzHeight(:,:,:)

    ! Locals
    integer                             :: nLon, nLat, nLev
    type(struct_hco), pointer           :: hco
    real(kind=4)                        :: latitude
    real(kind=4)                        :: gzH, b1, b2, A2, A3
    integer                             :: lonIndex, latIndex, lvlIndex

    ! gzHeight comes from external `vgd_levels` which does not know the
    ! mpi shifted indexes
    nLon = ubound(gzHeight, 1)
    nLat = ubound(gzHeight, 2)
    nLev = ubound(gzHeight, 3)
    allocate(alt(nLon, nLat, nLev))

    hco => gsv_getHco(statevector)

    do lonIndex = 1, nLon
      do latIndex = 1, nLat
        do lvlIndex = 1, nLev
          ! explicit shift of indexes
          latitude = hco%lat2d_4( lonIndex+statevector%myLonBeg-1,&
                                  latIndex+statevector%myLatBeg-1)
          gzH = gzHeight(lonIndex, latIndex, lvlIndex)
          ! gzH(alt) = g0 * (1 + b1*alt + b2*alt**2)
          b1 = -2.0/ec_wgs_a*(1.0+ec_wgs_f+ec_wgs_m-2*ec_wgs_f*latitude**2)
          b2 = 3.0/ec_wgs_a**2
          ! reversed series coefficients (Abramowitz and Stegun 3.6.25)
          A2 = -b1/2.0
          A3 = b1**2/2.0 - b2/3.0
          alt(lonIndex, latIndex, lvlIndex) = gzH + A2*gzH**2 + A3*gzH**3
        end do
      end do
    end do

  end function gz2alt_r4

  !---------------------------------------------------------
  ! calcHeight_gsv_nl_vcode2100x_r8
  !---------------------------------------------------------
  subroutine calcHeight_gsv_nl_vcode2100x_r8(statevector, Z_T, Z_M)
    !
    ! :Purpose: Retrieve heights for GEM-H statevector, return height values 
    !           in pointer arguments.
    !           real(8) version
    !
    implicit none

    ! Arguments
    type(struct_gsv),  intent(in)    :: statevector
    real(8), pointer,  intent(inout) :: Z_T(:,:,:,:)
    real(8), pointer,  intent(inout) :: Z_M(:,:,:,:)

    ! Locals
    integer ::  numStep, stepIndex, status
    real(kind=8), pointer   :: Hsfc(:,:), GZHeight_out(:,:,:)

    call msg('calcHeight_gsv_nl_vcode2100x_r8 (czp)', 'START', verb_opt=4)

    Hsfc => gsv_getHeightSfc(statevector)
    numStep = statevector%numStep

    do stepIndex = 1, numStep

      ! Z_M
      nullify(GZHeight_out)
      status = vgd_levels(statevector%vco%vgrid, &
                          ip1_list=statevector%vco%ip1_M, &
                          levels=GZHeight_out, &
                          sfc_field=Hsfc, &
                          in_log=.false.)
      if( status .ne. VGD_OK ) then
          call utl_abort('calcHeight_gsv_nl_vcode2100x_r8 (czp): ERROR with vgd_levels')
      end if
      Z_M(:,:,:,stepIndex) = gz2alt_r8(statevector, GZHeight_out)
      deallocate(GZHeight_out)

      ! Z_T
      nullify(GZHeight_out)
      status = vgd_levels(statevector%vco%vgrid, &
                          ip1_list=statevector%vco%ip1_T, &
                          levels=GZHeight_out, &
                          sfc_field=Hsfc, &
                          in_log=.false.)
      if( status .ne. VGD_OK ) then
          call utl_abort('calcHeight_gsv_nl_vcode2100x_r8 (czp): ERROR with vgd_levels')
      end if
      Z_T(:,:,:,stepIndex) = gz2alt_r8(statevector, GZHeight_out)
      deallocate(GZHeight_out)
    end do

    call msg('calcHeight_gsv_nl_vcode2100x_r8 (czp)', 'END', verb_opt=4)
  end subroutine calcHeight_gsv_nl_vcode2100x_r8

  !---------------------------------------------------------
  ! gz2alt_r8
  !---------------------------------------------------------
  function gz2alt_r8(statevector, gzHeight) result(alt)
    !
    ! :Purpose: Iterative conversion of geopotential height to geometric
    !           altitude.  (solution proposed by J. Aparicio)
    !           real(8) version.
    !
    implicit none

    ! Arguments
    real(kind=8), allocatable           :: alt(:,:,:)
    type(struct_gsv),      intent(in)   :: statevector
    real(kind=8), pointer, intent(in)   :: gzHeight(:,:,:)

    ! Locals
    integer                             :: nLon, nLat, nLev
    type(struct_hco), pointer           :: hco
    real(kind=8)                        :: latitude
    real(kind=8)                        :: gzH, b1, b2, A2, A3
    integer                             :: lonIndex, latIndex, lvlIndex

    ! gzHeight comes from external `vgd_levels` which does not know the
    ! mpi shifted indexes
    nLon = ubound(gzHeight, 1)
    nLat = ubound(gzHeight, 2)
    nLev = ubound(gzHeight, 3)
    allocate(alt(nLon, nLat, nLev))

    hco => gsv_getHco(statevector)

    do lonIndex = 1, nLon
      do latIndex = 1, nLat
        do lvlIndex = 1, nLev
          ! explicit shift of indexes
          latitude = hco%lat2d_4( lonIndex+statevector%myLonBeg-1,&
                                  latIndex+statevector%myLatBeg-1)
          gzH = gzHeight(lonIndex, latIndex, lvlIndex)
          ! gzH(alt) = g0 * (1 + b1*alt + b2*alt**2)
          b1 = -2.0D0/ec_wgs_a*(1.0D0+ec_wgs_f+ec_wgs_m-2*ec_wgs_f*latitude**2)
          b2 = 3.0D0/ec_wgs_a**2
          ! reversed series coefficients (Abramowitz and Stegun 3.6.25)
          A2 = -b1/2.0D0
          A3 = b1**2/2.0D0 - b2/3.0D0
          alt(lonIndex, latIndex, lvlIndex) = gzH + A2*gzH**2 + A3*gzH**3
        end do
      end do
    end do

  end function gz2alt_r8

  !---------------------------------------------------------
  ! calcHeight_gsv_nl_vcode500x
  !---------------------------------------------------------
  subroutine calcHeight_gsv_nl_vcode500x( statevector, &
                                          PTin_r4_opt, PMin_r4_opt, &
                                          PTin_r8_opt, PMin_r8_opt, &
                                          ZTout_r4_opt, ZMout_r4_opt, &
                                          ZTout_r8_opt, ZMout_r8_opt)
    !
    ! :Purpose: Compute heights for GEM-P statevector, return height values 
    !           in pointer arguments.
    !           Assumptions:
    !           1) nlev_T = nlev_M+1 (for vcode=5002)
    !           2) alt_T(nlev_T) = alt_M(nlev_M), both at the surface
    !           3) a thermo level exists at the top, higher than the highest
    !              momentum level
    !           4) the placement of the thermo levels means that alt_T is the
    !              average of 2 nearest alt_M (according to Ron and Claude)
    !
    implicit none

    ! Arguments
    type(struct_gsv),           intent(in)    :: statevector
    real(4), pointer, optional, intent(in)    :: PTin_r4_opt(:,:,:,:)
    real(4), pointer, optional, intent(in)    :: PMin_r4_opt(:,:,:,:)
    real(8), pointer, optional, intent(in)    :: PTin_r8_opt(:,:,:,:)
    real(8), pointer, optional, intent(in)    :: PMin_r8_opt(:,:,:,:)
    real(4), pointer, optional, intent(inout) :: ZTout_r4_opt(:,:,:,:)
    real(4), pointer, optional, intent(inout) :: ZMout_r4_opt(:,:,:,:)
    real(8), pointer, optional, intent(inout) :: ZTout_r8_opt(:,:,:,:)
    real(8), pointer, optional, intent(inout) :: ZMout_r8_opt(:,:,:,:)

    ! Locals
    integer ::  lev_M,lev_T,nlev_M,nlev_T,status,Vcode
    integer ::  numStep, stepIndex, latIndex,lonIndex
    real(4) ::  lat_4, heightSfcOffset_T_r4, heightSfcOffset_M_r4
    real(8) ::  delThick, ratioP
    real(8) ::  ScaleFactorBottom, ScaleFactorTop
    real(8) ::  P_M, P_M1, P_Mm1, P_T
    real(8) ::  hu, tt, Pr, cmp, h0, Rgh, P0, dh, rMt
    real(8) ::  sLat, cLat, lat_8

    real(8), allocatable :: tv(:), height_T(:), height_M(:)

    real(4), pointer     :: height_T_ptr_r4(:,:,:,:)
    real(4), pointer     :: height_M_ptr_r4(:,:,:,:)
    real(4), pointer     :: hu_ptr_r4(:,:,:,:),tt_ptr_r4(:,:,:,:)
    real(4), pointer     :: P_T_ptr_r4(:,:,:,:),P_M_ptr_r4(:,:,:,:)
    real(4), pointer     :: P0_ptr_r4(:,:,:,:)

    real(8), pointer     :: height_T_ptr_r8(:,:,:,:)
    real(8), pointer     :: height_M_ptr_r8(:,:,:,:)
    real(8), pointer     :: P_T_ptr_r8(:,:,:,:),P_M_ptr_r8(:,:,:,:)
    real(8), pointer     :: P0_ptr_r8(:,:,:,:)
    real(8), pointer     :: hu_ptr_r8(:,:,:,:),tt_ptr_r8(:,:,:,:)
    real(8), pointer     :: HeightSfc_ptr_r8(:,:)

    call msg('calcHeight_gsv_nl_vcode500x (czp)', 'START', verb_opt=4)

    nlev_T = gsv_getNumLev(statevector,'TH')
    nlev_M = gsv_getNumLev(statevector,'MM')
    Vcode = gsv_getVco(statevector)%vcode
    numStep = statevector%numStep

    allocate(tv(nlev_T))

    if (Vcode == 5002 .and. nlev_T /= nlev_M+1) then
      call utl_abort('calcHeight_gsv_nl_vcode500x (czp): nlev_T is not equal to nlev_M+1!')
    end if
    if (Vcode == 5005 .and. nlev_T /= nlev_M) then
      call utl_abort('calcHeight_gsv_nl_vcode500x (czp): nlev_T is not equal to nlev_M!')
    end if

    if (Vcode == 5005) then
      status = vgd_get( statevector%vco%vgrid, &
                        key='DHM - height of the diagnostic level (m)', &
                        value=heightSfcOffset_M_r4)
      status = vgd_get( statevector%vco%vgrid, &
                        key='DHT - height of the diagnostic level (t)', &
                        value=heightSfcOffset_T_r4)
      call msg('calcHeight_gsv_nl_vcode500x (czp)', &
           'height offset for near-sfc momentum level is:'//str(heightSfcOffset_M_r4)//' meters'&
           //new_line('')//'height offset for near-sfc thermo level is:'//str(heightSfcOffset_T_r4)//' meters', &
           verb_opt=2, mpiAll_opt=.false.)
      if ( .not.statevector%addHeightSfcOffset ) then
        call msg('calcHeight_gsv_nl_vcode500x (czp)', new_line('') &
             //'--------------------------------------------------------------------------'//new_line('')&
             //'BUT HEIGHT OFFSET REMOVED FOR DIAGNOSTIC LEVELS FOR BACKWARD COMPATIBILITY'//new_line('')&
             //'--------------------------------------------------------------------------', &
             verb_opt=2, mpiAll_opt=.false.)
      end if
    end if

    allocate(height_T(nlev_T))
    allocate(height_M(nlev_M))

    if ( statevector%dataKind == 4 ) then
      if ( .not. (present(PTin_r4_opt) .and. present(PMin_r4_opt))) then
        call utl_abort('calcHeight_gsv_nl_vcode500x (czp): dataKind=4: P{T,M}in_r4_opt expected')
      end if
      P_T_ptr_r4 => PTin_r4_opt
      P_M_ptr_r4 => PMin_r4_opt

      if ( .not. (present(ZTout_r4_opt) .and. present(ZMout_r4_opt))) then
        call utl_abort('calcHeight_gsv_nl_vcode500x (czp): dataKind=4: Z{T,M}out_r4_opt expected')
      end if
      height_M_ptr_r4 => ZMout_r4_opt 
      height_T_ptr_r4 => ZTout_r4_opt

      ! initialize the height pointer to zero
      height_M_ptr_r4(:,:,:,:) = 0.0
      height_T_ptr_r4(:,:,:,:) = 0.0

      call gsv_getField(statevector,hu_ptr_r4,'HU')
      call gsv_getField(statevector,tt_ptr_r4,'TT')
      call gsv_getField(statevector,P0_ptr_r4,'P0')

    else ! datakind = 8
      if ( .not. (present(PTin_r8_opt) .and. present(PMin_r8_opt))) then
        call utl_abort('calcHeight_gsv_nl_vcode500x (czp): dataKind=8: P{T,M}in_r8_opt expected')
      end if
      P_T_ptr_r8 => PTin_r8_opt
      P_M_ptr_r8 => PMin_r8_opt

      if ( .not. (present(ZTout_r8_opt) .and. present(ZMout_r8_opt))) then
        call utl_abort('calcHeight_gsv_nl_vcode500x (czp): dataKind=8: Z{T,M}out_r8_opt expected')
      end if
      height_M_ptr_r8 => ZMout_r8_opt 
      height_T_ptr_r8 => ZTout_r8_opt

      ! initialize the height pointer to zero
      height_M_ptr_r8(:,:,:,:) = 0.0d0
      height_T_ptr_r8(:,:,:,:) = 0.0d0

      call gsv_getField(statevector,hu_ptr_r8,'HU')
      call gsv_getField(statevector,tt_ptr_r8,'TT')
      call gsv_getField(statevector,P0_ptr_r8,'P0')
    end if

    HeightSfc_ptr_r8 => gsv_getHeightSfc(statevector)

    ! compute virtual temperature on thermo levels (corrected of compressibility)
    do_computeHeight_gsv_nl : do stepIndex = 1, numStep
      do latIndex = statevector%myLatBeg, statevector%myLatEnd
        do lonIndex = statevector%myLonBeg, statevector%myLonEnd

          height_T(:) = 0.0D0
          height_M(:) = 0.0D0

          ! latitude
          lat_4 = statevector%hco%lat2d_4(lonIndex,latIndex)
          lat_8 = real(lat_4,8)
          sLat = sin(lat_8)
          cLat = cos(lat_8)

          do lev_T = 1, nlev_T
            if ( statevector%dataKind == 4 ) then
              hu = real(hu_ptr_r4(lonIndex,latIndex,lev_T,stepIndex),8)
              tt = real(tt_ptr_r4(lonIndex,latIndex,lev_T,stepIndex),8)
              Pr = real(P_T_ptr_r4(lonIndex,latIndex,lev_T,stepIndex),8)
            else
              hu = hu_ptr_r8(lonIndex,latIndex,lev_T,stepIndex)
              tt = tt_ptr_r8(lonIndex,latIndex,lev_T,stepIndex)
              Pr = P_T_ptr_r8(lonIndex,latIndex,lev_T,stepIndex)
            end if
            cmp = gpscompressibility(Pr,tt,hu)

            tv(lev_T) = phf_fotvt8(tt,hu) * cmp
          end do

          rMT = HeightSfc_ptr_r8(lonIndex,latIndex)

          ! compute altitude on bottom momentum level
          if (Vcode == 5002) then
            height_M(nlev_M) = rMT
          else if (Vcode == 5005) then
            height_M(nlev_M) = rMT + heightSfcOffset_M_r4
          end if

          ! compute altitude on 2nd momentum level
          if (nlev_M > 1) then
            if ( statevector%dataKind == 4 ) then
              P_M = real(P_M_ptr_r4(lonIndex,latIndex,nlev_M-1,stepIndex), 8)
              P0  = real(P0_ptr_r4(lonIndex,latIndex,1,stepIndex), 8)
            else
              P_M = P_M_ptr_r8(lonIndex,latIndex,nlev_M-1,stepIndex)
              P0  = P0_ptr_r8(lonIndex,latIndex,1,stepIndex)
            end if

            ratioP  = log( P_M / P0 )

            ! Gravity acceleration
            h0  = rMT
            Rgh = phf_gravityalt(sLat,h0)
            dh  = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(nlev_T-1) * ratioP
            Rgh = phf_gravityalt(sLat, h0+0.5D0*dh)

            delThick = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(nlev_T-1) * ratioP
            height_M(nlev_M-1) = rMT + delThick
          end if

          ! compute altitude on rest of momentum levels
          do lev_M = nlev_M-2, 1, -1
            if ( statevector%dataKind == 4 ) then
              P_M = real(P_M_ptr_r4(lonIndex,latIndex,lev_M,stepIndex),8)
              P_M1 = real(P_M_ptr_r4(lonIndex,latIndex,lev_M+1,stepIndex),8)
            else
              P_M = P_M_ptr_r8(lonIndex,latIndex,lev_M,stepIndex)
              P_M1 = P_M_ptr_r8(lonIndex,latIndex,lev_M+1,stepIndex)
            end if

            ratioP  = log( P_M / P_M1 )

            if (Vcode == 5002) then
              lev_T = lev_M + 1
            else if (Vcode == 5005) then
              lev_T = lev_M
            end if

            ! Gravity acceleration
            h0  = height_M(lev_M+1)
            Rgh = phf_gravityalt(sLat,h0)
            dh  = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(lev_T) * ratioP
            Rgh = phf_gravityalt(sLat, h0+0.5D0*dh)

            delThick   = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(lev_T) * ratioP
            height_M(lev_M) = height_M(lev_M+1) + delThick
          end do

          ! compute Altitude on thermo levels
          if_computeHeight_gsv_nl_vcodes : if (Vcode == 5002) then
            height_T(nlev_T) = height_M(nlev_M)

            do lev_T = 2, nlev_T-1
              lev_M = lev_T ! momentum level just below thermo level being computed

              if ( statevector%dataKind == 4 ) then
                P_T   = real(P_T_ptr_r4(&
                              lonIndex,latIndex,lev_T,stepIndex), 8)
                P_M   = real(P_M_ptr_r4(&
                              lonIndex,latIndex,lev_M,stepIndex), 8)
                P_Mm1 = real(P_M_ptr_r4(&
                              lonIndex,latIndex,lev_M-1,stepIndex), 8)
              else
                P_T   = P_T_ptr_r8(lonIndex,latIndex,lev_T  ,stepIndex)
                P_M   = P_M_ptr_r8(lonIndex,latIndex,lev_M  ,stepIndex)
                P_Mm1 = P_M_ptr_r8(lonIndex,latIndex,lev_M-1,stepIndex)
              end if

              ScaleFactorBottom = log( P_T / P_Mm1 ) / log( P_M / P_Mm1 )
              ScaleFactorTop    = 1 - ScaleFactorBottom
              height_T(lev_T) = ScaleFactorBottom * height_M(lev_M) &
                + ScaleFactorTop * height_M(lev_M-1)
            end do

            ! compute altitude on top thermo level
            if ( statevector%dataKind == 4 ) then
              P_T = real(P_T_ptr_r4(lonIndex,latIndex,1,stepIndex),8)
              P_M = real(P_M_ptr_r4(lonIndex,latIndex,1,stepIndex),8)
            else
              P_T = P_T_ptr_r8(lonIndex,latIndex,1,stepIndex)
              P_M = P_M_ptr_r8(lonIndex,latIndex,1,stepIndex)
            end if

            ratioP = log( P_T / P_M )

            ! Gravity acceleration
            h0  = height_M(1)
            Rgh = phf_gravityalt(sLat, h0)
            dh  = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(1) * ratioP
            Rgh = phf_gravityalt(sLat, h0+0.5D0*dh)

            delThick   = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(1) * ratioP
            height_T(1) = height_M(1) + delThick

          else if (Vcode == 5005) then if_computeHeight_gsv_nl_vcodes
            height_T(nlev_T) = rMT + heightSfcOffset_T_r4

            do lev_T = 1, nlev_T-2
              lev_M = lev_T + 1  ! momentum level just below thermo level being computed
              if ( statevector%dataKind == 4 ) then
                P_T   = real(P_T_ptr_r4(&
                              lonIndex,latIndex,lev_T,stepIndex), 8)
                P_M   = real(P_M_ptr_r4(&
                              lonIndex,latIndex,lev_M,stepIndex), 8)
                P_Mm1 = real(P_M_ptr_r4(&
                              lonIndex,latIndex,lev_M-1,stepIndex), 8)
              else
                P_T   = P_T_ptr_r8(lonIndex,latIndex,lev_T  ,stepIndex)
                P_M   = P_M_ptr_r8(lonIndex,latIndex,lev_M  ,stepIndex)
                P_Mm1 = P_M_ptr_r8(lonIndex,latIndex,lev_M-1,stepIndex)
              end if

              ScaleFactorBottom = log( P_T / P_Mm1 ) / log( P_M / P_Mm1 )
              ScaleFactorTop    = 1 - ScaleFactorBottom
              height_T(lev_T) = ScaleFactorBottom * height_M(lev_M) &
                + ScaleFactorTop * height_M(lev_M-1)
            end do

            ! compute altitude on next to bottom thermo level
            if (nlev_T > 1) then
              if ( statevector%dataKind == 4 ) then
                P_T = real(P_T_ptr_r4(&
                            lonIndex,latIndex,nlev_T-1,stepIndex), 8)
                P0  = real(P0_ptr_r4(&
                            lonIndex,latIndex,1,stepIndex), 8)
              else
                P_T = P_T_ptr_r8(lonIndex,latIndex,nlev_T-1,stepIndex)
                P0  = P0_ptr_r8(lonIndex,latIndex,1,stepIndex)
              end if

              ratioP = log( P_T / P0 )

              h0  = rMT
              Rgh = phf_gravityalt(sLat,h0)
              dh  = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(nlev_T-1) * ratioP
              Rgh = phf_gravityalt(sLat, h0+0.5D0*dh)

              delThick =  (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(nlev_T-1) * &
                          ratioP
              height_T(nlev_T-1) = rMT + delThick
            end if
          end if if_computeHeight_gsv_nl_vcodes

          ! fill the height array
          if ( statevector%dataKind == 4 ) then
            do lev_T = 1, nlev_T
              height_T_ptr_r4(lonIndex,latIndex,lev_T,stepIndex) &
                  = real(height_T(lev_T),4)
            end do
            do lev_M = 1, nlev_M
              height_M_ptr_r4(lonIndex,latIndex,lev_M,stepIndex) &
                  = real(height_M(lev_M),4)
            end do
          else
            height_T_ptr_r8(lonIndex,latIndex,1:nlev_T,stepIndex) &
                = height_T(1:nlev_T)
            height_M_ptr_r8(lonIndex,latIndex,1:nlev_M,stepIndex) &
                = height_M(1:nlev_M)
          end if

          ! remove the height offset for the diagnostic levels for backward compatibility only
          if ( .not. statevector%addHeightSfcOffset ) then
            if ( statevector%dataKind == 4 ) then
              height_T_ptr_r4(lonIndex,latIndex,nlev_T,stepIndex) = &
                  real(rMT,4)
              height_M_ptr_r4(lonIndex,latIndex,nlev_M,stepIndex) = &
                  real(rMT,4)
            else
              height_T_ptr_r8(lonIndex,latIndex,nlev_T,stepIndex) = rMT
              height_M_ptr_r8(lonIndex,latIndex,nlev_M,stepIndex) = rMT
            end if
          end if

        end do
      end do
    end do do_computeHeight_gsv_nl

    deallocate(height_M)
    deallocate(height_T)
    deallocate(tv)

    call msg('calcHeight_gsv_nl_vcode500x (czp)', 'statevector%addHeightSfcOffset='&
         //str(statevector%addHeightSfcOffset), verb_opt=2)
    call msg('calcHeight_gsv_nl_vcode500x (czp)', 'END', verb_opt=4)
  end subroutine calcHeight_gsv_nl_vcode500x

  !---------------------------------------------------------
  ! calcHeight_gsv_tl
  !---------------------------------------------------------
  subroutine calcHeight_gsv_tl(statevector,statevectorRef)
    !
    ! :Purpose: Tangent height computation on statevector.
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout) :: statevector
    type(struct_gsv), intent(in)    :: statevectorRef

    ! Locals
    integer :: Vcode

    call utl_tmg_start(173,'low-level--czp_calcHeight_tl')
    call msg('calcHeight_gsv_tl (czp)', 'START', verb_opt=2)

    Vcode = gsv_getVco(statevectorRef)%vcode
    if (Vcode == 5005 .or. Vcode == 5002) then
      if ( .not. gsv_varExist(statevector,'P_*')  ) then
        call utl_abort('calcHeight_gsv_tl (czp): for vcode 500x, variables P_T and P_M must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'Z_*')  ) then
        call utl_abort('calcHeight_gsv_tl (czp): for vcode 500x, variables Z_T and Z_M must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'TT')  ) then
        call utl_abort('calcHeight_gsv_tl (czp): for vcode 500x, variable TT must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'HU')  ) then
        call utl_abort('calcHeight_gsv_tl (czp): for vcode 500x, variable HU must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'P0')  ) then
        call utl_abort('calcHeight_gsv_tl (czp): for vcode 500x, variable P0 must be allocated in gridstatevector')
      end if
      call calcHeight_gsv_tl_vcode500x
    else if (Vcode == 21001) then
      ! Development notes (@mad001)
      !   probably some some gsv_varExist(statevector,.) needed for GEM-H
      call calcHeight_gsv_tl_vcode2100x
    end if

    call msg('calcHeight_gsv_tl (czp)', 'END', verb_opt=2)
    call utl_tmg_stop(173)

    contains
      !---------------------------------------------------------
      ! calcHeight_gsv_tl_vcode2100x
      !---------------------------------------------------------
      subroutine calcHeight_gsv_tl_vcode2100x
        implicit none

        call utl_abort('calcHeight_gsv_tl (czp): vcode 21001 not implemented yet')

      end subroutine calcHeight_gsv_tl_vcode2100x

      !---------------------------------------------------------
      ! calcHeight_gsv_tl_vcode500x
      !---------------------------------------------------------
      subroutine calcHeight_gsv_tl_vcode500x
        implicit none

        ! Locals
        integer ::  lev_M,lev_T,nlev_M,nlev_T,Vcode_anl
        integer ::  numStep,stepIndex, latIndex,lonIndex
        real(8) :: ScaleFactorBottom, ScaleFactorTop
        real(8), allocatable :: delThick(:,:,:,:)
        real(8), pointer     :: height_T_ptr(:,:,:,:),height_M_ptr(:,:,:,:)
        real(8), pointer     :: P_T(:,:,:,:), P_M(:,:,:,:)
        real(pre_incrReal), pointer ::  delHeight_M_ptr_r48(:,:,:,:)
        real(pre_incrReal), pointer ::  delHeight_T_ptr_r48(:,:,:,:)
        real(pre_incrReal), pointer ::  delTT_r48(:,:,:,:), delHU_r48(:,:,:,:)
        real(pre_incrReal), pointer ::  delP0_r48(:,:,:,:)
        real(pre_incrReal), pointer ::  delP_T_r48(:,:,:,:), delP_M_r48(:,:,:,:)

        call msg('calcHeight_gsv_tl_vcode500x (czp)', 'START', verb_opt=4)
        Vcode_anl = gsv_getVco(statevectorRef)%vcode

        nlev_T = gsv_getNumLev(statevectorRef,'TH')
        nlev_M = gsv_getNumLev(statevectorRef,'MM')
        numStep = statevectorRef%numstep

        allocate(delThick(statevectorRef%myLonBeg:statevectorRef%myLonEnd, &
                          statevectorRef%myLatBeg:statevectorRef%myLatEnd, &
                          nlev_T,numStep))

        ! generate the height coefficients on the grid
        call calcHeightCoeff_gsv(statevectorRef)

        ! loop over all lat/lon/step

        call gsv_getField(statevectorRef,height_M_ptr,'Z_M')
        call gsv_getField(statevectorRef,height_T_ptr,'Z_T')
        call gsv_getField(statevectorRef,P_T,'P_T')
        call gsv_getField(statevectorRef,P_M,'P_M')
        call gsv_getField(statevector,delHeight_M_ptr_r48,'Z_M')
        call gsv_getField(statevector,delHeight_T_ptr_r48,'Z_T')
        call gsv_getField(statevector,delTT_r48,'TT')
        call gsv_getField(statevector,delHU_r48,'HU')
        call gsv_getField(statevector,delP0_r48,'P0')
        call gsv_getField(statevector,delP_T_r48,'P_T')
        call gsv_getField(statevector,delP_M_r48,'P_M')
        ! ensure increment at sfc is zero (fixed height level)
        delHeight_M_ptr_r48(:,:,nlev_M,:) = 0.0d0
        delHeight_T_ptr_r48(:,:,nlev_T,:) = 0.0d0

        if_computeHeight_gsv_tl_vcodes : if(Vcode_anl == 5002) then

          ! compute increment to thickness for each layer between the two momentum levels
          do stepIndex = 1, numStep
            do lev_T = 2, (nlev_T-1)
              do latIndex = statevectorRef%myLatBeg, statevectorRef%myLatEnd
                do lonIndex = statevectorRef%myLonBeg, statevectorRef%myLonEnd
                  delThick(lonIndex,latIndex,lev_T,stepIndex) = &
                      coeff_M_TT_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                      delTT_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                      coeff_M_HU_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                      delHU_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                      coeff_M_P0_delPM_gsv(lonIndex,latIndex,lev_T,stepIndex) *&
                      ( delP_M_r48(lonIndex,latIndex,lev_T  ,stepIndex) / &
                        P_M(lonIndex,latIndex,lev_T  ,stepIndex) - &
                        delP_M_r48(lonIndex,latIndex,lev_T-1,stepIndex) / &
                        P_M(lonIndex,latIndex,lev_T-1,stepIndex) ) + &
                      coeff_M_P0_dP_delPT_gsv(&
                                lonIndex,latIndex,lev_T,stepIndex) * &
                      delP_T_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                      coeff_M_P0_dP_delP0_gsv(&
                                lonIndex,latIndex,lev_T,stepIndex) * &
                      delP0_r48(lonIndex,latIndex,1,stepIndex)
                end do
              end do
            end do
          end do

          ! compute height increment on momentum levels above the surface
          do stepIndex = 1, numStep
            do lev_M = (nlev_M-1), 1, -1
              do latIndex = statevectorRef%myLatBeg, statevectorRef%myLatEnd
                do lonIndex = statevectorRef%myLonBeg, statevectorRef%myLonEnd
                  lev_T = lev_M + 1 ! thermo level just below momentum level being computed
                  delHeight_M_ptr_r48(lonIndex,latIndex,lev_M,stepIndex) =  &
                       delHeight_M_ptr_r48(&
                                  lonIndex,latIndex,lev_M+1,stepIndex) + &
                       delThick(lonIndex,latIndex,lev_T,stepIndex)
                end do
              end do
            end do
          end do

          ! compute height increment on thermo levels using weighted average of height increment of momentum levels
          do stepIndex = 1, numStep
            do lev_T = 1, (nlev_T-1)
              do latIndex = statevectorRef%myLatBeg, statevectorRef%myLatEnd
                do lonIndex = statevectorRef%myLonBeg, statevectorRef%myLonEnd
                  if ( lev_T == 1) then
                    ! compute height increment for top thermo level (from top momentum level)
                    delHeight_T_ptr_r48(lonIndex,latIndex,1,stepIndex) = &
                    delHeight_M_ptr_r48(lonIndex,latIndex,1,stepIndex) + &
                    coeff_T_TT_gsv(lonIndex,latIndex,stepIndex) * &
                    delTT_r48(lonIndex,latIndex,1,stepIndex) + &
                    coeff_T_HU_gsv(lonIndex,latIndex,stepIndex) * &
                    delHU_r48(lonIndex,latIndex,1,stepIndex) + &
                    coeff_T_P0_delP1_gsv(lonIndex,latIndex,stepIndex) * &
                    ( delP_M_r48(lonIndex,latIndex,1,stepIndex) / &
                      P_M(lonIndex,latIndex,1,stepIndex) - &
                      delP_T_r48(lonIndex,latIndex,1,stepIndex) / &
                      P_T(lonIndex,latIndex,1,stepIndex) ) + &
                    coeff_T_P0_dP_delPT_gsv(lonIndex,latIndex,stepIndex) * &
                    delP_T_r48(lonIndex,latIndex,1,stepIndex) + &
                    coeff_T_P0_dP_delP0_gsv(lonIndex,latIndex,stepIndex) * &
                    delP0_r48(lonIndex,latIndex,1,stepIndex)
                  else
                    lev_M = lev_T ! momentum level just below thermo level being computed
                    ScaleFactorBottom = &
                        (height_T_ptr(lonIndex,latIndex,lev_T,stepIndex) - &
                          height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex)) / &
                        (height_M_ptr(lonIndex,latIndex,lev_M,stepIndex) - &
                          height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex))
                    ScaleFactorTop    = 1 - ScaleFactorBottom
                    delHeight_T_ptr_r48(lonIndex,latIndex,lev_T,stepIndex) = &
                        ScaleFactorBottom * &
                        delHeight_M_ptr_r48(&
                                  lonIndex,latIndex,lev_M  ,stepIndex) + &
                        ScaleFactorTop * &
                        delHeight_M_ptr_r48(&
                                  lonIndex,latIndex,lev_M-1,stepIndex)
                  end if
                end do
              end do
            end do
          end do

        else if(Vcode_anl == 5005) then if_computeHeight_gsv_tl_vcodes

          ! compute increment to thickness for each layer between the two momentum levels
          do stepIndex = 1, numStep
            do lev_T = 1, (nlev_T-1)
              do latIndex = statevectorRef%myLatBeg, statevectorRef%myLatEnd
                do lonIndex = statevectorRef%myLonBeg, statevectorRef%myLonEnd
                  delThick(lonIndex,latIndex,lev_T,stepIndex) = &
                      coeff_M_TT_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                      delTT_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                      coeff_M_HU_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                      delHU_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                      coeff_M_P0_delPM_gsv(lonIndex,latIndex,lev_T,stepIndex) *&
                      ( delP_M_r48(lonIndex,latIndex,lev_T+1,stepIndex) / &
                        P_M(lonIndex,latIndex,lev_T+1,stepIndex) - &
                        delP_M_r48(lonIndex,latIndex,lev_T  ,stepIndex) / &
                        P_M(lonIndex,latIndex,lev_T  ,stepIndex) ) + &
                      coeff_M_P0_dP_delPT_gsv(&
                                      lonIndex,latIndex,lev_T,stepIndex) * &
                      delP_T_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                      coeff_M_P0_dP_delP0_gsv(&
                                      lonIndex,latIndex,lev_T,stepIndex) * &
                      delP0_r48(lonIndex,latIndex,1,stepIndex)
                end do
              end do
            end do
          end do

          ! compute height increment on momentum levels above the surface
          do stepIndex = 1, numStep
            do lev_M = (nlev_M-1), 1, -1
              do latIndex = statevectorRef%myLatBeg, statevectorRef%myLatEnd
                do lonIndex = statevectorRef%myLonBeg, statevectorRef%myLonEnd
                  lev_T = lev_M ! thermo level just below momentum level being computed
                  delHeight_M_ptr_r48(lonIndex,latIndex,lev_M,stepIndex) = &
                  delHeight_M_ptr_r48(lonIndex,latIndex,lev_M+1,stepIndex) + &
                  delThick(lonIndex,latIndex,lev_T,stepIndex)
                end do
              end do
            end do
          end do

          ! compute height increment on thermo levels using weighted average of height increment of momentum levels
          do stepIndex = 1, numStep
            do lev_T = 1, (nlev_T-1)
              do latIndex = statevectorRef%myLatBeg, statevectorRef%myLatEnd
                do lonIndex = statevectorRef%myLonBeg, statevectorRef%myLonEnd
                  lev_M = lev_T + 1 ! momentum level just below thermo level being computed
                  ScaleFactorBottom = &
                      (height_T_ptr(lonIndex,latIndex,lev_T,stepIndex) - &
                        height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex)) / &
                      (height_M_ptr(lonIndex,latIndex,lev_M,stepIndex) - &
                        height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex))
                  ScaleFactorTop    = 1 - ScaleFactorBottom
                  delHeight_T_ptr_r48(lonIndex,latIndex,lev_T,stepIndex) = &
                      ScaleFactorBottom * &
                      delHeight_M_ptr_r48(lonIndex,latIndex,lev_M  ,stepIndex) + &
                      ScaleFactorTop * &
                      delHeight_M_ptr_r48(lonIndex,latIndex,lev_M-1,stepIndex)
                end do
              end do
            end do
          end do

        end if if_computeHeight_gsv_tl_vcodes

        deallocate(delThick)

        call msg('calcHeight_gsv_tl_vcode500x (czp)', 'END', verb_opt=4)
      end subroutine calcHeight_gsv_tl_vcode500x

  end subroutine calcHeight_gsv_tl

  !---------------------------------------------------------
  ! calcHeight_gsv_ad
  !---------------------------------------------------------
  subroutine calcHeight_gsv_ad(statevector,statevectorRef)
    !
    ! :Purpose: Adjoint of height computation on statevector.
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout) :: statevector
    type(struct_gsv), intent(in)    :: statevectorRef

    ! Locals
    integer :: Vcode

    call utl_tmg_start(174,'low-level--czp_calcHeight_ad')
    call msg('calcHeight_gsv_ad', 'START', verb_opt=2)

    Vcode = gsv_getVco(statevectorRef)%vcode
    if (Vcode == 5005 .or. Vcode == 5002) then
      if ( .not. gsv_varExist(statevector,'P_*')  ) then
        call utl_abort('calcHeight_gsv_ad (czp): for vcode 500x, variables P_M and P_T must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'Z_*')  ) then
        call utl_abort('calcHeight_gsv_ad (czp): for vcode 500x, variables Z_M and Z_T must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'TT')  ) then
        call utl_abort('calcHeight_gsv_ad (czp): for vcode 500x, variable TT must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'HU')  ) then
        call utl_abort('calcHeight_gsv_ad (czp): for vcode 500x, variable HU must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'P0')  ) then
        call utl_abort('calcHeight_gsv_ad (czp): for vcode 500x, variable P0 must be allocated in gridstatevector')
      end if
      call calcHeight_gsv_ad_vcode500x
    else if (Vcode == 21001) then
      ! Development notes (@mad001)
      !   probably some some gsv_varExist(statevector,.) needed for GEM-H
      call calcHeight_gsv_ad_vcode2100x
    end if

    call msg('calcHeight_gsv_ad', 'END', verb_opt=2)
    call utl_tmg_stop(174)

    contains
      !---------------------------------------------------------
      ! calcHeight_gsv_ad_vcode2100x
      !---------------------------------------------------------
      subroutine calcHeight_gsv_ad_vcode2100x
        implicit none

        call utl_abort('calcHeight_gsv_ad (czp): vcode 21001 not implemented yet')

      end subroutine calcHeight_gsv_ad_vcode2100x

      !---------------------------------------------------------
      ! calcHeight_gsv_ad
      !---------------------------------------------------------
      subroutine calcHeight_gsv_ad_vcode500x
        implicit none

        ! Locals
        integer ::  lev_M,lev_T,nlev_M,nlev_T
        integer ::  numStep,stepIndex,latIndex,lonIndex
        real(8) ::  ScaleFactorBottom, ScaleFactorTop
        real(8), allocatable :: delThick(:,:,:,:)
        real(8), pointer     :: height_M_ptr(:,:,:,:),height_T_ptr(:,:,:,:)
        real(8), allocatable :: delHeight_M(:,:,:,:),delHeight_T(:,:,:,:)
        real(8), pointer     :: P_M(:,:,:,:),P_T(:,:,:,:)
        real(pre_incrReal), pointer :: delHeight_M_ptr_r48(:,:,:,:)
        real(pre_incrReal), pointer :: delHeight_T_ptr_r48(:,:,:,:)
        real(pre_incrReal), pointer :: delTT_r48(:,:,:,:),delHU_r48(:,:,:,:)
        real(pre_incrReal), pointer :: delP0_r48(:,:,:,:)
        real(pre_incrReal), pointer :: delP_M_r48(:,:,:,:),delP_T_r48(:,:,:,:)

        call msg('calcHeight_gsv_ad_vcode500x (czp)', 'START', verb_opt=4)

        nlev_T = gsv_getNumLev(statevectorRef,'TH')
        nlev_M = gsv_getNumLev(statevectorRef,'MM')
        numStep = statevectorRef%numstep

        allocate(delHeight_M(statevectorRef%myLonBeg:statevectorRef%myLonEnd, &
                             statevectorRef%myLatBeg:statevectorRef%myLatEnd, &
                             nlev_M,numStep))
        allocate(delHeight_T(statevectorRef%myLonBeg:statevectorRef%myLonEnd, &
                             statevectorRef%myLatBeg:statevectorRef%myLatEnd, &
                             nlev_T,numStep))
        allocate(delThick(statevectorRef%myLonBeg:statevectorRef%myLonEnd, &
                          statevectorRef%myLatBeg:statevectorRef%myLatEnd, &
                          0:nlev_T,numStep))

        ! generate the height coefficients on the grid
        call calcHeightCoeff_gsv(statevectorRef)

        ! loop over all lat/lon/step
        call gsv_getField(statevectorRef,height_M_ptr,'Z_M')
        call gsv_getField(statevectorRef,height_T_ptr,'Z_T')
        call gsv_getField(statevectorRef,P_T,'P_T')
        call gsv_getField(statevectorRef,P_M,'P_M')

        call gsv_getField(statevector,delHeight_M_ptr_r48,'Z_M')
        call gsv_getField(statevector,delHeight_T_ptr_r48,'Z_T')
        call gsv_getField(statevector,delTT_r48,'TT')
        call gsv_getField(statevector,delHU_r48,'HU')
        call gsv_getField(statevector,delP0_r48,'P0')
        call gsv_getField(statevector,delP_T_r48,'P_T')
        call gsv_getField(statevector,delP_M_r48,'P_M')
        delHeight_M(:,:,:,:) = delHeight_M_ptr_r48(:,:,:,:)
        delHeight_T(:,:,:,:) = delHeight_T_ptr_r48(:,:,:,:)

        if_computeHeight_gsv_ad_vcodes : if(Vcode == 5002) then

          ! adjoint of compute height increment on thermo levels by simple averaging
          do stepIndex = 1, numStep
            do lev_T = 1, (nlev_T-1)
              do latIndex = statevectorRef%myLatBeg, statevectorRef%myLatEnd
                do lonIndex = statevectorRef%myLonBeg, statevectorRef%myLonEnd
                  lev_M = lev_T ! momentum level just below thermo level being computed

                  ! adjoint of compute height increment on top thermo level
                  ! (from top momentum level)
                  if (lev_T == 1) then
                    delHeight_M(lonIndex,latIndex,1,stepIndex)  =  &
                        delHeight_M(lonIndex,latIndex,1,stepIndex) + &
                        delHeight_T(lonIndex,latIndex,1,stepIndex)

                    delTT_r48(lonIndex,latIndex,1,stepIndex) =  &
                        delTT_r48(lonIndex,latIndex,1,stepIndex) + &
                        coeff_T_TT_gsv(lonIndex,latIndex,stepIndex) * &
                        delHeight_T(lonIndex,latIndex,1,stepIndex)

                    delHU_r48(lonIndex,latIndex,1,stepIndex) =  &
                        delHU_r48(lonIndex,latIndex,1,stepIndex) + &
                        coeff_T_HU_gsv   (lonIndex,latIndex,stepIndex) * &
                        delHeight_T(lonIndex,latIndex,1,stepIndex)

                    delP_M_r48(lonIndex,latIndex,1,stepIndex) =  &
                        delP_M_r48(lonIndex,latIndex,1,stepIndex) + &
                        coeff_T_P0_delP1_gsv(lonIndex,latIndex,stepIndex) / &
                        P_M(lonIndex,latIndex,1,stepIndex) * &
                        delHeight_T(lonIndex,latIndex,1,stepIndex)

                    delP_T_r48(lonIndex,latIndex,1,stepIndex) =  &
                        delP_T_r48(lonIndex,latIndex,1,stepIndex) - &
                        coeff_T_P0_delP1_gsv(lonIndex,latIndex,stepIndex) / &
                        P_T(lonIndex,latIndex,1,stepIndex) * &
                        delHeight_T(lonIndex,latIndex,1,stepIndex)

                    delP_T_r48(lonIndex,latIndex,1,stepIndex) =  &
                        delP_T_r48(lonIndex,latIndex,1,stepIndex) + &
                        coeff_T_P0_dP_delPT_gsv(lonIndex,latIndex,stepIndex) * &
                        delHeight_T(lonIndex,latIndex,1,stepIndex)

                    delP0_r48(lonIndex,latIndex,1,stepIndex) =  &
                        delP0_r48(lonIndex,latIndex,1,stepIndex) + &
                        coeff_T_P0_dp_delP0_gsv(lonIndex,latIndex,stepIndex) * &
                        delHeight_T(lonIndex,latIndex,1,stepIndex)
                  else
                    ScaleFactorBottom =  &
                        (height_T_ptr(lonIndex,latIndex,lev_T,stepIndex) - &
                          height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex)) / &
                        (height_M_ptr(lonIndex,latIndex,lev_M,stepIndex) - &
                          height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex))
                    ScaleFactorTop    = 1 - ScaleFactorBottom

                    delHeight_M(lonIndex,latIndex,lev_M-1,stepIndex) =  &
                        delHeight_M(lonIndex,latIndex,lev_M-1,stepIndex) + &
                        ScaleFactorTop * &
                        delHeight_T(lonIndex,latIndex,lev_T,stepIndex)

                    delHeight_M(lonIndex,latIndex,lev_M,stepIndex) =  &
                        delHeight_M(lonIndex,latIndex,lev_M  ,stepIndex) + &
                        ScaleFactorBottom * &
                        delHeight_T(lonIndex,latIndex,lev_T,stepIndex)
                  end if
                end do
              end do
            end do
          end do

          ! adjoint of compute height increment on momentum levels above the surface
          delThick(:,:,0:1,:) = 0.0d0
          do stepIndex = 1, numStep
            do lev_M = 1, (nlev_M-1)
              do latIndex = statevectorRef%myLatBeg, statevectorRef%myLatEnd
                do lonIndex = statevectorRef%myLonBeg, statevectorRef%myLonEnd
                  lev_T = lev_M + 1 ! thermo level just below momentum level being computed
                  delThick(lonIndex,latIndex,lev_T,stepIndex) =  &
                       delThick(lonIndex,latIndex,lev_T-1,stepIndex) + &
                       delHeight_M(lonIndex,latIndex,lev_M  ,stepIndex)
                end do
              end do
            end do
          end do

          ! adjoint of compute increment to thickness for each layer between the two momentum levels
          do stepIndex = 1, numStep
            do lev_T = 2, nlev_T-1
              do latIndex = statevectorRef%myLatBeg, statevectorRef%myLatEnd
                do lonIndex = statevectorRef%myLonBeg, statevectorRef%myLonEnd
                  delTT_r48(lonIndex,latIndex,lev_T,stepIndex) =  &
                      delTT_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                      coeff_M_TT_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                      delThick(lonIndex,latIndex,lev_T,stepIndex)

                  delHU_r48(lonIndex,latIndex,lev_T,stepIndex) =  &
                      delHU_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                      coeff_M_HU_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                      delThick(lonIndex,latIndex,lev_T,stepIndex)

                  delP_M_r48(lonIndex,latIndex,lev_T,stepIndex)=  &
                      delP_M_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                      coeff_M_P0_delPM_gsv(lonIndex,latIndex,lev_T,stepIndex) /&
                      P_M(lonIndex,latIndex,lev_T,stepIndex) * &
                      delThick(lonIndex,latIndex,lev_T,stepIndex)

                  delP_M_r48(lonIndex,latIndex,lev_T-1,stepIndex) =  &
                      delP_M_r48(lonIndex,latIndex,lev_T-1,stepIndex) - &
                      coeff_M_P0_delPM_gsv(lonIndex,latIndex,lev_T,stepIndex) /&
                      P_M(lonIndex,latIndex,lev_T-1,stepIndex) * &
                      delThick(lonIndex,latIndex,lev_T,stepIndex)

                  delP_T_r48(lonIndex,latIndex,lev_T,stepIndex) =  &
                      delP_T_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                      coeff_M_P0_dP_delPT_gsv(&
                                    lonIndex,latIndex,lev_T,stepIndex) * &
                      delThick(lonIndex,latIndex,lev_T,stepIndex)

                  delP0_r48(lonIndex,latIndex,1,stepIndex) =  &
                      delP0_r48(lonIndex,latIndex,1,stepIndex) + &
                      coeff_M_P0_dP_delP0_gsv(&
                                    lonIndex,latIndex,lev_T,stepIndex) * &
                      delThick(lonIndex,latIndex,lev_T,stepIndex)
                end do
              end do
            end do
          end do

        else if(Vcode == 5005) then if_computeHeight_gsv_ad_vcodes

          ! adjoint of compute height increment on thermo levels by simple averaging
          do stepIndex = 1, numStep
            do lev_T = 1, (nlev_T-1)
              do latIndex = statevectorRef%myLatBeg, statevectorRef%myLatEnd
                do lonIndex = statevectorRef%myLonBeg, statevectorRef%myLonEnd
                  lev_M = lev_T+1 ! momentum level just below thermo level being computed
                  ScaleFactorBottom = &
                      (height_T_ptr(lonIndex,latIndex,lev_T,stepIndex) - &
                        height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex)) / &
                      (height_M_ptr(lonIndex,latIndex,lev_M,stepIndex) -&
                        height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex))
                  ScaleFactorTop = 1 - ScaleFactorBottom
                  delHeight_M(lonIndex,latIndex,lev_M-1,stepIndex) = &
                      delHeight_M(lonIndex,latIndex,lev_M-1,stepIndex) + &
                      ScaleFactorTop * &
                      delHeight_T(lonIndex,latIndex,lev_T,stepIndex)
                  delHeight_M(lonIndex,latIndex,lev_M,stepIndex) = &
                      delHeight_M(lonIndex,latIndex,lev_M  ,stepIndex) + &
                      ScaleFactorBottom * &
                      delHeight_T(lonIndex,latIndex,lev_T,stepIndex)
                end do
              end do
            end do
          end do

          ! adjoint of compute height increment on momentum levels
          delThick(:,:,0,:) = 0.0d0
          do stepIndex = 1, numStep
            do lev_M = 1, (nlev_M-1)
              do latIndex = statevectorRef%myLatBeg, statevectorRef%myLatEnd
                do lonIndex = statevectorRef%myLonBeg, statevectorRef%myLonEnd
                  lev_T = lev_M ! thermo level just below momentum level being computed
                  delThick(lonIndex,latIndex,lev_T,stepIndex) = &
                      delThick(lonIndex,latIndex,lev_T-1,stepIndex) + &
                      delHeight_M (lonIndex,latIndex,lev_M,stepIndex)
                end do
              end do
            end do
          end do

          do stepIndex = 1, numStep
            do lev_T = 1, nlev_T-1
              do latIndex = statevectorRef%myLatBeg, statevectorRef%myLatEnd
                do lonIndex = statevectorRef%myLonBeg, statevectorRef%myLonEnd
                  delTT_r48(lonIndex,latIndex,lev_T,stepIndex) =  &
                      delTT_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                      coeff_M_TT_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                      delThick(lonIndex,latIndex,lev_T,stepIndex)

                  delHU_r48(lonIndex,latIndex,lev_T,stepIndex) =  &
                      delHU_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                      coeff_M_HU_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                      delThick(lonIndex,latIndex,lev_T,stepIndex)

                  delP_M_r48(lonIndex,latIndex,lev_T+1,stepIndex) =  &
                      delP_M_r48(lonIndex,latIndex,lev_T+1,stepIndex) + &
                      coeff_M_P0_delPM_gsv(lonIndex,latIndex,lev_T,stepIndex) /&
                      P_M(lonIndex,latIndex,lev_T+1,stepIndex) * &
                      delThick(lonIndex,latIndex,lev_T,stepIndex)

                  delP_M_r48(lonIndex,latIndex,lev_T  ,stepIndex) =  &
                      delP_M_r48(lonIndex,latIndex,lev_T  ,stepIndex) - &
                      coeff_M_P0_delPM_gsv(lonIndex,latIndex,lev_T,stepIndex) /&
                      P_M(lonIndex,latIndex,lev_T  ,stepIndex) * &
                      delThick(lonIndex,latIndex,lev_T,stepIndex)

                  delP_T_r48(lonIndex,latIndex,lev_T  ,stepIndex) =  &
                      delP_T_r48(lonIndex,latIndex,lev_T  ,stepIndex) + &
                      coeff_M_P0_dP_delPT_gsv( lonIndex,latIndex,lev_T,stepIndex) * &
                      delThick(lonIndex,latIndex,lev_T,stepIndex)

                  delP0_r48(lonIndex,latIndex,1,stepIndex)     =  &
                      delP0_r48(lonIndex,latIndex,1,stepIndex) + &
                      coeff_M_P0_dP_delP0_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                      delThick(lonIndex,latIndex,lev_T,stepIndex)
                end do
              end do
            end do
          end do

        end if if_computeHeight_gsv_ad_vcodes

        deallocate(delThick)
        deallocate(delHeight_M)
        deallocate(delHeight_T)

        call msg('calcHeight_gsv_ad_vcode500x (czp)', 'END', verb_opt=4)
      end subroutine calcHeight_gsv_ad_vcode500x

  end subroutine calcHeight_gsv_ad

  !---------------------------------------------------------
  ! calcPressure_gsv_nl
  !---------------------------------------------------------
  subroutine calcPressure_gsv_nl(statevector, Ps_in_hPa_opt)
    !
    ! :Purpose: Pressure computation, values stored in statevector.
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout) :: statevector
    logical, optional, intent(in)   :: Ps_in_hPa_opt            ! If true, conversion from hPa to mbar will be done for surface pressure

    ! Locals
    integer :: Vcode

    real(4), pointer :: ptr_ZT_r4(:,:,:,:), ptr_ZM_r4(:,:,:,:)
    real(8), pointer :: ptr_ZT_r8(:,:,:,:), ptr_ZM_r8(:,:,:,:)
    real(4), pointer :: ptr_PT_r4(:,:,:,:), ptr_PM_r4(:,:,:,:)
    real(8), pointer :: ptr_PT_r8(:,:,:,:), ptr_PM_r8(:,:,:,:)

    call utl_tmg_start(177,'low-level--czp_calcPressure_nl')
    call msg('calcPressure_gsv_nl (czp)', 'START', verb_opt=2)

    Vcode = gsv_getVco(statevector)%vcode
    if (Vcode == 5005 .or. Vcode == 5002) then
      if ( gsv_getDataKind(statevector) == 4 ) then
        call gsv_getField(statevector, ptr_PT_r4, 'P_T')
        call gsv_getField(statevector, ptr_PM_r4, 'P_M')
        call calcPressure_gsv_nl_vcode500x_r4(statevector, &
                                              ptr_PT_r4, ptr_PM_r4, &
                                              Ps_in_hPa_opt=Ps_in_hPa_opt)
      else
        call gsv_getField(statevector, ptr_PT_r8, 'P_T')
        call gsv_getField(statevector, ptr_PM_r8, 'P_M')
        call calcPressure_gsv_nl_vcode500x_r8(statevector, &
                                              ptr_PT_r8, ptr_PM_r8, &
                                              Ps_in_hPa_opt=Ps_in_hPa_opt)
      end if
    else if (Vcode == 21001) then
      ! Development notes (@mad001)
      !   probably some some gsv_varExist(statevector,.) needed for GEM-H
      if ( gsv_getDataKind(statevector) == 4 ) then
        call gsv_getField(statevector, ptr_ZT_r4, 'Z_T')
        call gsv_getField(statevector, ptr_ZM_r4, 'Z_M')
        call gsv_getField(statevector, ptr_PT_r4, 'P_T')
        call gsv_getField(statevector, ptr_PM_r4, 'P_M')
        call calcPressure_gsv_nl_vcode2100x(statevector, &
                                            ZTin_r4_opt=ptr_ZT_r4, &
                                            ZMin_r4_opt=ptr_ZM_r4, &
                                            PTout_r4_opt=ptr_PT_r4, &
                                            PMout_r4_opt=ptr_PM_r4)
      else
        call gsv_getField(statevector, ptr_ZT_r8, 'Z_T')
        call gsv_getField(statevector, ptr_ZM_r8, 'Z_M')
        call gsv_getField(statevector, ptr_PT_r8, 'P_T')
        call gsv_getField(statevector, ptr_PM_r8, 'P_M')
        call calcPressure_gsv_nl_vcode2100x(statevector, &
                                            ZTin_r8_opt=ptr_ZT_r8, &
                                            ZMin_r8_opt=ptr_ZM_r8, &
                                            PTout_r8_opt=ptr_PT_r8, &
                                            PMout_r8_opt=ptr_PM_r8)
      end if
    end if

    if ( gsv_getDataKind(statevector) == 4 ) then
      call msg('calcPressure_gsv_nl (czp)', &
             new_line('')//'P_M = '&
           //str(ptr_PM_r4( statevector%myLonBeg,statevector%myLatBeg,:,1), vertical_opt=.false.) &
           //new_line('')//'P_T = '&
           //str(ptr_PT_r4( statevector%myLonBeg,statevector%myLatBeg,:,1), vertical_opt=.false.), &
           verb_opt=2)
    else
      call msg('calcPressure_gsv_nl (czp)', &
             new_line('')//'P_M = '&
           //str(ptr_PM_r8( statevector%myLonBeg,statevector%myLatBeg,:,1), vertical_opt=.false.) &
           //new_line('')//'P_T = '&
           //str(ptr_PT_r8( statevector%myLonBeg,statevector%myLatBeg,:,1), vertical_opt=.false.), &
           verb_opt=2)
    end if

    call msg('calcPressure_gsv_nl (czp)', 'END', verb_opt=2)
    call utl_tmg_stop(177)
  end subroutine calcPressure_gsv_nl

  !---------------------------------------------------------
  ! czp_calcReturnPressure_gsv_nl
  !---------------------------------------------------------
  subroutine czp_calcReturnPressure_gsv_nl( statevector, &
                                            ZTin_r4_opt, ZMin_r4_opt, &
                                            ZTin_r8_opt, ZMin_r8_opt, &
                                            PTout_r4_opt, PMout_r4_opt, & 
                                            PTout_r8_opt, PMout_r8_opt, &
                                            Ps_in_hPa_opt)
    !
    ! :Purpose: Compute or retrieve pressures and return values in pointer arguments.
    !           Proceeds to vcode dispatching.
    !
    implicit none

    ! Arguments
    type(struct_gsv),           intent(in)    :: statevector
    real(4), optional, pointer, intent(in)    :: ZTin_r4_opt(:,:,:,:)
    real(4), optional, pointer, intent(in)    :: ZMin_r4_opt(:,:,:,:)
    real(8), optional, pointer, intent(in)    :: ZTin_r8_opt(:,:,:,:)
    real(8), optional, pointer, intent(in)    :: ZMin_r8_opt(:,:,:,:)
    real(4), optional, pointer, intent(inout) :: PTout_r4_opt(:,:,:,:)
    real(4), optional, pointer, intent(inout) :: PMout_r4_opt(:,:,:,:)
    real(8), optional, pointer, intent(inout) :: PTout_r8_opt(:,:,:,:)
    real(8), optional, pointer, intent(inout) :: PMout_r8_opt(:,:,:,:)
    logical, optional,          intent(in)    :: Ps_in_hPa_opt            ! If true, conversion from hPa to mbar will be done for surface pressure

    ! Locals
    integer :: Vcode

    call utl_tmg_start(177,'low-level--czp_calcPressure_nl')
    call msg('czp_calcReturnPressure_gsv_nl', 'START', verb_opt=2)

    Vcode = gsv_getVco(statevector)%vcode
    if (Vcode == 5005 .or. Vcode == 5002) then
      if ( gsv_getDataKind(statevector) == 4 ) then
        if ( .not. (present(PTout_r4_opt) .and. present(PMout_r4_opt))) then
          call utl_abort('czp_calcReturnPressure_gsv_nl: dataKind=4: P{T,M}out_r4_opt expected')
        end if
        call calcPressure_gsv_nl_vcode500x_r4(statevector, PTout_r4_opt, PMout_r4_opt, &
                                              Ps_in_hPa_opt)
      else
        if ( .not. (present(PTout_r8_opt) .and. present(PMout_r8_opt))) then
          call utl_abort('czp_calcReturnPressure_gsv_nl: dataKind=8: P{T,M}out_r8_opt expected')
        end if
        call calcPressure_gsv_nl_vcode500x_r8(statevector, PTout_r8_opt, PMout_r8_opt, &
                                              Ps_in_hPa_opt)
      end if
    else if (Vcode == 21001) then
      if ( gsv_getDataKind(statevector) == 4 ) then
        if ( .not. (present(ZTin_r4_opt) .and. present(ZMin_r4_opt))) then
          call utl_abort('czp_calcReturnPressure_gsv_nl: dataKind=4: Z{T,M}out_r4_opt expected')
        end if
        if ( .not. (present(PTout_r4_opt) .and. present(PMout_r4_opt))) then
          call utl_abort('czp_calcReturnPressure_gsv_nl: dataKind=4: P{T,M}out_r4_opt expected')
        end if
        call calcPressure_gsv_nl_vcode2100x(statevector, &
                                            ZTin_r4_opt=ZTin_r4_opt, &
                                            ZMin_r4_opt=ZMin_r4_opt, &
                                            PTout_r4_opt=PTout_r4_opt, &
                                            PMout_r4_opt=PMout_r4_opt)
      else ! datakind = 8
        if ( .not. (present(ZTin_r8_opt) .and. present(ZMin_r8_opt))) then
          call utl_abort('czp_calcReturnPressure_gsv_nl: dataKind=8: Z{T,M}out_r8_opt expected')
        end if
        if ( .not. (present(PTout_r8_opt) .and. present(PMout_r8_opt))) then
          call utl_abort('czp_calcReturnPressure_gsv_nl: dataKind=8: P{T,M}out_r8_opt expected')
        end if
        call calcPressure_gsv_nl_vcode2100x(statevector, &
                                            ZTin_r8_opt=ZTin_r8_opt, &
                                            ZMin_r8_opt=ZMin_r8_opt, &
                                            PTout_r8_opt=PTout_r8_opt, &
                                            PMout_r8_opt=PMout_r8_opt)
      end if
    end if

    call msg('czp_calcReturnPressure_gsv_nl', 'END', verb_opt=2)
    call utl_tmg_stop(177)
  end subroutine czp_calcReturnPressure_gsv_nl
  
  !---------------------------------------------------------
  ! calcPressure_gsv_nl_vcode2100x
  !---------------------------------------------------------
  subroutine calcPressure_gsv_nl_vcode2100x(statevector, &
                                            ZTin_r4_opt, ZMin_r4_opt, &
                                            ZTin_r8_opt, ZMin_r8_opt, &
                                            PTout_r4_opt, PMout_r4_opt, &
                                            PTout_r8_opt, PMout_r8_opt)
    !
    ! :Purpose: Compute pressure and return values in pointer arguments.
    !           GEM-H statevector input.
    !
    implicit none

    ! Arguments
    type(struct_gsv),           intent(in)    :: statevector
    real(4), pointer, optional, intent(in)    :: ZTin_r4_opt(:,:,:,:)
    real(4), pointer, optional, intent(in)    :: ZMin_r4_opt(:,:,:,:)
    real(8), pointer, optional, intent(in)    :: ZTin_r8_opt(:,:,:,:)
    real(8), pointer, optional, intent(in)    :: ZMin_r8_opt(:,:,:,:)
    real(4), pointer, optional, intent(inout) :: PTout_r4_opt(:,:,:,:)
    real(4), pointer, optional, intent(inout) :: PMout_r4_opt(:,:,:,:)
    real(8), pointer, optional, intent(inout) :: PTout_r8_opt(:,:,:,:)
    real(8), pointer, optional, intent(inout) :: PMout_r8_opt(:,:,:,:)

    ! Locals
    integer ::  stepIndex, latIndex, lonIndex, numStep
    integer ::  lev_M,lev_T,nlev_M,nlev_T,status

    real(4) ::  heightSfcOffset_T_r4, heightSfcOffset_M_r4
    real(4) ::  lat_4
    real(8) ::  hu, tt, cmp, Rgh, P0, dh, tv0, rMt, Z_T, Z_M, Z_M1, logP
    real(8) ::  sLat, cLat, lat_8
    real(8) ::  ScaleFactorBottom

    real(8), allocatable :: tv(:),  pressure_T(:), pressure_M(:)

    real(4), pointer     :: height_T_ptr_r4(:,:,:,:)
    real(4), pointer     :: height_M_ptr_r4(:,:,:,:)
    real(4), pointer     :: hu_ptr_r4(:,:,:,:),tt_ptr_r4(:,:,:,:)
    real(4), pointer     :: P_T_ptr_r4(:,:,:,:),P_M_ptr_r4(:,:,:,:)
    real(4), pointer     :: P0_ptr_r4(:,:,:,:)

    real(8), pointer     :: height_T_ptr_r8(:,:,:,:)
    real(8), pointer     :: height_M_ptr_r8(:,:,:,:)
    real(8), pointer     :: P_T_ptr_r8(:,:,:,:),P_M_ptr_r8(:,:,:,:)
    real(8), pointer     :: P0_ptr_r8(:,:,:,:)
    real(8), pointer     :: hu_ptr_r8(:,:,:,:),tt_ptr_r8(:,:,:,:)
    real(8), pointer     :: HeightSfc_ptr_r8(:,:)

    call msg('calcPressure_gsv_nl_vcode2100x (czp)', 'START', verb_opt=4)
    
    nlev_T = gsv_getNumLev(statevector,'TH')
    nlev_M = gsv_getNumLev(statevector,'MM')
    numStep = statevector%numStep

    allocate(tv(nlev_T))

    if (nlev_T /= nlev_M) then
      call utl_abort('calcPressure_gsv_nl_vcode2100x: nlev_T is not equal to nlev_M!')
    end if

    status = vgd_get( statevector%vco%vgrid, &
                      key='DHM - height of the diagnostic level (m)', &
                      value=heightSfcOffset_M_r4)
    status = vgd_get( statevector%vco%vgrid, &
                      key='DHT - height of the diagnostic level (t)', &
                      value=heightSfcOffset_T_r4)
    call msg('calcPressure_gsv_nl_vcode2100x (czp)', &
         'height offset for near-sfc momentum level is:'//str(heightSfcOffset_M_r4)//' meters'&
         //new_line('')//'height offset for near-sfc thermo level is:'//str(heightSfcOffset_T_r4)//' meters', &
         verb_opt=2, mpiAll_opt=.false.)
    if ( .not.statevector%addHeightSfcOffset ) then
      call msg('calcPressure_gsv_nl_vcode2100x (czp)', new_line('') &
             //'--------------------------------------------------------------------------'//new_line('')&
             //'BUT HEIGHT OFFSET REMOVED FOR DIAGNOSTIC LEVELS FOR BACKWARD COMPATIBILITY'//new_line('')&
             //'--------------------------------------------------------------------------', &
             verb_opt=2, mpiAll_opt=.false.)
    end if

    allocate(pressure_T(nlev_T))
    allocate(pressure_M(nlev_M))

    if ( statevector%dataKind == 4 ) then
      if ( .not. (present(ZTin_r4_opt) .and. present(ZMin_r4_opt))) then
        call utl_abort('calcPressure_gsv_nl_vcode2100x (czp): dataKind=4: Z{T,M}in_r4_opt expected')
      end if
      height_T_ptr_r4 => ZTin_r4_opt
      height_M_ptr_r4 => ZMin_r4_opt

      if ( .not. (present(PTout_r4_opt) .and. present(PMout_r4_opt))) then
        call utl_abort('calcPressure_gsv_nl_vcode2100x (czp): dataKind=4: P{T,M}out_r4_opt expected')
      end if
      P_M_ptr_r4 => PMout_r4_opt
      P_T_ptr_r4 => PTout_r4_opt
      call gsv_getField(statevector,hu_ptr_r4,'HU')
      call gsv_getField(statevector,tt_ptr_r4,'TT')
      call gsv_getField(statevector,P0_ptr_r4,'P0')

      ! initialize the pressure pointer to zero
      P_M_ptr_r4(:,:,:,:) = 0.0
      P_T_ptr_r4(:,:,:,:) = 0.0
    else ! datakind = 8
      if ( .not. (present(ZTin_r8_opt) .and. present(ZMin_r8_opt))) then
        call utl_abort('calcPressure_gsv_nl_vcode2100x (czp): dataKind=4: Z{T,M}in_r8_opt expected')
      end if
      height_T_ptr_r8 => ZTin_r8_opt
      height_M_ptr_r8 => ZMin_r8_opt

      if ( .not. (present(PTout_r8_opt) .and. present(PMout_r8_opt))) then
        call utl_abort('calcPressure_gsv_nl_vcode2100x (czp): dataKind=8: P{T,M}_r8_opt expected')
      end if
      P_M_ptr_r8 => PMout_r8_opt
      P_T_ptr_r8 => PTout_r8_opt
      call gsv_getField(statevector,hu_ptr_r8,'HU')
      call gsv_getField(statevector,tt_ptr_r8,'TT')
      call gsv_getField(statevector,P0_ptr_r8,'P0')

      ! initialize the pressure pointer to zero
      P_M_ptr_r8(:,:,:,:) = 0.0d0
      P_T_ptr_r8(:,:,:,:) = 0.0d0
    end if
    HeightSfc_ptr_r8 => gsv_getHeightSfc(statevector)

    ! Development notes (@mad001)
    !   if feasible, consider reusing the same code for both 
    !   `calcPressure_{gsv,col}_nl_vcode2100x`
    do_computePressure_gsv_nl: do stepIndex = 1, numStep
      do latIndex = statevector%myLatBeg, statevector%myLatEnd
        do lonIndex = statevector%myLonBeg, statevector%myLonEnd

          pressure_T(:) = 0.0D0
          pressure_M(:) = 0.0D0

          ! latitude
          lat_4 = statevector%hco%lat2d_4(lonIndex,latIndex)
          lat_8 = real(lat_4,8)
          sLat = sin(lat_8)
          cLat = cos(lat_8)

          if ( statevector%dataKind == 4 ) then
            P0 = real(P0_ptr_r4(lonIndex,latIndex,1, stepIndex),8)
          else
            P0 = P0_ptr_r8(lonIndex,latIndex,1, stepIndex)
          end if

          rMT = HeightSfc_ptr_r8(lonIndex,latIndex)

          ! compute pressure on diagnostic levels
          if ( statevector%dataKind == 4 ) then
            hu = real(hu_ptr_r4(lonIndex,latIndex,nlev_T,stepIndex),8)
            tt = real(tt_ptr_r4(lonIndex,latIndex,nlev_T,stepIndex),8)
          else
            hu = hu_ptr_r8(lonIndex,latIndex,nlev_T,stepIndex)
            tt = tt_ptr_r8(lonIndex,latIndex,nlev_T,stepIndex)
          end if
          tv0 = phf_fotvt8(tt,hu)

          ! thermo diagnostic level
          if ( statevector%dataKind == 4 ) then
            Z_T = real(height_T_ptr_r4(lonIndex,latIndex,nlev_T,stepIndex),8)
          else
            Z_T = height_T_ptr_r8(lonIndex,latIndex,nlev_T,stepIndex)
          end if
          cmp = gpscompressibility(P0,tt,hu) 
          tv(nlev_T) = tv0*cmp
          dh = Z_T - rMT
          Rgh = phf_gravityalt(sLat, rMT+0.5D0*dh)
          pressure_T(nlev_T) = P0*exp(-Rgh*dh/MPC_RGAS_DRY_AIR_R8/tv(nlev_T))

          ! momentum diagnostic level
          if ( statevector%dataKind == 4 ) then
            Z_M = real(height_M_ptr_r4(lonIndex,latIndex,nlev_M,stepIndex),8)
          else
            Z_M = height_M_ptr_r8(lonIndex,latIndex,nlev_M,stepIndex)
          end if
          dh = Z_M - rMT
          Rgh = phf_gravityalt(sLat, rMT+0.5D0*dh)
          pressure_M(nlev_M) = P0*exp(-Rgh*dh/MPC_RGAS_DRY_AIR_R8/tv(nlev_T))

          ! compute pressure on all levels above except the last 
          do lev_M = nlev_M-1, 1, -1
            lev_T = lev_M ! thermo level just below
            if ( statevector%dataKind == 4 ) then
              hu = real(hu_ptr_r4(lonIndex,latIndex,lev_T,stepIndex),8)
              tt = real(tt_ptr_r4(lonIndex,latIndex,lev_T,stepIndex),8)
              Z_M = real(height_M_ptr_r4(lonIndex,latIndex,lev_M,stepIndex),8)
              Z_M1 = real(height_M_ptr_r4(lonIndex,latIndex,lev_M+1,stepIndex),8)
              Z_T = real(height_T_ptr_r4(lonIndex,latIndex,lev_T,stepIndex),8)
            else
              hu = hu_ptr_r8(lonIndex,latIndex,lev_T,stepIndex)
              tt = tt_ptr_r8(lonIndex,latIndex,lev_T,stepIndex)
              Z_M = height_M_ptr_r8(lonIndex,latIndex,lev_M,stepIndex)
              Z_M1 = height_M_ptr_r8(lonIndex,latIndex,lev_M+1,stepIndex)
              Z_T = height_T_ptr_r8(lonIndex,latIndex,lev_T,stepIndex)
            end if
            tv0 = phf_fotvt8(tt,hu)
            dh = Z_M - Z_M1
            Rgh = phf_gravityalt(sLat, Z_M1+0.5D0*dh)

            ! approximation of tv from pressure on previous momentum level
            cmp = gpscompressibility(pressure_M(lev_M+1),tt,hu) 
            tv(lev_T) = tv0*cmp
            pressure_M(lev_M) = pressure_M(lev_M+1) * &
                                exp(-Rgh*dh/MPC_RGAS_DRY_AIR_R8/tv(lev_T))
            ! first interpolation of thermo pressure
            scaleFactorBottom = (Z_T-Z_M1)/(Z_M-Z_M1)
            logP = (1.0D0-scaleFactorBottom)*log(pressure_M(lev_M+1)) + &
                                  scaleFactorBottom*log(pressure_M(lev_M))
            pressure_T(lev_T) = exp(logP)

            ! second iteration on tv
            cmp = gpscompressibility(pressure_T(lev_T),tt,hu)
            tv(lev_T) = tv0*cmp
            pressure_M(lev_M) = pressure_M(lev_M+1) * &
                                exp(-Rgh*dh/MPC_RGAS_DRY_AIR_R8/tv(lev_T))

            ! second iteration interpolation of thermo pressure
            logP = (1.0D0-scaleFactorBottom)*log(pressure_M(lev_M+1)) + &
                                  scaleFactorBottom*log(pressure_M(lev_M))
            pressure_T(lev_T) = exp(logP)

          end do

          ! fill the height array
          if ( statevector%dataKind == 4 ) then
            do lev_T = 1, nlev_T
              P_T_ptr_r4(lonIndex,latIndex,lev_T,stepIndex) = &
                  real(pressure_T(lev_T),4)
            end do
            do lev_M = 1, nlev_M
            P_M_ptr_r4(lonIndex,latIndex,lev_M,stepIndex) = &
                  real(pressure_M(lev_M),4)
            end do
          else
            P_T_ptr_r8(lonIndex,latIndex,1:nlev_T,stepIndex)=pressure_T(1:nlev_T)
            P_M_ptr_r8(lonIndex,latIndex,1:nlev_M,stepIndex)=pressure_M(1:nlev_M)
          end if

        end do ! lonIndex
      end do ! latIndex
    end do do_computePressure_gsv_nl

    deallocate(pressure_T)
    deallocate(pressure_M)
    deallocate(tv)

    call msg('calcPressure_gsv_nl_vcode2100x (czp)', 'END', verb_opt=4)
  end subroutine calcPressure_gsv_nl_vcode2100x

  !---------------------------------------------------------
  ! calcPressure_gsv_nl_vcode500x_r8
  !---------------------------------------------------------
  subroutine calcPressure_gsv_nl_vcode500x_r8(statevector, P_T, P_M, Ps_in_hPa_opt)
    !
    ! :Purpose: Pressure retrieval for GEM-P real(8) statevector, values
    !           values returned in pointers.
    !
    implicit none

    ! Arguments
    type(struct_gsv),           intent(in)    :: statevector
    real(8),           pointer, intent(inout) :: P_T(:,:,:,:)
    real(8),           pointer, intent(inout) :: P_M(:,:,:,:)
    logical, optional,          intent(in)    :: Ps_in_hPa_opt            ! If true, conversion from hPa to mbar will be done for surface pressure

    ! Locals
    real(kind=8), allocatable   :: Psfc(:,:)
    real(kind=8), pointer       :: Pressure_out(:,:,:)
    real(kind=8), pointer       :: field_Psfc(:,:,:,:)
    integer                     :: status, stepIndex, numStep

    call msg('calcPressure_gsv_nl_vcode500x_r8 (czp)', 'START', verb_opt=4)

    allocate(Psfc(statevector%myLonBeg:statevector%myLonEnd, &
                  statevector%myLatBeg:statevector%myLatEnd))

    call gsv_getField(statevector,field_Psfc,'P0')
    numStep = statevector%numStep

    do stepIndex = 1, numStep
      Psfc(:,:) = field_Psfc(:,:,1,stepIndex)
      if ( present(Ps_in_hPa_opt) ) then
        if ( Ps_in_hPa_opt ) Psfc = Psfc * mpc_pa_per_mbar_r8
      end if

      ! P_M
      nullify(Pressure_out)
      status = vgd_levels(statevector%vco%vgrid, &
                          ip1_list=statevector%vco%ip1_M, &
                          levels=Pressure_out, &
                          sfc_field=Psfc, &
                          in_log=.false.)
      if( status .ne. VGD_OK ) then
        call utl_abort('calcPressure_gsv_nl_vcode500x_r8 (czp): ERROR with vgd_levels')
      end if
      P_M(:,:,:,stepIndex) = Pressure_out(:,:,:)
      deallocate(Pressure_out)

      ! P_T
      nullify(Pressure_out)
      status = vgd_levels(statevector%vco%vgrid, &
                          ip1_list=statevector%vco%ip1_T, &
                          levels=Pressure_out, &
                          sfc_field=Psfc, &
                          in_log=.false.)
      if( status .ne. VGD_OK ) then
        call utl_abort('calcPressure_gsv_nl_vcode500x_r8 (czp): ERROR with vgd_levels')
      end if
      P_T(:,:,:,stepIndex) = Pressure_out(:,:,:)
      deallocate(Pressure_out)

    end do

    deallocate(Psfc)

    call msg('calcPressure_gsv_nl_vcode500x_r8 (czp)', 'END', verb_opt=4)
  end subroutine calcPressure_gsv_nl_vcode500x_r8

  !---------------------------------------------------------
  ! calcPressure_gsv_nl_vcode500x_r4
  !---------------------------------------------------------
  subroutine calcPressure_gsv_nl_vcode500x_r4(statevector, P_T, P_M, Ps_in_hPa_opt)
    !
    ! :Purpose: Pressure retrieval for GEM-P real(4) statevector, values
    !           values returned in pointers.
    !
    implicit none

    ! Arguments
    type(struct_gsv),           intent(in)    :: statevector
    real(4),           pointer, intent(inout) :: P_T(:,:,:,:)
    real(4),           pointer, intent(inout) :: P_M(:,:,:,:)
    logical, optional,          intent(in)    :: Ps_in_hPa_opt            ! If true, conversion from hPa to mbar will be done for surface pressure

    ! Locals
    real(kind=4), allocatable   :: Psfc(:,:)
    real(kind=4), pointer       :: Pressure_out(:,:,:)
    real(kind=4), pointer       :: field_Psfc(:,:,:,:)
    integer                     :: status, stepIndex, numStep

    call msg('calcPressure_gsv_nl_vcode500x_r4 (czp)', 'START', verb_opt=4)

    allocate(Psfc(statevector%myLonBeg:statevector%myLonEnd, &
                  statevector%myLatBeg:statevector%myLatEnd))

    call gsv_getField(statevector,field_Psfc,'P0')
    numStep = statevector%numStep

    do stepIndex = 1, numStep
      Psfc(:,:) = field_Psfc(:,:,1,stepIndex)
      if ( present(Ps_in_hPa_opt) ) then
        if ( Ps_in_hPa_opt ) Psfc = Psfc * mpc_pa_per_mbar_r4
      end if

      ! P_M
      nullify(Pressure_out)
      status = vgd_levels(statevector%vco%vgrid, &
                          ip1_list=statevector%vco%ip1_M, &
                          levels=Pressure_out, &
                          sfc_field=Psfc, &
                          in_log=.false.)
      if( status .ne. VGD_OK ) then
        call utl_abort('calcPressure_gsv_nl_vcode500x_r4 (czp): ERROR with vgd_levels')
      end if
      P_M(:,:,:,stepIndex) = Pressure_out(:,:,:)
      deallocate(Pressure_out)

      ! P_T
      nullify(Pressure_out)
      status = vgd_levels(statevector%vco%vgrid, &
                          ip1_list=statevector%vco%ip1_T, &
                          levels=Pressure_out, &
                          sfc_field=Psfc, &
                          in_log=.false.)
      if( status .ne. VGD_OK ) then
        call utl_abort('calcPressure_gsv_nl_vcode500x_r4 (czp): ERROR with vgd_levels')
      end if
      P_T(:,:,:,stepIndex) = Pressure_out(:,:,:)
      deallocate(Pressure_out)

    end do

    deallocate(Psfc)

    call msg('calcPressure_gsv_nl_vcode500x_r4 (czp)', 'START', verb_opt=4)
  end subroutine calcPressure_gsv_nl_vcode500x_r4

  !---------------------------------------------------------
  ! calcPressure_gsv_tl
  !---------------------------------------------------------
  subroutine calcPressure_gsv_tl( statevector, statevectorRef)
    !
    !:Purpose: Tangent of pressure computation.
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout) :: statevector      ! statevector that will contain the P_T/P_M increments
    type(struct_gsv), intent(in)    :: statevectorRef   ! statevector containing needed reference fields

    ! Locals
    integer :: Vcode

    call utl_tmg_start(178,'low-level--czp_calcPressure_tl')
    call msg('calcPressure_gsv_tl (czp)', 'START', verb_opt=2)

    Vcode = gsv_getVco(statevectorRef)%vcode
    if (Vcode == 5005 .or. Vcode == 5002) then
      if ( .not. gsv_varExist(statevector,'P_*')  ) then
        call utl_abort('calcPressure_gsv_tl (czp): for vcode 500x, variables P_T and P_M must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'P0')  ) then
        call utl_abort('calcPressure_gsv_tl (czp): for vcode 500x, variable P0 must be allocated in gridstatevector')
      end if
      call calcPressure_gsv_tl_vcode500x
    else if (Vcode == 21001) then
      ! Development notes (@mad001)
      !   probably some some gsv_varExist(statevector,.) needed for GEM-H
      call calcPressure_gsv_tl_vcode2100x
    end if

    call msg('calcPressure_gsv_tl (czp)', 'END', verb_opt=2)
    call utl_tmg_stop(178)

    contains

      !---------------------------------------------------------
      ! calcPressure_gsv_tl_vcode2100x
      !---------------------------------------------------------
      subroutine calcPressure_gsv_tl_vcode2100x
        implicit none

        call utl_abort('calcPressure_gsv_tl (czp): vcode 21001 not implemented yet')

      end subroutine calcPressure_gsv_tl_vcode2100x

      !---------------------------------------------------------
      ! calcPressure_gsv_tl_vcode500x
      !---------------------------------------------------------
      subroutine calcPressure_gsv_tl_vcode500x
        implicit none

        ! Locals
        real(8), allocatable  :: Psfc(:,:)
        real(4), pointer      :: delPsfc_r4(:,:,:,:)
        real(8), pointer      :: delPsfc_r8(:,:,:,:)
        real(8), pointer      :: field_Psfc(:,:,:,:)
        real(4), pointer      :: delP_T_r4(:,:,:,:)
        real(8), pointer      :: delP_T_r8(:,:,:,:)
        real(4), pointer      :: delP_M_r4(:,:,:,:)
        real(8), pointer      :: delP_M_r8(:,:,:,:)
        real(8), pointer      :: dP_dPsfc_T(:,:,:)
        real(8), pointer      :: dP_dPsfc_M(:,:,:)
        integer               :: status, stepIndex,lonIndex,latIndex
        integer               :: lev_M, lev_T, nlev_T, nlev_M, numStep

        call msg('calcPressure_gsv_tl_vcode500x (czp)', 'START', verb_opt=4)

        nullify(dP_dPsfc_T)
        nullify(dP_dPsfc_M)
        nullify(delPsfc_r4,delPsfc_r8)
        nullify(delP_T_r4,delP_T_r8)
        nullify(delP_M_r4,delP_M_r8)

        if (gsv_getDataKind(statevector) == 4) then
          call gsv_getField(statevector,delP_T_r4,'P_T')
          call gsv_getField(statevector,delP_M_r4,'P_M')
          call gsv_getField(statevector,delPsfc_r4,'P0')
        else
          call gsv_getField(statevector,delP_T_r8,'P_T')
          call gsv_getField(statevector,delP_M_r8,'P_M')
          call gsv_getField(statevector,delPsfc_r8,'P0')
        end if

        nullify(field_Psfc)
        call gsv_getField(statevectorRef,field_Psfc,'P0')

        nlev_T = gsv_getNumLev(statevector,'TH')
        nlev_M = gsv_getNumLev(statevector,'MM')
        numStep = statevector%numstep

        allocate(Psfc(statevector%myLonBeg:statevector%myLonEnd, &
                      statevector%myLatBeg:statevector%myLatEnd))

        do stepIndex = 1, numStep

          Psfc(:,:) = field_Psfc(:,:,1,stepIndex)

          ! dP_dPsfc_M
          nullify(dP_dPsfc_M)
          status = vgd_dpidpis(statevector%vco%vgrid, &
                               statevector%vco%ip1_M, &
                               dP_dPsfc_M, &
                               Psfc)
          if( status .ne. VGD_OK ) then
              call utl_abort('calcPressure_gsv_tl (czp): ERROR with vgd_dpidpis')
          end if
          ! calculate delP_M
          if (gsv_getDataKind(statevector) == 4) then
            do lev_M = 1, nlev_M
              do latIndex = statevector%myLatBeg, statevector%myLatEnd
                do lonIndex = statevector%myLonBeg, statevector%myLonEnd
                  delP_M_r4(lonIndex,latIndex,lev_M,stepIndex) =  &
                       dP_dPsfc_M(lonIndex-statevector%myLonBeg+1,&
                                  latIndex-statevector%myLatBeg+1,lev_M) * &
                       delPsfc_r4(lonIndex,latIndex,1,stepIndex)
                end do
              end do
            end do
          else
            do lev_M = 1, nlev_M
              do latIndex = statevector%myLatBeg, statevector%myLatEnd
                do lonIndex = statevector%myLonBeg, statevector%myLonEnd
                  delP_M_r8(lonIndex,latIndex,lev_M,stepIndex) =  &
                       dP_dPsfc_M(lonIndex-statevector%myLonBeg+1,&
                                  latIndex-statevector%myLatBeg+1,lev_M) * &
                       delPsfc_r8(lonIndex,latIndex,1,stepIndex)
                end do
              end do
            end do
          end if
          deallocate(dP_dPsfc_M)

          ! dP_dPsfc_T
          nullify(dP_dPsfc_T)
          status = vgd_dpidpis(statevector%vco%vgrid, &
                               statevector%vco%ip1_T, &
                               dP_dPsfc_T, &
                               Psfc)
          if( status .ne. VGD_OK ) then
              call utl_abort('calcPressure_gsv_tl (czp): ERROR with vgd_dpidpis')
          end if
          ! calculate delP_T
          if (gsv_getDataKind(statevector) == 4) then
            do lev_T = 1, nlev_T
              do latIndex = statevector%myLatBeg, statevector%myLatEnd
                do lonIndex = statevector%myLonBeg, statevector%myLonEnd
                  delP_T_r4(lonIndex,latIndex,lev_T,stepIndex) =  &
                       dP_dPsfc_T(lonIndex-statevector%myLonBeg+1,&
                                  latIndex-statevector%myLatBeg+1,lev_T) * &
                       delPsfc_r4(lonIndex,latIndex,1,stepIndex)
                end do
              end do
            end do
          else
            do lev_T = 1, nlev_T
              do latIndex = statevector%myLatBeg, statevector%myLatEnd
                do lonIndex = statevector%myLonBeg, statevector%myLonEnd
                  delP_T_r8(lonIndex,latIndex,lev_T,stepIndex) =  &
                       dP_dPsfc_T(lonIndex-statevector%myLonBeg+1,&
                                  latIndex-statevector%myLatBeg+1,lev_T) * &
                       delPsfc_r8(lonIndex,latIndex,1,stepIndex)
                end do
              end do
            end do
          end if
          deallocate(dP_dPsfc_T)

        end do

        deallocate(Psfc)

        call msg('calcPressure_gsv_tl_vcode500x (czp)', 'END', verb_opt=4)
      end subroutine calcPressure_gsv_tl_vcode500x

  end subroutine calcPressure_gsv_tl

  !---------------------------------------------------------
  ! calcPressure_gsv_ad
  !---------------------------------------------------------
  subroutine calcPressure_gsv_ad( statevector, statevectorRef)
    !
    !:Purpose: Adjoint of pressure computation.
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout) :: statevector    ! statevector that will contain increment of P_T/P_M
    type(struct_gsv), intent(in)    :: statevectorRef ! statevector containing needed reference fields

    ! Locals
    integer :: Vcode

    call utl_tmg_start(179,'low-level--czp_calcPressure_ad')
    call msg('calcPressure_gsv_ad (czp)', 'START', verb_opt=2)

    Vcode = gsv_getVco(statevectorRef)%vcode
    if (Vcode == 5005 .or. Vcode == 5002) then
      if ( .not. gsv_varExist(statevector,'P_*')  ) then
        call utl_abort('calcPressure_gsv_ad (czp): for vcode 500x, variables P_M and P_T must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'P0')  ) then
        call utl_abort('calcPressure_gsv_ad (czp): for vcode 500x, variable P0 must be allocated in gridstatevector')
      end if
      call calcPressure_gsv_ad_vcode500x
    else if (Vcode == 21001) then
      ! Development notes (@mad001)
      !   probably some some gsv_varExist(statevector,.) needed for GEM-H
      call calcPressure_gsv_ad_vcode2100x
    end if

    call msg('calcPressure_gsv_ad (czp)', 'END', verb_opt=2)
    call utl_tmg_stop(179)

    contains

      !---------------------------------------------------------
      ! calcPressure_gsv_ad_vcode2100x
      !---------------------------------------------------------
      subroutine calcPressure_gsv_ad_vcode2100x
        implicit none

        call utl_abort('calcPressure_gsv_ad (czp): vcode 21001 not implemented yet')

      end subroutine calcPressure_gsv_ad_vcode2100x

      !---------------------------------------------------------
      ! calcPressure_gsv_ad_vcode500x
      !---------------------------------------------------------
      subroutine calcPressure_gsv_ad_vcode500x
        implicit none

        ! Locals
        real(8), allocatable     :: Psfc(:,:)
        real(4), pointer         :: delPsfc_r4(:,:,:,:)
        real(8), pointer         :: delPsfc_r8(:,:,:,:)
        real(8), pointer         :: field_Psfc(:,:,:,:)
        real(4), pointer         :: delP_T_r4(:,:,:,:)
        real(8), pointer         :: delP_T_r8(:,:,:,:)
        real(4), pointer         :: delP_M_r4(:,:,:,:)
        real(8), pointer         :: delP_M_r8(:,:,:,:)
        real(8), pointer         :: dP_dPsfc_T(:,:,:)
        real(8), pointer         :: dP_dPsfc_M(:,:,:)
        integer                  :: status, stepIndex,lonIndex,latIndex
        integer                  :: lev_M, lev_T, nlev_T, nlev_M, numStep

        call msg('calcPressure_gsv_ad_vcode500x (czp)', 'START', verb_opt=4)

        nullify(delPsfc_r4, delPsfc_r8)
        nullify(field_Psfc)
        nullify(delP_T_r4, delP_T_r8)
        nullify(delP_M_r4, delP_M_r8)
        nullify(dP_dPsfc_T)
        nullify(dP_dPsfc_M)

        if (gsv_getDataKind(statevector) == 4) then
          call gsv_getField(statevector,delP_T_r4,'P_T')
          call gsv_getField(statevector,delP_M_r4,'P_M')
          call gsv_getField(statevector,delPsfc_r4,'P0')
        else
          call gsv_getField(statevector,delP_T_r8,'P_T')
          call gsv_getField(statevector,delP_M_r8,'P_M')
          call gsv_getField(statevector,delPsfc_r8,'P0')
        end if
        call gsv_getField(statevectorRef,field_Psfc,'P0')

        nlev_T = gsv_getNumLev(statevector,'TH')
        nlev_M = gsv_getNumLev(statevector,'MM')
        numStep = statevector%numstep

        allocate(Psfc(statevector%myLonBeg:statevector%myLonEnd, &
                      statevector%myLatBeg:statevector%myLatEnd))

        do stepIndex = 1, numStep

          Psfc(:,:) = field_Psfc(:,:,1,stepIndex)

          ! dP_dPsfc_M
          nullify(dP_dPsfc_M)
          status = vgd_dpidpis(statevector%vco%vgrid, &
                               statevector%vco%ip1_M, &
                               dP_dPsfc_M, &
                               Psfc)
          if( status .ne. VGD_OK ) then
            call utl_abort('calcPressure_gsv_ad (czp): ERROR with vgd_dpidpis')
          end if
          ! calculate delP_M
          if (gsv_getDataKind(statevector) == 4) then
            do lev_M = 1, nlev_M
              do latIndex = statevector%myLatBeg, statevector%myLatEnd
                do lonIndex = statevector%myLonBeg, statevector%myLonEnd
                  delPsfc_r4(lonIndex,latIndex,1,stepIndex) =  &
                       delPsfc_r4(lonIndex,latIndex,1,stepIndex) + &
                       dP_dPsfc_M(lonIndex-statevector%myLonBeg+1,&
                                  latIndex-statevector%myLatBeg+1,lev_M) * &
                       delP_M_r4(lonIndex,latIndex,lev_M,stepIndex)
                end do
              end do
            end do
          else
            do lev_M = 1, nlev_M
              do latIndex = statevector%myLatBeg, statevector%myLatEnd
                do lonIndex = statevector%myLonBeg, statevector%myLonEnd
                  delPsfc_r8(lonIndex,latIndex,1,stepIndex) =  &
                       delPsfc_r8(lonIndex,latIndex,1,stepIndex) + &
                       dP_dPsfc_M(lonIndex-statevector%myLonBeg+1,&
                                  latIndex-statevector%myLatBeg+1,lev_M) * &
                       delP_M_r8(lonIndex,latIndex,lev_M,stepIndex)
                end do
              end do
            end do
          end if
          deallocate(dP_dPsfc_M)

          ! dP_dPsfc_T
          nullify(dP_dPsfc_T)
          status = vgd_dpidpis(statevector%vco%vgrid, &
                               statevector%vco%ip1_T, &
                               dP_dPsfc_T, &
                               Psfc)
          if( status .ne. VGD_OK ) then
              call utl_abort('calcPressure_gsv_ad (czp): ERROR with vgd_dpidpis')
          end if
          ! calculate delP_T
          if (gsv_getDataKind(statevector) == 4) then
            do lev_T = 1, nlev_T
              do latIndex = statevector%myLatBeg, statevector%myLatEnd
                do lonIndex = statevector%myLonBeg, statevector%myLonEnd
                  delPsfc_r4(lonIndex,latIndex,1,stepIndex) =  &
                       delPsfc_r4(lonIndex,latIndex,1,stepIndex) + &
                       dP_dPsfc_T(lonIndex-statevector%myLonBeg+1,&
                                  latIndex-statevector%myLatBeg+1,lev_T) * &
                       delP_T_r4(lonIndex,latIndex,lev_T,stepIndex)
                end do
              end do
            end do
          else
            do lev_T = 1, nlev_T
              do latIndex = statevector%myLatBeg, statevector%myLatEnd
                do lonIndex = statevector%myLonBeg, statevector%myLonEnd
                  delPsfc_r8(lonIndex,latIndex,1,stepIndex) =  &
                       delPsfc_r8(lonIndex,latIndex,1,stepIndex) + &
                       dP_dPsfc_T(lonIndex-statevector%myLonBeg+1,&
                                  latIndex-statevector%myLatBeg+1,lev_T) * &
                       delP_T_r8(lonIndex,latIndex,lev_T,stepIndex)
                end do
              end do
            end do
          end if
          deallocate(dP_dPsfc_T)

        end do

        deallocate(Psfc)

        call msg('calcPressure_gsv_ad_vcode500x (czp)', 'END', verb_opt=4)
      end subroutine calcPressure_gsv_ad_vcode500x

  end subroutine calcPressure_gsv_ad

  !---------------------------------------------------------------------
  ! subroutines operating on struct_columnData
  !---------------------------------------------------------------------

  !---------------------------------------------------------
  ! calcZandP_col_nl
  !---------------------------------------------------------
  subroutine calcZandP_col_nl(column)
    !
    ! :Purpose: Compute pressure and height in the column in proper order
    !           depending on the vgrid kind
    !
    implicit none

    ! Arguments
    type(struct_columnData), intent(inout) :: column  ! column that will contain the Z_*/P_* fields

    ! Locals
    integer   :: Vcode

    call msg('calcZandP_col_nl (czp)', 'START', verb_opt=2)

    Vcode = column%vco%vcode
    if (Vcode == 5002 .or. Vcode == 5005) then
      ! if P_T, P_M not allocated : do nothing
      if (col_varExist(column,'P_*')) then
        call calcPressure_col_nl(column)
        if (col_varExist(column,'Z_*')) then
          call calcHeight_col_nl(column)
        end if
      end if
    else if (Vcode == 21001) then
      ! if Z_T, Z_M not allocated : do nothing
      if (col_varExist(column,'Z_*')) then
        call calcHeight_col_nl(column)
        if (col_varExist(column,'P_*')) then
          call calcPressure_col_nl(column)
        end if
      end if
    end if
  
    call msg('calcZandP_col_nl (czp)', 'END', verb_opt=2)
  end subroutine calcZandP_col_nl

  !---------------------------------------------------------
  ! calcZandP_col_tl
  !---------------------------------------------------------
  subroutine calcZandP_col_tl(columnInc, columnIncRef)
    !
    ! :Purpose: Compute pressure and height increment in the column in proper
    !           order depending on the vgrid kind
    !
    implicit none

    ! Arguments
    type(struct_columnData), intent(inout) :: columnInc    ! column that will contain the Z_*/P_* increments
    type(struct_columnData), intent(in)    :: columnIncRef ! column containing needed reference fields

    ! Locals
    integer   :: Vcode

    call msg('calcZandP_col_tl (czp)', 'START', verb_opt=2)

    Vcode = columnInc%vco%vcode
    if (Vcode == 5002 .or. Vcode == 5005) then
      ! if P_T, P_M not allocated : do nothing
      if (col_varExist(columnInc,'P_*')) then
        call calcPressure_col_tl( columnInc, columnIncRef)
        if (col_varExist(columnInc,'Z_*')) then
          call calcHeight_col_tl(columnInc, columnIncRef)
        end if
      end if
    else if (Vcode == 21001) then
      ! if Z_T, Z_M not allocated : do nothing
      if (col_varExist(columnInc,'Z_*')) then
        call calcHeight_col_tl(columnInc, columnIncRef)

        if (col_varExist(columnInc,'P_*')) then
          call calcPressure_col_tl( columnInc, columnIncRef)
        end if
      end if
    end if

    call msg('calcZandP_col_tl (czp)', 'END', verb_opt=2)
  end subroutine calcZandP_col_tl

  !---------------------------------------------------------
  ! calcZandP_col_ad
  !---------------------------------------------------------
  subroutine calcZandP_col_ad(columnInc, columnIncRef)
    !
    ! :Purpose: Adjoint of pressure and height increment computation in the
    !           column in proper order depending on the vgrid kind
    !
    implicit none

    ! Arguments
    type(struct_columnData), intent(inout) :: columnInc    ! column that will contain the Z_*/P_* increments
    type(struct_columnData), intent(in)    :: columnIncRef ! column containing needed reference fields

    ! Locals
    integer   :: Vcode


    call msg('calcZandP_col_ad (czp)', 'START', verb_opt=2)

    Vcode = columnInc%vco%vcode
    if (Vcode == 5002 .or. Vcode == 5005) then
      ! if Z_T, Z_M not allocated : do nothing
      if (col_varExist(columnInc,'Z_*')) then
        call calcHeight_col_ad(columnInc, columnIncRef)
        if (col_varExist(columnInc,'P_*')) then
          call calcPressure_col_ad( columnInc, columnIncRef)
        end if
      end if
    else if (Vcode == 21001) then
      ! if P_T, P_M not allocated : do nothing
      if (col_varExist(columnInc,'P_*')) then
        call calcPressure_col_ad( columnInc, columnIncRef)
        if (col_varExist(columnInc,'Z_*')) then
          call calcHeight_col_ad(columnInc, columnIncRef)
        end if
      end if
    end if

    call msg('calcZandP_col_ad (czp)', 'END', verb_opt=2)
  end subroutine calcZandP_col_ad

  !---------------------------------------------------------
  ! calcHeight_col_nl
  !---------------------------------------------------------
  subroutine calcHeight_col_nl(column)
    !
    ! :Purpose: Compute or retrieve heights on the column.
    !
    implicit none

    ! Arguments
    type(struct_columnData), intent(inout) :: column  ! column that will contain the Z_M/Z_T fields

    ! Locals
    real(8), pointer  ::  Z_T(:,:), Z_M(:,:)

    call msg('calcHeight_col_nl (czp)', 'START', verb_opt=2)

    Z_T => col_getAllColumns(column, 'Z_T')
    Z_M => col_getAllColumns(column, 'Z_M')
    call czp_calcReturnHeight_col_nl(column, Z_T, Z_M)

    call msg('calcHeight_col_nl (czp)', &
           new_line('')//'Z_M = '//str(col_getColumn(column,1,'Z_M')) &
         //new_line('')//'Z_T = '//str(col_getColumn(column,1,'Z_T')), &
         verb_opt=2)

    call msg('calcHeight_col_nl (czp)', 'END', verb_opt=2)
  end subroutine calcHeight_col_nl

  !---------------------------------------------------------
  ! czp_calcReturnHeight_col_nl
  !---------------------------------------------------------
  subroutine czp_calcReturnHeight_col_nl(column, Z_T, Z_M)
    !
    ! :Purpose: Return heights on the column.
    !
    implicit none

    ! Arguments
    type(struct_columnData),  intent(in)    :: column  ! reference column containing temperature and geopotential
    real(8), pointer,         intent(inout) :: Z_T(:,:) ! output pointer to computed column height values on thermodynamic levels
    real(8), pointer,         intent(inout) :: Z_M(:,:) ! output pointer to computed column height values on momentum levels

    ! Locals
    integer :: vcode

    call msg('czp_calcReturnHeight_col_nl (czp)', 'START', verb_opt=2)

    Vcode = col_getVco(column)%vcode
    if (Vcode == 5005 .or. Vcode == 5002) then
      if ( .not. (col_varExist(column,'P0') .and. col_varExist(column,'TT') .and. &
                  col_varExist(column,'HU'))  ) then
        call utl_abort('czp_calcReturnHeight_col_nl: for vcode 500x, variables P0, TT and HU must be allocated in column')
      end if
      call calcHeight_col_nl_vcode500x(column, Z_T, Z_M)
    else if (Vcode == 21001) then
      call calcHeight_col_nl_vcode2100x(column, Z_T, Z_M)

    end if

    call msg('czp_calcReturnHeight_col_nl (czp)', 'END', verb_opt=2)
  end subroutine czp_calcReturnHeight_col_nl

  !---------------------------------------------------------
  ! calcHeight_col_nl_vcode2100x
  !---------------------------------------------------------
  subroutine calcHeight_col_nl_vcode2100x(column, Z_T, Z_M)
    !
    ! :Purpose: Return heights on a GEM-H column.
    !
    implicit none

    ! Arguments
    type(struct_columnData),  intent(in)    :: column  ! reference column containing temperature and geopotential
    real(8), pointer,         intent(inout) :: Z_T(:,:) ! output pointer to computed column height values on thermodynamic levels
    real(8), pointer,         intent(inout) :: Z_M(:,:) ! output pointer to computed column height values on momentum levels

    ! Locals
    real(8), allocatable  :: hSfc(:,:)
    real(8), pointer      :: hPtr(:,:,:)
    integer :: numCol, colIndex, stat

    call msg('calcHeight_col_nl_vcode2100x (czp)', 'START', verb_opt=4)
    if ( col_getNumCol(column) <= 0 ) then
      call msg('calcHeight_col_nl_vcode2100x (czp)',&
           'END (number of columns <= 0)', verb_opt=2)
      return
    end if

    numCol = col_getNumCol(column)
    allocate(hSfc(1, numCol))
    do colIndex = 1, numCol
      hSfc(1,colIndex) = col_getHeight(column,1,colIndex, 'SF')
    end do

    ! momentum levels
    stat=vgd_levels(column%vco%vgrid, ip1_list=column%vco%ip1_M, &
                    levels=hPtr, sfc_field=hSfc, in_log=.false.)
    if ( stat .ne. VGD_OK ) then
      call utl_abort('calcHeight_col_nl_vcode2100x (czp): ERROR with vgd_levels')
    end if
    Z_M(:,:) = hPtr(1,:,:)
    if (associated(hPtr)) deallocate(hPtr)

    ! thermo levels
    stat=vgd_levels(column%vco%vgrid, ip1_list=column%vco%ip1_T, &
                    levels=hPtr, sfc_field=hSfc, in_log=.false.)
    if ( stat .ne. VGD_OK ) then
      call utl_abort('calcHeight_col_nl_vcode2100x (czp): ERROR with vgd_levels')
    end if
    Z_T(:,:) = hPtr(1,:,:)
    if (associated(hPtr)) deallocate(hPtr)

    deallocate(hSfc)
    call msg('calcHeight_col_nl_vcode2100x (czp)', 'END', verb_opt=4)
  end subroutine calcHeight_col_nl_vcode2100x

  !---------------------------------------------------------
  ! calcHeight_col_nl_vcode500x
  !---------------------------------------------------------
  subroutine calcHeight_col_nl_vcode500x(column, Z_T, Z_M)
    !
    ! :Purpose: Compute heights for GEM-P columns, return height values 
    !           in pointer arguments.
    !
    implicit none

    ! Arguments
    type(struct_columnData),  intent(in)    :: column  ! reference column containing temperature and geopotential
    real(8), pointer,         intent(inout) :: Z_T(:,:) ! output pointer to computed column height values on thermodynamic levels
    real(8), pointer,         intent(inout) :: Z_M(:,:) ! output pointer to computed column height values on momentum levels

    ! Developement notes (@mad001)
    !   Null subroutine, no computation needed at time of writing.
    !   The code is traversed because of `calcZandP_nl` call in `cvt` (agnostic if
    !   dealing with GEM-P or GEM-H), but the results for heights are not used 
    !   at this time.
    !   We keep that stub however for future functionalities.
    call msg('calcHeight_col_nl_vcode500x (czp)', 'END (nothing done)', verb_opt=4)
    return

    ! to prevent 'variable not used' remark
    if (.false.) then
      call col_zero(column)
      Z_T = 0.0
      Z_M = 0.0
    end if

  end subroutine calcHeight_col_nl_vcode500x


  !---------------------------------------------------------
  ! calcHeight_col_tl
  !---------------------------------------------------------
  subroutine calcHeight_col_tl(columnInc,columnIncRef)
    !
    ! :Purpose: Tangent height computation on the column.
    !
    implicit none

    ! Arguments
    type(struct_columnData), intent(inout)  :: columnInc
    type(struct_columnData), intent(in)     :: columnIncRef

    ! Locals
    integer :: Vcode

    call utl_tmg_start(173,'low-level--czp_calcHeight_tl')
    call msg('calcHeight_col_tl (czp)', 'START', verb_opt=2)

    Vcode = col_getVco(columnIncRef)%vcode
    if (Vcode == 5005 .or. Vcode == 5002) then
      if ( .not. col_varExist(columnInc,'P_*')  ) then
        call utl_abort('calcHeight_col_tl (czp): for vcode 500x, variables P_M and P_T must be allocated in column')
      end if
      if ( .not. col_varExist(columnInc,'Z_*')  ) then
        call utl_abort('calcHeight_col_tl (czp): for vcode 500x, variables Z_M and Z_T must be allocated in column')
      end if
      if ( .not. col_varExist(columnInc,'TT')  ) then
        call utl_abort('calcHeight_col_tl (czp): for vcode 500x, variable TT must be allocated in column')
      end if
      if ( .not. col_varExist(columnInc,'HU')  ) then
        call utl_abort('calcHeight_col_tl (czp): for vcode 500x, variable HU must be allocated in column')
      end if
      if ( .not. col_varExist(columnInc,'P0')  ) then
        call utl_abort('calcHeight_col_tl (czp): for vcode 500x, variable P0 must be allocated in column')
      end if
      call calcHeight_col_tl_vcode500x
    else if (Vcode == 21001) then
      ! Development notes (@mad001)
      !   probably some some gsv_varExist(statevector,.) needed for GEM-H
      call calcHeight_col_tl_vcode2100x
    end if

    call msg('calcHeight_col_tl (czp)', 'END', verb_opt=2)
    call utl_tmg_stop(173)

    contains
      !---------------------------------------------------------
      ! calcHeight_col_tl_vcode2100x
      !---------------------------------------------------------
      subroutine calcHeight_col_tl_vcode2100x
        implicit none

        call utl_abort('calcHeight_col_tl: vcode 21001 not implemented yet')

      end subroutine calcHeight_col_tl_vcode2100x

      !---------------------------------------------------------
      ! calcHeight_col_tl_vcode500x
      !---------------------------------------------------------
      subroutine calcHeight_col_tl_vcode500x
        implicit none

        ! Locals
        integer :: lev_M,lev_T,nlev_M,nlev_T,colIndex,numColumns
        real(8) :: ScaleFactorBottom, ScaleFactorTop
        real(8), allocatable :: delThick(:,:)
        real(8), pointer     :: height_T_ptr(:,:),height_M_ptr(:,:)
        real(8), pointer     :: P_T(:,:), P_M(:,:)
        real(8), pointer  :: delHeight_M_ptr(:,:),delHeight_T_ptr(:,:)
        real(8), pointer  :: delTT(:,:),delHU(:,:),delP0(:,:)
        real(8), pointer  :: delP_T(:,:), delP_M(:,:)

        call msg('calcHeight_col_tl_vcode500x (czp)', 'START', verb_opt=4)

        nlev_T = col_getNumLev(columnIncRef,'TH')
        nlev_M = col_getNumLev(columnIncRef,'MM')

        numColumns = col_getNumCol(columnInc)

        allocate(delThick(nlev_T,numColumns))
        delThick(:,:) = 0.0d0

        ! generate the height coefficients on the grid
        call calcHeightCoeff_col(columnIncRef)

        ! loop over all lat/lon/step

        height_M_ptr => col_getAllColumns(columnIncRef,'Z_M')
        height_T_ptr => col_getAllColumns(columnIncRef,'Z_T')
        P_M          => col_getAllColumns(columnIncRef,'P_M')
        P_T          => col_getAllColumns(columnIncRef,'P_T')

        delHeight_M_ptr => col_getAllColumns(columnInc,'Z_M')
        delHeight_T_ptr => col_getAllColumns(columnInc,'Z_T')
        delTT           => col_getAllColumns(columnInc,'TT')
        delHU           => col_getAllColumns(columnInc,'HU')
        delP0           => col_getAllColumns(columnInc,'P0')
        delP_M          => col_getAllColumns(columnInc,'P_M')
        delP_T          => col_getAllColumns(columnInc,'P_T')

        ! ensure increment at sfc is zero (fixed height level)
        delHeight_M_ptr(nlev_M,:) = 0.0d0
        delHeight_T_ptr(nlev_T,:) = 0.0d0

        if_computeHeight_col_tl_vcodes : if (Vcode == 5002) then

          ! compute increment to thickness for each layer between the two momentum levels
          do colIndex = 1, numColumns
            do lev_T = 2, (nlev_T-1)
              delThick(lev_T,colIndex) =  &
                    coeff_M_TT_col(lev_T,colIndex) * delTT(lev_T,colIndex) + &
                    coeff_M_HU_col(lev_T,colIndex) * delHU(lev_T,colIndex) + &
                    coeff_M_P0_delPM_col(lev_T,colIndex) * &
                    ( delP_M(lev_T  ,colIndex) / P_M(lev_T  ,colIndex) - &
                      delP_M(lev_T-1,colIndex) / P_M(lev_T-1,colIndex) ) + &
                    coeff_M_P0_dP_delPT_col(lev_T,colIndex) * &
                    delP_T(lev_T,colIndex) + &
                    coeff_M_P0_dP_delP0_col(lev_T,colIndex) * delP0(1,colIndex)
            end do
          end do

          ! compute height increment on momentum levels above the surface
          do colIndex = 1, numColumns
            do lev_M = (nlev_M-1), 1, -1
              lev_T = lev_M + 1 ! thermo level just below momentum level being computed
              delHeight_M_ptr(lev_M,colIndex) =  &
                   delHeight_M_ptr(lev_M+1,colIndex) + delThick(lev_T,colIndex)
            end do
          end do

          ! compute height increment on thermo levels using weighted average of height increment of momentum levels
          do colIndex = 1, numColumns
            do lev_T = 1, (nlev_T-1)
              if ( lev_T == 1) then
                ! compute height increment for top thermo level (from top momentum level)
                delHeight_T_ptr(1,colIndex) = delHeight_M_ptr(1,colIndex) +  &
                     coeff_T_TT_col(colIndex) * delTT(1,colIndex) + &
                     coeff_T_HU_col(colIndex) * delHU(1,colIndex) + &
                     coeff_T_P0_delP1_col(colIndex) * &
                     ( delP_M(1,colIndex) / P_M(1,colIndex) - &
                       delP_T(1,colIndex) / P_T(1,colIndex) ) + &
                     coeff_T_P0_dP_delPT_col(colIndex) * delP_T(1,colIndex) + &
                     coeff_T_P0_dP_delP0_col(colIndex) * delP0(1,colIndex)
              else
                lev_M = lev_T ! momentum level just below thermo level being computed
                ScaleFactorBottom = (height_T_ptr(lev_T,colIndex) - &
                    height_M_ptr(lev_M-1,colIndex)) / &
                    (height_M_ptr(lev_M,colIndex) - height_M_ptr(lev_M-1,colIndex))
                ScaleFactorTop    = 1 - ScaleFactorBottom
                delHeight_T_ptr(lev_T,colIndex) =  &
                     ScaleFactorBottom * delHeight_M_ptr(lev_M  ,colIndex) + &
                     ScaleFactorTop * delHeight_M_ptr(lev_M-1,colIndex)
              end if
            end do
          end do

        else if(Vcode == 5005) then if_computeHeight_col_tl_vcodes

          ! compute increment to thickness for each layer between the two momentum levels
          do colIndex = 1, numColumns
            do lev_T = 1, (nlev_T-1)
              delThick(lev_T,colIndex) =  &
                   coeff_M_TT_col(lev_T,colIndex) * delTT(lev_T,colIndex) + &
                   coeff_M_HU_col(lev_T,colIndex) * delHU(lev_T,colIndex) + &
                   coeff_M_P0_delPM_col(lev_T,colIndex) * &
                   ( delP_M(lev_T+1,colIndex) / P_M(lev_T+1,colIndex) - &
                     delP_M(lev_T  ,colIndex) / P_M(lev_T  ,colIndex) ) + &
                   coeff_M_P0_dP_delPT_col(lev_T,colIndex) * &
                   delP_T(lev_T,colIndex) + &
                   coeff_M_P0_dP_delP0_col(lev_T,colIndex) * delP0(1,colIndex)
            end do
          end do

          ! compute height increment on momentum levels above the surface
          do colIndex = 1, numColumns
            do lev_M = (nlev_M-1), 1, -1
              lev_T = lev_M ! thermo level just below momentum level being computed
              delHeight_M_ptr(lev_M,colIndex) =  &
                   delHeight_M_ptr(lev_M+1,colIndex) + delThick(lev_T,colIndex)
            end do
          end do

          ! compute height increment on thermo levels using weighted average of height increment of momentum levels
          do colIndex = 1, numColumns
            do lev_T = 1, (nlev_T-1)
              lev_M = lev_T + 1 ! momentum level just below thermo level being computed
              ScaleFactorBottom =  &
                   (height_T_ptr(lev_T,colIndex) - height_M_ptr(lev_M-1,colIndex)) / &
                   (height_M_ptr(lev_M,colIndex) - height_M_ptr(lev_M-1,colIndex))
              ScaleFactorTop    = 1 - ScaleFactorBottom
              delHeight_T_ptr(lev_T,colIndex) =  &
                   ScaleFactorBottom * delHeight_M_ptr(lev_M  ,colIndex) + &
                   ScaleFactorTop * delHeight_M_ptr(lev_M-1,colIndex)
            end do
          end do
        end if if_computeHeight_col_tl_vcodes

        deallocate(delThick)

        call msg('calcHeight_col_tl_vcode500x (czp)', 'END', verb_opt=4)
      end subroutine calcHeight_col_tl_vcode500x

  end subroutine calcHeight_col_tl

  !---------------------------------------------------------
  ! calcHeight_col_ad
  !---------------------------------------------------------
  subroutine calcHeight_col_ad(columnInc,columnIncRef)
    !
    ! :Purpose: Adjoint of height computation on the column.
    !
    !
    implicit none

    ! Arguments
    type(struct_columnData), intent(inout) :: columnInc
    type(struct_columnData), intent(in) :: columnIncRef

    ! Locals
    integer :: Vcode

    call utl_tmg_start(174,'low-level--czp_calcHeight_ad')
    call msg('calcHeight_col_ad (czp)', 'START', verb_opt=2)

    Vcode = col_getVco(columnIncRef)%vcode
    if (Vcode == 5005 .or. Vcode == 5002) then
      if ( .not. col_varExist(columnInc,'P_*')  ) then
        call utl_abort('calcHeight_col_ad (czp): for vcode 500x, variables P_M and P_T must be allocated in column')
      end if
      if ( .not. col_varExist(columnInc,'Z_*')  ) then
        call utl_abort('calcHeight_col_ad (czp): for vcode 500x, variables Z_M and Z_T must be allocated in column')
      end if
      if ( .not. col_varExist(columnInc,'TT')  ) then
        call utl_abort('calcHeight_col_ad (czp): for vcode 500x, variable TT must be allocated in column')
      end if
      if ( .not. col_varExist(columnInc,'HU')  ) then
        call utl_abort('calcHeight_col_ad (czp): for vcode 500x, variable HU must be allocated in column')
      end if
      if ( .not. col_varExist(columnInc,'P0')  ) then
        call utl_abort('calcHeight_col_ad (czp): for vcode 500x, variable P0 must be allocated in column')
      end if
      call calcHeight_col_ad_vcode500x
    else if (Vcode == 21001) then
      ! Development notes (@mad001)
      !   probably some some gsv_varExist(statevector,.) needed for GEM-H
      call calcHeight_col_ad_vcode2100x
    end if

    call msg('calcHeight_col_ad (czp)', 'END', verb_opt=2)
    call utl_tmg_stop(174)

    contains
      !---------------------------------------------------------
      ! calcHeight_col_ad_vcode2100x
      !---------------------------------------------------------
      subroutine calcHeight_col_ad_vcode2100x
        implicit none

        call utl_abort('calcHeight_col_ad (czp): vcode 21001 not implemented yet')

      end subroutine calcHeight_col_ad_vcode2100x

      !---------------------------------------------------------
      ! calcHeight_col_ad_vcode500x
      !---------------------------------------------------------
      subroutine calcHeight_col_ad_vcode500x
        implicit none

        ! Locals
        integer :: lev_M,lev_T,nlev_M,nlev_T,numColumns,colIndex
        real(8) :: ScaleFactorBottom, ScaleFactorTop
        real(8), allocatable :: delThick(:,:)
        real(8), pointer     :: height_M_ptr(:,:),height_T_ptr(:,:)
        real(8), allocatable :: delHeight_M(:,:),delHeight_T(:,:)
        real(8), pointer     :: P_M(:,:),P_T(:,:)
        real(8), pointer     :: delHeight_M_ptr(:,:),delHeight_T_ptr(:,:)
        real(8), pointer     :: delTT(:,:),delHU(:,:),delP0(:,:)
        real(8), pointer     :: delP_M(:,:),delP_T(:,:)

        call msg('calcHeight_col_ad_vcode500x (czp)', 'START', verb_opt=4)

        nlev_T = col_getNumLev(columnIncRef,'TH')
        nlev_M = col_getNumLev(columnIncRef,'MM')
        numColumns = col_getNumCol(columnIncRef)

        allocate(delHeight_M(nlev_M,numColumns))
        allocate(delHeight_T(nlev_T,numColumns))
        allocate(delThick(0:nlev_T,numColumns))

        ! generate the height coefficients on the grid
        call calcHeightCoeff_col(columnIncRef)

        height_M_ptr => col_getAllColumns(columnIncRef,'Z_M')
        height_T_ptr => col_getAllColumns(columnIncRef,'Z_T')
        P_M          => col_getAllColumns(columnIncRef,'P_M')
        P_T          => col_getAllColumns(columnIncRef,'P_T')

        delHeight_M_ptr => col_getAllColumns(columnInc,'Z_M')
        delHeight_T_ptr => col_getAllColumns(columnInc,'Z_T')
        delTT           => col_getAllColumns(columnInc,'TT')
        delHU           => col_getAllColumns(columnInc,'HU')
        delP0           => col_getAllColumns(columnInc,'P0')
        delP_M          => col_getAllColumns(columnInc,'P_M')
        delP_T          => col_getAllColumns(columnInc,'P_T')

        delHeight_M(:,:) = delHeight_M_ptr(:,:)
        delHeight_T(:,:) = delHeight_T_ptr(:,:)

        if_computeHeight_col_ad_vcodes : if(Vcode == 5002) then

          ! adjoint of compute height increment on thermo levels by simple averaging
          do colIndex = 1, numColumns
            do lev_T = 1, (nlev_T-1)
              lev_M = lev_T ! momentum level just below thermo level being computed

              ! adjoint of compute height increment on top thermo level (from top momentum level)
              if (lev_T == 1) then
                delHeight_M(1,colIndex)  =  &
                     delHeight_M(1,colIndex) + &
                     delHeight_T(1,colIndex)

                delTT(1,colIndex) =  &
                     delTT(1,colIndex) + &
                     coeff_T_TT_col(colIndex) * delHeight_T(1,colIndex)

                delHU(1,colIndex) =  &
                     delHU(1,colIndex) + &
                     coeff_T_HU_col   (colIndex) * delHeight_T(1,colIndex)

                delP_M(1,colIndex) =  &
                     delP_M(1,colIndex) + &
                     coeff_T_P0_delP1_col(colIndex) / P_M(1,colIndex) * &
                     delHeight_T(1,colIndex)

                delP_T(1,colIndex) =  &
                     delP_T(1,colIndex) - &
                     coeff_T_P0_delP1_col(colIndex) / P_T(1,colIndex) * &
                     delHeight_T(1,colIndex)

                delP_T(1,colIndex) =  &
                     delP_T(1,colIndex) + &
                     coeff_T_P0_dP_delPT_col(colIndex) * delHeight_T(1,colIndex)

                delP0(1,colIndex) =  &
                     delP0(1,colIndex) + &
                     coeff_T_P0_dp_delP0_col(colIndex) * delHeight_T(1,colIndex)
              else
                ScaleFactorBottom =  &
                    (height_T_ptr(lev_T,colIndex) - &
                      height_M_ptr(lev_M-1,colIndex)) / &
                    (height_M_ptr(lev_M,colIndex) - &
                      height_M_ptr(lev_M-1,colIndex))
                ScaleFactorTop    = 1 - ScaleFactorBottom

                delHeight_M(lev_M-1,colIndex) =  &
                     delHeight_M(lev_M-1,colIndex) + &
                     ScaleFactorTop * delHeight_T(lev_T,colIndex)

                delHeight_M(lev_M,colIndex) =  &
                     delHeight_M(lev_M  ,colIndex) + &
                     ScaleFactorBottom * delHeight_T(lev_T,colIndex)
              end if
            end do
          end do

          ! adjoint of compute height increment on momentum levels above the surface
          delThick(0:1,:) = 0.0d0
          do colIndex = 1, numColumns
            do lev_M = 1, (nlev_M-1)
              lev_T = lev_M + 1 ! thermo level just below momentum level being computed
              delThick(lev_T,colIndex) =  &
                   delThick(lev_T-1,colIndex) + &
                   delHeight_M(lev_M  ,colIndex)
            end do
          end do

          ! adjoint of compute increment to thickness for each layer between the two momentum levels
          do colIndex = 1, numColumns
            do lev_T = 2, nlev_T-1
              delTT(lev_T,colIndex) =  &
                  delTT            (lev_T,colIndex) + &
                  coeff_M_TT_col   (lev_T,colIndex) * delThick(lev_T,colIndex)

              delHU(lev_T,colIndex) =  &
                  delHU            (lev_T,colIndex) + &
                  coeff_M_HU_col   (lev_T,colIndex) * delThick(lev_T,colIndex)

              delP_M(lev_T,colIndex)=  &
                  delP_M(lev_T,colIndex) + &
                  coeff_M_P0_delPM_col(lev_T,colIndex) / P_M(lev_T,colIndex) * &
                  delThick(lev_T,colIndex)

              delP_M(lev_T-1,colIndex) =  &
                  delP_M(lev_T-1,colIndex) - &
                  coeff_M_P0_delPM_col(lev_T,colIndex) / P_M(lev_T-1,colIndex)*&
                  delThick(lev_T,colIndex)

              delP_T(lev_T,colIndex) =  &
                  delP_T(lev_T,colIndex) + &
                  coeff_M_P0_dP_delPT_col(lev_T,colIndex) * &
                  delThick(lev_T,colIndex)

              delP0(1,colIndex) =  &
                  delP0(1,colIndex) + &
                  coeff_M_P0_dP_delP0_col(lev_T,colIndex) * &
                  delThick(lev_T,colIndex)
            end do
          end do

        else if(Vcode == 5005) then if_computeHeight_col_ad_vcodes

          ! adjoint of compute height increment on thermo levels by simple averaging
          do colIndex = 1, numColumns
            do lev_T = 1, (nlev_T-1)
              lev_M = lev_T+1 ! momentum level just below thermo level being computed
              ScaleFactorBottom = (height_T_ptr(lev_T,colIndex) - &
                    height_M_ptr(lev_M-1,colIndex)) / &
                  (height_M_ptr(lev_M,colIndex) - &
                    height_M_ptr(lev_M-1,colIndex))
              ScaleFactorTop    = 1 - ScaleFactorBottom
              delHeight_M(lev_M-1,colIndex) = delHeight_M(lev_M-1,colIndex) + &
                  ScaleFactorTop * delHeight_T(lev_T  ,colIndex)
              delHeight_M(lev_M,colIndex)   = delHeight_M(lev_M  ,colIndex) + &
                  ScaleFactorBottom * delHeight_T(lev_T  ,colIndex)
            end do
          end do

          ! adjoint of compute height increment on momentum levels
          delThick(0,:) = 0.0d0
          do colIndex = 1, numColumns
            do lev_M = 1, (nlev_M-1)
              lev_T = lev_M ! thermo level just below momentum level being computed
              delThick(lev_T,colIndex) = delThick(lev_T-1,colIndex) + &
                                         delHeight_M (lev_M  ,colIndex)
            end do
          end do

          do colIndex = 1, numColumns
            do lev_T = 1, nlev_T-1
              delTT(lev_T,colIndex) =  &
                  delTT(lev_T,colIndex) + &
                  coeff_M_TT_col(lev_T,colIndex) * delThick(lev_T,colIndex)

              delHU(lev_T,colIndex) =  &
                  delHU(lev_T,colIndex) + &
                  coeff_M_HU_col(lev_T,colIndex) * delThick(lev_T,colIndex)

              delP_M(lev_T+1,colIndex) =  &
                  delP_M(lev_T+1,colIndex) + &
                  coeff_M_P0_delPM_col(lev_T,colIndex) / P_M(lev_T+1,colIndex)*&
                  delThick(lev_T,colIndex)

              delP_M(lev_T  ,colIndex) =  &
                  delP_M(lev_T  ,colIndex) - &
                  coeff_M_P0_delPM_col(lev_T,colIndex) / P_M(lev_T  ,colIndex)*&
                  delThick(lev_T,colIndex)

              delP_T(lev_T  ,colIndex) =  &
                  delP_T(lev_T  ,colIndex) + &
                  coeff_M_P0_dP_delPT_col(lev_T,colIndex) * &
                  delThick(lev_T,colIndex)

              delP0(1,colIndex)     =  &
                  delP0(1,colIndex) + &
                  coeff_M_P0_dP_delP0_col(lev_T,colIndex) * &
                  delThick(lev_T,colIndex)
            end do
          end do
        end if if_computeHeight_col_ad_vcodes

        deallocate(delThick)
        deallocate(delHeight_M)
        deallocate(delHeight_T)

        call msg('calcHeight_col_ad_vcode500x (czp)', 'END', verb_opt=4)
      end subroutine calcHeight_col_ad_vcode500x

  end subroutine calcHeight_col_ad

  !---------------------------------------------------------
  ! calcPressure_col_nl
  !---------------------------------------------------------
  subroutine calcPressure_col_nl(column)
    !
    ! :Purpose: Pressure computation on the column.
    !
    implicit none

    ! Arguments
    type(struct_columnData), intent(inout)  :: column

    ! Locals
    real(8), pointer  ::  P_T(:,:), P_M(:,:)

    call msg('calcPressure_col_nl (czp)', 'START', verb_opt=2)

    P_T => col_getAllColumns(column, 'P_T')
    P_M => col_getAllColumns(column, 'P_M')
    call czp_calcReturnPressure_col_nl(column, P_T, P_M)

    call msg('calcPressure_col_nl (czp)', &
           new_line('')//'P_M = '//str(col_getColumn(column,1,'P_M')) &
         //new_line('')//'P_T = '//str(col_getColumn(column,1,'P_T')), &
         verb_opt=2)

    call msg('calcPressure_col_nl (czp)', 'END', verb_opt=2)
  end subroutine calcPressure_col_nl

  !---------------------------------------------------------
  ! czp_calcReturnPressure_col_nl
  !---------------------------------------------------------
  subroutine czp_calcReturnPressure_col_nl(column, P_T, P_M)
    !
    ! :Purpose: Pressure computation on the column, return values in pointers.
    !
    implicit none

    ! Arguments
    type(struct_columnData),  intent(in)    :: column ! reference column
    real(8), pointer,         intent(inout) :: P_T(:,:) ! output pointer to computed column pressure values on thermodynamic levels
    real(8), pointer,         intent(inout) :: P_M(:,:) ! output pointer to computed column pressure values on momentum levels

    ! Locals
    integer :: Vcode

    call msg('czp_calcReturnPressure_col_nl (czp)', 'START', verb_opt=2)

    Vcode = col_getVco(column)%vcode
    if (Vcode == 5005 .or. Vcode == 5002) then
      if ( .not. col_varExist(column,'P0')  ) then
        call utl_abort('czp_calcReturnPressure_col_nl (czp): for vcode 500x, variable P0 must be allocated in column')
      end if
      call calcPressure_col_nl_vcode500x(column, P_T, P_M)
    else if (Vcode == 21001) then
      if ( .not. (col_varExist(column,'P0') .and. col_varExist(column,'TT') .and. &
                  col_varExist(column,'HU'))  ) then
        call utl_abort('czp_calcReturnPressure_col_nl (czp): for vcode 2100x, variables P0, TT and HU must be allocated in column')
      end if
      call calcPressure_col_nl_vcode2100x(column, P_T, P_M)
    end if

    call msg('czp_calcReturnPressure_col_nl (czp)', 'END', verb_opt=2)
  end subroutine czp_calcReturnPressure_col_nl

  !---------------------------------------------------------
  ! calcPressure_col_nl_vcode2100x
  !---------------------------------------------------------
  subroutine calcPressure_col_nl_vcode2100x(column, P_T, P_M)
    !
    ! :Purpose: Compute pressure and return values in pointer arguments.
    !           GEM-H column input.
    !
    implicit none
    ! Development notes (@mad001)
    !   if feasible, consider reusing the same code for both 
    !   `calcPressure_{gsv,col}_nl_vcode2100x`
    !   (@mab001) Also should remove need for `Z_T/M` to be allocated in column
    !   and use local array instead.

    ! Arguments
    type(struct_columnData),  intent(in)    :: column ! reference column
    real(8), pointer,         intent(inout) :: P_T(:,:) ! output pointer to computed column pressure values on thermodynamic levels
    real(8), pointer,         intent(inout) :: P_M(:,:) ! output pointer to computed column pressure values on momentum levels

    ! Locals
    real(8), allocatable  :: tv(:)
    integer :: numCol, nLev_T, nLev_M
    integer :: colIndex, lev_T, lev_M
    real(8) :: lat, sLat, cLat
    real(8) :: P0, rMT, hu, tt, tv0, cmp, dh, Rgh
    real(8) :: scaleFactorBottom, logP
    real(8) :: Z_T, Z_M, Z_M1

    call msg('calcPressure_col_nl_vcode2100x (czp)', 'START', verb_opt=4)

    numCol = col_getNumCol(column)
    nLev_M = col_getNumLev(column, 'MM')
    nLev_T = col_getNumLev(column, 'TH')

    allocate(tv(nLev_T))

    do_onAllcolumns: do colIndex = lbound(column%lat,1), ubound(column%lat,1)
      ! column%lat populated in innovation_mod from obsSpaceData latitudes
      lat = column%lat(colIndex)
      sLat = sin(lat)
      cLat = cos(lat)

      ! surface values
      P0  = col_getElem(  column, 1, colIndex, 'P0') ! surface pressure
      rMT = col_getHeight(column, 1, colIndex, 'SF') ! surface height

      ! compute pressure on diagnostic (nLev_{T,M}) levels
      hu = col_getElem(column, nLev_T, colIndex, 'HU')
      tt = col_getElem(column, nLev_T, colIndex, 'TT')
      tv0 = phf_fotvt8(tt,hu)

      ! thermo diagnostic level
      Z_T = col_getHeight(column, nLev_T, colIndex, 'TH')
      cmp = gpscompressibility(P0,tt,hu)
      tv(nlev_T) = tv0*cmp
      dh = Z_T - rMT
      Rgh = phf_gravityalt(sLat, rMT+0.5D0*dh)
      P_T(colIndex, nlev_T) = P0*exp(-Rgh*dh/MPC_RGAS_DRY_AIR_R8/tv(nlev_T))

      ! momentum diagnostic level
      Z_M = col_getHeight(column,nLev_M,colIndex,'MM')
      dh = Z_M - rMT
      Rgh = phf_gravityalt(sLat, rMT+0.5D0*dh)
      P_M(colIndex, nlev_M) = P0*exp(-Rgh*dh/MPC_RGAS_DRY_AIR_R8/tv(nlev_T))

      call msg('calcPressure_col_nl_vcode2100x (czp)', &
           'Column index '//str(colIndex) //': lat='//str(lat)&
           //'   Surface: height='//str(rMT)//', P0='//str(P0) &
           //'   TH_diag_lvl ('//str(nLev_T)//'): height='//str(Z_T)&
                                  //', P_T='//str(P_T(colIndex, nlev_T)) &
           //'   MM_diag_lvl ('//str(nLev_M)//'): height='//str(Z_M)&
                                  //', P_M='//str(P_M(colIndex, nlev_M)), &
           verb_opt=6)

      ! compute pressure on all levels above except the last
      do lev_M = nlev_M-1, 1, -1
        lev_T = lev_M ! thermo level just below
        hu   = col_getElem(  column, lev_T,   colIndex, 'HU')
        tt   = col_getElem(  column, lev_T,   colIndex, 'TT')
        Z_M  = col_getHeight(column, lev_M,   colIndex, 'MM')
        Z_M1 = col_getHeight(column, lev_M+1, colIndex, 'MM')
        Z_T  = col_getHeight(column, lev_T,   colIndex, 'TH')

        tv0 = phf_fotvt8(tt,hu)
        dh = Z_M - Z_M1
        Rgh = phf_gravityalt(sLat, Z_M1+0.5D0*dh)

        ! approximation of tv from pressure on previous momentum level
        cmp = gpscompressibility(P_M(colIndex, lev_M+1),tt,hu)
        tv(lev_T) = tv0*cmp
        P_M(colIndex, lev_M) = P_M(colIndex, lev_M+1) * &
                            exp(-Rgh*dh/MPC_RGAS_DRY_AIR_R8/tv(lev_T))
        ! first interpolation of thermo pressure
        scaleFactorBottom = (Z_T-Z_M1)/(Z_M-Z_M1)
        logP = (1.0D0-scaleFactorBottom)*log(P_M(colIndex, lev_M+1)) + &
                              scaleFactorBottom*log(P_M(colIndex, lev_M))
        P_T(colIndex, lev_T) = exp(logP)

        ! second iteration on tv
        cmp = gpscompressibility(P_T(colIndex, lev_T),tt,hu)
        tv(lev_T) = tv0*cmp
        P_M(colIndex, lev_M) = P_M(colIndex, lev_M+1) * &
                            exp(-Rgh*dh/MPC_RGAS_DRY_AIR_R8/tv(lev_T))

        ! second iteration interpolation of thermo pressure
        logP = (1.0D0-scaleFactorBottom)*log(P_M(colIndex, lev_M+1)) + &
                              scaleFactorBottom*log(P_M(colIndex, lev_M))
        P_T(colIndex, lev_T) = exp(logP)

      end do
    end do do_onAllColumns

    deallocate(tv)

    call msg('calcPressure_col_nl_vcode2100x (czp)', 'END', verb_opt=4)
  end subroutine calcPressure_col_nl_vcode2100x

  !---------------------------------------------------------
  ! calcPressure_col_nl_vcode500x
  !---------------------------------------------------------
  subroutine calcPressure_col_nl_vcode500x(column, P_T, P_M)
    implicit none

    ! Arguments
    type(struct_columnData),  intent(in)    :: column ! reference column
    real(8), pointer,         intent(inout) :: P_T(:,:) ! output pointer to computed column pressure values on thermodynamic levels
    real(8), pointer,         intent(inout) :: P_M(:,:) ! output pointer to computed column pressure values on momentum levels

    ! Locals
    real(kind=8), allocatable :: Psfc(:,:)
    real(kind=8), pointer     :: zppobs(:,:,:)
    integer :: headerIndex, status

    call msg('calcPressure_col_nl_vcode500x (czp)', 'START', verb_opt=4)
    if ( col_getNumCol(column) <= 0 ) then
      call msg('calcPressure_col_nl_vcode500x (czp)',&
           'END (number of columns <= 0)', verb_opt=2)
      return
    end if

    if (.not.col_varExist(column,'P0')) then
      call utl_abort('calcPressure_col_nl (czp): P0 must be present as an analysis variable!')
    end if

    allocate(Psfc(1,col_getNumCol(column)))
    do headerIndex = 1,col_getNumCol(column)
      Psfc(1,headerIndex) = col_getElem(column,1,headerIndex,'P0')
    end do

    status=vgd_levels(column%vco%vgrid,ip1_list=column%vco%ip1_M,  &
                      levels=zppobs,sfc_field=Psfc,in_log=.false.)
    if(status.ne.VGD_OK) call utl_abort('ERROR with vgd_levels')
    P_M(:,:) = zppobs(1,:,:)
    if (associated(zppobs))  deallocate(zppobs)

    status=vgd_levels(column%vco%vgrid,ip1_list=column%vco%ip1_T,  &
                      levels=zppobs,sfc_field=Psfc,in_log=.false.)
    if(status.ne.VGD_OK) call utl_abort('ERROR with vgd_levels')
    P_T(:,:) = zppobs(1,:,:)
    if (associated(zppobs)) deallocate(zppobs)

    deallocate(Psfc)

    call msg('calcPressure_col_nl_vcode500x (czp)', 'END', verb_opt=4)
  end subroutine calcPressure_col_nl_vcode500x


  !---------------------------------------------------------
  ! calcPressure_col_tl
  !---------------------------------------------------------
  subroutine calcPressure_col_tl( columnInc, columnIncRef)
    !
    !:Purpose: Tangent pressure computation on the column.
    !
    implicit none

    ! Arguments
    type(struct_columnData), intent(inout) :: columnInc    ! column that will contain the P_T/P_M increments
    type(struct_columnData), intent(in)    :: columnIncRef ! column containing needed reference fields

    ! Locals
    integer :: Vcode

    call msg('calcPressure_col_tl (czp)', 'START', verb_opt=2)

    Vcode = col_getVco(columnIncRef)%vcode
    if (Vcode == 5005 .or. Vcode == 5002) then
      if ( .not. col_varExist(columnInc,'P_*')  ) then
        call utl_abort('calcPressure_col_tl (czp): for vcode 500x, variables P_M and P_T must be allocated in column')
      end if
      if ( .not. col_varExist(columnInc,'P0')  ) then
        call utl_abort('calcPressure_col_tl (czp): for vcode 500x, variable P0 must be allocated in column')
      end if
      call calcPressure_col_tl_vcode500x
    else if (Vcode == 21001) then
      ! Development notes (@mad001)
      !   probably some some gsv_varExist(statevector,.) needed for GEM-H
      call calcPressure_col_tl_vcode2100x
    end if

    call msg('calcPressure_col_tl (czp)', 'END', verb_opt=2)

    contains
      !---------------------------------------------------------
      ! calcPressure_col_tl_vcode2100x
      !---------------------------------------------------------
      subroutine calcPressure_col_tl_vcode2100x
        implicit none

        call utl_abort('calcPressure_col_tl (czp): vcode 21001 not implemented yet')

      end subroutine calcPressure_col_tl_vcode2100x

      !---------------------------------------------------------
      ! calcPressure_col_tl_vcode500x
      !---------------------------------------------------------
      subroutine calcPressure_col_tl_vcode500x
        implicit none

        ! Locals
        real(8)          :: Psfc
        real(8), pointer :: delPsfc(:,:), PsfcRef(:,:)
        real(8), pointer :: delP_T(:,:), delP_M(:,:)
        real(8), pointer :: dP_dPsfc_T(:), dP_dPsfc_M(:)
        integer          :: status, colIndex
        integer          :: lev_M, lev_T, nlev_T, nlev_M, numColumns

        call msg('calcPressure_col_tl_vcode500x (czp)', 'START', verb_opt=4)

        nullify(dP_dPsfc_T)
        nullify(dP_dPsfc_M)
        nullify(delPsfc)
        nullify(delP_T)
        nullify(delP_M)

        delP_M  => col_getAllColumns(columnInc,'P_M')
        delP_T  => col_getAllColumns(columnInc,'P_T')
        delPsfc => col_getAllColumns(columnInc,'P0')
        PsfcRef => col_getAllColumns(columnIncRef,'P0')

        nlev_T = col_getNumLev(columnInc,'TH')
        nlev_M = col_getNumLev(columnInc,'MM')
        numColumns = col_getNumCol(columnInc)

        do colIndex = 1, numColumns

          Psfc = PsfcRef(1,colIndex)

          ! dP_dPsfc_M
          nullify(dP_dPsfc_M)
          status = vgd_dpidpis(columnInc%vco%vgrid, &
                               columnInc%vco%ip1_M, &
                               dP_dPsfc_M, &
                               Psfc)
          if( status .ne. VGD_OK ) then
              call utl_abort('calcPressure_col_tl (czp): ERROR with vgd_dpidpis')
          end if
          ! calculate delP_M
          do lev_M = 1, nlev_M
            delP_M(lev_M,colIndex) = dP_dPsfc_M(lev_M) * delPsfc(1,colIndex)
          end do
          deallocate(dP_dPsfc_M)

          ! dP_dPsfc_T
          nullify(dP_dPsfc_T)
          status = vgd_dpidpis(columnInc%vco%vgrid, &
                               columnInc%vco%ip1_T, &
                               dP_dPsfc_T, &
                               Psfc)
          if( status .ne. VGD_OK ) then
              call utl_abort('calcPressure_col_tl (czp): ERROR with vgd_dpidpis')
          end if
          ! calculate delP_T
          do lev_T = 1, nlev_T
            delP_T(lev_T,colIndex) = dP_dPsfc_T(lev_T) * delPsfc(1,colIndex)
          end do
          deallocate(dP_dPsfc_T)

        end do

        call msg('calcPressure_col_tl_vcode500x (czp)', 'END', verb_opt=4)
      end subroutine calcPressure_col_tl_vcode500x

  end subroutine calcPressure_col_tl

  !---------------------------------------------------------
  ! calcPressure_col_ad
  !---------------------------------------------------------
  subroutine calcPressure_col_ad( columnInc, columnIncRef)
    !
    !:Purpose: Adjoint of pressure computation on the column.
    !
    implicit none

    ! Arguments
    type(struct_columnData), intent(inout) :: columnInc    ! column that will contain increments of P_M/P_T
    type(struct_columnData), intent(in)    :: columnIncRef ! column containing needed reference fields

    ! Locals
    integer :: Vcode

    call msg('calcPressure_col_ad (czp)', 'START', verb_opt=2)

    Vcode = col_getVco(columnIncRef)%vcode
    if (Vcode == 5005 .or. Vcode == 5002) then
      if ( .not. col_varExist(columnInc,'P_*')  ) then
        call utl_abort('calcPressure_col_ad (czp): for vcode 500x, variables P_M and P_T must be allocated in column')
      end if
      if ( .not. col_varExist(columnInc,'P0')  ) then
        call utl_abort('calcPressure_col_ad (czp): for vcode 500x, variable P0 must be allocated in column')
      end if
      call calcPressure_col_ad_vcode500x
    else if (Vcode == 21001) then
      ! Development notes (@mad001)
      !   probably some some gsv_varExist(statevector,.) needed for GEM-H
      call calcPressure_col_ad_vcode2100x
    end if

    call msg('calcPressure_col_ad (czp)', 'END', verb_opt=2)

    contains
      !---------------------------------------------------------
      ! calcPressure_col_ad_vcode2100x
      !---------------------------------------------------------
      subroutine calcPressure_col_ad_vcode2100x
        implicit none

        call utl_abort('calcPressure_col_ad (czp): vcode 21001 not implemented yet')

      end subroutine calcPressure_col_ad_vcode2100x

      !---------------------------------------------------------
      ! calcPressure_col_ad_vcode500x
      !---------------------------------------------------------
      subroutine calcPressure_col_ad_vcode500x
        implicit none

        ! Locals
        real(8)          :: Psfc
        real(8), pointer :: delPsfc(:,:), PsfcRef(:,:)
        real(8), pointer :: delP_T(:,:), delP_M(:,:)
        real(8), pointer :: dP_dPsfc_T(:), dP_dPsfc_M(:)
        integer          :: status, colIndex
        integer          :: lev_M, lev_T, nlev_T, nlev_M, numColumns

        call msg('calcPressure_col_ad_vcode500x (czp)', 'START', verb_opt=4)

        nullify(delPsfc)
        nullify(PsfcRef)
        nullify(delP_T)
        nullify(delP_M)
        nullify(dP_dPsfc_T)
        nullify(dP_dPsfc_M)

        delP_M  => col_getAllColumns(columnInc,'P_M')
        delP_T  => col_getAllColumns(columnInc,'P_T')
        delPsfc => col_getAllColumns(columnInc,'P0')
        PsfcRef => col_getAllColumns(columnIncRef,'P0')

        nlev_T = col_getNumLev(columnInc,'TH')
        nlev_M = col_getNumLev(columnInc,'MM')
        numColumns = col_getNumCol(columnInc)

        do colIndex = 1, numColumns

          Psfc = PsfcRef(1,colIndex)

          ! dP_dPsfc_M
          nullify(dP_dPsfc_M)
          status = vgd_dpidpis(columnInc%vco%vgrid, &
                               columnInc%vco%ip1_M, &
                               dP_dPsfc_M, &
                               Psfc)
          if( status .ne. VGD_OK ) then
              call utl_abort('calcPressure_col_ad (czp): ERROR with vgd_dpidpis')
          end if
          ! calculate delP_M
          do lev_M = 1, nlev_M
            delPsfc(1,colIndex) = delPsfc(1,colIndex) + &
                 dP_dPsfc_M(lev_M) * delP_M(lev_M,colIndex)
          end do
          deallocate(dP_dPsfc_M)

          ! dP_dPsfc_T
          nullify(dP_dPsfc_T)
          status = vgd_dpidpis(columnInc%vco%vgrid, &
                               columnInc%vco%ip1_T, &
                               dP_dPsfc_T, &
                               Psfc)
          if( status .ne. VGD_OK ) then
            call utl_abort('calcPressure_col_ad (czp): ERROR with vgd_dpidpis')
          end if
          ! calculate delP_T
          do lev_T = 1, nlev_T
            delPsfc(1,colIndex) = &
                delPsfc(1,colIndex) + dP_dPsfc_T(lev_T) * delP_T(lev_T,colIndex)
          end do
          deallocate(dP_dPsfc_T)

        end do
        call msg('calcPressure_col_ad_vcode500x (czp)', 'END', verb_opt=4)
      end subroutine calcPressure_col_ad_vcode500x

  end subroutine calcPressure_col_ad

  !---------------------------------------------------------------------
  ! subroutines wrapping vgd_levels for 3D fields queries
  !---------------------------------------------------------------------
  
  !---------------------------------------------------------
  ! czp_fetch3DField
  !---------------------------------------------------------

  !---------------------------------------------------------------------
  ! subroutines wrapping vgd_levels for profile queries
  !---------------------------------------------------------------------

  !---------------------------------------------------------
  ! czp_fetchProfile
  !---------------------------------------------------------
  subroutine czp_fetchProfile(vco, sfcValue, profM_opt, profT_opt)
    !
    ! :Purpose: Return pressure profile for both momentum and thermodynamic levels
    !
    implicit none

    ! Arguments
    type(struct_vco),   intent(in)            :: vco          ! Vertical descriptor
    real(8),           intent(in)             :: sfcValue     ! Surface field reference for coordinate
    real(8), optional, intent(inout), pointer :: profM_opt(:) ! Momemtum levels profile
    real(8), optional, intent(inout), pointer :: profT_opt(:) ! Thermodynamic levels profile

    ! Locals
    integer :: status

    if (present(profM_opt)) then
      nullify(profM_opt)
      status = vgd_levels(vco%vgrid, ip1_list=vco%ip1_M, &
                          levels=profM_opt, sfc_field=sfcValue, in_log=.false.)
      if ( status .ne. VGD_OK ) then
        call utl_abort('calcPressure_prof (czp):  ERROR with vgd_levels (momentum levels)')
      end if
    end if

    if (present(profT_opt)) then
      nullify(profT_opt)
      status = vgd_levels(vco%vgrid, ip1_list=vco%ip1_T, &
                          levels=profT_opt, sfc_field=sfcValue, in_log=.false.)
      if ( status .ne. VGD_OK ) then
        call utl_abort('calcPressure_prof (czp):  ERROR with vgd_levels (thermodynamic levels)')
      end if
    end if
  end subroutine czp_fetchProfile


  !---------------------------------------------------------------------
  ! helper private functions and subroutines
  !---------------------------------------------------------------------

  !---------------------------------------------------------
  ! calcHeightCoeff_gsv
  !---------------------------------------------------------
  subroutine calcHeightCoeff_gsv(statevector)
    !
    ! :Purpose: Calculating the coefficients of height for
    !           czp_calcHeight_tl/czp_calcHeight_ad
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(in) :: statevector

    ! Locals
    integer :: lev_T,nlev_M,nlev_T,numStep,stepIndex,latIndex,lonIndex,Vcode
    real(8) :: hu,tt,Pr,height_T,cmp,cmp_TT,cmp_HU,cmp_P0_1,cmp_P0_2,ratioP1
    real(4) :: lat_4
    real(8) :: Rgh, sLat, lat_8
    real(8), pointer :: hu_ptr(:,:,:,:),tt_ptr(:,:,:,:)
    real(8), pointer :: P_T_ptr(:,:,:,:),P_M_ptr(:,:,:,:)
    real(8), pointer :: height_T_ptr(:,:,:,:)
    type(struct_vco), pointer :: vco

    logical, save :: firstTimeHeightCoeff_gsv = .true.

    if ( .not. firstTimeHeightCoeff_gsv ) return

    call msg('calcHeightCoeff_gsv (czp)', 'START', verb_opt=2)

    ! initialize and save coefficients for increased efficiency
    ! (assumes no relinearization)
    firstTimeHeightCoeff_gsv = .false.

    vco => gsv_getVco(statevector)
    Vcode = vco%vcode

    nlev_T = gsv_getNumLev(statevector,'TH')
    nlev_M = gsv_getNumLev(statevector,'MM')
    numStep = statevector%numstep

    ! saved arrays
    allocate(coeff_M_TT_gsv(&
        statevector%myLonBeg:statevector%myLonEnd, &
        statevector%myLatBeg:statevector%myLatEnd, &
        nlev_T,numStep))
    allocate(coeff_M_HU_gsv(&
        statevector%myLonBeg:statevector%myLonEnd, &
        statevector%myLatBeg:statevector%myLatEnd, &
        nlev_T,numStep))
    allocate(coeff_M_P0_delPM_gsv(&
        statevector%myLonBeg:statevector%myLonEnd, &
        statevector%myLatBeg:statevector%myLatEnd, &
        nlev_T,numStep))
    allocate(coeff_M_P0_dP_delPT_gsv(&
        statevector%myLonBeg:statevector%myLonEnd, &
        statevector%myLatBeg:statevector%myLatEnd, &
        nlev_T,numStep))
    allocate(coeff_M_P0_dP_delP0_gsv(&
        statevector%myLonBeg:statevector%myLonEnd, &
        statevector%myLatBeg:statevector%myLatEnd, &
        nlev_T,numStep))

    allocate(coeff_T_TT_gsv(&
        statevector%myLonBeg:statevector%myLonEnd, &
        statevector%myLatBeg:statevector%myLatEnd, &
        numStep))
    allocate(coeff_T_HU_gsv(&
        statevector%myLonBeg:statevector%myLonEnd, &
        statevector%myLatBeg:statevector%myLatEnd, &
        numStep))
    allocate(coeff_T_P0_delP1_gsv(&
        statevector%myLonBeg:statevector%myLonEnd, &
        statevector%myLatBeg:statevector%myLatEnd, &
        numStep))
    allocate(coeff_T_P0_dP_delPT_gsv(&
        statevector%myLonBeg:statevector%myLonEnd, &
        statevector%myLatBeg:statevector%myLatEnd, &
        numStep))
    allocate(coeff_T_P0_dP_delP0_gsv(&
        statevector%myLonBeg:statevector%myLonEnd, &
        statevector%myLatBeg:statevector%myLatEnd, &
        numStep))

    coeff_M_TT_gsv(:,:,:,:) = 0.0D0
    coeff_M_HU_gsv(:,:,:,:) = 0.0D0

    coeff_M_P0_delPM_gsv(:,:,:,:) = 0.0D0

    coeff_M_P0_dP_delPT_gsv(:,:,:,:) = 0.0D0
    coeff_M_P0_dP_delP0_gsv(:,:,:,:) = 0.0D0

    coeff_T_TT_gsv(:,:,:) = 0.0D0
    coeff_T_HU_gsv(:,:,:) = 0.0D0

    coeff_T_P0_delP1_gsv(:,:,:) = 0.0D0

    coeff_T_P0_dP_delPT_gsv(:,:,:) = 0.0D0
    coeff_T_P0_dP_delP0_gsv(:,:,:) = 0.0D0

    call gsv_getField(statevector,hu_ptr,'HU')
    call gsv_getField(statevector,tt_ptr,'TT')
    call gsv_getField(statevector,P_T_ptr,'P_T')
    call gsv_getField(statevector,P_M_ptr,'P_M')
    call gsv_getField(statevector,height_T_ptr,'Z_T')

    if_calcHeightCoeff_gsv_vcodes : if (Vcode == 5002) then

      do stepIndex = 1, numStep
        do lev_T = 1, (nlev_T-1)
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd
              if ( lev_T == 1 ) then
                ! compute height coefficients on only the top thermo level
                ratioP1 = log( P_M_ptr(lonIndex,latIndex,1,stepIndex) / &
                               P_T_ptr(lonIndex,latIndex,1,stepIndex) )
                hu = max(hu_ptr(lonIndex,latIndex,1,stepIndex),&
                                MPC_MINIMUM_HU_R8)
                tt = tt_ptr(lonIndex,latIndex,1,stepIndex)
                Pr = P_T_ptr(lonIndex,latIndex,1,stepIndex)
                height_T = height_T_ptr(lonIndex,latIndex,1,stepIndex)

                cmp = gpscompressibility(Pr,tt,hu)
                cmp_TT = gpscompressibility_TT(Pr,tt,hu)
                cmp_HU = gpscompressibility_HU(Pr,tt,hu)
                cmp_P0_1 = gpscompressibility_P0_1(Pr,tt,hu,1.0d0)
                cmp_P0_2 = gpscompressibility_P0_2(Pr,tt,hu)

                ! Gravity acceleration
                lat_4 = statevector%hco%lat2d_4(lonIndex,latIndex)
                lat_8 = real(lat_4,8)
                sLat = sin(lat_8)
                Rgh = phf_gravityalt(sLat, height_T)

                coeff_T_TT_gsv(lonIndex,latIndex,stepIndex) = &
                    (MPC_RGAS_DRY_AIR_R8 / Rgh) * (phf_fottva(hu,1.0D0) * cmp + &
                    phf_fotvt8(tt,hu) * cmp_TT) * ratioP1

                coeff_T_HU_gsv(lonIndex,latIndex,stepIndex) = &
                    (MPC_RGAS_DRY_AIR_R8 / Rgh) * (phf_folnqva(hu,tt,1.0d0) / &
                    hu * cmp + phf_fotvt8(tt,hu) * cmp_HU) * ratioP1

                coeff_T_P0_delP1_gsv(lonIndex,latIndex,stepIndex) = &
                    (MPC_RGAS_DRY_AIR_R8 / Rgh) * phf_fotvt8(tt,hu) * cmp

                coeff_T_P0_dP_delPT_gsv(lonIndex,latIndex,stepIndex) = &
                    (MPC_RGAS_DRY_AIR_R8 / Rgh) * phf_fotvt8(tt,hu) * cmp_P0_1 * &
                    ratioP1

                coeff_T_P0_dP_delP0_gsv(lonIndex,latIndex,stepIndex) = &
                    (MPC_RGAS_DRY_AIR_R8 / Rgh) * phf_fotvt8(tt,hu) * cmp_P0_2 * &
                    ratioP1
              else
                ! compute height coefficients on momentum levels
                ratioP1 = log( P_M_ptr(lonIndex,latIndex,lev_T  ,stepIndex) / &
                               P_M_ptr(lonIndex,latIndex,lev_T-1,stepIndex) )
                hu = max( hu_ptr(lonIndex,latIndex,lev_T,stepIndex),&
                          MPC_MINIMUM_HU_R8)
                tt = tt_ptr(lonIndex,latIndex,lev_T,stepIndex)
                Pr = P_T_ptr(lonIndex,latIndex,lev_T,stepIndex)
                height_T = height_T_ptr(lonIndex,latIndex,lev_T,stepIndex)

                cmp = gpscompressibility(Pr,tt,hu)
                cmp_TT = gpscompressibility_TT(Pr,tt,hu)
                cmp_HU = gpscompressibility_HU(Pr,tt,hu)
                cmp_P0_1 = gpscompressibility_P0_1(Pr,tt,hu,1.0d0)
                cmp_P0_2 = gpscompressibility_P0_2(Pr,tt,hu)

                ! Gravity acceleration
                lat_4 = statevector%hco%lat2d_4(lonIndex,latIndex)
                lat_8 = real(lat_4,8)
                sLat = sin(lat_8)
                Rgh = phf_gravityalt(sLat, height_T)

                coeff_M_TT_gsv(lonIndex,latIndex,lev_T,stepIndex) = &
                    (MPC_RGAS_DRY_AIR_R8 / Rgh) * (phf_fottva(hu,1.0D0) * cmp + &
                    phf_fotvt8(tt,hu) * cmp_TT) * ratioP1

                coeff_M_HU_gsv(lonIndex,latIndex,lev_T,stepIndex) = &
                    (MPC_RGAS_DRY_AIR_R8 / RgH) * (phf_folnqva(hu,tt,1.0d0) / hu * &
                    cmp + phf_fotvt8(tt,hu) * cmp_HU) * ratioP1

                coeff_M_P0_delPM_gsv   (lonIndex,latIndex,lev_T,stepIndex) = &
                    (MPC_RGAS_DRY_AIR_R8 / Rgh) * phf_fotvt8(tt,hu) * cmp

                coeff_M_P0_dP_delPT_gsv(lonIndex,latIndex,lev_T,stepIndex) = &
                    (MPC_RGAS_DRY_AIR_R8 / Rgh) * phf_fotvt8(tt,hu) * cmp_P0_1 * &
                    ratioP1

                coeff_M_P0_dP_delP0_gsv(lonIndex,latIndex,lev_T,stepIndex) = &
                    (MPC_RGAS_DRY_AIR_R8 / Rgh) * phf_fotvt8(tt,hu) * cmp_P0_2 * &
                    ratioP1
              end if
            end do
          end do
        end do
      end do

    else if (Vcode == 5005) then if_calcHeightCoeff_gsv_vcodes

      do stepIndex = 1, numStep
        do lev_T = 1, (nlev_T-1)
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd
              ! compute height coefficients on momentum levels
              ratioP1 = log( P_M_ptr(lonIndex,latIndex,lev_T+1,stepIndex) / &
                             P_M_ptr(lonIndex,latIndex,lev_T  ,stepIndex) )
              hu = max( hu_ptr(lonIndex,latIndex,lev_T,stepIndex),&
                        MPC_MINIMUM_HU_R8)
              tt = tt_ptr(lonIndex,latIndex,lev_T,stepIndex)
              Pr = P_T_ptr(lonIndex,latIndex,lev_T,stepIndex)
              height_T = height_T_ptr(lonIndex,latIndex,lev_T,stepIndex)

              cmp = gpscompressibility(Pr,tt,hu)
              cmp_TT = gpscompressibility_TT(Pr,tt,hu)
              cmp_HU = gpscompressibility_HU(Pr,tt,hu)
              cmp_P0_1 = gpscompressibility_P0_1(Pr,tt,hu,1.0d0)
              cmp_P0_2 = gpscompressibility_P0_2(Pr,tt,hu)

              ! Gravity acceleration
              lat_4 = statevector%hco%lat2d_4(lonIndex,latIndex)
              lat_8 = real(lat_4,8)
              sLat = sin(lat_8)
              Rgh = phf_gravityalt(sLat, height_T)

              coeff_M_TT_gsv(lonIndex,latIndex,lev_T,stepIndex) = &
                  (MPC_RGAS_DRY_AIR_R8 / Rgh) * (phf_fottva(hu,1.0D0) * cmp + &
                  phf_fotvt8(tt,hu) * cmp_TT) * ratioP1

              coeff_M_HU_gsv(lonIndex,latIndex,lev_T,stepIndex) = &
                  (MPC_RGAS_DRY_AIR_R8 / RgH) * (phf_folnqva(hu,tt,1.0d0) / hu * &
                  cmp + phf_fotvt8(tt,hu) * cmp_HU) * ratioP1

              coeff_M_P0_delPM_gsv(lonIndex,latIndex,lev_T,stepIndex) = &
                  (MPC_RGAS_DRY_AIR_R8 / Rgh) * phf_fotvt8(tt,hu) * cmp

              coeff_M_P0_dP_delPT_gsv(lonIndex,latIndex,lev_T,stepIndex) = &
                  (MPC_RGAS_DRY_AIR_R8 / Rgh) * phf_fotvt8(tt,hu) * cmp_P0_1 * &
                  ratioP1

              coeff_M_P0_dP_delP0_gsv(lonIndex,latIndex,lev_T,stepIndex) = &
                  (MPC_RGAS_DRY_AIR_R8 / Rgh) * phf_fotvt8(tt,hu) * cmp_P0_2 * &
                  ratioP1
            end do
          end do
        end do
      end do

    else if_calcHeightCoeff_gsv_vcodes
      call utl_abort('calcHeightCoeff_gsv (czp): only vcode 5002 and 5005 implemented')

    end if if_calcHeightCoeff_gsv_vcodes

    call msg('calcHeightCoeff_gsv (czp)', 'END', verb_opt=2)
  end subroutine calcHeightCoeff_gsv

  !---------------------------------------------------------
  ! calcHeightCoeff_col
  !---------------------------------------------------------
  subroutine calcHeightCoeff_col(column)
    !
    ! :Purpose: Calculating the coefficients of height for
    !           czp_calcHeight_tl/czp_calcHeight_ad
    !
    implicit none

    ! Arguments
    type(struct_columnData), intent(in) :: column

    ! Locals
    integer :: lev_T,nlev_M,nlev_T,numColumns,colIndex,Vcode
    real(8) :: hu,tt,Pr,height_T,cmp,cmp_TT,cmp_HU,cmp_P0_1,cmp_P0_2,ratioP1
    real(8) :: Rgh, sLat, lat_8
    real(8), pointer :: hu_ptr(:,:),tt_ptr(:,:)
    real(8), pointer :: P_T_ptr(:,:),P_M_ptr(:,:)
    real(8), pointer :: height_T_ptr(:,:)
    type(struct_vco), pointer :: vco

    logical, save :: firstTimeHeightCoeff_col = .true.

    if ( .not. firstTimeHeightCoeff_col ) return

    call msg('calcHeightCoeff_col (czp)', 'START', verb_opt=2)

    ! initialize and save coefficients for increased efficiency
    ! (assumes no relinearization)
    firstTimeHeightCoeff_col = .false.

    vco => col_getVco(column)
    Vcode = vco%vcode

    nlev_T = col_getNumLev(column,'TH')
    nlev_M = col_getNumLev(column,'MM')
    numColumns = col_getNumCol(column)

    ! saved arrays
    allocate(coeff_M_TT_col         (nlev_T,numColumns))
    allocate(coeff_M_HU_col         (nlev_T,numColumns))
    allocate(coeff_M_P0_delPM_col   (nlev_T,numColumns))
    allocate(coeff_M_P0_dP_delPT_col(nlev_T,numColumns))
    allocate(coeff_M_P0_dP_delP0_col(nlev_T,numColumns))

    allocate(coeff_T_TT_col         (numColumns))
    allocate(coeff_T_HU_col         (numColumns))
    allocate(coeff_T_P0_delP1_col   (numColumns))
    allocate(coeff_T_P0_dP_delPT_col(numColumns))
    allocate(coeff_T_P0_dP_delP0_col(numColumns))

    coeff_M_TT_col(:,:) = 0.0D0
    coeff_M_HU_col(:,:) = 0.0D0

    coeff_M_P0_delPM_col(:,:) = 0.0D0

    coeff_M_P0_dP_delPT_col(:,:) = 0.0D0
    coeff_M_P0_dP_delP0_col(:,:) = 0.0D0

    coeff_T_TT_col(:) = 0.0D0
    coeff_T_HU_col(:) = 0.0D0

    coeff_T_P0_delP1_col(:) = 0.0D0

    coeff_T_P0_dP_delPT_col(:) = 0.0D0
    coeff_T_P0_dP_delP0_col(:) = 0.0D0

    hu_ptr       => col_getAllColumns(column,'HU')
    tt_ptr       => col_getAllColumns(column,'TT')
    P_T_ptr      => col_getAllColumns(column,'P_T')
    P_M_ptr      => col_getAllColumns(column,'P_M')
    height_T_ptr => col_getAllColumns(column,'Z_T')

    if_calcHeightCoeff_col_vcodes : if (Vcode == 5002) then

      do colIndex = 1, numColumns
        do lev_T = 1, (nlev_T-1)
          if ( lev_T == 1 ) then
            ! compute height coefficients on only the top thermo level
            ratioP1 = log( P_M_ptr(1,colIndex) / &
                           P_T_ptr(1,colIndex) )
            hu = max(hu_ptr(1,colIndex),MPC_MINIMUM_HU_R8)
            tt = tt_ptr(1,colIndex)
            Pr = P_T_ptr(1,colIndex)
            height_T = height_T_ptr(1,colIndex)

            cmp = gpscompressibility(Pr,tt,hu)
            cmp_TT = gpscompressibility_TT(Pr,tt,hu)
            cmp_HU = gpscompressibility_HU(Pr,tt,hu)
            cmp_P0_1 = gpscompressibility_P0_1(Pr,tt,hu,1.0d0)
            cmp_P0_2 = gpscompressibility_P0_2(Pr,tt,hu)

            ! Gravity acceleration
            lat_8 = column%lat(colIndex)
            sLat = sin(lat_8)
            Rgh = phf_gravityalt(sLat, height_T)

            coeff_T_TT_col(colIndex) = &
                (MPC_RGAS_DRY_AIR_R8 / Rgh) * (phf_fottva(hu,1.0D0) * cmp + &
                phf_fotvt8(tt,hu) * cmp_TT) * ratioP1

            coeff_T_HU_col(colIndex) = &
                (MPC_RGAS_DRY_AIR_R8 / Rgh) * (phf_folnqva(hu,tt,1.0d0) / hu * &
                cmp + phf_fotvt8(tt,hu) * cmp_HU) * ratioP1

            coeff_T_P0_delP1_col   (colIndex) = &
                (MPC_RGAS_DRY_AIR_R8 / Rgh) * phf_fotvt8(tt,hu) * cmp

            coeff_T_P0_dP_delPT_col(colIndex) = &
                (MPC_RGAS_DRY_AIR_R8 / Rgh) * phf_fotvt8(tt,hu) * cmp_P0_1 * ratioP1

            coeff_T_P0_dP_delP0_col(colIndex) = &
                (MPC_RGAS_DRY_AIR_R8 / Rgh) * phf_fotvt8(tt,hu) * cmp_P0_2 * ratioP1
          else
            ! compute height coefficients on momentum levels
            ratioP1 = log( P_M_ptr(lev_T  ,colIndex) / &
                           P_M_ptr(lev_T-1,colIndex) )
            hu = max(hu_ptr(lev_T,colIndex),MPC_MINIMUM_HU_R8)
            tt = tt_ptr(lev_T,colIndex)
            Pr = P_T_ptr(lev_T,colIndex)
            height_T = height_T_ptr(lev_T,colIndex)

            cmp = gpscompressibility(Pr,tt,hu)
            cmp_TT = gpscompressibility_TT(Pr,tt,hu)
            cmp_HU = gpscompressibility_HU(Pr,tt,hu)
            cmp_P0_1 = gpscompressibility_P0_1(Pr,tt,hu,1.0d0)
            cmp_P0_2 = gpscompressibility_P0_2(Pr,tt,hu)

            ! Gravity acceleration
            lat_8 = column%lat(colIndex)
            sLat = sin(lat_8)
            Rgh = phf_gravityalt(sLat, height_T)

            coeff_M_TT_col(lev_T,colIndex) = &
                (MPC_RGAS_DRY_AIR_R8 / Rgh) * (phf_fottva(hu,1.0D0) * cmp + &
                phf_fotvt8(tt,hu) * cmp_TT) * ratioP1

            coeff_M_HU_col(lev_T,colIndex) = &
                (MPC_RGAS_DRY_AIR_R8 / RgH) * (phf_folnqva(hu,tt,1.0d0) / hu * &
                cmp + phf_fotvt8(tt,hu) * cmp_HU) * ratioP1

            coeff_M_P0_delPM_col(lev_T,colIndex) = &
                (MPC_RGAS_DRY_AIR_R8 / Rgh) * phf_fotvt8(tt,hu) * cmp

            coeff_M_P0_dP_delPT_col(lev_T,colIndex) = &
                (MPC_RGAS_DRY_AIR_R8 / Rgh) * phf_fotvt8(tt,hu) * cmp_P0_1 * ratioP1

            coeff_M_P0_dP_delP0_col(lev_T,colIndex) = &
                (MPC_RGAS_DRY_AIR_R8 / Rgh) * phf_fotvt8(tt,hu) * cmp_P0_2 * ratioP1
          end if
        end do
      end do

    else if (Vcode == 5005) then if_calcHeightCoeff_col_vcodes

      do colIndex = 1, numColumns
        do lev_T = 1, (nlev_T-1)
          ! compute height coefficients on momentum levels
          ratioP1 = log( P_M_ptr(lev_T+1,colIndex) / &
                         P_M_ptr(lev_T  ,colIndex) )
          hu = max(hu_ptr(lev_T,colIndex),MPC_MINIMUM_HU_R8)
          tt = tt_ptr(lev_T,colIndex)
          Pr = P_T_ptr(lev_T,colIndex)
          height_T = height_T_ptr(lev_T,colIndex)

          cmp = gpscompressibility(Pr,tt,hu)
          cmp_TT = gpscompressibility_TT(Pr,tt,hu)
          cmp_HU = gpscompressibility_HU(Pr,tt,hu)
          cmp_P0_1 = gpscompressibility_P0_1(Pr,tt,hu,1.0d0)
          cmp_P0_2 = gpscompressibility_P0_2(Pr,tt,hu)

          ! Gravity acceleration
          lat_8 = column%lat(colIndex)
          sLat = sin(lat_8)
          Rgh = phf_gravityalt(sLat, height_T)

          coeff_M_TT_col(lev_T,colIndex) = &
              (MPC_RGAS_DRY_AIR_R8 / Rgh) * (phf_fottva(hu,1.0D0) * cmp + &
              phf_fotvt8(tt,hu) * cmp_TT) * ratioP1

          coeff_M_HU_col(lev_T,colIndex) = &
              (MPC_RGAS_DRY_AIR_R8 / RgH) * (phf_folnqva(hu,tt,1.0d0) / hu * cmp + &
              phf_fotvt8(tt,hu) * cmp_HU) * ratioP1

          coeff_M_P0_delPM_col   (lev_T,colIndex) = &
              (MPC_RGAS_DRY_AIR_R8 / Rgh) * phf_fotvt8(tt,hu) * cmp

          coeff_M_P0_dP_delPT_col(lev_T,colIndex) = &
              (MPC_RGAS_DRY_AIR_R8 / Rgh) * phf_fotvt8(tt,hu) * cmp_P0_1 * ratioP1

          coeff_M_P0_dP_delP0_col(lev_T,colIndex) = &
              (MPC_RGAS_DRY_AIR_R8 / Rgh) * phf_fotvt8(tt,hu) * cmp_P0_2 * ratioP1
        end do
      end do

    end if if_calcHeightCoeff_col_vcodes

    call msg('calcHeightCoeff_col (czp)', 'END', verb_opt=2)
  end subroutine calcHeightCoeff_col

  !---------------------------------------------------------
  ! gpscompressibility
  !---------------------------------------------------------
  function gpscompressibility(p,t,q)
    implicit none

    ! Arguments
    real(8), intent(in)  :: p
    real(8), intent(in)  :: t
    real(8), intent(in)  :: q

    ! Locals
    real(8)              :: gpscompressibility

    real(8), parameter   :: a0= 1.58123D-6
    real(8), parameter   :: a1=-2.9331D-8
    real(8), parameter   :: a2= 1.1043D-10
    real(8), parameter   :: b0= 5.707D-6
    real(8), parameter   :: b1=-2.051D-8
    real(8), parameter   :: c0= 1.9898D-4
    real(8), parameter   :: c1=-2.376D-6
    real(8), parameter   :: d = 1.83D-11
    real(8), parameter   :: e =-0.765D-8

    real(8)         :: x,tc,pt,tc2,x2

    if ( t <= 0 ) call utl_abort('gpscompressibility: t <= 0')

    x  = gps_p_wa * q / (1.D0 + gps_p_wb * q)
    ! Estimate, from CIPM, Picard (2008)
    tc = t - MPC_K_C_DEGREE_OFFSET_R8
    pt = p / t
    tc2= tc * tc
    x2 = x * x
    gpscompressibility = 1.D0 - pt * &
        (a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2) + pt*pt*(d+e*x2)
  end function gpscompressibility

  !---------------------------------------------------------
  ! gpscompressibility_TT
  !---------------------------------------------------------
  function gpscompressibility_TT(p,t,q)
    implicit none

    ! Arguments
    real(8), intent(in)  :: p
    real(8), intent(in)  :: t
    real(8), intent(in)  :: q

    ! Locals
    real(8)              :: gpscompressibility_TT

    real(8), parameter   :: a0= 1.58123D-6
    real(8), parameter   :: a1=-2.9331D-8
    real(8), parameter   :: a2= 1.1043D-10
    real(8), parameter   :: b0= 5.707D-6
    real(8), parameter   :: b1=-2.051D-8
    real(8), parameter   :: c0= 1.9898D-4
    real(8), parameter   :: c1=-2.376D-6
    real(8), parameter   :: d = 1.83D-11
    real(8), parameter   :: e =-0.765D-8

    real(8)         :: x,tc,pt,tc2,x2
    real(8)         :: d_x,d_tc,d_pt,d_tc2,d_x2

    if ( t <= 0 ) call utl_abort('gpscompressibility_TT: t <= 0')

    x  = gps_p_wa * q / (1.D0 + gps_p_wb * q)
    d_x  = 0.0D0
    ! Estimate, from CIPM, Picard (2008)
    tc = t - MPC_K_C_DEGREE_OFFSET_R8
    d_tc = 1.D0
    pt = p / t
    d_pt = - p / t**2
    tc2= tc * tc
    d_tc2= 2 * tc * d_tc
    x2 = x * x
    d_x2 = 2 * x * d_x
    gpscompressibility_TT = &
        -d_pt * (a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2) - pt * &
        (a1*d_tc + a2*d_tc2 + b1*d_tc*x + (b0+b1*tc)*d_x + c1*d_tc*x2 + &
        (c0+c1*tc)*d_x2) + 2*pt*d_pt*(d+e*x2) + pt*pt*e*d_x2
  end function gpscompressibility_TT

  !---------------------------------------------------------
  ! gpscompressibility_HU
  !---------------------------------------------------------
  function gpscompressibility_HU(p,t,q)
    implicit none

    ! Arguments
    real(8), intent(in)  :: p
    real(8), intent(in)  :: t
    real(8), intent(in)  :: q

    ! Locals
    real(8)              :: gpscompressibility_HU

    real(8), parameter   :: a0= 1.58123D-6
    real(8), parameter   :: a1=-2.9331D-8
    real(8), parameter   :: a2= 1.1043D-10
    real(8), parameter   :: b0= 5.707D-6
    real(8), parameter   :: b1=-2.051D-8
    real(8), parameter   :: c0= 1.9898D-4
    real(8), parameter   :: c1=-2.376D-6
    real(8), parameter   :: d = 1.83D-11
    real(8), parameter   :: e =-0.765D-8

    real(8)         :: x,tc,pt,tc2,x2
    real(8)         :: d_x,d_tc,d_pt,d_tc2,d_x2

    if ( t <= 0 ) call utl_abort('gpscompressibility_HU: t <= 0')

    x  = gps_p_wa * q / (1.D0+gps_p_wb*q)
    d_x  = gps_p_wa * (1.0D0 / (1.D0+gps_p_wb*q) - q / (1.D0+gps_p_wb*q)**2 * gps_p_wb * 1.0D0)
    ! Estimate, from CIPM, Picard (2008)
    tc = t - MPC_K_C_DEGREE_OFFSET_R8
    d_tc = 0.0D0
    pt = p / t
    d_pt = 0.0D0
    tc2= tc * tc
    d_tc2= 2 * tc * d_tc
    x2 = x * x
    d_x2 = 2 * x * d_x
    gpscompressibility_HU = &
        -d_pt * (a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2) - &
        pt * (a1*d_tc + a2*d_tc2 + b1*d_tc*x + (b0+b1*tc)*d_x + c1*d_tc*x2 + &
        (c0+c1*tc)*d_x2) + 2*pt*d_pt*(d+e*x2) + pt*pt*e*d_x2
  end function gpscompressibility_HU

  !---------------------------------------------------------
  ! gpscompressibility_P0_1
  !---------------------------------------------------------
  function gpscompressibility_P0_1(p,t,q,dpdp0)
    ! gpscompressibility_P0_1 has dpdp0 dependency
    implicit none

    ! Arguments
    real(8), intent(in)  :: p
    real(8), intent(in)  :: t
    real(8), intent(in)  :: q
    real(8), intent(in)  :: dpdp0

    ! Locals
    real(8)              :: gpscompressibility_P0_1

    real(8), parameter   :: a0= 1.58123D-6
    real(8), parameter   :: a1=-2.9331D-8
    real(8), parameter   :: a2= 1.1043D-10
    real(8), parameter   :: b0= 5.707D-6
    real(8), parameter   :: b1=-2.051D-8
    real(8), parameter   :: c0= 1.9898D-4
    real(8), parameter   :: c1=-2.376D-6
    real(8), parameter   :: d = 1.83D-11
    real(8), parameter   :: e =-0.765D-8

    real(8)         :: x,tc,pt,tc2,x2
    real(8)         :: d_x,d_tc,d_pt,d_tc2,d_x2

    if ( t <= 0 ) call utl_abort('gpscompressibility_P0_1: t <= 0')

    x  = gps_p_wa * q / (1.D0+gps_p_wb*q)
    d_x  = 0.0D0
    ! Estimate, from CIPM, Picard (2008)
    tc = t - MPC_K_C_DEGREE_OFFSET_R8
    d_tc = 0.0D0
    pt = p / t
    d_pt = dpdp0 / t
    tc2= tc * tc
    d_tc2= 2 * tc * d_tc
    x2 = x * x
    d_x2 = 2 * x * d_x
    gpscompressibility_P0_1 = &
        -d_pt * (a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2) + &
        2*pt*d_pt*(d+e*x2)
  end function gpscompressibility_P0_1

  !---------------------------------------------------------
  ! gpscompressibility_P0_2
  !---------------------------------------------------------
  function gpscompressibility_P0_2(p,t,q)
    ! gpscompressibility_P0_2 has NO dpdp0 dependency
    implicit none

    ! Arguments
    real(8), intent(in)  :: p
    real(8), intent(in)  :: t
    real(8), intent(in)  :: q

    ! Locals
    real(8)              :: gpscompressibility_P0_2

    real(8), parameter   :: a0= 1.58123D-6
    real(8), parameter   :: a1=-2.9331D-8
    real(8), parameter   :: a2= 1.1043D-10
    real(8), parameter   :: b0= 5.707D-6
    real(8), parameter   :: b1=-2.051D-8
    real(8), parameter   :: c0= 1.9898D-4
    real(8), parameter   :: c1=-2.376D-6
    real(8), parameter   :: d = 1.83D-11
    real(8), parameter   :: e =-0.765D-8

    real(8)         :: x,tc,pt,tc2,x2
    real(8)         :: d_x,d_tc,d_tc2,d_x2

    if ( t <= 0 ) call utl_abort('gpscompressibility_P0_2: t <= 0')

    x  = gps_p_wa * q / (1.D0+gps_p_wb*q)
    d_x  = 0.0D0
    ! Estimate, from CIPM, Picard (2008)
    tc = t - MPC_K_C_DEGREE_OFFSET_R8
    d_tc = 0.0D0
    pt = p / t
    tc2= tc * tc
    d_tc2= 2 * tc * d_tc
    x2 = x * x
    d_x2 = 2 * x * d_x
    gpscompressibility_P0_2 = &
        -pt * (a1*d_tc + a2*d_tc2 + b1*d_tc*x + (b0+b1*tc)*d_x + c1*d_tc*x2 &
              + (c0+c1*tc)*d_x2) + pt*pt*e*d_x2
  end function gpscompressibility_P0_2

end module calcHeightAndPressure_mod

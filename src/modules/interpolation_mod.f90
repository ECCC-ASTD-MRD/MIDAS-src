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

module interpolation_mod
  ! MODULE interpolation_mod (prefix='int' category='4. Data Object transformations')
  !
  ! :Purpose: The grid-point state vector interpolation.
  !
  use mpi, only : mpi_status_size ! this is the mpi library module
  use midasMpi_mod
  use gridstatevector_mod
  use columnData_mod
  use varNameList_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use mathPhysConstants_mod
  use utilities_mod
  implicit none
  save
  private

  ! public subroutines and functions
  public :: int_interp_gsv
  public :: int_hInterp_gsv, int_hInterp_gsv_r4
  public :: int_vInterp_gsv, int_vInterp_gsv_r4
  public :: int_tInterp_gsv
  public :: int_vInterp_col
  public :: int_sint, int_uvint, int_ezgdef, int_cxgaig

  ! module interfaces
  ! -----------------

  interface int_sint
    module procedure int_sint_gsv
    module procedure int_sint_r4_2d
    module procedure int_sint_r8_2d
  end interface int_sint

  interface int_uvint
    module procedure int_uvint_gsv
    module procedure int_uvint_r4_2d
    module procedure int_uvint_r8_2d
  end interface int_uvint

contains

  !--------------------------------------------------------------------------
  ! int_interp_gsv
  !--------------------------------------------------------------------------
  subroutine int_interp_gsv(statevector_in,statevector_out,          &
                             PsfcReference_opt, PsfcReference_r4_opt, &
                             checkModelTop_opt)
    !
    ! :Purpose: high-level interpolation subroutine that proceed with
    !           horizontal and vertical interpolation
    !
    implicit none

    ! Arguments
    type(struct_gsv),       intent(in)    :: statevector_in               ! statevector that will contain the interpolated fields
    type(struct_gsv),       intent(inout) :: statevector_out              ! Reference statevector providing the horizontal and vertical structure

    real(8),      optional, intent(in)    :: PsfcReference_opt(:,:,:)     ! Provides a surface pressure field to be used instead of the first P0 level
    real(4),      optional, intent(in)    :: PsfcReference_r4_opt(:,:,:)  ! Provides a surface pressure field to be used instead of the first P0 level
    
    logical,      optional, intent(in)    :: checkModelTop_opt            ! If true, model top consistency will be checked in vertical interpolation
    
    ! Locals
    logical :: checkModelTop

    character(len=4), pointer :: varNamesToInterpolate(:)

    type(struct_gsv) :: statevector_in_varsLevs, statevector_in_varsLevs_hInterp
    type(struct_gsv) :: statevector_in_hInterp

    !
    !- Error traps
    !
    if (.not.gsv_isAllocated(statevector_in)) then
      call utl_abort('int_interp_gsv: gridStateVector_in not yet allocated')
    end if
    if (.not.gsv_isAllocated(statevector_out)) then
      call utl_abort('int_interp_gsv: gridStateVector_out not yet allocated')
    end if

    !
    !- Do the interpolation of statevector_in onto the grid of statevector_out
    !
    nullify(varNamesToInterpolate)
    call gsv_varNamesList(varNamesToInterpolate, statevector_in)

    !- Horizontal interpolation
    call gsv_allocate(statevector_in_VarsLevs, statevector_in%numstep, &
                      statevector_in%hco, statevector_in%vco,          &
                      mpi_local_opt=statevector_in%mpi_local, mpi_distribution_opt='VarsLevs',  &
                      dataKind_opt=statevector_in%dataKind,                                     &
                      allocHeightSfc_opt=statevector_in%heightSfcPresent, &
                      varNames_opt=varNamesToInterpolate, &
                      hInterpolateDegree_opt=statevector_out%hInterpolateDegree, &
                      hExtrapolateDegree_opt=statevector_out%hExtrapolateDegree )

    call gsv_transposeTilesToVarsLevs( statevector_in, statevector_in_VarsLevs )

    call gsv_allocate(statevector_in_VarsLevs_hInterp, statevector_in%numstep, &
                      statevector_out%hco, statevector_in%vco,  &
                      mpi_local_opt=statevector_out%mpi_local, mpi_distribution_opt='VarsLevs', &
                      dataKind_opt=statevector_out%dataKind,                                    &
                      allocHeightSfc_opt=statevector_out%heightSfcPresent, &
                      varNames_opt=varNamesToInterpolate, &
                      hInterpolateDegree_opt=statevector_out%hInterpolateDegree, &
                      hExtrapolateDegree_opt=statevector_out%hExtrapolateDegree )

    if (statevector_in_VarsLevs%dataKind == 4) then
      call int_hInterp_gsv_r4(statevector_in_VarsLevs, statevector_in_VarsLevs_hInterp)
    else
      call int_hInterp_gsv(statevector_in_VarsLevs, statevector_in_VarsLevs_hInterp)
    end if
    call gsv_deallocate(statevector_in_VarsLevs)

    call gsv_allocate(statevector_in_hInterp, statevector_in%numstep, &
                      statevector_out%hco, statevector_in%vco,      &
                      mpi_local_opt=statevector_out%mpi_local, mpi_distribution_opt='Tiles', &
                      dataKind_opt=statevector_out%dataKind,                                 &
                      allocHeightSfc_opt=statevector_out%heightSfcPresent, &
                      varNames_opt=varNamesToInterpolate )

    call gsv_transposeVarsLevsToTiles( statevector_in_varsLevs_hInterp, statevector_in_hInterp )
    call gsv_deallocate(statevector_in_varsLevs_hInterp)

    !- Vertical interpolation
    
    ! the default is to ensure that the top of the output grid is ~equal or lower than the top of the input grid 
    if ( present(checkModelTop_opt) ) then
      checkModelTop = checkModelTop_opt
    else
      checkModelTop = .true.
    end if
    
    if (statevector_in_VarsLevs%dataKind == 4) then
      call int_vInterp_gsv_r4(statevector_in_hInterp,statevector_out,PsfcReference_opt=PsfcReference_r4_opt, &
                               checkModelTop_opt=checkModelTop)
    else 
      call int_vInterp_gsv(statevector_in_hInterp,statevector_out,PsfcReference_opt=PsfcReference_opt, &
                            checkModelTop_opt=checkModelTop)
    end if

    call gsv_deallocate(statevector_in_hInterp)
    nullify(varNamesToInterpolate)

  end subroutine int_interp_gsv

  !--------------------------------------------------------------------------
  ! int_hInterp_gsv
  !--------------------------------------------------------------------------
  subroutine int_hInterp_gsv(statevector_in,statevector_out)
    !
    ! :Purpose: Horizontal interpolation
    !
    implicit none

    ! Arguments:
    type(struct_gsv),       intent(in)    :: statevector_in   ! statevector that will contain the horizontally interpolated fields
    type(struct_gsv),       intent(inout) :: statevector_out  ! Reference statevector providing the horizontal structure

    ! Locals:
    integer :: varIndex, levIndex, nlev, stepIndex, ierr, kIndex
    character(len=4) :: varName
    character(len=12):: interpolationDegree, extrapolationDegree

    if ( hco_equal(statevector_in%hco,statevector_out%hco) ) then
      write(*,*) 'int_hInterp_gsv: The input and output statevectors are already on same horizontal grids'
      call gsv_copy(statevector_in, statevector_out)
      return
    end if

    if ( .not. vco_equal(statevector_in%vco, statevector_out%vco) ) then
      call utl_abort('int_hInterp_gsv: The input and output statevectors are not on the same vertical levels.')
    end if

    if ( statevector_in%dataKind /= 8 .or. statevector_out%dataKind /= 8 ) then
      call utl_abort('int_hInterp_gsv: Incorrect value for dataKind. Only compatible with dataKind=4')
    end if

    ! set the interpolation degree
    interpolationDegree = statevector_out%hInterpolateDegree
    extrapolationDegree = statevector_out%hExtrapolateDegree

    if ( .not.statevector_in%mpi_local .and. .not.statevector_out%mpi_local ) then

      write(*,*) 'int_hInterp_gsv: before interpolation (no mpi)'

      step_loop: do stepIndex = 1, statevector_out%numStep
        ! copy over some time related parameters
        statevector_out%deet                      = statevector_in%deet
        statevector_out%dateOriginList(stepIndex) = statevector_in%dateOriginList(stepIndex)
        statevector_out%npasList(stepIndex)       = statevector_in%npasList(stepIndex)
        statevector_out%ip2List(stepIndex)        = statevector_in%ip2List(stepIndex)
        statevector_out%etiket                    = statevector_in%etiket

        ! Do horizontal interpolation for mpi global statevectors
        var_loop: do varIndex = 1, vnl_numvarmax
          varName = vnl_varNameList(varIndex)
          if ( .not. gsv_varExist(statevector_in,varName) ) cycle var_loop
          if ( trim(varName) == 'VV' ) cycle var_loop

          nlev = statevector_out%varNumLev(varIndex)

          ! horizontal interpolation

          if ( trim(varName) == 'UU' ) then
            ! interpolate both UV components and keep both UU and VV
            do levIndex = 1, nlev
              ierr = int_uvint( statevector_out, statevector_in, 'BOTH', levIndex, stepIndex, &
                                interpDegree=trim(interpolationDegree), extrapDegree_opt=trim(extrapolationDegree) )
            end do
          else
            ! interpolate scalar variable
            do levIndex = 1, nlev
              ierr = int_sint( statevector_out, statevector_in, varName, levIndex, stepIndex, &
                               interpDegree=trim(interpolationDegree), extrapDegree_opt=trim(extrapolationDegree) )
            end do
          end if
        end do var_loop

      end do step_loop

    else

      write(*,*) 'int_hInterp_gsv: before interpolation (with mpi)'

      if ( statevector_in%mpi_distribution /= 'VarsLevs' .or.   &
          statevector_out%mpi_distribution /= 'VarsLevs' ) then
        call utl_abort('int_hInterp_gsv: The input or output statevector is not distributed by VarsLevs.')
      end if

      do stepIndex = 1, statevector_out%numStep
        k_loop: do kIndex = statevector_in%mykBeg, statevector_in%mykEnd
          varName = gsv_getVarNameFromK(statevector_in,kIndex)
          if ( .not. gsv_varExist(statevector_in,varName) ) cycle k_loop

          ! horizontal interpolation

          if ( trim(varName) == 'UU' ) then
            ! interpolate both UV components and keep UU in main vector
            ierr = int_uvint( statevector_out, statevector_in, 'UU', kIndex, stepIndex, &
                              interpDegree=trim(interpolationDegree), extrapDegree_opt=trim(extrapolationDegree) )
          else if ( trim(varName) == 'VV' ) then
            ! interpolate both UV components and keep VV in main vector
            ierr = int_uvint( statevector_out, statevector_in, 'VV', kIndex, stepIndex, &
                              interpDegree=trim(interpolationDegree), extrapDegree_opt=trim(extrapolationDegree) )
          else
            ! interpolate scalar variable
            ierr = int_sint( statevector_out, statevector_in, 'ALL', kIndex, stepIndex, &
                             interpDegree=trim(interpolationDegree), extrapDegree_opt=trim(extrapolationDegree) )
          end if
        end do k_loop

      end do ! stepIndex

    end if

    if ( gsv_isAssocHeightSfc(statevector_in) .and. gsv_isAssocHeightSfc(statevector_out) ) then
      write(*,*) 'int_hInterp_gsv: interpolating surface height'
      ierr = int_sint( statevector_out, statevector_in, 'ZSFC', 1, 1, &
                       interpDegree=trim(interpolationDegree), extrapDegree_opt=trim(extrapolationDegree) )
    end if

  end subroutine int_hInterp_gsv

  !--------------------------------------------------------------------------
  ! int_hInterp_gsv_r4
  !--------------------------------------------------------------------------
  subroutine int_hInterp_gsv_r4(statevector_in,statevector_out)
    !
    ! :Purpose: Horizontal interpolation
    !
    implicit none

    ! Arguments:
    type(struct_gsv),       intent(in)    :: statevector_in   ! statevector that will contain the horizontally interpolated fields
    type(struct_gsv),       intent(inout) :: statevector_out  ! Reference statevector providing the horizontal structure

    ! Locals:
    integer :: varIndex, levIndex, nlev, stepIndex, ierr, kIndex
    character(len=4) :: varName
    character(len=12):: interpolationDegree, extrapolationDegree

    if ( hco_equal(statevector_in%hco,statevector_out%hco) ) then
      write(*,*) 'int_hInterp_gsv_r4: The input and output statevectors are already on same horizontal grids'
      call gsv_copy(statevector_in, statevector_out)
      return
    end if

    if ( .not. vco_equal(statevector_in%vco, statevector_out%vco) ) then
      call utl_abort('int_hInterp_gsv_r4: The input and output statevectors are not on the same vertical levels.')
    end if

    if ( statevector_in%dataKind /= 4 .or. statevector_out%dataKind /= 4 ) then
      call utl_abort('int_hInterp_gsv_r4: Incorrect value for dataKind. Only compatible with dataKind=4')
    end if

    ! set the interpolation degree
    interpolationDegree = statevector_out%hInterpolateDegree
    extrapolationDegree = statevector_out%hExtrapolateDegree

    if ( .not.statevector_in%mpi_local .and. .not.statevector_out%mpi_local ) then

      write(*,*) 'int_hInterp_gsv_r4: before interpolation (no mpi)'

      step_loop: do stepIndex = 1, statevector_out%numStep
        ! copy over some time related parameters
        statevector_out%deet                      = statevector_in%deet
        statevector_out%dateOriginList(stepIndex) = statevector_in%dateOriginList(stepIndex)
        statevector_out%npasList(stepIndex)       = statevector_in%npasList(stepIndex)
        statevector_out%ip2List(stepIndex)        = statevector_in%ip2List(stepIndex)
        statevector_out%etiket                    = statevector_in%etiket

        ! Do horizontal interpolation for mpi global statevectors
        var_loop: do varIndex = 1, vnl_numvarmax
          varName = vnl_varNameList(varIndex)
          if ( .not. gsv_varExist(statevector_in,varName) ) cycle var_loop
          if ( trim(varName) == 'VV' ) cycle var_loop

          nlev = statevector_out%varNumLev(varIndex)

          ! horizontal interpolation

          if ( trim(varName) == 'UU' ) then
            ! interpolate both UV components and keep both UU and VV
            do levIndex = 1, nlev
              ierr = int_uvint( statevector_out, statevector_in, 'BOTH', levIndex, stepIndex, &
                                interpDegree=trim(interpolationDegree), extrapDegree_opt=trim(extrapolationDegree) )
            end do
          else
            ! interpolate scalar variable
            do levIndex = 1, nlev
              ierr = int_sint( statevector_out, statevector_in, varName, levIndex, stepIndex, &
                               interpDegree=trim(InterpolationDegree), extrapDegree_opt=trim(extrapolationDegree) )
            end do
          end if
        end do var_loop

      end do step_loop

    else

      write(*,*) 'int_hInterp_gsv_r4: before interpolation (with mpi)'

      if ( statevector_in%mpi_distribution /= 'VarsLevs' .or.   &
          statevector_out%mpi_distribution /= 'VarsLevs' ) then
        call utl_abort('int_hInterp_gsv_r4: The input or output statevector is not distributed by VarsLevs.')
      end if

      step2_loop: do stepIndex = 1, statevector_out%numStep
        k_loop: do kIndex = statevector_in%mykBeg, statevector_in%mykEnd
          varName = gsv_getVarNameFromK(statevector_in,kIndex)
          if ( .not. gsv_varExist(statevector_in,varName) ) cycle k_loop

          ! horizontal interpolation

          if ( trim(varName) == 'UU' ) then
            ! interpolate both UV components and keep UU
            ierr = int_uvint( statevector_out, statevector_in, 'UU', kIndex, stepIndex, &
                              interpDegree=trim(interpolationDegree), extrapDegree_opt=trim(extrapolationDegree) )
          else if ( trim(varName) == 'VV' ) then
            ! interpolate both UV components and keep VV
            ierr = int_uvint( statevector_out, statevector_in, 'VV', kIndex, stepIndex, &
                              interpDegree=trim(interpolationDegree), extrapDegree_opt=trim(extrapolationDegree) )
          else
            ! interpolate scalar variable
            ierr = int_sint( statevector_out, statevector_in, 'ALL', kIndex, stepIndex, &
                             interpDegree=trim(InterpolationDegree), extrapDegree_opt=trim(extrapolationDegree) )
          end if
        end do k_loop

      end do step2_loop

    end if

    if ( gsv_isAssocHeightSfc(statevector_in) .and. gsv_isAssocHeightSfc(statevector_out)) then
      write(*,*) 'int_hInterp_gsv_r4: interpolating surface height'
      ierr = int_sint( statevector_out, statevector_in, 'ZSFC', 1, 1, &
                       interpDegree=trim(InterpolationDegree), extrapDegree_opt=trim(extrapolationDegree) )
    end if

  end subroutine int_hInterp_gsv_r4

  !--------------------------------------------------------------------------
  ! int_vInterp_gsv
  !--------------------------------------------------------------------------
  subroutine int_vInterp_gsv(statevector_in,statevector_out,Ps_in_hPa_opt, &
                              PsfcReference_opt,checkModelTop_opt)
    !
    ! :Purpose: Vertical interpolation
    !
    ! :Namelist parameters:
    !         :vInterpCopyLowestLevel:  if true, will overwrite values at the lowest
    !                                   levels to avoid extrapolation
    !
    implicit none

    ! Arguments:
    type(struct_gsv),  intent(in)     :: statevector_in           ! statevector that will contain the vertically interpolated fields
    type(struct_gsv),  intent(inout)  :: statevector_out          ! Reference statevector providing the vertical structure
    logical, optional, intent(in)     :: Ps_in_hPa_opt            ! If true, conversion from hPa to mbar will be done
    real(8), optional, intent(in)     :: PsfcReference_opt(:,:,:) ! Provides a surface pressure field to be used instead of the first P0 level
    logical, optional, intent(in)     :: checkModelTop_opt        ! Model top consistency will be checked prior to interpolation if true

    ! Locals:
    logical :: checkModelTop, vInterpCopyLowestLevel

    integer :: nlev_out, nlev_in, levIndex_out, levIndex_in, latIndex, lonIndex
    integer :: latIndex2, lonIndex2, varIndex, stepIndex
    integer :: status, nulnam, fnom, ierr, fclos
    real(8) :: zwb, zwt

    real(8)           :: psfc_in(statevector_in%myLonBeg:statevector_in%myLonEnd, &
                                 statevector_in%myLatBeg:statevector_in%myLatEnd)
    real(8), pointer  :: pres_out(:,:,:), pres_in(:,:,:), field_out(:,:,:,:), field_in(:,:,:,:)
    real(8), pointer  :: heightSfcIn(:,:), heightSfcOut(:,:)

    character(len=4) :: varName

    type(struct_vco), pointer :: vco_in, vco_out

    NAMELIST /NAMINT/vInterpCopyLowestLevel

    vInterpCopyLowestLevel = .false.
    if ( .not. utl_isNamelistPresent('NAMINT','./flnml') ) then
      if ( mmpi_myid == 0 ) then
        write(*,*) 'int_vInterp_gsv: namint is missing in the namelist.'
        write(*,*) '                     The default values will be taken.'
      end if
    else
      ! Read namelist NAMINT
      nulnam=0
      ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam,nml=namint,iostat=ierr)
      if (ierr.ne.0) call utl_abort('int_vInterp_gsv: Error reading namelist NAMINT')
      if (mmpi_myid.eq.0) write(*,nml=namint)
      ierr=fclos(nulnam)
    end if

    if ( vco_equal(statevector_in%vco, statevector_out%vco) ) then
      write(*,*) 'int_vInterp_gsv: The input and output statevectors are already on same vertical levels'
      call gsv_copy(statevector_in, statevector_out)
      return
    end if

    if ( .not. hco_equal(statevector_in%hco, statevector_out%hco) ) then
      call utl_abort('int_vInterp_gsv: The input and output statevectors are not on the same horizontal grid.')
    end if

    if ( statevector_in%dataKind /= 8 .or. statevector_out%dataKind /= 8 ) then
      call utl_abort('int_vInterp_gsv: Incorrect value for dataKind. Only compatible with dataKind=8')
    end if

    if (gsv_isAssocHeightSfc(statevector_in) .and. gsv_isAssocHeightSfc(statevector_out) ) then
      heightSfcIn => gsv_getHeightSfc(statevector_in)
      heightSfcOut => gsv_getHeightSfc(statevector_out)
      heightSfcOut(:,:) = heightSfcIn(:,:)
    end if

    vco_in => gsv_getVco(statevector_in)
    vco_out => gsv_getVco(statevector_out)

    ! the default is to ensure that the top of the output grid is ~equal or lower than the top of the input grid 
    if ( present(checkModelTop_opt) ) then
      checkModelTop = checkModelTop_opt
    else
      checkModelTop = .true.
    end if
    if (checkModelTop) then
      call vco_ensureCompatibleTops(vco_in, vco_out)
    end if

    step_loop: do stepIndex = 1, statevector_out%numStep

      ! copy over some time related and other parameters
      statevector_out%deet                      = statevector_in%deet
      statevector_out%dateOriginList(stepIndex) = statevector_in%dateOriginList(stepIndex)
      statevector_out%npasList(stepIndex)       = statevector_in%npasList(stepIndex)
      statevector_out%ip2List(stepIndex)        = statevector_in%ip2List(stepIndex)
      statevector_out%etiket                    = statevector_in%etiket
      statevector_out%onPhysicsGrid(:)          = statevector_in%onPhysicsGrid(:)
      statevector_out%hco_physics              => statevector_in%hco_physics

      if ( present(PsfcReference_opt) ) then
        psfc_in(:,:) = PsfcReference_opt(:,:,stepIndex)
      else
        call gsv_getField(statevector_in,field_in,'P0')
        psfc_in(:,:) = field_in(:,:,1,stepIndex)
      end if
      if ( present(Ps_in_hPa_opt) ) then
        if ( Ps_in_hPa_opt ) psfc_in = psfc_in * mpc_pa_per_mbar_r8
      end if

      var_loop: do varIndex = 1, vnl_numvarmax
        varName = vnl_varNameList(varIndex)
        if ( .not. gsv_varExist(statevector_in,varName) ) cycle var_loop

        nlev_in  = statevector_in%varNumLev(varIndex)
        nlev_out = statevector_out%varNumLev(varIndex)

        call gsv_getField(statevector_in ,field_in,varName)
        call gsv_getField(statevector_out,field_out,varName)

        ! for 2D fields and variables on "other" levels, just copy and cycle to next variable
        if ( (nlev_in == 1 .and. nlev_out == 1) .or. &
             (vnl_varLevelFromVarname(varName) == 'OT') ) then
          field_out(:,:,:,stepIndex) = field_in(:,:,:,stepIndex)
          cycle var_loop
        end if

        nullify(pres_out,pres_in)

        ! define pressures on input and output levels
        if (vnl_varLevelFromVarname(varName) == 'TH') then
          status=vgd_levels(vco_in%vgrid,           &
                            ip1_list=vco_in%ip1_T,  &
                            levels=pres_in,         &
                            sfc_field=psfc_in,      &
                            in_log=.false.)
          status=vgd_levels(vco_out%vgrid,           &
                            ip1_list=vco_out%ip1_T,  &
                            levels=pres_out,         &
                            sfc_field=psfc_in,       &
                            in_log=.false.)
        else
          status=vgd_levels(vco_in%vgrid,           &
                            ip1_list=vco_in%ip1_M,  &
                            levels=pres_in,         &
                            sfc_field=psfc_in,      &
                            in_log=.false.)
          status=vgd_levels(vco_out%vgrid,           &
                            ip1_list=vco_out%ip1_M,  &
                            levels=pres_out,         &
                            sfc_field=psfc_in,       &
                            in_log=.false.)
        end if

        ! do the vertical interpolation
        field_out(:,:,:,stepIndex) = 0.0d0
        !$OMP PARALLEL DO PRIVATE(latIndex,latIndex2,lonIndex,lonIndex2,levIndex_in,levIndex_out,zwb,zwt)
        do latIndex = statevector_out%myLatBeg, statevector_out%myLatEnd
          latIndex2 = latIndex - statevector_out%myLatBeg + 1
          do lonIndex = statevector_out%myLonBeg, statevector_out%myLonEnd
            lonIndex2 = lonIndex - statevector_out%myLonBeg + 1
            levIndex_in = 1
            do levIndex_out = 1, nlev_out
              levIndex_in = levIndex_in + 1
              do while(pres_out(lonIndex2,latIndex2,levIndex_out).gt.pres_in(lonIndex2,latIndex2,levIndex_in)  &
                       .and.levIndex_in.lt.nlev_in)
                levIndex_in = levIndex_in + 1
              end do
              levIndex_in = levIndex_in - 1
              zwb = log(pres_out(lonIndex2,latIndex2,levIndex_out)/pres_in(lonIndex2,latIndex2,levIndex_in))  &
                   /log(pres_in(lonIndex2,latIndex2,levIndex_in+1)/pres_in(lonIndex2,latIndex2,levIndex_in))
              zwt = 1.d0 - zwb
              field_out(lonIndex,latIndex,levIndex_out,stepIndex) =   &
                                 zwb*field_in(lonIndex,latIndex,levIndex_in+1,stepIndex) &
                               + zwt*field_in(lonIndex,latIndex,levIndex_in,stepIndex)
            end do
          end do
        end do
        !$OMP END PARALLEL DO

        ! overwrite values at the lowest levels to avoid extrapolation
        if (vInterpCopyLowestLevel) then
          field_out(:,:,nlev_out,stepIndex) = field_in(:,:,nlev_in,stepIndex)
        end if

        deallocate(pres_out)
        deallocate(pres_in)

      end do var_loop

    end do step_loop

  end subroutine int_vInterp_gsv

  !--------------------------------------------------------------------------
  ! int_vInterp_gsv_r4
  !--------------------------------------------------------------------------
  subroutine int_vInterp_gsv_r4(statevector_in,statevector_out,Ps_in_hPa_opt, &
                                 PsfcReference_opt,checkModelTop_opt)
    !
    ! :Purpose: Vertical interpolation of pressure defined fields
    !
    ! :Namelist parameters:
    !         :vInterpCopyLowestLevel:  if true, will overwrite values at the lowest
    !                                   levels to avoid extrapolation
    !
    implicit none

    ! Arguments:
    type(struct_gsv),  intent(in)     :: statevector_in           ! Reference statevector providing the vertical structure
    type(struct_gsv),  intent(inout)  :: statevector_out          ! statevector that will contain the vertically interpolated fields
    logical, optional, intent(in)     :: Ps_in_hPa_opt            ! If true, conversion from hPa to mbar will be done
    real(4), optional, intent(in)     :: PsfcReference_opt(:,:,:) ! Provides a surface pressure field to be used instead of the first P0 level
    logical, optional, intent(in)     :: checkModelTop_opt        ! Model top consistency will be checked prior to interpolation if true

    ! Locals:
    logical :: checkModelTop, vInterpCopyLowestLevel

    integer :: nlev_out, nlev_in, levIndex_out, levIndex_in, latIndex, lonIndex
    integer :: latIndex2, lonIndex2, varIndex, stepIndex
    integer :: status, nulnam, fnom, ierr, fclos
    real(4) :: zwb, zwt

    real(4)           :: psfc_in(statevector_in%myLonBeg:statevector_in%myLonEnd, &
                                 statevector_in%myLatBeg:statevector_in%myLatEnd)
    real(4), pointer  :: pres_out(:,:,:), pres_in(:,:,:), field_out(:,:,:,:), field_in(:,:,:,:)
    real(8), pointer  :: heightSfcIn(:,:), heightSfcOut(:,:)

    character(len=4) :: varName

    type(struct_vco), pointer :: vco_in, vco_out

    NAMELIST /NAMINT/vInterpCopyLowestLevel

    vInterpCopyLowestLevel = .false.
    if ( .not. utl_isNamelistPresent('NAMINT','./flnml') ) then
      if ( mmpi_myid == 0 ) then
        write(*,*) 'int_vInterp_gsv_r4: namint is missing in the namelist.'
        write(*,*) '                     The default values will be taken.'
      end if
    else
      ! Read namelist NAMINT
      nulnam=0
      ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam,nml=namint,iostat=ierr)
      if (ierr.ne.0) call utl_abort('int_vInterp_gsv_r4: Error reading namelist NAMINT')
      if (mmpi_myid.eq.0) write(*,nml=namint)
      ierr=fclos(nulnam)
    end if

    if ( vco_equal(statevector_in%vco, statevector_out%vco) ) then
      write(*,*) 'int_vInterp_gsv_r4: The input and output statevectors are already on same vertical levels'
      call gsv_copy(statevector_in, statevector_out)
      return
    end if

    if ( .not. hco_equal(statevector_in%hco, statevector_out%hco) ) then
      call utl_abort('int_vInterp_gsv_r4: The input and output statevectors are not on the same horizontal grid.')
    end if

    if ( statevector_in%dataKind /= 4 .or. statevector_out%dataKind /= 4 ) then
      call utl_abort('int_vInterp_gsv_r4: Incorrect value for dataKind. Only compatible with dataKind=4')
    end if

    if ( gsv_isAssocHeightSfc(statevector_in) .and. gsv_isAssocHeightSfc(statevector_out) ) then
      heightSfcIn => gsv_getHeightSfc(statevector_in)
      heightSfcOut => gsv_getHeightSfc(statevector_out)
      HeightSfcOut(:,:) = HeightSfcIn(:,:)
    end if

    vco_in => gsv_getVco(statevector_in)
    vco_out => gsv_getVco(statevector_out)

    ! the default is to ensure that the top of the output grid is ~equal or lower than the top of the input grid 
    if ( present(checkModelTop_opt) ) then
      checkModelTop = checkModelTop_opt
    else
      checkModelTop = .true.
    end if
    if (checkModelTop) then
      write(*,*) 'int_vInterp_gsv_r4: Checking that that the top of the destination grid is not higher than the top of the source grid.'
      call vco_ensureCompatibleTops(vco_in, vco_out)
    end if

    step_loop: do stepIndex = 1, statevector_out%numStep

      ! copy over some time related and other parameters
      statevector_out%deet                      = statevector_in%deet
      statevector_out%dateOriginList(stepIndex) = statevector_in%dateOriginList(stepIndex)
      statevector_out%npasList(stepIndex)       = statevector_in%npasList(stepIndex)
      statevector_out%ip2List(stepIndex)        = statevector_in%ip2List(stepIndex)
      statevector_out%etiket                    = statevector_in%etiket
      statevector_out%onPhysicsGrid(:)          = statevector_in%onPhysicsGrid(:)
      statevector_out%hco_physics              => statevector_in%hco_physics

      if ( present(PsfcReference_opt) ) then
        psfc_in(:,:) = PsfcReference_opt(:,:,stepIndex)
      else
        call gsv_getField(statevector_in,field_in,'P0')
        psfc_in(:,:) = field_in(:,:,1,stepIndex)
      end if
      if ( present(Ps_in_hPa_opt) ) then
        if ( Ps_in_hPa_opt ) psfc_in = psfc_in * mpc_pa_per_mbar_r4
      end if

      var_loop: do varIndex = 1, vnl_numvarmax
        varName = vnl_varNameList(varIndex)
        if ( .not. gsv_varExist(statevector_in,varName) ) cycle var_loop

        nlev_in  = statevector_in%varNumLev(varIndex)
        nlev_out = statevector_out%varNumLev(varIndex)

        call gsv_getField(statevector_in ,field_in,varName)
        call gsv_getField(statevector_out,field_out,varName)

        ! for 2D fields and variables on "other" levels, just copy and cycle to next variable
        if ( (nlev_in == 1 .and. nlev_out == 1) .or. &
             (vnl_varLevelFromVarname(varName) == 'OT') ) then
          field_out(:,:,:,stepIndex) = field_in(:,:,:,stepIndex)
          cycle var_loop
        end if

        nullify(pres_out,pres_in)

        ! define pressures on input and output levels
        if (vnl_varLevelFromVarname(varName) == 'TH') then
          status=vgd_levels(vco_in%vgrid,           &
                            ip1_list=vco_in%ip1_T,  &
                            levels=pres_in,         &
                            sfc_field=psfc_in,      &
                            in_log=.false.)
          status=vgd_levels(vco_out%vgrid,           &
                            ip1_list=vco_out%ip1_T,  &
                            levels=pres_out,         &
                            sfc_field=psfc_in,       &
                            in_log=.false.)
        else
          status=vgd_levels(vco_in%vgrid,           &
                            ip1_list=vco_in%ip1_M,  &
                            levels=pres_in,         &
                            sfc_field=psfc_in,      &
                            in_log=.false.)
          status=vgd_levels(vco_out%vgrid,           &
                            ip1_list=vco_out%ip1_M,  &
                            levels=pres_out,         &
                            sfc_field=psfc_in,       &
                            in_log=.false.)
        end if

        ! do the vertical interpolation
        field_out(:,:,:,stepIndex) = 0.0
        !$OMP PARALLEL DO PRIVATE(latIndex,latIndex2,lonIndex,lonIndex2,levIndex_in,levIndex_out,zwb,zwt)
        do latIndex = statevector_out%myLatBeg, statevector_out%myLatEnd
          latIndex2 = latIndex - statevector_out%myLatBeg + 1
          do lonIndex = statevector_out%myLonBeg, statevector_out%myLonEnd
            lonIndex2 = lonIndex - statevector_out%myLonBeg + 1
            levIndex_in = 1
            do levIndex_out = 1, nlev_out
              levIndex_in = levIndex_in + 1
              do while(pres_out(lonIndex2,latIndex2,levIndex_out).gt.pres_in(lonIndex2,latIndex2,levIndex_in)  &
                       .and.levIndex_in.lt.nlev_in)
                levIndex_in = levIndex_in + 1
              end do
              levIndex_in = levIndex_in - 1
              zwb = log(pres_out(lonIndex2,latIndex2,levIndex_out)/pres_in(lonIndex2,latIndex2,levIndex_in))  &
                   /log(pres_in(lonIndex2,latIndex2,levIndex_in+1)/pres_in(lonIndex2,latIndex2,levIndex_in))
              zwt = 1.0 - zwb
              field_out(lonIndex,latIndex,levIndex_out,stepIndex) =   &
                                 zwb*field_in(lonIndex,latIndex,levIndex_in+1,stepIndex) &
                               + zwt*field_in(lonIndex,latIndex,levIndex_in,stepIndex)
            end do
          end do
        end do
        !$OMP END PARALLEL DO

        ! overwrite values at the lowest levels to avoid extrapolation
        if (vInterpCopyLowestLevel) then
          field_out(:,:,nlev_out,stepIndex) = field_in(:,:,nlev_in,stepIndex)
        end if

        deallocate(pres_out)
        deallocate(pres_in)

      end do var_loop

    end do step_loop

  end subroutine int_vInterp_gsv_r4

  !--------------------------------------------------------------------------
  ! int_tInterp_gsv
  !--------------------------------------------------------------------------
  subroutine int_tInterp_gsv(statevector_in,statevector_out)
    !
    ! :Purpose: Time interpolation from statevector with low temporal resolution
    !           to statevector with high temporal resolution.
    !
    implicit none

    ! Arguments:
    type(struct_gsv),  intent(in)     :: statevector_in           ! statevector that will contain the temporal interpolated fields
    type(struct_gsv),  intent(inout)  :: statevector_out          ! Reference statevector providing the temporal structure

    ! Locals:
    integer :: kIndex, latIndex, lonIndex
    integer :: stepIndexIn1, stepIndexIn2, stepIndexOut, numStepIn, numStepOut
    integer :: lon1, lon2, lat1, lat2, k1, k2
    integer :: dateStampIn, dateStampOut 
    real(8) :: weight1, weight2
    real(8) :: deltaHour, deltaHourInOut

    real(4), pointer  :: gdIn_r4(:,:,:,:), gdOut_r4(:,:,:,:)
    real(8), pointer  :: gdIn_r8(:,:,:,:), gdOut_r8(:,:,:,:)

    write(*,*) 'int_tInterp_gsv: STARTING'

    if ( .not. vco_equal(statevector_in%vco, statevector_out%vco) ) then
      call utl_abort('int_tInterp_gsv: The input and output statevectors are not on the same vertical levels')
    end if

    if ( .not. hco_equal(statevector_in%hco, statevector_out%hco) ) then
      call utl_abort('int_tInterp_gsv: The input and output statevectors are not on the same horizontal grid.')
    end if

    if ( statevector_in%numStep > statevector_out%numStep ) then
      write(*,*) 'int_tInterp_gsv: numStep_out is less than numStep_in, calling gsv_copy.'
      call gsv_copy(statevector_in, statevector_out, allowTimeMismatch_opt=.true.)
      return
    else if ( statevector_in%numStep == statevector_out%numStep ) then
      write(*,*) 'int_tInterp_gsv: numStep_out is equal to numStep_in, calling gsv_copy.'
      call gsv_copy(statevector_in, statevector_out, allowVarMismatch_opt=.true.)
      return
    end if

    lon1 = statevector_in%myLonBeg
    lon2 = statevector_in%myLonEnd
    lat1 = statevector_in%myLatBeg
    lat2 = statevector_in%myLatEnd
    k1 = statevector_in%mykBeg
    k2 = statevector_in%mykEnd

    numStepIn = statevector_in%numStep
    numStepOut = statevector_out%numStep
    if ( mmpi_myid == 0 ) then
      write(*,*) 'int_tInterp_gsv: numStepIn=', numStepIn, ',numStepOut=',numStepOut
    end if

    ! compute positive deltaHour between two first stepIndex of statevector_in (input temporal grid). 
    ! If numStepIn == 1, no time interpolation needed (weights are set to zero).
    if ( numStepIn > 1 ) then
      call difdatr(statevector_in%dateStampList(1), statevector_in%dateStampList(2), deltaHour)
      deltaHour = abs(deltaHour)
    else
      deltaHour = 0.0d0
    end if

    do stepIndexOut = 1, numStepOut
      if ( numStepIn == 1 ) then
        stepIndexIn1 = 1
        stepIndexIn2 = 1
        weight1 = 1.0d0 
        weight2 = 0.0d0 
      else
        ! find staevector_in%dateStamp on the left and right
        if ( statevector_in%dateStampList(numStepIn) == statevector_out%dateStampList(stepIndexOut) ) then
          stepIndexIn2 = numStepIn
        else
          stepInLoop: do stepIndexIn2 = 1, numStepIn
            dateStampIn = statevector_in%dateStampList(stepIndexIn2)
            dateStampOut = statevector_out%dateStampList(stepIndexOut)
            call difdatr(dateStampIn, dateStampOut, deltaHourInOut)
            if ( deltaHourInOut > 0.0d0 ) exit stepInLoop
          end do stepInLoop
        end if
        stepIndexIn1 = stepIndexIn2 - 1

        ! compute deltaHour between left stepIndex of statevector_in and statevector_out
        dateStampIn = statevector_in%dateStampList(stepIndexIn1)
        dateStampOut = statevector_out%dateStampList(stepIndexOut)
        call difdatr(dateStampOut, dateStampIn, deltaHourInOut)
        if ( deltaHourInOut < 0.0d0 ) then 
          call utl_abort('int_tInterp_gsv: deltaHourInOut should be greater or equal to 0')
        end if

        ! compute the interpolation weights for left stepIndex of statevector_in (weight1) and right stepIndex of statevector_in (weight2) 
        weight1 = 1.0d0 - deltaHourInOut / deltaHour
        weight2 = deltaHourInOut / deltaHour
      end if

      if ( mmpi_myid == 0 ) then
        write(*,*) 'int_tInterp_gsv: for stepIndexOut=', stepIndexOut, &
                   ',stepIndexIn1=', stepIndexIn1, ',stepIndexIn2=', stepIndexIn2, &
                   ',weight1=', weight1, ',weight2=', weight2, &
                   ',deltaHourInOut/deltaHour=', deltaHourInOut,'/',deltaHour
      end if

      if ( statevector_in%dataKind == 4 .and. statevector_out%dataKind == 4 ) then
        call gsv_getField(statevector_in, gdIn_r4)
        call gsv_getField(statevector_out, gdOut_r4)
        !$OMP PARALLEL DO PRIVATE (latIndex,kIndex,lonIndex)
        do kIndex = k1, k2
          do latIndex = lat1, lat2
            do lonIndex = lon1, lon2
              gdOut_r4(lonIndex,latIndex,kIndex,stepIndexOut) =  &
                real(weight1,4) * gdIn_r4(lonIndex,latIndex,kIndex,stepIndexIn1) + &
                real(weight2,4) * gdIn_r4(lonIndex,latIndex,kIndex,stepIndexIn2)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      else if ( statevector_in%dataKind == 4 .and. statevector_out%dataKind == 8 ) then
        call gsv_getField(statevector_in, gdIn_r4)
        call gsv_getField(statevector_out, gdOut_r8)
        !$OMP PARALLEL DO PRIVATE (latIndex,kIndex,lonIndex)
        do kIndex = k1, k2
          do latIndex = lat1, lat2
            do lonIndex = lon1, lon2
              gdOut_r8(lonIndex,latIndex,kIndex,stepIndexOut) =  &
                weight1 * real(gdIn_r4(lonIndex,latIndex,kIndex,stepIndexIn1),8) + &
                weight2 * real(gdIn_r4(lonIndex,latIndex,kIndex,stepIndexIn2),8)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      else if ( statevector_in%dataKind == 8 .and. statevector_out%dataKind == 4 ) then
        call gsv_getField(statevector_in, gdIn_r8)
        call gsv_getField(statevector_out, gdOut_r4)
        !$OMP PARALLEL DO PRIVATE (latIndex,kIndex,lonIndex)
        do kIndex = k1, k2
          do latIndex = lat1, lat2
            do lonIndex = lon1, lon2
              gdOut_r4(lonIndex,latIndex,kIndex,stepIndexOut) =  &
                real(weight1 * gdIn_r8(lonIndex,latIndex,kIndex,stepIndexIn1),4) + &
                real(weight2 * gdIn_r8(lonIndex,latIndex,kIndex,stepIndexIn2),4)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      else if ( statevector_in%dataKind == 8 .and. statevector_out%dataKind == 8 ) then
        call gsv_getField(statevector_in, gdIn_r8)
        call gsv_getField(statevector_out, gdOut_r8)
        !$OMP PARALLEL DO PRIVATE (latIndex,kIndex,lonIndex)
        do kIndex = k1, k2
          do latIndex = lat1, lat2
            do lonIndex = lon1, lon2
              gdOut_r8(lonIndex,latIndex,kIndex,stepIndexOut) =  &
                weight1 * gdIn_r8(lonIndex,latIndex,kIndex,stepIndexIn1) + &
                weight2 * gdIn_r8(lonIndex,latIndex,kIndex,stepIndexIn2)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end if

    end do

    write(*,*) 'int_tInterp_gsv: END'

  end subroutine int_tInterp_gsv

  !--------------------------------------------------------------------------
  ! int_vInterp_col
  !--------------------------------------------------------------------------
  subroutine int_vInterp_col(column_in,column_out,varName,useColumnPressure_opt)
    !
    ! :Purpose: Vertical interpolation of a columData object
    !
    implicit none

    ! Arguments:
    type(struct_columnData),  intent(in)    :: column_in              ! Reference columnData providing the vertical structure
    type(struct_columnData),  intent(inout) :: column_out             ! columnData that will contain the vertically interpolated column
    character(len=*),         intent(in)    :: varName
    logical, optional,        intent(in)    :: useColumnPressure_opt  ! if .true. use P_* instead of the pressure provided by vgd_levels 

    ! Locals:
    real(8), pointer :: pres_in(:), pres_out(:)
    real(8), pointer :: pres1D_in(:)  , pres1D_out(:)
    real(8), pointer :: pres3D_in(:,:,:), pres3D_out(:,:,:)
    character(len=4) :: varLevel
    real(8)          :: zwb, zwt
    integer          :: levIndex_out, levIndex_in, columnIndex, status
    logical          :: vInterp, useColumnPressure

    integer, allocatable, target :: THlevelWanted(:), MMlevelWanted(:)
    integer, pointer :: levelWanted(:)
    real(8), allocatable :: psfc_in(:,:)

    varLevel = vnl_varLevelFromVarname(varName)

    if ( present(useColumnPressure_opt) ) then
      useColumnPressure = useColumnPressure_opt
    else
      useColumnPressure = .true.
    end if

    nullify(pres3D_in)
    nullify(pres3D_out)

    vInterp = .true.
    if ( .not. col_varExist(column_in,'P0' ) ) then
      write(*,*)
      write(*,*) 'int_vInterp_col: P0 is missing. Vertical interpolation WILL NOT BE PERFORMED'
      vInterp = .false.
    else if ( col_getNumLev(column_in ,'TH') <= 1 .or. &
              col_getNumLev(column_in ,'MM') <= 1 ) then
      vInterp = .false.
      write(*,*)
      write(*,*) 'int_vInterp_col: The input backgrounds are 2D. Vertical interpolation WILL NOT BE PERFORMED'
    end if

    if_vInterp: if (vInterp) then
      if ( .not. useColumnPressure ) then

        ! read psfc_in to use in vgd_levels
        allocate(psfc_in(1,col_getNumCol(column_in)))
        do columnIndex = 1,col_getNumCol(column_in)
          psfc_in(1,columnIndex) = col_getElem(column_out,1,columnIndex,'P0')
        end do

        ! Compute pressure
        if ( varLevel == 'TH' ) then
          ! pres3D_in
          status = vgd_levels(column_in%vco%vgrid,ip1_list=column_in%vco%ip1_T,  &
                              levels=pres3D_in,sfc_field=psfc_in,in_log=.false.)
          if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_levels')

          ! pres3D_out
          status = vgd_levels(column_out%vco%vgrid,ip1_list=column_out%vco%ip1_T,  &
                              levels=pres3D_out,sfc_field=psfc_in,in_log=.false.)
          if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_levels')

        else if ( varLevel == 'MM' ) then
          ! pres3D_in
          status = vgd_levels(column_in%vco%vgrid,ip1_list=column_in%vco%ip1_M,  &
                              levels=pres3D_in,sfc_field=psfc_in,in_log=.false.)
          if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_levels')

          ! pres3D_out
          status = vgd_levels(column_out%vco%vgrid,ip1_list=column_out%vco%ip1_M,  &
                              levels=pres3D_out,sfc_field=psfc_in,in_log=.false.)
          if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_levels')

        else
          call utl_abort('int_vInterp_col: only varLevel TH/MM is allowed')
        end if  ! varLevel
      end if    ! useColumnPressure

      do columnIndex = 1, col_getNumCol(column_out)

        ! pres1D_in
        if ( .not. useColumnPressure ) then
          pres1D_in => pres3D_in(1,columnIndex,:)
        else
          if ( varLevel == 'TH' ) then
            pres1D_in => col_getColumn(column_in,columnIndex,'P_T')
          else if ( varLevel == 'MM' ) then
            pres1D_in => col_getColumn(column_in,columnIndex,'P_M')
          else
            call utl_abort('int_vInterp_col: only varLevel TH/MM is allowed')
          end if
        end if

        ! pres1D_out
        if ( .not. useColumnPressure ) then
          pres1D_out => pres3D_out(1,columnIndex,:)
        else
          if ( varLevel == 'TH' ) then
            pres1D_out => col_getColumn(column_out,columnIndex,'P_T')
          else if ( varLevel == 'MM' ) then
            pres1D_out => col_getColumn(column_out,columnIndex,'P_M')
          else
            call utl_abort('int_vInterp_col: only varLevel TH/MM is allowed')
          end if
        end if

        pres_in  => col_getColumn(column_in ,columnIndex,varName)
        pres_out => col_getColumn(column_out,columnIndex,varName)

        if ( mmpi_myid == 0 .and. columnIndex == 1 .and. &
             (trim(varName) == 'P_T' ) ) then

          write(*,*) 'useColumnPressure=', useColumnPressure

          write(*,*) 'int_vInterp_col, COLUMN_IN(1):'
          write(*,*) trim(varName),':'
          write(*,*) pres_in(:)

          write(*,*) 'int_vInterp_col, COLUMN_OUT(1):'
          write(*,*) trim(varName),':'
          write(*,*) pres_out(:)
          write(*,*)
        end if

        levIndex_in = 1
        do levIndex_out = 1, col_getNumLev(column_out,varLevel)
          levIndex_in = levIndex_in + 1
          do while( pres1D_out(levIndex_out) .gt. pres1D_in(levIndex_in) .and. &
               levIndex_in .lt. col_getNumLev(column_in,varLevel) )
            levIndex_in = levIndex_in + 1
          end do
          levIndex_in = levIndex_in - 1
          zwb = log(pres1D_out(levIndex_out)/pres1D_in(levIndex_in))/  &
               log(pres1D_in(levIndex_in+1)/pres1D_in(levIndex_in))
          zwt = 1. - zwb
          if (  useColumnPressure .and. &
              (trim(varName) == 'P_T' .or. trim(varName) == 'P_M' ) ) then
            ! do nothing, i.e. use the pressures from column_in
          else if ( .not. useColumnPressure .and. &
              (trim(varName) == 'P_T' .or. trim(varName) == 'P_M' ) ) then
            pres_out(levIndex_out) = exp(zwb*log(pres_in(levIndex_in+1)) + zwt*log(pres_in(levIndex_in)))
          else
            pres_out(levIndex_out) = zwb*pres_in(levIndex_in+1) + zwt*pres_in(levIndex_in)
          end if
        end do
      end do

      if ( .not. useColumnPressure ) then
        deallocate(pres3D_in)
        deallocate(pres3D_out)
        deallocate(psfc_in)
      end if

    else if_vInterp

      if (column_out%vco%nlev_T > 0 .and. column_out%vco%nlev_M > 0) then

        ! Find which levels in column_in matches column_out
        allocate(THlevelWanted(column_out%vco%nlev_T))
        allocate(MMlevelWanted(column_out%vco%nlev_M))

        call vco_levelMatchingList( THlevelWanted, MMlevelWanted, & ! OUT
                                    column_out%vco, column_in%vco ) ! IN

        if ( any(THlevelWanted == -1) .or. any(MMlevelWanted == -1) ) then
          call utl_abort('int_vInterp_col: column_out is not a subsets of column_in!')
        end if

        ! Transfer the corresponding data
        do columnIndex = 1, col_getNumCol(column_out)
          pres_in  => col_getColumn(column_in ,columnIndex,varName)
          pres_out => col_getColumn(column_out,columnIndex,varName)
          if (vnl_varLevelFromVarname(varName) == 'TH') then
            levelWanted => THlevelWanted
          else
            levelWanted => MMlevelWanted
          end if
          do levIndex_out = 1, col_getNumLev(column_out,varLevel)
            pres_out(levIndex_out) = pres_in(levelWanted(levIndex_out))
          end do
        end do

        deallocate(THlevelWanted)
        deallocate(MMlevelWanted)

      else if (column_out%vco%nlev_depth > 0) then
        write(*,*) 'vco_levelMatchingList: no MM and TH levels, but depth levels exist'
        if (any(column_out%vco%depths(:) /= column_in%vco%depths(:))) then
          call utl_abort('int_vInterp_col: some depth levels not equal')
        else
          ! copy over depth levels
          do columnIndex = 1, col_getNumCol(column_out)
            pres_in  => col_getColumn(column_in ,columnIndex,varName)
            pres_out => col_getColumn(column_out,columnIndex,varName)
            do levIndex_out = 1, col_getNumLev(column_out,varLevel)
              pres_out(levIndex_out) = pres_in(levIndex_out)
            end do
          end do
        end if

      end if

    end if if_vInterp

  end subroutine int_vInterp_col



  subroutine int_setezopt(interpDegree, extrapDegree_opt)
    implicit none

    ! arguments
    character(len=*) :: interpDegree
    character(len=*), optional :: extrapDegree_opt

    ! locals
    character(len=12) :: extrapDegree
    integer           :: ierr, ezsetopt, ezsetval

    if ( trim(interpDegree) /= 'LINEAR' .and. &
         trim(interpDegree) /= 'CUBIC' .and. &
         trim(interpDegree) /= 'NEAREST' ) then
      write(*,*) 'int_setezopt: interpDegree = ', trim(interpDegree)
      call utl_abort('int_setezopt: invalid interpolation degree')
    end if

    if ( present(extrapDegree_opt) ) then
      extrapDegree = extrapDegree_opt
    else
      extrapDegree = 'VALUE'
    end if

    ierr = ezsetopt('INTERP_DEGREE', interpDegree)
    if ( trim(extrapDegree) == 'VALUE' ) then
      ierr = ezsetval('EXTRAP_VALUE', 0.0)
    end if
    ierr = ezsetopt('EXTRAP_DEGREE', extrapDegree)

  end subroutine int_setezopt


  function int_sint_gsv(stateVectorOut, stateVectorIn, varName, levIndex, stepIndex, &
                        interpDegree, extrapDegree_opt) result(ierr)
    implicit none

    ! arguments
    type(struct_gsv) :: stateVectorOut
    type(struct_gsv) :: stateVectorIn
    character(len=*) :: varName
    integer          :: levIndex
    integer          :: stepIndex
    character(len=*)           :: interpDegree
    character(len=*), optional :: extrapDegree_opt
    integer :: ierr

    ! locals
    real(4), pointer :: zout4(:,:,:,:), zin4(:,:,:,:)
    real(8), pointer :: zout8(:,:,:,:), zin8(:,:,:,:)
    real(8), pointer :: heightSfcOut(:,:), heightSfcIn(:,:)
    integer :: ezsint, ezdefset

    write(*,*) 'int_sint_gsv: starting'

    ! check if special interpolation is required
    if (stateVectorIn%hco%initialized .and. stateVectorOut%hco%initialized) then
      if (stateVectorIn%hco%grtyp == 'Y') then
        ! for now, only comptable for real(4)
        if(stateVectorOut%dataKind /= 4 .or. stateVectorIn%dataKind /= 4) then
          call utl_abort('int_sint_gsv: cloudToGrid only implemented for real(4)')
        end if
        ierr = int_sintCloudToGrid_gsv(stateVectorOut, stateVectorIn, varName, levIndex, stepIndex)
        return
      end if
    end if

    ! do the standard interpolation

    ierr = ezdefset(stateVectorOut%hco%EZscintID, stateVectorIn%hco%EZscintID)
    call int_setezopt(interpDegree, extrapDegree_opt)   

    if (trim(varName) == 'ZSFC') then

      heightSfcIn  => gsv_getHeightSfc(stateVectorIn)
      heightSfcOut => gsv_getHeightSfc(stateVectorOut)

      ! allocate real(4) buffers and copy to/from for interpolation
      allocate(zin4(stateVectorIn%hco%ni,stateVectorIn%hco%nj,1,1))
      allocate(zout4(stateVectorOut%hco%ni,stateVectorOut%hco%nj,1,1))
      zin4(:,:,1,1) = heightSfcIn(:,:)
      ierr = ezsint(zout4(:,:,1,1),zin4(:,:,1,1))
      heightSfcOut(:,:) = zout4(:,:,1,1)
      deallocate(zin4,zout4)

    else if (stateVectorOut%dataKind == 4 .and. stateVectorIn%dataKind == 4) then

      if (trim(varName) == 'ALL') then
        call gsv_getField(stateVectorOut, zout4)
        call gsv_getField(stateVectorIn,  zin4)
     else
        call gsv_getField(stateVectorOut, zout4, varName)
        call gsv_getField(stateVectorIn,  zin4,  varName)
      end if

      ierr = ezsint(zout4(:,:,levIndex,stepIndex),zin4(:,:,levIndex,stepIndex))

    else if (stateVectorOut%dataKind == 8 .and. stateVectorIn%dataKind == 8) then

      if (trim(varName) == 'ALL') then
        call gsv_getField(stateVectorOut, zout8)
        call gsv_getField(stateVectorIn,  zin8)
      else
        call gsv_getField(stateVectorOut, zout8, varName)
        call gsv_getField(stateVectorIn,  zin8,  varName)
      end if

      ! allocate real(4) buffers and copy to/from for interpolation
      allocate(zin4(stateVectorIn%hco%ni,stateVectorIn%hco%nj,1,1))
      allocate(zout4(stateVectorOut%hco%ni,stateVectorOut%hco%nj,1,1))
      zin4(:,:,1,1) = zin8(:,:,levIndex,stepIndex)
      ierr = ezsint(zout4(:,:,1,1),zin4(:,:,1,1))
      zout8(:,:,levIndex,stepIndex) = zout4(:,:,1,1)
      deallocate(zin4,zout4)

    else

      call utl_abort('int_sint_gsv: not implemented for mixed dataKind')

    end if

  end function int_sint_gsv


  function int_sint_r4_2d(zout4, zin4, hco_out, hco_in, interpDegree, extrapDegree_opt) result(ierr)
    implicit none

    ! arguments
    real(4) :: zout4(:,:), zin4(:,:)
    type(struct_hco), pointer :: hco_out, hco_in
    integer :: ierr
    character(len=*)           :: interpDegree
    character(len=*), optional :: extrapDegree_opt

    ! locals
    integer :: ezsint, ezdefset

    ! check if special interpolation is required
    if (hco_in%initialized .and. hco_out%initialized) then
      ierr = ezdefset(hco_out%EZscintID, hco_in%EZscintID)
      if (hco_in%grtyp == 'Y') then
        call utl_abort('int_sint_r4_2d: not implemented for cloud to grid')
        return
      end if
    end if

    ! do the standard interpolation
    call int_setezopt(interpDegree, extrapDegree_opt)   
    ierr = ezsint(zout4,zin4)

  end function int_sint_r4_2d


  function int_sintCloudToGrid_gsv(stateVectorGrid, stateVectorCloud, varName, levIndex, stepIndex) result(ierr)
    implicit none

    ! arguments
    type(struct_gsv)          :: stateVectorGrid, stateVectorCloud
    character(len=*)          :: varName
    integer                   :: levIndex, stepIndex
    integer                   :: ierr

    ! locals
    integer :: gdxyfll
    integer :: niCloud, njCloud, niGrid, njGrid
    integer :: top, bottom, left, right, np, lonIndexCloud, latIndexCloud, k, l, m, lonIndexGrid, latIndexGrid
    integer :: in(8), jn(8), ngp, nfill, nhole, nextrap0, nextrap1
    integer, allocatable :: icount(:,:), maskGrid(:,:), maskCloud(:,:)
    real(4), pointer     :: fieldCloud_4d(:,:,:,:), fieldGrid_4d(:,:,:,:)
    real(4), pointer     :: fieldCloud(:,:), fieldGrid(:,:)
    real(4), allocatable :: xCloud(:,:), yCloud(:,:)
    real(4) :: extrapValue

    extrapValue = 0.0

    niCloud = stateVectorCloud%hco%ni
    njCloud = stateVectorCloud%hco%nj
    niGrid  = stateVectorGrid%hco%ni
    njGrid  = stateVectorGrid%hco%nj

    if (trim(varName) == 'ALL') then
      call gsv_getField(stateVectorGrid,  fieldGrid_4d)
      call gsv_getField(stateVectorCloud, fieldCloud_4d)
    else
      call gsv_getField(stateVectorGrid,  fieldGrid_4d,  varName)
      call gsv_getField(stateVectorCloud, fieldCloud_4d, varName)
    end if
    fieldGrid  => fieldGrid_4d(:,:,levIndex,stepIndex)
    fieldCloud => fieldCloud_4d(:,:,levIndex,stepIndex)

    allocate(xCloud(niCloud, njCloud))
    allocate(yCloud(niCloud, njCloud))
    allocate(icount(niGrid, njGrid))
    allocate(maskCloud(niCloud, njCloud))
    allocate(maskGrid(niGrid, njGrid))

    ! set masks based on oceanMasks (if present)
    if (stateVectorCloud%oceanMask%maskPresent) then
      maskCloud(:,:) = 0
      if (trim(varName) == 'ALL') then
        where(stateVectorCloud%oceanMask%mask(:,:,gsv_getLevFromK(stateVectorCloud,levIndex))) maskCloud(:,:) = 1
      else
        where(stateVectorCloud%oceanMask%mask(:,:,levIndex)) maskCloud(:,:) = 1
      end if
    else
      maskCloud(:,:) = 1
    end if
    maskGrid(:,:)  = 1

    ! Calcul des pos. x-y des eclairs sur la grille modele
    ierr = gdxyfll(stateVectorGrid%hco%EZscintID, xCloud, yCloud, &
                   stateVectorCloud%hco%lat2d_4*MPC_DEGREES_PER_RADIAN_R4, &
                   stateVectorCloud%hco%lon2d_4*MPC_DEGREES_PER_RADIAN_R4, &
                   stateVectorCloud%hco%ni*stateVectorCloud%hco%nj)

    write(*,*) 'xCloud, yCloud, lonCloud, latCloud = ', xCloud(1,1), yCloud(1,1),  &
               stateVectorCloud%hco%lon2d_4(1,1)*MPC_DEGREES_PER_RADIAN_R4, &
               stateVectorCloud%hco%lat2d_4(1,1)*MPC_DEGREES_PER_RADIAN_R4

    write(*,*) 'shape(fieldCloud) = ', shape(fieldCloud), niCloud, njCloud

    fieldGrid(:,:) = 0.0
    icount(:,:) = 0

    do latIndexCloud = 1, njCloud
      do lonIndexCloud = 1, niCloud
        lonIndexGrid = nint(xCloud(lonIndexCloud,latIndexCloud))
        latIndexGrid = nint(yCloud(lonIndexCloud,latIndexCloud))
        if ( lonIndexGrid >= 1 .and. latIndexGrid >= 1 .and. &
             lonIndexGrid <= niGrid .and. latIndexGrid <= njGrid .and. &
             maskCloud(lonIndexCloud,latIndexCloud) == 1 ) then
          fieldGrid(lonIndexGrid,latIndexGrid) = fieldGrid(lonIndexGrid,latIndexGrid) +  &
                                                 fieldCloud(lonIndexCloud,latIndexCloud)
          icount(lonIndexGrid,latIndexGrid) = icount(lonIndexGrid,latIndexGrid) + 1
        end if
      end do
    end do

    do latIndexGrid = 1, njGrid
      do lonIndexGrid = 1, niGrid
        if(icount(lonIndexGrid,latIndexGrid) > 0) then
          fieldGrid(lonIndexGrid,latIndexGrid) = fieldGrid(lonIndexGrid,latIndexGrid)/ &
                                                 real(icount(lonIndexGrid,latIndexGrid))
        end if
      end do
    end do

    ! Now do something for grid points that don't have any value assigned
    nfill = 0
    nhole = 0
    nextrap0 = 0
    nextrap1 = 0
    do latIndexGrid = 1, njGrid
      do lonIndexGrid = 1, niGrid
        if (icount(lonIndexGrid,latIndexGrid) == 0) then
          nhole = nhole + 1

          top    = latIndexGrid + 1
          bottom = latIndexGrid - 1
          left   = lonIndexGrid - 1
          right  = lonIndexGrid + 1

          np = 0
          l = bottom
          do k = left, right
            np = np + 1
            in(np) = k
            jn(np) = l
          end do
          k = right
          do l = bottom + 1, top
            np = np + 1
            in(np) = k
            jn(np) = l
          end do
          l = top
          do k = right - 1, left, -1
            np = np + 1
            in(np) = k
            jn(np) = l
          end do
          k = left
          do l = top - 1, bottom + 1, -1
            np = np + 1
            in(np) = k
            jn(np) = l
          end do

          ngp = 0
          do m = 1, np

            if (in(m) >= 1 .and. in(m) <= niGrid .and.    &
                jn(m) >= 1 .and. jn(m) <= njGrid) then
              if ( icount(in(m),jn(m)) > 0 ) then

                fieldGrid(lonIndexGrid,latIndexGrid) = fieldGrid(lonIndexGrid,latIndexGrid) +  &
                                                       fieldGrid(in(m),jn(m))
                ngp = ngp + 1

              end if
            end if

          end do

          if (ngp /= 0) then
            fieldGrid(lonIndexGrid,latIndexGrid) = fieldGrid(lonIndexGrid,latIndexGrid)/real(ngp)
            nfill = nfill + 1
          else
            fieldGrid(lonIndexGrid,latIndexGrid) = extrapValue
            if (maskGrid(lonIndexGrid,latIndexGrid) == 1) then
              nextrap1 = nextrap1 + 1
            else if (maskGrid(lonIndexGrid,latIndexGrid) == 0) then
              nextrap0 = nextrap0 + 1
            else
              write(*,*) 'Expecting 0 or 1 for the mask field at grid point: ',lonIndexGrid,latIndexGrid
              call utl_abort('int_sintCloudToGrid_gsv')
            end if
          end if

        end if
      end do
    end do

    write(*,*) 'Total number of grid points: ', niGrid*njGrid
    write(*,*) 'Number of grid points that were not covered '// &
               'by the cloud of points: ', nhole
    write(*,*) 'Number of grid points filled by a neighbor: ', nfill
    write(*,*) 'Number of grid points with extrapolated value in masked area: ', nextrap0
    write(*,*) 'Number of grid points with extrapolated value in visible area: ', nextrap1

    deallocate(xCloud)
    deallocate(yCloud)
    deallocate(icount)
    deallocate(maskCloud)
    deallocate(maskGrid)

    ierr = 0

  end function int_sintCloudToGrid_gsv


  function int_sint_r8_2d(zout8, zin8, hco_out, hco_in, interpDegree, extrapDegree_opt) result(ierr)
    implicit none

    ! arguments
    real(8) :: zout8(:,:), zin8(:,:)
    type(struct_hco), pointer :: hco_out, hco_in
    integer :: ierr
    character(len=*)           :: interpDegree
    character(len=*), optional :: extrapDegree_opt

    ! locals
    integer :: nii, nji, nio, njo     
    integer :: jk1, jk2
    real(4), allocatable :: bufferi4(:,:), buffero4(:,:)
    integer :: ezsint, ezdefset

    ! check if special interpolation is required
    if (hco_in%initialized .and. hco_out%initialized) then
      ierr = ezdefset(hco_out%EZscintID, hco_in%EZscintID)
      if (hco_in%grtyp == 'Y') then
        call utl_abort('int_sint_r8_2d: cloudToGrid not implemented')
        return
      end if
    end if

    ! do the standard interpolation

    call int_setezopt(interpDegree, extrapDegree_opt)   

    nii = size(zin8,1)
    nji = size(zin8,2)

    nio = size(zout8,1)
    njo = size(zout8,2)

    allocate(bufferi4(nii,nji))
    allocate(buffero4(nio,njo))

    do jk2 = 1,nji
      do jk1 = 1,nii
        bufferi4(jk1,jk2) = zin8(jk1,jk2)
      end do
    end do

    ierr = ezsint(buffero4,bufferi4)

    do jk2 = 1,njo
      do jk1 = 1,nio
        zout8(jk1,jk2) = buffero4(jk1,jk2)
      end do
    end do

    deallocate(bufferi4)
    deallocate(buffero4)

  end function int_sint_r8_2d


  function int_uvint_gsv(stateVectorOut, stateVectorIn, varName, levIndex, stepIndex, &
                        interpDegree, extrapDegree_opt) result(ierr)
    implicit none

    ! arguments
    type(struct_gsv) :: stateVectorOut
    type(struct_gsv) :: stateVectorIn
    character(len=*) :: varName
    integer          :: levIndex
    integer          :: stepIndex
    character(len=*)           :: interpDegree
    character(len=*), optional :: extrapDegree_opt
    integer :: ierr

    ! locals
    real(4), pointer :: UUout4(:,:,:,:), VVout4(:,:,:,:), UUin4(:,:,:,:), VVin4(:,:,:,:)
    real(8), pointer :: UUout8(:,:,:,:), VVout8(:,:,:,:), UUin8(:,:,:,:), VVin8(:,:,:,:)
    real(4), pointer :: UVout4(:,:,:), UVin4(:,:,:)
    real(8), pointer :: UVout8(:,:,:), UVin8(:,:,:)
    integer :: ezuvint, ezdefset

    write(*,*) 'int_uvint_gsv: starting'

    ! check if special interpolation is required
    if (stateVectorIn%hco%initialized .and. stateVectorOut%hco%initialized) then
      if (stateVectorIn%hco%grtyp == 'Y') then
        call utl_abort('int_uvint_gsv: cloudToGrid not implemented')
      end if
    end if

    ! do the standard interpolation

    ierr = ezdefset(stateVectorOut%hco%EZscintID, stateVectorIn%hco%EZscintID)
    call int_setezopt(interpDegree, extrapDegree_opt)   

    if (stateVectorOut%dataKind == 4 .and. stateVectorIn%dataKind == 4) then

      if (trim(varName) == 'BOTH') then
        call gsv_getField(stateVectorOut, UUout4, 'UU')
        call gsv_getField(stateVectorOut, VVout4, 'VV')
        call gsv_getField(stateVectorIn,  UUin4,  'UU')
        call gsv_getField(stateVectorIn,  VVin4,  'VV')
        ierr = ezuvint(UUout4(:,:,levIndex,stepIndex),VVout4(:,:,levIndex,stepIndex), &
                       UUin4(:,:,levIndex,stepIndex), VVin4(:,:,levIndex,stepIndex))
      else if (trim(varName) == 'UU') then
        call gsv_getField  (stateVectorIn,  UUin4)
        call gsv_getField  (stateVectorOut, UUout4)
        call gsv_getFieldUV(stateVectorIn,  UVin4,  levIndex)
        call gsv_getFieldUV(stateVectorOut, UVout4, levIndex)
        ierr = ezuvint(UUout4(:,:,levIndex,stepIndex),UVout4(:,:,stepIndex), &
                       UUin4(:,:,levIndex,stepIndex), UVin4(:,:,stepIndex))
      else if (trim(varName) == 'VV') then
        call gsv_getField  (stateVectorIn,  VVin4)
        call gsv_getField  (stateVectorOut, VVout4)
        call gsv_getFieldUV(stateVectorIn,  UVin4,  levIndex)
        call gsv_getFieldUV(stateVectorOut, UVout4, levIndex)
        ierr = ezuvint(UVout4(:,:,stepIndex),VVout4(:,:,levIndex,stepIndex), &
                       UVin4(:,:,stepIndex), VVin4(:,:,levIndex,stepIndex))
      else
        write(*,*) 'varName = ', trim(varName)
        call utl_abort('int_uvint_gsv: unexpected varName')
      end if

    else if (stateVectorOut%dataKind == 8 .and. stateVectorIn%dataKind == 8) then

      ! allocate real(4) buffers for copying to/from for interpolation
      allocate(UUin4(stateVectorIn%hco%ni,stateVectorIn%hco%nj,1,1))
      allocate(VVin4(stateVectorIn%hco%ni,stateVectorIn%hco%nj,1,1))
      allocate(UUout4(stateVectorOut%hco%ni,stateVectorOut%hco%nj,1,1))
      allocate(VVout4(stateVectorOut%hco%ni,stateVectorOut%hco%nj,1,1))

      if (trim(varName) == 'BOTH') then
        call gsv_getField(stateVectorOut, UUout8, 'UU')
        call gsv_getField(stateVectorOut, VVout8, 'VV')
        call gsv_getField(stateVectorIn,  UUin8,  'UU')
        call gsv_getField(stateVectorIn,  VVin8,  'VV')
        UUin4(:,:,1,1) = UUin8(:,:,levIndex,stepIndex)
        VVin4(:,:,1,1) = VVin8(:,:,levIndex,stepIndex)
        ierr = ezuvint(UUout4(:,:,1,1),VVout4(:,:,1,1),UUin4(:,:,1,1),VVin4(:,:,1,1))
        UUout8(:,:,levIndex,stepIndex) = UUout4(:,:,1,1)
        VVout8(:,:,levIndex,stepIndex) = VVout4(:,:,1,1)
      else if (trim(varName) == 'UU') then
        call gsv_getField  (stateVectorIn,  UUin8)
        call gsv_getField  (stateVectorOut, UUout8)
        call gsv_getFieldUV(stateVectorIn,  UVin8,  levIndex)
        call gsv_getFieldUV(stateVectorOut, UVout8, levIndex)
        UUin4(:,:,1,1) = UUin8(:,:,levIndex,stepIndex)
        VVin4(:,:,1,1) = UVin8(:,:,stepIndex)
        ierr = ezuvint(UUout4(:,:,1,1),VVout4(:,:,1,1),UUin4(:,:,1,1),VVin4(:,:,1,1))
        UUout8(:,:,levIndex,stepIndex) = UUout4(:,:,1,1)
        UVout8(:,:,stepIndex)          = VVout4(:,:,1,1)
      else if (trim(varName) == 'VV') then
        call gsv_getField  (stateVectorIn,  VVin8)
        call gsv_getField  (stateVectorOut, VVout8)
        call gsv_getFieldUV(stateVectorIn,  UVin8,  levIndex)
        call gsv_getFieldUV(stateVectorOut, UVout8, levIndex)
        UUin4(:,:,1,1) = UVin8(:,:,stepIndex)
        VVin4(:,:,1,1) = VVin8(:,:,levIndex,stepIndex)
        ierr = ezuvint(UUout4(:,:,1,1),VVout4(:,:,1,1),UUin4(:,:,1,1),VVin4(:,:,1,1))
        UVout8(:,:,stepIndex)          = UUout4(:,:,1,1)
        VVout8(:,:,levIndex,stepIndex) = VVout4(:,:,1,1)
      else
        write(*,*) 'varName = ', trim(varName)
        call utl_abort('int_uvint_gsv: unexpected varName')
      end if

      deallocate(UUin4,VVin4,UUout4,VVout4)

    else

      call utl_abort('int_uvint_gsv: not implemented for mixed dataKind')

    end if

  end function int_uvint_gsv


  function int_uvint_r4_2d(uuout, vvout, uuin, vvin, hco_out, hco_in, interpDegree, extrapDegree_opt) result(ierr)
    implicit none

    ! arguments
    real(4) :: uuout(:,:), vvout(:,:)
    real(4) :: uuin(:,:) , vvin(:,:)
    type(struct_hco), pointer :: hco_out, hco_in
    character(len=*)           :: interpDegree
    character(len=*), optional :: extrapDegree_opt
    integer :: ierr

    ! locals
    integer :: ezuvint, ezdefset

    ! check if special interpolation is required
    if (hco_in%initialized .and. hco_out%initialized) then
      ierr = ezdefset(hco_out%EZscintID, hco_in%EZscintID)
      if (hco_in%grtyp == 'Y') then
        call utl_abort('int_uvint_r4_2d: cloudToGrid not implemented')
        return
      end if
    end if

    ! do the standard interpolation
    call int_setezopt(interpDegree, extrapDegree_opt)   
    ierr = ezuvint(uuout, vvout, uuin, vvin)

  end function int_uvint_r4_2d


  function int_uvint_r8_2d(uuout, vvout, uuin, vvin, hco_out, hco_in, interpDegree, extrapDegree_opt) result(ierr)
    implicit none

    ! arguments
    real(8) :: uuout(:,:), vvout(:,:)
    real(8) :: uuin(:,:) , vvin(:,:)
    type(struct_hco), pointer :: hco_out, hco_in
    character(len=*)           :: interpDegree
    character(len=*), optional :: extrapDegree_opt
    integer :: ierr

    ! locals
    integer :: nio, njo, nii, nji
    integer :: jk1, jk2
    real, allocatable :: bufuuout4(:,:), bufvvout4(:,:)
    real, allocatable :: bufuuin4(:,:), bufvvin4(:,:)
    integer :: ezuvint, ezdefset

    ! check if special interpolation is required
    if (hco_in%initialized .and. hco_out%initialized) then
      ierr = ezdefset(hco_out%EZscintID, hco_in%EZscintID)
      if (hco_in%grtyp == 'Y') then
        call utl_abort('int_uvint_r8_2d: cloudToGrid not implemented')
        return
      end if
    end if

    ! do the standard interpolation

    call int_setezopt(interpDegree, extrapDegree_opt)   

    nii = size(uuin,1)
    nji = size(uuin,2)

    nio = size(uuout,1)
    njo = size(uuout,2)

    allocate(bufuuout4(nio,njo))
    allocate(bufvvout4(nio,njo))
    allocate(bufuuin4(nii,nji))
    allocate(bufvvin4(nii,nji))

    do jk2 = 1,nji
      do jk1 = 1,nii
        bufuuin4(jk1,jk2) = uuin(jk1,jk2)
        bufvvin4(jk1,jk2) = vvin(jk1,jk2)
      end do
    end do

    ierr = ezuvint(bufuuout4, bufvvout4, bufuuin4, bufvvin4)

    do jk2 = 1,njo
      do jk1 = 1,nio
        uuout(jk1,jk2) = bufuuout4(jk1,jk2)
        vvout(jk1,jk2) = bufvvout4(jk1,jk2)
      end do
    end do

    deallocate(bufuuin4)
    deallocate(bufvvin4)
    deallocate(bufuuout4)
    deallocate(bufvvout4)

  end function int_uvint_r8_2d


  function int_ezgdef(ni, nj, grtyp, grtypref, ig1, ig2, ig3, ig4, ax, ay) result(vezgdef)
    implicit none

    integer :: vezgdef

    integer :: ni, nj, ig1, ig2, ig3, ig4
    real(8) :: ax(*), ay(*)
    character(len=*) :: grtyp, grtypref

    integer :: ier2,jk,ilenx,ileny
    real, allocatable :: bufax4(:), bufay4(:)

    integer :: ezgdef

    if      (grtyp .eq. 'Y') then
       ilenx=max(1,ni*nj)
       ileny=ilenx
    else if (grtyp .eq. 'Z') then
       ilenx=max(1,ni)
       ileny=max(1,nj)
    else
       call utl_abort('VEZGDEF: Grid type not supported')
    end if

    allocate(bufax4(ilenx))
    allocate(bufay4(ileny))

    do jk = 1,ilenx
       bufax4(jk) = ax(jk)
    end do
    do jk = 1,ileny
       bufay4(jk) = ay(jk)
    end do

    ier2 = ezgdef(ni, nj, grtyp, grtypref, ig1, ig2, ig3, ig4, &
         bufax4, bufay4)

    deallocate(bufax4)
    deallocate(bufay4)

    vezgdef=ier2

  end function int_ezgdef


  subroutine int_cxgaig(grtyp, ig1, ig2, ig3, ig4, xlat0, xlon0, dlat, dlon) 
    implicit none

    integer :: ig1, ig2, ig3, ig4   
    real(8) :: xlat0, xlon0, dlat, dlon 
    character(len=*) :: grtyp 

    real(4) :: xlat04, xlon04, dlat4, dlon4

    xlat04=xlat0
    xlon04=xlon0
    dlat4=dlat
    dlon4=dlon

    call cxgaig(grtyp, ig1, ig2, ig3, ig4, xlat04, xlon04, dlat4, dlon4)

  end subroutine int_cxgaig

end module interpolation_mod

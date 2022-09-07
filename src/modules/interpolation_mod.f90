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
  use mpi_mod
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

    real(8), pointer :: field_in_r8_ptr(:,:,:,:), field_out_r8_ptr(:,:,:,:)
    real(8), pointer :: fieldUU_in_r8_ptr(:,:,:,:), fieldUU_out_r8_ptr(:,:,:,:)
    real(8), pointer :: fieldVV_in_r8_ptr(:,:,:,:), fieldVV_out_r8_ptr(:,:,:,:)
    real(8), pointer :: fieldUV_in_r8_ptr(:,:,:), fieldUV_out_r8_ptr(:,:,:)

    character(len=4) :: varName
    character(len=12):: interpolationDegree, extrapolationDegree

    integer :: ezdefset

    real(8), pointer :: heightSfcIn(:,:), heightSfcOut(:,:)

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
          ierr = ezdefset(statevector_out%hco%EZscintID, statevector_in%hco%EZscintID)

          if ( trim(varName) == 'UU' ) then
            ! interpolate both UV components and keep both UU and VV
            call gsv_getField(statevector_in,fieldUU_in_r8_ptr,'UU')
            call gsv_getField(statevector_out,fieldUU_out_r8_ptr,'UU')
            call gsv_getField(statevector_in,fieldVV_in_r8_ptr,'VV')
            call gsv_getField(statevector_out,fieldVV_out_r8_ptr,'VV')
            do levIndex = 1, nlev
              ierr = utl_ezuvint( fieldUU_out_r8_ptr(:,:,levIndex,stepIndex), fieldVV_out_r8_ptr(:,:,levIndex,stepIndex), &
                                  fieldUU_in_r8_ptr(:,:,levIndex,stepIndex),  fieldVV_in_r8_ptr(:,:,levIndex,stepIndex),  & 
                                  interpDegree=trim(interpolationDegree), extrapDegree_opt=trim(extrapolationDegree) )
            end do
          else
            ! interpolate scalar variable
            call gsv_getField(statevector_in, field_in_r8_ptr, varName)
            call gsv_getField(statevector_out, field_out_r8_ptr, varName)
            do levIndex = 1, nlev
              ierr = utl_ezsint( field_out_r8_ptr(:,:,levIndex,stepIndex), field_in_r8_ptr(:,:,levIndex,stepIndex),  &
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
          ierr = ezdefset(statevector_out%hco%EZscintID, statevector_in%hco%EZscintID)

          if ( trim(varName) == 'UU' ) then
            ! interpolate both UV components and keep UU in main vector
            call gsv_getField(statevector_in,fieldUU_in_r8_ptr)
            call gsv_getField(statevector_out,fieldUU_out_r8_ptr)
            call gsv_getFieldUV(statevector_in,fieldUV_in_r8_ptr,kIndex)
            call gsv_getFieldUV(statevector_out,fieldUV_out_r8_ptr,kIndex)
            ierr = utl_ezuvint( fieldUU_out_r8_ptr(:,:,kIndex,stepIndex), fieldUV_out_r8_ptr(:,:,stepIndex), &
                                fieldUU_in_r8_ptr(:,:,kIndex,stepIndex),  fieldUV_in_r8_ptr(:,:,stepIndex), &
                                interpDegree=trim(interpolationDegree), extrapDegree_opt=trim(extrapolationDegree) ) 
          else if ( trim(varName) == 'VV' ) then
            ! interpolate both UV components and keep VV in main vector
            call gsv_getFieldUV(statevector_in,fieldUV_in_r8_ptr,kIndex)
            call gsv_getFieldUV(statevector_out,fieldUV_out_r8_ptr,kIndex)
            call gsv_getField(statevector_in,fieldVV_in_r8_ptr)
            call gsv_getField(statevector_out,fieldVV_out_r8_ptr)
            ierr = utl_ezuvint( fieldUV_out_r8_ptr(:,:,stepIndex), fieldVV_out_r8_ptr(:,:,kIndex,stepIndex), &
                                fieldUV_in_r8_ptr(:,:,stepIndex),  fieldVV_in_r8_ptr(:,:,kIndex,stepIndex), &
                                interpDegree=trim(interpolationDegree), extrapDegree_opt=trim(extrapolationDegree) ) 
          else
            ! interpolate scalar variable
            call gsv_getField(statevector_in,field_in_r8_ptr)
            call gsv_getField(statevector_out,field_out_r8_ptr)
            ierr = utl_ezsint( field_out_r8_ptr(:,:,kIndex,stepIndex), field_in_r8_ptr(:,:,kIndex,stepIndex), &
                               interpDegree=trim(interpolationDegree), extrapDegree_opt=trim(extrapolationDegree) )
          end if
        end do k_loop

      end do ! stepIndex

    end if

    if ( gsv_isAssocHeightSfc(statevector_in) .and. gsv_isAssocHeightSfc(statevector_out) ) then
      heightSfcIn => gsv_getHeightSfc(statevector_in)
      heightSfcOut => gsv_getHeightSfc(statevector_out)
      write(*,*) 'int_hInterp_gsv: interpolating surface height'
      ierr = ezdefset(statevector_out%hco%EZscintID, statevector_in%hco%EZscintID)
      ierr = utl_ezsint( heightSfcOut(:,:), heightSfcIn(:,:), &
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

    real(4), pointer :: field_in_r4_ptr(:,:,:,:), field_out_r4_ptr(:,:,:,:)
    real(4), pointer :: fieldUU_in_r4_ptr(:,:,:,:), fieldUU_out_r4_ptr(:,:,:,:)
    real(4), pointer :: fieldVV_in_r4_ptr(:,:,:,:), fieldVV_out_r4_ptr(:,:,:,:)
    real(4), pointer :: fieldUV_in_r4_ptr(:,:,:), fieldUV_out_r4_ptr(:,:,:)

    character(len=4) :: varName
    character(len=12):: interpolationDegree, extrapolationDegree
    integer :: ezdefset

    real(8), pointer :: heightSfcIn(:,:), heightSfcOut(:,:)

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
          ierr = ezdefset(statevector_out%hco%EZscintID, statevector_in%hco%EZscintID)

          if ( trim(varName) == 'UU' ) then
            ! interpolate both UV components and keep both UU and VV
            call gsv_getField(statevector_in,fieldUU_in_r4_ptr,'UU')
            call gsv_getField(statevector_out,fieldUU_out_r4_ptr,'UU')
            call gsv_getField(statevector_in,fieldVV_in_r4_ptr,'VV')
            call gsv_getField(statevector_out,fieldVV_out_r4_ptr,'VV')
            do levIndex = 1, nlev
              ierr = utl_ezuvint( fieldUU_out_r4_ptr(:,:,levIndex,stepIndex), fieldVV_out_r4_ptr(:,:,levIndex,stepIndex),   &
                                  fieldUU_in_r4_ptr(:,:,levIndex,stepIndex),  fieldVV_in_r4_ptr(:,:,levIndex,stepIndex),    &
                                  interpDegree=trim(interpolationDegree), extrapDegree_opt=trim(extrapolationDegree) )
            end do
          else
            ! interpolate scalar variable
            call gsv_getField(statevector_in, field_in_r4_ptr, varName)
            call gsv_getField(statevector_out, field_out_r4_ptr, varName)
            do levIndex = 1, nlev
              ierr = utl_ezsint( field_out_r4_ptr(:,:,levIndex,stepIndex), field_in_r4_ptr(:,:,levIndex,stepIndex),  &
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
          ierr = ezdefset(statevector_out%hco%EZscintID, statevector_in%hco%EZscintID)

          if ( trim(varName) == 'UU' ) then
            ! interpolate both UV components and keep UU
            call gsv_getField(statevector_in,fieldUU_in_r4_ptr)
            call gsv_getField(statevector_out,fieldUU_out_r4_ptr)
            call gsv_getFieldUV(statevector_in,fieldUV_in_r4_ptr,kIndex)
            call gsv_getFieldUV(statevector_out,fieldUV_out_r4_ptr,kIndex)
            ierr = utl_ezuvint( fieldUU_out_r4_ptr(:,:,kIndex,stepIndex), fieldUV_out_r4_ptr(:,:,stepIndex),   &
                                fieldUU_in_r4_ptr(:,:,kIndex,stepIndex),  fieldUV_in_r4_ptr(:,:,stepIndex), &
                                interpDegree=trim(InterpolationDegree), extrapDegree_opt=trim(extrapolationDegree) ) 
          else if ( trim(varName) == 'VV' ) then
            ! interpolate both UV components and keep VV
            call gsv_getFieldUV(statevector_in,fieldUV_in_r4_ptr,kIndex)
            call gsv_getFieldUV(statevector_out,fieldUV_out_r4_ptr,kIndex)
            call gsv_getField(statevector_in,fieldVV_in_r4_ptr)
            call gsv_getField(statevector_out,fieldVV_out_r4_ptr)
            ierr = utl_ezuvint( fieldUV_out_r4_ptr(:,:,stepIndex), fieldVV_out_r4_ptr(:,:,kIndex,stepIndex),   &
                                fieldUV_in_r4_ptr(:,:,stepIndex),  fieldVV_in_r4_ptr(:,:,kIndex,stepIndex),  &
                                interpDegree=trim(InterpolationDegree), extrapDegree_opt=trim(extrapolationDegree) ) 
          else
            ! interpolate scalar variable
            call gsv_getField(statevector_in,field_in_r4_ptr)
            call gsv_getField(statevector_out,field_out_r4_ptr)
            ierr = utl_ezsint( field_out_r4_ptr(:,:,kIndex,stepIndex), field_in_r4_ptr(:,:,kIndex,stepIndex),  &
                               interpDegree=trim(InterpolationDegree), extrapDegree_opt=trim(extrapolationDegree) )
          end if
        end do k_loop

      end do step2_loop

    end if

    if ( gsv_isAssocHeightSfc(statevector_in) .and. gsv_isAssocHeightSfc(statevector_out)) then
      heightSfcIn => gsv_getHeightSfc(statevector_in)
      heightSfcOut => gsv_getHeightSfc(statevector_out)
      write(*,*) 'int_hInterp_gsv_r4: interpolating surface height'
      ierr = ezdefset(statevector_out%hco%EZscintID, statevector_in%hco%EZscintID)
      ierr = utl_ezsint( heightSfcOut(:,:), heightSfcIn(:,:), &
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
      if ( mpi_myid == 0 ) then
        write(*,*) 'int_vInterp_gsv: namint is missing in the namelist.'
        write(*,*) '                     The default values will be taken.'
      end if
    else
      ! Read namelist NAMINT
      nulnam=0
      ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam,nml=namint,iostat=ierr)
      if (ierr.ne.0) call utl_abort('int_vInterp_gsv: Error reading namelist NAMINT')
      if (mpi_myid.eq.0) write(*,nml=namint)
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
      if ( mpi_myid == 0 ) then
        write(*,*) 'int_vInterp_gsv_r4: namint is missing in the namelist.'
        write(*,*) '                     The default values will be taken.'
      end if
    else
      ! Read namelist NAMINT
      nulnam=0
      ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam,nml=namint,iostat=ierr)
      if (ierr.ne.0) call utl_abort('int_vInterp_gsv_r4: Error reading namelist NAMINT')
      if (mpi_myid.eq.0) write(*,nml=namint)
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
    if ( mpi_myid == 0 ) then
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

      if ( mpi_myid == 0 ) then
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

        if ( mpi_myid == 0 .and. columnIndex == 1 .and. &
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


end module interpolation_mod

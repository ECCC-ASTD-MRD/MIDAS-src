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
  use midasMpi_mod
  use gridstatevector_mod
  use columnData_mod
  use calcHeightAndPressure_mod
  use varNameList_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use mathPhysConstants_mod
  use utilities_mod
  use kdtree2_mod
  implicit none
  save
  private

  ! public subroutines and functions
  public :: int_interp_gsv
  public :: int_hInterp_gsv
  public :: int_vInterp_gsv, int_vInterp_gsv_r4
  public :: int_tInterp_gsv
  public :: int_vInterp_col
  public :: int_hInterpScalar, int_ezgdef, int_cxgaig

  ! module interfaces
  ! -----------------

  interface int_hInterpScalar
    module procedure int_hInterpScalar_gsv
    module procedure int_hInterpScalar_r4_2d
    module procedure int_hInterpScalar_r8_2d
  end interface int_hInterpScalar

  interface int_hInterpUV
    module procedure int_hInterpUV_gsv
    module procedure int_hInterpUV_r4_2d
    module procedure int_hInterpUV_r8_2d
  end interface int_hInterpUV

  ! Namelist variables
  ! ------------------
  logical :: vInterpCopyLowestLevel ! overwrite values at lowest level to avoid extrapolation
  logical :: checkCloudToGridUnassigned ! abort if unmasked points not assigned from cloudToGrid interp
  integer :: maxBoxSize             ! max size used to fill values for cloudToGrid interpolation

contains

  !--------------------------------------------------------------------------
  ! int_readNml (private subroutine)
  !--------------------------------------------------------------------------
  subroutine int_readNml()
    !
    ! :Purpose: Read the namelist block NAMINT.
    !
    ! :Namelist parameters:
    !         :vInterpCopyLowestLevel:  if true, will overwrite values at the lowest
    !                                   levels to avoid extrapolation
    !         :maxBoxSize: maximum size in terms of grid points used for filling
    !                      in undefined values from neighbours when doing interpolation
    !                      from a cloud of points to a grid
    !         :checkCloudToGridUnassigned: abort if any unmasked points are not assigned
    !                                      after doing cloudToGrid interpolation
    !
    implicit none

    ! Locals:
    integer :: ierr, nulnam, fnom, fclos
    logical, save :: alreadyRead = .false.

    NAMELIST /NAMINT/ vInterpCopyLowestLevel, checkCloudToGridUnassigned, maxBoxSize

    if (alreadyRead) then
      return
    else
      alreadyRead = .true.
    end if

    ! default values
    vInterpCopyLowestLevel = .false.
    checkCloudToGridUnassigned = .true.
    maxBoxSize = 1

    if ( .not. utl_isNamelistPresent('NAMINT','./flnml') ) then
      if (mmpi_myid == 0) then
        write(*,*) 'int_readNml: namint is missing in the namelist.'
        write(*,*) '             The default values will be taken.'
      end if
    else
      ! Read namelist NAMINT
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam,nml=namint,iostat=ierr)
      if (ierr /= 0) call utl_abort('int_readlNml: Error reading namelist NAMINT')
      if (mmpi_myid == 0) write(*,nml=namint)
      ierr = fclos(nulnam)
    end if

  end subroutine int_readNml

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

    real(8),      optional, intent(in)    :: PsfcReference_opt(:,:,:,:)     ! Provides a surface pressure field to be used instead of the first P0 level
    real(4),      optional, intent(in)    :: PsfcReference_r4_opt(:,:,:,:)  ! Provides a surface pressure field to be used instead of the first P0 level
    
    logical,      optional, intent(in)    :: checkModelTop_opt            ! If true, model top consistency will be checked in vertical interpolation
    
    ! Locals:
    logical :: checkModelTop
    character(len=4), pointer :: varNamesToInterpolate(:)
    type(struct_gsv) :: statevector_in_varsLevs, statevector_in_varsLevs_hInterp
    type(struct_gsv) :: statevector_in_hInterp

    write(*,*) 'int_interpolate: START'

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

    call int_hInterp_gsv(statevector_in_VarsLevs, statevector_in_VarsLevs_hInterp)
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
      call int_vInterp_gsv_r4(statevector_in_hInterp,statevector_out,&
                              PsfcReference_opt=PsfcReference_r4_opt, &
                              checkModelTop_opt=checkModelTop)
    else 
      call int_vInterp_gsv( statevector_in_hInterp,statevector_out,&
                            PsfcReference_opt=PsfcReference_opt, &
                            checkModelTop_opt=checkModelTop)
    end if

    call gsv_deallocate(statevector_in_hInterp)
    nullify(varNamesToInterpolate)

    write(*,*) 'int_interpolate: END'
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
    type(struct_gsv),       intent(inout) :: statevector_in   ! statevector that will contain the horizontally interpolated fields
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
              ierr = int_hInterpUV( statevector_out, statevector_in, 'BOTH', levIndex, stepIndex, &
                                    interpDegree=trim(interpolationDegree), &
                                    extrapDegree_opt=trim(extrapolationDegree) )
            end do
          else
            ! interpolate scalar variable
            do levIndex = 1, nlev
              ierr = int_hInterpScalar( statevector_out, statevector_in, varName, levIndex, stepIndex, &
                                        interpDegree=trim(interpolationDegree), &
                                        extrapDegree_opt=trim(extrapolationDegree) )
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
            ierr = int_hInterpUV( statevector_out, statevector_in, 'UU', kIndex, stepIndex, &
                                  interpDegree=trim(interpolationDegree), &
                                  extrapDegree_opt=trim(extrapolationDegree) )
          else if ( trim(varName) == 'VV' ) then
            ! interpolate both UV components and keep VV in main vector
            ierr = int_hInterpUV( statevector_out, statevector_in, 'VV', kIndex, stepIndex, &
                                  interpDegree=trim(interpolationDegree), &
                                  extrapDegree_opt=trim(extrapolationDegree) )
          else
            ! interpolate scalar variable
            ierr = int_hInterpScalar( statevector_out, statevector_in, 'ALL', kIndex, stepIndex, &
                                      interpDegree=trim(interpolationDegree), &
                                      extrapDegree_opt=trim(extrapolationDegree) )
          end if
        end do k_loop

      end do ! stepIndex

    end if

    if ( gsv_isAssocHeightSfc(statevector_in) .and. gsv_isAssocHeightSfc(statevector_out) ) then
      write(*,*) 'int_hInterp_gsv: interpolating surface height'
      ierr = int_hInterpScalar( statevector_out, statevector_in, 'ZSFC', 1, 1, &
                                interpDegree=trim(interpolationDegree), &
                                extrapDegree_opt=trim(extrapolationDegree) )
    end if

  end subroutine int_hInterp_gsv

  !--------------------------------------------------------------------------
  ! int_vInterp_gsv
  !--------------------------------------------------------------------------
  subroutine int_vInterp_gsv(statevector_in,statevector_out,Ps_in_hPa_opt, &
                              PsfcReference_opt,checkModelTop_opt)
    !
    ! :Purpose: Vertical interpolation
    !
    implicit none

    ! Arguments:
    type(struct_gsv),          intent(in)     :: statevector_in             ! statevector that will contain the vertically interpolated fields
    type(struct_gsv),          intent(inout)  :: statevector_out            ! Reference statevector providing the vertical structure
    logical, optional,         intent(in)     :: Ps_in_hPa_opt              ! If true, conversion from hPa to mbar will be done for surface pressure
    real(8), optional, target, intent(in)     :: PsfcReference_opt(:,:,:,:) ! Provides a surface pressure field to be used instead of the first P0 level
    logical, optional,         intent(in)     :: checkModelTop_opt          ! Model top consistency will be checked prior to interpolation if true

    ! Locals:
    logical :: checkModelTop

    integer :: vcode_in, vcode_out
    integer :: nlev_out, nlev_in
    integer :: varIndex, stepIndex
    integer :: status, nulnam, fnom, ierr, fclos

    real(8), pointer  :: hLikeT_in(:,:,:,:), hLikeM_in(:,:,:,:)   ! abstract height dimensioned coordinate
    real(8), pointer  :: hLikeT_out(:,:,:,:), hLikeM_out(:,:,:,:) ! abstract height dimensioned coordinate
    real(8), pointer  :: field_in(:,:,:,:), field_out(:,:,:,:)
    real(8), pointer  :: PsfcRef(:,:,:,:)
    real(8), pointer  :: heightSfcIn(:,:), heightSfcOut(:,:)

    character(len=4) :: varName

    type(struct_vco), pointer :: vco_in, vco_out

    ! read the namelist
    call int_readNml()

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

    vco_in  => gsv_getVco(statevector_in)
    vco_out => gsv_getVco(statevector_out)
    vcode_in  = vco_in%vcode
    vcode_out = vco_out%vcode

    ! the default is to ensure that the top of the output grid is ~equal or lower than the top of the input grid 
    if ( present(checkModelTop_opt) ) then
      checkModelTop = checkModelTop_opt
    else
      checkModelTop = .true.
    end if
    if (checkModelTop) then
      call vco_ensureCompatibleTops(vco_in, vco_out)
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
        field_out(:,:,:,:) = field_in(:,:,:,:)
        cycle var_loop
      end if

      ! allocating for temporary pressure references
      allocate(hLikeT_in( statevector_in%myLonBeg:statevector_in%myLonEnd, &
                          statevector_in%myLatBeg:statevector_in%myLatEnd, &
                          gsv_getNumLev(statevector_in,'TH'), statevector_in%numStep))
      allocate(hLikeM_in( statevector_in%myLonBeg:statevector_in%myLonEnd, &
                          statevector_in%myLatBeg:statevector_in%myLatEnd, &
                          gsv_getNumLev(statevector_in,'MM'), statevector_in%numStep))
      allocate(hLikeT_out(statevector_out%myLonBeg:statevector_out%myLonEnd, &
                          statevector_out%myLatBeg:statevector_out%myLatEnd, &
                          gsv_getNumLev(statevector_out,'TH'), statevector_out%numStep))
      allocate(hLikeM_out(statevector_out%myLonBeg:statevector_out%myLonEnd, &
                          statevector_out%myLatBeg:statevector_out%myLatEnd, &
                          gsv_getNumLev(statevector_out,'MM'), statevector_out%numStep))

      ! define pressures on input and output levels
      if ( present(PsfcReference_opt) ) then
        PsfcRef => PsfcReference_opt
      else
        call gsv_getField(statevector_in,PsfcRef,'P0')
      end if

      ! call czp to compute height or pressure
      ! convert if necessary to height-like abstract coordinate for interpolation
      if ( vcode_in==5002 .or. vcode_in==5005 ) then
        call czp_calcReturnPressure_gsv_nl( statevector_in, &
                                            PT_r8=hLikeT_in, PM_r8=hLikeM_in, &
                                            PsfcRef_r8_opt=PsfcRef, &
                                            Ps_in_hPa_opt=Ps_in_hPa_opt)
        write(*,*) 'int_vInterp_gsv converting input coordinates to height-like vcode=', vcode_in
        call logP(hLikeM_in)
        call logP(hLikeT_in)

      else if ( vcode_in==21001 ) then
        call czp_calcReturnHeight_gsv_nl( statevector_in, &
                                          ZT_r8=hLikeT_in, ZM_r8=hLikeM_in)
      end if

      if ( vcode_out==5002 .or. vcode_out==5005 ) then
        call czp_calcReturnPressure_gsv_nl( statevector_out, &
                                            PT_r8=hLikeT_out, PM_r8=hLikeM_out, &
                                            PsfcRef_r8_opt=PsfcRef, &
                                            Ps_in_hPa_opt=Ps_in_hPa_opt)
        write(*,*) 'int_vInterp_gsv converting output coordinates to height-like vcode=', vcode_out
        call logP(hLikeM_out)
        call logP(hLikeT_out)
      else if ( vcode_out==21001 ) then
        call czp_calcReturnHeight_gsv_nl( statevector_out, &
                                          ZT_r8=hLikeT_out, ZM_r8=hLikeM_out)
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
  
        ! do the vertical interpolation
        field_out(:,:,:,stepIndex) = 0.0d0
        if (vnl_varLevelFromVarname(varName) == 'TH') then
          call hLike_interpolation(hLikeT_in, hLikeT_out)
        else
          call hLike_interpolation(hLikeM_in, hLikeM_out)
        end if


        ! overwrite values at the lowest levels to avoid extrapolation
        if (vInterpCopyLowestLevel) then
          field_out(:,:,nlev_out,stepIndex) = field_in(:,:,nlev_in,stepIndex)
        end if

      end do step_loop

      deallocate(hLikeT_in, hLikeM_in, hLikeT_out, hLikeM_out)
    end do var_loop

    contains

      subroutine hLike_interpolation(hLike_in, hLike_out)
        !
        ! :Purpose: Proceed to actual interpolation in H-logP representation
        !
        implicit none

        ! Arguments:
        real(8), pointer, intent(in)  :: hLike_in(:,:,:,:), hLike_out(:,:,:,:)   ! abstract height dimensioned coordinate

        ! Locals:
        integer :: latIndex, lonIndex, levIndex_out, levIndex_in
        real(8) :: zwb, zwt

        !$OMP PARALLEL DO PRIVATE(latIndex,lonIndex,levIndex_in,levIndex_out,zwb,zwt)
        do latIndex = statevector_out%myLatBeg, statevector_out%myLatEnd
          do lonIndex = statevector_out%myLonBeg, statevector_out%myLonEnd
            levIndex_in = 1
            do levIndex_out = 1, nlev_out
              levIndex_in = levIndex_in + 1
              do while(hLike_out(lonIndex,latIndex,levIndex_out,stepIndex) &
                        .gt.hLike_in(lonIndex,latIndex,levIndex_in,stepIndex)  &
                       .and.levIndex_in.lt.nlev_in)
                levIndex_in = levIndex_in + 1
              end do
              levIndex_in = levIndex_in - 1
              zwb = (hLike_out(lonIndex,latIndex,levIndex_out,stepIndex) &
                      - hLike_in(lonIndex,latIndex,levIndex_in,stepIndex)) &
                   /(hLike_in(lonIndex,latIndex,levIndex_in+1,stepIndex) &
                      - hLike_in(lonIndex,latIndex,levIndex_in,stepIndex))
              zwt = 1.d0 - zwb
              field_out(lonIndex,latIndex,levIndex_out,stepIndex) =   &
                                 zwb*field_in(lonIndex,latIndex,levIndex_in+1,stepIndex) &
                               + zwt*field_in(lonIndex,latIndex,levIndex_in,stepIndex)
            end do
          end do
        end do
        !$OMP END PARALLEL DO

      end subroutine hLike_interpolation

  end subroutine int_vInterp_gsv

  !--------------------------------------------------------------------------
  ! logP
  !--------------------------------------------------------------------------
  subroutine logP(presInLogOut)
    !
    ! :Purpose: compute log of pressurce field
    !
    implicit none

    ! Arguments:
    real(8), intent(inout) :: presInLogOut(:,:,:,:)

    ! Locals:
    integer :: latIdx, lonIdx, levIdx, stepIdx

    !$OMP PARALLEL DO PRIVATE(lonIdx,latIdx,levIdx,stepIdx)
    do lonIdx = lbound(presInLogOut,1), ubound(presInLogOut,1)
      do latIdx = lbound(presInLogOut,2), ubound(presInLogOut,2)
        do levIdx = lbound(presInLogOut,3), ubound(presInLogOut,3)
          do stepIdx = lbound(presInLogOut,4), ubound(presInLogOut,4)
            presInLogOut(lonIdx,latIdx,levIdx,stepIdx) = &
                 log(presInLogOut(lonIdx,latIdx,levIdx,stepIdx))
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine logP

  !--------------------------------------------------------------------------
  ! int_vInterp_gsv_r4
  !--------------------------------------------------------------------------
  subroutine int_vInterp_gsv_r4(statevector_in,statevector_out,Ps_in_hPa_opt, &
                                 PsfcReference_opt,checkModelTop_opt)
    !
    ! :Purpose: Vertical interpolation
    !
    !
    implicit none

    ! Arguments:
    type(struct_gsv),          intent(in)     :: statevector_in             ! statevector that will contain the vertically interpolated fields
    type(struct_gsv),          intent(inout)  :: statevector_out            ! Reference statevector providing the vertical structure
    logical, optional,         intent(in)     :: Ps_in_hPa_opt              ! If true, conversion from hPa to mbar will be done for surface pressure
    real(4), optional, target, intent(in)     :: PsfcReference_opt(:,:,:,:) ! Provides a surface pressure field to be used instead of the first P0 level
    logical, optional,         intent(in)     :: checkModelTop_opt          ! Model top consistency will be checked prior to interpolation if true

    ! Locals:
    logical :: checkModelTop

    integer :: vcode_in, vcode_out
    integer :: nlev_out, nlev_in
    integer :: varIndex, stepIndex
    integer :: status, nulnam, fnom, ierr, fclos

    real(4), pointer  :: hLikeT_in(:,:,:,:), hLikeM_in(:,:,:,:)   ! abstract height dimensioned coordinate
    real(4), pointer  :: hLikeT_out(:,:,:,:), hLikeM_out(:,:,:,:) ! abstract height dimensioned coordinate
    real(4), pointer  :: field_in(:,:,:,:), field_out(:,:,:,:)
    real(4), pointer  :: PsfcRef(:,:,:,:)
    real(8), pointer  :: heightSfcIn(:,:), heightSfcOut(:,:)

    character(len=4) :: varName

    type(struct_vco), pointer :: vco_in, vco_out

    ! read the namelist
    call int_readNml()

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

    if (gsv_isAssocHeightSfc(statevector_in) .and. gsv_isAssocHeightSfc(statevector_out) ) then
      heightSfcIn => gsv_getHeightSfc(statevector_in)
      heightSfcOut => gsv_getHeightSfc(statevector_out)
      heightSfcOut(:,:) = heightSfcIn(:,:)
    end if

    vco_in  => gsv_getVco(statevector_in)
    vco_out => gsv_getVco(statevector_out)
    vcode_in  = vco_in%vcode
    vcode_out = vco_out%vcode

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
        field_out(:,:,:,:) = field_in(:,:,:,:)
        cycle var_loop
      end if

      ! allocating for temporary pressure references
      allocate(hLikeT_in( statevector_in%myLonBeg:statevector_in%myLonEnd, &
                          statevector_in%myLatBeg:statevector_in%myLatEnd, &
                          gsv_getNumLev(statevector_in,'TH'), statevector_in%numStep))
      allocate(hLikeM_in( statevector_in%myLonBeg:statevector_in%myLonEnd, &
                          statevector_in%myLatBeg:statevector_in%myLatEnd, &
                          gsv_getNumLev(statevector_in,'MM'), statevector_in%numStep))
      allocate(hLikeT_out(statevector_out%myLonBeg:statevector_out%myLonEnd, &
                          statevector_out%myLatBeg:statevector_out%myLatEnd, &
                          gsv_getNumLev(statevector_out,'TH'), statevector_out%numStep))
      allocate(hLikeM_out(statevector_out%myLonBeg:statevector_out%myLonEnd, &
                          statevector_out%myLatBeg:statevector_out%myLatEnd, &
                          gsv_getNumLev(statevector_out,'MM'), statevector_out%numStep))

      ! define pressures on input and output levels
      if ( present(PsfcReference_opt) ) then
        PsfcRef => PsfcReference_opt
      else
        call gsv_getField(statevector_in,PsfcRef,'P0')
      end if

      ! call czp to compute height or pressure
      ! convert if necessary to height-like abstract coordinate for interpolation
      if ( vcode_in==5002 .or. vcode_in==5005 ) then
        call czp_calcReturnPressure_gsv_nl( statevector_in, &
                                            PT_r4=hLikeT_in, PM_r4=hLikeM_in, &
                                            PsfcRef_r4_opt=PsfcRef, &
                                            Ps_in_hPa_opt=Ps_in_hPa_opt)
        write(*,*) 'int_vInterp_gsv_r4 converting input coordinates to height-like vcode=', vcode_in
        call logP_r4(hLikeM_in)
        call logP_r4(hLikeT_in)

      else if ( vcode_in==21001 ) then
        call czp_calcReturnHeight_gsv_nl( statevector_in, &
                                          ZT_r4=hLikeT_in, ZM_r4=hLikeM_in)
      end if

      if ( vcode_out==5002 .or. vcode_out==5005 ) then
        call czp_calcReturnPressure_gsv_nl( statevector_out, &
                                            PT_r4=hLikeT_out, PM_r4=hLikeM_out, &
                                            PsfcRef_r4_opt=PsfcRef, &
                                            Ps_in_hPa_opt=Ps_in_hPa_opt)
        write(*,*) 'int_vInterp_gsv_r4 converting output coordinates to height-like vcode=', vcode_out
        call logP_r4(hLikeM_out)
        call logP_r4(hLikeT_out)
      else if ( vcode_out==21001 ) then
        call czp_calcReturnHeight_gsv_nl( statevector_out, &
                                          ZT_r4=hLikeT_out, ZM_r4=hLikeM_out)
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
  
        ! do the vertical interpolation
        field_out(:,:,:,stepIndex) = 0.0d0
        if (vnl_varLevelFromVarname(varName) == 'TH') then
          call hLike_interpolation_r4(hLikeT_in, hLikeT_out)
        else
          call hLike_interpolation_r4(hLikeM_in, hLikeM_out)
        end if


        ! overwrite values at the lowest levels to avoid extrapolation
        if (vInterpCopyLowestLevel) then
          field_out(:,:,nlev_out,stepIndex) = field_in(:,:,nlev_in,stepIndex)
        end if

      end do step_loop

      deallocate(hLikeT_in, hLikeM_in, hLikeT_out, hLikeM_out)
    end do var_loop

    contains

      subroutine hLike_interpolation_r4(hLike_in, hLike_out)
        !
        ! :Purpose: Proceed to actual interpolation in H-logP representation
        !
        implicit none

        ! Arguments:
        real(4), pointer, intent(in)  :: hLike_in(:,:,:,:), hLike_out(:,:,:,:)   ! abstract height dimensioned coordinate

        ! Locals:
        integer :: latIndex, lonIndex, levIndex_out, levIndex_in
        real(4) :: zwb, zwt

        !$OMP PARALLEL DO PRIVATE(latIndex,lonIndex,levIndex_in,levIndex_out,zwb,zwt)
        do latIndex = statevector_out%myLatBeg, statevector_out%myLatEnd
          do lonIndex = statevector_out%myLonBeg, statevector_out%myLonEnd
            levIndex_in = 1
            do levIndex_out = 1, nlev_out
              levIndex_in = levIndex_in + 1
              do while(hLike_out(lonIndex,latIndex,levIndex_out,stepIndex) &
                        .gt.hLike_in(lonIndex,latIndex,levIndex_in,stepIndex)  &
                       .and.levIndex_in.lt.nlev_in)
                levIndex_in = levIndex_in + 1
              end do
              levIndex_in = levIndex_in - 1
              zwb = (hLike_out(lonIndex,latIndex,levIndex_out,stepIndex) &
                      - hLike_in(lonIndex,latIndex,levIndex_in,stepIndex)) &
                   /(hLike_in(lonIndex,latIndex,levIndex_in+1,stepIndex) &
                      - hLike_in(lonIndex,latIndex,levIndex_in,stepIndex))
              zwt = 1.d0 - zwb
              field_out(lonIndex,latIndex,levIndex_out,stepIndex) =   &
                                 zwb*field_in(lonIndex,latIndex,levIndex_in+1,stepIndex) &
                               + zwt*field_in(lonIndex,latIndex,levIndex_in,stepIndex)
            end do
          end do
        end do
        !$OMP END PARALLEL DO

      end subroutine hLike_interpolation_r4

  end subroutine int_vInterp_gsv_r4

  !--------------------------------------------------------------------------
  ! logP_r4
  !--------------------------------------------------------------------------
  subroutine logP_r4(presInLogOut)
    !
    ! :Purpose: compute log of pressurce field
    !
    implicit none

    ! Arguments:
    real(4), intent(inout) :: presInLogOut(:,:,:,:)

    ! Locals:
    integer :: latIdx, lonIdx, levIdx, stepIdx

    !$OMP PARALLEL DO PRIVATE(lonIdx,latIdx,levIdx,stepIdx)
    do lonIdx = lbound(presInLogOut,1), ubound(presInLogOut,1)
      do latIdx = lbound(presInLogOut,2), ubound(presInLogOut,2)
        do levIdx = lbound(presInLogOut,3), ubound(presInLogOut,3)
          do stepIdx = lbound(presInLogOut,4), ubound(presInLogOut,4)
            presInLogOut(lonIdx,latIdx,levIdx,stepIdx) = &
                 log(presInLogOut(lonIdx,latIdx,levIdx,stepIdx))
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine logP_r4

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

    ! read the namelist
    call int_readNml()

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

  !--------------------------------------------------------------------------
  ! int_setezopt
  !--------------------------------------------------------------------------
  subroutine int_setezopt(interpDegree, extrapDegree_opt)
    !
    ! :Purpose: Wrapper subroutine for rmnlib routine setezopt.
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in)           :: interpDegree
    character(len=*), intent(in), optional :: extrapDegree_opt

    ! Locals:
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

  !--------------------------------------------------------------------------
  ! int_hInterpScalar_gsv
  !--------------------------------------------------------------------------
  function int_hInterpScalar_gsv(stateVectorOut, stateVectorIn, varName, levIndex, stepIndex, &
                        interpDegree, extrapDegree_opt) result(ierr)
    !
    ! :Purpose: Horizontal interpolation of 2D scalar field that use stateVector
    !           objects for input and output. Accessed through int_hInterpScalar.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: stateVectorOut
    type(struct_gsv), intent(inout) :: stateVectorIn
    character(len=*), intent(in)    :: varName
    integer         , intent(in)    :: levIndex
    integer         , intent(in)    :: stepIndex
    character(len=*), intent(in)           :: interpDegree
    character(len=*), intent(in), optional :: extrapDegree_opt
    integer :: ierr

    ! Locals:
    real(4), pointer :: fieldOut_r4(:,:,:,:), fieldIn_r4(:,:,:,:)
    real(8), pointer :: fieldOut_r8(:,:,:,:), fieldIn_r8(:,:,:,:)
    real(8), pointer :: heightSfcOut(:,:), heightSfcIn(:,:)
    integer :: ezsint, ezdefset

    ! read the namelist
    call int_readNml()

    ! check if special interpolation is required
    if (stateVectorIn%hco%initialized .and. stateVectorOut%hco%initialized) then
      if (stateVectorIn%hco%grtyp == 'Y') then
        ! for now, only comptable for real(4)
        if(stateVectorOut%dataKind /= 4 .or. stateVectorIn%dataKind /= 4) then
          call utl_abort('int_hInterpScalar_gsv: cloudToGrid only implemented for real(4)')
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
      allocate(fieldIn_r4(stateVectorIn%hco%ni,stateVectorIn%hco%nj,1,1))
      allocate(fieldOut_r4(stateVectorOut%hco%ni,stateVectorOut%hco%nj,1,1))
      fieldIn_r4(:,:,1,1) = heightSfcIn(:,:)
      ierr = ezsint(fieldOut_r4(:,:,1,1),fieldIn_r4(:,:,1,1))
      heightSfcOut(:,:) = fieldOut_r4(:,:,1,1)
      deallocate(fieldIn_r4,fieldOut_r4)

    else if (stateVectorOut%dataKind == 4 .and. stateVectorIn%dataKind == 4) then

      if (trim(varName) == 'ALL') then
        call gsv_getField(stateVectorOut, fieldOut_r4)
        call gsv_getField(stateVectorIn,  fieldIn_r4)
     else
        call gsv_getField(stateVectorOut, fieldOut_r4, varName)
        call gsv_getField(stateVectorIn,  fieldIn_r4,  varName)
      end if

      ierr = ezsint(fieldOut_r4(:,:,levIndex,stepIndex),fieldIn_r4(:,:,levIndex,stepIndex))

    else if (stateVectorOut%dataKind == 8 .and. stateVectorIn%dataKind == 8) then

      if (trim(varName) == 'ALL') then
        call gsv_getField(stateVectorOut, fieldOut_r8)
        call gsv_getField(stateVectorIn,  fieldIn_r8)
      else
        call gsv_getField(stateVectorOut, fieldOut_r8, varName)
        call gsv_getField(stateVectorIn,  fieldIn_r8,  varName)
      end if

      ! allocate real(4) buffers and copy to/from for interpolation
      allocate(fieldIn_r4(stateVectorIn%hco%ni,stateVectorIn%hco%nj,1,1))
      allocate(fieldOut_r4(stateVectorOut%hco%ni,stateVectorOut%hco%nj,1,1))
      fieldIn_r4(:,:,1,1) = fieldIn_r8(:,:,levIndex,stepIndex)
      ierr = ezsint(fieldOut_r4(:,:,1,1),fieldIn_r4(:,:,1,1))
      fieldOut_r8(:,:,levIndex,stepIndex) = fieldOut_r4(:,:,1,1)
      deallocate(fieldIn_r4,fieldOut_r4)

    else

      call utl_abort('int_hInterpScalar_gsv: not implemented for mixed dataKind')

    end if

  end function int_hInterpScalar_gsv

  !--------------------------------------------------------------------------
  ! int_hInterpScalar_r4_2d
  !--------------------------------------------------------------------------
  function int_hInterpScalar_r4_2d(fieldOut_r4, fieldIn_r4, interpDegree, extrapDegree_opt) result(ierr)
    !
    ! :Purpose: Horizontal interpolation of 2D scalar field that use real(4) arrays
    !           for input and output. Accessed through int_hInterpScalar.
    !
    implicit none

    ! Arguments:
    real(4), intent(inout) :: fieldOut_r4(:,:)
    real(4), intent(in)    :: fieldIn_r4(:,:)
    integer :: ierr
    character(len=*)           :: interpDegree
    character(len=*), optional :: extrapDegree_opt

    ! Locals:
    integer :: ezsint

    ! read the namelist
    call int_readNml()

    ! do the standard interpolation
    call int_setezopt(interpDegree, extrapDegree_opt)   
    ierr = ezsint(fieldOut_r4,fieldIn_r4)

  end function int_hInterpScalar_r4_2d

  !--------------------------------------------------------------------------
  ! int_sintCloudToGrid_gsv
  !--------------------------------------------------------------------------
  function int_sintCloudToGrid_gsv(stateVectorGrid, stateVectorCloud, varName, levIndex, stepIndex) result(ierr)
    !
    ! :Purpose: Perform horizontal interpolation for 1 level and time step (and variable)
    !           in the case where the input data is a cloud of points (i.e. a Y grid) and
    !           the output is on a regular grid. Accessed through int_hInterpScalar.
    !
    ! :Note:  When varName=='ALL', the argument levIndex is actually kIndex
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: stateVectorGrid
    type(struct_gsv), intent(inout) :: stateVectorCloud
    character(len=*), intent(in)    :: varName
    integer,          intent(in)    :: levIndex
    integer,          intent(in)    :: stepIndex
    integer                         :: ierr

    ! Locals:
    integer :: gdxyfll, omp_get_thread_num
    integer :: niCloud, njCloud, niGrid, njGrid, myThreadNum
    integer :: top, bottom, left, right, numBoxIndexes, lonIndexCloud, latIndexCloud
    integer :: boxSize, lonBoxIndex, latBoxIndex, boxIndex, lonIndexGrid, latIndexGrid
    integer :: lonBoxIndexes(100), latBoxIndexes(100), ngp, nfill(mmpi_numThread), nhole(mmpi_numThread), nextrap0, nextrap1
    integer, allocatable :: numFilledByAvg(:,:), filledByInterp(:,:), maskGrid(:,:), maskCloud(:,:)
    real(4), pointer     :: fieldCloud_4d(:,:,:,:), fieldGrid_4d(:,:,:,:)
    real(4), pointer     :: fieldCloud(:,:), fieldGrid(:,:)
    real(4), allocatable :: fieldGrid_tmp(:,:)
    real(4), allocatable :: xCloud(:,:), yCloud(:,:)
    real(8) :: extrapValue
    integer, parameter        :: maxNumLocalGridPointsSearch = 200000
    integer                   :: numLocalGridPointsFound, gridIndex
    type(kdtree2), pointer    :: tree => null()
    real(kdkind), allocatable :: positionArray(:,:)
    type(kdtree2_result)      :: searchResults(maxNumLocalGridPointsSearch)
    real(kdkind)              :: searchRadiusSquared
    real(kdkind)              :: refPosition(3)

    call utl_tmg_start(176, 'low-level--int_sintCloudToGrid_gsv')

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
    allocate(numFilledByAvg(niGrid, njGrid))
    allocate(filledByInterp(niGrid, njGrid))
    allocate(maskCloud(niCloud, njCloud))
    allocate(maskGrid(niGrid, njGrid))
    allocate(fieldGrid_tmp(niGrid,njGrid))

    ! set masks based on oceanMasks (if present)
    if (stateVectorCloud%oceanMask%maskPresent) then
      maskCloud(:,:) = 0
      if (trim(varName) == 'ALL') then
        ! when varName==ALL, the argument levIndex is actually kIndex
        where(stateVectorCloud%oceanMask%mask(:,:,gsv_getLevFromK(stateVectorCloud,levIndex))) maskCloud(:,:) = 1
      else
        where(stateVectorCloud%oceanMask%mask(:,:,levIndex)) maskCloud(:,:) = 1
      end if
    else
      maskCloud(:,:) = 1
    end if
    if (stateVectorGrid%oceanMask%maskPresent) then
      maskGrid(:,:) = 0
      if (trim(varName) == 'ALL') then
        ! when varName==ALL, the argument levIndex is actually kIndex
        where(stateVectorGrid%oceanMask%mask(:,:,gsv_getLevFromK(stateVectorGrid,levIndex))) maskGrid(:,:) = 1
      else
        where(stateVectorGrid%oceanMask%mask(:,:,levIndex)) maskGrid(:,:) = 1
      end if
    else
      maskGrid(:,:)  = 1
    end if

    ! Calcul des pos. x-y des eclairs sur la grille modele
    ierr = gdxyfll(stateVectorGrid%hco%EZscintID, xCloud, yCloud, &
                   stateVectorCloud%hco%lat2d_4*MPC_DEGREES_PER_RADIAN_R4, &
                   stateVectorCloud%hco%lon2d_4*MPC_DEGREES_PER_RADIAN_R4, &
                   stateVectorCloud%hco%ni*stateVectorCloud%hco%nj)

    ! Average values of cloud points neighbouring a grid location
    fieldGrid(:,:) = 0.0
    numFilledByAvg(:,:) = 0
    do latIndexCloud = 1, njCloud
      do lonIndexCloud = 1, niCloud
        lonIndexGrid = nint(xCloud(lonIndexCloud,latIndexCloud))
        latIndexGrid = nint(yCloud(lonIndexCloud,latIndexCloud))
        if ( lonIndexGrid >= 1 .and. latIndexGrid >= 1 .and. &
             lonIndexGrid <= niGrid .and. latIndexGrid <= njGrid ) then
          if ( maskCloud(lonIndexCloud,latIndexCloud) == 1 ) then
            fieldGrid(lonIndexGrid,latIndexGrid) = fieldGrid(lonIndexGrid,latIndexGrid) +  &
                                                   fieldCloud(lonIndexCloud,latIndexCloud)
            numFilledByAvg(lonIndexGrid,latIndexGrid) = numFilledByAvg(lonIndexGrid,latIndexGrid) + 1
          end if
        end if
      end do
    end do

    do latIndexGrid = 1, njGrid
      do lonIndexGrid = 1, niGrid
        if(numFilledByAvg(lonIndexGrid,latIndexGrid) > 0) then
          fieldGrid(lonIndexGrid,latIndexGrid) = fieldGrid(lonIndexGrid,latIndexGrid)/ &
                                                 real(numFilledByAvg(lonIndexGrid,latIndexGrid))
        end if
      end do
    end do

    ! Now do something for grid points that don't have any value assigned
    fieldGrid_tmp(:,:) = fieldGrid(:,:)
    nfill(:) = 0
    nhole(:) = 0
    boxSizeLoop: do boxSize = 1, maxBoxSize
      filledByInterp(:,:) = 0
      !$OMP PARALLEL DO PRIVATE(latIndexGrid,lonIndexGrid,myThreadNum,top,bottom,left,right,numBoxIndexes,lonBoxIndex,latBoxIndex,lonBoxIndexes,latBoxIndexes,ngp,boxIndex)
      do latIndexGrid = 1, njGrid
        myThreadNum = 1 + omp_get_thread_num()
        do lonIndexGrid = 1, niGrid
          if (numFilledByAvg(lonIndexGrid,latIndexGrid) > 0) cycle

          if (boxSize == 1) nhole(myThreadNum) = nhole(myThreadNum) + 1

          top    = latIndexGrid + boxSize
          bottom = latIndexGrid - boxSize
          left   = lonIndexGrid - boxSize
          right  = lonIndexGrid + boxSize

          numBoxIndexes = 0
          latBoxIndex = bottom
          do lonBoxIndex = left, right
            numBoxIndexes = numBoxIndexes + 1
            lonBoxIndexes(numBoxIndexes) = lonBoxIndex
            latBoxIndexes(numBoxIndexes) = latBoxIndex
          end do
          lonBoxIndex = right
          do latBoxIndex = bottom + 1, top
            numBoxIndexes = numBoxIndexes + 1
            lonBoxIndexes(numBoxIndexes) = lonBoxIndex
            latBoxIndexes(numBoxIndexes) = latBoxIndex
          end do
          latBoxIndex = top
          do lonBoxIndex = right - 1, left, -1
            numBoxIndexes = numBoxIndexes + 1
            lonBoxIndexes(numBoxIndexes) = lonBoxIndex
            latBoxIndexes(numBoxIndexes) = latBoxIndex
          end do
          lonBoxIndex = left
          do latBoxIndex = top - 1, bottom + 1, -1
            numBoxIndexes = numBoxIndexes + 1
            lonBoxIndexes(numBoxIndexes) = lonBoxIndex
            latBoxIndexes(numBoxIndexes) = latBoxIndex
          end do

          ngp = 0
          do boxIndex = 1, numBoxIndexes

            if (lonBoxIndexes(boxIndex) >= 1 .and. lonBoxIndexes(boxIndex) <= niGrid .and.    &
                latBoxIndexes(boxIndex) >= 1 .and. latBoxIndexes(boxIndex) <= njGrid) then
              if ( numFilledByAvg(lonBoxIndexes(boxIndex),latBoxIndexes(boxIndex)) > 0 ) then

                fieldGrid_tmp(lonIndexGrid,latIndexGrid) = &
                     fieldGrid_tmp(lonIndexGrid,latIndexGrid) +  &
                     fieldGrid(lonBoxIndexes(boxIndex),latBoxIndexes(boxIndex))
                ngp = ngp + 1

              end if
            end if

          end do

          if (ngp /= 0) then
            ! mark the grid point as being filled by interpolation on the grid
            filledByInterp(lonIndexGrid,latIndexGrid) = 1
            fieldGrid_tmp(lonIndexGrid,latIndexGrid) = fieldGrid_tmp(lonIndexGrid,latIndexGrid)/real(ngp)
            nfill(myThreadNum) = nfill(myThreadNum) + 1
          end if

        end do
      end do
      !$OMP END PARALLEL DO

      fieldGrid(:,:) = fieldGrid_tmp(:,:)
      numFilledByAvg(:,:) = numFilledByAvg(:,:) + filledByInterp(:,:)
    end do boxSizeLoop

    ! find any remaining grid points that are likely water and assign value
    if ( checkCloudToGridUnassigned .and. &
         stateVectorCloud%oceanMask%maskPresent .and. &
         .not. stateVectorGrid%oceanMask%maskPresent ) then

      ! create a kdtree to index all cloud locations
      allocate(positionArray(3,niCloud*njCloud))
      gridIndex = 0
      do lonIndexCloud = 1, niCloud
        do latIndexCloud = 1, njCloud
          gridIndex = gridIndex + 1
          positionArray(:,gridIndex) = &
               kdtree2_3dPosition(real(stateVectorCloud%hco%lon2d_4(lonIndexCloud,latIndexCloud),8), &
                                  real(stateVectorCloud%hco%lat2d_4(lonIndexCloud,latIndexCloud),8))
        end do
      end do
      tree => kdtree2_create(positionArray, sort=.true., rearrange=.true.) 
      
      ! assign grid mask to value of mask at nearest cloud location
      do latIndexGrid = 1, njGrid
        do lonIndexGrid = 1, niGrid
          if (numFilledByAvg(lonIndexGrid,latIndexGrid) > 0) cycle

          ! find nearest cloud location
          refPosition(:) = kdtree2_3dPosition(real(stateVectorGrid%hco%lon2d_4(lonIndexGrid,latIndexGrid),8), &
                                              real(stateVectorGrid%hco%lat2d_4(lonIndexGrid,latIndexGrid),8))

          searchRadiusSquared = (50.0D0*1000.0D0)**2
          call kdtree2_r_nearest(tp=tree, qv=refPosition, r2=searchRadiusSquared, &
                                 nfound=numLocalGridPointsFound, &
                                 nalloc=maxNumLocalGridPointsSearch, results=searchResults)
          if (numLocalGridPointsFound > maxNumLocalGridPointsSearch) then
            call utl_abort('The parameter maxNumLocalGridPointsSearch must be increased')
          end if
          if (numLocalGridPointsFound == 0) then
            ! very far from a cloud location, assume it is land
            maskGrid(lonIndexGrid,latIndexGrid) = 0
          else
            gridIndex = searchResults(1)%idx
            lonIndexCloud = 1 + ((gridIndex-1)/njCloud)
            latIndexCloud = gridIndex - ((gridIndex-1)/njCloud)*njCloud
            if(lonIndexCloud < 1 .or. lonIndexCloud > niCloud) call utl_abort('lonIndexCloud wrong')
            if(latIndexCloud < 1 .or. latIndexCloud > njCloud) call utl_abort('latIndexCloud wrong')
            maskGrid(lonIndexGrid,latIndexGrid) = maskCloud(lonIndexCloud,latIndexCloud)
          end if

        end do
      end do

      deallocate(positionArray)
      call kdtree2_destroy(tree)

    end if

    ! fill in remaining grid points
    if (count(numFilledByAvg(:,:) > 0) > 0) then

      ! compute spatial mean of assigned grid values
      extrapValue = 0.0d0
      do latIndexGrid = 1, njGrid
        do lonIndexGrid = 1, niGrid
          if (numFilledByAvg(lonIndexGrid,latIndexGrid) == 0) cycle
          extrapValue = extrapValue + real(fieldGrid(lonIndexGrid,latIndexGrid),8)
        end do
      end do
      extrapValue = extrapValue / real(count(numFilledByAvg(:,:) > 0),8)

      ! set unassigned grid points and count them
      nextrap0 = 0
      nextrap1 = 0
      do latIndexGrid = 1, njGrid
        do lonIndexGrid = 1, niGrid

          if (numFilledByAvg(lonIndexGrid,latIndexGrid) > 0) cycle

          fieldGrid(lonIndexGrid,latIndexGrid) = extrapValue

          if (maskGrid(lonIndexGrid,latIndexGrid) == 1) then
            nextrap1 = nextrap1 + 1
          else if (maskGrid(lonIndexGrid,latIndexGrid) == 0) then
            nextrap0 = nextrap0 + 1
          else
            write(*,*) 'Expecting 0 or 1 for the mask field at grid point: ', &
                       lonIndexGrid,latIndexGrid,maskGrid(lonIndexGrid,latIndexGrid)
            call utl_abort('int_sintCloudToGrid_gsv')
          end if

        end do
      end do
    end if

    if ( checkCloudToGridUnassigned ) then
      write(*,*) 'Total number of grid points:                                   ', niGrid*njGrid
      write(*,*) 'Number of grid points not covered by the cloud of points:      ', sum(nhole(:))
      write(*,*) 'Number of grid points filled by neighbours:                    ', sum(nfill(:))
      write(*,*) 'Number of grid points with extrapolated value in masked area:  ', nextrap0
      write(*,*) 'Number of grid points with extrapolated value in visible area: ', nextrap1

      if ( nextrap1 > 0 ) then
        call utl_abort('int_sintCloudToGrid_gsv: Values at some unmasked grid points were not assigned')
      end if
    end if

    deallocate(fieldGrid_tmp)
    deallocate(xCloud)
    deallocate(yCloud)
    deallocate(numFilledByAvg)
    deallocate(filledByInterp)
    deallocate(maskCloud)
    deallocate(maskGrid)

    ierr = 0

    call utl_tmg_stop(176)

  end function int_sintCloudToGrid_gsv

  !--------------------------------------------------------------------------
  ! int_hInterpScalar_r8_2d
  !--------------------------------------------------------------------------
  function int_hInterpScalar_r8_2d(fieldOut_r8, fieldIn_r8, interpDegree, extrapDegree_opt) result(ierr)
    !
    ! :Purpose: Horizontal interpolation of 2D scalar field that use real(8) arrays
    !           for input and output. Accessed through int_hInterpScalar.
    !
    implicit none

    ! Arguments:
    real(8), intent(inout) :: fieldOut_r8(:,:)
    real(8), intent(in)    :: fieldIn_r8(:,:)
    integer :: ierr
    character(len=*)           :: interpDegree
    character(len=*), optional :: extrapDegree_opt

    ! Locals:
    integer :: nii, nji, nio, njo     
    integer :: jk1, jk2
    real(4), allocatable :: bufferi4(:,:), buffero4(:,:)
    integer :: ezsint

    ! read the namelist
    call int_readNml()

    ! do the standard interpolation

    call int_setezopt(interpDegree, extrapDegree_opt)   

    nii = size(fieldIn_r8,1)
    nji = size(fieldIn_r8,2)

    nio = size(fieldOut_r8,1)
    njo = size(fieldOut_r8,2)

    allocate(bufferi4(nii,nji))
    allocate(buffero4(nio,njo))

    do jk2 = 1,nji
      do jk1 = 1,nii
        bufferi4(jk1,jk2) = fieldIn_r8(jk1,jk2)
      end do
    end do

    ierr = ezsint(buffero4,bufferi4)

    do jk2 = 1,njo
      do jk1 = 1,nio
        fieldOut_r8(jk1,jk2) = buffero4(jk1,jk2)
      end do
    end do

    deallocate(bufferi4)
    deallocate(buffero4)

  end function int_hInterpScalar_r8_2d

  !--------------------------------------------------------------------------
  ! int_hInterpUV_gsv
  !--------------------------------------------------------------------------
  function int_hInterpUV_gsv(stateVectorOut, stateVectorIn, varName, levIndex, stepIndex, &
                             interpDegree, extrapDegree_opt) result(ierr)
    !
    ! :Purpose: Horizontal interpolation of 2D vector field that use stateVector objects
    !           for input and output. Accessed through int_hInterpUV.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: stateVectorOut
    type(struct_gsv), intent(inout) :: stateVectorIn
    character(len=*), intent(in)    :: varName
    integer         , intent(in)    :: levIndex
    integer         , intent(in)    :: stepIndex
    character(len=*), intent(in)           :: interpDegree
    character(len=*), intent(in), optional :: extrapDegree_opt
    integer :: ierr

    ! Locals:
    real(4), pointer :: UUout4(:,:,:,:), VVout4(:,:,:,:), UUin4(:,:,:,:), VVin4(:,:,:,:)
    real(8), pointer :: UUout8(:,:,:,:), VVout8(:,:,:,:), UUin8(:,:,:,:), VVin8(:,:,:,:)
    real(4), pointer :: UVout4(:,:,:), UVin4(:,:,:)
    real(8), pointer :: UVout8(:,:,:), UVin8(:,:,:)
    integer :: ezuvint, ezdefset

    write(*,*) 'int_hInterpUV_gsv: starting'

    ! read the namelist
    call int_readNml()

    ! check if special interpolation is required
    if (stateVectorIn%hco%initialized .and. stateVectorOut%hco%initialized) then
      if (stateVectorIn%hco%grtyp == 'Y') then
        call utl_abort('int_hInterpUV_gsv: cloudToGrid not implemented')
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
        call utl_abort('int_hInterpUV_gsv: unexpected varName')
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
        call utl_abort('int_hInterpUV_gsv: unexpected varName')
      end if

      deallocate(UUin4,VVin4,UUout4,VVout4)

    else

      call utl_abort('int_hInterpUV_gsv: not implemented for mixed dataKind')

    end if

  end function int_hInterpUV_gsv

  !--------------------------------------------------------------------------
  ! int_hInterpUV_r4_2d
  !--------------------------------------------------------------------------
  function int_hInterpUV_r4_2d(uuout, vvout, uuin, vvin, interpDegree, extrapDegree_opt) result(ierr)
    !
    ! :Purpose: Horizontal interpolation of 2D vector field that use real(4) arrays
    !           for input and output. Accessed through int_hInterpUV.
    !
    implicit none

    ! Arguments:
    real(4), intent(inout) :: uuout(:,:)
    real(4), intent(inout) :: vvout(:,:)
    real(4), intent(in)    :: uuin(:,:)
    real(4), intent(in)    :: vvin(:,:)
    character(len=*), intent(in)           :: interpDegree
    character(len=*), intent(in), optional :: extrapDegree_opt
    integer :: ierr

    ! Locals:
    integer :: ezuvint

    ! read the namelist
    call int_readNml()

    ! do the standard interpolation
    call int_setezopt(interpDegree, extrapDegree_opt)   
    ierr = ezuvint(uuout, vvout, uuin, vvin)

  end function int_hInterpUV_r4_2d

  !--------------------------------------------------------------------------
  ! int_hInterpUV_r8_2d
  !--------------------------------------------------------------------------
  function int_hInterpUV_r8_2d(uuout, vvout, uuin, vvin, interpDegree, extrapDegree_opt) result(ierr)
    !
    ! :Purpose: Horizontal interpolation of 2D vector field that use real(8) arrays
    !           for input and output. Accessed through int_hInterpUV.
    !
    implicit none

    ! Arguments:
    real(8), intent(inout) :: uuout(:,:)
    real(8), intent(inout) :: vvout(:,:)
    real(8), intent(in)    :: uuin(:,:)
    real(8), intent(in)    :: vvin(:,:)
    character(len=*), intent(in)           :: interpDegree
    character(len=*), intent(in), optional :: extrapDegree_opt
    integer :: ierr

    ! Locals:
    integer :: nio, njo, nii, nji
    integer :: jk1, jk2
    real, allocatable :: bufuuout4(:,:), bufvvout4(:,:)
    real, allocatable :: bufuuin4(:,:), bufvvin4(:,:)
    integer :: ezuvint

    ! read the namelist
    call int_readNml()

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

  end function int_hInterpUV_r8_2d

  !--------------------------------------------------------------------------
  ! int_ezgdef
  !--------------------------------------------------------------------------
  function int_ezgdef(ni, nj, grtyp, grtypref, ig1, ig2, ig3, ig4, ax, ay) result(vezgdef)
    !
    ! :Purpose: Subroutine wrapper for rmnlib procedure ezgdef.
    !
    implicit none

    ! Arguments:
    integer :: vezgdef
    integer, intent(in) :: ni
    integer, intent(in) :: nj
    integer, intent(in) :: ig1
    integer, intent(in) :: ig2
    integer, intent(in) :: ig3
    integer, intent(in) :: ig4
    real(8), intent(in) :: ax(:)
    real(8), intent(in) :: ay(:)
    character(len=*), intent(in) :: grtyp
    character(len=*), intent(in) :: grtypref

    ! Locals:
    integer :: ier2, jk, ilenx, ileny
    real(4), allocatable :: bufax4(:), bufay4(:)
    integer :: ezgdef

    if (grtyp .eq. 'Y') then
      ilenx = max(1,ni*nj)
      ileny = ilenx
    else if (grtyp .eq. 'Z') then
      ilenx = max(1,ni)
      ileny = max(1,nj)
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

    vezgdef = ier2

  end function int_ezgdef

  !--------------------------------------------------------------------------
  ! int_cxgaig
  !--------------------------------------------------------------------------
  subroutine int_cxgaig(grtyp, ig1, ig2, ig3, ig4, xlat0, xlon0, dlat, dlon) 
    !
    ! :Purpose: Subroutine wrapper for rmnlib procedure cxgaig.
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: ig1
    integer, intent(in) :: ig2
    integer, intent(in) :: ig3
    integer, intent(in) :: ig4   
    real(8), intent(in) :: xlat0
    real(8), intent(in) :: xlon0
    real(8), intent(in) :: dlat
    real(8), intent(in) :: dlon 
    character(len=*), intent(in) :: grtyp 

    ! Locals:
    real(4) :: xlat04, xlon04, dlat4, dlon4

    xlat04 = xlat0
    xlon04 = xlon0
    dlat4 = dlat
    dlon4 = dlon

    call cxgaig(grtyp, ig1, ig2, ig3, ig4, xlat04, xlon04, dlat4, dlon4)

  end subroutine int_cxgaig

end module interpolation_mod

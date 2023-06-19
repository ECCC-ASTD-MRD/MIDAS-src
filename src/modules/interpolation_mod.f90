
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
  use message_mod
  use kdTree2_mod
  implicit none
  save
  private

  ! public subroutines and functions
  public :: int_interp_gsv
  public :: int_hInterp_gsv
  public :: int_vInterp_gsv
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
    call msg('int_readNml', 'START', verb_opt=2)

    ! default values
    vInterpCopyLowestLevel = .false.
    checkCloudToGridUnassigned = .true.
    maxBoxSize = 1

    if ( .not. utl_isNamelistPresent('NAMINT','./flnml') ) then
      call msg('int_readNml', 'namint is missing in the namelist.'&
           //'The default values will be taken.', mpiAll_opt=.false.)
    else
      ! Read namelist NAMINT
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam,nml=namint,iostat=ierr)
      if (ierr /= 0) call utl_abort('int_readlNml: Error reading namelist NAMINT')
      if (mmpi_myid == 0) write(*,nml=namint)
      ierr = fclos(nulnam)
    end if

    call msg('int_readNml', 'END', verb_opt=2)
  end subroutine int_readNml

  !--------------------------------------------------------------------------
  ! int_interp_gsv
  !--------------------------------------------------------------------------
  subroutine int_interp_gsv(statevector_in, statevector_out, statevectorRef_opt, &
                            checkModelTop_opt)
    !
    ! :Purpose: high-level interpolation subroutine that proceed with
    !           horizontal and vertical interpolation
    !
    implicit none

    ! Arguments:
    type(struct_gsv),           intent(in)    :: statevector_in     ! Statevector input
    type(struct_gsv),           intent(inout) :: statevector_out    ! Statevector output (with target grids)
    type(struct_gsv), optional, intent(in)    :: statevectorRef_opt ! Reference statevector with fields P0, TT and HU
    logical,          optional, intent(in)    :: checkModelTop_opt  ! If true, model top consistency checked
    
    ! Locals:
    logical :: checkModelTop
    character(len=4), pointer :: varNamesToInterpolate(:)
    type(struct_gsv) :: statevector_in_varsLevs, statevector_in_varsLevs_hInterp
    type(struct_gsv) :: statevector_in_hInterp

    call msg('int_interp_gsv', 'START', verb_opt=2)

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
                      dataKind_opt=gsv_getDataKind(statevector_in), &
                      allocHeightSfc_opt=statevector_in%heightSfcPresent, &
                      varNames_opt=varNamesToInterpolate, &
                      hInterpolateDegree_opt=statevector_out%hInterpolateDegree, &
                      hExtrapolateDegree_opt=statevector_out%hExtrapolateDegree )

    call gsv_transposeTilesToVarsLevs( statevector_in, statevector_in_VarsLevs )

    call gsv_allocate(statevector_in_VarsLevs_hInterp, statevector_in%numstep, &
                      statevector_out%hco, statevector_in%vco,  &
                      mpi_local_opt=statevector_out%mpi_local, mpi_distribution_opt='VarsLevs', &
                      dataKind_opt=gsv_getDataKind(statevector_out), &
                      allocHeightSfc_opt=statevector_out%heightSfcPresent, &
                      varNames_opt=varNamesToInterpolate, &
                      hInterpolateDegree_opt=statevector_out%hInterpolateDegree, &
                      hExtrapolateDegree_opt=statevector_out%hExtrapolateDegree )

    call int_hInterp_gsv(statevector_in_VarsLevs, statevector_in_VarsLevs_hInterp)
    call gsv_deallocate(statevector_in_VarsLevs)

    call gsv_allocate(statevector_in_hInterp, statevector_in%numstep, &
                      statevector_out%hco, statevector_in%vco,      &
                      mpi_local_opt=statevector_out%mpi_local, mpi_distribution_opt='Tiles', &
                      dataKind_opt=gsv_getDataKind(statevector_out), &
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
    
    call int_vInterp_gsv(statevector_in_hInterp,statevector_out,&
                         statevectorRef_opt=statevectorRef_opt, &
                         checkModelTop_opt=checkModelTop)

    call gsv_deallocate(statevector_in_hInterp)
    nullify(varNamesToInterpolate)

    call msg('int_interp_gsv', 'END', verb_opt=2)
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
    type(struct_gsv), intent(inout) :: statevector_in  ! Statevector input
    type(struct_gsv), intent(inout) :: statevector_out ! Statevector with target horiz and vert grids and result

    ! Locals:
    integer :: varIndex, levIndex, nlev, stepIndex, ierr, kIndex
    character(len=4) :: varName
    character(len=12):: interpolationDegree, extrapolationDegree

    call msg('int_hInterp_gsv', 'START', verb_opt=2)

    if ( hco_equal(statevector_in%hco,statevector_out%hco) ) then
      call msg('int_hInterp_gsv', 'The input and output statevectors are already on same horizontal grids')
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

      call msg('int_hInterp_gsv', 'before interpolation (no mpi)')

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

      call msg('int_hInterp_gsv', 'before interpolation (with mpi)')

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
      call msg('int_hInterp_gsv','interpolating surface height')
      ierr = int_hInterpScalar( statevector_out, statevector_in, 'ZSFC', 1, 1, &
                                interpDegree=trim(interpolationDegree), &
                                extrapDegree_opt=trim(extrapolationDegree) )
    end if

    call msg('int_hInterp_gsv', 'END', verb_opt=2)
  end subroutine int_hInterp_gsv

  !--------------------------------------------------------------------------
  ! int_vInterp_gsv
  !--------------------------------------------------------------------------
  subroutine int_vInterp_gsv(statevector_in,statevector_out,statevectorRef_opt, &
                             Ps_in_hPa_opt,checkModelTop_opt)
    !
    ! :Purpose: Vertical interpolation.
    !           Interpolation coordinates are either height or log(P)
    !           based on target vertical descriptor vcode.
    !           * When target vertical coordinates are pressure (interpolating
    !             to GEM-P), log(P) is used;
    !           * When they are height (interpolating to GEM-H), height is used.
    !           Call to ``calcHeightAndPressure_mod`` will return required
    !           interpolation coordinates (log is then applied when required).
    !
    implicit none

    ! Arguments:
    type(struct_gsv),           target, intent(in)    :: statevector_in     ! Statevector input
    type(struct_gsv),                   intent(inout) :: statevector_out    ! Statevector with target horiz/vert grids and result
    type(struct_gsv), optional, target, intent(in)    :: statevectorRef_opt ! Reference statevector with P0, TT and HU
    logical,          optional,         intent(in)    :: Ps_in_hPa_opt      ! If true, surface pressure is in hPa, not Pa
    logical,          optional,         intent(in)    :: checkModelTop_opt  ! Model top consistency will be checked

    call msg('int_vInterp_gsv', 'START', verb_opt=2)

    ! read the namelist
    call int_readNml()

    if ( gsv_getDataKind(statevector_in) == 8 & 
         .and. gsv_getDataKind(statevector_out) == 8 ) then 
      call vInterp_gsv_r8(statevector_in,statevector_out,statevectorRef_opt, &
                          Ps_in_hPa_opt,checkModelTop_opt)
    else if ( gsv_getDataKind(statevector_in) == 4 &
         .and. gsv_getDataKind(statevector_out) == 4 ) then 
      call vInterp_gsv_r4(statevector_in,statevector_out,statevectorRef_opt, &
                          Ps_in_hPa_opt,checkModelTop_opt)
    else
      call utl_abort('int_vInterp_gsv: input and output statevectors must be of same dataKind')
    end if

    call msg('int_vInterp_gsv', 'END', verb_opt=2)

  end subroutine int_vInterp_gsv

  !--------------------------------------------------------------------------
  ! vInterp_gsv_r8
  !--------------------------------------------------------------------------
  subroutine vInterp_gsv_r8(statevector_in,statevector_out,statevectorRef_opt, &
                            Ps_in_hPa_opt,checkModelTop_opt)
    !
    ! :Purpose: Vertical interpolation, ``real(8)`` version.
    !
    implicit none

    ! Arguments:
    type(struct_gsv),           target, intent(in)    :: statevector_in     ! Statevector input
    type(struct_gsv),                   intent(inout) :: statevector_out    ! Statevector with target horiz/vert grids and result
    type(struct_gsv), optional, target, intent(in)    :: statevectorRef_opt ! Reference statevector with P0, TT and HU
    logical,          optional,         intent(in)    :: Ps_in_hPa_opt      ! If true, surface pressure is in hPa, not Pa
    logical,          optional,         intent(in)    :: checkModelTop_opt  ! Model top consistency will be checked

    ! Locals:
    logical :: checkModelTop, hLikeCalc
    integer :: vcode_in, vcode_out
    integer :: nlev_out, nlev_in
    integer :: varIndex, stepIndex
    type(struct_gsv), pointer   :: statevectorRef
    type(struct_gsv)            :: statevectorRef_out
    real(8), pointer  :: hLikeT_in(:,:,:,:), hLikeM_in(:,:,:,:)   ! abstract height dimensioned coordinate
    real(8), pointer  :: hLikeT_out(:,:,:,:), hLikeM_out(:,:,:,:) ! abstract height dimensioned coordinate
    real(8), pointer  :: field_in(:,:,:,:), field_out(:,:,:,:)
    real(8), pointer  :: heightSfcIn(:,:), heightSfcOut(:,:)
    real(8), pointer  :: tmpCoord_T(:,:,:,:), tmpCoord_M(:,:,:,:)
    character(len=4) :: varName
    type(struct_vco), pointer :: vco_in, vco_out

    call msg('vInterp_gsv_r8', 'START', verb_opt=4)

    vco_in  => gsv_getVco(statevector_in)
    vco_out => gsv_getVco(statevector_out)
    vcode_in  = vco_in%vcode
    vcode_out = vco_out%vcode
    nullify(hLikeT_in,hLikeM_in,hLikeT_out,hLikeM_out)

    hLikeCalc = .true.

    if (present(statevectorRef_opt)) then
      if ( .not. vco_equal(gsv_getVco(statevectorRef_opt), gsv_getVco(statevector_in))) then
        call utl_abort('vInterp_gsv_r8: reference must have input vertical structure')
      end if
      statevectorRef => statevectorRef_opt
    else
      statevectorRef => statevector_in
    end if

    if ( vco_equal(statevector_in%vco, statevector_out%vco) ) then
      call msg('vInterp_gsv_r8', 'The input and output statevectors are already on same vertical levels')
      call gsv_copy(statevector_in, statevector_out)
      return
    end if

    if ( .not. hco_equal(statevector_in%hco, statevector_out%hco) ) then
      call utl_abort('vInterp_gsv_r8: The input and output statevectors are not on the same horizontal grid.')
    end if

    if ( gsv_getDataKind(statevector_in) /= 8 .or. gsv_getDataKind(statevector_out) /= 8 ) then
      call utl_abort('vInterp_gsv_r8: Incorrect value for dataKind. Only compatible with dataKind=8')
    end if

    if (gsv_isAssocHeightSfc(statevector_in) .and. gsv_isAssocHeightSfc(statevector_out) ) then
      heightSfcIn => gsv_getHeightSfc(statevector_in)
      heightSfcOut => gsv_getHeightSfc(statevector_out)
      heightSfcOut(:,:) = heightSfcIn(:,:)
    end if

    ! the default is to ensure that the top of the output grid is ~equal or lower than the top of the input grid 
    if ( present(checkModelTop_opt) ) then
      checkModelTop = checkModelTop_opt
    else
      checkModelTop = .true.
    end if
    if (checkModelTop) then
      call msg('vInterp_gsv_r8', ' Checking that that the top of the destination grid is not higher than the top of the source grid.')
      if ( vcode_in == 21001 .or. vcode_out == 21001 ) then
        call msg('vInterp_gsv_r8', 'bypassing top check, '&
             //'vcode_in='//str(vcode_in)//', vcode_out='//str(vcode_out))
        ! Development notes (@mad001)
        !   we should consider having a new criterion that works for GEM-H as well
      else
        call czp_ensureCompatibleTops(vco_in, vco_out)
      end if
    end if

    ! create reference state on output vco with reference surface fields
    call gsv_allocate(statevectorRef_out, statevectorRef%numstep, &
                      statevectorRef%hco, statevector_out%vco,      &
                      mpi_local_opt=statevectorRef%mpi_local, mpi_distribution_opt='Tiles', &
                      dataKind_opt=gsv_getDataKind(statevectorRef), &
                      allocHeightSfc_opt=statevectorRef%heightSfcPresent, &
                      varNames_opt=(/'P0','P0LS'/) )
    call gsv_copy(stateVectorRef, stateVectorRef_out, allowVcoMismatch_opt=.true., &
                  allowVarMismatch_opt=.true.)

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

      hLikeCalcIf: if (hLikeCalc) then
        hLikeCalc = .false.

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
        allocate(tmpCoord_T(statevector_in%myLonBeg:statevector_in%myLonEnd, &
                            statevector_in%myLatBeg:statevector_in%myLatEnd, &
                            gsv_getNumLev(statevector_in,'TH'), statevector_in%numStep))
        allocate(tmpCoord_M(statevector_in%myLonBeg:statevector_in%myLonEnd, &
                            statevector_in%myLatBeg:statevector_in%myLatEnd, &
                            gsv_getNumLev(statevector_in,'MM'), statevector_in%numStep))

        ! output grid GEM-P interpolation in log-pressure
        if ( vcode_out==5002 .or. vcode_out==5005 .or. vcode_out==5100 ) then
          call czp_calcReturnPressure_gsv_nl( statevectorRef_out, &
                                              PTout_r8_opt=hLikeT_out, &
                                              PMout_r8_opt=hLikeM_out, &
                                              Ps_in_hPa_opt=Ps_in_hPa_opt)

          if ( vcode_in==5002 .or. vcode_in==5005 .or. vcode_in==5100 ) then
            call czp_calcReturnPressure_gsv_nl( statevectorRef, &
                                                PTout_r8_opt=hLikeT_in, &
                                                PMout_r8_opt=hLikeM_in, &
                                                Ps_in_hPa_opt=Ps_in_hPa_opt)
          else if ( vcode_in==21001 ) then
            call czp_calcReturnHeight_gsv_nl( statevectorRef, &
                                              ZTout_r8_opt=tmpCoord_T, &
                                              ZMout_r8_opt=tmpCoord_M)
            call czp_calcReturnPressure_gsv_nl( statevectorRef, &
                                                ZTin_r8_opt=tmpCoord_T, &
                                                ZMin_r8_opt=tmpCoord_M, &
                                                PTout_r8_opt=hLikeT_in, &
                                                PMout_r8_opt=hLikeM_in, &
                                                Ps_in_hPa_opt=Ps_in_hPa_opt)
          end if

          call msg('vInterp_gsv_r8','converting pressure coordinates to height-like, '&
                   //'vcode_in='//str(vcode_in)//', vcode_out='//str(vcode_out))
          call logP_r8(hLikeM_out)
          call logP_r8(hLikeT_out)
          call logP_r8(hLikeM_in)
          call logP_r8(hLikeT_in)

          ! output grid GEM-H interpolation in height
        else if ( vcode_out==21001 ) then
          call czp_calcReturnHeight_gsv_nl( statevectorRef_out, &
                                            ZTout_r8_opt=hLikeT_out, &
                                            ZMout_r8_opt=hLikeM_out)
          if ( vcode_in==21001 ) then
            call czp_calcReturnHeight_gsv_nl( statevectorRef, &
                                              ZTout_r8_opt=hLikeT_in, &
                                              ZMout_r8_opt=hLikeM_in)
          else if ( vcode_in==5002 .or. vcode_in==5005 .or. vcode_in==5100 ) then
            call czp_calcReturnPressure_gsv_nl( statevectorRef, &
                                                PTout_r8_opt=tmpCoord_T, &
                                                PMout_r8_opt=tmpCoord_M, &
                                                Ps_in_hPa_opt=Ps_in_hPa_opt)
            call czp_calcReturnHeight_gsv_nl( statevectorRef, &
                                              PTin_r8_opt=tmpCoord_T, &
                                              PMin_r8_opt=tmpCoord_M, &
                                              ZTout_r8_opt=hLikeT_in, &
                                              ZMout_r8_opt=hLikeM_in)
          end if
        end if
        deallocate(tmpCoord_T)
        deallocate(tmpCoord_M)

      end if hLikeCalcIf

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
          call hLike_interpolation_r8(hLikeT_in, hLikeT_out)
        else
          call hLike_interpolation_r8(hLikeM_in, hLikeM_out)
        end if

        ! overwrite values at the lowest levels to avoid extrapolation
        if (vInterpCopyLowestLevel) then
          field_out(:,:,nlev_out,stepIndex) = field_in(:,:,nlev_in,stepIndex)
        end if

      end do step_loop

    end do var_loop

    if (associated(hLikeT_in)) then
      deallocate(hLikeT_in, hLikeM_in, hLikeT_out, hLikeM_out)
    end if
    call gsv_deallocate(statevectorRef_out)

    call msg('vInterp_gsv_r8', 'END', verb_opt=4)

    contains

      subroutine hLike_interpolation_r8(hLike_in, hLike_out)
        !
        ! :Purpose: Proceed to actual interpolation in H-logP representation
        !
        implicit none

        ! Arguments:
        real(8), pointer, intent(in) :: hLike_in(:,:,:,:)  ! abstract height dimensioned input coordinate
        real(8), pointer, intent(in) :: hLike_out(:,:,:,:) ! abstract height dimensioned target coordinate

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

      end subroutine hLike_interpolation_r8

  end subroutine vInterp_gsv_r8

  !--------------------------------------------------------------------------
  ! logP_r8
  !--------------------------------------------------------------------------
  subroutine logP_r8(presInLogOut)
    !
    ! :Purpose: compute log of pressurce field, real(8) version
    !
    implicit none

    ! Arguments:
    real(8), intent(inout) :: presInLogOut(:,:,:,:)

    ! Locals:
    integer :: latIdx, lonIdx, levIdx, stepIdx

    !$OMP PARALLEL DO PRIVATE(lonIdx,latIdx,levIdx,stepIdx)
    do stepIdx = lbound(presInLogOut,4), ubound(presInLogOut,4)
      do levIdx = lbound(presInLogOut,3), ubound(presInLogOut,3)
        do latIdx = lbound(presInLogOut,2), ubound(presInLogOut,2)
          do lonIdx = lbound(presInLogOut,1), ubound(presInLogOut,1)
            presInLogOut(lonIdx,latIdx,levIdx,stepIdx) = &
                 log(presInLogOut(lonIdx,latIdx,levIdx,stepIdx))
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine logP_r8

  !--------------------------------------------------------------------------
  ! vInterp_gsv_r4
  !--------------------------------------------------------------------------
  subroutine vInterp_gsv_r4(statevector_in,statevector_out,statevectorRef_opt, &
                            Ps_in_hPa_opt,checkModelTop_opt)
    !
    ! :Purpose: Vertical interpolation, ``real(4)`` version.
    !
    !
    implicit none

    ! Arguments:
    type(struct_gsv),           target, intent(in)    :: statevector_in     ! Statevector input
    type(struct_gsv),                   intent(inout) :: statevector_out    ! Statevector with the target horiz/vert grids and result
    type(struct_gsv), optional, target, intent(in)    :: statevectorRef_opt ! Reference statevector with P0, TT and HU
    logical,          optional,         intent(in)    :: Ps_in_hPa_opt      ! If true, surface pressure in in hPa, not Pa
    logical,          optional,         intent(in)    :: checkModelTop_opt  ! Model top consistency will be checked

    ! Locals:
    logical :: checkModelTop, hLikeCalc
    integer :: vcode_in, vcode_out
    integer :: nlev_out, nlev_in
    integer :: varIndex, stepIndex
    type(struct_gsv), pointer   :: statevectorRef
    type(struct_gsv)            :: statevectorRef_out
    real(4), pointer  :: hLikeT_in(:,:,:,:), hLikeM_in(:,:,:,:)   ! abstract height dimensioned coordinate
    real(4), pointer  :: hLikeT_out(:,:,:,:), hLikeM_out(:,:,:,:) ! abstract height dimensioned coordinate
    real(4), pointer  :: field_in(:,:,:,:), field_out(:,:,:,:)
    real(8), pointer  :: heightSfcIn(:,:), heightSfcOut(:,:)
    real(4), pointer  :: tmpCoord_T(:,:,:,:), tmpCoord_M(:,:,:,:)
    character(len=4) :: varName
    type(struct_vco), pointer :: vco_in, vco_out

    call msg('vInterp_gsv_r4', 'START', verb_opt=4)

    vco_in  => gsv_getVco(statevector_in)
    vco_out => gsv_getVco(statevector_out)
    vcode_in  = vco_in%vcode
    vcode_out = vco_out%vcode
    nullify(hLikeT_in,hLikeM_in,hLikeT_out,hLikeM_out)

    hLikeCalc = .true.

    if (present(statevectorRef_opt)) then
      if ( .not. vco_equal(gsv_getVco(statevectorRef_opt), gsv_getVco(statevector_in))) then
        call utl_abort('vInterp_gsv_r4: reference must have input vertical structure')
      end if
      statevectorRef => statevectorRef_opt
    else
      statevectorRef => statevector_in
    end if

    if ( vco_equal(statevector_in%vco, statevector_out%vco) ) then
      call msg('vInterp_gsv_r4', 'The input and output statevectors are already on same vertical levels')
      call gsv_copy(statevector_in, statevector_out)
      return
    end if

    if ( .not. hco_equal(statevector_in%hco, statevector_out%hco) ) then
      call utl_abort('vInterp_gsv_r4: The input and output statevectors are not on the same horizontal grid.')
    end if

    ! DBGmad remove?
    if ( gsv_getDataKind(statevector_in) /= 4 .or. gsv_getDataKind(statevector_out) /= 4 ) then
      call utl_abort('vInterp_gsv_r4: Incorrect value for dataKind. Only compatible with dataKind=4')
    end if

    if (gsv_isAssocHeightSfc(statevector_in) .and. gsv_isAssocHeightSfc(statevector_out) ) then
      heightSfcIn => gsv_getHeightSfc(statevector_in)
      heightSfcOut => gsv_getHeightSfc(statevector_out)
      heightSfcOut(:,:) = heightSfcIn(:,:)
    end if

    ! DBGmad move to int_vInterp_gsv?
    ! the default is to ensure that the top of the output grid is ~equal or lower than the top of the input grid 
    if ( present(checkModelTop_opt) ) then
      checkModelTop = checkModelTop_opt
    else
      checkModelTop = .true.
    end if
    if (checkModelTop) then
      call msg('vInterp_gsv_r4', ' Checking that that the top of the destination grid is not higher than the top of the source grid.')
      if ( vcode_in == 21001 .or. vcode_out == 21001 ) then
        call msg('vInterp_gsv_r4', 'bypassing top check, '&
             //'vcode_in='//str(vcode_in)//', vcode_out='//str(vcode_out))
        ! Development notes (@mad001)
        !   we should consider having a new criterion that works for GEM-H as well
      else
        call czp_ensureCompatibleTops(vco_in, vco_out)
      end if
    end if

    ! create reference state on output vco with reference surface fields
    call gsv_allocate(statevectorRef_out, statevectorRef%numstep, &
                      statevectorRef%hco, statevector_out%vco,      &
                      mpi_local_opt=statevectorRef%mpi_local,  &
                      mpi_distribution_opt=statevectorRef%mpi_distribution, &
                      dataKind_opt=gsv_getDataKind(statevectorRef), &
                      allocHeightSfc_opt=statevectorRef%heightSfcPresent, &
                      varNames_opt=(/'P0','P0LS'/) )
    call gsv_copy(stateVectorRef, stateVectorRef_out, allowVcoMismatch_opt=.true., &
                  allowVarMismatch_opt=.true.)

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

      hLikeCalcIf: if (hLikeCalc) then
        hLikeCalc = .false.

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
        allocate(tmpCoord_T(statevector_in%myLonBeg:statevector_in%myLonEnd, &
                            statevector_in%myLatBeg:statevector_in%myLatEnd, &
                            gsv_getNumLev(statevector_in,'TH'), statevector_in%numStep))
        allocate(tmpCoord_M(statevector_in%myLonBeg:statevector_in%myLonEnd, &
                            statevector_in%myLatBeg:statevector_in%myLatEnd, &
                          gsv_getNumLev(statevector_in,'MM'), statevector_in%numStep))

        ! output grid GEM-P interpolation in log-pressure
        if ( vcode_out==5002 .or. vcode_out==5005 .or. vcode_out==5100 ) then
          call czp_calcReturnPressure_gsv_nl( statevectorRef_out, &
                                              PTout_r4_opt=hLikeT_out, &
                                              PMout_r4_opt=hLikeM_out, &
                                              Ps_in_hPa_opt=Ps_in_hPa_opt)

          if ( vcode_in==5002 .or. vcode_in==5005 .or. vcode_in==5100 ) then
            call czp_calcReturnPressure_gsv_nl( statevectorRef, &
                                                PTout_r4_opt=hLikeT_in, &
                                                PMout_r4_opt=hLikeM_in, &
                                                Ps_in_hPa_opt=Ps_in_hPa_opt)
          else if ( vcode_in==21001 ) then
            call czp_calcReturnHeight_gsv_nl( statevectorRef, &
                                              ZTout_r4_opt=tmpCoord_T, &
                                              ZMout_r4_opt=tmpCoord_M)
            call czp_calcReturnPressure_gsv_nl( statevectorRef, &
                                                ZTin_r4_opt=tmpCoord_T, &
                                                ZMin_r4_opt=tmpCoord_M, &
                                                PTout_r4_opt=hLikeT_in, &
                                                PMout_r4_opt=hLikeM_in, &
                                                Ps_in_hPa_opt=Ps_in_hPa_opt)
          end if

          call msg('vInterp_gsv_r4','converting pressure coordinates to height-like, '&
                   //'vcode_in='//str(vcode_in)//', vcode_out='//str(vcode_out))
          call logP_r4(hLikeM_out)
          call logP_r4(hLikeT_out)
          call logP_r4(hLikeM_in)
          call logP_r4(hLikeT_in)

          ! output grid GEM-H interpolation in height
        else if ( vcode_out==21001 ) then
          call czp_calcReturnHeight_gsv_nl( statevectorRef_out, &
                                            ZTout_r4_opt=hLikeT_out, &
                                            ZMout_r4_opt=hLikeM_out)
          if ( vcode_in==21001 ) then
            call czp_calcReturnHeight_gsv_nl( statevectorRef, &
                                              ZTout_r4_opt=hLikeT_in, &
                                              ZMout_r4_opt=hLikeM_in)
          else if ( vcode_in==5002 .or. vcode_in==5005 .or. vcode_in==5100 ) then
            call czp_calcReturnPressure_gsv_nl( statevectorRef, &
                                                PTout_r4_opt=tmpCoord_T, &
                                                PMout_r4_opt=tmpCoord_M, &
                                                Ps_in_hPa_opt=Ps_in_hPa_opt)
            call czp_calcReturnHeight_gsv_nl( statevectorRef, &
                                              PTin_r4_opt=tmpCoord_T, &
                                              PMin_r4_opt=tmpCoord_M, &
                                              ZTout_r4_opt=hLikeT_in, &
                                              ZMout_r4_opt=hLikeM_in)
          end if
        end if
        deallocate(tmpCoord_T)
        deallocate(tmpCoord_M)

      end if hLikeCalcIf

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

    end do var_loop

    if (associated(hLikeT_in)) then
      deallocate(hLikeT_in, hLikeM_in, hLikeT_out, hLikeM_out)
    end if
    call gsv_deallocate(statevectorRef_out)

    call msg('vInterp_gsv_r4', 'END', verb_opt=4)

    contains

      subroutine hLike_interpolation_r4(hLike_in, hLike_out)
        !
        ! :Purpose: Proceed to actual interpolation in H-logP representation
        !
        implicit none

        ! Arguments:
        real(4), pointer, intent(in)  :: hLike_in(:,:,:,:)  ! abstract height dimensioned input coordinate
        real(4), pointer, intent(in)  :: hLike_out(:,:,:,:) ! abstract height dimensioned target coordinate

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

  end subroutine vInterp_gsv_r4

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
    do stepIdx = lbound(presInLogOut,4), ubound(presInLogOut,4)
      do levIdx = lbound(presInLogOut,3), ubound(presInLogOut,3)
        do latIdx = lbound(presInLogOut,2), ubound(presInLogOut,2)
          do lonIdx = lbound(presInLogOut,1), ubound(presInLogOut,1)
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
    type(struct_gsv),  intent(in)    :: statevector_in  ! Statevector input
    type(struct_gsv),  intent(inout) :: statevector_out ! Statevector with target temporal structure and results

    ! Locals:
    integer :: kIndex, latIndex, lonIndex
    integer :: stepIndexIn1, stepIndexIn2, stepIndexOut, numStepIn, numStepOut
    integer :: lon1, lon2, lat1, lat2, k1, k2
    integer :: dateStampIn, dateStampOut 
    real(8) :: weight1, weight2
    real(8) :: deltaHour, deltaHourInOut
    real(4), pointer  :: gdIn_r4(:,:,:,:), gdOut_r4(:,:,:,:)
    real(8), pointer  :: gdIn_r8(:,:,:,:), gdOut_r8(:,:,:,:)

    call msg('int_tInterp_gsv', 'START', verb_opt=2)

    ! read the namelist
    call int_readNml()

    if ( .not. vco_equal(statevector_in%vco, statevector_out%vco) ) then
      call utl_abort('int_tInterp_gsv: The input and output statevectors are not on the same vertical levels')
    end if

    if ( .not. hco_equal(statevector_in%hco, statevector_out%hco) ) then
      call utl_abort('int_tInterp_gsv: The input and output statevectors are not on the same horizontal grid.')
    end if

    if ( statevector_in%numStep > statevector_out%numStep ) then
      call msg('int_tInterp_gsv', 'numStep_out is less than numStep_in, calling gsv_copy.')
      call gsv_copy(statevector_in, statevector_out, allowTimeMismatch_opt=.true.)
      return
    else if ( statevector_in%numStep == statevector_out%numStep ) then
      call msg('int_tInterp_gsv', 'numStep_out is equal to numStep_in, calling gsv_copy.')
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
    call msg('int_tInterp_gsv', 'numStepIn='//str(numStepIn)&
         //', numStepOut='//str(numStepOut), mpiAll_opt=.false.)

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
        deltaHourInOut = 0.0d0
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

      call msg('int_tInterp_gsv', 'for stepIndexOut='//str(stepIndexOut) &
           //', stepIndexIn1='//str(stepIndexIn1)//', stepIndexIn2='//str(stepIndexIn2) &
           //', weight1='//str(weight1)//', weight2='//str(weight2) &
           //', deltaHourInOut/deltaHour='//str(deltaHourInOut)//'/'//str(deltaHour), &
           mpiAll_opt=.false.)

      if ( gsv_getDataKind(statevector_in) == 4 .and. gsv_getDataKind(statevector_out) == 4 ) then
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
      else if ( gsv_getDataKind(statevector_in) == 4 .and. gsv_getDataKind(statevector_out) == 8 ) then
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
      else if ( gsv_getDataKind(statevector_in) == 8 .and. gsv_getDataKind(statevector_out) == 4 ) then
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
      else if ( gsv_getDataKind(statevector_in) == 8 .and. gsv_getDataKind(statevector_out) == 8 ) then
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

    call msg('int_tInterp_gsv', 'END', verb_opt=2)

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
    type(struct_columnData),  intent(in)    :: column_in              ! ColumnData input
    type(struct_columnData),  intent(inout) :: column_out             ! columnData with the vert structure and results
    character(len=*),         intent(in)    :: varName                ! variable name to be interpolated
    logical, optional,        intent(in)    :: useColumnPressure_opt  ! if .true. use P_* instead of the pressure provided by calcHeightAndPressure_mod

    ! Locals:
    real(8), pointer :: varInterp_in(:), varInterp_out(:)
    real(8), pointer :: coordRef_in(:)  , coordRef_out(:)
    real(8), pointer :: ZT_in(:,:), ZM_in(:,:), PT_in(:,:), PM_in(:,:)
    real(8), pointer :: PT_out(:,:), PM_out(:,:)
    character(len=4) :: varLevel
    real(8)          :: zwb, zwt
    integer          :: levIndex_out, levIndex_in, columnIndex
    integer          :: vcode_in, vcode_out
    integer          :: nLevIn_T, nLevIn_M, nLevOut_T, nLevOut_M
    logical          :: vInterp, useColumnPressure
    integer, allocatable, target :: THlevelWanted(:), MMlevelWanted(:)
    integer, pointer :: levelWanted(:)

    call msg('int_vInterp_col', 'START', verb_opt=2)
    varLevel = vnl_varLevelFromVarname(varName)

    if ( present(useColumnPressure_opt) ) then
      useColumnPressure = useColumnPressure_opt
    else
      useColumnPressure = .true.
    end if

    call msg('int_vInterp_col', varName//' ('//varLevel &
         //'), useColumnPressure='//str(useColumnPressure), verb_opt=3)

    vInterp = .true.
    if ( .not. col_varExist(column_in,'P0' ) ) then
      call msg('int_vInterp_col', 'P0 is missing. Vertical interpolation WILL NOT BE PERFORMED')
      vInterp = .false.
    else if ( col_getNumLev(column_in ,'TH') <= 1 .or. &
              col_getNumLev(column_in ,'MM') <= 1 ) then
      vInterp = .false.
      call msg('int_vInterp_col', 'The input backgrounds are 2D. Vertical interpolation WILL NOT BE PERFORMED')
    end if

    vcode_in  = column_in%vco%vcode
    vcode_out = column_out%vco%vcode
    nLevIn_T  = col_getNumLev(column_in,  'TH')
    nLevIn_M  = col_getNumLev(column_in,  'MM')
    nLevOut_T = col_getNumLev(column_out, 'TH')
    nLevOut_M = col_getNumLev(column_out, 'MM')

    if_vInterp: if (vInterp) then
      if ( .not. useColumnPressure ) then

        call msg('int_vInterp_col', 'vcode_in='//str(vcode_in)//', '&
             //'vcode_out='//str(vcode_out), verb_opt=3)

        ! Compute interpolation variables (pressures and or heights)
        allocate(PT_in( col_getNumCol(column_in),  nLevIn_T))
        allocate(PM_in( col_getNumCol(column_in),  nLevIn_M))
        allocate(ZT_in( col_getNumCol(column_in),  nLevIn_T))
        allocate(ZM_in( col_getNumCol(column_in),  nLevIn_M))
        allocate(PT_out(col_getNumCol(column_out), nLevOut_T))
        allocate(PM_out(col_getNumCol(column_out), nLevOut_M))

        call czp_calcReturnPressure_col_nl(column_in, PT_in, PM_in)
        call czp_calcReturnHeight_col_nl(column_in, ZT_in, ZM_in)
        call czp_calcReturnPressure_col_nl(column_out, PT_out, PM_out)

      end if    ! useColumnPressure

      do columnIndex = 1, col_getNumCol(column_out)

        ! coordRef_{in,out}
        if ( varLevel == 'TH' ) then
          if ( .not. useColumnPressure ) then
            coordRef_in  => PT_in(columnIndex,:)
            coordRef_out => PT_out(columnIndex,:)
          else
            coordRef_in  => col_getColumn(column_in,columnIndex,'P_T')
            coordRef_out => col_getColumn(column_out,columnIndex,'P_T')
          end if
        else if ( varLevel == 'MM' ) then
          if ( .not. useColumnPressure ) then
            coordRef_in  => PM_in(columnIndex,:)
            coordRef_out => PM_out(columnIndex,:)
          else
            coordRef_in  => col_getColumn(column_in,columnIndex,'P_M')
            coordRef_out => col_getColumn(column_out,columnIndex,'P_M')
          end if
        else
          call utl_abort('int_vInterp_col: only varLevel TH/MM is allowed')
        end if

        varInterp_in  => col_getColumn(column_in ,columnIndex,varName)
        varInterp_out => col_getColumn(column_out,columnIndex,varName)

        if ( columnIndex == 1 .and. (trim(varName) == 'P_T' ) ) then
          call msg('int_vInterp_col', 'useColumnPressure='//str(useColumnPressure) &
               //', '//trim(varName)//':' &
               //new_line('')//'COLUMN_IN(1):'//str(varInterp_in(:))&
               //new_line('')//'COLUMN_OUT(1):'//str(varInterp_out(:)), &
               mpiAll_opt=.false.)
        end if

        ! actual interpolation
        !   Development notes (@mad001)
        !     Potential issue with GEM-H height based interpolation
        !     we should consider to convert to pure logP interpolation
        !     as we did for int_vInterp_gsv
        !     see also #466 https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/issues/466#note_497052
        levIndex_in = 1
        do levIndex_out = 1, col_getNumLev(column_out,varLevel)
          levIndex_in = levIndex_in + 1
          do while( coordRef_out(levIndex_out) .gt. coordRef_in(levIndex_in) .and. &
               levIndex_in .lt. col_getNumLev(column_in,varLevel) )
            levIndex_in = levIndex_in + 1
          end do
          levIndex_in = levIndex_in - 1
          zwb = log(coordRef_out(levIndex_out)/coordRef_in(levIndex_in))/  &
               log(coordRef_in(levIndex_in+1)/coordRef_in(levIndex_in))
          zwt = 1. - zwb
          if (  useColumnPressure .and. &
              (trim(varName) == 'P_T' .or. trim(varName) == 'P_M' ) ) then
            ! do nothing, i.e. use the pressures from column_in
          else if ( .not. useColumnPressure .and. &
              (trim(varName) == 'P_T' .or. trim(varName) == 'P_M' ) ) then
            varInterp_out(levIndex_out) = exp(zwb*log(varInterp_in(levIndex_in+1)) + zwt*log(varInterp_in(levIndex_in)))
          else
            varInterp_out(levIndex_out) = zwb*varInterp_in(levIndex_in+1) + zwt*varInterp_in(levIndex_in)
          end if
        end do
      end do

      call msg('int_vInterp_col', trim(varName)//' (in): '//str(varInterp_in)&
           //new_line('')//trim(varName)//' (out): '//str(varInterp_out), verb_opt=3)

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
          varInterp_in  => col_getColumn(column_in ,columnIndex,varName)
          varInterp_out => col_getColumn(column_out,columnIndex,varName)
          if (vnl_varLevelFromVarname(varName) == 'TH') then
            levelWanted => THlevelWanted
          else
            levelWanted => MMlevelWanted
          end if
          do levIndex_out = 1, col_getNumLev(column_out,varLevel)
            varInterp_out(levIndex_out) = varInterp_in(levelWanted(levIndex_out))
          end do
        end do

        deallocate(THlevelWanted)
        deallocate(MMlevelWanted)

      else if (column_out%vco%nlev_depth > 0) then
        call msg('int_vInterp_col', 'vco_levelMatchingList: no MM and TH levels, but depth levels exist')
        if (any(column_out%vco%depths(:) /= column_in%vco%depths(:))) then
          call utl_abort('int_vInterp_col: some depth levels not equal')
        else
          ! copy over depth levels
          do columnIndex = 1, col_getNumCol(column_out)
            varInterp_in  => col_getColumn(column_in ,columnIndex,varName)
            varInterp_out => col_getColumn(column_out,columnIndex,varName)
            do levIndex_out = 1, col_getNumLev(column_out,varLevel)
              varInterp_out(levIndex_out) = varInterp_in(levIndex_out)
            end do
          end do
        end if

      end if

    end if if_vInterp

    call msg('int_vInterp_col', 'END', verb_opt=2)
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
    character(len=*),           intent(in) :: interpDegree
    character(len=*), optional, intent(in) :: extrapDegree_opt

    ! Locals:
    character(len=12) :: extrapDegree
    integer           :: ierr, ezsetopt, ezsetval

    call msg('int_setezopt', 'START', verb_opt=4)
    if ( trim(interpDegree) /= 'LINEAR' .and. &
         trim(interpDegree) /= 'CUBIC' .and. &
         trim(interpDegree) /= 'NEAREST' ) then

      call utl_abort('int_setezopt: invalid interpolation degree = '//trim(interpDegree))
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

    call msg('int_setezopt', 'END', verb_opt=4)
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
    type(struct_gsv),           intent(inout) :: stateVectorOut
    type(struct_gsv),           intent(inout) :: stateVectorIn
    character(len=*),           intent(in)    :: varName
    integer         ,           intent(in)    :: levIndex
    integer         ,           intent(in)    :: stepIndex
    character(len=*),           intent(in)    :: interpDegree
    character(len=*), optional, intent(in)    :: extrapDegree_opt
    ! Result:
    integer :: ierr

    ! Locals:
    real(4), pointer :: fieldOut_r4(:,:,:,:), fieldIn_r4(:,:,:,:)
    real(8), pointer :: fieldOut_r8(:,:,:,:), fieldIn_r8(:,:,:,:)
    real(8), pointer :: heightSfcOut(:,:), heightSfcIn(:,:)
    integer :: ezsint, ezdefset

    call msg('int_hInterpScalar_gsv', 'START', verb_opt=4)
    ! read the namelist
    call int_readNml()

    ! check if special interpolation is required
    if (stateVectorIn%hco%initialized .and. stateVectorOut%hco%initialized) then
      if (stateVectorIn%hco%grtyp == 'Y') then
        ! for now, only comptable for real(4)
        if( gsv_getDataKind(stateVectorOut) /= 4 .or. gsv_getDataKind(stateVectorIn) /= 4) then
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

    else if ( gsv_getDataKind(stateVectorOut) == 4 .and. gsv_getDataKind(stateVectorIn) == 4) then

      if (trim(varName) == 'ALL') then
        call gsv_getField(stateVectorOut, fieldOut_r4)
        call gsv_getField(stateVectorIn,  fieldIn_r4)
     else
        call gsv_getField(stateVectorOut, fieldOut_r4, varName)
        call gsv_getField(stateVectorIn,  fieldIn_r4,  varName)
      end if

      ierr = ezsint(fieldOut_r4(:,:,levIndex,stepIndex),fieldIn_r4(:,:,levIndex,stepIndex))

    else if ( gsv_getDataKind(stateVectorOut) == 8 .and. gsv_getDataKind(stateVectorIn) == 8) then

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

    call msg('int_hInterpScalar_gsv', 'END', verb_opt=4)
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
    real(4),                    intent(inout) :: fieldOut_r4(:,:)
    real(4),                    intent(in)    :: fieldIn_r4(:,:)
    character(len=*),           intent(in)    :: interpDegree
    character(len=*), optional, intent(in)    :: extrapDegree_opt
    ! Result:
    integer :: ierr

    ! Locals:
    integer :: ezsint

    call msg('int_hInterpScalar_r4_2d', 'START', verb_opt=4)
    ! read the namelist
    call int_readNml()

    ! do the standard interpolation
    call int_setezopt(interpDegree, extrapDegree_opt)   
    ierr = ezsint(fieldOut_r4,fieldIn_r4)

    call msg('int_hInterpScalar_r4_2d', 'END', verb_opt=4)
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
    ! Result:
    integer                         :: ierr

    ! Locals:
    integer :: gdxyfll, omp_get_thread_num
    integer :: niCloud, njCloud, niGrid, njGrid, myThreadNum
    integer :: top, bottom, left, right, numBoxIndexes, lonIndexCloud, latIndexCloud
    integer :: boxSize, lonBoxIndex, latBoxIndex, boxIndex, lonIndexGrid, latIndexGrid
    integer :: lonBoxIndexes(100), latBoxIndexes(100), ngp
    integer :: nfill(mmpi_numThread), nhole(mmpi_numThread), nextrap0, nextrap1
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
    call msg('int_sintCloudToGrid_gsv', 'START', verb_opt=2)

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
            call msg('int_sintCloudToGrid_gsv', &
                 'Expecting 0 or 1 for the mask field at grid point: ' &
                 //str(lonIndexGrid)//', '//str(latIndexGrid)//', ' &
                 //str(maskGrid(lonIndexGrid,latIndexGrid)))
            call utl_abort('int_sintCloudToGrid_gsv')
          end if

        end do
      end do
    end if

    if ( checkCloudToGridUnassigned ) then
      call msg('int_sintCloudToGrid_gsv', & 
             new_line('')//'Total number of grid points:                                   '//str(niGrid*njGrid) &
           //new_line('')//'Number of grid points not covered by the cloud of points:      '//str(sum(nhole(:))) &
           //new_line('')//'Number of grid points filled by neighbours:                    '//str(sum(nfill(:))) &
           //new_line('')//'Number of grid points with extrapolated value in masked area:  '//str(nextrap0) &
           //new_line('')//'Number of grid points with extrapolated value in visible area: '//str(nextrap1))

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

    call msg('int_sintCloudToGrid_gsv', 'END', verb_opt=2)
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
    real(8),                    intent(inout) :: fieldOut_r8(:,:)
    real(8),                    intent(in)    :: fieldIn_r8(:,:)
    character(len=*),           intent(in)    :: interpDegree
    character(len=*), optional, intent(in)    :: extrapDegree_opt
    ! Result:
    integer :: ierr

    ! Locals:
    integer :: nii, nji, nio, njo     
    integer :: jk1, jk2
    real(4), allocatable :: bufferi4(:,:), buffero4(:,:)
    integer :: ezsint

    call msg('int_hInterpScalar_r8_2d', 'START', verb_opt=4)
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

    call msg('int_hInterpScalar_r8_2d', 'END', verb_opt=4)
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
    type(struct_gsv),           intent(inout) :: stateVectorOut
    type(struct_gsv),           intent(inout) :: stateVectorIn
    character(len=*),           intent(in)    :: varName
    integer         ,           intent(in)    :: levIndex
    integer         ,           intent(in)    :: stepIndex
    character(len=*),           intent(in)    :: interpDegree
    character(len=*), optional, intent(in)    :: extrapDegree_opt
    ! Result:
    integer :: ierr

    ! Locals:
    real(4), pointer :: UUout4(:,:,:,:), VVout4(:,:,:,:), UUin4(:,:,:,:), VVin4(:,:,:,:)
    real(8), pointer :: UUout8(:,:,:,:), VVout8(:,:,:,:), UUin8(:,:,:,:), VVin8(:,:,:,:)
    real(4), pointer :: UVout4(:,:,:), UVin4(:,:,:)
    real(8), pointer :: UVout8(:,:,:), UVin8(:,:,:)
    integer :: ezuvint, ezdefset

    call msg('int_hInterpUV_gsv', 'START', verb_opt=4)

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

    if ( gsv_getDataKind(stateVectorOut) == 4 .and. gsv_getDataKind(stateVectorIn) == 4) then

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
        call utl_abort('int_hInterpUV_gsv: unexpected varName: '//trim(varName))
      end if

    else if ( gsv_getDataKind(stateVectorOut) == 8 .and. gsv_getDataKind(stateVectorIn) == 8) then

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
        call utl_abort('int_hInterpUV_gsv: unexpected varName: '//trim(varName))
      end if

      deallocate(UUin4,VVin4,UUout4,VVout4)

    else

      call utl_abort('int_hInterpUV_gsv: not implemented for mixed dataKind')

    end if

    call msg('int_hInterpUV_gsv', 'END', verb_opt=4)
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
    real(4),                    intent(inout) :: uuout(:,:)
    real(4),                    intent(inout) :: vvout(:,:)
    real(4),                    intent(in)    :: uuin(:,:)
    real(4),                    intent(in)    :: vvin(:,:)
    character(len=*),           intent(in)    :: interpDegree
    character(len=*), optional, intent(in)    :: extrapDegree_opt
    ! Result:
    integer :: ierr

    ! Locals:
    integer :: ezuvint

    call msg('int_hInterpUV_r4_2d', 'START', verb_opt=4)
    ! read the namelist
    call int_readNml()

    ! do the standard interpolation
    call int_setezopt(interpDegree, extrapDegree_opt)   
    ierr = ezuvint(uuout, vvout, uuin, vvin)

    call msg('int_hInterpUV_r4_2d', 'END', verb_opt=4)
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
    real(8),                    intent(inout) :: uuout(:,:)
    real(8),                    intent(inout) :: vvout(:,:)
    real(8),                    intent(in)    :: uuin(:,:)
    real(8),                    intent(in)    :: vvin(:,:)
    character(len=*),           intent(in)    :: interpDegree
    character(len=*), optional, intent(in)    :: extrapDegree_opt
    ! Result:
    integer :: ierr

    ! Locals:
    integer :: nio, njo, nii, nji
    integer :: jk1, jk2
    real, allocatable :: bufuuout4(:,:), bufvvout4(:,:)
    real, allocatable :: bufuuin4(:,:), bufvvin4(:,:)
    integer :: ezuvint

    call msg('int_hInterpUV_r8_2d', 'START', verb_opt=4)
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

    call msg('int_hInterpUV_r8_2d', 'END', verb_opt=4)
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
    integer,          intent(in) :: ni
    integer,          intent(in) :: nj
    integer,          intent(in) :: ig1
    integer,          intent(in) :: ig2
    integer,          intent(in) :: ig3
    integer,          intent(in) :: ig4
    real(8),          intent(in) :: ax(:)
    real(8),          intent(in) :: ay(:)
    character(len=*), intent(in) :: grtyp
    character(len=*), intent(in) :: grtypref
    ! Result:
    integer :: vezgdef

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
    integer,          intent(in) :: ig1
    integer,          intent(in) :: ig2
    integer,          intent(in) :: ig3
    integer,          intent(in) :: ig4   
    real(8),          intent(in) :: xlat0
    real(8),          intent(in) :: xlon0
    real(8),          intent(in) :: dlat
    real(8),          intent(in) :: dlon 
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

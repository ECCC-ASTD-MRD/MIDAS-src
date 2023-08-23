
module columnData_mod
  ! MODULE columnData_mod (prefix='col' category='6. High-level data objects')
  !
  !:Purpose:  A derived type and related procedures for storing and manipulating
  !           vertical columns of analysis variables on model or analysis grid
  !           levels. These columns are generally produced by horizontally
  !           interpolating a gridStateVector object to the observation
  !           locations.
  !
  use midasMpi_mod
  use varNameList_mod
  use verticalCoord_mod
  use mathPhysConstants_mod
  use utilities_mod

  implicit none
  save
  private

  ! public variables and types
  public :: col_rhumin, col_minValVarKindCH, struct_columnData

  ! public subroutines and functions
  public :: col_setup, col_allocate, col_deallocate
  public :: col_varExist, col_getOffsetFromVarno
  public :: col_getNumLev, col_getNumCol, col_getVarNameFromK
  public :: col_getPressure, col_getHeight, col_setHeightSfc
  public :: col_zero, col_getAllColumns, col_getColumn, col_getElem, col_getVco, col_setVco
  public :: col_getLevIndexFromVarLevIndex, col_add, col_copy

  type struct_columnData
    integer           :: nk, numCol
    logical           :: allocated=.false.
    real(8), pointer  :: all(:,:)
    real(8), pointer  :: heightSfc(:)
    real(8), pointer  :: oltv(:,:,:)    ! Tangent linear operator of virtual temperature
    integer, pointer  :: varOffset(:),varNumLev(:)
    logical           :: varExistList(vnl_numVarMax)
    type(struct_vco), pointer :: vco => null()
    logical           :: addHeightSfcOffset = .false.
    real(8), pointer  :: lat(:)
  end type struct_columnData

  real(8) :: col_rhumin
  logical :: varExistList(vnl_numvarmax)

  ! Minimum values for variables of CH kind
  real(8) :: col_minValVarKindCH(vnl_numVarMax)

  ! Namelist variables
  real(8) :: rhumin                         ! minimum humidity value imposed after interpolation to columns
  logical :: addHeightSfcOffset             ! choose to add non-zero height offset to diagnostic (sfc) levels
  real(8) :: minValVarKindCH(vnl_numVarMax) ! variable-dependent minimum value applied to chemistry variables

contains

  !--------------------------------------------------------------------------
  ! col_setup
  !--------------------------------------------------------------------------
  subroutine col_setup
    implicit none

    ! Locals:
    integer :: varIndex, loopIndex
    integer :: fnom,fclos,nulnam,ierr
    integer :: numVar3D, numVar2D, numVarOther

    ! Namelist variables (local)
    character(len=4) :: anlvar(vnl_numvarmax)           ! list of state variable names
    character(len=8) :: anltime_bin                     ! can be 'MIDDLE', 'FIRST' or 'LAST'
    logical          :: conversionVarKindCHtoMicrograms ! activate unit conversion for CH variables
    logical          :: abortOnMpiImbalance             ! choose to abort program when MPI imbalance is too large

    namelist /namstate/anlvar,rhumin,anltime_bin,addHeightSfcOffset,conversionVarKindCHtoMicrograms, &
                       minValVarKindCH, abortOnMpiImbalance

    if(mmpi_myid == 0) write(*,*) 'col_setup: List of known (valid) variable names'
    if(mmpi_myid == 0) write(*,*) 'col_setup: varNameList3D=',vnl_varNameList3D
    if(mmpi_myid == 0) write(*,*) 'col_setup: varNameList2D=',vnl_varNameList2D
    if(mmpi_myid == 0) write(*,*) 'col_setup: varNameList  =',vnl_varNameList

    ! Read NAMELIST NAMSTATE to find which fields are needed

    anlvar(:) = '    '
    rhumin = MPC_MINIMUM_HU_R8
    anltime_bin = 'MIDDLE'
    addHeightSfcOffset = .false.
    conversionVarKindCHtoMicrograms = .false.
    minValVarKindCH(:) = mpc_missingValue_r8
    abortOnMpiImbalance = .true.

    nulnam=0
    ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=namstate,iostat=ierr)
    if(ierr.ne.0) call utl_abort('col_setup: Error reading namelist')
    if(mmpi_myid == 0) write(*,nml=namstate)
    ierr=fclos(nulnam)

    col_rhumin = rhumin
    col_minValVarKindCH(:)=minValVarKindCH(:)

    if( varneed('Z_T') .or. varneed('Z_M') ) then
      call utl_abort('col_setup: height can not be specified as analysis variable in namelist!')
    end if
    if( varneed('P_T') .or. varneed('P_M') ) then
      call utl_abort('col_setup: pressure can not be specified as analysis variable in namelist!')
    end if

    numVar3D    = 0
    numVar2D    = 0
    numVarOther = 0

    do varIndex = 1, vnl_numvarmax3D
      if (varneed(vnl_varNameList3D(varIndex))) then
        varExistList(varIndex) = .true.
        numVar3D = numVar3D + 1
      else
        varExistList(varIndex) = .false.
      end if
    end do

    do varIndex = 1, vnl_numvarmax2D
      if (varneed(vnl_varNameList2D(varIndex))) then
        varExistList(varIndex+vnl_numvarmax3D) = .true.
        numVar2D = numVar2D + 1
      else
        varExistList(varIndex+vnl_numvarmax3D) = .false.
      end if
    end do

    do varIndex = 1, vnl_numvarmaxOther
      if (varneed(vnl_varNameListOther(varIndex))) then
        varExistList(varIndex+vnl_numvarmax3D+vnl_numvarmax2D) = .true.
        numVarOther = numVarOther + 1
      else
        varExistList(varIndex+vnl_numvarmax3D+vnl_numvarmax2D) = .false.
      end if
    end do

    ! Setup to assign min values to apply
    
    ! Check for input values only for variables of CH kind
    do varIndex = 1, vnl_numvarmax
      if ( trim(AnlVar(varIndex)) == '' ) exit
      if ( vnl_varKindFromVarname(AnlVar(varIndex)) == 'CH' ) then
        if ( minValVarKindCH(varIndex) < 0.99d0 * MPC_missingValue_R8 ) then
          if ( trim(AnlVar(varIndex)) == 'AF' .or. trim(AnlVar(varIndex)) == 'AC' ) then
            ! Set for particulate matter in micrograms/cm^3
            minValVarKindCH(varIndex) = MPC_MINIMUM_PM_R8
          else
            ! Set for concentrations in micrograms/kg
            minValVarKindCH(varIndex) = MPC_MINIMUM_CH_R8
          end if
        end if
      end if
    end do

    ! Assign min values to apply
    col_minValVarKindCH(:) = MPC_missingValue_R8
    do varIndex = 1, vnl_numvarmax
      if ( varExistList(varIndex) ) then
        do loopIndex = 1, vnl_numvarmax
          if ( trim(AnlVar(loopIndex)) == '' ) exit
          if ( trim(vnl_varNameList(varIndex)) == trim(AnlVar(loopIndex)) ) &
             col_minValVarKindCH(varIndex) = minValVarKindCH(loopIndex)
        end do
      end if 
    end do

    if(mmpi_myid == 0) write(*,*) 'col_setup: numVar3D (no Z_T/Z_M/P_T/P_M included), numVar2D, numVarOther = ', numVar3D, numVar2D, numVarOther
    if(mmpi_myid == 0) write(*,*) 'col_setup: varExistList (no Z_T/Z_M/P_T/P_M included) = ',varExistList

    contains

      logical function varneed(varName)
        character(len=*) :: varName
        integer :: jvar
 
        varneed = .false.
        do jvar = 1, vnl_numVarMax
          if (trim(varName) == trim(anlvar(jvar))) then
            varneed = .true.
          end if
        end do

      end function varneed

  end subroutine col_setup

  !--------------------------------------------------------------------------
  ! col_zero
  !--------------------------------------------------------------------------
  subroutine col_zero(column)
    implicit none

    ! Arguments:
    type(struct_columnData), intent(inout) :: column

    if (column%numCol > 0) then
      column%all(:,:) = 0.0d0
      column%heightSfc(:) = 0.0d0
    end if

  end subroutine col_zero

  !--------------------------------------------------------------------------
  ! col_allocate
  !--------------------------------------------------------------------------
  subroutine col_allocate(column, numCol, beSilent_opt, setToZero_opt, varNames_opt)
    implicit none

    ! Arguments:
    type(struct_columnData),    intent(inout) :: column
    integer,                    intent(in)    :: numCol
    logical,          optional, intent(in)    :: beSilent_opt
    logical,          optional, intent(in)    :: setToZero_opt
    character(len=*), optional, intent(in)    :: varNames_opt(:)

    ! Locals:
    integer :: iloc, varIndex, varIndex2, numVar
    logical :: beSilent, setToZero

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .false.
    end if

    if ( present(setToZero_opt) ) then
      setToZero = setToZero_opt
    else
      setToZero = .true.
    end if

    if ( present(varNames_opt) ) then      
      column%varExistList(:) = .false.
      numVar = size( varNames_opt ) 
      do varIndex2 = 1, numVar
        varIndex = vnl_varListIndex(varNames_opt(varIndex2))
        column%varExistList(varIndex) = .true.
      end do
    else
      ! set the variable list using the global ExistList
      column%varExistList(:) = varExistList(:)
    end if

    if ( column%varExistList(vnl_varListIndex('TT')) .and. &
         column%varExistList(vnl_varListIndex('HU')) .and. &
         column%varExistList(vnl_varListIndex('P0')) ) then
      if ( col_getNumLev(column,'TH') > 0 ) column%varExistList(vnl_varListIndex('Z_T')) = .true.
      if ( col_getNumLev(column,'MM') > 0 ) column%varExistList(vnl_varListIndex('Z_M')) = .true.
    end if

    if ( column%varExistList(vnl_varListIndex('P0')) ) then
      if ( col_getNumLev(column,'TH') > 0 ) column%varExistList(vnl_varListIndex('P_T')) = .true.
      if ( col_getNumLev(column,'MM') > 0 ) column%varExistList(vnl_varListIndex('P_M')) = .true.
    end if

    ! add P0LS to the varExistList if vcode=5100
    if (column%vco%vcode == 5100) then
      column%varExistList(vnl_varListIndex('P0LS')) = .true.
    end if

    column%numCol = numCol

    if(.not.column%vco%initialized) then
      call utl_abort('col_allocate: VerticalCoord has not been initialized!')
    end if

    allocate(column%varOffset(vnl_numvarmax))
    column%varOffset(:)=0
    allocate(column%varNumLev(vnl_numvarmax))
    column%varNumLev(:)=0

    iloc = 0
    do varIndex = 1, vnl_numvarmax3d
      if(column%varExistList(varIndex)) then
        column%varOffset(varIndex) = iloc
        column%varNumLev(varIndex) = col_getNumLev(column,vnl_varLevelFromVarname(vnl_varNameList(varIndex)))
        iloc = iloc + column%varNumLev(varIndex)
      end if
    end do
    do varIndex2 = 1, vnl_numvarmax2d
      varIndex = varIndex2+vnl_numvarmax3d
      if(column%varExistList(varIndex)) then
        column%varOffset(varIndex) = iloc
        column%varNumLev(varIndex) = 1
        iloc = iloc + 1
      end if
    end do
    do varIndex2 = 1, vnl_numvarmaxOther
      varIndex = varIndex2+vnl_numvarmax3d+vnl_numvarmax2d
      if(column%varExistList(varIndex)) then
        column%varOffset(varIndex) = iloc
        column%varNumLev(varIndex) = col_getNumLev(column,'OT',vnl_varNameListOther(varIndex2))
        iloc = iloc + column%varNumLev(varIndex)
      end if
    end do

    if (iloc == 0) then
      call utl_abort('col_allocate: Nothing to allocate')
    end if

    column%nk = iloc

    if(column%numCol.le.0) then
      if ( .not.beSilent ) write(*,*) 'col_allocate: number of columns is zero, not allocated'
    else         
      allocate(column%all(column%nk,column%numCol))
      if ( setToZero ) column%all(:,:)=0.0d0

      allocate(column%heightSfc(column%numCol))
      column%heightSfc(:)=0.0d0

      allocate(column%oltv(2,col_getNumLev(column,'TH'),numCol))
      if ( setToZero ) column%oltv(:,:,:)=0.0d0

      allocate(column%lat(numCol))
      if ( setToZero ) column%lat(:)=0.0d0
    end if
 
    if(mmpi_myid == 0 .and. .not.beSilent) write(*,*) 'col_allocate: column%nk = ', column%nk
    if(mmpi_myid == 0 .and. .not.beSilent) write(*,*) 'col_allocate: varOffset=',column%varOffset
    if(mmpi_myid == 0 .and. .not.beSilent) write(*,*) 'col_allocate: varNumLev=',column%varNumLev

    column%addHeightSfcOffset = addHeightSfcOffset

    column%allocated=.true.

  end subroutine col_allocate

  !--------------------------------------------------------------------------
  ! col_deallocate
  !--------------------------------------------------------------------------
  subroutine col_deallocate(column)
    implicit none

    ! Arguments:
    type(struct_columnData), intent(inout) :: column

    deallocate(column%varOffset)
    deallocate(column%varNumLev)

    if(column%numCol.gt.0) then
      deallocate(column%all)
      deallocate(column%heightSfc)
      deallocate(column%oltv)
    end if

    column%allocated=.false.

  end subroutine col_deallocate

  !--------------------------------------------------------------------------
  ! col_varExist
  !--------------------------------------------------------------------------
  recursive function col_varExist(column_opt,varName) result(varExist)
    implicit none

    ! Arguments:
    type(struct_columnData), optional, intent(in) :: column_opt
    character(len=*),                  intent(in) :: varName
    ! Result:
    logical                                       :: varExist 

    if (varName == 'Z_*') then
      varExist =  col_varExist(column_opt, 'Z_T') .and. &
                  col_varExist(column_opt, 'Z_M')
    else if (varName == 'P_*') then
      varExist =  col_varExist(column_opt, 'P_T') .and. &
                  col_varExist(column_opt, 'P_M')
    else
      if ( present(column_opt) ) then
        if ( column_opt%varExistList(vnl_varListIndex(varName)) ) then
          varExist = .true.
        else
          varExist = .false.
        end if
      else
        if ( varExistList(vnl_varListIndex(varName)) ) then
          varExist = .true.
        else
          varExist = .false.
        end if
      end if

      if (present(column_opt)) then
        varExist = column_opt % varExistList(vnl_varListIndex(varName))
      else
        varExist = varExistList(vnl_varListIndex(varName))
      end if
    end if
  
  end function col_varExist

  !--------------------------------------------------------------------------
  ! col_getOffsetFromVarno
  !--------------------------------------------------------------------------
  function col_getOffsetFromVarno(column,varnum,varNumberChm_opt,modelName_opt) result(offset)
    implicit none

    ! Arguments:
    type(struct_columnData),    intent(in) :: column
    integer,                    intent(in) :: varnum
    integer,          optional, intent(in) :: varNumberChm_opt
    character(len=*), optional, intent(in) :: modelName_opt
    ! Result:
    integer                 :: offset

    offset=column%varOffset(vnl_varListIndex(vnl_varnameFromVarnum(varnum,varNumberChm_opt=varNumberChm_opt,modelName_opt=modelName_opt)))

  end function col_getOffsetFromVarno

  !--------------------------------------------------------------------------
  ! col_getLevIndexFromVarLevIndex
  !--------------------------------------------------------------------------
  function col_getLevIndexFromVarLevIndex(column, varLevIndex) result(levIndex)
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in) :: column
    integer,                 intent(in) :: varLevIndex
    ! Result:
    integer                 :: varIndex

    ! Locals:
    integer                 :: levIndex

    do varIndex = 1, vnl_numvarmax
      if ( column%varExistList(varIndex) ) then
        if ( (varLevIndex >= (column%varOffset(varIndex) + 1)) .and.  &
            (varLevIndex <= (column%varOffset(varIndex) + column%varNumLev(varIndex))) ) then
          levIndex = varLevIndex - column%varOffset(varIndex)
          return
        end if
      end if
    end do

    write(*,*) 'col_getLevIndexFromVarLevIndex: varLevIndex out of range: ', varLevIndex
    call utl_abort('col_getLevIndexFromVarLevIndex')

  end function col_getLevIndexFromVarLevIndex

  !--------------------------------------------------------------------------
  ! col_getVarNameFromK
  !--------------------------------------------------------------------------
  function col_getVarNameFromK(column,kIndex) result(varName)
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in) :: column
    integer,                 intent(in) :: kIndex
    ! Result:
    character(len=4)    :: varName

    ! Locals:
    integer             :: varIndex

    do varIndex = 1, vnl_numvarmax
      if ( column%varExistList(varIndex) ) then
        if ( (kIndex >= (column%varOffset(varIndex) + 1)) .and.  &
            (kIndex <= (column%varOffset(varIndex) + column%varNumLev(varIndex))) ) then
          varName = vnl_varNameList(varIndex)
          return
        end if
      end if
    end do

    write(*,*) 'col_getVarNameFromK: kIndex out of range: ', kIndex
    call utl_abort('col_getVarNameFromK')

  end function col_getVarNameFromK

  !--------------------------------------------------------------------------
  ! col_getPressure
  !--------------------------------------------------------------------------
  function col_getPressure(column,ilev,headerIndex,varLevel) result(pressure)
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in) :: column
    integer,                 intent(in) :: ilev
    integer,                 intent(in) :: headerIndex
    character(len=*),        intent(in) :: varLevel
    ! Result:
    real(8)                             :: pressure

    ! Locals:
    integer                             :: ilev1

    if (varLevel == 'TH' .and. col_varExist(column,'P_T')) then
      ilev1 = 1 + column%varOffset(vnl_varListIndex('P_T'))
      pressure = column%all(ilev1+ilev-1,headerIndex)
    elseif (varLevel == 'MM' .and. col_varExist(column,'P_M') ) then
      ilev1 = 1 + column%varOffset(vnl_varListIndex('P_M'))
      pressure = column%all(ilev1+ilev-1,headerIndex)
    else
      call utl_abort('col_getPressure: Unknown variable type: ' // varLevel)
    end if

  end function col_getPressure
 
  !--------------------------------------------------------------------------
  ! col_getHeight
  !--------------------------------------------------------------------------
  function col_getHeight(column,ilev,headerIndex,varLevel) result(height)
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in) :: column
    integer,                 intent(in) :: ilev
    integer,                 intent(in) :: headerIndex
    character(len=*),        intent(in) :: varLevel
    ! Result:
    real(8)                             :: height

    ! Locals:
    integer                             :: ilev1

    if (varLevel == 'TH') then
      if (.not. col_varExist(column,'Z_T') ) then
        call utl_abort('col_getHeight: Z_T not found!')
      end if
      ilev1 = 1 + column%varOffset(vnl_varListIndex('Z_T'))
      height = column%all(ilev1+ilev-1,headerIndex)
    else if (varLevel == 'MM') then
      if (.not. col_varExist(column,'Z_M') ) then
        call utl_abort('col_getHeight: Z_M not found!')
      end if 
      ilev1 = 1 + column%varOffset(vnl_varListIndex('Z_M'))
      height = column%all(ilev1+ilev-1,headerIndex)
    else if (varLevel == 'SF' ) then
      height = column%heightSfc(headerIndex)
    else
      call utl_abort('col_getHeight: unknown varLevel! ' // varLevel)
    end if

  end function col_getHeight

  !--------------------------------------------------------------------------
  ! col_setHeightsSfc
  !--------------------------------------------------------------------------
  subroutine col_setHeightSfc(column,headerIndex,height)
    implicit none

    ! Arguments:
    type(struct_columnData), intent(inout) :: column
    integer,                 intent(in)    :: headerIndex
    real(8),                 intent(in)    :: height

    column%heightSfc(headerIndex) = height

  end subroutine col_setHeightSfc

  !--------------------------------------------------------------------------
  ! col_getAllColumns
  !--------------------------------------------------------------------------
  function col_getAllColumns(column,varName_opt) result(allColumns)
    implicit none

    ! Arguments:
    type(struct_columnData),    intent(in) :: column
    character(len=*), optional, intent(in) :: varName_opt
    ! Result:
    real(8), pointer                       :: allColumns(:,:)

    ! Locals:
    integer                                :: ilev1,ilev2

    if ( column%numCol > 0 ) then
      if(present(varName_opt)) then
        if ( col_varExist(column,varName_opt) ) then
          ilev1 = column%varOffset(vnl_varListIndex(varName_opt))+1
          ilev2 = ilev1 - 1 + column%varNumLev(vnl_varListIndex(varName_opt))
          allColumns => column%all(ilev1:ilev2,:)
        else
          call utl_abort('col_getAllColumns: Unknown variable name! ' // varName_opt)
        end if
      else
        allColumns => column%all(:,:)
      end if
    else
      allColumns => null()
    end if

  end function col_getAllColumns

  !--------------------------------------------------------------------------
  ! col_getColumn
  !--------------------------------------------------------------------------
  function col_getColumn(column,headerIndex,varName_opt) result(onecolumn)
    implicit none

    ! Arguments:
    type(struct_columnData),    intent(in) :: column
    integer,                    intent(in) :: headerIndex
    character(len=*), optional, intent(in) :: varName_opt
    ! Result:
    real(8), pointer                       :: onecolumn(:)

    ! Locals:
    integer                                :: ilev1,ilev2

    if(present(varName_opt)) then
      if(col_varExist(column,varName_opt)) then
        ilev1 = column%varOffset(vnl_varListIndex(varName_opt))+1
        ilev2 = ilev1 - 1 + column%varNumLev(vnl_varListIndex(varName_opt))
        onecolumn => column%all(ilev1:ilev2,headerIndex)
      else
        call utl_abort('col_getColumn: Unknown variable name! ' // varName_opt)
      end if
    else
      onecolumn => column%all(:,headerIndex)
    end if

  end function col_getColumn

  !--------------------------------------------------------------------------
  ! col_getElem
  !--------------------------------------------------------------------------
  function col_getElem(column,ilev,headerIndex,varName_opt) result(value)
    implicit none

    ! Arguments:
    type(struct_columnData),    intent(in) :: column
    integer,                    intent(in) :: ilev
    integer,                    intent(in) :: headerIndex
    character(len=*), optional, intent(in) :: varName_opt
    ! Result:
    real(8)                                :: value

    if(present(varName_opt)) then
      if(.not.col_varExist(column,varName_opt)) call utl_Abort('col_getElem: Unknown variable name! ' // varName_opt)
      value = column%all(column%varOffset(vnl_varListIndex(varName_opt))+ilev,headerIndex)
    else
      value = column%all(ilev,headerIndex)
    end if

  end function col_getElem

  !--------------------------------------------------------------------------
  ! col_getNumLev
  !--------------------------------------------------------------------------
  function col_getNumLev(column,varLevel,varName_opt) result(nlev)
    implicit none

    ! Arguments:
    type(struct_columnData),    intent(in) :: column
    character(len=*),           intent(in) :: varLevel
    character(len=*), optional, intent(in) :: varName_opt
    ! Result:
    integer                             :: nlev

    nlev = vco_getNumLev(column%vco,varLevel,varName_opt)

  end function col_getNumLev

  !--------------------------------------------------------------------------
  ! col_getNumCol
  !--------------------------------------------------------------------------
  function col_getNumCol(column) result(numColumn)
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in) :: column
    ! Result:
    integer                             :: numColumn

    numColumn = column%numCol

  end function col_getNumCol

  !--------------------------------------------------------------------------
  ! col_getVco
  !--------------------------------------------------------------------------
  function col_getVco(column) result(vco_ptr)
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in) :: column
    ! Result:
    type(struct_vco), pointer :: vco_ptr

    vco_ptr => column%vco

  end function col_getVco

  !--------------------------------------------------------------------------
  ! col_setVco
  !--------------------------------------------------------------------------
  subroutine col_setVco(column,vco_ptr)
    implicit none

    ! Arguments:
    type(struct_columnData),   intent(inout) :: column
    type(struct_vco), pointer, intent(in)    :: vco_ptr

    column%vco => vco_ptr

  end subroutine col_setVco

  !--------------------------------------------------------------------------
  ! col_add
  !--------------------------------------------------------------------------
  subroutine col_add(columnIn, columnInout, scaleFactor_opt)
    !
    ! :Purpose: Adds two columns
    !           columnInout = columnInout + scaleFactor_opt * columnIn
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)     :: columnIn           ! first operand 
    type(struct_columnData), intent(inout)  :: columnInout        ! second operand, will receive the result
    real(8), optional,       intent(in)     :: scaleFactor_opt    ! optional scaling of the second operand prior to the addition

    ! Locals:
    real(8), pointer                  :: ptrColInOut(:,:)
    real(8), pointer                  :: ptrColIn(:,:)

    if (columnInout%nk /= columnIn%nk) &
                call utl_abort('col_add: Number of levels between two column object are not the same')

    if (columnInout%numCol /= columnIn%numCol) &
                call utl_abort('col_add: Number of columns between two column object are not the same')

    ptrColInOut => col_getAllColumns(columnInout)
    ptrColIn => col_getAllColumns(columnIn)

    if (present(scaleFactor_opt)) then
      ptrColInOut(:,:) = ptrColInOut(:,:) + scaleFactor_opt * ptrColIn(:,:)
    else
      ptrColInOut(:,:) = ptrColInOut(:,:) + ptrColIn(:,:)
    end if

  end subroutine col_add

  !--------------------------------------------------------------------------
  ! col_copy
  !--------------------------------------------------------------------------
  subroutine col_copy(columnIn, columnOut)
    !
    ! :Purpose: Copy column object from columnIn to columnOut
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)   :: columnIn  ! Source column to be copied from
    type(struct_columnData), intent(inout)  :: columnOut ! Destination column to be copied into

    if (.not. columnIn%allocated) &
          call utl_abort('col_copy: columnIn is not allocated')

    if (.not. columnOut%allocated) &
          call utl_abort('col_copy: columnOut is not allocated')

    if (any(columnIn%varNumLev(:) /= columnOut%varNumLev(:))) &
          call utl_abort('col_copy: varNumLev in columnIn and columnOut are not equal')
   
    if (.not. vco_equal(col_getVco(columnIn), col_getVco(columnOut))) &
          call utl_abort('col_copy: Vco in columnIn and columnOut are not equal')
    
    !Copy Content
    columnOut%nk = columnIn%nk
    columnOut%numCol = columnIn%numCol
    columnOut%varExistList = columnIn%varExistList
    columnOut%addHeightSfcOffset = columnIn%addHeightSfcOffset
    columnOut%varExistList = columnIn%varExistList
    columnOut%all(:,:) =  columnIn%all(:,:)
    columnOut%heightSfc(:) = columnIn%heightSfc(:)
    columnOut%oltv(:,:,:) = columnIn%oltv(:,:,:)
    columnOut%lat(:) = columnIn%lat(:)

  end subroutine col_copy

end module columnData_mod

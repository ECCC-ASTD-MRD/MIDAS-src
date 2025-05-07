
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

  ! Public variables and types
  real(8), public, protected :: col_rhumin
  real(8), public, protected :: col_minValVarKindCH(vnl_numVarMax) ! Minimum values for variables of CH kind

  ! Public dervied type
  public :: struct_columnData

  ! Public subroutines and functions
  public :: col_setup, col_allocate, col_deallocate
  public :: col_varExist, col_getOffsetFromVarno, col_getOffsetFromVarName
  public :: col_getNumLev, col_getNumCol, col_getNumVarLev, col_getVarNameFromVarLev
  public :: col_addHeightSfcOffset
  public :: col_getPressure, col_getHeight, col_setHeightSfc, col_copyHeightSfc
  public :: col_zero, col_getAllColumns, col_getColumn, col_getElem
  public :: col_getLat, col_setLat, col_getOltv, col_setOltv
  public :: col_getVco, col_setVco
  public :: col_getLevFromVarLev, col_add, col_copy, col_copyLat

  type struct_columnData
    private
    integer                   :: numVarLev
    integer                   :: numCol
    logical                   :: allocated=.false.
    logical                   :: addHeightSfcOffset = .false.
    real(8),          pointer :: all(:,:)
    real(8),          pointer :: heightSfc(:)
    real(8),          pointer :: oltv(:,:,:)    ! Tangent linear operator of virtual temperature
    integer,          pointer :: varOffset(:)
    integer,          pointer :: varNumLev(:)
    logical                   :: varExistList(vnl_numVarMax)
    type(struct_vco), pointer :: vco => null()
    real(8),          pointer :: lat(:)
  end type struct_columnData

  logical :: varExistList(vnl_numvarmax)

  ! Namelist variables
  real(8) :: rhumin                         ! minimum humidity value imposed after interpolation to columns
  logical :: addHeightSfcOffset             ! choose to add non-zero height offset to diagnostic (sfc) levels
  real(8) :: minValVarKindCH(vnl_numVarMax) ! variable-dependent minimum value applied to chemistry variables

contains

  !--------------------------------------------------------------------------
  ! col_setup
  !--------------------------------------------------------------------------
  subroutine col_setup()
    !
    !:Purpose: Read the namelist and setup some module variables.
    !
    implicit none

    ! Locals:
    integer :: varIndex, loopIndex
    integer :: ierr
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

    call utl_tmg_start(181,'low-level--readNML')
    read(utl_flnml,nml=namstate,iostat=ierr)
    if(ierr /= 0) call utl_abort('col_setup: Error reading namelist')
    if(mmpi_myid == 0) write(*,nml=namstate)
    call utl_tmg_stop(181)

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
        !
        !:Purpose: Detemine if the given variable name is to be included for allocation.
        !
        implicit none

        character(len=*) :: varName
        integer :: jvar
 
        varneed = .false.
        NEED_LOOP: do jvar = 1, vnl_numVarMax
          if (trim(varName) == trim(anlvar(jvar))) then
            varneed = .true.
            exit NEED_LOOP
          end if
        end do NEED_LOOP

      end function varneed

  end subroutine col_setup

  !--------------------------------------------------------------------------
  ! col_zero
  !--------------------------------------------------------------------------
  subroutine col_zero(column)
    !
    !:Purpose: Set the data column components (variables and surface
    !          height) to zero.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(inout) :: column ! The `columnData` object

    if (column%numCol > 0) then
      column%all(:,:) = 0.0d0
    end if

  end subroutine col_zero

  !--------------------------------------------------------------------------
  ! col_allocate
  !--------------------------------------------------------------------------
  subroutine col_allocate(column, numCol, beSilent_opt, setToZero_opt, varNames_opt)
    !
    !:Purpose: Allocate contents of the `column` object. The user can either
    !          specify the list of variables to allocate or rely on the defaults
    !          as specified by the namelist `namstate`.
    !
    implicit none

    ! Arguments:
    type(struct_columnData),    intent(inout) :: column          ! The `columnData` object
    integer,                    intent(in)    :: numCol          ! Number of columns (header indexes)
    logical,          optional, intent(in)    :: beSilent_opt    ! Control over verbose output
    logical,          optional, intent(in)    :: setToZero_opt   ! Control to set contents to zero after allocating
    character(len=*), optional, intent(in)    :: varNames_opt(:) ! List of variable names to use for allocation

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

    ! add MELS to the varExistList if vcode=21001 and SLEVE is active
    if (column%vco%vcode == 21001 .and. column%vco%sleveCoord) then
      column%varExistList(vnl_varListIndex('MELS')) = .true.
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
        if (column%varNumLev(varIndex) <= 0) then
          call utl_abort('col_allocate: Number of levels is invalid for varName = ' // &
                         trim(vnl_varNameList(varIndex)))
        end if
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
        if (column%varNumLev(varIndex) <= 0) then
          call utl_abort('col_allocate: Number of levels is invalid for varName = ' // &
                         trim(vnl_varNameListOther(varIndex2)))
        end if
      end if
    end do

    if (iloc == 0) then
      call utl_abort('col_allocate: Nothing to allocate')
    end if

    column%numVarLev = iloc

    if(column%numCol.le.0) then
      if ( .not.beSilent ) write(*,*) 'col_allocate: number of columns is zero, not allocated'
    else         
      allocate(column%all(column%numVarLev,column%numCol))
      if ( setToZero ) column%all(:,:)=0.0d0

      allocate(column%heightSfc(column%numCol))
      column%heightSfc(:)=0.0d0

      allocate(column%oltv(2,col_getNumLev(column,'TH'),numCol))
      if ( setToZero ) column%oltv(:,:,:)=0.0d0

      allocate(column%lat(numCol))
      if ( setToZero ) column%lat(:)=0.0d0
    end if
 
    if(mmpi_myid == 0 .and. .not.beSilent) write(*,*) 'col_allocate: column%numVarLev = ', column%numVarLev
    if(mmpi_myid == 0 .and. .not.beSilent) write(*,*) 'col_allocate: varOffset=',column%varOffset
    if(mmpi_myid == 0 .and. .not.beSilent) write(*,*) 'col_allocate: varNumLev=',column%varNumLev

    column%addHeightSfcOffset = addHeightSfcOffset

    column%allocated=.true.

  end subroutine col_allocate

  !--------------------------------------------------------------------------
  ! col_deallocate
  !--------------------------------------------------------------------------
  subroutine col_deallocate(column)
    !
    !:Purpose: Deallocate the data contents of the `column` object.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(inout) :: column ! The `columnData` object

    deallocate(column%varOffset)
    deallocate(column%varNumLev)

    if(column%numCol > 0) then
      deallocate(column%all)
      deallocate(column%heightSfc)
      deallocate(column%oltv)
      deallocate(column%lat)
    end if

    column%allocated=.false.

  end subroutine col_deallocate

  !--------------------------------------------------------------------------
  ! col_varExist
  !--------------------------------------------------------------------------
  recursive function col_varExist(column_opt,varName) result(varExist)
    !
    !:Purpose: Check if the supplied variable name is part of this
    !          `column` object. Wildcard names are possible for height
    !          and pressure to include checks on these variables on
    !          both momentum and thermodynamic levels (`Z_*`, `P_*`).
    !
    implicit none

    ! Arguments:
    type(struct_columnData), optional, intent(in) :: column_opt ! The `columnData` object
    character(len=*),                  intent(in) :: varName    ! The variable name
    ! Result:
    logical                                       :: varExist   ! Indicates if the variable exists in object

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
    !
    !:Purpose: Return the "offset" for a given observation variable "number"
    !          (that is, the bufr code for an observation element) within the
    !          "varsLevs" list of variables and levels.
    !
    implicit none

    ! Arguments:
    type(struct_columnData),    intent(in) :: column           ! The `columnData` object
    integer,                    intent(in) :: varnum           ! Observation variable "number" (bufr code)
    integer,          optional, intent(in) :: varNumberChm_opt ! More specific obs code for chemistry obs
    character(len=*), optional, intent(in) :: modelName_opt    ! The "model name", only needed for some obs types.
    ! Result:
    integer                                :: offset           ! The returned offset value within "varsLevs"

    offset=column%varOffset(vnl_varListIndex(vnl_varnameFromVarnum(varnum,varNumberChm_opt=varNumberChm_opt,modelName_opt=modelName_opt)))

  end function col_getOffsetFromVarno

  !--------------------------------------------------------------------------
  ! col_getOffsetFromVarName
  !--------------------------------------------------------------------------
  function col_getOffsetFromVarName(column,varName) result(offset)
    !
    !:Purpose: Return the "offset" for a given variable name within the
    !          "varsLevs" list of variables and levels.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in) :: column   ! The `columnData` object
    character(len=*),        intent(in) :: varName  ! Variable name
    ! Result:
    integer                             :: offset   ! The returned offset value within "varsLevs"

    offset=column%varOffset(vnl_varListIndex(trim(varName)))

  end function col_getOffsetFromVarName

  !--------------------------------------------------------------------------
  ! col_getLevFromVarLev
  !--------------------------------------------------------------------------
  function col_getLevFromVarLev(column, varLevIndex) result(levIndex)
    !
    !:Purpose: Return the level index for a given value of the "varsLevs" index.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in) :: column      ! The `columnData` object
    integer,                 intent(in) :: varLevIndex ! The index into the "varsLevs" array
    ! Result:
    integer                             :: levIndex    ! The returned level index

    ! Locals:
    integer :: varIndex

    do varIndex = 1, vnl_numvarmax
      if ( column%varExistList(varIndex) ) then
        if ( (varLevIndex >= (column%varOffset(varIndex) + 1)) .and.  &
            (varLevIndex <= (column%varOffset(varIndex) + column%varNumLev(varIndex))) ) then
          levIndex = varLevIndex - column%varOffset(varIndex)
          return
        end if
      end if
    end do

    write(*,*) 'col_getLevFromVarLev: varLevIndex out of range: ', varLevIndex
    call utl_abort('col_getLevFromVarLev')

  end function col_getLevFromVarLev

  !--------------------------------------------------------------------------
  ! col_getVarNameFromVarLev
  !--------------------------------------------------------------------------
  function col_getVarNameFromVarLev(column,varLevIndex) result(varName)
    !
    !:Purpose: Return the variable name for a given value of the
    !          "varsLevs" index.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in) :: column  ! The `columnData` object
    integer,                 intent(in) :: varLevIndex  ! The index into "varsLevs" array
    ! Result:
    character(len=4)                    :: varName ! The returned variable name

    ! Locals:
    integer             :: varIndex

    do varIndex = 1, vnl_numvarmax
      if ( column%varExistList(varIndex) ) then
        if ( (varLevIndex >= (column%varOffset(varIndex) + 1)) .and.  &
            (varLevIndex <= (column%varOffset(varIndex) + column%varNumLev(varIndex))) ) then
          varName = vnl_varNameList(varIndex)
          return
        end if
      end if
    end do

    write(*,*) 'col_getVarNameFromVarLev: varLevIndex out of range: ', varLevIndex
    call utl_abort('col_getVarNameFromVarLev')

  end function col_getVarNameFromVarLev

  !--------------------------------------------------------------------------
  ! col_getPressure
  !--------------------------------------------------------------------------
  function col_getPressure(column,ilev,headerIndex,varLevel) result(pressure)
    !
    !:Purpose: Return the pressure for a given level index, header/column
    !          index and type of levels.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in) :: column      ! The `columnData` object
    integer,                 intent(in) :: ilev        ! The level index
    integer,                 intent(in) :: headerIndex ! The column/header index
    character(len=*),        intent(in) :: varLevel    ! The type of vertical level
    ! Result:
    real(8)                             :: pressure    ! The returned pressure value

    ! Locals:
    integer                             :: ilev1

    if (headerIndex > column%numCol .or. headerIndex < 1) then
      write(*,*) 'headerIndex = ', headerIndex
      call utl_abort('col_getPressure: headerIndex out of range')
    end if

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
    !
    !:Purpose: Return the height for a given level index, header/column
    !          index and type of levels.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in) :: column      ! The `columnData` object
    integer,                 intent(in) :: ilev        ! The level index
    integer,                 intent(in) :: headerIndex ! The column/header index
    character(len=*),        intent(in) :: varLevel    ! The type of vertical level
    ! Result:
    real(8)                             :: height      ! The returned height value

    ! Locals:
    integer                             :: ilev1

    if (headerIndex > column%numCol .or. headerIndex < 1) then
      write(*,*) 'headerIndex = ', headerIndex
      call utl_abort('col_getHeight: headerIndex out of range')
    end if

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
    !
    !:Purpose: Set the height of the surface for a given header/column index.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(inout) :: column      ! The `columnData` object
    integer,                 intent(in)    :: headerIndex ! The column/header index
    real(8),                 intent(in)    :: height      ! The supplied height value

    if (headerIndex > column%numCol .or. headerIndex < 1) then
      write(*,*) 'headerIndex = ', headerIndex
      call utl_abort('col_setHeightSfc: headerIndex out of range')
    end if

    column%heightSfc(headerIndex) = height

  end subroutine col_setHeightSfc

  !--------------------------------------------------------------------------
  ! col_getOltv
  !--------------------------------------------------------------------------
  function col_getOltv(column, varIndex, levIndex, headerIndex) result(value)
    !
    !:Purpose: Get the value of `oltv` for a given level index and
    !          header/column index. The second argument selects between
    !          either temperature or humidity.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in) :: column      ! The `columnData` object
    integer,                 intent(in) :: varIndex    ! Indicates temperature or humidity (1 or 2)
    integer,                 intent(in) :: levIndex    ! The level index
    integer,                 intent(in) :: headerIndex ! The column/header index
    ! Result:
    real(8)                             :: value       ! The returned value of "oltv"

    if (headerIndex > column%numCol .or. headerIndex < 1) then
      write(*,*) 'headerIndex = ', headerIndex
      call utl_abort('col_getOltv: headerIndex out of range')
    end if

    value = column%oltv(varIndex, levIndex, headerIndex)

  end function col_getOltv
  
  !--------------------------------------------------------------------------
  ! col_setOltv
  !--------------------------------------------------------------------------
  subroutine col_setOltv(column, varIndex, levIndex, headerIndex, value)
    !
    !:Purpose: Set the value of `oltv` for a given level index and
    !          header/column index. The second argument selects between
    !          either temperature or humidity.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(inout) :: column      ! The `columnData` object
    integer,                 intent(in)    :: varIndex    ! Indicates temperature or humidity (1 or 2)
    integer,                 intent(in)    :: levIndex    ! The level index
    integer,                 intent(in)    :: headerIndex ! The column/header index
    real(8),                 intent(in)    :: value       ! The value of "oltv" for setting

    if (headerIndex > column%numCol .or. headerIndex < 1) then
      write(*,*) 'headerIndex = ', headerIndex
      call utl_abort('col_setOltv: headerIndex out of range')
    end if

    column%oltv(varIndex, levIndex, headerIndex) = value

  end subroutine col_setOltv
  
  !--------------------------------------------------------------------------
  ! col_getAllColumns
  !--------------------------------------------------------------------------
  function col_getAllColumns(column,varName_opt) result(allColumns)
    !
    !:Purpose: Return a pointer to either a portion of the main `column`
    !          object data array for a given variable or to the entire array.
    !
    implicit none

    ! Arguments:
    type(struct_columnData),    intent(in) :: column      ! The `columnData` object
    character(len=*), optional, intent(in) :: varName_opt ! The variable name
    ! Result:
    real(8), pointer                       :: allColumns(:,:) ! Resulting pointer to complete array

    ! Locals:
    integer :: ilev1, ilev2

    if ( column%numCol > 0 ) then
      if (present(varName_opt)) then
        if ( col_varExist(column,varName_opt) ) then
          ilev1 = column%varOffset(vnl_varListIndex(varName_opt))+1
          ilev2 = ilev1 - 1 + column%varNumLev(vnl_varListIndex(varName_opt))
          allColumns => column%all(ilev1:ilev2,:)
        else
          call utl_abort('col_getAllColumns: Unknown variable name! ' // varName_opt)
        end if
      else
        ! No variable name specified, return full columns
        ilev1 = 1
        ilev2 = column%numVarLev
        allColumns => column%all(ilev1:ilev2,:)
      end if
    else
      allColumns => null()
    end if

  end function col_getAllColumns

  !--------------------------------------------------------------------------
  ! col_getColumn
  !--------------------------------------------------------------------------
  function col_getColumn(column,headerIndex,varName_opt) result(onecolumn)
    !
    !:Purpose: For a single header/column index, return a pointer to either
    !          a portion of the main `column` object data array for a given
    !          variable or to the entire array.
    !
    implicit none

    ! Arguments:
    type(struct_columnData),    intent(in) :: column       ! The `columnData` object
    integer,                    intent(in) :: headerIndex  ! The column/header index
    character(len=*), optional, intent(in) :: varName_opt  ! The variable name
    ! Result:
    real(8), pointer                       :: onecolumn(:) ! The return array pointer

    ! Locals:
    integer                                :: ilev1,ilev2

    if (headerIndex > column%numCol .or. headerIndex < 1) then
      write(*,*) 'headerIndex = ', headerIndex
      call utl_abort('col_getColumn: headerIndex out of range')
    end if

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
    !
    !:Purpose: Return a single value from the main `column` data array for
    !          a given level index, header/column index and variable name.
    !          If the variable name is not given then the level index is
    !          actually used as the index into the complete "varsLevs"
    !          array of all variables and levels.
    !
    implicit none

    ! Arguments:
    type(struct_columnData),    intent(in) :: column      ! The `columnData` object
    integer,                    intent(in) :: ilev        ! The level index
    integer,                    intent(in) :: headerIndex ! The column/header index
    character(len=*), optional, intent(in) :: varName_opt ! The variable name
    ! Result:
    real(8)                                :: value       ! The returned value

    if (headerIndex > column%numCol .or. headerIndex < 1) then
      write(*,*) 'headerIndex = ', headerIndex
      call utl_abort('col_getElem: headerIndex out of range')
    end if

    if(present(varName_opt)) then
      if(.not.col_varExist(column,varName_opt)) then
        call utl_Abort('col_getElem: Unknown variable name! ' // varName_opt)
      end if
      value = column%all(column%varOffset(vnl_varListIndex(varName_opt))+ilev,headerIndex)
    else
      if (ilev > column%numVarLev .or. ilev < 1) then
        write(*,*) 'varsLevs index = ', ilev
        call utl_abort('col_getElem: varsLevs index out of range')
      end if
      value = column%all(ilev,headerIndex)
    end if

  end function col_getElem

  !--------------------------------------------------------------------------
  ! col_getLat
  !--------------------------------------------------------------------------
  function col_getLat(column,headerIndex) result(value)
    !
    !:Purpose: Return the latitude for a given header/column index.
    !
    implicit none

    ! Arguments:
    type(struct_columnData),    intent(in) :: column      ! The `columnData` object
    integer,                    intent(in) :: headerIndex ! The column/header index
    ! Result:
    real(8)                                :: value       ! The returned latitude value

    if (headerIndex > column%numCol .or. headerIndex < 1) then
      write(*,*) 'headerIndex = ', headerIndex
      call utl_abort('col_getLat: headerIndex out of range')
    end if

    value = column%lat(headerIndex)

  end function col_getLat

  !--------------------------------------------------------------------------
  ! col_setLat
  !--------------------------------------------------------------------------
  subroutine col_setLat(column, headerIndex, value)
    !
    !:Purpose: Set the latitude for a given header/column index.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(inout) :: column      ! The `columnData` object
    integer,                 intent(in)    :: headerIndex ! The column/header index
    real(8),                 intent(in)    :: value       ! The latitude value for setting

    if (headerIndex > column%numCol .or. headerIndex < 1) then
      write(*,*) 'headerIndex = ', headerIndex
      call utl_abort('col_setLat: headerIndex out of range')
    end if

    column%lat(headerIndex) = value

  end subroutine col_setLat

  !--------------------------------------------------------------------------
  ! col_getNumLev
  !--------------------------------------------------------------------------
  function col_getNumLev(column,varLevel,varName_opt) result(nlev)
    !
    !:Purpose: Return the number of levels for a given type of vertical
    !          level and (optionally) the variable name.
    !
    implicit none

    ! Arguments:
    type(struct_columnData),    intent(in) :: column      ! The `columnData` object
    character(len=*),           intent(in) :: varLevel    ! The type of vertical level
    character(len=*), optional, intent(in) :: varName_opt ! The variable name
    ! Result:
    integer                                :: nlev        ! The returned number of levels

    nlev = vco_getNumLev(column%vco,varLevel,varName_opt)

  end function col_getNumLev

  !--------------------------------------------------------------------------
  ! col_getNumCol
  !--------------------------------------------------------------------------
  function col_getNumCol(column) result(numColumn)
    !
    !:Purpose: Return the number of columns.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in) :: column    ! The `columnData` object
    ! Result:
    integer                             :: numColumn ! The returned number of columns

    numColumn = column%numCol

  end function col_getNumCol

  !--------------------------------------------------------------------------
  ! col_getNumVarLev
  !--------------------------------------------------------------------------
  function col_getNumVarLev(column) result(numVarLev)
    !
    !:Purpose: Return the number of variables x levels.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in) :: column ! The `columnData` object
    ! Result:
    integer                             :: numVarLev   ! The returned number of varsLevs

    numVarLev = column%numVarLev

  end function col_getNumVarLev

  !--------------------------------------------------------------------------
  ! col_addHeightSfcOffset
  !--------------------------------------------------------------------------
  function col_addHeightSfcOffset(column) result(addHeightSfcOffset)
    !
    !:Purpose: Return the value of addHeightSfcOffset
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in) :: column             ! The `columnData` object
    ! Result:
    logical                             :: addHeightSfcOffset ! The returned value of addHeightSfcOffset

    addHeightSfcOffset = column%addHeightSfcOffset

  end function col_addHeightSfcOffset

  !--------------------------------------------------------------------------
  ! col_getVco
  !--------------------------------------------------------------------------
  function col_getVco(column) result(vco_ptr)
    !
    !:Purpose: Return a pointer to the vertical coordinate inside
    !          the `column` object.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in) :: column ! The `columnData` object
    ! Result:
    type(struct_vco), pointer           :: vco_ptr ! The returned `vco` object pointer

    vco_ptr => column%vco

  end function col_getVco

  !--------------------------------------------------------------------------
  ! col_setVco
  !--------------------------------------------------------------------------
  subroutine col_setVco(column,vco_ptr)
    !
    !:Purpose: Set the pointer to the vertical coordinate inside
    !          the `column` object.
    !
    implicit none

    ! Arguments:
    type(struct_columnData),   intent(inout) :: column  ! The `columnData` object
    type(struct_vco), pointer, intent(in)    :: vco_ptr ! The `vco` object pointer for setting

    column%vco => vco_ptr

  end subroutine col_setVco

  !--------------------------------------------------------------------------
  ! col_add
  !--------------------------------------------------------------------------
  subroutine col_add(columnIn, columnInout, scaleFactor_opt)
    !
    ! :Purpose: Adds two columns:
    !
    !           columnInout = columnInout + scaleFactor_opt * columnIn
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)     :: columnIn        ! First operand 
    type(struct_columnData), intent(inout)  :: columnInout     ! Second operand, will receive the result
    real(8), optional,       intent(in)     :: scaleFactor_opt ! Optional scaling of second operand before addition

    ! Locals:
    real(8), pointer                        :: ptrColInOut(:,:)
    real(8), pointer                        :: ptrColIn(:,:)

    if (columnInout%numVarLev /= columnIn%numVarLev) then
      call utl_abort('col_add: Number of levels in columnIn and columnInout are not equal')
    end if

    if (columnInout%numCol /= columnIn%numCol) then
      call utl_abort('col_add: Number of columns in columnIn and columnInout are not equal')
    end if

    if (.not. columnIn%allocated) then
      call utl_abort('col_add: columnIn is not allocated')
    end if

    if (.not. columnInout%allocated) then
      call utl_abort('col_add: columnInout is not allocated')
    end if

    if (any(columnIn%varNumLev(:) /= columnInout%varNumLev(:))) then
      call utl_abort('col_add: varNumLev in columnIn and columnInout are not equal')
    end if
   
    if (.not. vco_equal(col_getVco(columnIn), col_getVco(columnInout))) then
      call utl_abort('col_add: Vco in columnIn and columnInout are not equal')
    end if

    ptrColInOut => columnInout%all
    ptrColIn => columnIn%all

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
    ! :Purpose: Copy column object from columnIn to columnOut.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)     :: columnIn  ! Source column to be copied from
    type(struct_columnData), intent(inout)  :: columnOut ! Destination column to be copied into

    if (columnOut%numVarLev /= columnIn%numVarLev) then
      call utl_abort('col_copy: Number of levels in columnIn and columnOut are not equal')
    end if

    if (columnOut%numCol /= columnIn%numCol) then
      call utl_abort('col_copy: Number of columns in columnIn and columnOut are not equal')
    end if

    if (.not. columnIn%allocated) then
      call utl_abort('col_copy: columnIn is not allocated')
    end if

    if (.not. columnOut%allocated) then
      call utl_abort('col_copy: columnOut is not allocated')
    end if

    if (any(columnIn%varNumLev(:) /= columnOut%varNumLev(:))) then
      call utl_abort('col_copy: varNumLev in columnIn and columnOut are not equal')
    end if
   
    if (.not. vco_equal(col_getVco(columnIn), col_getVco(columnOut))) then
      call utl_abort('col_copy: Vco in columnIn and columnOut are not equal')
    end if
    
    ! Copy content
    columnOut%addHeightSfcOffset = columnIn%addHeightSfcOffset
    columnOut%all(:,:) =  columnIn%all(:,:)
    columnOut%heightSfc(:) = columnIn%heightSfc(:)
    columnOut%oltv(:,:,:) = columnIn%oltv(:,:,:)
    columnOut%lat(:) = columnIn%lat(:)

  end subroutine col_copy

  !--------------------------------------------------------------------------
  ! col_copyLat
  !--------------------------------------------------------------------------
  subroutine col_copyLat(columnIn, columnOut)
    !
    ! :Purpose: Copy only latitude variable from columnIn to columnOut.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)     :: columnIn  ! Source column to be copied from
    type(struct_columnData), intent(inout)  :: columnOut ! Destination column to be copied into

    if (columnOut%numCol /= columnIn%numCol) then
      call utl_abort('col_copyLat: Number of columns in columnIn and columnOut are not equal')
    end if

    if (.not. columnIn%allocated) then
      call utl_abort('col_copyLat: columnIn is not allocated')
    end if

    if (.not. columnOut%allocated) then
      call utl_abort('col_copyLat: columnOut is not allocated')
    end if

    ! Copy latitude
    columnOut%lat(:) = columnIn%lat(:)

  end subroutine col_copyLat

  !--------------------------------------------------------------------------
  ! col_copyHeightSfc
  !--------------------------------------------------------------------------
  subroutine col_copyHeightSfc(columnIn, columnOut)
    !
    ! :Purpose: Copy only surface height variable from columnIn to columnOut.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)     :: columnIn  ! Source column to be copied from
    type(struct_columnData), intent(inout)  :: columnOut ! Destination column to be copied into

    if (columnOut%numCol /= columnIn%numCol) then
      call utl_abort('col_copyHeightSfc: Number of columns in columnIn and columnOut are not equal')
    end if

    if (.not. columnIn%allocated) then
      call utl_abort('col_copyHeightSfc: columnIn is not allocated')
    end if

    if (.not. columnOut%allocated) then
      call utl_abort('col_copyHeightSfc: columnOut is not allocated')
    end if

    ! Copy latitude
    columnOut%heightSfc(:) = columnIn%heightSfc(:)

  end subroutine col_copyHeightSfc

end module columnData_mod

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

module columnData_mod
  ! MODULE columnData_mod (prefix='col' category='2. High-level data objects')
  !
  ! :Purpose: A derived type and related procedures for storing and manipulating
  !           vertical columns of analysis variables on model or analysis grid
  !           levels. These columns are general produced by horizontally
  !           interpolating a gridStateVector object to the observation
  !           locations.
  !
  use mpi_mod
  use earthConstants_mod
  use varNameList_mod
  use verticalCoord_mod
  use mathPhysConstants_mod
  use utilities_mod

  implicit none
  save
  private

  ! public variables and types
  public :: col_rhumin, struct_columnData

  ! public subroutines and functions
  public :: col_setup, col_allocate, col_deallocate
  public :: col_varExist, col_getOffsetFromVarno
  public :: col_getNumLev, col_getNumCol, col_getVarNameFromK
  public :: col_getPressure, col_getPressureDeriv, col_calcPressure, col_vintProf, col_getHeight, col_setHeightSfc
  public :: col_zero, col_getAllColumns, col_getColumn, col_getElem, col_getVco, col_setVco

  type struct_columnData
    integer           :: nk, numCol
    logical           :: allocated=.false.
    logical           :: mpi_local
    real(8), pointer  :: all(:,:)
    real(8), pointer  :: HeightSfc(:,:)
    real(8), pointer  :: oltv(:,:,:)    ! Tangent linear operator of virtual temperature
    integer, pointer  :: varOffset(:),varNumLev(:)
    logical           :: varExistList(vnl_numVarMax)
    type(struct_vco), pointer :: vco => null()
    logical           :: addHeightSfcOffset = .false.
  end type struct_columnData

  real(8) :: rhumin, col_rhumin
  logical varExistList(vnl_numvarmax)
  logical :: addHeightSfcOffset ! controls adding non-zero height offset to diag levels

contains


  subroutine col_setup
    implicit none
    integer :: jvar
    integer :: fnom,fclos,nulnam,ierr
    integer :: numVar3D, numVar2D
    character(len=4) :: anlvar(vnl_numvarmax)
    character(len=8) :: anltime_bin
    logical :: unitConversion_varKindCH
    namelist /namstate/anlvar,rhumin,anltime_bin,addHeightSfcOffset,unitConversion_varKindCH

    if(mpi_myid == 0) write(*,*) 'col_setup: List of known (valid) variable names'
    if(mpi_myid == 0) write(*,*) 'col_setup: varNameList3D=',vnl_varNameList3D
    if(mpi_myid == 0) write(*,*) 'col_setup: varNameList2D=',vnl_varNameList2D
    if(mpi_myid == 0) write(*,*) 'col_setup: varNameList  =',vnl_varNameList

    ! Read NAMELIST NAMSTATE to find which fields are needed

    anlvar(:) = '    '
    rhumin = MPC_MINIMUM_HU_R8
    anltime_bin = 'MIDDLE'
    addHeightSfcOffset = .false.
    unitConversion_varKindCH = .false.

    nulnam=0
    ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=namstate,iostat=ierr)
    if(ierr.ne.0) call utl_abort('col_setup: Error reading namelist')
    if(mpi_myid == 0) write(*,nml=namstate)
    ierr=fclos(nulnam)

    col_rhumin = rhumin

    if( varneed('Z_T') .or. varneed('Z_M') ) call utl_abort('col_setup: height can not be specified as analysis variable in namelist!')
    if( varneed('P_T') .or. varneed('P_M') ) call utl_abort('col_setup: pressure can not be specified as analysis variable in namelist!')

    numVar3D = 0
    numVar2D = 0

    do jvar = 1, vnl_numvarmax3D
      if (varneed(vnl_varNameList3D(jvar))) then
        varExistList(jvar) = .true.
        numVar3D = numVar3D + 1
      else
        varExistList(jvar) = .false.
      endif
    enddo

    do jvar = 1, vnl_numvarmax2D
      if (varneed(vnl_varNameList2D(jvar))) then
        varExistList(jvar+vnl_numvarmax3D) = .true.
        numVar2D = numVar2D + 1
      else
        varExistList(jvar+vnl_numvarmax3D) = .false.
      endif
    enddo

    if(mpi_myid == 0) write(*,*) 'col_setup: numVar3D (no Z_T/Z_M/P_T/P_M included), numVar2D = ', numVar3D, numVar2D
    if(mpi_myid == 0) write(*,*) 'col_setup: varExistList (no Z_T/Z_M/P_T/P_M included) = ',varExistList

    if(mpi_myid == 0) WRITE(*,*)' DIMENSIONS OF MODEL STATE ARRAYS:'
    if(mpi_myid == 0) WRITE(*,FMT=9120) numVar3D, numVar2D
 9120 FORMAT(4X,'  numVar3D =',I6,', numVar2D    =',I6)

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


  subroutine col_zero(column)
    implicit none
    type(struct_columnData) :: column

    if(column%numCol.gt.0) then
      column%all(:,:) = 0.0d0
      column%HeightSfc(:,:) = 0.0d0
    endif

  end subroutine col_zero


  subroutine col_allocate(column, numCol, mpiLocal_opt, beSilent_opt, setToZero_opt)
    implicit none

    ! arguments
    type(struct_columnData) :: column
    integer, intent(in)     :: numCol
    logical, optional       :: mpiLocal_opt
    logical, optional       :: beSilent_opt
    logical, optional       :: setToZero_opt

    ! locals
    integer :: iloc, jvar, jvar2
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

    ! set the variable list using the global ExistList
    column%varExistList(:) = varExistList(:)

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

    column%numCol = numCol
    if(present(mpiLocal_opt)) then
      column%mpi_local=mpiLocal_opt
    else
      column%mpi_local=.true.
      if (mpi_myid == 0 .and. .not.beSilent) write(*,*) 'col_allocate: assuming columnData is mpi-local'
    endif

    if(.not.column%vco%initialized) then
      call utl_abort('col_allocate: VerticalCoord has not been initialized!')
    endif

    allocate(column%varOffset(vnl_numvarmax))
    column%varOffset(:)=0
    allocate(column%varNumLev(vnl_numvarmax))
    column%varNumLev(:)=0

    iloc=0
    do jvar = 1, vnl_numvarmax3d
      if(column%varExistList(jvar)) then
        column%varOffset(jvar)=iloc
        column%varNumLev(jvar)=col_getNumLev(column,vnl_varLevelFromVarname(vnl_varNameList(jvar)))
        iloc = iloc + column%varNumLev(jvar)
      endif
    enddo
    do jvar2 = 1, vnl_numvarmax2d
      jvar=jvar2+vnl_numvarmax3d
      if(column%varExistList(jvar)) then
        column%varOffset(jvar)=iloc
        column%varNumLev(jvar)=1
        iloc = iloc + 1
      endif
    enddo

    if (iloc == 0) then
      call utl_abort('col_allocate: Nothing to allocate')
    end if

    column%nk = iloc

    if(column%numCol.le.0) then
      if ( .not.beSilent ) write(*,*) 'col_allocate: number of columns is zero, not allocated'
    else         
      allocate(column%all(column%nk,column%numCol))
      if ( setToZero ) column%all(:,:)=0.0d0

      allocate(column%HeightSfc(1,column%numCol))
      column%HeightSfc(:,:)=0.0d0

      allocate(column%oltv(2,col_getNumLev(column,'TH'),numCol))
      if ( setToZero ) column%oltv(:,:,:)=0.0d0

    endif
 
    if(mpi_myid == 0 .and. .not.beSilent) write(*,*) 'col_allocate: column%nk = ', column%nk
    if(mpi_myid == 0 .and. .not.beSilent) write(*,*) 'col_allocate: varOffset=',column%varOffset
    if(mpi_myid == 0 .and. .not.beSilent) write(*,*) 'col_allocate: varNumLev=',column%varNumLev

    column%addHeightSfcOffset = addHeightSfcOffset

    column%allocated=.true.

  end subroutine col_allocate


  subroutine col_deallocate(column)
    implicit none

    type(struct_columnData) :: column

    deallocate(column%varOffset)
    deallocate(column%varNumLev)

    if(column%numCol.gt.0) then
      deallocate(column%all)
      deallocate(column%HeightSfc)
      deallocate(column%oltv)
    endif

    column%allocated=.false.

  end subroutine col_deallocate


  function col_varExist(column_opt,varName) result(varExist)
    implicit none
    type(struct_columnData), optional :: column_opt
    character(len=*), intent(in)      :: varName
    logical                           :: varExist 

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

  end function col_varExist


  function col_getOffsetFromVarno(column,varnum,varNumberChm_opt) result(offset)
    implicit none
    type(struct_columnData) :: column
    integer, intent(in)     :: varnum
    integer, intent(in), optional :: varNumberChm_opt
    integer                 :: offset

    offset=column%varOffset(vnl_varListIndex(vnl_varnameFromVarnum(varnum,varNumberChm_opt=varNumberChm_opt)))

  end function col_getOffsetFromVarno

  function col_getVarNameFromK(column,kIndex) result(varName)
    implicit none
    type(struct_columnData) :: column
    integer, intent(in) :: kIndex
    character(len=4)    :: varName
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

  subroutine col_calcPressure(column, beSilent_opt)
    implicit none
    type(struct_columnData), intent(inout) :: column
    logical, optional :: beSilent_opt

    real(kind=8), allocatable :: Psfc(:,:),zppobs2(:,:)
    real(kind=8), pointer     :: zppobs1(:,:,:) => null()
    integer :: headerIndex, status, ilev1, ilev2
    logical                   :: beSilent

    if ( col_getNumCol(column) <= 0 ) return

    if (.not.col_varExist(column,'P0')) then
      call utl_abort('col_calcPressure: P0 must be present as an analysis variable!')
    endif

    allocate(Psfc(1,col_getNumCol(column)))
    do headerIndex = 1,col_getNumCol(column)
      Psfc(1,headerIndex) = col_getElem(column,1,headerIndex,'P0')
    enddo

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .false.
    end if

    if ( .not.beSilent ) write(*,*) 'col_calcPressure: computing pressure on staggered or UNstaggered levels'

    status=vgd_levels(column%vco%vgrid,ip1_list=column%vco%ip1_M,  &
                      levels=zppobs1,sfc_field=Psfc,in_log=.false.)
    if(status.ne.VGD_OK) call utl_abort('ERROR with vgd_levels')
    allocate(zppobs2(col_getNumLev(column,'MM'),col_getNumCol(column)))
    zppobs2 = transpose(zppobs1(1,:,:))
    ilev1 = 1 + column%varOffset(vnl_varListIndex('P_M'))
    ilev2 = ilev1 - 1 + column%varNumLev(vnl_varListIndex('P_M'))
    column%all(ilev1:ilev2,:) = zppobs2(:,:)
    if (associated(zppobs1))  deallocate(zppobs1)
    deallocate(zppobs2)

    status=vgd_levels(column%vco%vgrid,ip1_list=column%vco%ip1_T,  &
                      levels=zppobs1,sfc_field=Psfc,in_log=.false.)
    if(status.ne.VGD_OK) call utl_abort('ERROR with vgd_levels')
    allocate(zppobs2(col_getNumLev(column,'TH'),col_getNumCol(column)))
    zppobs2 = transpose(zppobs1(1,:,:))
    ilev1 = 1 + column%varOffset(vnl_varListIndex('P_T'))
    ilev2 = ilev1 - 1 + column%varNumLev(vnl_varListIndex('P_T'))
    column%all(ilev1:ilev2,:) = zppobs2(:,:)
    if (associated(zppobs1)) deallocate(zppobs1)
    deallocate(zppobs2)

    deallocate(Psfc)

  end subroutine col_calcPressure

  subroutine col_vintprof(column_in,column_out,varName,useColumnPressure_opt)
    implicit none
    type(struct_columnData), intent(inout) :: column_out
    type(struct_columnData), intent(in) :: column_in
    character(len=*) :: varName
    logical, optional :: useColumnPressure_opt

    real(8), pointer :: column_ptr_in(:), column_ptr_out(:)
    real(8), pointer :: pres1Dptr_in(:)  , pres1Dptr_out(:)
    real(8), pointer :: pres3Dptr_in(:,:,:), pres3Dptr_out(:,:,:)
    character(len=2) :: varLevel
    real(8)          :: zwb, zwt
    integer          :: jlevo, jlevi, columnIndex, status
    logical          :: vInterp, useColumnPressure

    integer, allocatable, target :: THlevelWanted(:), MMlevelWanted(:)
    integer, pointer :: levelWanted(:)
    real(8), allocatable :: Psfc(:,:)

    varLevel = vnl_varLevelFromVarname(varName)

    if ( present(useColumnPressure_opt) ) then
      useColumnPressure = useColumnPressure_opt
    else
      useColumnPressure = .true.
    end if

    nullify(pres3Dptr_in)
    nullify(pres3Dptr_out)

    vInterp = .true.
    if ( .not. col_varExist(column_in,'P0' ) ) then
      write(*,*)
      write(*,*) 'col_vintprof: P0 is missing. Vertical interpolation WILL NOT BE PERFORMED'
      vInterp = .false.
    else if ( col_getNumLev(column_in ,'TH') <= 1 .or. &
              col_getNumLev(column_in ,'MM') <= 1 ) then
      vInterp = .false.
      write(*,*)
      write(*,*) 'col_vintprof: The input backgrounds are 2D. Vertical interpolation WILL NOT BE PERFORMED'
    end if

    if (vInterp) then
      if ( .not. useColumnPressure ) then

        ! read Psfc to use in vgd_levels
        allocate(Psfc(1,col_getNumCol(column_in)))
        do columnIndex = 1,col_getNumCol(column_in)
          Psfc(1,columnIndex) = col_getElem(column_out,1,columnIndex,'P0')
        end do

        ! Compute pressure
        if ( varLevel == 'TH' ) then
          ! pres3Dptr_in 
          status = vgd_levels(column_in%vco%vgrid,ip1_list=column_in%vco%ip1_T,  &
                              levels=pres3Dptr_in,sfc_field=Psfc,in_log=.false.)
          if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_levels')

          ! pres3Dptr_out
          status = vgd_levels(column_out%vco%vgrid,ip1_list=column_out%vco%ip1_T,  &
                              levels=pres3Dptr_out,sfc_field=Psfc,in_log=.false.)
          if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_levels')

        else if ( varLevel == 'MM' ) then
          ! pres3Dptr_in 
          status = vgd_levels(column_in%vco%vgrid,ip1_list=column_in%vco%ip1_M,  &
                              levels=pres3Dptr_in,sfc_field=Psfc,in_log=.false.)
          if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_levels')

          ! pres3Dptr_out
          status = vgd_levels(column_out%vco%vgrid,ip1_list=column_out%vco%ip1_M,  &
                              levels=pres3Dptr_out,sfc_field=Psfc,in_log=.false.)
          if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_levels')

        else
          call utl_abort('col_vintprof: only varLevel TH/MM is allowed')
        end if  ! varLevel
      end if    ! useColumnPressure 

      do columnIndex = 1, col_getNumCol(column_out)

        ! pres1Dptr_in 
        if ( .not. useColumnPressure ) then
          pres1Dptr_in => pres3Dptr_in(1,columnIndex,:)
        else
          if ( varLevel == 'TH' ) then
            pres1Dptr_in => col_getColumn(column_in,columnIndex,'P_T')
          else if ( varLevel == 'MM' ) then
            pres1Dptr_in => col_getColumn(column_in,columnIndex,'P_M')
          else
            call utl_abort('col_vintprof: only varLevel TH/MM is allowed')
          end if
        end if

        ! pres1Dptr_out 
        if ( .not. useColumnPressure ) then
          pres1Dptr_out => pres3Dptr_out(1,columnIndex,:)
        else
          if ( varLevel == 'TH' ) then
            pres1Dptr_out => col_getColumn(column_out,columnIndex,'P_T')
          else if ( varLevel == 'MM' ) then
            pres1Dptr_out => col_getColumn(column_out,columnIndex,'P_M')
          else
            call utl_abort('col_vintprof: only varLevel TH/MM is allowed')
          end if
        end if

        column_ptr_in  => col_getColumn(column_in ,columnIndex,varName)
        column_ptr_out => col_getColumn(column_out,columnIndex,varName)

        if ( mpi_myid == 0 .and. columnIndex == 1 .and. &
             (trim(varName) == 'P_T' ) ) then
             !(trim(varName) == 'P_T' .or. trim(varName) == 'P_M') ) then

          write(*,*) 'useColumnPressure=', useColumnPressure

          write(*,*) 'col_vintprof, COLUMN_IN(1):'
          write(*,*) trim(varName),':'
          write(*,*) column_ptr_in(:)

          write(*,*) 'col_vintprof, COLUMN_OUT(1):'
          write(*,*) trim(varName),':'
          write(*,*) column_ptr_out(:)
          write(*,*)
        end if

        jlevi = 1
        do jlevo = 1, col_getNumLev(column_out,varLevel)
          jlevi = jlevi + 1
          do while( pres1Dptr_out(jlevo) .gt. pres1Dptr_in(jlevi) .and. &
               jlevi .lt. col_getNumLev(column_in,varLevel) )
            jlevi = jlevi + 1
          enddo
          jlevi = jlevi - 1
          zwb = log(pres1Dptr_out(jlevo)/pres1Dptr_in(jlevi))/  &
               log(pres1Dptr_in(jlevi+1)/pres1Dptr_in(jlevi))
          zwt = 1. - zwb
          if (  useColumnPressure .and. &
              (trim(varName) == 'P_T' .or. trim(varName) == 'P_M' ) ) then
            ! do nothing, i.e. use the pressures from column_in
          else if ( .not. useColumnPressure .and. &
              (trim(varName) == 'P_T' .or. trim(varName) == 'P_M' ) ) then
            column_ptr_out(jlevo) = exp(zwb*log(column_ptr_in(jlevi+1)) + zwt*log(column_ptr_in(jlevi)))
            !column_ptr_out(jlevo) = zwb*column_ptr_in(jlevi+1) + zwt*column_ptr_in(jlevi)
          else
            column_ptr_out(jlevo) = zwb*column_ptr_in(jlevi+1) + zwt*column_ptr_in(jlevi)
          end if
        enddo
      enddo

      if ( .not. useColumnPressure ) then
        deallocate(pres3Dptr_in)
        deallocate(pres3Dptr_out)
        deallocate(Psfc)
      end if
      
    else

      ! Find which levels in column_in matches column_out
      allocate(THlevelWanted(column_out%vco%nlev_T))
      allocate(MMlevelWanted(column_out%vco%nlev_M))

      call vco_levelMatchingList( THlevelWanted, MMlevelWanted, & ! OUT
                                  column_out%vco, column_in%vco ) ! IN

      if ( any(THlevelWanted == -1) .or. any(MMlevelWanted == -1) ) then
        call utl_abort('col_vintprof: column_out is not a subsets of column_in!')
      end if

      ! Transfer the corresponding data
      do columnIndex = 1, col_getNumCol(column_out)
        column_ptr_in  => col_getColumn(column_in ,columnIndex,varName)
        column_ptr_out => col_getColumn(column_out,columnIndex,varName)
        if (vnl_varLevelFromVarname(varName) == 'TH') then
          levelWanted => THlevelWanted
        else
          levelWanted => MMlevelWanted
        end if
        do jlevo = 1, col_getNumLev(column_out,varLevel)
          column_ptr_out(jlevo) = column_ptr_in(levelWanted(jlevo))
        end do
      end do

      deallocate(THlevelWanted)
      deallocate(MMlevelWanted)

    end if

  end subroutine col_vintprof


  function col_getPressure(column,ilev,headerIndex,varLevel) result(pressure)
    implicit none
    type(struct_columnData), intent(in) :: column
    integer, intent(in)                 :: ilev,headerIndex
    character(len=*), intent(in)        :: varLevel
    real(8)                             :: pressure
    integer                             :: ilev1

    if (varLevel == 'TH' .and. col_varExist(column,'P_T')) then
      ilev1 = 1 + column%varOffset(vnl_varListIndex('P_T'))
      pressure = column%all(ilev1+ilev-1,headerIndex)
    elseif (varLevel == 'MM' .and. col_varExist(column,'P_M') ) then
      ilev1 = 1 + column%varOffset(vnl_varListIndex('P_M'))
      pressure = column%all(ilev1+ilev-1,headerIndex)
    else
      call utl_abort('col_getPressure: Unknown variable type: ' // varLevel)
    endif

  end function col_getPressure
 

  function col_getPressureDeriv(column,ilev,headerIndex,varLevel) result(dP_dPsfc)
    implicit none
    type(struct_columnData), intent(in) :: column
    integer, intent(in)                 :: ilev,headerIndex
    character(len=*), intent(in)        :: varLevel
    real(8)                             :: dP_dPsfc, Psfc
    real(8), pointer                    :: dP_dPsfc_col(:) => null()
    integer :: status

    Psfc = col_getElem(column,1,headerIndex,'P0')

    if (varLevel == 'TH') then
      status = vgd_dpidpis(column%vco%vgrid,column%vco%ip1_T,dP_dPsfc_col,Psfc)
    elseif (varLevel == 'MM' ) then
      status = vgd_dpidpis(column%vco%vgrid,column%vco%ip1_M,dP_dPsfc_col,Psfc)
    else
      call utl_abort('col_getPressureDeriv: Unknown variable type: ' // varLevel)
    endif

    dP_dPsfc = dP_dPsfc_col(ilev)

    deallocate(dP_dPsfc_col)

  end function col_getPressureDeriv


  function col_getHeight(column,ilev,headerIndex,varLevel) result(height)
    implicit none
    type(struct_columnData), intent(in) :: column
    integer, intent(in)                 :: ilev,headerIndex
    character(len=*), intent(in)        :: varLevel
    real(8)                             :: height
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
      height = column%HeightSfc(1,headerIndex)
    else
      call utl_abort('col_getHeight: unknown varLevel! ' // varLevel)
    end if

  end function col_getHeight


  subroutine col_setHeightSfc(column,headerIndex,height)
    implicit none
    type(struct_columnData)             :: column
    integer, intent(in)                 :: headerIndex
    real(8), intent(in)                 :: height

    column%HeightSfc(1,headerIndex) = height

  end subroutine col_setHeightSfc


  function col_getAllColumns(column,varName_opt) result(allColumns)
    implicit none
    type(struct_columnData), intent(in)    :: column
    character(len=*), intent(in), optional :: varName_opt
    real(8), pointer                       :: allColumns(:,:)
    integer                                :: ilev1,ilev2

    if ( column%numCol > 0 ) then
      if(present(varName_opt)) then
        if ( col_varExist(column,varName_opt) ) then
          ilev1 = column%varOffset(vnl_varListIndex(varName_opt))+1
          ilev2 = ilev1 - 1 + column%varNumLev(vnl_varListIndex(varName_opt))
          allColumns => column%all(ilev1:ilev2,:)
        else
          call utl_abort('col_getAllColumns: Unknown variable name! ' // varName_opt)
        endif
      else
        allColumns => column%all(:,:)
      endif
    else
      allColumns => null()
    end if

  end function col_getAllColumns


  function col_getColumn(column,headerIndex,varName_opt) result(onecolumn)
    implicit none
    type(struct_columnData), intent(in)    :: column
    integer, intent(in)                    :: headerIndex
    character(len=*), intent(in), optional :: varName_opt
    real(8), pointer                       :: onecolumn(:)
    integer                                :: ilev1,ilev2

    if(present(varName_opt)) then
      if(col_varExist(column,varName_opt)) then
        ilev1 = column%varOffset(vnl_varListIndex(varName_opt))+1
        ilev2 = ilev1 - 1 + column%varNumLev(vnl_varListIndex(varName_opt))
        onecolumn => column%all(ilev1:ilev2,headerIndex)
      else
        call utl_abort('col_getColumn: Unknown variable name! ' // varName_opt)
      endif
    else
      onecolumn => column%all(:,headerIndex)
    endif

  end function col_getColumn


  function col_getElem(column,ilev,headerIndex,varName_opt) result(value)
    implicit none
    type(struct_columnData), intent(in)    :: column
    integer, intent(in)                    :: ilev
    integer, intent(in)                    :: headerIndex
    character(len=*), intent(in), optional :: varName_opt
    real(8)                                :: value

    if(present(varName_opt)) then
      if(.not.col_varExist(column,varName_opt)) call utl_Abort('col_getElem: Unknown variable name! ' // varName_opt)
      value = column%all(column%varOffset(vnl_varListIndex(varName_opt))+ilev,headerIndex)
    else
      value = column%all(ilev,headerIndex)
    endif

  end function col_getElem


  function col_getNumLev(column,varLevel) result(nlev)
    implicit none
    type(struct_columnData), intent(in) :: column
    character(len=*), intent(in)        :: varLevel
    integer                             :: nlev

    nlev = vco_getNumLev(column%vco,varLevel)

  end function col_getNumLev


  function col_getNumCol(column) result(numColumn)
    implicit none
    type(struct_columnData), intent(in) :: column
    integer                             :: numColumn

    numColumn = column%numCol

  end function col_getNumCol


  function col_getVco(column) result(vco_ptr)
    implicit none
    type(struct_columnData)   :: column
    type(struct_vco), pointer :: vco_ptr

    vco_ptr => column%vco

  end function col_getVco


  subroutine col_setVco(column,vco_ptr)
    implicit none
    type(struct_columnData)   :: column
    type(struct_vco), pointer :: vco_ptr

    column%vco => vco_ptr

  end subroutine col_setVco

end module columnData_mod

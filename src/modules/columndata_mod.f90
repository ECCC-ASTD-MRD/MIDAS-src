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
!! MODULE columnData (prefix="col" category='2. High-level data objects')
!!
!! *Purpose*: A derived type and related procedures for storing and manipulating
!!            vertical columns of analysis variables on model or analysis grid levels.
!!            These columns are general produced by horizontally interpolating
!!            a gridStateVector object to the observation locations.
!!
!--------------------------------------------------------------------------
module columnData_mod
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
  public :: col_getNumLev, col_getNumCol
  public :: col_getPressure, col_getPressureDeriv, col_calcPressure, col_vintProf, col_getHeight, col_getGZsfc, col_setGZsfc
  public :: col_zero, col_fillmvo, col_getAllColumns, col_getColumn, col_getElem, col_getVco, col_setVco

  type struct_columnData
    integer           :: numCol
    logical           :: allocated=.false.
    logical           :: mpi_local
    real(8), pointer  :: all(:,:)
    real(8), pointer  :: gz_T(:,:),gz_M(:,:),gz_sfc(:)
    real(8), pointer  :: pressure_T(:,:),pressure_M(:,:)
    real(8), pointer  :: dP_dPsfc_T(:,:),dP_dPsfc_M(:,:)
    real(8), pointer  :: oltv(:,:,:)    ! Tangent linear operator of virtual temperature
    integer, pointer  :: varOffset(:),varNumLev(:)
    type(struct_vco), pointer :: vco => null()
    logical           :: addGZsfcOffset = .false.
  end type struct_columnData

  real(8) :: rhumin, col_rhumin
  logical nmvoexist(vnl_numvarmax)
  integer :: nvo3d,nvo2d
  logical :: AddGZSfcOffset ! controls adding non-zero GZ offset to diag levels

contains


  subroutine col_setup
    implicit none
    integer :: jvar, ipos
    integer :: fnom,fclos,nulnam,ierr
    character(len=4) :: anlvar(vnl_numvarmax)
    character(len=8) :: anltime_bin
    namelist /namstate/anlvar,rhumin,anltime_bin,AddGZSfcOffset

    if(mpi_myid == 0) write(*,*) 'col_setup: List of known (valid) variable names'
    if(mpi_myid == 0) write(*,*) 'col_setup: varNameList3D=',vnl_varNameList3D
    if(mpi_myid == 0) write(*,*) 'col_setup: varNameList2D=',vnl_varNameList2D
    if(mpi_myid == 0) write(*,*) 'col_setup: varNameList  =',vnl_varNameList

    ! Read NAMELIST NAMSTATE to find which fields are needed

    anlvar(:) = '    '
    rhumin = MPC_MINIMUM_HU_R8
    anltime_bin = 'MIDDLE'
    AddGZSfcOffset = .false.

    nulnam=0
    ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=namstate,iostat=ierr)
    if(ierr.ne.0) call utl_abort('col_setup: Error reading namelist')
    if(mpi_myid == 0) write(*,nml=namstate)
    ierr=fclos(nulnam)

    col_rhumin = rhumin

    if(varneed('GZ')) call utl_abort('col_setup: GZ can no longer be included as a variable in columnData!')

    nvo3d  = 0
    nvo2d  = 0

    do jvar = 1, vnl_numvarmax3D
      if (varneed(vnl_varNameList3D(jvar))) then
        nmvoexist(jvar) = .true.
        nvo3d = nvo3d + 1
      else
        nmvoexist(jvar) = .false.
      endif
    enddo

    do jvar = 1, vnl_numvarmax2D
      if (varneed(vnl_varNameList2D(jvar))) then
        nmvoexist(jvar+vnl_numvarmax3D) = .true.
        nvo2d = nvo2d + 1
      else
        nmvoexist(jvar+vnl_numvarmax3D) = .false.
      endif
    enddo

    if(mpi_myid == 0) write(*,*) 'col_setup: nvo3d,nvo2d=',nvo3d,nvo2d
    if(mpi_myid == 0) write(*,*) 'col_setup: nmvoexist =',nmvoexist

    if(mpi_myid == 0) WRITE(*,*)' DIMENSIONS OF MODEL STATE ARRAYS:'
    if(mpi_myid == 0) WRITE(*,FMT=9120) NVO3D,NVO2D
 9120 FORMAT(4X,'  NVO3D =',I6,' NVO2D    =',I6)

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
      column%pressure_M(:,:) = 0.0d0
      column%pressure_T(:,:) = 0.0d0
      column%dP_dPsfc_T(:,:) = 0.0d0
      column%dP_dPsfc_M(:,:) = 0.0d0
      column%gz_M(:,:) = 0.0d0
      column%gz_T(:,:) = 0.0d0
      column%gz_sfc(:) = 0.0d0
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
    integer :: nkgdimo, iloc, jvar, jvar2
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
      if(nmvoexist(jvar)) then
        column%varOffset(jvar)=iloc
        column%varNumLev(jvar)=col_getNumLev(column,vnl_varLevelFromVarname(vnl_varNameList(jvar)))
        iloc = iloc + column%varNumLev(jvar)
      endif
    enddo
    do jvar2 = 1, vnl_numvarmax2d
      jvar=jvar2+vnl_numvarmax3d
      if(nmvoexist(jvar)) then
        column%varOffset(jvar)=iloc
        column%varNumLev(jvar)=1
        iloc = iloc + 1
      endif
    enddo
    nkgdimo=iloc

    if(column%numCol.le.0) then
      if ( .not.beSilent ) write(*,*) 'col_allocate: number of columns is zero, not allocated'
    else         
      allocate(column%all(nkgdimo,column%numCol))
      if ( setToZero ) column%all(:,:)=0.0d0

      allocate(column%gz_T(col_getNumLev(column,'TH'),column%numCol))
      allocate(column%gz_M(col_getNumLev(column,'MM'),column%numCol))
      allocate(column%gz_sfc(column%numCol))
      if ( setToZero ) column%gz_T(:,:)=0.0d0
      if ( setToZero ) column%gz_M(:,:)=0.0d0
      column%gz_sfc(:)=0.0d0

      allocate(column%pressure_T(col_getNumLev(column,'TH'),column%numCol))
      allocate(column%pressure_M(col_getNumLev(column,'MM'),column%numCol))
      if ( setToZero ) column%pressure_T(:,:)=0.0d0
      if ( setToZero ) column%pressure_M(:,:)=0.0d0

      allocate(column%dP_dPsfc_T(col_getNumLev(column,'TH'),column%numCol))
      allocate(column%dP_dPsfc_M(col_getNumLev(column,'MM'),column%numCol))
      if ( setToZero ) column%dP_dPsfc_T(:,:)=0.0d0
      if ( setToZero ) column%dP_dPsfc_M(:,:)=0.0d0

      allocate(column%oltv(2,col_getNumLev(column,'TH'),numCol))
      if ( setToZero ) column%oltv(:,:,:)=0.0d0

    endif
 
    if(mpi_myid == 0 .and. .not.beSilent) write(*,*) 'col_allocate: nkgdimo = ',nkgdimo
    if(mpi_myid == 0 .and. .not.beSilent) write(*,*) 'col_allocate: varOffset=',column%varOffset
    if(mpi_myid == 0 .and. .not.beSilent) write(*,*) 'col_allocate: varNumLev=',column%varNumLev

    column%addGZsfcOffset = addGZsfcOffset

    column%allocated=.true.

  end subroutine col_allocate


  subroutine col_deallocate(column)
    implicit none

    type(struct_columnData) :: column

    deallocate(column%varOffset)
    deallocate(column%varNumLev)

    if(column%numCol.gt.0) then
      deallocate(column%all)
      deallocate(column%gz_T)
      deallocate(column%gz_M)
      deallocate(column%gz_sfc)
      deallocate(column%pressure_T)
      deallocate(column%pressure_M)
      deallocate(column%dP_dPsfc_T)
      deallocate(column%dP_dPsfc_M)
      deallocate(column%oltv)
    endif

    column%allocated=.false.

  end subroutine col_deallocate


  subroutine col_fillmvo(columnghr,pvar,varName,varLevel_opt)
    !
    !**s/r fillmvo - Fill in values for a complete set of columns at once
    !
    !Arguments
    !
    !       input:
    !             columnghr          : HR or BG column object
    !             VARNAME (character*4): NOMVAR of the state variable
    !             PVAR(knlev,knobs)   : Variable to transfer in COMMVO(G)(HR)
    implicit none

    type(struct_columnData) :: columnghr
    real(8)          :: pvar(:,:)
    character(len=*) :: varName
    character(len=*), optional :: varLevel_opt

    integer jobs, jlev, status, vcode, nlev_M, nlev_T
    real(8), pointer :: column_ptr(:)
    real :: gz_sfcOffset_r4

    status = vgd_get(columnghr%vco%vgrid,key='ig_1 - vertical coord code',value=vcode)
    nlev_M = col_getNumLev(columnghr,'MM')
    nlev_T = col_getNumLev(columnghr,'TH')

    ! Pressure
    select case(trim(varName))
    case('PRES')
      if(present(varLevel_opt)) then
        if(varLevel_opt == 'TH') then
          do jobs = 1, columnghr%numCol
            do jlev = 1, nlev_T
              columnghr%pressure_T(jlev,jobs) = pvar(jlev,jobs)
            enddo
          enddo
        elseif(varLevel_opt == 'MM') then
          do jobs = 1, columnghr%numCol
            do jlev = 1, nlev_M
              columnghr%pressure_M(jlev,jobs) = pvar(jlev,jobs)
            enddo
          enddo
        else
          call utl_abort('col_fillmvo: must specify varLevel TH or MM for Pressure! ' // varLevel_opt)
        endif
      else
        call utl_abort('col_fillmvo: must specify varLevel for Pressure!')
      endif

    ! Height
    case('GZ')
      if( present( varLevel_opt ) ) then
        if( varLevel_opt == 'TH') then
          if (vcode == 5005 .and. AddGZSfcOffset) then
            status = vgd_get(columnghr%vco%vgrid,key='DHT - height of the diagnostic level (t)',value=gz_sfcOffset_r4)
          else
            gz_sfcOffset_r4 = 0.0
          endif
          if(mpi_myid == 0) write(*,*) 'col_fillmvo: GZ offset for near-sfc thermo level is:   ', gz_sfcOffset_r4, ' metres'
          do jobs = 1, columnghr%numCol
            do jlev = 1, nlev_T
              columnghr%gz_T(jlev,jobs) = pvar(jlev,jobs)
              if ( jlev == nlev_T ) then
                ! set the surface GZ
                columnghr%gz_sfc(jobs) = columnghr%gz_T(nlev_T,jobs)
                ! adjust heights of near-surface GZ
                columnghr%gz_T(nlev_T,jobs) = columnghr%gz_T(nlev_T,jobs) + real(gz_sfcOffset_r4,8) * RG
              end if
            end do
          end do
        elseif(varLevel_opt == 'MM') then
          if( vcode == 5005 .and. AddGZSfcOffset ) then
            status = vgd_get( columnghr%vco%vgrid,key = 'DHM - height of the diagnostic level (m)', value = gz_sfcOffset_r4 )
          else
            gz_sfcOffset_r4 = 0.0
          endif
          if( mpi_myid == 0 ) write(*,*) 'col_fillmvo: GZ offset for near-sfc momentum level is: ', gz_sfcOffset_r4, ' metres'
          do jobs = 1, columnghr%numCol
            do jlev = 1, nlev_M
              columnghr%gz_M( jlev, jobs ) = pvar( jlev, jobs )
              if ( jlev == nlev_M ) then
                ! set the surface GZ
                columnghr%gz_sfc(jobs) = columnghr%gz_M( nlev_M, jobs )
                ! adjust heights of near-surface GZ
                columnghr%gz_M( nlev_M, jobs ) = columnghr%gz_M( nlev_M, jobs ) + real( gz_sfcOffset_r4,8 ) * RG
              end if 
            end do
          end do
        else
          call utl_abort('col_fillmvo: must specify varLevel TH or MM for GZ! ' // varLevel_opt )
        end if
      else
        call utl_abort('col_fillmvo: must specify varLevel for GZ!')
      end if

    ! All the other variables that are stored in column%all
    case default
      if ( col_varExist( trim(varName) )) then
        do jobs = 1, columnghr%numCol
          column_ptr => col_getColumn( columnghr, jobs, varName )
          do jlev = 1, col_getNumLev( columnghr, vnl_varLevelFromVarname( varName ) )
            column_ptr(jlev) = pvar( jlev, jobs ) 
          end do
        end do

      ! Unknown variable name
      else
        call utl_abort('col_fillmvo: Unknown variable name: ' // varName)
      endif

    end select

  END SUBROUTINE col_fillmvo


  function col_varExist(varName) result(varExist)
    implicit none
    character(len=*), intent(in) :: varName
    logical                      :: varExist 

    if(trim(varName) == 'GZ' .or. trim(varName) == 'PRES') then
      ! pressure and height always available
      varExist = .true.
    elseif(nmvoexist(vnl_varListIndex(varName))) then
      varExist = .true.
    else
      varExist = .false.
    endif

  end function col_varExist


  function col_getOffsetFromVarno(column,varnum,varNumberChm_opt) result(offset)
    !
    !   Revisions:
    !             Y.J. Rochon (ARQI), Jan. 2015
    !             - Added optional varCHnumber
    !          
    implicit none
    type(struct_columnData) :: column
    integer, intent(in)     :: varnum
    integer, intent(in), optional :: varNumberChm_opt
    integer                 :: offset

    offset=column%varOffset(vnl_varListIndex(vnl_varnameFromVarnum(varnum,varNumberChm_opt=varNumberChm_opt)))

  end function col_getOffsetFromVarno


  subroutine col_calcPressure(column, beSilent_opt)
    implicit none
    type(struct_columnData), intent(inout) :: column
    logical, optional :: beSilent_opt

    real(kind=8), allocatable :: Psfc(:,:),zppobs2(:,:)
    real(kind=8), pointer     :: zppobs1(:,:,:) => null()
    real(kind=8), pointer     :: dP_dPsfc(:,:,:) => null()
    integer :: jobs, status
    logical                   :: beSilent

    if ( col_getNumCol(column) <= 0 ) return

    if (.not.col_varExist('P0')) then
      call utl_abort('col_calcPressure: P0 must be present as an analysis variable!')
    endif

    allocate(Psfc(1,col_getNumCol(column)))
    do jobs = 1,col_getNumCol(column)
      Psfc(1,jobs) = col_getElem(column,1,jobs,'P0')
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
    call col_fillmvo(column,zppobs2,'PRES','MM')
    if (associated(zppobs1))  deallocate(zppobs1)
    deallocate(zppobs2)

    status=vgd_levels(column%vco%vgrid,ip1_list=column%vco%ip1_T,  &
                      levels=zppobs1,sfc_field=Psfc,in_log=.false.)
    if(status.ne.VGD_OK) call utl_abort('ERROR with vgd_levels')
    allocate(zppobs2(col_getNumLev(column,'TH'),col_getNumCol(column)))
    zppobs2 = transpose(zppobs1(1,:,:))
    call col_fillmvo(column,zppobs2,'PRES','TH')
    if (associated(zppobs1)) deallocate(zppobs1)
    deallocate(zppobs2)

    if ( .not.beSilent ) write(*,*) 'col_calcPressure: computing derivate of pressure wrt surface pressure'

    status = vgd_dpidpis(column%vco%vgrid,column%vco%ip1_M,dP_dPsfc,Psfc)
    allocate(zppobs2(col_getNumLev(column,'MM'),col_getNumCol(column)))
    zppobs2 = transpose(dP_dPsfc(1,:,:))
    column%dP_dPsfc_M(:,:) = zppobs2(:,:)
    deallocate(dP_dPsfc)
    deallocate(zppobs2)

    status = vgd_dpidpis(column%vco%vgrid,column%vco%ip1_T,dP_dPsfc,Psfc)
    allocate(zppobs2(col_getNumLev(column,'TH'),col_getNumCol(column)))
    zppobs2 = transpose(dP_dPsfc(1,:,:))
    column%dP_dPsfc_T(:,:) = zppobs2(:,:)
    deallocate(dP_dPsfc)
    deallocate(zppobs2)

    deallocate(Psfc)

  end subroutine col_calcPressure

  !--------------------------------------------------------------------------
  ! col_vintprof
  !--------------------------------------------------------------------------
  subroutine col_vintprof(column_in,column_out,varName)
    implicit none
    type(struct_columnData), intent(inout) :: column_out
    type(struct_columnData), intent(in) :: column_in
    character(len=*) :: varName

    real(kind=8), pointer :: column_ptr_in(:),column_ptr_out(:)
    character(len=2) :: varLevel
    real(kind=8)     :: zwb,zwt
    integer          :: jlevo,jlevi,jprof
    logical          :: vInterp

    integer, allocatable, target :: THlevelWanted(:), MMlevelWanted(:)
    integer, pointer :: levelWanted(:)

    varLevel = vnl_varLevelFromVarname(varName)

    vInterp = .true.
    if ( .not. col_varExist('P0' ) ) then
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
      do jprof = 1, col_getNumCol(column_out)
        column_ptr_in  => col_getColumn(column_in ,jprof,varName)
        column_ptr_out => col_getColumn(column_out,jprof,varName)
        jlevi = 1
        do jlevo = 1, col_getNumLev(column_out,varLevel)
          jlevi = jlevi + 1
          do while(col_getPressure(column_out,jlevo,jprof,varLevel) .gt.  &
               col_getPressure(column_in ,jlevi,jprof,varLevel) .and. &
               jlevi .lt. col_getNumLev(column_in,varLevel) )
            jlevi = jlevi + 1
          enddo
          jlevi = jlevi - 1
          zwb = log(col_getPressure(column_out,jlevo,jprof,varLevel)/col_getPressure(column_in,jlevi,jprof,varLevel))/  &
               log(col_getPressure(column_in,jlevi+1,jprof,varLevel)/col_getPressure(column_in,jlevi,jprof,varLevel))
          zwt = 1. - zwb
          column_ptr_out(jlevo) = zwb*column_ptr_in(jlevi+1) + zwt*column_ptr_in(jlevi)
        enddo
      enddo
      
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
      do jprof = 1, col_getNumCol(column_out)
        column_ptr_in  => col_getColumn(column_in ,jprof,varName)
        column_ptr_out => col_getColumn(column_out,jprof,varName)
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

  !--------------------------------------------------------------------------
  ! col_getPressure
  !--------------------------------------------------------------------------
  function col_getPressure(column,ilev,headerIndex,varLevel) result(pressure)
    implicit none
    type(struct_columnData), intent(in) :: column
    integer, intent(in)                 :: ilev,headerIndex
    character(len=*), intent(in)        :: varLevel
    real(8)                             :: pressure

    if (varLevel == 'TH') then
      pressure = column%pressure_t(ilev,headerIndex)
    elseif (varLevel == 'MM' ) then
      pressure = column%pressure_m(ilev,headerIndex)
    else
      call utl_abort('col_getPressure: Unknown variable type: ' // varLevel)
    endif

  end function col_getPressure
 

  function col_getPressureDeriv(column,ilev,headerIndex,varLevel) result(dP_dPsfc)
    implicit none
    type(struct_columnData), intent(in) :: column
    integer, intent(in)                 :: ilev,headerIndex
    character(len=*), intent(in)        :: varLevel
    real(8)                             :: dP_dPsfc

    if (varLevel == 'TH') then
      dP_dPsfc = column%dP_dPsfc_t(ilev,headerIndex)
    elseif (varLevel == 'MM' ) then
      dP_dPsfc = column%dP_dPsfc_m(ilev,headerIndex)
    else
      call utl_abort('col_getPressureDeriv: Unknown variable type: ' // varLevel)
    endif

  end function col_getPressureDeriv


  function col_getHeight(column,ilev,headerIndex,varLevel) result(height)
    implicit none
    type(struct_columnData), intent(in) :: column
    integer, intent(in)                 :: ilev,headerIndex
    character(len=*), intent(in)        :: varLevel
    real(8)                             :: height
    integer                             :: ilev1

    if (varLevel == 'TH') then
      height = column%gz_t(ilev,headerIndex)
    elseif (varLevel == 'MM' ) then
      height = column%gz_m(ilev,headerIndex)
    elseif (varLevel == 'SF' ) then
      height = column%gz_sfc(headerIndex)
    else
      call utl_abort('col_getHeight: unknown varLevel! ' // varLevel)
    endif

  end function col_getHeight


  function col_getGZsfc(column,headerIndex) result(height)
    implicit none
    type(struct_columnData), intent(in) :: column
    integer, intent(in)                 :: headerIndex
    real(8)                             :: height

    height = column%gz_sfc(headerIndex)

  end function col_getGZsfc


  subroutine col_setGZsfc(column,headerIndex,height)
    implicit none
    type(struct_columnData)             :: column
    integer, intent(in)                 :: headerIndex
    real(8), intent(in)                 :: height

    column%gz_sfc(headerIndex) = height

  end subroutine col_setGZsfc


  function col_getAllColumns(column,varName_opt) result(allColumns)
    implicit none
    type(struct_columnData), intent(in)    :: column
    character(len=*), intent(in), optional :: varName_opt
    real(8), pointer                       :: allColumns(:,:)
    integer                                :: ilev1,ilev2

    if ( column%numCol > 0 ) then
      if(present(varName_opt)) then
        if(trim(varName_opt) == 'GZ') then
          call utl_abort('col_getAllColumns: Cannot call this for GZ!')
        elseif(col_varExist(varName_opt)) then
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


  function col_getColumn(column,headerIndex,varName_opt,varLevel_opt) result(onecolumn)
    implicit none
    type(struct_columnData), intent(in)    :: column
    integer, intent(in)                    :: headerIndex
    character(len=*), intent(in), optional :: varName_opt
    character(len=*), intent(in), optional :: varLevel_opt
    real(8), pointer                       :: onecolumn(:)
    integer                                :: ilev1,ilev2

    if(present(varName_opt)) then
      if(col_varExist(varName_opt)) then

        select case(trim(varName_opt))
        case('PRES')
          if(present(varLevel_opt)) then 
            if(varLevel_opt == 'TH') then
              onecolumn => column%pressure_T(:,headerIndex)
            elseif(varLevel_opt == 'MM') then
              onecolumn => column%pressure_M(:,headerIndex)
            else
              call utl_abort('col_getColumn: varLevel must MM or TH for Pressure! ' // varLevel_opt)         
            endif
          else
            call utl_abort('col_getColumn: varLevel must be specified for Pressure!')
          endif

        case('GZ')
          if(present(varLevel_opt)) then 
            if(varLevel_opt == 'TH') then
              onecolumn => column%gz_T(:,headerIndex)
            elseif(varLevel_opt == 'MM') then
              onecolumn => column%gz_M(:,headerIndex)
            else
              call utl_abort('col_getColumn: varLevel must MM or TH for Height! ' // varLevel_opt)         
            endif
          else
            call utl_abort('col_getColumn: varLevel must be specified for Height!')         
          endif

        case default ! all other variable names
          ilev1 = column%varOffset(vnl_varListIndex(varName_opt))+1
          ilev2 = ilev1 - 1 + column%varNumLev(vnl_varListIndex(varName_opt))
          onecolumn => column%all(ilev1:ilev2,headerIndex)
        end select

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
      if(.not.col_varExist(varName_opt)) call utl_Abort('col_getElem: Unknown variable name! ' // varName_opt)
      if(trim(varName_opt) == 'GZ') call utl_Abort('col_getElem: cannot call for GZ!')
      if(trim(varName_opt) == 'PRES') call utl_Abort('col_getElem: cannot call for Pressure!')
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

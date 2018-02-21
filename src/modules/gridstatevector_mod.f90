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
!! MODULE gridStateVector (prefix="gsv")
!!
!! *Purpose*: The grid-point state vector and related information.
!!
!--------------------------------------------------------------------------
module gridStateVector_mod
  use mpivar_mod
  use earthConstants_mod
  use varNameList_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use mathPhysConstants_mod
  use timeCoord_mod
  use utilities_mod
  implicit none
  save
  private

  ! public structure definition
  public :: struct_gsv

  ! public subroutines and functions
  public :: gsv_setup, gsv_allocate, gsv_deallocate, gsv_zero, gsv_3dto4d, gsv_3dto4dAdj
  public :: gsv_getOffsetFromVarName, gsv_getLevFromK, gsv_getVarNameFromK, gsv_hPad
  public :: gsv_writeToFileMpi, gsv_writeToFile, gsv_readFromFile, gsv_readTrials, gsv_readFile
  public :: gsv_hInterpolate, gsv_hInterpolate_r4, gsv_vInterpolate, gsv_vInterpolate_r4
  public :: gsv_getField_r8, gsv_getField3D_r8, gsv_getField_r4, gsv_getField3D_r4
  public :: gsv_getField_i2, gsv_getField3D_i2, gsv_convertToInteger
  public :: gsv_getFieldUV_r8, gsv_getFieldUV_r4, gsv_getGZsfc
  public :: gsv_getIntOffset, gsv_getIntMultFactor
  public :: gsv_getDateStamp, gsv_getNumLev, gsv_getNumLevFromVarName
  public :: gsv_add, gsv_power, gsv_scale, gsv_scaleVertical, gsv_copy, gsv_stddev
  public :: gsv_getVco, gsv_getHco
  public :: gsv_horizSubSample, gsv_interpolateAndAdd, gsv_interpolate
  public :: gsv_varKindExist, gsv_varExist
  public :: gsv_multEnergyNorm, gsv_dotProduct

  ! public entities accessed through inheritance
  public :: struct_vco, vco_SetupFromFile
  public :: vnl_varnameFromVarnum, vnl_varLevelFromVarnum, vnl_varLevelFromVarname
  public :: vnl_numvarmax2d, vnl_numvarmax3d,vnl_numvarmax
  public :: vnl_varNameList2d, vnl_varNameList3d, vnl_varNameList
  public :: vgd_get,vgd_levels,vgd_ok,vgd_dpidpis,vgd_write

  type struct_gsv
    ! These are the main data storage arrays
    real(8), pointer    :: gd_r8(:,:,:,:) => null()
    real(8), pointer    :: gd3d_r8(:,:,:) => null()
    real(4), pointer    :: gd_r4(:,:,:,:) => null()
    real(4), pointer    :: gd3d_r4(:,:,:) => null()
    integer(2), pointer :: gd_i2(:,:,:,:) => null()
    integer(2), pointer :: gd3d_i2(:,:,:) => null()
    logical             :: gzSfcPresent = .false.
    real(8), pointer    :: gzSfc(:,:) => null()  ! surface GZ, if VarsLevs then only on proc 0
    ! These are used when distribution is VarLevs to keep corresponding UV
    ! components together on each mpi task to facilitate horizontal interpolation
    logical             :: UVComponentPresent = .false.
    integer             :: myUVkBeg, myUVkEnd
    real(8), pointer    :: gdUV_r8(:,:,:,:) => null()
    real(4), pointer    :: gdUV_r4(:,:,:,:) => null()
    integer(2), pointer :: gdUV_i2(:,:,:,:) => null()
    ! These are needed for converting between real4 <-> integer2
    real(8), pointer    :: intOffset(:,:) => null()
    real(8), pointer    :: intMultFactor(:,:) => null()
    ! All the remaining extra information
    integer             :: dataKind = 8
    integer             :: ni, nj, nk, numStep, anltime
    integer             :: latPerPE, latPerPEmax, myLatBeg, myLatEnd
    integer             :: lonPerPE, lonPerPEmax, myLonBeg, myLonEnd
    integer             :: mykCount, mykBeg, mykEnd
    integer, pointer    :: allLatBeg(:), allLatEnd(:), allLatPerPE(:)
    integer, pointer    :: allLonBeg(:), allLonEnd(:), allLonPerPE(:)
    integer, pointer    :: allkCount(:), allkBeg(:), allkEnd(:)
    integer, pointer    :: dateStampList(:) => null()
    integer, pointer    :: dateStamp3d
    logical             :: allocated=.false.
    type(struct_vco), pointer :: vco => null()
    type(struct_hco), pointer :: hco => null()
    integer, pointer    :: varOffset(:), varNumLev(:)
    logical             :: mpi_local=.false.
    character(len=8)    :: mpi_distribution='None'  ! or "Tiles" or "VarsLevs"
    integer             :: horizSubSample
    logical             :: varExistList(vnl_numVarMax)
    character(len=12)   :: hInterpolateDegree='Empty' ! or "CUBIC" or "NEAREST"
  end type struct_gsv  

  logical :: varExistList(vnl_numVarMax)
  character(len=8) :: ANLTIME_BIN
  integer :: get_max_rss
  real(8) :: rhumin

  contains

  !--------------------------------------------------------------------------
  ! gsv_getOffsetFromVarName
  !--------------------------------------------------------------------------
  function gsv_getOffsetFromVarName(statevector,varName) result(offset)
    implicit none
    type(struct_gsv)             :: statevector
    character(len=*), intent(in) :: varName
    integer                      :: offset

    offset=statevector%varOffset(vnl_varListIndex(varName))

  end function gsv_getOffsetFromVarName

  !--------------------------------------------------------------------------
  ! gsv_getVarNameFromK
  !--------------------------------------------------------------------------
  function gsv_getVarNameFromK(statevector,kIndex) result(varName)
    implicit none
    type(struct_gsv)    :: statevector
    integer, intent(in) :: kIndex
    character(len=4)    :: varName
    integer             :: varIndex

    do varIndex = 1, vnl_numvarmax
      if ( statevector%varExistList(varIndex) ) then
        if ( (kIndex >= (statevector%varOffset(varIndex) + 1)) .and.  &
            (kIndex <= (statevector%varOffset(varIndex) + statevector%varNumLev(varIndex))) ) then
          varName = vnl_varNameList(varIndex)
          return
        end if
      end if
    end do

    write(*,*) 'gsv_getVarNameFromK: kIndex out of range: ', kIndex
    call utl_abort('gsv_getVarNameFromK')

  end function gsv_getVarNameFromK

  !--------------------------------------------------------------------------
  ! gsv_getLevFromK
  !--------------------------------------------------------------------------
  function gsv_getLevFromK(statevector,kIndex) result(levIndex)
    implicit none
    type(struct_gsv)    :: statevector
    integer, intent(in) :: kIndex
    integer             :: levIndex
    integer             :: varIndex

    do varIndex = 1, vnl_numvarmax
      if ( statevector%varExistList(varIndex) ) then
        if ( (kIndex >= (statevector%varOffset(varIndex) + 1)) .and.  &
            (kIndex <= (statevector%varOffset(varIndex) + statevector%varNumLev(varIndex))) ) then
          levIndex = kIndex - statevector%varOffset(varIndex)
          return
        end if
      end if
    end do

    write(*,*) 'gsv_getLevFromK: kIndex out of range: ', kIndex
    call utl_abort('gsv_getLevFromK')

  end function gsv_getLevFromK

  !--------------------------------------------------------------------------
  ! gsv_varExist
  !--------------------------------------------------------------------------
  function gsv_varExist(statevector_opt,varName) result(varExist)
    implicit none
    type(struct_gsv), optional   :: statevector_opt
    character(len=*), intent(in) :: varName
    logical                      :: varExist 

    if ( present(statevector_opt) ) then
      if ( statevector_opt%varExistList(vnl_varListIndex(varName)) ) then
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

  end function gsv_varExist

  !--------------------------------------------------------------------------
  ! gsv_getNumLev
  !--------------------------------------------------------------------------
  function gsv_getNumLev(statevector,varLevel) result(nlev)
    implicit none
    type(struct_gsv), intent(in)  :: statevector
    character(len=*), intent(in)  :: varLevel
    integer                       :: nlev

    nlev = vco_getNumLev(statevector%vco,varLevel)

  end function gsv_getNumLev

  !--------------------------------------------------------------------------
  ! gsv_getNumLevFromVarName
  !--------------------------------------------------------------------------
  function gsv_getNumLevFromVarName(statevector,varName) result(nlev)
    implicit none
    type(struct_gsv), intent(in)  :: statevector
    character(len=*), intent(in)  :: varName
    integer                       :: nlev

    nlev = statevector%varNumLev(vnl_varListIndex(varName))

  end function gsv_getNumLevFromVarName

  !--------------------------------------------------------------------------
  ! gsv_setup
  !--------------------------------------------------------------------------
  SUBROUTINE gsv_setup
    implicit none
    integer :: varIndex, fnom, fclos, nulnam, ierr
    CHARACTER(len=4) :: ANLVAR(VNL_NUMVARMAX)
    logical :: AddGZSfcOffset = .false. ! controls adding non-zero GZ offset to diag levels
    NAMELIST /NAMSTATE/ANLVAR,rhumin,ANLTIME_BIN,AddGZSfcOffset

    if (mpi_myid.eq.0) write(*,*) 'gsv_setup: List of known (valid) variable names'
    if (mpi_myid.eq.0) write(*,*) 'gsv_setup: varNameList3D=',vnl_varNameList3D
    if (mpi_myid.eq.0) write(*,*) 'gsv_setup: varNameList2D=',vnl_varNameList2D

    ! Read namelist NAMSTATE to find which fields are needed

    ANLVAR(1:vnl_numvarmax) = '    '
    ANLTIME_BIN = 'MIDDLE'
    rhumin = MPC_MINIMUM_HU_R8

    nulnam=0
    ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=namstate,iostat=ierr)
    if (ierr.ne.0) call utl_abort('gsv_setup: Error reading namelist')
    if (mpi_myid.eq.0) write(*,nml=namstate)
    ierr=fclos(nulnam)

    if (varneed('GZ')) call utl_abort('gsv_setup: GZ can no longer be included as a variable in gridStateVector!')

    do varIndex = 1, vnl_numvarmax3D
      if (varneed(vnl_varNameList3D(varIndex))) then
        varExistList(varIndex) = .true.
      else
        varExistList(varIndex) = .false.
      end if
    end do

    do varIndex = 1, vnl_numvarmax2D
      if (varneed(vnl_varNameList2D(varIndex))) then
        varExistList(varIndex+vnl_numVarMax3D) = .true.
      else
        varExistList(varIndex+vnl_numVarMax3D) = .false.
      end if
    end do

    if (mpi_myid.eq.0) write(*,*) 'gsv_setup: global varExistList =',varExistList

    ! Check value for ANLTIME_BIN
    if (ANLTIME_BIN .ne. 'MIDDLE' .and. ANLTIME_BIN .ne. 'FIRST' .and.  ANLTIME_BIN .ne. 'LAST') then
      call utl_abort('gsv_setup: Problem setting ANLTIME_BIN. Verify NAMSTATE namelist. Aborting!')
    end if

    return

    contains

      logical function varneed(varName)
        character(len=*) :: varName
        integer :: varIndex
 
        varneed=.false.
        do varIndex=1,VNL_NUMVARMAX
          if (trim(varName) == trim(anlvar(varIndex))) then
            varneed=.true.
          end if
        end do

      end function varneed

  end subroutine gsv_setup

  !--------------------------------------------------------------------------
  ! gsv_allocate
  !--------------------------------------------------------------------------
  subroutine gsv_allocate(statevector, numStep, hco_ptr, vco_ptr, dateStamp_opt, dateStampList_opt,  &
                          mpi_local_opt, mpi_distribution_opt, horizSubSample_opt,                   &
                          varNames_opt, dataKind_opt, allocGZsfc_opt, hInterpolateDegree_opt)
    implicit none

    ! arguments
    type(struct_gsv)           :: statevector
    integer, intent(in)        :: numStep
    type(struct_hco), pointer  :: hco_ptr
    type(struct_vco), pointer  :: vco_ptr
    integer, optional          :: dateStamp_opt
    integer, optional          :: dateStampList_opt(:)
    logical, optional          :: mpi_local_opt
    character(len=*), optional :: mpi_distribution_opt
    integer, optional          :: horizSubSample_opt
    character(len=*), optional :: varNames_opt(:)  ! allow specification of variables
    integer, optional          :: dataKind_opt
    logical, optional          :: allocGZsfc_opt
    character(len=*), optional :: hInterpolateDegree_opt

    integer :: ierr,iloc,varIndex,varIndex2,stepIndex,lon1,lat1,k1,kIndex

    ! set the horizontal and vertical coordinates
    call gsv_sethco(statevector,hco_ptr)
    call gsv_setvco(statevector,vco_ptr)

    if (.not.statevector%vco%initialized) then
       call utl_abort('statevector_allocate: VerticalCoord has not been initialized!')
    end if

    if ( statevector%allocated ) then
      if (mpi_myid.eq.0) write(*,*) 'gridStateVector already allocated! Deallocating first.'
      call gsv_deallocate(statevector)
    end if

    if ( present(dataKind_opt) ) statevector%dataKind = dataKind_opt

    if ( present(varNames_opt) ) then      
      statevector%varExistList(:) = .false.
      do varIndex2 = 1, size(varNames_opt)
        varIndex = vnl_varListIndex(varNames_opt(varIndex2))
        statevector%varExistList(varIndex) = .true.
      end do
    else
      ! use the global variable list
      statevector%varExistList(:) = varExistList(:)
    end if

    if ( present(horizSubSample_opt) ) then
      ! user has chosen a coarser grid than specified in hco
      statevector%horizSubSample = horizSubSample_opt
    else
      ! default is no sub-sampling
      statevector%horizSubSample = 1
    end if

    if ( present(hInterpolateDegree_opt) ) then
      ! set the horizontal interpolation degree
      statevector%hInterpolateDegree = trim(hInterpolateDegree_opt)
    else
      ! default is linear horizontal interpolation
      statevector%hInterpolateDegree = 'LINEAR'
    end if

    ! compute the number of global grid points for a given subSample level
    statevector%ni = ceiling(real(statevector%hco%ni,8) / real(statevector%horizSubSample,8))
    statevector%nj = ceiling(real(statevector%hco%nj,8) / real(statevector%horizSubSample,8))

    if ( statevector%ni * statevector%horizSubSample /= statevector%hco%ni ) then
      write(*,*) 'gsv_allocate: number of longitudes is not evenly divisible at this subSample level'
      write(*,*) 'gsv_allocate: ni, horizSubSample = ', statevector%ni, statevector%horizSubSample
      call utl_abort('gsv_allocate')
    end if

    if ( statevector%nj * statevector%horizSubSample /= statevector%hco%nj ) then
      write(*,*) 'gsv_allocate: number of latitudes is not evenly divisible at this subSample level'
      write(*,*) 'gsv_allocate: nj, horizSubSample = ', statevector%nj, statevector%horizSubSample
      call utl_abort('gsv_allocate')
    end if

    statevector%numStep=numStep

    if ( present(mpi_local_opt) ) then
      statevector%mpi_local = mpi_local_opt
    else
      statevector%mpi_local = .false.
    end if

    if ( present(mpi_distribution_opt) ) then
      if ( trim(mpi_distribution_opt) .ne. 'Tiles'    .and. &
           trim(mpi_distribution_opt) .ne. 'VarsLevs' .and. &
           trim(mpi_distribution_opt) .ne. 'None' ) then
        call utl_abort('gsv_allocate: Unknown value of mpi_distribution: ' // trim(mpi_distribution_opt))
      end if
      statevector%mpi_distribution = mpi_distribution_opt
    else
      if ( statevector%mpi_local ) then
        statevector%mpi_distribution = 'Tiles'
      else
        statevector%mpi_distribution = 'None'
      end if
    end if

    ! determine lat/lon index ranges
    if( statevector%mpi_distribution == 'Tiles' ) then
      call mpivar_setup_latbands(statevector%nj,  &
                                 statevector%latPerPE, statevector%latPerPEmax, &
                                 statevector%myLatBeg, statevector%myLatEnd)
      call mpivar_setup_lonbands(statevector%ni,  &
                                 statevector%lonPerPE, statevector%lonPerPEmax, &
                                 statevector%myLonBeg, statevector%myLonEnd)
    else
      statevector%latPerPE    = statevector%nj
      statevector%latPerPEmax = statevector%nj
      statevector%myLatBeg    = 1
      statevector%myLatEnd    = statevector%nj
      statevector%lonPerPE    = statevector%ni
      statevector%lonPerPEmax = statevector%ni
      statevector%myLonBeg    = 1
      statevector%myLonEnd    = statevector%ni
    end if

    allocate(statevector%varOffset(vnl_numvarmax))
    statevector%varOffset(:)=0
    allocate(statevector%varNumLev(vnl_numvarmax))
    statevector%varNumLev(:)=0

    iloc=0
    do varIndex = 1, vnl_numvarmax3d
      if ( statevector%varExistList(varIndex) ) then
          statevector%varOffset(varIndex)=iloc
          statevector%varNumLev(varIndex)=gsv_getNumLev(statevector,vnl_varLevelFromVarname(vnl_varNameList(varIndex)))
          iloc = iloc + statevector%varNumLev(varIndex)
        end if
    end do
    do varIndex2 = 1, vnl_numvarmax2d
        varIndex=varIndex2+vnl_numvarmax3d
      if ( statevector%varExistList(varIndex) ) then
          statevector%varOffset(varIndex)=iloc
          statevector%varNumLev(varIndex)=1
          iloc = iloc + 1
        end if
    end do
    statevector%nk=iloc

    ! determine range of values for the "k" index (vars+levels)
    if ( statevector%mpi_distribution == 'VarsLevs' ) then
      call mpivar_setup_varslevels(statevector%nk, statevector%mykBeg, &
                                   statevector%mykEnd, statevector%mykCount)
    else
      statevector%mykCount = statevector%nk
      statevector%mykBeg = 1
      statevector%mykEnd = statevector%nk
    end if

    statevector%UVComponentPresent = .false.
    if ( statevector%mpi_distribution == 'VarsLevs' ) then
      statevector%myUVkBeg = -1
      do kIndex = statevector%mykBeg, statevector%mykEnd
        if ( gsv_getVarNameFromK(statevector,kIndex) == 'UU' .or.  &
            gsv_getVarNameFromK(statevector,kIndex) == 'VV' ) then
          statevector%UVComponentPresent = .true.
          if ( statevector%myUVkBeg == -1 ) statevector%myUVkBeg = kIndex
          statevector%myUVkEnd = kIndex
        end if
      end do
      if (statevector%UVComponentPresent) then
        write(*,*) 'gsv_allocate: UV component present on this mpi task in k range = ', &
                   statevector%myUVkBeg, statevector%myUVkEnd
      end if
    end if

    if ( statevector%mpi_local ) then
      allocate(statevector%allLonBeg(mpi_npex))
      CALL rpn_comm_allgather(statevector%myLonBeg,1,"mpi_integer",       &
                              statevector%allLonBeg,1,"mpi_integer","EW",ierr)
      allocate(statevector%allLonEnd(mpi_npex))
      CALL rpn_comm_allgather(statevector%myLonEnd,1,"mpi_integer",       &
                              statevector%allLonEnd,1,"mpi_integer","EW",ierr)
      allocate(statevector%allLonPerPE(mpi_npex))
      CALL rpn_comm_allgather(statevector%lonPerPE,1,"mpi_integer",       &
                              statevector%allLonPerPE,1,"mpi_integer","EW",ierr)
  
      allocate(statevector%allLatBeg(mpi_npey))
      CALL rpn_comm_allgather(statevector%myLatBeg,1,"mpi_integer",       &
                              statevector%allLatBeg,1,"mpi_integer","NS",ierr)
      allocate(statevector%allLatEnd(mpi_npey))
      CALL rpn_comm_allgather(statevector%myLatEnd,1,"mpi_integer",       &
                              statevector%allLatEnd,1,"mpi_integer","NS",ierr)
      allocate(statevector%allLatPerPE(mpi_npey))
      CALL rpn_comm_allgather(statevector%LatPerPE,1,"mpi_integer",       &
                              statevector%allLatPerPE,1,"mpi_integer","NS",ierr)

      allocate(statevector%allkCount(mpi_nprocs))
      CALL rpn_comm_allgather(statevector%mykCount,1,"mpi_integer",       &
                              statevector%allkCount,1,"mpi_integer","GRID",ierr)
      allocate(statevector%allkBeg(mpi_nprocs))
      CALL rpn_comm_allgather(statevector%mykBeg,1,"mpi_integer",       &
                              statevector%allkBeg,1,"mpi_integer","GRID",ierr)
      allocate(statevector%allkEnd(mpi_nprocs))
      CALL rpn_comm_allgather(statevector%mykEnd,1,"mpi_integer",       &
                              statevector%allkEnd,1,"mpi_integer","GRID",ierr)
    end if

    select case (ANLTIME_BIN)
    case ("FIRST")
       statevector%anltime=1
    case ("MIDDLE")
       statevector%anltime=nint((real(numStep,8)+1.0d0)/2.0d0)
    case ("LAST")
       statevector%anltime=numStep
    end select          

    if (present(dateStamp_opt) .and. present(dateStampList_opt)) then
      call utl_abort('gsv_allocate: Either dateStamp or dateStampList should be presented but not both')
    elseif (present(dateStampList_opt)) then
      allocate(statevector%dateStampList(numStep))
      do stepIndex = 1, numStep
       statevector%dateStampList(stepIndex)= dateStampList_opt(stepIndex)
      end do
      statevector%dateStamp3d => statevector%dateStampList(statevector%anltime)
    elseif (present(dateStamp_opt)) then
      allocate(statevector%dateStampList(numStep))
      call tim_getstamplist(statevector%dateStampList,numStep,dateStamp_opt)
      statevector%dateStamp3d => statevector%dateStampList(statevector%anltime)
    else
      nullify(statevector%dateStamplist)
    end if

    if (statevector%dataKind==8) then
      allocate(statevector%gd_r8(statevector%myLonBeg:statevector%myLonEnd,  &
                                 statevector%myLatBeg:statevector%myLatEnd,  &
                                 statevector%mykBeg:statevector%mykEnd,numStep),stat=ierr)
      if (statevector%UVComponentPresent) then
        allocate(statevector%gdUV_r8(statevector%myLonBeg:statevector%myLonEnd,  &
                                     statevector%myLatBeg:statevector%myLatEnd,  &
                                     statevector%myUVkBeg:statevector%myUVkEnd,numStep),stat=ierr)
      end if
    elseif (statevector%dataKind==4) then
      allocate(statevector%gd_r4(statevector%myLonBeg:statevector%myLonEnd,  &
                                 statevector%myLatBeg:statevector%myLatEnd,  &
                                 statevector%mykBeg:statevector%mykEnd,numStep),stat=ierr)
      if (statevector%UVComponentPresent) then
        allocate(statevector%gdUV_r4(statevector%myLonBeg:statevector%myLonEnd,  &
                                     statevector%myLatBeg:statevector%myLatEnd,  &
                                     statevector%myUVkBeg:statevector%myUVkEnd,numStep),stat=ierr)
      end if
    elseif (statevector%dataKind==2) then
      allocate(statevector%gd_i2(statevector%myLonBeg:statevector%myLonEnd,  &
                                 statevector%myLatBeg:statevector%myLatEnd,  &
                                 statevector%mykBeg:statevector%mykEnd,numStep),stat=ierr)
      if (statevector%UVComponentPresent) then
        allocate(statevector%gdUV_i2(statevector%myLonBeg:statevector%myLonEnd,  &
                                     statevector%myLatBeg:statevector%myLatEnd,  &
                                     statevector%myUVkBeg:statevector%myUVkEnd,numStep),stat=ierr)
      end if
    else
      call utl_abort('gsv_allocate: unknown value of datakind')
    end if
    if (ierr.ne.0) then
      write(*,*) 'gridStateVector: Problem allocating memory! id=1 ',ierr
      call utl_abort('gsv_allocate')
    end if

    if ( present(allocGZsfc_opt) ) then
      if ( allocGZsfc_opt ) then
        ! if VarsLevs, then only proc 0 allocates surface GZ, otherwise all procs do
        statevector%gzSfcPresent = .true.
        if ( ( statevector%mpi_distribution == 'VarsLevs' .and. mpi_myid == 0 ) .or. &
             statevector%mpi_distribution /= 'VarsLevs' ) then
          write(*,*) 'gsv_allocate: allocating gzSfc on this mpi task'
          allocate(statevector%gzSfc(statevector%myLonBeg:statevector%myLonEnd,  &
                                     statevector%myLatBeg:statevector%myLatEnd))
          statevector%gzSfc(:,:) = 0.0d0
        end if
      end if
    end if

    lon1=statevector%myLonBeg
    lat1=statevector%myLatBeg
    k1=statevector%mykBeg
    if (statevector%dataKind==8) then
      statevector%gd3d_r8(lon1:,lat1:,k1:) => statevector%gd_r8(:,:,:,statevector%anltime)
    elseif (statevector%dataKind==4) then
      statevector%gd3d_r4(lon1:,lat1:,k1:) => statevector%gd_r4(:,:,:,statevector%anltime)
    elseif (statevector%dataKind==2) then
      statevector%gd3d_i2(lon1:,lat1:,k1:) => statevector%gd_i2(:,:,:,statevector%anltime)
    end if

    statevector%allocated=.true.

  end subroutine gsv_allocate

  !--------------------------------------------------------------------------
  ! gsv_convertToInteger
  !--------------------------------------------------------------------------
  subroutine gsv_convertToInteger(statevector)
    implicit none
    type(struct_gsv) :: statevector
    integer          :: stepIndex,lonIndex,kIndex,latIndex,lat1,lat2,lon1,lon2,k1,k2
    integer(2)       :: int2value = 1
    real(8)          :: real8value = 1.0d0

    if (.not.statevector%allocated) then
      call utl_abort('gsv_convertToInteger: gridStateVector not yet allocated!')
    end if

    lon1=statevector%myLonBeg
    lon2=statevector%myLonEnd
    lat1=statevector%myLatBeg
    lat2=statevector%myLatEnd
    k1=statevector%mykBeg
    k2=statevector%mykEnd

    if ( .not. associated(statevector%intOffset) ) then
      allocate(statevector%intOffset(k1:k2, statevector%numStep))
    else
      call utl_abort('gsv_convertToInteger: intOffset already allocated!')
    end if
    if ( .not. associated(statevector%intMultFactor) ) then
      allocate(statevector%intMultFactor(k1:k2, statevector%numStep))
    else
      call utl_abort('gsv_convertToInteger: intMultFactor already allocated!')
    end if

    ! only implement for real(4) -> integer(2), for now
    if ( statevector%dataKind /= 4 ) then
      call utl_abort('gsv_convertToInteger: unknown or invalid value of datakind')
    end if

    if ( .not. associated(statevector%gd_i2) ) then
      allocate(statevector%gd_i2(lon1:lon2,lat1:lat2,k1:k2,statevector%numstep))
    else
      call utl_abort('gsv_convertToInteger: gd_i2 already allocated!')
    end if

    do stepIndex = 1, statevector%numStep
      do kIndex = k1, k2
        statevector%intOffset(kIndex,stepIndex) = 0.5d0 * ( real( minval(statevector%gd_r4(:,:,kIndex,stepIndex)), 8 ) +  &
                                                    real( maxval(statevector%gd_r4(:,:,kIndex,stepIndex)), 8 ) )
        statevector%intMultFactor(kIndex,stepIndex) = 0.5d0 * ( real( maxval(statevector%gd_r4(:,:,kIndex,stepIndex)), 8 ) -  &
                                                        real( minval(statevector%gd_r4(:,:,kIndex,stepIndex)), 8 ) ) /  &
                                                      real( huge(int2value) - 1, 8 )

        if ( abs(statevector%intMultFactor(kIndex,stepIndex)) > epsilon(real8value) ) then
          statevector%gd_i2(:,:,kIndex,stepIndex) = nint( ( real(statevector%gd_r4(:,:,kIndex,stepIndex),8) -  &
                                                    statevector%intOffset(kIndex,stepIndex) ) /  &
                                                  statevector%intMultFactor(kIndex,stepIndex), 2 )
        else
          statevector%gd_i2(:,:,kIndex,stepIndex) = nint( ( real(statevector%gd_r4(:,:,kIndex,stepIndex),8) -  &
                                                    statevector%intOffset(kIndex,stepIndex) ), 2 )
        end if
      end do
    end do
    deallocate(statevector%gd_r4)

    statevector%dataKind = 2

  end subroutine gsv_convertToInteger

  !--------------------------------------------------------------------------
  ! gsv_zero
  !--------------------------------------------------------------------------
  subroutine gsv_zero(statevector)
    implicit none
    type(struct_gsv) :: statevector
    integer          :: stepIndex,lonIndex,kIndex,latIndex,lat1,lat2,lon1,lon2,k1,k2

    if (.not.statevector%allocated) then
      call utl_abort('gsv_zero: gridStateVector not yet allocated! Aborting.')
    end if

    lon1=statevector%myLonBeg
    lon2=statevector%myLonEnd
    lat1=statevector%myLatBeg
    lat2=statevector%myLatEnd
    k1=statevector%mykBeg
    k2=statevector%mykEnd

    if ( associated(statevector%gzSfc) ) statevector%gzSfc(:,:) = 0.0d0

    if ( statevector%dataKind==8 ) then

!$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex)    
       do kIndex = k1, k2
          do stepIndex = 1, statevector%numStep
             do latIndex = lat1, lat2
                do lonIndex = lon1, lon2
                   statevector%gd_r8(lonIndex,latIndex,kIndex,stepIndex) = 0.0d0
                end do
             end do
          end do
       end do
!$OMP END PARALLEL DO

    elseif ( statevector%dataKind==4 ) then

!$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex)
       do kIndex = k1, k2
          do stepIndex = 1, statevector%numStep
             do latIndex = lat1, lat2
                do lonIndex = lon1, lon2
                   statevector%gd_r4(lonIndex,latIndex,kIndex,stepIndex) = 0.0
                end do
             end do
          end do
       end do
!$OMP END PARALLEL DO

    else
       call utl_abort('gsv_zero: unknown value of datakind')
    end if
    
  end subroutine gsv_zero

  !--------------------------------------------------------------------------
  ! GSV_interpolateAndAdd
  !--------------------------------------------------------------------------
  subroutine gsv_interpolateAndAdd(statevector_in,statevector_inout,scaleFactor_opt, &
                                   PsfcReference_opt)
    implicit none
    type(struct_gsv)  :: statevector_in, statevector_inout

    real(8), optional :: scaleFactor_opt
    real(8), optional :: PsfcReference_opt(:,:,:)

    type(struct_gsv) :: statevector_in_hvInterp

    character(len=4), pointer :: varNamesToInterpolate(:)

    !
    !- Error traps
    !
    if (.not.statevector_in%allocated) then
      call utl_abort('gsv_interpolateAndAdd: gridStateVector_in not yet allocated! Aborting.')
    end if
    if (.not.statevector_inout%allocated) then
      call utl_abort('gsv_interpolateAndAdd: gridStateVector_inout not yet allocated! Aborting.')
    end if

    nullify(varNamesToInterpolate)
    call vnl_varNamesFromExistList(varNamesToInterpolate, statevector_in%varExistlist(:))

    !
    !- Do the interpolation of statevector_in onto the grid of statevector_inout
    !
    call gsv_allocate(statevector_in_hvInterp, statevector_inout%numstep,                       &
                      statevector_inout%hco, statevector_inout%vco,                             &
                      mpi_local_opt=statevector_inout%mpi_local, mpi_distribution_opt='Tiles',  &
                      dataKind_opt=statevector_inout%dataKind,                                  &
                      allocGZsfc_opt=statevector_inout%gzSfcPresent,                            &
                      varNames_opt=varNamesToInterpolate)

    call gsv_interpolate(statevector_in,statevector_in_hvInterp,PsfcReference_opt=PsfcReference_opt)

    !
    !- Do the summation
    !
    call gsv_add(statevector_in_hvInterp,statevector_inout,scaleFactor_opt=scaleFactor_opt)

    call gsv_deallocate(statevector_in_hvInterp)
    deallocate(varNamesToInterpolate)

  end subroutine gsv_interpolateAndAdd

  !--------------------------------------------------------------------------
  ! GSV_interpolate
  !--------------------------------------------------------------------------
  subroutine gsv_interpolate(statevector_in,statevector_out, &
                             PsfcReference_opt, PsfcReference_r4_opt)
    implicit none
    type(struct_gsv)  :: statevector_in,statevector_out

    real(8), optional :: PsfcReference_opt(:,:,:)
    real(4), optional :: PsfcReference_r4_opt(:,:,:)

    type(struct_gsv) :: statevector_in_varLevs, statevector_in_varLevs_hInterp
    type(struct_gsv) :: statevector_in_hInterp

    character(len=4), pointer :: varNamesToInterpolate(:)

    !
    !- Error traps
    !
    if (.not.statevector_in%allocated) then
      call utl_abort('gsv_interpolate: gridStateVector_in not yet allocated! Aborting.')
    end if
    if (.not.statevector_out%allocated) then
      call utl_abort('gsv_interpolate: gridStateVector_out not yet allocated! Aborting.')
    end if

    !
    !- Do the interpolation of statevector_in onto the grid of statevector_out
    !
    nullify(varNamesToInterpolate)
    call vnl_varNamesFromExistList(varNamesToInterpolate, statevector_in%varExistlist(:))

    !- Horizontal interpolation
    call gsv_allocate(statevector_in_VarLevs, statevector_in%numstep, &
                      statevector_in%hco, statevector_in%vco,          &
                      mpi_local_opt=statevector_in%mpi_local, mpi_distribution_opt='VarsLevs',  &
                      dataKind_opt=statevector_in%dataKind,                                     &
                      allocGZsfc_opt=statevector_in%gzSfcPresent, &
                      varNames_opt=varNamesToInterpolate)

    call transposeLatLonToVarsLevs( statevector_in, statevector_in_VarLevs )

    call gsv_allocate(statevector_in_VarLevs_hInterp, statevector_in%numstep, &
                      statevector_out%hco, statevector_in%vco,  &
                      mpi_local_opt=statevector_out%mpi_local, mpi_distribution_opt='VarsLevs', &
                      dataKind_opt=statevector_out%dataKind,                                    &
                      allocGZsfc_opt=statevector_out%gzSfcPresent, &
                      varNames_opt=varNamesToInterpolate)

    if (statevector_in_VarLevs%dataKind == 4) then
      call gsv_hInterpolate_r4(statevector_in_VarLevs, statevector_in_VarLevs_hInterp)
    else
      call gsv_hInterpolate(statevector_in_VarLevs, statevector_in_VarLevs_hInterp)
    end if
    call gsv_deallocate(statevector_in_VarLevs)

    call gsv_allocate(statevector_in_hInterp, statevector_in%numstep, &
                      statevector_out%hco, statevector_in%vco,      &
                      mpi_local_opt=statevector_out%mpi_local, mpi_distribution_opt='Tiles', &
                      dataKind_opt=statevector_out%dataKind,                                 &
                      allocGZsfc_opt=statevector_out%gzSfcPresent, &
                      varNames_opt=varNamesToInterpolate)

    call transposeVarsLevsToLatLon( statevector_in_varLevs_hInterp, statevector_in_hInterp )
    call gsv_deallocate(statevector_in_varLevs_hInterp)

    !- Vertical interpolation
    if (statevector_in_VarLevs%dataKind == 4) then
      call gsv_vInterpolate_r4(statevector_in_hInterp,statevector_out,PsfcReference_opt=PsfcReference_r4_opt)
    else 
      call gsv_vInterpolate(statevector_in_hInterp,statevector_out,PsfcReference_opt=PsfcReference_opt)
    end if

    call gsv_deallocate(statevector_in_hInterp)
    nullify(varNamesToInterpolate)

  end subroutine gsv_interpolate

  !--------------------------------------------------------------------------
  ! GSV_add
  !--------------------------------------------------------------------------
  subroutine gsv_add(statevector_in,statevector_inout,scaleFactor_opt)
    implicit none
    type(struct_gsv)  :: statevector_in,statevector_inout
    real(8), optional :: scaleFactor_opt

    integer           :: stepIndex,lonIndex,kIndex,latIndex,lon1,lon2,lat1,lat2,k1,k2

    if (.not.statevector_in%allocated) then
      call utl_abort('gsv_add: gridStateVector_in not yet allocated! Aborting.')
    end if
    if (.not.statevector_inout%allocated) then
      call utl_abort('gsv_add: gridStateVector_inout not yet allocated! Aborting.')
    end if

    lon1=statevector_in%myLonBeg
    lon2=statevector_in%myLonEnd
    lat1=statevector_in%myLatBeg
    lat2=statevector_in%myLatEnd
    k1=statevector_in%mykBeg
    k2=statevector_in%mykEnd

    if ( statevector_inout%dataKind == 8 .and. statevector_in%dataKind == 8 ) then

      if (present(scaleFactor_opt)) then
        do stepIndex = 1, statevector_inout%numStep
!$OMP PARALLEL DO PRIVATE (latIndex,kIndex,lonIndex)    
          do kIndex = k1, k2
            do latIndex = lat1, lat2
              do lonIndex = lon1, lon2
                statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex) = statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex) +  &
                                                  scaleFactor_opt * statevector_in%gd_r8(lonIndex,latIndex,kIndex,stepIndex)
              end do
            end do
          end do
!$OMP END PARALLEL DO
        end do
      else
        do stepIndex = 1, statevector_inout%numStep
!$OMP PARALLEL DO PRIVATE (latIndex,kIndex,lonIndex)    
          do kIndex = k1, k2
            do latIndex = lat1, lat2
              do lonIndex = lon1, lon2
                statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex) = statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex) +  &
                                                                statevector_in%gd_r8(lonIndex,latIndex,kIndex,stepIndex)
              end do
            end do
          end do
!$OMP END PARALLEL DO
        end do
      end if

    elseif ( statevector_inout%dataKind == 4 .and. statevector_in%dataKind == 4 ) then

      if (present(scaleFactor_opt)) then
!$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex)    
        do kIndex = k1, k2
          do stepIndex = 1, statevector_inout%numStep
            do latIndex = lat1, lat2
              do lonIndex = lon1, lon2
                statevector_inout%gd_r4(lonIndex,latIndex,kIndex,stepIndex) = statevector_inout%gd_r4(lonIndex,latIndex,kIndex,stepIndex) +  &
                                          real(scaleFactor_opt,4) * statevector_in%gd_r4(lonIndex,latIndex,kIndex,stepIndex)
              end do
            end do
          end do
        end do
!$OMP END PARALLEL DO
      else
!$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex)    
        do kIndex = k1, k2
          do stepIndex = 1, statevector_inout%numStep
            do latIndex = lat1, lat2
              do lonIndex = lon1, lon2
                statevector_inout%gd_r4(lonIndex,latIndex,kIndex,stepIndex) = statevector_inout%gd_r4(lonIndex,latIndex,kIndex,stepIndex) +  &
                                                                statevector_in%gd_r4(lonIndex,latIndex,kIndex,stepIndex)
              end do
            end do
          end do
        end do
!$OMP END PARALLEL DO
      end if

    else
      call utl_abort('gsv_add: Data type must be the same for both statevectors')
    end if

  end subroutine gsv_add

  !--------------------------------------------------------------------------
  ! GSV_copy
  !--------------------------------------------------------------------------
  subroutine gsv_copy(statevector_in,statevector_out)
    implicit none
    type(struct_gsv)  :: statevector_in,statevector_out
    integer           :: stepIndex,lonIndex,kIndex,latIndex,lon1,lon2,lat1,lat2,k1,k2

    if (.not.statevector_in%allocated) then
      call utl_abort('gsv_copy: gridStateVector_in not yet allocated! Aborting.')
    end if
    if (.not.statevector_out%allocated) then
      call utl_abort('gsv_copy: gridStateVector_out not yet allocated! Aborting.')
    end if

    lon1=statevector_in%myLonBeg
    lon2=statevector_in%myLonEnd
    lat1=statevector_in%myLatBeg
    lat2=statevector_in%myLatEnd
    k1=statevector_in%mykBeg
    k2=statevector_in%mykEnd

    if ( associated(statevector_in%gzSfc) .and. associated(statevector_out%gzSfc) ) then
      statevector_out%gzSfc(:,:) = statevector_in%gzSfc(:,:)
    end if

    if ( statevector_out%dataKind == 8 .and. statevector_in%dataKind == 8 ) then

!$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex)    
      do kIndex = k1, k2
        do stepIndex = 1, statevector_out%numStep
          do latIndex = lat1, lat2
            do lonIndex = lon1, lon2
              statevector_out%gd_r8(lonIndex,latIndex,kIndex,stepIndex) =  &
                statevector_in%gd_r8(lonIndex,latIndex,kIndex,stepIndex)
            end do
          end do
        end do
      end do
!$OMP END PARALLEL DO

    elseif ( statevector_out%dataKind == 4 .and. statevector_in%dataKind == 4 ) then

!$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex)    
      do kIndex = k1, k2
        do stepIndex = 1, statevector_out%numStep
          do latIndex = lat1, lat2
            do lonIndex = lon1, lon2
              statevector_out%gd_r4(lonIndex,latIndex,kIndex,stepIndex) =  &
                statevector_in%gd_r4(lonIndex,latIndex,kIndex,stepIndex)
            end do
          end do
        end do
      end do
!$OMP END PARALLEL DO

    else
      call utl_abort('gsv_copy: Data type must be the same for both statevectors')
    end if

    if ( associated(statevector_in%gzSfc) .and. associated(statevector_out%gzSfc) ) then
      statevector_out%gzSfc(:,:) = statevector_in%gzSfc(:,:)
    end if

  end subroutine gsv_copy

  !--------------------------------------------------------------------------
  ! GSV_hPad
  !--------------------------------------------------------------------------
  subroutine gsv_hPad(statevector_in,statevector_out)
    implicit none
    type(struct_gsv)  :: statevector_in,statevector_out

    integer :: stepIndex,lonIndex,kIndex,latIndex
    integer :: lonBeg_in, lonEnd_in, latBeg_in, latEnd_in, kBeg, kEnd

    real(8) :: paddingValue_r8
    real(4) :: paddingValue_r4

    if (.not.statevector_in%allocated) then
      call utl_abort('gsv_hPad: gridStateVector_in not yet allocated! Aborting.')
    end if
    if (.not.statevector_out%allocated) then
      call utl_abort('gsv_hPad: gridStateVector_out not yet allocated! Aborting.')
    end if
    if (statevector_in%mpi_local .or. statevector_out%mpi_local) then
       call utl_abort('gsv_hPad: both gridStateVectors must be NO MPI! Aborting.')
    end if

    lonBeg_in=statevector_in%myLonBeg
    lonEnd_in=statevector_in%myLonEnd
    latBeg_in=statevector_in%myLatBeg
    latEnd_in=statevector_in%myLatEnd
    kBeg=statevector_in%mykBeg
    kEnd=statevector_in%mykEnd

    if (lonBeg_in > statevector_out%myLonBeg .or. &
        lonEnd_in > statevector_out%myLonEnd .or. &
        latBeg_in > statevector_out%myLatBeg .or. &
        latEnd_in > statevector_out%myLatEnd ) then
      call utl_abort('gsv_hPad: StateVector_out is SMALLER than StateVector_in! Aborting.')
    end if
    if ( kBeg /= statevector_out%mykBeg .or. kEnd /= statevector_out%mykEnd) then
      call utl_abort('gsv_hPad: Vertical levels are not compatible! Aborting.')
    end if

    if ( statevector_out%dataKind == 8 .and. statevector_in%dataKind == 8 ) then

      if ( associated(statevector_in%gzSfc) .and. associated(statevector_out%gzSfc) ) then
        statevector_out%gzSfc(:,:) = 0.d0
        statevector_out%gzSfc(lonBeg_in:lonEnd_in,latBeg_in:latEnd_in) = statevector_in%gzSfc(:,:)
      end if

!$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex,paddingValue_r8)    
      do kIndex = kBeg, kEnd
        if (trim(gsv_getVarNameFromK(statevector_out,kIndex)) == 'P0') then
          paddingValue_r8 = 1000.d0 ! 1000 hPa
        else
          paddingValue_r8 = 0.d0
        end if
        do stepIndex = 1, statevector_out%numStep
          statevector_out%gd_r8(:,:,kIndex,stepIndex) = paddingValue_r8
          do latIndex = latBeg_in, latEnd_in
            do lonIndex = lonBeg_in, lonEnd_in
              statevector_out%gd_r8(lonIndex,latIndex,kIndex,stepIndex) = statevector_in%gd_r8(lonIndex,latIndex,kIndex,stepIndex)
            end do
          end do
        end do
      end do
!$OMP END PARALLEL DO

    elseif ( statevector_out%dataKind == 4 .and. statevector_in%dataKind == 4 ) then

      if ( associated(statevector_in%gzSfc) .and. associated(statevector_out%gzSfc) ) then
        statevector_out%gzSfc(:,:) = 0.0
        statevector_out%gzSfc(lonBeg_in:lonEnd_in,latBeg_in:latEnd_in) = statevector_in%gzSfc(:,:)
      end if

!$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex,paddingValue_r4)    
      do kIndex = kBeg, kEnd
        if (trim(gsv_getVarNameFromK(statevector_out,kIndex)) == 'P0') then
          paddingValue_r4 = 1000.0 ! 1000 hPa
        else
          paddingValue_r4 = 0.0
        end if
        do stepIndex = 1, statevector_out%numStep
          statevector_out%gd_r4(:,:,kIndex,stepIndex) = paddingValue_r4
          do latIndex = latBeg_in, latEnd_in
            do lonIndex = lonBeg_in, lonEnd_in
              statevector_out%gd_r4(lonIndex,latIndex,kIndex,stepIndex) = statevector_in%gd_r4(lonIndex,latIndex,kIndex,stepIndex)
            end do
          end do
        end do
      end do
!$OMP END PARALLEL DO

    else
      call utl_abort('gsv_hPad: Data type must be the same for both statevectors')
    end if

  end subroutine gsv_hPad

  !--------------------------------------------------------------------------
  ! GSV_power
  !--------------------------------------------------------------------------
  subroutine gsv_power(statevector_inout,power,scaleFactor_opt)
    implicit none
    type(struct_gsv)    :: statevector_inout
    real(8), intent(in) :: power
    integer          :: stepIndex,lonIndex,kIndex,latIndex,lon1,lon2,lat1,lat2,k1,k2
    real(8), optional :: scaleFactor_opt

    if (.not.statevector_inout%allocated) then
      call utl_abort('gsv_power: gridStateVector_inout not yet allocated! Aborting.')
    end if

    lon1=statevector_inout%myLonBeg
    lon2=statevector_inout%myLonEnd
    lat1=statevector_inout%myLatBeg
    lat2=statevector_inout%myLatEnd
    k1=statevector_inout%mykBeg
    k2=statevector_inout%mykEnd

    if ( statevector_inout%dataKind == 8 ) then

      if (present(scaleFactor_opt)) then
!$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex)    
        do kIndex = k1, k2
          do stepIndex = 1, statevector_inout%numStep
            do latIndex = lat1, lat2
              do lonIndex = lon1, lon2
                statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex) = &
                     scaleFactor_opt * (statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex))**power
              end do
            end do
          end do
        end do
!$OMP END PARALLEL DO
      else
!$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex)    
        do kIndex = k1, k2
          do stepIndex = 1, statevector_inout%numStep
            do latIndex = lat1, lat2
              do lonIndex = lon1, lon2
                statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex) = &
                     (statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex))**power
              end do
            end do
          end do
        end do
!$OMP END PARALLEL DO
      end if

    elseif ( statevector_inout%dataKind == 4 ) then

      if (present(scaleFactor_opt)) then
!$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex)    
        do kIndex = k1, k2
          do stepIndex = 1, statevector_inout%numStep
            do latIndex = lat1, lat2
              do lonIndex = lon1, lon2
                statevector_inout%gd_r4(lonIndex,latIndex,kIndex,stepIndex) = &
                     real(scaleFactor_opt,4) * (statevector_inout%gd_r4(lonIndex,latIndex,kIndex,stepIndex))**power
              end do
            end do
          end do
        end do
!$OMP END PARALLEL DO
      else
!$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex)    
        do kIndex = k1, k2
          do stepIndex = 1, statevector_inout%numStep
            do latIndex = lat1, lat2
              do lonIndex = lon1, lon2
                statevector_inout%gd_r4(lonIndex,latIndex,kIndex,stepIndex) = &
                     (statevector_inout%gd_r4(lonIndex,latIndex,kIndex,stepIndex))**power
              end do
            end do
          end do
        end do
!$OMP END PARALLEL DO
      end if

    end if

  end subroutine gsv_power

  !--------------------------------------------------------------------------
  ! GSV_stddev
  !--------------------------------------------------------------------------
  subroutine gsv_stddev(statevector_in,stddev)
    implicit none
    type(struct_gsv) :: statevector_in
    real(8)          :: stddev(:)
    integer          :: stepIndex,lonIndex,kIndex,latIndex,lon1,lon2,lat1,lat2,k1,k2

    if (.not.statevector_in%allocated) then
      call utl_abort('gsv_stddev: gridStateVector_in not yet allocated! Aborting.')
    end if

    lon1=statevector_in%myLonBeg
    lon2=statevector_in%myLonEnd
    lat1=statevector_in%myLatBeg
    lat2=statevector_in%myLatEnd
    k1=statevector_in%mykBeg
    k2=statevector_in%mykEnd

    stddev(:) = 0.0d0
    do kIndex = k1, k2
      do latIndex = lat1, lat2
        do stepIndex = 1, statevector_in%numStep
          do lonIndex = lon1, lon2
            if ( statevector_in%dataKind == 8 ) then
              stddev(kIndex) = stddev(kIndex) + statevector_in%gd_r8(lonIndex,latIndex,kIndex,stepIndex)**2
            else
              stddev(kIndex) = stddev(kIndex) + real(statevector_in%gd_r4(lonIndex,latIndex,kIndex,stepIndex),8)**2
            end if
          end do
        end do
      end do
      stddev(kIndex) = stddev(kIndex) / real(statevector_in%numStep,8)
      call mpi_allreduce_sumreal8scalar(stddev(kIndex),"GRID")
      stddev(kIndex) = stddev(kIndex) / real(statevector_in%ni * statevector_in%nj,8)
      if (stddev(kIndex).gt.0.0d0) stddev(kIndex) = sqrt(stddev(kIndex))
    end do

  end subroutine gsv_stddev

  !--------------------------------------------------------------------------
  ! GSV_scale
  !--------------------------------------------------------------------------
  subroutine gsv_scale(statevector_inout,scaleFactor)
    implicit none
    type(struct_gsv) :: statevector_inout
    integer          :: stepIndex,lonIndex,kIndex,latIndex,lon1,lon2,lat1,lat2,k1,k2
    real(8)          :: scaleFactor

    if (.not.statevector_inout%allocated) then
      call utl_abort('gsv_scale: gridStateVector_inout not yet allocated! Aborting.')
    end if

    lon1=statevector_inout%myLonBeg
    lon2=statevector_inout%myLonEnd
    lat1=statevector_inout%myLatBeg
    lat2=statevector_inout%myLatEnd
    k1=statevector_inout%mykBeg
    k2=statevector_inout%mykEnd

    if ( statevector_inout%dataKind == 8 ) then

!$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex)    
      do kIndex = k1, k2
        do stepIndex = 1, statevector_inout%numStep
          do latIndex = lat1, lat2
            do lonIndex = lon1, lon2
              statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex) = &
                   scaleFactor * statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex)
            end do
          end do
        end do
      end do
!$OMP END PARALLEL DO

    else

!$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex)    
      do kIndex = k1, k2
        do stepIndex = 1, statevector_inout%numStep
          do latIndex = lat1, lat2
            do lonIndex = lon1, lon2
              statevector_inout%gd_r4(lonIndex,latIndex,kIndex,stepIndex) = &
                   real(scaleFactor,4) * statevector_inout%gd_r4(lonIndex,latIndex,kIndex,stepIndex)
            end do
          end do
        end do
      end do
!$OMP END PARALLEL DO

    end if

  end subroutine gsv_scale


  !--------------------------------------------------------------------------
  ! GSV_scaleVertical
  !--------------------------------------------------------------------------
  subroutine gsv_scaleVertical(statevector_inout,scaleFactor)
    implicit none
    type(struct_gsv) :: statevector_inout
    integer          :: stepIndex,lonIndex,kIndex,latIndex,lon1,lon2,lat1,lat2,k1,k2
    real(8)          :: scaleFactor(:)

    if (.not.statevector_inout%allocated) then
      call utl_abort('gsv_scaleVertical: gridStateVector_inout not yet allocated! Aborting.')
    end if

    lon1=statevector_inout%myLonBeg
    lon2=statevector_inout%myLonEnd
    lat1=statevector_inout%myLatBeg
    lat2=statevector_inout%myLatEnd
    k1=statevector_inout%mykBeg
    k2=statevector_inout%mykEnd

    if ( statevector_inout%dataKind == 8 ) then

!$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex)    
      do kIndex = k1, k2
        do stepIndex = 1, statevector_inout%numStep
          do latIndex = lat1, lat2
            do lonIndex = lon1, lon2
              statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex) = &
                   scaleFactor(kIndex) * statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex)
            end do
          end do
        end do
      end do
!$OMP END PARALLEL DO

    else

!$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex)    
      do kIndex = k1, k2
        do stepIndex = 1, statevector_inout%numStep
          do latIndex = lat1, lat2
            do lonIndex = lon1, lon2
              statevector_inout%gd_r4(lonIndex,latIndex,kIndex,stepIndex) = &
                   real(scaleFactor(kIndex),4) * statevector_inout%gd_r4(lonIndex,latIndex,kIndex,stepIndex)
            end do
          end do
        end do
      end do
!$OMP END PARALLEL DO

    end if

  end subroutine gsv_scaleVertical

  !--------------------------------------------------------------------------
  ! gsv_3dto4d
  !--------------------------------------------------------------------------
  subroutine gsv_3dto4d(statevector_inout)
    implicit none
    type(struct_gsv) :: statevector_inout
    integer          :: stepIndex,lonIndex,kIndex,latIndex,lon1,lon2,lat1,lat2,k1,k2

    if (.not.statevector_inout%allocated) then
      call utl_abort('gsv_3dto4d: statevector not yet allocated! Aborting.')
    end if

    lon1=statevector_inout%myLonBeg
    lon2=statevector_inout%myLonEnd
    lat1=statevector_inout%myLatBeg
    lat2=statevector_inout%myLatEnd
    k1=statevector_inout%mykBeg
    k2=statevector_inout%mykEnd

    if (statevector_inout%numStep.eq.1) return

    if ( statevector_inout%dataKind == 8 ) then

!$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex)    
      do kIndex = k1, k2
        do stepIndex = 1, statevector_inout%numStep
          if (stepIndex.ne.statevector_inout%anltime) then
            do latIndex = lat1, lat2
              do lonIndex = lon1, lon2
                statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex) = statevector_inout%gd3d_r8(lonIndex,latIndex,kIndex)
              end do
            end do
          end if
        end do
      end do
!$OMP END PARALLEL DO

    else

!$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex)    
      do kIndex = k1, k2
        do stepIndex = 1, statevector_inout%numStep
          if (stepIndex.ne.statevector_inout%anltime) then
            do latIndex = lat1, lat2
              do lonIndex = lon1, lon2
                statevector_inout%gd_r4(lonIndex,latIndex,kIndex,stepIndex) = statevector_inout%gd3d_r4(lonIndex,latIndex,kIndex)
              end do
            end do
          end if
        end do
      end do
!$OMP END PARALLEL DO

    end if

  end subroutine gsv_3dto4d

  !--------------------------------------------------------------------------
  ! gsv_3dto4dAdj
  !--------------------------------------------------------------------------
  subroutine gsv_3dto4dAdj(statevector_inout)
    implicit none
    type(struct_gsv) :: statevector_inout
    integer          :: stepIndex,lonIndex,kIndex,latIndex,lon1,lon2,lat1,lat2,k1,k2
    real(4), allocatable :: gd2d_tmp_r4(:,:)
    real(8), allocatable :: gd2d_tmp(:,:)

    if (.not.statevector_inout%allocated) then
      call utl_abort('gsv_3dto4dAdj: statevector not yet allocated! Aborting.')
    end if

    lon1=statevector_inout%myLonBeg
    lon2=statevector_inout%myLonEnd
    lat1=statevector_inout%myLatBeg
    lat2=statevector_inout%myLatEnd
    k1=statevector_inout%mykBeg
    k2=statevector_inout%mykEnd

    if (statevector_inout%numStep.eq.1) return

    if ( statevector_inout%dataKind == 8 ) then

      allocate(gd2d_tmp(lon1:lon2,lat1:lat2))
!$OMP PARALLEL DO PRIVATE (latIndex,kIndex,lonIndex,stepIndex,gd2d_tmp)
      do kIndex = k1, k2
        gd2d_tmp(:,:) = 0.0d0
        do stepIndex = 1, statevector_inout%numStep
          do latIndex = lat1, lat2
            do lonIndex = lon1, lon2
              gd2d_tmp(lonIndex,latIndex) = gd2d_tmp(lonIndex,latIndex) +   &
                                    statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex)
            end do
          end do
        end do
        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2
            statevector_inout%gd3d_r8(lonIndex,latIndex,kIndex) = gd2d_tmp(lonIndex,latIndex)
          end do
        end do
      end do
!$OMP END PARALLEL DO
      deallocate(gd2d_tmp)

    else

      allocate(gd2d_tmp_r4(lon1:lon2,lat1:lat2))
!$OMP PARALLEL DO PRIVATE (latIndex,kIndex,lonIndex,stepIndex,gd2d_tmp)
      do kIndex = k1, k2
        gd2d_tmp_r4(:,:) = 0.0d0
        do stepIndex = 1, statevector_inout%numStep
          do latIndex = lat1, lat2
            do lonIndex = lon1, lon2
              gd2d_tmp_r4(lonIndex,latIndex) = gd2d_tmp_r4(lonIndex,latIndex) +   &
                                       statevector_inout%gd_r4(lonIndex,latIndex,kIndex,stepIndex)
            end do
          end do
        end do
        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2
            statevector_inout%gd3d_r4(lonIndex,latIndex,kIndex) = gd2d_tmp_r4(lonIndex,latIndex)
          end do
        end do
      end do
!$OMP END PARALLEL DO
      deallocate(gd2d_tmp_r4)

    end if

  end subroutine gsv_3dto4dAdj

  !--------------------------------------------------------------------------
  ! gsv_deallocate
  !--------------------------------------------------------------------------
  subroutine gsv_deallocate(statevector)
    implicit none

    type(struct_gsv) :: statevector
    integer        :: ierr

    if (.not.statevector%allocated) then
      call utl_abort('gsv_deallocate: gridStateVector not yet allocated! Aborting.')
    end if

    statevector%allocated=.false.

    if ( statevector%mpi_local ) then
      deallocate(statevector%allLonBeg)
      deallocate(statevector%allLonEnd)
      deallocate(statevector%allLonPerPE)
      deallocate(statevector%allLatBeg)
      deallocate(statevector%allLatEnd)
      deallocate(statevector%allLatPerPE)
      deallocate(statevector%allkBeg)
      deallocate(statevector%allkEnd)
      deallocate(statevector%allkCount)
    end if

    if (statevector%dataKind==8) then
      deallocate(statevector%gd_r8,stat=ierr)
      nullify(statevector%gd_r8)
      if (statevector%UVComponentPresent) then 
        deallocate(statevector%gdUV_r8)
        nullify(statevector%gdUV_r8)
      end if
    elseif (statevector%dataKind==4) then
      deallocate(statevector%gd_r4,stat=ierr)
      nullify(statevector%gd_r4)
      if (statevector%UVComponentPresent) then 
        deallocate(statevector%gdUV_r4)
        nullify(statevector%gdUV_r4)
      end if
    elseif (statevector%dataKind==2) then
      deallocate(statevector%gd_i2,stat=ierr)
      nullify(statevector%gd_i2)
      if (statevector%UVComponentPresent) then 
        deallocate(statevector%gdUV_i2)
        nullify(statevector%gdUV_i2)
      end if
    end if
    if (ierr.ne.0) then
      write(*,*) 'gsv_deallocate: Problem detected. IERR =',ierr
    end if

    statevector%gzSfcPresent = .false.
    if ( associated(statevector%gzSfc) ) then
      deallocate(statevector%gzSfc)
      nullify(statevector%gzSfc)
    end if

    if ( associated(statevector%dateStampList) ) then
      deallocate(statevector%dateStampList)
      nullify(statevector%dateStampList)
    end if
    deallocate(statevector%varOffset)
    deallocate(statevector%varNumLev)

  end subroutine GSV_deallocate

  !--------------------------------------------------------------------------
  ! gsv_getField_r8
  !--------------------------------------------------------------------------
  function gsv_getField_r8(statevector,varName_opt) result(field)
    implicit none
    type(struct_gsv), intent(in)           :: statevector
    character(len=*), intent(in), optional :: varName_opt
    real(8),pointer                        :: field(:,:,:,:)
    integer                                :: ilev1,ilev2,lon1,lat1,k1

    lon1 = statevector%myLonBeg
    lat1 = statevector%myLatBeg
    k1 = statevector%mykBeg

    if (.not. associated(statevector%gd_r8)) call utl_abort('gsv_getField_r8: data with type r8 not allocated')

    if (present(varName_opt)) then
      if (statevector%mpi_distribution == 'VarsLevs') then
        call utl_abort('gsv_getField_r8: cannot specify a varName for VarsLevs mpi distribution')
      end if
      if (gsv_varExist(statevector,varName_opt)) then
        ilev1 = 1 + statevector%varOffset(vnl_varListIndex(varName_opt))
        ilev2 = ilev1 - 1 + statevector%varNumLev(vnl_varListIndex(varName_opt))
        field(lon1:,lat1:,1:,1:) => statevector%gd_r8(:,:,ilev1:ilev2,:)
      else
        call utl_abort('gsv_getField_r8: Unknown variable name! ' // varName_opt)
      end if
    else
      field(lon1:,lat1:,k1:,1:) => statevector%gd_r8(:,:,:,:)
    end if

  end function gsv_getField_r8

  !--------------------------------------------------------------------------
  ! gsv_getField3D_r8
  !--------------------------------------------------------------------------
  function gsv_getField3D_r8(statevector,varName_opt,stepIndex_opt) result(field3D)
    implicit none
    type(struct_gsv), intent(in)           :: statevector
    character(len=*), intent(in), optional :: varName_opt
    integer, intent(in), optional          :: stepIndex_opt
    real(8),pointer                        :: field3D(:,:,:)
    integer                                :: ilev1,ilev2,lon1,lat1,k1

    lon1=statevector%myLonBeg
    lat1=statevector%myLatBeg
    k1=statevector%mykBeg

    if (.not. associated(statevector%gd3d_r8)) call utl_abort('gsv_getField3D_r8: data with type r8 not allocated')

    if (present(varName_opt)) then
      if (statevector%mpi_distribution == 'VarsLevs') then
        call utl_abort('gsv_getField3D_r8: cannot specify a varName for VarsLevs mpi distribution')
      end if
      if (gsv_varExist(statevector,varName_opt)) then
        ilev1 = 1 + statevector%varOffset(vnl_varListIndex(varName_opt))
        ilev2 = ilev1 - 1 + statevector%varNumLev(vnl_varListIndex(varName_opt))
        if (present(stepIndex_opt)) then
          field3D(lon1:,lat1:,1:) => statevector%gd_r8(:,:,ilev1:ilev2,stepIndex_opt)
        else
          field3D(lon1:,lat1:,1:) => statevector%gd3d_r8(:,:,ilev1:ilev2)
        end if
      else
        call utl_abort('gsv_getField3D_r8: Unknown variable name! ' // varName_opt)
      end if
    else
      if (present(stepIndex_opt)) then
        field3D(lon1:,lat1:,k1:) => statevector%gd_r8(:,:,:,stepIndex_opt)
      else
        field3D(lon1:,lat1:,k1:) => statevector%gd3d_r8(:,:,:)
      end if
    end if

  end function gsv_getField3D_r8

  !--------------------------------------------------------------------------
  ! gsv_getField_r4
  !--------------------------------------------------------------------------
  function gsv_getField_r4(statevector,varName_opt) result(field)
    implicit none
    type(struct_gsv), intent(in)           :: statevector
    character(len=*), intent(in), optional :: varName_opt
    real(4),pointer                        :: field(:,:,:,:)
    integer                                :: ilev1,ilev2,lon1,lat1,k1

    lon1=statevector%myLonBeg
    lat1=statevector%myLatBeg
    k1=statevector%mykBeg

    if (.not. associated(statevector%gd_r4)) call utl_abort('gsv_getField_r4: data with type r4 not allocated')

    if (present(varName_opt)) then
      if (statevector%mpi_distribution == 'VarsLevs') then
        call utl_abort('gsv_getField_r4: cannot specify a varName for VarsLevs mpi distribution')
      end if
      if (gsv_varExist(statevector,varName_opt)) then
        ilev1 = 1 + statevector%varOffset(vnl_varListIndex(varName_opt))
        ilev2 = ilev1 - 1 + statevector%varNumLev(vnl_varListIndex(varName_opt))
        field(lon1:,lat1:,1:,1:) => statevector%gd_r4(:,:,ilev1:ilev2,:)
      else
        call utl_abort('gsv_getField_r4: Unknown variable name! ' // varName_opt)
      end if
    else
      field(lon1:,lat1:,k1:,1:) => statevector%gd_r4(:,:,:,:)
    end if

  end function gsv_getField_r4

  !--------------------------------------------------------------------------
  ! gsv_getField3D_r4
  !--------------------------------------------------------------------------
  function gsv_getField3D_r4(statevector,varName_opt,stepIndex_opt) result(field3D)
    implicit none
    type(struct_gsv), intent(in)           :: statevector
    character(len=*), intent(in), optional :: varName_opt
    integer, intent(in), optional          :: stepIndex_opt
    real(4),pointer                        :: field3D(:,:,:)
    integer                                :: ilev1,ilev2,lon1,lat1,k1

    lon1=statevector%myLonBeg
    lat1=statevector%myLatBeg
    k1=statevector%mykBeg

    if (.not. associated(statevector%gd3d_r4)) call utl_abort('gsv_getField3D_r4: data with type r4 not allocated')

    if (present(varName_opt)) then
      if (statevector%mpi_distribution == 'VarsLevs') then
        call utl_abort('gsv_getField3D_r4: cannot specify a varName for VarsLevs mpi distribution')
      end if
      if (gsv_varExist(statevector,varName_opt)) then
        ilev1 = 1 + statevector%varOffset(vnl_varListIndex(varName_opt))
        ilev2 = ilev1 - 1 + statevector%varNumLev(vnl_varListIndex(varName_opt))
        if (present(stepIndex_opt)) then
          field3D(lon1:,lat1:,1:) => statevector%gd_r4(:,:,ilev1:ilev2,stepIndex_opt)
        else
          field3D(lon1:,lat1:,1:) => statevector%gd3d_r4(:,:,ilev1:ilev2)
        end if
      else
        call utl_abort('gsv_getField3D_r4: Unknown variable name! ' // varName_opt)
      end if
    else
      if (present(stepIndex_opt)) then
        field3D(lon1:,lat1:,k1:) => statevector%gd_r4(:,:,:,stepIndex_opt)
      else
        field3D(lon1:,lat1:,k1:) => statevector%gd3d_r4(:,:,:)
      end if
    end if

  end function gsv_getField3D_r4

  !--------------------------------------------------------------------------
  ! gsv_getField_i2
  !--------------------------------------------------------------------------
  function gsv_getField_i2(statevector,varName_opt) result(field)
    implicit none
    type(struct_gsv), intent(in)           :: statevector
    character(len=*), intent(in), optional :: varName_opt
    integer(2),pointer                     :: field(:,:,:,:)
    integer                                :: ilev1,ilev2,lon1,lat1,k1

    lon1=statevector%myLonBeg
    lat1=statevector%myLatBeg
    k1=statevector%mykBeg

    if (.not. associated(statevector%gd_i2)) call utl_abort('gsv_getField_i2: data with type i2 not allocated')

    if (present(varName_opt)) then
      if (statevector%mpi_distribution == 'VarsLevs') then
        call utl_abort('gsv_getField_i2: cannot specify a varName for VarsLevs mpi distribution')
      end if
      if (gsv_varExist(statevector,varName_opt)) then
        ilev1 = 1 + statevector%varOffset(vnl_varListIndex(varName_opt))
        ilev2 = ilev1 - 1 + statevector%varNumLev(vnl_varListIndex(varName_opt))
        field(lon1:,lat1:,1:,1:) => statevector%gd_i2(:,:,ilev1:ilev2,:)
      else
        call utl_abort('gsv_getField_i2: Unknown variable name! ' // varName_opt)
      end if
    else
      field(lon1:,lat1:,k1:,1:) => statevector%gd_i2(:,:,:,:)
    end if

  end function gsv_getField_i2

  !--------------------------------------------------------------------------
  ! gsv_getField3D_i2
  !--------------------------------------------------------------------------
  function gsv_getField3D_i2(statevector,varName_opt,stepIndex_opt) result(field3D)
    implicit none
    type(struct_gsv), intent(in)           :: statevector
    character(len=*), intent(in), optional :: varName_opt
    integer, intent(in), optional          :: stepIndex_opt
    integer(2),pointer                     :: field3D(:,:,:)
    integer                                :: ilev1,ilev2,lon1,lat1,k1

    lon1=statevector%myLonBeg
    lat1=statevector%myLatBeg
    k1=statevector%mykBeg

    if (.not. associated(statevector%gd3d_i2)) call utl_abort('gsv_getField3D_i2: data with type i2 not allocated')

    if (present(varName_opt)) then
      if (statevector%mpi_distribution == 'VarsLevs') then
        call utl_abort('gsv_getField3D_i2: cannot specify a varName for VarsLevs mpi distribution')
      end if
      if (gsv_varExist(statevector,varName_opt)) then
        ilev1 = 1 + statevector%varOffset(vnl_varListIndex(varName_opt))
        ilev2 = ilev1 - 1 + statevector%varNumLev(vnl_varListIndex(varName_opt))
        if (present(stepIndex_opt)) then
          field3D(lon1:,lat1:,1:) => statevector%gd_i2(:,:,ilev1:ilev2,stepIndex_opt)
        else
          field3D(lon1:,lat1:,1:) => statevector%gd3d_i2(:,:,ilev1:ilev2)
        end if
      else
        call utl_abort('gsv_getField3D_i2: Unknown variable name! ' // varName_opt)
      end if
    else
      if (present(stepIndex_opt)) then
        field3D(lon1:,lat1:,k1:) => statevector%gd_i2(:,:,:,stepIndex_opt)
      else
        field3D(lon1:,lat1:,k1:) => statevector%gd3d_i2(:,:,:)
      end if
    end if

  end function gsv_getField3D_i2

  !--------------------------------------------------------------------------
  ! gsv_getFieldUV_r8
  !--------------------------------------------------------------------------
  function gsv_getFieldUV_r8(statevector) result(field)
    implicit none
    type(struct_gsv), intent(in)           :: statevector
    real(8),pointer                        :: field(:,:,:,:)
    integer                                :: lon1,lat1,k1

    lon1 = statevector%myLonBeg
    lat1 = statevector%myLatBeg
    k1 = statevector%mykBeg

    if (.not. associated(statevector%gd_r8)) call utl_abort('gsv_getFieldUV_r8: data with type r8 not allocated')

    field(lon1:,lat1:,k1:,1:) => statevector%gdUV_r8(:,:,:,:)

  end function gsv_getFieldUV_r8

  !--------------------------------------------------------------------------
  ! gsv_getFieldUV_r4
  !--------------------------------------------------------------------------
  function gsv_getFieldUV_r4(statevector) result(field)
    implicit none
    type(struct_gsv), intent(in)           :: statevector
    real(4),pointer                        :: field(:,:,:,:)
    integer                                :: lon1,lat1,k1

    lon1=statevector%myLonBeg
    lat1=statevector%myLatBeg
    k1=statevector%mykBeg

    if (.not. associated(statevector%gd_r4)) call utl_abort('gsv_getFieldUV_r4: data with type r4 not allocated')

    field(lon1:,lat1:,k1:,1:) => statevector%gdUV_r4(:,:,:,:)

  end function gsv_getFieldUV_r4

  !--------------------------------------------------------------------------
  ! gsv_getFieldUV_i2
  !--------------------------------------------------------------------------
  function gsv_getFieldUV_i2(statevector) result(field)
    implicit none
    type(struct_gsv), intent(in)           :: statevector
    integer(2),pointer                        :: field(:,:,:,:)
    integer                                :: lon1,lat1,k1

    lon1=statevector%myLonBeg
    lat1=statevector%myLatBeg
    k1=statevector%mykBeg

    if (.not. associated(statevector%gd_i2)) call utl_abort('gsv_getFieldUV_i2: data with type i2 not allocated')

    field(lon1:,lat1:,k1:,1:) => statevector%gdUV_i2(:,:,:,:)

  end function gsv_getFieldUV_i2

  !--------------------------------------------------------------------------
  ! gsv_getGZsfc
  !--------------------------------------------------------------------------
  function gsv_getGZsfc(statevector) result(field)
    implicit none
    type(struct_gsv), intent(in)           :: statevector
    real(8),pointer                        :: field(:,:)
    integer                                :: lon1,lat1

    lon1 = statevector%myLonBeg
    lat1 = statevector%myLatBeg

    if ( .not. statevector%gzSfcPresent ) call utl_abort('gsv_getGZsfc: data not allocated')

    if ( associated(statevector%gzSfc) ) then
      field(lon1:,lat1:) => statevector%gzSfc(:,:)
    else
      nullify(field)
    end if

  end function gsv_getGZsfc

  !--------------------------------------------------------------------------
  ! gsv_getIntOffset
  !--------------------------------------------------------------------------
  function gsv_getIntOffset(statevector,varName_opt) result(offset)
    implicit none
    type(struct_gsv), intent(in) :: statevector
    character(len=*), optional   :: varName_opt
    real(8), pointer             :: offset(:,:)
    integer                      :: k1, ilev1, ilev2

    k1=statevector%mykBeg

    if (.not. associated(statevector%intOffset)) call utl_abort('gsv_getIntOffset: intOffset not allocated')

    if (present(varName_opt)) then
      if (statevector%mpi_distribution == 'VarsLevs') then
        call utl_abort('gsv_getIntOffset: cannot specify a varName for VarsLevs mpi distribution')
      end if
      if (gsv_varExist(statevector,varName_opt)) then
        ilev1 = 1 + statevector%varOffset(vnl_varListIndex(varName_opt))
        ilev2 = ilev1 - 1 + statevector%varNumLev(vnl_varListIndex(varName_opt))
        offset(1:,1:) => statevector%intOffset(ilev1:ilev2,:)
      else
        call utl_abort('gsv_getIntOffset: Unknown variable name! ' // varName_opt)
      end if
    else
      offset(k1:,1:) => statevector%intOffset(:,:)
    end if

  end function gsv_getIntOffset

  !--------------------------------------------------------------------------
  ! gsv_getIntMultFactor
  !--------------------------------------------------------------------------
  function gsv_getIntMultFactor(statevector, varName_opt) result(multFactor)
    implicit none
    type(struct_gsv), intent(in) :: statevector
    character(len=*), optional   :: varName_opt
    real(8), pointer             :: multFactor(:,:)
    integer                      :: k1, ilev1, ilev2

    k1=statevector%mykBeg

    if (.not. associated(statevector%intMultFactor)) call utl_abort('gsv_getIntMultFactor: intMultFactor not allocated')

    if (present(varName_opt)) then
      if (statevector%mpi_distribution == 'VarsLevs') then
        call utl_abort('gsv_getIntMultFactor: cannot specify a varName for VarsLevs mpi distribution')
      end if
      if (gsv_varExist(statevector,varName_opt)) then
        ilev1 = 1 + statevector%varOffset(vnl_varListIndex(varName_opt))
        ilev2 = ilev1 - 1 + statevector%varNumLev(vnl_varListIndex(varName_opt))
        multFactor(1:,1:) => statevector%intMultFactor(ilev1:ilev2,:)
      else
        call utl_abort('gsv_getIntMultFactor: Unknown variable name! ' // varName_opt)
      end if
    else
      multFactor(k1:,1:) => statevector%intMultFactor(:,:)
    end if

  end function gsv_getIntMultFactor

  !--------------------------------------------------------------------------
  ! gsv_getDateStamp
  !--------------------------------------------------------------------------
  function gsv_getDateStamp(statevector,stepIndex_opt) result(dateStamp)
    implicit none
    type(struct_gsv), intent(in)   :: statevector
    integer, intent(in), optional  :: stepIndex_opt 
    integer                        :: dateStamp

    if (associated(statevector%dateStampList)) then
      if (present(stepIndex_opt)) then
        if (stepIndex_opt.gt.0.and.stepIndex_opt.le.statevector%numStep) then
          dateStamp=statevector%dateStampList(stepIndex_opt)
        else
          write(*,*) 'gsv_getDateStamp: requested step is out of range! Step,numStep=',stepIndex_opt,statevector%numStep
          call utl_abort('gsv_getDateStamp')
        end if    
      else
        dateStamp=statevector%dateStamp3D
      end if
    else
      call utl_abort('gsv_getDateStamp: dateStampList was not created during allocation!')
    end if

  end function gsv_getDateStamp

  !--------------------------------------------------------------------------
  ! gsv_getVco
  !--------------------------------------------------------------------------
  function gsv_getVco(statevector) result(vco_ptr)
    implicit none
    type(struct_gsv)          :: statevector
    type(struct_vco), pointer :: vco_ptr

    vco_ptr => statevector%vco

  end function gsv_getVco

  !--------------------------------------------------------------------------
  ! gsv_setVcO
  !--------------------------------------------------------------------------
  subroutine gsv_setVco(statevector,vco_ptr)
    implicit none
    type(struct_gsv) :: statevector
    type(struct_vco), pointer  :: vco_ptr

    statevector%vco => vco_ptr

  end subroutine gsv_setVco

  !--------------------------------------------------------------------------
  ! gsv_getHco
  !--------------------------------------------------------------------------
  function gsv_getHco(statevector) result(hco_ptr)
    implicit none
    type(struct_gsv)          :: statevector
    type(struct_hco), pointer :: hco_ptr

    hco_ptr => statevector%hco

  end function gsv_getHco

  !--------------------------------------------------------------------------
  ! gsv_setHco
  !--------------------------------------------------------------------------
  subroutine gsv_setHco(statevector,hco_ptr)
    implicit none
    type(struct_gsv)          :: statevector
    type(struct_hco), pointer :: hco_ptr

    statevector%hco => hco_ptr

  end subroutine gsv_setHco

  !--------------------------------------------------------------------------
  ! gsv_readFromFile
  !--------------------------------------------------------------------------
  subroutine gsv_readFromFile(statevector_out, fileName, etiket_in, typvar_in, stepIndex_opt,  &
                              unitConversion_opt, HUcontainsLQ_opt, PsfcReference_opt, readGZsfc_opt)
    implicit none
    ! Note this routine currently only works correctly for reading FULL FIELDS,
    ! not increments or perturbations... because of the HU -> LQ conversion

    ! arguments
    type(struct_gsv)              :: statevector_out

    character(len=*), intent(in)  :: fileName
    character(len=*), intent(in)  :: etiket_in
    character(len=*), intent(in)  :: typvar_in

    integer, optional             :: stepIndex_opt

    logical, optional             :: unitConversion_opt
    logical, optional             :: HUcontainsLQ_opt
    logical, optional             :: readGZsfc_opt

    real(8), optional             :: PsfcReference_opt(:,:)

    ! locals
    type(struct_gsv) :: statevector_file, statevector_tiles, statevector_hinterp, statevector_vinterp

    real(4), pointer     :: field3d_r4_ptr(:,:,:)
    real(8), pointer     :: field_in_ptr(:,:,:,:), field_out_ptr(:,:,:,:)
    real(8), allocatable :: PsfcReference3D(:,:,:)
    real(4) :: factor_r4

    integer :: kIndex, lonIndex, latIndex, varIndex, stepIndex, levIndex

    character(len=4) :: varName
    character(len=4), pointer :: varNamesToRead(:)

    logical :: doHorizInterp, doVertInterp, unitConversion, HUcontainsLQ
    logical :: readGZsfc, readSubsetOfLevels

    type(struct_vco), pointer :: vco_file
    type(struct_hco), pointer :: hco_file

    nullify(vco_file, hco_file)
    nullify(field3d_r4_ptr, field_in_ptr, field_out_ptr)

    write(*,*) ''
    write(*,*) 'gsv_readFromFile: START'
    call tmg_start(7,'READFROMFILE')

    if ( present(stepIndex_opt) ) then
      stepIndex = stepIndex_opt
    else
      stepIndex = statevector_out%anltime
    end if

    if ( present(unitConversion_opt) ) then
      unitConversion = unitConversion_opt
    else
      unitConversion = .true.
    end if

    if ( present(HUcontainsLQ_opt) ) then
      HUcontainsLQ = HUcontainsLQ_opt
    else
      HUcontainsLQ = .true.
    end if

    if ( present(readGZsfc_opt) ) then
      readGZsfc = readGZsfc_opt
    else
      readGZsfc = .false.
    end if

    nullify(varNamesToRead)
    call vnl_varNamesFromExistList(varNamesToRead, statevector_out%varExistlist(:))

    ! set up vertical and horizontal coordinate for input file
    call vco_SetupFromFile(vco_file,trim(fileName),beSilent_opt=.true.)
    readSubsetOfLevels = vco_subsetOrNot(statevector_out%vco, vco_file)
    if ( readSubsetOfLevels ) then
      ! use the output vertical grid provided to read only a subset of the verical levels
      write(*,*)
      write(*,*) 'gsv_readFromFile: read only a subset of the verical levels'
      call vco_deallocate(vco_file)
      vco_file => statevector_out%vco
    else
      write(*,*)
      write(*,*) 'gsv_readFromFile: all the vertical levels will be read'
    end if

    varName=gsv_getVarNameFromK(statevector_out,1)
    call hco_SetupFromFile(hco_file,trim(fileName), ' ',gridName_opt='FILEGRID',varName_opt=varName)

    ! test if horizontal and/or vertical interpolation needed for statevector grid
    if (readSubsetOfLevels) then
      doVertInterp = .false.
    else
      doVertInterp = .not.vco_equal(vco_file,statevector_out%vco)
    end if
    doHorizInterp = .not.hco_equal(hco_file,statevector_out%hco)
    write(*,*) 'gsv_readFromFile: doVertInterp = ', doVertInterp, ', doHorizInterp = ', doHorizInterp

    !-- 1.0 Read the file, distributed over mpi task with respect to variables/levels

    ! initialize single precision 3D working copy of statevector for reading file
    call gsv_allocate(statevector_file, 1, hco_file, vco_file,                &
                      dateStamp_opt=statevector_out%datestamplist(stepIndex), &
                      mpi_local_opt=.true., mpi_distribution_opt='VarsLevs',  &
                      dataKind_opt=4, allocGZsfc_opt=readGZsfc,               &
                      varNames_opt=varNamesToRead)

    call gsv_readFile(statevector_file, filename, etiket_in, typvar_in,  &
                      readGZsfc_opt=readGZsfc)

    !-- 2.0 Horizontal Interpolation

    ! initialize single precision 3D working copy of statevector for horizontal interpolation result
    call gsv_allocate(statevector_hinterp, 1, statevector_out%hco, vco_file,  &
                      dateStamp_opt=statevector_out%datestamplist(stepIndex), &
                      mpi_local_opt=.true., mpi_distribution_opt='VarsLevs',  &
                      dataKind_opt=4, allocGZsfc_opt=readGZsfc,               &
                      varNames_opt=varNamesToRead)

    call gsv_hInterpolate_r4(statevector_file, statevector_hinterp)

    call gsv_deallocate(statevector_file)

    !-- 3.0 Scale factor and unit conversion

    field3d_r4_ptr => gsv_getField3D_r4(statevector_hinterp)
    K_LOOP: do kIndex = statevector_hinterp%mykBeg, statevector_hinterp%mykEnd
      varName = gsv_getVarNameFromK(statevector_hinterp,kIndex)

      if ( unitConversion ) then
        if ( trim(varName) == 'UU' .or. trim(varName) == 'VV') then
          factor_r4 = MPC_M_PER_S_PER_KNOT_R4 ! knots -> m/s
        else if ( trim(varName) == 'P0' ) then
          factor_r4 = 100.0 ! hPa -> Pa
        else
          factor_r4 = 1.0 ! no conversion
        end if
      else
        factor_r4 = 1.0
      end if
      field3d_r4_ptr(:,:,kIndex) = factor_r4 * field3d_r4_ptr(:,:,kIndex)

    end do K_LOOP

    !-- 4.0 MPI communication from vars/levels to lat/lon tiles

    ! initialize double precision 3D working copy of statevector for mpi communication result
    call gsv_allocate(statevector_tiles, 1, statevector_out%hco, vco_file,    &
                      dateStamp_opt=statevector_out%datestamplist(stepIndex), &
                      mpi_local_opt=.true., mpi_distribution_opt='Tiles',     &
                      dataKind_opt=8, allocGZsfc_opt=readGZsfc,               &
                      varNames_opt=varNamesToRead)

    call transposeVarsLevsToLatLon(statevector_hinterp, statevector_tiles)

    call gsv_deallocate(statevector_hinterp)

    !-- 5.0 Vertical interpolation

    ! initialize double precision 3D working copy of statevector for mpi communication result
    call gsv_allocate(statevector_vinterp, 1, statevector_out%hco, statevector_out%vco, &
                      dateStamp_opt=statevector_out%datestamplist(stepIndex), &
                      mpi_local_opt=.true., mpi_distribution_opt='Tiles', dataKind_opt=8, &
                      allocGZsfc_opt=readGZsfc, varNames_opt=varNamesToRead)

    if (present(PsfcReference_opt) ) then
      allocate(PsfcReference3D(statevector_tiles%myLonBeg:statevector_tiles%myLonEnd,statevector_tiles%myLatBeg:statevector_tiles%myLatEnd,1))
      PsfcReference3D(:,:,1) = PsfcReference_opt(:,:)
      call gsv_vInterpolate(statevector_tiles,statevector_vinterp,PsfcReference_opt=PsfcReference3D)
      deallocate(PsfcReference3D)
    else
      call gsv_vInterpolate(statevector_tiles,statevector_vinterp)
    end if

    call gsv_deallocate(statevector_tiles)

    !-- 6.0 Copy result to output statevector and convert HU to LQ
    ! THIS SHOULD BE REPLACED MOSTLY BY CALL TO gsv_copy

    if ( associated(statevector_vinterp%gzSfc) .and. associated(statevector_out%gzSfc) ) then
      statevector_out%gzSfc(:,:) = statevector_vinterp%gzSfc(:,:)
    end if

    do varIndex = 1, vnl_numvarmax
      if ( .not. gsv_varExist(statevector_out,vnl_varNameList(varIndex)) ) cycle
      field_in_ptr => gsv_getField_r8(statevector_vinterp, vnl_varNameList(varIndex))
      field_out_ptr => gsv_getField_r8(statevector_out, vnl_varNameList(varIndex))

      if ( trim(vnl_varNameList(varIndex)) == 'HU' .and. HUcontainsLQ ) then
!$OMP PARALLEL DO PRIVATE (latIndex,levIndex,lonIndex)
        do latIndex = statevector_out%myLatBeg, statevector_out%myLatEnd
          do levIndex = 1, statevector_out%varNumLev(varIndex)
            do lonIndex = statevector_out%myLonBeg, statevector_out%myLonEnd
              field_out_ptr(lonIndex,latIndex,levIndex,stepIndex) = log(max(field_in_ptr(lonIndex,latIndex,levIndex,1),MPC_MINIMUM_HU_R8))
            end do
          end do
        end do
!$OMP END PARALLEL DO
      else
!$OMP PARALLEL DO PRIVATE (latIndex,levIndex,lonIndex)
        do latIndex = statevector_out%myLatBeg, statevector_out%myLatEnd
          do levIndex = 1, statevector_out%varNumLev(varIndex)
            do lonIndex = statevector_out%myLonBeg, statevector_out%myLonEnd
              field_out_ptr(lonIndex,latIndex,levIndex,stepIndex) = field_in_ptr(lonIndex,latIndex,levIndex,1)
            end do
          end do
        end do
!$OMP END PARALLEL DO
      end if

    end do ! varIndex

    call gsv_deallocate(statevector_vinterp)
    deallocate(varNamesToRead)

    call tmg_stop(7)
    write(*,*) 'gsv_readFromFile: END'

  end subroutine gsv_readFromFile

  !--------------------------------------------------------------------------
  ! gsv_readFile
  !--------------------------------------------------------------------------
  subroutine gsv_readFile(statevector, filename, etiket_in, typvar_in, readGZsfc_opt)
    implicit none

    ! arguments
    type(struct_gsv)              :: statevector
    character(len=*), intent(in)  :: fileName
    character(len=*), intent(in)  :: etiket_in
    character(len=*), intent(in)  :: typvar_in

    logical, optional             :: readGZsfc_opt

    ! locals
    integer :: nulfile, ierr, ip1, ni_file, nj_file, nk_file, kIndex, stepIndex, ikey, levIndex
    integer :: ni_var, nj_var, nk_var
    integer :: fnom, fstouv, fclos, fstfrm, fstlir, fstinf
    integer :: fstprm, EZscintID_var, ezdefset, ezqkdef

    integer :: dateo_var, deet_var, npas_var, nbits_var, datyp_var
    integer :: ip1_var, ip2_var, ip3_var, swa_var, lng_var, dltf_var, ubc_var
    integer :: extra1_var, extra2_var, extra3_var
    integer :: ig1_var, ig2_var, ig3_var, ig4_var

    character(len=4 ) :: nomvar_var
    character(len=2 ) :: typvar_var
    character(len=1 ) :: grtyp_var
    character(len=12) :: etiket_var

    real(4), pointer :: field_r4_ptr(:,:,:,:)
    real(4), pointer :: fieldUV_r4_ptr(:,:,:,:)
    real(4), pointer :: gd2d_file_r4(:,:)
    real(4), allocatable :: gd2d_var_r4(:,:)

    character(len=4)  :: varName
    character(len=2)  :: varLevel

    type(struct_vco), pointer :: vco_file
    type(struct_hco), pointer :: hco_file

    nullify(gd2d_file_r4)

    vco_file => gsv_getVco(statevector)

    if ( statevector%mpi_distribution /= 'VarsLevs' .and. &
         statevector%mpi_local ) then
      call utl_abort('gsv_readFile: statevector must have ' //   &
                     'complete horizontal fields on each mpi task.')
    end if

    !- Open input field
    nulfile = 0
    write(*,*) 'gsv_readFile: file name = ',trim(fileName)
    ierr = fnom(nulfile,trim(fileName),'RND+OLD+R/O',0)
       
    if ( ierr >= 0 ) then
      write(*,*)'gsv_readFile: file opened with unit number ',nulfile
      ierr  =  fstouv(nulfile,'RND+OLD')
    else
      call utl_abort('gsv_readFile: problem opening input file, aborting!')
    end if

    if (nulfile == 0 ) then
      call utl_abort('gsv_readFile: unit number for input file not valid!')
    end if

    ! Read surface GZ if requested
    if ( present(readGZsfc_opt) ) then
      if ( readGZsfc_opt .and. associated( statevector%GZsfc ) ) then
        write(*,*) 'gsv_readFile: reading the surface GZ'
        varName = 'GZ'
        ip1 = statevector%vco%ip1_sfc
        ikey = fstinf(nulfile, ni_file, nj_file, nk_file,  &
                      -1, etiket_in, &
                      -1, -1, -1, typvar_in, varName)
        allocate(gd2d_file_r4(ni_file,nj_file))
        gd2d_file_r4(:,:) = 0.0d0
        ierr=fstlir(gd2d_file_r4(:,:),nulfile,ni_file, nj_file, nk_file,  &
                    -1,etiket_in,ip1,-1,-1,  &
                    typvar_in,varName)
        if (ierr.lt.0)then
          write(*,*) varName,ip1
          call utl_abort('gsv_readFile: Problem with reading surface GZ from file')
        end if
        statevector%GZsfc(:,:) = real(gd2d_file_r4(1:statevector%hco%ni,1:statevector%hco%nj),8)*10.0d0*RG
        deallocate(gd2d_file_r4)
        nullify(gd2d_file_r4)
      end if
    end if

    nullify(hco_file)
    if (statevector%hco%global) then
      ! In global mode, allow for possibility that input is Z grid equivalent to output Gaussian grid
      if ( statevector%mykCount > 0 ) then
        call hco_SetupFromFile(hco_file, filename, ' ', 'INPUTFILE')
      end if 
    else
      ! In LAM mode, force the input file dimensions to be always identical to the input statevector dimensions
      hco_file => statevector%hco
      ni_file=statevector%ni
      nj_file=statevector%nj
    end if
    allocate(gd2d_file_r4(hco_file%ni,hco_file%nj))
    gd2d_file_r4(:,:) = 0.0

    ! Read all other fields needed for this MPI task
    field_r4_ptr => gsv_getField_r4(statevector)
    if (statevector%mpi_distribution == 'VarsLevs') then
      fieldUV_r4_ptr => gsv_getFieldUV_r4(statevector)
    end if
    do stepIndex = 1, statevector%numStep
      K_LOOP: do kIndex = statevector%mykBeg, statevector%mykEnd
        varName = gsv_getVarNameFromK(statevector,kIndex)
        levIndex = gsv_getLevFromK(statevector,kIndex)

        if (.not.gsv_varExist(statevector,varName)) cycle K_LOOP

        ! do not try to read diagnostic variables
        if ( trim(vnl_varTypeFromVarname(varName)) == 'DIAG') cycle K_LOOP

        varLevel = vnl_varLevelFromVarname(varName)
        if (varLevel == 'MM') then
          ip1 = vco_file%ip1_M(levIndex)
        elseif (varLevel == 'TH') then
          ip1 = vco_file%ip1_T(levIndex)
        elseif (varLevel == 'SF') then
          ip1 = -1
        else
          call utl_abort('gsv_readFile: unknown varLevel')
        end if

        ! Make sure that the input variable has the same grid size than hco_file   
        ikey = fstinf(nulfile, ni_var, nj_var, nk_var,         &
                      statevector%datestamplist(1), etiket_in, &
                      -1, -1, -1, typvar_in, varName)

        if ( ni_var == hco_file%ni .and. nj_var == hco_file%nj ) then
          ierr=fstlir(gd2d_file_r4(:,:),nulfile,ni_file, nj_file, nk_file,  &
                     statevector%datestamplist(stepIndex),etiket_in,ip1,-1,-1,  &
                     typvar_in,varName)
        else
          ! Special cases for variables that are on a different horizontal grid in LAM (e.g. TG)
          write(*,*)
          write(*,*) 'gsv_readFile: variable on a different horizontal grid = ',trim(varName)
          if (statevector%hco%global) then
            call utl_abort('gsv_readFile: This is not allowed in global mode!')
          end if
          ierr = fstprm( ikey,                                                               & ! IN
                         dateo_var, deet_var, npas_var, ni_var, nj_var, nk_var, nbits_var,   & ! OUT
                         datyp_var, ip1_var, ip2_var, ip3_var, typvar_var, nomvar_var,       & ! OUT
                         etiket_var, grtyp_var, ig1_var, ig2_var, ig3_var, ig4_var, swa_var, & ! OUT
                         lng_var, dltf_var, ubc_var, extra1_var, extra2_var, extra3_var )      ! OUT

          EZscintID_var  = ezqkdef( ni_var, nj_var, grtyp_var, ig1_var, ig2_var, ig3_var, ig4_var, nulfile ) ! IN

          allocate(gd2d_var_r4(ni_var,nj_var))
          gd2d_var_r4(:,:) = 0.0

          ierr=fstlir(gd2d_var_r4(:,:),nulfile,ni_var, nj_var, nk_var,  &
                     statevector%datestamplist(stepIndex),etiket_in,ip1,-1,-1,  &
                     typvar_in,varName)

          ierr = ezdefset(hco_file%EZscintID,EZscintID_var)
          ierr = utl_ezsint( gd2d_file_r4, gd2d_var_r4, interpDegree_opt='NEAREST', extrapDegree_opt='NEUTRAL' )

          deallocate(gd2d_var_r4)
        end if

        field_r4_ptr(:,:,kIndex,stepIndex) = gd2d_file_r4(1:statevector%hco%ni,1:statevector%hco%nj)

        if (ierr.lt.0)then
          if (varName == 'HU') then
            ! HU variable not found in file, try reading LQ
            ierr=fstlir(gd2d_file_r4(:,:),nulfile, ni_file, nj_file, nk_file,  &
                        statevector%datestamplist(stepIndex),etiket_in,ip1,-1,-1,  &
                        typvar_in,'LQ')
            field_r4_ptr(:,:,kIndex,stepIndex) = gd2d_file_r4(1:statevector%hco%ni,1:statevector%hco%nj)
            if (ierr.lt.0)then
              write(*,*) 'LQ',ip1,statevector%datestamplist(stepIndex)
              call utl_abort('gsv_readFile: Problem with reading file')
            end if
          else
            write(*,*) varName,ip1,statevector%datestamplist(stepIndex)
            call utl_abort('gsv_readFile: Problem with reading file')
          end if
        end if

        ! When mpi distribution could put UU on a different mpi task than VV
        ! here we re-read the corresponding UV component and store it
        if (statevector%mpi_distribution == 'VarsLevs') then
          if (varName == 'UU') then
            ierr=fstlir(gd2d_file_r4(:,:),nulfile, ni_file, nj_file, nk_file,  &
                        statevector%datestamplist(stepIndex),etiket_in,ip1,-1,-1,  &
                        typvar_in,'VV')
            fieldUV_r4_ptr(:,:,kIndex,stepIndex) = gd2d_file_r4(1:statevector%hco%ni,1:statevector%hco%nj)
          elseif (varName == 'VV') then
            ierr=fstlir(gd2d_file_r4(:,:),nulfile, ni_file, nj_file, nk_file,  &
                        statevector%datestamplist(stepIndex),etiket_in,ip1,-1,-1,  &
                        typvar_in,'UU')
            fieldUV_r4_ptr(:,:,kIndex,stepIndex) = gd2d_file_r4(1:statevector%hco%ni,1:statevector%hco%nj)
          end if
        end if

      end do K_LOOP
    end do

    if (statevector%hco%global .and. statevector%mykCount > 0) call hco_deallocate(hco_file)

    ierr = fstfrm(nulfile)
    ierr = fclos(nulfile)        
    if ( associated(gd2d_file_r4) ) deallocate(gd2d_file_r4)

  end subroutine gsv_readFile

  !--------------------------------------------------------------------------
  ! transposeVarsLevsToLatLon
  !--------------------------------------------------------------------------
  subroutine transposeVarsLevsToLatLon(statevector_in, statevector_out)
    ! Transpose the data from mpi_distribution=VarsLevs to Tiles
    implicit none

    ! arguments
    type(struct_gsv) :: statevector_in
    type(struct_gsv) :: statevector_out

    ! locals
    integer :: youridx, youridy, yourid, nsize, maxkcount, ierr
    integer :: sendrecvKind, inKind, outKind
    integer, allocatable :: displs(:), nsizes(:)
    real(4), allocatable :: gd_send_r4(:,:,:,:,:), gd_recv_r4(:,:,:,:,:)
    real(8), allocatable :: gd_send_r8(:,:,:,:,:), gd_recv_r8(:,:,:,:,:)
    real(4), pointer     :: field_in_r4_ptr(:,:,:,:), field_out_r4_ptr(:,:,:,:)
    real(8), pointer     :: field_in_r8_ptr(:,:,:,:), field_out_r8_ptr(:,:,:,:)
    real(8), allocatable :: gd_send_GZ(:,:,:), gd_recv_GZ(:,:)
    real(8), pointer     :: field_GZ_in_ptr(:,:), field_GZ_out_ptr(:,:)

    if ( statevector_in%mpi_distribution /= 'VarsLevs' ) then
      call utl_abort('transposeVarsLevsToLatLon: input statevector must have VarsLevs mpi distribution') 
    end if

    if ( statevector_out%mpi_distribution /= 'Tiles' ) then
      call utl_abort('transposeVarsLevsToLatLon: out statevector must have Tiles mpi distribution')
    end if

    inKind = statevector_in%dataKind
    outKind = statevector_out%dataKind
    if ( inKind == 4 .or. outKind == 4 ) then
      sendrecvKind = 4
    else
      sendrecvKind = 8
    end if

    maxkCount = maxval(statevector_in%allkCount(:))
    if( sendrecvKind == 4 ) then
      allocate(gd_send_r4(statevector_out%lonPerPEmax,statevector_out%latPerPEmax,maxkCount,statevector_out%numStep,mpi_nprocs))
      allocate(gd_recv_r4(statevector_out%lonPerPEmax,statevector_out%latPerPEmax,maxkCount,statevector_out%numStep,mpi_nprocs))
      gd_send_r4(:,:,:,:,:) = 0.0
      gd_recv_r4(:,:,:,:,:) = 0.0
    else
      allocate(gd_send_r8(statevector_out%lonPerPEmax,statevector_out%latPerPEmax,maxkCount,statevector_out%numStep,mpi_nprocs))
      allocate(gd_recv_r8(statevector_out%lonPerPEmax,statevector_out%latPerPEmax,maxkCount,statevector_out%numStep,mpi_nprocs))
      gd_send_r8(:,:,:,:,:) = 0.0
      gd_recv_r8(:,:,:,:,:) = 0.0
    end if

    if ( sendrecvKind == 4 .and. inKind == 4 ) then
      field_in_r4_ptr => gsv_getField_r4(statevector_in)
!$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
      do youridy = 0, (mpi_npey-1)
        do youridx = 0, (mpi_npex-1)
          yourid = youridx + youridy*mpi_npex
          gd_send_r4(1:statevector_out%allLonPerPE(youridx+1),  &
                     1:statevector_out%allLatPerPE(youridy+1),  &
                     1:statevector_in%mykCount,:,yourid+1) =  &
            field_in_r4_ptr(statevector_out%allLonBeg(youridx+1):statevector_out%allLonEnd(youridx+1),  &
                            statevector_out%allLatBeg(youridy+1):statevector_out%allLatEnd(youridy+1),  &
                            statevector_in%mykBeg:statevector_in%mykEnd,:)
        end do
      end do
!$OMP END PARALLEL DO
    elseif ( sendrecvKind == 4 .and. inKind == 8 ) then
      field_in_r8_ptr => gsv_getField_r8(statevector_in)
!$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
      do youridy = 0, (mpi_npey-1)
        do youridx = 0, (mpi_npex-1)
          yourid = youridx + youridy*mpi_npex
          gd_send_r4(1:statevector_out%allLonPerPE(youridx+1),  &
                     1:statevector_out%allLatPerPE(youridy+1),  &
                     1:statevector_in%mykCount,:,yourid+1) =  &
            real(field_in_r8_ptr(statevector_out%allLonBeg(youridx+1):statevector_out%allLonEnd(youridx+1),  &
                                 statevector_out%allLatBeg(youridy+1):statevector_out%allLatEnd(youridy+1),  &
                                 statevector_in%mykBeg:statevector_in%mykEnd,:),4)
        end do
      end do
!$OMP END PARALLEL DO
    elseif ( sendrecvKind == 8 .and. inKind == 4 ) then
      field_in_r4_ptr => gsv_getField_r4(statevector_in)
!$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
      do youridy = 0, (mpi_npey-1)
        do youridx = 0, (mpi_npex-1)
          yourid = youridx + youridy*mpi_npex
          gd_send_r8(1:statevector_out%allLonPerPE(youridx+1),  &
                     1:statevector_out%allLatPerPE(youridy+1),  &
                     1:statevector_in%mykCount,:,yourid+1) =  &
            real(field_in_r4_ptr(statevector_out%allLonBeg(youridx+1):statevector_out%allLonEnd(youridx+1),  &
                                 statevector_out%allLatBeg(youridy+1):statevector_out%allLatEnd(youridy+1),  &
                                 statevector_in%mykBeg:statevector_in%mykEnd,:),8)
        end do
      end do
!$OMP END PARALLEL DO
    elseif ( sendrecvKind == 8 .and. inKind == 8 ) then
      field_in_r8_ptr => gsv_getField_r8(statevector_in)
!$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
      do youridy = 0, (mpi_npey-1)
        do youridx = 0, (mpi_npex-1)
          yourid = youridx + youridy*mpi_npex
          gd_send_r8(1:statevector_out%allLonPerPE(youridx+1),  &
                     1:statevector_out%allLatPerPE(youridy+1),  &
                     1:statevector_in%mykCount,:,yourid+1) =  &
            field_in_r8_ptr(statevector_out%allLonBeg(youridx+1):statevector_out%allLonEnd(youridx+1),  &
                            statevector_out%allLatBeg(youridy+1):statevector_out%allLatEnd(youridy+1),  &
                            statevector_in%mykBeg:statevector_in%mykEnd,:)
        end do
      end do
!$OMP END PARALLEL DO
    end if

    nsize = statevector_out%lonPerPEmax * statevector_out%latPerPEmax * statevector_out%numStep * maxkCount
    if(mpi_nprocs > 1) then
      if( sendrecvKind == 4 ) then
        call rpn_comm_alltoall(gd_send_r4, nsize, "mpi_real4",  &
                               gd_recv_r4, nsize, "mpi_real4", "GRID", ierr)
      else
        call rpn_comm_alltoall(gd_send_r8, nsize, "mpi_real8",  &
                               gd_recv_r8, nsize, "mpi_real8", "GRID", ierr)
      end if
    else
      if ( sendrecvKind == 4 ) then
         gd_recv_r4(:,:,:,:,1) = gd_send_r4(:,:,:,:,1)
      else
         gd_recv_r8(:,:,:,:,1) = gd_send_r8(:,:,:,:,1)
      end if
    end if

    if ( sendrecvKind == 4 .and. outKind == 4 ) then
      field_out_r4_ptr => gsv_getField_r4(statevector_out)
!$OMP PARALLEL DO PRIVATE(yourid)
      do yourid = 0, (mpi_nprocs-1)
        field_out_r4_ptr(statevector_out%myLonBeg:statevector_out%myLonEnd, &
                         statevector_out%myLatBeg:statevector_out%myLatEnd, &
                         statevector_in%allkBeg(yourid+1):statevector_in%allkEnd(yourid+1), :) =   &
          gd_recv_r4(1:statevector_out%lonPerPE,  &
                     1:statevector_out%latPerPE,  &
                     1:statevector_in%allkCount(yourid+1), :, yourid+1)
      end do
!$OMP END PARALLEL DO
    elseif ( sendrecvKind == 4 .and. outKind == 8 ) then
      field_out_r8_ptr => gsv_getField_r8(statevector_out)
!$OMP PARALLEL DO PRIVATE(yourid)
      do yourid = 0, (mpi_nprocs-1)
        field_out_r8_ptr(statevector_out%myLonBeg:statevector_out%myLonEnd, &
                         statevector_out%myLatBeg:statevector_out%myLatEnd, &
                         statevector_in%allkBeg(yourid+1):statevector_in%allkEnd(yourid+1), :) =   &
          real(gd_recv_r4(1:statevector_out%lonPerPE,  &
                          1:statevector_out%latPerPE,  &
                          1:statevector_in%allkCount(yourid+1), :, yourid+1),8)
      end do
!$OMP END PARALLEL DO
    elseif ( sendrecvKind == 8 .and. outKind == 4 ) then
      field_out_r4_ptr => gsv_getField_r4(statevector_out)
!$OMP PARALLEL DO PRIVATE(yourid)
      do yourid = 0, (mpi_nprocs-1)
        field_out_r4_ptr(statevector_out%myLonBeg:statevector_out%myLonEnd, &
                         statevector_out%myLatBeg:statevector_out%myLatEnd, &
                         statevector_in%allkBeg(yourid+1):statevector_in%allkEnd(yourid+1), :) =   &
          real(gd_recv_r8(1:statevector_out%lonPerPE,  &
                          1:statevector_out%latPerPE,  &
                          1:statevector_in%allkCount(yourid+1), :, yourid+1),4)
      end do
!$OMP END PARALLEL DO
    elseif ( sendrecvKind == 8 .and. outKind == 8 ) then
      field_out_r8_ptr => gsv_getField_r8(statevector_out)
!$OMP PARALLEL DO PRIVATE(yourid)
      do yourid = 0, (mpi_nprocs-1)
        field_out_r8_ptr(statevector_out%myLonBeg:statevector_out%myLonEnd, &
                         statevector_out%myLatBeg:statevector_out%myLatEnd, &
                         statevector_in%allkBeg(yourid+1):statevector_in%allkEnd(yourid+1),:) =   &
          gd_recv_r8(1:statevector_out%lonPerPE,  &
                     1:statevector_out%latPerPE,  &
                     1:statevector_in%allkCount(yourid+1), :, yourid+1)
      end do
!$OMP END PARALLEL DO
    end if

    if ( sendrecvKind == 4 ) then
      deallocate(gd_send_r4)
      deallocate(gd_recv_r4)
    else
      deallocate(gd_send_r8)
      deallocate(gd_recv_r8)
    end if

    if ( statevector_in%gzSfcPresent .and. statevector_out%gzSfcPresent ) then
      allocate(gd_send_GZ(statevector_out%lonPerPEmax,statevector_out%latPerPEmax,mpi_nprocs))
      allocate(gd_recv_GZ(statevector_out%lonPerPEmax,statevector_out%latPerPEmax))
      gd_send_GZ(:,:,:) = 0.0d0
      gd_recv_GZ(:,:) = 0.0d0
      field_GZ_in_ptr => gsv_getGZsfc(statevector_in)
      field_GZ_out_ptr => gsv_getGZsfc(statevector_out)

      if ( mpi_myid == 0 ) then
!$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
        do youridy = 0, (mpi_npey-1)
          do youridx = 0, (mpi_npex-1)
            yourid = youridx + youridy*mpi_npex
            gd_send_GZ(1:statevector_out%allLonPerPE(youridx+1),  &
                       1:statevector_out%allLatPerPE(youridy+1), yourid+1) =  &
              field_GZ_in_ptr(statevector_out%allLonBeg(youridx+1):statevector_out%allLonEnd(youridx+1),  &
                              statevector_out%allLatBeg(youridy+1):statevector_out%allLatEnd(youridy+1))
          end do
        end do
!$OMP END PARALLEL DO
      end if

      nsize = statevector_out%lonPerPEmax * statevector_out%latPerPEmax
      allocate(displs(mpi_nprocs))
      allocate(nsizes(mpi_nprocs))
      do yourid = 0, (mpi_nprocs-1)
        displs(yourid+1) = yourid*nsize
        nsizes(yourid+1) = nsize
      end do
      call rpn_comm_scatterv(gd_send_GZ, nsizes, displs, 'mpi_double_precision', &
                             gd_recv_GZ, nsize, 'mpi_double_precision', &
                             0, 'GRID', ierr)

      field_GZ_out_ptr(statevector_out%myLonBeg:statevector_out%myLonEnd, &
                       statevector_out%myLatBeg:statevector_out%myLatEnd) =   &
        gd_recv_GZ(1:statevector_out%lonPerPE,  &
                   1:statevector_out%latPerPE)

      deallocate(displs)
      deallocate(nsizes)
      deallocate(gd_recv_GZ)
      deallocate(gd_send_GZ)
    end if

  end subroutine transposeVarsLevsToLatLon

  !--------------------------------------------------------------------------
  ! transposeLatLonToVarsLevs
  !--------------------------------------------------------------------------
  subroutine transposeLatLonToVarsLevs(statevector_in, statevector_out)
    ! Transpose the data from mpi_distribution=Tiles to VarsLevs
    implicit none

    ! arguments
    type(struct_gsv) :: statevector_in
    type(struct_gsv) :: statevector_out

    ! locals
    integer :: youridx, youridy, yourid, nsize, maxkcount, ierr
    integer :: sendrecvKind, inKind, outKind
    real(4), allocatable :: gd_send_r4(:,:,:,:,:), gd_recv_r4(:,:,:,:,:)
    real(8), allocatable :: gd_send_r8(:,:,:,:,:), gd_recv_r8(:,:,:,:,:)
    real(4), pointer     :: field_in_r4_ptr(:,:,:,:), field_out_r4_ptr(:,:,:,:)
    real(8), pointer     :: field_in_r8_ptr(:,:,:,:), field_out_r8_ptr(:,:,:,:)
    real(8), pointer     :: field_GZ_in_ptr(:,:), field_GZ_out_ptr(:,:)
    real(8), allocatable :: gd_recv_GZ(:,:,:)

    if ( statevector_in%mpi_distribution /= 'Tiles' ) then
      call utl_abort('transposeLatLonToVarsLevs: input statevector must have Tiles mpi distribution') 
    end if

    if ( statevector_out%mpi_distribution /= 'VarsLevs' ) then
      call utl_abort('transposeLatLonToVarsLevs: out statevector must have VarsLevs mpi distribution')
    end if

    inKind = statevector_in%dataKind
    outKind = statevector_out%dataKind
    if ( inKind == 4 .or. outKind == 4 ) then
      sendrecvKind = 4
    else
      sendrecvKind = 8
    end if

    maxkCount = maxval(statevector_out%allkCount(:))
    if( sendrecvKind == 4 ) then
      allocate(gd_send_r4(statevector_in%lonPerPEmax,statevector_in%latPerPEmax,maxkCount,statevector_in%numStep,mpi_nprocs))
      allocate(gd_recv_r4(statevector_in%lonPerPEmax,statevector_in%latPerPEmax,maxkCount,statevector_in%numStep,mpi_nprocs))
      gd_send_r4(:,:,:,:,:) = 0.0
      gd_recv_r4(:,:,:,:,:) = 0.0
    else
      allocate(gd_send_r8(statevector_in%lonPerPEmax,statevector_in%latPerPEmax,maxkCount,statevector_in%numStep,mpi_nprocs))
      allocate(gd_recv_r8(statevector_in%lonPerPEmax,statevector_in%latPerPEmax,maxkCount,statevector_in%numStep,mpi_nprocs))
      gd_send_r8(:,:,:,:,:) = 0.0
      gd_recv_r8(:,:,:,:,:) = 0.0
    end if

    if ( sendrecvKind == 8 .and. inKind == 8 ) then
      field_in_r8_ptr => gsv_getField_r8(statevector_in)
!$OMP PARALLEL DO PRIVATE(yourid)
      do yourid = 0, (mpi_nprocs-1)
        gd_send_r8(1:statevector_in%lonPerPE, &
                   1:statevector_in%latPerPE, &
                   1:statevector_out%allkCount(yourid+1), :, yourid+1) =  &
            field_in_r8_ptr(statevector_in%myLonBeg:statevector_in%myLonEnd, &
                            statevector_in%myLatBeg:statevector_in%myLatEnd, &
                            statevector_out%allkBeg(yourid+1):statevector_out%allkEnd(yourid+1), :)
      end do
!$OMP END PARALLEL DO
    else
      call utl_abort('transposeLatLonToLevsVars: not compatible yet with these data types')
    end if

    nsize = statevector_in%lonPerPEmax * statevector_in%latPerPEmax * statevector_in%numStep * maxkCount
    if(mpi_nprocs > 1) then
      if( sendrecvKind == 4 ) then
        call rpn_comm_alltoall(gd_send_r4, nsize, "mpi_real4",  &
                               gd_recv_r4, nsize, "mpi_real4", "GRID", ierr)
      else
        call rpn_comm_alltoall(gd_send_r8, nsize, "mpi_real8",  &
                               gd_recv_r8, nsize, "mpi_real8", "GRID", ierr)
      end if
    else
      if ( sendrecvKind == 4 ) then
         gd_recv_r4(:,:,:,:,1) = gd_send_r4(:,:,:,:,1)
      else
         gd_recv_r8(:,:,:,:,1) = gd_send_r8(:,:,:,:,1)
      end if
    end if

    if ( sendrecvKind == 8 .and. outKind == 8 ) then
      field_out_r8_ptr => gsv_getField_r8(statevector_out)
!$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
      do youridy = 0, (mpi_npey-1)
        do youridx = 0, (mpi_npex-1)
          yourid = youridx + youridy*mpi_npex
          field_out_r8_ptr(statevector_in%allLonBeg(youridx+1):statevector_in%allLonEnd(youridx+1),  &
                           statevector_in%allLatBeg(youridy+1):statevector_in%allLatEnd(youridy+1),  &
                           statevector_out%mykBeg:statevector_out%mykEnd, :) = &
              gd_recv_r8(1:statevector_in%allLonPerPE(youridx+1),  &
                         1:statevector_in%allLatPerPE(youridy+1),  &
                         1:statevector_out%mykCount, :, yourid+1)
        end do
      end do
!$OMP END PARALLEL DO
    else
      call utl_abort('transposeLatLonToLevsVars: not compatible yet with these data types')
    end if

    if ( sendrecvKind == 4 ) then
      deallocate(gd_send_r4)
      deallocate(gd_recv_r4)
    else
      deallocate(gd_send_r8)
      deallocate(gd_recv_r8)
    end if

    if ( statevector_in%gzSfcPresent .and. statevector_out%gzSfcPresent ) then
      allocate(gd_recv_GZ(statevector_out%lonPerPEmax,statevector_out%latPerPEmax,mpi_nprocs))
      field_GZ_in_ptr => gsv_getGZsfc(statevector_in)
      field_GZ_out_ptr => gsv_getGZsfc(statevector_out)

      nsize = statevector_out%lonPerPEmax * statevector_out%latPerPEmax
      call rpn_comm_gather(field_GZ_in_ptr, nsize, 'mpi_double_precision',  &
                           gd_recv_GZ,      nsize, 'mpi_double_precision', 0, 'GRID', ierr )

      if ( mpi_myid == 0 ) then
!$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
        do youridy = 0, (mpi_npey-1)
          do youridx = 0, (mpi_npex-1)
            yourid = youridx + youridy*mpi_npex
            field_GZ_out_ptr(statevector_in%allLonBeg(youridx+1):statevector_in%allLonEnd(youridx+1),  &
                             statevector_in%allLatBeg(youridy+1):statevector_in%allLatEnd(youridy+1)) = &
                gd_recv_GZ(1:statevector_in%allLonPerPE(youridx+1),  &
                           1:statevector_in%allLatPerPE(youridy+1), yourid+1)
          end do
        end do
!$OMP END PARALLEL DO
      end if

      deallocate(gd_recv_GZ)
    end if

  end subroutine transposeLatLonToVarsLevs

  !--------------------------------------------------------------------------
  ! gsv_hInterpolate
  !--------------------------------------------------------------------------
  subroutine gsv_hInterpolate(statevector_in,statevector_out)
    ! s/r gsv_hInterpolate  - Horizontal interpolation of pressure defined fields
    implicit none

    ! arguments
    type(struct_gsv) :: statevector_in
    type(struct_gsv) :: statevector_out

    ! locals
    integer :: varIndex, levIndex, nlev, stepIndex, ierr, kIndex

    real(8), pointer :: field_in_r8_ptr(:,:,:,:), field_out_r8_ptr(:,:,:,:)
    real(8), pointer :: fieldUU_in_r8_ptr(:,:,:,:), fieldUU_out_r8_ptr(:,:,:,:)
    real(8), pointer :: fieldVV_in_r8_ptr(:,:,:,:), fieldVV_out_r8_ptr(:,:,:,:)
    real(8), allocatable :: field2dUU_out_r8(:,:), field2dVV_out_r8(:,:)

    character(len=4) :: varName
    character(len=12):: interpolationDegree

    integer :: ezdefset

    if ( hco_equal(statevector_in%hco,statevector_out%hco) ) then
      write(*,*) 'gsv_hInterpolate: The input and output statevectors are already on same horizontal grids'
      call gsv_copy(statevector_in, statevector_out)
      return
    end if

    if ( .not. vco_equal(statevector_in%vco, statevector_out%vco) ) then
      call utl_abort('gsv_hInterpolate: The input and output statevectors are not on the same vertical levels.')
    end if

    if ( statevector_in%dataKind /= 8 .or. statevector_out%dataKind /= 8 ) then
      call utl_abort('gsv_hInterpolate: Incorrect value for dataKind. Only compatible with dataKind=4')
    end if

    ! set the interpolation degree
    interpolationDegree = statevector_out%hInterpolateDegree

    if ( .not.statevector_in%mpi_local .and. .not.statevector_out%mpi_local ) then

      write(*,*) 'gsv_hInterpolate: before interpolation (no mpi)'

      do stepIndex = 1, statevector_out%numStep
        ! Do horizontal interpolation for mpi global statevectors
        VAR_LOOP: do varIndex = 1, vnl_numvarmax
          varName = vnl_varNameList(varIndex)
          if ( .not. gsv_varExist(statevector_in,varName) ) cycle VAR_LOOP
          if ( trim(varName) == 'VV' ) cycle VAR_LOOP

          nlev = statevector_out%varNumLev(varIndex)

          ! horizontal interpolation
          ierr = ezdefset(statevector_out%hco%EZscintID, statevector_in%hco%EZscintID)

          if ( trim(varName) == 'UU' ) then
            ! interpolate both UV components and keep both UU and VV
            fieldUU_in_r8_ptr => gsv_getField_r8(statevector_in,'UU')
            fieldUU_out_r8_ptr => gsv_getField_r8(statevector_out,'UU')
            fieldVV_in_r8_ptr => gsv_getField_r8(statevector_in,'VV')
            fieldVV_out_r8_ptr => gsv_getField_r8(statevector_out,'VV')
            do levIndex = 1, nlev
              ierr = utl_ezuvint( fieldUU_out_r8_ptr(:,:,levIndex,stepIndex), fieldVV_out_r8_ptr(:,:,levIndex,stepIndex), &
                                  fieldUU_in_r8_ptr(:,:,levIndex,stepIndex),  fieldVV_in_r8_ptr(:,:,levIndex,stepIndex),  & 
                                  interpDegree_opt=trim(interpolationDegree) ) 
            end do
          else
            ! interpolate scalar variable
            field_in_r8_ptr => gsv_getField_r8(statevector_in, varName)
            field_out_r8_ptr => gsv_getField_r8(statevector_out, varName)
            do levIndex = 1, nlev
              ierr = utl_ezsint( field_out_r8_ptr(:,:,levIndex,stepIndex), field_in_r8_ptr(:,:,levIndex,stepIndex),  &
                                 interpDegree_opt=trim(interpolationDegree) )
            end do
          end if
        end do VAR_LOOP

      end do ! stepIndex

    else

      write(*,*) 'gsv_hInterpolate: before interpolation (with mpi)'

      if ( statevector_in%mpi_distribution /= 'VarsLevs' .or.   &
          statevector_out%mpi_distribution /= 'VarsLevs' ) then
        call utl_abort('gsv_hInterpolate: The input or output statevector is not distributed by VarsLevs.')
      end if

      do stepIndex = 1, statevector_out%numStep
        K_LOOP: do kIndex = statevector_in%mykBeg, statevector_in%mykEnd
          varName = gsv_getVarNameFromK(statevector_in,kIndex)
          if ( .not. gsv_varExist(statevector_in,varName) ) cycle K_LOOP

          ! horizontal interpolation
          ierr = ezdefset(statevector_out%hco%EZscintID, statevector_in%hco%EZscintID)

          if ( trim(varName) == 'UU' ) then
            ! interpolate both UV components and keep UU
            fieldUU_in_r8_ptr => gsv_getField_r8(statevector_in)
            fieldUU_out_r8_ptr => gsv_getField_r8(statevector_out)
            fieldVV_in_r8_ptr => gsv_getFieldUV_r8(statevector_in)
            allocate( field2dVV_out_r8(statevector_out%ni, statevector_out%nj) )
            ierr = utl_ezuvint( fieldUU_out_r8_ptr(:,:,kIndex,stepIndex), field2dVV_out_r8(:,:),   &
                                fieldUU_in_r8_ptr(:,:,kIndex,stepIndex),  fieldVV_in_r8_ptr(:,:,kIndex,stepIndex), &
                                interpDegree_opt=trim(interpolationDegree) ) 
            deallocate( field2dVV_out_r8 )
          elseif ( trim(varName) == 'VV' ) then
            ! interpolate both UV components and keep VV
            fieldUU_in_r8_ptr => gsv_getFieldUV_r8(statevector_in)
            allocate( field2dUU_out_r8(statevector_out%ni, statevector_out%nj) )
            fieldVV_in_r8_ptr => gsv_getField_r8(statevector_in)
            fieldVV_out_r8_ptr => gsv_getField_r8(statevector_out)
            ierr = utl_ezuvint( field2dUU_out_r8(:,:), fieldVV_out_r8_ptr(:,:,kIndex,stepIndex),   &
                                fieldUU_in_r8_ptr(:,:,kIndex,stepIndex), fieldVV_in_r8_ptr(:,:,kIndex,stepIndex), &
                                interpDegree_opt=trim(interpolationDegree) ) 
            deallocate( field2dUU_out_r8 )
          else
            ! interpolate scalar variable
            field_in_r8_ptr => gsv_getField_r8(statevector_in)
            field_out_r8_ptr => gsv_getField_r8(statevector_out)
            ierr = utl_ezsint( field_out_r8_ptr(:,:,kIndex,stepIndex), field_in_r8_ptr(:,:,kIndex,stepIndex), &
                               interpDegree_opt=trim(interpolationDegree) )
          end if
        end do K_LOOP

      end do ! stepIndex

    end if

    if ( associated(statevector_in%gzSfc) .and. associated(statevector_out%gzSfc) ) then
      write(*,*) 'gsv_hInterpolate: interpolating surface GZ'
      ierr = ezdefset(statevector_out%hco%EZscintID, statevector_in%hco%EZscintID)
      ierr = utl_ezsint( statevector_out%GZsfc(:,:), statevector_in%GZsfc(:,:), &
                         interpDegree_opt=trim(interpolationDegree) )
    end if

  end subroutine gsv_hInterpolate

  !--------------------------------------------------------------------------
  ! gsv_hInterpolate_r4
  !--------------------------------------------------------------------------
  subroutine gsv_hInterpolate_r4(statevector_in,statevector_out)
    ! s/r gsv_hInterpolate_r4  - Horizontal interpolation of pressure defined fields
    implicit none

    ! arguments
    type(struct_gsv) :: statevector_in
    type(struct_gsv) :: statevector_out

    ! locals
    integer :: varIndex, levIndex, nlev, stepIndex, ierr, kIndex

    real(4), pointer :: field_in_r4_ptr(:,:,:,:), field_out_r4_ptr(:,:,:,:)
    real(4), pointer :: fieldUU_in_r4_ptr(:,:,:,:), fieldUU_out_r4_ptr(:,:,:,:)
    real(4), pointer :: fieldVV_in_r4_ptr(:,:,:,:), fieldVV_out_r4_ptr(:,:,:,:)
    real(4), allocatable :: field2dUU_out_r4(:,:), field2dVV_out_r4(:,:)

    character(len=4) :: varName
    character(len=12):: interpolationDegree
    integer :: ezdefset

    if ( hco_equal(statevector_in%hco,statevector_out%hco) ) then
      write(*,*) 'gsv_hInterpolate_r4: The input and output statevectors are already on same horizontal grids'
      call gsv_copy(statevector_in, statevector_out)
      return
    end if

    if ( .not. vco_equal(statevector_in%vco, statevector_out%vco) ) then
      call utl_abort('gsv_hInterpolate_r4: The input and output statevectors are not on the same vertical levels.')
    end if

    if ( statevector_in%dataKind /= 4 .or. statevector_out%dataKind /= 4 ) then
      call utl_abort('gsv_hInterpolate_r4: Incorrect value for dataKind. Only compatible with dataKind=4')
    end if

    ! set the interpolation degree
    interpolationDegree = statevector_out%hInterpolateDegree

    if ( .not.statevector_in%mpi_local .and. .not.statevector_out%mpi_local ) then

      write(*,*) 'gsv_hInterpolate_r4: before interpolation (no mpi)'

      do stepIndex = 1, statevector_out%numStep
        ! Do horizontal interpolation for mpi global statevectors
        VAR_LOOP: do varIndex = 1, vnl_numvarmax
          varName = vnl_varNameList(varIndex)
          if ( .not. gsv_varExist(statevector_in,varName) ) cycle VAR_LOOP
          if ( trim(varName) == 'VV' ) cycle VAR_LOOP

          nlev = statevector_out%varNumLev(varIndex)

          ! horizontal interpolation
          ierr = ezdefset(statevector_out%hco%EZscintID, statevector_in%hco%EZscintID)

          if ( trim(varName) == 'UU' ) then
            ! interpolate both UV components and keep both UU and VV
            fieldUU_in_r4_ptr => gsv_getField_r4(statevector_in,'UU')
            fieldUU_out_r4_ptr => gsv_getField_r4(statevector_out,'UU')
            fieldVV_in_r4_ptr => gsv_getField_r4(statevector_in,'VV')
            fieldVV_out_r4_ptr => gsv_getField_r4(statevector_out,'VV')
            do levIndex = 1, nlev
              ierr = utl_ezuvint( fieldUU_out_r4_ptr(:,:,levIndex,stepIndex), fieldVV_out_r4_ptr(:,:,levIndex,stepIndex),   &
                                  fieldUU_in_r4_ptr(:,:,levIndex,stepIndex),  fieldVV_in_r4_ptr(:,:,levIndex,stepIndex),    &
                                  interpDegree_opt=trim(InterpolationDegree) ) 
            end do
          else
            ! interpolate scalar variable
            field_in_r4_ptr => gsv_getField_r4(statevector_in, varName)
            field_out_r4_ptr => gsv_getField_r4(statevector_out, varName)
            do levIndex = 1, nlev
              ierr = utl_ezsint( field_out_r4_ptr(:,:,levIndex,stepIndex), field_in_r4_ptr(:,:,levIndex,stepIndex),  &
                                 interpDegree_opt=trim(InterpolationDegree) )
            end do
          end if
        end do VAR_LOOP

      end do ! stepIndex

    else

      write(*,*) 'gsv_hInterpolate_r4: before interpolation (with mpi)'

      if ( statevector_in%mpi_distribution /= 'VarsLevs' .or.   &
          statevector_out%mpi_distribution /= 'VarsLevs' ) then
        call utl_abort('gsv_hInterpolate_r4: The input or output statevector is not distributed by VarsLevs.')
      end if

      do stepIndex = 1, statevector_out%numStep
        K_LOOP: do kIndex = statevector_in%mykBeg, statevector_in%mykEnd
          varName = gsv_getVarNameFromK(statevector_in,kIndex)
          if ( .not. gsv_varExist(statevector_in,varName) ) cycle K_LOOP

          ! horizontal interpolation
          ierr = ezdefset(statevector_out%hco%EZscintID, statevector_in%hco%EZscintID)

          if ( trim(varName) == 'UU' ) then
            ! interpolate both UV components and keep UU
            fieldUU_in_r4_ptr => gsv_getField_r4(statevector_in)
            fieldUU_out_r4_ptr => gsv_getField_r4(statevector_out)
            fieldVV_in_r4_ptr => gsv_getFieldUV_r4(statevector_in)
            allocate( field2dVV_out_r4(statevector_out%ni, statevector_out%nj) )
            ierr = utl_ezuvint( fieldUU_out_r4_ptr(:,:,kIndex,stepIndex), field2dVV_out_r4(:,:),   &
                                fieldUU_in_r4_ptr(:,:,kIndex,stepIndex),  fieldVV_in_r4_ptr(:,:,kIndex,stepIndex), &
                                interpDegree_opt=trim(InterpolationDegree) ) 
            deallocate( field2dVV_out_r4 )
          elseif ( trim(varName) == 'VV' ) then
            ! interpolate both UV components and keep VV
            fieldUU_in_r4_ptr => gsv_getFieldUV_r4(statevector_in)
            allocate( field2dUU_out_r4(statevector_out%ni, statevector_out%nj) )
            fieldVV_in_r4_ptr => gsv_getField_r4(statevector_in)
            fieldVV_out_r4_ptr => gsv_getField_r4(statevector_out)
            ierr = utl_ezuvint( field2dUU_out_r4(:,:), fieldVV_out_r4_ptr(:,:,kIndex,stepIndex),   &
                                fieldUU_in_r4_ptr(:,:,kIndex,stepIndex),  fieldVV_in_r4_ptr(:,:,kIndex,stepIndex),  &
                                interpDegree_opt=trim(InterpolationDegree) ) 
            deallocate( field2dUU_out_r4 )
          else
            ! interpolate scalar variable
            field_in_r4_ptr => gsv_getField_r4(statevector_in)
            field_out_r4_ptr => gsv_getField_r4(statevector_out)
            ierr = utl_ezsint( field_out_r4_ptr(:,:,kIndex,stepIndex), field_in_r4_ptr(:,:,kIndex,stepIndex),  &
                               interpDegree_opt=trim(InterpolationDegree) )
          end if
        end do K_LOOP

      end do ! stepIndex

    end if

    if ( associated(statevector_in%gzSfc) .and. associated(statevector_out%gzSfc) ) then
      write(*,*) 'gsv_hInterpolate_r4: interpolating surface GZ'
      ierr = ezdefset(statevector_out%hco%EZscintID, statevector_in%hco%EZscintID)
      ierr = utl_ezsint( statevector_out%GZsfc(:,:), statevector_in%GZsfc(:,:),  &
                         interpDegree_opt=trim(InterpolationDegree) )
    end if

  end subroutine gsv_hInterpolate_r4

  !--------------------------------------------------------------------------
  ! gsv_vInterpolate
  !--------------------------------------------------------------------------
  subroutine gsv_vInterpolate(statevector_in,statevector_out,Ps_in_hPa_opt,PsfcReference_opt)
    ! s/r gsv_vInterpolate  - Vertical interpolation of pressure defined fields
    implicit none

    ! arguments
    type(struct_gsv)  :: statevector_in
    type(struct_gsv)  :: statevector_out
    logical, optional :: Ps_in_hPa_opt
    real(8), optional :: PsfcReference_opt(:,:,:)

    ! locals
    character(len=4) :: varName
    type(struct_vco), pointer :: vco_in, vco_out
    integer :: nlev_out, nlev_in, levIndex_out, levIndex_in, latIndex, lonIndex
    integer :: status, latIndex2, lonIndex2, varIndex, stepIndex
    real(8) :: zwb, zwt
    real(8), pointer  :: pres_out(:,:,:), pres_in(:,:,:), field_out(:,:,:,:), field_in(:,:,:,:)
    real(8) :: psfc_in(statevector_in%myLonBeg:statevector_in%myLonEnd, &
                       statevector_in%myLatBeg:statevector_in%myLatEnd)

    if ( vco_equal(statevector_in%vco, statevector_out%vco) ) then
      write(*,*) 'gsv_vInterpolate: The input and output statevectors are already on same vertical levels'
      call gsv_copy(statevector_in, statevector_out)
      return
    end if

    if ( .not. hco_equal(statevector_in%hco, statevector_out%hco) ) then
      call utl_abort('gsv_vInterpolate: The input and output statevectors are not on the same horizontal grid.')
    end if

    if ( statevector_in%dataKind /= 8 .or. statevector_out%dataKind /= 8 ) then
      call utl_abort('gsv_vInterpolate: Incorrect value for dataKind. Only compatible with dataKind=8')
    end if

    if ( associated(statevector_in%gzSfc) .and. associated(statevector_out%gzSfc) ) then
      statevector_out%gzSfc(:,:) = statevector_in%gzSfc(:,:)
    end if

    vco_in => gsv_getVco(statevector_in)
    vco_out => gsv_getVco(statevector_out)

    do stepIndex = 1, statevector_out%numStep

      if ( present(PsfcReference_opt) ) then
        psfc_in(:,:) = PsfcReference_opt(:,:,stepIndex)
      else
        field_in => gsv_getField_r8(statevector_in,'P0')
        psfc_in(:,:) = field_in(:,:,1,stepIndex)
      end if
      if ( present(Ps_in_hPa_opt) ) then
        if ( Ps_in_hPa_opt ) psfc_in = psfc_in * MPC_PA_PER_MBAR_R8
      end if

      VAR_LOOP: do varIndex = 1, vnl_numvarmax
        varName = vnl_varNameList(varIndex)
        if ( .not. gsv_varExist(statevector_in,varName) ) cycle VAR_LOOP

        nlev_in  = statevector_in%varNumLev(varIndex)
        nlev_out = statevector_out%varNumLev(varIndex)

        field_in  => gsv_getField_r8(statevector_in ,varName)
        field_out => gsv_getField_r8(statevector_out,varName)

        ! for 2D fields, just copy and cycle to next variable
        if ( nlev_in == 1 .and. nlev_out == 1 ) then
          field_out(:,:,:,stepIndex) = field_in(:,:,:,stepIndex)
          cycle VAR_LOOP
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

        deallocate(pres_out)
        deallocate(pres_in)

      end do VAR_LOOP

    end do ! stepIndex

  end subroutine gsv_vInterpolate

  !--------------------------------------------------------------------------
  ! gsv_vInterpolate_r4
  !--------------------------------------------------------------------------
  subroutine gsv_vInterpolate_r4(statevector_in,statevector_out,Ps_in_hPa_opt,PsfcReference_opt)
    ! s/r gsv_vInterpolate_r4  - Vertical interpolation of pressure defined fields
    implicit none

    ! arguments
    type(struct_gsv)  :: statevector_in
    type(struct_gsv)  :: statevector_out
    logical, optional :: Ps_in_hPa_opt
    real(4), optional :: PsfcReference_opt(:,:,:)

    ! locals
    character(len=4) :: varName
    type(struct_vco), pointer :: vco_in, vco_out
    integer :: nlev_out, nlev_in, levIndex_out, levIndex_in, latIndex, lonIndex
    integer :: status, latIndex2, lonIndex2, varIndex, stepIndex
    real(4) :: zwb, zwt
    real(4), pointer  :: pres_out(:,:,:), pres_in(:,:,:), field_out(:,:,:,:), field_in(:,:,:,:)
    real(4) :: psfc_in(statevector_in%myLonBeg:statevector_in%myLonEnd, &
                       statevector_in%myLatBeg:statevector_in%myLatEnd)

    if ( vco_equal(statevector_in%vco, statevector_out%vco) ) then
      write(*,*) 'gsv_vInterpolate_r4: The input and output statevectors are already on same vertical levels'
      call gsv_copy(statevector_in, statevector_out)
      return
    end if

    if ( .not. hco_equal(statevector_in%hco, statevector_out%hco) ) then
      call utl_abort('gsv_vInterpolate_r4: The input and output statevectors are not on the same horizontal grid.')
    end if

    if ( statevector_in%dataKind /= 4 .or. statevector_out%dataKind /= 4 ) then
      call utl_abort('gsv_vInterpolate_r4: Incorrect value for dataKind. Only compatible with dataKind=4')
    end if

    if ( associated(statevector_in%gzSfc) .and. associated(statevector_out%gzSfc) ) then
      statevector_out%gzSfc(:,:) = statevector_in%gzSfc(:,:)
    end if

    vco_in => gsv_getVco(statevector_in)
    vco_out => gsv_getVco(statevector_out)

    do stepIndex = 1, statevector_out%numStep

      if ( present(PsfcReference_opt) ) then
        psfc_in(:,:) = PsfcReference_opt(:,:,stepIndex)
      else
        field_in => gsv_getField_r4(statevector_in,'P0')
        psfc_in(:,:) = field_in(:,:,1,stepIndex)
      end if
      if ( present(Ps_in_hPa_opt) ) then
        if ( Ps_in_hPa_opt ) psfc_in = psfc_in * MPC_PA_PER_MBAR_R4
      end if

      VAR_LOOP: do varIndex = 1, vnl_numvarmax
        varName = vnl_varNameList(varIndex)
        if ( .not. gsv_varExist(statevector_in,varName) ) cycle VAR_LOOP

        nlev_in  = statevector_in%varNumLev(varIndex)
        nlev_out = statevector_out%varNumLev(varIndex)

        field_in  => gsv_getField_r4(statevector_in ,varName)
        field_out => gsv_getField_r4(statevector_out,varName)

        ! for 2D fields, just copy and cycle to next variable
        if ( nlev_in == 1 .and. nlev_out == 1 ) then
          field_out(:,:,:,stepIndex) = field_in(:,:,:,stepIndex)
          cycle VAR_LOOP
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

        deallocate(pres_out)
        deallocate(pres_in)

      end do VAR_LOOP

    end do ! stepIndex

  end subroutine gsv_vInterpolate_r4

  !--------------------------------------------------------------------------
  ! gsv_writeToFileMpi
  !--------------------------------------------------------------------------
  subroutine gsv_writeToFileMpi(statevector, fileName, etiket_in, scaleFactor_opt, ip3_opt, &
       stepIndex_opt, typvar_opt, HUcontainsLQ_opt, unitConversion_opt)
    implicit none

    ! arguments
    type(struct_gsv)             :: statevector
    character(len=*), intent(in) :: fileName
    character(len=*), intent(in) :: etiket_in
    real(8), optional,intent(in) :: scaleFactor_opt
    integer, optional,intent(in) :: ip3_opt, stepIndex_opt
    character(len=*), optional, intent(in) :: typvar_opt
    logical, optional,intent(in) :: HUcontainsLQ_opt
    logical, optional,intent(in) :: unitConversion_opt

    ! locals
    integer :: fclos, fnom, fstouv, fstfrm
    integer :: nulfile, stepIndex
    real(4), allocatable :: work2d_r4(:,:), gd_send_r4(:,:,:), gd_recv_r4(:,:,:)
    real(4)   :: factor_r4, work_r4
    integer   :: ierr, fstecr
    integer   :: var, k, kgdim
    integer :: ni, nj, nk
    integer :: dateo, npak, ji, jj, levIndex, levIndex_last, nlev, varIndex
    integer :: ip1, ip2, ip3, deet, npas, datyp
    integer :: ig1 ,ig2 ,ig3 ,ig4
    integer :: batchnum, levIndex2, yourid, nsize, youridy, youridx
    integer :: writeLevPE(500)
    character(len=1)  :: grtyp
    character(len=4)  :: nomvar
    character(len=2)  :: typvar
    character(len=12) :: etiket
    character(len=3)  :: charprocid
    character(len=128) :: fileName_full
    logical :: iDoWriting, unitConversion
    real(8), pointer :: field_r8(:,:,:,:)
    real(4), pointer :: field_r4(:,:,:,:)
    real(8), allocatable :: work3d(:,:,:)

    write(*,*) 'gsv_writeToFileMpi: START'

    if ( .not. statevector%mpi_local ) then
      call utl_abort('gsv_writeToFileMpi: cannot use this subroutine for mpiglobal data!')
    end if

    call tmg_start(5,'WRITETOFILE')

    if (present(ip3_opt)) then
      ip3 = ip3_opt
    else
      ip3 = 0
    end if

    if (present(unitConversion_opt)) then
      unitConversion = unitConversion_opt
    else
      unitConversion = .true.
    end if

    ! if step index not specified, choose anltime (usually center of window)
    if (present(stepIndex_opt)) then
      stepIndex = stepIndex_opt
    else
      stepIndex = statevector%anltime
    end if

    ! initialization of parameters for writing to file
    npak   = -32
    dateo  = statevector%dateStampList(stepIndex)
    deet   = 0
    npas   = 0
    ni     = statevector%ni
    nj     = statevector%nj
    nk     = 1
    ip2    = 0
    if ( present(typvar_opt) ) then
      typvar = trim(typvar_opt)
    else
      typvar = 'R'
    end if
    etiket = trim(Etiket_in)
    grtyp  = statevector%hco%grtyp
    ig1    = statevector%hco%ig1
    ig2    = statevector%hco%ig2
    ig3    = statevector%hco%ig3
    ig4    = statevector%hco%ig4
    datyp  = 1

    nlev = max(gsv_getNumLev(statevector,'MM'),gsv_getNumLev(statevector,'TH'))
    ! for each level determine which processor should do the reading
    do levIndex = 1, nlev
      writeLevPE(levIndex) = mod(levIndex-1,mpi_nprocs)
    end do
    iDoWriting = (mpi_myid+1).le.nlev

    !
    !- Write the global StateVector
    !
    if (iDoWriting) then

      !- Open output field
      nulfile = 0
      write(charProcId,'(i3.3)') mpi_myid
      fileName_full = trim(fileName)//'_proc'//charProcId
      write(*,*) 'gsv_writeToFileMpi: file name = ',trim(fileName_full)
      ierr = fnom(nulfile,trim(fileName_full),'RND+APPEND',0)
       
      if ( ierr >= 0 ) then
        write(*,*)'gsv_writeToFileMpi: file opened with unit number ',nulfile
        ierr  =  fstouv(nulfile,'RND')
      else
        call utl_abort('gsv_writeToFileMpi: problem opening output file, aborting!')
      end if

      if (nulfile == 0 ) then
        call utl_abort('gsv_writeToFileMpi: unit number for output file not valid!')
      end if

      !- Write TicTacToc
      if (mpi_myid == 0) call WriteTicTacToc(statevector,nulfile) ! IN

    end if

    allocate(work2d_r4(statevector%ni,statevector%nj))
    allocate(gd_send_r4(statevector%lonPerPEmax,statevector%latPerPEmax,mpi_nprocs))
    allocate(gd_recv_r4(statevector%lonPerPEmax,statevector%latPerPEmax,mpi_nprocs))

    do varIndex = 1, vnl_numvarmax 
 
      if (gsv_varExist(statevector,vnl_varNameList(varIndex)) ) then

        nlev = statevector%varNumLev(varIndex)

        allocate(work3d(statevector%myLonBeg:statevector%myLonEnd, &
                        statevector%myLatBeg:statevector%myLatEnd,nlev))
        if ( statevector%dataKind == 8 ) then
          field_r8 => gsv_getField_r8(statevector,vnl_varNameList(varIndex))
          work3d(:,:,:) = field_r8(statevector%myLonBeg:statevector%myLonEnd, &
                                   statevector%myLatBeg:statevector%myLatEnd,1:nlev,stepIndex)
        else
          field_r4 => gsv_getField_r4(statevector,vnl_varNameList(varIndex))
          work3d(:,:,:) = real(field_r4(statevector%myLonBeg:statevector%myLonEnd, &
                                        statevector%myLatBeg:statevector%myLatEnd,1:nlev,stepIndex), 8)
        end if

        do batchnum = 1, ceiling(dble(nlev)/dble(mpi_nprocs))

          levIndex      = (batchnum-1)*mpi_nprocs + mpi_myid + 1
          levIndex_last = min(nlev,batchnum*mpi_nprocs)

!$OMP PARALLEL DO PRIVATE(levIndex2,yourid)
          do levIndex2 = 1+(batchnum-1)*mpi_nprocs, levIndex_last
            yourid = writeLevPE(levIndex2)
            gd_send_r4(1:statevector%lonPerPE,  &
                       1:statevector%latPerPE, yourid+1) =  &
                real(work3d(statevector%myLonBeg:statevector%myLonEnd, &
                            statevector%myLatBeg:statevector%myLatEnd, levIndex2),4)
          end do
!$OMP END PARALLEL DO

          nsize = statevector%lonPerPEmax * statevector%latPerPEmax
          if(mpi_nprocs.gt.1) then
            call rpn_comm_alltoall(gd_send_r4,nsize,"mpi_real4",  &
                                   gd_recv_r4,nsize,"mpi_real4","GRID",ierr)
          else
            gd_recv_r4(:,:,1) = gd_send_r4(:,:,1)
          end if

!$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
          do youridy = 0, (mpi_npey-1)
            do youridx = 0, (mpi_npex-1)
              yourid = youridx + youridy*mpi_npex
                work2d_r4(statevector%allLonBeg(youridx+1):statevector%allLonEnd(youridx+1),  &
                          statevector%allLatBeg(youridy+1):statevector%allLatEnd(youridy+1)) = &
                  gd_recv_r4(1:statevector%allLonPerPE(youridx+1),  &
                             1:statevector%allLatPerPE(youridy+1), yourid+1)
            end do
          end do
!$OMP END PARALLEL DO

          ! now do writing
          if (iDoWriting.and.levIndex.le.nlev) then

            ! Set the ip1 value
            if (vnl_varLevelFromVarname(vnl_varNameList(varIndex)) == 'MM') then
               ip1 = statevector%vco%ip1_M(levIndex)
            elseif (vnl_varLevelFromVarname(vnl_varNameList(varIndex)) == 'TH') then
               ip1 = statevector%vco%ip1_T(levIndex)
            elseif (vnl_varLevelFromVarname(vnl_varNameList(varIndex)) == 'SF') then
               ip1 = 0
            else
               write(*,*) 'gsv_writeToFileMpi: unknown type of vertical level: ',  &
                          vnl_varLevelFromVarname(vnl_varNameList(varIndex))
               call utl_abort('gsv_writeToFileMpi')
            end if

            ! Set the output variable name
            nomvar = trim(vnl_varNameList(varIndex))
            if ( trim(nomvar) == 'HU' .and. present(HUcontainsLQ_opt) ) then
               if ( HUcontainsLQ_opt ) nomvar = 'LQ'
            end if

            ! Set the conversion factor
            if ( unitConversion ) then
              if ( trim(nomvar) == 'UU' .or. trim(nomvar) == 'VV') then
                factor_r4 = MPC_KNOTS_PER_M_PER_S_R4 ! m/s -> knots
              else if ( trim(nomvar) == 'P0' ) then
                factor_r4 = 0.01 ! Pa -> hPa
              else
                factor_r4 = 1.0
              end if
            else
              factor_r4 = 1.0
            end if

            if (present(scaleFactor_opt)) factor_r4 = factor_r4 * real(scaleFactor_opt,4)

            !- Scale
            work2d_r4(:,:) = factor_r4 * work2d_r4(:,:)

            !- Writing to file
            ierr = fstecr(work2d_r4, work_r4, npak, nulfile, dateo, deet, npas, ni, nj, &
                          nk, ip1, ip2, ip3, typvar, nomvar, etiket, grtyp,      &
                          ig1, ig2, ig3, ig4, datyp, .false.)

          end if ! iDoWriting

        end do ! batchnum

        deallocate(work3d)

      end if ! varExist

    end do ! varIndex

    deallocate(work2d_r4)
    deallocate(gd_send_r4)
    deallocate(gd_recv_r4)

    if (iDoWriting) then
      ierr = fstfrm(nulfile)
      ierr = fclos(nulfile)        
    end if
       
    call tmg_stop(5)
    write(*,*) 'gsv_writeToFileMpi: END'

  end subroutine gsv_writeToFileMpi

  !--------------------------------------------------------------------------
  ! gsv_writeToFile
  !--------------------------------------------------------------------------
  subroutine gsv_writeToFile(statevector, fileName, etiket_in, scaleFactor_opt, ip3_opt, &
       stepIndex_opt, typvar_opt, HUcontainsLQ_opt, unitConversion_opt, writeGZsfc_opt, numBits_opt)
    implicit none

    ! arguments
    type(struct_gsv)             :: statevector
    character(len=*), intent(in) :: fileName
    character(len=*), intent(in) :: etiket_in
    real(8), optional,intent(in) :: scaleFactor_opt
    integer, optional,intent(in) :: ip3_opt, stepIndex_opt
    character(len=*), optional, intent(in) :: typvar_opt
    logical, optional,intent(in) :: HUcontainsLQ_opt
    logical, optional,intent(in) :: unitConversion_opt
    logical, optional,intent(in) :: writeGZsfc_opt
    integer, optional,intent(in) :: numBits_opt

    ! locals
    integer :: fclos, fnom, fstouv, fstfrm
    integer :: nulfile, stepIndex
    real(4), allocatable :: work2d_r4(:,:), gd_send_r4(:,:), gd_recv_r4(:,:,:)
    real(4)   :: factor_r4, work_r4
    integer   :: ierr, fstecr
    integer   :: var, k, kgdim
    integer :: ni, nj, nk
    integer :: dateo, npak, ji, jj, levIndex, levIndex_last, nlev, varIndex
    integer :: ip1, ip2, ip3, deet, npas, datyp
    integer :: ig1 ,ig2 ,ig3 ,ig4
    integer :: batchnum, levIndex2, yourid, nsize, youridy, youridx
    character(len=1)  :: grtyp
    character(len=4)  :: nomvar
    character(len=2)  :: typvar
    character(len=12) :: etiket
    logical :: iDoWriting, unitConversion
    real(8), pointer :: field_r8(:,:,:,:)
    real(4), pointer :: field_r4(:,:,:,:)

    write(*,*) 'gsv_writeToFile: START'

    if ( .not. statevector%mpi_local ) then
      write(*,*) 'gsv_writeToFile: writing statevector that is already mpiglobal!'
    end if

    call tmg_start(5,'WRITETOFILE')

    if (present(ip3_opt)) then
      ip3 = ip3_opt
    else
      ip3 = 0
    end if

    if (present(unitConversion_opt)) then
      unitConversion = unitConversion_opt
    else
      unitConversion = .true.
    end if

    ! if step index not specified, choose anltime (usually center of window)
    if (present(stepIndex_opt)) then
      stepIndex = stepIndex_opt
    else
      stepIndex = statevector%anltime
    end if

    if ( present(numBits_opt) ) then
      npak = -numBits_opt
    else
      npak = -32
    end if

    ! initialization of parameters for writing to file
    dateo  = statevector%dateStampList(stepIndex)
    deet   = 0
    npas   = 0
    ni     = statevector%ni
    nj     = statevector%nj
    nk     = 1
    ip2    = 0
    if ( present(typvar_opt) ) then
      typvar = trim(typvar_opt)
    else
      typvar = 'R'
    end if
    etiket = trim(Etiket_in)
    grtyp  = statevector%hco%grtyp
    ig1    = statevector%hco%ig1
    ig2    = statevector%hco%ig2
    ig3    = statevector%hco%ig3
    ig4    = statevector%hco%ig4
    datyp  = 134

    nlev = max(gsv_getNumLev(statevector,'MM'),gsv_getNumLev(statevector,'TH'))

    ! only proc 0 does writing or each proc when data is global 
    ! (assuming only called for proc with global data)
    iDoWriting = (mpi_myid == 0) .or. (.not. statevector%mpi_local)

    !
    !- Write the global StateVector
    !
    if (iDoWriting) then

      !- Open output field
      nulfile = 0
      write(*,*) 'gsv_writeToFile: file name = ',trim(fileName)
      ierr = fnom(nulfile,trim(fileName),'RND+APPEND',0)
       
      if ( ierr >= 0 ) then
        write(*,*)'gsv_writeToFile: file opened with unit number ',nulfile
        ierr  =  fstouv(nulfile,'RND')
      else
        call utl_abort('gsv_writeToFile: problem opening output file, aborting!')
      end if

      if (nulfile == 0 ) then
        call utl_abort('gsv_writeToFile: unit number for output file not valid!')
      end if

      !- Write TicTacToc
      if (mpi_myid == 0) call WriteTicTacToc(statevector,nulfile) ! IN

    end if

    allocate(gd_send_r4(statevector%lonPerPEmax,statevector%latPerPEmax))
    if( mpi_myid == 0 .or. (.not. statevector%mpi_local) ) then
      allocate(gd_recv_r4(statevector%lonPerPEmax,statevector%latPerPEmax,mpi_nprocs))
      allocate(work2d_r4(statevector%ni,statevector%nj))
    else
      allocate(gd_recv_r4(1,1,1))
      allocate(work2d_r4(1,1))
    end if

    ! Write surface GZ, if requested
    if ( present(writeGZsfc_opt) ) then
      if ( writeGZsfc_opt .and. associated(statevector%GZsfc) ) then
        write(*,*) 'gsv_writeToFile: writing surface GZ'

        ! MPI communication
        gd_send_r4(1:statevector%lonPerPE, &
                   1:statevector%latPerPE) =  &
             real(statevector%GZsfc(statevector%myLonBeg:statevector%myLonEnd, &
                                    statevector%myLatBeg:statevector%myLatEnd),4)
        if( (mpi_nprocs > 1) .and. (statevector%mpi_local) ) then
          nsize = statevector%lonPerPEmax * statevector%latPerPEmax
          call rpn_comm_gather(gd_send_r4, nsize, "mpi_real4",  &
                               gd_recv_r4, nsize, "mpi_real4", 0, "GRID", ierr )
        else
          ! just copy when either nprocs is 1 or data is global
          gd_recv_r4(:,:,1) = gd_send_r4(:,:)
        end if
        if ( mpi_myid == 0 .and. statevector%mpi_local ) then
!$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
          do youridy = 0, (mpi_npey-1)
            do youridx = 0, (mpi_npex-1)
              yourid = youridx + youridy*mpi_npex
                work2d_r4(statevector%allLonBeg(youridx+1):statevector%allLonEnd(youridx+1),  &
                          statevector%allLatBeg(youridy+1):statevector%allLatEnd(youridy+1)) = &
                  gd_recv_r4(1:statevector%allLonPerPE(youridx+1),  &
                             1:statevector%allLatPerPE(youridy+1),yourid+1)
            end do
          end do
!$OMP END PARALLEL DO
        elseif ( .not. statevector%mpi_local ) then
          work2d_r4(:,:) = gd_recv_r4(:,:,1)
        end if

        ! now do writing
        if (iDoWriting) then
          ip1 = statevector%vco%ip1_sfc
          nomvar = 'GZ'

          !- Scale
          factor_r4 = real(1.0d0/(10.0d0 * RG),4)
          work2d_r4(:,:) = factor_r4 * work2d_r4(:,:)

          !- Writing to file
          ierr = fstecr(work2d_r4, work_r4, npak, nulfile, dateo, deet, npas, ni, nj, &
                        nk, ip1, ip2, ip3, typvar, nomvar, etiket, grtyp,      &
                        ig1, ig2, ig3, ig4, datyp, .false.)
        end if ! iDoWriting

      end if
    end if

    do varIndex = 1, vnl_numvarmax 
 
      if (gsv_varExist(statevector,vnl_varNameList(varIndex)) ) then

        nlev = statevector%varNumLev(varIndex)

        do levIndex = 1, nlev

          if ( statevector%dataKind == 8 ) then
            field_r8 => gsv_getField_r8(statevector,vnl_varNameList(varIndex))
            gd_send_r4(1:statevector%lonPerPE,  &
                       1:statevector%latPerPE) =  &
                real(field_r8(statevector%myLonBeg:statevector%myLonEnd, &
                              statevector%myLatBeg:statevector%myLatEnd,levIndex,stepIndex),4)
          else
            field_r4 => gsv_getField_r4(statevector,vnl_varNameList(varIndex))
            gd_send_r4(1:statevector%lonPerPE,  &
                       1:statevector%latPerPE) =  &
                field_r4(statevector%myLonBeg:statevector%myLonEnd, &
                         statevector%myLatBeg:statevector%myLatEnd,levIndex,stepIndex)
          end if

          nsize = statevector%lonPerPEmax*statevector%latPerPEmax
          if( (mpi_nprocs > 1) .and. (statevector%mpi_local) ) then
            call rpn_comm_gather(gd_send_r4, nsize, "mpi_real4",  &
                                 gd_recv_r4, nsize, "mpi_real4", 0, "GRID", ierr )
          else
            ! just copy when either nprocs is 1 or data is global
            gd_recv_r4(:,:,1) = gd_send_r4(:,:)
          end if

          if ( mpi_myid == 0 .and. statevector%mpi_local ) then
!$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
            do youridy = 0, (mpi_npey-1)
              do youridx = 0, (mpi_npex-1)
                yourid = youridx + youridy*mpi_npex
                  work2d_r4(statevector%allLonBeg(youridx+1):statevector%allLonEnd(youridx+1),  &
                            statevector%allLatBeg(youridy+1):statevector%allLatEnd(youridy+1)) = &
                    gd_recv_r4(1:statevector%allLonPerPE(youridx+1),  &
                               1:statevector%allLatPerPE(youridy+1), yourid+1)
              end do
            end do
!$OMP END PARALLEL DO
          elseif ( .not. statevector%mpi_local ) then
            work2d_r4(:,:) = gd_recv_r4(:,:,1)
          end if

          ! now do writing
          if (iDoWriting) then

            ! Set the ip1 value
            if (vnl_varLevelFromVarname(vnl_varNameList(varIndex)) == 'MM') then
               ip1 = statevector%vco%ip1_M(levIndex)
            elseif (vnl_varLevelFromVarname(vnl_varNameList(varIndex)) == 'TH') then
               ip1 = statevector%vco%ip1_T(levIndex)
            elseif (vnl_varLevelFromVarname(vnl_varNameList(varIndex)) == 'SF') then
               ip1 = 0
            else
               write(*,*) 'gsv_writeToFile: unknown type of vertical level: ',  &
                          vnl_varLevelFromVarname(vnl_varNameList(varIndex))
               call utl_abort('gsv_writeToFile')
            end if

            ! Set the output variable name
            nomvar = trim(vnl_varNameList(varIndex))
            if ( trim(nomvar) == 'HU' .and. present(HUcontainsLQ_opt) ) then
               if ( HUcontainsLQ_opt ) nomvar = 'LQ'
            end if

            ! Set the conversion factor
            if ( unitConversion ) then
              if ( trim(nomvar) == 'UU' .or. trim(nomvar) == 'VV') then
                factor_r4 = MPC_KNOTS_PER_M_PER_S_R4 ! m/s -> knots
              else if ( trim(nomvar) == 'P0' ) then
                factor_r4 = 0.01 ! Pa -> hPa
              else
                factor_r4 = 1.0
              end if
            else
              factor_r4 = 1.0
            end if

            if (present(scaleFactor_opt)) factor_r4 = factor_r4 * real(scaleFactor_opt,4)

            !- Scale
            work2d_r4(:,:) = factor_r4 * work2d_r4(:,:)

            !- Writing to file
            ierr = fstecr(work2d_r4, work_r4, npak, nulfile, dateo, deet, npas, ni, nj, &
                          nk, ip1, ip2, ip3, typvar, nomvar, etiket, grtyp,      &
                          ig1, ig2, ig3, ig4, datyp, .false.)

          end if ! iDoWriting

        end do ! levIndex

      end if ! varExist

    end do ! varIndex

    deallocate(work2d_r4)
    deallocate(gd_send_r4)
    deallocate(gd_recv_r4)

    if (iDoWriting) then
      ierr = fstfrm(nulfile)
      ierr = fclos(nulfile)        
    end if
       
    call tmg_stop(5)
    write(*,*) 'gsv_writeToFile: END'

  end subroutine gsv_writeToFile

  !--------------------------------------------------------------------------
  ! WriteTicTacToc
  !--------------------------------------------------------------------------
  subroutine WriteTicTacToc(statevector,iun)
    use vGrid_Descriptors , only: vgrid_descriptor, vgd_write, VGD_OK
    use MathPhysConstants_mod, only : MPC_DEGREES_PER_RADIAN_R8
    implicit none

    type(struct_gsv)    :: statevector
    integer, intent(in) :: iun

    integer :: ier

    integer :: dateo, npak, status, fstecr
    integer :: ip1,ip2,ip3,deet,npas,datyp,ig1,ig2,ig3,ig4
    integer :: ig1_tictac,ig2_tictac,ig3_tictac,ig4_tictac

    character(len=1)  :: grtyp
    character(len=2)  :: typvar
    character(len=12) :: etiket

    !
    !- 1.  Writing Tic-Tac
    !
    if ( statevector % hco % grtyp == 'Z' ) then
       npak     = -32
       deet     =  0
       ip1      =  statevector%hco%ig1
       ip2      =  statevector%hco%ig2
       ip3      =  0
       npas     =  0
       datyp    =  1
       grtyp    = 'E'
       typvar   = 'X'
       etiket   = 'TICTICTACTAC'
       dateo =  0

       call cxgaig ( grtyp,                                          & ! IN
                     ig1_tictac, ig2_tictac, ig3_tictac, ig4_tictac, & ! OUT
                     real(statevector%hco%xlat1), real(statevector%hco%xlon1),   & ! IN
                     real(statevector%hco%xlat2), real(statevector%hco%xlon2)  )   ! IN

       ig1      =  ig1_tictac
       ig2      =  ig2_tictac
       ig3      =  ig3_tictac
       ig4      =  ig4_tictac

       ier = utl_fstecr(statevector%hco%lon*MPC_DEGREES_PER_RADIAN_R8, npak, &
                        iun, dateo, deet, npas, statevector%ni, 1, 1, ip1,    &
                        ip2, ip3, typvar, '>>', etiket, grtyp, ig1,          &
                        ig2, ig3, ig4, datyp, .true.)

       ier = utl_fstecr(statevector%hco%lat*MPC_DEGREES_PER_RADIAN_R8, npak, &
                        iun, dateo, deet, npas, 1, statevector%nj, 1, ip1,    &
                        ip2, ip3, typvar, '^^', etiket, grtyp, ig1,          &
                        ig2, ig3, ig4, datyp, .true.)

    else if ( statevector % hco % grtyp == 'U' ) then
      etiket   = 'TICTICTACTAC'
      npak     = -32
      ier = fstecr(statevector%hco%tictacU, statevector%hco%tictacU, npak, iun, 0, 0, 0, size(statevector%hco%tictacU), 1, 1  , &
                   statevector%hco%ig1, statevector%hco%ig2,  statevector%hco%ig3, 'X', '^>', etiket, &
                   'F', 1, 0, 0, 0, 5, .false.)

    end if

    !
    !- Writing Toc-Toc
    !
    if ( trim(statevector%vco%setupType) == 'FromFile' ) then 
       status = vgd_write(statevector%vco%vgrid,iun,'fst')
       if ( status /= VGD_OK ) then
          call utl_abort('WriteTicTacToc: ERROR with vgd_write')
       end if
    end if

  end subroutine WriteTicTacToc

  !--------------------------------------------------------------------------
  ! gsv_horizSubSample
  !--------------------------------------------------------------------------
  subroutine gsv_horizSubSample(statevector_in,statevector_out,horizSubSample)
    implicit none
    type(struct_gsv)    :: statevector_in, statevector_out
    integer, intent(in) :: horizSubSample
    integer             :: relativeFactor, middleStep
    real(8)             :: ratio_r8
    integer             :: stepIndex, lonIndex, latIndex, ilon_in, ilon_out, ilat_in, ilat_out
    integer             :: ilon_in1, ilon_in2, lonIndex_in, ilat_in1, ilat_in2, latIndex_in

    if ( statevector_out%allocated ) then
      call gsv_deallocate(statevector_out)
    end if

    ! allocate the output statevector
    if ( associated(statevector_in%dateStampList) ) then
      middleStep = nint((real(statevector_in%numStep,8)+1.0d0)/2.0d0)
      call gsv_allocate(statevector_out,  &
                        statevector_in%numStep, statevector_in%hco, statevector_in%vco, &
                        datestamp_opt = statevector_in%dateStampList(middleStep),  &
                        mpi_local_opt = statevector_in%mpi_local,  &
                        horizSubSample_opt = horizSubSample)
    else
      call gsv_allocate(statevector_out,  &
                        statevector_in%numStep, statevector_in%hco, statevector_in%vco, &
                        mpi_local_opt = statevector_in%mpi_local, &
                        horizSubSample_opt = horizSubSample)
    end if

    if ( statevector_out%horizSubSample == statevector_in%horizSubSample ) then
      if ( mpi_myid == 0 ) write(*,*) 'gsv_horizSubSample: already at the selected subsample level: ', &
                                     statevector_out%horizSubSample
      call gsv_copy(statevector_in,statevector_out)
      return
    end if

    if ( statevector_out%horizSubSample > statevector_in%horizSubSample ) then

      ! simple averaging onto a coarser grid
      if ( mpi_myid == 0 ) write(*,*) 'gsv_horizSubSample: increasing subsample level from ',  &
                                     statevector_in%horizSubSample, ' to ',  &
                                     statevector_out%horizSubSample

      ratio_r8 = real(statevector_out%horizSubSample,8)/real(statevector_in%horizSubSample,8)
      if ( abs(ratio_r8 - real(nint(ratio_r8),8)) > 1.0d-5 ) then
        write(*,*) 'gsv_horizSubSample: original subsample level=', statevector_in%horizSubSample
        write(*,*) 'gsv_horizSubSample: new      subsample level=', statevector_out%horizSubSample
        call utl_abort('gsv_horizSubSample: relative change of subsample level not an integer!')
      end if
      relativeFactor = nint(ratio_r8)

      call gsv_zero(statevector_out)
!$OMP PARALLEL DO PRIVATE(stepIndex, latIndex, ilat_out, ilat_in1, ilat_in2, latIndex_in, &
!$OMP                     lonIndex, ilon_out, ilon_in1, ilon_in2, lonIndex_in)
      do stepIndex = 1, statevector_out%numStep

        do latIndex = 1, statevector_out%latPerPE
          ilat_out = latIndex + statevector_out%myLatBeg - 1
          ilat_in1 = relativeFactor*(latIndex-1) + statevector_in%myLatBeg
          ilat_in2 = ilat_in1 + relativeFactor - 1
          do latIndex_in = ilat_in1, ilat_in2

            do lonIndex = 1, statevector_out%lonPerPE
              ilon_out = lonIndex + statevector_out%myLonBeg - 1
              ilon_in1 = relativeFactor*(lonIndex-1) + statevector_in%myLonBeg
              ilon_in2 = ilon_in1 + relativeFactor - 1
              do lonIndex_in = ilon_in1, ilon_in2

                statevector_out%gd_r8(ilon_out, ilat_out, :, stepIndex) =  &
                  statevector_out%gd_r8(ilon_out, ilat_out, :, stepIndex) +  &
                  statevector_in%gd_r8(lonIndex_in, ilat_in, :, stepIndex)

              end do ! lonIndex_in
            end do ! lonIndex

          end do ! latIndex_in
        end do ! latIndex

      end do ! stepIndex
!$OMP END PARALLEL DO 

    else

      ! interpolate to a finer grid
      if ( mpi_myid == 0 ) write(*,*) 'gsv_horizSubSample: decreasing subsample level from ',  &
                                     statevector_in%horizSubSample, ' to ',  &
                                     statevector_out%horizSubSample

      ratio_r8 = real(statevector_in%horizSubSample)/real(statevector_out%horizSubSample)
      if ( abs(ratio_r8 - real(nint(ratio_r8),8)) > 1.0d-5 ) then
        write(*,*) 'gsv_horizSubSample: original subsample level=', statevector_in%horizSubSample
        write(*,*) 'gsv_horizSubSample: new      subsample level=', statevector_out%horizSubSample
        call utl_abort('gsv_horizSubSample: relative change of subsample level not an integer!')
      end if
      relativeFactor = nint(ratio_r8)

      ! copy each input grid point to multiple output grid points
!$OMP PARALLEL DO PRIVATE(stepIndex, latIndex, ilat_out, ilat_in, lonIndex, ilon_out, ilon_in)
      do stepIndex = 1, statevector_out%numStep

        do latIndex = 1, statevector_out%latPerPE
          ilat_out = latIndex + statevector_out%myLatBeg - 1
          ilat_in  = floor(real(latIndex-1,8)/real(relativeFactor,8)) + statevector_in%myLatBeg
          do lonIndex = 1, statevector_out%lonPerPE
            ilon_out = lonIndex + statevector_in%myLonBeg - 1
            ilon_in  = floor(real(lonIndex-1,8)/real(relativeFactor,8)) + statevector_in%myLonBeg
            statevector_out%gd_r8(ilon_out, ilat_out, :, stepIndex) =  &
              statevector_in%gd_r8(ilon_in, ilat_in, :, stepIndex)
          end do ! lonIndex
        end do ! latIndex

      end do ! stepIndex
!$OMP END PARALLEL DO 

    end if

  end subroutine gsv_horizSubSample

  !--------------------------------------------------------------------------
  ! gsv_readTrials
  !--------------------------------------------------------------------------
  subroutine gsv_readTrials(hco_in, vco_in,statevector_trial, HUcontainsLQ_opt, &
                            hInterpolateDegree_opt)
  !
  ! Author: Y. Rochon, Feb 2017 (addition recommended by Mark Buehner)
  !         Bulk of content originally in vtr_setupTrials.
  ! 
  !---------------------------------------------------------------------------
    implicit none

    type(struct_hco), pointer :: hco_in
    type(struct_vco), pointer :: vco_in
    type(struct_gsv)          :: statevector_trial
    logical, optional         :: HUcontainsLQ_opt
    character(len=*),optional :: hInterpolateDegree_opt

    integer              :: fnom, fstouv, fclos, fstfrm, fstinf
    integer              :: ierr, ikey, stepIndex, trialIndex, numTrials, nulTrial
    integer              :: ni_file, nj_file, nk_file, dateStamp
    integer, allocatable :: dateStampList(:)
    character(len=2)     :: fileNumber
    character(len=30)    :: fileName
    logical              :: fileExists, HUcontainsLQ

    if ( present(HUcontainsLQ_opt) ) then
      HUcontainsLQ = HUcontainsLQ_opt
    else
      HUcontainsLQ = .true.
    end if

    ! initialize statevector_trial
    call gsv_allocate(statevector_trial, tim_nstepobsinc, hco_in, vco_in,     &
                      dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                      allocGZsfc_opt=.true., hInterpolateDegree_opt=hInterpolateDegree_opt)

    ! initialize list of dates for the 4D analysis increment
    allocate(datestamplist(tim_nStepObsInc))
    call tim_getstamplist(datestamplist,tim_nStepObsInc,tim_getDatestamp())

    ! check existence of all trial files
    numTrials = 0
    do 
      write(fileNumber,'(I2.2)') numTrials+1
      fileName = './trlm_' // trim(fileNumber)
      inquire(file=trim(fileName),exist=fileExists)
      if (fileExists) then
        numTrials = numTrials + 1
      elseif ( (.not. fileExists) .and. numTrials > 0 ) then
        exit  
      elseif ( (.not. fileExists) .and. numTrials == 0 ) then
        call utl_abort('gsv_readTrials: No trial files found')
      end if
    end do

    if (numTrials.ne.tim_nstepobs) then
      write(*,*) 'numTrials, tim_nstepobs = ',numTrials, tim_nstepobs
      call utl_abort('gsv_readTrials: numTrials /= tim_nstepobs')
    end if

    ! loop over times for which increment is computed
    do stepIndex = 1, tim_nstepobsinc
      dateStamp = dateStampList(stepIndex)
      if (mpi_myid.eq.0) write(*,*) 'gsv_readTrials: reading background for time step: ',stepIndex, dateStamp

      ! identify which trial file corresponds with current datestamp
      ikey = 0
      do trialIndex = 1, numTrials
        write(fileNumber,'(I2.2)') trialIndex
        fileName = './trlm_' // trim(fileNumber)
        nulTrial = 0
        ierr = fnom(nulTrial,trim(fileName),'RND+OLD+R/O',0)
        ierr = fstouv(nulTrial,'RND+OLD')
        ikey = fstinf(nulTrial, ni_file, nj_file, nk_file,  &
               dateStamp, ' ', -1, -1, -1, ' ', 'P0')
        ierr = fstfrm(nulTrial)
        ierr = fclos(nulTrial)
        if (ikey.gt.0) exit
      end do

      if (ikey.le.0) then 
        write(*,*) 'stepIndex, dateStamp = ', stepIndex, dateStamp
        call utl_abort('gsv_readTrials: trial file not found for this increment timestep')
      end if

      ! read the trial file for this timestep
      call gsv_readFromFile(statevector_trial, fileName, ' ', 'P', stepIndex,  &
                            HUcontainsLQ_opt=HUcontainsLQ, readGZsfc_opt=.true.)

    end do ! stepIndex

  end subroutine gsv_readTrials

  !--------------------------------------------------------------------------
  ! gsv_varKindExist
  !--------------------------------------------------------------------------
  function gsv_varKindExist(varKind) result(KindFound)
  !
  ! Author: M. Sitwell Sept 2015
  !
  ! Purpose: Checks if any of the variables to be assimilated (i.e. specified in the
  !          namelist NAMSTATE) are part of the specified variable kind
  ! 
  ! IN
  ! 
  !    varKind      Variable kind (e.g. MT or CH)
  !  
  ! OUT
  !
  !   KindFound     Logical indicating if var kind found 
  !
  !-------------------------------------------------------------------------

    implicit none
    character(len=*) :: varKind
    logical :: KindFound
    
    integer :: varIndex

    KindFound = .false.

    VARLIST: do varIndex = 1, vnl_numvarmax
       if (gsv_varExist(varName=vnl_varNameList(varIndex))) then
          if (vnl_varKindFromVarname(vnl_varNameList(varIndex)).eq.varKind) then
             KindFound = .true.
             exit VARLIST 
          end if
       end if
    end do VARLIST

  end function gsv_varKindExist

  !--------------------------------------------------------------------------
  ! gsv_multEnergyNorm
  !--------------------------------------------------------------------------
  SUBROUTINE gsv_multEnergyNorm(statevector_inout, statevector_ref )
    implicit none
    type(struct_gsv)     :: statevector_inout, statevector_ref
    integer              :: jstep, jlon, jlev, jlat, jlon2, jlat2, status, nLev_M, nLev_T
    real(8)              :: scaleFactor, scaleFactorConst, scaleFactorLat, scaleFactorLon, scaleFactorLev
    real(8)              :: pfac, tfac, qfac
    real(8)              :: sumScale , sumeu, sumev, sumep, sumet, sumeq
    real(8), pointer     :: field_UU(:,:,:,:), field_VV(:,:,:,:), field_T(:,:,:,:), field_LQ(:,:,:,:)
    real(8), pointer     :: field_Psfc(:,:,:,:), field_TG(:,:,:,:),Psfc_ptr(:,:,:)
    real(8), pointer     :: Press_T(:,:,:) 
    real(8), pointer     :: Press_M(:,:,:)
    real(8), allocatable :: Psfc_ref(:,:)
    real(8), parameter   :: T_r = 280.0D0
    real(8), parameter   :: Psfc_r = 100000.0D0 ! unit Pa

    if (mpi_myid == 0) write(*,*) 'gsv_multEnergyNorm: START'
    nullify(Press_T,Press_M)

    ! the factors for TT, HU and Ps (for wind is 1)
    tfac = MPC_CP_DRY_AIR_R8/T_r                                 ! temperature factor (c_p/T_r)
    qfac = MPC_HEAT_CONDENS_WATER_R8**2/(MPC_CP_DRY_AIR_R8*T_r)  ! humidity factor ( (l_p*l_p)/(c_p*T_r) )
    pfac = MPC_RGAS_DRY_AIR_R8*T_r/(Psfc_r**2)                   ! surface pressure factor (R*T_r/Psfc_r^2)

    if (.not.statevector_inout%allocated) then
      call utl_abort('gsv_multEnergyNorm: gridStateVector_inout not yet allocated! Aborting.')
    end if

    nLev_M = gsv_getNumLev(statevector_inout,'MM')
    nLev_T = gsv_getNumLev(statevector_inout,'TH')

    ! compute 3D log pressure fields
    Psfc_ptr => gsv_getField3D_r8(statevector_ref,'P0')
    allocate(Psfc_ref(statevector_inout%lonPerPEmax,statevector_inout%latPerPEmax))
    Psfc_ref(:,:) =  &
                  Psfc_ptr(statevector_inout%myLonBeg:statevector_inout%myLonEnd,  &
                  statevector_inout%myLatBeg:statevector_inout%myLatEnd, 1)
    status = vgd_levels(statevector_inout%vco%vgrid, &
                        ip1_list=statevector_inout%vco%ip1_T,  &
                        levels=Press_T,   &
                        sfc_field=Psfc_ref,      &
                        in_log=.false.)
    status = vgd_levels(statevector_inout%vco%vgrid, &
                        ip1_list=statevector_inout%vco%ip1_M,  &
                        levels=Press_M,   &
                        sfc_field=Psfc_ref,      &
                        in_log=.false.)
    ! dlat * dlon
    scaleFactorConst = statevector_inout%hco%dlat*statevector_inout%hco%dlon

    ! for wind components
    field_UU => gsv_getField_r8(statevector_inout,'UU')
    field_VV => gsv_getField_r8(statevector_inout,'VV')
    sumeu = 0.0D0
    sumev = 0.0D0
    sumScale = 0.0D0  
    do jlev = 1, nLev_M
      do jstep = 1, statevector_inout%numStep
        do jlat = statevector_inout%myLatBeg, statevector_inout%myLatEnd
          jlat2 = jlat - statevector_inout%myLatBeg + 1
          scaleFactorLat = cos(statevector_inout%hco%lat(jlat))

          do jlon = statevector_inout%myLonBeg, statevector_inout%myLonEnd
            jlon2 = jlon - statevector_inout%myLonBeg + 1
            ! do all thermo levels for which there is a momentum level above and below
            if ( jlev  ==  nLev_M) then
              scaleFactorLev = Press_M(jlon2, jlat2, nLev_M)-Press_T(jlon2, jlat2, nLev_T-1) 
            else if ( Press_T(jlon2, jlat2, jlev) < 10000.0D0) then 
              scaleFactorLev = 0.0D0
            else
              scaleFactorLev = Press_T(jlon2, jlat2, jlev+1) -  Press_T(jlon2, jlat2, jlev)
            end if

              scaleFactor = scaleFactorConst * scaleFactorLat * scaleFactorLev
              sumScale = sumScale + scaleFactor

              sumeu = sumeu + &
                      0.5 * field_UU(jlon,jlat,jlev,jstep) * field_UU(jlon,jlat,jlev,jstep) * scaleFactor
              sumev = sumev + &
                      0.5 * field_VV(jlon,jlat,jlev,jstep) * field_VV(jlon,jlat,jlev,jstep) * scaleFactor

              field_UU(jlon,jlat,jlev,jstep) = &
                   field_UU(jlon,jlat,jlev,jstep) * 0.5 * scaleFactor
              field_VV(jlon,jlat,jlev,jstep) = &
                   field_VV(jlon,jlat,jlev,jstep) * 0.5 * scaleFactor
          end do !jlon
        end do !jlat
      end do ! jstep
    end do ! jlev

    call mpi_allreduce_sumreal8scalar(sumeu,"GRID")
    call mpi_allreduce_sumreal8scalar(sumev,"GRID")
    call mpi_allreduce_sumreal8scalar(sumScale,"GRID")

    sumeu = sumeu/sumScale
    sumev = sumev/sumScale

    if (mpi_myid == 0)  write(*,*) 'energy for UU=', sumeu
    if (mpi_myid == 0)  write(*,*) 'energy for VV=', sumev

    field_UU(:,:,:,:) = field_UU(:,:,:,:)/sumScale
    field_VV(:,:,:,:) = field_VV(:,:,:,:)/sumScale

    ! for Temperature
    field_T => gsv_getField_r8(statevector_inout,'TT')
    sumScale = 0.0D0
    sumet = 0.0D0

    do jlev = 1, nLev_T
      do jstep = 1, statevector_inout%numStep
        do jlat = statevector_inout%myLatBeg, statevector_inout%myLatEnd
          jlat2 = jlat - statevector_inout%myLatBeg + 1
          scaleFactorLat = cos(statevector_inout%hco%lat(jlat))

          ! do all thermo levels for which there is a momentum level above and below
          do jlon = statevector_inout%myLonBeg, statevector_inout%myLonEnd
            jlon2 = jlon - statevector_inout%myLonBeg + 1

            if (jlev == nLev_T) then  !surface
              scaleFactorLev =  Press_T(jlon2, jlat2, nLev_T)-Press_T(jlon2, jlat2, nLev_T-1) 
            else if (jlev == 1)  then  ! top
              scaleFactorLev = 0.0D0
            else if ( Press_M(jlon2, jlat2, jlev-1) < 10000.0D0) then 
              scaleFactorLev = 0.0D0
            else
              scaleFactorLev = Press_M(jlon2, jlat2, jlev ) - Press_M(jlon2, jlat2, jlev-1)
            end if
              scaleFactor = scaleFactorConst * scaleFactorLat * scaleFactorLev
              sumet = sumet + &
                   0.5 * tfac * field_T(jlon,jlat,jlev,jstep) * field_T(jlon,jlat,jlev,jstep) * scaleFactor
              sumScale = sumScale + scaleFactor
              field_T(jlon,jlat,jlev,jstep) = &
                           field_T(jlon,jlat,jlev,jstep) * 0.5 * tfac * scaleFactor
          end do
        end do
      end do ! jstep
    end do ! jlev
    call mpi_allreduce_sumreal8scalar(sumet,"GRID")
    call mpi_allreduce_sumreal8scalar(sumScale,"GRID")
    sumet = sumet/sumScale
    if (mpi_myid == 0)  write(*,*) 'energy for TT=', sumet
    field_T(:,:,:,:) = field_T(:,:,:,:)/sumScale

    ! humidity (set to zero, for now)
    field_LQ => gsv_getField_r8(statevector_inout,'HU')
    sumScale = 0.0D0
    sumeq = 0.0D0

    do jlev = 1, nLev_T
      do jstep = 1, statevector_inout%numStep
        do jlat = statevector_inout%myLatBeg, statevector_inout%myLatEnd
          jlat2 = jlat - statevector_inout%myLatBeg + 1
          scaleFactorLat = cos(statevector_inout%hco%lat(jlat))
          ! do all thermo levels for which there is a momentum level above and below
          do jlon = statevector_inout%myLonBeg, statevector_inout%myLonEnd
            jlon2 = jlon - statevector_inout%myLonBeg + 1

            if ( jlev == nLev_T) then !surface
              scaleFactorLev =  Press_T(jlon2, jlat2, nLev_T) - Press_T(jlon2, jlat2, nLev_T-1) 
            else if (jlev == 1)  then  ! top
              scaleFactorLev = 0.0D0
            else if ( Press_M(jlon2, jlat2, jlev-1) < 10000.0D0) then 
              scaleFactorLev = 0.0D0
            else
              scaleFactorLev = Press_M(jlon2, jlat2, jlev ) - Press_M(jlon2, jlat2, jlev-1)
            end if

            scaleFactor = scaleFactorConst * scaleFactorLat * scaleFactorLev
            sumScale = sumScale + scaleFactor

            sumeq = sumeq + 0.5 * qfac * &
                    field_LQ(jlon,jlat,jlev,jstep) * field_LQ(jlon,jlat,jlev,jstep) * scaleFactor

            field_LQ(jlon,jlat,jlev,jstep) = &
                       field_LQ(jlon,jlat,jlev,jstep) * 0.5 * scaleFactor * qfac * 0.0

          end do
        end do
      end do ! jstep
    end do ! jlat
    call mpi_allreduce_sumreal8scalar(sumScale,"GRID")
    field_LQ(:,:,:,:) = field_LQ(:,:,:,:)/sumScale*0.0

    ! surface pressure
    field_Psfc => gsv_getField_r8(statevector_inout,'P0')
    sumScale = 0.0D0
    sumep = 0.0

    do jstep = 1, statevector_inout%numStep
      do jlat = statevector_inout%myLatBeg, statevector_inout%myLatEnd
        scaleFactorLat = cos(statevector_inout%hco%lat(jlat))
        do jlon = statevector_inout%myLonBeg, statevector_inout%myLonEnd
          scaleFactor = scaleFactorConst * scaleFactorLat
          sumScale = sumScale + scaleFactor
          sumep = sumep + 0.5 * pfac * &
                  field_Psfc(jlon,jlat,1,jstep) * field_Psfc(jlon,jlat,1,jstep) * scaleFactor
          field_Psfc(jlon,jlat,1,jstep) = &
            field_Psfc(jlon,jlat,1,jstep) * 0.5 * scaleFactor * pfac
        end do
      end do ! jlat
    end do ! jstep

    call mpi_allreduce_sumreal8scalar(sumep,"GRID")
    call mpi_allreduce_sumreal8scalar(sumScale,"GRID")

    sumep = sumep/sumScale

    if (mpi_myid == 0)  write(*,*) 'energy for Ps=', sumep

    field_Psfc(:,:,:,:) =  field_Psfc(:,:,:,:)/sumScale

    ! skin temperature (set to zero for now)
    field_TG => gsv_getField_r8(statevector_inout,'TG')
    sumScale = 0.0D0
    do jstep = 1, statevector_inout%numStep
      do jlat = statevector_inout%myLatBeg, statevector_inout%myLatEnd
        scaleFactorLat = cos(statevector_inout%hco%lat(jlat))
        do jlon = statevector_inout%myLonBeg, statevector_inout%myLonEnd
          scaleFactor = scaleFactorConst * scaleFactorLat
          sumScale = sumScale + scaleFactor
          field_TG(jlon,jlat,1,jstep) = &
                  field_TG(jlon,jlat,1,jstep) * 0.5 * scaleFactor * 0.0
        end do
      end do ! jlat
    end do ! jstep

    call mpi_allreduce_sumreal8scalar(sumScale,"GRID")

    field_TG(:,:,:,:) = field_TG(:,:,:,:)/sumScale * 0.0
    if (mpi_myid == 0) write(*,*) 'energy for total=', sumeu + sumev + sumet + sumep

    deallocate(Press_T,Press_M)
    deallocate(Psfc_ref)

    if (mpi_myid == 0) write(*,*) 'gsv_multEnergyNorm: END'

  END SUBROUTINE gsv_multEnergyNorm

  !--------------------------------------------------------------------------
  ! gsv_dotProduct
  !--------------------------------------------------------------------------
  SUBROUTINE gsv_dotProduct(statevector_a,statevector_b,dotsum)
    implicit none

    type(struct_gsv) :: statevector_a,statevector_b
    real(8)          :: dotsum
    integer          :: jstep,jlon,jlev,jlat,lon1,lon2,lat1,lat2
    integer          :: k1,k2

    if (.not.statevector_a%allocated) then
      call utl_abort('gsv_dotProduct: gridStateVector_in not yet allocated! Aborting.')
    end if
    if (.not.statevector_b%allocated) then
      call utl_abort('gsv_dotProduct: gridStateVector_inout not yet allocated! Aborting.')
    end if

    lon1 = statevector_a%myLonBeg
    lon2 = statevector_a%myLonEnd
    lat1 = statevector_a%myLatBeg
    lat2 = statevector_a%myLatEnd
    k1 = statevector_a%mykBeg
    k2 = statevector_a%mykEnd

    dotsum = 0.0D0
    do jstep = 1, statevector_a%numStep
      do jlev = k1,k2
        do jlat = lat1, lat2
          do jlon = lon1, lon2
            dotsum = dotsum + statevector_a%gd_r8(jlon,jlat,jlev,jstep) * &
                              statevector_b%gd_r8(jlon,jlat,jlev,jstep)
          end do 
        end do
      end do
    end do

    call mpi_allreduce_sumreal8scalar(dotsum,"GRID")

  END SUBROUTINE gsv_dotProduct

end module gridStateVector_mod

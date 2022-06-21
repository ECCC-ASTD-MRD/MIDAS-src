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

module gridStateVector_mod
  ! MODULE gridStateVector_mod (prefix='gsv' category='6. High-level data objects')
  !
  ! :Purpose: The grid-point state vector and related information.
  !
  use mpi, only : mpi_status_size ! this is the mpi library module
  use codePrecision_mod
  use midasMpi_mod
  use earthConstants_mod
  use varNameList_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use oceanMask_mod
  use mathPhysConstants_mod
  use timeCoord_mod
  use utilities_mod
  use message_mod
  use physicsFunctions_mod
  implicit none
  save
  private

  ! public structure definition
  public :: gsv_rhumin, struct_gsv
  public :: gsv_minValVarKindCH

  ! public subroutines and functions
  public :: gsv_setup, gsv_allocate, gsv_deallocate, gsv_zero, gsv_3dto4d, gsv_3dto4dAdj
  public :: gsv_getOffsetFromVarName, gsv_getLevFromK, gsv_getVarNameFromK, gsv_getMpiIdFromK, gsv_hPad
  public :: gsv_modifyVarName, gsv_modifyDate
  public :: gsv_transposeTilesToStep, gsv_transposeStepToTiles, gsv_transposeTilesToMpiGlobal
  public :: gsv_transposeTilesToVarsLevs, gsv_transposeTilesToVarsLevsAd
  public :: gsv_transposeVarsLevsToTiles
  public :: gsv_getField, gsv_getFieldUV
  public :: gsv_getHeightSfc, gsv_isAssocHeightSfc
  public :: gsv_getDateStamp, gsv_getNumLev, gsv_getNumLevFromVarName
  public :: gsv_add, gsv_power, gsv_scale, gsv_scaleVertical, gsv_copy, gsv_copy4Dto3D
  public :: gsv_copyHeightSfc, gsv_copyMask
  public :: gsv_getVco, gsv_getHco, gsv_getHco_physics, gsv_getDataKind, gsv_getNumK
  public :: gsv_horizSubSample
  public :: gsv_varKindExist, gsv_varExist, gsv_varNamesList
  public :: gsv_dotProduct, gsv_schurProduct
  public :: gsv_field3d_hbilin, gsv_smoothHorizontal
  public :: gsv_communicateTimeParams, gsv_resetTimeParams, gsv_getInfo, gsv_isInitialized
  public :: gsv_applyMaskLAM, gsv_containsNonZeroValues
  public :: gsv_isAllocated
  public :: gsv_transposesteptovarslevs

  ! public module variables
  public :: gsv_conversionVarKindCHtoMicrograms

  interface gsv_getField
    module procedure gsv_getFieldWrapper_r4
    module procedure gsv_getFieldWrapper_r8
    module procedure gsv_getField3D_r4
    module procedure gsv_getField3D_r8
  end interface gsv_getField

  interface gsv_getFieldUV
    module procedure gsv_getFieldUVWrapper_r4
    module procedure gsv_getFieldUVWrapper_r8
  end interface gsv_getFieldUV

  type struct_gdUV
    real(8), pointer :: r8(:,:,:) => null()
    real(4), pointer :: r4(:,:,:) => null()
  end type struct_gdUV

  type struct_gsv
    ! This is the derived type of the statevector object

    ! These are the main data storage arrays
    logical, private          :: allocated=.false.
    real(8), pointer, private :: gd_r8(:,:,:,:) => null()
    real(8), pointer, private :: gd3d_r8(:,:,:) => null()
    real(4), pointer, private :: gd_r4(:,:,:,:) => null()
    real(4), pointer, private :: gd3d_r4(:,:,:) => null()
    type(struct_ocm)    :: oceanMask
    logical             :: heightSfcPresent = .false.
    real(8), pointer, private :: heightSfc(:,:) => null()  ! for VarsLevs, heightSfc only on proc 0

    ! These are used when distribution is VarLevs to keep corresponding UV
    ! components together on each mpi task to facilitate horizontal interpolation
    logical             :: UVComponentPresent = .false.  ! wind component present on this mpi task
    logical             :: extraUVallocated = .false.    ! extra winds (gdUV) are allocated
    integer             :: myUVkBeg, myUVkEnd, myUVkCount
    type(struct_gdUV), pointer, private :: gdUV(:) => null()

    ! All the remaining extra information
    integer             :: dataKind = 8 ! default value
    integer             :: ni, nj, nk, numStep, anltime
    integer             :: latPerPE, latPerPEmax, myLatBeg, myLatEnd
    integer             :: lonPerPE, lonPerPEmax, myLonBeg, myLonEnd
    integer             :: mykCount, mykBeg, mykEnd
    integer, pointer    :: allLatBeg(:), allLatEnd(:), allLatPerPE(:)
    integer, pointer    :: allLonBeg(:), allLonEnd(:), allLonPerPE(:)
    integer, pointer    :: allkCount(:), allkBeg(:), allkEnd(:)
    integer, pointer    :: allUVkCount(:), allUVkBeg(:), allUVkEnd(:)
    integer, pointer    :: dateStampList(:) => null()
    integer, pointer    :: dateStamp3d
    integer, pointer    :: dateOriginList(:)
    integer, pointer    :: npasList(:), ip2List(:)
    integer             :: deet
    character(len=12)   :: etiket
    type(struct_vco), pointer :: vco => null()
    type(struct_hco), pointer :: hco => null()
    type(struct_hco), pointer :: hco_physics => null()
    integer, pointer    :: varOffset(:), varNumLev(:)
    logical, pointer    :: onPhysicsGrid(:)
    logical             :: mpi_local=.false.
    character(len=8)    :: mpi_distribution='None'  ! or 'Tiles' or 'VarsLevs'
    integer             :: horizSubSample
    logical             :: varExistList(vnl_numVarMax)
    character(len=12)   :: hInterpolateDegree='UNSPECIFIED' ! or 'LINEAR' or 'CUBIC' or 'NEAREST'
    character(len=12)   :: hExtrapolateDegree='MAXIMUM' ! or 'VALUE' or 'MINIMUM' or 'NEUTRAL'
    logical             :: addHeightSfcOffset = .false.
  end type struct_gsv

  logical :: varExistList(vnl_numVarMax)
  character(len=8) :: ANLTIME_BIN
  integer, external :: get_max_rss
  real(8) :: rhumin, gsv_rhumin
  logical :: addHeightSfcOffset ! controls adding non-zero height offset to diag levels
  logical :: abortOnMpiImbalance

  ! Min values imposed for input trial and output analysis (and related increment)
  ! for variables of CH kind of the AnlVar list.
  real(8) :: minValVarKindCH(vnl_numVarMax), gsv_minValVarKindCH(vnl_numVarMax)
  ! Logical to turn on unit conversion for variables of CH kind of the AnlVar list
  ! when unitConversion=.true.
  logical :: gsv_conversionVarKindCHtoMicrograms
     
  ! arrays used for transpose VarsLevs <-> Tiles
  real(4), allocatable :: gd_send_varsLevs_r4(:,:,:,:), gd_recv_varsLevs_r4(:,:,:,:)
  real(8), allocatable :: gd_send_varsLevs_r8(:,:,:,:), gd_recv_varsLevs_r8(:,:,:,:)

  ! initialized 
  logical :: initialized = .false.

  contains

  !--------------------------------------------------------------------------
  ! gsv_getOffsetFromVarName
  !--------------------------------------------------------------------------
  function gsv_getOffsetFromVarName(statevector,varName) result(offset)
    !
    ! :Purpose: Returns the offset for the given variable provided it exists
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)  :: statevector
    character(len=*), intent(in)  :: varName
    integer                       :: offset

    ! Locals:
    integer :: varIndex

    varIndex = vnl_varListIndex(varName)
    if (.not. statevector%varExistList(varIndex)) then
      call utl_abort('gsv_getOffsetFromVarName: specified varName does not exist in stateVector: ' // trim(varName))
    end if
    offset=statevector%varOffset(varIndex)

  end function gsv_getOffsetFromVarName

  !--------------------------------------------------------------------------
  ! gsv_getVarNameFromK
  !--------------------------------------------------------------------------
  function gsv_getVarNameFromK(statevector,kIndex) result(varName)
    !
    ! :Purpose: Returns the variable name from a given kIndex
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(in)  :: statevector
    integer,          intent(in)  :: kIndex
    character(len=4)              :: varName

    ! Locals:
    integer :: varIndex

    do varIndex = 1, vnl_numvarmax
      if (statevector%varExistList(varIndex)) then
        if ((kIndex >= (statevector%varOffset(varIndex) + 1)) .and.  &
            (kIndex <= (statevector%varOffset(varIndex) + statevector%varNumLev(varIndex)))) then
          varName = vnl_varNameList(varIndex)
          return
        end if
      end if
    end do

    call utl_abort('gsv_getVarNameFromK: kIndex out of range: '//str(kIndex))

  end function gsv_getVarNameFromK

  !--------------------------------------------------------------------------
  ! gsv_getLevFromK
  !--------------------------------------------------------------------------
  function gsv_getLevFromK(statevector,kIndex) result(levIndex)
    !
    ! :Purpose: Returns level index from a given kIndex
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(in) :: statevector
    integer,          intent(in) :: kIndex
    integer                      :: levIndex

    ! Locals
    integer             :: varIndex

    do varIndex = 1, vnl_numvarmax
      if (statevector%varExistList(varIndex)) then
        if ((kIndex >= (statevector%varOffset(varIndex) + 1)) .and.  &
            (kIndex <= (statevector%varOffset(varIndex) + statevector%varNumLev(varIndex)))) then
          levIndex = kIndex - statevector%varOffset(varIndex)
          return
        end if
      end if
    end do

    call utl_abort('gsv_getLevFromK: kIndex out of range: '//str(kIndex))

  end function gsv_getLevFromK

  !--------------------------------------------------------------------------
  ! gsv_getMpiIdFromK
  !--------------------------------------------------------------------------
  function gsv_getMpiIdFromK(statevector,kIndex) result(MpiId)
    !
    ! :Purpose: Returns MPI id from the given kIndex
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in) :: statevector
    integer,          intent(in) :: kIndex
    integer                      :: MpiId
    
    ! Locals:
    integer             :: procIndex

    do procIndex = 1, mmpi_nprocs
      if ((kIndex >= statevector%allKBeg(procIndex)) .and.  &
          (kIndex <= statevector%allKEnd(procIndex))) then
          MpiId = procIndex - 1
          return
      end if
    end do

    call utl_abort('gsv_getMpiIdFromK: kIndex out of range: '//str(kIndex))

  end function gsv_getMpiIdFromK

  !--------------------------------------------------------------------------
  ! gsv_varExist
  !--------------------------------------------------------------------------
  recursive function gsv_varExist(statevector_opt,varName) result(varExist)
    !
    ! :Purpose: Boolean fonction returning .true. if the queried variable
    !           exists in the statevector if provided or in the global variable
    !           list otherwise.
    !           For 'Z_*' and 'P_*' variables, the statevector argument is
    !           mandatory.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), optional, intent(in)  :: statevector_opt
    character(len=*),           intent(in)  :: varName
    logical                                 :: varExist 

    if (varName == 'Z_*') then
      varExist =  gsv_varExist(statevector_opt, 'Z_T') .and. &
                  gsv_varExist(statevector_opt, 'Z_M')
    else if (varName == 'P_*') then
      varExist =  gsv_varExist(statevector_opt, 'P_T') .and. &
                  gsv_varExist(statevector_opt, 'P_M')
    else
      if (present(statevector_opt)) then
        if (statevector_opt%varExistList(vnl_varListIndex(varName))) then
          varExist = .true.
        else
          varExist = .false.
        end if
      else
        if (varExistList(vnl_varListIndex(varName))) then
          varExist = .true.
        else
          varExist = .false.
        end if
      end if
    end if

  end function gsv_varExist

  !--------------------------------------------------------------------------
  ! gsv_varNamesList
  !--------------------------------------------------------------------------
  subroutine gsv_varNamesList(varNames,statevector)
    !
    ! :Purpose: Lists all variables present in the statevector 
    !
    implicit none
    
    ! Arguments:
    character(len=4), pointer,  intent(inout) :: varNames(:)
    type(struct_gsv), optional, intent(in)    :: statevector
    
    ! Locals:
    integer :: varLevIndex, varNumberIndex, varIndex, numFound
    character(len=4) :: varName

    if (associated(varNames)) then
      call utl_abort('gsv_varNamesList: varNames must be NULL pointer on input')
    end if
 
    !
    !- 1. How many variables do we have?
    !
    numFound = 0
    if (present(statevector)) then
      do varIndex = 1, vnl_numvarmax
        if (gsv_varExist(statevector,vnl_varNameList(varIndex))) numFound = numFound + 1
      end do
    else
      do varIndex = 1, vnl_numvarmax
        if (varExistList(varIndex)) numFound = numFound + 1
      end do
    end if

    !
    !- 2. List the variables
    !
    allocate(varNames(numFound))
    varNames(:) = ''

    varNumberIndex = 0
    if (present(statevector)) then
      !- 2.1 List the variables based on the varLevIndex ordering
      do varLevIndex = 1, statevector%nk
        varName = gsv_getVarNameFromK(statevector,varLevIndex)
        if (.not. ANY(varNames(:) == varName)) then
          varNumberIndex = varNumberIndex + 1
          varNames(varNumberIndex) = varName
        end if
      end do
    else
      !- 2.2 List the variables based on the varnamelist_mod ordering
      do varIndex = 1, vnl_numvarmax
        if (varExistList(varIndex)) then
          varName = vnl_varNameList(varIndex)
          if (.not. ANY(varNames(:) == varName)) then
            varNumberIndex = varNumberIndex + 1
            varNames(varNumberIndex) = varName
          end if
        end if
      end do
    end if

  end subroutine gsv_varNamesList

  !--------------------------------------------------------------------------
  ! gsv_getNumLev
  !--------------------------------------------------------------------------
  function gsv_getNumLev(statevector,varLevel,varName_opt) result(nlev)
    !
    ! :Purpose: Returns the number of levels for a given type of variable;
    !           varLevel can be one of 'TH', 'MM', 'SF', 'SFMM', 'SFTH', 'DP',
    !           'SFDP' or 'OT'.
    !
    implicit none

    ! Arguments:
    type(struct_gsv),           intent(in) :: statevector ! Input statevector
    character(len=*),           intent(in) :: varLevel    ! Variable type in 'TH', 'MM', 'SF', 'SFMM', 'SFTH', 'DP', 'SFDP' or 'OT'
    character(len=*), optional, intent(in) :: varName_opt ! Variable name when varLevel='OT'

    ! Locals:
    integer                       :: nlev

    nlev = vco_getNumLev(statevector%vco,varLevel,varName_opt)

  end function gsv_getNumLev

  !--------------------------------------------------------------------------
  ! gsv_getNumK
  !--------------------------------------------------------------------------
  function gsv_getNumK(statevector) result(numK)
    !
    ! :Purpose: Returns the number of k indexes on the current MPI process
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)  :: statevector
    integer                       :: numK

    numK = 1 + statevector%mykEnd - statevector%mykBeg

  end function gsv_getNumK

  !--------------------------------------------------------------------------
  ! gsv_getDataKind
  !--------------------------------------------------------------------------
  function gsv_getDataKind(statevector) result(dataKind)
    !
    ! :Purpose: Returns the real kind (4 or 8 bytes floating point value) of 
    !           the input statevector
    !
    implicit none

    ! arguments
    type(struct_gsv), intent(in)  :: statevector
    integer                       :: dataKind

    dataKind = statevector%dataKind

  end function gsv_getDataKind

  !--------------------------------------------------------------------------
  ! gsv_getNumLevFromVarName
  !--------------------------------------------------------------------------
  function gsv_getNumLevFromVarName(statevector,varName) result(nlev)
    !
    ! :Purpose: Returns the number of levels for a given variable
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)  :: statevector
    character(len=*), intent(in)  :: varName
    integer                       :: nlev

    nlev = statevector%varNumLev(vnl_varListIndex(varName))

  end function gsv_getNumLevFromVarName

  !--------------------------------------------------------------------------
  ! gsv_setup
  !--------------------------------------------------------------------------
  subroutine gsv_setup
    !
    ! :Purpose: Initialises the gridstatevector module global structure.
    !
    ! :Namelist parameters:
    !         :anlvar:          Analysis variable (chemistry analysis only)
    !
    !         :rhumin:          Minimum relative humidity value
    !
    !         :anlTime_bin:     Analysis time reference ('MIDDLE', 'FIRST' or 'LAST')
    !
    !         :addHeightSfcOffset:
    !                           Global statevector height offset  
    !
    !         :conversionVarKindCHtoMicrograms:
    !                           If .true. will apply some unit conversion when
    !                           writing to file (chemistry analysis only)
    !
    !         :minValVarKindCH: Minimal values imposed for input trial and 
    !                           output analysis and related increment for 
    !                           variables of CH kind (chemistry analysis only)
    !                           
    !         :abortOnMpiImbalance:
    !                           If .true., will abort when MPI topology is 
    !                           inappropriate
    !
    implicit none

    ! Locals:
    logical           :: conversionVarKindCHtoMicrograms
    integer           :: varIndex, fnom, fclos, nulnam, ierr, loopIndex
    character(len=4)  :: anlvar(vnl_numVarMax)

    NAMELIST /NAMSTATE/anlvar, rhumin, anlTime_bin, addHeightSfcOffset, &
                       conversionVarKindCHtoMicrograms, minValVarKindCH, &
                       abortOnMpiImbalance

    if (initialized) return

    call msg('gsv_setup', 'List of known (valid) variable names', mpiAll_opt=.false.)
    call msg('gsv_setup', 'varNameList3D   ='//str(vnl_varNameList3D(:)), mpiAll_opt=.false.)
    call msg('gsv_setup', 'varNameList3D   ='//str(vnl_varNameList3D(:)), mpiAll_opt=.false.)
    call msg('gsv_setup', 'varNameListOther='//str(vnl_varNameListOther(:)), mpiAll_opt=.false.)

    ! Read namelist NAMSTATE to find which fields are needed

    anlvar(:) = '    '
    rhumin = mpc_minimum_hu_r8
    anltime_bin = 'MIDDLE'
    addHeightSfcOffset = .false.
    conversionVarKindCHtoMicrograms = .false.
    minValVarKindCH(:) = mpc_missingValue_r8
    abortOnMpiImbalance = .true.

    nulnam=0
    ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=namstate,iostat=ierr)
    if (ierr.ne.0) call utl_abort('gsv_setup: Error reading namelist NAMSTATE')
    if (mmpi_myid.eq.0) write(*,nml=namstate)
    ierr=fclos(nulnam)

    gsv_rhumin = rhumin
    gsv_conversionVarKindCHtoMicrograms = conversionVarKindCHtoMicrograms

    if (varneed('Z_T') .or. varneed('Z_M')) call utl_abort('gsv_setup: height can not be specified as analysis variable in namelist!')
    if (varneed('P_T') .or. varneed('P_M')) call utl_abort('gsv_setup: pressure can not be specified as analysis variable in namelist!')

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

    do varIndex = 1, vnl_numvarmaxOther
      if (varneed(vnl_varNameListOther(varIndex))) then
        varExistList(varIndex+vnl_numVarMax3D+vnl_numVarMax2D) = .true.
      else
        varExistList(varIndex+vnl_numVarMax3D+vnl_numVarMax2D) = .false.
      end if
    end do

    ! Setup to assign min values to apply
    
    ! Check for input values only for variables of CH kind
    do varIndex = 1, vnl_numvarmax
      if (trim(AnlVar(varIndex)) == '') exit
      if (vnl_varKindFromVarname(AnlVar(varIndex)) == 'CH') then
        if (minValVarKindCH(varIndex) < 0.99d0 * mpc_missingValue_r8) then
          if (trim(AnlVar(varIndex)) == 'AF' .or. trim(AnlVar(varIndex)) == 'AC') then
            ! Set for particulate matter in micrograms/cm^3
            minValVarKindCH(varIndex) = mpc_minimum_pm_r8
          else
            ! Set for concentrations in micrograms/kg
            minValVarKindCH(varIndex) = mpc_minimum_ch_r8
          end if
        end if
      end if
    end do

    ! Assign min values to apply
    gsv_minValVarKindCH(:) = mpc_missingValue_r8
    do varIndex = 1, vnl_numvarmax
      if (varExistList(varIndex)) then
        do loopIndex = 1, vnl_numvarmax
          if (trim(AnlVar(loopIndex)) == '') exit
          if (trim(vnl_varNameList(varIndex)) == trim(AnlVar(loopIndex))) &
             gsv_minValVarKindCH(varIndex) = minValVarKindCH(loopIndex)
        end do
      end if
    end do

    call msg('gsv_setup','global varExistList ='//str(varExistList), mpiAll_opt=.false.)

    ! Check value for ANLTIME_BIN
    if (ANLTIME_BIN .ne. 'MIDDLE' .and. ANLTIME_BIN .ne. 'FIRST' .and.  ANLTIME_BIN .ne. 'LAST') then
      call utl_abort('gsv_setup: Problem setting ANLTIME_BIN. Verify NAMSTATE namelist')
    end if

    initialized = .true.

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
  ! gsv_isInitialized
  !--------------------------------------------------------------------------
  function gsv_isInitialized() result(gsvInitialized)
    !
    ! :Purpose: To verify gsv_setup has already run.
    !
    implicit none

    ! Argument:
    logical :: gsvInitialized

    gsvInitialized = .false.

    if (initialized) gsvInitialized = .true.

  end function gsv_isInitialized

  !--------------------------------------------------------------------------
  ! gsv_isAllocated
  !--------------------------------------------------------------------------
  function gsv_isAllocated(stateVector) result(isAllocated)
    !
    ! :Purpose: To verify if a stateVector is allocated.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)  :: stateVector
    logical                       :: isAllocated

    isAllocated = stateVector%allocated

  end function gsv_isAllocated

  !--------------------------------------------------------------------------
  ! gsv_allocate
  !--------------------------------------------------------------------------
  subroutine gsv_allocate(statevector, numStep, hco_ptr, vco_ptr, dateStamp_opt, dateStampList_opt,  &
                          mpi_local_opt, mpi_distribution_opt, horizSubSample_opt,                   &
                          varNames_opt, dataKind_opt, allocHeightSfc_opt, hInterpolateDegree_opt,    &
                          hExtrapolateDegree_opt, allocHeight_opt, allocPressure_opt, beSilent_opt)
    !
    ! :Purpose: Allocates the struct_gsv memory, sets horizontal and vertical 
    !           coordinates, sets some options and MPI distribution 
    !           configurations
    !
    implicit none

    ! Arguments:
    type(struct_gsv),           intent(inout) :: statevector            ! statevector to be allocated
    integer,                    intent(in)    :: numStep                ! number of time steps
    type(struct_hco), pointer,  intent(in)    :: hco_ptr                ! horizontal structure
    type(struct_vco), pointer,  intent(in)    :: vco_ptr                ! vertical structure
    integer,          optional, intent(in)    :: dateStamp_opt          ! reference datestamp
    integer,          optional, intent(in)    :: dateStampList_opt(:)   ! explicit datestamp list
    logical,          optional, intent(in)    :: mpi_local_opt          ! if .false. no MPI distribution will be used 
    character(len=*), optional, intent(in)    :: mpi_distribution_opt   ! MPI distribution strategy in {'Tiles', 'VarsLevs', 'None'} defaults to 'Tiles'
    integer,          optional, intent(in)    :: horizSubSample_opt     ! horizontal subsampling factor (to get a coarser grid)
    character(len=*), optional, intent(in)    :: varNames_opt(:)        ! allow specification of variables
    integer,          optional, intent(in)    :: dataKind_opt           ! real kind (4 or 8 bytes; defaults to 8)
    logical,          optional, intent(in)    :: allocHeightSfc_opt     ! toggle allocation of surface height field
    logical,          optional, intent(in)    :: allocHeight_opt        ! force the allocation of 'Z_T' and 'Z_M'
    logical,          optional, intent(in)    :: allocPressure_opt      ! force the allocation of 'P_T' and 'P_M'
    character(len=*), optional, intent(in)    :: hInterpolateDegree_opt ! set the horizontal interpolation degree
    character(len=*), optional, intent(in)    :: hExtrapolateDegree_opt ! set the horizontal extrapolation degree
    logical,          optional, intent(in)    :: beSilent_opt           ! limit outputs to listing

    ! Locals:
    integer :: ierr,iloc,varIndex,varIndex2,stepIndex,lon1,lat1,k1,kIndex,kIndex2,levUV
    character(len=4) :: UVname
    logical :: beSilent, allocPressure, allocHeight
    integer :: verbLevel

    call utl_tmg_start(168, 'low-level--gsv_allocate')

    call utl_tmg_start(168, 'low-level--gsv_allocate')

    if (.not. initialized) then
      call msg('gsv_allocate','gsv_setup must be called first to be able to use this module. Call it now')
      call gsv_setup
    end if

    if (present(beSilent_opt)) then
      beSilent = beSilent_opt
    else
      beSilent = .true.
    end if

    ! set the horizontal and vertical coordinates
    statevector%hco => hco_ptr
    statevector%vco => vco_ptr

    if (.not.statevector%vco%initialized) then
       call utl_abort('statevector_allocate: VerticalCoord has not been initialized!')
    end if

    if (statevector%allocated) then
      call msg('gsv_allocate', 'gridStateVector already allocated! Deallocating first.', mpiAll_opt=.false.)
      call gsv_deallocate(statevector)
    end if

    if (present(dataKind_opt)) statevector%dataKind = dataKind_opt

    if (present(varNames_opt)) then
      if (present(allocHeight_opt)) call utl_abort('gsv_allocate: to allocate Z_T/Z_M, set them in varNames_opt')
      if (present(allocPressure_opt)) call utl_abort('gsv_allocate: to allocate P_T/P_M, set them in varNames_opt')
    end if

    if (present(varNames_opt)) then      
      statevector%varExistList(:) = .false.
      do varIndex2 = 1, size(varNames_opt)
        varIndex = vnl_varListIndex(varNames_opt(varIndex2))
        statevector%varExistList(varIndex) = .true.
      end do
    else
      ! use the global variable list and set allocHeight/allocPressure if needed
      statevector%varExistList(:) = varExistList(:)

      if (present(allocHeight_opt)) then
        allocHeight = allocHeight_opt
      else
        if (statevector%varExistList(vnl_varListIndex('TT ')) .and. &
            statevector%varExistList(vnl_varListIndex('HU ')) .and. &
            statevector%varExistList(vnl_varListIndex('P0 '))) then
          allocHeight = .true.
        else
          allocHeight = .false.
        end if
      end if

      if (present(allocPressure_opt)) then
        allocPressure = allocPressure_opt
      else
        if (statevector%varExistList(vnl_varListIndex('P0 '))) then
          allocPressure = .true.
        else
          allocPressure = .false.
        end if
      end if

      ! add Z_T/Z_M and P_T/P_M to the varExistList
      if (allocHeight) then
        if (gsv_getNumLev(statevector,'TH') > 0) statevector%varExistList(vnl_varListIndex('Z_T ')) = .true.
        if (gsv_getNumLev(statevector,'MM') > 0) statevector%varExistList(vnl_varListIndex('Z_M ')) = .true.
      end if
      if (allocPressure) then
        if (gsv_getNumLev(statevector,'TH') > 0) statevector%varExistList(vnl_varListIndex('P_T ')) = .true.
        if (gsv_getNumLev(statevector,'MM') > 0) statevector%varExistList(vnl_varListIndex('P_M ')) = .true.
      end if
    end if

    if (present(horizSubSample_opt)) then
      ! user has chosen a coarser grid than specified in hco
      statevector%horizSubSample = horizSubSample_opt
    else
      ! default is no sub-sampling
      statevector%horizSubSample = 1
    end if

    if (present(hInterpolateDegree_opt)) then
      ! set the horizontal interpolation degree (intentionally no default value)
      statevector%hInterpolateDegree = trim(hInterpolateDegree_opt)
    end if

    if (present(hExtrapolateDegree_opt)) then
      ! set the horizontal extrapolation degree (intentionally no default value)
      statevector%hExtrapolateDegree = trim(hExtrapolateDegree_opt)
    end if

    ! compute the number of global grid points for a given subSample level
    statevector%ni = ceiling(real(statevector%hco%ni,8) / real(statevector%horizSubSample,8))
    statevector%nj = ceiling(real(statevector%hco%nj,8) / real(statevector%horizSubSample,8))

    if (statevector%ni * statevector%horizSubSample /= statevector%hco%ni) then
      call msg('gsv_allocate',' number of longitudes is not evenly divisible at this subSample level'&
               //' ni='//str(statevector%ni)// ', horizSubSample = '&
               //str(statevector%horizSubSample))
      call utl_abort('gsv_allocate')
    end if

    if (statevector%nj * statevector%horizSubSample /= statevector%hco%nj) then
      call msg('gsv_allocate','number of latitudes is not evenly divisible at this subSample level'&
               //' nj='//str(statevector%nj)//', horizSubSample = '&
               //str(statevector%horizSubSample))
      call utl_abort('gsv_allocate')
    end if

    statevector%numStep=numStep

    if (present(mpi_local_opt)) then
      statevector%mpi_local = mpi_local_opt
    else
      statevector%mpi_local = .false.
    end if

    if (present(mpi_distribution_opt)) then
      if (trim(mpi_distribution_opt) .ne. 'Tiles'    .and. &
          trim(mpi_distribution_opt) .ne. 'VarsLevs' .and. &
          trim(mpi_distribution_opt) .ne. 'None') then
        call utl_abort('gsv_allocate: Unknown value of mpi_distribution: ' // trim(mpi_distribution_opt))
      end if
      statevector%mpi_distribution = mpi_distribution_opt
    else
      if (statevector%mpi_local) then
        statevector%mpi_distribution = 'Tiles'
      else
        statevector%mpi_distribution = 'None'
      end if
    end if

    ! determine lat/lon index ranges
    if (stateVector%mpi_distribution == 'Tiles') then
      call mmpi_setup_latbands(statevector%nj,  &
                               statevector%latPerPE, statevector%latPerPEmax, &
                               statevector%myLatBeg, statevector%myLatEnd)
      call mmpi_setup_lonbands(statevector%ni,  &
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
    allocate(statevector%onPhysicsGrid(vnl_numvarmax))
    statevector%onPhysicsGrid(:) = .false.

    iloc=0
    if (present(varNames_opt)) then

      do varIndex2 = 1, size(varNames_opt)
        varIndex = vnl_varListIndex(varNames_opt(varIndex2))
        statevector%varOffset(varIndex)=iloc
        statevector%varNumLev(varIndex)=  &
             gsv_getNumLev(statevector, &
                           vnl_varLevelFromVarname(vnl_varNameList(varIndex)),  &
                           vnl_varNameList(varIndex))
        iloc = iloc + statevector%varNumLev(varIndex)
      end do

    else

      do varIndex = 1, vnl_numvarmax3d
        if (statevector%varExistList(varIndex)) then
          statevector%varOffset(varIndex)=iloc
          statevector%varNumLev(varIndex)=  &
               gsv_getNumLev(statevector,  &
                             vnl_varLevelFromVarname(vnl_varNameList(varIndex)))
          iloc = iloc + statevector%varNumLev(varIndex)
        end if
      end do
      do varIndex2 = 1, vnl_numvarmax2d
        varIndex=varIndex2+vnl_numvarmax3d
        if (statevector%varExistList(varIndex)) then
          statevector%varOffset(varIndex)=iloc
          statevector%varNumLev(varIndex)=1
          iloc = iloc + 1
        end if
      end do
      do varIndex2 = 1, vnl_numvarmaxOther
        varIndex=varIndex2+vnl_numvarmax3d+vnl_numvarmax2d
        if (statevector%varExistList(varIndex)) then
          statevector%varOffset(varIndex)=iloc
          statevector%varNumLev(varIndex)=  &
               gsv_getNumLev(statevector,  &
                              vnl_varLevelFromVarname(vnl_varNameList(varIndex)), &
                              vnl_varNameList(varIndex))
          iloc = iloc + statevector%varNumLev(varIndex)
        end if
      end do

    end if

    if (iloc == 0) then
      call utl_abort('gsv_allocate:  Nothing to allocate')
    end if

    statevector%nk=iloc

    if (beSilent) then
      verbLevel = msg_NEVER
    else
      verbLevel = 2
    end if
    call msg('gsv_allocate', 'statevector%nk = '//str(statevector%nk)&
             //new_line('')//'varOffset='//str(statevector%varOffset)&
             //new_line('')//'varNumLev='//str(statevector%varNumLev),&
             verb_opt=verbLevel, mpiAll_opt=.false.)

    ! determine range of values for the 'k' index (vars+levels)
    if (statevector%mpi_distribution == 'VarsLevs') then
      call mmpi_setup_varslevels(statevector%nk, statevector%mykBeg, &
                                 statevector%mykEnd, statevector%mykCount)
    else
      statevector%mykCount = statevector%nk
      statevector%mykBeg = 1
      statevector%mykEnd = statevector%nk
    end if

    ! determine if a wind component exists on this mpi task
    statevector%UVComponentPresent = .false.
    statevector%myUVkCount = 0
    statevector%myUVkBeg = 0
    statevector%myUVkEnd = -1
    do kIndex = statevector%mykBeg, statevector%mykEnd
      if (gsv_getVarNameFromK(statevector,kIndex) == 'UU' .or.  &
           gsv_getVarNameFromK(statevector,kIndex) == 'VV') then
        statevector%UVComponentPresent = .true.
        if (statevector%myUVkBeg == 0) statevector%myUVkBeg = kIndex
        statevector%myUVkEnd = kIndex
        statevector%myUVkCount = statevector%myUVkCount + 1
      end if
    end do

    ! determine if a separate complementary wind component needed, which
    ! is the case when mpi distribution could mean that both components 
    ! are not available on same mpi task, or the statevector was allocated
    ! with only one of the components
    statevector%extraUVallocated = .false.
    if (statevector%mpi_distribution == 'VarsLevs' .or.   &
         (gsv_varExist(statevector,'UU') .and. .not. gsv_varExist(statevector,'VV')) .or. &
         (gsv_varExist(statevector,'VV') .and. .not. gsv_varExist(statevector,'UU'))) then
      statevector%extraUVallocated = statevector%UVComponentPresent
    end if

    if (statevector%mpi_local) then
      allocate(statevector%allLonBeg(mmpi_npex))
      CALL rpn_comm_allgather(statevector%myLonBeg,1,'mpi_integer',       &
                              statevector%allLonBeg,1,'mpi_integer','EW',ierr)
      allocate(statevector%allLonEnd(mmpi_npex))
      CALL rpn_comm_allgather(statevector%myLonEnd,1,'mpi_integer',       &
                              statevector%allLonEnd,1,'mpi_integer','EW',ierr)
      allocate(statevector%allLonPerPE(mmpi_npex))
      CALL rpn_comm_allgather(statevector%lonPerPE,1,'mpi_integer',       &
                              statevector%allLonPerPE,1,'mpi_integer','EW',ierr)
  
      allocate(statevector%allLatBeg(mmpi_npey))
      CALL rpn_comm_allgather(statevector%myLatBeg,1,'mpi_integer',       &
                              statevector%allLatBeg,1,'mpi_integer','NS',ierr)
      allocate(statevector%allLatEnd(mmpi_npey))
      CALL rpn_comm_allgather(statevector%myLatEnd,1,'mpi_integer',       &
                              statevector%allLatEnd,1,'mpi_integer','NS',ierr)
      allocate(statevector%allLatPerPE(mmpi_npey))
      CALL rpn_comm_allgather(statevector%LatPerPE,1,'mpi_integer',       &
                              statevector%allLatPerPE,1,'mpi_integer','NS',ierr)

      call gsv_checkMpiDistribution(stateVector)

      allocate(statevector%allkCount(mmpi_nprocs))
      CALL rpn_comm_allgather(statevector%mykCount,1,'mpi_integer',       &
                              statevector%allkCount,1,'mpi_integer','grid',ierr)
      allocate(statevector%allkBeg(mmpi_nprocs))
      CALL rpn_comm_allgather(statevector%mykBeg,1,'mpi_integer',       &
                              statevector%allkBeg,1,'mpi_integer','grid',ierr)
      allocate(statevector%allkEnd(mmpi_nprocs))
      CALL rpn_comm_allgather(statevector%mykEnd,1,'mpi_integer',       &
                              statevector%allkEnd,1,'mpi_integer','grid',ierr)

      allocate(statevector%allUVkCount(mmpi_nprocs))
      CALL rpn_comm_allgather(statevector%myUVkCount,1,'mpi_integer',       &
                              statevector%allUVkCount,1,'mpi_integer','grid',ierr)
      allocate(statevector%allUVkBeg(mmpi_nprocs))
      CALL rpn_comm_allgather(statevector%myUVkBeg,1,'mpi_integer',       &
                              statevector%allUVkBeg,1,'mpi_integer','grid',ierr)
      allocate(statevector%allUVkEnd(mmpi_nprocs))
      CALL rpn_comm_allgather(statevector%myUVkEnd,1,'mpi_integer',       &
                              statevector%allUVkEnd,1,'mpi_integer','grid',ierr)
    end if

    select case (ANLTIME_BIN)
    case ('FIRST')
       statevector%anltime=1
    case ('MIDDLE')
       statevector%anltime=nint((real(numStep,8)+1.0d0)/2.0d0)
    case ('LAST')
       statevector%anltime=numStep
    case default
      call utl_abort('gsv_allocate: unsupported value for ANLTIME_BIN = '//trim(ANLTIME_BIN))
    end select

    if (present(dateStamp_opt) .and. present(dateStampList_opt)) then
      call utl_abort('gsv_allocate: Either dateStamp or dateStampList should be presented but not both')
    else if (present(dateStampList_opt)) then
      allocate(statevector%dateStampList(numStep))
      do stepIndex = 1, numStep
        statevector%dateStampList(stepIndex)= dateStampList_opt(stepIndex)
      end do
      statevector%dateStamp3d => statevector%dateStampList(statevector%anltime)
    else if (present(dateStamp_opt)) then
      allocate(statevector%dateStampList(numStep))
      if (numStep == 1) then
        statevector%dateStampList(1) = dateStamp_opt
      else
        call tim_getstamplist(statevector%dateStampList,numStep,dateStamp_opt)
      end if
      statevector%dateStamp3d => statevector%dateStampList(statevector%anltime)
    else
      nullify(statevector%dateStamplist)
    end if
    allocate(statevector%dateOriginList(numStep))
    allocate(statevector%npasList(numStep))
    allocate(statevector%ip2List(numStep))

    call gsv_resetTimeParams(statevector)

    if (statevector%dataKind == 8) then
      allocate(statevector%gd_r8(statevector%myLonBeg:statevector%myLonEnd,  &
                                 statevector%myLatBeg:statevector%myLatEnd,  &
                                 statevector%mykBeg:statevector%mykEnd,numStep),stat=ierr)
      if (statevector%UVComponentPresent) then
        allocate(statevector%gdUV(statevector%myUVkBeg:statevector%myUVkEnd))
        if (statevector%extraUVallocated) then
          do kIndex = statevector%myUVkBeg, statevector%myUVkEnd
            allocate(statevector%gdUV(kIndex)%r8(statevector%myLonBeg:statevector%myLonEnd,  &
                                                 statevector%myLatBeg:statevector%myLatEnd,  &
                                                 numStep))
            statevector%gdUV(kIndex)%r8(:,:,:) = 0.0d0
          end do
        else
          ! in this case, both components available on each mpi task, so just point to it
          do kIndex = statevector%myUVkBeg, statevector%myUVkEnd
            levUV = gsv_getLevFromK(statevector, kIndex)
            UVname = complementaryUVname(gsv_getVarNameFromK(statevector,kIndex))
            kIndex2 = levUV + gsv_getOffsetFromVarName(statevector,UVname)
            lon1 = statevector%myLonBeg
            lat1 = statevector%myLatBeg
            statevector%gdUV(kIndex)%r8(lon1:,lat1:,1:) => statevector%gd_r8(:,:,kIndex2,:)
          end do
        end if
      end if
    else if (statevector%dataKind == 4) then
      allocate(statevector%gd_r4(statevector%myLonBeg:statevector%myLonEnd,  &
                                 statevector%myLatBeg:statevector%myLatEnd,  &
                                 statevector%mykBeg:statevector%mykEnd,numStep),stat=ierr)
      if (statevector%UVComponentPresent) then
        allocate(statevector%gdUV(statevector%myUVkBeg:statevector%myUVkEnd))
        if (statevector%extraUVallocated) then
          do kIndex = statevector%myUVkBeg, statevector%myUVkEnd
            allocate(statevector%gdUV(kIndex)%r4(statevector%myLonBeg:statevector%myLonEnd,  &
                                                    statevector%myLatBeg:statevector%myLatEnd,  &
                                                    numStep))
            statevector%gdUV(kIndex)%r4(:,:,:) = 0.0
          end do
        else
          ! in this case, both components available on each mpi task, so just point to it
          do kIndex = statevector%myUVkBeg, statevector%myUVkEnd
            levUV = gsv_getLevFromK(statevector, kIndex)
            UVname = complementaryUVname(gsv_getVarNameFromK(statevector,kIndex))
            kIndex2 = levUV + gsv_getOffsetFromVarName(statevector,UVname)
            lon1 = statevector%myLonBeg
            lat1 = statevector%myLatBeg
            statevector%gdUV(kIndex)%r4(lon1:,lat1:,1:) => statevector%gd_r4(:,:,kIndex2,:)
          end do
        end if
      end if
    else
      call utl_abort('gsv_allocate: unknown value of datakind')
    end if
    if (ierr.ne.0) then
      call utl_abort('gsv_allocate: Problem allocating memory! id=1 '//str(ierr))
    end if

    if (present(allocHeightSfc_opt)) then
      if (allocHeightSfc_opt) then
        ! if VarsLevs, then only proc 0 allocates surface height, otherwise all procs do
        statevector%heightSfcPresent = .true.
        if ((statevector%mpi_distribution == 'VarsLevs' .and. mmpi_myid == 0) .or. &
             statevector%mpi_distribution /= 'VarsLevs') then
          allocate(statevector%HeightSfc(statevector%myLonBeg:statevector%myLonEnd,  &
                                     statevector%myLatBeg:statevector%myLatEnd))
          statevector%HeightSfc(:,:) = 0.0d0
        end if
      end if
    end if

    lon1=statevector%myLonBeg
    lat1=statevector%myLatBeg
    k1=statevector%mykBeg
    if (statevector%dataKind == 8) then
      statevector%gd3d_r8(lon1:,lat1:,k1:) => statevector%gd_r8(:,:,:,statevector%anltime)
    else if (statevector%dataKind == 4) then
      statevector%gd3d_r4(lon1:,lat1:,k1:) => statevector%gd_r4(:,:,:,statevector%anltime)
    end if

    statevector%addHeightSfcOffset = addHeightSfcOffset

    statevector%allocated=.true.

    call utl_tmg_stop(168)

  end subroutine gsv_allocate

  !--------------------------------------------------------------------------
  ! gsv_communicateTimeParams
  !--------------------------------------------------------------------------
  subroutine gsv_communicateTimeParams(statevector)
    !
    ! :Purpose: Ensures all mpi tasks have certain time and other parameters
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: statevector

    ! Locals:
    integer :: deet, ierr
    integer :: ip2List(statevector%numStep), npasList(statevector%numStep)
    integer :: dateOriginList(statevector%numStep)
    logical :: onPhysicsGrid(vnl_numVarMax)

    call rpn_comm_allreduce(statevector%deet, deet, 1,  &
                            'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)
    statevector%deet = deet
    call rpn_comm_allreduce(statevector%ip2List, ip2List, statevector%numStep,  &
                            'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)
    statevector%ip2List(:) = ip2List(:)
    call rpn_comm_allreduce(statevector%npasList, npasList, statevector%numStep,  &
                            'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)
    statevector%npasList(:) = npasList(:)
    call rpn_comm_allreduce(statevector%dateOriginList, dateOriginList, statevector%numStep,  &
                            'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)
    statevector%dateOriginList(:) = dateOriginList(:)

    call rpn_comm_allreduce(statevector%onPhysicsGrid(:), onPhysicsGrid(:), size(onPhysicsGrid),  &
                            'MPI_LOGICAL', 'MPI_LOR', 'GRID', ierr)
    statevector%onPhysicsGrid(:) = onPhysicsGrid(:)

    call msg('gsv_communicateTimeParams', 'deet = '//str(deet) &
         //new_line('')//'ip2List = '//str(ip2List(:)) &
         //new_line('')//'npasList = '//str(npasList(:)) &
         //new_line('')//'dateOriginList = '//str(dateOriginList(:)))

  end subroutine gsv_communicateTimeParams

  !--------------------------------------------------------------------------
  ! gsv_resetTimeParams
  !--------------------------------------------------------------------------
  subroutine gsv_resetTimeParams(statevector)
    !
    ! :Purpose: Resets certain time parameters to "missing" values
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: statevector

    statevector%dateOriginList(:) =  mpc_missingValue_int
    statevector%npasList(:)       =  mpc_missingValue_int
    statevector%ip2List(:)        =  mpc_missingValue_int
    statevector%deet              =  mpc_missingValue_int
    statevector%etiket            =  "UNDEFINED"

  end subroutine gsv_resetTimeParams

  !--------------------------------------------------------------------------
  ! gsv_checkMpiDistribution
  !--------------------------------------------------------------------------
  subroutine gsv_checkMpiDistribution(stateVector)
    !
    ! :Purpose: Checks the distribution of latitude and longitude gridpoints
    !           over the mpi tasks. If the variation in the number of grid
    !           points in either direction is too large, other mpi topologies
    !           will be suggested in the listing and the program could
    !           potentially abort. The printing to the listing is limited
    !           to only the first 5 calls.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in) :: statevector

    ! Locals:
    integer       :: npex, npey
    integer       :: lonPerPEmin, lonPerPEmax, latPerPEmin, latPerPEmax
    integer, save :: numCalls = 0

    ! check if distribution of gridpoints over mpi tasks is very uneven
    if (maxval(statevector%allLonPerPE) > 2*minval(statevector%allLonPerPE) .or. &
        maxval(statevector%allLatPerPE) > 2*minval(statevector%allLatPerPE)) then
      numCalls = numCalls + 1
      if (mmpi_myid == 0 .and. (numCalls <= 5)) then
        call msg('gsv_checkMpiDistribution', & 
            new_line('')//'=============================================================' &
          //new_line('')//'WARNING: bad choice of mpi topology!' &
          //new_line('')//'   mpi x, y dimensions = '//str(mmpi_npex)//', '//str(mmpi_npey) &
          //new_line('')//'   min(lonPerPE) = '//str(minval(statevector%allLonPerPE)) &
          //new_line('')//'   max(lonPerPE) = '//str(maxval(statevector%allLonPerPE)) &
          //new_line('')//'   min(latPerPE) = '//str(minval(statevector%allLatPerPE)) &
          //new_line('')//'   max(latPerPE) = '//str(maxval(statevector%allLatPerPE)) &
          //new_line(''), mpiAll_opt=.false.)

        ! make suggestions for mpi x diminension
        if (maxval(statevector%allLonPerPE) > 2*minval(statevector%allLonPerPE)) then
          call msg('gsv_checkMpiDistribution', & 
                   'Please choose a value of mpi x dimension that gives a smaller ' &
                   //'difference between min and max of lonPerPE. Here are some options:', &
                   mpiAll_opt=.false.)
          do npex = 1, 2*mmpi_npex
            lonPerPEmin = floor(real(stateVector%ni)/real(npex))
            lonPerPEmax = stateVector%ni - (npex - 1) * lonPerPEmin
            if (lonPerPEmax < 2*lonPerPEmin) then
              call msg('gsv_checkMpiDistribution','mpi x dimension = '//str(npex) &
                   //', difference between min and max lonPerPE = ' &
                   //str(lonPerPEmax - lonPerPEmin), mpiAll_opt=.false.)
            end if
          end do
        end if

        ! make suggestions for mpi y dimension
        if (maxval(statevector%allLatPerPE) > 2*minval(statevector%allLatPerPE)) then
          call msg('gsv_checkMpiDistribution',&
               'Please choose a value of mpi y dimension that gives a smaller ' &
               //'difference between min and max of latPerPE. Here are some options:',&
               mpiAll_opt=.false.)
          do npey = 1, 2*mmpi_npey
            latPerPEmin = floor(real(stateVector%nj)/real(npey))
            latPerPEmax = stateVector%nj - (npey - 1) * latPerPEmin
            if (latPerPEmax < 2*latPerPEmin) then
              call msg('gsv_checkMpiDistribution','mpi y dimension = '//str(npey)  &
                //', difference between min and max latPerPE = ' &
                //str(latPerPEmax - latPerPEmin), mpiAll_opt=.false.)
            end if
          end do
        end if

        call msg('gsv_checkMpiDistribution', new_line('')//'=============================================================', mpiAll_opt=.false.)
      else
        ! After 5 calls, just give a short message
        call msg('gsv_checkMpiDistribution','WARNING: bad choice of mpi topology!', mpiAll_opt=.false.)
      end if

      if (abortOnMpiImbalance) call utl_abort('gsv_checkMpiDistribution: Please choose a better mpi topology')
    end if

  end subroutine gsv_checkMpiDistribution
    
  !--------------------------------------------------------------------------
  ! gsv_complementaryUVname
  !--------------------------------------------------------------------------
  function complementaryUVname(UV_in) result(UV_out)
    !
    ! :Purpose: Returns the other wind component name
    !           UU -> VV, VV -> UU
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in)  :: UV_in
    character(len=4)              :: UV_out

    ! return UU on VV input, and vice versa
    if (trim(UV_in) == 'UU') then
      UV_out = 'VV  '
    else if (trim(UV_in) == 'VV') then
      UV_out = 'UU  '
    else
      call utl_abort('complementaryUVname: invalid input, UV_in = '//trim(UV_in))
    end if
  end function complementaryUVname

  !--------------------------------------------------------------------------
  ! gsv_modifyDate
  !--------------------------------------------------------------------------
  subroutine gsv_modifyDate(statevector, dateStamp, modifyDateOrigin_opt)
    !
    ! :Purpose: Modifies a statevector reference date
    !
    implicit none
  
    ! Arguments
    type(struct_gsv),  intent(inout) :: statevector
    integer,           intent(in)    :: dateStamp
    logical, optional, intent(in)    :: modifyDateOrigin_opt

    if (statevector%numStep == 1) then
      statevector%dateStampList(1) = dateStamp
      if(present(modifyDateOrigin_opt)) statevector%dateOriginList(1) = dateStamp 
    else
      call tim_getstamplist(statevector%dateStampList, statevector%numStep, dateStamp)
      if(present(modifyDateOrigin_opt)) call tim_getstamplist(statevector%dateOriginList, statevector%numStep, dateStamp)
    end if
    statevector%dateStamp3d => statevector%dateStampList(statevector%anltime)

  end subroutine gsv_modifyDate

  !--------------------------------------------------------------------------
  ! gsv_modifyVarName
  !--------------------------------------------------------------------------
  subroutine gsv_modifyVarName(statevector, oldVarName, newVarName) 
    !
    ! :Purpose: Replaces a variable with a variable of the same vertical level type
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: statevector
    character(len=*), intent(in)    :: oldVarName
    character(len=*), intent(in)    :: newVarName
    
    ! Locals:
    integer :: varIndex_oldVarName, varIndex_newVarName

    ! Test the compatibility of the modifications
    if (.not. gsv_varExist(statevector,oldVarName)) then
      call utl_abort('gsv_modifyVarName: the varName to replace does not exist '//trim(oldVarName))
    end if
    if (gsv_varExist(statevector,newVarName)) then
      call utl_abort('gsv_modifyVarName: the varName to add already exist '//trim(oldVarName))
    end if
    if (vnl_varLevelFromVarname(newVarName) /= vnl_varLevelFromVarname(oldVarName)) then
      call utl_abort('gsv_modifyVarName: the level type are different')
    end if

    ! Find  varIndex_oldVarName & varIndex_newVarName
    varIndex_oldVarName = vnl_varListIndex(oldVarName)
    varIndex_newVarName = vnl_varListIndex(newVarName)

    ! Change the ExistList
    statevector%varExistList(varIndex_oldVarName) = .false.
    statevector%varExistList(varIndex_newVarName) = .true.

    ! Change the offset
    statevector%varOffset(varIndex_newVarName) = statevector%varOffset(varIndex_oldVarName)
    statevector%varOffset(varIndex_oldVarName) = 0

    ! Change the number of levels
    statevector%varNumLev(varIndex_newVarName) = statevector%varNumLev(varIndex_oldVarName)
    statevector%varNumLev(varIndex_oldVarName) = 0

  end subroutine gsv_modifyVarName

  !--------------------------------------------------------------------------
  ! gsv_zero
  !--------------------------------------------------------------------------
  subroutine gsv_zero(statevector)
    !
    ! :Purpose: Zeros all struct_gsv arrays
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: statevector

    ! Locals:
    integer :: stepIndex,lonIndex,kIndex,latIndex,lat1,lat2,lon1,lon2,k1,k2,k1UV,k2UV

    if (.not.statevector%allocated) then
      call utl_abort('gsv_zero: gridStateVector not yet allocated')
    end if

    lon1=statevector%myLonBeg
    lon2=statevector%myLonEnd
    lat1=statevector%myLatBeg
    lat2=statevector%myLatEnd
    k1=statevector%mykBeg
    k2=statevector%mykEnd
    k1UV = statevector%myUVkBeg
    k2UV = statevector%myUVkEnd

    if (associated(statevector%HeightSfc)) statevector%HeightSfc(:,:) = 0.0d0

    if (statevector%dataKind == 8) then

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

      if (statevector%extraUVallocated) then
        !$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex)    
        do kIndex = k1UV, k2UV
          do stepIndex = 1, statevector%numStep
            do latIndex = lat1, lat2
              do lonIndex = lon1, lon2
                statevector%gdUV(kIndex)%r8(lonIndex,latIndex,stepIndex) = 0.0d0
              end do
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end if

    else if (statevector%dataKind == 4) then

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

      if (statevector%extraUVallocated) then
        !$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex)    
        do kIndex = k1UV, k2UV
          do stepIndex = 1, statevector%numStep
            do latIndex = lat1, lat2
              do lonIndex = lon1, lon2
                statevector%gdUV(kIndex)%r4(lonIndex,latIndex,stepIndex) = 0.0
              end do
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end if

    else
      call utl_abort('gsv_zero: unknown value of datakind')
    end if
    
  end subroutine gsv_zero

  !--------------------------------------------------------------------------
  ! gsv_add
  !--------------------------------------------------------------------------
  subroutine gsv_add(statevector_in,statevector_inout,scaleFactor_opt)
    !
    ! :Purpose: Adds two statevectors
    !           statevector_inout = statevector_inout + scaleFactor_opt * statevector_in
    !
    implicit none

    ! Arguments:
    type(struct_gsv),  intent(in)     :: statevector_in     ! first operand 
    type(struct_gsv),  intent(inout)  :: statevector_inout  ! second operand, will receive the result
    real(8), optional, intent(in)     :: scaleFactor_opt    ! optional scaling of the second operand prior to the addition

    ! Locals:
    integer           :: stepIndex,lonIndex,kIndex,latIndex,lon1,lon2,lat1,lat2,k1,k2

    if (.not.statevector_in%allocated) then
      call utl_abort('gsv_add: gridStateVector_in not yet allocated')
    end if
    if (.not.statevector_inout%allocated) then
      call utl_abort('gsv_add: gridStateVector_inout not yet allocated')
    end if

    lon1=statevector_in%myLonBeg
    lon2=statevector_in%myLonEnd
    lat1=statevector_in%myLatBeg
    lat2=statevector_in%myLatEnd
    k1=statevector_in%mykBeg
    k2=statevector_in%mykEnd

    if (statevector_inout%dataKind == 8 .and. statevector_in%dataKind == 8) then

      if (present(scaleFactor_opt)) then
        do stepIndex = 1, statevector_inout%numStep
          !$OMP PARALLEL DO PRIVATE (latIndex,kIndex,lonIndex)    
          do kIndex = k1, k2
            do latIndex = lat1, lat2
              do lonIndex = lon1, lon2
                statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex) = &
                     statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex) +  &
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
                statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex) = &
                     statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex) +  &
                     statevector_in%gd_r8(lonIndex,latIndex,kIndex,stepIndex)
              end do
            end do
          end do
          !$OMP END PARALLEL DO
        end do
      end if

    else if (statevector_inout%dataKind == 4 .and. statevector_in%dataKind == 4) then

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
  ! gsv_schurProduct
  !--------------------------------------------------------------------------
  subroutine gsv_schurProduct(statevector_in,statevector_inout)
    !
    ! :Purpose: Applies the Schur product of two statevector
    !           statevector_inout(i,j,k,l) = statevector_inout(i,j,k,l) * statevector_in(i,j,k,l) 
    !
    implicit none

    ! Arguments:
    type(struct_gsv),  intent(in)     :: statevector_in     ! first operand 
    type(struct_gsv),  intent(inout)  :: statevector_inout  ! second operand, will receive the result

    ! Locals:
    integer :: stepIndex,lonIndex,kIndex,latIndex,lon1,lon2,lat1,lat2,k1,k2

    if (.not.statevector_in%allocated) then
      call utl_abort('gsv_schurProduct: gridStateVector_in not yet allocated')
    end if
    if (.not.statevector_inout%allocated) then
      call utl_abort('gsv_schurProduct: gridStateVector_inout not yet allocated')
    end if

    lon1=statevector_in%myLonBeg
    lon2=statevector_in%myLonEnd
    lat1=statevector_in%myLatBeg
    lat2=statevector_in%myLatEnd
    k1=statevector_in%mykBeg
    k2=statevector_in%mykEnd

    if (statevector_inout%dataKind == 8 .and. statevector_in%dataKind == 8) then

      do stepIndex = 1, statevector_inout%numStep
        !$OMP PARALLEL DO PRIVATE (latIndex,kIndex,lonIndex)    
        do kIndex = k1, k2
          do latIndex = lat1, lat2
            do lonIndex = lon1, lon2
              statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex) = &
                   statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex) *  &
                   statevector_in%gd_r8(lonIndex,latIndex,kIndex,stepIndex)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do
      

    else if (statevector_inout%dataKind == 4 .and. statevector_in%dataKind == 4) then

      do stepIndex = 1, statevector_inout%numStep
        !$OMP PARALLEL DO PRIVATE (latIndex,kIndex,lonIndex)    
        do kIndex = k1, k2
          do latIndex = lat1, lat2
            do lonIndex = lon1, lon2
              statevector_inout%gd_r4(lonIndex,latIndex,kIndex,stepIndex) = &
                   statevector_inout%gd_r4(lonIndex,latIndex,kIndex,stepIndex) *  &
                   statevector_in%gd_r4(lonIndex,latIndex,kIndex,stepIndex)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do

    else
      call utl_abort('gsv_schurProduct: Data type must be the same for both statevectors')
    end if

  end subroutine gsv_schurProduct

  !--------------------------------------------------------------------------
  ! gsv_copy
  !--------------------------------------------------------------------------
  subroutine gsv_copy(statevector_in, statevector_out, stepIndexOut_opt, &
                      allowTimeMismatch_opt, allowVarMismatch_opt)
    !
    ! :Purpose: Copies a statevector
    !
    implicit none

    ! Arguments:
    type(struct_gsv),  intent(in)    :: statevector_in
    type(struct_gsv),  intent(inout) :: statevector_out
    integer, optional, intent(in)    :: stepIndexOut_opt
    logical, optional, intent(in)    :: allowTimeMismatch_opt
    logical, optional, intent(in)    :: allowVarMismatch_opt

    ! Locals:
    logical            :: timeMismatch, allowVarMismatch, varMismatch

    integer :: stepIndex, lonIndex, kIndex, latIndex, levIndex, varIndex, numCommonVar 
    integer :: lon1, lon2, lat1, lat2, k1, k2, step1, step2, stepIn, nlev_in

    real(4), pointer :: field_out_r4(:,:,:,:), field_in_r4(:,:,:,:)
    real(8), pointer :: field_out_r8(:,:,:,:), field_in_r8(:,:,:,:)

    character(len=4), allocatable :: varNameListCommon(:)
    character(len=4)              :: varName
    character(len=10)             :: gsvCopyType 
    character(len=4), pointer     :: varNamesList_in(:), varNamesList_out(:)

    if (present(allowVarMismatch_opt)) then
      allowVarMismatch = allowVarMismatch_opt
    else
      allowVarMismatch = .false.
    end if
    varMismatch = .false.

    timeMismatch = .false.
    if (present(allowTimeMismatch_opt)) then
      if (allowTimeMismatch_opt) then
        if (statevector_in%numStep < statevector_out%numStep) then
          call utl_abort('gsv_copy: numStep_in less than numStep_out, which is not allowed')
        end if
        if (statevector_in%numStep /= statevector_out%numStep) then
          timeMismatch = .true.
        end if
      else
        if (statevector_in%numStep /= statevector_out%numStep) then
          call utl_abort('gsv_copy: numStep_in not equal to numStep_out')
        end if
      end if
    end if

    if (present(stepIndexOut_opt) .and. present(allowTimeMismatch_opt)) then
      call utl_abort('gsv_copy: Cannot specify both stepIndexOut_opt ' //  &
                     'and allowTimeMismatch_opt in the same call')
    end if

    if (.not.statevector_in%allocated) then
      call utl_abort('gsv_copy: gridStateVector_in not yet allocated')
    end if
    if (.not.statevector_out%allocated) then
      call utl_abort('gsv_copy: gridStateVector_out not yet allocated')
    end if

    if (statevector_in%mpi_distribution == 'VarsLevs') allowVarMismatch = .false.

    nullify(varNamesList_in)
    nullify(varNamesList_out)
    call gsv_varNamesList(varNamesList_in,statevector_in)
    call gsv_varNamesList(varNamesList_out,statevector_out)

    if (size(varNamesList_in(:)) /= size(varNamesList_out(:))) then
      varMismatch = .true.
    else 
      if (all(varNamesList_in(:) == varNamesList_out(:))) then
        varMismatch = .false.
      else
        varMismatch = .true.
      end if 
    end if
    deallocate(varNamesList_out)
    deallocate(varNamesList_in)

    ! if varMismatch and allowVarMismatch -> copy by varName, else copy by kIndex
    if (varMismatch .and. allowVarMismatch) then 
      gsvCopyType = 'VarName'
    else if (.not. varMismatch) then
      gsvCopyType = 'kIndex'
    else 
      call utl_abort('gsv_copy: varMismatch and allowVarMismatch do not agree! Aborting.')
    end if

    call msg('gsv_copy', 'gsvCopyType='//gsvCopyType &
         //', timeMismatch='//str(timeMismatch) &
         //', varMismatch='//str(varMismatch) &
         //', allowVarMismatch='//str(allowVarMismatch), verb_opt=2)

    ! build list of common variables and see if there is a mismatch
    allocate(varNameListCommon(vnl_numvarmax))
    varNameListCommon(:) = '    '
    if (varMismatch) then
      numCommonVar = 0
      do varIndex = 1, vnl_numvarmax
        varName = vnl_varNameList(varIndex)
        if (gsv_varExist(statevector_in,varName) .and. gsv_varExist(statevector_out,varName)) then
          numCommonVar = numCommonVar + 1
          varNameListCommon(numCommonVar) = varName 
        end if 
      end do
    end if

    lon1 = statevector_in%myLonBeg
    lon2 = statevector_in%myLonEnd
    lat1 = statevector_in%myLatBeg
    lat2 = statevector_in%myLatEnd
    k1 = statevector_in%mykBeg
    k2 = statevector_in%mykEnd
    ! If stepIndexOut_opt present then copy from step 1 to stepIndexOut_opt
    if (present(stepIndexOut_opt)) then
      step1 = stepIndexOut_opt
      step2 = stepIndexOut_opt
    else
      step1 = 1
      step2 = statevector_out%numStep
    end if

    ! copy over some time related parameters
    statevector_out%deet   = statevector_in%deet
    statevector_out%etiket = statevector_in%etiket
    do stepIndex = step1, step2
      if (present(stepIndexOut_opt)) then
        stepIn = 1
      else if(timeMismatch) then
        stepIn_Loop0: do stepIn = 1, statevector_in%numStep
          if (statevector_in%dateStampList(stepIn) ==  &
              statevector_out%dateStampList(stepIndex)) exit stepIn_loop0
        end do stepIn_Loop0
      else
        stepIn = stepIndex
      end if
      statevector_out%dateOriginList(stepIndex) = statevector_in%dateOriginList(stepIn)
      statevector_out%npasList(stepIndex)       = statevector_in%npasList(stepIn)
      statevector_out%ip2List(stepIndex)        = statevector_in%ip2List(stepIn)
    end do

    if (associated(statevector_in%HeightSfc) .and. associated(statevector_out%HeightSfc)) then
      statevector_out%HeightSfc(:,:) = statevector_in%HeightSfc(:,:)
    end if

    if (statevector_out%dataKind == 8 .and. statevector_in%dataKind == 8) then

      if (trim(gsvCopyType) == 'kIndex') then
        !$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex,stepIn)
        do kIndex = k1, k2
          do stepIndex = step1, step2
            if (present(stepIndexOut_opt)) then
              stepIn = 1
            else if(timeMismatch) then
              stepIn_Loop: do stepIn = 1, statevector_in%numStep
                if (statevector_in%dateStampList(stepIn) ==  &
                    statevector_out%dateStampList(stepIndex)) exit stepIn_loop
              end do stepIn_Loop
            else
              stepIn = stepIndex
            end if
            do latIndex = lat1, lat2
              do lonIndex = lon1, lon2
                statevector_out%gd_r8(lonIndex,latIndex,kIndex,stepIndex) =  &
                  statevector_in%gd_r8(lonIndex,latIndex,kIndex,stepIn)
              end do
            end do
          end do
        end do
        !$OMP END PARALLEL DO

      else
        do varIndex = 1, numCommonVar
          varName = varNameListCommon(varIndex)

          nlev_in = gsv_getNumLevFromVarName(statevector_in,varName)

          call gsv_getField(statevector_in ,field_in_r8, varName)
          call gsv_getField(statevector_out,field_out_r8, varName)

          !$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,levIndex,lonIndex,stepIn)
          do stepIndex = step1, step2
            if (present(stepIndexOut_opt)) then
              stepIn = 1
            else if(timeMismatch) then
              stepIn_Loop2: do stepIn = 1, statevector_in%numStep
                if (statevector_in%dateStampList(stepIn) ==  &
                    statevector_out%dateStampList(stepIndex)) exit stepIn_loop2
              end do stepIn_Loop2
            else
              stepIn = stepIndex
            end if
            do levIndex = 1, nlev_in
              do latIndex = lat1, lat2
                do lonIndex = lon1, lon2
                  field_out_r8(lonIndex,latIndex,levIndex,stepIndex) =  &
                    field_in_r8(lonIndex,latIndex,levIndex,stepIn)
                end do
              end do
            end do
          end do
          !$OMP END PARALLEL DO

        end do
      end if

    else if (statevector_out%dataKind == 4 .and. statevector_in%dataKind == 4) then

      if (trim(gsvCopyType) == 'kIndex') then
        !$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex,stepIn)
        do kIndex = k1, k2
          do stepIndex = step1, step2
            if (present(stepIndexOut_opt)) then
              stepIn = 1
            else if(timeMismatch) then
              stepIn_Loop3: do stepIn = 1, statevector_in%numStep
                if (statevector_in%dateStampList(stepIn) ==  &
                    statevector_out%dateStampList(stepIndex)) exit stepIn_loop3
              end do stepIn_Loop3
            else
              stepIn = stepIndex
            end if
            do latIndex = lat1, lat2
              do lonIndex = lon1, lon2
                statevector_out%gd_r4(lonIndex,latIndex,kIndex,stepIndex) =  &
                  statevector_in%gd_r4(lonIndex,latIndex,kIndex,stepIn)
              end do
            end do
          end do
        end do
        !$OMP END PARALLEL DO

      else
        do varIndex = 1, numCommonVar
          varName = varNameListCommon(varIndex)

          nlev_in = gsv_getNumLevFromVarName(statevector_in,varName)

          call gsv_getField(statevector_in ,field_in_r4, varName)
          call gsv_getField(statevector_out,field_out_r4, varName)

          !$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,levIndex,lonIndex,stepIn)
          do stepIndex = step1, step2
            if (present(stepIndexOut_opt)) then
              stepIn = 1
            else if(timeMismatch) then
              stepIn_Loop4: do stepIn = 1, statevector_in%numStep
                if (statevector_in%dateStampList(stepIn) ==  &
                    statevector_out%dateStampList(stepIndex)) exit stepIn_loop4
              end do stepIn_Loop4
            else
              stepIn = stepIndex
            end if
            do levIndex = 1, nlev_in
              do latIndex = lat1, lat2
                do lonIndex = lon1, lon2
                  field_out_r4(lonIndex,latIndex,levIndex,stepIndex) =  &
                    field_in_r4(lonIndex,latIndex,levIndex,stepIn)
                end do
              end do
            end do
          end do
          !$OMP END PARALLEL DO
        end do
      end if

    else if (statevector_out%dataKind == 4 .and. statevector_in%dataKind == 8) then

      if (trim(gsvCopyType) == 'kIndex') then
        !$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex,stepIn)
        do kIndex = k1, k2
          do stepIndex = step1, step2
            if (present(stepIndexOut_opt)) then
              stepIn = 1
            else if(timeMismatch) then
              stepIn_Loop5: do stepIn = 1, statevector_in%numStep
                if (statevector_in%dateStampList(stepIn) ==  &
                    statevector_out%dateStampList(stepIndex)) exit stepIn_loop5
              end do stepIn_Loop5
            else
              stepIn = stepIndex
            end if
            do latIndex = lat1, lat2
              do lonIndex = lon1, lon2
                statevector_out%gd_r4(lonIndex,latIndex,kIndex,stepIndex) =  &
                  real(statevector_in%gd_r8(lonIndex,latIndex,kIndex,stepIn),4)
              end do
            end do
          end do
        end do
        !$OMP END PARALLEL DO

      else
        do varIndex = 1, numCommonVar
          varName = varNameListCommon(varIndex)

          nlev_in = gsv_getNumLevFromVarName(statevector_in,varName)

          call gsv_getField(statevector_in ,field_in_r8, varName)
          call gsv_getField(statevector_out,field_out_r4, varName)

          !$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,levIndex,lonIndex,stepIn)
          do stepIndex = step1, step2
            if (present(stepIndexOut_opt)) then
              stepIn = 1
            else if(timeMismatch) then
              stepIn_Loop6: do stepIn = 1, statevector_in%numStep
                if (statevector_in%dateStampList(stepIn) ==  &
                    statevector_out%dateStampList(stepIndex)) exit stepIn_loop6
              end do stepIn_Loop6
            else
              stepIn = stepIndex
            end if
            do levIndex = 1, nlev_in
              do latIndex = lat1, lat2
                do lonIndex = lon1, lon2
                  field_out_r4(lonIndex,latIndex,levIndex,stepIndex) =  &
                    real(field_in_r8(lonIndex,latIndex,levIndex,stepIn),4)
                end do
              end do
            end do
          end do
          !$OMP END PARALLEL DO

        end do
      end if

    else if (statevector_out%dataKind == 8 .and. statevector_in%dataKind == 4) then

      if (trim(gsvCopyType) == 'kIndex') then
        !$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex,stepIn)
        do kIndex = k1, k2
          do stepIndex = step1, step2
            if (present(stepIndexOut_opt)) then
              stepIn = 1
            else if(timeMismatch) then
              stepIn_Loop7: do stepIn = 1, statevector_in%numStep
                if (statevector_in%dateStampList(stepIn) ==  &
                    statevector_out%dateStampList(stepIndex)) exit stepIn_loop7
              end do stepIn_Loop7
            else
              stepIn = stepIndex
            end if
            do latIndex = lat1, lat2
              do lonIndex = lon1, lon2
                statevector_out%gd_r8(lonIndex,latIndex,kIndex,stepIndex) =  &
                  real(statevector_in%gd_r4(lonIndex,latIndex,kIndex,stepIn),8)
              end do
            end do
          end do
        end do
        !$OMP END PARALLEL DO

      else
        do varIndex = 1, numCommonVar
          varName = varNameListCommon(varIndex)

          nlev_in = gsv_getNumLevFromVarName(statevector_in,varName)

          call gsv_getField(statevector_in ,field_in_r4, varName)
          call gsv_getField(statevector_out,field_out_r8, varName)

          !$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,levIndex,lonIndex,stepIn)
          do stepIndex = step1, step2
            if (present(stepIndexOut_opt)) then
              stepIn = 1
            else if(timeMismatch) then
              stepIn_Loop8: do stepIn = 1, statevector_in%numStep
                if (statevector_in%dateStampList(stepIn) ==  &
                    statevector_out%dateStampList(stepIndex)) exit stepIn_loop8
              end do stepIn_Loop8
            else
              stepIn = stepIndex
            end if
            do levIndex = 1, nlev_in
              do latIndex = lat1, lat2
                do lonIndex = lon1, lon2
                  field_out_r8(lonIndex,latIndex,levIndex,stepIndex) =  &
                    real(field_in_r4(lonIndex,latIndex,levIndex,stepIn),8)
                end do
              end do
            end do
          end do
          !$OMP END PARALLEL DO

        end do
      end if

    else
      call utl_abort('gsv_copy: Unknown data types')
    end if

    deallocate(varNameListCommon)

    ! Copy mask if it exists
    call gsv_copyMask(statevector_in, statevector_out)

  end subroutine gsv_copy

  !--------------------------------------------------------------------------
  ! gsv_copyMask
  !--------------------------------------------------------------------------
  subroutine gsv_copyMask(statevector_in,statevector_out)
    !
    ! :Purpose: Copy ocean mask, if it exists.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)     :: statevector_in
    type(struct_gsv), intent(inout)  :: statevector_out

    ! Copy mask if it exists
    call ocm_copyMask(statevector_in%oceanMask, statevector_out%oceanMask)

  end subroutine gsv_copyMask

  !--------------------------------------------------------------------------
  ! gsv_copy4Dto3D
  !--------------------------------------------------------------------------
  subroutine gsv_copy4Dto3D(statevector_in,statevector_out)
    !
    ! :Purpose: Copies contents of a 4D statevector into a 3D statevector
    !           object by extracting the middle time step.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)     :: statevector_in
    type(struct_gsv), intent(inout)  :: statevector_out

    ! Locals:
    integer :: middleStepIndex, lonIndex, latIndex, kIndex
    integer :: lon1, lon2, lat1, lat2, k1, k2, numStepIn, numStepOut

    if (.not.statevector_in%allocated) then
      call utl_abort('gsv_copy4Dto3D: gridStateVector_in not yet allocated')
    end if
    if (.not.statevector_out%allocated) then
      call utl_abort('gsv_copy4Dto3D: gridStateVector_out not yet allocated')
    end if

    lon1 = statevector_in%myLonBeg
    lon2 = statevector_in%myLonEnd
    lat1 = statevector_in%myLatBeg
    lat2 = statevector_in%myLatEnd
    k1 = statevector_in%mykBeg
    k2 = statevector_in%mykEnd
    numStepIn  =  statevector_in%numStep
    numStepOut =  statevector_out%numStep

    if (numStepOut /= 1) call utl_abort('gsv_copy4Dto3D: output statevector must have only 1 timestep')
    if (numStepIn == 1) then
      call msg('gsv_copy4Dto3D', 'WARNING: input statevector only has 1 timestep, will simply copy.')
    end if
    middleStepIndex = (numStepIn + 1) / 2

    if (associated(statevector_in%HeightSfc) .and. associated(statevector_out%HeightSfc)) then
      statevector_out%HeightSfc(:,:) = statevector_in%HeightSfc(:,:)
    end if

    if (statevector_out%dataKind == 8 .and. statevector_in%dataKind == 8) then

      !$OMP PARALLEL DO PRIVATE (kIndex,latIndex,lonIndex)
      do kIndex = k1, k2
        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2
            statevector_out%gd_r8(lonIndex,latIndex,kIndex,1) =  &
              statevector_in%gd_r8(lonIndex,latIndex,kIndex,middleStepIndex)
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    else if (statevector_out%dataKind == 4 .and. statevector_in%dataKind == 4) then

      !$OMP PARALLEL DO PRIVATE (kIndex,latIndex,lonIndex)
      do kIndex = k1, k2
        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2
            statevector_out%gd_r4(lonIndex,latIndex,kIndex,1) =  &
              statevector_in%gd_r4(lonIndex,latIndex,kIndex,middleStepIndex)
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    else
      call utl_abort('gsv_copy4Dto3D: Data type must be the same for both statevectors')
    end if

  end subroutine gsv_copy4Dto3D

  !--------------------------------------------------------------------------
  ! gsv_copyHeightSfc
  !--------------------------------------------------------------------------
  subroutine gsv_copyHeightSfc(statevector_in,statevector_out)
    !
    ! :Purpose: Copies HeightSfc data from one statevector to another
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)    :: statevector_in
    type(struct_gsv), intent(inout) :: statevector_out

    if (.not.statevector_in%allocated) then
      call utl_abort('gsv_copyHeightSfc: gridStateVector_in not yet allocated')
    end if
    if (.not.statevector_out%allocated) then
      call utl_abort('gsv_copyHeightSfc: gridStateVector_out not yet allocated')
    end if

    if (.not. associated(statevector_in%HeightSfc)) then
      call utl_abort('gsv_copyHeightSfc: HeightSfc in gridStateVector_in not allocated')
    end if
    if (.not. associated(statevector_out%HeightSfc)) then
      call utl_abort('gsv_copyHeightSfc: HeightSfc in gridStateVector_out not allocated')
    end if

    statevector_out%HeightSfc(:,:) = statevector_in%HeightSfc(:,:)

  end subroutine gsv_copyHeightSfc

  !--------------------------------------------------------------------------
  ! gsv_hPad
  !--------------------------------------------------------------------------
  subroutine gsv_hPad(statevector_in,statevector_out)
    !
    ! :Purpose: Copies a statevector to a horizontally larger one and pad with a
    !           predefined value of 0 (or 1000 for 'P0'). 
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)     :: statevector_in
    type(struct_gsv), intent(inout)  :: statevector_out

    ! Locals:
    integer :: stepIndex,lonIndex,kIndex,latIndex
    integer :: lonBeg_in, lonEnd_in, latBeg_in, latEnd_in, kBeg, kEnd

    real(8) :: paddingValue_r8
    real(4) :: paddingValue_r4

    if (.not.statevector_in%allocated) then
      call utl_abort('gsv_hPad: gridStateVector_in not yet allocated')
    end if
    if (.not.statevector_out%allocated) then
      call utl_abort('gsv_hPad: gridStateVector_out not yet allocated')
    end if
    if (statevector_in%mpi_local .or. statevector_out%mpi_local) then
       call utl_abort('gsv_hPad: both gridStateVectors must be NO MPI')
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
        latEnd_in > statevector_out%myLatEnd) then
      call utl_abort('gsv_hPad: StateVector_out is SMALLER than StateVector_in')
    end if
    if (kBeg /= statevector_out%mykBeg .or. kEnd /= statevector_out%mykEnd) then
      call utl_abort('gsv_hPad: Vertical levels are not compatible')
    end if

    ! copy over some time related parameters
    statevector_out%deet   = statevector_in%deet
    statevector_out%etiket = statevector_in%etiket
    do stepIndex = 1, statevector_out%numStep
      statevector_out%dateOriginList(stepIndex) = statevector_in%dateOriginList(stepIndex)
      statevector_out%npasList(stepIndex)       = statevector_in%npasList(stepIndex)
      statevector_out%ip2List(stepIndex)        = statevector_in%ip2List(stepIndex)
    end do

    if (statevector_out%dataKind == 8 .and. statevector_in%dataKind == 8) then

      if (associated(statevector_in%HeightSfc) .and. associated(statevector_out%HeightSfc)) then
        statevector_out%HeightSfc(:,:) = 0.d0
        statevector_out%HeightSfc(lonBeg_in:lonEnd_in,latBeg_in:latEnd_in) = &
             statevector_in%HeightSfc(:,:)
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
              statevector_out%gd_r8(lonIndex,latIndex,kIndex,stepIndex) = &
                   statevector_in%gd_r8(lonIndex,latIndex,kIndex,stepIndex)
            end do
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    else if (statevector_out%dataKind == 4 .and. statevector_in%dataKind == 4) then

      if (associated(statevector_in%HeightSfc) .and. associated(statevector_out%HeightSfc)) then
        statevector_out%HeightSfc(:,:) = 0.0
        statevector_out%HeightSfc(lonBeg_in:lonEnd_in,latBeg_in:latEnd_in) = statevector_in%HeightSfc(:,:)
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
  ! gsv_power
  !--------------------------------------------------------------------------
  subroutine gsv_power(statevector_inout,power,scaleFactor_opt)
    !
    ! :Purpose: Applies the power function
    !           statevector_inout(i,j,k,l) = scaleFactor_opt * (statevector_inout(i,j,k,l)**power)
    !
    implicit none

    ! Arguments:
    type(struct_gsv),  intent(inout)  :: statevector_inout
    real(8),           intent(in)     :: power
    real(8), optional, intent(in)     :: scaleFactor_opt    ! optional scaling applied on the power of the operand

    ! Locals:
    integer :: stepIndex,lonIndex,kIndex,latIndex,lon1,lon2,lat1,lat2,k1,k2

    if (.not.statevector_inout%allocated) then
      call utl_abort('gsv_power: gridStateVector_inout not yet allocated')
    end if

    lon1=statevector_inout%myLonBeg
    lon2=statevector_inout%myLonEnd
    lat1=statevector_inout%myLatBeg
    lat2=statevector_inout%myLatEnd
    k1=statevector_inout%mykBeg
    k2=statevector_inout%mykEnd

    if (statevector_inout%dataKind == 8) then

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

    else if (statevector_inout%dataKind == 4) then

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
  ! gsv_scale
  !--------------------------------------------------------------------------
  subroutine gsv_scale(statevector_inout,scaleFactor)
    !
    ! :Purpose: Applies scaling factor to a statevector
    !           statevector_inout = scaleFactor * statevector_inout
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: statevector_inout
    real(8),          intent(in)    :: scaleFactor

    ! Locals:
    integer :: stepIndex,lonIndex,kIndex,latIndex,lon1,lon2,lat1,lat2,k1,k2

    if (.not.statevector_inout%allocated) then
      call utl_abort('gsv_scale: gridStateVector_inout not yet allocated')
    end if

    lon1=statevector_inout%myLonBeg
    lon2=statevector_inout%myLonEnd
    lat1=statevector_inout%myLatBeg
    lat2=statevector_inout%myLatEnd
    k1=statevector_inout%mykBeg
    k2=statevector_inout%mykEnd

    if (statevector_inout%dataKind == 8) then

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
  ! gsv_scaleVertical
  !--------------------------------------------------------------------------
  subroutine gsv_scaleVertical(statevector_inout,scaleFactor)
    !
    ! :Purpose: Applies a specific scaling to each level
    !           statevector_inout(:,:,k,:) = scaleFactor(k) * statevector_inout(:,:,k,:)
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: statevector_inout
    real(8),          intent(in)    :: scaleFactor(:)

    ! Locals:
    integer          :: stepIndex,lonIndex,kIndex,latIndex,lon1,lon2,lat1,lat2,k1,k2,levIndex
    character(len=4) :: varLevel

    if (.not.statevector_inout%allocated) then
      call utl_abort('gsv_scaleVertical: gridStateVector_inout not yet allocated')
    end if

    lon1=statevector_inout%myLonBeg
    lon2=statevector_inout%myLonEnd
    lat1=statevector_inout%myLatBeg
    lat2=statevector_inout%myLatEnd
    k1=statevector_inout%mykBeg
    k2=statevector_inout%mykEnd

    if (statevector_inout%dataKind == 8) then

      !$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex,levIndex,varLevel)
      do kIndex = k1, k2
        varLevel = vnl_varLevelFromVarname(gsv_getVarNameFromK(statevector_inout, kIndex))
        if (varLevel == 'SF' .or. varLevel == 'SFMM' .or. varLevel == 'SFTH') then
          ! use lowest momentum level for surface variables
          levIndex = gsv_getNumLev(statevector_inout, 'MM')
        else
          levIndex = gsv_getLevFromK(statevector_inout, kIndex)
        end if
        do stepIndex = 1, statevector_inout%numStep
          do latIndex = lat1, lat2
            do lonIndex = lon1, lon2
              statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex) = &
                   scaleFactor(levIndex) * statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex)
            end do
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    else

      !$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex,levIndex,varLevel)
      do kIndex = k1, k2
        varLevel = vnl_varLevelFromVarname(gsv_getVarNameFromK(statevector_inout, kIndex))
        if (varLevel == 'SF' .or. varLevel == 'SFMM' .or. varLevel == 'SFTH') then
          ! use lowest momentum level for surface variables
          levIndex = gsv_getNumLev(statevector_inout, 'MM')
        else
          levIndex = gsv_getLevFromK(statevector_inout, kIndex)
        end if
        do stepIndex = 1, statevector_inout%numStep
          do latIndex = lat1, lat2
            do lonIndex = lon1, lon2
              statevector_inout%gd_r4(lonIndex,latIndex,kIndex,stepIndex) = &
                   real(scaleFactor(levIndex),4) * statevector_inout%gd_r4(lonIndex,latIndex,kIndex,stepIndex)
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
    !
    ! :Purpose: Copies the 3D data array to all time steps of the 4D array of the
    !           same statevector
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: statevector_inout

    ! Locals:
    integer :: stepIndex,lonIndex,kIndex,latIndex,lon1,lon2,lat1,lat2,k1,k2

    if (.not.statevector_inout%allocated) then
      call utl_abort('gsv_3dto4d: statevector not yet allocated')
    end if

    lon1=statevector_inout%myLonBeg
    lon2=statevector_inout%myLonEnd
    lat1=statevector_inout%myLatBeg
    lat2=statevector_inout%myLatEnd
    k1=statevector_inout%mykBeg
    k2=statevector_inout%mykEnd

    if (statevector_inout%numStep.eq.1) return

    if (statevector_inout%dataKind == 8) then

      !$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,kIndex,lonIndex)    
      do kIndex = k1, k2
        do stepIndex = 1, statevector_inout%numStep
          if (stepIndex.ne.statevector_inout%anltime) then
            do latIndex = lat1, lat2
              do lonIndex = lon1, lon2
                statevector_inout%gd_r8(lonIndex,latIndex,kIndex,stepIndex) = &
                     statevector_inout%gd3d_r8(lonIndex,latIndex,kIndex)
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
                statevector_inout%gd_r4(lonIndex,latIndex,kIndex,stepIndex) = &
                     statevector_inout%gd3d_r4(lonIndex,latIndex,kIndex)
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
    !
    ! :Purpose: Adjoint code of the 3dto4d copy to all time steps 
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: statevector_inout

    ! Locals:
    integer :: stepIndex,lonIndex,kIndex,latIndex,lon1,lon2,lat1,lat2,k1,k2

    real(4), allocatable :: gd2d_tmp_r4(:,:)
    real(8), allocatable :: gd2d_tmp(:,:)

    if (.not.statevector_inout%allocated) then
      call utl_abort('gsv_3dto4dAdj: statevector not yet allocated')
    end if

    lon1=statevector_inout%myLonBeg
    lon2=statevector_inout%myLonEnd
    lat1=statevector_inout%myLatBeg
    lat2=statevector_inout%myLatEnd
    k1=statevector_inout%mykBeg
    k2=statevector_inout%mykEnd

    if (statevector_inout%numStep.eq.1) return

    if (statevector_inout%dataKind == 8) then

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
            statevector_inout%gd3d_r8(lonIndex,latIndex,kIndex) = &
                 gd2d_tmp(lonIndex,latIndex)
          end do
        end do
      end do
      !$OMP END PARALLEL DO
      deallocate(gd2d_tmp)

    else

      allocate(gd2d_tmp_r4(lon1:lon2,lat1:lat2))
      !$OMP PARALLEL DO PRIVATE (latIndex,kIndex,lonIndex,stepIndex,gd2d_tmp_r4)
      do kIndex = k1, k2
        gd2d_tmp_r4(:,:) = 0.0
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
            statevector_inout%gd3d_r4(lonIndex,latIndex,kIndex) = &
                 gd2d_tmp_r4(lonIndex,latIndex)
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
    !
    ! :Purpose: Deallocates the struct_gsv memory structure
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: statevector

    ! Locals:
    integer :: ierr, kIndex

    if (.not.statevector%allocated) then
      call utl_abort('gsv_deallocate: gridStateVector not yet allocated')
    end if

    statevector%allocated=.false.

    if (statevector%mpi_local) then
      deallocate(statevector%allLonBeg)
      deallocate(statevector%allLonEnd)
      deallocate(statevector%allLonPerPE)
      deallocate(statevector%allLatBeg)
      deallocate(statevector%allLatEnd)
      deallocate(statevector%allLatPerPE)
      deallocate(statevector%allkBeg)
      deallocate(statevector%allkEnd)
      deallocate(statevector%allkCount)
      deallocate(statevector%allUVkBeg)
      deallocate(statevector%allUVkEnd)
      deallocate(statevector%allUVkCount)
    end if

    if (statevector%dataKind == 8) then
      deallocate(statevector%gd_r8,stat=ierr)
      nullify(statevector%gd_r8)
      if (statevector%UVComponentPresent) then 
        do kIndex = statevector%myUVkBeg, statevector%myUVkEnd
          if (statevector%extraUVallocated) deallocate(statevector%gdUV(kIndex)%r8)
          nullify(statevector%gdUV(kIndex)%r8)
        end do
        deallocate(statevector%gdUV)
        nullify(statevector%gdUV)
      end if
    else if (statevector%dataKind == 4) then
      deallocate(statevector%gd_r4,stat=ierr)
      nullify(statevector%gd_r4)
      if (statevector%UVComponentPresent) then 
        do kIndex = statevector%myUVkBeg, statevector%myUVkEnd
          if (statevector%extraUVallocated) deallocate(statevector%gdUV(kIndex)%r4)
          nullify(statevector%gdUV(kIndex)%r4)
        end do
        deallocate(statevector%gdUV)
        nullify(statevector%gdUV)
      end if
    end if
    if (ierr.ne.0) then
      write(*,*) 'gsv_deallocate: Problem detected. IERR =',ierr
    end if

    statevector%heightSfcPresent = .false.
    if (associated(statevector%HeightSfc)) then
      deallocate(statevector%HeightSfc)
      nullify(statevector%HeightSfc)
    end if

    if (associated(statevector%dateStampList)) then
      deallocate(statevector%dateStampList)
      nullify(statevector%dateStampList)
    end if
    deallocate(statevector%dateOriginList)
    deallocate(statevector%npasList)
    deallocate(statevector%ip2List)
    deallocate(statevector%varOffset)
    deallocate(statevector%varNumLev)

    call ocm_deallocate(statevector%oceanMask)

  end subroutine gsv_deallocate

  !--------------------------------------------------------------------------
  ! gsv_getField main routine and wrappers for r4 and r8
  !--------------------------------------------------------------------------
  subroutine gsv_getFieldWrapper_r4(statevector,field_r4,varName_opt)
    !
    ! :Purpose: Returns a pointer to the 4D data array. 
    !           Wrapper for the kind 4 real.
    !
    implicit none

    ! Arguments:
    type(struct_gsv),           intent(in)    :: statevector
    real(4), pointer,           intent(inout) :: field_r4(:,:,:,:)
    character(len=*), optional, intent(in)    :: varName_opt

    if (statevector%dataKind /= 4) call utl_abort('gsv_getFieldWrapper_r4: wrong dataKind')
    call gsv_getField_r48(statevector,field_r4=field_r4,varName_opt=varName_opt)

  end subroutine gsv_getFieldWrapper_r4

  subroutine gsv_getFieldWrapper_r8(statevector,field_r8,varName_opt)
    !
    ! :Purpose: Returns a pointer to the 4D data array. 
    !           Wrapper for the kind 8 real.
    !
    implicit none

    ! Arguments:
    type(struct_gsv),           intent(in)    :: statevector
    real(8), pointer,           intent(inout) :: field_r8(:,:,:,:)
    character(len=*), optional, intent(in)    :: varName_opt

    if (statevector%dataKind /= 8) call utl_abort('gsv_getFieldWrapper_r8: wrong dataKind')
    call gsv_getField_r48(statevector,field_r8=field_r8,varName_opt=varName_opt)

  end subroutine gsv_getFieldWrapper_r8

  subroutine gsv_getField_r48(statevector,field_r4,field_r8,varName_opt)
    !
    ! :Purpose: Returns a pointer to the 4D data array. 
    !           Wrapper pairing the proper real data kind
    !
    implicit none

    ! Arguments:
    type(struct_gsv),           intent(in)    :: statevector
    character(len=*), optional, intent(in)    :: varName_opt
    real(4), pointer, optional, intent(inout) :: field_r4(:,:,:,:)
    real(8), pointer, optional, intent(inout) :: field_r8(:,:,:,:)

    ! Locals:
    integer                                :: ilev1, ilev2, lon1, lat1, k1

    lon1 = statevector%myLonBeg
    lat1 = statevector%myLatBeg
    k1 = statevector%mykBeg

    if (present(varName_opt)) then
      if (statevector%mpi_distribution == 'VarsLevs') then
        call utl_abort('gsv_getField_r48: cannot specify a varName for VarsLevs mpi distribution')
      end if
      if (gsv_varExist(statevector,varName_opt)) then
        ilev1 = 1 + statevector%varOffset(vnl_varListIndex(varName_opt))
        ilev2 = ilev1 - 1 + statevector%varNumLev(vnl_varListIndex(varName_opt))
        if (gsv_getDataKind(statevector) == 4) then
          field_r4(lon1:,lat1:,1:,1:) => statevector%gd_r4(:,:,ilev1:ilev2,:)
        else
          field_r8(lon1:,lat1:,1:,1:) => statevector%gd_r8(:,:,ilev1:ilev2,:)
        end if
      else
        call utl_abort('gsv_getField_r48: Unknown variable name! ' // varName_opt)
      end if
    else
      if (gsv_getDataKind(statevector) == 4) then
        field_r4(lon1:,lat1:,k1:,1:) => statevector%gd_r4(:,:,:,:)
      else
        field_r8(lon1:,lat1:,k1:,1:) => statevector%gd_r8(:,:,:,:)
      end if
    end if

  end subroutine gsv_getField_r48

  !--------------------------------------------------------------------------
  ! gsv_getField3D_r8
  !--------------------------------------------------------------------------
  subroutine gsv_getField3D_r8(statevector,field3D,varName_opt,stepIndex_opt)
    !
    ! :Purpose: Returns a pointer to the 3D data array. 
    !           Wrapper for the kind 8 real.
    !
    implicit none

    ! Arguments:
    type(struct_gsv),           intent(in)    :: statevector
    real(8), pointer,           intent(inout) :: field3D(:,:,:)
    character(len=*), optional, intent(in)    :: varName_opt
    integer,          optional, intent(in)    :: stepIndex_opt

    ! Locals:
    integer                                :: ilev1,ilev2,lon1,lat1,k1

    lon1=statevector%myLonBeg
    lat1=statevector%myLatBeg
    k1=statevector%mykBeg

    if (.not. associated(statevector%gd3d_r8)) then
      call utl_abort('gsv_getField3D_r8: data with type r8 not allocated')
    end if

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

  end subroutine gsv_getField3D_r8

  !--------------------------------------------------------------------------
  ! gsv_getField3D_r4
  !--------------------------------------------------------------------------
  subroutine gsv_getField3D_r4(statevector,field3d,varName_opt,stepIndex_opt)
    !
    ! :Purpose: Returns a pointer to the 3D data array. 
    !           Wrapper for the kind 4 real.
    !
    implicit none

    ! Arguments:
    type(struct_gsv),           intent(in)    :: statevector
    real(4), pointer,           intent(inout) :: field3d(:,:,:)
    character(len=*), optional, intent(in)    :: varName_opt
    integer,          optional, intent(in)    :: stepIndex_opt

    ! Locals:
    integer                                :: ilev1,ilev2,lon1,lat1,k1

    lon1=statevector%myLonBeg
    lat1=statevector%myLatBeg
    k1=statevector%mykBeg

    if (.not. associated(statevector%gd3d_r4)) then
      call utl_abort('gsv_getField3D_r4: data with type r4 not allocated')
    end if

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

  end subroutine gsv_getField3D_r4

  !--------------------------------------------------------------------------
  ! gsv_getFieldUV main routine and wrappers for r4 and r8
  !--------------------------------------------------------------------------
  subroutine gsv_getFieldUVWrapper_r4(statevector,field_r4,kIndex)
    !
    ! :Purpose: Returns a pointer to the UV data array. 
    !           Wrapper for the kind 4 real.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)    :: statevector
    real(4), pointer, intent(inout) :: field_r4(:,:,:)
    integer,          intent(in)    :: kIndex

    call gsv_getFieldUV_r48(statevector,field_r4=field_r4,kIndex=kIndex)
    
  end subroutine gsv_getFieldUVWrapper_r4
  
  subroutine gsv_getFieldUVWrapper_r8(statevector,field_r8,kIndex)
    !
    ! :Purpose: Returns a pointer to the UV data array. 
    !           Wrapper for the kind 8 real.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)    :: statevector
    real(8), pointer, intent(inout) :: field_r8(:,:,:)
    integer,          intent(in)    :: kIndex

    call gsv_getFieldUV_r48(statevector,field_r8=field_r8,kIndex=kIndex)
    
  end subroutine gsv_getFieldUVWrapper_r8
  
  subroutine gsv_getFieldUV_r48(statevector,field_r4,field_r8,kIndex)
    !
    ! :Purpose: Returns a pointer to the UV data array. 
    !           Wrapper pairing the proper real data kind
    !
    implicit none

    ! Arguments:
    type(struct_gsv),           intent(in)    :: statevector
    integer,                    intent(in)    :: kIndex
    real(4), optional, pointer, intent(inout) :: field_r4(:,:,:)
    real(8), optional, pointer, intent(inout) :: field_r8(:,:,:)

    ! Locals:
    integer                                :: lon1,lat1

    lon1 = statevector%myLonBeg
    lat1 = statevector%myLatBeg

    if (gsv_getDataKind(statevector) == 4) then
      if (.not. associated(statevector%gdUV(kIndex)%r4)) call utl_abort('gsv_getFieldUV_r48: data with type r4 not allocated')
      field_r4(lon1:,lat1:,1:) => statevector%gdUV(kIndex)%r4(:,:,:)
    else
      if (.not. associated(statevector%gdUV(kIndex)%r8)) call utl_abort('gsv_getFieldUV_r48: data with type r8 not allocated')
      field_r8(lon1:,lat1:,1:) => statevector%gdUV(kIndex)%r8(:,:,:)
    end if

  end subroutine gsv_getFieldUV_r48

  !--------------------------------------------------------------------------
  ! gsv_isAssocHeightSfc
  !--------------------------------------------------------------------------
  function gsv_isAssocHeightSfc(statevector) result(isAssociated)
    !
    ! :Purpose: Returns .true. if HeightSfc is associated
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(in)        :: statevector
    logical                             :: isAssociated

    isAssociated = .false.
    if (associated(statevector%HeightSfc)) then
      isAssociated = .true.
    end if

  end function gsv_isAssocHeightSfc

  !--------------------------------------------------------------------------
  ! gsv_getHeightSfc
  !--------------------------------------------------------------------------
  function gsv_getHeightSfc(statevector) result(field)
    !
    ! :Purpose: Returns an access pointer to HeightSfc
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(in)           :: statevector
    real(8),pointer                        :: field(:,:)

    ! Locals:
    integer                                :: lon1,lat1

    lon1 = statevector%myLonBeg
    lat1 = statevector%myLatBeg

    if (.not. statevector%heightSfcPresent) call utl_abort('gsv_getHeightSfc: data not allocated')

    if (associated(statevector%HeightSfc)) then
      field(lon1:,lat1:) => statevector%HeightSfc(:,:)
    else
      nullify(field)
    end if

  end function gsv_getHeightSfc

  !--------------------------------------------------------------------------
  ! gsv_getDateStamp
  !--------------------------------------------------------------------------
  function gsv_getDateStamp(statevector,stepIndex_opt) result(dateStamp)
    !
    ! :Purpose: Returns the reference datestamp (or the datestamp of a specified
    !           time step.
    !
    implicit none

    ! Arguments
    type(struct_gsv),  intent(in)  :: statevector
    integer, optional, intent(in)  :: stepIndex_opt 
    integer                        :: dateStamp

    if (associated(statevector%dateStampList)) then
      if (present(stepIndex_opt)) then
        if (stepIndex_opt.gt.0.and.stepIndex_opt.le.statevector%numStep) then
          dateStamp=statevector%dateStampList(stepIndex_opt)
        else
          call utl_abort('gsv_getDateStamp: requested step is out of range! Step=' &
               //str(stepIndex_opt)//',numStep='//str(statevector%numStep))
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
    !
    ! :Purpose: Returns an access pointer to the statevector vco 
    !           (vertical coordinate structure)
    !
    implicit none
    type(struct_gsv), intent(in) :: statevector
    type(struct_vco), pointer    :: vco_ptr

    vco_ptr => statevector%vco

  end function gsv_getVco

  !--------------------------------------------------------------------------
  ! gsv_getHco
  !--------------------------------------------------------------------------
  function gsv_getHco(statevector) result(hco_ptr)
    !
    ! :Purpose: Returns an access pointer to the statevector hco 
    !           (horizontal coordinate structure)
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in) :: statevector
    type(struct_hco), pointer    :: hco_ptr

    hco_ptr => statevector%hco

  end function gsv_getHco

  !--------------------------------------------------------------------------
  ! gsv_getHco_physics
  !--------------------------------------------------------------------------
  function gsv_getHco_physics(statevector) result(hco_ptr)
    !
    ! :Purpose: Returns an access pointer to the statevector physics hco 
    !           (horizontal coordinate structure)
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in) :: statevector
    type(struct_hco), pointer    :: hco_ptr

    hco_ptr => statevector%hco_physics

  end function gsv_getHco_physics

  !--------------------------------------------------------------------------
  ! gsv_transposeVarsLevsToTiles
  !--------------------------------------------------------------------------
  subroutine gsv_transposeVarsLevsToTiles(statevector_in, statevector_out)
    !
    ! :Purpose: Transposes the data from mpi_distribution=VarsLevs to Tiles
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)    :: statevector_in
    type(struct_gsv), intent(inout) :: statevector_out

    ! Locals:
    integer :: youridx, youridy, yourid, nsize, maxkcount, ierr
    integer :: sendrecvKind, inKind, outKind, stepIndex
    integer, allocatable :: displs(:), nsizes(:)
    real(4), pointer     :: field_in_r4_ptr(:,:,:,:), field_out_r4_ptr(:,:,:,:)
    real(8), pointer     :: field_in_r8_ptr(:,:,:,:), field_out_r8_ptr(:,:,:,:)
    real(8), allocatable :: gd_send_height(:,:,:), gd_recv_height(:,:)
    real(8), pointer     :: field_height_in_ptr(:,:), field_height_out_ptr(:,:)

    if (statevector_in%mpi_distribution /= 'VarsLevs') then
      call utl_abort('gsv_transposeVarsLevsToTiles: input statevector must have VarsLevs mpi distribution') 
    end if

    if (statevector_out%mpi_distribution /= 'Tiles') then
      call utl_abort('gsv_transposeVarsLevsToTiles: out statevector must have Tiles mpi distribution')
    end if

    call utl_tmg_start(164,'low-level--gsv_varsLevsToTiles')

    inKind = statevector_in%dataKind
    outKind = statevector_out%dataKind
    if (inKind == 4 .or. outKind == 4) then
      sendrecvKind = 4
    else
      sendrecvKind = 8
    end if

    maxkCount = maxval(statevector_in%allkCount(:))
    if (sendrecvKind == 4) then
      call utl_reAllocate(gd_send_varsLevs_r4, statevector_out%lonPerPEmax, statevector_out%latPerPEmax, &
                           maxkCount, mmpi_nprocs)
      call utl_reAllocate(gd_recv_varsLevs_r4, statevector_out%lonPerPEmax, statevector_out%latPerPEmax, &
                          maxkCount, mmpi_nprocs)
    else
      call utl_reAllocate(gd_send_varsLevs_r8, statevector_out%lonPerPEmax, statevector_out%latPerPEmax, &
                          maxkCount, mmpi_nprocs)
      call utl_reAllocate(gd_recv_varsLevs_r8, statevector_out%lonPerPEmax, statevector_out%latPerPEmax, &
                          maxkCount, mmpi_nprocs)
    end if

    do stepIndex = 1, statevector_out%numStep

      if (sendrecvKind == 4 .and. inKind == 4) then
        call gsv_getField(statevector_in,field_in_r4_ptr)
        !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
        do youridy = 0, (mmpi_npey-1)
          do youridx = 0, (mmpi_npex-1)
            yourid = youridx + youridy*mmpi_npex
            gd_send_varsLevs_r4(1:statevector_out%allLonPerPE(youridx+1),  &
                                1:statevector_out%allLatPerPE(youridy+1),  &
                                1:statevector_in%mykCount, yourid+1) =  &
              field_in_r4_ptr(statevector_out%allLonBeg(youridx+1):statevector_out%allLonEnd(youridx+1),  &
                              statevector_out%allLatBeg(youridy+1):statevector_out%allLatEnd(youridy+1),  &
                              statevector_in%mykBeg:statevector_in%mykEnd, stepIndex)
          end do
        end do
        !$OMP END PARALLEL DO
      else if (sendrecvKind == 4 .and. inKind == 8) then
        call gsv_getField(statevector_in,field_in_r8_ptr)
        !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
        do youridy = 0, (mmpi_npey-1)
          do youridx = 0, (mmpi_npex-1)
            yourid = youridx + youridy*mmpi_npex
            gd_send_varsLevs_r4(1:statevector_out%allLonPerPE(youridx+1),  &
                                1:statevector_out%allLatPerPE(youridy+1),  &
                                1:statevector_in%mykCount, yourid+1) =  &
              real(field_in_r8_ptr(statevector_out%allLonBeg(youridx+1):statevector_out%allLonEnd(youridx+1),  &
                                   statevector_out%allLatBeg(youridy+1):statevector_out%allLatEnd(youridy+1),  &
                                   statevector_in%mykBeg:statevector_in%mykEnd, stepIndex),4)
          end do
        end do
        !$OMP END PARALLEL DO
      else if (sendrecvKind == 8 .and. inKind == 4) then
        call gsv_getField(statevector_in,field_in_r4_ptr)
        !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
        do youridy = 0, (mmpi_npey-1)
          do youridx = 0, (mmpi_npex-1)
            yourid = youridx + youridy*mmpi_npex
            gd_send_varsLevs_r8(1:statevector_out%allLonPerPE(youridx+1),  &
                                1:statevector_out%allLatPerPE(youridy+1),  &
                                1:statevector_in%mykCount, yourid+1) =  &
              real(field_in_r4_ptr(statevector_out%allLonBeg(youridx+1):statevector_out%allLonEnd(youridx+1),  &
                                   statevector_out%allLatBeg(youridy+1):statevector_out%allLatEnd(youridy+1),  &
                                   statevector_in%mykBeg:statevector_in%mykEnd, stepIndex),8)
          end do
        end do
        !$OMP END PARALLEL DO
      else if (sendrecvKind == 8 .and. inKind == 8) then
        call gsv_getField(statevector_in,field_in_r8_ptr)
        !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
        do youridy = 0, (mmpi_npey-1)
          do youridx = 0, (mmpi_npex-1)
            yourid = youridx + youridy*mmpi_npex
            gd_send_varsLevs_r8(1:statevector_out%allLonPerPE(youridx+1),  &
                                1:statevector_out%allLatPerPE(youridy+1),  &
                                1:statevector_in%mykCount, yourid+1) =  &
              field_in_r8_ptr(statevector_out%allLonBeg(youridx+1):statevector_out%allLonEnd(youridx+1),  &
                              statevector_out%allLatBeg(youridy+1):statevector_out%allLatEnd(youridy+1),  &
                              statevector_in%mykBeg:statevector_in%mykEnd, stepIndex)
          end do
        end do
        !$OMP END PARALLEL DO
      end if

      nsize = statevector_out%lonPerPEmax * statevector_out%latPerPEmax * maxkCount
      if (mmpi_nprocs > 1) then
        if (sendrecvKind == 4) then
          call rpn_comm_alltoall(gd_send_varsLevs_r4, nsize, 'mpi_real4',  &
                                 gd_recv_varsLevs_r4, nsize, 'mpi_real4', 'grid', ierr)
        else
          call rpn_comm_alltoall(gd_send_varsLevs_r8, nsize, 'mpi_real8',  &
                                 gd_recv_varsLevs_r8, nsize, 'mpi_real8', 'grid', ierr)
        end if
      else
        if (sendrecvKind == 4) then
          gd_recv_varsLevs_r4(:,:,:,1) = gd_send_varsLevs_r4(:,:,:,1)
        else
          gd_recv_varsLevs_r8(:,:,:,1) = gd_send_varsLevs_r8(:,:,:,1)
        end if
      end if

      if (sendrecvKind == 4 .and. outKind == 4) then
        call gsv_getField(statevector_out,field_out_r4_ptr)
        !$OMP PARALLEL DO PRIVATE(yourid)
        do yourid = 0, (mmpi_nprocs-1)
          field_out_r4_ptr(statevector_out%myLonBeg:statevector_out%myLonEnd, &
                           statevector_out%myLatBeg:statevector_out%myLatEnd, &
                           statevector_in%allkBeg(yourid+1):statevector_in%allkEnd(yourid+1), stepIndex) =   &
            gd_recv_varsLevs_r4(1:statevector_out%lonPerPE,  &
                                1:statevector_out%latPerPE,  &
                                1:statevector_in%allkCount(yourid+1), yourid+1)
        end do
        !$OMP END PARALLEL DO
      else if (sendrecvKind == 4 .and. outKind == 8) then
        call gsv_getField(statevector_out,field_out_r8_ptr)
        !$OMP PARALLEL DO PRIVATE(yourid)
        do yourid = 0, (mmpi_nprocs-1)
          field_out_r8_ptr(statevector_out%myLonBeg:statevector_out%myLonEnd, &
                           statevector_out%myLatBeg:statevector_out%myLatEnd, &
                           statevector_in%allkBeg(yourid+1):statevector_in%allkEnd(yourid+1), stepIndex) =   &
            real(gd_recv_varsLevs_r4(1:statevector_out%lonPerPE,  &
                                     1:statevector_out%latPerPE,  &
                                     1:statevector_in%allkCount(yourid+1), yourid+1),8)
        end do
        !$OMP END PARALLEL DO
      else if (sendrecvKind == 8 .and. outKind == 4) then
        call gsv_getField(statevector_out,field_out_r4_ptr)
        !$OMP PARALLEL DO PRIVATE(yourid)
        do yourid = 0, (mmpi_nprocs-1)
          field_out_r4_ptr(statevector_out%myLonBeg:statevector_out%myLonEnd, &
                           statevector_out%myLatBeg:statevector_out%myLatEnd, &
                           statevector_in%allkBeg(yourid+1):statevector_in%allkEnd(yourid+1), stepIndex) =   &
            real(gd_recv_varsLevs_r8(1:statevector_out%lonPerPE,  &
                                     1:statevector_out%latPerPE,  &
                                     1:statevector_in%allkCount(yourid+1), yourid+1),4)
        end do
        !$OMP END PARALLEL DO
      else if (sendrecvKind == 8 .and. outKind == 8) then
        call gsv_getField(statevector_out,field_out_r8_ptr)
        !$OMP PARALLEL DO PRIVATE(yourid)
        do yourid = 0, (mmpi_nprocs-1)
          field_out_r8_ptr(statevector_out%myLonBeg:statevector_out%myLonEnd, &
                           statevector_out%myLatBeg:statevector_out%myLatEnd, &
                           statevector_in%allkBeg(yourid+1):statevector_in%allkEnd(yourid+1), stepIndex) =   &
            gd_recv_varsLevs_r8(1:statevector_out%lonPerPE,  &
                                1:statevector_out%latPerPE,  &
                                1:statevector_in%allkCount(yourid+1), yourid+1)
        end do
        !$OMP END PARALLEL DO
      end if

    end do ! stepIndex

    if (statevector_in%heightSfcPresent .and. statevector_out%heightSfcPresent) then
      allocate(gd_send_height(statevector_out%lonPerPEmax,statevector_out%latPerPEmax,mmpi_nprocs))
      allocate(gd_recv_height(statevector_out%lonPerPEmax,statevector_out%latPerPEmax))
      gd_send_height(:,:,:) = 0.0d0
      gd_recv_height(:,:) = 0.0d0
      field_height_in_ptr => gsv_getHeightSfc(statevector_in)
      field_height_out_ptr => gsv_getHeightSfc(statevector_out)

      if (mmpi_myid == 0) then
        !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
        do youridy = 0, (mmpi_npey-1)
          do youridx = 0, (mmpi_npex-1)
            yourid = youridx + youridy*mmpi_npex
            gd_send_height(1:statevector_out%allLonPerPE(youridx+1),  &
                           1:statevector_out%allLatPerPE(youridy+1), yourid+1) =  &
              field_height_in_ptr(statevector_out%allLonBeg(youridx+1):statevector_out%allLonEnd(youridx+1),  &
                                  statevector_out%allLatBeg(youridy+1):statevector_out%allLatEnd(youridy+1))
          end do
        end do
        !$OMP END PARALLEL DO
      end if

      nsize = statevector_out%lonPerPEmax * statevector_out%latPerPEmax
      allocate(displs(mmpi_nprocs))
      allocate(nsizes(mmpi_nprocs))
      do yourid = 0, (mmpi_nprocs-1)
        displs(yourid+1) = yourid*nsize
        nsizes(yourid+1) = nsize
      end do
      call rpn_comm_scatterv(gd_send_height, nsizes, displs, 'mpi_double_precision', &
                             gd_recv_height, nsize, 'mpi_double_precision', &
                             0, 'grid', ierr)

      field_height_out_ptr(statevector_out%myLonBeg:statevector_out%myLonEnd, &
                           statevector_out%myLatBeg:statevector_out%myLatEnd) =   &
        gd_recv_height(1:statevector_out%lonPerPE,  &
                       1:statevector_out%latPerPE)

      deallocate(displs)
      deallocate(nsizes)
      deallocate(gd_recv_height)
      deallocate(gd_send_height)
    end if

    ! Copy over the mask, if it exists
    call gsv_copyMask(statevector_in, statevector_out)

    ! Copy metadata
    if (associated(statevector_in%dateStampList)) statevector_out%dateStampList(:) = statevector_in%dateStampList(:)
    if (associated(statevector_in%dateOriginList)) statevector_out%dateOriginList(:) = statevector_in%dateOriginList(:)
    if (associated(statevector_in%npasList))      statevector_out%npasList(:) = statevector_in%npasList(:)
    if (associated(statevector_in%ip2List))       statevector_out%ip2List(:) = statevector_in%ip2List(:)
    statevector_out%deet = statevector_in%deet
    statevector_out%etiket = statevector_in%etiket

    call utl_tmg_stop(164)

  end subroutine gsv_transposeVarsLevsToTiles

  !--------------------------------------------------------------------------
  ! gsv_transposeTilesToVarsLevs
  !--------------------------------------------------------------------------
  subroutine gsv_transposeTilesToVarsLevs(statevector_in, statevector_out)
    !
    !:Purpose: Transposes the data from mpi_distribution=Tiles to VarsLevs
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)    :: statevector_in
    type(struct_gsv), intent(inout) :: statevector_out

    ! Locals:
    integer :: youridx, youridy, yourid, nsize, maxkcount, ierr, mpiTagUU, mpiTagVV
    integer :: sendrecvKind, inKind, outKind, kIndexUU, kIndexVV, MpiIdUU, MpiIdVV
    integer :: levUV, stepIndex, numSend, numRecv
    integer :: requestIdSend(stateVector_out%nk), requestIdRecv(stateVector_out%nk)
    integer :: mpiStatuses(mpi_status_size,stateVector_out%nk)
    real(4), pointer     :: field_in_r4_ptr(:,:,:,:), field_out_r4_ptr(:,:,:,:)
    real(8), pointer     :: field_in_r8_ptr(:,:,:,:), field_out_r8_ptr(:,:,:,:)
    real(8), pointer     :: field_height_in_ptr(:,:), field_height_out_ptr(:,:)
    real(8), allocatable :: gd_send_height(:,:), gd_recv_height(:,:,:)

    call msg('gsv_transposeTilesToVarsLevs','START', verb_opt=2)
    call msg_memUsage('gsv_transposeTilesToVarsLevs')

    if (statevector_in%mpi_distribution /= 'Tiles') then
      call utl_abort('gsv_transposeTilesToVarsLevs: input statevector must have Tiles mpi distribution') 
    end if

    if (statevector_out%mpi_distribution /= 'VarsLevs') then
      call utl_abort('gsv_transposeTilesToVarsLevs: out statevector must have VarsLevs mpi distribution')
    end if

    inKind = statevector_in%dataKind
    outKind = statevector_out%dataKind
    if (inKind == 4 .or. outKind == 4) then
      sendrecvKind = 4
    else
      sendrecvKind = 8
    end if

    if (sendrecvKind == 4) then
      call utl_tmg_start(165,'low-level--gsv_tilesToVarsLevs_r4')
    else
      call utl_tmg_start(166,'low-level--gsv_tilesToVarsLevs_r8')
    end if

    maxkCount = maxval(statevector_out%allkCount(:))
    if (sendrecvKind == 4) then
      call utl_reAllocate(gd_send_varsLevs_r4, statevector_in%lonPerPEmax, statevector_in%latPerPEmax, &
                          maxkCount, mmpi_nprocs)
      call utl_reAllocate(gd_recv_varsLevs_r4, statevector_in%lonPerPEmax, statevector_in%latPerPEmax, &
                          maxkCount, mmpi_nprocs)
    else
      call utl_reAllocate(gd_send_varsLevs_r8, statevector_in%lonPerPEmax, statevector_in%latPerPEmax, &
                          maxkCount, mmpi_nprocs)
      call utl_reAllocate(gd_recv_varsLevs_r8, statevector_in%lonPerPEmax, statevector_in%latPerPEmax, &
                          maxkCount, mmpi_nprocs)
    end if

    do stepIndex = 1, statevector_in%numStep

      if (sendrecvKind == 8 .and. inKind == 8) then
        call gsv_getField(statevector_in,field_in_r8_ptr)
        !$OMP PARALLEL DO PRIVATE(yourid)
        do yourid = 0, (mmpi_nprocs-1)
          gd_send_varsLevs_r8(1:statevector_in%lonPerPE, &
                              1:statevector_in%latPerPE, &
                              1:statevector_out%allkCount(yourid+1), yourid+1) =  &
              field_in_r8_ptr(statevector_in%myLonBeg:statevector_in%myLonEnd, &
                              statevector_in%myLatBeg:statevector_in%myLatEnd, &
                              statevector_out%allkBeg(yourid+1):statevector_out%allkEnd(yourid+1), stepIndex)
        end do
        !$OMP END PARALLEL DO
      else if (sendrecvKind == 4 .and. inKind == 4) then
        call gsv_getField(statevector_in,field_in_r4_ptr)
        !$OMP PARALLEL DO PRIVATE(yourid)
        do yourid = 0, (mmpi_nprocs-1)
          gd_send_varsLevs_r4(1:statevector_in%lonPerPE, &
                              1:statevector_in%latPerPE, &
                              1:statevector_out%allkCount(yourid+1), yourid+1) =  &
              field_in_r4_ptr(statevector_in%myLonBeg:statevector_in%myLonEnd, &
                              statevector_in%myLatBeg:statevector_in%myLatEnd, &
                              statevector_out%allkBeg(yourid+1):statevector_out%allkEnd(yourid+1), stepIndex)
        end do
        !$OMP END PARALLEL DO
      else if (sendrecvKind == 4 .and. inKind == 8) then
        call gsv_getField(statevector_in,field_in_r8_ptr)
        !$OMP PARALLEL DO PRIVATE(yourid)
        do yourid = 0, (mmpi_nprocs-1)
          gd_send_varsLevs_r4(1:statevector_in%lonPerPE, &
                              1:statevector_in%latPerPE, &
                              1:statevector_out%allkCount(yourid+1), yourid+1) =  &
              real(field_in_r8_ptr(statevector_in%myLonBeg:statevector_in%myLonEnd, &
                                   statevector_in%myLatBeg:statevector_in%myLatEnd, &
                                   statevector_out%allkBeg(yourid+1):statevector_out%allkEnd(yourid+1), stepIndex),4)
        end do
        !$OMP END PARALLEL DO
      else
        call utl_abort('gsv_transposeTilesToVarsLevs: Incompatible mix of real 4 and 8 before alltoall mpi comm')
      end if

      nsize = statevector_in%lonPerPEmax * statevector_in%latPerPEmax * maxkCount
      if (mmpi_nprocs > 1) then
        if (sendrecvKind == 4) then
          call rpn_comm_alltoall(gd_send_varsLevs_r4, nsize, 'mpi_real4',  &
                                 gd_recv_varsLevs_r4, nsize, 'mpi_real4', 'grid', ierr)
        else
          call rpn_comm_alltoall(gd_send_varsLevs_r8, nsize, 'mpi_real8',  &
                                 gd_recv_varsLevs_r8, nsize, 'mpi_real8', 'grid', ierr)
        end if
      else
        if (sendrecvKind == 4) then
          gd_recv_varsLevs_r4(:,:,:,1) = gd_send_varsLevs_r4(:,:,:,1)
        else
          gd_recv_varsLevs_r8(:,:,:,1) = gd_send_varsLevs_r8(:,:,:,1)
        end if
      end if

      if (sendrecvKind == 8 .and. outKind == 8) then
        call gsv_getField(statevector_out,field_out_r8_ptr)
        !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
        do youridy = 0, (mmpi_npey-1)
          do youridx = 0, (mmpi_npex-1)
            yourid = youridx + youridy*mmpi_npex
            field_out_r8_ptr(statevector_in%allLonBeg(youridx+1):statevector_in%allLonEnd(youridx+1),  &
                             statevector_in%allLatBeg(youridy+1):statevector_in%allLatEnd(youridy+1),  &
                             statevector_out%mykBeg:statevector_out%mykEnd, stepIndex) = &
                gd_recv_varsLevs_r8(1:statevector_in%allLonPerPE(youridx+1),  &
                                    1:statevector_in%allLatPerPE(youridy+1),  &
                                    1:statevector_out%mykCount, yourid+1)
          end do
        end do
        !$OMP END PARALLEL DO
      else if (sendrecvKind == 4 .and. outKind == 4) then
        call gsv_getField(statevector_out,field_out_r4_ptr)
        !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
        do youridy = 0, (mmpi_npey-1)
          do youridx = 0, (mmpi_npex-1)
            yourid = youridx + youridy*mmpi_npex
            field_out_r4_ptr(statevector_in%allLonBeg(youridx+1):statevector_in%allLonEnd(youridx+1),  &
                             statevector_in%allLatBeg(youridy+1):statevector_in%allLatEnd(youridy+1),  &
                             statevector_out%mykBeg:statevector_out%mykEnd, stepIndex) = &
                gd_recv_varsLevs_r4(1:statevector_in%allLonPerPE(youridx+1),  &
                                    1:statevector_in%allLatPerPE(youridy+1),  &
                                    1:statevector_out%mykCount, yourid+1)
          end do
        end do
        !$OMP END PARALLEL DO
      else if (sendrecvKind == 4 .and. outKind == 8) then
        call gsv_getField(statevector_out,field_out_r8_ptr)
        !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
        do youridy = 0, (mmpi_npey-1)
          do youridx = 0, (mmpi_npex-1)
            yourid = youridx + youridy*mmpi_npex
            field_out_r8_ptr(statevector_in%allLonBeg(youridx+1):statevector_in%allLonEnd(youridx+1),  &
                             statevector_in%allLatBeg(youridy+1):statevector_in%allLatEnd(youridy+1),  &
                             statevector_out%mykBeg:statevector_out%mykEnd, stepIndex) = &
                real(gd_recv_varsLevs_r4(1:statevector_in%allLonPerPE(youridx+1),  &
                                         1:statevector_in%allLatPerPE(youridy+1),  &
                                         1:statevector_out%mykCount, yourid+1), 8)
          end do
        end do
        !$OMP END PARALLEL DO
      else
        call utl_abort('gsv_transposeTilesToVarsLevs: Incompatible mix of real 4 and 8 after alltoall mpi comm')
      end if

      ! send copy of wind component to task that has other component
      if (gsv_varExist(stateVector_out, 'UU') .and.  &
          gsv_varExist(stateVector_out, 'VV')) then

        numSend = 0
        numRecv = 0
        LOOP_KINDEX: do kIndexUU = 1, stateVector_out%nk
          if (gsv_getVarNameFromK(stateVector_out, kIndexUU) /= 'UU') cycle LOOP_KINDEX

          ! get k index for corresponding VV component
          levUV = gsv_getLevFromK(stateVector_out, kIndexUU)
          kIndexVV = levUV + gsv_getOffsetFromVarName(statevector_out,'VV')

          ! get Mpi task id for both U and V
          MpiIdUU = gsv_getMpiIdFromK(statevector_out,kIndexUU)
          MpiIdVV = gsv_getMpiIdFromK(statevector_out,kIndexVV)

          if (MpiIdUU == MpiIdVV .and. mmpi_myid == MpiIdUU) then
            if (outKind == 8) then
              statevector_out%gdUV(kIndexUU)%r8(:, :, stepIndex) =  statevector_out%gd_r8(:, :, kIndexVV, stepIndex)
              statevector_out%gdUV(kIndexVV)%r8(:, :, stepIndex) =  statevector_out%gd_r8(:, :, kIndexUU, stepIndex)
            else
              statevector_out%gdUV(kIndexUU)%r4(:, :, stepIndex) =  statevector_out%gd_r4(:, :, kIndexVV, stepIndex)
              statevector_out%gdUV(kIndexVV)%r4(:, :, stepIndex) =  statevector_out%gd_r4(:, :, kIndexUU, stepIndex)
            end if
            cycle LOOP_KINDEX
          end if

          ! need to do mpi communication to exchange wind components
          nsize = statevector_out%ni * statevector_out%nj
          mpiTagUU = kIndexUU
          mpiTagVV = kIndexVV
          if (mmpi_myid == MpiIdUU) then ! I have UU

            numRecv = numRecv + 1
            if (outKind == 8) then
              call mpi_irecv(statevector_out%gdUV(kIndexUU)%r8(:, :, stepIndex),  &
                             nsize, mmpi_datyp_real8, MpiIdVV, mpiTagVV,  &
                             mmpi_comm_grid, requestIdRecv(numRecv), ierr)
            else
              call mpi_irecv(statevector_out%gdUV(kIndexUU)%r4(:, :, stepIndex),  &
                             nsize, mmpi_datyp_real4, MpiIdVV, mpiTagVV,  &
                             mmpi_comm_grid, requestIdRecv(numRecv), ierr)
            end if

            numSend = numSend + 1
            if (outKind == 8) then
              call mpi_isend(statevector_out%gd_r8(:, :, kIndexUU, stepIndex),  &
                             nsize, mmpi_datyp_real8, MpiIdVV, mpiTagUU,  &
                             mmpi_comm_grid, requestIdSend(numSend), ierr)
            else
              call mpi_isend(statevector_out%gd_r4(:, :, kIndexUU, stepIndex),  &
                             nsize, mmpi_datyp_real4, MpiIdVV, mpiTagUU,  &
                             mmpi_comm_grid, requestIdSend(numSend), ierr)
            end if

          else if (mmpi_myid == MpiIdVV) then  ! I have VV

            numRecv = numRecv + 1
            if (outKind == 8) then
              call mpi_irecv(statevector_out%gdUV(kIndexVV)%r8(:, :, stepIndex),  &
                             nsize, mmpi_datyp_real8, MpiIdUU, mpiTagUU,  &
                             mmpi_comm_grid, requestIdRecv(numRecv), ierr)
            else
              call mpi_irecv(statevector_out%gdUV(kIndexVV)%r4(:, :, stepIndex),  &
                             nsize, mmpi_datyp_real4, MpiIdUU, mpiTagUU,  &
                             mmpi_comm_grid, requestIdRecv(numRecv), ierr)
            end if

            numSend = numSend + 1
            if (outKind == 8) then
              call mpi_isend(statevector_out%gd_r8(:, :, kIndexVV, stepIndex),  &
                             nsize, mmpi_datyp_real8, MpiIdUU, mpiTagVV,  &
                             mmpi_comm_grid, requestIdSend(numSend), ierr)
            else
              call mpi_isend(statevector_out%gd_r4(:, :, kIndexVV, stepIndex),  &
                             nsize, mmpi_datyp_real4, MpiIdUU, mpiTagVV,  &
                             mmpi_comm_grid, requestIdSend(numSend), ierr)
            end if

          end if

        end do LOOP_KINDEX

        if (numRecv > 0) then
          call mpi_waitAll(numRecv, requestIdRecv(1:numRecv), mpiStatuses(:,1:numRecv), ierr)
        end if

        if (numSend > 0) then
          call mpi_waitAll(numSend, requestIdSend(1:numSend), mpiStatuses(:,1:numSend), ierr)
        end if

      end if ! UU and VV exist

    end do ! stepIndex

    ! gather up surface height onto task 0
    if (statevector_in%heightSfcPresent .and. statevector_out%heightSfcPresent) then
      allocate(gd_send_height(statevector_in%lonPerPEmax,statevector_in%latPerPEmax))
      allocate(gd_recv_height(statevector_in%lonPerPEmax,statevector_in%latPerPEmax,mmpi_nprocs))
      field_height_in_ptr => gsv_getHeightSfc(statevector_in)
      field_height_out_ptr => gsv_getHeightSfc(statevector_out)

      gd_send_height(:,:) = 0.0D0
      gd_send_height(1:statevector_in%lonPerPE,1:statevector_in%latPerPE) = field_height_in_ptr(:,:)

      nsize = statevector_in%lonPerPEmax * statevector_in%latPerPEmax
      call rpn_comm_gather(gd_send_height, nsize, 'mpi_double_precision',  &
                           gd_recv_height, nsize, 'mpi_double_precision', 0, 'grid', ierr)

      if (mmpi_myid == 0) then
        !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
        do youridy = 0, (mmpi_npey-1)
          do youridx = 0, (mmpi_npex-1)
            yourid = youridx + youridy*mmpi_npex
            field_height_out_ptr(statevector_in%allLonBeg(youridx+1):statevector_in%allLonEnd(youridx+1),  &
                             statevector_in%allLatBeg(youridy+1):statevector_in%allLatEnd(youridy+1)) = &
                gd_recv_height(1:statevector_in%allLonPerPE(youridx+1),  &
                           1:statevector_in%allLatPerPE(youridy+1), yourid+1)
          end do
        end do
        !$OMP END PARALLEL DO
      end if

      deallocate(gd_send_height)
      deallocate(gd_recv_height)
    end if ! heightSfcPresent

    ! Copy over the mask, if it exists
    call gsv_copyMask(statevector_in, statevector_out)

    call msg('gsv_transposeTilesToVarsLevs','END', verb_opt=2)
    call msg_memUsage('gsv_transposeTilesToVarsLevs')

    if (sendrecvKind == 4) then
      call utl_tmg_stop(165)
    else
      call utl_tmg_stop(166)
    end if

  end subroutine gsv_transposeTilesToVarsLevs

  !--------------------------------------------------------------------------
  ! gsv_transposeTilesToVarsLevsAd
  !--------------------------------------------------------------------------
  subroutine gsv_transposeTilesToVarsLevsAd(statevector_in, statevector_out)
    !
    ! :Purpose: Adjoint of Transpose the data from mpi_distribution=Tiles to
    !           VarsLevs
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)    :: statevector_in
    type(struct_gsv), intent(inout) :: statevector_out

    ! Locals:
    integer :: youridx, youridy, yourid, nsize, maxkcount, ierr, mpiIdUU, mpiIdVV
    integer :: sendrecvKind, inKind, outKind, kIndexUU, kIndexVV, kIndex
    integer :: levUV, mpiTagUU, mpiTagVV, stepIndex, numSend, numRecv
    integer :: requestIdSend(stateVector_in%nk), requestIdRecv(stateVector_in%nk)
    integer :: mpiStatuses(mpi_status_size,stateVector_in%nk)
    integer, allocatable :: displs(:), nsizes(:)
    real(4), pointer     :: field_in_r4_ptr(:,:,:,:), field_out_r4_ptr(:,:,:,:)
    real(8), pointer     :: field_in_r8_ptr(:,:,:,:), field_out_r8_ptr(:,:,:,:)
    real(8), allocatable :: gd_send_height(:,:,:), gd_recv_height(:,:)
    real(8), pointer     :: field_height_in_ptr(:,:), field_height_out_ptr(:,:)
    real(4), allocatable :: gdUV_r4(:,:,:), gd_r4(:,:,:)
    real(8), allocatable :: gdUV_r8(:,:,:), gd_r8(:,:,:)

    call utl_tmg_start(167,'low-level--gsv_tilesToVarsLevsAD')

    if (statevector_in%mpi_distribution /= 'VarsLevs') then
      call utl_abort('gsv_transposeTilesToVarsLevsAd: input statevector must have VarsLevs mpi distribution') 
    end if

    if (statevector_out%mpi_distribution /= 'Tiles') then
      call utl_abort('gsv_transposeTilesToVarsLevsAd: out statevector must have Tiles mpi distribution')
    end if

    inKind = statevector_in%dataKind
    outKind = statevector_out%dataKind
    if (inKind == 4 .and. outKind == 4) then
      sendrecvKind = 4
    else if (inKind == 8 .and. outKind == 8) then
      sendrecvKind = 8
    else
      call utl_abort('gsv_transposeTilesToVarsLevsAd: input and output must have same dataKind')
    end if

    do stepIndex = 1, statevector_out%numStep

      ! adjoint of: send copy of wind component to task that has the other component
      if (gsv_varExist(stateVector_in, 'UU') .and.  &
          gsv_varExist(stateVector_in, 'VV')) then

        if (statevector_in%UVComponentPresent) then
          if (sendrecvKind == 4) then
            allocate(gdUV_r4(statevector_in%ni,statevector_in%nj, &
                             statevector_in%myUVkBeg:statevector_in%myUVkEnd))
            allocate(gd_r4(statevector_in%ni,statevector_in%nj, &
                           statevector_in%myUVkBeg:statevector_in%myUVkEnd))
          else
            allocate(gdUV_r8(statevector_in%ni,statevector_in%nj, &
                             statevector_in%myUVkBeg:statevector_in%myUVkEnd))
            allocate(gd_r8(statevector_in%ni,statevector_in%nj, &
                           statevector_in%myUVkBeg:statevector_in%myUVkEnd))
          end if
        end if

        numSend = 0
        numRecv = 0
        LOOP_KINDEX: do kIndexUU = 1, stateVector_in%nk
          if (gsv_getVarNameFromK(stateVector_in, kIndexUU) /= 'UU') cycle LOOP_KINDEX

          ! get k index for corresponding VV component
          levUV = gsv_getLevFromK(stateVector_in, kIndexUU)
          kIndexVV = levUV + gsv_getOffsetFromVarName(statevector_in,'VV')

          ! get Mpi task id for both U and V
          MpiIdUU = gsv_getMpiIdFromK(statevector_in,kIndexUU)
          MpiIdVV = gsv_getMpiIdFromK(statevector_in,kIndexVV)

          if (MpiIdUU == MpiIdVV .and.  mmpi_myid == MpiIdUU) then
            if (sendrecvKind == 4) then
              gd_r4(:, :, kIndexUU) = statevector_in%gdUV(kIndexVV)%r4(:, :, stepIndex)
              gd_r4(:, :, kIndexVV) = statevector_in%gdUV(kIndexUU)%r4(:, :, stepIndex)
            else
              gd_r8(:, :, kIndexUU) = statevector_in%gdUV(kIndexVV)%r8(:, :, stepIndex)
              gd_r8(:, :, kIndexVV) = statevector_in%gdUV(kIndexUU)%r8(:, :, stepIndex)
            end if
            cycle LOOP_KINDEX
          end if

          ! need to do mpi communication to exchange wind components
          nsize = statevector_in%ni * statevector_in%nj
          mpiTagUU = kIndexUU
          mpiTagVV = kIndexVV

          if (mmpi_myid == MpiIdUU) then ! I have UU

            numRecv = numRecv + 1
            if (sendrecvKind == 4) then
              call mpi_irecv(gd_r4(:, :, kIndexUU),  &
                             nsize, mmpi_datyp_real4, MpiIdVV, mpiTagVV,  &
                             mmpi_comm_grid, requestIdRecv(numRecv), ierr)
            else
              call mpi_irecv(gd_r8(:, :, kIndexUU),  &
                             nsize, mmpi_datyp_real8, MpiIdVV, mpiTagVV,  &
                             mmpi_comm_grid, requestIdRecv(numRecv), ierr)
            end if

            numSend = numSend + 1
            if (sendrecvKind == 4) then
              gdUV_r4(:, :, kIndexUU) = statevector_in%gdUV(kIndexUU)%r4(:, :, stepIndex)
              call mpi_isend(gdUV_r4(:, :, kIndexUU),  &
                             nsize, mmpi_datyp_real4, MpiIdVV, mpiTagUU,  &
                             mmpi_comm_grid, requestIdSend(numSend), ierr)
            else
              gdUV_r8(:, :, kIndexUU) = statevector_in%gdUV(kIndexUU)%r8(:, :, stepIndex)
              call mpi_isend(gdUV_r8(:, :, kIndexUU),  &
                             nsize, mmpi_datyp_real8, MpiIdVV, mpiTagUU,  &
                             mmpi_comm_grid, requestIdSend(numSend), ierr)
            end if

          else if (mmpi_myid == MpiIDVV) then ! I have VV

            numRecv = numRecv + 1
            if (sendrecvKind == 4) then
              call mpi_irecv(gd_r4(:, :, kIndexVV),  &
                             nsize, mmpi_datyp_real4, MpiIdUU, mpiTagUU,  &
                             mmpi_comm_grid, requestIdRecv(numRecv), ierr)
            else
              call mpi_irecv(gd_r8(:, :, kIndexVV),  &
                             nsize, mmpi_datyp_real8, MpiIdUU, mpiTagUU,  &
                             mmpi_comm_grid, requestIdRecv(numRecv), ierr)
            end if

            numSend = numSend + 1
            if (sendrecvKind == 4) then
              gdUV_r4(:, :, kIndexVV) = statevector_in%gdUV(kIndexVV)%r4(:, :, stepIndex)
              call mpi_isend(gdUV_r4(:, :, kIndexVV),  &
                             nsize, mmpi_datyp_real4, MpiIdUU, mpiTagVV,  &
                             mmpi_comm_grid, requestIdSend(numSend), ierr)
            else
              gdUV_r8(:, :, kIndexVV) = statevector_in%gdUV(kIndexVV)%r8(:, :, stepIndex)
              call mpi_isend(gdUV_r8(:, :, kIndexVV),  &
                             nsize, mmpi_datyp_real8, MpiIdUU, mpiTagVV,  &
                             mmpi_comm_grid, requestIdSend(numSend), ierr)
            end if

          end if

        end do LOOP_KINDEX

        if (numRecv > 0) then
          call mpi_waitAll(numRecv, requestIdRecv(1:numRecv), mpiStatuses(:,1:numRecv), ierr)
        end if

        if (numSend > 0) then
          call mpi_waitAll(numSend, requestIdSend(1:numSend), mpiStatuses(:,1:numSend), ierr)
        end if

        if (statevector_in%UVComponentPresent) then
          do kIndex = statevector_in%myUVkBeg, statevector_in%myUVkEnd
            if (sendrecvKind == 4) then
              statevector_in%gd_r4(:, :, kIndex, stepIndex) =   &
                   statevector_in%gd_r4(:, :, kIndex, stepIndex) +  &
                   gd_r4(:, :, kIndex)
            else
              statevector_in%gd_r8(:, :, kIndex, stepIndex) =   &
                   statevector_in%gd_r8(:, :, kIndex, stepIndex) +  &
                   gd_r8(:, :, kIndex)
            end if
          end do
          if (sendrecvKind == 4) then
            deallocate(gdUV_r4)
            deallocate(gd_r4)
          else
            deallocate(gdUV_r8)
            deallocate(gd_r8)
          end if
        end if

      end if ! UU and VV exist

      ! do allToAll to redistribute main variables from VarsLevs to Tiles
      maxkCount = maxval(statevector_in%allkCount(:))
      if (sendrecvKind == 4) then
        call utl_reAllocate(gd_send_varsLevs_r4, statevector_out%lonPerPEmax, statevector_out%latPerPEmax, &
                            maxkCount, mmpi_nprocs)
        call utl_reAllocate(gd_recv_varsLevs_r4, statevector_out%lonPerPEmax, statevector_out%latPerPEmax, &
                            maxkCount, mmpi_nprocs)
      else
        call utl_reAllocate(gd_send_varsLevs_r8, statevector_out%lonPerPEmax, statevector_out%latPerPEmax, &
                            maxkCount, mmpi_nprocs)
        call utl_reAllocate(gd_recv_varsLevs_r8, statevector_out%lonPerPEmax, statevector_out%latPerPEmax, &
                            maxkCount, mmpi_nprocs)
      end if

      if (sendrecvKind == 4 .and. inKind == 4) then
        call gsv_getField(statevector_in,field_in_r4_ptr)
        !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
        do youridy = 0, (mmpi_npey-1)
          do youridx = 0, (mmpi_npex-1)
            yourid = youridx + youridy*mmpi_npex
            gd_send_varsLevs_r4(1:statevector_out%allLonPerPE(youridx+1),  &
                                1:statevector_out%allLatPerPE(youridy+1),  &
                                1:statevector_in%mykCount, yourid+1) =  &
              field_in_r4_ptr(statevector_out%allLonBeg(youridx+1):statevector_out%allLonEnd(youridx+1),  &
                              statevector_out%allLatBeg(youridy+1):statevector_out%allLatEnd(youridy+1),  &
                              statevector_in%mykBeg:statevector_in%mykEnd, stepIndex)
          end do
        end do
        !$OMP END PARALLEL DO
      else if (sendrecvKind == 8 .and. inKind == 8) then
        call gsv_getField(statevector_in,field_in_r8_ptr)
        !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
        do youridy = 0, (mmpi_npey-1)
          do youridx = 0, (mmpi_npex-1)
            yourid = youridx + youridy*mmpi_npex
            gd_send_varsLevs_r8(1:statevector_out%allLonPerPE(youridx+1),  &
                                1:statevector_out%allLatPerPE(youridy+1),  &
                                1:statevector_in%mykCount, yourid+1) =  &
              field_in_r8_ptr(statevector_out%allLonBeg(youridx+1):statevector_out%allLonEnd(youridx+1),  &
                              statevector_out%allLatBeg(youridy+1):statevector_out%allLatEnd(youridy+1),  &
                              statevector_in%mykBeg:statevector_in%mykEnd, stepIndex)
          end do
        end do
        !$OMP END PARALLEL DO
      else if (sendrecvKind == 4 .and. inKind == 8) then
        call gsv_getField(statevector_in,field_in_r8_ptr)
        !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
        do youridy = 0, (mmpi_npey-1)
          do youridx = 0, (mmpi_npex-1)
            yourid = youridx + youridy*mmpi_npex
            gd_send_varsLevs_r4(1:statevector_out%allLonPerPE(youridx+1),  &
                                1:statevector_out%allLatPerPE(youridy+1),  &
                                1:statevector_in%mykCount, yourid+1) =  &
              real(field_in_r8_ptr(statevector_out%allLonBeg(youridx+1):statevector_out%allLonEnd(youridx+1),  &
                                   statevector_out%allLatBeg(youridy+1):statevector_out%allLatEnd(youridy+1),  &
                                   statevector_in%mykBeg:statevector_in%mykEnd, stepIndex),4)
          end do
        end do
        !$OMP END PARALLEL DO
      else
        call utl_abort('gsv_transposeTilesToVarsLevsAd: Incompatible mix of real 4 and 8 before alltoall mpi comm')
      end if

      nsize = statevector_out%lonPerPEmax * statevector_out%latPerPEmax * maxkCount
      if (mmpi_nprocs > 1) then
        if (sendrecvKind == 4) then
          call rpn_comm_alltoall(gd_send_varsLevs_r4, nsize, 'mpi_real4',  &
                                 gd_recv_varsLevs_r4, nsize, 'mpi_real4', 'grid', ierr)
        else
          call rpn_comm_alltoall(gd_send_varsLevs_r8, nsize, 'mpi_real8',  &
                                 gd_recv_varsLevs_r8, nsize, 'mpi_real8', 'grid', ierr)
        end if
      else
        if (sendrecvKind == 4) then
           gd_recv_varsLevs_r4(:,:,:,1) = gd_send_varsLevs_r4(:,:,:,1)
        else
           gd_recv_varsLevs_r8(:,:,:,1) = gd_send_varsLevs_r8(:,:,:,1)
        end if
      end if

      if (sendrecvKind == 4 .and. outKind == 4) then
        call gsv_getField(statevector_out,field_out_r4_ptr)
        !$OMP PARALLEL DO PRIVATE(yourid)
        do yourid = 0, (mmpi_nprocs-1)
          field_out_r4_ptr(statevector_out%myLonBeg:statevector_out%myLonEnd, &
                           statevector_out%myLatBeg:statevector_out%myLatEnd, &
                           statevector_in%allkBeg(yourid+1):statevector_in%allkEnd(yourid+1), stepIndex) =   &
            gd_recv_varsLevs_r4(1:statevector_out%lonPerPE,  &
                                1:statevector_out%latPerPE,  &
                                1:statevector_in%allkCount(yourid+1), yourid+1)
        end do
        !$OMP END PARALLEL DO
      else if (sendrecvKind == 8 .and. outKind == 8) then
        call gsv_getField(statevector_out,field_out_r8_ptr)
        !$OMP PARALLEL DO PRIVATE(yourid)
        do yourid = 0, (mmpi_nprocs-1)
          field_out_r8_ptr(statevector_out%myLonBeg:statevector_out%myLonEnd, &
                           statevector_out%myLatBeg:statevector_out%myLatEnd, &
                           statevector_in%allkBeg(yourid+1):statevector_in%allkEnd(yourid+1), stepIndex) =   &
            gd_recv_varsLevs_r8(1:statevector_out%lonPerPE,  &
                                1:statevector_out%latPerPE,  &
                                1:statevector_in%allkCount(yourid+1), yourid+1)
        end do
        !$OMP END PARALLEL DO
      else if (sendrecvKind == 4 .and. outKind == 8) then
        call gsv_getField(statevector_out,field_out_r8_ptr)
        !$OMP PARALLEL DO PRIVATE(yourid)
        do yourid = 0, (mmpi_nprocs-1)
          field_out_r8_ptr(statevector_out%myLonBeg:statevector_out%myLonEnd, &
                           statevector_out%myLatBeg:statevector_out%myLatEnd, &
                           statevector_in%allkBeg(yourid+1):statevector_in%allkEnd(yourid+1), stepIndex) =   &
            real(gd_recv_varsLevs_r4(1:statevector_out%lonPerPE,  &
                                     1:statevector_out%latPerPE,  &
                                     1:statevector_in%allkCount(yourid+1), yourid+1),8)
        end do
        !$OMP END PARALLEL DO
      else
        call utl_abort('gsv_transposeTilesToVarsLevsAd: Incompatible mix of real 4 and 8 after alltoall mpi comm')
      end if

    end do ! stepIndex

    ! scatter surface height from task 0 to all others, 1 tile per task
    if (statevector_in%heightSfcPresent .and. statevector_out%heightSfcPresent) then
      allocate(gd_send_height(statevector_out%lonPerPEmax,statevector_out%latPerPEmax,mmpi_nprocs))
      allocate(gd_recv_height(statevector_out%lonPerPEmax,statevector_out%latPerPEmax))
      gd_send_height(:,:,:) = 0.0d0
      gd_recv_height(:,:) = 0.0d0
      field_height_in_ptr => gsv_getHeightSfc(statevector_in)
      field_height_out_ptr => gsv_getHeightSfc(statevector_out)

      if (mmpi_myid == 0) then
        !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
        do youridy = 0, (mmpi_npey-1)
          do youridx = 0, (mmpi_npex-1)
            yourid = youridx + youridy*mmpi_npex
            gd_send_height(1:statevector_out%allLonPerPE(youridx+1),  &
                       1:statevector_out%allLatPerPE(youridy+1), yourid+1) =  &
              field_height_in_ptr(statevector_out%allLonBeg(youridx+1):statevector_out%allLonEnd(youridx+1),  &
                              statevector_out%allLatBeg(youridy+1):statevector_out%allLatEnd(youridy+1))
          end do
        end do
        !$OMP END PARALLEL DO
      end if

      nsize = statevector_out%lonPerPEmax * statevector_out%latPerPEmax
      allocate(displs(mmpi_nprocs))
      allocate(nsizes(mmpi_nprocs))
      do yourid = 0, (mmpi_nprocs-1)
        displs(yourid+1) = yourid*nsize
        nsizes(yourid+1) = nsize
      end do
      call rpn_comm_scatterv(gd_send_height, nsizes, displs, 'mpi_double_precision', &
                             gd_recv_height, nsize, 'mpi_double_precision', &
                             0, 'grid', ierr)

      field_height_out_ptr(statevector_out%myLonBeg:statevector_out%myLonEnd, &
                       statevector_out%myLatBeg:statevector_out%myLatEnd) =   &
        gd_recv_height(1:statevector_out%lonPerPE,  &
                   1:statevector_out%latPerPE)

      deallocate(displs)
      deallocate(nsizes)
      deallocate(gd_recv_height)
      deallocate(gd_send_height)
    end if

    ! Copy over the mask, if it exists
    call gsv_copyMask(statevector_in, statevector_out)

    call utl_tmg_stop(167)

  end subroutine gsv_transposeTilesToVarsLevsAd

  !--------------------------------------------------------------------------
  ! gsv_horizSubSample
  !--------------------------------------------------------------------------
  subroutine gsv_horizSubSample(statevector_in,statevector_out,horizSubSample)
    !
    ! :Purpose: Subsamples the horizontal statevector grid by an integral factor
    !           and transform accordingly the fields 
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)   :: statevector_in
    type(struct_gsv), intent(out)  :: statevector_out
    integer,          intent(in)   :: horizSubSample

    ! Locals:
    integer             :: relativeFactor, middleStep
    real(8)             :: ratio_r8
    integer             :: stepIndex, lonIndex, latIndex, ilon_in, ilon_out, ilat_in, ilat_out
    integer             :: ilon_in1, ilon_in2, lonIndex_in, ilat_in1, ilat_in2, latIndex_in
    character(len=4), pointer :: varNames(:)

    if (statevector_out%allocated) then
      call gsv_deallocate(statevector_out)
    end if

    nullify(varNames)
    call gsv_varNamesList(varNames, statevector_in)

    ! allocate the output statevector
    if (associated(statevector_in%dateStampList)) then
      middleStep = nint((real(statevector_in%numStep,8)+1.0d0)/2.0d0)
      call gsv_allocate(statevector_out,  &
                        statevector_in%numStep, statevector_in%hco, statevector_in%vco, &
                        datestamp_opt = statevector_in%dateStampList(middleStep),  &
                        mpi_local_opt = statevector_in%mpi_local,  &
                        horizSubSample_opt = horizSubSample, varNames_opt=varNames)

    else
      call gsv_allocate(statevector_out,  &
                        statevector_in%numStep, statevector_in%hco, statevector_in%vco, &
                        mpi_local_opt = statevector_in%mpi_local, &
                        horizSubSample_opt = horizSubSample, varNames_opt=varNames)
    end if

    deallocate(varNames)

    if (statevector_out%horizSubSample == statevector_in%horizSubSample) then
      call msg('gsv_horizSubSample', &
           'gsv_horizSubSample: already at the selected subsample level: ' &
           //str(statevector_out%horizSubSample), mpiAll_opt=.false.)
      call gsv_copy(statevector_in,statevector_out)
      return
    end if

    if (statevector_out%horizSubSample > statevector_in%horizSubSample) then

      ! simple averaging onto a coarser grid
      call msg('gsv_horizSubSample', 'increasing subsample level from ' &
           //str(statevector_in%horizSubSample)//' to ' &
           //str(statevector_out%horizSubSample), mpiAll_opt=.false.)

      ratio_r8 = real(statevector_out%horizSubSample,8)/real(statevector_in%horizSubSample,8)
      if (abs(ratio_r8 - real(nint(ratio_r8),8)) > 1.0d-5) then
        call msg('gsv_horizSubSample', &
                             'original subsample level='//str(statevector_in%horizSubSample) &
             //new_line('')//'new      subsample level='//str(statevector_out%horizSubSample))
        call utl_abort('gsv_horizSubSample: relative change of subsample level not an integer')
      end if
      relativeFactor = nint(ratio_r8)

      call gsv_zero(statevector_out)
      !$OMP PARALLEL DO PRIVATE(stepIndex, latIndex, ilat_out, ilat_in1, ilat_in2, &
      !$OMP latIndex_in, lonIndex, ilon_out, ilon_in1, ilon_in2, lonIndex_in)
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
      call msg('gsv_horizSubSample', 'decreasing subsample level from ' &
           //str(statevector_in%horizSubSample)//' to '  &
           //str(statevector_out%horizSubSample), mpiAll_opt=.false.)

      ratio_r8 = real(statevector_in%horizSubSample)/real(statevector_out%horizSubSample)
      if (abs(ratio_r8 - real(nint(ratio_r8),8)) > 1.0d-5) then
        call msg('gsv_horizSubSample', &
                             'original subsample level='//str(statevector_in%horizSubSample) &
             //new_line('')//'new      subsample level='//str(statevector_out%horizSubSample))
        call utl_abort('gsv_horizSubSample: relative change of subsample level not an integer')
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
  ! gsv_transposeStepToVarsLevs
  !--------------------------------------------------------------------------
  subroutine gsv_transposeStepToVarsLevs(stateVector_1step_r4, &
                                         stateVector_VarsLevs, stepIndexBeg)
    !
    ! :Purpose: Transposes the data from a timestep MPI distribution (1 timestep
    !           per MPI task) to the `mpi_distribution='VarsLevs'` distribution.
    !
    ! :Comment: Step-wise distribution is mostly only used for file I/O.
    !           When in such implicit time distribution, it is necessery that
    !           `mpi_local=.false.` and `numStep=1`.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)    :: stateVector_1step_r4
    type(struct_gsv), intent(inout) :: stateVector_VarsLevs
    integer,          intent(in)    :: stepIndexBeg

    ! Locals:
    integer :: ierr, maxkCount, numStepInput, numkToSend, stepIndexInput
    integer :: nsize
    integer :: sendsizes(mmpi_nprocs), recvsizes(mmpi_nprocs), senddispls(mmpi_nprocs), recvdispls(mmpi_nprocs)
    integer :: kIndex, kIndex2, levUV, procIndex, stepIndex
    logical :: thisProcIsAsender(mmpi_nprocs)
    real(4), allocatable :: gd_send_r4(:,:,:), gd_recv_r4(:,:,:)
    real(4), pointer     :: field_in_r4(:,:,:,:), field_out_r4(:,:,:,:)

    call utl_tmg_start(162,'low-level--gsv_stepToVarsLevs')

    if (statevector_VarsLevs%mpi_distribution /= 'VarsLevs') then
      call utl_abort('gsv_transposeStepToVarsLevs: output statevector must have VarsLevs mpi distribution') 
    end if

    ! do mpi transpose to get 4D stateVector into VarsLevs form
    call rpn_comm_barrier('GRID',ierr)
    call msg('gsv_transposeStepToVarsLevs', 'START', verb_opt=2)
    call msg_memUsage('gsv_transposeStepToVarsLevs')

    ! first do surface height, just need to copy over on task 0, assuming task 0 read a file
    if (mmpi_myid == 0 .and. stateVector_VarsLevs%heightSfcPresent) then
      if (.not.stateVector_1step_r4%allocated) then
        call utl_abort('gsv_transposeStepToVarsLevs: Problem with HeightSfc')
      end if
      stateVector_VarsLevs%HeightSfc(:,:) = stateVector_1step_r4%HeightSfc(:,:)
    end if

    ! determine which tasks have something to send and let everyone know
    do procIndex = 1, mmpi_nprocs
      thisProcIsAsender(procIndex) = .false.
      if (mmpi_myid == (procIndex-1) .and. stateVector_1step_r4%allocated) then
        thisProcIsAsender(procIndex) = .true.
      end if
      call rpn_comm_bcast(thisProcIsAsender(procIndex), 1,  &
                          'MPI_LOGICAL', procIndex-1, 'GRID', ierr)
    end do

    numStepInput = 0
    do procIndex = 1, mmpi_nprocs
      if (thisProcIsAsender(procIndex)) numStepInput = numStepInput + 1
    end do
    call msg('gsv_transposeStepToVarsLevs', 'numStepInput = '//str(numStepInput))

    maxkCount = maxval(stateVector_VarsLevs%allkCount(:))
    numkToSend = min(mmpi_nprocs,stateVector_VarsLevs%nk)
    allocate(gd_recv_r4(stateVector_VarsLevs%ni,stateVector_VarsLevs%nj,numStepInput))
    gd_recv_r4(:,:,:) = 0.0
    if (stateVector_1step_r4%allocated) then
      allocate(gd_send_r4(stateVector_VarsLevs%ni,stateVector_VarsLevs%nj,numkToSend))
    else
      allocate(gd_send_r4(1,1,1))
    end if
    gd_send_r4(:,:,:) = 0.0

    call gsv_getField(stateVector_VarsLevs,field_out_r4)

    ! prepare for alltoallv
    nsize = stateVector_VarsLevs%ni * stateVector_VarsLevs%nj

    ! only send the data from tasks with data, same amount to all
    sendsizes(:) = 0
    if (stateVector_1step_r4%allocated) then
      do procIndex = 1, numkToSend
        sendsizes(procIndex) = nsize
      end do
    end if
    senddispls(1) = 0
    do procIndex = 2, mmpi_nprocs
      senddispls(procIndex) = senddispls(procIndex-1) + sendsizes(procIndex-1)
    end do

    ! all tasks recv only from those with data
    recvsizes(:) = 0
    if ((1+mmpi_myid) <= numkToSend) then
      do procIndex = 1, mmpi_nprocs
        if (thisProcIsAsender(procIndex)) then
          recvsizes(procIndex) = nsize
        end if
      end do
    end if
    recvdispls(1) = 0
    do procIndex = 2, mmpi_nprocs
      recvdispls(procIndex) = recvdispls(procIndex-1) + recvsizes(procIndex-1)
    end do

    ! loop to send (at most) 1 level to (at most) all other mpi tasks
    do kIndex = 1, maxkCount

      ! prepare the complete 1 timestep for sending on all tasks that read something
      if (stateVector_1step_r4%allocated) then

        call gsv_getField(stateVector_1step_r4,field_in_r4)
        !$OMP PARALLEL DO PRIVATE(procIndex,kIndex2)
        do procIndex = 1, mmpi_nprocs
          ! compute kIndex value being sent
          kIndex2 = kIndex + stateVector_VarsLevs%allkBeg(procIndex) - 1
          if (kIndex2 <= stateVector_VarsLevs%allkEnd(procIndex)) then
            if(procIndex > numkToSend) then
              call utl_abort('gsv_transposeStepToVarsLevs: ERROR with numkToSend? '&
                   //'procIndex='//str(procIndex) &
                   //', numkToSend = '//str(numkToSend))
            end if
            gd_send_r4(:,:,procIndex) = field_in_r4(:,:,kIndex2,1)
          end if
        end do
        !$OMP END PARALLEL DO

      end if

      call mpi_alltoallv(gd_send_r4, sendsizes, senddispls, mmpi_datyp_real4,  &
                         gd_recv_r4, recvsizes, recvdispls, mmpi_datyp_real4, mmpi_comm_grid, ierr)

      stepIndex = stepIndexBeg - 1
      stepIndexInput = 0
      do procIndex = 1, mmpi_nprocs
        ! skip if this task has nothing to send
        if (.not. thisProcIsAsender(procIndex)) cycle

        stepIndex = stepIndex + 1
        if (stepIndex > stateVector_VarsLevs%numStep) then
          call utl_abort('gsv_transposeStepToVarsLevs: stepIndex > numStep')
        end if

        ! all tasks copy the received step data into correct slot
        kIndex2 = kIndex + stateVector_VarsLevs%mykBeg - 1
        if (kIndex2 <= stateVector_VarsLevs%mykEnd) then
          stepIndexInput = stepIndexInput + 1
          field_out_r4(:,:,kIndex2,stepIndex) = gd_recv_r4(:,:,stepIndexInput)
        end if
      end do

    end do ! kIndex

    ! also do extra transpose for complementary wind components when wind component present
    ! this should really be changed to only send the data needed for UV, not total k
    if (gsv_varExist(stateVector_VarsLevs, 'UU') .or. gsv_varExist(stateVector_VarsLevs, 'VV')) then
      do kIndex = 1, maxkCount

        ! prepare the complete 1 timestep for sending on all tasks that read something
        if (stateVector_1step_r4%allocated) then

          !$OMP PARALLEL DO PRIVATE(procIndex,kIndex2,levUV)
          ! loop over all tasks we are sending to
          do procIndex = 1, mmpi_nprocs
            kIndex2 = kIndex + stateVector_VarsLevs%allkBeg(procIndex) - 1
            if (kIndex2 <= stateVector_VarsLevs%allkEnd(procIndex) .and. &
                kIndex2 >= stateVector_1step_r4%myUVkBeg .and.           &
                kIndex2 <= stateVector_1step_r4%myUVkEnd) then
              gd_send_r4(:,:,procIndex) = stateVector_1step_r4%gdUV(kIndex2)%r4(:,:,1)
            end if
          end do
          !$OMP END PARALLEL DO

        end if

        call mpi_alltoallv(gd_send_r4, sendsizes, senddispls, mmpi_datyp_real4,  &
                           gd_recv_r4, recvsizes, recvdispls, mmpi_datyp_real4, mmpi_comm_grid, ierr)

        stepIndex = stepIndexBeg - 1
        stepIndexInput = 0
        ! loop over all tasks from which we receive something
        do procIndex = 1, mmpi_nprocs
          ! skip if this task did not send anything
          if (.not. thisProcIsAsender(procIndex)) cycle

          stepIndex = stepIndex + 1
          if (stepIndex > stateVector_VarsLevs%numStep) then
            call utl_abort('gsv_transposeStepToVarsLevs: stepIndex > numStep')
          end if

          ! all tasks copy the received step data into correct slot
          kIndex2 = kIndex + stateVector_VarsLevs%mykBeg - 1
          if (kIndex2 <= stateVector_VarsLevs%mykEnd) then
            stepIndexInput = stepIndexInput + 1
            if (kIndex2 >= stateVector_VarsLevs%myUVkBeg .and.  &
                kIndex2 <= stateVector_VarsLevs%myUVkEnd) then
              statevector_varsLevs%gdUV(kIndex2)%r4(:,:,stepIndex) = gd_recv_r4(:,:,stepIndexInput)
            end if
          end if
        end do

      end do ! kIndex

    end if ! do complementary wind component transpose

    deallocate(gd_send_r4)
    deallocate(gd_recv_r4)

    ! Copy the mask if it is present
    if (stateVector_1step_r4%allocated) then
      call gsv_copyMask(stateVector_1step_r4, stateVector_varsLevs)
    end if
    call ocm_communicateMask(stateVector_varsLevs%oceanMask)

    call msg_memUsage('gsv_transposeStepToVarsLevs')
    call msg('gsv_transposeStepToVarsLevs','END', verb_opt=2)

    call utl_tmg_stop(162)

  end subroutine gsv_transposeStepToVarsLevs

  !--------------------------------------------------------------------------
  ! gsv_transposeStepToTiles
  !--------------------------------------------------------------------------
  subroutine gsv_transposeStepToTiles(stateVector_1step, stateVector_tiles, stepIndexBeg)
    !
    ! :Purpose: Transposes the data from a timestep MPI distribution (1 timestep
    !           per MPI task) to the `mpi_distribution='Tiles'` distribution 
    !           (4D lat-lon tiles).
    !
    ! :Comment: Step-wise distribution is mostly only used for file I/O.
    !           When in such implicit time distribution, it is necessery that
    !           `mpi_local=.false.` and `numStep=1`.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)    :: stateVector_1step
    type(struct_gsv), intent(inout) :: stateVector_tiles
    integer,          intent(in)    :: stepIndexBeg

    ! Locals:
    integer :: ierr, yourid, youridx, youridy, nsize, numStepInput, stepCount
    integer :: displs(mmpi_nprocs), nsizes(mmpi_nprocs)
    integer :: senddispls(mmpi_nprocs), sendsizes(mmpi_nprocs)
    integer :: recvdispls(mmpi_nprocs), recvsizes(mmpi_nprocs)
    integer :: kIndex, procIndex, stepIndex, indexBeg, indexEnd
    integer :: inKind, outKind, sendrecvKind, inKindLocal
    logical :: thisProcIsAsender(mmpi_nprocs), allZero, allZero_mpiglobal
    real(8), allocatable :: gd_send_height(:,:,:), gd_recv_height(:,:)
    real(4), allocatable :: gd_send_1d_r4(:), gd_recv_3d_r4(:,:,:)
    real(8), allocatable :: gd_send_1d_r8(:), gd_recv_3d_r8(:,:,:)
    real(4), pointer     :: field_in_r4_ptr(:,:,:,:), field_out_r4_ptr(:,:,:,:)
    real(8), pointer     :: field_in_r8_ptr(:,:,:,:), field_out_r8_ptr(:,:,:,:)

    call rpn_comm_barrier('GRID',ierr)

    call utl_tmg_start(163,'low-level--gsv_stepToTiles')

    if (statevector_tiles%mpi_distribution /= 'Tiles') then
      call utl_abort('gsv_transposeStepToTiles: output statevector must have Tiles mpi distribution')
    end if

    call msg('gsv_transposeStepToTiles','START', verb_opt=2)

    ! determine which tasks have something to send and let everyone know
    do procIndex = 1, mmpi_nprocs
      thisProcIsAsender(procIndex) = .false.
      if (mmpi_myid == (procIndex-1) .and. stateVector_1step%allocated) then
        thisProcIsAsender(procIndex) = .true.
      end if
      call rpn_comm_bcast(thisProcIsAsender(procIndex), 1,  &
                          'MPI_LOGICAL', procIndex-1, 'GRID', ierr)
    end do

    ! only send the data from tasks with data, same amount to all
    sendsizes(:) = 0
    if (stateVector_1step%allocated) then
      do youridy = 0, (mmpi_npey-1)
        do youridx = 0, (mmpi_npex-1)
          yourid = youridx + youridy*mmpi_npex
          nsize = stateVector_tiles%allLonPerPE(youridx+1) * stateVector_tiles%allLatPerPE(youridy+1)
          sendsizes(yourid+1) = nsize
        end do
      end do
    end if
    senddispls(:) = 0
    do yourid = 1, (mmpi_nprocs-1)
      senddispls(yourid+1) = senddispls(yourid) + sendsizes(yourid)
    end do

    ! all tasks recv, but only from those with data
    recvsizes(:) = 0
    nsize = stateVector_tiles%lonPerPE * stateVector_tiles%latPerPE
    do yourid = 0, (mmpi_nprocs-1) ! recv from this task
      if (thisProcIsAsender(yourid+1)) then
        recvsizes(yourid+1) = nsize
      end if
    end do
    recvdispls(:) = 0
    do yourid = 1, (mmpi_nprocs-1)
      recvdispls(yourid+1) = recvdispls(yourid) + recvsizes(yourid)
    end do

    numStepInput = 0
    do yourid = 0, (mmpi_nprocs-1)
      if (thisProcIsAsender(yourid+1)) numStepInput = numStepInput + 1
    end do

    if (stateVector_1step%allocated) then
      inKindLocal = stateVector_1step%dataKind
    else
      inKindLocal = -1
    end if
    call rpn_comm_allreduce(inKindLocal, inKind, 1,  &
                            'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)
    outKind = stateVector_tiles%dataKind
    if (inKind == 4 .or. outKind == 4) then
      sendrecvKind = 4
    else
      sendrecvKind = 8
    end if

    ! allocate arrays used for mpi communication of 1 level/variable at a time
    if (sendrecvKind == 4) then
      allocate(gd_recv_3d_r4(stateVector_tiles%lonPerPE,stateVector_tiles%latPerPE,numStepInput))
      gd_recv_3d_r4(:,:,:) = 0.0
      if (stateVector_1step%allocated) then
        allocate(gd_send_1d_r4(stateVector_tiles%ni*stateVector_tiles%nj))
      else
        allocate(gd_send_1d_r4(1))
      end if
      gd_send_1d_r4(:) = 0.0
    else
      allocate(gd_recv_3d_r8(stateVector_tiles%lonPerPE,stateVector_tiles%latPerPE,numStepInput))
      gd_recv_3d_r8(:,:,:) = 0.0d0
      if (stateVector_1step%allocated) then
        allocate(gd_send_1d_r8(stateVector_tiles%ni*stateVector_tiles%nj))
      else
        allocate(gd_send_1d_r8(1))
      end if
      gd_send_1d_r8(:) = 0.0d0
    end if

    if (stateVector_tiles%dataKind == 4) then
      call gsv_getField(stateVector_tiles,field_out_r4_ptr)
    else if (stateVector_tiles%dataKind == 8) then
      call gsv_getField(stateVector_tiles,field_out_r8_ptr)
    else
      call utl_abort('gsv_transposeStepToTiles: stateVector_tiles%dataKind not real 4 or 8')
    end if

    kIndex_Loop: do kIndex = 1, stateVector_tiles%nk

      ! determine if there is data to send for this kIndex
      if (stateVector_1step%allocated) then
        if (stateVector_1step%dataKind == 4) then
          call gsv_getField(stateVector_1step,field_in_r4_ptr)
          allZero = (maxval(abs(field_in_r4_ptr(:, :, kIndex, 1))) == 0.0)
        else if (stateVector_1step%dataKind == 8) then
          call gsv_getField(stateVector_1step,field_in_r8_ptr)
          allZero = (maxval(abs(field_in_r8_ptr(:, :, kIndex, 1))) == 0.0D0)
        end if
      else
        allZero = .true.
      end if
      call rpn_comm_allReduce(allZero,allZero_mpiglobal,1,'mpi_logical','mpi_land','GRID',ierr)
      if (allZero_mpiglobal) then
        ! Field equal to zero, skipping this kIndex to save time
        cycle kIndex_Loop
      end if

      ! prepare the complete 1 timestep for sending on all tasks that have read something
      if (stateVector_1step%allocated) then
        !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid,nsize,indexBeg,indexEnd)
        do youridy = 0, (mmpi_npey-1)
          do youridx = 0, (mmpi_npex-1)
            yourid = youridx + youridy*mmpi_npex
            nsize = stateVector_tiles%allLonPerPE(youridx+1) *  &
                    stateVector_tiles%allLatPerPE(youridy+1)
            indexBeg = senddispls(yourid+1) + 1
            indexEnd = senddispls(yourid+1) + nsize
            if (sendrecvKind == 4 .and. stateVector_1step%dataKind == 4) then
              gd_send_1d_r4(indexBeg:indexEnd) =  &
                        reshape(field_in_r4_ptr(stateVector_tiles%allLonBeg(youridx+1):stateVector_tiles%allLonEnd(youridx+1), &
                                                 stateVector_tiles%allLatBeg(youridy+1):stateVector_tiles%allLatEnd(youridy+1), &
                                                 kIndex, 1), (/ nsize /))
            else if (sendrecvKind == 4 .and. stateVector_1step%dataKind == 8) then
              gd_send_1d_r4(indexBeg:indexEnd) =  &
                   real(reshape(field_in_r8_ptr(stateVector_tiles%allLonBeg(youridx+1):stateVector_tiles%allLonEnd(youridx+1), &
                                                  stateVector_tiles%allLatBeg(youridy+1):stateVector_tiles%allLatEnd(youridy+1), &
                                                  kIndex, 1), (/ nsize /)), 4)
            else if (sendrecvKind == 8 .and. stateVector_1step%dataKind == 8) then
              gd_send_1d_r8(indexBeg:indexEnd) =  &
                        reshape(field_in_r8_ptr(stateVector_tiles%allLonBeg(youridx+1):stateVector_tiles%allLonEnd(youridx+1), &
                                                 stateVector_tiles%allLatBeg(youridy+1):stateVector_tiles%allLatEnd(youridy+1), &
                                                 kIndex, 1), (/ nsize /))
            else
              call utl_abort('gsv_stepToTiles: unexpected combination of real kinds')
            end if
          end do
        end do
        !$OMP END PARALLEL DO
      end if

      if (sendrecvKind == 4) then
        call mpi_alltoallv(gd_send_1d_r4, sendsizes, senddispls, mmpi_datyp_real4, &
                           gd_recv_3d_r4, recvsizes, recvdispls, mmpi_datyp_real4, &
                           mmpi_comm_grid, ierr)
      else if (sendrecvKind == 8) then
        call mpi_alltoallv(gd_send_1d_r8, sendsizes, senddispls, mmpi_datyp_real8, &
                           gd_recv_3d_r8, recvsizes, recvdispls, mmpi_datyp_real8, &
                           mmpi_comm_grid, ierr)
      end if

      stepIndex = stepIndexBeg - 1
      stepCount = 0
      do procIndex = 1, mmpi_nprocs

        ! skip if this task had nothing to send
        if (.not. thisProcIsAsender(procIndex)) cycle

        stepCount = stepCount + 1
        stepIndex = stepIndex + 1
        if (stepIndex > stateVector_tiles%numStep) then
          call utl_abort('gsv_transposeStepToTiles: stepIndex > numStep')
        end if
        if (sendrecvKind == 4 .and. stateVector_tiles%dataKind == 4) then
          field_out_r4_ptr(stateVector_tiles%myLonBeg:stateVector_tiles%myLonEnd,  &
                           stateVector_tiles%myLatBeg:stateVector_tiles%myLatEnd,  &
                           kIndex, stepIndex) =   &
              gd_recv_3d_r4(1:stateVector_tiles%lonPerPE,1:stateVector_tiles%latPerPE,stepCount)
        else if (sendrecvKind == 4 .and. stateVector_tiles%dataKind == 8) then
          field_out_r8_ptr(stateVector_tiles%myLonBeg:stateVector_tiles%myLonEnd,  &
                           stateVector_tiles%myLatBeg:stateVector_tiles%myLatEnd,  &
                           kIndex, stepIndex) =   &
              real(gd_recv_3d_r4(1:stateVector_tiles%lonPerPE,1:stateVector_tiles%latPerPE,stepCount),8)
        else if (sendrecvKind == 8 .and. stateVector_tiles%dataKind == 8) then
          field_out_r8_ptr(stateVector_tiles%myLonBeg:stateVector_tiles%myLonEnd,  &
                           stateVector_tiles%myLatBeg:stateVector_tiles%myLatEnd,  &
                           kIndex, stepIndex) =   &
              gd_recv_3d_r8(1:stateVector_tiles%lonPerPE,1:stateVector_tiles%latPerPE,stepCount)
        end if

      end do ! procIndex
    end do kIndex_Loop ! kIndex

    if (sendrecvKind == 4) then
      deallocate(gd_recv_3d_r4)
      deallocate(gd_send_1d_r4)
    else if (sendrecvKind == 8) then
      deallocate(gd_recv_3d_r8)
      deallocate(gd_send_1d_r8)
    end if

    ! now send HeightSfc from task 0 to all others
    if (stateVector_tiles%heightSfcPresent) then

      allocate(gd_recv_height(stateVector_tiles%lonPerPEmax,stateVector_tiles%latPerPEmax))
      gd_recv_height(:,:) = 0.0d0

      ! prepare data to send from task 0
      if (mmpi_myid == 0) then
        if (.not. stateVector_1step%allocated) then
          call utl_abort('gsv_transposeStepToVarsLevs: Problem with HeightSfc')
        end if

        allocate(gd_send_height(stateVector_tiles%lonPerPEmax,stateVector_tiles%latPerPEmax,mmpi_nprocs))
        gd_send_height(:,:,:) = 0.0d0

        !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
        do youridy = 0, (mmpi_npey-1)
          do youridx = 0, (mmpi_npex-1)
            yourid = youridx + youridy*mmpi_npex
            gd_send_height(1:stateVector_tiles%allLonPerPE(youridx+1),  &
                    1:stateVector_tiles%allLatPerPE(youridy+1), yourid+1) =  &
                real(stateVector_1step%HeightSfc(&
                     stateVector_tiles%allLonBeg(youridx+1):stateVector_tiles%allLonEnd(youridx+1), &
                     stateVector_tiles%allLatBeg(youridy+1):stateVector_tiles%allLatEnd(youridy+1)), 8)
          end do
        end do
        !$OMP END PARALLEL DO

      else
        allocate(gd_send_height(1,1,1))
      end if

      ! distribute from task 0 to all tasks
      nsize = stateVector_tiles%lonPerPEmax * stateVector_tiles%latPerPEmax
      do procIndex = 1, mmpi_nprocs
        displs(procIndex) = (procIndex-1)*nsize
        nsizes(procIndex) = nsize
      end do
      call rpn_comm_scatterv(gd_send_height, nsizes, displs, 'mpi_real8', &
                             gd_recv_height, nsize, 'mpi_real8', &
                             0, 'grid', ierr)

      stateVector_tiles%HeightSfc(&
                stateVector_tiles%myLonBeg:stateVector_tiles%myLonEnd,    &
                stateVector_tiles%myLatBeg:stateVector_tiles%myLatEnd) =  &
          gd_recv_height(1:stateVector_tiles%lonPerPE,1:stateVector_tiles%latPerPE)

      deallocate(gd_recv_height)
      deallocate(gd_send_height)

    end if ! heightSfcPresent

    ! Copy the mask if it is present
    if (stateVector_1step%allocated) then
      call gsv_copyMask(stateVector_1step, stateVector_tiles)
    end if
    call ocm_communicateMask(stateVector_tiles%oceanMask)

    call msg('gsv_transposeStepToTiles','END', verb_opt=2)

    call utl_tmg_stop(163)

  end subroutine gsv_transposeStepToTiles

  !--------------------------------------------------------------------------
  ! gsv_transposeTilesToStep
  !--------------------------------------------------------------------------
  subroutine gsv_transposeTilesToStep(stateVector_1step, stateVector_tiles, stepIndexBeg)
    !
    ! :Purpose: Transposes the data from a `mpi_distribution='Tiles'` distribution 
    !           (4D lat-lon tiles) to a timestep MPI distribution (1 timestep per
    !           MPI task)
    !
    ! :Comment: Step-wise distribution is mostly only used for file I/O.
    !           When in such implicit time distribution, it is necessery that
    !           `mpi_local=.false.` and `numStep=1`.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout)  :: stateVector_1step
    type(struct_gsv), intent(in)     :: stateVector_tiles
    integer,          intent(in)     :: stepIndexBeg

    ! Locals:
    integer :: ierr, yourid, youridx, youridy, nsize
    integer :: kIndex, procIndex, stepIndex, numStepOutput, stepCount
    integer :: inKind, outKind, sendrecvKind, outKindLocal
    logical :: thisProcIsAreceiver(mmpi_nprocs)
    integer :: sendsizes(mmpi_nprocs), recvsizes(mmpi_nprocs), senddispls(mmpi_nprocs), recvdispls(mmpi_nprocs)
    real(4), allocatable :: gd_send_r4(:,:,:), gd_recv_r4(:,:,:)
    real(8), allocatable :: gd_send_r8(:,:,:), gd_recv_r8(:,:,:)
    real(8), allocatable :: gd_send_height(:,:), gd_recv_height(:,:,:)
    real(4), pointer     :: field_in_r4_ptr(:,:,:,:), field_out_r4_ptr(:,:,:,:)
    real(8), pointer     :: field_in_r8_ptr(:,:,:,:), field_out_r8_ptr(:,:,:,:)

    if (statevector_tiles%mpi_distribution /= 'Tiles') then
      call utl_abort('gsv_transposeTilesToStep: input statevector must have Tiles mpi distribution')
    end if

    call rpn_comm_barrier('GRID',ierr)
    call msg('gsv_transposeTilesToStep', 'START', verb_opt=2)

    ! determine which tasks have something to receive and let everyone know
    do procIndex = 1, mmpi_nprocs
      thisProcIsAreceiver(procIndex) = .false.
      if (mmpi_myid == (procIndex-1) .and. stateVector_1step%allocated) then
        thisProcIsAreceiver(procIndex) = .true.
      end if
      call rpn_comm_bcast(thisProcIsAreceiver(procIndex), 1,  &
                          'MPI_LOGICAL', procIndex-1, 'GRID', ierr)
    end do

    numStepOutput = 0
    do procIndex = 1, mmpi_nprocs
      if (thisProcIsAreceiver(procIndex)) numStepOutput = numStepOutput + 1
    end do
    call msg('gsv_transposeTilesToStep','numStepOutput = '//str(numStepOutput))

    ! size of each message
    nsize = stateVector_tiles%lonPerPEmax * stateVector_tiles%latPerPEmax

    ! only recv the data on tasks that need data
    recvsizes(:) = 0
    if (stateVector_1step%allocated) then
      do procIndex = 1, mmpi_nprocs
        recvsizes(procIndex) = nsize
      end do
    end if
    recvdispls(1) = 0
    do procIndex = 2, mmpi_nprocs
      recvdispls(procIndex) = recvdispls(procIndex-1) + recvsizes(procIndex-1)
    end do

    ! all tasks send only to those that need data
    sendsizes(:) = 0
    do procIndex = 1, mmpi_nprocs
      if (thisProcIsAreceiver(procIndex)) then
        sendsizes(procIndex) = nsize
      end if
    end do
    senddispls(1) = 0
    do procIndex = 2, mmpi_nprocs
      senddispls(procIndex) = senddispls(procIndex-1) + sendsizes(procIndex-1)
    end do

    if (stateVector_1step%allocated) then
      outKindLocal = stateVector_1step%dataKind
    else
      outKindLocal = -1
    end if
    call rpn_comm_allreduce(outKindLocal, outKind, 1,  &
                            'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)
    inKind = stateVector_tiles%dataKind
    if (inKind == 4 .or. outKind == 4) then
      sendrecvKind = 4
    else
      sendrecvKind = 8
    end if

    ! allocate arrays used for mpi communication of 1 level/variable at a time
    if (sendrecvKind == 4) then
      allocate(gd_send_r4(stateVector_tiles%lonPerPEmax,stateVector_tiles%latPerPEmax,numStepOutput))
      gd_send_r4(:,:,:) = 0.0
      if (stateVector_1step%allocated) then
        allocate(gd_recv_r4(stateVector_tiles%lonPerPEmax,stateVector_tiles%latPerPEmax,mmpi_nprocs))
      else
        allocate(gd_recv_r4(1,1,1))
      end if
      gd_recv_r4(:,:,:) = 0.0
    else
      allocate(gd_send_r8(stateVector_tiles%lonPerPEmax,stateVector_tiles%latPerPEmax,numStepOutput))
      gd_send_r8(:,:,:) = 0.0d0
      if (stateVector_1step%allocated) then
        allocate(gd_recv_r8(stateVector_tiles%lonPerPEmax,stateVector_tiles%latPerPEmax,mmpi_nprocs))
      else
        allocate(gd_recv_r8(1,1,1))
      end if
      gd_recv_r8(:,:,:) = 0.0d0
    end if

    if (stateVector_tiles%dataKind == 4) then
      call gsv_getField(stateVector_tiles,field_in_r4_ptr)
    else
      call gsv_getField(stateVector_tiles,field_in_r8_ptr)
    end if

    do kIndex = 1, stateVector_tiles%nk

      ! prepare data for distribution from all tasks to only those that need it
      stepIndex = stepIndexBeg - 1
      stepCount = 0

      do procIndex = 1, mmpi_nprocs

        ! skip if this task has nothing to receive
        if (.not. thisProcIsAreceiver(procIndex)) cycle

        stepCount = stepCount + 1
        stepIndex = stepIndex + 1
        if (stepIndex > stateVector_tiles%numStep) then
          call utl_abort('gsv_transposeTilesToStep: stepIndex > numStep')
        end if

        if (sendrecvKind == 4 .and. stateVector_tiles%dataKind == 4) then
          gd_send_r4(1:stateVector_tiles%lonPerPE,1:stateVector_tiles%latPerPE,stepCount) = &
               field_in_r4_ptr(stateVector_tiles%myLonBeg:stateVector_tiles%myLonEnd,  &
                               stateVector_tiles%myLatBeg:stateVector_tiles%myLatEnd,  &
                               kIndex, stepIndex)
        else if (sendrecvKind == 4 .and. stateVector_tiles%dataKind == 8) then
          gd_send_r4(1:stateVector_tiles%lonPerPE,1:stateVector_tiles%latPerPE,stepCount) = &
               real(field_in_r8_ptr(stateVector_tiles%myLonBeg:stateVector_tiles%myLonEnd,  &
                                    stateVector_tiles%myLatBeg:stateVector_tiles%myLatEnd,  &
                                    kIndex, stepIndex),4)
        else if (sendrecvKind == 8 .and. stateVector_tiles%dataKind == 8) then
          gd_send_r8(1:stateVector_tiles%lonPerPE,1:stateVector_tiles%latPerPE,stepCount) = &
               field_in_r8_ptr(stateVector_tiles%myLonBeg:stateVector_tiles%myLonEnd,  &
                               stateVector_tiles%myLatBeg:stateVector_tiles%myLatEnd,  &
                               kIndex, stepIndex)
        else
          call utl_abort('gsv_transposeTilesToStep: unexpected combination of real kinds')
        end if

      end do ! procIndex

      if (sendrecvKind == 4) then
        call mpi_alltoallv(gd_send_r4, sendsizes, senddispls, mmpi_datyp_real4, &
                           gd_recv_r4, recvsizes, recvdispls, mmpi_datyp_real4, &
                           mmpi_comm_grid, ierr)
      else if (sendrecvKind == 8) then
        call mpi_alltoallv(gd_send_r8, sendsizes, senddispls, mmpi_datyp_real8, &
                           gd_recv_r8, recvsizes, recvdispls, mmpi_datyp_real8, &
                           mmpi_comm_grid, ierr)
      end if

      ! copy over the complete 1 timestep received
      if (stateVector_1step%allocated) then

        if (sendrecvKind == 4 .and. stateVector_1step%dataKind == 4) then

          call gsv_getField(stateVector_1step,field_out_r4_ptr)
          !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
          do youridy = 0, (mmpi_npey-1)
            do youridx = 0, (mmpi_npex-1)
              yourid = youridx + youridy*mmpi_npex
              field_out_r4_ptr(stateVector_tiles%allLonBeg(youridx+1):stateVector_tiles%allLonEnd(youridx+1), &
                               stateVector_tiles%allLatBeg(youridy+1):stateVector_tiles%allLatEnd(youridy+1), &
                               kIndex, 1) = &
                    gd_recv_r4(1:stateVector_tiles%allLonPerPE(youridx+1),  &
                               1:stateVector_tiles%allLatPerPE(youridy+1), yourid+1)
            end do
          end do
          !$OMP END PARALLEL DO

        else if (sendrecvKind == 4 .and. stateVector_1step%dataKind == 8) then

          call gsv_getField(stateVector_1step,field_out_r8_ptr)
          !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
          do youridy = 0, (mmpi_npey-1)
            do youridx = 0, (mmpi_npex-1)
              yourid = youridx + youridy*mmpi_npex
              field_out_r8_ptr(stateVector_tiles%allLonBeg(youridx+1):stateVector_tiles%allLonEnd(youridx+1), &
                               stateVector_tiles%allLatBeg(youridy+1):stateVector_tiles%allLatEnd(youridy+1), &
                               kIndex, 1) = &
               real(gd_recv_r4(1:stateVector_tiles%allLonPerPE(youridx+1),  &
                               1:stateVector_tiles%allLatPerPE(youridy+1), yourid+1), 8)
            end do
          end do
          !$OMP END PARALLEL DO

        else if (sendrecvKind == 8 .and. stateVector_1step%dataKind == 8) then

          call gsv_getField(stateVector_1step,field_out_r8_ptr)
          !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
          do youridy = 0, (mmpi_npey-1)
            do youridx = 0, (mmpi_npex-1)
              yourid = youridx + youridy*mmpi_npex
              field_out_r8_ptr(stateVector_tiles%allLonBeg(youridx+1):stateVector_tiles%allLonEnd(youridx+1), &
                               stateVector_tiles%allLatBeg(youridy+1):stateVector_tiles%allLatEnd(youridy+1), &
                               kIndex, 1) = &
                    gd_recv_r8(1:stateVector_tiles%allLonPerPE(youridx+1),  &
                               1:stateVector_tiles%allLatPerPE(youridy+1), yourid+1)
            end do
          end do
          !$OMP END PARALLEL DO

        else
          call utl_abort('gsv_transposeTilesToStep: unexpected combination of real kinds')
        end if

      end if

    end do ! kIndex

    if (sendrecvKind == 4) then
      deallocate(gd_recv_r4)
      deallocate(gd_send_r4)
    else if (sendrecvKind == 8) then
      deallocate(gd_recv_r8)
      deallocate(gd_send_r8)
    end if

    ! now gather the same HeightSfc onto each task that is a receiver
    if (stateVector_tiles%heightSfcPresent) then

      allocate(gd_send_height(stateVector_tiles%lonPerPEmax,stateVector_tiles%latPerPEmax))
      gd_send_height(:,:) = 0.0d0
      if (stateVector_1step%allocated) then
        allocate(gd_recv_height(stateVector_tiles%lonPerPEmax,stateVector_tiles%latPerPEmax,mmpi_nprocs))
      else
        allocate(gd_recv_height(1,1,1))
      end if
      gd_recv_height(:,:,:) = 0.0d0

      ! prepare tile to send on each task
      gd_send_height(1:stateVector_tiles%lonPerPE,1:stateVector_tiles%latPerPE) = &
          stateVector_tiles%HeightSfc(stateVector_tiles%myLonBeg:stateVector_tiles%myLonEnd,    &
                                  stateVector_tiles%myLatBeg:stateVector_tiles%myLatEnd)

      ! gather from all tasks onto each task with a receiving statevector
      do procIndex = 1, mmpi_nprocs

        ! skip if this task has nothing to receive
        if (.not. thisProcIsAreceiver(procIndex)) cycle

        nsize = stateVector_tiles%lonPerPEmax * stateVector_tiles%latPerPEmax
        call rpn_comm_gather(gd_send_height, nsize, 'mpi_real8', &
                             gd_recv_height, nsize, 'mpi_real8', &
                             procIndex-1, 'grid', ierr)

        ! copy over the complete 1 timestep received
        if (mmpi_myid == procIndex-1) then
          !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
          do youridy = 0, (mmpi_npey-1)
            do youridx = 0, (mmpi_npex-1)
              yourid = youridx + youridy*mmpi_npex
              stateVector_1step%HeightSfc(&
                   stateVector_tiles%allLonBeg(youridx+1):stateVector_tiles%allLonEnd(youridx+1), &
                   stateVector_tiles%allLatBeg(youridy+1):stateVector_tiles%allLatEnd(youridy+1)) = &
                   gd_recv_height(1:stateVector_tiles%allLonPerPE(youridx+1),  &
                              1:stateVector_tiles%allLatPerPE(youridy+1), yourid+1)
            end do
          end do
          !$OMP END PARALLEL DO
        end if

      end do ! procIndex

      deallocate(gd_recv_height)
      deallocate(gd_send_height)

    end if ! heightSfcPresent

    ! Copy mask if it exists on mpi task with step data allocated
    if (stateVector_1step%allocated) then
      call gsv_copyMask(stateVector_tiles, stateVector_1step)
    end if

    call msg('gsv_transposeTilesToStep', 'END', verb_opt=2)

  end subroutine gsv_transposeTilesToStep

  !--------------------------------------------------------------------------
  ! gsv_transposeTilesToMpiGlobal
  !--------------------------------------------------------------------------
  subroutine gsv_transposeTilesToMpiGlobal(stateVector_mpiGlobal, stateVector_tiles)
    !
    !:Purpose: Does MPI transpose (allGather) from `mpi_distribution='Tiles'` 
    !          (4D lat-lon tiles) to global 4D stateVector on each MPI task 
    !          where it is allocated.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout)  :: stateVector_mpiGlobal
    type(struct_gsv), intent(in)     :: stateVector_tiles

    ! Locals:
    integer :: ierr, yourid, youridx, youridy, nsize
    integer :: kIndex, stepIndex, numStep

    real(4), allocatable :: gd_send_r4(:,:), gd_recv_r4(:,:,:)
    real(8), allocatable :: gd_send_r8(:,:), gd_recv_r8(:,:,:)
    real(4), pointer     :: field_out_r4(:,:,:,:), field_in_r4(:,:,:,:)
    real(8), pointer     :: field_out_r8(:,:,:,:), field_in_r8(:,:,:,:)

    if (stateVector_tiles%mpi_distribution /= 'Tiles') then
      call utl_abort('gsv_transposeTilesToMpiGlobal: input statevector must have Tiles mpi distribution')
    end if

    if (stateVector_mpiGlobal%allocated) then
      if (stateVector_mpiGlobal%numStep /= stateVector_tiles%numStep) then
        call utl_abort('gsv_transposeTilesToMpiGlobal: input and output ' // &
                       'stateVectors must have same numStep')
      end if
    end if

    numStep = stateVector_tiles%numStep

    call rpn_comm_barrier('GRID',ierr)
    call msg('gsv_transposeTilesToMpiGlobal', 'START', verb_opt=2)

    ! size of each message
    nsize = stateVector_tiles%lonPerPEmax * stateVector_tiles%latPerPEmax

    ! allocate arrays used for mpi communication of 1 level/variable at a time
    allocate(gd_send_r4(stateVector_tiles%lonPerPEmax,  &
                        stateVector_tiles%latPerPEmax))
    gd_send_r4(:,:) = 0.0
    allocate(gd_recv_r4(stateVector_tiles%lonPerPEmax,  &
                        stateVector_tiles%latPerPEmax, mmpi_nprocs))
    gd_recv_r4(:,:,:) = 0.0

    if (stateVector_tiles%dataKind == 4) then
      call gsv_getField(stateVector_tiles,field_in_r4)
    else
      call gsv_getField(stateVector_tiles,field_in_r8)
    end if
    if (stateVector_mpiGlobal%allocated) then
      if (stateVector_mpiGlobal%dataKind == 4) then
        call gsv_getField(stateVector_mpiGlobal,field_out_r4)
      else
        call gsv_getField(stateVector_mpiGlobal,field_out_r8)
      end if
    end if

    ! do allGather for 1 2D field/stepIndex at a time
    do stepIndex = 1, numStep
      do kIndex = 1, stateVector_tiles%nk
        if (stateVector_tiles%dataKind == 4) then
          gd_send_r4(1:stateVector_tiles%lonPerPE,1:stateVector_tiles%latPerPE) =  &
               field_in_r4(stateVector_tiles%myLonBeg:stateVector_tiles%myLonEnd,  &
                           stateVector_tiles%myLatBeg:stateVector_tiles%myLatEnd,  &
                           kIndex, stepIndex)
        else
          gd_send_r4(1:stateVector_tiles%lonPerPE,1:stateVector_tiles%latPerPE) =      &
               real(field_in_r8(stateVector_tiles%myLonBeg:stateVector_tiles%myLonEnd, &
                                stateVector_tiles%myLatBeg:stateVector_tiles%myLatEnd, &
                                kIndex, stepIndex), 4)
        end if

        call rpn_comm_allgather(gd_send_r4, nsize, 'mpi_real4',  &
                                gd_recv_r4, nsize, 'mpi_real4', 'grid', ierr)

        ! copy over the complete 2D field for 1 stepIndex received
        if (stateVector_mpiGlobal%allocated) then

          !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
          do youridy = 0, (mmpi_npey-1)
            do youridx = 0, (mmpi_npex-1)
              yourid = youridx + youridy*mmpi_npex
              if (stateVector_mpiGlobal%dataKind == 4) then
                field_out_r4(stateVector_tiles%allLonBeg(youridx+1):stateVector_tiles%allLonEnd(youridx+1), &
                             stateVector_tiles%allLatBeg(youridy+1):stateVector_tiles%allLatEnd(youridy+1), &
                             kIndex, stepIndex) = &
                     gd_recv_r4(1:stateVector_tiles%allLonPerPE(youridx+1),  &
                                1:stateVector_tiles%allLatPerPE(youridy+1), yourid+1)
              else
                field_out_r8(stateVector_tiles%allLonBeg(youridx+1):stateVector_tiles%allLonEnd(youridx+1), &
                             stateVector_tiles%allLatBeg(youridy+1):stateVector_tiles%allLatEnd(youridy+1), &
                             kIndex, stepIndex) = &
                     real(gd_recv_r4(1:stateVector_tiles%allLonPerPE(youridx+1),  &
                                      1:stateVector_tiles%allLatPerPE(youridy+1), yourid+1), 4)
              end if
            end do
          end do
          !$OMP END PARALLEL DO

        end if

      end do ! kIndex
    end do ! stepIndex

    deallocate(gd_recv_r4)
    deallocate(gd_send_r4)

    ! now gather the same HeightSfc onto each task that is a receiver
    if (stateVector_tiles%heightSfcPresent) then

      allocate(gd_send_r8(stateVector_tiles%lonPerPEmax,stateVector_tiles%latPerPEmax))
      gd_send_r8(:,:) = 0.0d0
      allocate(gd_recv_r8(stateVector_tiles%lonPerPEmax,stateVector_tiles%latPerPEmax,mmpi_nprocs))
      gd_recv_r8(:,:,:) = 0.0d0

      ! prepare tile to send on each task
      gd_send_r8(1:stateVector_tiles%lonPerPE,1:stateVector_tiles%latPerPE) = &
          stateVector_tiles%HeightSfc(stateVector_tiles%myLonBeg:stateVector_tiles%myLonEnd,    &
                                      stateVector_tiles%myLatBeg:stateVector_tiles%myLatEnd)

      call rpn_comm_allGather(gd_send_r8, nsize, 'mpi_real8',  &
                              gd_recv_r8, nsize, 'mpi_real8', 'grid', ierr)

      if (stateVector_mpiGlobal%allocated) then

        !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
        do youridy = 0, (mmpi_npey-1)
          do youridx = 0, (mmpi_npex-1)
            yourid = youridx + youridy*mmpi_npex
            stateVector_mpiGlobal%HeightSfc(&
                   stateVector_tiles%allLonBeg(youridx+1):stateVector_tiles%allLonEnd(youridx+1), &
                   stateVector_tiles%allLatBeg(youridy+1):stateVector_tiles%allLatEnd(youridy+1)) = &
                   gd_recv_r8(1:stateVector_tiles%allLonPerPE(youridx+1),  &
                              1:stateVector_tiles%allLatPerPE(youridy+1), yourid+1)
          end do
        end do
        !$OMP END PARALLEL DO

      end if

      deallocate(gd_recv_r8)
      deallocate(gd_send_r8)

    end if ! heightSfcPresent

    ! Copy over the mask, if it exists
    if (stateVector_mpiGlobal%allocated) then
      call gsv_copyMask(stateVector_tiles, stateVector_mpiGlobal)
    end if

    call msg('gsv_transposeTilesToMpiGlobal', 'END', verb_opt=2)

  end subroutine gsv_transposeTilesToMpiGlobal

  !--------------------------------------------------------------------------
  ! gsv_varKindExist
  !--------------------------------------------------------------------------
  function gsv_varKindExist(varKind) result(KindFound)
    !
    ! :Purpose: To check whether any of the variables to be assimilated
    !           (i.e. specified in the namelist NAMSTATE) are part of the
    !           specified variable kind
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: varKind    ! Variable kind (e.g. MT or CH)
    logical                      :: KindFound  ! Logical indicating whether var kind found 

    ! Locals:
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
  ! gsv_dotProduct
  !--------------------------------------------------------------------------
  subroutine gsv_dotProduct(stateVector_a,stateVector_b,dotsum)
    !
    ! :Purpose: Computes the dot product of two statevectors
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)  :: stateVector_a,stateVector_b
    real(8),          intent(out) :: dotsum

    ! Locals:
    integer          :: jstep,jlon,jlev,jlat,lon1,lon2,lat1,lat2
    integer          :: k1,k2

    if (.not.stateVector_a%allocated) then
      call utl_abort('gsv_dotProduct: gridStateVector_in not yet allocated')
    end if
    if (.not.statevector_b%allocated) then
      call utl_abort('gsv_dotProduct: gridStateVector_inout not yet allocated')
    end if

    lon1 = stateVector_a%myLonBeg
    lon2 = stateVector_a%myLonEnd
    lat1 = stateVector_a%myLatBeg
    lat2 = stateVector_a%myLatEnd
    k1 = stateVector_a%mykBeg
    k2 = stateVector_a%mykEnd

    dotsum = 0.0D0
    do jstep = 1, stateVector_a%numStep
      do jlev = k1,k2
        do jlat = lat1, lat2
          do jlon = lon1, lon2
            dotsum = dotsum + stateVector_a%gd_r8(jlon,jlat,jlev,jstep) * &
                              stateVector_b%gd_r8(jlon,jlat,jlev,jstep)
          end do 
        end do
      end do
    end do

    call mmpi_allreduce_sumreal8scalar(dotsum,'grid')

  end subroutine gsv_dotProduct

  !--------------------------------------------------------------------------
  ! gsv_field3d_hbilin
  !--------------------------------------------------------------------------
  subroutine gsv_field3d_hbilin(field,nlong,nlat,nlev,xlong,xlat,vlev, &
                                fieldout,nlongout,nlatout,nlevout,xlongout, &
                                xlatout,vlevout)
    !
    ! :Purpose: Horizontal bilinear interpolation from a 3D regular gridded field
    !           to another 3D regular gridded field.
    !
    !           This version can be used with fields that are not part of the
    !           background state, such as climatologies.
    !
    !           This version does not depend on gridstatevector data
    !           types/structures.
    !
    implicit none

    ! Arguments:
    real(8), intent(in) :: field(nlong,nlat,nlev) ! 3D field
    integer, intent(in) :: nlong ! number of latitudes
    integer, intent(in) :: nlat  ! number of longitudes
    integer, intent(in) :: nlev  ! number of vertical levels
    real(8), intent(in) :: xlong(nlong) ! longitudes (radians)
    real(8), intent(in) :: xlat(nlat)   ! latitudes (radians)
    real(8), intent(in) :: vlev(nlev)   ! vertical levels of input field (in pressure)
    real(8), intent(out) :: fieldout(nlongout,nlatout,nlevout) ! 3D field
    integer, intent(in) :: nlongout ! number or latitudes
    integer, intent(in) :: nlatout  ! number of target longitudes
    integer, intent(in) :: nlevout  ! Number of target vertical levels
    real(8), intent(in) :: xlongout(nlongout) ! target longitudes (radians) 
    real(8), intent(in) :: xlatout(nlatout)   ! target of target latitudes (radians)
    real(8), intent(in) :: vlevout(nlevout)   ! Target vertical levels (in pressure)

    ! Locals:
    real(8) :: lnvlev(nlev),lnvlevout(nlevout),plong2
    integer :: ilev,ilon,ilat,ilatp1,i,j,ilongout,ilatout
    logical :: same_vlev

    real(8) :: DLDX, DLDY, DLDP, DLW1, DLW2, DLW3, DLW4

    ! Check if vertical interpolation needed
    if (nlev /= nlevout) then
      same_vlev = .false.
    else
      if (any(abs(vlev-vlevout) > 0.01*vlev)) then
        same_vlev = .false.
      else
        same_vlev = .true.
      end if
    end if 
    
    ! Find near lat/long grid points
    
    do ilongout = 1, nlongout

      if (nlongout > 1) then    
        plong2 = xlongout(ilongout)
        if (plong2 < 0.0) plong2 = 2.D0*mpc_pi_r8 + plong2
        do ilon = 2,nlong
          if  (xlong(ilon-1) < xlong(ilon)) then
            if (plong2 >= xlong(ilon-1) .and. plong2 <= xlong(ilon)) exit
          else 
            ! Assumes this is a transition between 360 to 0 (if it exists). Skip over.
          end if
        end do
        ilon = ilon-1
      else
        ilon = 1
      end if
      
      do ilatout = 1, nlatout
         
        do ilat = 2, nlat
          if (xlatout(ilatout) <= xlat(ilat)) exit
        end do
        ilat = min(ilat-1,nlat)
        ilatp1 = min(ilat+1,nlat)
    
        ! Set lat/long interpolation weights
    
        if (nlongout > 1) then
          DLDX = (xlongout(ilongout) - xlong(ilon))/(xlong(ilon+1)-xlong(ilon))
          if (ilat < nlat) then
            DLDY = (xlatout(ilatout) - xlat(ilat))/(xlat(ilat+1)-xlat(ilat))
          else
            DLDY = (xlatout(ilatout) - xlat(ilat))/(xlat(ilat)-xlat(ilat-1))
          end if

          DLW1 = (1.d0-DLDX) * (1.d0-DLDY)
          DLW2 =       DLDX  * (1.d0-DLDY)
          DLW3 = (1.d0-DLDX) *       DLDY
          DLW4 =       DLDX  *       DLDY
        else
          if (ilat < nlat) then
            DLDY = (xlatout(ilatout) - xlat(ilat))/(xlat(ilat+1)-xlat(ilat))
          else
            DLDY = (xlatout(ilatout) - xlat(ilat))/(xlat(ilat)-xlat(ilat-1))
          end if

          DLW1 = (1.d0-DLDY)
          DLW3 = DLDY        
        end if
        
        ! Set vertical interpolation weights (assumes pressure vertical coordinate)
    
        if (.not.same_vlev) then
          lnvlevout(:) = log(vlevout(:))    
          lnvlev(:) = log(vlev(:))    
          
          ilev = 1
          do i = 1, nlevout
            do j = ilev, nlev          
               if (lnvlevout(i) < lnvlev(j)) exit    ! assumes both lnvlevout and lnvlev increase with increasing index value
            end do
            ilev = j-1
            if (ilev < 1) then
               ilev = 1
            else if (ilev >= nlev) then
               ilev = nlev-1
            end if
       
            DLDP = (lnvlev(ilev+1)-lnvlevout(i))/(lnvlev(ilev+1)-lnvlev(ilev))
            
            fieldout(ilongout,ilatout,i) = DLDP* (DLW1 * field(ilon,ilat,ilev) &
                           + DLW2 * field(ilon+1,ilat,ilev) &
                           + DLW3 * field(ilon,ilatp1,ilev) &
                           + DLW4 * field(ilon+1,ilatp1,ilev)) &
             + (1.d0-DLDP)* (DLW1 * field(ilon,ilat,ilev+1) &
                           + DLW2 * field(ilon+1,ilat,ilev+1) &
                           + DLW3 * field(ilon,ilatp1,ilev+1) &
                           + DLW4 * field(ilon+1,ilatp1,ilev+1))                               
          end do
        else if (nlongout > 1) then
          do ilev = 1, nlevout           
            fieldout(ilongout,ilatout,ilev) = DLW1 * field(ilon,ilat,ilev) &
                           + DLW2 * field(ilon+1,ilat,ilev) &
                           + DLW3 * field(ilon,ilatp1,ilev) &
                           + DLW4 * field(ilon+1,ilatp1,ilev)
          end do
        else 
          do ilev = 1, nlevout           
            fieldout(ilongout,ilatout,ilev) = DLW1 * field(ilon,ilat,ilev) &
                           + DLW3 * field(ilon,ilatp1,ilev) 
          end do
        end if 
      end do
    end do
        
  end subroutine gsv_field3d_hbilin   

  !--------------------------------------------------------------------------
  ! gsv_smoothHorizontal
  !--------------------------------------------------------------------------
  subroutine gsv_smoothHorizontal(stateVector_inout, horizontalScale, maskNegatives_opt, &
       varName_opt, binInteger_opt, binReal_opt, binRealThreshold_opt)
    !
    ! :Purpose: To apply a horizontal smoothing to all of the fields according
    !           to the specified horizontal length scale
    !
    implicit none

    ! Arguments:
    type(struct_gsv), target,   intent(inout) :: stateVector_inout
    real(8),                    intent(in)    :: horizontalScale
    logical, optional,          intent(in)    :: maskNegatives_opt
    character(len=*), optional, intent(in)    :: varName_opt
    real(4), optional, pointer, intent(in)    :: binInteger_opt(:,:,:)
    real(8), optional, pointer, intent(in)    :: binReal_opt(:,:)
    real(8), optional,          intent(in)    :: binRealThreshold_opt

    ! Locals:
    type(struct_gsv), pointer :: stateVector
    type(struct_gsv), target  :: stateVector_varsLevs

    integer :: latIndex, lonIndex, stepIndex, kIndex
    integer :: latIndex2, lonIndex2, maxDeltaIndex, minDeltaIndex, count
    integer :: latBeg, latEnd, lonBeg, lonEnd
    integer :: myBinInteger

    real(8), allocatable :: smoothedField(:,:)
    real(8) :: lat1_r8, lon1_r8, lat2_r8, lon2_r8, distance
    real(8) :: binRealThreshold, myBinReal

    real(4), pointer :: binInteger(:,:,:)
    real(8), pointer :: binReal(:,:)

    logical :: maskNegatives, binIntegerTest, binRealTest
    
    call utl_tmg_start(169, 'low-level--gsv_smoothHorizontal')

    if (horizontalScale <= 0.0d0) then
      call msg('gsv_smoothHorizontal', 'specified scale <= 0, returning')
      return
    end if

    if (present(binInteger_opt)) then
      binIntegerTest= .true.
      binInteger => binInteger_opt
    else
      binIntegerTest = .false.
    end if

    if (present(binReal_opt)) then
      binRealTest= .true.
      binReal => binReal_opt
      if (present(binRealThreshold_opt)) then
        binRealThreshold = binRealThreshold_opt
      else
        call utl_abort('gsv_smoothHorizontal: a binRealThreshold_opt value must be provided with binReal_opt')
      end if
    else
      binRealTest = .false.
    end if

    if (stateVector_inout%mpi_distribution /= 'VarsLevs' .and. stateVector_inout%mpi_local) then
    
      call gsv_allocate(statevector_varsLevs, statevector_inout%numStep, statevector_inout%hco, &
                        statevector_inout%vco, dataKind_opt=statevector_inout%dataKind,         &
                        mpi_local_opt=.true., mpi_distribution_opt='VarsLevs', &
                        allocHeight_opt=.false., allocPressure_opt=.false.)
      call gsv_transposeTilesToVarsLevs(statevector_inout, statevector_varsLevs)
      stateVector => stateVector_varsLevs
      
    else
    
      stateVector => stateVector_inout
      
    end if

    if (present(maskNegatives_opt)) then
      maskNegatives = maskNegatives_opt
    else
      maskNegatives = .false.
    end if

    allocate(smoothedField(stateVector%ni,stateVector%nj))

    ! figure out the maximum possible number of grid points to search
    if (stateVector%hco%dlat > 0.0d0) then
      maxDeltaIndex = ceiling(1.5d0 * horizontalScale / (ec_ra * max(stateVector%hco%dlat,stateVector%hco%dlon)))
      write(*,*) 'gsv_smoothHorizontal: maxDistance, maxDeltaIndex = ', horizontalScale / 1000.0d0, 'km', &
           maxDeltaIndex, max(stateVector%hco%dlat,stateVector%hco%dlon)
    else if(stateVector%hco%dlat == 0.0d0) then
      maxDeltaIndex = ceiling(1.5d0 * horizontalScale / statevector%hco%minGridSpacing)
      write(*,*) 'gsv_smoothHorizontal: maxDistance: ', horizontalScale / 1000.0d0, ' km,', &
           ' maxDeltaIndex: ', maxDeltaIndex, ', minGridSpacing: ', statevector%hco%minGridSpacing
    else
      call utl_abort('gsv_smoothHorizontal: cannot compute a value for maxDeltaIndex')
    end if

    ! apply a simple footprint operator type of averaging within specified radius
    do stepIndex = 1, stateVector%numStep
      do kIndex = stateVector%mykBeg, stateVector%mykEnd
      
        if (present(varName_opt)) then
          if (gsv_getVarNameFromK(stateVector,kIndex) /= trim(varName_opt)) cycle
        end if
	
        smoothedField(:,:) = 0.0d0
        do latIndex = 1, stateVector%nj
          do lonIndex = 1, stateVector%ni
	  
            lat1_r8 = stateVector%hco%lat2d_4(lonIndex,latIndex)
            lon1_r8 = stateVector%hco%lon2d_4(lonIndex,latIndex)
            count = 0
            latBeg = max(1, latIndex - maxDeltaIndex)
            latEnd = min(stateVector%nj, latIndex + maxDeltaIndex)
            lonBeg = max(1, lonIndex - maxDeltaIndex)
            lonEnd = min(stateVector%ni, lonIndex + maxDeltaIndex)
            if (binIntegerTest) then
              myBinInteger = int(binInteger(lonIndex,latIndex,1)) 
            end if
            
	    if (binRealTest) then
              myBinReal = binReal(lonIndex,latIndex)
            end if
	    
            do latIndex2 = latBeg, latEnd
              do lonIndex2 = lonBeg, lonEnd

                ! skip negative value if it should be masked
                if (maskNegatives .and. stateVector%dataKind == 8) then
                  if (stateVector%gd_r8(lonIndex2,latIndex2,kIndex,stepIndex) < 0.0d0) cycle
                else if (maskNegatives .and. statevector%dataKind == 4) then
                  if (stateVector%gd_r4(lonIndex2,latIndex2,kIndex,stepIndex) < 0.0) cycle
                end if

                ! skip value if it is beyond the specified distance
                lat2_r8 = stateVector%hco%lat2d_4(lonIndex2,latIndex2)
                lon2_r8 = stateVector%hco%lon2d_4(lonIndex2,latIndex2)
                distance = phf_calcDistanceFast(lat2_r8, lon2_r8, lat1_r8, lon1_r8)
                if (distance > horizontalScale) cycle

                ! skip value if it lie in a different bin
                if (binIntegerTest) then
                  if (int(binInteger(lonIndex2,latIndex2,1)) /=  myBinInteger .or. & 
                      int(binInteger(lonIndex2,latIndex2,1)) == -1) cycle
                end if
		
                if (binRealTest) then
                  if (abs(binReal(lonIndex2,latIndex2) - myBinReal) > binRealThreshold) cycle
                end if

                count = count + 1
                if (stateVector%dataKind == 8) then
                  smoothedField(lonIndex,latIndex) = smoothedField(lonIndex,latIndex) + &
                                                     stateVector%gd_r8(lonIndex2,latIndex2,kIndex,stepIndex)
                else
                  smoothedField(lonIndex,latIndex) = smoothedField(lonIndex,latIndex) + &
                                                     real(stateVector%gd_r4(lonIndex2,latIndex2,kIndex,stepIndex),8)
                end if
            
	      end do
            end do
            
	    if (count > 0) then
              smoothedField(lonIndex,latIndex) = smoothedField(lonIndex,latIndex) / real(count,8)
            else
              if (maskNegatives) smoothedField(lonIndex,latIndex) = mpc_missingValue_r8
            end if
	    
          end do
        end do
	
        if (stateVector%dataKind == 8) then
          stateVector%gd_r8(:,:,kIndex,stepIndex) = smoothedField(:,:)
        else
          stateVector%gd_r4(:,:,kIndex,stepIndex) = real(smoothedField(:,:),4)
        end if
      end do
    end do

    deallocate(smoothedField)

    if (stateVector_inout%mpi_distribution /= 'VarsLevs' .and. stateVector_inout%mpi_local) then
      call gsv_transposeVarsLevsToTiles(statevector_varsLevs, statevector_inout)
      call gsv_deallocate(statevector_varsLevs)
    end if
    
    call utl_tmg_stop(169)

  end subroutine gsv_smoothHorizontal

  !--------------------------------------------------------------------------
  ! gsv_getInfo
  !--------------------------------------------------------------------------
  ! DBGmad : TODO use a `gsv_str()` string representation?
  subroutine gsv_getInfo(stateVector, message)
    !:Purpose: Writes out grid state vector parameters
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in) :: stateVector
    character(len=*) :: message

    write(*,*)
    write(*,*) message
    write(*,*) '------------------- START -------------------'
    write(*,*) 'heightSfcPresent = ',stateVector%heightSfcPresent
    write(*,*) 'UVComponentPresent = ',stateVector%UVComponentPresent
    write(*,*) 'extraUVallocated = ',stateVector%extraUVallocated
    write(*,*) 'myUVkBeg = ',stateVector%myUVkBeg
    write(*,*) 'myUVkEnd = ',stateVector%myUVkEnd
    write(*,*) 'myUVkCount = ',stateVector%myUVkCount
    write(*,*) 'dataKind = ',stateVector%dataKind
    write(*,*) 'ni = ',stateVector%ni
    write(*,*) 'nj = ',stateVector%nj
    write(*,*) 'nk = ',stateVector%nk
    write(*,*) 'numStep = ',stateVector%numStep
    write(*,*) 'anltime = ',stateVector%anltime
    write(*,*) 'latPerPE = ',stateVector%latPerPE
    write(*,*) 'latPerPEmax = ',stateVector%latPerPEmax
    write(*,*) 'myLatBeg = ',stateVector%myLatBeg
    write(*,*) 'myLatEnd = ',stateVector%myLatEnd
    write(*,*) 'lonPerPE = ',stateVector%lonPerPE
    write(*,*) 'lonPerPEmax = ',stateVector%lonPerPEmax
    write(*,*) 'myLonBeg = ',stateVector%myLonBeg
    write(*,*) 'myLonEnd = ',stateVector%myLonEnd
    write(*,*) 'mykCount = ',stateVector%mykCount
    write(*,*) 'mykBeg = ',stateVector%mykBeg
    write(*,*) 'mykEnd = ',stateVector%mykEnd
    if (associated(stateVector%allLatBeg)) write(*,*) 'allLatBeg = ',stateVector%allLatBeg
    if (associated(stateVector%allLatEnd)) write(*,*) 'allLatEnd = ',stateVector%allLatEnd
    if (associated(stateVector%allLatPerPE)) write(*,*) 'allLatPerPE = ',stateVector%allLatPerPE
    if (associated(stateVector%allLonBeg)) write(*,*) 'allLonBeg = ',stateVector%allLonBeg
    if (associated(stateVector%allLonEnd)) write(*,*) 'allLonEnd = ',stateVector%allLonEnd
    if (associated(stateVector%allLonPerPE)) write(*,*) 'allLonPerPE = ',stateVector%allLonPerPE
    if (associated(stateVector%allkCount)) write(*,*) 'allkCount = ',stateVector%allkCount
    if (associated(stateVector%allkBeg)) write(*,*) 'allkBeg = ',stateVector%allkBeg
    if (associated(stateVector%allkEnd)) write(*,*) 'allkEnd = ',stateVector%mykEnd
    if (associated(stateVector%allUVkCount)) write(*,*) 'allUVkCount = ',stateVector%allUVkCount
    if (associated(stateVector%allUVkBeg)) write(*,*) 'allUVkBeg = ',stateVector%allUVkBeg
    if (associated(stateVector%allUVkEnd)) write(*,*) 'allUVkEnd = ',stateVector%myUVkEnd
    if (associated(stateVector%dateStampList)) write(*,*) 'dateStampList = ',stateVector%dateStampList
    if (associated(stateVector%dateStamp3d)) write(*,*) 'dateStamp3d = ',stateVector%dateStamp3d
    if (associated(stateVector%dateOriginList)) write(*,*) 'dateOriginList = ',stateVector%dateOriginList
    if (associated(stateVector%npasList)) write(*,*) 'npasList = ',stateVector%npasList
    if (associated(stateVector%ip2List)) write(*,*) 'ip2List = ',stateVector%ip2List
    write(*,*) 'deet = ',stateVector%deet
    write(*,*) 'etiket = ',stateVector%etiket
    write(*,*) 'allocated = ',stateVector%allocated
    if (associated(stateVector%varOffset)) write(*,*) 'varOffset = ',stateVector%varOffset
    if (associated(stateVector%varNumLev)) write(*,*) 'varNumLev = ',stateVector%varNumLev
    write(*,*) 'mpi_local = ',stateVector%mpi_local
    write(*,*) 'mpi_distribution = ',stateVector%mpi_distribution
    write(*,*) 'horizSubSample = ',stateVector%horizSubSample
    write(*,*) 'varExistList = ',stateVector%varExistList
    write(*,*) 'hInterpolateDegree = ',stateVector%hInterpolateDegree
    write(*,*) 'hExtrapolateDegree = ',stateVector%hExtrapolateDegree
    write(*,*) 'addHeightSfcOffset = ',stateVector%addHeightSfcOffset
    write(*,*) '-------------------- END --------------------'

  end subroutine gsv_getInfo

  !--------------------------------------------------------------------------
  ! gsv_applyMaskLAM
  !--------------------------------------------------------------------------
  subroutine gsv_applyMaskLAM(statevector_inout, maskLAM)
    !:Purpose: To apply a mask to a state vector for LAM grid
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout)  :: statevector_inout
    type(struct_gsv), intent(in)     :: maskLAM

    ! Locals
    real(4), pointer :: increment_r4(:,:,:,:)
    real(8), pointer :: increment_r8(:,:,:,:)
    real(pre_incrReal), pointer :: analIncMask(:,:,:)
    integer :: latIndex, kIndex, lonIndex, stepIndex

    call msg('gsv_applyMaskLAM','START', verb_opt=2)
    
    call gsv_getField(maskLAM,analIncMask)

    if (statevector_inout%dataKind == 4) then
      ! apply mask using increment_r4
      call gsv_getField(statevector_inout,increment_r4)
      do stepIndex = 1, statevector_inout%numStep
        !$OMP PARALLEL DO PRIVATE (latIndex,kIndex,lonIndex)
        do kIndex = 1, statevector_inout%nk
          do latIndex =  statevector_inout%myLatBeg,  statevector_inout%myLatEnd
            do lonIndex =  statevector_inout%myLonBeg,  statevector_inout%myLonEnd
              increment_r4(lonIndex,latIndex,kIndex,stepIndex) =      &
                   increment_r4(lonIndex,latIndex,kIndex,stepIndex) * &
                   analIncMask(lonIndex,latIndex,1)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do
    else
      ! apply mask using increment_r8
      call gsv_getField(statevector_inout,increment_r8)
      do stepIndex = 1, statevector_inout%numStep
        !$OMP PARALLEL DO PRIVATE (latIndex,kIndex,lonIndex)
        do kIndex = 1, statevector_inout%nk
          do latIndex =  statevector_inout%myLatBeg,  statevector_inout%myLatEnd
            do lonIndex =  statevector_inout%myLonBeg,  statevector_inout%myLonEnd
              increment_r8(lonIndex,latIndex,kIndex,stepIndex) =      &
                   increment_r8(lonIndex,latIndex,kIndex,stepIndex) * &
                   analIncMask(lonIndex,latIndex,1)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do
    end if

    call msg('gsv_applyMaskLAM','END', verb_opt=2)

  end subroutine gsv_applyMaskLAM

  !--------------------------------------------------------------------------
  ! gsv_containsNonZeroValues
  !--------------------------------------------------------------------------
  function gsv_containsNonZeroValues(stateVector) result(stateVectorHasNonZeroValue)
    !:Purpose: To check if stateVector has any non-zero value.
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(in) :: stateVector
    logical                      :: stateVectorHasNonZeroValue

    ! Locals
    real(4), pointer             :: field_r4_ptr(:,:,:,:)
    real(8), pointer             :: field_r8_ptr(:,:,:,:)
    logical                      :: allZero, allZero_mpiglobal
    integer                      :: ierr

    if (.not. stateVector%allocated) then
      stateVectorHasNonZeroValue = .false.
      return
    end if

    if (stateVector%dataKind == 4) then
      call gsv_getField(stateVector,field_r4_ptr)
      allZero = (maxval(abs(field_r4_ptr(:,:,:,:))) == 0.0)
    else if (stateVector%dataKind == 8) then
      call gsv_getField(stateVector,field_r8_ptr)
      allZero = (maxval(abs(field_r8_ptr(:,:,:,:))) == 0.0D0)
    end if

    call rpn_comm_allReduce(allZero,allZero_mpiglobal,1,'mpi_logical','mpi_land','GRID',ierr)
    stateVectorHasNonZeroValue = .not. allZero_mpiglobal

  end function gsv_containsNonZeroValues

end module gridStateVector_mod

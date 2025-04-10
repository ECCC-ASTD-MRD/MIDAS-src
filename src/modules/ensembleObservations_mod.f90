
MODULE ensembleObservations_mod
  ! MODULE ensembleObservations_mod (prefix='eob' category='6. High-level data objects')
  !
  !:Purpose:  Store and manipulate ensemble of quanitites in observation space.
  !           This module uses the kdtree2 module for efficiently finding the
  !           nearest observations within the local volume.
  !
  use kdTree2_mod
  use message_mod
  use columnData_mod
  use tovs_mod
  use rttov_types, only: rttov_transmission, rttov_profile
  use parkind1, only: jpim, jprb
  use midasMpi_mod
  use oceanMask_mod
  use verticalCoord_mod
  use obsSpaceData_mod
  use randomNumber_mod
  use mathPhysConstants_mod
  use utilities_mod
  use earthConstants_mod
  use bufr_mod
  use codePrecision_mod
  use codtyp_mod
  use obsfamilylist_mod
  use varnamelist_mod
  use localizationFunction_mod
  use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
  use mpi   ! This is the standard mpi library

  implicit none
  save
  private

  ! Public variables
  logical, public, protected :: eob_simObsAssim

  ! Public derived type
  public :: struct_eob

  ! Public procedures
  public :: eob_init, eob_allocate, eob_deallocate, eob_allGather, eob_getLocalBodyIndices
  public :: eob_setYb, eob_setYa, eob_setDeterYb, eob_setLatLonObsCod, eob_setObsErrInv
  public :: eob_setHPHT, eob_calcAndRemoveMeanYb, eob_setVertLocation, eob_setAssFlag, eob_copy, eob_zero
  public :: eob_calcRandPert, eob_setSigiSigo, eob_setTypeVertCoord, eob_setSimObsVal
  public :: eob_backgroundCheck, eob_huberNorm, eob_rejectRadNearSfc, eob_setMeanOMP
  public :: eob_removeObsNearLand, eob_readFromFiles, eob_writeToFiles

  integer, parameter   :: maxNumLocalObsSearch = 500000
  logical              :: psvObsAssim
  integer              :: numSimObsFam
  integer              :: numPsvObsFam
  integer              :: numSimCodTyp(ofl_numFamily), numPsvCodTyp(ofl_numFamily)
  integer              :: numSimVarNum(vnl_numvarmax), numPsvVarNum(vnl_numvarmax)
  integer, allocatable :: simCodTyp(:,:), psvCodTyp(:,:)
  integer, allocatable :: simVarNum(:,:), psvVarNum(:,:)

  type struct_eob
    logical                       :: allocated      = .false.
    logical                       :: mpiGlobal      = .false.
    logical                       :: meanRemoved    = .false.
    integer                       :: numMembers       ! number of ensemble members
    integer                       :: firstMember      ! starting value of member index
    integer                       :: numObs           ! number of observations
    integer                       :: fileMemberIndex1 = 1 ! first member number in ensemble set
    character(len=20)             :: typeVertCoord = 'undefined' ! 'logPressure' or 'depth'
    type(struct_obs), pointer     :: obsSpaceData     ! pointer to obsSpaceData object
    real(8), allocatable          :: lat(:), lon(:)   ! lat/lon of observation
    real(8), allocatable          :: vertLocation(:)  ! in ln(pres) or meters, used for localization
    real(8), allocatable          :: obsErrInv(:)     ! inverse of obs error variances
    real(8), allocatable          :: obsErrInv_sim(:) ! like obsErrInv, used when simulating observations
    integer                       :: Yb_window = -1   ! handle to shared data "window" for Yb_r4
    type(c_ptr)                   :: Yb_baseptr                 ! pointer used to allocated shared memory
    real(4), pointer              :: Yb_r4(:,:) => null()       ! background ensemble perturbation in obs space
    integer                       :: Ya_window = -1             ! handle to shared data "window" for Ya_r4
    type(c_ptr)                   :: Ya_baseptr                 ! pointer used to allocated shared memory
    real(4), pointer              :: Ya_r4(:,:) => null()       ! analysis ensemble perturbation in obs space
    type(c_ptr)                   :: randPert_baseptr           ! pointer used to allocated shared memory
    real(4), pointer              :: randPert_r4(:,:) => null() ! unbiased random perturbations with covariance equal to R
    integer                       :: randPert_window = -1       ! handle to shared data "window" for Ya_r4
    real(8), allocatable          :: meanYb(:)        ! ensemble mean background state in obs space
    real(8), allocatable          :: deterYb(:)       ! deterministic background state in obs space
    real(8), allocatable          :: obsValue(:)      ! the observed value
    integer, allocatable          :: assFlag(:)       ! assimilation flag
    integer, allocatable          :: codTyp(:)        ! observation code type
  end type struct_eob

  type(kdtree2), pointer :: tree => null()

  type struct_eobOutput
    real(8) :: distMax(codtyp_maxNumber)
    real(8) :: distMean(codtyp_maxNumber)
    real(8) :: locFun(codtyp_maxNumber)
    real(8) :: obsErr(codtyp_maxNumber)
    real(8) :: ensSpread(codtyp_maxNumber)
    real(8) :: trace(codtyp_maxNumber)
    real(8) :: vertTop(codtyp_maxNumber)
    real(8) :: vertBottom(codtyp_maxNumber)
  end type struct_eobOutput

  ! namelist variables
  character(len=2)  :: simObsFamily(ofl_numFamily) ! observation families for simulation
  character(len=2)  :: psvObsFamily(ofl_numFamily) ! observation families for passive assimilation
  character(len=codtyp_name_length) :: simCodTypName(ofl_numFamily,codtyp_maxNumber) ! codtyp names for sim. obs families
  character(len=codtyp_name_length) :: psvCodTypName(ofl_numFamily,codtyp_maxNumber) ! codtyp names for psv. obs families
  character(len=4) :: simVarName(ofl_numFamily,vnl_numvarmax) ! varName(s) for sim. obs families
  character(len=4) :: psvVarName(ofl_numFamily,vnl_numvarmax) ! varName(s) for psv. obs families
  namelist /NAMENSOBS/simObsFamily, psvObsFamily, simCodTypName, psvCodTypName, simVarName, psvVarName


CONTAINS

  !--------------------------------------------------------------------------
  ! eob_init
  !--------------------------------------------------------------------------
  subroutine eob_init()
    !
    !: Purpose: This subroutine reads the namelist section NAMENSOBS for this module.
    !
    implicit none

    ! Locals:
    integer :: ierr, obsfamIndex, codtypIndex, varnumIndex
    logical, save :: eob_initialized = .false.

    if (eob_initialized) return
    write(*,*) 'eob_init: starting'
    eob_initialized = .true.

    ! default values for namelist variables
    simObsFamily(:)    = ''
    simCodTypName(:,:) = ''
    simVarName(:,:)    = ''
    psvObsFamily(:)    = ''
    psvCodTypName(:,:) = ''
    psvVarName(:,:)    = ''

    ! for tracking the number of non-empty chars in namelist variable arrays;
    ! these are used in loops in various subroutines
    numSimObsFam = 0
    numPsvObsFam = 0
    numSimCodTyp(:) = 0
    numPsvCodTyp(:) = 0
    numSimVarNum(:) = 0
    numPsvVarNum(:) = 0

    ! read namelist
    if (utl_isNamelistPresent('namensobs','./flnml')) then
      call utl_tmg_start(181,'low-level--readNML')
      read(utl_flnml, nml=namensobs, iostat=ierr)
      if (ierr /= 0) call utl_abort('eob_init: Error reading namelist namensobs')
      call utl_tmg_stop(181)
      do obsfamIndex = 1, ofl_numFamily
        if (trim(simObsFamily(obsfamIndex)) /= '') then
          numSimObsFam = numSimObsFam + 1
        end if
        if (trim(psvObsFamily(obsfamIndex)) /= '') then
          numPsvObsFam = numPsvObsFam + 1
        end if
      end do

      do obsfamIndex = 1, ofl_numFamily
        ! Simulation functionality section
        if (trim(simObsFamily(obsfamIndex)) /= '') then
          ! simulated observation family specified for current
          ! obsfamIndex; check to see if any codtyp names specified
          do codtypIndex = 1, codtyp_maxNumber
            if (trim(simCodTypName(obsfamIndex,codtypIndex)) /= '') then
              numSimCodTyp(obsfamIndex) = numSimCodTyp(obsfamIndex) + 1
              if (.not. allocated(simCodTyp)) then
                allocate(simCodTyp(numSimObsFam,codtyp_maxNumber))
                simCodTyp(:,:) = -999
              end if
              ! store CodTyp for simulated obs family
              simCodTyp(obsfamIndex,codtypIndex) = codtyp_get_codtyp(simCodTypName(obsfamIndex,codtypIndex))
            end if
          end do ! codtypIndex
          ! also check to see if any varnames specified
          do varnumIndex = 1, vnl_numvarmax
            if (trim(simVarName(obsfamIndex,varnumIndex)) /= '') then
              numSimVarNum(obsfamIndex) = numSimVarNum(obsfamIndex) + 1
              if (.not. allocated(simVarNum)) then
                allocate(simVarNum(numSimObsFam,vnl_numvarmax))
                simVarNum(:,:) = -999
              end if
              ! store VarNum for simulated obs family
              simVarNum(obsfamIndex,varnumIndex) = vnl_varnumFromVarName(simVarName(obsfamIndex,varnumIndex))
            end if
          end do ! varnumIndex
        end if ! simObsFamily

        ! Passive functionality section
        if (trim(psvObsFamily(obsfamIndex)) /= '')  then
          ! passive observation family specified for current
          ! obsfamIndex; check to see if any codtyp names specified
          do codtypIndex = 1, codtyp_maxNumber
            if (trim(psvCodTypName(obsfamIndex,codtypIndex)) /= '') then
              numPsvCodTyp(obsfamIndex) = numPsvCodTyp(obsfamIndex) + 1
              if (.not. allocated(psvCodTyp)) then
                allocate(psvCodTyp(numPsvObsFam,codtyp_maxNumber))
                psvCodTyp(:,:) = -999
              end if
              ! store CodTyp for passive obs family
              psvCodTyp(obsfamIndex,codtypIndex) = codtyp_get_codtyp(psvCodTypName(obsfamIndex,codtypIndex))
            end if
          end do ! codtypIndex
          ! also check to see if any varnames specified
          do varnumIndex = 1, vnl_numvarmax
            if (trim(psvVarName(obsfamIndex,varnumIndex)) /= '') then
              numPsvVarNum(obsfamIndex) = numPsvVarNum(obsfamIndex) + 1
              if (.not. allocated(PsvVarNum)) then
                allocate(psvVarNum(numPsvObsFam,vnl_numvarmax))
                psvVarNum(:,:) = -999
              end if
              ! store VarNum for passive obs family
              psvVarNum(obsfamIndex,varnumIndex) = vnl_varnumFromVarName(psvVarName(obsfamIndex,varnumIndex))
            end if
          end do ! varnumIndex
        end if ! psvObsFamily

      end do ! obsFamIndex
    else
      write(*,*)
      write(*,*) 'eob_init: namensobs is missing in the namelist. The default value will be taken.'
    end if

    eob_simObsAssim = numSimObsFam > 0
    psvObsAssim = numPsvObsFam > 0

    write(*,*) 'eob_init: eob_simObsAssim = ', eob_simObsAssim
    write(*,*) 'eob_init: psvObsAssim     = ', psvObsAssim

  end subroutine eob_init

  !--------------------------------------------------------------------------
  ! eob_allocate
  !--------------------------------------------------------------------------
  subroutine eob_allocate(ensObs, numMembers, numObs, obsSpaceData, &
                          fileMemberIndex1_opt, allocYb_opt, includeMember0_opt)
    !
    ! :Purpose: Allocate an ensObs object
    !
    implicit none

    ! Arguments:
    type(struct_eob)        , intent(inout) :: ensObs               ! eob object to be allocated
    integer                 , intent(in)    :: numMembers           ! number of ensemble members
    integer                 , intent(in)    :: numObs               ! number of observations
    type(struct_obs), target, intent(in)    :: obsSpaceData         ! obs object
    integer, optional       , intent(in)    :: fileMemberIndex1_opt ! first member index of this set (default is 1)
    logical, optional       , intent(in)    :: allocYb_opt          ! request that Yb_r4 is allocated (default is true)
    logical, optional       , intent(in)    :: includeMember0_opt   ! choose to start at member 0 instead of 1

    ! Locals:
    logical                                 :: allocYb, includeMember0

    if (ensObs%allocated) then
      write(*,*) 'eob_allocate: this object is already allocated, deallocating first.'
      call eob_deallocate(ensObs)
    end if

    if (present(fileMemberIndex1_opt)) then
      ensObs%fileMemberIndex1 = fileMemberIndex1_opt
    end if

    if (present(allocYb_opt)) then
      allocYb = allocYb_opt
    else
      allocYb = .true.
    end if

    if (present(includeMember0_opt)) then
      includeMember0 = includeMember0_opt
    else
      includeMember0 = .false.
    end if

    call eob_init()

    ensObs%obsSpaceData  => obsSpaceData
    ensObs%numMembers    = numMembers
    ensObs%numObs        = numObs
    if (includeMember0) then
      write(*,*) 'eob_allocate: NOTE: memberIndex starts from 0, not 1'
      ensObs%firstMember = 0
    else
      ensObs%firstMember = 1
    end if

    allocate(ensObs%lat(ensObs%numObs))
    allocate(ensObs%lon(ensObs%numObs))
    allocate(ensObs%vertLocation(ensObs%numObs))
    allocate(ensObs%obsValue(ensObs%numObs))
    allocate(ensObs%obsErrInv(ensObs%numObs))
    if (eob_simObsAssim) allocate(ensObs%obsErrInv_sim(ensObs%numObs))
    if (allocYb) allocate(ensObs%Yb_r4(ensObs%firstMember:ensObs%numMembers,ensObs%numObs))
    allocate(ensObs%meanYb(ensObs%numObs))
    allocate(ensObs%deterYb(ensObs%numObs))
    allocate(ensObs%assFlag(ensObs%numObs))
    allocate(ensObs%codTyp(ensObs%numObs))

    ensObs%allocated = .true.

    call msg_memUsage('eob_allocate')

  end subroutine eob_allocate

  !--------------------------------------------------------------------------
  ! eob_deallocate
  !--------------------------------------------------------------------------
  subroutine eob_deallocate(ensObs)
    implicit none

    ! Arguments:
    type(struct_eob), intent(inout) :: ensObs ! eob object to be deallocated

    ! Locals:
    integer :: ierr

    if (.not. ensObs%allocated) return

    deallocate(ensObs%lat)
    deallocate(ensObs%lon)
    deallocate(ensObs%vertLocation)
    deallocate(ensObs%obsValue)
    deallocate(ensObs%obsErrInv)
    if (allocated(ensObs%obsErrInv_sim)) deallocate(ensObs%obsErrInv_sim)
    deallocate(ensObs%meanYb)
    deallocate(ensObs%deterYb)
    deallocate(ensObs%assFlag)
    deallocate(ensObs%codTyp)

    if (ensObs%Yb_window /= -1) then
      ! using shared memory
      call mpi_win_free(ensObs%Yb_window, ierr)
      ensObs%Yb_window = -1
    else if (associated(ensObs%Yb_r4)) then
      deallocate(ensObs%Yb_r4)
    end if
    nullify(ensObs%Yb_r4)

    if (ensObs%Ya_window /= -1) then
      ! using shared memory
      call mpi_win_free(ensObs%Ya_window, ierr)
      ensObs%Ya_window = -1
    else if (associated(ensObs%Ya_r4)) then
      deallocate(ensObs%Ya_r4)
    end if
    nullify(ensObs%Ya_r4)

    if (ensObs%randPert_window /= -1) then
      ! using shared memory
      call mpi_win_free(ensObs%randPert_window, ierr)
      ensObs%randPert_window = -1
    else if (associated(ensObs%randPert_r4)) then
      deallocate(ensObs%randPert_r4)
    end if
    nullify(ensObs%randPert_r4)

    ensObs%allocated = .false.

  end subroutine eob_deallocate

  !--------------------------------------------------------------------------
  ! eob_zero
  !--------------------------------------------------------------------------
  subroutine eob_zero(ensObs)
    !
    ! :Purpose: Initialize an ensObs object to zero
    !
    implicit none

    ! Arguments:
    type(struct_eob), intent(inout) :: ensObs ! eob object to be modified

    if ( .not.ensObs%allocated ) then
      call utl_abort('eob_zero: this object is not allocated')
    end if

    ensObs%lat(:)           = 0.0d0
    ensObs%lon(:)           = 0.0d0
    ensObs%vertLocation(:)  = 0.0d0
    ensObs%obsValue(:)      = 0.0d0
    ensObs%obsErrInv(:)     = 0.0d0
    if (allocated(ensObs%obsErrInv_sim)) ensObs%obsErrInv_sim(:) = 0.0
    ensObs%Yb_r4(:,:)       = 0.0
    if (associated(ensObs%Ya_r4)) ensObs%Ya_r4(:,:) = 0.0
    if (associated(ensObs%randPert_r4)) ensObs%randPert_r4(:,:) = 0.0
    ensObs%meanYb(:)        = 0.0d0
    ensObs%deterYb(:)       = 0.0d0
    ensObs%assFlag(:)       = 0
    ensObs%codTyp(:)        = 0

    call msg_memUsage('eob_zero')

  end subroutine eob_zero

  !--------------------------------------------------------------------------
  ! eob_setTypeVertCoord
  !--------------------------------------------------------------------------
  subroutine eob_setTypeVertCoord(ensObs, typeVertCoord)
    !
    ! :Purpose: Set the type of vertical coordinate ('logPressure' or 'depth').
    !
    implicit none

    ! Arguments:
    type(struct_eob), intent(inout) :: ensObs        ! eob object to be modified
    character(len=*), intent(in)    :: typeVertCoord ! type of vertical coordinate

    if ( trim(typeVertCoord) /= 'logPressure' .and. &
         trim(typeVertCoord) /= 'depth' ) then
      write(*,*) 'eob_setTypeVertCoord: typeVertCoord = ', trim(typeVertCoord)
      call utl_abort('eob_setTypeVertCoord: Unknown type of vertical coordinate')
    end if

    ensObs%typeVertCoord = typeVertCoord

  end subroutine eob_setTypeVertCoord

  !--------------------------------------------------------------------------
  ! eob_clean (private routine)
  !--------------------------------------------------------------------------
  subroutine eob_clean(ensObs,ensObsClean)
    !
    ! :Purpose: Remove all obs from the ensObs object that are not
    !           flagged for assimilation. Put the cleaned result in the
    !           locally created output object.
    !
    implicit none

    ! Arguments:
    type(struct_eob), intent(inout) :: ensObs      ! input eob object to be cleaned
    type(struct_eob), intent(out)   :: ensObsClean ! resulting cleaned eob object

    ! Locals:
    integer :: obsIndex, obsCleanIndex, numObsClean

    if (ensObs%mpiglobal) then
      call utl_abort('eob_clean: Cannot be called with an mpiglobal object')
    end if

    if (ensObs%firstMember==0) then
      call utl_abort('eob_clean: Cannot be called when firstMember is 0')
    end if

    call eob_setAssFlag(ensObs)

    numObsClean = 0
    do obsIndex = 1, ensObs%numObs
      if (ensObs%assFlag(obsIndex) == 1) numObsClean = numObsClean + 1
    end do

    write(*,*) 'eob_clean: reducing numObs from ', ensObs%numObs, ' to ', numObsClean
    call eob_allocate(ensObsClean, ensObs%numMembers, numObsClean, ensObs%obsSpaceData)
    if (associated(ensObs%Ya_r4)) then
      allocate(ensObsClean%Ya_r4(ensObs%numMembers,numObsClean))
    end if
    if (associated(ensObs%randPert_r4)) then
      allocate(ensObsClean%randPert_r4(ensObs%numMembers,numObsClean))
    end if

    obsCleanIndex = 0
    do obsIndex = 1, ensObs%numObs
      if (ensObs%assFlag(obsIndex) == 1) then
        obsCleanIndex = obsCleanIndex + 1
        ensObsClean%lat(obsCleanIndex)           = ensObs%lat(obsIndex)
        ensObsClean%lon(obsCleanIndex)           = ensObs%lon(obsIndex)
        ensObsClean%vertLocation(obsCleanIndex)  = ensObs%vertLocation(obsIndex)
        ensObsClean%obsErrInv(obsCleanIndex)     = ensObs%obsErrInv(obsIndex)
        if (allocated(ensObs%obsErrInv_sim)) then
          ensObsClean%obsErrInv_sim(obsCleanIndex) = ensObs%obsErrInv_sim(obsIndex)
        end if
        ensObsClean%Yb_r4(:,obsCleanIndex)       = ensObs%Yb_r4(:,obsIndex)
        if (associated(ensObs%Ya_r4)) then
          ensObsClean%Ya_r4(:,obsCleanIndex) = ensObs%Ya_r4(:,obsIndex)
        end if
        if (associated(ensObs%randPert_r4)) then
          ensObsClean%randPert_r4(:,obsCleanIndex) = ensObs%randPert_r4(:,obsIndex)
        end if
        ensObsClean%meanYb(obsCleanIndex)        = ensObs%meanYb(obsIndex)
        ensObsClean%deterYb(obsCleanIndex)       = ensObs%deterYb(obsIndex)
        ensObsClean%obsValue(obsCleanIndex)      = ensObs%obsValue(obsIndex)
        ensObsClean%assFlag(obsCleanIndex)       = ensObs%assFlag(obsIndex)
        ensObsClean%codTyp(obsCleanIndex)        = ensObs%codTyp(obsIndex)
      end if
    end do

  end subroutine eob_clean

  !--------------------------------------------------------------------------
  ! eob_copy
  !--------------------------------------------------------------------------
  subroutine eob_copy(ensObsIn,ensObsOut)
    implicit none

    ! Arguments:
    type(struct_eob), intent(in)    :: ensObsIn   ! input eob object
    type(struct_eob), intent(inout) :: ensObsOut  ! output eob object

    if (ensObsIn%mpiglobal) then
      call utl_abort('eob_copy: Cannot be called with an mpiglobal object')
    end if
    if (ensObsIn%firstMember==0) then
      call utl_abort('eob_copy: Cannot be called when firstMember is 0')
    end if

    ensObsOut%lat(:)           = ensObsIn%lat(:)
    ensObsOut%lon(:)           = ensObsIn%lon(:)
    ensObsOut%vertLocation(:)  = ensObsIn%vertLocation(:)
    ensObsOut%obsErrInv(:)     = ensObsIn%obsErrInv(:)
    if (allocated(ensObsIn%obsErrInv_sim)) then
      ensObsOut%obsErrInv(:) = ensObsIn%obsErrInv_sim(:)
    end if
    ensObsOut%Yb_r4(:,:)       = ensObsIn%Yb_r4(:,:)
    if (associated(ensObsIn%Ya_r4)) then
      allocate(ensObsOut%Ya_r4(ensObsIn%numMembers,ensObsIn%numObs))
      ensObsOut%Ya_r4(:,:) = ensObsIn%Ya_r4(:,:)
    end if
    if (associated(ensObsIn%randPert_r4)) then
      allocate(ensObsOut%randPert_r4(ensObsIn%numMembers,ensObsIn%numObs))
      ensObsOut%randPert_r4(:,:) = ensObsIn%randPert_r4(:,:)
    end if
    ensObsOut%meanYb(:)        = ensObsIn%meanYb(:)
    ensObsOut%deterYb(:)       = ensObsIn%deterYb(:)
    ensObsOut%obsValue(:)      = ensObsIn%obsValue(:)
    ensObsOut%assFlag(:)       = ensObsIn%assFlag(:)
    ensObsOut%codTyp(:)        = ensObsIn%codTyp(:)
    ensObsOut%typeVertCoord    = ensObsIn%typeVertCoord

  end subroutine eob_copy

  !--------------------------------------------------------------------------
  ! eob_allGather
  !--------------------------------------------------------------------------
  subroutine eob_allGather(ensObs,ensObs_mpiglobal)
    !
    ! :Purpose: Collect obs information distributed over all mpi tasks and
    !           make it available on all mpi tasks. The output ensObs object
    !           will be allocated within this subroutine.
    !
    implicit none

    ! Arguments:
    type(struct_eob), intent(inout)  :: ensObs           ! input mpi local eob object
    type(struct_eob), intent(out)    :: ensObs_mpiglobal ! output mpi global eob object (all obs on each mpi task)

    ! Locals:
    type(struct_eob) :: ensObsClean
    integer :: ierr, nsize, procIndex, memberIndex, numObs_mpiglobal
    integer :: allNumObs(mmpi_nprocs), displs(mmpi_nprocs)
    integer :: dispUnit, arrayShape(2)
    integer, allocatable :: requestId(:)
    integer(kind=mpi_address_kind) :: windowSize
    real(4), allocatable :: tempBuffer(:,:)

    write(*,*) 'eob_allGather: starting'
    call msg_memUsage('eob_allGather')

    if (ensObs%mpiglobal) then
      call utl_abort('eob_allGather: Input eob object is already mpiglobal')
    end if
    if (ensObs%firstMember == 0) then
      call utl_abort('eob_allGather: Cannot be called when firstMember is 0')
    end if

    call utl_tmg_start(10,'--Observations')
    call utl_tmg_start(24,'----Eob_AllGather')

    ! refresh assimilation flag and then clean ensObs before communicating and writing
    call eob_setAssFlag(ensObs)
    call eob_clean(ensObs,ensObsClean)

    call rpn_comm_allgather(ensObsClean%numObs, 1, 'mpi_integer',  &
                            allNumObs, 1, 'mpi_integer', &
                            'GRID', ierr)
    numObs_mpiglobal = sum(allNumObs(:))

    if (ensObs_mpiglobal%allocated) then
      call utl_abort('eob_allGather: output ensObs object must not be already allocated')
    end if
    call eob_allocate(ensObs_mpiglobal, ensObsClean%numMembers, numObs_mpiglobal, ensObsClean%obsSpaceData, &
                      fileMemberIndex1_opt=ensObs%fileMemberIndex1, allocYb_opt=.false.)
    ensObs_mpiglobal%mpiglobal = .true.

    ! Allocate large mpi global arrays that will be put in shared memory
    arrayShape=(/ ensObsClean%numMembers,numObs_mpiglobal /)
    if (mmpi_myidHost == 0) then
      windowSize = int(product(arrayShape),MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
    else
      windowSize = 0_MPI_ADDRESS_KIND
    end if
    dispUnit = 1
    write(*,*) 'eob_allGather: allocate shared memory arrays with size = ', windowSize

    ! Shared memory for Yb_r4
    call mpi_win_allocate_shared(windowSize, dispUnit, MPI_INFO_NULL, &
                                 mmpi_mpicomm_shared, ensObs_mpiglobal%Yb_baseptr, &
                                 ensObs_mpiglobal%Yb_window, ierr)
    if (mmpi_myidHost /= 0) then
      call mpi_win_shared_query(ensObs_mpiglobal%Yb_window, 0, windowSize, dispUnit, &
                                ensObs_mpiglobal%Yb_baseptr, ierr)
    end if
    call c_f_pointer(ensObs_mpiglobal%Yb_baseptr, ensObs_mpiglobal%Yb_r4, arrayShape)

    ! Share memory for Ya_r4
    if (associated(ensObsClean%Ya_r4)) then
      call mpi_win_allocate_shared(windowSize, dispUnit, MPI_INFO_NULL, &
                                   mmpi_mpicomm_shared, ensObs_mpiglobal%Ya_baseptr, &
                                   ensObs_mpiglobal%Ya_window, ierr)
      if (mmpi_myidHost /= 0) then
        call mpi_win_shared_query(ensObs_mpiglobal%Ya_window, 0, windowSize, dispUnit, &
                                  ensObs_mpiglobal%Ya_baseptr, ierr)
      end if
      call c_f_pointer(ensObs_mpiglobal%Ya_baseptr, ensObs_mpiglobal%Ya_r4, arrayShape)
    end if

    ! Share memory for randPert_r4
    if (associated(ensObsClean%randPert_r4)) then
      call mpi_win_allocate_shared(windowSize, dispUnit, MPI_INFO_NULL, &
                                   mmpi_mpicomm_shared, ensObs_mpiglobal%randPert_baseptr, &
                                   ensObs_mpiglobal%randPert_window, ierr)
      if (mmpi_myidHost /= 0) then
        call mpi_win_shared_query(ensObs_mpiglobal%randPert_window, 0, windowSize, dispUnit, &
                                  ensObs_mpiglobal%randPert_baseptr, ierr)
      end if
      call c_f_pointer(ensObs_mpiglobal%randPert_baseptr, ensObs_mpiglobal%randPert_r4, arrayShape)
    end if

    ! The baseptr now associated with Fortran pointer for access to global shared data
    if (mmpi_myid == 0) then
      allocate(tempBuffer(ensObsClean%numMembers,numObs_mpiglobal))
    else
      allocate(tempBuffer(ensObsClean%numMembers,1))
    end if

    ensObs_mpiglobal%typeVertCoord = ensObsClean%typeVertCoord

    if (mmpi_myid == 0) then
      displs(1) = 0
      do procIndex = 2, mmpi_nprocs
        displs(procIndex) = displs(procIndex-1) + allNumObs(procIndex-1)
      end do
    else
      displs(:) = 0
    end if

    call rpn_comm_gatherv(ensObsClean%lat, ensObsClean%numObs, 'mpi_real8', &
                          ensObs_mpiglobal%lat, allNumObs, displs, 'mpi_real8',  &
                          0, 'GRID', ierr)
    call rpn_comm_gatherv(ensObsClean%lon, ensObsClean%numObs, 'mpi_real8', &
                          ensObs_mpiglobal%lon, allNumObs, displs, 'mpi_real8',  &
                          0, 'GRID', ierr)
    call rpn_comm_gatherv(ensObsClean%vertLocation, ensObsClean%numObs, 'mpi_real8', &
                          ensObs_mpiglobal%vertLocation, allNumObs, displs, 'mpi_real8',  &
                          0, 'GRID', ierr)
    call rpn_comm_gatherv(ensObsClean%obsValue, ensObsClean%numObs, 'mpi_real8', &
                          ensObs_mpiglobal%obsValue, allNumObs, displs, 'mpi_real8',  &
                          0, 'GRID', ierr)
    call rpn_comm_gatherv(ensObsClean%obsErrInv, ensObsClean%numObs, 'mpi_real8', &
                          ensObs_mpiglobal%obsErrInv, allNumObs, displs, 'mpi_real8',  &
                          0, 'GRID', ierr)
    call rpn_comm_gatherv(ensObsClean%meanYb, ensObsClean%numObs, 'mpi_real8', &
                          ensObs_mpiglobal%meanYb, allNumObs, displs, 'mpi_real8',  &
                          0, 'GRID', ierr)
    call rpn_comm_gatherv(ensObsClean%deterYb, ensObsClean%numObs, 'mpi_real8', &
                          ensObs_mpiglobal%deterYb, allNumObs, displs, 'mpi_real8',  &
                          0, 'GRID', ierr)
    call rpn_comm_gatherv(ensObsClean%assFlag, ensObsClean%numObs, 'mpi_integer', &
                          ensObs_mpiglobal%assFlag, allNumObs, displs, 'mpi_integer',  &
                          0, 'GRID', ierr)
    call rpn_comm_gatherv(ensObsClean%codTyp, ensObsClean%numObs, 'mpi_integer', &
                          ensObs_mpiglobal%codTyp, allNumObs, displs, 'mpi_integer',  &
                          0, 'GRID', ierr)
    if (allocated(ensObsClean%obsErrInv_sim)) then
      call rpn_comm_gatherv(ensObsClean%obsErrInv_sim, ensObsClean%numObs, 'mpi_real8', &
                            ensObs_mpiglobal%obsErrInv_sim, allNumObs, displs, 'mpi_real8',  &
                            0, 'GRID', ierr)
    end if

    do memberIndex = 1, ensObsClean%numMembers
      call rpn_comm_gatherv(ensObsClean%Yb_r4(memberIndex,:), ensObsClean%numObs, 'mpi_real4', &
                            tempBuffer(memberIndex,:), allNumObs, displs, 'mpi_real4',  &
                            0, 'GRID', ierr)
    end do
    if (mmpi_myid == 0) then
      ensObs_mpiglobal%Yb_r4(:,:) = tempBuffer(:,:)
    end if
    call mpi_win_fence(0, ensObs_mpiglobal%Yb_window, ierr)

    if (associated(ensObsClean%Ya_r4)) then
      do memberIndex = 1, ensObsClean%numMembers
        call rpn_comm_gatherv(ensObsClean%Ya_r4(memberIndex,:), ensObsClean%numObs, 'mpi_real4', &
                              tempBuffer(memberIndex,:), allNumObs, displs, 'mpi_real4',  &
                              0, 'GRID', ierr)
      end do
      if (mmpi_myid == 0) then
        ensObs_mpiglobal%Ya_r4(:,:) = tempBuffer(:,:)
      end if
      call mpi_win_fence(0, ensObs_mpiglobal%Ya_window, ierr)
    end if

    if (associated(ensObsClean%randPert_r4)) then
      do memberIndex = 1, ensObsClean%numMembers
        call rpn_comm_gatherv(ensObsClean%randPert_r4(memberIndex,:), ensObsClean%numObs, 'mpi_real4', &
                              tempBuffer(memberIndex,:), allNumObs, displs, 'mpi_real4',  &
                              0, 'GRID', ierr)
      end do
      if (mmpi_myid == 0) then
        ensObs_mpiglobal%randPert_r4(:,:) = tempBuffer(:,:)
      end if
      call mpi_win_fence(0, ensObs_mpiglobal%randPert_window, ierr)
    end if

    call eob_deallocate(ensObsClean)

    deallocate(tempBuffer)

    call rpn_comm_bcast(ensObs_mpiglobal%lat, ensObs_mpiglobal%numObs, 'mpi_real8',  &
                        0, 'GRID', ierr)
    call rpn_comm_bcast(ensObs_mpiglobal%lon, ensObs_mpiglobal%numObs, 'mpi_real8',  &
                        0, 'GRID', ierr)
    call rpn_comm_bcast(ensObs_mpiglobal%vertLocation, ensObs_mpiglobal%numObs, 'mpi_real8',  &
                        0, 'GRID', ierr)
    call rpn_comm_bcast(ensObs_mpiglobal%obsValue, ensObs_mpiglobal%numObs, 'mpi_real8',  &
                        0, 'GRID', ierr)
    call rpn_comm_bcast(ensObs_mpiglobal%obsErrInv, ensObs_mpiglobal%numObs, 'mpi_real8',  &
                        0, 'GRID', ierr)
    if (allocated(ensObs_mpiglobal%obsErrInv_sim)) then
      call rpn_comm_bcast(ensObs_mpiglobal%obsErrInv_sim, ensObs_mpiglobal%numObs, 'mpi_real8',  &
                          0, 'GRID', ierr)
    end if
    call rpn_comm_bcast(ensObs_mpiglobal%meanYb, ensObs_mpiglobal%numObs, 'mpi_real8',  &
                        0, 'GRID', ierr)
    call rpn_comm_bcast(ensObs_mpiglobal%deterYb, ensObs_mpiglobal%numObs, 'mpi_real8',  &
                        0, 'GRID', ierr)
    call mmpi_bcast(ensObs_mpiglobal%assFlag, ensObs_mpiglobal%numObs)
    call mmpi_bcast(ensObs_mpiglobal%codTyp,  ensObs_mpiglobal%numObs)

    ! For shared memory we only need to send data to the other node masters
    nsize = ensObs_mpiglobal%numMembers*ensObs_mpiglobal%numObs
    if (mmpi_myid == 0) then
      write(*,*) 'eob_allGather: send data from task 0 to all node masters: ', mmpi_nodeMasters(:)
      allocate(requestId(size(mmpi_nodeMasters)-1))
      do procIndex = 2, size(mmpi_nodeMasters)
        call mpi_isend(ensObs_mpiglobal%Yb_r4(:,:),  &
                       nsize, mmpi_datyp_real4, mmpi_nodeMasters(procIndex), 1,  &
                       mmpi_comm_grid, requestId(procIndex-1), ierr)
      end do
      if (size(mmpi_nodeMasters) > 1) then
        call mpi_waitAll(size(mmpi_nodeMasters)-1, requestId(:), MPI_STATUSES_IGNORE, ierr)
      end if
      if (associated(ensObs_mpiglobal%Ya_r4)) then
        do procIndex = 2, size(mmpi_nodeMasters)
          call mpi_isend(ensObs_mpiglobal%Ya_r4(:,:),  &
                         nsize, mmpi_datyp_real4, mmpi_nodeMasters(procIndex), 2,  &
                         mmpi_comm_grid, requestId(procIndex-1), ierr)
        end do
        if (size(mmpi_nodeMasters) > 1) then
          call mpi_waitAll(size(mmpi_nodeMasters)-1, requestId(:), MPI_STATUSES_IGNORE, ierr)
        end if
      end if
      if (associated(ensObs_mpiglobal%randPert_r4)) then
        do procIndex = 2, size(mmpi_nodeMasters)
          call mpi_isend(ensObs_mpiglobal%randPert_r4(:,:),  &
                         nsize, mmpi_datyp_real4, mmpi_nodeMasters(procIndex), 3,  &
                         mmpi_comm_grid, requestId(procIndex-1), ierr)
        end do
        if (size(mmpi_nodeMasters) > 1) then
          call mpi_waitAll(size(mmpi_nodeMasters)-1, requestId(:), MPI_STATUSES_IGNORE, ierr)
        end if
      end if
    else if(mmpi_myidHost == 0) then
      write(*,*) 'eob_allGather: I am a nodeMaster, receive data for shared array'
      call mpi_recv(ensObs_mpiglobal%Yb_r4(:,:),  &
                    nsize, mmpi_datyp_real4, MPI_ANY_SOURCE, 1,  &
                    mmpi_comm_grid, MPI_STATUS_IGNORE, ierr)
      if (associated(ensObs_mpiglobal%Ya_r4)) then
        call mpi_recv(ensObs_mpiglobal%Ya_r4(:,:),  &
                      nsize, mmpi_datyp_real4, MPI_ANY_SOURCE, 2,  &
                      mmpi_comm_grid, MPI_STATUS_IGNORE, ierr)
      end if
      if (associated(ensObs_mpiglobal%randPert_r4)) then
        call mpi_recv(ensObs_mpiglobal%randPert_r4(:,:),  &
                      nsize, mmpi_datyp_real4, MPI_ANY_SOURCE, 3,  &
                      mmpi_comm_grid, MPI_STATUS_IGNORE, ierr)
      end if
    else
      write(*,*) 'eob_allGather: I am not a nodeMaster, do nothing'
    end if

    call mpi_win_fence(0, ensObs_mpiglobal%Yb_window, ierr)
    if (associated(ensObs_mpiglobal%Ya_r4)) then
      call mpi_win_fence(0, ensObs_mpiglobal%Ya_window, ierr)
    end if
    if (associated(ensObs_mpiglobal%randPert_r4)) then
      call mpi_win_fence(0, ensObs_mpiglobal%randPert_window, ierr)
    end if

    write(*,*) 'eob_allGather: Yb_r4(1,1), Yb_r4(end,0.5*end) = ', &
               ensObs_mpiglobal%Yb_r4(1,1), &
               ensObs_mpiglobal%Yb_r4(ensObs_mpiglobal%numMembers,nint(0.5*ensObs_mpiglobal%numObs))

    if (associated(ensObs_mpiglobal%randPert_r4)) then
      write(*,*) 'eob_allGather: randPert_r4(1,1), randPert_r4(end,0.5*end) = ', &
                 ensObs_mpiglobal%randPert_r4(1,1), &
                 ensObs_mpiglobal%randPert_r4(ensObs_mpiglobal%numMembers,nint(0.5*ensObs_mpiglobal%numObs))
    end if

    write(*,*) 'eob_allGather: total number of obs to be assimilated =', sum(ensObs_mpiglobal%assFlag(:))

    call utl_tmg_stop(24)
    call utl_tmg_stop(10)

    write(*,*) 'eob_allGather: finished'
    call msg_memUsage('eob_allGather')

  end subroutine eob_allGather

  !--------------------------------------------------------------------------
  ! eob_writeToFiles
  !--------------------------------------------------------------------------
  subroutine eob_writeToFiles(ensObs, outputFilenamePrefix, writeObsInfo, &
                              numGroupsToDivideMembers_opt, &
                              maxNumMembersPerGroup_opt)
    !
    ! :Purpose: Write the contents of an ensObs mpi local object to files
    !
    implicit none

    ! Arguments:
    type(struct_eob),  intent(in) :: ensObs                       ! input eob object
    character(len=*),  intent(in) :: outputFilenamePrefix         ! prefix for output files
    logical,           intent(in) :: writeObsInfo                 ! choose to write info common to all members
    integer, optional, intent(in) :: numGroupsToDivideMembers_opt ! number of groups to divide members in files
    integer, optional, intent(in) :: maxNumMembersPerGroup_opt    ! maximum number of members per group

    ! Locals:
    integer :: unitNum, ierr, obsIndex, memberIndex
    integer :: obsVcoCode(ensObs%numObs), obsAssFlag(ensObs%numObs)
    integer :: obsFlag(ensObs%numObs)
    integer, allocatable :: memberIndexArray(:)
    character(len=40) :: fileName
    character(len=4)  :: myidxStr, myidyStr
    character(len=30) :: fileNameExtention
    integer :: fnom, fclos
    logical :: fileExists

    if (.not. ensObs%allocated) then
      call utl_abort('eob_writeToFiles: this object is not allocated')
    end if
    if (ensObs%firstMember == 0) then
      call utl_abort('eob_writeToFiles: Cannot be called when firstMember is 0')
    end if

    call obs_extractObsIntBodyColumn(obsVcoCode, ensObs%obsSpaceData, OBS_VCO)
    call obs_extractObsIntBodyColumn(obsAssFlag, ensObs%obsSpaceData, OBS_ASS)
    call obs_extractObsIntBodyColumn(obsFlag, ensObs%obsSpaceData, OBS_FLG)

    write(myidxStr,'(I4.4)') (mmpi_myidx + 1)
    write(myidyStr,'(I4.4)') (mmpi_myidy + 1)
    fileNameExtention = trim(myidxStr) // '_' // trim(myidyStr)

    ! write observation info to a file
    if (writeObsInfo) then
      fileName = 'eob_obsInfo_' // trim(fileNameExtention)
      write(*,*) 'eob_writeToFiles: writing ',trim(filename)
      inquire(file=trim(fileName),exist=fileExists)
      if ( fileExists ) then
        call utl_abort('eob_writeToFiles: file should not exist')
      end if

      unitNum = 0
      ierr = fnom(unitNum, fileName, 'FTN+SEQ+UNF+R/W', 0)
      write(unitNum) ensObs%numMembers, ensObs%numObs
      write(unitNum) (ensObs%lat(obsIndex), obsIndex = 1, ensObs%numObs)
      write(unitNum) (ensObs%lon(obsIndex), obsIndex = 1, ensObs%numObs)
      write(unitNum) (obsVcoCode(obsIndex), obsIndex = 1, ensObs%numObs)
      write(unitNum) (ensObs%obsValue(obsIndex), obsIndex = 1, ensObs%numObs)
      write(unitNum) (obsAssFlag(obsIndex), obsIndex = 1, ensObs%numObs)
      write(unitNum) (obsFlag(obsIndex), obsIndex = 1, ensObs%numObs)
      ierr = fclos(unitNum)
    end if

    ! get memberIndex in the full ensemble set
    allocate(memberIndexArray(ensObs%numMembers))
    call getMemberIndexInFullEnsSet(ensObs, memberIndexArray, &
                                    numGroupsToDivideMembers_opt=numGroupsToDivideMembers_opt, &
                                    maxNumMembersPerGroup_opt=maxNumMembersPerGroup_opt)

    ! Open file and write ensObs%Yb for all the members to one file
    fileName = trim(outputFilenamePrefix) // '_' // trim(fileNameExtention)
    write(*,*) 'eob_writeToFiles: writing ',trim(filename)
    inquire(file=trim(fileName),exist=fileExists)
    if (fileExists) then
      call utl_abort('eob_writeToFiles: file should not exist')
    end if

    unitNum = 0
    ierr = fnom(unitNum, fileName, 'FTN+SEQ+UNF+R/W', 0)
    write(unitNum) ensObs%numMembers
    write(unitNum) (memberIndexArray(memberIndex), memberIndex = 1, ensObs%numMembers)
    do memberIndex = 1, ensObs%numMembers
      if (mmpi_myid == 0) then
        write(*,*) 'eob_writeToFiles: fileMemberIndex1=', ensObs%fileMemberIndex1, &
                   ', memberIndex=', memberIndex, &
                   ', memberIndex in full ensemble set=', memberIndexArray(memberIndex)
      end if
      write(unitNum) (ensObs%Yb_r4(memberIndex,obsIndex), obsIndex = 1, ensObs%numObs)
    end do
    ierr = fclos(unitNum)

    deallocate(memberIndexArray)

  end subroutine eob_writeToFiles

  !--------------------------------------------------------------------------
  ! eob_readFromFiles
  !--------------------------------------------------------------------------
  subroutine eob_readFromFiles(ensObs, numMembersToRead, inputFilenamePrefix, &
                               readObsInfo)
    !
    ! :Purpose: Read mpi local ensObs%Yb object from file. Several files in separate subdirectories
    !           can be read. Some examples of path+filename are:
    !           ensObs_0001/eob_HX_0001_0001
    !           ensObs_0002/eob_HX_0001_0001
    !
    implicit none

    ! Arguments:
    type(struct_eob), intent(inout) :: ensObs               ! eob object to be filled
    integer         , intent(in)    :: numMembersToRead     ! number of members to read
    character(len=*), intent(in)    :: inputFilenamePrefix  ! prefix of filename to read
    logical,          intent(in)    :: readObsInfo          ! read info common to all members

    ! Locals:
    real(8) :: latFromFile(ensObs%numObs), lonFromFile(ensObs%numObs)
    real(8) :: obsValueFromFile(ensObs%numObs)
    integer :: obsVcoCode(ensObs%numObs), obsVcoCodeFromFile(ensObs%numObs)
    integer :: obsFlag(ensObs%numObs), assFlagFrom1File(ensObs%numObs)
    integer :: assFlagFromAllFiles(ensObs%numObs)
    integer :: unitNum, ierr, memberIndex, obsIndex, numObsFromFile
    integer :: numMembersFromFile, fnom, fclos
    integer :: fileIndex, numMembersAlreadyRead
    integer, allocatable :: memberIndexFromFile(:)
    logical :: fileExists
    character(len=256) :: fileName
    character(len=100) :: fileBaseName
    character(len=4)   :: myidxStr, myidyStr
    character(len=3)   :: fileIndexStr
    character(len=30)  :: fileNameExtention

    if ( .not. ensObs%allocated ) then
      call utl_abort('eob_readFromFiles: this object is not allocated')
    end if
    if (ensObs%firstMember == 0) then
      call utl_abort('eob_readFromFiles: Cannot be called when firstMember is 0')
    end if

    call obs_extractObsIntBodyColumn(obsVcoCode, ensObs%obsSpaceData, OBS_VCO)

    write(myidxStr,'(I4.4)') (mmpi_myidx + 1)
    write(myidyStr,'(I4.4)') (mmpi_myidy + 1)
    fileNameExtention = trim(myidxStr) // '_' // trim(myidyStr)

    if (readObsInfo) then
      ! loop on all file from different directories, read obsInfo and check they match ensObs
      fileBaseName = 'eob_obsInfo_' // trim(fileNameExtention)

      fileIndex = 0
      numMembersAlreadyRead = 0
      do while (numMembersAlreadyRead < numMembersToRead)
        fileIndex = fileIndex + 1
        write(fileIndexStr,'(i3.3)') fileIndex
        fileName = './ensObs_' // fileIndexStr // '/' // fileBaseName

        write(*,*) 'eob_readFromFiles: reading ',trim(fileName)
        inquire(file=trim(fileName),exist=fileExists)
        if (.not. fileExists) then
          write(*,*) 'fileName=', fileName
          call utl_abort('eob_readFromFiles: file does not exist')
        end if

        unitNum = 0
        ierr = fnom(unitNum,trim(fileName),'FTN+SEQ+UNF',0)
        read(unitNum) numMembersFromFile, numObsFromFile
        if (ensObs%numObs /= numObsFromFile) then
          call utl_abort('eob_readFromFiles: ensObs%numObs does not match with that of file')
        end if

        read(unitNum) (latFromFile(obsIndex), obsIndex = 1, ensObs%numObs)
        read(unitNum) (lonFromFile(obsIndex), obsIndex = 1, ensObs%numObs)
        read(unitNum) (obsVcoCodeFromFile(obsIndex), obsIndex = 1, ensObs%numObs)
        read(unitNum) (obsValueFromFile(obsIndex), obsIndex = 1, ensObs%numObs)

        if (maxval(abs(latFromFile(:) - ensObs%lat(:))) > 1.0d-5 .or. &
            maxval(abs(lonFromFile(:) - ensObs%lon(:))) > 1.0d-5 .or. &
            maxval(abs(obsValueFromFile(:) - ensObs%obsValue(:))) > 1.0d-7 .or. &
            .not. all(obsVcoCodeFromFile(:) == obsVcoCode(:))) then

          call utl_abort('eob_readFromFiles: obsInfo file do not match ensObs')
        end if

        ! Read assimilation flag for of all files and apply a "logical or" to get the value
        !   to put in obsSpaceData. Read obs flag only on the first file.
        read(unitNum) (assFlagFrom1File(obsIndex), obsIndex = 1, ensObs%numObs)
        if (numMembersAlreadyRead == 0) then
          read(unitNum) (obsFlag(obsIndex), obsIndex = 1, ensObs%numObs)
        end if

        if (numMembersAlreadyRead == 0) assFlagFromAllFiles(:) = assFlagFrom1File(:)
        where (assFlagFrom1File(:) == obs_notAssimilated .and. numMembersAlreadyRead > 0)
          assFlagFromAllFiles(:) = obs_notAssimilated
        end where

        ierr = fclos(unitNum)

        numMembersAlreadyRead = numMembersAlreadyRead + numMembersFromFile
      end do

      ! update assimilation flag in obsSpaceData
      do obsIndex = 1, ensObs%numObs
        ! skip this obs it is already set to be assimilated
        if (assFlagFromAllFiles(obsIndex) == obs_assimilated) cycle

        call obs_bodySet_i(ensObs%obsSpaceData, OBS_ASS, obsIndex, obs_notAssimilated)
        call obs_bodySet_i(ensObs%obsSpaceData, OBS_FLG, obsIndex, obsFlag(obsIndex))
      end do
    end if

    ! loop on all files from different directories to read ensObs%Yb for all members
    fileBaseName = trim(inputFilenamePrefix) // '_' // trim(fileNameExtention)

    fileIndex = 0
    numMembersAlreadyRead = 0
    do while (numMembersAlreadyRead < numMembersToRead)
      fileIndex = fileIndex + 1
      write(fileIndexStr,'(i3.3)') fileIndex
      fileName = './ensObs_' // fileIndexStr // '/' // fileBaseName

      write(*,*) 'eob_readFromFiles: reading ',trim(fileName)
      inquire(file=trim(fileName),exist=fileExists)
      if (.not. fileExists) then
        write(*,*) 'fileName=', fileName
        call utl_abort('eob_readFromFiles: file does not exist')
      end if

      unitNum = 0
      ierr = fnom(unitNum,trim(fileName),'FTN+SEQ+UNF',0)
      read(unitNum) numMembersFromFile
      allocate(memberIndexFromFile(numMembersFromFile))
      read(unitNum) (memberIndexFromFile(memberIndex), memberIndex = 1, numMembersFromFile)
      do memberIndex = 1, numMembersFromFile
        read(unitNum) (ensObs%Yb_r4(memberIndexFromFile(memberIndex),obsIndex), obsIndex = 1, ensObs%numObs)
      end do
      ierr = fclos(unitNum)

      deallocate(memberIndexFromFile)

      numMembersAlreadyRead = numMembersAlreadyRead + numMembersFromFile
    end do

  end subroutine eob_readFromFiles

  !--------------------------------------------------------------------------
  ! eob_getLocalBodyIndices
  !--------------------------------------------------------------------------
  function eob_getLocalBodyIndices(ensObs,localBodyIndices,locFunSelected,lat,lon,vertLocation,  &
                                   hLocalize,vLocalize,numObsFound, &
                                   maxNumLocalObsPerType, localSelectionOutput, localObsSorting) &
                                   result(numObsSelected)
    !
    ! :Purpose: Return a list of values of bodyIndex for all observations within
    !           the local volume around the specified lat/lon used for assimilation
    !           (as defined by h/vLocalize). The kdtree2 module is used to efficiently
    !           perform this task. The kdtree itself is constructed on the first call.
    !
    implicit none

    ! Arguments:
    type(struct_eob), intent(in)  :: ensObs                ! input eob object
    integer         , intent(out) :: localBodyIndices(:)   ! body indexes of selected local obs
    real(8)         , intent(out) :: locFunSelected(:)     ! values of localization function of selected local obs
    real(8)         , intent(in)  :: lat                   ! reference location lat
    real(8)         , intent(in)  :: lon                   ! reference location lon
    real(8)         , intent(in)  :: vertLocation          ! reference location vertical position
    real(8)         , intent(in)  :: hLocalize             ! horizontal localization distance
    real(8)         , intent(in)  :: vLocalize             ! vertical localization distance
    integer         , intent(out) :: numObsFound           ! total number of local obs within the local volume
    logical         , intent(in)  :: localSelectionOutput  ! output information about the selection of observations
    integer         , intent(in)  :: maxNumLocalObsPerType ! maximum number of each obs type assimilated locally
    character(len=*), intent(in)  :: localObsSorting       ! sort method: HORIZONTAL distance, LOCFUN, or MINTRACE
    ! Result:
    integer                       :: numObsSelected        ! number of local obs up to the array size

    ! Locals:
    integer :: bodyIndex, numObsFoundSearch, maxNumLocalObs, localObsIndex, sortedIndex
    type(struct_eobOutput)            :: eobOut ! variables for the output information
    real(kdkind), allocatable         :: positionArray(:,:)
    type(kdtree2_result), allocatable :: searchResults(:)        ! kdtree results
    real(kdkind)                      :: maxRadius
    real(kdkind)                      :: refPosition(3)
    integer                           :: numObsSelectedByCodtyp(codtyp_maxNumber)    ! number of selected observations for each type
    integer                           :: numObsNotSelectedByCodtyp(codtyp_maxNumber) ! number of not selected observations for each type
    integer                           :: codTyp
    real(8)                           :: hDistance, vDistance
    real(8), allocatable              :: sortValue(:)         ! values used to sort the kdtree results, choice of values is given by localObsSorting
    integer, allocatable              :: sortIndex(:)         ! sorted indices of the kdtree results
    real(8), allocatable              :: locFun(:)            ! localization function of observations inside the localization volume

    ! create the kdtree on the first call
    if (.not. associated(tree)) then
      write(*,*) 'eob_getLocalBodyIndices: start creating kdtree'
      call msg_memUsage('eob_getLocalBodyIndices')
      allocate(positionArray(3,ensObs%numObs))
      do bodyIndex = 1, ensObs%numObs
        positionArray(1,bodyIndex) = ec_ra * sin(ensObs%lon(bodyIndex)) * cos(ensObs%lat(bodyIndex))
        positionArray(2,bodyIndex) = ec_ra * cos(ensObs%lon(bodyIndex)) * cos(ensObs%lat(bodyIndex))
        positionArray(3,bodyIndex) = ec_ra *                              sin(ensObs%lat(bodyIndex))
      end do
      tree => kdtree2_create(positionArray, rearrange=.true.)
      write(*,*) 'eob_getLocalBodyIndices: done creating kdtree'
      call msg_memUsage('eob_getLocalBodyIndices')
    end if

    ! do the search
    maxNumLocalObs = size(localBodyIndices)
    maxRadius = hLocalize**2
    refPosition(1) = ec_ra * sin(lon) * cos(lat)
    refPosition(2) = ec_ra * cos(lon) * cos(lat)
    refPosition(3) = ec_ra *            sin(lat)
    allocate(searchResults(maxNumLocalObsSearch))
    call kdtree2_r_nearest(tp=tree, qv=refPosition, r2=maxRadius, nfound=numObsFoundSearch,&
                           nalloc=maxNumLocalObsSearch, results=searchResults)
    if (numObsFoundSearch > maxNumLocalObsSearch) then
      call utl_abort('eob_getLocalBodyIndices: the parameter maxNumLocalObsSearch must be increased')
    end if

    ! searchResults is an array of maxNumLocalObsSearch observations.
    ! Its first numObsFoundSearch observations are all the ones inside the localization column
    ! and are sorted by horizontal distance.

    ! sortValue are the values used for sorting the numObsFoundSearch obs of the localization column.
    ! sortIndex are the indices of the searchResults array.
    !
    ! An observation with index sortIndex(i) has sort value sortValue(i), localization function locFun(i)
    ! and is found in searchResults(sortIndex(i)) and in ensObs(searchResults(sortIndex(i))%idx).

    ! create values to sort the observations
    allocate(sortValue(numObsFoundSearch))
    allocate(sortIndex(numObsFoundSearch))
    allocate(locFun(numObsFoundSearch))
    ! loop on searchResults
    do localObsIndex = 1, numObsFoundSearch
      bodyIndex           = searchResults(localObsIndex)%idx
      hDistance           = sqrt(searchResults(localObsIndex)%dis)
      vDistance           = abs( vertLocation - ensObs%vertLocation(bodyIndex) )
      locFun(localObsindex) = lfn_Response(hDistance,hLocalize)
      if ( vLocalize > 0.0d0 .and. vertLocation /= MPC_missingValue_R8 ) then
         locFun(localObsIndex) = locFun(localObsindex) * lfn_Response(vDistance,vLocalize)
      endif
      sortIndex(localObsIndex) = localObsIndex
      if ( trim(localObsSorting) == 'HORIZONTAL' ) then
        sortValue(localObsIndex) = hDistance
      else if ( trim(localObsSorting) == 'LOCFUN' ) then
        sortValue(localObsIndex) = locFun(localObsindex)
      else if ( trim(localObsSorting) == 'MINTRACE' ) then
        sortValue(localObsIndex) = locFun(localObsIndex) * &
                                   sum(ensObs%Yb_r4(:,bodyIndex)**2) * ensObs%obsErrInv(bodyIndex)
      else
        call utl_abort('eob_getLocalBodyIndices: localObsSorting has unknown value:'//trim(localObsSorting))
      endif
    end do
    ! sort observations by sortValue
    if (numObsFoundSearch > 0) then
      call utl_heapsort1d(sortValue,sortIndex)
      if ( trim(localObsSorting) == 'LOCFUN' .or. trim(localObsSorting) == 'MINTRACE' ) then
        ! In these cases, we select observations by decreasing order of their sort value
        sortValue = sortValue(size(sortValue):1:-1)
        sortIndex = sortIndex(size(sortIndex):1:-1)
      endif
    endif

    if(localSelectionOutput) then
      call eob_zeroOutput(eobOut)
    endif

    if ( vLocalize > 0.0d0 .and. vertLocation /= MPC_missingValue_R8 ) then
      ! copy search results to output vectors, only those within vertical localization distance
      numObsFound      = 0
      numObsSelected   = 0
      numObsSelectedByCodtyp(:)    = 0
      numObsNotSelectedByCodtyp(:) = 0
      ! loop on all the observations inside the horizontal radius column
      do sortedIndex = 1, numObsFoundSearch
        localObsIndex = sortIndex(sortedIndex)
        bodyIndex     = searchResults(localObsIndex)%idx
        ! if the observation is inside the localization volume
        if (locFun(localObsIndex) > 0.0d0 .and. ensObs%assFlag(bodyIndex) == 1) then
          numObsFound  = numObsFound + 1
          codTyp       = ensObs%codTyp(bodyIndex)
          ! select the observation if the maximum numbers have not been reached
          if (numObsSelectedByCodtyp(codTyp) < maxNumLocalObsPerType .and. &
              numObsSelected < maxNumLocalObs) then
            numObsSelected                    = numObsSelected + 1
            localBodyIndices(numObsSelected)  = bodyIndex
            locFunSelected(numObsSelected)    = locFun(localObsIndex)
            numObsSelectedByCodtyp(codTyp)    = numObsSelectedByCodtyp(codTyp) + 1
            if(localSelectionOutput) then
              call eob_updateOutput(eobOut, ensObs, searchResults, codTyp, &
                                    localObsIndex, bodyIndex, vertLocation, locFun)
            end if
          else
            numObsNotSelectedByCodtyp(codTyp) = numObsNotSelectedByCodtyp(codTyp) + 1
          end if
        end if
      end do
    else
      ! no vertical location, so just copy results
      numObsFound      = 0
      numObsSelected   = 0
      numObsSelectedByCodtyp(:)    = 0
      numObsNotSelectedByCodtyp(:) = 0
      do sortedIndex = 1, numObsFoundSearch
        localObsIndex = sortIndex(sortedIndex)
        if (ensObs%assFlag(searchResults(localObsIndex)%idx) == 1) then
          numObsFound = numObsFound + 1
          if (numObsSelected < maxNumLocalObs) then
            numObsSelected = numObsSelected + 1
            localBodyIndices(numObsSelected) = searchResults(localObsIndex)%idx
            locFunSelected(numObsSelected)   = locFun(localObsIndex)
          end if
        end if
      end do
    end if

    ! send information about this grid point to an output text file
    if(localSelectionOutput) then
      call eob_writeOutput(eobOut, lat, lon, vertlocation, &
                           numObsSelectedByCodtyp, numObsNotSelectedByCodtyp, &
                           numObsFound, numObsSelected)
    end if

    deallocate(locFun)
    deallocate(sortIndex)
    deallocate(sortValue)
    deallocate(searchResults)

  end function eob_getLocalBodyIndices

  !--------------------------------------------------------------------------
  ! eob_zeroOutput
  !--------------------------------------------------------------------------
  subroutine eob_zeroOutput(eobOut)
    !
    !:Purpose:  Initialize the members of a struct_eobOutput variable to zero
    !
    implicit none

    ! Arguments:
    type(struct_eobOutput), intent(inout) :: eobOut ! eobOutput object

    eobOut%distMax(:)    = 0.0d0
    eobOut%distMean(:)   = 0.0d0
    eobOut%vertTop(:)    = 0.0d0
    eobOut%vertBottom(:) = 0.0d0
    eobOut%locFun(:)     = 0.0d0
    eobOut%trace(:)      = 0.0d0
    eobOut%obsErr(:)     = 0.0d0
    eobOut%ensSpread(:)  = 0.0d0

  end subroutine eob_zeroOutput

  !--------------------------------------------------------------------------
  ! eob_updateOutput
  !--------------------------------------------------------------------------
  subroutine eob_updateOutput(eobOut, ensObs, searchResults, codTyp, &
                              localObsIndex, bodyIndex, vertLocation, locFun)
    !
    !:Purpose:  Update the values of a struct_eobOutput variable
    !
    implicit none

    ! Arguments:
    type(struct_eobOutput), intent(inout) :: eobOut           ! eobOutput object
    type(struct_eob),       intent(in)    :: ensObs           ! eob object
    type(kdtree2_result),   intent(in)    :: searchResults(:) ! kdtree results
    integer,                intent(in)    :: codTyp           ! eobOut index
    integer,                intent(in)    :: localObsIndex    ! searchResults index
    integer,                intent(in)    :: bodyIndex        ! ensObs index
    real(8),                intent(in)    :: vertLocation     ! grid point vert location
    real(8),                intent(in)    :: locFun(:)        ! Localization function

    eobOut%distMax(codTyp)    = max(eobOut%distMax(codTyp), &
                                     sqrt(searchResults(localObsIndex)%dis))
    eobOut%vertBottom(codTyp) = max(eobOut%vertBottom(codTyp), &
                                     vertLocation-ensObs%vertLocation(bodyIndex))
    eobOut%vertTop(codTyp)    = min(eobOut%vertTop(codTyp), &
                                     vertLocation-ensObs%vertLocation(bodyIndex))
    eobOut%locFun(codTyp)     = eobOut%locFun(codTyp) + &
                                locFun(localObsIndex)
    eobOut%distMean(codTyp)   = eobOut%distMean(codTyp) + &
                                sqrt(searchResults(localObsIndex)%dis)
    eobOut%trace(codTyp)      = eobOut%trace(codTyp) + &
                                locFun(localObsIndex) * &
                                sum(ensObs%Yb_r4(:,bodyIndex)**2) * &
                                ensObs%obsErrInv(bodyIndex)
    eobOut%obsErr(codTyp)     = eobOut%obsErr(codTyp) + &
                                1.0d0 / max(ensObs%obsErrInv(bodyIndex),1.0d-30)
    eobOut%ensSpread(codTyp)  = eobOut%ensSpread(codTyp) + &
                                sum(ensObs%Yb_r4(:,bodyIndex)**2)

  end subroutine eob_updateOutput

  !--------------------------------------------------------------------------
  ! eob_writeOutput
  !--------------------------------------------------------------------------
  subroutine eob_writeOutput(eobOut, lat, lon, vertlocation, &
                             numObsSelectedByCodtyp, numObsNotSelectedByCodtyp, &
                             numObsFound, numObsSelected)
    !
    !:Purpose:  Create a text file eob_glbi_XXXXX with XXXX the mpi task
    !           and append the variable struct_eobOutput to the text file.
    !
    implicit none

    ! Arguments:
    type(struct_eobOutput), intent(inout) :: eobOut                       ! eobOutput object to write
    real(8),                intent(in)    :: lat, lon, vertlocation       ! LETKF gridpoint location
    integer,                intent(in)    :: numObsSelectedByCodtyp(:)    ! number of selected observations for each type
    integer,                intent(in)    :: numObsNotSelectedByCodtyp(:) ! number of not selected observations for each type
    integer,                intent(in)    :: numObsFound                  ! number of obs within the local volume
    integer,                intent(in)    :: numObsSelected               ! number of obs selected

    ! Locals:
    character(len=50)                  :: outfilename          ! filename for the output
    integer                            :: fclos, funit, ierr
    integer                            :: codTyp
    logical, save                      :: firstcall = .true.   ! firstcall to create the ouput and its header

    write(outfilename, '(I5.5)') mmpi_myid ! we assume there are less than 100 000 mpi tasks...
    outfilename = './eob_glbi_'//trim(adjustl(outfilename))
    if (firstcall) then
      firstcall = .false.
      call utl_open_asciifile(outfilename,funit)
      write(funit,'(12(A,1X))') 'gdp', &
            'lat', &
            'lon', &
            'lnp', &
            'dist', &
            'vbottom', &
            'vtop', &
            'numtot', &
            'numsel', &
            'meandist', &
            'meanlfn', &
            'sumtrace'
      write(funit,'(12(A,1X))') 'typ', &
            'codtyp', &
            'numsel', &
            'numrej', &
            'dist', &
            'vbottom', &
            'vtop', &
            'meandist', &
            'meanlfn', &
            'sumtrace', &
            'meanerr', &
            'meanspd'
    else
      call utl_open_asciifile(outfilename,funit)
    endif

    ! generic information about the gridpoint
    write(funit,'(A3,1X,2(F7.3,1X),4(F12.4,1X),I7,1X,I7,3(1X,ES10.2))') 'gdp', &
          lat*MPC_DEGREES_PER_RADIAN_R8, &
          lon*MPC_DEGREES_PER_RADIAN_R8, &
          vertlocation, &
          maxval(eobOut%distMax), &
          maxval(eobOut%vertBottom), &
          minval(eobOut%vertTop), &
          numObsFound, &
          numObsSelected, &
          sum(eobOut%distMean) / real(max(1,numObsSelected)), &
          sum(eobOut%locFun)   / real(max(1,numObsSelected)), &
          sum(eobOut%trace)
    ! followed by more detailed information for each observation type
    do codTyp = 1, size(numObsSelectedByCodtyp)
      if ((numObsSelectedByCodtyp(codTyp) + numObsNotSelectedByCodtyp(codTyp)) /= 0) then
        eobOut%distMean(codTyp)  = eobOut%distMean(codTyp) /  &
                                   real(max(1,numObsSelectedByCodtyp(codTyp)))
        eobOut%locFun(codTyp)    = eobOut%locFun(codTyp) /    &
                                   real(max(1,numObsSelectedByCodtyp(codTyp)))
        eobOut%obsErr(codTyp)    = eobOut%obsErr(codTyp) /    &
                                   real(max(1,numObsSelectedByCodtyp(codTyp)))
        eobOut%ensSpread(codTyp) = eobOut%ensSpread(codTyp) / &
                                   real(max(1,numObsSelectedByCodtyp(codTyp)))
        write(funit,'(*(A3,1X,I3,2(1X,I7),1X,F12.2,2(1X,F7.2),1X,F12.2,4(1X,ES9.2)))') 'typ', &
              codTyp, &
              numObsSelectedByCodtyp(codTyp), &
              numObsNotSelectedByCodtyp(codTyp), &
              eobOut%distMax(codTyp), &
              eobOut%vertBottom(codTyp), &
              eobOut%vertTop(codTyp), &
              eobOut%distMean(codTyp), &
              eobOut%locFun(codTyp), &
              eobOut%trace(codTyp), &
              eobOut%obsErr(codTyp), &
              eobOut%ensSpread(codTyp)
      end if
    end do

    ierr = fclos(funit)

  end subroutine eob_writeOutput

  !--------------------------------------------------------------------------
  ! eob_setLatLonObsCod
  !--------------------------------------------------------------------------
  subroutine eob_setLatLonObsCod(ensObs)
    implicit none

    ! Arguments:
    type(struct_eob), intent(inout) :: ensObs ! eob object to modify

    call obs_extractObsRealHeaderColumn(ensObs%lat, ensObs%obsSpaceData, OBS_LAT)
    call obs_extractObsRealHeaderColumn(ensObs%lon, ensObs%obsSpaceData, OBS_LON)
    call obs_extractObsRealBodyColumn(ensObs%obsValue, ensObs%obsSpaceData, OBS_VAR)
    call obs_extractObsIntHeaderColumn(ensObs%codTyp, ensObs%obsSpaceData, OBS_ITY)

  end subroutine eob_setLatLonObsCod

  !--------------------------------------------------------------------------
  ! eob_setobsErrInv
  !--------------------------------------------------------------------------
  subroutine eob_setobsErrInv(ensObs)
    implicit none

    ! Arguments:
    type(struct_eob), intent(inout) :: ensObs ! eob object to modify

    ! Locals:
    integer :: obsIndex

    call obs_extractObsRealBodyColumn(ensObs%obsErrInv, ensObs%obsSpaceData, OBS_OER)
    do obsIndex = 1, ensObs%numObs
      if(ensObs%obsErrInv(obsIndex) > 0.0d0) then
        ensObs%obsErrInv(obsIndex) = 1.0d0/(ensObs%obsErrInv(obsIndex)**2)
      else
        ensObs%obsErrInv(obsIndex) = 0.0d0
      end if
    end do

    ! read namelist if necessary and calculate obs error inverse for
    ! passive and simulated observations
    call eob_init()
    if (psvObsAssim) call eob_setPsvObsErrInv(ensObs)
    if (eob_simObsAssim) call eob_setSimObsErrInv(ensObs)

  end subroutine eob_setobsErrInv

  !--------------------------------------------------------------------------
  ! eob_setPsvObsErrInv
  !--------------------------------------------------------------------------
  subroutine eob_setPsvObsErrInv(ensObs)
    !
    !:Purpose:  Updates the inverse of the observation error variance
    !           for passive osbervations and stores this in ensObs%obsErrInv.
    !           This is done assuming that ensObs%obsErrInv was already set.
    !
    implicit none

    ! Arguments:
    type(struct_eob),  intent(inout) :: ensObs ! eob object to modify

    ! Locals:
    integer       :: obsIndex, headerIndex
    integer       :: codtyp, varnum, obsfamIndex
    character(2)  :: obsfamCurrent
    logical       :: psvFlag

    do obsIndex = 1, ensObs%numObs
      psvFlag = .false.
      headerIndex = obs_bodyElem_i(ensObs%obsSpaceData, OBS_HIND, obsIndex)
      obsfamCurrent = obs_getFamily(ensObs%obsSpaceData, headerIndex_opt=headerIndex)
      ! update obs error inverse to 0 if current observation is passive
      if (ANY(psvObsFamily == obsfamCurrent)) then
        obsfamIndex = utl_findloc(psvObsFamily(:), obsfamCurrent)
        if ((numPsvCodTyp(obsfamIndex) > 0) .and. (numPsvVarNum(obsfamIndex) > 0)) then
          ! at least 1 codtyp AND varnum specified for current obs family so
          ! see if current observation matches any of those codtypes AND
          ! any of those varnums
          codtyp = obs_headElem_i(ensObs%obsSpaceData, OBS_ITY, headerIndex)
          varnum = obs_bodyElem_i(ensObs%obsSpaceData, OBS_VNM, obsIndex)
          if (ANY(psvCodTyp(obsfamIndex,:) == codtyp) .and. ANY(psvVarNum(obsfamIndex,:) == varnum)) then
            ensObs%obsErrInv(obsIndex) = 0.0d0
            psvFlag = .true.
          end if
        else if (numPsvCodTyp(obsfamIndex) > 0) then
          ! at least 1 codtype is specified for current obs family so
          ! see if current observation matches any of those codtypes
          codtyp = obs_headElem_i(ensObs%obsSpaceData, OBS_ITY, headerIndex)
          if (ANY(psvCodTyp(obsfamIndex,:) == codtyp)) then
            ensObs%obsErrInv(obsIndex) = 0.0d0
            psvFlag = .true.
          end if
        else if (numPsvVarNum(obsfamIndex) > 0) then
          ! at least 1 varnum is specified for current obs family so
          ! see if current observation matches any of those varnums
          varnum = obs_bodyElem_i(ensObs%obsSpaceData, OBS_VNM, obsIndex)
          if (ANY(psvVarNum(obsfamIndex,:) == varnum)) then
            ensObs%obsErrInv(obsIndex) = 0.0d0
            psvFlag = .true.
          end if
        else
          ! passive observation family doesn't include any codtype or
          ! any varnum, so set error inverse to 0 irrespective of current
          ! observation codtype and varnum
          ensObs%obsErrInv(obsIndex) = 0.0d0
          psvFlag = .true.
        end if
        ! set OBS_FLG to indicate passive observation
        if (psvFlag) call obs_bodySet_i(ensObs%obsSpaceData,OBS_FLG,obsIndex, &
                                        ibset(obs_bodyElem_i(ensObs%obsSpaceData,OBS_FLG,obsIndex),25))
      end if
    end do

  end subroutine eob_setPsvObsErrInv

  !--------------------------------------------------------------------------
  ! eob_setSimObsErrInv
  !--------------------------------------------------------------------------
  subroutine eob_setSimObsErrInv(ensObs)
    !
    !:Purpose:  Computes the inverse of the observation error variance if
    !           simulating any observations. Stores this in
    !           ensObs%obsErrInv_sim.
    !
    implicit none

    ! Arguments:
    type(struct_eob),  intent(inout) :: ensObs ! eob object to modify

    ! Locals:
    integer       :: obsIndex, headerIndex
    integer       :: codtyp, varnum, obsfamIndex
    character(2)  :: obsfamCurrent
    logical       :: simFlag

    write(*,*) 'eob_setSimObsErrInv: starting'

    ! set to copy of regular obs error inverse
    ensObs%obsErrInv_sim(:) = ensObs%obsErrInv(:)

    ! loop through all observations, and update the obs error inverse
    ! to a value of 0 for simulated observations
    do obsIndex = 1, ensObs%numObs
      simFlag = .false.
      headerIndex = obs_bodyElem_i(ensObs%obsSpaceData, OBS_HIND, obsIndex)
      obsfamCurrent = obs_getFamily(ensObs%obsSpaceData, headerIndex_opt=headerIndex)
      ! update obs error inverse for mean update to 0 if current observation is simulated
      if (ANY(simObsFamily == obsfamCurrent)) then
        obsfamIndex = utl_findloc(simObsFamily(:), obsfamCurrent)
        if ((numSimCodTyp(obsfamIndex) > 0) .and. (numSimVarNum(obsfamIndex) > 0)) then
          ! at least 1 codtyp AND varnum specified for current obs family so
          ! see if current observation matches any of those codtypes AND
          ! any of those varnums
          codtyp = obs_headElem_i(ensObs%obsSpaceData, OBS_ITY, headerIndex)
          varnum = obs_bodyElem_i(ensObs%obsSpaceData, OBS_VNM, obsIndex)
          if (ANY(simCodTyp(obsfamIndex,:) == codtyp) .and. ANY(simVarNum(obsfamIndex,:) == varnum)) then
            ensObs%obsErrInv_sim(obsIndex) = 0.0d0
            simFlag = .true.
          end if
        else if (numSimCodTyp(obsfamIndex) > 0) then
          ! at least 1 codtype is specified for current obs family so
          ! see if current observation matches any of those codtypes
          codtyp = obs_headElem_i(ensObs%obsSpaceData, OBS_ITY, headerIndex)
          if (ANY(simCodTyp(obsfamIndex,:) == codtyp)) then
            ensObs%obsErrInv_sim(obsIndex) = 0.0d0
            simFlag = .true.
          end if
        else if (numSimVarNum(obsfamIndex) > 0) then
          ! at least 1 varnum is specified for current obs family so
          ! see if current observation matches any of those varnums
          varnum = obs_bodyElem_i(ensObs%obsSpaceData, OBS_VNM, obsIndex)
          if (ANY(simVarNum(obsfamIndex,:) == varnum)) then
            ensObs%obsErrInv(obsIndex) = 0.0d0
            simFlag = .true.
          end if
        else
          ! simulated observation family doesn't include any codtype or
          ! any varnum, so set error inverse to 0 irrespective of current
          ! observation codtype and varnum
          ensObs%obsErrInv_sim(obsIndex) = 0.0d0
          simFlag = .true.
        end if
        ! set OBS_FLG to indicate simulated observation
        if (simFlag) call obs_bodySet_i(ensObs%obsSpaceData,OBS_FLG,obsIndex, &
                                        ibset(obs_bodyElem_i(ensObs%obsSpaceData,OBS_FLG,obsIndex),24))
      end if
    end do

  end subroutine eob_setSimObsErrInv

  !--------------------------------------------------------------------------
  ! eob_setVertLocation
  !--------------------------------------------------------------------------
  subroutine eob_setVertLocation(ensObs, columnMeanTrl)
    !
    ! :Purpose: Set the vertical location value for each observation that
    !           will be used when doing vertical localization. For
    !           radiance observations, the level of the maximum value
    !           of the derivative of transmission is used. This value
    !           is also written in obsSpaceData in the OBS_ZHA column.
    !
    implicit none

    ! Arguments:
    type(struct_eob)       , intent(inout) :: ensObs        ! eob object to modify
    type(struct_columnData), intent(in)    :: columnMeanTrl ! column object to supply vertical level info

    ! Locals:
    integer          :: obsIndex, headerIndex, channelIndex, numTovsLevels, sensorIndex
    integer          :: levIndex, levIndexBelow, levIndexAbove, nLev_M
    integer          :: varNumber(ensObs%numObs), obsVcoCode(ensObs%numObs), codType(ensObs%numObs)
    real(8)          :: obsHeight, interpFactor, obsPPP(ensObs%numObs)
    real(8), pointer :: sfcPres_ptr(:,:), presM_ptr(:,:), heightM_ptr(:,:)
    type(rttov_profile), pointer :: profiles(:)
    type(struct_vco),    pointer :: vco_ptr
    logical          :: verbose = .false.

    call eob_setAssFlag(ensObs)

    call obs_extractObsRealBodyColumn(obsPPP, ensObs%obsSpaceData, OBS_PPP)
    call obs_extractObsIntBodyColumn(varNumber, ensObs%obsSpaceData, OBS_VNM)
    call obs_extractObsIntBodyColumn(obsVcoCode, ensObs%obsSpaceData, OBS_VCO)
    call obs_extractObsIntHeaderColumn(codType, ensObs%obsSpaceData, OBS_ITY)

    if (ensObs%typeVertCoord == 'logPressure') then

      presM_ptr   => col_getAllColumns(columnMeanTrl,'P_M')
      heightM_ptr => col_getAllColumns(columnMeanTrl,'Z_M')
      sfcPres_ptr => col_getAllColumns(columnMeanTrl,'P0')
      nLev_M = col_getNumLev(columnMeanTrl,'MM')

      call tvs_getProfile(profiles,'nl')

    end if

    OBS_LOOP: do obsIndex = 1, ensObs%numObs
      headerIndex = obs_bodyElem_i(ensObs%obsSpaceData,OBS_HIND,obsIndex)

      if( varNumber(obsIndex) == BUFR_NETS .or. varNumber(obsIndex) == BUFR_NEPS   .or.  &
          varNumber(obsIndex) == BUFR_NEUS .or. varNumber(obsIndex) == BUFR_NEVS   .or.  &
          varNumber(obsIndex) == BUFR_NESS .or. varNumber(obsIndex) == BUFR_NEPN   .or.  &
          varNumber(obsIndex) == BUFR_VIS  .or. varNumber(obsIndex) == BUFR_LOGVIS .or.  &
          varNumber(obsIndex) == BUFR_GUST .or.  &
          varNumber(obsIndex) == BUFR_radarPrecip .or. varNumber(obsIndex) == BUFR_logRadarPrecip ) then

        ! all surface observations
        if (ensObs%typeVertCoord == 'logPressure') then
          ensObs%vertLocation(obsIndex) = log(sfcPres_ptr(1,headerIndex))
        else if (ensObs%typeVertCoord == 'depth') then
          ensObs%vertLocation(obsIndex) = 0.0D0
        else
          call utl_abort('eob_setVertLocation: unknown typeVertCoord:' // trim(ensObs%typeVertCoord))
        end if

      else if (varNumber(obsIndex) == BUFR_NEZD) then

        ! ZTD observation, try 0.7*Psfc (i.e. ~700hPa when Psfc=1000hPa)
        if (ensObs%typeVertCoord == 'logPressure') then
          ensObs%vertLocation(obsIndex) = log(0.7D0 * sfcPres_ptr(1,headerIndex))
        else
          call utl_abort('eob_setVertLocation: ZTD obs only compatible with logPressure coordinate')
        end if

      else if (obsPPP(obsIndex) > 0.0d0 .and. obsVcoCode(obsIndex)==2) then

        ! all pressure level observations
        if (ensObs%typeVertCoord == 'logPressure') then
          ensObs%vertLocation(obsIndex) = log(obsPPP(obsIndex))
        else
          call utl_abort('eob_setVertLocation: pressure obs only compatible with logPressure coordinate')
        end if

      else if(obsVcoCode(obsIndex) == 1 .and. .not.bufr_isOceanObs(varNumber(obsIndex))) then

        if (ensObs%typeVertCoord /= 'logPressure') then
          ! skip this obs if it will not be assimilated
          if (ensObs%assFlag(obsIndex) == 0) cycle OBS_LOOP
          ! otherwise, abort
          write(*,*) 'eob_setVertLocation: varNum = ', varNumber(obsIndex)
          call utl_abort('eob_setVertLocation: height level obs only compatible with logPressure coordinate')
        end if

        ! all height level observations (not including surface obs)
        obsHeight = obsPPP(obsIndex)

        ! find level just below the observation
        levIndexBelow = 0
        LEV_LOOP: do levIndex = 1, nLev_M
          if (obsHeight > heightM_ptr(levIndex,headerIndex)) then
            levIndexBelow = levIndex
            exit LEV_LOOP
          end if
        end do LEV_LOOP

        ! set the log pressure for observation
        if (levIndexBelow == 1) then
          ! above top level, use top level pressure
          ensObs%vertLocation(obsIndex) = log(presM_ptr(1,headerIndex))
        else if (levIndexBelow == 0) then
          ! below bottom level, use surface pressure
          ensObs%vertLocation(obsIndex) = log(sfcPres_ptr(1,headerIndex))
        else
          ! interpolate
          levIndexAbove = levIndexBelow - 1
          interpFactor = ( obsHeight                              - heightM_ptr(levIndexBelow,headerIndex) ) /  &
                         ( heightM_ptr(levIndexAbove,headerIndex) - heightM_ptr(levIndexBelow,headerIndex) )
          ensObs%vertLocation(obsIndex) = interpFactor           * log(presM_ptr(levIndexAbove,headerIndex)) +  &
                                          (1.0D0 - interpFactor) * log(presM_ptr(levIndexBelow,headerIndex))
        end if

      else if(tvs_isIdBurpTovs(codType(obsIndex))) then

        if (ensObs%typeVertCoord /= 'logPressure') then
          call utl_abort('eob_setVertLocation: radiance obs only compatible with logPressure coordinate')
        end if

        sensorIndex = tvs_lsensor(headerIndex)
        numTovsLevels   = size(tvs_transmission(headerIndex)%tau_levels,1)
        channelIndex = nint(obsPPP(obsIndex))
        channelIndex = max(0,min(channelIndex,tvs_maxChannelNumber+1))
        channelIndex = channelIndex - tvs_channelOffset(sensorIndex)
        channelIndex = utl_findloc(tvs_ichan(:,sensorIndex), channelIndex)
        if (channelIndex > 0 .and. ensObs%assFlag(obsIndex)==1) then
          call max_transmission(tvs_transmission(headerIndex), numTovsLevels, &
                                channelIndex, profiles(headerIndex)%p, ensObs%vertLocation(obsIndex))
          if(mmpi_myid == 0 .and. verbose) then
            write(*,*) 'eob_setVertLocation for tovs: ', codType(obsIndex), &
                       obsPPP(obsIndex), 0.01*exp(ensObs%vertLocation(obsIndex))
          end if
        else
          ensObs%vertLocation(obsIndex) = log(500.0D2)
        end if

      else if (varNumber(obsIndex) == BUFR_SST) then

        if (ensObs%typeVertCoord /= 'depth') then
          call utl_abort('eob_setVertLocation: SST obs only compatible with ocean depth coordinate')
        end if

        ! SST observations
        vco_ptr => col_getVco(columnMeanTrl)
        ensObs%vertLocation(obsIndex) = minval(vco_ptr%depths(:))

      else if(ensObs%assFlag(obsIndex)==1) then

        write(*,*) 'eob_setLatLonPresObs: ERROR! cannot compute pressure for this observation: ',  &
                   obsPPP(obsIndex), varNumber(obsIndex), obsVcoCode(obsIndex)
        call utl_abort('eob_setVertLocation')

      end if

      ! write the value into obsSpaceData for later diagnostics
      if (ensObs%assFlag(obsIndex)==1) then
        call obs_bodySet_r(ensObs%obsSpaceData, OBS_ZHA, obsIndex,  &
                           ensObs%vertLocation(obsIndex))
      end if

    end do OBS_LOOP

    nullify(profiles)

  end subroutine eob_setVertLocation

  !--------------------------------------------------------------------------
  ! eob_setAssFlag
  !--------------------------------------------------------------------------
  subroutine eob_setAssFlag(ensObs)
    implicit none

    ! Arguments:
    type(struct_eob), intent(inout) :: ensObs  ! eob object to modify

    call obs_extractObsIntBodyColumn(ensObs%assFlag, ensObs%obsSpaceData, OBS_ASS)

  end subroutine eob_setAssFlag

  !--------------------------------------------------------------------------
  ! eob_setYb
  !--------------------------------------------------------------------------
  subroutine eob_setYb(ensObs, memberIndex)
    implicit none

    ! Arguments:
    type(struct_eob), intent(inout) :: ensObs      ! eob object to modify
    integer         , intent(in)    :: memberIndex ! index of member to set

    if (memberIndex < ensObs%firstMember .or. memberIndex > ensObs%numMembers) then
      write(*,*) 'memberIndex = ', memberIndex
      call utl_abort('eob_setYb: memberIndex value is invalid')
    end if

    ! get the Y-HX value from obsSpaceData
    call obs_extractObsRealBodyColumn_r4(ensObs%Yb_r4(memberIndex,:), ensObs%obsSpaceData, OBS_OMP)

    ! now compute HX = Y - (Y-HX)
    ensObs%Yb_r4(memberIndex,:) = ensObs%obsValue(:) - ensObs%Yb_r4(memberIndex,:)

  end subroutine eob_setYb

  !--------------------------------------------------------------------------
  ! eob_setYa (like eob_setYb but for the analysis)
  !--------------------------------------------------------------------------
  subroutine eob_setYa(ensObs, memberIndex, obsColumnName)
    implicit none

    ! Arguments:
    type(struct_eob), intent(inout)  :: ensObs        ! eob object to modify
    integer         , intent(in)     :: memberIndex   ! index of member to set
    integer         , intent(in)     :: obsColumnName ! name of obs column to get the value

    if ( .not. associated(ensObs%Ya_r4) ) then
      call utl_abort('eob_setYa: ensObs%Ya_r4 must be allocated and it is not')
    end if

    ! get the Y-HX value from obsSpaceData
    call obs_extractObsRealBodyColumn_r4(ensObs%Ya_r4(memberIndex,:), ensObs%obsSpaceData, obsColumnName)

    ! now compute HX = Y - (Y-HX)
    ensObs%Ya_r4(memberIndex,:) = ensObs%obsValue(:) - ensObs%Ya_r4(memberIndex,:)

  end subroutine eob_setYa

  !--------------------------------------------------------------------------
  ! eob_setSimObsVal
  !--------------------------------------------------------------------------
  subroutine eob_setSimObsVal(ensObs)
    !
    ! :Purpose: Set the observed value for simulated observations to
    !           the background ensemble mean in observation space.
    !
    implicit none

    ! Arguments:
    type(struct_eob) , intent(inout)  :: ensObs ! eob object to modify

    ! Locals:
    integer       :: obsIndex, headerIndex
    integer       :: codtyp, varnum, obsfamIndex
    character(2)  :: obsfamCurrent

    if (.not. eob_simObsAssim) return
    ! Loop through observations and set y to mean(H(x)) if y is in obs family of interest
    do obsIndex = 1, ensObs%numObs
      headerIndex = obs_bodyElem_i(ensObs%obsSpaceData, OBS_HIND, obsIndex)
      obsfamCurrent = obs_getFamily(ensObs%obsSpaceData, headerIndex_opt=headerIndex)
      if (ANY(simObsFamily == obsfamCurrent)) then
        obsfamIndex = utl_findloc(simObsFamily(:), obsfamCurrent)
        if ((numSimCodTyp(obsfamIndex) > 0) .and. (numSimVarNum(obsfamIndex) > 0)) then
          ! at least 1 codtyp AND varnum specified for current obs family so
          ! see if current observation matches any of those codtypes AND
          ! any of those varnums
          codtyp = obs_headElem_i(ensObs%obsSpaceData, OBS_ITY, headerIndex)
          varnum = obs_bodyElem_i(ensObs%obsSpaceData, OBS_VNM, obsIndex)
          if (ANY(simCodTyp(obsfamIndex,:) == codtyp) .and. ANY(simVarNum(obsfamIndex,:) == varnum)) then
            ensObs%obsvalue(obsIndex) = ensObs%meanYb(obsIndex)
          end if
        else if (numSimCodTyp(obsfamIndex) > 0) then
          ! at least 1 codtype is specified for current obs family so
          ! see if current observation matches any of those codtypes
          codtyp = obs_headElem_i(ensObs%obsSpaceData, OBS_ITY, headerIndex)
          if (ANY(simCodTyp(obsfamIndex,:) == codtyp)) then
            ensObs%obsvalue(obsIndex) = ensObs%meanYb(obsIndex)
          end if
        else if (numSimVarNum(obsfamIndex) > 0) then
          ! at least 1 varnum is specified for current obs family so
          ! see if current observation matches any of those varnums
          varnum = obs_bodyElem_i(ensObs%obsSpaceData, OBS_VNM, obsIndex)
          if (ANY(simVarNum(obsfamIndex,:) == varnum)) then
            ensObs%obsvalue(obsIndex) = ensObs%meanYb(obsIndex)
          end if
        else
          ! simulated observation family doesn't include any codtype
          ! so set irrespective of current observation's codtype
          ensObs%obsvalue(obsIndex) = ensObs%meanYb(obsIndex)
        end if
      end if
    end do

  end subroutine eob_setSimObsVal

  !--------------------------------------------------------------------------
  ! eob_setDeterYb
  !--------------------------------------------------------------------------
  subroutine eob_setDeterYb(ensObs)
    implicit none

    ! Arguments:
    type(struct_eob), intent(inout) :: ensObs ! eob object to modify

    ! get the Y-HX value from obsSpaceData
    call obs_extractObsRealBodyColumn(ensObs%DeterYb(:), ensObs%obsSpaceData, OBS_OMP)

    ! now compute HX = Y - (Y-HX)
    ensObs%DeterYb(:) = ensObs%obsValue(:) - ensObs%DeterYb(:)

  end subroutine eob_setDeterYb

  !--------------------------------------------------------------------------
  ! eob_calcAndRemoveMeanYb
  !--------------------------------------------------------------------------
  subroutine eob_calcAndRemoveMeanYb(ensObs)
    implicit none

    ! Arguments:
    type(struct_eob), intent(inout) :: ensObs ! eob object to modify

    ! Locals:
    integer :: obsIndex

    do obsIndex = 1, ensObs%numObs
      ensObs%meanYb(obsIndex) = sum(ensObs%Yb_r4(:,obsIndex)) / ensObs%numMembers
      ensObs%Yb_r4(:,obsIndex) = ensObs%Yb_r4(:,obsIndex) - ensObs%meanYb(obsIndex)
    end do

    ensObs%meanRemoved = .true.

  end subroutine eob_calcAndRemoveMeanYb

  !--------------------------------------------------------------------------
  ! eob_calcRandPert
  !--------------------------------------------------------------------------
  subroutine eob_calcRandPert(ensObs, randomSeed)
    implicit none

    ! Arguments:
    type(struct_eob), intent(inout) :: ensObs     ! eob object to modify
    integer         , intent(in)    :: randomSeed ! value of random seed to use

    ! Locals:
    integer :: obsIndex, memberIndex
    real(4) :: meanRandPert, sigObs

    if (ensObs%mpiglobal) then
      call utl_abort('eob_calcRandPert: Input eob object must not already be mpiglobal')
    end if

    call rng_setup(abs(randomSeed))
    if ( associated(ensObs%randPert_r4) ) then
      call utl_abort('eob_calcRandPert: ensObs%randPert_r4 must not be already allocated')
    end if

    allocate(ensObs%randPert_r4(ensObs%numMembers,ensObs%numObs))

    do obsIndex = 1, ensObs%numObs
      sigObs = obs_bodyElem_r(ensObs%obsSpaceData, OBS_OER, obsIndex)
      do memberIndex = 1, ensObs%numMembers
        ensObs%randPert_r4(memberIndex,obsIndex) = sigObs * rng_gaussian()
      end do

      meanRandPert = sum(ensObs%randPert_r4(:,obsIndex)) / real(ensObs%numMembers,4)
      ensObs%randPert_r4(:,obsIndex) = ensObs%randPert_r4(:,obsIndex) - meanRandPert
    end do

  end subroutine eob_calcRandPert

  !--------------------------------------------------------------------------
  ! eob_setMeanOMP
  !--------------------------------------------------------------------------
  subroutine eob_setMeanOMP(ensObs)
    implicit none

    ! Arguments:
    type(struct_eob), intent(inout) :: ensObs ! eob object to modify

    ! Locals:
    integer :: obsIndex

    do obsIndex = 1, ensObs%numObs
      if (obs_bodyElem_i(ensObs%obsSpaceData, OBS_ASS, obsIndex) == obs_notAssimilated) cycle
      call obs_bodySet_r(ensObs%obsSpaceData, OBS_OMP, obsIndex,  &
                         ensObs%obsValue(obsIndex)-ensObs%meanYb(obsIndex))
    end do

  end subroutine eob_setMeanOMP

  !--------------------------------------------------------------------------
  ! eob_setHPHT
  !--------------------------------------------------------------------------
  subroutine eob_setHPHT(ensObs)
    implicit none

    ! Arguments:
    type(struct_eob), intent(inout) :: ensObs ! eob object to modify

    ! Locals:
    integer :: obsIndex, memberIndex
    real(8) :: hpht

    do obsIndex = 1, ensObs%numObs
      if (obs_bodyElem_i(ensObs%obsSpaceData, OBS_ASS, obsIndex) == obs_notAssimilated) cycle
      hpht = 0.0d0
      do memberIndex = 1, ensObs%numMembers
        hpht = hpht + ensObs%Yb_r4(memberIndex,obsIndex)**2 / ensObs%numMembers
      end do
      if (hpht > 0.0D0) then
        hpht = sqrt(hpht)
      else
        hpht = 0.0D0
      end if
      call obs_bodySet_r(ensObs%obsSpaceData, OBS_HPHT, obsIndex, hpht)

    end do

  end subroutine eob_setHPHT

  !--------------------------------------------------------------------------
  ! eob_backgroundCheck
  !--------------------------------------------------------------------------
  subroutine eob_backgroundCheck(ensObs)
    !
    ! :Purpose: Apply additional background using the ensemble spread.
    !
    implicit none

    ! Arguments::
    type(struct_eob), intent(inout) :: ensObs ! eob object to modify

    ! Locals:
    integer :: bodyIndexBeg, bodyIndexEnd, headerIndex, ivar, bodyIndex
    integer :: numRejected, numRejectedMpiGlobal, ierr, windCount
    real    :: sigo, sigb, omp, sig, reject_limit
    logical :: reject_wind

    numRejected = 0
    reject_limit = 5.0
    do headerIndex = 1, obs_numheader(ensObs%obsSpaceData)
      bodyIndexBeg = obs_headElem_i(ensObs%obsSpaceData, OBS_RLN, headerIndex)
      bodyIndexEnd = obs_headElem_i(ensObs%obsSpaceData, OBS_NLV, headerIndex) + bodyIndexBeg - 1
      reject_wind = .false.

      do bodyIndex = bodyIndexBeg, bodyIndexEnd
        if (obs_bodyElem_i(ensObs%obsSpaceData, OBS_ASS, bodyIndex) == obs_notAssimilated) cycle
        sigo = obs_bodyElem_r(ensObs%obsSpaceData, OBS_OER, bodyIndex)
        sigb = obs_bodyElem_r(ensObs%obsSpaceData, OBS_HPHT, bodyIndex)
        ! cut off at reject_limit standard deviations
        sig = reject_limit*(sigo**2 + sigb**2)**0.5
        omp = abs(obs_bodyElem_r(ensObs%obsSpaceData, OBS_OMP, bodyIndex))
        if (omp > sig) then
          call obs_bodySet_i(ensObs%obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
          call obs_bodySet_i(ensObs%obsSpaceData, OBS_FLG, bodyIndex,  &
                             IBSET(obs_bodyElem_i(ensObs%obsSpaceData, OBS_FLG, bodyIndex),9))
          ivar = obs_bodyElem_i(ensObs%obsSpaceData, OBS_VNM, bodyIndex)
          reject_wind = bufr_isWindComponent(ivar)
          write(*,*) 'eob_backgroundCheck: headerIndex, omp, sig, ivar, reject_wind', &
                                           headerIndex, omp, sig, ivar, reject_wind
          numRejected = numRejected + 1
        end if
      end do

      ! take the same reject decision for both components of the
      ! wind vector (and any other winds for that station).
      ! N.B.: This seems to assume only one level per header, this is generally
      ! the case for wind observations currently, since radiosondes are 4D
      if (reject_wind) then
        ! first count how many wind observations we have for this station
        windCount = 0
        do bodyIndex = bodyIndexBeg, bodyIndexEnd
          if (bufr_isWindComponent(obs_bodyElem_i(ensObs%obsSpaceData, OBS_VNM, bodyIndex))) then
            windCount = windCount + 1
          end if
        end do
        if (windCount > 2) then
          write(*,*) 'eob_backgroundCheck: WARNING'
          write(*,*) 'Station ',headerIndex,' has ',windCount,' wind observations '
          write(*,*) 'Perhaps old radiosonde format - std dev not changed for other wind component'
        else
          do bodyIndex = bodyIndexBeg, bodyIndexEnd
            if (obs_bodyElem_i(ensObs%obsSpaceData, OBS_ASS, bodyIndex) == obs_notassimilated) cycle
            ivar = obs_bodyElem_i(ensObs%obsSpaceData, OBS_VNM, bodyIndex)
            if (bufr_isWindComponent(ivar) .and.  &
                obs_bodyElem_i(ensObs%obsSpaceData, OBS_ASS, bodyIndex) == obs_assimilated) then
              call obs_bodySet_i(ensObs%obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
              call obs_bodySet_i(ensObs%obsSpaceData, OBS_FLG, bodyIndex,  &
                                 IBSET(obs_bodyElem_i(ensObs%obsSpaceData, OBS_FLG, bodyIndex),9))
              write(*,*) 'eob_backgroundCheck: other wind component headerIndex, ivar', &
                                               headerIndex, ivar
              numRejected = numRejected + 1
            end if
          end do
        end if
      end if
    end do

    call rpn_comm_allreduce(numRejected,numRejectedMpiGlobal,1,'mpi_integer','mpi_sum','GRID',ierr)
    write(*,*)
    write(*,*) 'eob_backgroundCheck: number of observations rejected (local) =', numRejected
    write(*,*) 'eob_backgroundCheck: number of observations rejected (global)=', numRejectedMpiGlobal

  end subroutine eob_backgroundCheck

  !--------------------------------------------------------------------------
  ! eob_removeObsNearLand
  !--------------------------------------------------------------------------
  subroutine eob_removeObsNearLand(ensObs, oceanMask, minDistanceToLand)
    !
    ! :Purpose: Reject observations that are close to land as determined by
    !           the argument "minDistanceToLand". In the case of a depth-
    !           varying land mask, the first level (should be the surface) is
    !           used.
    !
    implicit none

    ! Arguments:
    type(struct_eob), intent(inout) :: ensObs            ! eob object to modify
    type(struct_ocm), intent(in)    :: oceanMask         ! ocm object to define coast line
    real(8),          intent(in)    :: minDistanceToLand ! minimum distance from land

    ! Locals:
    integer :: headerIndex, bodyIndex, bodyIndexBeg, bodyIndexEnd, levIndex
    integer :: numRejected, numRejectedMpiGlobal, ierr
    real(8) :: obsLon, obsLat

    write(*,*) 'eob_removeObsNearLand: starting'

    numRejected = 0

    HEADER_LOOP: do headerIndex = 1, obs_numheader(ensObs%obsSpaceData)
      bodyIndexBeg = obs_headElem_i(ensObs%obsSpaceData, OBS_RLN, headerIndex)
      bodyIndexEnd = obs_headElem_i(ensObs%obsSpaceData, OBS_NLV, headerIndex) + bodyIndexBeg - 1

      levIndex = 1
      obsLat = obs_headElem_r(ensObs%obsSpaceData, OBS_LAT, headerIndex)
      obsLon = obs_headElem_r(ensObs%obsSpaceData, OBS_LON, headerIndex)

      ! skip this obs if it is far from land
      if (ocm_farFromLand(oceanMask, levIndex, obsLon, obsLat, minDistanceToLand)) then
        cycle HEADER_LOOP
      end if

      ! otherwise it is rejected
      BODY_LOOP: do bodyIndex = bodyIndexBeg, bodyIndexEnd
        if (obs_bodyElem_i(ensObs%obsSpaceData, OBS_ASS, bodyIndex) == obs_notAssimilated) cycle BODY_LOOP

        call obs_bodySet_i(ensObs%obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
        call obs_bodySet_i(ensObs%obsSpaceData, OBS_FLG, bodyIndex,  &
                           IBSET(obs_bodyElem_i(ensObs%obsSpaceData, OBS_FLG, bodyIndex),9))
        numRejected = numRejected + 1
      end do BODY_LOOP
    end do HEADER_LOOP

    call rpn_comm_allreduce(numRejected,numRejectedMpiGlobal,1,'mpi_integer','mpi_sum','GRID',ierr)
    write(*,*)
    write(*,*) 'eob_removeObsNearLand: number of observations rejected (local) =', numRejected
    write(*,*) 'eob_removeObsNearLand: number of observations rejected (global)=', numRejectedMpiGlobal

    write(*,*) 'eob_removeObsNearLand: finished'

  end subroutine eob_removeObsNearLand

  !--------------------------------------------------------------------------
  ! eob_setSigiSigo
  !--------------------------------------------------------------------------
  subroutine eob_setSigiSigo(ensObs)
    !
    ! :Purpose: Apply huber norm quality control procedure. This modifies
    !           the OBS_OER value, but before that its value is copied into
    !           OBS_SIGO and also OBS_SIGI computed
    !
    implicit none

    ! Arguments:
    type(struct_eob), intent(inout) :: ensObs ! eob object to modify

    ! Locals:
    integer           :: bodyIndex
    real(pre_obsReal) :: sigo, sigb, sigi

    ! Set 'sigi' and 'sigo' before oer is modified by Huber norm
    do bodyIndex = 1, obs_numbody(ensObs%obsSpaceData)
      sigb = obs_bodyElem_r(ensObs%obsSpaceData, OBS_HPHT, bodyIndex)
      sigo = obs_bodyElem_r(ensObs%obsSpaceData, OBS_OER, bodyIndex)
      if (obs_bodyElem_i(ensObs%obsSpaceData, OBS_ASS, bodyIndex) == obs_assimilated) then
        sigi = (sigo**2 + sigb**2)**0.5
        call obs_bodySet_r(ensObs%obsSpaceData, OBS_SIGI, bodyIndex, sigi)
        call obs_bodySet_r(ensObs%obsSpaceData, OBS_SIGO, bodyIndex, sigo)
      end if
    end do

  end subroutine eob_setSigiSigo

  !--------------------------------------------------------------------------
  ! eob_huberNorm
  !--------------------------------------------------------------------------
  subroutine eob_huberNorm(ensObs)
    !
    ! :Purpose: Apply huber norm quality control procedure. This modifies
    !           the OBS_OER value.
    !
    implicit none

    ! Arguments:
    type(struct_eob), intent(inout) :: ensObs ! eob object to modify

    ! Locals:
    integer           :: huberCount, huberCountMpiGlobal, ivar, windCount, ierr
    integer           :: bodyIndex, bodyIndexBeg, bodyIndexEnd, headerIndex
    real(pre_obsReal) :: c_limit, sig, sigo, sigb, omp, sigo_hub, sigo_hub_wind
    logical           :: reject_wind

    c_limit = 2.0
    huberCount = 0
    do headerIndex = 1, obs_numheader(ensObs%obsSpaceData)
      bodyIndexBeg = obs_headElem_i(ensObs%obsSpaceData, OBS_RLN, headerIndex)
      bodyIndexEnd = obs_headElem_i(ensObs%obsSpaceData, OBS_NLV, headerIndex) + bodyIndexBeg - 1
      reject_wind = .false.
      sigo_hub_wind = 0.0
      do bodyIndex = bodyIndexBeg, bodyIndexEnd
        sigo = obs_bodyElem_r(ensObs%obsSpaceData, OBS_OER, bodyIndex)
        sigb = obs_bodyElem_r(ensObs%obsSpaceData, OBS_HPHT, bodyIndex)
        ! cut off at reject_limit standard deviations
        sig = c_limit*(sigo**2 + sigb**2)**0.5
        omp = abs(obs_bodyElem_r(ensObs%obsSpaceData, OBS_OMP, bodyIndex))
        if (omp > sig) then
          ! redefining the observational error such that the innovation
          ! is at exactly c_limit standard deviations.
          huberCount = huberCount + 1
          sigo_hub = ((omp/c_limit)**2-sigb**2)**0.5
          ivar = obs_bodyElem_i(ensObs%obsSpaceData, OBS_VNM, bodyIndex)
          if (bufr_isWindComponent(ivar)) then
            ! the other wind components will be changed at the same time (next if)
            reject_wind = .true.
            call obs_bodySet_r(ensObs%obsSpaceData, OBS_OER, bodyIndex, sigo_hub)
            ! this is for the special case both components of the wind innovation
            ! suggest using the Huber norm.
            if (sigo_hub > sigo_hub_wind) then
              sigo_hub_wind = sigo_hub
            end if
          else
            call obs_bodySet_r(ensObs%obsSpaceData, OBS_OER, bodyIndex, sigo_hub)
          end if
        end if
      end do

      ! Give the same inflated observation error to both wind components.
      ! this part assumes that modern sondes are used (at most two
      ! wind observations per station.
      if (reject_wind) then
        ! first count how many wind observations we have for this station
        windCount = 0
        do bodyIndex = bodyIndexBeg, bodyIndexEnd
          if (bufr_isWindComponent(obs_bodyElem_i(ensObs%obsSpaceData, OBS_VNM, bodyIndex))) then
            windCount = windCount + 1
          end if
        end do
        if (windCount > 2) then
          write(*,*) 'Warning Hubernorm'
          write(*,*) 'Station ',headerIndex,' has ',windCount,' wind observations '
          write(*,*) 'Perhaps old radiosonde format - std dev not changed for other wind component'
        else
          huberCount = huberCount + 1
          do bodyIndex = bodyIndexBeg, bodyIndexEnd
            if (bufr_isWindComponent(obs_bodyElem_i(ensObs%obsSpaceData, OBS_VNM, bodyIndex))) then
              call obs_bodySet_r(ensObs%obsSpaceData, OBS_OER, bodyIndex, sigo_hub_wind)
            end if
          end do
        end if
      end if
    end do

    call rpn_comm_allreduce(huberCount, huberCountMpiGlobal, 1, 'mpi_integer', 'mpi_sum', 'GRID', ierr)
    write(*,*)
    write(*,*) 'eob_huberNorm: number of obs with increased error stddev (local) = ', huberCount
    write(*,*) 'eob_huberNorm: number of obs with increased error stddev (global)= ', huberCountMpiGlobal

  end subroutine eob_huberNorm

  !--------------------------------------------------------------------------
  ! eob_rejectRadNearSfc
  !--------------------------------------------------------------------------
  subroutine eob_rejectRadNearSfc(ensObs)
    !
    ! :Purpose: Reject all radiance observations with peak sensitivity
    !           too close to the surface.
    !
    implicit none

    ! Arguments:
    type(struct_eob), intent(inout) :: ensObs ! eob object to modify

    ! Locals:
    integer :: acceptCount, rejectCount, acceptCountMpiGlobal, rejectCountMpiGlobal
    integer :: varNumber, bodyIndex, ierr
    ! reject lower than 975 hPa
    real, parameter    :: logPresRadianceLimit = 11.4876E0

    acceptCount = 0
    rejectCount = 0
    do bodyIndex = 1, obs_numbody(ensObs%obsSpaceData)
      varNumber = obs_bodyElem_i(ensObs%obsSpaceData, OBS_VNM, bodyIndex)
      if (varNumber == bufr_nbt1 .or.  &
          varNumber == bufr_nbt2 .or.  &
          varNumber == bufr_nbt3) then
        if (ensObs%vertLocation(bodyIndex) < logPresRadianceLimit) then
          acceptCount = acceptCount + 1
        else
          rejectCount = rejectCount + 1
          call obs_bodySet_i(ensObs%obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
          call obs_bodySet_i(ensObs%obsSpaceData, OBS_FLG, bodyIndex,  &
                             IBSET(obs_bodyElem_i(ensObs%obsSpaceData, OBS_FLG, bodyIndex),9))
        end if
      end if
    end do

    call rpn_comm_allreduce(acceptCount, acceptCountMpiGlobal, 1, 'mpi_integer', 'mpi_sum', 'GRID', ierr)
    call rpn_comm_allreduce(rejectCount, rejectCountMpiGlobal, 1, 'mpi_integer', 'mpi_sum', 'GRID', ierr)
    write(*,*)
    write(*,*) 'eob_rejectRadNearSfc: Number of accepted, rejected observations (local) : ',  &
               acceptCount, rejectCount
    write(*,*) 'eob_rejectRadNearSfc: Number of accepted, rejected observations (global): ',  &
               acceptCountMpiGlobal, rejectCountMpiGlobal

  end subroutine eob_rejectRadNearSfc

  !--------------------------------------------------------------------------
  ! getMemberIndexInFullEnsSet (private routine)
  !--------------------------------------------------------------------------
  subroutine getMemberIndexInFullEnsSet(ensObs, memberIndexArray, &
                                        numGroupsToDivideMembers_opt, &
                                        maxNumMembersPerGroup_opt)
    !
    ! :Purpose: get memberIndex array corresponding to the full ensemble set. This
    !           is useful when ensObs is a subset of full ensemble members.
    !           If first member in ensObs is member 6, to get the full ensemble set equivalent of ensObs members:
    !           a) When members are not grouped (numGroupsToDivideMembers=1), all members are offset by
    !              memberIndexOffset (e.g. 6, 7, ..., 6+ensObs%numMermbers)
    !           b) When members are grouped (numGroupsToDivideMembers/=1), members within each group
    !              are offset by memberIndexOffset but there is increment of maxNumMembersPerGroup_opt
    !              to jump to the next group (e.g. if maxNumMembersPerGroup_opt=10, for first group 6, 7, 8, 9,
    !              for second group 6+10, 7+10, 8+10, 9+10, and so on)
    !
    implicit none

    ! Arguments:
    type(struct_eob),  intent(in)    :: ensObs                       ! input eob object
    integer,           intent(inout) :: memberIndexArray(:)          ! resulting list of member indexes
    integer, optional, intent(in)    :: numGroupsToDivideMembers_opt ! number of groups to divide the members
    integer, optional, intent(in)    :: maxNumMembersPerGroup_opt    ! maximum number of members per group

    ! Locals:
    integer :: memberIndex, groupIndex, memberIndexOffset, memberIndexInGroup
    integer :: numGroupsToDivideMembers, numMembersPerGroup

    if (present(numGroupsToDivideMembers_opt)) then
      numGroupsToDivideMembers = numGroupsToDivideMembers_opt
    else
      numGroupsToDivideMembers = 1
    end if

    memberIndexOffset = ensObs%fileMemberIndex1

    if (numGroupsToDivideMembers == 1) then
      do memberIndex = 1, ensObs%numMembers
        memberIndexArray(memberIndex) = memberIndex + memberIndexOffset - 1
      end do
    else
      if (.not. present(maxNumMembersPerGroup_opt)) then
        call utl_abort('getMemberIndexInFullEnsSet: maxNumMembersPerGroup_opt input argument missing')
      end if

      ! divide members into groups
      numMembersPerGroup = ensObs%numMembers / numGroupsToDivideMembers
      if (numMembersPerGroup > maxNumMembersPerGroup_opt) then
        call utl_abort('getMemberIndexInFullEnsSet: numMembersPerGroup > maxNumMembersPerGroup_opt')
      end if

      memberIndex = 0
      do groupIndex = 1, numGroupsToDivideMembers
        do memberIndexInGroup = 1, numMembersPerGroup
          memberIndex = memberIndex + 1
          memberIndexArray(memberIndex) = (groupIndex - 1) * maxNumMembersPerGroup_opt + &
                                          memberIndexInGroup + memberIndexOffset - 1
        end do
      end do

    end if

  end subroutine getMemberIndexInFullEnsSet

  !--------------------------------------------------------------------------
  ! max_transmission (private routine)
  !--------------------------------------------------------------------------
  subroutine max_transmission(transmission, numLevels, transIndex, rttovPres, maxLnP)
    !
    ! :Purpose: Determine the height in log pressure where we find the maximum
    !           value of the first derivative of transmission with respect to
    !           log pressure
    !
    implicit none

    ! Arguments:
    type(rttov_transmission), intent(in)  :: transmission ! transmission (rttov type)
    integer(kind=jpim)      , intent(in)  :: numLevels    ! number of RTTOV levels
    integer                 , intent(in)  :: transIndex   ! index of transmission%tau_levels
    real(kind=jprb), pointer, intent(in)  :: rttovPres(:) ! pressure of RTTOV levels
    real(8)                 , intent(out) :: maxLnP       ! computed log pressure of maximum

    ! Locals:
    integer :: levIndex
    real(8) :: lnPres(numLevels), avgPres(numLevels-1)
    real(8) :: diffTau, derivTau(numLevels), maxDeriv
    integer :: nAvgLev, maxIndex

    nAvgLev = numLevels - 1
    lnPres(:) = log(rttovPres(:)*MPC_PA_PER_MBAR_R8)
    ! calculate the first derivative of transmission with respect to log pressure
    ! and find the level index for its maximum
    maxDeriv = -0.1d0
    derivTau(1) = 0.0d0
    maxIndex = numLevels
    do levIndex = 2, numLevels
      avgPres(levIndex-1) = 0.5d0*(lnPres(levIndex)+lnPres(levIndex-1))
      diffTau = transmission%tau_levels(levIndex-1,transIndex) - transmission%tau_levels(levIndex,transIndex)
      derivTau(levIndex) = diffTau / (lnPres(levIndex)-lnPres(levIndex-1))
      if (derivTau(levIndex)>maxDeriv) then
        maxDeriv = derivTau(levIndex)
        maxIndex = levIndex
      end if
    end do

    ! get the height in log pressure for the level index (maxIndex) found above
    if (maxIndex==1) maxIndex = maxIndex + 1
    if ((maxIndex==2).or.(maxIndex==numLevels)) then
      maxLnP = avgPres(maxIndex-1)
    else
      call get_peak(maxIndex,nAvgLev,avgPres,derivTau,maxLnP)
    end if

  end subroutine max_transmission

  !--------------------------------------------------------------------------
  ! get_peak (private routine)
  !--------------------------------------------------------------------------
  subroutine get_peak(maxIndex,nlev,lnp,deriv,maxLnP)
    !
    ! :Purpose: Do quadratic interpolation to find pressure of peak transmission.
    !
    implicit none

    ! Arguments:
    integer, intent(in)    :: maxIndex      ! level index of the max value
    integer, intent(in)    :: nlev          ! number of levels for dimensioning arguments
    real(8), intent(in)    :: lnp(nlev)     ! natural log of pressure
    real(8), intent(in)    :: deriv(nlev+1) ! derivative with respect to log pressure
    real(8), intent(inout) :: maxLnP        ! computed log pressure of maximum

    ! Locals:
    external :: dgesv
    integer, parameter :: N=3
    integer :: info
    integer, parameter :: lda=N, ldb=N, nrhs=1
    integer :: ipiv(N)
    real(8) :: A(lda,N),B(ldb,nrhs)
    integer :: index1, index2

    index2 = 0
    do index1=maxIndex-1,maxIndex+1
      index2 = index2 + 1
      A(index2,1) = lnp(index1-1)*lnp(index1-1)
      A(index2,2) = lnp(index1-1)
      A(index2,3) = 1.0d0
      B(index2,1) = deriv(index1)
    end do

    call dgesv(N,nrhs,A,lda,ipiv,B,ldb,info)

    if (info==0) then
      maxLnP = -0.5*(B(2,1)/B(1,1))
    else
      maxLnP = lnp(maxIndex-1)
    end if

  end subroutine get_peak

end module ensembleObservations_mod

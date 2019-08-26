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

module gridVariableTransforms_mod
  ! MODULE gridVariableTransforms (prefix='gvt' category='3. High-level transformations')
  !
  ! :Purpose: To store various functions for variable transforms using inputs
  !           from gridStateVector(s). Outputs are also placed in a
  !           gridStateVector.
  !
  use mpivar_mod
  use mathPhysConstants_mod
  use earthConstants_mod
  use timeCoord_mod
  use gridStateVector_mod
  use ensembleStateVector_mod
  use lamSpectralTransform_mod
  use globalSpectralTransform_mod
  use analysisGrid_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use utilities_mod
  use varNameList_mod
  use tt2phi_mod
  implicit none
  save
  private

  ! public procedures
  public :: gvt_setup, gvt_transform, gvt_getStateVectorTrial

  logical                   :: huTrialsInitialized  = .false.
  logical                   :: heightTrialsInitialized  = .false.
  type(struct_hco), pointer :: hco_anl => null()
  type(struct_vco), pointer :: vco_anl => null()

  type(struct_gsv), target :: stateVectorTrialHU
  type(struct_gsv), target :: stateVectorTrialHeight

  ! module interfaces
  interface gvt_transform
    module procedure gvt_transform_gsv
    module procedure gvt_transform_ens
  end interface gvt_transform

CONTAINS

  subroutine gvt_setup(hco_in,vco_in)
    ! 
    ! :Purpose: To set up a variable transformation object
    !
    ! :Arguments:
    !           :hco_in: horizontal coordinate object input 
    !           :vco_in: vertical   coordinate object input
    !
    implicit none

    type(struct_hco), pointer :: hco_in
    type(struct_vco), pointer :: vco_in
    
    if (huTrialsInitialized) return
    if (heightTrialsInitialized) return

    write(*,*) 'gvt_setup: starting'

    hco_anl => hco_in
    vco_anl => vco_in

    write(*,*) 'gvt_setup: done'

  end subroutine gvt_setup


  subroutine gvt_setupTrials(varName)

    implicit none

    type(struct_gsv) :: statevector_noZnoP
    character(len=*), intent(in) :: varName

    select case ( trim(varName) )
    case ('HU')
      ! initialize stateVectorTrialHU on analysis grid
      call gsv_allocate(stateVectorTrialHU, tim_nstepobsinc, hco_anl, vco_anl,   &
                        dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                        allocHeightSfc_opt=.true., hInterpolateDegree_opt='LINEAR', &
                        varNames_opt=(/'HU','P0'/) )

      ! read trial files using default horizontal interpolation degree
      call gsv_readTrials( stateVectorTrialHU )  ! IN/OUT

      huTrialsInitialized = .true.
    case ('height')
      ! initialize stateVectorTrialHeight on analysis grid
      call gsv_allocate(stateVectorTrialHeight, tim_nstepobsinc, hco_anl, vco_anl,   &
                        dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                        allocHeightSfc_opt=.true., hInterpolateDegree_opt='LINEAR', &
                        varNames_opt=(/'Z_T','Z_M','P_T','P_M','TT','HU','P0'/))

      ! initialize statevector_noZnoP on analysis grid
      call gsv_allocate(statevector_noZnoP, tim_nstepobsinc, hco_anl, vco_anl, &
                        dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                        allocHeightSfc_opt=.true., hInterpolateDegree_opt='LINEAR', &
                        varNames_opt=(/'TT','HU','P0'/))
      write(*,*) 'gvt_setupTrials: statevector_noZnoP allocated'

      ! read trial files using default horizontal interpolation degree
      call gsv_readTrials( statevector_noZnoP )  ! IN/OUT

      ! copy the statevectors
      call gsv_copy( statevector_noZnoP, stateVectorTrialHeight, allowMismatch_opt=.true. )

      call gsv_deallocate(statevector_noZnoP)

      ! do height/P calculation of the grid
      call PsfcToP_nl( stateVectorTrialHeight )
      call tt2phi( stateVectorTrialHeight )

      heightTrialsInitialized = .true.
    case default
      call utl_abort('gvt_setupTrials: unknown variable ='//trim(varName))
    end select

  end subroutine gvt_setupTrials


  subroutine gvt_transform_gsv(statevector, transform, statevectorOut_opt,  &
                               stateVectorRef_opt)
    implicit none
   
    ! Arguments
    type(struct_gsv)           :: statevector
    character(len=*), intent(in) :: transform
    type(struct_gsv), optional :: statevectorOut_opt
    type(struct_gsv), optional :: statevectorRef_opt


    select case(trim(transform))

    case ('UVtoVortDiv')
      if (present(statevectorOut_opt)) then
        call utl_abort('gvt_transform: for UVtoVortDiv, the option statevectorOut_opt is not yet available')
      end if
      call UVtoVortDiv_gsv(statevector)
    case ('VortDivToPsiChi')
      if ( .not. gsv_varExist(statevector,'QR') .or. .not. gsv_varExist(statevector,'DD') ) then
        call utl_abort('gvt_transform: for VortDivToPsiChi, variables QR and DD must be allocated in gridstatevector')
      end if
      if (present(statevectorOut_opt)) then
        call utl_abort('gvt_transform: for VortDivToPsiChi, the option statevectorOut_opt is not yet available')
      end if
      call VortDivToPsiChi_gsv(statevector)
    case ('UVtoPsiChi')
      if ( .not. gsv_varExist(statevector,'PP') .or. .not. gsv_varExist(statevector,'CC') ) then
        call utl_abort('gvt_transform: for UVToPsiChi, variables PP and CC must be allocated in gridstatevector')
      end if
      if (present(statevectorOut_opt)) then
        call utl_abort('gvt_transform: for UVToPsiChi, the option statevectorOut_opt is not yet available')
      end if
      call UVtoPsiChi_gsv(statevector)
    case ('LQtoHU')
      if ( .not. gsv_varExist(statevector,'HU') ) then
        call utl_abort('gvt_transform: for LQtoHU, variable HU must be allocated in gridstatevector')
      end if
      if (present(statevectorOut_opt)) then
        call utl_abort('gvt_transform: for LQtoHU, the option statevectorOut_opt is not yet available')
      end if
      call LQtoHU(statevector)
    case ('HUtoLQ')
      if ( .not. gsv_varExist(statevector,'HU') ) then
        call utl_abort('gvt_transform: for HUtoLQ, variable HU must be allocated in gridstatevector')
      end if
      if (present(statevectorOut_opt)) then
        call utl_abort('gvt_transform: for HUtoLQ, the option statevectorOut_opt is not yet available')
      end if
      call HUtoLQ_gsv(statevector)
    case ('LQtoHU_tlm')
      if ( .not. gsv_varExist(statevector,'HU') ) then
        call utl_abort('gvt_transform: for LQtoHU_tlm, variable HU must be allocated in gridstatevector')
      end if
      if (present(statevectorOut_opt)) then
        call utl_abort('gvt_transform: for LQtoHU_tlm, the option statevectorOut_opt is not yet available')
      end if
      call LQtoHU_tlm(statevector, stateVectorRef_opt)
    case ('HUtoLQ_tlm')
      if ( .not. gsv_varExist(statevector,'HU') ) then
        call utl_abort('gvt_transform: for HUtoLQ_tlm, variable HU must be allocated in gridstatevector')
      end if
      if (present(statevectorOut_opt)) then
        call utl_abort('gvt_transform: for HUtoLQ_ad, the option statevectorOut_opt is not yet available')
      end if
      call HUtoLQ_tlm(statevector, stateVectorRef_opt)
    case ('LQtoHU_ad')
      if ( .not. gsv_varExist(statevector,'HU') ) then
        call utl_abort('gvt_transform: for LQtoHU_ad, variable HU must be allocated in gridstatevector')
      end if
      if (present(statevectorOut_opt)) then
        call utl_abort('gvt_transform: for LQtoHU_ad, the option statevectorOut_opt is not yet available')
      end if
      call LQtoHU_tlm(statevector, stateVectorRef_opt) ! self-adjoint
    case ('HUtoLQ_ad')
      if ( .not. gsv_varExist(statevector,'HU') ) then
        call utl_abort('gvt_transform: for HUtoLQ_ad, variable HU must be allocated in gridstatevector')
      end if
      if (present(statevectorOut_opt)) then
        call utl_abort('gvt_transform: for HUtoLQ_ad, the option statevectorOut_opt is not yet available')
      end if
      call HUtoLQ_tlm(statevector, stateVectorRef_opt) ! self-adjoint

    case ('LVIStoVIS')
      if (present(statevectorOut_opt)) then
        if ( .not. gsv_varExist(statevector,'LVIS')) then
          call utl_abort('gvt_transform: for LVIStoVIS, variable LVIS must be allocated in gridstatevector IN')
        end if
        if ( .not. gsv_varExist(statevectorOut_opt,'VIS')) then
          call utl_abort('gvt_transform: for LVIStoVIS, variable VIS must be allocated in gridstatevector OUT')
        end if
        call LVIStoVIS(statevector, statevectorOut_opt=statevectorOut_opt)
      else
        if ( .not. gsv_varExist(statevector,'VIS') .or. .not. gsv_varExist(statevector,'LVIS') ) then
          call utl_abort('gvt_transform: for LVIStoVIS, variables LVIS and VIS must be allocated in gridstatevector')
        end if
        call LVIStoVIS(statevector)
      end if
    case ('TTHUtoHeight_nl')
      if ( .not. gsv_varExist(statevector,'TT')  ) then
        call utl_abort('gvt_transform: for TTHUtoHeight_nl, variable TT must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'HU')  ) then
        call utl_abort('gvt_transform: for TTHUtoHeight_nl, variable HU must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'P0')  ) then
        call utl_abort('gvt_transform: for TTHUtoHeight_nl, variable P0 must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'Z_T')  ) then
        call utl_abort('gvt_transform: for TTHUtoHeight_nl, variable Z_T must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'Z_M')  ) then
        call utl_abort('gvt_transform: for TTHUtoHeight_nl, variable Z_M must be allocated in gridstatevector')
      end if
      call TTHUtoHeight_nl(statevector)

    case ('TTHUtoHeight_tl')
      if ( .not. gsv_varExist(statevector,'TT')  ) then
        call utl_abort('gvt_transform: for TTHUtoHeight_tl, variable TT must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'HU')  ) then
        call utl_abort('gvt_transform: for TTHUtoHeight_tl, variable HU must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'P0')  ) then
        call utl_abort('gvt_transform: for TTHUtoHeight_tl, variable P0 must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'Z_T')  ) then
        call utl_abort('gvt_transform: for TTHUtoHeight_tl, variable Z_T must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'Z_M')  ) then
        call utl_abort('gvt_transform: for TTHUtoHeight_tl, variable Z_M must be allocated in gridstatevector')
      end if
      call TTHUtoHeight_tl(statevector)

    case ('TTHUtoHeight_ad')
      if ( .not. gsv_varExist(statevector,'TT')  ) then
        call utl_abort('gvt_transform: for TTHUtoHeight_ad, variable TT must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'HU')  ) then
        call utl_abort('gvt_transform: for TTHUtoHeight_ad, variable HU must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'P0')  ) then
        call utl_abort('gvt_transform: for TTHUtoHeight_ad, variable P0 must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'Z_T')  ) then
        call utl_abort('gvt_transform: for TTHUtoHeight_ad, variable Z_T must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'Z_M')  ) then
        call utl_abort('gvt_transform: for TTHUtoHeight_ad, variable Z_M must be allocated in gridstatevector')
      end if
      call TTHUtoHeight_ad(statevector)

    case ('PsfcToP_nl')
      if ( .not. gsv_varExist(statevector,'P_T')  ) then
        call utl_abort('gvt_transform: for PsfcToP_nl, variable P_T must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'P_M')  ) then
        call utl_abort('gvt_transform: for PsfcToP_nl, variable P_M must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'P0')  ) then
        call utl_abort('gvt_transform: for PsfcToP_nl, variable P0 must be allocated in gridstatevector')
      end if
      call PsfcToP_nl(statevector)

    case ('PsfcToP_tl')
      if ( .not. gsv_varExist(statevector,'P_T')  ) then
        call utl_abort('gvt_transform: for PsfcToP_tl, variable P_T must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'P_M')  ) then
        call utl_abort('gvt_transform: for PsfcToP_tl, variable P_M must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'P0')  ) then
        call utl_abort('gvt_transform: for PsfcToP_tl, variable P0 must be allocated in gridstatevector')
      end if
      call PsfcToP_tl(statevector)

    case ('PsfcToP_ad')
      if ( .not. gsv_varExist(statevector,'P_T')  ) then
        call utl_abort('gvt_transform: for PsfcToP_ad, variable P_T must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'P_M')  ) then
        call utl_abort('gvt_transform: for PsfcToP_ad, variable P_M must be allocated in gridstatevector')
      end if
      if ( .not. gsv_varExist(statevector,'P0')  ) then
        call utl_abort('gvt_transform: for PsfcToP_ad, variable P0 must be allocated in gridstatevector')
      end if
      call PsfcToP_ad(statevector)

    case default
      write(*,*)
      write(*,*) 'Unsupported function : ', trim(transform)
      call utl_abort('gvt_transform')
    end select

  end subroutine gvt_transform_gsv


  subroutine gvt_transform_ens(ens,transform)
    implicit none
   
    type(struct_ens) :: ens
 
    character(len=*), intent(in) :: transform

    select case(trim(transform))
    case ('HUtoLQ')
      call HUtoLQ_ens(ens)
    case ('UVtoPsiChi')
      call UVtoPsiChi_ens(ens)
    case ('UVtoVortDiv')
      call UVtoVortDiv_ens(ens)
    case default
      call utl_abort('gvt_transform_ens: Unsupported function '//trim(transform))
    end select

  end subroutine gvt_transform_ens

  !--------------------------------------------------------------------------
  ! gvt_getStateVectorTrial
  !--------------------------------------------------------------------------
  function gvt_getStateVectorTrial(varName) result(statevector_ptr)
    implicit none

    ! arguments
    character(len=*), intent(in) :: varName

    ! local
    type(struct_gsv), pointer  :: statevector_ptr

    select case ( trim(varName) )
    case ('HU')
      if ( .not. huTrialsInitialized ) call gvt_setupTrials('HU')
      statevector_ptr => stateVectorTrialHU
      huTrialsInitialized = .true.

    case ('height')
      if ( .not. heightTrialsInitialized ) call gvt_setupTrials('height')
      statevector_ptr => stateVectorTrialHeight
      heightTrialsInitialized = .true.

    case default
      call utl_abort('gvt_getStateVectorTrial: unknown variable ='//trim(varName))
    end select

  end function gvt_getStateVectorTrial

  !--------------------------------------------------------------------------
  ! LQtoHU
  !--------------------------------------------------------------------------
  subroutine LQtoHU(statevector)
    implicit none

    type(struct_gsv) :: statevector
    integer :: lonIndex,latIndex,levIndex,stepIndex

    real(8), pointer :: hu_ptr(:,:,:,:), lq_ptr(:,:,:,:)

    hu_ptr => gsv_getField_r8(statevector,'HU')
    lq_ptr => gsv_getField_r8(statevector,'HU')

    do stepIndex = 1, statevector%numStep
      !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
      do levIndex = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname('HU'))
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do lonIndex = statevector%myLonBeg, statevector%myLonEnd
            hu_ptr(lonIndex,latIndex,levIndex,stepIndex) = exp(lq_ptr(lonIndex,latIndex,levIndex,stepIndex))
          end do
        end do
      end do
      !$OMP END PARALLEL DO
    end do

  end subroutine LQtoHU


  subroutine HUtoLQ_gsv(statevector)
    implicit none

    type(struct_gsv)    :: statevector
    integer :: lonIndex,latIndex,levIndex,stepIndex

    real(4), pointer :: hu_ptr_r4(:,:,:,:), lq_ptr_r4(:,:,:,:)
    real(8), pointer :: hu_ptr_r8(:,:,:,:), lq_ptr_r8(:,:,:,:)

    if ( statevector%dataKind == 8 ) then

      hu_ptr_r8   => gsv_getField_r8(statevector,'HU')
      lq_ptr_r8   => gsv_getField_r8(statevector,'HU')

      do stepIndex = 1, statevector%numStep
        !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
        do levIndex = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname('HU'))
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd
              lq_ptr_r8(lonIndex,latIndex,levIndex,stepIndex) = log(max(hu_ptr_r8(lonIndex,latIndex,levIndex,stepIndex),gsv_rhumin))
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do

    else

      hu_ptr_r4   => gsv_getField_r4(statevector,'HU')
      lq_ptr_r4   => gsv_getField_r4(statevector,'HU')

      do stepIndex = 1, statevector%numStep
        !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
        do levIndex = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname('HU'))
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd
              lq_ptr_r4(lonIndex,latIndex,levIndex,stepIndex) = log(max(hu_ptr_r4(lonIndex,latIndex,levIndex,stepIndex),real(gsv_rhumin,4)))
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do

    end if

  end subroutine HUtoLQ_gsv


  subroutine HUtoLQ_ens(ens)
    implicit none

    type(struct_ens) :: ens

    integer :: lonIndex, latIndex, levIndex, stepIndex, memberIndex
    integer :: myLatBeg, myLatEnd, myLonBeg, myLonEnd

    character(len=4) :: varName

    real(4), pointer :: ptr4d_r4(:,:,:,:)

    call ens_getLatLonBounds(ens, myLonBeg, myLonEnd, myLatBeg, myLatEnd)

    do levIndex = 1, ens_getNumK(ens)

      varName = ens_getVarNameFromK(ens,levIndex)
      if (varName /= 'HU') cycle

      ptr4d_r4 => ens_getOneLev_r4(ens,levIndex)
    
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          do stepIndex = 1, ens_getNumStep(ens)
            do memberIndex = 1, ens_getNumMembers(ens)
              ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex) = log(max(ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex),real(gsv_rhumin,4)))
            end do
          end do
        end do
      end do

    end do

  end subroutine HUtoLQ_ens


  subroutine LQtoHU_tlm(statevector, stateVectorRef_opt)
    implicit none

    ! Arguments
    type(struct_gsv)           :: statevector
    type(struct_gsv), optional :: statevectorRef_opt

    ! Locals
    integer :: lonIndex,latIndex,levIndex,stepIndex
    real(8), pointer :: hu_ptr(:,:,:,:), lq_ptr(:,:,:,:), hu_trial(:,:,:,:)

    if (present(statevectorRef_opt)) then
      hu_trial => gsv_getField_r8(stateVectorRef_opt,'HU')
    else
      if ( .not. huTrialsInitialized ) call gvt_setupTrials('HU')
      hu_trial => gsv_getField_r8(stateVectorTrialHU,'HU')
    end if

    hu_ptr   => gsv_getField_r8(statevector      ,'HU')
    lq_ptr   => gsv_getField_r8(statevector      ,'HU')

    do stepIndex = 1, statevector%numStep
      !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
      do levIndex = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname('HU'))
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do lonIndex = statevector%myLonBeg, statevector%myLonEnd       
            hu_ptr(lonIndex,latIndex,levIndex,stepIndex) =  lq_ptr(lonIndex,latIndex,levIndex,stepIndex)*max(hu_trial(lonIndex,latIndex,levIndex,stepIndex),MPC_MINIMUM_HU_R8)
          end do
        end do
      end do
      !$OMP END PARALLEL DO
    end do

  end subroutine LQtoHU_tlm


  subroutine HUtoLQ_tlm(statevector, stateVectorRef_opt)
    implicit none

    ! Arguments
    type(struct_gsv)           :: statevector
    type(struct_gsv), optional :: statevectorRef_opt

    ! Locals
    integer :: lonIndex,latIndex,levIndex,stepIndex
    real(8), pointer :: hu_ptr(:,:,:,:), lq_ptr(:,:,:,:), hu_trial(:,:,:,:)

    if (present(statevectorRef_opt)) then
      hu_trial => gsv_getField_r8(stateVectorRef_opt,'HU')
    else
      if ( .not. huTrialsInitialized ) call gvt_setupTrials('HU')
      hu_trial => gsv_getField_r8(stateVectorTrialHU,'HU')
    end if

    hu_ptr   => gsv_getField_r8(statevector      ,'HU')
    lq_ptr   => gsv_getField_r8(statevector      ,'HU')

    do stepIndex = 1, statevector%numStep
      !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
      do levIndex = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname('HU'))
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do lonIndex = statevector%myLonBeg, statevector%myLonEnd
            lq_ptr(lonIndex,latIndex,levIndex,stepIndex) = hu_ptr(lonIndex,latIndex,levIndex,stepIndex)/max(hu_trial(lonIndex,latIndex,levIndex,stepIndex),MPC_MINIMUM_HU_R8)
          end do
        end do
      end do
      !$OMP END PARALLEL DO
    end do

  end subroutine HUtoLQ_tlm


  subroutine LVIStoVIS(statevector_in, statevectorOut_opt)
    implicit none

    type(struct_gsv) :: statevector_in
    type(struct_gsv), optional :: statevectorOut_opt
    integer :: lonIndex,latIndex,levIndex,stepIndex

    real(4), pointer :: vis_ptr_r4(:,:,:,:), lvis_ptr_r4(:,:,:,:)
    real(8), pointer :: vis_ptr_r8(:,:,:,:), lvis_ptr_r8(:,:,:,:)

    if ( statevector_in%dataKind == 8 ) then

      if (present(statevectorOut_opt)) then
        vis_ptr_r8  => gsv_getField_r8(statevectorOut_opt,'VIS')
      else
        vis_ptr_r8  => gsv_getField_r8(statevector_in,'VIS')
      end if
      lvis_ptr_r8 => gsv_getField_r8(statevector_in,'LVIS')
      
      do stepIndex = 1, statevector_in%numStep
        !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
        do levIndex = 1, gsv_getNumLev(statevector_in,vnl_varLevelFromVarname('LVIS'))
          do latIndex = statevector_in%myLatBeg, statevector_in%myLatEnd
            do lonIndex = statevector_in%myLonBeg, statevector_in%myLonEnd
              vis_ptr_r8(lonIndex,latIndex,levIndex,stepIndex) = exp(lvis_ptr_r8(lonIndex,latIndex,levIndex,stepIndex))
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do

    else

      if (present(statevectorOut_opt)) then
        vis_ptr_r4  => gsv_getField_r4(statevectorOut_opt,'VIS')
      else
        vis_ptr_r4  => gsv_getField_r4(statevector_in,'VIS')
      end if
      lvis_ptr_r4 => gsv_getField_r4(statevector_in,'LVIS')

      do stepIndex = 1, statevector_in%numStep
        !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
        do levIndex = 1, gsv_getNumLev(statevector_in,vnl_varLevelFromVarname('LVIS'))
          do latIndex = statevector_in%myLatBeg, statevector_in%myLatEnd
            do lonIndex = statevector_in%myLonBeg, statevector_in%myLonEnd
              vis_ptr_r4(lonIndex,latIndex,levIndex,stepIndex) = exp(lvis_ptr_r4(lonIndex,latIndex,levIndex,stepIndex))
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do

    end if

  end subroutine LVIStoVIS

  ! TTHUtoHeight_nl
  !--------------------------------------------------------------------------
  subroutine TTHUtoHeight_nl(statevector)
    implicit none

    type(struct_gsv)    :: statevector

    call tt2phi(statevector)

  end subroutine TTHUtoHeight_nl

  !--------------------------------------------------------------------------
  ! TTHUtoHeight_tl
  !--------------------------------------------------------------------------
  subroutine TTHUtoHeight_tl(statevector)
    implicit none

    type(struct_gsv)    :: statevector

    if ( .not. heightTrialsInitialized ) call gvt_setupTrials('height')

    call tt2phi_tl(statevector, stateVectorTrialHeight)

  end subroutine TTHUtoHeight_tl

  !--------------------------------------------------------------------------
  ! TTHUtoHeight_ad
  !--------------------------------------------------------------------------
  subroutine TTHUtoHeight_ad(statevector)
    implicit none

    type(struct_gsv)    :: statevector

    if ( .not. heightTrialsInitialized ) call gvt_setupTrials('height')

    call tt2phi_ad(statevector,stateVectorTrialHeight)

  end subroutine TTHUtoHeight_ad

  !--------------------------------------------------------------------------
  ! PsfcToP_nl
  !--------------------------------------------------------------------------
  subroutine PsfcToP_nl(statevector)
    implicit none

    type(struct_gsv)    :: statevector

    if ( statevector%dataKind == 8 ) then
      call calcpressure_nl_r8(statevector)
    else if ( statevector%dataKind == 4 ) then
      call calcpressure_nl_r4(statevector)
    end if

  end subroutine PsfcToP_nl

  !--------------------------------------------------------------------------
  ! PsfcToP_tl
  !--------------------------------------------------------------------------
  subroutine PsfcToP_tl(statevector)
    implicit none

    type(struct_gsv)    :: statevector

    if ( .not. heightTrialsInitialized ) call gvt_setupTrials('height')

    call calcpressure_tl(statevector,stateVectorTrialHeight)

  end subroutine PsfcToP_tl

  !--------------------------------------------------------------------------
  ! PsfcToP_ad
  !--------------------------------------------------------------------------
  subroutine PsfcToP_ad(statevector)
    implicit none

    type(struct_gsv)    :: statevector

    if ( .not. heightTrialsInitialized ) call gvt_setupTrials('height')

    call calcpressure_ad(statevector,stateVectorTrialHeight)

  end subroutine PsfcToP_ad

  !--------------------------------------------------------------------------
  ! UVtoVortDiv_gsv
  !--------------------------------------------------------------------------
  subroutine UVtoVortDiv_gsv(statevector)
    implicit none
   
    type(struct_gsv) :: statevector
    integer :: stepIndex
 
    real(8), pointer :: uu_ptr(:,:,:,:), vv_ptr(:,:,:,:)
    real(8), pointer :: qr_ptr(:,:,:,:), dd_ptr(:,:,:,:)

    integer :: nlev_M

    uu_ptr => gsv_getField_r8(statevector,'UU')
    vv_ptr => gsv_getField_r8(statevector,'VV')
    if (gsv_varExist(statevector,'QR') .and. &
        gsv_varExist(statevector,'DD')) then
      qr_ptr => gsv_getField_r8(statevector,'QR')
      dd_ptr => gsv_getField_r8(statevector,'DD')
    else
      qr_ptr => gsv_getField_r8(statevector,'UU')
      dd_ptr => gsv_getField_r8(statevector,'VV')
    end if

    nlev_M =  gsv_getNumLev  (statevector,'MM')
    
    if (  statevector%hco%global ) then
      call utl_abort('UVtoVortDiv_gsv: global mode not available')
    else
      do stepIndex = 1, statevector%numStep
        call agd_UVToVortDiv(qr_ptr(:,:,:,stepIndex), dd_ptr(:,:,:,stepIndex), & ! OUT
                             uu_ptr(:,:,:,stepIndex), vv_ptr(:,:,:,stepIndex), & ! IN
                             nlev_M )                                            ! IN          
      end do
    end if

  end subroutine UVtoVortDiv_gsv


  subroutine vortDivToPsiChi_gsv(statevector)
    implicit none
   
    type(struct_gsv) :: statevector
    integer :: stepIndex

    logical, save :: firstTime = .true.

    real(8), pointer :: qr_ptr(:,:,:,:), dd_ptr(:,:,:,:)
    real(8), pointer :: pp_ptr(:,:,:,:), cc_ptr(:,:,:,:)

    type(struct_lst), save :: lst_lapl   ! Spectral transform Parameters for Vort/Div -> Psi/Chi
    integer :: nlev_M
    integer :: nTrunc

    nTrunc = max(statevector%hco%ni,statevector%hco%nj)/4 + 1

    if (gsv_varExist(statevector,'QR') .and. &
        gsv_varExist(statevector,'DD')) then
      qr_ptr => gsv_getField_r8(statevector,'QR')
      dd_ptr => gsv_getField_r8(statevector,'DD')
    else
      qr_ptr => gsv_getField_r8(statevector,'UU')
      dd_ptr => gsv_getField_r8(statevector,'VV')
    end if

    if (gsv_varExist(statevector,'PP') .and. &
        gsv_varExist(statevector,'CC')) then
      pp_ptr => gsv_getField_r8(statevector,'PP')
      cc_ptr => gsv_getField_r8(statevector,'CC')
    else
      pp_ptr => gsv_getField_r8(statevector,'UU')
      cc_ptr => gsv_getField_r8(statevector,'VV')
    end if

    nlev_M =  gsv_getNumLev  (statevector,'MM')

    if ( statevector%hco%global ) then
      call utl_abort('vortDivToPsiChi_gsv: global mode not available')
    else
      if (firstTime) then
        call lst_Setup(lst_lapl,                                & ! OUT
                       statevector%hco%ni, statevector%hco%nj,  & ! IN
                       statevector%hco%dlon, nTrunc,            & ! IN
                       'LatLonMN', nlev_M )                       ! IN
        firstTime = .false.
      end if

      do stepIndex = 1, statevector%numStep
          
        pp_ptr(:,:,:,stepIndex) = qr_ptr(:,:,:,stepIndex)
        cc_ptr(:,:,:,stepIndex) = dd_ptr(:,:,:,stepIndex)
         
        ! Vort -> Psi
        call lst_Laplacian(lst_lapl,                & ! IN
                           pp_ptr(:,:,:,stepIndex), & ! INOUT
                           'Inverse', nlev_M )        ! IN

        ! Div -> Chi
        call lst_Laplacian(lst_lapl,                & ! IN
                           cc_ptr(:,:,:,stepIndex), & ! INOUT
                           'Inverse', nlev_M )        ! IN

      end do

    end if

  end subroutine VortDivToPsiChi_gsv


  subroutine UVtoPsiChi_gsv(statevector)
    implicit none
   
    type(struct_gsv) :: statevector

    integer :: stepIndex 
    real(8), pointer :: uu_ptr(:,:,:,:), vv_ptr(:,:,:,:)
    real(8), pointer :: psi_ptr(:,:,:,:), chi_ptr(:,:,:,:)
    real(8), allocatable :: gridState(:,:,:), spectralState(:,:,:)
    real(8) :: dla2
    integer :: nlev_M, levIndex
    integer :: ila_mpiglobal, ila_mpilocal

    ! spectral transform configuration (saved)
    integer, save :: gstID = -1
    integer, save :: nla_mpilocal, maxMyNla, ntrunc
    integer, save :: mymBeg,mymEnd,mymSkip,mymCount
    integer, save :: mynBeg,mynEnd,mynSkip,mynCount
    integer, save, pointer :: ilaList_mpiglobal(:), ilaList_mpilocal(:)

    write(*,*) 'UVtoPsiChi_gsv: starting'
    call flush(6)

    if ( .not. statevector%hco%global ) then

      call UVtoVortDiv_gsv    (statevector)
      call vortDivToPsiChi_gsv(statevector)

    else

      uu_ptr  => gsv_getField_r8(statevector,'UU')
      vv_ptr  => gsv_getField_r8(statevector,'VV')
      psi_ptr => gsv_getField_r8(statevector,'PP')
      chi_ptr => gsv_getField_r8(statevector,'CC')
      nlev_M = gsv_getNumLev(statevector,'MM')

      if ( gstID < 0 ) then
        !ntrunc = statevector%nj
        ntrunc = 180
        gstID = gst_setup(statevector%ni,statevector%nj,ntrunc,2*nlev_M)
        call mpivar_setup_m(ntrunc,mymBeg,mymEnd,mymSkip,mymCount)
        call mpivar_setup_n(ntrunc,mynBeg,mynEnd,mynSkip,mynCount)
        call gst_ilaList_mpiglobal(ilaList_mpiglobal,nla_mpilocal,maxMyNla,gstID,mymBeg,mymEnd,mymSkip,mynBeg,mynEnd,mynSkip)
        call gst_ilaList_mpilocal(ilaList_mpilocal,gstID,mymBeg,mymEnd,mymSkip,mynBeg,mynEnd,mynSkip)
      end if

      dla2   = dble(ra)*dble(ra)
      allocate(gridState(statevector%lonPerPE,statevector%latPerPE,2*nlev_M))
      allocate(spectralState(nla_mpilocal,2,2*nlev_M))

      do stepIndex = 1, statevector%numStep

        gridState(:,:,1:nlev_M)            = uu_ptr(:,:,:,stepIndex)
        gridState(:,:,(nlev_M+1):2*nlev_M) = vv_ptr(:,:,:,stepIndex)

        call gst_setID(gstID)
        call gst_gdsp(spectralState,gridState,nlev_M)

        do levIndex = 1, 2*nlev_M
          do ila_mpilocal = 1, nla_mpilocal
            ila_mpiglobal = ilaList_mpiglobal(ila_mpilocal)
            spectralState(ila_mpilocal,:,levIndex) =   &
                             spectralState(ila_mpilocal,:,levIndex) *   &
                             dla2*gst_getR1snp1(ila_mpiglobal)
          enddo
        enddo

        call gst_speree(spectralState,gridState)

        psi_ptr(:,:,:,stepIndex) = gridState(:,:,1:nlev_M)
        chi_ptr(:,:,:,stepIndex) = gridState(:,:,(nlev_M+1):2*nlev_M)

      end do

      write(*,*) 'deallocate'; call flush(6)
      deallocate(gridState)
      deallocate(spectralState)

    end if

    write(*,*) 'UVtoPsiChi_gsv: finished'
    call flush(6)

  end subroutine UVtoPsiChi_gsv
  

  subroutine UVtoPsiChi_ens(ens)
    implicit none
   
    type(struct_ens) :: ens

    type(struct_hco), pointer :: hco_ens => null()
    type(struct_gsv) :: gridStateVector_oneMember

    integer :: memberIndex

    write(*,*)
    write(*,*) 'gvt_UVtoPsiChi_ens: starting'

    hco_ens => ens_getHco(ens)

    if (hco_ens%global ) then
      call utl_abort('gvt_UVtoPsiChi_ens: global mode not yet available')
    end if

    !
    !- 1.  Create a working stateVector
    !
    call gsv_allocate(gridStateVector_oneMember, 1, hco_ens, ens_getVco(ens), &
                      varNames_opt=(/'UU','VV'/), datestamp_opt=tim_getDatestamp(), &
                      mpi_local_opt=.true., dataKind_opt=8)

    !
    !- 2.  Loop on members
    !
    do memberIndex = 1, ens_getNumMembers(ens)

      !- 2.1 Copy to a stateVector
      call ens_copyMember(ens, gridStateVector_oneMember, memberIndex)

      !- 2.2 Do the transform
      call UVtoPsiChi_gsv(gridStateVector_oneMember)

      !- 2.3 Put the result back in the input ensembleStateVector
      call ens_insertMember(ens, gridStateVector_oneMember, memberIndex)
    end do

    !
    !- 3. Cleaning
    !
    call gsv_deallocate(gridStateVector_oneMember)

    write(*,*) 'gvt_UVtoPsiChi_ens: finished'

  end subroutine UVtoPsiChi_ens


  subroutine UVtoVortDiv_ens(ens)
    implicit none
   
    type(struct_ens) :: ens

    type(struct_hco), pointer :: hco_ens => null()
    type(struct_gsv) :: gridStateVector_oneMember

    integer :: memberIndex

    write(*,*)
    write(*,*) 'gvt_UVtoVortDiv_ens: starting'

    hco_ens => ens_getHco(ens)

    if (hco_ens%global ) then
      call utl_abort('gvt_UVtoVortDiv_ens: global mode not yet available')
    end if

    !
    !- 1.  Create a working stateVector
    !
    call gsv_allocate(gridStateVector_oneMember, 1, hco_ens, ens_getVco(ens), &
                      varNames_opt=(/'UU','VV'/), datestamp_opt=tim_getDatestamp(), &
                      mpi_local_opt=.true., dataKind_opt=8)

    !
    !- 2.  Loop on members
    !
    do memberIndex = 1, ens_getNumMembers(ens)

      !- 2.1 Copy to a stateVector
      call ens_copyMember(ens, gridStateVector_oneMember, memberIndex)

      !- 2.2 Do the transform
      call UVtoVortDiv_gsv(gridStateVector_oneMember)

      !- 2.3 Put the result back in the input ensembleStateVector
      call ens_insertMember(ens, gridStateVector_oneMember, memberIndex)
    end do

    !
    !- 3. Cleaning
    !
    call gsv_deallocate(gridStateVector_oneMember)

    write(*,*) 'gvt_UVtoVortDiv_ens: finished'

  end subroutine UVtoVortDiv_ens

  !--------------------------------------------------------------------------
  ! calcpressure_nl_r8
  !--------------------------------------------------------------------------
  subroutine calcPressure_nl_r8(statevector, beSilent_opt)
    !
    !:Purpose: double-precision calculation of the pressure on the grid.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: statevector ! statevector that will contain P_T/P_M
    logical, optional               :: beSilent_opt

    ! Locals:
    real(kind=8), allocatable   :: Psfc(:,:)
    real(kind=8), pointer       :: Pressure_out(:,:,:) 
    real(kind=8), pointer       :: field_Psfc(:,:,:,:)
    integer                     :: status, stepIndex, numStep
    logical                     :: beSilent
    real(8), pointer            :: P_T(:,:,:,:)
    real(8), pointer            :: P_M(:,:,:,:)

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .true.
    end if

    if ( .not. beSilent ) write(*,*) 'calcPressure_nl_r8: computing pressure on staggered or UNstaggered levels'

    if ( .not. gsv_varExist(statevector,'P_T') .or. .not. gsv_varExist(statevector,'P_M') .or. .not. gsv_varExist(statevector,'P0')) then
      call utl_abort('calcPressure_nl_r8: P_T/P_M/P0 do not exist in statevector!')
    end if

    nullify(P_T)
    nullify(P_M)

    P_T => gsv_getField_r8(statevector,'P_T')
    P_M => gsv_getField_r8(statevector,'P_M')

    allocate(Psfc(statevector%myLonBeg:statevector%myLonEnd, &
                  statevector%myLatBeg:statevector%myLatEnd))
    field_Psfc => gsv_getField_r8(statevector,'P0')
    numStep = statevector%numStep

    do stepIndex = 1, numStep

      Psfc(:,:) = field_Psfc(:,:,1,stepIndex)

      ! P_M
      nullify(Pressure_out)
      status = vgd_levels(statevector%vco%vgrid, &
                        ip1_list=statevector%vco%ip1_M, &
                        levels=Pressure_out, &
                        sfc_field=Psfc, &
                        in_log=.false.)
      if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_levels')
      P_M(:,:,:,stepIndex) = Pressure_out(:,:,:)
      deallocate(Pressure_out)

      ! P_T
      nullify(Pressure_out)
      status = vgd_levels(statevector%vco%vgrid, &
                        ip1_list=statevector%vco%ip1_T, &
                        levels=Pressure_out, &
                        sfc_field=Psfc, &
                        in_log=.false.)
      if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_levels')
      P_T(:,:,:,stepIndex) = Pressure_out(:,:,:)
      deallocate(Pressure_out)

      if ( .not. beSilent .and. stepIndex == 1 ) then
        write(*,*) 'stepIndex=',stepIndex, ',P_M='
        write(*,*) P_M(statevector%myLonBeg,statevector%myLatBeg,:,stepIndex)
        write(*,*) 'stepIndex=',stepIndex, ',P_T='
        write(*,*) P_T(statevector%myLonBeg,statevector%myLatBeg,:,stepIndex)
      end if

    end do

    deallocate(Psfc)

    if ( .not. beSilent ) write(*,*) 'calcPressure_nl_r8: END'

  end subroutine calcPressure_nl_r8

  !--------------------------------------------------------------------------
  ! calcpressure_nl_r4
  !--------------------------------------------------------------------------
  subroutine calcPressure_nl_r4(statevector, beSilent_opt)
    !
    !:Purpose: single-precision calculation of the pressure on the grid.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: statevector ! statevector that will contain P_T/P_M
    logical, optional               :: beSilent_opt

    ! Locals:
    real(kind=4), allocatable   :: Psfc(:,:)
    real(kind=4), pointer       :: Pressure_out(:,:,:) 
    real(kind=4), pointer       :: field_Psfc(:,:,:,:)
    integer                     :: status, stepIndex, numStep
    logical                     :: beSilent
    real(4), pointer            :: P_T(:,:,:,:)
    real(4), pointer            :: P_M(:,:,:,:)

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .true.
    end if

    if ( .not. beSilent ) write(*,*) 'calcPressure_nl_r4: computing pressure on staggered or UNstaggered levels'

    if ( .not. gsv_varExist(statevector,'P_T') .or. .not. gsv_varExist(statevector,'P_M') .or. .not. gsv_varExist(statevector,'P0')) then
      call utl_abort('calcPressure_nl_r4: P_T/P_M/P0 do not exist in statevector!')
    end if

    nullify(P_T)
    nullify(P_M)

    P_T => gsv_getField_r4(statevector,'P_T')
    P_M => gsv_getField_r4(statevector,'P_M')

    allocate(Psfc(statevector%myLonBeg:statevector%myLonEnd, &
                  statevector%myLatBeg:statevector%myLatEnd))
    field_Psfc => gsv_getField_r4(statevector,'P0')
    numStep = statevector%numStep

    do stepIndex = 1, numStep

      Psfc(:,:) = field_Psfc(:,:,1,stepIndex)

      ! P_T
      nullify(Pressure_out)
      status = vgd_levels(statevector%vco%vgrid, &
                        ip1_list=statevector%vco%ip1_M, &
                        levels=Pressure_out, &
                        sfc_field=Psfc, &
                        in_log=.false.)
      if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_levels')
      P_M(:,:,:,stepIndex) = Pressure_out(:,:,:)
      deallocate(Pressure_out)

      ! P_M
      nullify(Pressure_out)
      status = vgd_levels(statevector%vco%vgrid, &
                        ip1_list=statevector%vco%ip1_T, &
                        levels=Pressure_out, &
                        sfc_field=Psfc, &
                        in_log=.false.)
      if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_levels')
      P_T(:,:,:,stepIndex) = Pressure_out(:,:,:)
      deallocate(Pressure_out)

      if ( .not. beSilent .and. stepIndex == 1 ) then
        write(*,*) 'stepIndex=',stepIndex, ',P_M='
        write(*,*) P_M(statevector%myLonBeg,statevector%myLatBeg,:,stepIndex)
        write(*,*) 'stepIndex=',stepIndex, ',P_T='
        write(*,*) P_T(statevector%myLonBeg,statevector%myLatBeg,:,stepIndex)
      end if

    end do

    deallocate(Psfc)

  end subroutine calcPressure_nl_r4

  !--------------------------------------------------------------------------
  ! calcpressure_tl
  !--------------------------------------------------------------------------
  subroutine calcPressure_tl(statevector, statevector_trial, beSilent_opt)
    !
    !:Purpose: calculation of the pressure increment on the grid.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: statevector       ! statevector that will contain the P_T/P_M increments
    type(struct_gsv), intent(in)    :: statevector_trial ! statevector that has the Psfc
    logical, optional               :: beSilent_opt

    ! Locals:
    real(8), allocatable  :: Psfc(:,:)
    real(8), pointer      :: delPsfc(:,:,:,:)
    real(8), pointer      :: field_Psfc(:,:,:,:)
    real(8), pointer      :: delP_T(:,:,:,:)
    real(8), pointer      :: delP_M(:,:,:,:)
    real(8), pointer      :: dP_dPsfc_T(:,:,:)
    real(8), pointer      :: dP_dPsfc_M(:,:,:)
    integer               :: status, stepIndex,lonIndex,latIndex
    integer               :: lev_M, lev_T, nlev_T, nlev_M, numStep
    logical               :: beSilent

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .true.
    end if

    if ( .not. beSilent ) write(*,*) 'calcPressure_tl: computing delP_T/delP_M on the gridstatevector'

    nullify(delPsfc)
    nullify(field_Psfc)
    nullify(delP_T)
    nullify(delP_M)
    nullify(dP_dPsfc_T)
    nullify(dP_dPsfc_M)

    delP_T => gsv_getField_r8(statevector,'P_T')
    delP_M => gsv_getField_r8(statevector,'P_M')
    delPsfc => gsv_getField_r8(statevector,'P0')
    field_Psfc => gsv_getField_r8(statevector_trial,'P0')

    nlev_T = gsv_getNumLev(statevector,'TH')
    nlev_M = gsv_getNumLev(statevector,'MM')
    numStep = statevector%numstep

    allocate(Psfc(statevector%myLonBeg:statevector%myLonEnd, &
                  statevector%myLatBeg:statevector%myLatEnd))

    do stepIndex = 1, numStep

      Psfc(:,:) = field_Psfc(:,:,1,stepIndex)

      ! dP_dPsfc_M
      nullify(dP_dPsfc_M)
      status = vgd_dpidpis(statevector%vco%vgrid, &
                           statevector%vco%ip1_M, &
                           dP_dPsfc_M, &
                           Psfc)
      if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_dpidpis')
      ! calculate delP_M
      do lev_M = 1, nlev_M
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do lonIndex = statevector%myLonBeg, statevector%myLonEnd
            delP_M(lonIndex,latIndex,lev_M,stepIndex) = dP_dPsfc_M(lonIndex-statevector%myLonBeg+1,latIndex-statevector%myLatBeg+1,lev_M) * delPsfc(lonIndex,latIndex,1,stepIndex)
          end do
        end do
      end do
      deallocate(dP_dPsfc_M)

      ! dP_dPsfc_T
      nullify(dP_dPsfc_T)
      status = vgd_dpidpis(statevector%vco%vgrid, &
                           statevector%vco%ip1_T, &
                           dP_dPsfc_T, &
                           Psfc)
      if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_dpidpis')
      ! calculate delP_T
      do lev_T = 1, nlev_T
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do lonIndex = statevector%myLonBeg, statevector%myLonEnd
            delP_T(lonIndex,latIndex,lev_T,stepIndex) = dP_dPsfc_T(lonIndex-statevector%myLonBeg+1,latIndex-statevector%myLatBeg+1,lev_T) * delPsfc(lonIndex,latIndex,1,stepIndex)
          end do
        end do
      end do
      deallocate(dP_dPsfc_T)

    end do

    deallocate(Psfc)

  end subroutine calcPressure_tl

  !--------------------------------------------------------------------------
  ! calcpressure_ad
  !--------------------------------------------------------------------------
  subroutine calcPressure_ad(statevector, statevector_trial, beSilent_opt)
    !
    !:Purpose: adjoint of calculation of the pressure on the grid.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: statevector       ! statevector that will contain ncrement of Psfc.
    type(struct_gsv), intent(in)    :: statevector_trial ! statevector that has the Psfc
    logical, optional               :: beSilent_opt

    ! Locals:
    real(kind=8), allocatable   :: Psfc(:,:)
    real(kind=8), pointer       :: delPsfc(:,:,:,:)
    real(kind=8), pointer       :: field_Psfc(:,:,:,:)
    real(8), pointer            :: delP_T(:,:,:,:)
    real(8), pointer            :: delP_M(:,:,:,:)
    real(8), pointer            :: dP_dPsfc_T(:,:,:)
    real(8), pointer            :: dP_dPsfc_M(:,:,:)
    integer                     :: status, stepIndex,lonIndex,latIndex
    integer                     :: lev_M, lev_T, nlev_T, nlev_M, numStep
    logical                     :: beSilent

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .true.
    end if

    if ( .not. beSilent ) write(*,*) 'calcPressure_ad: computing delP_T/delP_M on the gridstatevector'

    nullify(delPsfc)
    nullify(field_Psfc)
    nullify(delP_T)
    nullify(delP_M)
    nullify(dP_dPsfc_T)
    nullify(dP_dPsfc_M)

    delP_T => gsv_getField_r8(statevector,'P_T')
    delP_M => gsv_getField_r8(statevector,'P_M')
    delPsfc => gsv_getField_r8(statevector,'P0')
    field_Psfc => gsv_getField_r8(statevector_trial,'P0')

    nlev_T = gsv_getNumLev(statevector,'TH')
    nlev_M = gsv_getNumLev(statevector,'MM')
    numStep = statevector%numstep

    allocate(Psfc(statevector%myLonBeg:statevector%myLonEnd, &
                  statevector%myLatBeg:statevector%myLatEnd))

    do stepIndex = 1, numStep

      Psfc(:,:) = field_Psfc(:,:,1,stepIndex)

      ! dP_dPsfc_M
      nullify(dP_dPsfc_M)
      status = vgd_dpidpis(statevector%vco%vgrid, &
                           statevector%vco%ip1_M, &
                           dP_dPsfc_M, &
                           Psfc)
      if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_dpidpis')
      ! calculate delP_M
      do lev_M = 1, nlev_M
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do lonIndex = statevector%myLonBeg, statevector%myLonEnd
            delPsfc(lonIndex,latIndex,1,stepIndex) = delPsfc(lonIndex,latIndex,1,stepIndex) + dP_dPsfc_M(lonIndex-statevector%myLonBeg+1,latIndex-statevector%myLatBeg+1,lev_M) * delP_M(lonIndex,latIndex,lev_M,stepIndex)
          end do
        end do
      end do
      deallocate(dP_dPsfc_M)

      ! dP_dPsfc_T
      nullify(dP_dPsfc_T)
      status = vgd_dpidpis(statevector%vco%vgrid, &
                           statevector%vco%ip1_T, &
                           dP_dPsfc_T, &
                           Psfc)
      if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_dpidpis')
      ! calculate delP_T
      do lev_T = 1, nlev_T
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do lonIndex = statevector%myLonBeg, statevector%myLonEnd
            delPsfc(lonIndex,latIndex,1,stepIndex) = delPsfc(lonIndex,latIndex,1,stepIndex) + dP_dPsfc_T(lonIndex-statevector%myLonBeg+1,latIndex-statevector%myLatBeg+1,lev_T) * delP_T(lonIndex,latIndex,lev_T,stepIndex)
          end do
        end do
      end do
      deallocate(dP_dPsfc_T)

    end do

    deallocate(Psfc)

  end subroutine calcPressure_ad

end module gridVariableTransforms_mod
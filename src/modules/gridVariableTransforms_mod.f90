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
  use mpi_mod
  use mpivar_mod
  use codePrecision_mod
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
  use calcHeightAndPressure_mod
  use utilities_mod
  use humiditylimits_mod
  
  implicit none
  save
  private

  ! public procedures
  public :: gvt_setup, gvt_transform, gvt_getStateVectorTrial
  public :: gvt_setupRefFromTrialFiles, gvt_setupRefFromStateVector

  logical :: varKindCHTrialsInitialized(vnl_numVarMax)  = .false.

  type(struct_hco), pointer :: hco_anl => null()
  type(struct_vco), pointer :: vco_anl => null()
  type(struct_hco), pointer :: hco_trl => null()
  type(struct_vco), pointer :: vco_trl => null()

  type(struct_gsv), target :: stateVectorRefHU
  type(struct_gsv), target :: stateVectorTrialvarKindCH(vnl_numVarMax)
  type(struct_gsv), target :: stateVectorRefHeight

  ! module interfaces
  interface gvt_transform
    module procedure gvt_transform_gsv
    module procedure gvt_transform_ens
  end interface gvt_transform

CONTAINS

  !--------------------------------------------------------------------------
  ! gvt_setup
  !--------------------------------------------------------------------------
  subroutine gvt_setup(hco_in,hco_core,vco_in)
    ! 
    ! :Purpose: To set up a variable transformation object
    !
    implicit none

    ! Arguments
    type(struct_hco), pointer, intent(inout) :: hco_in
    type(struct_hco), pointer, intent(inout) :: hco_core
    type(struct_vco), pointer, intent(inout) :: vco_in
    
    if ( gsv_containsNonZeroValues(stateVectorRefHU) ) return
    if ( gsv_containsNonZeroValues(stateVectorRefHeight) ) return
    if ( any(varKindCHTrialsInitialized(:)) ) return

    write(*,*) 'gvt_setup: starting'

    hco_anl => hco_in
    vco_anl => vco_in

    call agd_setupFromHco(hco_anl,hco_core)

    write(*,*) 'gvt_setup: done'

  end subroutine gvt_setup

  !--------------------------------------------------------------------------
  ! gvt_setupRefFromTrialFiles
  !--------------------------------------------------------------------------
  subroutine gvt_setupRefFromTrialFiles(varName, varKind_opt)
    ! 
    ! :Purpose: Initialise reference statevector from file 
    !
    ! :Arguments:
    !   :varKind_opt: optional variable "kind" argument presently used to 
    !                 initialise the reference state in a chemical assimilation
    !                 context.
    !
    implicit none

    ! Arguments
    character(len=*), intent(in) :: varName ! reference variable/type used
    character(len=*), optional, intent(in) :: varKind_opt ! additional variable/type information mandatory for some initialization

    ! Locals
    type(struct_gsv) :: statevector_noZnoP
    integer :: varIndex
    
    select case ( trim(varName) )
    case ('HU')
      ! initialize stateVectorRefHU on analysis grid
      call gsv_allocate(stateVectorRefHU, tim_nstepobsinc, hco_anl, vco_anl,   &
                        dateStamp_opt=tim_getDateStamp(), &
                        mpi_local_opt=.true., allocHeightSfc_opt=.true., &
                        hInterpolateDegree_opt='LINEAR', &
                        varNames_opt=(/'HU','P0'/) )

      ! read trial files using default horizontal interpolation degree
      call gsv_readTrials( stateVectorRefHU )  ! IN/OUT

    case ('height')
      if ( .not. stateVectorRefHeight%allocated ) then
        ! initialize stateVectorRefHeight on analysis grid
        call gsv_allocate(stateVectorRefHeight, tim_nstepobsinc, hco_anl, &
                          vco_anl, dateStamp_opt=tim_getDateStamp(), &
                          mpi_local_opt=.true., allocHeightSfc_opt=.true., &
                          hInterpolateDegree_opt='LINEAR', &
                          varNames_opt=&
                              (/'Z_T','Z_M','P_T','P_M','TT','HU','P0'/) )
        call gsv_zero(stateVectorRefHeight)
      end if

      ! initialize statevector_noZnoP on analysis grid
      call gsv_allocate(statevector_noZnoP, tim_nstepobsinc, hco_anl, vco_anl, &
                        dateStamp_opt=tim_getDateStamp(), &
                        mpi_local_opt=.true., allocHeightSfc_opt=.true., &
                        hInterpolateDegree_opt='LINEAR', &
                        varNames_opt=(/'TT','HU','P0'/))
      write(*,*) 'gvt_setupRefFromTrialFiles: statevector_noZnoP allocated'

      ! read trial files using default horizontal interpolation degree
      call gsv_readTrials( statevector_noZnoP )  ! IN/OUT

      ! copy the statevectors
      call gsv_copy(statevector_noZnoP, stateVectorRefHeight, &
                    allowVarMismatch_opt=.true. )

      call gsv_deallocate(statevector_noZnoP)

      ! do height/P calculation of the grid
      call czp_calcZandP_nl( stateVectorRefHeight )

    case default
      if ( present(varKind_opt) ) then
        if (varKind_opt == 'CH' .and. &
            vnl_varKindFromVarname(varName) == varKind_opt ) then
        
          varIndex = vnl_varListIndex(varName)
          
          ! initialize stateVectorTrialvarKindCH(varIndex) on analysis grid
          
          call gsv_allocate(stateVectorTrialvarKindCH(varIndex), &
                            tim_nstepobsinc, hco_anl, vco_anl,   &
                            dateStamp_opt=tim_getDateStamp(), &
                            mpi_local_opt=.true., allocHeightSfc_opt=.true., &
                            hInterpolateDegree_opt='LINEAR', &
                            varNames_opt=(/trim(varName),'P0'/) )

          ! read trial files using default horizontal interpolation degree
          call gsv_readTrials( stateVectorTrialvarKindCH(varIndex) )  ! IN/OUT

          varKindCHTrialsInitialized(varIndex) = .true.
           
        else
          call utl_abort('gvt_setupRefFromTrialFiles: unknown variable ='//trim(varName))
        end if 
      else
        call utl_abort('gvt_setupRefFromTrialFiles: unknown variable ='//trim(varName))
      end if
    end select

  end subroutine gvt_setupRefFromTrialFiles

  !--------------------------------------------------------------------------
  ! gvt_transform_gsv
  !--------------------------------------------------------------------------
  subroutine gvt_transform_gsv( statevector, transform, statevectorOut_opt,  &
                                stateVectorRef_opt, varName_opt, &
                                allowOverWrite_opt )
    ! 
    ! :Purpose: Top-level switch routine for transformations on the grid.
    !
    implicit none
   
    ! Arguments
    type(struct_gsv), intent(inout)           :: statevector ! statevector operand of the transformation
    character(len=*), intent(in)              :: transform ! string identifying the requested transformation
    type(struct_gsv), optional, intent(inout) :: statevectorOut_opt
    type(struct_gsv), optional, intent(in)    :: statevectorRef_opt ! reference statevector necessary for some transformation
    logical, optional, intent(in)             :: allowOverWrite_opt
    character(len=*), optional, intent(in)    :: varName_opt ! additional variable/type information mandatory for some transformation

    ! check stateVector and statevectorOut_opt are on the same grid
    if ( present(stateVectorRef_opt) ) then
      if ( .not. hco_equal(stateVectorRef_opt%hco, stateVector%hco) .or. &
           .not. vco_equal(stateVectorRef_opt%vco, stateVector%vco) ) then
        call utl_abort('gvt_transform: stateVectorRef_opt not on same grid as stateVector')
      end if
    end if

    select case(trim(transform))

    case ('AllTransformedToModel') ! Do all transformed variables: LPRtoPR, LVIStoVIS
      if ( gsv_varExist(statevector,'LPR') ) then
        write(*,*) 'gvt_transform: calling LPRtoPR stateVector transformation'
        call LPRtoPR_gsv(statevector, statevectorOut_opt=statevectorOut_opt,  &
                         allowOverWrite_opt=allowOverWrite_opt)
      end if
      if ( gsv_varExist(statevector,'LVIS') ) then
        write(*,*) 'gvt_transform: calling LVIStoVIS stateVector transformation'
        call LVIStoVIS(statevector, statevectorOut_opt=statevectorOut_opt,  &
                       allowOverWrite_opt=allowOverWrite_opt)
      end if

    case ('UVtoVortDiv')
      if (present(statevectorOut_opt)) then
        call utl_abort('gvt_transform: for UVtoVortDiv, the option statevectorOut_opt is not yet available')
      end if
      call UVtoVortDiv_gsv(statevector)

    case ('VortDivToPsiChi')
      if (.not. gsv_varExist(statevector,'QR') .or. &
          .not. gsv_varExist(statevector,'DD') ) then
        call utl_abort('gvt_transform: for VortDivToPsiChi, variables QR and DD must be allocated in gridstatevector')
      end if
      if (present(statevectorOut_opt)) then
        call utl_abort('gvt_transform: for VortDivToPsiChi, the option statevectorOut_opt is not yet available')
      end if
      call VortDivToPsiChi_gsv(statevector)

    case ('UVtoPsiChi')
      if (.not. gsv_varExist(statevector,'PP') .or. &
          .not. gsv_varExist(statevector,'CC') ) then
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
        call utl_abort('gvt_transform: for HUtoLQ_tlm, the option statevectorOut_opt is not yet available')
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

    case ('LPRoPR')
      if ( .not. gsv_varExist(statevector,'PR') ) then
        call utl_abort('gvt_transform: for LPRtoPR, variable PR must be allocated in gridstatevector')
      end if
      call LPRtoPR_gsv(statevector, statevectorOut_opt,  &
                       allowOverWrite_opt=allowOverWrite_opt)

    case ('PRtoLPR')
      if ( .not. gsv_varExist(statevector,'PR') ) then
        call utl_abort('gvt_transform: for PRtoLPR, variable PR must be allocated in gridstatevector')
      end if
      call PRtoLPR_gsv(statevector, statevectorOut_opt)

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
        if (.not. gsv_varExist(statevector,'VIS') .or. &
            .not. gsv_varExist(statevector,'LVIS') ) then
          call utl_abort('gvt_transform: for LVIStoVIS, variables LVIS and VIS must be allocated in gridstatevector')
        end if
        call LVIStoVIS(statevector)
      end if

    case ('ZandP_nl')
      if ( present(statevectorOut_opt) ) then
        call utl_abort('gvt_transform: for ZandP_nl, the option statevectorOut_opt is not yet available')
      end if
      call czp_calcZandP_nl(stateVector)

    case ('ZandP_tl')
      if ( present(statevectorOut_opt) ) then
        call utl_abort('gvt_transform: for ZandP_tl, the option statevectorOut_opt is not yet available')
      end if
      call ZandP_tl(stateVector)

    case ('ZandP_ad')
      if ( present(statevectorOut_opt) ) then
        call utl_abort('gvt_transform: for ZandP_ad, the option statevectorOut_opt is not yet available')
      end if
      call ZandP_ad(stateVector)

    case ('expCH_tlm')
      if ( .not. gsv_varKindExist('CH')  ) then
        call utl_abort('gvt_transform: for expCH_tlm, variables of CH kind must be allocated in gridstatevector')
      else if ( .not.present(varName_opt) ) then
        call utl_abort('gvt_transform: for expCH_tlm, missing variable name')
      else  if ( .not. gsv_varExist(statevector,trim(varName_opt)) ) then
        call utl_abort('gvt_transform: for expCH_tlm, variable '//trim(varName_opt)//' must be allocated in gridstatevector')
      else if ( vnl_varKindFromVarname(trim(varName_opt)) /= 'CH' ) then
        call utl_abort('gvt_transform: Invalid kind of varName for expCH_tlm')
      else if ( .not. present(statevectorRef_opt)) then
        call utl_abort('gvt_transform: for expCH_tlm, the option statevectorRef_opt must be present')
      else if (present(statevectorOut_opt)) then
        call utl_abort('gvt_transform: for expCH_tlm, the option statevectorOut_opt is not yet available')
      end if
      call expCH_tlm(statevector, varName_opt, stateVectorRef_opt)

   case ('expCH_ad')
      if ( .not. gsv_varKindExist('CH')  ) then
        call utl_abort('gvt_transform: for expCH_ad, variables of CH kind must be allocated in gridstatevector')
      else if ( .not.present(varName_opt) ) then
        call utl_abort('gvt_transform: for expCH_ad missing variable name')
      else  if ( .not. gsv_varExist(statevector,trim(varName_opt)) ) then
        call utl_abort('gvt_transform: for expCH_ad, variable '//trim(varName_opt)//' must be allocated in gridstatevector')
      else if ( vnl_varKindFromVarname(trim(varName_opt)) /= 'CH' ) then
        call utl_abort('gvt_transform: Invalid kind of varName for expCH_ad')
      else if ( .not. present(statevectorRef_opt)) then
        call utl_abort('gvt_transform: for expCH_ad, the option statevectorRef_opt must be present')
      else if (present(statevectorOut_opt)) then
        call utl_abort('gvt_transform: for expCH_ad, the option statevectorOut_opt is not yet available')
      end if
      call expCH_tlm(statevector, varName_opt, stateVectorRef_opt) ! self-adjoint

    case ('CH_bounds')
      if ( .not. gsv_varKindExist('CH')  ) then
        call utl_abort('gvt_transform: for CH_bounds, variables of CH kind must be allocated in gridstatevector')
      end if
      call CH_bounds(statevector)

    case ('GLtoLG')
      if ( .not. present(statevectorRef_opt)) then
        call utl_abort('gvt_transform: for GLtoLG, the option statevectorRef_opt must be present')
      end if
      call GLtoLG(statevector, stateVectorRef_opt)

    case default
      write(*,*)
      write(*,*) 'Unsupported function : ', trim(transform)
      call utl_abort('gvt_transform')
    end select

  end subroutine gvt_transform_gsv

  !--------------------------------------------------------------------------
  ! gvt_transform_ens
  !--------------------------------------------------------------------------
    subroutine gvt_transform_ens( ens,transform, allowOverWrite_opt, &
                                  varName_opt, huMinValue_opt)
    ! 
    ! :Purpose: Top-level switch routine for ensemble transformations on the 
    !           grid
    !
    implicit none
   
    ! Arguments
    type(struct_ens), intent(inout) :: ens ! operand (ensemble of statevector)
    character(len=*), intent(in)    :: transform ! string identifying the requested transformation
    logical, optional, intent(in)   :: allowOverWrite_opt
    character(len=*), optional, intent(in)  :: varName_opt ! additional variable/type information mandatory for some transformation
    real(8), optional, intent(in)   :: huMinValue_opt ! HU min value for 'HUtoLQ' transformation

    select case(trim(transform))
    case ('AllTransformedToModel') ! Do all transformed variables: LPRtoPR
      if ( ens_varExist(ens,'LPR') ) then
        write(*,*) 'gvt_transform: calling LPRtoPR ensemble transformation'
        call LPRtoPR_ens(ens, allowOverWrite_opt=allowOverWrite_opt)
      end if
    case ('HUtoLQ')
      call HUtoLQ_ens(ens, huMinValue_opt=huMinValue_opt)
    case ('LPRtoPR')
      call LPRtoPR_ens(ens, allowOverWrite_opt=allowOverWrite_opt)
    case ('UVtoPsiChi')
      call UVtoPsiChi_ens(ens)
    case ('UVtoVortDiv')
      call UVtoVortDiv_ens(ens)
    case ('logCH')
      if ( .not.present(varName_opt) ) then
        call utl_abort('gvt_transform: for logCH missing variable name')
      else if ( vnl_varKindFromVarname(trim(varName_opt)) /= 'CH' ) then
        call utl_abort('gvt_transform: Invalid kind of varName for logCH')
      end if 
      call logCH_ens(ens,varName_opt)
    case default
      call utl_abort('gvt_transform_ens: Unsupported function '//trim(transform))
    end select

  end subroutine gvt_transform_ens

  !--------------------------------------------------------------------------
  ! gvt_getStateVectorTrial
  !--------------------------------------------------------------------------
  function gvt_getStateVectorTrial(varName) result(statevector_ptr)
    ! 
    ! :Purpose: Returns a pointer to requested reference statevector 
    !
    implicit none

    ! Arguments
    character(len=*), intent(in) :: varName ! reference variable/type requested 

    ! Locals
    type(struct_gsv), pointer  :: statevector_ptr

    select case ( trim(varName) )
    case ('HU')
      if ( .not. gsv_containsNonZeroValues(stateVectorRefHU) ) then
        call utl_abort('gvt_getStateVectorTrial: do trials to stateVectorRefHU transform at higher level')
      end if
      statevector_ptr => stateVectorRefHU

    case ('height')
      if ( .not. gsv_containsNonZeroValues(stateVectorRefHeight) ) then
        call utl_abort('gvt_getStateVectorTrial: do trials to stateVectorRefHeight transform at higher level')
      end if
      statevector_ptr => stateVectorRefHeight

    case default
      call utl_abort('gvt_getStateVectorTrial: unknown variable ='//trim(varName))
    end select

  end function gvt_getStateVectorTrial

  !--------------------------------------------------------------------------
  ! gvt_setupRefFromStateVector
  !--------------------------------------------------------------------------
  subroutine gvt_setupRefFromStateVector( stateVectorOnTrlGrid, varName, &
                                          applyLimitOnHU_opt )
    !
    !:Purpose: computing the reference stateVector on the analysis grid at each 
    !          outer-loop iteration. The calculation is skipped if stateVectorRef* is 
    !          initialized (gsv_containsNonZeroValue(stateVectorRef*)=.true.).
    !          The input stateVector is the high spatial/temporal resolution
    !          statevector used for reading the trials and should contain 
    !          TT/HU/P0 if stateVectorRefHeight is asked for.
    !
    implicit none

    ! Arguments
    type(struct_gsv),  intent(in) :: stateVectorOnTrlGrid ! high spatial/temporal resolution statevector 
    character(len=*),  intent(in) :: varName ! reference variable/type used
    logical, intent(in), optional :: applyLimitOnHU_opt

    ! Locals
    type(struct_gsv)         :: stateVectorLowResTime
    type(struct_gsv)         :: stateVectorLowResTimeSpace
    type(struct_gsv), target :: stateVectorRefHUTT

    logical :: allocHeightSfc
    character(len=4), pointer :: varNames(:)

    if ( mpi_myid == 0 ) write(*,*) 'gvt_setupRefFromStateVector: START'

    if ( .not. associated(hco_trl) ) hco_trl => gsv_getHco(stateVectorOnTrlGrid)
    if ( .not. associated(vco_trl) ) vco_trl => gsv_getVco(stateVectorOnTrlGrid)

    select case ( trim(varName) )
    case ('HU')
      if ( .not. present(applyLimitOnHU_opt) ) then
        call utl_abort('gvt_setupRefFromStateVector: applyLimitOnHU_opt for RefHU missing')
      end if

      if ( .not. stateVectorRefHU%allocated ) then
        call gsv_allocate(stateVectorRefHU, tim_nstepobsinc, hco_anl, vco_anl,&
                          dateStamp_opt=tim_getDateStamp(), &
                          mpi_local_opt=.true., allocHeightSfc_opt=.true., &
                          hInterpolateDegree_opt='LINEAR', &
                          varNames_opt=(/'HU','P0'/) )
      else
        call gsv_zero( stateVectorRefHU )
      end if

      allocHeightSfc = ( vco_trl%Vcode /= 0 )

      call gsv_allocate(stateVectorRefHUTT, tim_nstepobsinc, hco_anl, vco_anl,&
                        dateStamp_opt=tim_getDateStamp(),&
                        mpi_local_opt=.true., allocHeightSfc_opt=.true., &
                        hInterpolateDegree_opt='LINEAR', &
                        varNames_opt=(/'HU','TT','P0'/) )

      ! First, degrade the time steps
      call gsv_allocate(stateVectorLowResTime, tim_nstepobsinc, hco_trl, &
                        vco_trl, dataKind_opt=pre_incrReal, &
                        dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true.,&
                        allocHeightSfc_opt=allocHeightSfc, &
                        hInterpolateDegree_opt='LINEAR', & 
                        allocHeight_opt=.false., allocPressure_opt=.false. )
      call gsv_copy( stateVectorOnTrlGrid, stateVectorLowResTime, &
                     allowTimeMismatch_opt=.true., allowVarMismatch_opt=.true. )

      ! Second, interpolate to the low-resolution spatial grid.
      nullify(varNames)
      call gsv_varNamesList(varNames, stateVectorLowResTime)
      call gsv_allocate(stateVectorLowResTimeSpace, tim_nstepobsinc, hco_anl, &
                        vco_anl, dataKind_opt=pre_incrReal, &
                        dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true.,&
                        allocHeightSfc_opt=.true., &
                        hInterpolateDegree_opt='LINEAR', varNames_opt=varNames)
      call gsv_interpolate(stateVectorLowResTime, stateVectorLowResTimeSpace)

      ! Now copy only P0, HU, and TT to create reference stateVector.
      call gsv_copy(stateVectorLowResTimeSpace, stateVectorRefHUTT, &
                    allowTimeMismatch_opt=.false., allowVarMismatch_opt=.true.)
      call gsv_deallocate(stateVectorLowResTimeSpace)
      call gsv_deallocate(stateVectorLowResTime)

      ! Impose humidity limits on stateVectorRefHUTT
      if ( applyLimitOnHU_opt ) then
        write(*,*) 'var: impose limits on stateVectorRefHUTT'
        call qlim_saturationLimit(stateVectorRefHUTT)
        call qlim_rttovLimit(stateVectorRefHUTT)
      end if

      call gsv_copy(stateVectorRefHUTT, stateVectorRefHU, &
                    allowTimeMismatch_opt=.false., allowVarMismatch_opt=.true.)
      call gsv_deallocate(stateVectorRefHUTT)

    case ('height')
      if ( .not. stateVectorRefHeight%allocated ) then
        call gsv_allocate( stateVectorRefHeight, tim_nstepobsinc, hco_anl, &
                           vco_anl, dateStamp_opt=tim_getDateStamp(), &
                           mpi_local_opt=.true., allocHeightSfc_opt=.true., &
                           hInterpolateDegree_opt='LINEAR', &
                           varNames_opt=&
                              (/'Z_T','Z_M','P_T','P_M','TT','HU','P0'/) )
      else
        call gsv_zero( stateVectorRefHeight )
      end if

      ! First, degrade the time steps
      call gsv_allocate(stateVectorLowResTime, tim_nstepobsinc, hco_trl, &
                        vco_trl, dateStamp_opt=tim_getDateStamp(), &
                        mpi_local_opt=.true., allocHeightSfc_opt=.true., &
                        hInterpolateDegree_opt='LINEAR', &
                        varNames_opt=(/'TT','HU','P0'/) )
      call gsv_copy(stateVectorOnTrlGrid, stateVectorLowResTime, &
                    allowTimeMismatch_opt=.true., allowVarMismatch_opt=.true.)

      ! Second, interpolate to the low-resolution spatial grid.
      nullify(varNames)
      call gsv_varNamesList(varNames, stateVectorLowResTime)
      call gsv_allocate(stateVectorLowResTimeSpace, tim_nstepobsinc, hco_anl,&
                        vco_anl, dateStamp_opt=tim_getDateStamp(), &
                        mpi_local_opt=.true., allocHeightSfc_opt=.true., &
                        hInterpolateDegree_opt='LINEAR', varNames_opt=varNames)
      call gsv_interpolate(stateVectorLowResTime, stateVectorLowResTimeSpace)

      ! Now copy to create final stateVector height.
      call gsv_copy(stateVectorLowResTimeSpace, stateVectorRefHeight, &
                    allowTimeMismatch_opt=.false., allowVarMismatch_opt=.true.)
      call gsv_deallocate(stateVectorLowResTimeSpace)
      call gsv_deallocate(stateVectorLowResTime)

      ! do height/P calculation of the grid
      call czp_calcZandP_nl(stateVectorRefHeight)

    end select

    if ( mpi_myid == 0 ) write(*,*) 'gvt_setupRefFromStateVector: END'

  end subroutine gvt_setupRefFromStateVector

  !--------------------------------------------------------------------------
  ! LQtoHU
  !--------------------------------------------------------------------------
  subroutine LQtoHU(statevector)
    ! 
    ! :Purpose: Specific humidity logarithm exponentiation. 
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout) :: statevector

    ! Locals
    integer :: lonIndex,latIndex,levIndex,stepIndex

    real(8), pointer :: hu_ptr(:,:,:,:), lq_ptr(:,:,:,:)

    call gsv_getField(statevector,hu_ptr,'HU')
    call gsv_getField(statevector,lq_ptr,'HU')

    do stepIndex = 1, statevector%numStep
      !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
      do levIndex = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname('HU'))
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do lonIndex = statevector%myLonBeg, statevector%myLonEnd
            hu_ptr(lonIndex,latIndex,levIndex,stepIndex) = &
                exp(lq_ptr(lonIndex,latIndex,levIndex,stepIndex))
          end do
        end do
      end do
      !$OMP END PARALLEL DO
    end do

  end subroutine LQtoHU

  !--------------------------------------------------------------------------
  ! HUtoLQ_gsv
  !--------------------------------------------------------------------------
  subroutine HUtoLQ_gsv(statevector)
    ! 
    ! :Purpose: Logarithmic transformation of specific humidity 
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout) :: statevector

    ! Locals
    integer :: lonIndex,latIndex,levIndex,stepIndex

    real(4), pointer :: hu_ptr_r4(:,:,:,:), lq_ptr_r4(:,:,:,:)
    real(8), pointer :: hu_ptr_r8(:,:,:,:), lq_ptr_r8(:,:,:,:)

    if ( gsv_getDataKind(statevector) == 8 ) then

      call gsv_getField(statevector,hu_ptr_r8,'HU')
      call gsv_getField(statevector,lq_ptr_r8,'HU')

      do stepIndex = 1, statevector%numStep
        !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
        do levIndex = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname('HU'))
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd
              lq_ptr_r8(lonIndex,latIndex,levIndex,stepIndex) = &
                  log(max(hu_ptr_r8(lonIndex,latIndex,levIndex,stepIndex),&
                          gsv_rhumin) )
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do

    else

      call gsv_getField(statevector,hu_ptr_r4,'HU')
      call gsv_getField(statevector,lq_ptr_r4,'HU')

      do stepIndex = 1, statevector%numStep
        !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
        do levIndex = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname('HU'))
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd
              lq_ptr_r4(lonIndex,latIndex,levIndex,stepIndex) = &
                  log(max(hu_ptr_r4(lonIndex,latIndex,levIndex,stepIndex),&
                          real(gsv_rhumin,4)) )
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do

    end if

  end subroutine HUtoLQ_gsv

  !--------------------------------------------------------------------------
  ! HUtoLQ_ens
  !--------------------------------------------------------------------------
  subroutine HUtoLQ_ens(ens, huMinValue_opt)
    ! 
    ! :Purpose: Specific humidity logarithm exponentiation (ensemble processing) 
    !
    implicit none

    ! Arguments
    type(struct_ens), intent(inout) :: ens
    real(8), optional, intent(in)   :: huMinValue_opt

    ! Locals
    integer :: lonIndex, latIndex, levIndex, stepIndex, memberIndex
    integer :: myLatBeg, myLatEnd, myLonBeg, myLonEnd
    character(len=4) :: varName
    real(4), pointer :: ptr4d_r4(:,:,:,:)
    real(4)          :: huMinValue

    if (present(huMinValue_opt)) then
      huMinValue = huMinValue_opt
    else
      huMinValue = MPC_MINIMUM_HU_R4
    end if

    call ens_getLatLonBounds(ens, myLonBeg, myLonEnd, myLatBeg, myLatEnd)

    do levIndex = 1, ens_getNumK(ens)

      varName = ens_getVarNameFromK(ens,levIndex)
      if (varName /= 'HU') cycle

      ptr4d_r4 => ens_getOneLev_r4(ens,levIndex)
    
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          do stepIndex = 1, ens_getNumStep(ens)
            do memberIndex = 1, ens_getNumMembers(ens)
              ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex) = &
                   log(max( ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex),&
                            huMinValue) )
            end do
          end do
        end do
      end do

    end do

  end subroutine HUtoLQ_ens

  !--------------------------------------------------------------------------
  ! LQtoHU_tlm
  !--------------------------------------------------------------------------
  subroutine LQtoHU_tlm(statevector, stateVectorRef_opt)
    ! 
    ! :Purpose: Tangent linear of exponentiation transformation of specific 
    !           humidity in logarithmic representation
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout)         :: statevector
    type(struct_gsv), optional, intent(in)  :: statevectorRef_opt

    ! Locals
    integer :: lonIndex,latIndex,levIndex,stepIndex
    real(8), pointer :: hu_ptr_r8(:,:,:,:), lq_ptr_r8(:,:,:,:)
    real(8), pointer :: hu_trial(:,:,:,:)
    real(4), pointer :: hu_ptr_r4(:,:,:,:), lq_ptr_r4(:,:,:,:)

    if ( present(statevectorRef_opt) ) then
      call gsv_getField(stateVectorRef_opt,hu_trial,'HU')
    else
      if ( .not. gsv_containsNonZeroValues(stateVectorRefHU) ) then
        call utl_abort('LQtoHU_tlm: do trials to stateVectorRefHU transform at higher level')
      end if
      call gsv_getField(stateVectorRefHU,hu_trial,'HU')
    end if

    if ( gsv_getDataKind(statevector) == 4 ) then
      call gsv_getField(statevector,hu_ptr_r4,'HU')
      call gsv_getField(statevector,lq_ptr_r4,'HU')

      do stepIndex = 1, statevector%numStep
        !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
        do levIndex = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname('HU'))
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd       
              hu_ptr_r4(lonIndex,latIndex,levIndex,stepIndex) =  &
                   lq_ptr_r4(lonIndex,latIndex,levIndex,stepIndex)*  &
                   max( hu_trial(lonIndex,latIndex,levIndex,stepIndex),&
                        MPC_MINIMUM_HU_R4)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do
    else
      call gsv_getField(statevector,hu_ptr_r8,'HU')
      call gsv_getField(statevector,lq_ptr_r8,'HU')

      do stepIndex = 1, statevector%numStep
        !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
        do levIndex = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname('HU'))
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd       
              hu_ptr_r8(lonIndex,latIndex,levIndex,stepIndex) =  &
                   lq_ptr_r8(lonIndex,latIndex,levIndex,stepIndex)*  &
                   max( hu_trial(lonIndex,latIndex,levIndex,stepIndex),&
                        MPC_MINIMUM_HU_R8)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do
    end if

  end subroutine LQtoHU_tlm

  !--------------------------------------------------------------------------
  ! HUtoLQ_tlm
  !--------------------------------------------------------------------------
  subroutine HUtoLQ_tlm(statevector, stateVectorRef_opt)
    ! 
    ! :Purpose: Tangent linear of logarithmic transformation of specific 
    !           humidity
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout)         :: statevector
    type(struct_gsv), optional, intent(in)  :: statevectorRef_opt

    ! Locals
    integer           :: lonIndex,latIndex,levIndex,stepIndex
    real(8), pointer  :: hu_ptr_r8(:,:,:,:), lq_ptr_r8(:,:,:,:)
    real(8), pointer  :: hu_trial(:,:,:,:)
    real(4), pointer  :: hu_ptr_r4(:,:,:,:), lq_ptr_r4(:,:,:,:)

    if ( present(statevectorRef_opt) ) then
      call gsv_getField(stateVectorRef_opt,hu_trial,'HU')
    else
      if ( .not. gsv_containsNonZeroValues(stateVectorRefHU) ) then
        call utl_abort('HUtoLQ_tlm: do trials to stateVectorRefHU transform at higher level')
      end if
      call gsv_getField(stateVectorRefHU,hu_trial,'HU')
    end if

    if ( gsv_getDataKind(statevector) == 4 ) then
      call gsv_getField(statevector,hu_ptr_r4,'HU')
      call gsv_getField(statevector,lq_ptr_r4,'HU')

      do stepIndex = 1, statevector%numStep
        !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
        do levIndex = 1,gsv_getNumLev(statevector,vnl_varLevelFromVarname('HU'))
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd
              lq_ptr_r4(lonIndex,latIndex,levIndex,stepIndex) =  &
                   hu_ptr_r4(lonIndex,latIndex,levIndex,stepIndex) / &
                   max( hu_trial(lonIndex,latIndex,levIndex,stepIndex),&
                        MPC_MINIMUM_HU_R8)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do
    else
      call gsv_getField(statevector,hu_ptr_r8,'HU')
      call gsv_getField(statevector,lq_ptr_r8,'HU')

      do stepIndex = 1, statevector%numStep
        !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
        do levIndex = 1,gsv_getNumLev(statevector,vnl_varLevelFromVarname('HU'))
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd
              lq_ptr_r8(lonIndex,latIndex,levIndex,stepIndex) =  &
                   hu_ptr_r8(lonIndex,latIndex,levIndex,stepIndex) / &
                   max( hu_trial(lonIndex,latIndex,levIndex,stepIndex),&
                        MPC_MINIMUM_HU_R8)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do
    end if

  end subroutine HUtoLQ_tlm

  !--------------------------------------------------------------------------
  ! LPRtoPR_gsv
  !--------------------------------------------------------------------------
  subroutine LPRtoPR_gsv(stateVector, statevectorOut_opt, allowOverWrite_opt)
    ! 
    ! :Purpose: Quantity of precipitation logarithm exponentiation 
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout)           :: stateVector
    type(struct_gsv), optional, intent(inout) :: statevectorOut_opt
    logical, optional, intent(in)             :: allowOverWrite_opt

    ! Locals
    integer :: lonIndex,latIndex,levIndex,stepIndex
    logical :: overWriteNeeded
    real(4), pointer :: pr_ptr_r4(:,:,:,:), lpr_ptr_r4(:,:,:,:)
    real(8), pointer :: pr_ptr_r8(:,:,:,:), lpr_ptr_r8(:,:,:,:)

    ! Check if overWrite of PR is needed, but not allowed
    overWriteNeeded = .false.
    if (present(statevectorOut_opt)) then
      if (.not. gsv_varExist(stateVectorOut_opt,'PR')) then
        overWriteNeeded = .true.
      end if
    else
      if (.not. gsv_varExist(stateVector,'PR')) then
        overWriteNeeded = .true.
      end if
    end if
    if (overWriteNeeded) then
      if(present(allowOverWrite_opt)) then
        if (.not.allowOverWrite_opt) then
          call utl_abort('LPRtoPR_gsv: allowOverWrite_opt is false, but PR not present in stateVector')
        end if
      else
        call utl_abort('LPRtoPR_gsv: allowOverWrite_opt not specified, but PR not present in stateVector')
      end if
    end if

    if ( gsv_getDataKind(stateVector) == 8 ) then

      if (present(statevectorOut_opt)) then
        if (gsv_varExist(stateVectorOut_opt,'PR')) then
          call gsv_getField(statevectorOut_opt,pr_ptr_r8,'PR')
        else
          call gsv_getField(statevectorOut_opt,pr_ptr_r8,'LPR')
        end if
      else
        if (gsv_varExist(stateVector,'PR')) then
          call gsv_getField(stateVector,pr_ptr_r8,'PR')
        else
          call gsv_getField(stateVector,pr_ptr_r8,'LPR')
        end if
      end if
      call gsv_getField(stateVector,lpr_ptr_r8,'LPR')
      
      do stepIndex = 1, stateVector%numStep
        !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
        do levIndex = 1, gsv_getNumLev(stateVector,vnl_varLevelFromVarname('LPR'))
          do latIndex = stateVector%myLatBeg, stateVector%myLatEnd
            do lonIndex = stateVector%myLonBeg, stateVector%myLonEnd
              if (lpr_ptr_r8(lonIndex,latIndex,levIndex,stepIndex) > &
                  0.99D0*MPC_missingValue_R8) then

                pr_ptr_r8(lonIndex,latIndex,levIndex,stepIndex) =   &
                     max( 0.0D0, exp( lpr_ptr_r8(lonIndex,latIndex,levIndex,&
                                      stepIndex)) - &
                          MPC_MINIMUM_PR_R8)
              else
                pr_ptr_r8(lonIndex,latIndex,levIndex,stepIndex) = &
                    MPC_missingValue_R8
              end if
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do

    else

      if (present(statevectorOut_opt)) then
        if (gsv_varExist(stateVectorOut_opt,'PR')) then
          call gsv_getField(statevectorOut_opt,pr_ptr_r4,'PR')
        else
          call gsv_getField(statevectorOut_opt,pr_ptr_r4,'LPR')
        end if
      else
        if (gsv_varExist(stateVector,'PR')) then
          call gsv_getField(stateVector,pr_ptr_r4,'PR')
        else
          call gsv_getField(stateVector,pr_ptr_r4,'LPR')
        end if
      end if
      call gsv_getField(stateVector,lpr_ptr_r4,'LPR')

      do stepIndex = 1, stateVector%numStep
        !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
        do levIndex = 1,gsv_getNumLev(stateVector,vnl_varLevelFromVarname('LPR'))
          do latIndex = stateVector%myLatBeg, stateVector%myLatEnd
            do lonIndex = stateVector%myLonBeg, stateVector%myLonEnd
              if (lpr_ptr_r4(lonIndex,latIndex,levIndex,stepIndex) > &
                  0.99*MPC_missingValue_R4) then

                pr_ptr_r4(lonIndex,latIndex,levIndex,stepIndex) =   &
                     max( 0.0, exp( lpr_ptr_r4(lonIndex,latIndex,levIndex,&
                                    stepIndex)) - &
                          MPC_MINIMUM_PR_R4)
              else
                pr_ptr_r4(lonIndex,latIndex,levIndex,stepIndex) = &
                    MPC_missingValue_R4
              end if
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do

    end if

    ! Change the variable name from LPR to PR if PR does not exist
    if (present(statevectorOut_opt)) then
      if (.not. gsv_varExist(stateVectorOut_opt,'PR')) then
        call gsv_modifyVarName(stateVectorOut_opt,'LPR','PR')
      end if
    else
      if (.not. gsv_varExist(stateVector,'PR')) then
        call gsv_modifyVarName(stateVector,'LPR','PR')
      end if
    end if

  end subroutine LPRtoPR_gsv

  !--------------------------------------------------------------------------
  ! LPRtoPR_ens
  !--------------------------------------------------------------------------
  subroutine LPRtoPR_ens(ens, allowOverWrite_opt)
    ! 
    ! :Purpose: Quantity of precipitation logarithm exponentiation 
    !           (ensemble processing)
    !
    implicit none

    ! Arguments
    type(struct_ens), intent(inout) :: ens
    logical, optional, intent(in)   :: allowOverWrite_opt

    ! Locals
    integer :: lonIndex, latIndex, levIndex, stepIndex, memberIndex
    integer :: myLatBeg, myLatEnd, myLonBeg, myLonEnd, kIndexLPR, kIndexPR
    logical :: overWriteNeeded
    real(4), pointer :: PR_ptr_r4(:,:,:,:), LPR_ptr_r4(:,:,:,:)

    ! Check if overWrite of PR is needed, but not allowed
    overWriteNeeded = .false.
    if (.not. ens_varExist(ens,'PR')) then
      overWriteNeeded = .true.
    end if
    if (overWriteNeeded) then
      if(present(allowOverWrite_opt)) then
        if (.not.allowOverWrite_opt) then
          call utl_abort('LPRtoPR_ens: allowOverWrite_opt is false, but PR not present in ensemble')
        end if
      else
        call utl_abort('LPRtoPR_ens: allowOverWrite_opt not specified, but PR not present in ensemble')
      end if
    end if

    call ens_getLatLonBounds(ens, myLonBeg, myLonEnd, myLatBeg, myLatEnd)

    levIndex = 1
    kIndexLPR = ens_getKFromLevVarName(ens, levIndex, 'LPR')
    if (ens_varExist(ens,'PR')) then
      kIndexPR = ens_getKFromLevVarName(ens, levIndex, 'PR')
    else
      kIndexPR = ens_getKFromLevVarName(ens, levIndex, 'LPR')
    end if

    LPR_ptr_r4 => ens_getOneLev_r4(ens,kIndexLPR)
    PR_ptr_r4  => ens_getOneLev_r4(ens,kIndexPR)

    do latIndex = myLatBeg, myLatEnd
      do lonIndex = myLonBeg, myLonEnd
        do stepIndex = 1, ens_getNumStep(ens)
          do memberIndex = 1, ens_getNumMembers(ens)
              PR_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) =   &
                   max( 0.0, exp( LPR_ptr_r4(memberIndex,stepIndex,lonIndex,&
                                  latIndex)) - &
                        MPC_MINIMUM_PR_R4)
              ! maybe also impose minimum value of LPR
              !LPR_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) =   &
              !     max(log(MPC_MINIMUM_PR_R4), LPR_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex))
          end do
        end do
      end do
    end do

    ! Change the variable name from LPR to PR if PR does not exist
    if (.not. ens_varExist(ens,'PR')) then
      call ens_modifyVarName(ens,'LPR','PR')
    end if

  end subroutine LPRtoPR_ens

  !--------------------------------------------------------------------------
  ! PRtoLPR_gsv
  !--------------------------------------------------------------------------
  subroutine PRtoLPR_gsv(stateVector, statevectorOut_opt)
    ! 
    ! :Purpose: Logarithmic transformation of quantity of precipitation
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout)           :: stateVector
    type(struct_gsv), optional, intent(inout) :: statevectorOut_opt

    ! Locals
    integer :: lonIndex,latIndex,levIndex,stepIndex
    real(4), pointer :: pr_ptr_r4(:,:,:,:), lpr_ptr_r4(:,:,:,:)
    real(8), pointer :: pr_ptr_r8(:,:,:,:), lpr_ptr_r8(:,:,:,:)

    if ( gsv_getDataKind(stateVector) == 8 ) then

      if (present(statevectorOut_opt)) then
        call gsv_getField(statevectorOut_opt,lpr_ptr_r8,'LPR')
      else
        call gsv_getField(stateVector,lpr_ptr_r8,'LPR')
      end if
      call gsv_getField(stateVector,pr_ptr_r8,'PR')
      
      do stepIndex = 1, stateVector%numStep
        !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
        do levIndex = 1, gsv_getNumLev(stateVector,vnl_varLevelFromVarname('PR'))
          do latIndex = stateVector%myLatBeg, stateVector%myLatEnd
            do lonIndex = stateVector%myLonBeg, stateVector%myLonEnd
              if (pr_ptr_r8(lonIndex,latIndex,levIndex,stepIndex) > -1.0D0) then
                lpr_ptr_r8(lonIndex,latIndex,levIndex,stepIndex) =  &
                     log( MPC_MINIMUM_PR_R8 + &
                          max(0.0d0,pr_ptr_r8(lonIndex,latIndex,levIndex,&
                                              stepIndex)))
              else
                lpr_ptr_r8(lonIndex,latIndex,levIndex,stepIndex) = &
                    MPC_missingValue_R8
              end if
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do

    else

      if (present(statevectorOut_opt)) then
        call gsv_getField(statevectorOut_opt,lpr_ptr_r4,'LPR')
      else
        call gsv_getField(stateVector,lpr_ptr_r4,'LPR')
      end if
      call gsv_getField(stateVector,pr_ptr_r4,'PR')

      do stepIndex = 1, stateVector%numStep
        !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
        do levIndex = 1,gsv_getNumLev(stateVector,vnl_varLevelFromVarname('PR'))
          do latIndex = stateVector%myLatBeg, stateVector%myLatEnd
            do lonIndex = stateVector%myLonBeg, stateVector%myLonEnd
              if (pr_ptr_r4(lonIndex,latIndex,levIndex,stepIndex) > -1.0) then
                lpr_ptr_r4(lonIndex,latIndex,levIndex,stepIndex) =  &
                     log( MPC_MINIMUM_PR_R4 + &
                          max(0.0,pr_ptr_r4(lonIndex,latIndex,levIndex,&
                                            stepIndex)))
              else
                lpr_ptr_r4(lonIndex,latIndex,levIndex,stepIndex) = &
                    MPC_missingValue_R4
              end if
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do

    end if

  end subroutine PRtoLPR_gsv

  !--------------------------------------------------------------------------
  ! LVIStoVIS
  !--------------------------------------------------------------------------
  subroutine LVIStoVIS(stateVector, statevectorOut_opt, allowOverWrite_opt)
    ! 
    ! :Purpose: Visibility logarithm exponentiation 
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout)           :: stateVector
    type(struct_gsv), optional, intent(inout) :: statevectorOut_opt
    logical, optional, intent(in)             :: allowOverWrite_opt

    ! Locals
    integer :: lonIndex,latIndex,levIndex,stepIndex
    logical :: overWriteNeeded
    real(4), pointer :: vis_ptr_r4(:,:,:,:), lvis_ptr_r4(:,:,:,:)
    real(8), pointer :: vis_ptr_r8(:,:,:,:), lvis_ptr_r8(:,:,:,:)

    ! Check if overWrite of PR is needed, but not allowed
    overWriteNeeded = .false.
    if (present(statevectorOut_opt)) then
      if (.not. gsv_varExist(stateVectorOut_opt,'PR')) then
        overWriteNeeded = .true.
      end if
    else
      if (.not. gsv_varExist(stateVector,'PR')) then
        overWriteNeeded = .true.
      end if
    end if
    if (overWriteNeeded) then
      if(present(allowOverWrite_opt)) then
        if (.not.allowOverWrite_opt) then
          call utl_abort('LVIStoVIS_gsv: allowOverWrite_opt is false, but PR not present in stateVector')
        end if
      else
        call utl_abort('LVIStoVIS_gsv: allowOverWrite_opt not specified, but PR not present in stateVector')
      end if
    end if

    if ( gsv_getDataKind(stateVector) == 8 ) then

      if (present(statevectorOut_opt)) then
        if (gsv_varExist(stateVectorOut_opt,'VIS')) then
          call gsv_getField(statevectorOut_opt,vis_ptr_r8,'VIS')
        else
          call gsv_getField(statevectorOut_opt,vis_ptr_r8,'LVIS')
        end if
      else
        if (gsv_varExist(stateVector,'VIS')) then
          call gsv_getField(stateVector,vis_ptr_r8,'VIS')
        else
          call gsv_getField(stateVector,vis_ptr_r8,'LVIS')
        end if
      end if
      call gsv_getField(stateVector,lvis_ptr_r8,'LVIS')
      
      do stepIndex = 1, stateVector%numStep
        !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
        do levIndex = 1, gsv_getNumLev( stateVector,&
                                        vnl_varLevelFromVarname('LVIS'))
          do latIndex = stateVector%myLatBeg, stateVector%myLatEnd
            do lonIndex = stateVector%myLonBeg, stateVector%myLonEnd
              vis_ptr_r8(lonIndex,latIndex,levIndex,stepIndex) = &
                  exp(lvis_ptr_r8(lonIndex,latIndex,levIndex,stepIndex))
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do

    else

      if (present(statevectorOut_opt)) then
        if (gsv_varExist(stateVectorOut_opt,'VIS')) then
          call gsv_getField(statevectorOut_opt,vis_ptr_r4,'VIS')
        else
          call gsv_getField(statevectorOut_opt,vis_ptr_r4,'LVIS')
        end if
      else
        if (gsv_varExist(stateVector,'VIS')) then
          call gsv_getField(stateVector,vis_ptr_r4,'VIS')
        else
          call gsv_getField(stateVector,vis_ptr_r4,'LVIS')
        end if
      end if
      call gsv_getField(stateVector,lvis_ptr_r4,'LVIS')

      do stepIndex = 1, stateVector%numStep
        !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
        do levIndex = 1, gsv_getNumLev( stateVector,&
                                        vnl_varLevelFromVarname('LVIS'))
          do latIndex = stateVector%myLatBeg, stateVector%myLatEnd
            do lonIndex = stateVector%myLonBeg, stateVector%myLonEnd
              vis_ptr_r4(lonIndex,latIndex,levIndex,stepIndex) = &
                  exp(lvis_ptr_r4(lonIndex,latIndex,levIndex,stepIndex))
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do

    end if

    ! Change the variable name from LVIS to VIS if VIS does not exist
    if (present(statevectorOut_opt)) then
      if (.not. gsv_varExist(stateVectorOut_opt,'VIS')) then
        call gsv_modifyVarName(stateVectorOut_opt,'LVIS','VIS')
      end if
    else
      if (.not. gsv_varExist(stateVector,'VIS')) then
        call gsv_modifyVarName(stateVector,'LVIS','VIS')
      end if
    end if

  end subroutine LVIStoVIS


  !--------------------------------------------------------------------------
  ! ZandP_tl
  !--------------------------------------------------------------------------
  subroutine ZandP_tl(stateVector)
    ! 
    ! :Purpose: Tangeant linear of height and pressure computation.  
    !           The computation order depends on the native model 
    !           representation (height or pressure based).
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout) :: stateVector

    call czp_calcZandP_tl(statevector, stateVectorRefHeight)

  end subroutine ZandP_tl

  !--------------------------------------------------------------------------
  ! ZandP_ad
  !--------------------------------------------------------------------------
  subroutine ZandP_ad(stateVector)
    ! 
    ! :Purpose: Adjoint of the tangeant linear computation of both height and 
    !           pressure. The computation order depends on the native model 
    !           representation (height or pressure based).
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout) :: stateVector

    call czp_calcZandP_ad(statevector, stateVectorRefHeight)

  end subroutine ZandP_ad

  !--------------------------------------------------------------------------
  ! UVtoVortDiv_gsv
  !--------------------------------------------------------------------------
  subroutine UVtoVortDiv_gsv(statevector)
    ! 
    ! :Purpose: Wind components to relative vorticity and divergence transformation.
    !
    implicit none
   
    ! Arguments
    type(struct_gsv), intent(inout) :: statevector

    ! Locals
    real(8), pointer  :: uu_ptr(:,:,:,:), vv_ptr(:,:,:,:)
    real(8), pointer  :: qr_ptr(:,:,:,:), dd_ptr(:,:,:,:)
    integer           :: stepIndex
    integer           :: nlev_M

    call gsv_getField(statevector,uu_ptr,'UU')
    call gsv_getField(statevector,vv_ptr,'VV')
    if (gsv_varExist(statevector,'QR') .and. &
        gsv_varExist(statevector,'DD')) then
      call gsv_getField(statevector,qr_ptr,'QR')
      call gsv_getField(statevector,dd_ptr,'DD')
    else
      call gsv_getField(statevector,qr_ptr,'UU')
      call gsv_getField(statevector,dd_ptr,'VV')
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

  !--------------------------------------------------------------------------
  ! vortDivtoPsiChi_gsv
  !--------------------------------------------------------------------------
  subroutine vortDivToPsiChi_gsv(statevector)
    ! 
    ! :Purpose: Relative vorticity and divergence to stream function and 
    !           velocity potential transformation.
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout) :: statevector

    ! Locals
    integer           :: stepIndex
    logical, save     :: firstTime = .true.
    real(8), pointer  :: qr_ptr(:,:,:,:), dd_ptr(:,:,:,:)
    real(8), pointer  :: pp_ptr(:,:,:,:), cc_ptr(:,:,:,:)

    type(struct_lst), save :: lst_lapl   ! Spectral transform Parameters for Vort/Div -> Psi/Chi
    integer :: nlev_M
    integer :: nTrunc

    nTrunc = max(statevector%hco%ni,statevector%hco%nj)/4 + 1

    if (gsv_varExist(statevector,'QR') .and. &
        gsv_varExist(statevector,'DD')) then
      call gsv_getField(statevector,qr_ptr,'QR')
      call gsv_getField(statevector,dd_ptr,'DD')
    else
      call gsv_getField(statevector,qr_ptr,'UU')
      call gsv_getField(statevector,dd_ptr,'VV')
    end if

    if (gsv_varExist(statevector,'PP') .and. &
        gsv_varExist(statevector,'CC')) then
      call gsv_getField(statevector,pp_ptr,'PP')
      call gsv_getField(statevector,cc_ptr,'CC')
    else
      call gsv_getField(statevector,pp_ptr,'UU')
      call gsv_getField(statevector,cc_ptr,'VV')
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

  !--------------------------------------------------------------------------
  ! UVtoPsiChi_gsv
  !--------------------------------------------------------------------------
  subroutine UVtoPsiChi_gsv(statevector)
    ! 
    ! :Purpose: Wind components to stream function and velocity potential 
    !           transformation.
    !
    implicit none
   
    ! Arguments
    type(struct_gsv), intent(inout) :: statevector

    ! Locals
    integer               :: stepIndex 
    real(8), pointer      :: uu_ptr(:,:,:,:), vv_ptr(:,:,:,:)
    real(8), pointer      :: psi_ptr(:,:,:,:), chi_ptr(:,:,:,:)
    real(8), allocatable  :: gridState(:,:,:), spectralState(:,:,:)
    real(8)               :: dla2
    integer               :: nlev_M, levIndex
    integer               :: ila_mpiglobal, ila_mpilocal

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

      call gsv_getField(statevector,uu_ptr,'UU')
      call gsv_getField(statevector,vv_ptr,'VV')
      call gsv_getField(statevector,psi_ptr,'PP')
      call gsv_getField(statevector,chi_ptr,'CC')
      nlev_M = gsv_getNumLev(statevector,'MM')

      if ( gstID < 0 ) then
        !ntrunc = statevector%nj
        ntrunc = 180
        gstID = gst_setup(statevector%ni,statevector%nj,ntrunc,2*nlev_M)
        call mpivar_setup_m(ntrunc,mymBeg,mymEnd,mymSkip,mymCount)
        call mpivar_setup_n(ntrunc,mynBeg,mynEnd,mynSkip,mynCount)
        call gst_ilaList_mpiglobal( ilaList_mpiglobal,nla_mpilocal,maxMyNla,&
                                    gstID,mymBeg,mymEnd,mymSkip,mynBeg,mynEnd,&
                                    mynSkip)
        call gst_ilaList_mpilocal(ilaList_mpilocal,gstID,mymBeg,mymEnd,mymSkip,&
                                  mynBeg,mynEnd,mynSkip)
      end if

      dla2   = ec_ra * ec_ra
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
  
  !--------------------------------------------------------------------------
  ! UVtoPsiChi_ens
  !--------------------------------------------------------------------------
  subroutine UVtoPsiChi_ens(ens)
    ! 
    ! :Purpose: Wind components to stream function and velocity potential 
    !           transformation (ensemble processing).
    !
    implicit none
   
    ! Arguments
    type(struct_ens), intent(inout) :: ens

    ! Locals
    type(struct_hco), pointer :: hco_ens => null()
    type(struct_gsv)          :: gridStateVector_oneMember
    integer                   :: memberIndex

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
                      varNames_opt=(/'UU','VV'/), &
                      datestamp_opt=tim_getDatestamp(), &
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

  !--------------------------------------------------------------------------
  ! UVtoVorDiv_ens
  !--------------------------------------------------------------------------
  subroutine UVtoVortDiv_ens(ens)
    ! 
    ! :Purpose: Wind components to relative vorticity and divergence transformation
    !           (ensemble processing)
    !
    implicit none
   
    ! Arguments
    type(struct_ens), intent(inout) :: ens

    ! Locals
    type(struct_hco), pointer :: hco_ens => null()
    type(struct_gsv)          :: gridStateVector_oneMember
    integer                   :: memberIndex

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
                      varNames_opt=(/'UU','VV'/), &
                      datestamp_opt=tim_getDatestamp(), &
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
  ! logCH_ens
  !--------------------------------------------------------------------------
  subroutine logCH_ens(ens,varName)
    ! 
    ! :Purpose: Logarithmic transformation of a chemical species concentration
    !           (ensemble processing)
    !
    implicit none

    ! Arguments
    type(struct_ens), intent(inout) :: ens
    character(len=4), intent(in)    :: varName

    ! Locals
    integer           :: lonIndex, latIndex, levIndex, stepIndex, memberIndex
    integer           :: myLatBeg, myLatEnd, myLonBeg, myLonEnd
    character(len=4)  :: varName_ens
    real(4)           :: minVal
    real(4), pointer  :: ptr4d_r4(:,:,:,:)

    call ens_getLatLonBounds(ens, myLonBeg, myLonEnd, myLatBeg, myLatEnd)

    do levIndex = 1, ens_getNumK(ens)

      varName_ens = ens_getVarNameFromK(ens,levIndex)
      if ( trim(varName_ens) /= trim(varName) ) cycle

      ptr4d_r4 => ens_getOneLev_r4(ens,levIndex)
      minVal=real(gsv_minValVarKindCH(vnl_varListIndex(varName)),4)
    
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          do stepIndex = 1, ens_getNumStep(ens)
            do memberIndex = 1, ens_getNumMembers(ens)
              ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex) = &
                  log(max(ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex),&
                          minVal))
            end do
          end do
        end do
      end do

    end do

  end subroutine logCH_ens

  !--------------------------------------------------------------------------
  ! expCH_tlm
  !--------------------------------------------------------------------------
  subroutine expCH_tlm(statevector, varName, stateVectorRef_opt)
    ! 
    ! :Purpose: Tangent linear of exponentiation of chemical species
    !           concentration.
    !           Transform d[log(x)] to dx where x = 'stateVectorRef_opt',
    !           the input 'statevector' component is d[log(x)] and
    !           the output 'statevector' component is dx.
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout)         :: statevector
    character(len=*), intent(in)            :: varName
    type(struct_gsv), optional, intent(in)  :: statevectorRef_opt
    
    ! Locals
    integer           :: lonIndex,latIndex,levIndex,stepIndex,varIndex
    real(8), pointer  :: var_ptr(:,:,:,:), logVar_ptr(:,:,:,:)
    real(8), pointer  :: var_trial(:,:,:,:)
    real(8)           :: minVal

    if ( present(statevectorRef_opt) ) then
       call gsv_getField(stateVectorRef_opt,var_trial,trim(varName))
    else
      varIndex = vnl_varListIndex(varName)
      if ( .not. varKindCHTrialsInitialized(varIndex) ) then
        call utl_abort('expCH_tlm: do trials to stateVectorRefChem transform at higher level')
      end if
      call gsv_getField(stateVectorTrialvarKindCH(varIndex),var_trial,&
                        trim(varName))
    end if

    call gsv_getField(statevector,var_ptr,trim(varName))
    call gsv_getField(statevector,logVar_ptr,trim(varName))

    minVal=gsv_minValVarKindCH(vnl_varListIndex(varName))

    do stepIndex = 1, statevector%numStep
      !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
      do levIndex = 1, gsv_getNumLev( statevector,&
                                      vnl_varLevelFromVarname(trim(varName)))
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do lonIndex = statevector%myLonBeg, statevector%myLonEnd       
            var_ptr(lonIndex,latIndex,levIndex,stepIndex) =  &
                logVar_ptr(lonIndex,latIndex,levIndex,stepIndex) &
                * max(var_trial(lonIndex,latIndex,levIndex,stepIndex),minVal)
          end do
        end do
      end do
      !$OMP END PARALLEL DO
    end do

  end subroutine expCH_tlm

  !--------------------------------------------------------------------------
  ! CH_bounds
  !--------------------------------------------------------------------------
  subroutine CH_bounds(statevector)
    ! 
    ! :Purpose: Impose boundary values to variables of CH kind.
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout) :: statevector
    
    ! Locals
    integer           :: varIndex,lonIndex,latIndex,levIndex,stepIndex
    real(8), pointer  :: var_ptr(:,:,:,:)
    real(8)           :: minVal
    character(len=4)  :: varName

    do varIndex = 1, vnl_numvarmax
      varName = vnl_varNameList(varIndex)
      if ( .not.gsv_varExist(varName=trim(varName)) ) cycle
      if ( vnl_varKindFromVarname(trim(varName)) /= 'CH' ) cycle

      call gsv_getField(statevector,var_ptr,trim(varName))

      minVal=gsv_minValVarKindCH(vnl_varListIndex(varName))

      do stepIndex = 1, statevector%numStep
        !$OMP PARALLEL DO PRIVATE(lonIndex,latIndex,levIndex)
        do levIndex = 1, gsv_getNumLev( statevector,&
                                        vnl_varLevelFromVarname(trim(varName)))
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd       
              var_ptr(lonIndex,latIndex,levIndex,stepIndex) = &
                  max(var_ptr(lonIndex,latIndex,levIndex,stepIndex),minVal)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do

    end do

  end subroutine CH_bounds

  !--------------------------------------------------------------------------
  ! GLtoLG
  !--------------------------------------------------------------------------
  subroutine GLtoLG(stateVector, stateVectorRef)
    !
    ! :Purpose: Solve laplaces equation at a subset of gridpoints
    !           subject to the boundary conditions imposed by the
    !           surrounding points.  Uses the method of sequential 
    !           relaxation or Liebmann relaxation, which converges
    !           more rapidly than the simultaneous relaxation (see
    !           numerical weather analysis and prediction by P. D. 
    !           Thompson, 1961, pp92-98). NOTE: this subroutine
    !           currently uses the oceanMask only for the first level.
    !           Therefore, if it is applied to 3D masked fields that
    !           have level-dependent masks, the code will abort so
    !           that the user can make the necessary adjustments to
    !           ensure the correct behaviour.
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout) :: statevector
    type(struct_gsv), intent(in)    :: stateVectorRef

    ! Locals
    type(struct_gsv) :: statevector_analysis_1step_r8
    type(struct_gsv) :: statevector_trial_1step_r8

    real(8), pointer  :: LGAnal_ptr(:,:,:,:), LGTrial_ptr(:,:,:,:)
    real(8), pointer  :: GL_ptr(:,:,:,:)
    real(8)           :: alpha, factor, correc, rms, maxAbsCorr, basic

    integer :: lonIndex, latIndex, levIndex, stepIndex
    integer :: ipass, numPass, numCorrect

    logical :: orca12

    ! abort if 3D mask is present, since we may not handle this situation correctly
    if ( stateVector%oceanMask%nLev > 1 ) then
      call utl_abort('gvt_GLtoLG: 3D mask present - this case not properly handled')
    end if

    ! allocate statevector for single time steps
    if ( mpi_myid < stateVector%numStep ) then
      call gsv_allocate(stateVector_analysis_1step_r8, 1, stateVector%hco, &
                        stateVector%vco, mpi_local_opt=.false., &
                        dataKind_opt=8, varNames_opt=(/'GL','LG'/) )
      call gsv_allocate(stateVector_trial_1step_r8, 1, stateVectorRef%hco, &
                        stateVectorRef%vco, mpi_local_opt=.false., &
                        dataKind_opt=8, varNames_opt=(/'GL','LG'/) )
    end if

    call gsv_transposeTilesToStep(stateVector_analysis_1step_r8, stateVector, 1)
    call gsv_transposeTilesToStep(stateVector_trial_1step_r8, stateVectorRef, 1)

    if ( stateVector_analysis_1step_r8%allocated ) then

      call gsv_getField(stateVector_analysis_1step_r8,     GL_ptr,'GL')
      call gsv_getField(stateVector_analysis_1step_r8, LGAnal_ptr,'LG')
      call gsv_getField(stateVector_trial_1step_r8,   LGTrial_ptr,'LG')

      alpha = 1.975d0
      if ( stateVector%ni == 4322 .and. stateVector%nj == 3059 ) then
        orca12 = .true.
        numPass = 1000
        factor = alpha*0.88d0
      else
        orca12 = .false.
        numPass = 500
        factor = alpha
      end if

      write(*,*) 'gvt_GLtoLG Liebmann relaxation'
      write(*,*) 'gvt_GLtoLG Number of free points: ',  count(.not. stateVector%oceanMask%mask)
      write(*,*) 'gvt_GLtoLG Number of fixed points: ', count(      stateVector%oceanMask%mask)
      write(*,*) 'gvt_GLtoLG Total number of grid points: ', stateVector%ni*stateVector%nj
      write(*,*) 'gvt_GLtoLG Total number of iterations: ', numPass

      do stepIndex = 1, statevector%numStep
        write(*,*) 'gvt_GLtoLG stepIndex = ',stepIndex
        do levIndex = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname('LG'))

          write(*,*) 'gvt_GLtoLG levIndex = ',levIndex

          ! Initialisation
          do latIndex = 1, stateVector%nj
            do lonIndex = 1, stateVector%ni
              if ( stateVector%oceanMask%mask(lonIndex,latIndex,1) ) then
                LGAnal_ptr(lonIndex,latIndex,levIndex,stepIndex) = &
                    GL_ptr(lonIndex,latIndex,levIndex,stepIndex)
              else
                LGAnal_ptr(lonIndex,latIndex,levIndex,stepIndex) = &
                    LGTrial_ptr(lonIndex,latIndex,levIndex,stepIndex)
              end if
            end do
          end do

          do ipass=1,numPass

            rms = 0.0d0
            numCorrect = 0
            maxAbsCorr = 0.0d0
            do latIndex = 2, stateVector%nj-1
              do lonIndex = 2, stateVector%ni-1
                if ( .not. stateVector%oceanMask%mask(lonIndex,latIndex,1) ) then
                  basic = ( LGAnal_ptr(lonIndex+1,latIndex,levIndex,stepIndex)+&
                            LGAnal_ptr(lonIndex-1,latIndex,levIndex,stepIndex)+&
                            LGAnal_ptr(lonIndex,latIndex+1,levIndex,stepIndex)+&
                            LGAnal_ptr(lonIndex,latIndex-1,levIndex,stepIndex)&
                            ) / 4.0d0
                  correc = factor*(basic - LGAnal_ptr(lonIndex,latIndex,&
                                                      levIndex,stepIndex))
                  LGAnal_ptr(lonIndex,latIndex,levIndex,stepIndex) = &
                      LGAnal_ptr(lonIndex,latIndex,levIndex,stepIndex) + correc
                  rms = rms + correc*correc
                  numCorrect = numCorrect + 1
                  if( abs(correc) > maxAbsCorr ) then
                    maxAbsCorr = abs(correc)
                  end if
                end if
              end do
            end do

            if( orca12 ) then
              ! Periodicity in the X direction
              lonIndex = 1
              do latIndex = 2, stateVector%nj-1
                if ( .not. stateVector%oceanMask%mask(lonIndex,latIndex,1) ) then
                  basic = ( LGAnal_ptr(lonIndex+1, latIndex, levIndex,&
                                        stepIndex) + &
                            LGAnal_ptr( stateVector%ni-2, latIndex, levIndex,&
                                        stepIndex) + &
                            LGAnal_ptr( lonIndex, latIndex+1, levIndex, &
                                        stepIndex) + &
                            LGAnal_ptr( lonIndex, latIndex-1, levIndex, &
                                        stepIndex) ) / 4.0d0
                  correc = factor*(basic - LGAnal_ptr(lonIndex,latIndex,&
                                                      levIndex,stepIndex))
                  LGAnal_ptr(lonIndex,latIndex,levIndex,stepIndex) = &
                      LGAnal_ptr(lonIndex,latIndex,levIndex,stepIndex) + correc
                  rms = rms + correc*correc
                  numCorrect = numCorrect + 1
                  if( abs(correc) > maxAbsCorr ) then
                    maxAbsCorr = abs(correc)
                  end if
                end if
              end do
              lonIndex = stateVector%ni
              do latIndex = 2, stateVector%nj-1
                if ( .not. stateVector%oceanMask%mask(lonIndex,latIndex,1) ) then
                  basic = ( LGAnal_ptr(3, latIndex,  levIndex,stepIndex) + &
                            LGAnal_ptr(lonIndex-1,latIndex,levIndex,stepIndex)+&
                            LGAnal_ptr(lonIndex,latIndex+1,levIndex,stepIndex)+&
                            LGAnal_ptr(lonIndex,latIndex-1,levIndex,stepIndex)&
                            ) / 4.0d0
                  correc = factor*(basic - LGAnal_ptr(lonIndex,latIndex,&
                                                      levIndex,stepIndex))
                  LGAnal_ptr(lonIndex,latIndex,levIndex,stepIndex) = &
                      LGAnal_ptr(lonIndex,latIndex,levIndex,stepIndex) + correc
                  rms = rms + correc*correc
                  numCorrect = numCorrect + 1
                  if( abs(correc) > maxAbsCorr ) then
                    maxAbsCorr = abs(correc)
                  end if
                end if
              end do
              ! North fold
              latIndex = stateVector%nj
              do lonIndex = 2, stateVector%ni-1
                if ( .not. stateVector%oceanMask%mask(lonIndex,latIndex,1) ) then
                  basic = ( LGAnal_ptr(lonIndex+1,latIndex,levIndex,stepIndex)+&
                            LGAnal_ptr(lonIndex-1,latIndex,levIndex,stepIndex)+&
                            LGAnal_ptr( stateVector%ni+2-lonIndex, latIndex-2,&
                                        levIndex,stepIndex) + &
                            LGAnal_ptr( stateVector%ni+2-lonIndex, latIndex-1,&
                                        levIndex,stepIndex) ) / 4.0d0
                  correc = factor*(basic - LGAnal_ptr(lonIndex,latIndex, &
                                                      levIndex,stepIndex))
                  LGAnal_ptr(lonIndex,latIndex,levIndex,stepIndex) = &
                      LGAnal_ptr(lonIndex,latIndex,levIndex,stepIndex) + correc
                  rms = rms + correc*correc
                  numCorrect = numCorrect + 1
                  if( abs(correc) > maxAbsCorr ) then
                    maxAbsCorr = abs(correc)
                  end if
                end if
              end do
            end if

          end do

          if( numCorrect > 0 ) rms = sqrt(rms/real(numCorrect))

          write(*,*) 'gvt_GLtoLG number of points corrected = ',numCorrect
          write(*,*) 'gvt_GLtoLG RMS correction during last iteration: ',rms
          write(*,*) 'gvt_GLtoLG MAX absolute correction during last iteration: ', maxAbsCorr
          write(*,*) 'gvt_GLtoLG Field min value = ',minval(LGAnal_ptr(:,:,levIndex,stepIndex))
          write(*,*) 'gvt_GLtoLG Field max value = ',maxval(LGAnal_ptr(:,:,levIndex,stepIndex))

          if( maxAbsCorr > 1.0 ) then
            call utl_abort('gvt_GLtoLG Unstable algorithm !')
          end if

        end do
      end do

      !
      !- Impose limits [0,1] on sea ice concentration analysis
      !
      LGAnal_ptr(:,:,:,:) = min(LGAnal_ptr(:,:,:,:), 1.0d0)
      LGAnal_ptr(:,:,:,:) = max(LGAnal_ptr(:,:,:,:), 0.0d0)

    end if

    call gsv_transposeStepToTiles(stateVector_analysis_1step_r8,stateVector,1)

    if ( mpi_myid < stateVector%numStep ) then
      call gsv_deallocate(stateVector_analysis_1step_r8)
      call gsv_deallocate(stateVector_trial_1step_r8)
    end if

  end subroutine GLtoLG

end module gridVariableTransforms_mod

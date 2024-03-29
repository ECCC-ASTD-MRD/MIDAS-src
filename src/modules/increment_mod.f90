
module increment_mod
  ! MODULE increment_mod (prefix='inc' category='1. High-level functionality')
  !
  !:Purpose:  To add a 4D increment to a given 4D background/reference state and
  !           to output the results
  !
  use codePrecision_mod
  use midasMpi_mod
  use timeCoord_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use interpolation_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use humidityLimits_mod
  use utilities_mod
  use message_mod
  use gridVariableTransforms_mod
  use BMatrix_mod
  use varNamelist_mod
  implicit none
  save
  private

  ! public procedures
  public :: inc_computeHighResAnalysis, inc_writeIncAndAnalHighRes, inc_getIncrement, inc_writeIncrement
  public :: inc_writeAnalysis, inc_analPostProcessing

  ! namelist variables
  integer           :: writeNumBits          ! number of bits to use when writing analysis and high-res increment
  logical           :: writeHiresIncrement   ! choose to write the high-res increment to a file
  logical           :: imposeRttovHuLimits   ! choose to impose "rttov" humidity limits to analysis
  logical           :: useAnalIncMask        ! for LAM only, choose to apply scale factor from a mask file to the increment  
  character(len=12) :: etiket_anlm           ! 'etiket' used when writing the analysis
  character(len=12) :: etiket_rehm           ! 'etiket' used when writing the high-res increment
  character(len=12) :: etiket_rebm           ! 'etiket' used when writing the low-res increment
  character(len=12) :: hInterpolationDegree  ! type of interpolation to use: 'LINEAR' or 'CUBIC'
  logical           :: applyLiebmann         ! choose to apply Liebmann extrapolation to SST and/or sea ice
  logical           :: SSTSpread             ! choose to apply spatial spreading of the SST increment onto land
  integer           :: SSTSpreadMaxBoxSize   ! control the amount of SST increment spreading
  character(len=10) :: SSTSubgrid            ! select subgrid on which to apply spreading: 'Yin' or 'Yan'

CONTAINS

  !--------------------------------------------------------------------------
  ! readNameList
  !--------------------------------------------------------------------------
  subroutine readNameList
    !
    ! :Purpose: Reading NAMINC namelist by any subroutines in increment_mod module.
    !
    implicit none

    ! Locals:
    integer :: nulnam, ierr
    integer, external :: fnom, fclos
    logical, save :: nmlAlreadyRead = .false.

    NAMELIST /NAMINC/ writeHiresIncrement, etiket_rehm, etiket_anlm, &
         etiket_rebm, writeNumBits, imposeRttovHuLimits, hInterpolationDegree, &
         useAnalIncMask, applyLiebmann, SSTSpread, SSTSpreadMaxBoxSize, SSTSubgrid

    if ( .not. nmlAlreadyRead ) then
      nmlAlreadyRead = .true.

      !- Setting default values
      writeHiresIncrement = .true.
      imposeRttovHuLimits = .true.
      useAnalIncMask      = .false.
      etiket_rehm = 'INCREMENT'
      etiket_rebm = 'INCREMENT'
      etiket_anlm = 'ANALYSIS'
      writeNumBits = 16
      hInterpolationDegree = 'LINEAR'
      applyLiebmann = .false.
      SSTSpread  = .false.
      SSTSpreadMaxBoxSize = 0
      SSTsubgrid = '   '

      if ( .not. utl_isNamelistPresent('NAMINC','./flnml') ) then
        call msg('readNameList (inc)', &
             'NAMINC is missing in the namelist. The default values will be taken.', &
             mpiAll_opt=.false.)
      else
        ! Reading the namelist
        nulnam = 0
        ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
        read(nulnam, nml=naminc, iostat=ierr)
        if ( ierr /= 0) call utl_abort('readNameList: Error reading namelist')
        ierr = fclos(nulnam)
      end if
      if ( mmpi_myid == 0 ) write(*,nml=naminc)
    end if

  end subroutine readNameList

  !--------------------------------------------------------------------------
  ! inc_computeHighResAnalysis
  !--------------------------------------------------------------------------
  subroutine inc_computeHighResAnalysis( statevectorIncLowRes,                            & ! IN
                                         stateVectorUpdateHighRes, stateVectorPsfcHighRes)  ! OUT
    !
    ! :Purpose: Computing high-resolution analysis on the trial grid.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)    :: statevectorIncLowRes
    type(struct_gsv), intent(inout) :: stateVectorUpdateHighRes
    type(struct_gsv), intent(out)   :: stateVectorPsfcHighRes

    ! Locals:
    type(struct_gsv) :: statevectorRef, statevectorRefPsfc
    type(struct_gsv) :: statevectorIncRefLowRes
    type(struct_gsv) :: statevectorTrlRefVars
    type(struct_gsv) :: statevectorTrlLowResTime
    type(struct_gsv) :: statevectorTrlLowResVert
    type(struct_gsv) :: statevector_mask
    type(struct_gsv) :: stateVectorHighRes

    type(struct_vco), pointer :: vco_trl, vco_inc
    type(struct_hco), pointer :: hco_trl, hco_inc

    real(pre_incrReal), pointer :: PsfcIncPrior(:,:,:,:), PsfcIncMasked(:,:,:,:)
    real(pre_incrReal), pointer :: analIncMask(:,:,:)

    integer              :: stepIndex, numStep_inc, numStep_trl
    integer, allocatable :: dateStampList(:)

    character(len=4), allocatable :: varNamesRef(:), varNamesPsfc(:)
    logical  :: allocHeightSfc, refStateNeeded

    call msg('inc_computeHighResAnalysis', 'START', verb_opt=2)
    call msg_memUsage('inc_computeHighResAnalysis')

    call utl_tmg_start(80,'--Increment')
    call utl_tmg_start(81,'----ComputeHighResAnalysis')

    ! Set/Read values for the namelist NAMINC
    call readNameList

    if ( gsv_isAllocated(stateVectorPsfcHighRes) ) then
      call utl_abort('inc_computeHighResAnalysis: '&
                      //'stateVectorPsfcHighRes should not be allocated yet')
    end if

    ! If P0 present it is infered the statevector is a 3D atmospheric state
    refStateNeeded =  gsv_varExist(varName='P0')

    ! Setup timeCoord module (date read from trial file)
    numStep_inc = tim_nstepobsinc ! low-res time
    numStep_trl = tim_nstepobs    ! high-res time
    allocate(dateStampList(numStep_inc))
    call tim_getstamplist(dateStampList,numStep_inc,tim_getDatestamp())

    hco_trl => gsv_getHco(stateVectorUpdateHighRes)
    vco_trl => gsv_getVco(stateVectorUpdateHighRes)
    hco_inc => gsv_getHco(statevectorIncLowRes)
    vco_inc => gsv_getVco(statevectorIncLowRes)

    if (vco_trl%Vcode == 0 .or. .not. gsv_varExist(varName='P0')) then
      allocHeightSfc = .false.
    else
      allocHeightSfc = stateVectorUpdateHighRes%heightSfcPresent
    end if

    ! Read the analysis mask (in LAM mode only) - N.B. different from land/sea mask!!!
    if (.not. hco_trl%global .and. useAnalIncMask) then
      call gio_getMaskLAM(statevector_mask, hco_trl, vco_trl, hInterpolationDegree)
    end if

    ref_building: if ( refStateNeeded ) then
      ! Build a reference variables
      !
      ! - Allocate the target reference statevector for vertical interpolation
      ! that will be done in inc_interpolateAndAdd.
      ! The reference needs to have the vertical structure of the input which is the
      ! increment, but since horizontal interpolation is done first, it needs the
      ! trial horizontal structure.
      if (vco_trl%Vcode == 5100) then
        allocate(varNamesRef(4))
        varNamesRef = (/'P0', 'P0LS', 'TT', 'HU'/)
        allocate(varNamesPsfc(2))
        varNamesPsfc = (/'P0','P0LS'/)
      else
        allocate(varNamesRef(3))
        varNamesRef = (/'P0', 'TT', 'HU'/)
        allocate(varNamesPsfc(1))
        varNamesPsfc = (/'P0'/)
      end if
      call gsv_allocate( statevectorRef, numStep_inc, hco_trl, vco_inc, &
                         dataKind_opt=pre_incrReal, &
                         dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true.,  &
                         varNames_opt=varNamesRef, allocHeightSfc_opt=allocHeightSfc, &
                         hInterpolateDegree_opt=hInterpolationDegree )

      ! - Restriction of increment to only reference state variables
      call gsv_allocate( statevectorIncRefLowRes, numStep_inc, hco_inc, vco_inc,  &
                         dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true.,  &
                         dataKind_opt=pre_incrReal, varNames_opt=varNamesRef, &
                         allocHeightSfc_opt=allocHeightSfc)
      call gsv_copy(statevectorIncLowRes, statevectorIncRefLowRes, &
                    allowVarMismatch_opt=.true.)

      ! - Horizontal *only* interpolation of analysis increment (ref state variables)
      !   using int_interp_gsv since it also does needed mpi transpositions
      call msg('inc_computeHighResAnalysis', &
               'horizontal interpolation of the Psfc increment', mpiAll_opt=.false.)
      call int_interp_gsv(statevectorIncRefLowRes,statevectorRef)
      call gsv_deallocate(statevectorIncRefLowRes)
      !- Apply mask to increment in LAM
      if (.not. hco_trl%global .and. useAnalIncMask) then
        call gsv_getField(statevectorRef,PsfcIncPrior,'P0')
        call gsv_getField(statevectorRef,PsfcIncMasked,'P0')
        call gsv_getField(statevector_mask,analIncMask)
        do stepIndex = 1, statevectorRef%numStep
          PsfcIncMasked(:,:,1,stepIndex) = PsfcIncPrior(:,:,1,stepIndex) * analIncMask(:,:,1)
        end do
        if (gsv_varExist(statevectorRef,'P0LS')) then
          call gsv_getField(statevectorRef,PsfcIncPrior,'P0LS')
          call gsv_getField(statevectorRef,PsfcIncMasked,'P0LS')
          do stepIndex = 1, statevectorRef%numStep
            PsfcIncMasked(:,:,1,stepIndex) = PsfcIncPrior(:,:,1,stepIndex) * analIncMask(:,:,1)
          end do
        end if
      end if

      ! - Compute analysis of reference state variables to use for vertical interpolation
      !   of increment
      call msg('inc_computeHighResAnalysis', &
           'Computing reference variables analysis to use for interpolation of increment', &
           mpiAll_opt=.false.)

      !   - build a restriction to reference state variables of trial
      call gsv_allocate(statevectorTrlRefVars, numStep_trl, &
                        hco_trl, vco_trl, dateStamp_opt=tim_getDateStamp(), &
                        mpi_local_opt=.true., dataKind_opt=pre_incrReal, &
                        varNames_opt=varNamesRef, allocHeightSfc_opt=allocHeightSfc )
      call gsv_copy(stateVectorUpdateHighRes, statevectorTrlRefVars,&
                    allowVarMismatch_opt=.true.)
      !   - bring trial to low time res
      call gsv_allocate(statevectorTrlLowResTime, numStep_inc, hco_trl, vco_trl,  &
                        dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true.,  &
                        dataKind_opt=pre_incrReal, varNames_opt=varNamesRef, &
                        allocHeightSfc_opt=allocHeightSfc)
      call gsv_copy(statevectorTrlRefVars, statevectorTrlLowResTime, &
                    allowTimeMismatch_opt=.true.)
      !   - bring trial to low res (increment) vertical structure
      !     (no vertical interpolation reference needed as it is a full state,
      !     not a incremental state)
      call gsv_allocate(statevectorTrlLowResVert, numStep_inc, hco_trl, vco_inc,  &
                        dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true.,  &
                        dataKind_opt=pre_incrReal, varNames_opt=varNamesRef, &
                        allocHeightSfc_opt=allocHeightSfc)
      call int_vInterp_gsv(statevectorTrlLowResTime, statevectorTrlLowResVert)
      call gsv_deallocate(statevectorTrlLowResTime)
      call gsv_deallocate(statevectorTrlRefVars)
      !   - build the analysis reference state for the increment vertical interpolation.
      !     At that stage, statevectorRef contains the increment on the trial
      !     horizontal grid, but analysis vertical structure; consistent for addition.
      call gsv_add(statevectorTrlLowResVert, statevectorRef)
      call gsv_deallocate(statevectorTrlLowResVert)

      ! - Time interpolation to get high-res Psfc analysis increment to output
      call msg('inc_computeHighResAnalysis', &
               'Time interpolation to get high-res Psfc analysis increment', &
               mpiAll_opt=.false.)

      !   - Restriction to P0 only + vco set on vco_trl for later consistency
      call gsv_allocate(statevectorRefPsfc, numStep_inc, hco_trl, vco_trl, &
                        dataKind_opt=pre_incrReal, &
                        dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true.,  &
                        varNames_opt=varNamesPsfc, allocHeightSfc_opt=allocHeightSfc, &
                        hInterpolateDegree_opt=hInterpolationDegree )
      call gsv_copy(statevectorRef, statevectorRefPsfc, allowVarMismatch_opt=.true., &
                    allowVcoMismatch_opt=.true.)
      !   - high res time interpolation
      call gsv_allocate(stateVectorPsfcHighRes, numStep_trl, hco_trl, vco_trl, &
                        dataKind_opt=pre_incrReal, &
                        dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true.,  &
                        varNames_opt=varNamesPsfc, allocHeightSfc_opt=allocHeightSfc, &
                        hInterpolateDegree_opt=hInterpolationDegree)
      call int_tInterp_gsv(statevectorRefPsfc, stateVectorPsfcHighRes)
    else ref_building
      call msg('inc_computeHighResAnalysis','P0 not present, not computing reference', &
               verb_opt=3)
    end if ref_building

    ! Compute the analysis
    call msg('inc_computeHighResAnalysis', 'compute the analysis', mpiAll_opt=.false.)

    ! Copy of trial to be passed to inc_interpolateAndAdd
    call gsv_allocate( stateVectorHighRes, numStep_trl, hco_trl, vco_trl, &
                       dataKind_opt=stateVectorUpdateHighRes%dataKind, &
                       dateStamp_opt=tim_getDateStamp(), &
                       mpi_local_opt=stateVectorUpdateHighRes%mpi_local, &
                       allocHeightSfc_opt=stateVectorUpdateHighRes%heightSfcPresent, &
                       hInterpolateDegree_opt=hInterpolationDegree, &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_copy( stateVectorUpdateHighRes,stateVectorHighRes, &
                   allowVarMismatch_opt=.true.)

    ! Interpolate low-res increments to high-res and add to the initial state
    if (.not. hco_trl%global .and. useAnalIncMask) then
      if ( refStateNeeded ) then
        call inc_interpolateAndAdd(statevectorIncLowRes, stateVectorHighRes, &
                                   statevectorRef_opt=statevectorRef, &
                                   statevectorMaskLAM_opt=statevector_mask)
      else
        call inc_interpolateAndAdd(statevectorIncLowRes, stateVectorHighRes, &
                                   statevectorMaskLAM_opt=statevector_mask)
      end if
    else
      if ( refStateNeeded ) then
        call inc_interpolateAndAdd(statevectorIncLowRes, stateVectorHighRes, &
                                   statevectorRef_opt=statevectorRef)
      else
        call inc_interpolateAndAdd(statevectorIncLowRes, stateVectorHighRes)
      end if
    end if

    call gsv_copy( stateVectorHighRes, stateVectorUpdateHighRes, &
                   allowVarMismatch_opt=.true.)
    call gsv_deallocate(stateVectorHighRes)

    if ( gsv_isAllocated(statevectorRef) ) call gsv_deallocate(statevectorRef)
    if ( gsv_isAllocated(statevector_mask) ) call gsv_deallocate(statevector_mask)

    call utl_tmg_stop(81)
    call utl_tmg_stop(80)

    call msg('inc_computeHighResAnalysis', 'END', verb_opt=2)

  end subroutine inc_computeHighResAnalysis

  !--------------------------------------------------------------------------
  ! inc_analPostProcessing
  !--------------------------------------------------------------------------
  subroutine inc_analPostProcessing (stateVectorPsfcHighRes, stateVectorUpdateHighRes, & ! IN
                                     stateVectorTrial, stateVectorPsfc, stateVectorAnal) ! OUT
    !
    ! :Purpose: Post processing of the high resolution analysis including degrading
    !           temporal resolution, variable transform, and humidity clipping.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout)  :: stateVectorPsfcHighRes
    type(struct_gsv), intent(in)     :: stateVectorUpdateHighRes
    type(struct_gsv), intent(out)    :: stateVectorTrial
    type(struct_gsv), intent(out)    :: stateVectorPsfc
    type(struct_gsv), intent(out)    :: stateVectorAnal

    ! Locals:
    type(struct_vco), pointer :: vco_trl => null()
    type(struct_hco), pointer :: hco_trl => null()
    character(len=4), allocatable :: varNamesPsfc(:)

    real(pre_incrReal), pointer :: oceanIce_ptr(:,:,:,:)

    logical  :: allocHeightSfc

    call msg('inc_analPostProcessing', 'START', verb_opt=2)
    call msg_memUsage('inc_analPostProcessing')

    call utl_tmg_start(80,'--Increment')
    call utl_tmg_start(82,'----AnalPostProcessing')

    ! Re-read trials to make stateVectorTrial with degraded timesteps available
    hco_trl => gsv_getHco(stateVectorUpdateHighRes)
    vco_trl => gsv_getVco(stateVectorUpdateHighRes)
    allocHeightSfc = ( vco_trl%Vcode /= 0 )

    call gsv_allocate(stateVectorTrial, tim_nstepobsinc, hco_trl, vco_trl,  &
                      dataKind_opt = pre_incrReal, &
                      dateStamp_opt = tim_getDateStamp(), mpi_local_opt = .true., &
                      allocHeightSfc_opt = allocHeightSfc, hInterpolateDegree_opt = 'LINEAR', &
                      allocHeight_opt = .false., allocPressure_opt = .false.)
    call gsv_zero(stateVectorTrial)
    call gio_readTrials(stateVectorTrial)
    call msg_memUsage('inc_analPostProcessing')

    if (vco_trl%Vcode == 0 .or. .not. gsv_varExist(varName='P0')) then
      allocHeightSfc = .false.
    else
      allocHeightSfc = stateVectorTrial%heightSfcPresent
    end if

    if (vco_trl%Vcode == 5100) then
      allocate(varNamesPsfc(2))
      varNamesPsfc = (/'P0','P0LS'/)
    else
      allocate(varNamesPsfc(1))
      varNamesPsfc = (/'P0'/)
    end if

    ! Degrade timesteps stateVector high-res Psfc, and high-res analysis
    if (gsv_varExist(varName='P0')) then
      call gsv_allocate(stateVectorPsfc, tim_nstepobsinc, hco_trl, vco_trl, &
                        dataKind_opt = pre_incrReal, &
                        dateStamp_opt = tim_getDateStamp(), mpi_local_opt = .true.,  &
                        varNames_opt = varNamesPsfc, allocHeightSfc_opt = allocHeightSfc, &
                        hInterpolateDegree_opt = hInterpolationDegree)
      call gsv_copy(stateVectorPsfcHighRes, stateVectorPsfc, &
                    allowTimeMismatch_opt = .true., allowVarMismatch_opt=.true.)
      call gsv_deallocate(stateVectorPsfcHighRes)
    end if

    call gsv_allocate(stateVectorAnal, tim_nstepobsinc, hco_trl, vco_trl,  &
                      dataKind_opt = pre_incrReal, &
                      dateStamp_opt = tim_getDateStamp(), mpi_local_opt = .true., &
                      allocHeightSfc_opt = allocHeightSfc, hInterpolateDegree_opt = 'LINEAR', &
                      allocHeight_opt = .false., allocPressure_opt = .false.)
    call gsv_copy(stateVectorUpdateHighRes, stateVectorAnal, &
                  allowVarMismatch_opt = .true., allowTimeMismatch_opt = .true.)

    if (gsv_varExist(stateVectorAnal, 'GL')) then
      ! Impose limits [0,1] on sea ice concentration analysis
      call gsv_getField(stateVectorAnal, oceanIce_ptr, 'GL')
      oceanIce_ptr(:,:,:,:) = min(oceanIce_ptr(:,:,:,:), 1.0d0)
      oceanIce_ptr(:,:,:,:) = max(oceanIce_ptr(:,:,:,:), 0.0d0)
    end if

    if (applyLiebmann) then

      ! Start the variable transformations
      if (gsv_varExist(stateVectorAnal, 'GL')) then
        if (gsv_varExist(stateVectorAnal, 'LG')) then
          ! Compute the continuous sea ice concentration field (LG)
          call gvt_transform(stateVectorAnal, 'oceanIceContinuous', stateVectorRef_opt = stateVectorTrial, varName_opt = 'LG')
        end if
      end if

      ! Start the variable transformations
      if(gsv_varExist(stateVectorAnal,'TM')) then
        ! Compute the continuous SST field (TM)
        call gvt_transform(stateVectorAnal, 'oceanIceContinuous', stateVectorRef_opt = stateVectorTrial, varName_opt = 'TM')
      end if

    end if

    if (SSTSpread) then
      if(gsv_varExist(stateVectorAnal,'TM')) then
        ! Compute SST spreading (TM) on neigbouring land points
        call gvt_transform(stateVectorAnal, 'SSTSpread', varName_opt = 'TM', maxBoxSize_opt = SSTSpreadMaxBoxSize, subgrid_opt = SSTSubgrid)
      end if
    end if

    ! Convert all transformed variables into model variables (e.g. LVIS->VIS, LPR->PR) for analysis
    call gvt_transform(stateVectorAnal, 'AllTransformedToModel', allowOverWrite_opt = .true.)

    ! Impose limits on humidity analysis
    call msg('inc_analPostProcessing', 'calling qlim_saturationLimit')
    call qlim_saturationLimit(stateVectorAnal)
    if (imposeRttovHuLimits) call qlim_rttovLimit(stateVectorAnal)

    if (gsv_varKindExist('CH')) then
      ! Apply boundaries to analysis of CH kind variables as needed.
      call msg('inc_analPostProcessing', &
           'applying minimum values to analysis for variables of CH kind')
      call gvt_transform(stateVectorAnal, 'CH_bounds')
    end if

    call utl_tmg_stop(82)
    call utl_tmg_stop(80)

    call msg('inc_analPostProcessing', 'END', verb_opt=2)

  end subroutine inc_analPostProcessing

  !--------------------------------------------------------------------------
  ! inc_writeIncAndAnalHighRes
  !--------------------------------------------------------------------------
  subroutine inc_writeIncAndAnalHighRes( stateVectorTrial, stateVectorPsfc, &
                                         stateVectorAnal )
    !
    ! :Purpose: Write the high-resolution analysis increments to the rehm file.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), target, intent(inout) :: stateVectorTrial
    type(struct_gsv),         intent(inout) :: stateVectorPsfc
    type(struct_gsv),         intent(inout) :: stateVectorAnal

    ! Locals:
    type(struct_gsv) :: stateVectorIncHighRes
    type(struct_gsv) :: stateVector_1step_r4, stateVectorPsfc_1step_r4
    character(len=4), allocatable :: varNamesPsfc(:)

    type(struct_vco), pointer :: vco_trl => null()
    type(struct_hco), pointer :: hco_trl => null()

    integer              :: stepIndex, stepIndexBeg, stepIndexEnd, stepIndexToWrite, numStep
    integer              :: dateStamp, numBatch, batchIndex, procToWrite
    integer, allocatable :: dateStampList(:)

    character(len=256)  :: incFileName, anlFileName
    character(len=4)    :: coffset
    character(len=4), pointer :: varNames(:)

    real(8)             :: deltaHours
    real(8), pointer    :: HeightSfc_increment(:,:), HeightSfc_trial(:,:)

    logical  :: allocHeightSfc, writeHeightSfc

    call msg('inc_writeIncAndAnalHighRes', 'START', verb_opt=2)
    call msg_memUsage('inc_writeIncAndAnalHighRes')

    call utl_tmg_start(80,'--Increment')
    call utl_tmg_start(83,'----WriteIncAndAnalHighRes')

    ! Setup timeCoord module (date read from trial file)
    numStep = tim_nstepobsinc
    allocate(dateStampList(numStep))
    call tim_getstamplist(dateStampList,numStep,tim_getDatestamp())

    hco_trl => gsv_getHco(stateVectorTrial)
    vco_trl => gsv_getVco(stateVectorTrial)
    if (vco_trl%Vcode == 0 .or. .not. gsv_varExist(varName='P0')) then
      allocHeightSfc = .false.
    else
      allocHeightSfc = stateVectorTrial%heightSfcPresent
    end if
    writeHeightSfc = allocHeightSfc

    if (vco_trl%Vcode == 5100) then
      allocate(varNamesPsfc(2))
      varNamesPsfc = (/'P0','P0LS'/)
    else
      allocate(varNamesPsfc(1))
      varNamesPsfc = (/'P0'/)
    end if

    ! Convert all transformed variables into model variables (e.g. LVIS->VIS, LPR->PR) for original trial
    call gvt_transform(stateVectorTrial,   'AllTransformedToModel',allowOverWrite_opt=.true.)

    ! Copy the surface height from trial into stateVectorPsfc and stateVectorAnal
    if ( allocHeightSfc ) then
      HeightSfc_trial     => gsv_getHeightSfc(stateVectorTrial)
      HeightSfc_increment => gsv_getHeightSfc(stateVectorPsfc)
      HeightSfc_increment(:,:) = HeightSfc_trial(:,:)

      nullify(HeightSfc_increment)
      HeightSfc_increment => gsv_getHeightSfc(stateVectorAnal)
      HeightSfc_increment(:,:) = HeightSfc_trial(:,:)
    end if

    ! reAllocate incHighRes with the names of the model variables (e.g. VIS, PR)
    nullify(varNames)
    call gsv_varNamesList(varNames, stateVectorAnal)
    call gsv_allocate(stateVectorIncHighRes, numStep, hco_trl, vco_trl, &
                      dataKind_opt=pre_incrReal,  &
                      dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                      allocHeightSfc_opt=allocHeightSfc, hInterpolateDegree_opt=hInterpolationDegree, &
                      varNames_opt=varNames )

    ! Recompute increments and write out the files in parallel with respect to time steps
    call gsv_copy(stateVectorAnal, stateVectorIncHighRes)
    call gsv_add(stateVectorTrial, stateVectorIncHighRes, -1.0d0)

    ! figure out number of batches of time steps for writing
    numBatch = ceiling(real(stateVectorAnal%numStep) / real(mmpi_nprocs))
    call msg('inc_writeIncAndAnalHighRes', &
         'writing will be done by number of batches = '//str(numBatch))

    batch_loop: do batchIndex = 1, numBatch

      stepIndexBeg = 1 + (batchIndex - 1) * mmpi_nprocs
      stepIndexEnd = min(stateVectorAnal%numStep, stepIndexBeg + mmpi_nprocs - 1)
      call msg('inc_writeIncAndAnalHighRes', 'batchIndex = '//str(batchIndex) &
           //', stepIndexBeg = '//str(stepIndexBeg) &
           //', stepIndexEnd = '//str(stepIndexEnd))

      ! figure out which time step I will write, if any (-1 if none)
      stepIndexToWrite = -1
      do stepIndex = stepIndexBeg, stepIndexEnd
        procToWrite = nint( real(stepIndex - stepIndexBeg) * real(mmpi_nprocs) / real(stepIndexEnd - stepIndexBeg + 1) )
        if ( procToWrite == mmpi_myid ) stepIndexToWrite = stepIndex
        call msg('inc_writeIncAndAnalHighRes', ' stepIndex = '//str(stepIndex) &
             //', procToWrite = '//str(procToWrite), mpiAll_opt=.false.)
      end do

      ! determine date and allocate stateVector for storing just 1 time step, if I do writing
      if ( stepIndexToWrite /= -1 ) then

        dateStamp = stateVectorAnal%dateStampList(stepIndexToWrite)
        call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
        if(nint(deltaHours*60.0d0).lt.0) then
          write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
        else
          write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
        endif

        call gsv_allocate( stateVector_1step_r4, 1, stateVectorAnal%hco, stateVectorAnal%vco, &
                           dateStamp_opt=dateStamp, mpi_local_opt=.false., dataKind_opt=4,        &
                           allocHeightSfc_opt=allocHeightSfc, &
                           varNames_opt=varNames )
        call gsv_allocate( stateVectorPsfc_1step_r4, 1, stateVectorAnal%hco, stateVectorAnal%vco, &
                           dateStamp_opt=dateStamp, mpi_local_opt=.false., dataKind_opt=4,        &
                           varNames_opt=varNamesPsfc, allocHeightSfc_opt=allocHeightSfc )
        call gsv_zero(stateVectorPsfc_1step_r4)
      end if

      ! transpose ANALYSIS data from Tiles to Steps
      call gsv_transposeTilesToStep(stateVector_1step_r4, stateVectorAnal, stepIndexBeg)

      ! write the ANALYSIS file for one timestep on all tasks with data
      if ( stepIndexToWrite /= -1 ) then
        call msg_memUsage('inc_writeIncAndAnalHighRes')
        anlFileName = './anlm_' // trim(coffset) // 'm'
        call gio_writeToFile( stateVector_1step_r4, trim(anlFileName), etiket_anlm,  &
                              typvar_opt='A', writeHeightSfc_opt=writeHeightSfc, &
                              numBits_opt=writeNumBits, containsFullField_opt=.true. )
      end if

      if (writeHiresIncrement) then

        ! transpose INCREMENT data from Tiles to Steps
        call gsv_transposeTilesToStep(stateVector_1step_r4, stateVectorIncHighRes, stepIndexBeg)

        ! write the INCREMENT file for one timestep on all tasks with data
        if ( stepIndexToWrite /= -1 ) then
          call msg_memUsage('inc_writeIncAndAnalHighRes')
          incFileName = './rehm_' // trim(coffset) // 'm'
          call gio_writeToFile( stateVector_1step_r4, trim(incFileName), etiket_rehm,  &
                                typvar_opt='R', &
                                numBits_opt=writeNumBits, containsFullField_opt=.false. )
        end if

        if (gsv_varExist(varName='P0')) then

          ! transpose ANALYSIS PSFC AND height_SFC ONLY data from Tiles to Steps
          call gsv_transposeTilesToStep(stateVectorPsfc_1step_r4, stateVectorPsfc, stepIndexBeg)

          ! Also write analysis value of Psfc and surface height to increment file
          if ( stepIndexToWrite /= -1 ) then
            call msg_memUsage('inc_writeIncAndAnalHighRes')
            incFileName = './rehm_' // trim(coffset) // 'm'
            call gio_writeToFile( stateVectorPsfc_1step_r4, trim(incFileName), etiket_rehm,  &
                                  typvar_opt='A', writeHeightSfc_opt=writeHeightSfc, &
                                  numBits_opt=writeNumBits, containsFullField_opt=.true. )
          end if

        end if

      end if ! writeHiresIncrement

      if ( stepIndexToWrite /= -1 ) then
        call gsv_deallocate(stateVector_1step_r4)
        call gsv_deallocate(stateVectorPsfc_1step_r4)
      end if

    end do batch_loop

    call gsv_deallocate(stateVectorAnal)
    if (gsv_varExist(varName='P0')) then
      call gsv_deallocate(stateVectorPsfc)
    end if
    call gsv_deallocate(stateVectorIncHighRes)
    call gsv_deallocate(stateVectorTrial)

    call utl_tmg_stop(83)
    call utl_tmg_stop(80)

    call msg('inc_writeIncAndAnalHighRes', 'END', verb_opt=2)

  end subroutine inc_writeIncAndAnalHighRes

  !--------------------------------------------------------------------------
  ! inc_getIncrement
  !--------------------------------------------------------------------------
  subroutine inc_getIncrement(incr_cv,statevector_incr,nvadim_mpilocal)
    !
    ! :Purpose: Get true analysis increment from control vector.
    !
    implicit none

    ! Arguments:
    real(8),          intent(in)    :: incr_cv(:)
    type(struct_gsv), intent(inout) :: statevector_incr
    integer,          intent(in)    :: nvadim_mpilocal

    call utl_tmg_start(80,'--Increment')
    call utl_tmg_start(84,'----GetIncrement')

    ! compute increment from control vector (multiply by B^1/2)
    call bmat_sqrtB( incr_cv, nvadim_mpilocal, statevector_incr )

    ! Compute new diagnotics based on NAMSTATE
    if ( gsv_varExist(statevector_incr,'QR') .and. gsv_varExist(statevector_incr,'DD') ) then
       call msg('inc_getIncrement', 'User is asking for Vort-Div analysis increment')
       call gvt_transform( statevector_incr, & ! INOUT
                           'UVtoVortDiv' )     ! IN
       if ( gsv_varExist(statevector_incr,'PP') .and. gsv_varExist(statevector_incr,'CC') ) then
          call msg('inc_getIncrement', 'User is asking for Psi-Chi analysis increment')
          call gvt_transform( statevector_incr, & ! INOUT
                              'VortDivToPsiChi')  ! IN
       end if
    end if

    call utl_tmg_stop(84)
    call utl_tmg_stop(80)

  end subroutine inc_getIncrement

  !--------------------------------------------------------------------------
  ! inc_writeIncrement
  !--------------------------------------------------------------------------
  subroutine inc_writeIncrement( stateVector_incr, &
                                 ip3ForWriteToFile_opt )
    !
    ! :Purpose: Write the low-resolution analysis increments to the rebm file.
    !
    implicit none

    ! Arguments:
    type(struct_gsv),          intent(in) :: stateVector_incr
    integer,         optional, intent(in) :: ip3ForWriteToFile_opt

    ! Locals:
    integer              :: stepIndex, dateStamp
    real(8)              :: deltaHours
    character(len=4)     :: coffset
    character(len=30)    :: fileName

    call msg('inc_writeIncrement', 'START', verb_opt=2)

    call utl_tmg_start(80,'--Increment')
    call utl_tmg_start(85,'----WriteIncrement')

    ! loop over times for which increment is computed
    do stepIndex = 1, tim_nstepobsinc
      if (gsv_isAllocated(statevector_incr)) then
        dateStamp = gsv_getDateStamp(stateVector_incr,stepIndex)
        call msg('inc_writeIncrement', 'writing increment for time step: ' &
             //str(stepIndex)//', dateStamp = '//str(dateStamp), mpiAll_opt=.false.)

        ! write the increment file for this time step
        call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
        if(nint(deltaHours*60.0d0).lt.0) then
          write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
        else
          write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
        end if
        fileName = './rebm_' // trim(coffset) // 'm'
        call gio_writeToFile( stateVector_incr, fileName, etiket_rebm, scaleFactor_opt=1.0d0, &
                              ip3_opt=ip3ForWriteToFile_opt, stepIndex_opt=stepIndex, &
                              containsFullField_opt=.false. )
      end if
    end do

    call utl_tmg_stop(85)
    call utl_tmg_stop(80)

    call msg('inc_writeIncrement', 'END', verb_opt=2)


  end subroutine inc_writeIncrement

  !--------------------------------------------------------------------------
  ! inc_writeAnalysis
  !--------------------------------------------------------------------------
  subroutine inc_writeAnalysis(statevector_anal)
    ! :Purpose: To write to output standard file the analysid from statevector strucure (1Dvar case)
    !           to output the results
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in) :: statevector_anal

    ! Locals:
    integer              :: stepIndex, dateStamp
    real(8)              :: deltaHours
    character(len=4)     :: coffset
    character(len=30)    :: fileName

    call msg('inc_writeAnalysis', 'START', verb_opt=2)

    call utl_tmg_start(80,'--Increment')
    call utl_tmg_start(86,'----WriteAnalysis')

    !
    !- Set/Read values for the namelist NAMINC
    !
    call readNameList()

    ! loop over times for which increment is computed
    do stepIndex = 1, tim_nstepobsinc
      if (gsv_isAllocated(statevector_anal)) then
        dateStamp = gsv_getDateStamp(statevector_anal,stepIndex)
        call msg('inc_writeAnalysis', 'writing analysis for time step: ' &
             //str(stepIndex)//', dateStamp = '//str(dateStamp), mpiAll_opt=.false.)

        ! write the increment file for this time step
        call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
        if(nint(deltaHours*60.0d0).lt.0) then
          write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
        else
          write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
        end if
        fileName = './anlm_' // trim(coffset) // 'm'
        call gio_writeToFile( statevector_anal, fileName, etiket_anlm, scaleFactor_opt = 1.0d0, &
             ip3_opt = 0, stepIndex_opt = stepIndex, containsFullField_opt=.true. )
      end if
    end do

    call utl_tmg_stop(86)
    call utl_tmg_stop(80)

    call msg('inc_writeAnalysis', 'END', verb_opt=2)
  end subroutine inc_writeAnalysis

  !--------------------------------------------------------------------------
  ! inc_interpolateAndAdd
  !--------------------------------------------------------------------------
  subroutine inc_interpolateAndAdd(statevector_in,statevector_inout, &
                                   statevectorRef_opt, statevectorMaskLAM_opt, &
                                   scaleFactor_opt)
    !
    ! :Purpose: Interpolate the low-resolution increments to trial grid and add to
    !           get the high-resolution analysis.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), target,   intent(in)    :: statevector_in
    type(struct_gsv),           intent(inout) :: statevector_inout
    type(struct_gsv), optional, intent(in)    :: statevectorRef_opt         ! Reference statevector providing optional fields (P0, TT, HU)
    type(struct_gsv), optional, intent(in)    :: statevectorMaskLAM_opt
    real(8),          optional, intent(in)    :: scaleFactor_opt

    ! Locals:
    type(struct_gsv)          :: statevector_in_hvInterp
    type(struct_gsv)          :: statevector_in_hvtInterp
    type(struct_gsv), target  :: statevector_in_withP0LS
    type(struct_gsv), pointer :: statevector_in_ptr

    character(len=4), pointer :: varNamesToInterpolate(:), varNamesToInterpolate_withP0LS(:)

    call msg('inc_interpolateAndAdd', 'START', verb_opt=2)

    ! Error traps
    if ( .not. gsv_isAllocated(statevector_in) ) then
      call utl_abort('inc_interpolateAndAdd: gridStateVector_in not yet allocated! Aborting.')
    end if
    if ( .not. gsv_isAllocated(statevector_inout) ) then
      call utl_abort('inc_interpolateAndAdd: gridStateVector_inout not yet allocated! Aborting.')
    end if
    if ( present(statevectorRef_opt) ) then
      if ( statevector_in%numstep /= statevectorRef_opt%numstep) then
        call utl_abort('inc_interpolateAndAdd: statevector_in and statevectorRef_opt numstep inconsistent')
      end if
    end if

    nullify(varNamesToInterpolate)
    call vnl_varNamesFromExistList(varNamesToInterpolate, statevector_in%varExistlist(:))

    ! Check if P0LS needs to be added to the increment
    if (statevector_inout%vco%vcode == 5100 .and. &
         .not.gsv_varExist(statevector_in, 'P0LS')) then
      varNamesToInterpolate_withP0LS => vnl_addToVarNames(varNamesToInterpolate,'P0LS')
      call gsv_allocate(statevector_in_withP0LS, statevector_in%numstep,                       &
                        statevector_in%hco, statevector_in%vco,                                &
                        dateStamp_opt=tim_getDateStamp(),                                      &
                        mpi_local_opt=statevector_in%mpi_local, mpi_distribution_opt='Tiles',  &
                        dataKind_opt=statevector_in%dataKind,                                  &
                        allocHeightSfc_opt=statevector_in%heightSfcPresent,                    &
                        varNames_opt=varNamesToInterpolate_withP0LS,                           &
                        hInterpolateDegree_opt=statevector_in%hInterpolateDegree)
      call gsv_zero(statevector_in_withP0LS)
      call gsv_copy(statevector_in, statevector_in_withP0LS, &
                    allowVarMismatch_opt=.true.)
      statevector_in_ptr => statevector_in_withP0LS
      deallocate(varNamesToInterpolate)
      varNamesToInterpolate => varNamesToInterpolate_withP0LS
    else
      statevector_in_ptr => statevector_in
    end if
    ! Do the spatial interpolation of statevector_in onto the grid of statevector_inout
    call gsv_allocate(statevector_in_hvInterp, statevector_in%numstep,                          &
                      statevector_inout%hco, statevector_inout%vco,                             &
                      dateStamp_opt=tim_getDateStamp(),                                         &
                      mpi_local_opt=statevector_inout%mpi_local, mpi_distribution_opt='Tiles',  &
                      dataKind_opt=statevector_inout%dataKind,                                  &
                      allocHeightSfc_opt=statevector_inout%heightSfcPresent,                    &
                      varNames_opt=varNamesToInterpolate,                                       &
                      hInterpolateDegree_opt=statevector_inout%hInterpolateDegree,              &
                      hExtrapolateDegree_opt='VALUE' )

    call int_interp_gsv(statevector_in, statevector_in_hvInterp, &
                        statevectorRef_opt=statevectorRef_opt)

    ! Do the time interpolation
    call gsv_allocate(statevector_in_hvtInterp, statevector_inout%numstep,                      &
                      statevector_inout%hco, statevector_inout%vco,                             &
                      dateStamp_opt=tim_getDateStamp(),                                         &
                      mpi_local_opt=statevector_inout%mpi_local, mpi_distribution_opt='Tiles',  &
                      dataKind_opt=statevector_inout%dataKind,                                  &
                      allocHeightSfc_opt=statevector_inout%heightSfcPresent,                    &
                      varNames_opt=varNamesToInterpolate,                                       &
                      hInterpolateDegree_opt=statevector_inout%hInterpolateDegree,              &
                      hExtrapolateDegree_opt='VALUE' )
    call int_tInterp_gsv(statevector_in_hvInterp, statevector_in_hvtInterp)
    call gsv_deallocate(statevector_in_hvInterp)

    ! Masking
    if (present(statevectorMaskLAM_opt)) then
      call gsv_applyMaskLAM(statevector_in_hvtInterp,statevectorMaskLAM_opt)
    end if

    ! Do the summation
    call gsv_add(statevector_in_hvtInterp,statevector_inout,scaleFactor_opt=scaleFactor_opt)

    call gsv_deallocate(statevector_in_hvtInterp)
    deallocate(varNamesToInterpolate)

    call msg('inc_interpolateAndAdd', 'END', verb_opt=2)

  end subroutine inc_interpolateAndAdd

END MODULE increment_mod

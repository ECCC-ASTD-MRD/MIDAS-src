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

module increment_mod
  ! MODULE increment_mod (prefix='inc' category='1. High-level functionality')
  !
  ! :Purpose: To add a 4D increment to a given 4D background/reference state and
  !           to output the results
  !
  use codePrecision_mod
  use mpi_mod
  use mpivar_mod
  use timeCoord_mod
  use gridStateVector_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use humidityLimits_mod
  use utilities_mod
  use gridVariableTransforms_mod
  use BMatrix_mod
  use varNamelist_mod
  implicit none
  save
  private

  ! public procedures
  public :: inc_computeHighResAnalysis, inc_writeIncrementHighRes, inc_getIncrement, inc_writeIncrement
  public :: inc_writeAnalysis, inc_analPostProcessing

  integer, external :: get_max_rss

  ! namelist variables
  integer  :: writeNumBits
  logical  :: writeHiresIncrement
  logical  :: imposeRttovHuLimits, useAnalIncMask
  character(len=12) :: etiket_anlm, etiket_rehm, etiket_rebm
  character(len=12) :: hInterpolationDegree

CONTAINS

  !--------------------------------------------------------------------------
  ! readNameList
  !--------------------------------------------------------------------------
  subroutine readNameList
    !
    ! :Purpose: Reading NAMINC namelist by any subroutines in increment_mod module.
    !
    implicit none

    integer :: nulnam, ierr
    integer, external :: fnom, fclos
    logical, save :: nmlAlreadyRead = .false.
    NAMELIST /NAMINC/ writeHiresIncrement, etiket_rehm, etiket_anlm, &
         etiket_rebm, writeNumBits, imposeRttovHuLimits, hInterpolationDegree, &
         useAnalIncMask

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

      if ( .not. utl_isNamelistPresent('NAMINC','./flnml') ) then
        if ( mpi_myid == 0 ) then
          write(*,*) 'NAMINC is missing in the namelist. The default values will be taken.'
        end if

      else
        ! Reading the namelist
        nulnam = 0
        ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
        read(nulnam, nml=naminc, iostat=ierr)
        if ( ierr /= 0) call utl_abort('readNameList: Error reading namelist')
        ierr = fclos(nulnam)
      end if
      if ( mpi_myid == 0 ) write(*,nml=naminc)
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
    type(struct_gsv), intent(in) :: statevectorIncLowRes
    type(struct_gsv), intent(inout) :: stateVectorUpdateHighRes
    type(struct_gsv), intent(inout) :: stateVectorPsfcHighRes

    ! Locals:
    type(struct_gsv) :: statevectorPsfcLowRes
    type(struct_gsv) :: statevectorPsfcLowResTime
    type(struct_gsv) :: statevector_mask
    type(struct_gsv) :: statevectorPsfc
    type(struct_gsv) :: stateVectorHighRes

    type(struct_vco), pointer :: vco_trl => null()
    type(struct_hco), pointer :: hco_trl => null()

    integer              :: stepIndex, numStep
    integer, allocatable :: dateStampList(:)
    integer              :: get_max_rss

    real(pre_incrReal), pointer :: PsfcTrial(:,:,:,:), PsfcAnalysis(:,:,:,:)
    real(pre_incrReal), pointer :: PsfcIncrement(:,:,:,:)
    real(pre_incrReal), pointer :: PsfcIncLowResFrom3Dgsv(:,:,:,:), PsfcIncLowRes(:,:,:,:)
    real(pre_incrReal), pointer :: analIncMask(:,:,:)

    logical  :: allocHeightSfc

    write(*,*) 'inc_computeHighResAnalysis: STARTING'
    write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

    ! Set/Read values for the namelist NAMINC
    call readNameList

    ! Setup timeCoord module (date read from trial file)
    numStep = tim_nstepobsinc
    allocate(dateStampList(numStep))
    call tim_getstamplist(dateStampList,numStep,tim_getDatestamp())

    hco_trl => gsv_getHco(stateVectorUpdateHighRes)
    vco_trl => gsv_getVco(stateVectorUpdateHighRes)
    if (vco_trl%Vcode == 0 .or. .not. gsv_varExist(varName='P0')) then
      allocHeightSfc = .false.
    else
      allocHeightSfc = stateVectorUpdateHighRes%heightSfcPresent
    end if

    ! Read the analysis mask (in LAM mode only) - N.B. different from land/sea mask!!!
    if (.not. hco_trl%global .and. useAnalIncMask) then
      call gsv_getMaskLAM(statevector_mask, hco_trl, vco_trl, hInterpolationDegree)
    end if

    ! Get the increment of Psfc
    if ( gsv_varExist(varName='P0') ) then
      call gsv_allocate( statevectorPsfc, numStep, hco_trl, vco_trl, &
                         dataKind_opt=pre_incrReal, &
                         dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true.,  &
                         varNames_opt=(/'P0'/), allocHeightSfc_opt=allocHeightSfc, &
                         hInterpolateDegree_opt=hInterpolationDegree )

      if( mpi_myid == 0 ) write(*,*) ''
      if( mpi_myid == 0 ) write(*,*) 'inc_computeHighResAnalysis: horizontal interpolation of the Psfc increment'

      ! Extract Psfc inc at low resolution
      call gsv_allocate( statevectorPsfcLowRes, numStep,  &
                         statevectorIncLowRes%hco, vco_trl,  &
                         dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true.,  &
                         dataKind_opt=pre_incrReal, &
                         varNames_opt=(/'P0'/) )
      call gsv_getField(statevectorPsfcLowRes,PsfcIncLowRes,'P0')
      call gsv_getField(statevectorIncLowRes,PsfcIncLowResFrom3Dgsv,'P0')
      PsfcIncLowRes(:,:,1,:) = PsfcIncLowResFrom3Dgsv(:,:,1,:)

      ! Spatial interpolation of Psfc analysis increment
      call gsv_interpolate(statevectorPsfcLowRes,statevectorPsfc)
      call gsv_deallocate(statevectorPsfcLowRes)

      ! Compute analysis Psfc to use for interpolation of increment
      if( mpi_myid == 0 ) write(*,*) 'inc_computeHighResAnalysis: Computing Psfc analysis to use for interpolation of increment'
      call gsv_allocate(statevectorPsfcLowResTime, tim_nstepobsinc, hco_trl, vco_trl,  &
                        dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true.,  &
                        dataKind_opt=pre_incrReal, &
                        varNames_opt=(/'P0'/))
      call gsv_copy( stateVectorUpdateHighRes, statevectorPsfcLowResTime, &
                     allowVarMismatch_opt=.true., allowTimeMismatch_opt=.true. )

      call gsv_getField(statevectorPsfcLowResTime,PsfcTrial,'P0')
      call gsv_getField(stateVectorPsfc,PsfcIncrement,'P0')
      call gsv_getField(stateVectorPsfc,PsfcAnalysis,'P0')

      if (.not. hco_trl%global .and. useAnalIncMask) then
        call gsv_getField(statevector_mask,analIncMask)
        do stepIndex = 1, stateVectorPsfc%numStep
          PsfcAnalysis(:,:,1,stepIndex) = PsfcTrial(:,:,1,stepIndex) + &
                                          PsfcIncrement(:,:,1,stepIndex) * analIncMask(:,:,1)
        end do
      else
        PsfcAnalysis(:,:,1,:) = PsfcTrial(:,:,1,:) + PsfcIncrement(:,:,1,:)
      end if

      ! Time interpolation to get high-res Psfc analysis increment
      if( mpi_myid == 0 ) write(*,*) 'inc_computeHighResAnalysis: Time interpolation to get high-res Psfc analysis increment'
      if ( .not. gsv_isAllocated(stateVectorPsfcHighRes) ) then
        call gsv_allocate( stateVectorPsfcHighRes, tim_nstepobs, hco_trl, vco_trl, &
                           dataKind_opt=pre_incrReal, &
                           dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true.,  &
                           varNames_opt=(/'P0'/), allocHeightSfc_opt=allocHeightSfc, &
                           hInterpolateDegree_opt=hInterpolationDegree)
      else
        call gsv_zero( stateVectorPsfcHighRes )
      end if
      call gsv_tInterpolate(statevectorPsfc, stateVectorPsfcHighRes)
    end if

    ! Compute the analysis
    if( mpi_myid == 0 ) write(*,*) ''
    if( mpi_myid == 0 ) write(*,*) 'inc_computeHighResAnalysis: compute the analysis'
    call tmg_start(181,'INC_COMPUTEANL')

    ! Interpolate low-res increments to high-res and add to the initial state
    call gsv_allocate( stateVectorHighRes, tim_nstepobs, hco_trl, vco_trl, &
                       dataKind_opt=stateVectorUpdateHighRes%dataKind, &
                       dateStamp_opt=tim_getDateStamp(), &
                       mpi_local_opt=stateVectorUpdateHighRes%mpi_local, &
                       allocHeightSfc_opt=stateVectorUpdateHighRes%heightSfcPresent, &
                       hInterpolateDegree_opt=hInterpolationDegree, &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_copy( stateVectorUpdateHighRes,stateVectorHighRes, &
                   allowVarMismatch_opt=.true.)

    if (.not. hco_trl%global .and. useAnalIncMask) then
      if (gsv_varExist(varName='P0')) then
        call inc_interpolateAndAdd(statevectorIncLowRes, stateVectorHighRes, &
                                   PsfcReference_opt=PsfcAnalysis(:,:,1,:), statevectorMaskLAM_opt=statevector_mask)
      else
        call inc_interpolateAndAdd(statevectorIncLowRes, stateVectorHighRes, &
                                   statevectorMaskLAM_opt=statevector_mask)
      end if
    else
      if (gsv_varExist(varName='P0')) then
        call inc_interpolateAndAdd(statevectorIncLowRes, stateVectorHighRes, &
                                   PsfcReference_opt=PsfcAnalysis(:,:,1,:))
      else
        call inc_interpolateAndAdd(statevectorIncLowRes, stateVectorHighRes)
      end if
    end if
    call tmg_stop(181)

    call gsv_copy( stateVectorHighRes, stateVectorUpdateHighRes, &
                   allowVarMismatch_opt=.true.)
    call gsv_deallocate(stateVectorHighRes)

    if ( gsv_isAllocated(statevectorPsfc) ) call gsv_deallocate(statevectorPsfc)
    if ( gsv_isAllocated(statevector_mask) ) call gsv_deallocate(statevector_mask)

    write(*,*) 'inc_computeHighResAnalysis: END'

  end subroutine inc_computeHighResAnalysis

  !--------------------------------------------------------------------------
  ! inc_analPostProcessing
  !--------------------------------------------------------------------------
  subroutine inc_analPostProcessing (stateVectorPsfcHighRes, stateVectorUpdateHighRes, &  ! IN
                                     stateVectorTrial, stateVectorPsfc, stateVectorAnal ) ! OUT
    !
    ! :Purpose: Post processing of the high resolution analysis including degrading 
    !           temporal resolution, variable transform, and humidity clipping.
    !
    implicit none 

    ! Arguments:
    type(struct_gsv), intent(in)  :: stateVectorPsfcHighRes
    type(struct_gsv), intent(in)  :: stateVectorUpdateHighRes
    type(struct_gsv), intent(out) :: stateVectorTrial
    type(struct_gsv), intent(out) :: stateVectorPsfc
    type(struct_gsv), intent(out) :: stateVectorAnal

    ! Locals:
    type(struct_vco), pointer :: vco_trl => null()
    type(struct_hco), pointer :: hco_trl => null()

    real(pre_incrReal), pointer :: oceanIce_ptr(:,:,:,:)

    logical  :: allocHeightSfc

    write(*,*) 'inc_analPostProcessing: STARTING'
    write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

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
    call gsv_readTrials(stateVectorTrial)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if (vco_trl%Vcode == 0 .or. .not. gsv_varExist(varName='P0')) then
      allocHeightSfc = .false.
    else
      allocHeightSfc = stateVectorTrial%heightSfcPresent
    end if

    ! Degrade timesteps stateVector high-res Psfc, and high-res analysis
    if (gsv_varExist(varName='P0')) then
      call gsv_allocate(stateVectorPsfc, tim_nstepobsinc, hco_trl, vco_trl, &
                        dataKind_opt = pre_incrReal, &
                        dateStamp_opt = tim_getDateStamp(), mpi_local_opt = .true.,  &
                        varNames_opt = (/'P0'/), allocHeightSfc_opt = allocHeightSfc, &
                        hInterpolateDegree_opt = hInterpolationDegree)
      call gsv_copy(stateVectorPsfcHighRes, stateVectorPsfc, allowTimeMismatch_opt = .true.)
      call gsv_deallocate(stateVectorPsfcHighRes)
    end if

    call gsv_allocate(stateVectorAnal, tim_nstepobsinc, hco_trl, vco_trl,  &
                      dataKind_opt = pre_incrReal, &
                      dateStamp_opt = tim_getDateStamp(), mpi_local_opt = .true., &
                      allocHeightSfc_opt = allocHeightSfc, hInterpolateDegree_opt = 'LINEAR', &
                      allocHeight_opt = .false., allocPressure_opt = .false.)
    call gsv_copy(stateVectorUpdateHighRes, stateVectorAnal, &
                  allowVarMismatch_opt = .true., allowTimeMismatch_opt = .true.)

    ! Start the variable transformations
    if (gsv_varExist(stateVectorAnal, 'GL')) then
      ! Impose limits [0,1] on sea ice concentration analysis
      call gsv_getField(stateVectorAnal, oceanIce_ptr, 'GL')
      oceanIce_ptr(:,:,:,:) = min(oceanIce_ptr(:,:,:,:), 1.0d0)
      oceanIce_ptr(:,:,:,:) = max(oceanIce_ptr(:,:,:,:), 0.0d0)
      if (gsv_varExist(stateVectorAnal, 'LG' ) ) then
        ! Compute the continuous sea ice concentration field (LG)
        call gvt_transform( stateVectorAnal, 'oceanIceContinuous', stateVectorRef_opt = stateVectorTrial, varName_opt = 'LG' )
      end if
    end if

    ! Start the variable transformations
    if( gsv_varExist(stateVectorAnal,'TM') ) then
      call gsv_getField( stateVectorAnal, oceanIce_ptr, 'TM' )
      ! Compute the continuous SST field (TM)
      call gvt_transform( stateVectorAnal, 'oceanIceContinuous', stateVectorRef_opt = stateVectorTrial, varName_opt = 'TM' )
    end if

    ! Convert all transformed variables into model variables (e.g. LVIS->VIS, LPR->PR) for analysis
    call gvt_transform(stateVectorAnal, 'AllTransformedToModel', allowOverWrite_opt = .true.)

    ! Impose limits on humidity analysis
    call tmg_start(182, 'INC_QLIMITS')
    write(*,*) 'inc_analPostProcessing: calling qlim_saturationLimit'
    call qlim_saturationLimit(stateVectorAnal)
    if (imposeRttovHuLimits) call qlim_rttovLimit(stateVectorAnal)
    call tmg_stop(182)

    if (gsv_varKindExist('CH')) then
      ! Apply boundaries to analysis of CH kind variables as needed.
      write(*,*) 'inc_analPostProcessing: applying minimum values to analysis for variables of CH kind'
      call gvt_transform(stateVectorAnal, 'CH_bounds')
    end if

    write(*,*) 'inc_analPostProcessing: END'

  end subroutine inc_analPostProcessing

  !--------------------------------------------------------------------------
  ! inc_writeIncrementHighRes
  !--------------------------------------------------------------------------
  subroutine inc_writeIncrementHighRes( stateVectorTrial, stateVectorPsfc, &
                                        stateVectorAnal )
    !
    ! :Purpose: Write the high-resolution analysis increments to the rehm file.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout), target :: stateVectorTrial
    type(struct_gsv), intent(in) :: stateVectorPsfc
    type(struct_gsv), intent(in) :: stateVectorAnal

    ! Locals:
    type(struct_gsv) :: stateVectorIncHighRes
    type(struct_gsv) :: stateVector_1step_r4, stateVectorPsfc_1step_r4

    type(struct_vco), pointer :: vco_trl => null()
    type(struct_hco), pointer :: hco_trl => null()

    integer              :: stepIndex, stepIndexBeg, stepIndexEnd, stepIndexToWrite, numStep
    integer              :: dateStamp, numBatch, batchIndex, procToWrite
    integer, allocatable :: dateStampList(:)
    integer              :: get_max_rss

    character(len=256)  :: incFileName, anlFileName
    character(len=4)    :: coffset
    character(len=4), pointer :: varNames(:)

    real(8)             :: deltaHours
    real(8), pointer    :: HeightSfc_increment(:,:), HeightSfc_trial(:,:)

    logical  :: allocHeightSfc, writeHeightSfc

    write(*,*) 'inc_writeIncrementHighRes: STARTING'
    write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

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
    call tmg_start(183,'INC_WRITEANLMREHM')
    call gsv_copy(stateVectorAnal, stateVectorIncHighRes)
    call gsv_add(stateVectorTrial, stateVectorIncHighRes, -1.0d0)

    ! figure out number of batches of time steps for writing
    numBatch = ceiling(real(stateVectorAnal%numStep) / real(mpi_nprocs))
    write(*,*) 'inc_writeIncrementHighRes: writing will be done by number of batches = ', numBatch

    batch_loop: do batchIndex = 1, numBatch

      stepIndexBeg = 1 + (batchIndex - 1) * mpi_nprocs
      stepIndexEnd = min(stateVectorAnal%numStep, stepIndexBeg + mpi_nprocs - 1)
      write(*,*) 'inc_writeIncrementHighRes: batchIndex, stepIndexBeg/End = ', batchIndex, stepIndexBeg, stepIndexEnd

      ! figure out which time step I will write, if any (-1 if none)
      stepIndexToWrite = -1
      do stepIndex = stepIndexBeg, stepIndexEnd
        procToWrite = nint( real(stepIndex - stepIndexBeg) * real(mpi_nprocs) / real(stepIndexEnd - stepIndexBeg + 1) )
        if ( procToWrite == mpi_myid ) stepIndexToWrite = stepIndex
        if ( mpi_myid == 0 ) write(*,*) 'inc_writeIncrementHighRes: stepIndex, procToWrite = ', stepIndex, procToWrite
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
                           varNames_opt=(/'P0'/), allocHeightSfc_opt=allocHeightSfc )
      end if

      ! transpose ANALYSIS data from Tiles to Steps
      call gsv_transposeTilesToStep(stateVector_1step_r4, stateVectorAnal, stepIndexBeg)

      ! write the ANALYSIS file for one timestep on all tasks with data
      if ( stepIndexToWrite /= -1 ) then
        write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
        anlFileName = './anlm_' // trim(coffset) // 'm'
        call gsv_writeToFile( stateVector_1step_r4, trim(anlFileName), etiket_anlm,  &
                              typvar_opt='A', writeHeightSfc_opt=writeHeightSfc, &
                              numBits_opt=writeNumBits, containsFullField_opt=.true. )
      end if

      if (writeHiresIncrement) then

        ! transpose INCREMENT data from Tiles to Steps
        call gsv_transposeTilesToStep(stateVector_1step_r4, stateVectorIncHighRes, stepIndexBeg)

        ! write the INCREMENT file for one timestep on all tasks with data
        if ( stepIndexToWrite /= -1 ) then
          write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
          incFileName = './rehm_' // trim(coffset) // 'm'
          call gsv_writeToFile( stateVector_1step_r4, trim(incFileName), etiket_rehm,  &
                                typvar_opt='R', &
                                numBits_opt=writeNumBits, containsFullField_opt=.false. )
        end if

        if (gsv_varExist(varName='P0')) then

          ! transpose ANALYSIS PSFC AND height_SFC ONLY data from Tiles to Steps
          call gsv_transposeTilesToStep(stateVectorPsfc_1step_r4, stateVectorPsfc, stepIndexBeg)

          ! Also write analysis value of Psfc and surface height to increment file
          if ( stepIndexToWrite /= -1 ) then
            write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
            incFileName = './rehm_' // trim(coffset) // 'm'
            call gsv_writeToFile( stateVectorPsfc_1step_r4, trim(incFileName), etiket_rehm,  &
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
    call tmg_stop(183)

    call gsv_deallocate(stateVectorAnal)
    if (gsv_varExist(varName='P0')) then
      call gsv_deallocate(stateVectorPsfc)
    end if
    call gsv_deallocate(stateVectorIncHighRes)
    call gsv_deallocate(stateVectorTrial)

    write(*,*) 'inc_writeIncrementHighRes: END'

  end subroutine inc_writeIncrementHighRes

  !--------------------------------------------------------------------------
  ! inc_getIncrement
  !--------------------------------------------------------------------------
  subroutine inc_getIncrement(incr_cv,statevector_incr,nvadim_mpilocal)
    !
    ! :Purpose: Get true analysis increment from control vector.
    !
    implicit none

    ! arguments
    real(8)          :: incr_cv(:)
    type(struct_gsv) :: statevector_incr
    integer          :: nvadim_mpilocal

    ! compute increment from control vector (multiply by B^1/2)
    call bmat_sqrtB( incr_cv, nvadim_mpilocal, statevector_incr )

    ! Compute new diagnotics based on NAMSTATE
    if ( gsv_varExist(statevector_incr,'QR') .and. gsv_varExist(statevector_incr,'DD') ) then
       write(*,*)
       write(*,*) 'User is asking for Vort-Div analysis increment'
       call gvt_transform( statevector_incr, & ! INOUT
                           'UVtoVortDiv' )     ! IN
       if ( gsv_varExist(statevector_incr,'PP') .and. gsv_varExist(statevector_incr,'CC') ) then
          write(*,*)
          write(*,*) 'User is asking for Psi-Chi analysis increment'
          call gvt_transform( statevector_incr, & ! INOUT
                              'VortDivToPsiChi')  ! IN
       end if
    end if

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
    type(struct_gsv)          :: stateVector_incr
    integer,         optional :: ip3ForWriteToFile_opt

    ! Locals:
    integer              :: stepIndex, dateStamp
    real(8)              :: deltaHours
    character(len=4)     :: coffset
    character(len=30)    :: fileName

    if ( mpi_myid == 0 ) write(*,*) 'inc_writeIncrement: STARTING'

    ! loop over times for which increment is computed
    do stepIndex = 1, tim_nstepobsinc
      if (gsv_isAllocated(statevector_incr)) then
        dateStamp = gsv_getDateStamp(stateVector_incr,stepIndex)
        if ( mpi_myid == 0 ) write(*,*) 'inc_writeIncrement: writing increment for time step: ', &
                                         stepIndex, dateStamp

        ! write the increment file for this time step
        call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
        if(nint(deltaHours*60.0d0).lt.0) then
          write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
        else
          write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
        end if
        fileName = './rebm_' // trim(coffset) // 'm'
        call gsv_writeToFile( stateVector_incr, fileName, etiket_rebm, scaleFactor_opt=1.0d0, &
                              ip3_opt=ip3ForWriteToFile_opt, stepIndex_opt=stepIndex, &
                              containsFullField_opt=.false. )
      end if
    end do

    if ( mpi_myid == 0 ) write(*,*) 'inc_writeIncrement: END'

  end subroutine inc_writeIncrement

  !--------------------------------------------------------------------------
  ! inc_writeAnalysis
  !--------------------------------------------------------------------------
  subroutine inc_writeAnalysis(statevector_anal)
    ! :Purpose: To write to output standard file the analysid from statevector strucure (1Dvar case) 
    !           to output the results
    !
    implicit none

    ! arguments
    type(struct_gsv)     :: statevector_anal
    ! locals
    integer              :: stepIndex, dateStamp
    real(8)              :: deltaHours
    character(len=4)     :: coffset
    character(len=30)    :: fileName

    if(mpi_myid == 0) write(*,*) 'inc_writeAnalysis: STARTING'

    !
    !- Set/Read values for the namelist NAMINC
    !
    call readNameList()

    ! loop over times for which increment is computed
    do stepIndex = 1, tim_nstepobsinc
      if (gsv_isAllocated(statevector_anal)) then
        dateStamp = gsv_getDateStamp(statevector_anal,stepIndex)
        if(mpi_myid == 0) write(*,*) 'inc_writeAnalysis: writing analysis for time step: ',stepIndex, dateStamp

        ! write the increment file for this time step
        call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
        if(nint(deltaHours*60.0d0).lt.0) then
          write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
        else
          write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
        end if
        fileName = './anlm_' // trim(coffset) // 'm'
        call gsv_writeToFile( statevector_anal, fileName, etiket_anlm, scaleFactor_opt = 1.0d0, &
             ip3_opt = 0, stepIndex_opt = stepIndex, containsFullField_opt=.true. )
      end if
    end do

  end subroutine inc_writeAnalysis

  !--------------------------------------------------------------------------
  ! inc_interpolateAndAdd
  !--------------------------------------------------------------------------
  subroutine inc_interpolateAndAdd(statevector_in,statevector_inout,scaleFactor_opt, &
                                   PsfcReference_opt, statevectorMaskLAM_opt)
    !
    ! :Purpose: Interpolate the low-resolution increments to trial grid and add to 
    !           get the high-resolution analysis.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)    :: statevector_in
    type(struct_gsv), intent(inout) :: statevector_inout
    real(8), intent(in), optional   :: scaleFactor_opt
    real(pre_incrReal), intent(in), optional :: PsfcReference_opt(:,:,:)
    type(struct_gsv), intent(in), optional   :: statevectorMaskLAM_opt

    ! Locals:
    type(struct_gsv) :: statevector_in_hvInterp
    type(struct_gsv) :: statevector_in_hvtInterp

    real(4), allocatable        :: PsfcReference_r4(:,:,:)
    real(8), allocatable        :: PsfcReference_r8(:,:,:)

    character(len=4), pointer :: varNamesToInterpolate(:)

    write(*,*) 'inc_interpolateAndAdd: STARTING'

    ! Error traps
    if ( .not. gsv_isAllocated(statevector_in) ) then
      call utl_abort('inc_interpolateAndAdd: gridStateVector_in not yet allocated! Aborting.')
    end if
    if ( .not. gsv_isAllocated(statevector_inout) ) then
      call utl_abort('inc_interpolateAndAdd: gridStateVector_inout not yet allocated! Aborting.')
    end if
    if ( present(PsfcReference_opt) ) then
      if ( statevector_in%numstep /= size(PsfcReference_opt,3) ) then
        call utl_abort('inc_interpolateAndAdd: statevector_in%numstep /= numStep of Psfc input')
      end if
    end if

    nullify(varNamesToInterpolate)
    call vnl_varNamesFromExistList(varNamesToInterpolate, statevector_in%varExistlist(:))

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

    ! Need to copy input PsfcReference to both real4 and real8 to isolate gsv from pre_incrReal
    if (gsv_getDataKind(statevector_in)==4) then
      allocate(PsfcReference_r4(statevector_in_hvInterp%myLonBeg:statevector_in_hvInterp%myLonEnd, &
                                statevector_in_hvInterp%myLatBeg:statevector_in_hvInterp%myLatEnd, &
                                statevector_in_hvInterp%numStep))
      if (present(PsfcReference_opt)) then
        PsfcReference_r4(:,:,:) = PsfcReference_opt(:,:,:)
      end if
      call gsv_interpolate(statevector_in,statevector_in_hvInterp,PsfcReference_r4_opt=PsfcReference_r4)
      deallocate(PsfcReference_r4)
    else
      allocate(PsfcReference_r8(statevector_in_hvInterp%myLonBeg:statevector_in_hvInterp%myLonEnd, &
                                statevector_in_hvInterp%myLatBeg:statevector_in_hvInterp%myLatEnd, &
                                statevector_in_hvInterp%numStep))
      if (present(PsfcReference_opt)) then
        PsfcReference_r8(:,:,:) = PsfcReference_opt(:,:,:)
      end if
      call gsv_interpolate(statevector_in,statevector_in_hvInterp,PsfcReference_opt=PsfcReference_r8)
      deallocate(PsfcReference_r8)
    end if

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
    call gsv_tInterpolate(statevector_in_hvInterp, statevector_in_hvtInterp)
    call gsv_deallocate(statevector_in_hvInterp)

    ! Masking
    if (present(statevectorMaskLAM_opt)) then
      call gsv_applyMaskLAM(statevector_in_hvtInterp,statevectorMaskLAM_opt)
    end if

    ! Do the summation
    call gsv_add(statevector_in_hvtInterp,statevector_inout,scaleFactor_opt=scaleFactor_opt)

    call gsv_deallocate(statevector_in_hvtInterp)
    deallocate(varNamesToInterpolate)

    write(*,*) 'inc_interpolateAndAdd: END'

  end subroutine inc_interpolateAndAdd

END MODULE increment_mod

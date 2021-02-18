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
  public :: inc_computeAndWriteAnalysis, inc_getIncrement, inc_writeIncrement

  ! namelist variables
  integer  :: writeNumBits
  logical  :: writeHiresIncrement
  logical  :: imposeRttovHuLimits, useAnalIncMask
  character(len=12) :: etiket_anlm, etiket_rehm, etiket_rebm
  character(len=12) :: hInterpolationDegree

CONTAINS

  !--------------------------------------------------------------------------
  ! readNameList - called by the other public subroutines in this module
  !--------------------------------------------------------------------------
  subroutine readNameList
    implicit none

    integer :: nulnam, ierr
    integer, external :: fnom, fclos
    NAMELIST /NAMINC/ writeHiresIncrement, etiket_rehm, etiket_anlm, &
         etiket_rebm, writeNumBits, imposeRttovHuLimits, hInterpolationDegree, &
         useAnalIncMask


    !- Setting default values
    writeHiresIncrement = .true.
    imposeRttovHuLimits = .true.
    useAnalIncMask      = .false.
    etiket_rehm = 'INCREMENT'
    etiket_rebm = 'INCREMENT'
    etiket_anlm = 'ANALYSIS'
    writeNumBits = 16
    hInterpolationDegree = 'LINEAR'

    !- Read the namelist
    nulnam = 0
    ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    read(nulnam, nml=naminc, iostat=ierr)
    if ( ierr /= 0) call utl_abort('increment_mod readNameList: Error reading namelist')
    if ( mpi_myid == 0 ) write(*,nml=naminc)
    ierr = fclos(nulnam)

  end subroutine readNameList

  !--------------------------------------------------------------------------
  ! inc_computeAndWriteAnalysis
  !--------------------------------------------------------------------------
  subroutine inc_computeAndWriteAnalysis(statevector_incLowRes_opt, stateVectorTrial_opt)
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in), optional :: statevector_incLowRes_opt
    type(struct_gsv), intent(in), target, optional :: statevectorTrial_opt

    ! Locals:
    type(struct_gsv) :: statevector_incHighRes
    type(struct_gsv), pointer :: statevector_trial
    type(struct_gsv) :: statevector_analysis
    type(struct_gsv) :: statevector_PsfcLowRes, statevector_Psfc
    type(struct_gsv) :: statevector_1step_r4, statevector_Psfc_1step_r4
    type(struct_gsv) :: statevector_mask

    type(struct_vco), pointer :: vco_trl => null()
    type(struct_vco), pointer :: vco_inc => null()
    type(struct_hco), pointer :: hco_trl => null()

    integer              :: stepIndex, stepIndexBeg, stepIndexEnd, stepIndexToWrite, numStep
    integer              :: dateStamp, numBatch, batchIndex, procToWrite
    integer, allocatable :: dateStampList(:)
    integer              :: get_max_rss, latIndex, kIndex, lonIndex

    character(len=256)  :: trialFileName, incFileName, anlFileName
    character(len=4)    :: coffset
    character(len=4), pointer :: anlVar(:)
    character(len=4), pointer :: varNames(:)

    real(8)             :: deltaHours
    real(pre_incrReal), pointer :: PsfcTrial(:,:,:,:), PsfcAnalysis(:,:,:,:), analInc(:,:,:,:)
    real(pre_incrReal), pointer :: PsfcIncrement(:,:,:,:)
    real(pre_incrReal), pointer :: PsfcIncLowResFrom3Dgsv(:,:,:,:), PsfcIncLowRes(:,:,:,:)
    real(8), pointer            :: HeightSfc_increment(:,:), HeightSfc_trial(:,:)
    real(pre_incrReal), pointer :: GL_ptr(:,:,:,:)
    real(pre_incrReal), pointer :: analIncMask(:,:,:)
    real(8), allocatable        :: PsfcAnalysis_r8(:,:)

    logical  :: allocHeightSfc, writeHeightSfc, useIncLevelsOnly

    !
    !- Set/Read values for the namelist NAMINC
    !
    call readNameList()

    write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

    ! Setup timeCoord module (date read from trial file)
    numStep = tim_nstepobsinc
    allocate(dateStampList(numStep))
    call tim_getstamplist(dateStampList,numStep,tim_getDatestamp())

    !- Do we need to read all the vertical levels from the trial fields?
    if (present(statevector_incLowRes_opt)) then
      vco_inc => statevector_incLowRes_opt%vco
    else
      call difdatr(datestamplist(1),tim_getDatestamp(),deltaHours)
      if(nint(deltaHours*60.0d0).lt.0) then
        write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
      else
        write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
      endif
      incFileName = './rebm_' // trim(coffset) // 'm'

      if ( mpi_myid == 0 ) then
        call vco_setupFromFile(vco_inc, incFileName)
      end if
      call vco_mpiBcast(vco_inc)
    end if

    if (present(stateVectorTrial_opt)) then

      !- If stateVectorTrial is supplied, then just use it

      hco_trl => gsv_getHco(stateVectorTrial_opt)
      vco_trl => gsv_getVco(stateVectorTrial_opt)
      if (vco_trl%Vcode == 0 .or. .not. gsv_varExist(varName='P0')) then
        allocHeightSfc = .false.
      else
        allocHeightSfc = stateVectorTrial_opt%heightSfcPresent
      end if
      
      !- In some cases we need to just extract a subset of levels from the trials
      useIncLevelsOnly = vco_subsetOrNot(vco_inc, vco_trl)
      if ( useIncLevelsOnly ) then
        write(*,*) 'inc_computeAndWriteAnalysis: extract only the increment levels from the trials'
        allocate(statevector_trial)
        call gsv_allocate(statevector_trial, tim_nstepobsinc, hco_trl, vco_inc,   &
                          dataKind_opt=pre_incrReal, &
                          dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                          allocHeightSfc_opt=allocHeightSfc, hInterpolateDegree_opt=hInterpolationDegree, &
                          allocHeight_opt=.false., allocPressure_opt=.false.)
        call gsv_interpolate(stateVectorTrial_opt, statevector_trial)        
        vco_trl => gsv_getVco(statevector_trial)
      else
        write(*,*) 'inc_computeAndWriteAnalysis: use the supplied trials directly'
        statevector_trial => stateVectorTrial_opt
      end if

    else

      !- Initialize the trial state grid
      if (mpi_myid == 0) write(*,*) ''
      if (mpi_myid == 0) write(*,*) 'inc_computeAndWriteAnalysis: Set hco parameters for trial grid'
      trialFileName = './trlm_01'

      nullify(anlVar)
      call gsv_varNamesList(anlVar)
      call hco_setupFromFile( hco_trl, trim(trialFileName), ' ', varName_opt=anlVar(1))

      if ( mpi_myid == 0 ) then
        call vco_setupFromFile( vco_trl, trim(trialFileName) )
      end if
      call vco_mpiBcast(vco_trl)

      if (vco_trl%Vcode == 0 .or. .not. gsv_varExist(varName='P0')) then
        allocHeightSfc = .false.
      else
        allocHeightSfc = .true.
      end if

      useIncLevelsOnly = vco_subsetOrNot(vco_inc, vco_trl)
      if ( useIncLevelsOnly ) then
        ! Read only the increment levels
        write(*,*)
        write(*,*) 'inc_computeAndWriteAnalysis: only the increment levels will be read in the trials'
        call  vco_deallocate(vco_trl)
        vco_trl => vco_inc
      else
        ! Read them all
        write(*,*)
        write(*,*) 'inc_computeAndWriteAnalysis: all the vertical levels will be read in the trials'
        if (.not. present(statevector_incLowRes_opt)) then
          call vco_deallocate(vco_inc)
        end if
      end if

      write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

      !
      !- Read trial files
      !
      if(mpi_myid == 0) write(*,*) ''
      if(mpi_myid == 0) write(*,*) 'inc_computeAndWriteAnalysis: reading background state for all time steps'
      allocate(statevector_trial)
      call gsv_allocate(statevector_trial, tim_nstepobsinc, hco_trl, vco_trl,   &
                        dataKind_opt=pre_incrReal, &
                        dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                        allocHeightSfc_opt=allocHeightSfc, hInterpolateDegree_opt=hInterpolationDegree, &
                        allocHeight_opt=.false., allocPressure_opt=.false.)

      call tmg_start(180,'INC_READTRIALS')
      call gsv_readTrials(statevector_trial)
      call tmg_stop(180)

    end if ! present(stateVectorTrial_opt)

    !
    !- Read the analysis mask (in LAM mode only) - N.B. different from land/sea mask!!!
    !
    if (.not. hco_trl%global .and. useAnalIncMask) then
      call gsv_allocate(statevector_mask, 1, hco_trl, vco_trl, dateStamp_opt=-1, &
                        dataKind_opt=pre_incrReal, &
                        mpi_local_opt=.true., varNames_opt=(/'MSKC'/),           &
                        hInterpolateDegree_opt=hInterpolationDegree)
      call gsv_readFromFile(statevector_mask, './analinc_mask', ' ', ' ', unitConversion_opt=.false., &
                            vcoFileIn_opt=vco_trl)
    end if

    !
    !- Get the increment of Psfc
    !
    if (gsv_varExist(varName='P0')) then
      call gsv_allocate(statevector_Psfc, numStep, hco_trl, vco_trl, &
                        dataKind_opt=pre_incrReal, &
                        dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true.,  &
                        varNames_opt=(/'P0'/), allocHeightSfc_opt=allocHeightSfc, &
                        hInterpolateDegree_opt=hInterpolationDegree)

      if (present(statevector_incLowRes_opt)) then
        if(mpi_myid == 0) write(*,*) ''
        if(mpi_myid == 0) write(*,*) 'inc_computeAndWriteAnalysis: horiz interpolation of the Psfc increment'

        ! Extract Psfc inc at low resolution
        call gsv_allocate(statevector_PsfcLowRes, numStep,  &
                          statevector_incLowRes_opt%hco, vco_trl,  &
                          dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true.,  &
                          dataKind_opt=pre_incrReal, &
                          varNames_opt=(/'P0'/))
        call gsv_getField(statevector_PsfcLowRes,PsfcIncLowRes,'P0')
        call gsv_getField(statevector_incLowRes_opt,PsfcIncLowResFrom3Dgsv,'P0')
        PsfcIncLowRes(:,:,1,:) = PsfcIncLowResFrom3Dgsv(:,:,1,:)

        ! Interpolate
        call gsv_interpolate(statevector_PsfcLowRes,statevector_Psfc)
        call gsv_deallocate(statevector_PsfcLowRes)

      else
        !- Read from file
        do stepIndex = 1, numStep
          dateStamp = datestamplist(stepIndex)
          if(mpi_myid == 0) write(*,*) ''
          if(mpi_myid == 0) write(*,*) 'inc_computeAndWriteAnalysis: reading Psfc increment for time step: ',stepIndex, dateStamp

          call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
          if(nint(deltaHours*60.0d0).lt.0) then
            write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
          else
            write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
          endif
          incFileName = './rebm_' // trim(coffset) // 'm'

          call gsv_readFromFile( statevector_Psfc, trim(incFileName), ' ', ' ', stepIndex,  &
                                 containsFullField_opt=.false. )
        end do

      end if

      !
      !- Compute analysis Psfc to use for interpolation of increment
      !
      call gsv_getField(statevector_trial,PsfcTrial,'P0')
      call gsv_getField(statevector_Psfc,PsfcIncrement,'P0')
      call gsv_getField(statevector_Psfc,PsfcAnalysis,'P0')

      if (.not. hco_trl%global .and. useAnalIncMask) then
        call gsv_getField(statevector_mask,analIncMask)
        do stepIndex = 1, statevector_trial%numStep
          PsfcAnalysis(:,:,1,stepIndex) = PsfcTrial(:,:,1,stepIndex) + PsfcIncrement(:,:,1,stepIndex)*analIncMask(:,:,1)
        end do
      else
        PsfcAnalysis(:,:,1,:) = PsfcTrial(:,:,1,:) + PsfcIncrement(:,:,1,:)
      end if

      !
      !- Copy the surface height from trial into statevector_Psfc
      !
      HeightSfc_increment => gsv_getHeightSfc(statevector_Psfc)
      HeightSfc_trial     => gsv_getHeightSfc(statevector_trial)
      HeightSfc_increment(:,:) = HeightSfc_trial(:,:)
    end if
    writeHeightSfc = allocHeightSfc

    !
    !- Compute the analysis
    !
    if(mpi_myid == 0) write(*,*) ''
    if(mpi_myid == 0) write(*,*) 'inc_computeAndWriteAnalysis: compute the analysis'
    call tmg_start(181,'INC_COMPUTEANL')

    call gsv_allocate(statevector_incHighRes, numStep, hco_trl, vco_trl, &
                      dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                      dataKind_opt=pre_incrReal, allocHeightSfc_opt=allocHeightSfc, &
                      hInterpolateDegree_opt=hInterpolationDegree, &
                      allocHeight_opt=.false., allocPressure_opt=.false.)
    call gsv_allocate(statevector_analysis, numStep, hco_trl, vco_trl, &
                      dataKind_opt=pre_incrReal, &
                      dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                      allocHeightSfc_opt=allocHeightSfc, hInterpolateDegree_opt=hInterpolationDegree, &
                      allocHeight_opt=.false., allocPressure_opt=.false.)
    call gsv_copy(statevector_trial, statevector_analysis)

    if (present(statevector_incLowRes_opt)) then
      !- Interpolate and add the input increments
      if (.not. hco_trl%global .and. useAnalIncMask) then
        if (gsv_varExist(varName='P0')) then
          call inc_interpolateAndAdd(statevector_incLowRes_opt, statevector_analysis, &
                                     PsfcReference_opt=PsfcAnalysis(:,:,1,:), mask2d_opt=statevector_mask)
        else
          call inc_interpolateAndAdd(statevector_incLowRes_opt, statevector_analysis, &
                                     mask2d_opt=statevector_mask)
        end if
      else
        if (gsv_varExist(varName='P0')) then
          call inc_interpolateAndAdd(statevector_incLowRes_opt, statevector_analysis,&
                                     PsfcReference_opt=PsfcAnalysis(:,:,1,:))
        else
          call inc_interpolateAndAdd(statevector_incLowRes_opt, statevector_analysis)
        end if
      end if
    else
      !- Read the increments from files
      do stepIndex = 1, numStep
        dateStamp = datestamplist(stepIndex)
        if(mpi_myid == 0) write(*,*) ''
        if(mpi_myid == 0) write(*,*) 'inc_computeAndWriteAnalysis: reading increment for time step: ',stepIndex, dateStamp

        call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
        if(nint(deltaHours*60.0d0).lt.0) then
          write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
        else
          write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
        endif
        incFileName = './rebm_' // trim(coffset) // 'm'

        write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
        if (gsv_varExist(varName='P0')) then
          allocate(PsfcAnalysis_r8(statevector_Psfc%myLonBeg:statevector_Psfc%myLonEnd, &
                                   statevector_Psfc%myLatBeg:statevector_Psfc%myLatEnd))
          PsfcAnalysis_r8(:,:) = PsfcAnalysis(:,:,1,stepIndex)
          call gsv_readFromFile( statevector_incHighRes, trim(incFileName), ' ', ' ', stepIndex,  &
                                 PsfcReference_opt=PsfcAnalysis_r8,  &
                                 containsFullField_opt=.false. )
          deallocate(PsfcAnalysis_r8)
        else
          call gsv_readFromFile( statevector_incHighRes, trim(incFileName), ' ', ' ', stepIndex, &
                                 containsFullField_opt=.false. )
        end if
      end do

      if (.not. hco_trl%global .and. useAnalIncMask) then
        nullify(analIncMask)
        call gsv_getField(statevector_incHighRes,analInc)
        call gsv_getField(statevector_mask,analIncMask)
        do stepIndex = 1, statevector_incHighRes%numStep
          !$OMP PARALLEL DO PRIVATE (latIndex,kIndex,lonIndex)
          do kIndex = 1, statevector_incHighRes%nk
            do latIndex = statevector_incHighRes%myLatBeg, statevector_incHighRes%myLatEnd
              do lonIndex = statevector_incHighRes%myLonBeg, statevector_incHighRes%myLonEnd
                analInc(lonIndex,latIndex,kIndex,stepIndex) = analInc(lonIndex,latIndex,kIndex,stepIndex) * &
                                                              analIncMask(lonIndex,latIndex,1)
              end do
            end do
          end do
        end do
      end if

      call gsv_add(statevector_incHighRes, statevector_analysis)

    end if
    call tmg_stop(181)

    if( gsv_varExist(stateVector_analysis,'GL') ) then

      !
      !- Impose limits [0,1] on sea ice concentration analysis
      !
      call gsv_getField(stateVector_analysis,GL_ptr,'GL')
      GL_ptr(:,:,:,:) = min(GL_ptr(:,:,:,:), 1.0d0)
      GL_ptr(:,:,:,:) = max(GL_ptr(:,:,:,:), 0.0d0)

      if( gsv_varExist(stateVector_analysis,'LG') ) then

        !
        !- Compute the continuous sea ice concentration field (LG)
        !
        call gvt_transform(stateVector_analysis,'GLtoLG',stateVectorRef_opt=stateVector_trial)

      end if

    end if

    !
    !- Convert all transformed variables into model variables (e.g. LVIS->VIS, LPR->PR)
    !
    call gvt_transform(stateVector_analysis,'AllTransformedToModel',allowOverWrite_opt=.true.)
    call gvt_transform(stateVector_trial,   'AllTransformedToModel',allowOverWrite_opt=.true.)

    !- reAllocate incHighRes with the names of the model variables (e.g. VIS, PR)
    nullify(varNames)
    call gsv_varNamesList(varNames, stateVector_analysis)
    call gsv_deallocate( stateVector_incHighRes )
    call gsv_allocate(statevector_incHighRes, numStep, hco_trl, vco_trl, &
                      dataKind_opt=pre_incrReal,  &
                      dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                      allocHeightSfc_opt=allocHeightSfc, hInterpolateDegree_opt=hInterpolationDegree, &
                      varNames_opt=varNames )

    !
    !- Impose limits on humidity analysis
    !
    call tmg_start(182,'INC_QLIMITS')
    write(*,*) 'inc_computeAndWriteAnalysis: calling qlim_saturationLimit'
    call qlim_saturationLimit(statevector_analysis)
    if( imposeRttovHuLimits ) call qlim_rttovLimit(statevector_analysis)

    if (gsv_varKindExist('CH')) then
      !
      !- Apply boundaries to analysis of CH kind variables as needed.
      !
      write(*,*) 'inc_computeAndWriteAnalysis: applying minimum values to analysis for variables of CH kind'
      call gvt_transform(stateVector_analysis,'CH_bounds')
    end if
    
    !
    !- Recompute increments
    !
    call gsv_copy(statevector_analysis, statevector_incHighRes)
    call gsv_add(statevector_trial, statevector_incHighRes, -1.0d0)

    call tmg_stop(182)

    !
    !- Write out the files in parallel with respect to time steps
    !
    call tmg_start(183,'INC_WRITEANLMREHM')

    ! figure out number of batches of time steps for writing
    numBatch = ceiling(real(stateVector_analysis%numStep) / real(mpi_nprocs))
    write(*,*) 'inc_computeAndWriteAnalysis: writing will be done by number of batches = ', numBatch

    batch_loop: do batchIndex = 1, numBatch

      stepIndexBeg = 1 + (batchIndex - 1) * mpi_nprocs
      stepIndexEnd = min(stateVector_analysis%numStep, stepIndexBeg + mpi_nprocs - 1)
      write(*,*) 'inc_computeAndWriteAnalysis: batchIndex, stepIndexBeg/End = ', batchIndex, stepIndexBeg, stepIndexEnd

      ! figure out which time step I will write, if any (-1 if none)
      stepIndexToWrite = -1
      do stepIndex = stepIndexBeg, stepIndexEnd
        procToWrite = nint( real(stepIndex - stepIndexBeg) * real(mpi_nprocs) / real(stepIndexEnd - stepIndexBeg + 1) )
        if ( procToWrite == mpi_myid ) stepIndexToWrite = stepIndex
        if ( mpi_myid == 0 ) write(*,*) 'inc_computeAndWriteAnalysis: stepIndex, procToWrite = ', stepIndex, procToWrite
      end do

      ! determine date and allocate stateVector for storing just 1 time step, if I do writing
      if ( stepIndexToWrite /= -1 ) then

        dateStamp = stateVector_analysis%dateStampList(stepIndexToWrite)
        call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
        if(nint(deltaHours*60.0d0).lt.0) then
          write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
        else
          write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
        endif

        call gsv_allocate( stateVector_1step_r4, 1, stateVector_analysis%hco, stateVector_analysis%vco, &
                           dateStamp_opt=dateStamp, mpi_local_opt=.false., dataKind_opt=4,        &
                           allocHeightSfc_opt=allocHeightSfc, &
                           varNames_opt=varNames )
        call gsv_allocate( stateVector_Psfc_1step_r4, 1, stateVector_analysis%hco, stateVector_analysis%vco, &
                           dateStamp_opt=dateStamp, mpi_local_opt=.false., dataKind_opt=4,        &
                           varNames_opt=(/'P0'/), allocHeightSfc_opt=allocHeightSfc )
      end if

      ! transpose ANALYSIS data from Tiles to Steps
      call gsv_transposeTilesToStep(stateVector_1step_r4, stateVector_analysis, stepIndexBeg)

      ! write the ANALYSIS file for one timestep on all tasks with data
      if ( stepIndexToWrite /= -1 ) then
        write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
        anlFileName = './anlm_' // trim(coffset) // 'm'
        call gsv_writeToFile( statevector_1step_r4, trim(anlFileName), etiket_anlm,  &
                              typvar_opt='A', writeHeightSfc_opt=writeHeightSfc, &
                              numBits_opt=writeNumBits, containsFullField_opt=.true. )
      end if

      if (writeHiresIncrement) then

        ! transpose INCREMENT data from Tiles to Steps
        call gsv_transposeTilesToStep(stateVector_1step_r4, stateVector_incHighRes, stepIndexBeg)

        ! write the INCREMENT file for one timestep on all tasks with data
        if ( stepIndexToWrite /= -1 ) then
          write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
          incFileName = './rehm_' // trim(coffset) // 'm'
          call gsv_writeToFile( statevector_1step_r4, trim(incFileName), etiket_rehm,  &
                                typvar_opt='R', &
                                numBits_opt=writeNumBits, containsFullField_opt=.false. )
        end if

        if (gsv_varExist(varName='P0')) then

          ! transpose ANALYSIS PSFC AND height_SFC ONLY data from Tiles to Steps
          call gsv_transposeTilesToStep(stateVector_Psfc_1step_r4, stateVector_Psfc, stepIndexBeg)

          ! Also write analysis value of Psfc and surface height to increment file
          if ( stepIndexToWrite /= -1 ) then
            write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
            incFileName = './rehm_' // trim(coffset) // 'm'
            call gsv_writeToFile( statevector_Psfc_1step_r4, trim(incFileName), etiket_rehm,  &
                                  typvar_opt='A', writeHeightSfc_opt=writeHeightSfc, &
                                  numBits_opt=writeNumBits, containsFullField_opt=.true. )
          end if

        end if

      end if ! writeHiresIncrement

      if ( stepIndexToWrite /= -1 ) then
        call gsv_deallocate(stateVector_1step_r4)
        call gsv_deallocate(stateVector_Psfc_1step_r4)
      end if

    end do batch_loop
    call tmg_stop(183)

    call gsv_deallocate(statevector_analysis)
    if (gsv_varExist(varName='P0')) then
      call gsv_deallocate(statevector_Psfc)
    end if
    call gsv_deallocate(statevector_incHighRes)
    call gsv_deallocate(statevector_trial)

  end subroutine inc_computeAndWriteAnalysis

  !--------------------------------------------------------------------------
  ! inc_getIncrement
  !--------------------------------------------------------------------------
  subroutine inc_getIncrement(incr_cv,statevector_incr,nvadim_mpilocal)

    implicit none

    ! arguments
    real(8) :: incr_cv(:)
    type(struct_gsv) :: statevector_incr
    integer :: nvadim_mpilocal

    ! compute increment from control vector (multiply by B^1/2)
    call bmat_sqrtB(incr_cv, nvadim_mpilocal, statevector_incr)

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
  subroutine inc_writeIncrement(statevector_incr)

    implicit none

    ! arguments
    type(struct_gsv)     :: statevector_incr
    ! locals
    integer              :: stepIndex, dateStamp
    real(8)              :: deltaHours
    character(len=4)     :: coffset
    character(len=30)    :: fileName

    if(mpi_myid == 0) write(*,*) 'inc_writeIncrement: STARTING'

    !
    !- Set/Read values for the namelist NAMINC
    !
    call readNameList()

    ! loop over times for which increment is computed
    do stepIndex = 1, tim_nstepobsinc

      dateStamp = gsv_getDateStamp(statevector_incr,stepIndex)
      if(mpi_myid == 0) write(*,*) 'inc_writeIncrement: writing increment for time step: ',stepIndex, dateStamp

      ! write the increment file for this time step
      call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
      if(nint(deltaHours*60.0d0).lt.0) then
        write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
      else
        write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
      endif
      fileName = './rebm_' // trim(coffset) // 'm'
      call gsv_writeToFile( statevector_incr, fileName, etiket_rebm, 1.0d0, 0,  &
                            stepIndex, containsFullField_opt=.false. )

    enddo

  end subroutine inc_writeIncrement

  !--------------------------------------------------------------------------
  ! inc_interpolateAndAdd
  !--------------------------------------------------------------------------
  subroutine inc_interpolateAndAdd(statevector_in,statevector_inout,scaleFactor_opt, &
                                   PsfcReference_opt, mask2d_opt)
    implicit none
    type(struct_gsv)  :: statevector_in, statevector_inout

    real(8), optional :: scaleFactor_opt
    real(pre_incrReal), optional :: PsfcReference_opt(:,:,:)
    type(struct_gsv),   optional :: mask2d_opt

    type(struct_gsv) :: statevector_in_hvInterp

    real(pre_incrReal), pointer :: increment(:,:,:,:)
    real(pre_incrReal), pointer :: analIncMask(:,:,:)
    real(4), allocatable        :: PsfcReference_r4(:,:,:)
    real(8), allocatable        :: PsfcReference_r8(:,:,:)

    integer :: latIndex,kIndex,lonIndex, stepIndex

    character(len=4), pointer :: varNamesToInterpolate(:)

    !
    !- Error traps
    !
    if (.not.statevector_in%allocated) then
      call utl_abort('inc_interpolateAndAdd: gridStateVector_in not yet allocated! Aborting.')
    end if
    if (.not.statevector_inout%allocated) then
      call utl_abort('inc_interpolateAndAdd: gridStateVector_inout not yet allocated! Aborting.')
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

    !
    !- Masking
    !
    if (present(mask2d_opt)) then
      call gsv_getField(statevector_in_hvInterp,increment)
      call gsv_getField(mask2d_opt,analIncMask)
      do stepIndex = 1, statevector_in_hvInterp%numStep
        !$OMP PARALLEL DO PRIVATE (latIndex,kIndex,lonIndex)
        do kIndex = 1, statevector_in_hvInterp%nk
          do latIndex =  statevector_in_hvInterp%myLatBeg,  statevector_in_hvInterp%myLatEnd
            do lonIndex =  statevector_in_hvInterp%myLonBeg,  statevector_in_hvInterp%myLonEnd
              increment(lonIndex,latIndex,kIndex,stepIndex) =      &
                   increment(lonIndex,latIndex,kIndex,stepIndex) * &
                   analIncMask(lonIndex,latIndex,1)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do
    end if

    !
    !- Do the summation
    !
    call gsv_add(statevector_in_hvInterp,statevector_inout,scaleFactor_opt=scaleFactor_opt)

    call gsv_deallocate(statevector_in_hvInterp)
    deallocate(varNamesToInterpolate)

  end subroutine inc_interpolateAndAdd

END MODULE increment_mod

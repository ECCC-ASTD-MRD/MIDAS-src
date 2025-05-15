module var1DIdealize_mod
    ! MODULE var1DIdealize_mod (prefix='var1Di' category='4. Data Object transformations')
    !
    ! :Purpose: contains all 1Dvar-related methods.
    !
    use midasMpi_mod
    use codeprecision_mod
    use mathphysconstants_mod
    use columnData_mod
    use columnVariableTransforms_mod
    use controlVector_mod
    use gridStatevector_mod
    use horizontalCoord_mod
    use obsSpaceData_mod
    use timeCoord_mod
    use utilities_mod
    use verticalCoord_mod
    use randomNumber_mod
    use bMatrix1Dvar_mod
    use innovation_mod
    use var1D_mod
    use gridStateVectorFileIO_mod
    use interpolation_mod
    use varNameList_mod
    use increment_mod
    use rMatrix_mod
    use humidityLimits_mod
    use obsoperators_mod
    use obsErrors_mod
    use tovs_mod
    use surfaceEmissivity_mod

    implicit none
    save
    private

    ! Public procedures
    public :: var1Di_simulateBackgroundState, var1Di_simulateObservation
    public :: var1Di_estSigmaBObsSpace

  contains

  !--------------------------------------------------------------------------
  ! var1Di_simulateBackgroundState
  !--------------------------------------------------------------------------
  subroutine var1Di_simulateBackgroundState(columnTruthOnTrlLev, columnSimTrlOnTrlLev, &
                                            obsSpaceData, vco_anl, seed, inflateEmissErr)
    !
    !:Purpose: Simulate the background state by adding a perturbation from the reference state (Truth)
    !
    implicit none

    ! Arguments:
    type(struct_columnData), target, intent(inout) :: columnTruthOnTrlLev   ! The true state column
    type(struct_columnData), target, intent(out)   :: columnSimTrlOnTrlLev  ! Simulated background column
    type(struct_obs),                intent(in)    :: obsSpaceData          ! ObsSpacedata object
    type(struct_vco), pointer,       intent(in)    :: vco_anl               ! Analysis vertical coordinate structure
    integer,                         intent(in)    :: seed                  ! Seed to random number generator
    real(8),                         intent(in)    :: inflateEmissErr       ! Emissvity error inflation scale factor

    ! Locals:
    type(struct_columnData)         :: columnPertOnAnLev
    type(struct_columnData)         :: columnTruthOnAnlLev
    type(struct_columnData)         :: columnPertOnTrlLev
    real(8), allocatable            :: controlVector(:)
    integer                         :: cvIndex
    type(struct_gsv)                :: stateVectorPertOnAnLevTruth
    type(struct_gsv)                :: stateVectorPertOnTrlLevTruth
    type(struct_gsv)                :: stateVectorTrlOnTrlLevTruth
    type(struct_gsv)                :: stateVectorTrlOnTrlLevSim
    type(struct_gsv)                :: stateVectorTrlOnAnlLevTruth
    character(len=50)               :: prefixFileName
    logical                         :: containsFullField
    integer                         :: randomSeed

    allocate(controlVector(cvm_nvadim))
    ! Generate perturbation sampling following gaussian distribution with zero mean and one std

    if (seed == 0) then
      randomSeed = var1Di_randomSeed()
    else
      randomSeed = seed + mmpi_myid
    end if

    call rng_setup(abs(randomSeed))

    do cvIndex = 1, cvm_nvadim
      controlVector(cvIndex) = rng_gaussian()
    end do

    ! Compute (B^1/2)*Pert (column)
    call col_setVco(columnPertOnAnLev, vco_anl)
    call col_allocate(columnPertOnAnLev, col_getNumCol(columnTruthOnTrlLev), setToZero_opt=.true.)

    if (col_varExist(columnPertOnAnLev, 'EMMW') .and. inflateEmissErr /= MPC_missingValue_R8) then
      call bmat1D_sqrtB(controlVector, cvm_nvadim, columnPertOnAnLev, obsSpaceData, &
                        inflateEmissErr_opt = inflateEmissErr)
    else
      call bmat1D_sqrtB(controlVector, cvm_nvadim, columnPertOnAnLev, obsSpaceData)
    end if

    call col_setVco(columnPertOnTrlLev, col_getVco(columnTruthOnTrlLev))
    call col_allocate(columnPertOnTrlLev, col_getNumCol(columnTruthOnTrlLev), setToZero_opt=.true.)

    ! Interpolate (B^1/2)*Pert from analysis to trial level
    call var1Di_vInterpPertAnLev2TrlLev(columnPertOnAnLev, columnPertOnTrlLev, columnTruthOnTrlLev)

    call col_setVco(columnSimTrlOnTrlLev, col_getVco(columnTruthOnTrlLev))
    call col_allocate(columnSimTrlOnTrlLev, col_getNumCol(columnTruthOnTrlLev), &
                      setToZero_opt=.true.)
    call col_copy(columnTruthOnTrlLev, columnSimTrlOnTrlLev)

    ! Add the truth and (B^1/2)*Pert columns
    call col_add(columnPertOnTrlLev, columnSimTrlOnTrlLev)

    ! Compute the pressure levels
    call cvt_transform(columnSimTrlOnTrlLev, 'ZandP_nl')

    ! Restrict the simulated humidity background within physically reasonable values.
    call qlim_rttovLimit(columnSimTrlOnTrlLev)

    ! Restict simulate background surface emissivity within physical values
    if (col_varExist(columnSimTrlOnTrlLev, 'EMMW')) then
      call sse_emissivityRttovLimits(columnSimTrlOnTrlLev, columnRef_opt = columnTruthOnTrlLev)
    end if

    ! Interpolate the truth from trial to analysis increment levels
    call col_setVco(columnTruthOnAnlLev, vco_anl)
    call col_allocate(columnTruthOnAnlLev, col_getNumCol(columnTruthOnTrlLev), setToZero_opt=.true.)
    call inn_setupColumnsOnAnlIncLev(columnTruthOnTrlLev, columnTruthOnAnlLev)

    ! Write trial into standard files
    prefixFileName = 'SimTrialOnTrlLev'
    containsFullField = .true.
    call var1d_transferColumnToYGrid(stateVectorTrlOnTrlLevSim, obsSpaceData, columnSimTrlOnTrlLev, bmat1D_includeAnlVar)
    call var1Di_writeSimTrial(stateVectorTrlOnTrlLevSim, prefixFileName, 'ANALYSIS', containsFullField)

    prefixFileName = 'TruthOnTrlLev'
    containsFullField = .true.
    call var1d_transferColumnToYGrid(stateVectorTrlOnTrlLevTruth, obsSpaceData, columnTruthOnTrlLev, bmat1D_includeAnlVar)
    call var1Di_writeSimTrial(stateVectorTrlOnTrlLevTruth, prefixFileName, 'ANALYSIS', containsFullField)

    prefixFileName = 'TruthOnAnlLev'
    containsFullField = .true.
    call var1d_transferColumnToYGrid(stateVectorTrlOnAnlLevTruth, obsSpaceData, columnTruthOnAnlLev, bmat1D_includeAnlVar)
    call var1Di_writeSimTrial(stateVectorTrlOnAnlLevTruth, prefixFileName, 'ANALYSIS', containsFullField)

    prefixFileName = 'PertOnTrlLev'
    containsFullField = .false.
    call var1d_transferColumnToYGrid(stateVectorPertOnTrlLevTruth, obsSpaceData, columnPertOnTrlLev, bmat1D_includeAnlVar)
    call var1Di_writeSimTrial(stateVectorPertOnTrlLevTruth, prefixFileName, 'INCREMENT', containsFullField)

    prefixFileName = 'PertOnAnlLev'
    containsFullField = .false.
    call var1d_transferColumnToYGrid(stateVectorPertOnAnLevTruth, obsSpaceData, columnPertOnAnLev, bmat1D_includeAnlVar)
    call var1Di_writeSimTrial(stateVectorPertOnAnLevTruth, prefixFileName, 'INCREMENT', containsFullField)

    if (mmpi_myId ==0) then
      call gsv_deallocate(stateVectorPertOnTrlLevTruth)
      call gsv_deallocate(stateVectorPertOnAnLevTruth)
      call gsv_deallocate(stateVectorTrlOnAnlLevTruth)
    end if

    call col_deallocate(columnPertOnTrlLev)
    call col_deallocate(columnPertOnAnLev)
    call col_deallocate(columnTruthOnAnlLev)
    deallocate(controlVector)

  end subroutine var1Di_simulateBackgroundState

  !--------------------------------------------------------------------------
  ! var1Di_vInterpPertAnLev2TrlLev
  !--------------------------------------------------------------------------
  subroutine var1Di_vInterpPertAnLev2TrlLev(columnAnlLev, columnTrlLev, columnPresRef)
    !
    ! :Purpose: Vertically Interpolate the generated perturbation from analysis to trial level.
    !           An reference column on trial level is required to compute the pressure level, which
    !           is used for the vertical interpolation
    !
    implicit none

    ! arguments
    type(struct_columnData), intent(in)     :: columnAnlLev  ! Column data in analysis level
    type(struct_columnData), intent(inout)  :: columnTrlLev  ! Column data in trial level
    type(struct_columnData), intent(in)     :: columnPresRef ! Column data where sfc pressure variables
                                                             ! will be used for vertical interpolation

    ! locals:
    integer                    :: numColumns
    integer                    :: columnIndex, varIndex
    real(8), allocatable       :: pSfcRef(:,:)
    real(8), pointer           :: columnAnlLev_ptr(:), columnTrlLev_ptr(:)

    write(*,*) 'var1Di_vInterpPertAnLev2TrlLev: Starting'

    ! Check the column size
    if (.not. (col_getNumCol(columnAnlLev) == col_getNumCol(columnTrlLev) .and.    &
        col_getNumCol(columnAnlLev) == col_getNumCol(columnPresRef))) then
      write(*,*) 'Column size columnAnlLev, columnTrlLev and columnPresRef', col_getNumCol(columnAnlLev), &
                  col_getNumCol(columnTrlLev), col_getNumCol(columnPresRef)
      call utl_abort('var1Di_vInterpPertAnLev2TrlLev: The columnAnlLev, columnTrlLev and columnPresRef ' // &
                      'do not have equal number of columns')
    end if

    numColumns = col_getNumCol(columnAnlLev)
    write(*,*) 'var1Di_vInterpPertAnLev2TrlLev: Column size', numColumns

    ! Extract the surface pressure from the columnPresRef
    allocate(pSfcRef(1,numColumns))
    do columnIndex = 1, numColumns
      pSfcRef(1, columnIndex) = col_getElem(columnPresRef, 1, columnIndex, 'P0')
    end do

    ! Vertical Interpolation
    do varIndex = 1, vnl_numvarmax3D
      ! Check if varName is an analysis variable
      if (.not. varneed(vnl_varNameList3D(varIndex))) cycle
      if (.not. col_varExist(columnAnlLev, vnl_varNameList3D(varIndex)) ) cycle
      call int_vInterp_col(columnAnlLev, columnTrlLev, vnl_varNameList3D(varIndex), sfcPressureRef_opt=pSfcRef)
    end do

    ! Copy 2D surface variables
    do varIndex = 1, vnl_numvarmax2D
      if (.not. varneed(vnl_varNameList2D(varIndex))) cycle
      if (.not. col_varExist(columnAnlLev, vnl_varNameList2D(varIndex))) cycle
      if (col_getNumCol(columnAnlLev) > 0) then
        do columnIndex = 1, col_getNumCol(columnAnlLev)
          columnTrlLev_ptr => col_getColumn(columnTrlLev , columnIndex, vnl_varNameList2D(varIndex))
          columnAnlLev_ptr => col_getColumn(columnAnlLev, columnIndex, vnl_varNameList2D(varIndex))
          columnTrlLev_ptr(:) = columnAnlLev_ptr(:)
        end do
      end if
    end do

    ! Copy other variables
    do varIndex = 1, vnl_numvarmaxOther
      if (.not. varneed(vnl_varNameListOther(varIndex))) cycle
      if (.not. col_varExist(columnAnlLev, vnl_varNameListOther(varIndex))) cycle
      if (col_getNumCol(columnAnlLev) <= 0) cycle

      do columnIndex = 1, col_getNumCol(columnAnlLev)
        columnTrlLev_ptr => col_getColumn(columnTrlLev, columnIndex, vnl_varNameListOther(varIndex))
        columnAnlLev_ptr => col_getColumn(columnAnlLev, columnIndex, vnl_varNameListOther(varIndex))
        columnTrlLev_ptr(:) = columnAnlLev_ptr(:)
      end do
    end do

    deallocate(pSfcRef)
    write(*,*) 'var1Di_vInterpPertAnLev2TrlLev: Finished'
    contains

    logical function varneed(varName)
      implicit none
      ! Arguements:
      character(len=*) :: varName ! Variable Name

      ! Locals:
      integer          :: varIndex2

      varneed=.false.
      do varIndex2=1,bmat1D_numIncludeAnlVar
        if (trim(varName) == trim(bmat1D_includeAnlVar(varIndex2))) then
          varneed=.true.
       end if
      end do

    end function varneed
  end subroutine var1Di_vInterpPertAnLev2TrlLev

  !--------------------------------------------------------------------------
  ! var1Di_writeSimTrial
  !--------------------------------------------------------------------------
  subroutine var1Di_writeSimTrial(statevectorSim, prefixFileName, etiket, containsFullField)
    !
    ! :Purpose: Write the simulate background state from statevector strucure (1Dvar case)
    !           into output standard file
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)     :: statevectorSim     ! Statevector to be written in file
    character(len=*), intent(in)     :: prefixFileName     ! Prefix of the filename
    character(len=*), intent(in)     :: etiket             ! Etiket of the filename
    logical,          intent(in)     :: containsFullField  ! Logical for full field values or Perturbation/Increments

    ! Locals:
    integer              :: stepIndex, dateStamp
    real(8)              :: deltaHours
    character(len=4)     :: coffset
    character(len=100)   :: fileName

    if(mmpi_myid == 0) write(*,*) 'var1Di_writeSimTrial: STARTING'

    ! loop over times for which increment is computed
    do stepIndex = 1, tim_nstepobsinc
      if (gsv_isAllocated(statevectorSim)) then
        dateStamp = gsv_getDateStamp(statevectorSim,stepIndex)
        if (mmpi_myid == 0) write(*,*) 'var1Di_writeSimTrial: writing for time step: ',stepIndex, dateStamp

        ! write the increment file for this time step
        call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
        if (nint(deltaHours*60.0d0).lt.0) then
          write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
        else
          write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
        end if

        fileName = './'//trim(prefixFileName)//'_' // trim(coffset) // 'm'
        call gio_writeToFile( statevectorSim, fileName, trim(etiket), scaleFactor_opt = 1.0d0, &
                              ip3_opt = 0, stepIndex_opt = stepIndex, containsFullField_opt=containsFullField, &
                              numBits_opt=32 )
      end if
    end do

    if(mmpi_myid == 0) write(*,*) 'var1Di_writeSimTrial: Finished'
  end subroutine var1Di_writeSimTrial

  !--------------------------------------------------------------------------
  ! var1Di_simulateObservation
  !--------------------------------------------------------------------------
  subroutine var1Di_simulateObservation(columnTruthOnTrlLev, obsSpaceData, datestamp, simObsSeed, simEmissSeed, useSimObsErr)
    !
    !:Purpose: Simulate the observation (only TOVS obs) by adding a perturbation from the reference data
    !          Additional changes are needed to generalize for all observations (not just TOVS obs)
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)    :: columnTruthOnTrlLev ! True column state
    type(struct_obs),        intent(inout) :: obsSpaceData        ! ObsSpacedata object
    integer,                 intent(in)    :: datestamp           ! Date stamp
    integer,                 intent(in)    :: simObsSeed          ! Seed to random number to simulate observation
    integer,                 intent(in)    :: simEmissSeed        ! Seed to random number to simulate emissivity
    logical,                 intent(in)    :: useSimObsErr        ! Simulate Observation Error Covariance

    ! Locals:
    logical              :: bgckMode, beSilent
    integer              :: randomSeed
    integer              :: headerIndex, bodyIndex, obsIndex
    integer              :: bodyIndexBeg, bodyIndexEnd, idatyp, count, channelNumber, channelIndex
    real(8), allocatable :: pert(:), obsPert(:), list_OER(:)
    integer, allocatable :: list_chanNumber(:), list_bodyIndex(:), list_chanIndex(:)

    beSilent = .false.
    bgckMode = .false.

    write(*,*) 'var1Di_simulateObservation: Starting'

    ! Compute the Truth the observation space
    write(*,*) 'var1Di_simulateObservation: Computing the truth in Obs Space'

    if (.not. allocated(tvs_emissivity) .and. obs_columnActive_RB(obsSpaceData, OBS_SEM)) then
      call tvs_allocateEmissivity(tvs_maxChannelNumber)
    end if

    ! Prepare atmospheric profiles for all tovs observation points for use in rttov
    call tvs_fillProfiles(columnTruthOnTrlLev, obsSpaceData, datestamp, 'nl', beSilent)

    ! Compute radiance
    call tvs_rttov(obsSpaceData, bgckMode, beSilent)

     ! loop over all header indices of the 'TO' family
    call obs_set_current_header_list(obsSpaceData,'TO')

    ! Store the true state (Observation Space) into ObsSpaceData
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER

      ! process only radiance data to be assimilated?
      idatyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if (.not. tvs_isIdBurpTovs(idatyp)) then
        write(*,*) 'var1Di_simulateObservation: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
        cycle HEADER
      end if

      bodyIndexBeg = obs_headElem_i(obsspacedata, OBS_RLN, headerIndex)
      bodyIndexEnd = obs_headElem_i(obsspacedata, OBS_NLV, headerIndex) + bodyIndexBeg - 1

      do bodyIndex = bodyIndexBeg, bodyIndexEnd
        if (obs_bodyElem_i(obsspacedata, OBS_ASS, bodyIndex) == obs_assimilated) then
          call tvs_getChannelNumIndexFromPPP(obsSpaceData, headerIndex, bodyIndex, &
                                                channelNumber, channelIndex)
          call obs_bodySet_r(obsSpaceData, OBS_TRUO, bodyIndex, tvs_radiance(headerIndex)%bt(channelIndex))

          if (allocated(tvs_emissivity) .and. obs_columnActive_RB(obsSpaceData, OBS_SEM)) then
            call obs_bodySet_r(obsSpaceData, OBS_SEM, bodyIndex, tvs_emissivity(channelIndex, headerIndex))
          end if
        end if
      end do
    end do HEADER

     ! loop over all header indices of the 'TO' family
    call obs_set_current_header_list(obsSpaceData,'TO')

    if (simObsSeed == 0) then
      randomSeed = var1Di_randomSeed()
    else
      randomSeed = simObsSeed + mmpi_myid
    end if

    call rng_setup(abs(randomSeed))

    HEADER2: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER2

      ! process only radiance data to be assimilated?
      idatyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if (.not. tvs_isIdBurpTovs(idatyp)) then
        write(*,*) 'var1Di_simulateObservation: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
        cycle HEADER2
      end if

      bodyIndexBeg = obs_headElem_i(obsspacedata, OBS_RLN, headerIndex)
      bodyIndexEnd = obs_headElem_i(obsspacedata, OBS_NLV, headerIndex) + bodyIndexBeg - 1

      if (tvs_isIdBurpTovs(idatyp)) then

        allocate(pert(tvs_maxChannelNumber))
        allocate(obsPert(tvs_maxChannelNumber))
        allocate(list_OER(tvs_maxChannelNumber))
        allocate(list_chanNumber(tvs_maxChannelNumber))
        allocate(list_chanIndex(tvs_maxChannelNumber))
        allocate(list_bodyIndex(tvs_maxChannelNumber))

        ! Read the Sigma O from ObsSpaceData
        count = 0
        do bodyIndex = bodyIndexBeg, bodyIndexEnd
          if (obs_bodyElem_i(obsspacedata, OBS_ASS, bodyIndex) == obs_assimilated) then
            call tvs_getChannelNumIndexFromPPP(obsSpaceData, headerIndex, bodyIndex, &
                                                channelNumber, channelIndex)
            count = count + 1
            list_bodyIndex(count) = bodyIndex
            list_chanNumber(count) = channelNumber
            list_chanIndex(count) = channelIndex
            list_OER(count) = obs_bodyElem_r(obsspacedata, OBS_OER, bodyIndex)
            pert(count) = rng_gaussian()
          end if
        end do

        if (count > 0) then
          ! Compute Observation Perturbation
          call rmat_Rsqrt(tvs_lsensor(headerIndex), count, pert(1:count), obsPert(1:count), list_chanNumber(1:count),&
                          list_OER(1:count))

          ! Update the obs value in ObsSpacedata
          do obsIndex = 1, count
            call obs_bodySet_r(obsSpaceData, OBS_VAR, list_bodyIndex(obsIndex), tvs_radiance(headerIndex)%bt(list_chanIndex(obsIndex)) + obsPert(obsIndex))
          end do
        end if

        deallocate(pert)
        deallocate(obsPert)
        deallocate(list_OER)
        deallocate(list_chanNumber)
        deallocate(list_chanIndex)
        deallocate(list_bodyIndex)
      end if
    end do HEADER2


    if (useSimObsErr) then
      if (col_varExist(columnTruthOnTrlLev, 'EMMW')) then
        ! Compute R matrix based on the difference between simulated obseration and truth
        ! that is based on the true state and simulated emissivity
        ! (H(x_true, emiss_true)+ obs_error) - H(x_true, emiss_true)
        ! Estimate and update R-Matrix.
        call rmat_updateRmat(obsSpaceData, OBS_VAR, OBS_TRUO)
        call rmat_writeRCorrFile
      else
        ! Compute R matrix based on the difference between simulated obseration and truth
        ! that is based on the true state and simulated emissivity
        ! (H(x_true, emiss_true) + obs_error) - H(x_true, emiss_sim)

        ! Simulate emissivity
        call var1Di_simulateEmissivity(obsSpaceData, simEmissSeed)
        tvs_useSfcEmissObsSpace = .true.

        ! Prepare atmospheric profiles for all tovs observation points for use in rttov
        call tvs_fillProfiles(columnTruthOnTrlLev, obsSpaceData, datestamp, 'nl', beSilent)

        ! Compute radiance
        call tvs_rttov(obsSpaceData, bgckMode, beSilent)

        ! loop over all header indices of the 'TO' family
        call obs_set_current_header_list(obsSpaceData,'TO')

        ! Store the true state (Observation Space) into ObsSpaceData
        HEADER3: do
          headerIndex = obs_getHeaderIndex(obsSpaceData)
          if (headerIndex < 0) exit HEADER3

          ! process only radiance data to be assimilated?
          idatyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
          if (.not. tvs_isIdBurpTovs(idatyp)) then
            write(*,*) 'var1Di_simulateObservation: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
            cycle HEADER3
          end if

          bodyIndexBeg = obs_headElem_i(obsspacedata, OBS_RLN, headerIndex)
          bodyIndexEnd = obs_headElem_i(obsspacedata, OBS_NLV, headerIndex) + bodyIndexBeg - 1

          do bodyIndex = bodyIndexBeg, bodyIndexEnd
            if (obs_bodyElem_i(obsspacedata, OBS_ASS, bodyIndex) == obs_assimilated) then
              call tvs_getChannelNumIndexFromPPP(obsSpaceData, headerIndex, bodyIndex, &
                                                    channelNumber, channelIndex)
              call obs_bodySet_r(obsSpaceData, OBS_ETRU, bodyIndex, tvs_radiance(headerIndex)%bt(channelIndex))
            end if
          end do
        end do HEADER3

        ! Estimate and update R-Matrix.
        call rmat_updateRmat(obsSpaceData, OBS_VAR, OBS_ETRU)
        call rmat_writeRCorrFile
      end if
    end if

    ! Compute the Jacobian
    if (tvs_computeJacobian) then
      call tvs_fillProfiles(columnTruthOnTrlLev, obsSpaceData, datestamp, "nl", beSilent)

      call tvs_rttov_k(columnTruthOnTrlLev, obsSpaceData)
    end if

    write(*,*) 'var1Di_simulateObservation: Finished '
  end subroutine var1Di_simulateObservation

  !--------------------------------------------------------------------------
  ! var1Di_simulateEmissivity
  !--------------------------------------------------------------------------
  subroutine var1Di_simulateEmissivity(obsSpaceData, simEmissSeed)
    !
    !:Purpose: Simulate surface emissivity (Only works for AMSU-A Observations)
    !
    implicit none

    ! Arguments:
    type(struct_obs),        intent(inout) :: obsSpaceData  ! ObsSpaceData object
    integer,                 intent(in)    :: simEmissSeed  ! Seed to random number to simulate emissivity

    ! Locals:
    integer, allocatable         :: emissChanList(:), chanListCMat(:)
    real, allocatable            :: emissStdErr(:)
    character(len=23), parameter :: filenameCorrEmiss = 'Cmat_SfcEmiss_amsua.dat'
    real(8), allocatable         :: emissErrCMat(:,:)
    integer                      :: randomSeed, count, channelNumber, nchanCMat, emissNumChan
    integer                      :: headerIndex, bodyIndex, matchChanIndex, sensorIndex, obsIndex
    integer                      :: bodyIndexBeg, bodyIndexEnd, idatyp
    real(8), allocatable         :: pert(:), emissPert(:), list_EMER(:)
    integer, allocatable         :: list_chanNumber(:), list_bodyIndex(:)
    real(8)                      :: emissErrStdPerChan, emissDiffTmp
    real(8), allocatable         :: sfcEmissivityOriginal(:), sfcEmissivityUpdated(:)
    real(8), parameter           :: missingValueEmisAtlas = -1.0d0

    ! Read Surface Emissivity Error Stdev
    call sse_readEmissError(emissChanList, emissStdErr, emissNumChan)

    ! Read Surface Emissivity Error Correlation Matrix
    call sse_readCEmissMatrixByFileName(filenameCorrEmiss, emissErrCMat, chanListCMat, nchanCMat)

    ! Initialize random number generation for simulating surface emissivity
    if (simEmissSeed == 0) then
      randomSeed = var1Di_randomSeed()
    else
      randomSeed = simEmissSeed + mmpi_myid
    end if

    call rng_setup(abs(randomSeed))

    ! loop over all header indices of the 'TO' family
    call obs_set_current_header_list(obsSpaceData,'TO')

    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER

      ! process only radiance data to be assimilated?
      idatyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if (.not. tvs_isIdBurpTovs(idatyp)) then
        write(*,*) 'var1Di_simulateObservation: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
        cycle HEADER
      end if

      sensorIndex = tvs_lsensor(headerIndex)

      bodyIndexBeg = obs_headElem_i(obsspacedata, OBS_RLN, headerIndex)
      bodyIndexEnd = obs_headElem_i(obsspacedata, OBS_NLV, headerIndex) + bodyIndexBeg - 1

      allocate(pert(emissNumChan))
      allocate(emissPert(emissNumChan))
      allocate(list_EMER(emissNumChan))
      allocate(list_chanNumber(emissNumChan))
      allocate(list_bodyIndex(emissNumChan))

      count = 0
      do bodyIndex = bodyIndexBeg, bodyIndexEnd
        if (obs_bodyElem_i(obsspacedata, OBS_ASS, bodyIndex) == obs_assimilated) then
          count = count + 1
          channelNumber = obs_bodyElem_r(obsspacedata, OBS_PPP, bodyIndex)

          ! Match tvs_channelnumber with the channel index in emissivity error std file
          matchChanIndex = FINDLOC(emissChanList, channelNumber, dim=1)
          if (matchChanIndex == 0) then
            call utl_abort('var1Di_simulateEmissivity: Unable to find emissivity error for a channel')
          end if

            emissErrStdPerChan = emissStdErr(matchChanIndex)
            list_EMER(count) = emissErrStdPerChan
            list_bodyIndex(count) = bodyIndex
            list_chanNumber(count) = channelNumber - tvs_channelOffset(sensorIndex)
            pert(count) = rng_gaussian()
        end if
      end do

      if (count > 0) then
        ! Generate the emissivity errors
        call sse_emissErrMatSqrt(count, pert(1:count), emissPert(1:count), list_chanNumber(1:count), list_EMER(1:count), &
                       emissErrCMat, chanListCMat, nchanCMat)

        allocate(sfcEmissivityOriginal(count))
        allocate(sfcEmissivityUpdated(count))

        do obsIndex = 1, count
          sfcEmissivityOriginal(obsIndex) = obs_bodyElem_r(obsspacedata, OBS_SEM, list_bodyIndex(obsIndex))

          ! Generate the simulated emissivity
          sfcEmissivityUpdated(obsIndex) = sfcEmissivityOriginal(obsIndex) + emissPert(obsIndex)
        end do

        ! QC check: Simulated surface emissivity can only be between 0 and 1 or missing value
        if (any(sfcEmissivityOriginal(:) == missingValueEmisAtlas)) then
          ! Fill missing value if ref column also have missing value
          sfcEmissivityUpdated(:) = missingValueEmisAtlas
        else if(any(sfcEmissivityOriginal(:) < 0.0d0)) then
          ! Fill missing value if ref column have negative emissivity
          sfcEmissivityUpdated(:) = missingValueEmisAtlas
        else if(any(sfcEmissivityOriginal(:) >= 0.0d0) .and. any(sfcEmissivityUpdated(:) < 0.0d0)) then
          ! Limit simulated emissivity to zero if ref column emissivity is equal to zero or greater
          emissDiffTmp = minval(sfcEmissivityUpdated(:)) - 0.0d0
          sfcEmissivityUpdated(:) = sfcEmissivityUpdated(:) - emissDiffTmp * 1.01
        else if(any(sfcEmissivityOriginal(:) <= 1.0d0) .and. any(sfcEmissivityUpdated(:) > 1.0d0)) then
          ! Limit simulated emissivity to one
          emissDiffTmp = maxval(sfcEmissivityUpdated(:)) - 1.0d0
          sfcEmissivityUpdated(:) = sfcEmissivityUpdated(:) - emissDiffTmp * 1.01
        end if

        do obsIndex = 1, count
          call obs_bodySet_r(obsSpaceData, OBS_TSEM, list_bodyIndex(obsIndex), sfcEmissivityOriginal(obsIndex))
          call obs_bodySet_r(obsSpaceData, OBS_SEM, list_bodyIndex(obsIndex), sfcEmissivityUpdated(obsIndex))
          call obs_bodySet_r(obsSpaceData, OBS_EMER, list_bodyIndex(obsIndex), list_EMER(obsIndex))
        end do

        deallocate(sfcEmissivityOriginal)
        deallocate(sfcEmissivityUpdated)

      end if

      deallocate(pert)
      deallocate(emissPert)
      deallocate(list_EMER)
      deallocate(list_chanNumber)
      deallocate(list_bodyIndex)
    end do HEADER
  end subroutine var1Di_simulateEmissivity

  !--------------------------------------------------------------------------
  ! var1Di_estSigmaBObsSpace
  !--------------------------------------------------------------------------
  subroutine var1Di_estSigmaBObsSpace(columnTruthOnTrlLev, numSamplesHBHT, obsSpaceData, vco_anl, datestamp, inflateEmissErr)
    !
    ! :Purpose: Estimating background error STD in observations space
    !
    implicit none

    ! Arguments
    type(struct_columnData), target, intent(in)    :: columnTruthOnTrlLev   ! The true state column
    integer,                         intent(in)    :: numSamplesHBHT        ! Number of samples to compute HBHT
    type(struct_obs),                intent(inout) :: obsSpaceData          ! ObsSpaceData object
    type(struct_vco), pointer,       intent(in)    :: vco_anl               ! Analysis vertical coordinate structure
    integer,                         intent(in)    :: datestamp             ! Date stamp
    real(8),                         intent(in)    :: inflateEmissErr       ! Emissvity error inflation scale factor

    ! Locals:
    type(struct_columnData)         :: columnPertOnAnLev
    type(struct_columnData)         :: columnPertOnTrlLev
    type(struct_columnData)         :: columnSimTrlOnTrlLev
    integer                         :: cvIndex, obsCount
    real(8), allocatable            :: controlVector(:)
    real(8), allocatable            :: errHx(:,:)
    integer, allocatable            :: errHxBodyList(:)
    integer                         :: headerIndex, bodyIndex
    integer                         :: channelNumber, channelIndex
    integer                         :: bodyIndexBeg, bodyIndexEnd, idatyp, obsIndex
    real(8)                         :: meanErrHx, stddevErrHx
    logical                         :: bgckMode, beSilent, nonEmptyBodyColumn, nonEmptyBodyColumn_mpiglobal
    integer                         :: randomSeed, sampleIndex
    real                            :: columnValue

    if (.not. obs_columnActive_RB(obsSpaceData, OBS_TRUO)) then
      call utl_abort('var1Di_estSigmaBObsSpace: The truth in observation space must computed stored ' // &
                      'OBS_TRUO obsSpaceData column')
    end if

    ! Check if the ObsSpace column OBS_HPHT is empty
    nonEmptyBodyColumn = .false.
    HEADERCHCK: do headerIndex = 1, obs_numHeader(obsSpaceData)

      bodyIndexBeg = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
      bodyIndexEnd = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) + bodyIndexBeg - 1

      BODYCHCK: do bodyIndex = bodyIndexBeg, bodyIndexEnd
        columnValue = obs_bodyElem_r(obsSpaceData, OBS_HPHT, bodyIndex)

        if (columnValue /= MPC_missingValue_R8) then
          nonEmptyBodyColumn = .true.
          exit HEADERCHCK
        end if

      end do BODYCHCK
    end do HEADERCHCK

    call mmpi_allReduce(nonEmptyBodyColumn, nonEmptyBodyColumn_mpiglobal, mmpi_lor)

    if (nonEmptyBodyColumn_mpiglobal) then
      call utl_abort('var1Di_estSigmaBObsSpace: ObsSpace column OBS_HPHT is already being used elsewhere')
    end if

    ! Compute background Error in observation space

    allocate(errHx(numSamplesHBHT, obs_numbody_max(obsSpaceData)))
    allocate(errHxBodyList(obs_numbody_max(obsSpaceData)))

    randomSeed = var1Di_randomSeed()

    call rng_setup(abs(randomSeed))

    do sampleIndex = 1, numSamplesHBHT

      allocate(controlVector(cvm_nvadim))

      ! Generate perturbation sampling
      do cvIndex = 1, cvm_nvadim
        controlVector(cvIndex) = rng_gaussian()
      end do

      ! Compute (B^1/2)*Pert (column)
      call col_setVco(columnPertOnAnLev, vco_anl)
      call col_allocate(columnPertOnAnLev, obs_numheader(obsSpaceData), setToZero_opt=.true.)

      if (col_varExist(columnPertOnAnLev, 'EMMW') .and. inflateEmissErr /= MPC_missingValue_R8) then
        call bmat1D_sqrtB(controlVector, cvm_nvadim, columnPertOnAnLev, obsSpaceData, &
                        inflateEmissErr_opt = inflateEmissErr)
      else
        call bmat1D_sqrtB(controlVector, cvm_nvadim, columnPertOnAnLev, obsSpaceData)
      end if

      call col_setVco(columnPertOnTrlLev, col_getVco(columnTruthOnTrlLev))
      call col_allocate(columnPertOnTrlLev, col_getNumCol(columnTruthOnTrlLev), setToZero_opt=.true.)

      ! Interpolate (B^1/2)*Pert from analysis to trial level
      call var1Di_vInterpPertAnLev2TrlLev(columnPertOnAnLev, columnPertOnTrlLev, columnTruthOnTrlLev)

      call col_setVco(columnSimTrlOnTrlLev, col_getVco(columnTruthOnTrlLev))
      call col_allocate(columnSimTrlOnTrlLev, col_getNumCol(columnTruthOnTrlLev), setToZero_opt=.true.)
      call col_copy(columnTruthOnTrlLev, columnSimTrlOnTrlLev)

      ! Add the truth and (B^1/2)*Pert columns
      call col_add(columnPertOnTrlLev, columnSimTrlOnTrlLev)

      ! Compute the pressure levels
      call cvt_transform(columnSimTrlOnTrlLev, 'ZandP_nl')

      ! Restrict the simulated humidity background within physically reasonable values.
      call qlim_rttovLimit(columnSimTrlOnTrlLev)

      beSilent = .false.
      bgckMode = .false.

      ! Prepare atmospheric profiles for all tovs observation points for use in rttov
      call tvs_fillProfiles(columnSimTrlOnTrlLev, obsSpaceData, datestamp, "nl", beSilent)

      ! Compute radiance
      call tvs_rttov(obsSpaceData, bgckMode, beSilent)

      obsCount = 0

      ! loop over all header indices of the 'TO' family
      call obs_set_current_header_list(obsSpaceData,'TO')

      HEADER: do
        headerIndex = obs_getHeaderIndex(obsSpaceData)
        if (headerIndex < 0) exit HEADER

        ! process only radiance data to be assimilated
        idatyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
        if (.not. tvs_isIdBurpTovs(idatyp)) then
          write(*,*) 'var1Di_estSigmaBObsSpace: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
          cycle HEADER
        end if

        bodyIndexBeg = obs_headElem_i(obsspacedata, OBS_RLN, headerIndex)
        bodyIndexEnd = obs_headElem_i(obsspacedata, OBS_NLV, headerIndex) + bodyIndexBeg - 1

        do bodyIndex = bodyIndexBeg, bodyIndexEnd
          if (obs_bodyElem_i(obsspacedata, OBS_ASS, bodyIndex) == obs_assimilated) then
            call tvs_getChannelNumIndexFromPPP(obsSpaceData, headerIndex, bodyIndex, &
                                                channelNumber, channelIndex)
            obsCount = obsCount + 1

            errHx(sampleIndex, obsCount) = tvs_radiance(headerIndex)%bt(channelIndex) - &
                                    obs_bodyElem_r(obsspacedata, OBS_TRUO, bodyIndex)

            errHxBodyList(obsCount) = bodyIndex
          end if
        end do
      end do HEADER

      deallocate(controlVector)
      call col_deallocate(columnPertOnAnLev)
      call col_deallocate(columnPertOnTrlLev)
      call col_deallocate(columnSimTrlOnTrlLev)
    end do

    ! Compute the background error Stdev in observation space
    do obsIndex = 1, obsCount
      meanErrHx = sum(errHx(1:numSamplesHBHT, obsIndex)) / numSamplesHBHT
      stddevErrHx = sqrt(sum((errHx(1:numSamplesHBHT, obsIndex) - meanErrHx)**2) / numSamplesHBHT)
      call obs_bodySet_r(obsSpaceData, OBS_HPHT, errHxBodyList(obsIndex), stddevErrHx)
    end do

  end subroutine var1Di_estSigmaBObsSpace

  !--------------------------------------------------------------------------
  ! var1Di_randomSeed
  !--------------------------------------------------------------------------
  function var1Di_randomSeed() result(randomSeed)
    !
    ! :Purpose: Generate a random seed based on the date stamp and MPI ID
    !
    implicit none

    ! Results:
    integer           :: randomSeed  ! Generated random seed

    ! Locals:
    integer           :: dateStamp, timePrint, datePrint
    integer           :: ierr, imode
    integer, external :: newdate

    imode = -3 ! stamp to printable date and time: YYYYMMDD, HHMMSShh
    dateStamp = tim_getDateStamp()
    ierr = newdate(dateStamp, datePrint, timePrint, imode)
    timePrint = timePrint/1000000
    datePrint =  datePrint*100 + timePrint
    randomSeed = (datePrint - 100000000*(datePrint/100000000)) + mmpi_myid

  end function var1Di_randomSeed

end module var1DIdealize_mod

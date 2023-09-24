module var1DIdealize_mod
    ! MODULE var1DIdealize_mod (prefix='var1D' category='4. Data Object transformations')
    !
    ! :Purpose: contains all 1Dvar-related methods.
    !
    use columnData_mod
    use columnVariableTransforms_mod
    use controlVector_mod
    use gridStatevector_mod
    use horizontalCoord_mod
    use midasMpi_mod 
    use obsSpaceData_mod
    use timeCoord_mod
    use utilities_mod
    use verticalCoord_mod
    use codeprecision_mod
    use mathphysconstants_mod
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
    use tovsNL_mod
    use tovsLin_mod
  
  
    implicit none
    save
    private
  
    ! public procedures
    public :: var1DIdealize_simulateBackgroundState, var1DIdealize_simulateObservation

  contains

  !--------------------------------------------------------------------------
  ! var1DIdealize_simulateBackgroundState
  !--------------------------------------------------------------------------
  subroutine var1DIdealize_simulateBackgroundState(columnTruthOnTrlLev, columnSimTrlOnTrlLev, &
                                                   obsSpaceData, vco_anl, seed)
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
    
    ! locals:
    type(struct_columnData), target :: columnPertOnAnLev
    type(struct_columnData), target :: columnTruthOnAnlLev
    type(struct_columnData), target :: columnPertOnTrlLev
    real(8), allocatable            :: controlVector(:)
    integer                         :: cvIndex
    type(struct_gsv)                :: stateVectorPertOnAnLevTruth
    type(struct_gsv)                :: stateVectorPertOnTrlLevTruth
    type(struct_gsv)                :: stateVectorTrlOnTrlLevTruth
    type(struct_gsv)                :: stateVectorTrlOnTrlLevSim
    type(struct_gsv)                :: stateVectorTrlOnAnlLevTruth
    character(len=50)               :: prefixFileName
    logical                         :: containsFullField

    allocate(controlVector(cvm_nvadim))
    ! Generate perturbation sampling following gaussian distribution with zero mean and one std
    call rng_setup(abs(seed))
    do cvIndex = 1, cvm_nvadim
      controlVector(cvIndex) = rng_gaussian()
    end do

    ! Compute (B^1/2)*Pert (column)
    call col_setVco(columnPertOnAnLev, vco_anl)
    call col_allocate(columnPertOnAnLev, col_getNumCol(columnTruthOnTrlLev), setToZero_opt=.true.)
    call bmat1D_sqrtB(controlVector, cvm_nvadim, columnPertOnAnLev, obsSpaceData)

    call col_setVco(columnPertOnTrlLev, col_getVco(columnTruthOnTrlLev))
    call col_allocate(columnPertOnTrlLev, col_getNumCol(columnTruthOnTrlLev), setToZero_opt=.true.)
    
    ! Interpolate (B^1/2)*Pert from analysis to trial level
    call var1DIdealize_vInterpPertAnLev2TrlLev(columnPertOnAnLev, columnPertOnTrlLev, columnTruthOnTrlLev)

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

    ! Interpolate the truth from trial to analysis increment levels
    call col_setVco(columnTruthOnAnlLev, vco_anl)
    call col_allocate(columnTruthOnAnlLev, col_getNumCol(columnTruthOnTrlLev), setToZero_opt=.true.)
    call inn_setupColumnsOnAnlIncLev(columnTruthOnTrlLev, columnTruthOnAnlLev)

    ! Write trial into standard files
    prefixFileName = 'SimTrialOnTrlLev'
    containsFullField = .true.
    call var1d_transferColumnToYGrid(stateVectorTrlOnTrlLevSim, obsSpaceData, columnSimTrlOnTrlLev, bmat1D_includeAnlVar)
    call var1DIdealize_writeSimTrial(stateVectorTrlOnTrlLevSim, prefixFileName, 'ANALYSIS', containsFullField)

    prefixFileName = 'TruthOnTrlLev'
    containsFullField = .true.
    call var1d_transferColumnToYGrid(stateVectorTrlOnTrlLevTruth, obsSpaceData, columnTruthOnTrlLev, bmat1D_includeAnlVar)
    call var1DIdealize_writeSimTrial(stateVectorTrlOnTrlLevTruth, prefixFileName, 'ANALYSIS', containsFullField)

    prefixFileName = 'TruthOnAnlLev'
    containsFullField = .true.
    call var1d_transferColumnToYGrid(stateVectorTrlOnAnlLevTruth, obsSpaceData, columnTruthOnAnlLev, bmat1D_includeAnlVar)
    call var1DIdealize_writeSimTrial(stateVectorTrlOnAnlLevTruth, prefixFileName, 'ANALYSIS', containsFullField)

    prefixFileName = 'PertOnTrlLev'
    containsFullField = .false.
    call var1d_transferColumnToYGrid(stateVectorPertOnTrlLevTruth, obsSpaceData, columnPertOnTrlLev, bmat1D_includeAnlVar)
    call var1DIdealize_writeSimTrial(stateVectorPertOnTrlLevTruth, prefixFileName, 'INCREMENT', containsFullField)

    prefixFileName = 'PertOnAnlLev'
    containsFullField = .false.
    call var1d_transferColumnToYGrid(stateVectorPertOnAnLevTruth, obsSpaceData, columnPertOnAnLev, bmat1D_includeAnlVar)
    call var1DIdealize_writeSimTrial(stateVectorPertOnAnLevTruth, prefixFileName, 'INCREMENT', containsFullField)

    if (mmpi_myId ==0) then
      call gsv_deallocate(stateVectorPertOnTrlLevTruth)
      call gsv_deallocate(stateVectorPertOnAnLevTruth)
      call gsv_deallocate(stateVectorTrlOnAnlLevTruth)
    end if

    call col_deallocate(columnPertOnTrlLev)
    call col_deallocate(columnPertOnAnLev)
    call col_deallocate(columnTruthOnAnlLev)
    deallocate(controlVector)

  end subroutine var1DIdealize_simulateBackgroundState

  !--------------------------------------------------------------------------
  ! var1DIdealize_vInterpPertAnLev2TrlLev
  !--------------------------------------------------------------------------
  subroutine var1DIdealize_vInterpPertAnLev2TrlLev(columnAnlLev, columnTrlLev, columnPresRef)
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

    write(*,*) 'var1DIdealize_vInterpPertAnLev2TrlLev: Starting'

    ! Check the column size
    if (.not. (col_getNumCol(columnAnlLev) == col_getNumCol(columnTrlLev) .and.    &
        col_getNumCol(columnAnlLev) == col_getNumCol(columnPresRef))) then
      write(*,*) 'Column size columnAnlLev, columnTrlLev and columnPresRef', col_getNumCol(columnAnlLev), &
                  col_getNumCol(columnTrlLev), col_getNumCol(columnPresRef)
      call utl_abort('var1DIdealize_vInterpPertAnLev2TrlLev: The columnAnlLev, columnTrlLev and columnPresRef &
                                 do not have equal number of columns')
    end if

    numColumns = col_getNumCol(columnAnlLev)
    write(*,*) 'var1DIdealize_vInterpPertAnLev2TrlLev: Column size', numColumns

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
          columnTrlLev_ptr  => col_getColumn(columnTrlLev , columnIndex, vnl_varNameList2D(varIndex))
          columnAnlLev_ptr => col_getColumn(columnAnlLev, columnIndex, vnl_varNameList2D(varIndex))
          columnTrlLev_ptr(:) = columnAnlLev_ptr(:)
        end do
      end if
    end do

    deallocate(pSfcRef)
    write(*,*) 'var1DIdealize_vInterpPertAnLev2TrlLev: Finished'
    contains

    logical function varneed(varName)
      implicit none
      ! Arguements: 
      character(len=*) :: varName ! Variable Name

      ! Locals:
      integer          :: varIndex2

      varneed=.false.
      do varIndex2=1,VNL_NUMVARMAX
        if (trim(varName) == trim(bmat1D_includeAnlVar(varIndex2))) then
          varneed=.true.
       end if
      end do

    end function varneed
  end subroutine var1DIdealize_vInterpPertAnLev2TrlLev

  !--------------------------------------------------------------------------
  ! var1DIdealize_writeSimTrial
  !--------------------------------------------------------------------------
  subroutine var1DIdealize_writeSimTrial(statevectorSim, prefixFileName, etiket, containsFullField)
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
   
    if(mmpi_myid == 0) write(*,*) 'var1DIdealize_writeSimTrial: STARTING'

    ! loop over times for which increment is computed
    do stepIndex = 1, tim_nstepobsinc
      if (gsv_isAllocated(statevectorSim)) then
        dateStamp = gsv_getDateStamp(statevectorSim,stepIndex)
        if (mmpi_myid == 0) write(*,*) 'var1DIdealize_writeSimTrial: writing for time step: ',stepIndex, dateStamp

        ! write the increment file for this time step
        call difdatr(dateStamp,tim_getDatestamp(),deltaHours)
        if (nint(deltaHours*60.0d0).lt.0) then
          write(coffset,'(I4.3)') nint(deltaHours*60.0d0)
        else
          write(coffset,'(I3.3)') nint(deltaHours*60.0d0)
        end if

        fileName = './'//trim(prefixFileName)//'_' // trim(coffset) // 'm'
        call gio_writeToFile( statevectorSim, fileName, trim(etiket), scaleFactor_opt = 1.0d0, &
                              ip3_opt = 0, stepIndex_opt = stepIndex, containsFullField_opt=containsFullField )
      end if
    end do

    if(mmpi_myid == 0) write(*,*) 'var1DIdealize_writeSimTrial: Finished'
  end subroutine var1DIdealize_writeSimTrial

  !--------------------------------------------------------------------------
  ! var1DIdealize_simulateObservation
  !--------------------------------------------------------------------------
  subroutine var1DIdealize_simulateObservation(columnTrlOnTrlLevTruth, obsSpaceData, datestamp, seed, useSimObsErr, varMode)
    !
    !:Purpose: Simulate the observation (only TOVS obs) by adding a perturbation from the reference data
    !          Additional changes are needed to generalize for all observations (not just TOVS obs)
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)    :: columnTrlOnTrlLevTruth ! True column state
    type(struct_obs),        intent(inout) :: obsSpaceData           ! ObsSpacedata object
    integer,                 intent(in)    :: datestamp              ! Date stamp
    integer,                 intent(in)    :: seed                   ! Seed to random number generator 
    logical,                 intent(in)    :: useSimObsErr           ! Simulate Observation Error Covariance
    character(len=*),        intent(in)    :: varMode                ! Variational Mode

    ! Locals:
    logical              :: bgckMode, beSilent
    integer              :: tovsIndex
    integer              :: headerIndex, bodyIndex, obsIndex, sensorIdIndex
    integer              :: idata, idatend, idatyp, count, channelNumber, channelIndex, tmpObs_OER
    real(8), allocatable :: pert(:), obsPert(:), list_OER(:)
    integer, allocatable :: list_chanNumber(:), list_bodyIndex(:), list_chanIndex(:)
    real(8), allocatable :: obsErrStdev(:,:)
    type(rmat_matrix), allocatable, target :: estR(:) 
    
    beSilent = .false.
    bgckMode = .false.

    write(*,*) 'var1DIdealize_simulateObservation: Starting'
    
    ! Compute the Truth the observation space
    write(*,*) 'var1DIdealize_simulateObservation: Computing the truth in Obs Space'

    ! Prepare atmospheric profiles for all tovs observation points for use in rttov
    call tvs_fillProfiles(columnTrlOnTrlLevTruth, obsSpaceData, datestamp, "nl", beSilent)

    ! Compute radiance
    call tvs_rttov(obsSpaceData, bgckMode, beSilent)

     ! loop over all header indices of the 'TO' family
    call obs_set_current_header_list(obsSpaceData,'TO')

    ! Store the true state (Observation Space) into ObsSpaceData
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER

      ! process only radiance data to be assimilated?
      idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
      if ( .not. tvs_isIdBurpTovs(idatyp) ) then
        write(*,*) 'var1DIdealize_simulateObservation: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
        cycle HEADER
      end if

      tovsIndex = tvs_tovsIndex(headerIndex)
      if (tovsIndex == -1) cycle HEADER

      idata   = obs_headElem_i(obsspacedata, OBS_RLN, headerIndex)
      idatend = obs_headElem_i(obsspacedata, OBS_NLV, headerIndex) + idata - 1

      do bodyIndex = idata, idatend
        if (obs_bodyElem_i(obsspacedata, OBS_ASS, bodyIndex ) == obs_assimilated) then
          call tvs_getChannelNumIndexFromPPP(obsSpaceData, headerIndex, bodyIndex, &
                                                channelNumber, channelIndex)
          call obs_bodySet_r(obsSpaceData, OBS_TRUO, bodyIndex, tvs_radiance(tovsIndex)%bt(channelIndex))
        end if 
      end do
    end do HEADER

    ! Generate Simulated Observations
    write(*,*) 'var1DIdealize_simulateObservation: Use simulated Obs and Emissivity Errors, useSimObsErr ', useSimObsErr
    
    if(useSimObsErr) then
      ! Prepare atmospheric profiles for all tovs observation points for use in rttov
      call tvs_fillProfiles(columnTrlOnTrlLevTruth, obsSpaceData, datestamp, "nl", beSilent)

      ! Compute radiance
      call tvs_rttov(obsSpaceData, bgckMode, beSilent, SimSfcEmiss_opt = .True.)
    end if

     ! loop over all header indices of the 'TO' family
    call obs_set_current_header_list(obsSpaceData,'TO')

    call rng_setup(abs(seed + mmpi_myid))

    HEADER2: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER2

      ! process only radiance data to be assimilated?
      idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
      if (.not. tvs_isIdBurpTovs(idatyp)) then
        write(*,*) 'var1DIdealize_simulateObservation: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
        cycle HEADER2
      end if

      tovsIndex = tvs_tovsIndex(headerIndex)
      if (tovsIndex == -1) cycle HEADER2

      idata   = obs_headElem_i(obsspacedata, OBS_RLN, headerIndex)
      idatend = obs_headElem_i(obsspacedata, OBS_NLV, headerIndex) + idata - 1

      if (tvs_isIdBurpTovs(idatyp)) then

        allocate(pert(tvs_maxChannelNumber))
        allocate(obsPert(tvs_maxChannelNumber))
        allocate(list_OER(tvs_maxChannelNumber))
        allocate(list_chanNumber(tvs_maxChannelNumber))
        allocate(list_chanIndex(tvs_maxChannelNumber))
        allocate(list_bodyIndex(tvs_maxChannelNumber))

        ! Read the Sigma O from ObsSpaceData
        count = 0
        do bodyIndex = idata, idatend
          if (obs_bodyElem_i(obsspacedata, OBS_ASS, bodyIndex ) == obs_assimilated) then
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

        if (.not. count > 0) then 
          if (allocated(pert)) deallocate(pert)
          if (allocated(obsPert)) deallocate(obsPert)
          if (allocated(list_OER)) deallocate(list_OER)
          if (allocated(list_chanNumber)) deallocate(list_chanNumber)
          if (allocated(list_chanIndex)) deallocate(list_chanIndex)
          if (allocated(list_bodyIndex)) deallocate(list_bodyIndex)
          cycle HEADER2
        end if
        
        ! Compute Observation Perturbation
        call rmat_Rsqrt(tvs_lsensor(tvs_tovsIndex(headerIndex)), count, pert(1:count), obsPert(1:count), list_chanNumber(1:count),&
                         list_OER(1:count), tvs_tovsIndex(headerIndex))

        ! Update the obs value in ObsSpacedata
        do obsIndex = 1, count
          call obs_bodySet_r(obsSpaceData, OBS_VAR, list_bodyIndex(obsIndex), tvs_radiance(tovsIndex)%bt(list_chanIndex(obsIndex)) + obsPert(obsIndex))
        end do

        if (allocated(pert)) deallocate(pert)
        if (allocated(obsPert)) deallocate(obsPert)
        if (allocated(list_OER)) deallocate(list_OER)
        if (allocated(list_chanNumber)) deallocate(list_chanNumber)
        if (allocated(list_chanIndex)) deallocate(list_chanIndex)
        if (allocated(list_bodyIndex)) deallocate(list_bodyIndex)
      end if
    end do HEADER2
    
    ! Estimate and update R-Matrix.
    allocate(estR(tvs_nsensors))

    call rmat_estimateR(obsSpaceData, estR)
    call rmat_updateRmat(estR, obsSpaceData)
    call rmat_writeRCorrFile

    deallocate(estR)

    write(*,*) 'Finish var1DIdealize_simulateObservation'
  end subroutine var1DIdealize_simulateObservation

end module var1DIdealize_mod

  

   

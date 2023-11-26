module surfaceEmissivity_mod 
  ! MODULE surfaceEmissivity_mod  (prefix='sse ' category='7. Low-level data objects')
  !
  !:Purpose: Manipulate surface emissivity for idealized assimilation of 
  !          surface sensitive AMSU-A channels.
  ! 
  use obsSpaceData_mod
  use rttov_types, only : rttov_emissivity, rttov_chanprof
  use rttov_const, only : surftype_land, surftype_seaice
  use MathPhysConstants_mod
  use midasMpi_mod
  use randomNumber_mod
  use utilities_mod
  use columnData_mod
  use verticalCoord_mod
  use varNameList_mod

  implicit none
  save
  private

  ! public procedures  
  public :: sse_simulateEmissivity, sse_setupEmissivityfromState, sse_extractEmissivityCol
  public :: sse_updateBEnsMatEmissFromBHi

contains

  !--------------------------------------------------------------------------
  ! sse_simulateEmissivity
  !--------------------------------------------------------------------------
  subroutine sse_simulateEmissivity(obsSpaceData, sensorTovsIndexes, sfcEmissivity, sensorId, &
                                    nSensor, sensorList, tovsHeaderIndexList, tovsIndexList, &
                                    tovsChannelOffset, tovsIChan, tovsMaxChanneNum, instrumentName)
    !
    !:Purpose: Simulate surface emissivity (Only works for AMSU-A channel 1-5)
    !
    implicit none

    ! Arguments:
    type(struct_obs),        intent(inout) :: obsSpaceData           ! ObsSpaceData object
    integer,                 intent(in)    :: sensorTovsIndexes(:)   ! Sensor TOVS Indexes
    type (rttov_emissivity), intent(inout) :: sfcEmissivity(:)       ! Simulated Surf. Emissivity
    integer,                 intent(in)    :: sensorId               ! Sensor ID
    integer,                 intent(in)    :: nSensor                ! Number of individual sensors.(tvs_nsenors)
    integer,                 intent(in)    :: sensorList(:)          ! Sensor number for each profile (tvs_lsensor)
    integer,                 intent(in)    :: tovsHeaderIndexList(:) ! Observation position in obsSpaceData header for each profile (tvs_headerIndex)
    integer,                 intent(in)    :: tovsIndexList(:)       ! Index in TOVS structures for each observation header in obsSpaceData (tvs_tovsIndex)
    integer,                 intent(in)    :: tovsChannelOffset(:)   ! RTTOV channel mapping offset (tvs_channelOffset)
    integer,                 intent(in)    :: tovsIChan(:,:)         ! List of channels per instrument (tvs_ichan)
    integer,                 intent(in)    :: tovsMaxChanneNum       ! Max. value for channel number (tvs_maxChannelNumber)
    character(len=15),       intent(in)    :: instrumentName(:)      ! Satellite name (tvs_instrumentName)

    ! Locals:
    integer, allocatable :: emissChanList(:)
    real, allocatable    :: emissStdErr(:)
    integer              :: emissNumChan
    integer              :: idata, idatend
    real(8), allocatable :: pert(:), emissPert(:), list_EMER(:)
    integer, allocatable :: list_chanNumber(:), list_bodyIndex(:), list_chanIndex(:), list_btIndex(:)
    integer              :: headerIndex, bodyIndex, obsIndex, count, channelNumber, channelIndex, tovsIndex, matchChanIndex
    real(8)              :: emissErrStdPerChan
    character (len=64)   :: filenameCorrEmiss
    integer              :: btIndex, profileIndex, sensorIndex
    real(8), allocatable :: emissErrCMat(:,:)
    integer, allocatable :: chanListCMat(:)
    integer              :: nchanCMat
    integer              :: ierr

    ! Namelist variables:
    integer   :: simEmissSeed   ! Seed used to simulate surface emissivity
    integer   :: maxSfcSenChan  ! Highest peaking AMSUA channels to be considered to simulate surface meissivity

    namelist /namSimEmiss/ simEmissSeed, maxSfcSenChan

    ! Set default values of namelist variables
    simEmissSeed = MPC_missingValue_INT
    maxSfcSenChan = MPC_missingValue_INT

    write(*,*) 'sse_simulateEmissivity: STARTING'

    do sensorIndex = 1, nSensor
      if (trim(instrumentName(sensorIndex)) /= 'AMSUA') then
        write(*,*) 'Satellite Name: ', trim(instrumentName(sensorIndex))
        call utl_abort('sse_simulateEmissivity: simulating surface emissivity can only be used for AMSUA instrument')
      end if
    end do

    ! Read the namelist for directives
    call utl_tmg_start(181,'low-level--readNML')
    read(utl_flnml, nml = namSimEmiss, iostat = ierr)
    call utl_tmg_stop(181)
    if (ierr /= 0) call utl_abort('sse_simulateEmissivity: Error reading namelist')
    if (mmpi_myid == 0) write(*, nml = namSimEmiss)

    if (simEmissSeed == MPC_missingValue_INT .or. maxSfcSenChan == MPC_missingValue_INT) then 
      call utl_abort('sse_simulateEmissivity: namelist variable not set')
    end if

    ! Read Surface Emissivity Error Stdev
    call sse_readEmissError(emissChanList, emissStdErr, emissNumChan, maxSfcSenChan)

    filenameCorrEmiss = 'Cmat_SfcEmiss_amsua.dat'

    ! Read Surface Emissivity Error Correlation Matrix
    call sse_readCEmissMatrixByFileName(filenameCorrEmiss, emissErrCMat, chanListCMat, nchanCMat)
    
    ! Initialize random number generation for simulating surface emissivity
    call rng_setup(abs(simEmissSeed + (mmpi_myid * nSensor) + sensorId))

    ! The count of btIndex follows the same logic in tvs_countRadiances
    btIndex = 0

    PROFILE: do profileIndex = 1, size(sensorTovsIndexes)
      tovsIndex = sensorTovsIndexes(profileIndex)
      sensorIndex = sensorList(tovsIndex)
      headerIndex = tovsHeaderIndexList(tovsIndex)
      if (headerIndex < 0) cycle PROFILE

      idata = obs_headElem_i(obsspacedata, OBS_RLN, headerIndex)
      idatend = obs_headElem_i(obsspacedata, OBS_NLV, headerIndex) + idata - 1
      
      allocate(pert(emissNumChan))
      allocate(emissPert(emissNumChan))
      allocate(list_EMER(emissNumChan))
      allocate(list_chanNumber(emissNumChan))
      allocate(list_chanIndex(emissNumChan))
      allocate(list_bodyIndex(emissNumChan))
      allocate(list_btIndex(emissNumChan))

      count = 0
      do bodyIndex = idata, idatend
        if (obs_bodyElem_i(obsspacedata, OBS_ASS, bodyIndex) == obs_assimilated) then

          btIndex = btIndex + 1

          call sse_getChannelNumIndexFromPPP(obsSpaceData, headerIndex, bodyIndex, &
                                            channelNumber, channelIndex, tovsIndexList, &
                                            sensorList, tovsMaxChanneNum, tovsChannelOffset, &
                                            tovsIChan)

          if (channelNumber + tovsChannelOffset(sensorIndex) <= maxSfcSenChan) then
            count = count + 1

            ! Match tvs_channelnumber with the channel index in emissivity error std file
            matchChanIndex = FINDLOC(emissChanList, channelNumber + tovsChannelOffset(sensorIndex), dim=1)
            if (matchChanIndex == 0) then
              call utl_abort('sse_simulateEmissivity: Unable to find emissivity error for a channel')
            end if
            emissErrStdPerChan = emissStdErr(matchChanIndex)
            list_EMER(count) = emissErrStdPerChan
            list_bodyIndex(count) = bodyIndex
            list_chanNumber(count) = channelNumber
            list_chanIndex(count) = channelIndex
            pert(count) = rng_gaussian()
            list_btIndex(count) = btIndex
          end if
        end if
      end do

      if (count > 0 .and. tovsIndex > 0) then 
        ! Generate the emissivity errors
        call sse_emissErrMatSqrt(count, pert(1:count), emissPert(1:count), list_chanNumber(1:count), list_EMER(1:count), &
                       emissErrCMat, chanListCMat, nchanCMat)

        do obsIndex = 1, count
          ! Generate the simulated emissivity
          sfcEmissivity(list_btIndex(obsIndex))%emis_in = sfcEmissivity(list_btIndex(obsIndex))%emis_in + emissPert(obsIndex)

          !Store simulated emissivity error STDev. into ObsSpaceData
          call obs_bodySet_r(obsSpaceData, OBS_EMER, list_bodyIndex(obsIndex), list_EMER(obsIndex))
        
          ! QC check: Simulated surface emissivity can only be between 0 and 1
          if (sfcEmissivity(list_btIndex(obsIndex))%emis_in > 1.0d0) then
            sfcEmissivity(list_btIndex(obsIndex))%emis_in = 1.0d0
          else if (sfcEmissivity(list_btIndex(obsIndex))%emis_in < 0.0d0) then
            sfcEmissivity(list_btIndex(obsIndex))%emis_in = 0.0d0
          end if

        end do
      end if

      deallocate(pert)
      deallocate(emissPert)
      deallocate(list_EMER)
      deallocate(list_chanNumber)
      deallocate(list_chanIndex)
      deallocate(list_bodyIndex)
      deallocate(list_btIndex)

    end do PROFILE

    deallocate(emissErrCMat)
    deallocate(chanListCMat)
    deallocate(emissChanList)
    deallocate(emissStdErr)

    write(*,*) 'sse_simulateEmissivity: FINISHED'
  end subroutine sse_simulateEmissivity

  !--------------------------------------------------------------------------
  ! sse_readEmissError
  !--------------------------------------------------------------------------
  subroutine sse_readEmissError(chanList, emissStdErr, numChan, maxSfcSenChan)
    !
    !:Purpose: Reading emissivity error std for AMSU-A
    !
    implicit none

    ! Arguemnts:
    integer, allocatable, intent(out) :: chanList(:)    ! Channel List
    real, allocatable,    intent(out) :: emissStdErr(:) ! Surf. Emissivity Error Stdev.
    integer,              intent(out) :: numChan        ! Number of Channels
    integer,              intent(in)  :: maxSfcSenChan  ! High peaking surface sensitive channel

    ! Locals:
    integer, external      :: fnom, fclos
    integer                :: iunit, ierr
    character(len=20)      :: fileName
    character(len = 512)   :: tmpStr
    integer                :: chanIndex

    fileName = 'EmissErrStd.dat'
    iunit = 0
    ierr = fnom(iunit, fileName, 'SEQ+FMT', 0)
    if (ierr /= 0) call utl_abort('sse_readEmissError: Error reading &
                                   surface emissivity error stdev. file') 

    ! Read temporary strings 
    read(iunit, *) tmpStr
    read(iunit, '(A)') tmpStr
    read(iunit, '(A)') tmpStr
    read(iunit, '(A)') tmpStr

    ! Read number of channels
    read(iunit, *) numChan

    read(iunit, '(A)') tmpStr

    allocate(chanList(numChan))
    allocate(emissStdErr(numChan))

    ! Read Emissivity Error Stdev.
    do chanIndex = 1, numChan
      read(iunit,*) chanList(chanIndex), emissStdErr(chanIndex)
    end do
    read(iunit, '(A)') tmpStr
  
    ierr = fclos(iunit)
    if (ierr /= 0) call utl_abort ('sse_readEmissError')

    if (maxSfcSenChan > maxval(chanList)) then
      write(*,*) 'The selected highest peaking channel is not in the emissivity error std file'
      call utl_abort('sse_readEmissError')
    end if

  end subroutine sse_readEmissError

  !--------------------------------------------------------------------------
  ! sse_readCEmissMatrixByFileName
  !--------------------------------------------------------------------------
  subroutine sse_readCEmissMatrixByFileName(infile, corrMat, chanList, nchan)
    !
    !:Purpose:  Read surface emissivity error correlation file
    !
    implicit none

    ! Arguments:
    character (len=*),    intent(in)  :: infile       ! Name of input file
    real(8), allocatable, intent(out) :: corrMat(:,:) ! Emissivity Error Correlation matrix
    integer, allocatable, intent(out) :: chanList(:)  ! Channel list
    integer,              intent(out) :: nchan        ! Number of channels

    ! Locals:
    integer              :: chanIndex1, chanIndex2
    integer              :: iunit, ierr, count, readChanIndex
    integer, external    :: fnom,fclos
    real(8)              :: tmpCor
    integer, allocatable :: foundChanIndex(:)
 
    iunit = 0
    ierr = fnom(iunit, trim(infile), 'FTN+SEQ+R/O', 0)
    if (ierr /= 0) then
      write(*,*) 'Cannot open ' // trim(infile)
      call utl_abort('sse_readCEmissMatrixByFileName')
    end if

    write(*,*) 'sse_readCEmissMatrixByFileName: Reading ' // trim(infile)
    
    read(iunit,*) nchan
  
    allocate(foundChanIndex(nchan))
    allocate(corrMat(nchan, nchan))
    allocate(chanList(nchan))
    corrMat = 0.d0

    do chanIndex1 = 1, nchan
      corrMat(chanIndex1, chanIndex1) = 1.d0
    end do
    
    count = 0
    foundChanIndex(:) = -1
    do chanIndex1 = 1, nchan
      read(iunit,*) readChanIndex
      foundChanIndex(chanIndex1) = chanIndex1
      chanList(chanIndex1) = readChanIndex
      count = count + 1
    end do
    if (count /= nchan) then
      write(*,*) 'Warning: Missing information in ' // trim(infile)
    end if

    do
      read(iunit, *, iostat=ierr) chanIndex1, chanIndex2, tmpCor
      if (ierr /= 0) exit
      if (foundChanIndex(chanIndex1) /= -1 .and. foundChanIndex(chanIndex2) /= -1) then
        corrMat(foundChanIndex(chanIndex1), foundChanIndex(chanIndex2)) = tmpCor
        corrMat(foundChanIndex(chanIndex2), foundChanIndex(chanIndex1)) = tmpCor
      end if
    end do

    ierr= fclos(iunit)
    deallocate(foundChanIndex)

  end subroutine sse_readCEmissMatrixByFileName

  !--------------------------------------------------------------------------
  ! sse_emissErrMatSqrt
  !--------------------------------------------------------------------------
  subroutine sse_emissErrMatSqrt(nsubset, obsIn, obsOut, list_sub, list_oer, Cmat, chanList, nchan)
    !
    ! :Purpose: Apply the operator EmissivityErrorCovarianceMatrix**1/2 to obsIn
    !           result in obsOut for the subset of channels specified
    !           in list_sub
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: nsubset           ! Number of subset channels in Emissivity Error Covariance Matrix
    integer, intent(in)  :: list_sub(nsubset) ! List of subset channels in Emissivity Error Covariance Matrix
    real(8), intent(in)  :: list_oer(nsubset) ! List of Emissivity Error Stdev
    real(8), intent(in)  :: obsIn(nsubset)    ! Sampling Perturbation
    integer, intent(in)  :: nchan             ! Number of channels in correlation matrix
    integer, intent(in)  :: chanList(nchan)   ! Channel list in correlation matrix
    real(8), intent(in)  :: Cmat(nchan,nchan) ! Emissivity error correlation matrix
    real(8), intent(out) :: obsOut(nsubset)   ! Error Perturbation
    
    ! Locals:
    real(8)               :: alpha, beta, variance
    real(8), allocatable  :: emissErrMat(:,:)
    integer               :: foundChanIndex(nsubset)
    integer               :: chanIndex1, chanIndex2
  
    allocate(emissErrMat(nsubset, nsubset))

    foundChanIndex(:) = -1
    do chanIndex1 = 1, nsubset
      chanLoop: do chanIndex2 = 1, nchan
        if (list_sub(chanIndex1) == chanList(chanIndex2)) then
          foundChanIndex(chanIndex1) = chanIndex2
          exit chanLoop 
        end if
      end do chanLoop
    end do

    if (any(foundChanIndex == -1)) then
      write(*,*) 'Missing information for some channel !'
      write(*,*) list_sub(:)
      write(*,*) foundChanIndex(:)
      call utl_abort('sse_emissErrMatSqrt')
    end if

    do chanIndex2 = 1, nsubset
      do chanIndex1 = 1, nsubset
        variance = list_oer(chanIndex1) * list_oer(chanIndex2)
        emissErrMat(chanIndex1,chanIndex2) = variance * Cmat(foundChanIndex(chanIndex1), foundChanIndex(chanIndex2))
      end do
    end do

    ! Calculation of F**1/2
    call utl_matSqrt(emissErrMat, nsubset, 1.d0, .false.)

    alpha = 1.d0
    beta = 0.d0
    obsOut = 0.d0

    ! Optimized symetric matrix vector product from Lapack
    call dsymv('L', nsubset, alpha, emissErrMat, nsubset, obsIn, 1, beta, obsOut, 1)

    deallocate(emissErrMat)

  end subroutine sse_emissErrMatSqrt

  !--------------------------------------------------------------------------
  !  sse_setupEmissivityfromState
  !--------------------------------------------------------------------------
  subroutine sse_setupEmissivityfromState(emissivity_inout, obsSpaceData, bodyIndexFromBtIndex, & 
                                   chanprof, sensorTovsIndexes, tovsheaderIndex, tovsIndexList, &
                                   nSensor, sensorList, instrumentName, tovsMaxChanneNum, &
                                   tovsChannelOffset, tovsIChan, surftype, &
                                   columnIn_opt, columnInout_opt, emissivityTovsIndex_opt, originalEmissivity_opt)
    !
    !:Purpose: Setup the emissivity_inout for tangent linear and adjoint of computation of 
    !          radiance with rttov_tl and rttov_ad
    !
    implicit none
  

    ! Arguments:
    type(rttov_emissivity), pointer, intent(inout) :: emissivity_inout(:)
    type(struct_obs),                intent(in)    :: obsSpaceData           ! ObsSpaceData object
    integer,                         intent(in)    :: bodyIndexFromBtIndex(:)
    type(rttov_chanprof),            intent(in)    :: chanprof(:)
    integer,                         intent(in)    :: sensorTovsIndexes(:)
    integer,                         intent(in)    :: tovsheaderIndex(:)
    integer,                         intent(in)    :: tovsIndexList(:)        ! Tovs Index List (tvs_tovsIndex)
    integer,                         intent(in)    :: nSensor
    integer,                         intent(in)    :: sensorList(:)           ! Sensor number for each profile (tvs_lsensor)
    character(len=15),               intent(in)    :: instrumentName(:)      ! Satellite name (tvs_instrumentName)
    integer,                         intent(in)    :: tovsMaxChanneNum       ! Max. value for channel number (tvs_maxChannelNumber)
    integer,                         intent(in)    :: tovsChannelOffset(:)    ! RTTOV channel mapping offset (tvs_channelOffset)
    integer,                         intent(in)    :: tovsIChan(:,:)          ! List of channels per instrument (tvs_ichan)
    integer,                         intent(in)    :: surftype(:)
    type(struct_columnData), target, optional,        intent(inout)    :: columnInout_opt
    type(struct_columnData), target, optional,        intent(in)    :: columnIn_opt
    real(8), optional, target,                      intent(in)    :: emissivityTovsIndex_opt(:,:)
    real(8), optional, target,                     intent(in) :: originalEmissivity_opt(:)


    ! Locals:
    integer :: profileIndex, btIndex, sensorIndex, tovsIndex
    integer :: bodyIndex, headerIndex
    integer :: channelNumber, channelIndex
    integer :: btCount, numSourceFile
    real(8), pointer :: emissivityCol(:), emissivityTovsIndex(:,:), originalEmissivity(:)
    type(struct_columnData), pointer :: column
    character(len=2)  :: profType

    write(*,*) 'sse_setupEmissivityfromState: STARTING'

    do sensorIndex = 1, nSensor
      if (trim(instrumentName(sensorIndex)) /= 'AMSUA') then
        write(*,*) 'Satellite Name: ', trim(instrumentName(sensorIndex))
        call utl_abort('sse_setupEmissivityfromState: !! can only be used for AMSUA instrument')
      end if
    end do

    numSourceFile = 0

    ! Setup for the tangent linear computation for surface emissivity
    if (present(columnIn_opt)) then 
      column => columnIn_opt
      profType = 'tl'
      numSourceFile = numSourceFile + 1
    end if

    ! Setup for the adjoint computation for surface emissivity
    if (present(columnInout_opt)) then 
      column => columnInout_opt
      profType = 'ad'
      numSourceFile = numSourceFile + 1
    end if

    ! Extracting surface emissivity from background trial
    if (present(emissivityTovsIndex_opt)) then 
      emissivityTovsIndex => emissivityTovsIndex_opt
      profType = 'dt'
      numSourceFile = numSourceFile + 1

      if (present(originalEmissivity_opt)) then 
        originalEmissivity => originalEmissivity_opt
      else 
        call utl_abort('sse_setupEmissivityfromState: originalEmissivity_opt must be supplied when &
                       extracting surface emissivity from background trial')
      end if
    end if

    if (numSourceFile /= 1) then
      call utl_abort('sse_setupEmissivityfromState: Only one source of emissivity variable must be supplied')
    end if

    btCount = size(emissivity_inout)
    do btIndex = 1, btCount
      bodyIndex = bodyIndexFromBtIndex(btIndex)
      profileIndex = chanprof(btIndex)%prof
      tovsIndex  = sensorTovsIndexes(profileIndex)
      headerIndex = tovsheaderIndex(tovsIndex)

      call sse_getChannelNumIndexFromPPP(obsSpaceData, headerIndex, bodyIndex, &
                                              channelNumber, channelIndex, tovsIndexList, &
                                              sensorList, tovsMaxChanneNum, tovsChannelOffset, &
                                              tovsIChan)

      if (profType == 'tl') then
        emissivityCol => col_getColumn(column, headerIndex,'EMMW')
        emissivity_inout(btIndex)%emis_in = emissivityCol(channelNumber)
      else if (profType == 'ad') then
        emissivityCol => col_getColumn(column, headerIndex,'EMMW')
        emissivityCol(channelNumber) = emissivity_inout(btIndex)%emis_out
      else if (profType == 'dt') then
        if (surftype(tovsIndex) == surftype_land .and. &
            emissivityTovsIndex(tovsIndex, channelNumber) > 1.d0) then
        end if
        if (surftype(tovsIndex) == surftype_land .and. &
           emissivityTovsIndex(tovsIndex, channelNumber) > 0.d0 .and. &
           emissivityTovsIndex(tovsIndex, channelNumber) <= 1.d0) then

          emissivity_inout(btIndex)%emis_in = emissivityTovsIndex(tovsIndex, channelNumber)
        else
          emissivity_inout(btIndex)%emis_in = originalEmissivity(btIndex)
        end if
      end if
    end do
  
  end subroutine sse_setupEmissivityfromState

  !--------------------------------------------------------------------------
  !  sse_extractEmissivityCol
  !--------------------------------------------------------------------------
  subroutine sse_extractEmissivityCol(column, emissivityFromTrl, profileCount, sensorTovsIndexes, sensorHeaderIndexes, nobtov)
    !
    !:Purpose: Extract emissivity from column objects
    !
    implicit none

    ! Arguments:
    type(struct_columnData),         intent(in)      :: column             ! column structure
    real(8), allocatable,            intent(inout)   :: emissivityFromTrl(:, :)
    integer,                         intent(in)      :: profileCount
    integer,                         intent(in)      :: sensorTovsIndexes(:)
    integer,                         intent(in)      :: sensorHeaderIndexes(:)
    integer,                         intent(in)      :: nobtov
    
    ! Locals:
    integer :: emissTrlNumChan
    integer :: profileIndex, tovsIndex, headerIndex
 
    if (.not. col_varExist(column, 'EMMW')) return

    emissTrlNumChan = col_getNumLev(column, 'OT', varName_opt = 'EMMW')
    
    if (.not. allocated(emissivityFromTrl)) allocate(emissivityFromTrl(nobtov, emissTrlNumChan))

    do  profileIndex = 1 , profileCount 
      tovsIndex = sensorTovsIndexes(profileIndex)
      headerIndex = sensorHeaderIndexes(profileIndex)

      emissivityFromTrl(tovsIndex, :) = col_getColumn(column, headerIndex, 'EMMW')
    end do
  end subroutine sse_extractEmissivityCol

  !--------------------------------------------------------------------------
  !  sse_getChannelNumIndexFromPPP
  !--------------------------------------------------------------------------
  subroutine sse_getChannelNumIndexFromPPP(obsSpaceData, headerIndex, bodyIndex, &
                                           channelNumber, channelIndex, tovsIndexList, &
                                           sensorList, tovsMaxChanneNum, tovsChannelOffset, &
                                           tovsIChan)
    !
    !:Purpose: Get channel number/index from obs_ppp for TO observations. 
    !          Based on tvs_getChannelNumIndexFromPPP subroutine in tovsNL_mod module. 
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(in)  :: obsSpaceData            ! ObsSpaceData object
    integer,          intent(in)  :: headerIndex             ! Header Index
    integer,          intent(in)  :: bodyIndex               ! Body Index
    integer,          intent(in)  :: tovsIndexList(:)        ! Tovs Index List (tvs_tovsIndex)
    integer,          intent(in)  :: sensorList(:)           ! Sensor number for each profile (tvs_lsensor)
    integer,          intent(in)  :: tovsMaxChanneNum        ! Max. value for channel number (tvs_maxChannelNumber)
    integer,          intent(in)  :: tovsChannelOffset(:)    ! RTTOV channel mapping offset (tvs_channelOffset)
    integer,          intent(in)  :: tovsIChan(:,:)          ! List of channels per instrument (local)
    integer,          intent(out) :: channelNumber           ! Channel Number
    integer,          intent(out) :: channelIndex            ! Channel Index

    ! Locals:
    integer :: tovsIndex, sensorIndex

    tovsIndex = tovsIndexList(headerIndex)
    sensorIndex = sensorList(tovsIndex)

    channelNumber = nint(obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex))
    channelNumber = max(0, min(channelNumber, tovsMaxChanneNum + 1))
    channelNumber = channelNumber - tovsChannelOffset(sensorIndex)
    channelIndex = utl_findloc(tovsIChan(:, sensorIndex), channelNumber)

  end subroutine sse_getChannelNumIndexFromPPP

  !--------------------------------------------------------------------------
  !  sse_extractEmissivityCol
  !--------------------------------------------------------------------------
  subroutine sse_updateBEnsMatEmissFromBHi(BmatHiLand, BmatHiSea, BmatEns, latLandHi, latSeaHi, obsSpaceData, & 
                                           validHeaderIndex, validHeaderCount, bmat1D_includeAnlVar, vco_1DVarBmat, vco_1DVar)
    !
    !:Purpose: Extract surface emissivity elements in the static B-Matrix and 
    !          copy to the ensemble B-Matrix. It assumes that surface emissivity 
    !          errors have zero correlation with other state variables.
    !
    implicit none

    ! Arguments:
    real(8),                   intent(in)    :: BmatHiLand(:,:,:)       ! Static B-matrix over Land
    real(8),                   intent(in)    :: BmatHiSea(:,:,:)        ! Static B-matrix over Sea
    real(8),                   intent(inout) :: bMatEns(:,:,:)          ! Ensemble B-Matrix
    type(struct_obs),          intent(in)    :: obsSpaceData            ! ObsSpaceData object
    integer,                   intent(in)    :: validHeaderCount        ! Number of valid header
    integer,                   intent(in)    :: validHeaderIndex(:)     ! Valid header indexes
    character(len=4),          intent(in)    :: bmat1D_includeAnlVar(:) ! Analysis variables in B-Matrix
    type(struct_vco), pointer, intent(in)    :: vco_1DVarBmat           ! Vertical coordinate object of B-Matrix
    type(struct_vco), pointer, intent(in)    :: vco_1DVar               ! Vertical coordinate object of var1D
    real(4), allocatable,      intent(in)    :: latLandHi(:)            ! Latitude band for the B-Matrix over land
    real(4), allocatable,      intent(in)    :: latSeaHi(:)             ! Latitude band for the B-Matrix over sea
    
    ! Locals:
    integer :: columnIndex, headerIndex, varIndex, emissVarIndex
    real(8) :: latitude
    integer :: terrainType, landSea, surfaceType
    integer :: varListEmiss, nlevEmiss, latitudeBandIndex(1)
    
    ! Get number of levels in the surface emssivity control vector
    varListEmiss = vnl_varListIndexOther('EMMW')
    nlevEmiss = vco_1DVar%nlev_Other(varListEmiss)

    ! Determine the start varIndex for surface emissivity in the B-Matrix.
    emissVarIndex = 0 
    do varIndex = 1, size(bmat1D_includeAnlVar)
      select case(bmat1D_includeAnlVar(varIndex))
      case('TT')
        emissVarIndex = emissVarIndex + vco_1DVarBmat%nLev_T
      case('HU')
        emissVarIndex = emissVarIndex + vco_1DVarBmat%nLev_T
      case('UU','VV')
        emissVarIndex = emissVarIndex + vco_1DVarBmat%nlev_M
      case('TG')
        emissVarIndex = emissVarIndex + 1
      case('P0')
        emissVarIndex = emissVarIndex + 1
      case('EMMW')
        emissVarIndex = emissVarIndex + 1
      case default
        call utl_abort('sse_updateBEnsMatEmissFromBHi: unsupported variable ' // bmat1D_includeAnlVar(varIndex))
      end select
    end do
    
    ! Copy emissivity errors from static to ensemble B-Matrix
    do columnIndex = 1, validHeaderCount
      headerIndex = validHeaderIndex(columnIndex)
      latitude = obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex) !radian 

      terrainType = obs_headElem_i(obsSpaceData,OBS_TTYP,headerIndex)
      landSea = obs_headElem_i(obsSpaceData,OBS_STYP,headerIndex)

      if ( terrainType ==  0 ) then
        surfaceType  = surftype_seaice
      else
        surfaceType  = landSea
      end if
      if (surfaceType == 1) then !Sea
        latitudeBandIndex = minloc(abs(latitude - latSeaHi(:)))
        bMatEns(columnIndex, emissVarIndex:emissVarIndex + nlevEmiss - 1, emissVarIndex:emissVarIndex + nlevEmiss - 1) = &
            BmatHiSea(latitudeBandIndex(1), emissVarIndex:emissVarIndex + nlevEmiss -1, emissVarIndex:emissVarIndex + nlevEmiss -1)
      else
        latitudeBandIndex = minloc(abs(latitude - latLandHi(:)))
        bMatEns(columnIndex, emissVarIndex:emissVarIndex + nlevEmiss -1, emissVarIndex:emissVarIndex + nlevEmiss -1) = &
            BmatHiLand(latitudeBandIndex(1), emissVarIndex:emissVarIndex + nlevEmiss -1, emissVarIndex:emissVarIndex + nlevEmiss -1)
      end if
     end do

  end subroutine sse_updateBEnsMatEmissFromBHi
end module surfaceEmissivity_mod
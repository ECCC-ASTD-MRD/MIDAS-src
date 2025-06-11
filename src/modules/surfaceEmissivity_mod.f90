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

  ! Public procedures
  public :: sse_setupEmissivityfromState, sse_extractEmissivityCol
  public :: sse_updateBEnsMatEmissFromBHi, sse_emissivityRttovLimits
  public :: sse_readEmissError, sse_readCEmissMatrixByFileName, sse_emissErrMatSqrt
  public :: sse_emissFromObsSpace

contains

  !--------------------------------------------------------------------------
  ! sse_emissFromObsSpace
  !--------------------------------------------------------------------------
  subroutine sse_emissFromObsSpace(obsSpaceData, updatedEmissivity, bodyIndexFromBtIndex, chanprof, sensorHeaderIndexes)
    !
    !:Purpose: Reading emissivity from ObsSpaceData
    !
    implicit none

    ! Arguments:
    type(struct_obs),        intent(in)    :: obsSpaceData            ! ObsSpacedata object
    type(rttov_emissivity),  intent(inout) :: updatedEmissivity(:)    ! Update emissivity
    integer,                 intent(in)    :: bodyIndexFromBtIndex(:) ! Provides the bodyIndex in ObsSpaceData based on btIndex
    type(rttov_chanprof),    intent(in)    :: chanprof(:)             ! Profile and channel index object
    integer,                 intent(in)    :: sensorHeaderIndexes(:)  ! Sensor obsSpaceDate header Indexes

    ! Locals:
    integer :: btCount, profileCount
    integer :: profileIndex, btIndex
    integer :: bodyIndex, headerIndex
    real(8) :: emissObsSpace
    integer :: sfcType

    write(*,*) 'sse_emissFromObsSpace: Use emissivity from obsSpaceData'

    btCount = size(updatedEmissivity)
    profileCount = size(sensorHeaderIndexes)

    do profileIndex = 1, profileCount !loop on profiles

      headerIndex = sensorHeaderIndexes(profileIndex)

      do btIndex = 1, btCount !loop on channels
        if (chanprof(btIndex)%prof == profileIndex) then
          bodyIndex = bodyIndexFromBtIndex(btIndex)
          emissObsSpace = obs_bodyElem_r(obsspacedata, OBS_SEM, bodyIndex)
          sfcType = obs_headElem_i(obsspacedata, OBS_STYP, headerIndex)

          if (sfcType == 0 .and. emissObsSpace > 0.d0 .and. emissObsSpace <= 1.d0 ) then
            updatedEmissivity(btIndex)%emis_in = emissObsSpace
          end if
        end if
      end do
    end do

  end subroutine sse_emissFromObsSpace

  !--------------------------------------------------------------------------
  ! sse_readEmissError
  !--------------------------------------------------------------------------
  subroutine sse_readEmissError(chanList, emissStdErr, numChan)
    !
    !:Purpose: Reading emissivity error std for AMSU-A
    !
    implicit none

    ! Arguemnts:
    integer, allocatable, intent(out) :: chanList(:)    ! Channel List
    real, allocatable,    intent(out) :: emissStdErr(:) ! Surf. Emissivity Error Stdev.
    integer,              intent(out) :: numChan        ! Number of Channels

    ! Locals:
    integer, external      :: fnom, fclos
    integer                :: iunit, ierr
    character(len=20)      :: fileName
    character(len = 512)   :: tmpStr
    integer                :: chanIndex

    fileName = 'EmissErrStd.dat'
    iunit = 0
    ierr = fnom(iunit, fileName, 'SEQ+FMT', 0)
    if (ierr /= 0) call utl_abort('sse_readEmissError: Error reading surface emissivity error stdev. file')

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
                                   chanprof, sensorHeaderIndexes, &
                                   nSensor, sensorList, instrumentName, tovsMaxChanneNum, &
                                   tovsChannelOffset, tovsIChan, surftype, &
                                   columProfTl_opt, columProfAd_opt, emissivityProfDt_opt)
    !
    !:Purpose: Setup emissivity for the rttov direct, tangent linear and adjoint computation when emissivity is included
    !          as an analysis variable.  When profType = 'dt', the background surface surface emissivity state is updated
    !          in 'emissivity_inout' using 'emissivityProfDt_opt' as input. When profType = 'tl', the tangent linear surface
    !          emissivity is updated in 'emissivity_inout' from 'columProfTl_opt' as input. When profType = 'ad', the adjoint surface
    !          emissivity is updated in 'columProfAd_opt' from 'emissivity_inout' as input.
    !
    implicit none

    ! Arguments:
    type(rttov_emissivity), pointer,           intent(inout) :: emissivity_inout(:)         ! Updated Surface emissivity object
    type(struct_obs),                          intent(in)    :: obsSpaceData                ! ObsSpaceData object
    integer,                                   intent(in)    :: bodyIndexFromBtIndex(:)     ! Provides the bodyIndex in ObsSpaceData based on btIndex
    type(rttov_chanprof),                      intent(in)    :: chanprof(:)                 ! Profile and channel index object
    integer,                                   intent(in)    :: sensorHeaderIndexes(:)      ! Sensor obsSpaceData header Indexes
    integer,                                   intent(in)    :: nSensor                     ! Number of individual sensors.
    integer,                                   intent(in)    :: sensorList(:)               ! Sensor number for each profile (tvs_lsensor)
    character(len=15),                         intent(in)    :: instrumentName(:)           ! Satellite name (tvs_instrumentName)
    integer,                                   intent(in)    :: tovsMaxChanneNum            ! Max. value for channel number (tvs_maxChannelNumber)
    integer,                                   intent(in)    :: tovsChannelOffset(:)        ! RTTOV channel mapping offset (tvs_channelOffset)
    integer,                                   intent(in)    :: tovsIChan(:,:)              ! List of channels per instrument (tvs_ichan)
    integer,                                   intent(in)    :: surftype(:)                 ! Surface type of the profile (land/sea)
    type(struct_columnData), target, optional, intent(inout) :: columProfAd_opt             ! Column object for rttov adjoint computation
    type(struct_columnData), target, optional, intent(in)    :: columProfTl_opt             ! Column object for rttov tangent computation
    real(8), optional, target,                 intent(in)    :: emissivityProfDt_opt(:,:)   ! Column object for rttov direct computation

    ! Locals:
    integer                          :: profileIndex, btIndex, sensorIndex
    integer                          :: bodyIndex, headerIndex
    integer                          :: channelNumber, channelIndex
    integer                          :: btCount, numSourceFile
    real(8), pointer                 :: emissivityCol(:), emissivityTovsIndex(:,:)
    type(struct_columnData), pointer :: column
    character(len=2)                 :: profType

    write(*,*) 'sse_setupEmissivityfromState: STARTING'

    do sensorIndex = 1, nSensor
      if (trim(instrumentName(sensorIndex)) /= 'AMSUA') then
        write(*,*) 'Satellite Name: ', trim(instrumentName(sensorIndex))
        call utl_abort('sse_setupEmissivityfromState: !! can only be used for AMSUA instrument')
      end if
    end do

    numSourceFile = 0

    ! Setup for the tangent linear computation for surface emissivity
    if (present(columProfTl_opt)) then
      column => columProfTl_opt
      profType = 'tl'
      numSourceFile = numSourceFile + 1
    end if

    ! Setup for the adjoint computation for surface emissivity
    if (present(columProfAd_opt)) then
      column => columProfAd_opt
      profType = 'ad'
      numSourceFile = numSourceFile + 1
    end if

    ! Extracting surface emissivity from background trial
    if (present(emissivityProfDt_opt)) then
      emissivityTovsIndex => emissivityProfDt_opt
      profType = 'dt'
      numSourceFile = numSourceFile + 1
    end if

    if (numSourceFile /= 1) then
      call utl_abort('sse_setupEmissivityfromState: Only one source of emissivity variable must be supplied')
    end if

    btCount = size(emissivity_inout)
    do btIndex = 1, btCount
      bodyIndex = bodyIndexFromBtIndex(btIndex)
      profileIndex = chanprof(btIndex)%prof
      headerIndex = sensorHeaderIndexes(profileIndex)

      call sse_getChannelNumIndexFromPPP(obsSpaceData, headerIndex, bodyIndex, &
                                         channelNumber, channelIndex, &
                                         sensorList, tovsMaxChanneNum, tovsChannelOffset, &
                                         tovsIChan)

      if (profType == 'tl') then
        emissivityCol => col_getColumn(column, headerIndex,'EMMW')
        emissivity_inout(btIndex)%emis_in = emissivityCol(channelNumber)
      else if (profType == 'ad') then
        emissivityCol => col_getColumn(column, headerIndex,'EMMW')
        emissivityCol(channelNumber) = emissivity_inout(btIndex)%emis_out
      else if (profType == 'dt') then
        if (surftype(headerIndex) == surftype_land .and. &
           emissivityTovsIndex(headerIndex, channelNumber) > 0.d0 .and. &
           emissivityTovsIndex(headerIndex, channelNumber) <= 1.d0) then

          emissivity_inout(btIndex)%emis_in = emissivityTovsIndex(headerIndex, channelNumber)
        end if
      end if
    end do

  end subroutine sse_setupEmissivityfromState

  !--------------------------------------------------------------------------
  !  sse_extractEmissivityCol
  !--------------------------------------------------------------------------
  subroutine sse_extractEmissivityCol(column, emissivityFromTrl, profileCount, sensorHeaderIndexes, headerEnd)
    !
    !:Purpose: Extract emissivity from column objects
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)    :: column                   ! column structure
    real(8), allocatable,    intent(inout) :: emissivityFromTrl(:, :)  ! Extract emissivity
    integer,                 intent(in)    :: profileCount             ! Profile count
    integer,                 intent(in)    :: sensorHeaderIndexes(:)   ! Sensor header index
    integer,                 intent(in)    :: headerEnd                ! end headerIndex for TOVS observations

    ! Locals:
    integer :: emissTrlNumChan
    integer :: profileIndex, headerIndex
    real(8), pointer :: onecolumn(:)

    emissTrlNumChan = col_getNumLev(column, 'OT', varName_opt = 'EMMW')

    if (.not. allocated(emissivityFromTrl)) allocate(emissivityFromTrl(headerEnd, emissTrlNumChan))

    do profileIndex = 1 , profileCount
      headerIndex = sensorHeaderIndexes(profileIndex)
      oneColumn => col_getColumn(column, headerIndex, 'EMMW')
      emissivityFromTrl(headerIndex, :) = oneColumn(:)
    end do
  end subroutine sse_extractEmissivityCol

  !--------------------------------------------------------------------------
  !  sse_getChannelNumIndexFromPPP
  !--------------------------------------------------------------------------
  subroutine sse_getChannelNumIndexFromPPP(obsSpaceData, headerIndex, bodyIndex, &
                                           channelNumber, channelIndex, &
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
    integer,          intent(in)  :: sensorList(:)           ! Sensor number for each profile (tvs_lsensor)
    integer,          intent(in)  :: tovsMaxChanneNum        ! Max. value for channel number (tvs_maxChannelNumber)
    integer,          intent(in)  :: tovsChannelOffset(:)    ! RTTOV channel mapping offset (tvs_channelOffset)
    integer,          intent(in)  :: tovsIChan(:,:)          ! List of channels per instrument (local)
    integer,          intent(out) :: channelNumber           ! Channel Number
    integer,          intent(out) :: channelIndex            ! Channel Index

    ! Locals:
    integer :: sensorIndex

    sensorIndex = sensorList(headerIndex)

    channelNumber = nint(obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex))
    channelNumber = max(0, min(channelNumber, tovsMaxChanneNum + 1))
    channelNumber = channelNumber - tovsChannelOffset(sensorIndex)
    channelIndex = utl_findloc(tovsIChan(:, sensorIndex), channelNumber)

  end subroutine sse_getChannelNumIndexFromPPP

  !--------------------------------------------------------------------------
  !  sse_updateBEnsMatEmissFromBHi
  !--------------------------------------------------------------------------
  subroutine sse_updateBEnsMatEmissFromBHi(BmatHiLand, BmatHiSea, BmatEns, latLandHi, latSeaHi, obsSpaceData, &
                                           validHeaderIndex, validHeaderCount, bmat1D_includeAnlVar, column)
    !
    !:Purpose: Extract surface emissivity elements in the static B-Matrix and
    !          copy to the ensemble B-Matrix. It assumes that surface emissivity
    !          errors have zero correlation with other analysis variables.
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
    real(4), allocatable,      intent(in)    :: latLandHi(:)            ! Latitude band for the B-Matrix over land
    real(4), allocatable,      intent(in)    :: latSeaHi(:)             ! Latitude band for the B-Matrix over sea
    type(struct_columnData),   intent(in)    :: column                  ! Column object
    ! Locals:
    integer :: columnIndex, headerIndex, varIndex, emissVarIndex
    real(8) :: latitude
    integer :: terrainType, landSea, surfaceType
    integer :: nlevEmiss, latitudeBandIndex(1)

    if (.not. any(bmat1D_includeAnlVar(:) == 'EMMW')) then
      call utl_abort('sse_updateBEnsMatEmissFromBHi: bmat1D_includeAnlVar does not include EMMW')
    end if

    ! Get number of levels in the surface emssivity analysis variable
    nlevEmiss = col_getNumLev(column, 'OT', varname_opt = 'EMMW')

    ! Determine the start varIndex for surface emissivity in the B-Matrix.
    emissVarIndex = 0
    do varIndex = 1, size(bmat1D_includeAnlVar)
      emissVarIndex = emissVarIndex + &
          col_getNumLev(column, vnl_varLevelFromVarname(bmat1D_includeAnlVar(varIndex)), varname_opt = trim(bmat1D_includeAnlVar(varIndex)))
    end do

    ! Copy emissivity errors from static to ensemble B-Matrix
    do columnIndex = 1, validHeaderCount
      headerIndex = validHeaderIndex(columnIndex)
      latitude = obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex) !radian

      terrainType = obs_headElem_i(obsSpaceData,OBS_TTYP,headerIndex)
      landSea = obs_headElem_i(obsSpaceData,OBS_STYP,headerIndex)

      if (terrainType ==  0) then
        surfaceType  = surftype_seaice
      else
        surfaceType  = landSea
      end if

      if (surfaceType == 1) then !Sea
        latitudeBandIndex = minloc(abs(latitude - latSeaHi(:)))
        bMatEns(columnIndex,  emissVarIndex - nlevEmiss + 1: emissVarIndex, emissVarIndex - nlevEmiss + 1: emissVarIndex) = &
            BmatHiSea(latitudeBandIndex(1), emissVarIndex - nlevEmiss + 1: emissVarIndex, emissVarIndex - nlevEmiss + 1: emissVarIndex)
      else
        latitudeBandIndex = minloc(abs(latitude - latLandHi(:)))
        bMatEns(columnIndex,  emissVarIndex - nlevEmiss + 1: emissVarIndex, emissVarIndex - nlevEmiss + 1: emissVarIndex) = &
            BmatHiLand(latitudeBandIndex(1), emissVarIndex - nlevEmiss + 1: emissVarIndex, emissVarIndex - nlevEmiss + 1: emissVarIndex)
      end if
    end do

  end subroutine sse_updateBEnsMatEmissFromBHi

  !--------------------------------------------------------------------------
  ! sse_emissivityRttovLimits
  !--------------------------------------------------------------------------
  subroutine sse_emissivityRttovLimits(column, columnRef_opt)
    !
    !:Purpose: Restrict simulated background surface emissivity within physical
    !          values, between 0 and 1. Fill column%all(:,:) with missing emissivity
    !          values based on a reference column object, if present.
    !
    implicit none

    ! Arguments:
    type(struct_columnData),           intent(inout)  :: column         ! Perturbed column object
    type(struct_columnData), optional, intent(in)     :: columnRef_opt  ! Reference column object (non-perturbed)

    ! Locals:
    real(8), pointer     :: emissPtr(:,:), emissPtrRef(:,:)
    integer              :: numCol, numLev, numColRef, numLevRef
    integer              :: icol, ilev
    real(8), parameter   :: missingValue = -1.0d0
    real(8)              :: emissDiffTmp

    numCol = col_getNumCol(column)
    numLev = col_getNumLev(column,'OT', varname_opt = 'EMMW')
    emissPtr => col_getAllColumns(column, 'EMMW')

    if (present(columnRef_opt)) then
      numColRef = col_getNumCol(columnRef_opt)
      numLevRef = col_getNumLev(columnRef_opt, 'OT', varname_opt = 'EMMW')
      emissPtrRef => col_getAllColumns(columnRef_opt, 'EMMW')

      if (numCol /= numColRef .and. numLev /= numLevRef) then
        call utl_abort('sse_emissivityRttovLimits: number of column and ' // &
                        'levels between column and reference column are not the same')
      end if
    end if

    colLoop: do icol = 1, numCol
      if (present(columnRef_opt)) then
        ! Fill missing value if ref column also have missing value
        do ilev = 1, size(emissPtrRef, 1)
          if ( utl_isEqual(emissPtrRef(ilev, icol), missingValue) ) then
            emissPtr(:, icol) = missingValue
            cycle colLoop
          end if
        end do

        ! Fill missing value if ref column have negative emissivity
        if(any(emissPtrRef(:, icol) < 0.0d0)) then
          emissPtr(:, icol) = missingValue
        ! Limit simulated emissivity to zero if ref column emissivity is equal to zero or greater
        else if(any(emissPtrRef(:, icol) >= 0.0d0) .and. any(emissPtr(:, icol) < 0.0d0)) then
          emissDiffTmp = minval(emissPtr(:, icol)) - 0.0d0
          emissPtr(:, icol) = emissPtr(:, icol) - emissDiffTmp*1.01
        ! Limit simulated emissivity to one
        else if(any(emissPtrRef(:, icol) <= 1.0d0) .and. any(emissPtr(:, icol) > 1.0d0)) then
          emissDiffTmp = maxval(emissPtr(:, icol)) - 1.0d0
          emissPtr(:, icol) = emissPtr(:, icol) - emissDiffTmp*1.01
        end if
      else
        ! Limit simulated emissivity to one
        if(any(emissPtr(:, icol) < 0.0d0)) then
          emissDiffTmp = minval(emissPtr(:, icol)) - 0.0d0
          emissPtr(:, icol) = emissPtr(:, icol) - emissDiffTmp * 1.01
        ! Limit simulated emissivity to zero
        else if(any(emissPtr(:, icol) > 1.0d0)) then
          emissDiffTmp = maxval(emissPtr(:, icol)) - 1.0d0
          emissPtr(:, icol) = emissPtr(:, icol) - emissDiffTmp * 1.01
        end if
      end if
    end do colLoop

  end subroutine sse_emissivityRttovLimits

end module surfaceEmissivity_mod

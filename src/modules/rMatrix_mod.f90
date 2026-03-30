
module rMatrix_mod
  ! MODULE rMatrix_mod (prefix='rmat' category='2. B and R matrices')
  !
  !:Purpose:  Module to handle non-diagonal observation-error covariance
  !           matrices for assimilation of radiances.
  !
  use rmn_fnom
  use linearAlgebra_mod
  use midasMpi_mod
  use rttovInterfaces_mod
  use rttov_const, only: errorstatus_success
  use utilities_mod
  use runtimeInfo_mod
  use obsSpaceData_mod
  use tovs_mod
  use mathPhysConstants_mod

  implicit none
  private
  save

  ! Public variables
  ! This is a namelist variable
  logical, public, protected :: rmat_lnondiagr ! choose to use non-diagonal R matrix (i.e. non-zero correlations)

  ! Public subroutines
  public :: rmat_init, rmat_cleanup, rmat_readCMatrix, rmat_RsqrtInverseOneObs, rmat_RsqrtInverseAllObs, rmat_Rsqrt
  public :: rmat_updateRmat, rmat_writeRCorrFile, rmat_getRmatrix

  type rmat_matrix
    real(8), pointer :: Rmat(:,:) => null()
    integer, pointer :: listChans(:) => null()
    integer          :: nchans=0
  end type rmat_matrix

  type(rmat_matrix), target, allocatable :: Rcorr_inst(:) ! non diagonal Correlation matrices for each instrument
  type(rmat_matrix), target, allocatable :: R_tovs(:)     ! non diagonal R matrices used for the assimilation of all radiances

  ! Private namelist module variables
  real(8) :: rmat_estLatMax      ! Max latitude criteria for obs to be used in estimating R-Matrix
  real(8) :: rmat_estLatMin      ! Min latitude criteria for obs to be used in estimating R-Matrix
  integer :: rmat_estLandSeaExcl ! Land/Sea criteria to exclude obs that are used in estimating R-Matrix(0: land, 1: sea)
  real(8) :: rmat_estElevMax     ! Max sfc elevation criteria to exclude obs to be used in estimating R-Matrix (km)

contains

  subroutine rmat_init(nsensors, headerEnd)

    implicit none

    ! Arguments:
    integer, intent(in) :: nsensors
    integer, intent(in) :: headerEnd

    ! Locals:
    integer :: ierr
    namelist /NAMRMAT/rmat_lnonDiagR, rmat_estLatMax, rmat_estLatMin, rmat_estLandSeaExcl
    namelist /NAMRMAT/rmat_estElevMax

    ! Default value for parameter rmat_lnondiagr, don't use interchannel correlation by default
    rmat_lnonDiagR = .false.
    rmat_estLatMax = mpc_missingvalue_r8
    rmat_estLatMin = mpc_missingvalue_r8
    rmat_estLandSeaExcl = mpc_missingvalue_int
    rmat_estElevMax = mpc_missingvalue_r8

    ! Read the parameters from NAMRMAT
    call rti_tmg_start(181,'low-level--readNML')
    read(utl_flnml,nml=namrmat, iostat=ierr)
    if (ierr /= 0) call rti_abort('rmat_init: Error reading namelist')
    if (mmpi_myid == 0) write(*,nml=namrmat)
    call rti_tmg_stop(181)
    if (rmat_lnonDiagR) then
      allocate(Rcorr_inst(nsensors))
      allocate(R_tovs(headerEnd))
    end if

  end subroutine rmat_init


  subroutine rmat_cleanup()
    implicit none

    if (rmat_lnondiagr) then
      deallocate(Rcorr_inst)
      deallocate(R_tovs)
    end if

  end subroutine rmat_cleanup


  subroutine rmat_readCMatrix(instrument, sensor_id, ichan )
    implicit none

    ! Arguments:
    integer, intent(in) :: instrument(3)
    integer, intent(in) :: sensor_id
    integer, intent(in) :: ichan(:)

    ! Locals:
    character (len=64) :: filename
    integer :: err

    call rttov_coeffname (err, instrument, coeffname=filename, filetype="Cmat")

    if (err == errorstatus_success) then
       filename = trim(filename) // ".dat"
       call rmat_readCMatrixByFileName(filename,Rcorr_inst(sensor_id), ichan )
    else
      write(*,*) "Unknown instrument ",instrument(:)
      call rti_abort("rmat_read_C_matrix")
    end if

  end subroutine rmat_readCMatrix


  subroutine rmat_readCMatrixByFileName(infile,C,chanList_opt)
    implicit none

    ! Arguments:
    character (len=*), intent(in)    :: infile          ! name of input file
    type(rmat_matrix), intent(inout) :: C               ! correlation matrix structure
    integer, optional, intent(in)    :: chanList_opt(:) ! list of requested channels (if missing will read all file content)

    ! Locals:
    integer :: i, j, iu, ierr, count, ich, nchn, nch
    real(8) :: x
    integer, allocatable :: foundChanIndex(:)

    nchn = -1
    if (present(chanList_opt)) then
      nchn = size(chanList_opt)
    end if

    iu = 0
    ierr = fnom(iu,trim(infile),'FTN+SEQ+R/O',0)
    if (ierr /= 0) then
      write(*,*) "Cannot open " // trim(infile)
      call rti_abort("rmat_readCMatrixByFileName")
    end if

    write(*,*) "rmat_readCMatrixByFileName: Reading " // trim(infile)

    read(iu,*) nch
    if (nchn == -1) then
      nchn = nch
    else
      if(nchn > nch) then
        write(*,*) "Not enough channels in " // trim(infile)
        write(*,*) nchn, nch
        call rti_abort("rmat_readCMatrixByFileName")
      end if
    end if
    allocate(foundChanIndex(nch))

    C%nchans = nchn
    allocate(C%Rmat(nchn,nchn))
    allocate(C%listChans(nchn))
    C%Rmat(:,:) = 0.d0
    do i = 1,nchn
      C%Rmat(i,i) = 1.d0
    end do
    count = 0
    foundChanIndex(:) = -1
    do i = 1, nch
      read(iu,*) ich
      if (present(chanList_opt)) then
        bj:do j = 1, nchn
          if (ich == chanList_opt(j)) then
            count = count + 1
            foundChanIndex(i) = j
            C%listChans(count) = ich
            exit bj
          end if
        end do bj
      else
        foundChanIndex(i) = i
        C%listChans(i) = ich
        count = count + 1
      end if
    end do
    if (count /= nchn) then
      write(*,*) "Warning: Missing information in " // trim(infile)
      do j = 1, nchn
        write(*,*) j, chanList_opt(j)
      end do
      write(*,*) "Not important if there is no observation of this family"
    end if

    do
      read(iu,*,iostat=ierr) i, j, x
      if (ierr /= 0) exit
      if (foundChanIndex(i) /= -1 .and. foundChanIndex(j) /= -1) then
        C%Rmat(foundChanIndex(i),foundChanIndex(j)) = x
        C%Rmat(foundChanIndex(j),foundChanIndex(i)) = x
      end if
    end do

    ierr= fclos(iu)
    deallocate(foundChanIndex)

  end subroutine rmat_readCMatrixByFileName


  subroutine rmat_RsqrtInverseOneObs(sensor_id,nsubset,obsIn,obsOut,list_sub,list_oer,headerIndex)
    !
    ! :Purpose: Apply the operator R**-1/2 to obsIn
    !           result in obsOut for the subset of channels specified
    !           in list_sub
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: sensor_id
    integer, intent(in)  :: nsubset
    integer, intent(in)  :: list_sub(nsubset)
    real(8), intent(in)  :: list_oer(nsubset)
    real(8), intent(in)  :: obsIn(nsubset)
    real(8), intent(out) :: obsOut(nsubset)
    integer, intent(in)  :: headerIndex

    ! Locals:
    real (8) :: Rsub(nsubset,nsubset), alpha, beta, variance
    integer :: foundChanIndex(nsubset)
    integer :: i,j

    if (R_tovs(headerIndex)%nchans == 0) then

      if (sensor_id <= 0 .or. sensor_id > size(Rcorr_inst)) then
        write(*,*) "invalid sensor_id",sensor_id,size(Rcorr_inst)
        call rti_abort('rmat_RsqrtInverseOneObs')
      end if

      foundChanIndex(:) = -1
      do i=1,nsubset
        bj: do j = 1, Rcorr_inst(sensor_id)%nchans
          if (list_sub(i) == Rcorr_inst(sensor_id)%listChans(j)) then
            foundChanIndex(i) = j
            exit bj
          end if
        end do bj
      end do
      if (any(foundChanIndex == -1)) then
        write(*,*) "Missing information for some channel !"
        write(*,*) list_sub(:)
        write(*,*) foundChanIndex(:)
        call rti_abort('rmat_RsqrtInverseOneObs')
      end if
      R_tovs(headerIndex)%nchans = nsubset
      allocate(R_tovs(headerIndex)%listChans(nsubset))
      R_tovs(headerIndex)%listChans(1:nsubset) = list_sub(1:nsubset)
      do j=1,nsubset
        do i=1,nsubset
          variance = list_oer(i) * list_oer(j)
          Rsub(i,j) = variance * Rcorr_inst(sensor_id)%Rmat(foundChanIndex(i),foundChanIndex(j))
        end do
      end do
      ! Calculation of R**-1/2
      call rti_tmg_start(20,'----RmatMatSqrt')
      call linalg_matSqrt(Rsub,nsubset,-1.d0,.false.)
      call rti_tmg_stop(20)
      allocate(R_tovs(headerIndex)%Rmat(nsubset,nsubset))
      do j=1,nsubset
        do i=1,nsubset
          R_tovs(headerIndex)%Rmat(i,j) = Rsub(i,j)
        end do
      end do
    end if

    call rti_tmg_start(21,'----RmatMatMult')
    alpha = 1.d0
    beta = 0.d0
    obsOut = 0.d0
    ! Optimized symetric matrix vector product from Blas
    call dsymv("L", nsubset, alpha, R_tovs(headerIndex)%Rmat, nsubset, obsIn, 1, beta, obsOut, 1)
    call rti_tmg_stop(21)

  end subroutine rmat_RsqrtInverseOneObs


  !--------------------------------------------------------------------------
  ! rmat_RsqrtInverseAllObs
  !--------------------------------------------------------------------------
  subroutine rmat_RsqrtInverseAllObs( obsSpaceData, elem_dest_i, elem_src_i )
    !
    !:Purpose: To apply observation-error variances to ROBDATA8(k_src,*) and to
    !          store it in the elem_src_s of obsspacedata
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsspacedata
    integer,          intent(in)    :: elem_dest_i ! destination index
    integer,          intent(in)    :: elem_src_i  ! source index

    ! Locals:
    integer :: bodyIndex, headerIndex
    integer :: idata, idatend, idatyp, count, channelNumber, channelIndex
    real(8) :: obsIn( tvs_maxChannelNumber ), obsOut( tvs_maxChannelNumber )
    integer :: list_chan( tvs_maxChannelNumber )
    real(8) :: list_OER( tvs_maxChannelNumber )

    ! NOTE I tried using openMP on this loop, but it increased the cost from 4sec to 80sec!!!
    do headerIndex =1, obs_numHeader(obsspacedata)

      idata   = obs_headElem_i( obsspacedata, OBS_RLN, headerIndex )
      idatend = obs_headElem_i( obsspacedata, OBS_NLV, headerIndex ) + idata - 1
      idatyp  = obs_headElem_i( obsspacedata, OBS_ITY, headerIndex )

      if ( tvs_isIdBurpTovs(idatyp) .and. rmat_lnondiagr ) then

        count = 0

        do bodyIndex = idata, idatend

          if (obs_bodyElem_i( obsspacedata, OBS_ASS, bodyIndex ) == obs_assimilated ) then
            call tvs_getChannelNumIndexFromPPP( obsSpaceData, headerIndex, bodyIndex, &
                                                channelNumber, channelIndex )
            count = count + 1
            list_chan( count ) = channelNumber
            list_OER( count ) = obs_bodyElem_r( obsspacedata, OBS_OER, bodyIndex )
            obsIn( count ) = obs_bodyElem_r( obsspacedata, elem_src_i, bodyIndex )
          end if

        end do

        if (count > 0) then

          call rmat_RsqrtInverseOneObs(tvs_lsensor(headerIndex), count, obsIn(1:count), obsOut(1:count), list_chan(1:count), list_OER(1:count), headerIndex)

          count = 0
          do bodyIndex = idata, idatend
            if ( obs_bodyElem_i( obsspacedata, OBS_ASS, bodyIndex ) == obs_assimilated) then
              count = count + 1
              call obs_bodySet_r(obsspacedata, elem_dest_i, bodyIndex,obsOut(count))
            end if
          end do

        else

          do bodyIndex = idata, idatend
            call obs_bodySet_r(obsspacedata, elem_dest_i, bodyIndex, 0.d0)
          end do

        end if

      else

        do bodyIndex = idata, idatend
          if (obs_bodyElem_i( obsspacedata, OBS_ASS, bodyIndex ) == obs_assimilated) then
            call obs_bodySet_r( obsspacedata, elem_dest_i, bodyIndex, &
                 obs_bodyElem_r( obsspacedata, elem_src_i, bodyIndex) / obs_bodyElem_r( obsspacedata, OBS_OER, bodyIndex ))
          end if
        end do

      end if ! is it a radiance in non diagonal R mode ?

    end do !loop on header

  end subroutine rmat_RsqrtInverseAllObs

  !--------------------------------------------------------------------------
  ! rmat_Rsqrt
  !--------------------------------------------------------------------------
  subroutine rmat_Rsqrt(sensor_id, nsubset, obsIn, obsOut, list_sub, list_oer)
    !
    ! :Purpose: Apply the operator R**1/2 to obsIn
    !           result in obsOut for the subset of channels specified
    !           in list_sub
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: sensor_id          ! Sensor ID
    integer, intent(in)  :: nsubset            ! Number of subset channels in R-matrix
    real(8), intent(in)  :: obsIn(nsubset)     ! Sampling Perturbation
    real(8), intent(out) :: obsOut(nsubset)    ! Error Perturbation
    integer, intent(in)  :: list_sub(nsubset)  ! List of subset channels in R-matrix
    real(8), intent(in)  :: list_oer(nsubset)  ! List of Obs Error

    ! Locals:
    real(8)              :: alpha, beta, variance
    real(8), allocatable :: Rsub(:,:)
    integer              :: foundChanIndex(nsubset)
    integer              :: chanIndex1, chanIndex2

    allocate(Rsub(nsubset, nsubset))
    if (sensor_id <= 0 .or. sensor_id > size(Rcorr_inst)) then
      write(*,*) 'invalid sensor_id', sensor_id,size(Rcorr_inst)
      call rti_abort('rmat_Rsqrt')
    end if

    ! Check error correlation for all channels are available
    foundChanIndex(:) = -1
    do chanIndex1 = 1, nsubset
      chanLoop: do chanIndex2 = 1, Rcorr_inst(sensor_id)%nchans
        if (list_sub(chanIndex1) == Rcorr_inst(sensor_id)%listChans(chanIndex2)) then
          foundChanIndex(chanIndex1) = chanIndex2
          exit chanLoop
        end if
      end do chanLoop
    end do

    if (any(foundChanIndex == -1)) then
      write(*,*) 'Missing information for some channel !'
      write(*,*) list_sub(:)
      write(*,*) foundChanIndex(:)
      call rti_abort('rmat_Rsqrt')
    end if

    do chanIndex2 = 1, nsubset
      do chanIndex1 = 1, nsubset
        variance = list_oer(chanIndex1) * list_oer(chanIndex2)
        Rsub(chanIndex1, chanIndex2) = variance * Rcorr_inst(sensor_id)%Rmat(foundChanIndex(chanIndex1), foundChanIndex(chanIndex2))
      end do
    end do

    ! Calculation of R**1/2
    call linalg_matSqrt(Rsub, nsubset, 1.d0, .false.)

    alpha = 1.d0
    beta = 0.d0
    obsOut = 0.d0

    ! Optimized symetric matrix vector product from Lapack
    call dsymv('L', nsubset, alpha, Rsub, nsubset, obsIn, 1, beta, obsOut, 1)

    deallocate(Rsub)

  end subroutine rmat_Rsqrt


  !--------------------------------------------------------------------------
  ! rmat_getRmatrix
  !--------------------------------------------------------------------------
  subroutine rmat_getRmatrix(sensor_id, list_sub, list_oer, Rsub)
    !
    ! :Purpose: get R for a given instrument, channel and standard deviation list
    !           result in Rsub
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: sensor_id            ! Sensor ID
    integer, intent(in)  :: list_sub(:)          ! List of subset channels in R-matrix
    real(8), intent(in)  :: list_oer(:)          ! List of Obs Error
    real(8), intent(out) :: Rsub(:,:)            ! output R matrix

    ! Locals:
    real(8)   :: variance
    integer   :: foundChanIndex(size(list_sub))
    integer   :: chanIndex1, chanIndex2, nSubset

    nSubset = size(list_sub)
    if (sensor_id <= 0 .or. sensor_id > size(Rcorr_inst)) then
      write(*,*) 'invalid sensor_id', sensor_id,size(Rcorr_inst)
      call rti_abort('rmat_getRMatrix')
    end if

    ! Check error correlation for all channels are available
    foundChanIndex(:) = -1
    do chanIndex1 = 1, nSubset
      chanLoop: do chanIndex2 = 1, Rcorr_inst(sensor_id)%nchans
        if (list_sub(chanIndex1) == Rcorr_inst(sensor_id)%listChans(chanIndex2)) then
          foundChanIndex(chanIndex1) = chanIndex2
          exit chanLoop
        end if
      end do chanLoop
    end do

    if (any(foundChanIndex == -1)) then
      write(*,*) 'Missing information for some channel !'
      write(*,*) list_sub(:)
      write(*,*) foundChanIndex(:)
      call rti_abort('rmat_getRMatrix')
    end if

    do chanIndex2 = 1, nSubset
      do chanIndex1 = 1, nSubset
        variance = list_oer(chanIndex1) * list_oer(chanIndex2)
        Rsub(chanIndex1, chanIndex2) = variance * Rcorr_inst(sensor_id)%Rmat(foundChanIndex(chanIndex1), foundChanIndex(chanIndex2))
      end do
    end do

  end subroutine rmat_getRmatrix


  !--------------------------------------------------------------------------
  ! rmat_updateRmat
  !--------------------------------------------------------------------------
  subroutine rmat_updateRmat(obsSpaceData, obsSpaceIndexObs, obsSpaceIndexTrue)
    !
    !:Purpose: Estimate observation error std and their correlations
    !          based on truth and simulated obs
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData      ! ObsSpaceData object
    integer,          intent(in)    :: obsSpaceIndexObs  ! ObsSpace Index corresponding to observation
    integer,          intent(in)    :: obsSpaceIndexTrue ! ObsSpace Index corresponding to truth


    ! Locals:
    integer              :: headerIndex, bodyIndex, sensorIndex, taskIndex, chanIndex
    integer              :: channelNumber, channelIndex
    integer              :: idata, idatend, idatyp
    integer              :: assimChan, tag
    integer, allocatable :: headerCount(:), headerCountMpiGlobal(:), headerCountAllTasks(:,:)
    real(8), allocatable :: obsErrSum(:,:), obsErrSumAllTasks(:,:,:), meanObsErrMpiGlobal(:,:)
    real(8), allocatable :: obsErrSqrdSum(:,:)
    real(8), allocatable :: vector(:,:), localRmat(:,:)
    real(8), allocatable :: obsErrStdev(:,:)
    real(8)              :: obsErr, tmpObsErr
    integer              :: numchan, ichan, jchan
    integer, allocatable :: chanList(:,:), chanListAllTasks(:,:,:)

    type(rmat_matrix), target, allocatable  :: estR(:)
    type(rmat_matrix), target, allocatable  :: ObsErrSqrdMat(:)
    type(rmat_matrix), target, allocatable :: RCorr(:)

    allocate(headerCountMpiGlobal(tvs_nsensors))
    allocate(obsErrSum(tvs_nsensors, maxval(tvs_nchanMpiGlobal)))
    allocate(obsErrSumAllTasks(tvs_nsensors, maxval(tvs_nchanMpiGlobal), mmpi_nprocs))
    allocate(obsErrSqrdSum(tvs_nsensors, maxval(tvs_nchanMpiGlobal)))
    allocate(headerCount(tvs_nsensors))
    allocate(headerCountAllTasks(tvs_nsensors, mmpi_nprocs))
    allocate(meanObsErrMpiGlobal(tvs_nsensors, maxval(tvs_nchanMpiGlobal)))
    allocate(chanList(tvs_nsensors, maxval(tvs_nchanMpiGlobal)))
    allocate(chanListAllTasks(tvs_nsensors, maxval(tvs_nchanMpiGlobal), mmpi_nprocs))

    allocate(estR(tvs_nsensors))
    allocate(ObsErrSqrdMat(tvs_nsensors))

    do sensorIndex = 1, tvs_nsensors
      ! Allocate estimate R-Matrix
      estR(sensorIndex)%nchans = tvs_nchanMpiGlobal(sensorIndex)
      allocate(estR(sensorIndex)%Rmat(tvs_nchanMpiGlobal(sensorIndex), tvs_nchanMpiGlobal(sensorIndex)))
      allocate(estR(sensorIndex)%listChans(tvs_nchanMpiGlobal(sensorIndex)))
      estR(sensorIndex)%Rmat(:, :) = 0.0d0

      ObsErrSqrdMat(sensorIndex)%nchans = tvs_nchanMpiGlobal(sensorIndex)
      allocate(ObsErrSqrdMat(sensorIndex)%Rmat(tvs_nchanMpiGlobal(sensorIndex), tvs_nchanMpiGlobal(sensorIndex)))
      allocate(ObsErrSqrdMat(sensorIndex)%listChans(tvs_nchanMpiGlobal(sensorIndex)))
      ObsErrSqrdMat(sensorIndex)%Rmat(:, :) = 0.0d0
    end do

    ! Initialize values
    meanObsErrMpiGlobal(:,:) = 0.0d0
    headerCount(:) = 0
    headerCountMpiGlobal(:) = 0
    obsErrSum(:,:) = 0.0d0
    chanList(:,:) = MPC_missingValue_INT
    chanListAllTasks(:,:,:) = MPC_missingValue_INT

    ! Count the number of header profile that will be used to compute R-Matrix.
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER1: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER1

      idatyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if (.not. tvs_isIdBurpTovs(idatyp)) then
        write(*,*) 'rmat_updateRmat: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
        cycle HEADER1
      end if

      sensorIndex = tvs_lsensor(headerIndex)

      ! Exclude observations with land sea mask
      if (rmat_estLandSeaExcl /= mpc_missingvalue_int .and. &
          obs_headElem_i(obsSpaceData, OBS_STYP, headerIndex) == rmat_estLandSeaExcl) then
        cycle HEADER1
      end if

      ! Exclude observations with elevation greater than rmat_estElevMax
      if (.not. utl_isEqual(rmat_estElevMax, mpc_missingvalue_r8) .and. &
          obs_headElem_r(obsspacedata, OBS_ELEV, headerIndex) > rmat_estElevMax) then
        cycle HEADER1
      end if

      ! Max latitude criteria
      if (.not. utl_isEqual(rmat_estLatMax, mpc_missingvalue_r8) .and. &
          obs_headElem_r(obsspacedata, OBS_LAT, headerIndex) * MPC_DEGREES_PER_RADIAN_R8 > rmat_estLatMax) then
        cycle HEADER1
      end if

      ! Min latitude criteria
      if (.not. utl_isEqual(rmat_estLatMin, mpc_missingvalue_r8) .and. &
          obs_headElem_r(obsspacedata, OBS_LAT, headerIndex) * MPC_DEGREES_PER_RADIAN_R8 < rmat_estLatMin) then
        cycle HEADER1
      end if

      idata   = obs_headElem_i(obsspacedata, OBS_RLN, headerIndex)
      idatend = obs_headElem_i(obsspacedata, OBS_NLV, headerIndex) + idata - 1

      if (tvs_isIdBurpTovs(idatyp)) then
        ! Check if all required channels are flag as assimilated in order to be used to estimate R-Matrix.
        assimChan = 0
        do bodyIndex = idata, idatend
          if (obs_bodyElem_i( obsspacedata, OBS_ASS, bodyIndex) == obs_assimilated) assimChan = assimChan + 1
        end do
        if (assimChan /= tvs_nchanMpiGlobal(sensorIndex)) cycle HEADER1

        allocate(vector(1,tvs_nchanMpiGlobal(sensorIndex)))

        headerCount(sensorIndex) = headerCount(sensorIndex) + 1
        assimChan = 0
        do bodyIndex = idata, idatend
          if (obs_bodyElem_i( obsspacedata, OBS_ASS, bodyIndex) == obs_assimilated) then
            assimChan = assimChan + 1

            ! Append channel list into estR object
            if (headerCount(sensorIndex) == 1) then
              call tvs_getChannelNumIndexFromPPP(obsSpaceData, headerIndex, bodyIndex, &
                                                channelNumber, channelIndex)

              chanList(sensorIndex, assimChan) = channelNumber
            end if

            ! Compute the observation error based on difference between observation and truth
            obsErr = obs_bodyElem_r(obsspacedata, obsSpaceIndexObs, bodyIndex) - &
                    obs_bodyElem_r(obsspacedata, obsSpaceIndexTrue, bodyIndex)
            vector(1, assimChan) = obsErr

            ! Compute the sum of observation sum(x)
            obsErrSum(sensorIndex, assimChan) = obsErrSum(sensorIndex, assimChan) + obsErr

          end if
        end do

        ! Compute sum of the squared observation error matrix sum(x^2)
        ObsErrSqrdMat(sensorIndex)%Rmat(:,:) = ObsErrSqrdMat(sensorIndex)%Rmat(:,:) + &
                                               matmul(transpose(vector), vector)

        deallocate(vector)
      end if
    end do HEADER1

    ! Gather all count of profiles, and 'sum of observations errors' from all mpi tasks
    call mmpi_allGather(headerCount, headerCountAllTasks)
    call mmpi_allGather(obsErrSum, obsErrSumAllTasks)
    call mmpi_allGather(chanList, chanListAllTasks)
    deallocate(headerCount)

    Sensor: do sensorIndex = 1, tvs_nsensors
      do taskIndex = 1, mmpi_nprocs
        if (headerCountAllTasks(sensorIndex, taskIndex) > 0) then
          estR(sensorIndex)%listChans(:) = chanListAllTasks(sensorIndex, 1:tvs_nchanMpiGlobal(sensorIndex), taskIndex)
          cycle Sensor
        end if
      end do
    end do Sensor
    deallocate(chanListAllTasks)
    deallocate(chanList)

    ! Compute the mean of observation error and total profile count to all MPI tasks
    do sensorIndex = 1, tvs_nsensors
      headerCountMpiGlobal(sensorIndex) = sum(headerCountAllTasks(sensorIndex,:))
      do chanIndex = 1, tvs_nchanMpiGlobal(sensorIndex)
        meanObsErrMpiGlobal(sensorIndex, chanIndex) = sum(obsErrSumAllTasks(sensorIndex, chanIndex, :)) / &
                                                      (headerCountMpiGlobal(sensorIndex) -1)
      end do
    end do

    deallocate(obsErrSum)
    deallocate(obsErrSumAllTasks)

    ! Compute the R Matrix.
    do sensorIndex = 1, tvs_nsensors

      ! Send ObsErrSqrdMat from all MPI tasks to MPI task 0.
      if (mmpi_myid > 0) then
        tag = mmpi_myId
        call mmpi_send(ObsErrSqrdMat(sensorIndex)%Rmat(:,:), tag, procID = 0)
      end if

      if (mmpi_myid == 0) then

        ! Collect ObsErrSqrdMat from all MPI tasks
        do taskIndex = 1, mmpi_nprocs - 1
          allocate(localRmat(tvs_nchanMpiGlobal(sensorIndex), tvs_nchanMpiGlobal(sensorIndex)))
          tag = taskIndex
          call mmpi_recv(localRmat, tag, taskIndex)
          ObsErrSqrdMat(sensorIndex)%Rmat(:,:) = ObsErrSqrdMat(sensorIndex)%Rmat(:,:) + localRmat
          deallocate(localRmat)
        end do

        !Compute Rmatix = mean(x^2) + (mean(x))^2
        allocate(vector(1,tvs_nchanMpiGlobal(sensorIndex)))
        vector(1,:) = meanObsErrMpiGlobal(sensorIndex, 1:tvs_nchanMpiGlobal(sensorIndex))
        estR(sensorIndex)%Rmat(:,:) = ObsErrSqrdMat(sensorIndex)%Rmat(:,:) / (headerCountMpiGlobal(sensorIndex) -1) + &
                                      matmul(transpose(vector), vector)
        deallocate(vector)
      end if

      call mmpi_bcast(estR(sensorIndex)%Rmat)
    end do
    deallocate(ObsErrSqrdMat)
    deallocate(meanObsErrMpiGlobal)
    deallocate(headerCountAllTasks)
    deallocate(headerCountMpiGlobal)

    ! Allocate estimate R-Matrix
    allocate(RCorr(tvs_nsensors))
    do sensorIndex = 1, tvs_nsensors
      RCorr(sensorIndex)%nchans = tvs_nchanMpiGlobal(sensorIndex)
      allocate(RCorr(sensorIndex)%Rmat(tvs_nchanMpiGlobal(sensorIndex), tvs_nchanMpiGlobal(sensorIndex)))
      allocate(RCorr(sensorIndex)%listChans(tvs_nchanMpiGlobal(sensorIndex)))
    end do

    allocate(obsErrStdev(tvs_nsensors, maxval(tvs_nchanMpiGlobal)))

    ! Extract Sigma(o) and observation error correlation matrix from estimated Rmatrix
    do sensorIndex = 1, tvs_nsensors
      numchan = estR(sensorIndex)%nchans

      do ichan = 1, numchan
        if (estR(sensorIndex)%Rmat(ichan, ichan) < 0.0d0) then
          call rti_abort('rmat_updateRmat: Unable to take SQRT of negative values of estR%Rmat')
        end if
        obsErrStdev(sensorIndex, ichan) = SQRT(estR(sensorIndex)%Rmat(ichan, ichan))
      end do

      do ichan = 1, numchan
        do jchan = 1, numchan
          if ( utl_isEqual(obsErrStdev(sensorIndex, ichan), 0.0d0) .or. utl_isEqual(obsErrStdev(sensorIndex, jchan), 0.0d0) ) then
            RCorr(sensorIndex)%Rmat(ichan, jchan) = 0.0d0
          else
            RCorr(sensorIndex)%Rmat(ichan, jchan) = estR(sensorIndex)%Rmat(ichan, jchan) / &
                       (obsErrStdev(sensorIndex, ichan) * obsErrStdev(sensorIndex, jchan))
          end if
        end do
      end do
    end do

    ! Update simga(o) into ObsSpaceData
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER2: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER2

      idatyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if ( .not. tvs_isIdBurpTovs(idatyp) ) then
        write(*,*) 'rmat_updateRmat: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
        cycle HEADER2
      end if

      sensorIndex = tvs_lsensor(headerIndex)

      idata   = obs_headElem_i(obsspacedata, OBS_RLN, headerIndex)
      idatend = obs_headElem_i(obsspacedata, OBS_NLV, headerIndex) + idata - 1

      do bodyIndex = idata, idatend
        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) cycle

        call tvs_getChannelNumIndexFromPPP(obsSpaceData, headerIndex, bodyIndex, &
                                                channelNumber, channelIndex)

        ! Copy observation error from OBS_OER to OBS_OERI
        tmpObsErr = obs_bodyElem_r(obsSpaceData, OBS_OER, bodyIndex)
        call obs_bodySet_r(obsSpaceData, OBS_OERI, bodyIndex, tmpObsErr)

        ! Update observation error
        call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, obsErrStdev(sensorIndex, channelIndex))

      end do
    end do HEADER2

    ! Copy R-Matrix related content into Rcorr_inst object
    do sensorIndex = 1, tvs_nsensors
      Rcorr_inst(sensorIndex)%Rmat = RCorr(sensorIndex)%Rmat
      Rcorr_inst(sensorIndex)%listChans = estR(sensorIndex)%listChans
      Rcorr_inst(sensorIndex)%nchans = estR(sensorIndex)%nchans
    end do

    deallocate(obsErrStdev)
    deallocate(RCorr)
    deallocate(estR)

  end subroutine rmat_updateRmat

  !--------------------------------------------------------------------------
  ! rmat_writeRCorrFile
  !--------------------------------------------------------------------------
  subroutine rmat_writeRCorrFile
    !
    !:Purpose: Write the observation error correlation matrix into Cmat files
    !
    implicit none

    ! Locals:
    integer :: sensorIndex
    character (len=64) :: filename
    integer :: err, iunit, numChan, chanIndex1, chanIndex2

    SENSOR: do sensorIndex = 1, tvs_nsensors
      if (.not. tvs_isReallyPresentMpiGLobal(sensorIndex)) cycle SENSOR

      ! Construct correlation file name
      call rttov_coeffname (err, tvs_listSensors(:,sensorIndex), coeffname=filename, filetype='Cmat')
      if (err /= errorstatus_success) then
        write(*,*) 'Unknown instrument ', tvs_listSensors(:,sensorIndex)
        call rti_abort('rmat_writeRCorrFile')
      end if

      numChan = Rcorr_inst(sensorIndex)%nchans

      iunit = 0
      err = fnom(iunit, 'update_' // trim(filename),'FTN+SEQ+R/W',0)
      write(iunit,'(I3)') numChan

      do chanIndex1 = 1, numChan
        write(iunit, '(I3)') Rcorr_inst(sensorIndex)%listChans(chanIndex1)
      end do

      do chanIndex1 = 1, numChan
        do chanIndex2 = chanIndex1 + 1, numChan
          write(iunit, '(I3, I3, F9.5)') Rcorr_inst(sensorIndex)%listChans(chanIndex1), Rcorr_inst(sensorIndex)%listChans(chanIndex2), &
                       Rcorr_inst(sensorIndex)%Rmat(chanIndex1, chanIndex2)
        end do
      end do

      err = fclos(iunit)
    end do SENSOR
  end subroutine rmat_writeRCorrFile

end module rMatrix_mod

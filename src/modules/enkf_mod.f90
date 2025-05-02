module enkf_mod
  ! MODULE enkf_mod (prefix='enkf' category='1. High-level functionality')
  !
  !:Purpose:  Various routines that are useful for implementing
  !           an EnKF in MIDAS, including the LETKF.
  !
  use mpi                 ! this is the mpi library module
  use midasMpi_mod
  use utilities_mod
  use mathPhysConstants_mod
  use timeCoord_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use ensembleStateVector_mod
  use gridStateVector_mod
  use obsSpaceData_mod
  use tovsNL_mod
  use ensembleObservations_mod
  use gridVariableTransforms_mod
  use localizationFunction_mod
  use varNameList_mod
  use codePrecision_mod
  use codTyp_mod
  use calcHeightAndPressure_mod
  implicit none
  save
  private

  ! public types
  public :: struct_enkfInterpInfo

  ! public procedures
  public :: enkf_setupInterpInfo, enkf_LETKFanalyses, enkf_modifyAMSUBobsError
  public :: enkf_rejectHighLatIR, enkf_getModulatedState, enkf_setupModulationFactor

  ! for weight interpolation
  type struct_enkfInterpInfo
    integer              :: latLonStep
    integer, allocatable :: numIndexes(:,:)
    integer, allocatable :: lonIndexes(:,:,:)
    integer, allocatable :: latIndexes(:,:,:)
    real(8), allocatable :: interpWeights(:,:,:)
    integer              :: myLonBegHalo
    integer              :: myLonEndHalo
    integer              :: myLatBegHalo
    integer              :: myLatEndHalo
    integer              :: myLonBeg
    integer              :: myLonEnd
    integer              :: myLatBeg
    integer              :: myLatEnd
  end type struct_enkfInterpInfo

  integer, external :: get_max_rss

contains

  !----------------------------------------------------------------------
  ! enkf_LETKFanalyses
  !----------------------------------------------------------------------
  subroutine enkf_LETKFanalyses(algorithm, numSubEns, randomShuffleSubEns,  &
                                ensembleAnl, ensembleTrl, &
                                ensObs_mpiglobal, ensObsGain_mpiglobal, &
                                stateVectorMeanAnl, &
                                wInterpInfo, maxNumLocalObs, maxNumLocalObsPerType, &
                                hLocalize, hLocalizePressure, hLinearLoc, vLocalize,  &
                                mpiDistribution, numRetainedEigen, myNumLatLonSendFactor, &
                                localSelectionOutput, localObsSorting)
    !
    !:Purpose: Local subroutine containing the code for computing
    !          the LETKF analyses for all ensemble members, ensemble
    !          mean. Two MPI distributions are currently available for
    !          distributing the calculation of the LETKF weights at each
    !          grid point:
    !
    !           - 'MASTERWORKER': Dynamic load balancing with the first MPI task acting
    !             as the 'MASTER' with the sole job of assigning work to all of the other
    !             MPI tasks. When all of the LETKF weights are computed, then the MASTER
    !             send a signal to all workers that we are done and can move on to the next
    !             part of the algorithm.
    !           - 'ROUNDROBIN': Deal out all grid points where weights are computed to all
    !             mpi tasks without considering where they are needed for use. This should
    !             provide the best load balance in terms of the cost of covariance
    !             calculation which can vary depending on the local density of observations,
    !             but requires the most communication.
    !
    implicit none

    ! Arguments:
    character(len=*),            intent(in)    :: algorithm
    integer         ,            intent(in)    :: numSubEns
    logical         ,            intent(in)    :: randomShuffleSubEns
    type(struct_ens), pointer,   intent(inout) :: ensembleTrl
    type(struct_ens),            intent(inout) :: ensembleAnl
    type(struct_eob), target,    intent(in)    :: ensObs_mpiglobal
    type(struct_eob),            intent(in)    :: ensObsGain_mpiglobal
    type(struct_gsv),            intent(in)    :: stateVectorMeanAnl
    type(struct_enkfInterpInfo), intent(in)    :: wInterpInfo
    integer,                     intent(in)    :: maxNumLocalObs
    integer,                     intent(in)    :: maxNumLocalObsPerType
    real(8),                     intent(in)    :: hLocalize(:)
    real(8),                     intent(in)    :: hLocalizePressure(:)
    logical,                     intent(in)    :: hLinearLoc
    real(8),                     intent(in)    :: vLocalize
    character(len=*),            intent(in)    :: mpiDistribution
    integer,                     intent(in)    :: numRetainedEigen
    integer,                     intent(in)    :: myNumLatLonSendFactor
    logical,                     intent(in)    :: localSelectionOutput
    character(len=*),            intent(in)    :: localObsSorting

    ! Locals:
    character :: readySignal
    integer :: workerProcID, finishedSignal, assignmentTag, readyTag, numFinished
    integer :: mpiStatus(MPI_STATUS_SIZE)
    integer, allocatable :: waitStatusesSend(:,:), waitStatusesRecv(:,:)
    integer :: latLonIndex, nEns, nLev_M, nLev_depth, nLev_weights
    integer :: memberIndex, memberIndex1, memberIndex2, ierr
    integer :: procIndex, procIndexSend, latLonIndexMpiGlobal
    integer :: latIndex, lonIndex, stepIndex, varLevIndex, levIndex, levIndex2
    integer :: countMaxExceeded, maxCountMaxExceeded, numGridPointWeights
    integer :: myNumLatLonRecv, numLatLonMpiGlobal, myNumLatLonCalcMax, myNumLatLonSendMax
    integer :: sendTag, recvTag, nsize, numRecv, numSend
    integer :: myLonBeg, myLonEnd, myLatBeg, myLatEnd, numVarLev
    integer :: myLonBegHalo, myLonEndHalo, myLatBegHalo, myLatEndHalo
    integer :: nEnsGain, eigenVectorColumnIndex
    integer :: memberIndexInModEns
    integer :: requestIdRecvFinished(mmpi_nprocs-1), requestIdSendFinished(mmpi_nprocs-1)
    integer :: requestIdSignal

    integer, allocatable :: levFromK(:)
    integer, allocatable :: myLatIndexesRecv(:), myLonIndexesRecv(:)
    integer, allocatable :: latIndexesSendMpiGlobal(:), lonIndexesSendMpiGlobal(:)
    integer, allocatable :: numProcsSendMpiGlobal(:)
    integer, allocatable :: procIndexesSendMpiGlobal(:,:)
    integer, allocatable :: requestIdRecv(:), requestIdSend(:)
    integer, allocatable :: latLonTagMpiGlobal(:,:)
    integer, allocatable :: levIndex2FromVarLevIndex(:)

    real(4)              :: pert_r4, modulationFactor_r4
    real(8), allocatable :: weightsMembers(:,:,:,:), weightsMembersLatLon(:,:,:)
    real(8), allocatable :: weightsMean(:,:,:,:), weightsMeanLatLon(:,:,:)
    real(8), allocatable :: memberAnlPert(:)
    real(4), allocatable :: vertLocation_r4(:,:,:)
    real(4), pointer     :: meanTrl_ptr_r4(:,:,:,:), meanAnl_ptr_r4(:,:,:,:), meanInc_ptr_r4(:,:,:,:)
    real(4), pointer     :: memberTrl_ptr_r4(:,:,:,:), memberAnl_ptr_r4(:,:,:,:)
    real(8), allocatable :: weightsSendCombined(:,:,:), weightsRecvCombined(:,:,:)

    character(len=4)     :: varName
    character(len=2), allocatable :: varKindFromK(:)
    character(len=4), allocatable :: varLevelFromK(:)

    type(struct_hco), pointer :: hco_ens
    type(struct_vco), pointer :: vco_ens
    type(struct_gsv)          :: stateVectorMeanInc
    type(struct_gsv)          :: stateVectorMeanTrl

    logical :: procAlreadyFinished(mmpi_nprocs)
    logical :: useModulatedEns, masterIsMe

    call utl_tmg_start(131,'--LETKFanalysis')

    write(*,*)
    write(*,*) 'enkf_LETKFanalyses: starting'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    nEns       = ens_getNumMembers(ensembleAnl)
    useModulatedEns = (numRetainedEigen > 0)
    if ( useModulatedEns ) then
      nEnsGain   = nEns * numRetainedEigen
    else
      nEnsGain   = nEns
    end if

    !
    ! Set things up for the redistribution of work across mpi tasks
    !
    call enkf_LETKFsetupMpiDistribution(numLatLonMpiGlobal, myNumLatLonRecv,  &
                                        myLatIndexesRecv, myLonIndexesRecv,   &
                                        latIndexesSendMpiGlobal, lonIndexesSendMpiGlobal,   &
                                        procIndexesSendMpiGlobal, &
                                        numProcsSendMpiGlobal, wInterpInfo)

    write(*,*) 'enkf_LETKFanalyses: numLatLonMpiGlobal, maxval(numProcsSendMpiGlobal) = ', &
                                    numLatLonMpiGlobal, maxval(numProcsSendMpiGlobal)

    ! Compute maximum expected number of grid points where weights computed on each mpi task 
    if (trim(mpiDistribution) == 'ROUNDROBIN') then
      myNumLatLonCalcMax = ceiling(real(numLatLonMpiGlobal)/real(mmpi_nprocs))
    else if (trim(mpiDistribution) == 'MASTERWORKER') then
      myNumLatLonCalcMax = myNumLatLonSendFactor*ceiling(real(numLatLonMpiGlobal)/real(mmpi_nprocs))
    else
      write(*,*) 'mpiDistribution = ', trim(mpiDistribution)
      call utl_abort('enkf_LETKFanalyses: unknown mpiDistribution')
    end if
    myNumLatLonSendMax = ceiling(real(maxval(numProcsSendMpiGlobal))*real(myNumLatLonCalcMax))
    write(*,*) 'enkf_LETKFanalyses: myNumLatLonCalcMax, myNumLatLonSendMax = ',  &
                                    myNumLatLonCalcMax, myNumLatLonSendMax
    allocate(requestIdSend(myNumLatLonSendMax))
    allocate(requestIdRecv(myNumLatLonRecv))
    allocate(waitStatusesSend(MPI_STATUS_SIZE,myNumLatLonSendMax))
    allocate(waitStatusesRecv(MPI_STATUS_SIZE,myNumLatLonRecv))
    allocate(weightsRecvCombined(nEnsGain, nEns+1, myNumLatLonRecv))
    allocate(weightsSendCombined(nEnsGain, nEns+1, myNumLatLonCalcMax))
    requestIdRecv(:) = 0
    requestIdSend(:) = 0
    waitStatusesRecv(:,:) = 0
    waitStatusesSend(:,:) = 0
    weightsRecvCombined(:,:,:) = 0.0d0
    weightsSendCombined(:,:,:) = 0.0d0

    nLev_M     = ens_getNumLev(ensembleAnl, 'MM')
    nLev_depth = ens_getNumLev(ensembleAnl, 'DP')
    nLev_weights = max(nLev_M,nLev_depth)
    if ( useModulatedEns ) nLev_weights = 1
    hco_ens => ens_getHco(ensembleAnl)
    vco_ens => ens_getVco(ensembleAnl)
    myLonBeg = stateVectorMeanAnl%myLonBeg
    myLonEnd = stateVectorMeanAnl%myLonEnd
    myLatBeg = stateVectorMeanAnl%myLatBeg
    myLatEnd = stateVectorMeanAnl%myLatEnd
    numVarLev    = stateVectorMeanAnl%nk
    myLonBegHalo = wInterpInfo%myLonBegHalo
    myLonEndHalo = wInterpInfo%myLonEndHalo
    myLatBegHalo = wInterpInfo%myLatBegHalo
    myLatEndHalo = wInterpInfo%myLatEndHalo

    allocate(memberAnlPert(nEns))

    ! Weights for mean analysis
    allocate(weightsMean(nEnsGain,1,myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
    weightsMean(:,:,:,:) = 0.0d0
    allocate(weightsMeanLatLon(nEnsGain,1,myNumLatLonCalcMax))
    weightsMeanLatLon(:,:,:) = 0.0d0
    ! Weights for member analyses
    allocate(weightsMembers(nEnsGain,nEns,myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
    weightsMembers(:,:,:,:) = 0.0d0
    allocate(weightsMembersLatLon(nEnsGain,nEns,myNumLatLonCalcMax))
    weightsMembersLatLon(:,:,:) = 0.0d0

    call gsv_allocate( stateVectorMeanTrl, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_zero(stateVectorMeanTrl)
    call gsv_allocate( stateVectorMeanInc, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_zero(stateVectorMeanInc)

    call ens_computeMean(ensembleTrl)
    call ens_copyEnsMean(ensembleTrl, stateVectorMeanTrl)

    call lfn_Setup(LocFunctionWanted='FifthOrder')

    ! compute 3D field of vertical location needed for localization
    if (vLocalize > 0.0d0) then
      call enkf_computeVertLocation(vertLocation_r4,stateVectorMeanTrl)
    end if

    allocate(varLevelFromK(numVarLev))
    allocate(levFromK(numVarLev))
    allocate(varKindFromK(numVarLev))
    do varLevIndex = 1, numVarLev
      varName = gsv_getVarNameFromK(stateVectorMeanInc,varLevIndex)
      varLevelFromK(varLevIndex) = vnl_varLevelFromVarname(varName)
      levFromK(varLevIndex) = gsv_getLevFromK(stateVectorMeanInc,varLevIndex)
      varKindFromK(varLevIndex) = vnl_varKindFromVarname(varName)
    end do

    call utl_tmg_start(141,'----Barr')
    call rpn_comm_barrier('GRID',ierr)
    call utl_tmg_stop(141)

    ! get mpi global list of tags used for mpi send/recv
    call utl_tmg_start(142, '----GetGlobalTags')
    allocate(latLonTagMpiGlobal(stateVectorMeanAnl%ni,stateVectorMeanAnl%nj))
    call enkf_LETKFgetMpiGlobalTags(latLonTagMpiGlobal,myLatIndexesRecv,myLonIndexesRecv)
    call utl_tmg_stop(142)

    ! Compute the weights for ensemble mean and members
    countMaxExceeded = 0
    maxCountMaxExceeded = 0
    numGridPointWeights = 0
    LEV_LOOP: do levIndex = 1, nLev_weights
      write(*,*) 'Computing ensemble updates for vertical level = ', levIndex
      call utl_printTime(reset_opt = (levIndex==1))

      !
      ! First post all recv instructions for communication of weights
      !
      call utl_tmg_start(132,'----CommWeights')
      call utl_tmg_start(148,'----CommWeights-irecv')

      numSend = 0
      numRecv = 0
      do latLonIndex = 1, myNumLatLonRecv
        latIndex = myLatIndexesRecv(latLonIndex)
        lonIndex = myLonIndexesRecv(latLonIndex)
        recvTag = latLonTagMpiGlobal(lonIndex,latIndex) + (levIndex-1)*maxval(latLonTagMpiGlobal)
        recvTag = 1 + mod(recvTag-1, mmpi_maxTagValue - 10) 

        nsize = nEnsGain * (nEns+1)
        numRecv = numRecv + 1
        weightsRecvCombined(:,:,latLonIndex) = -999.8d0
        call mpi_irecv( weightsRecvCombined(:,:,latLonIndex),  &
                        nsize, mmpi_datyp_real8, mpi_any_source, recvTag,  &
                        mmpi_comm_grid, requestIdRecv(numRecv), ierr )
      end do
      call utl_tmg_stop(148)
      call utl_tmg_stop(132)

      ! Set tag values for sending signals
      assignmentTag = mmpi_maxTagValue - 1
      readyTag      = mmpi_maxTagValue - 2

      masterIsMe = (mmpi_myid == 0) .and. (trim(mpiDistribution) == 'MASTERWORKER')

      MASTER_WORKER: if (masterIsMe) then  ! I am the master

        call utl_tmg_start(132,'----CommWeights')
        call utl_tmg_start(145,'----CommWeights-signals')

        procAlreadyFinished(:) = .false.

        write(*,*) 'Send all assignments'
        call utl_printTime()

        ! Loop over all gridpoints where calculations need to be performed
        write(*,*) 'Start of loop over all global grid points where weights computed'
        do latLonIndexMpiGlobal = 1, numLatLonMpiGlobal

          ReadyLoop: do ! Loop until we get a valid read signal

            ! Determine which MPI task is ready for a new work assignment
            call MPI_RECV(readySignal, 1, MPI_CHARACTER, mpi_any_source, readyTag, &  
                          mmpi_comm_grid, mpiStatus, ierr)
            workerProcID = mpiStatus(MPI_SOURCE)

            if (readySignal == 'N') then
              procAlreadyFinished(workerProcID+1) = .true.
            end if

            ! Assign this MPI task the next gridpoint to be calculated
            ! NOTE: this assignment will be ignored when readySignal is 'N'
            call MPI_RSEND(latLonIndexMpiGlobal, 1, MPI_INTEGER, workerProcID, assignmentTag, &
                           mmpi_comm_grid, mpiStatus, ierr)

            if (readySignal == 'Y') then
              exit ReadyLoop
            end if

          end do ReadyLoop

        end do

        ! Now that all work is done, we need to inform all workers
        write(*,*) 'Finished all grid points, send *finished* signal to all mpi tasks'
        call utl_printTime()

        numFinished = 0
        do procIndex = 2, mmpi_nprocs

          ! Skip this proc if he already told us that he is finished
          if (procAlreadyFinished(procIndex)) cycle

          numFinished = numFinished + 1

          ! Wait for signal for every task that last assignment is complete
          call MPI_IRECV(readySignal, 1, MPI_CHARACTER, procIndex-1, readyTag, &  
                         mmpi_comm_grid, requestIdRecvFinished(numFinished), ierr)

          ! Tell this worker we are done
          finishedSignal = 0
          call MPI_ISEND(finishedSignal, 1, MPI_INTEGER, procIndex-1, assignmentTag, &
                         mmpi_comm_grid, requestIdSendFinished(numFinished), ierr)

        end do
        call mpi_waitAll(numFinished, requestIdSendFinished(1:numFinished), MPI_STATUSES_IGNORE, ierr)
        call mpi_waitAll(numFinished, requestIdRecvFinished(1:numFinished), MPI_STATUSES_IGNORE, ierr)

        write(*,*) 'All *finished* signals received and acknowledged'
        call utl_printTime()

        call utl_tmg_stop(145)
        call utl_tmg_stop(132)

      else  ! I am a worker

        ! Main loop over grid points for computing analysis weights
        latLonIndex = 0
        latLonIndexMpiGlobal = 0
        LATLON_LOOP: do

          if (trim(mpiDistribution) == 'MASTERWORKER') then

            call utl_tmg_start(132,'----CommWeights')
            call utl_tmg_start(145,'----CommWeights-signals')

            write(*,*) 'Requesting an assignment'
            call utl_printTime()

            if ( (latLonIndex+1) <= myNumLatLonCalcMax) then

              ! Post recv to obtain assignment or signal that we are done
              call MPI_IRECV(latLonIndexMpiGlobal, 1, MPI_INTEGER, 0, assignmentTag, &  
                             mmpi_comm_grid, requestIdSignal, ierr)

              ! Signal that I am ready for assignment
              readySignal = 'Y'
              call MPI_SEND(readySignal, 1, MPI_CHAR, 0, readyTag, &
                            mmpi_comm_grid, mpiStatus, ierr)

              call mpi_wait(requestIdSignal, MPI_STATUS_IGNORE, ierr)

            else

              ! Post recv to obtain assignment that WILL BE IGNORED
              call MPI_IRECV(latLonIndexMpiGlobal, 1, MPI_INTEGER, 0, assignmentTag, &  
                             mmpi_comm_grid, requestIdSignal, ierr)

              ! Signal that I reached my maximum number of assignments
              write(*,*) 'Reached maximum number of assignments, wait until others are finished'
              readySignal = 'N'
              call MPI_SEND(readySignal, 1, MPI_CHAR, 0, readyTag, &
                            mmpi_comm_grid, mpiStatus, ierr)

              call mpi_wait(requestIdSignal, MPI_STATUS_IGNORE, ierr)

              ! No assignment, so set latLonIndex to 0 to exit LATLON_LOOP
              latLonIndexMpiGlobal = 0

            end if

            write(*,*) 'Received assignment: ', latLonIndexMpiGlobal
            call utl_printTime()

            call utl_tmg_stop(145)
            call utl_tmg_stop(132)

            ! Check if we are done
            if (latLonIndexMpiGlobal == 0) then
              write(*,*) 'Received the *finished* signal, after computing weights for this many gridpoints: ', latLonIndex
              exit LATLON_LOOP
            end if

          else if (trim(mpiDistribution) == 'ROUNDROBIN') then

            ! Find the next value of latLonIndexMpiGlobal that I am responsible for
            do
              latLonIndexMpiGlobal = latLonIndexMpiGlobal + 1
              if (latLonIndexMpiGlobal > numLatLonMpiGlobal) then
                write(*,*) 'Finished processing this many gridpoints: ', latLonIndex
                exit LATLON_LOOP
              end if
              if (mod(latLonIndexMpiGlobal-1,mmpi_nprocs)==mmpi_myid) then
                exit
              end if
            end do              

          else

            write(*,*) 'mpiDistribution = ', trim(mpiDistribution)
            call utl_abort('enkf_LETKFanalyses: unknown value of mpiDistribution')

          end if

          ! Increment main loop index
          latLonIndex = latLonIndex + 1

          ! Check if latLonIndex is larger than expected for allocations
          if (latLonIndex > myNumLatLonCalcMax) then
            write(*,*) 'enkf_LETKFanalyses: We encountered more latLonIndex values on this mpi task than expected!'
            write(*,*) '                    You should probably increase the value of namelist variable: myNumLatLonSendFactor'
            call utl_abort('enkf_LETKFanalyses: increase value of myNumLatLonSendFactor')
          end if

          latIndex = latIndexesSendMpiGlobal(latLonIndexMpiGlobal)
          lonIndex = lonIndexesSendMpiGlobal(latLonIndexMpiGlobal)

          numGridPointWeights = numGridPointWeights + 1

          !
          ! Call main subroutine for computing the LETKF analysis weights for this lev/lat/lon
          !
          call enkf_LETKFcomputeWeights(weightsMeanLatLon(:,:,latLonIndex), &
                                        weightsMembersLatLon(:,:,latLonIndex), &
                                        algorithm, numSubEns, randomShuffleSubEns, &
                                        ensembleAnl, levIndex, latIndex, lonIndex, &
                                        hLocalize, hLocalizePressure, hLinearLoc, vLocalize, &
                                        vertLocation_r4, numRetainedEigen, &
                                        maxNumLocalObs, maxNumLocalObsPerType, &
                                        countMaxExceeded, maxCountMaxExceeded, &
                                        ensObs_mpiglobal, ensObsGain_mpiglobal, &
                                        localSelectionOutput,localObsSorting)

          !
          ! Now post all send instructions (each lat-lon may be sent to multiple tasks)
          !
          call utl_tmg_start(132,'----CommWeights')
          call utl_tmg_start(149,'----CommWeights-isend')
          latIndex = latIndexesSendMpiGlobal(latLonIndexMpiGlobal)
          lonIndex = lonIndexesSendMpiGlobal(latLonIndexMpiGlobal)
          weightsSendCombined(:,1:nEns,latLonIndex) = weightsMembersLatLon(:,:,latLonIndex)
          weightsSendCombined(:,(nEns+1),latLonIndex) = weightsMeanLatLon(:,1,latLonIndex)

          do procIndex = 1, numProcsSendMpiGlobal(latLonIndexMpiGlobal)
            sendTag = latLonTagMpiGlobal(lonIndex,latIndex) + (levIndex-1)*maxval(latLonTagMpiGlobal)
            sendTag = 1 + mod(sendTag-1, mmpi_maxTagValue - 10) 
            procIndexSend = procIndexesSendMpiGlobal(latLonIndexMpiGlobal, procIndex)

            nsize = nEnsGain * (nEns+1)
            numSend = numSend + 1
            if (numSend > myNumLatLonSendMax) then
              call utl_abort('numSend larger than allowed limit')
            end if
            call mpi_isend( weightsSendCombined(:,:,latLonIndex),  &
                            nsize, mmpi_datyp_real8, procIndexSend-1, sendTag,  &
                            mmpi_comm_grid, requestIdSend(numSend), ierr )
          end do
          call utl_tmg_stop(149)
          call utl_tmg_stop(132)

        end do LATLON_LOOP

      end if MASTER_WORKER

      !
      ! Wait for RECV communications to finish before continuing
      !
      call utl_tmg_start(132,'----CommWeights')
      call utl_tmg_start(147,'----CommWeights-waitRecv')
      if ( numRecv > 0 ) then
        call mpi_waitAll(numRecv, requestIdRecv(1:numRecv), waitStatusesRecv(:,1:numRecv), ierr)
        if (ierr == mpi_err_in_status) call utl_abort('Error code return by mpi_waitAll for IRECV')
      end if
      call utl_tmg_stop(147)

      call utl_tmg_start(148,'----CommWeights-irecv')
      ! Copy from combined to main separate arrays
      do latLonIndex = 1, myNumLatLonRecv
        latIndex = myLatIndexesRecv(latLonIndex)
        lonIndex = myLonIndexesRecv(latLonIndex)
        weightsMembers(:,:,lonIndex,latIndex) = weightsRecvCombined(:,1:nEns,latLonIndex)
        weightsMean(:,1,lonIndex,latIndex)    = weightsRecvCombined(:,(nEns+1),latLonIndex)
        if (any(weightsMembers(:,:,lonIndex,latIndex) < -999.0d0)) then
          write(*,*) 'latLonIndex, latIndex, lonIndex = ', latLonIndex, latIndex, lonIndex
          write(*,*) 'weightsMembers = ', weightsMembers(:,:,lonIndex,latIndex)
          call utl_abort('Invalid weight value received!!!')
        end if
        if (any(weightsMean(:,:,lonIndex,latIndex) < -999.0d0)) then
          write(*,*) 'latLonIndex, latIndex, lonIndex = ', latLonIndex, latIndex, lonIndex
          write(*,*) 'weightsMean = ', weightsMean(:,:,lonIndex,latIndex)
          call utl_abort('Invalid weight value received!!!')
        end if
      end do
      call utl_tmg_stop(148)
      call utl_tmg_stop(132)

      !
      ! Interpolate weights from coarse to full resolution
      !
      call utl_tmg_start(140,'----InterpolateWeights')
      if (wInterpInfo%latLonStep > 1) then
        call enkf_interpWeights(wInterpInfo, weightsMean)
        call enkf_interpWeights(wInterpInfo, weightsMembers)
      end if
      call utl_tmg_stop(140)

      call utl_tmg_start(143,'----ApplyWeights')

      !
      ! Apply the weights to compute the ensemble mean and members
      !
      call gsv_getField(stateVectorMeanInc,meanInc_ptr_r4)
      call gsv_getField(stateVectorMeanTrl,meanTrl_ptr_r4)
      call gsv_getField(stateVectorMeanAnl,meanAnl_ptr_r4)

      if (.not. allocated(levIndex2FromVarLevIndex)) then
        allocate(levIndex2FromVarLevIndex(numVarLev))
        do varLevIndex = 1, numVarLev
          ! Only treat varLevIndex values that correspond with current levIndex
          if (varLevelFromK(varLevIndex) == 'SF'   .or. varLevelFromK(varLevIndex) == 'SFMM' .or. &
              varLevelFromK(varLevIndex) == 'SFTH' .or. varLevelFromK(varLevIndex) == 'SS') then
            if (varKindFromK(varLevIndex) == 'OC') then
              levIndex2 = 1
            else
              levIndex2 = max(nLev_M,nLev_depth)
            end if
          else if (varLevelFromK(varLevIndex) == 'MM' .or. varLevelFromK(varLevIndex) == 'TH' .or. varLevelFromK(varLevIndex) == 'DP') then
            levIndex2 = levFromK(varLevIndex)
          else if (varLevelFromK(varLevIndex) == 'OT') then
            ! Most (all?) variables using the 'other' coordinate are surface
            levIndex2 = max(nLev_M,nLev_depth)
          else
            write(*,*) 'varLevel = ', varLevelFromK(varLevIndex)
            call utl_abort('enkf_LETKFanalyses: unknown varLevel')
          end if
          levIndex2FromVarLevIndex(varLevIndex) = levIndex2
        end do
      end if

      !$OMP PARALLEL DO PRIVATE(latIndex, lonIndex, varLevIndex, levIndex2, memberTrl_ptr_r4, memberAnl_ptr_r4), &
      !$OMP PRIVATE(memberAnlPert, stepIndex, memberIndex, memberIndex2, memberIndex1, eigenVectorColumnIndex, pert_r4), &
      !$OMP PRIVATE(memberIndexInModEns, modulationFactor_r4)
      do latIndex = myLatBeg, myLatEnd
        LON_LOOP5: do lonIndex = myLonBeg, myLonEnd

          ! skip this grid point if all weights zero (no nearby obs)
          if (all(weightsMean(:,1,lonIndex,latIndex) == 0.0d0)) cycle LON_LOOP5

          ! Compute the ensemble mean increment and analysis
          do varLevIndex = 1, numVarLev
            levIndex2 = levIndex2FromVarLevIndex(varLevIndex)
            if (levIndex2 /= levIndex .and. .not. useModulatedEns) cycle
            memberTrl_ptr_r4 => ens_getOneLev_r4(ensembleTrl,varLevIndex)
            do stepIndex = 1, tim_nstepobsinc
              ! mean increment
              if ( useModulatedEns ) then
                do eigenVectorColumnIndex = 1, numRetainedEigen
                  call getModulationFactor( stateVectorMeanInc%vco, levIndex2, &
                                            eigenVectorColumnIndex, numRetainedEigen, &
                                            nEns, vLocalize, &
                                            modulationFactor_r4 )

                  do memberIndex = 1, nEns
                    pert_r4 = modulationFactor_r4 * ( memberTrl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) -  &
                                        meanTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) )

                    ! Index of the modulated ensemble member corresponding to original
                    ! ensemble member index (memberIndex1) and eigenVectorColumnIndex.
                    memberIndexInModEns = (eigenVectorColumnIndex - 1) * nEns + memberIndex

                    meanInc_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) =  &
                        meanInc_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) +  &
                        weightsMean(memberIndexInModEns,1,lonIndex,latIndex) * pert_r4
                  end do
                end do
              else
                do memberIndex = 1, nEns
                  meanInc_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) =  &
                       meanInc_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) +  &
                       weightsMean(memberIndex,1,lonIndex,latIndex) *  &
                       (memberTrl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) -  &
                        meanTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex))
                end do
              end if

              ! mean analysis
              meanAnl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) =  &
                   meanTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) +  &
                   meanInc_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex)
            end do ! stepIndex
          end do ! varLevIndex

          ! Compute the ensemble member analyses
          call utl_tmg_start(144,'------ApplyWeightsMember')
          do varLevIndex = 1, numVarLev
            levIndex2 = levIndex2FromVarLevIndex(varLevIndex)
            if (levIndex2 /= levIndex .and. .not. useModulatedEns) cycle
            memberTrl_ptr_r4 => ens_getOneLev_r4(ensembleTrl,varLevIndex)
            memberAnl_ptr_r4 => ens_getOneLev_r4(ensembleAnl,varLevIndex)
            do stepIndex = 1, tim_nstepobsinc

              ! Compute analysis member perturbation
              memberAnlPert(:) = 0.0d0

              if ( useModulatedEns ) then
                do memberIndex2 = 1, nEns
                  do eigenVectorColumnIndex = 1, numRetainedEigen
                    call getModulationFactor( stateVectorMeanInc%vco, levIndex2, &
                          eigenVectorColumnIndex, numRetainedEigen, &
                          nEns, vLocalize, &
                          modulationFactor_r4 )

                    do memberIndex1 = 1, nEns
                      ! Compute background ensemble perturbations for the modulated ensemble (Xb_Mod)
                      pert_r4 = modulationFactor_r4 * ( memberTrl_ptr_r4(memberIndex1,stepIndex,lonIndex,latIndex) -  &
                                          meanTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) )

                      ! Index of the modulated ensemble member corresponding to original
                      ! ensemble member index (memberIndex1) and eigenVectorColumnIndex.
                      memberIndexInModEns = (eigenVectorColumnIndex - 1) * nEns + memberIndex1
                      
                      ! sum Xb_Mod * Wa over all modulated ensembles to get member perturbations for
                      !   original ensemble (memberIndex2)
                      memberAnlPert(memberIndex2) = memberAnlPert(memberIndex2) + &
                           weightsMembers(memberIndexInModEns,memberIndex2,lonIndex,latIndex) *  pert_r4
                    end do
                  end do

                  ! Compute final member perturbations by removing background original ensemble perturbations
                  memberAnlPert(memberIndex2) = (memberTrl_ptr_r4(memberIndex2,stepIndex,lonIndex,latIndex) -  &
                                                 meanTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex)) + &
                                                 memberAnlPert(memberIndex2)

                end do ! memberIndex2
              else
                do memberIndex2 = 1, nEns
                  do memberIndex1 = 1, nEns
                    memberAnlPert(memberIndex2) = memberAnlPert(memberIndex2) + &
                         weightsMembers(memberIndex1,memberIndex2,lonIndex,latIndex) *  &
                         (memberTrl_ptr_r4(memberIndex1,stepIndex,lonIndex,latIndex) -  &
                         meanTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex))
                  end do ! memberIndex1
                end do ! memberIndex2
              end if

              ! Add analysis member perturbation to mean analysis
              memberAnl_ptr_r4(:,stepIndex,lonIndex,latIndex) =  &
                   meanAnl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) + memberAnlPert(:)


            end do ! stepIndex
          end do ! varLevIndex
          call utl_tmg_stop(144)

        end do LON_LOOP5
      end do
      !$OMP END PARALLEL DO

      call utl_tmg_stop(143)

      !
      ! Wait for SEND communications to finish before continuing
      !
      call utl_tmg_start(132,'----CommWeights')
      call utl_tmg_start(146,'----CommWeights-waitSend')
      if ( numSend > 0 ) then
        call mpi_waitAll(numSend, requestIdSend(1:numSend), waitStatusesSend(:,1:numSend), ierr)
        if (ierr == mpi_err_in_status) call utl_abort('Error code return by mpi_waitAll for ISEND')
      end if
      call utl_tmg_stop(146)
      call utl_tmg_stop(132)
      ! Should be safe to reuse buffer used for mpi_isend and the call to mpi_waitAll
      weightsSendCombined(:,:,:) = -999.9d0

    end do LEV_LOOP

    if (countMaxExceeded > 0) then
      write(*,*) 'enkf_LETKFanalyses: WARNING: Found more local obs than specified max number at ', &
                 real(100*countMaxExceeded)/real(numGridPointWeights), '% of grid points.'
      write(*,*) '                      Maximum number found was ', maxCountMaxExceeded,  &
                 ' which is greater than specified number ', maxNumLocalObs
      write(*,*) '                      Therefore will keep closest obs only.'
    end if

    call utl_tmg_start(141,'----Barr')
    call rpn_comm_barrier('GRID',ierr)
    call utl_tmg_stop(141)

    call gsv_deallocate(stateVectorMeanInc)
    call gsv_deallocate(stateVectorMeanTrl)

    write(*,*) 'enkf_LETKFanalyses: done'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*)

    call utl_tmg_stop(131)

  end subroutine enkf_LETKFanalyses

  !----------------------------------------------------------------------
  ! enkf_LETKFcomputeWeights (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_LETKFcomputeWeights(weightsMeanLatLon, weightsMembersLatLon, &
                                      algorithm, numSubEns, randomShuffleSubEns, &
                                      ensembleAnl, levIndex, latIndex, lonIndex, &
                                      hLocalize, hLocalizePressure, hLinearLoc, vLocalize, &
                                      vertLocation_r4, numRetainedEigen, &
                                      maxNumLocalObs, maxNumLocalObsPerType, &
                                      countMaxExceeded, maxCountMaxExceeded, &
                                      ensObs_mpiglobal, ensObsGain_mpiglobal, &
                                      localSelectionOutput, localObsSorting)
    !
    !:Purpose:
    !
    implicit none

    ! Arguments:
    real(8),                  intent(inout) :: weightsMeanLatLon(:,:)
    real(8),                  intent(inout) :: weightsMembersLatLon(:,:)
    character(len=*),         intent(in)    :: algorithm
    integer         ,         intent(in)    :: numSubEns
    logical         ,         intent(in)    :: randomShuffleSubEns
    type(struct_ens),         intent(inout) :: ensembleAnl
    integer,                  intent(in)    :: levIndex
    integer,                  intent(in)    :: latIndex
    integer,                  intent(in)    :: lonIndex
    real(8),                  intent(in)    :: hLocalize(:)
    real(8),                  intent(in)    :: hLocalizePressure(:)
    logical,                  intent(in)    :: hLinearLoc
    real(8),                  intent(in)    :: vLocalize
    real(4),                  intent(in)    :: vertLocation_r4(:,:,:)
    integer,                  intent(in)    :: numRetainedEigen
    integer,                  intent(in)    :: maxNumLocalObs
    integer,                  intent(in)    :: maxNumLocalObsPerType
    integer,                  intent(inout) :: countMaxExceeded
    integer,                  intent(inout) :: maxCountMaxExceeded
    type(struct_eob), target, intent(in)    :: ensObs_mpiglobal
    type(struct_eob),         intent(in)    :: ensObsGain_mpiglobal
    logical,                  intent(in)    :: localSelectionOutput
    character(len=*),         intent(in)    :: localObsSorting

    ! Locals:
    type(struct_hco), pointer :: hco_ens

    integer :: nEns, nEnsGain, matrixRank, imode, dateStamp, ierr
    integer :: datePrint, timePrint, newDate, randomSeed
    integer :: numLocalObs, numLocalObsFound, localObsIndex
    integer :: bodyIndex, memberIndex, memberIndex1, memberIndex2
    integer :: subEnsIndex, subEnsIndex2, memberIndexCV, memberIndexCV1, memberIndexCV2
    integer :: nEnsIndependentPerSubEns, nEnsPerSubEns, nEnsPerSubEns_mod
    integer :: eigenVectorColumnIndex, memberIndexInModEns

    !logical :: hLocalizeIsConstant, useModulatedEns
    logical :: useModulatedEns

    real(8) :: hLoc, anlLat, anlLon, anlVertLocation
    real(8) :: tolerance

    real(8), allocatable, target, save :: YbTinvRYb_pert(:,:)
    real(8), allocatable, target, save :: YbTinvRCopy_pert(:,:)
 
    integer, allocatable, save :: localBodyIndices(:)
    integer, allocatable, save :: memberIndexSubEns(:,:), memberIndexSubEns_mod(:,:)
    integer, allocatable, save :: memberIndexSubEnsComp(:,:)
    integer, allocatable, save :: randomMemberIndexArray(:)

    real(8), allocatable, save :: locFun(:), PaSqrt_pert(:,:)
    real(8), pointer, save :: PaInv_mean(:,:), Pa_mean(:,:)
    real(8), pointer, save :: YbTinvR_mean(:,:), YbTinvRCopy_mean(:,:), YbTinvRYb_mean(:,:)
    real(8), pointer, save :: eigenValues_mean(:), eigenVectors_mean(:,:)
    
    real(8), allocatable, target, save :: PaInv_pert(:,:), Pa_pert(:,:)
    real(8), allocatable, target, save :: YbTinvR_pert(:,:)
    real(8), allocatable, target, save :: eigenValues_pert(:), eigenVectors_pert(:,:)

    real(8), allocatable, save :: YbTinvRYb_CV(:,:), YbTinvRYb_mod(:,:)
    real(8), allocatable, save :: eigenValues_CV(:), eigenVectors_CV(:,:)
    real(8), allocatable, save :: weightsTemp(:), weightsTemp2(:)

    logical, save :: firstCall = .true.

    hco_ens => ens_getHco(ensembleAnl)

    nEns       = ens_getNumMembers(ensembleAnl)
    useModulatedEns = (numRetainedEigen > 0)
    if ( useModulatedEns ) then
      nEnsGain   = nEns * numRetainedEigen
    else
      nEnsGain   = nEns
    end if

    if (firstCall) then
      allocate(YbTinvRYb_pert(nEnsGain,nEnsGain))
      allocate(YbTinvRCopy_pert(maxNumLocalObs,nEnsGain))
    end if

    ! Quantities needed for CVLETKF and CVLETKF-PERTOBS and CVLETKF-ME
    if ( trim(algorithm) == 'CVLETKF' .or. trim(algorithm) == 'CVLETKF-PERTOBS' .or. &
         trim(algorithm) == 'CVLETKF-ME' ) then
      nEnsPerSubEns = nEns / numSubEns
      if ( (nEnsPerSubEns * numSubEns) /= nEns ) then
        call utl_abort('enkf_LETKFanalyses: ensemble size not divisible by numSubEnsembles')
      end if
      if (numSubEns <= 1) then
        call utl_abort('enkf_LETKFanalyses: for CVLETKF(-PERTOBS)(-ME) algorithm, numSubEns must be greater than 1')
      end if
      if ( .not. useModulatedEns ) then
        nEnsIndependentPerSubEns = nEns - nEnsPerSubEns
      else
        nEnsPerSubEns_mod = nEnsPerSubEns * numRetainedEigen
        nEnsIndependentPerSubEns = nEnsGain - nEnsPerSubEns_mod
      end if
      if (firstCall) then
        allocate(YbTinvRYb_CV(nEnsIndependentPerSubEns,nEnsIndependentPerSubEns))
        allocate(eigenValues_CV(nEnsIndependentPerSubEns))
        allocate(eigenVectors_CV(nEnsIndependentPerSubEns,nEnsIndependentPerSubEns))
        allocate(memberIndexSubEns(nEnsPerSubEns,numSubEns))
        allocate(memberIndexSubEnsComp(nEnsIndependentPerSubEns,numSubEns))
        if ( useModulatedEns ) allocate(memberIndexSubEns_mod(nEnsPerSubEns_mod,numSubEns))
      end if
      if (.not.randomShuffleSubEns) then
        ! form subensembles with contiguous sequential groups of members
        do subEnsIndex = 1, numSubEns
          do memberIndex = 1, nEnsPerSubEns
            memberIndexSubEns(memberIndex,subEnsIndex) =  &
                (subEnsIndex-1)*nEnsPerSubEns + memberIndex
          end do
        end do
        if ( useModulatedEns ) then
          do subEnsIndex = 1, numSubEns
            memberIndex2 = 0
            do memberIndex = 1, nEnsPerSubEns
              do eigenVectorColumnIndex = 1, numRetainedEigen
                memberIndex2 = memberIndex2 + 1
                memberIndexInModEns = (eigenVectorColumnIndex - 1) * nEns + &
                                        memberIndex
                memberIndexSubEns_mod(memberIndex2,subEnsIndex) =  &
                     (subEnsIndex-1)*nEnsPerSubEns + memberIndexInModEns
              end do
            end do
          end do
        end if
      else
        ! compute random seed from the date for randomly forming subensembles
        imode = -3 ! stamp to printable date and time: YYYYMMDD, HHMMSShh
        dateStamp = tim_getDateStamp()
        ierr = newdate(dateStamp, datePrint, timePrint, imode)
        timePrint = timePrint/1000000
        datePrint =  datePrint*100 + timePrint
        ! Remove the century, keeping 2 digits of the year
        randomSeed = datePrint - 100000000*(datePrint/100000000)
        if (firstCall) allocate(randomMemberIndexArray(nEns))
        do memberIndex = 1, nEns
          randomMemberIndexArray(memberIndex) = memberIndex
        end do
        call utl_randomOrderInt(randomMemberIndexArray,randomSeed)
        if (firstCall) then
          write(*,*) 'enkf_LETKFcomputeWeights: seed for random shuffle of sub ens = ', randomSeed
          write(*,*) 'enkf_LETKFcomputeWeights: randomOrder = ', randomMemberIndexArray(:)
        end if
        do subEnsIndex = 1, numSubEns
          do memberIndex = 1, nEnsPerSubEns
            memberIndexSubEns(memberIndex,subEnsIndex) =  &
                 randomMemberIndexArray((subEnsIndex-1)*nEnsPerSubEns + memberIndex)
          end do
        end do
        if ( useModulatedEns ) then
          do subEnsIndex = 1, numSubEns
            memberIndex2 = 0
            do memberIndex = 1, nEnsPerSubEns
              do eigenVectorColumnIndex = 1, numRetainedEigen
                memberIndex2 = memberIndex2 + 1
                memberIndexSubEns_mod(memberIndex2,subEnsIndex) =  &
                      randomMemberIndexArray((subEnsIndex-1)*nEnsPerSubEns + memberIndex) + &
                      (eigenVectorColumnIndex - 1) * nEns
              end do
            end do
          end do
        end if        
      end if

      do subEnsIndex = 1, numSubEns
        memberIndex = 1
        do subEnsIndex2 = 1, numSubEns
          if (subEnsIndex2 == subEnsIndex) cycle
          
          if ( .not. useModulatedEns ) then
            memberIndexSubEnsComp(memberIndex:memberIndex+nEnsPerSubEns-1,subEnsIndex) =  &
              memberIndexSubEns(:,subEnsIndex2)
            memberIndex = memberIndex + nEnsPerSubEns
          else
            memberIndexSubEnsComp(memberIndex:memberIndex+nEnsPerSubEns_mod-1,subEnsIndex) =  &
              memberIndexSubEns_mod(:,subEnsIndex2)
            memberIndex = memberIndex + nEnsPerSubEns_mod
          end if
        end do
      end do

      if ( firstCall .and. mmpi_myid == 0 ) then
        write(*,*) 'nEns, numSubEns, nEnsPerSubEns, nEnsIndependentPerSubEns = ',  &
                  nEns, numSubEns, nEnsPerSubEns, nEnsIndependentPerSubEns
        do subEnsIndex = 1, numSubEns
          write(*,*) 'memberIndexSubEns = '
          write(*,*) memberIndexSubEns(:,subEnsIndex)
          if ( useModulatedEns ) then
            write(*,*) 'memberIndexSubEns_mod = '
            write(*,*) memberIndexSubEns_mod(:,subEnsIndex)
          end if
          write(*,*) 'memberIndexSubEnsComp = '
          write(*,*) memberIndexSubEnsComp(:,subEnsIndex)
        end do
      end if

    end if ! if CVLETKF(-PERTOBS)(-ME) algorithm

    if (firstCall) then
      allocate(localBodyIndices(maxNumLocalObs))
      allocate(locFun(maxNumLocalObs))
      !
      ! Compute gridded 3D ensemble weights
      !
      allocate(YbTinvR_pert(nEnsGain,maxNumLocalObs))
      if ( trim(algorithm) == 'CVLETKF-ME' .or. &
           trim(algorithm) == 'LETKF-Gain-ME' ) then
        allocate(YbTinvRYb_mod(nEnsGain,nEns))
      end if
      allocate(eigenValues_pert(nEnsGain))
      allocate(eigenVectors_pert(nEnsGain,nEnsGain))
      allocate(PaInv_pert(nEnsGain,nEnsGain))
      allocate(PaSqrt_pert(nEnsGain,nEnsGain))
      allocate(Pa_pert(nEnsGain,nEnsGain))

      if (eob_simObsAssim) then
        allocate(YbTinvR_mean(nEnsGain,maxNumLocalObs))
        allocate(YbTinvRCopy_mean(maxNumLocalObs,nEnsGain))
        allocate(YbTinvRYb_mean(nEnsGain,nEnsGain))
        allocate(eigenValues_mean(nEnsGain))
        allocate(eigenVectors_mean(nEnsGain,nEnsGain))
        allocate(PaInv_mean(nEnsGain,nEnsGain))
        allocate(Pa_mean(nEnsGain,nEnsGain))
      else
        YbTinvR_mean => YbTinvR_pert
        YbTinvRCopy_mean => YbTinvRCopy_pert
        YbTinvRYb_mean => YbTinvRYb_pert
        eigenValues_mean => eigenValues_pert
        eigenVectors_mean => eigenVectors_pert
        PaInv_mean => PaInv_pert
        Pa_mean => Pa_pert
      end if
    
      allocate(weightsTemp(nEnsGain))
      allocate(weightsTemp2(nEnsGain))
      weightsTemp(:) = 0.0d0
      weightsTemp2(:) = 0.0d0
    end if

    ! lat-lon of the grid point for which we are doing the analysis
    anlLat = hco_ens%lat2d_4(lonIndex,latIndex)
    anlLon = hco_ens%lon2d_4(lonIndex,latIndex)
    ! if there is vertical localization
    if (vLocalize > 0.0d0) then
      anlVertLocation = real(vertLocation_r4(lonIndex,latIndex,levIndex),8)
    end if

    ! Get list of nearby observations and localization functions to gridpoint.
    ! With modulated-ensembles, we get observations in entire column.
    call utl_tmg_start(133,'----GetLocalBodyIndices')
    if ( useModulatedEns ) anlVertLocation = MPC_missingValue_R8

    ! Find horizontal localization value for this vertical level
    call enkf_getLocalizationRadius(hLocalize, hLocalizePressure, anlVertLocation, hLinearLoc, hLoc)

    numLocalObs = eob_getLocalBodyIndices(ensObs_mpiglobal, localBodyIndices,     &
                                          locFun, anlLat, anlLon, anlVertLocation,  &
                                          hLoc, vLocalize, numLocalObsFound, &
                                          maxNumLocalObsPerType, localSelectionOutput, localObsSorting)
    if (numLocalObsFound > maxNumLocalObs) then
      countMaxExceeded = countMaxExceeded + 1
      maxCountMaxExceeded = max(maxCountMaxExceeded, numLocalObsFound)
    end if
    call utl_tmg_stop(133)

    call utl_tmg_start(134,'----CalculateWeights')

    ! Extract initial quantities YbTinvR and first term of PaInv (YbTinvR*Yb)
    do localObsIndex = 1, numLocalObs
      bodyIndex = localBodyIndices(localObsIndex)

      do memberIndex = 1, nEnsGain
        ! YbTinvR for updating ensemble perturbations
        YbTinvR_pert(memberIndex,localObsIndex) =  &
             ensObsGain_mpiglobal%Yb_r4(memberIndex, bodyIndex) * &
             locFun(localObsIndex) * ensObsGain_mpiglobal%obsErrInv(bodyIndex)
      end do
      if (eob_simObsAssim) then
        do memberIndex = 1, nEnsGain
          ! YbTinvR for the ensemble mean update for EDA observation simulation experiment
          YbTinvR_mean(memberIndex,localObsIndex) =  &
               ensObsGain_mpiglobal%Yb_r4(memberIndex, bodyIndex) * &
               locFun(localObsIndex) * ensObsGain_mpiglobal%obsErrInv_sim(bodyIndex)
        end do
      end if
    end do ! localObsIndex

    call utl_tmg_start(136,'------CalcYbTinvRYb')

    ! make copy of YbTinvR, and ensObsGain_mpiglobal%Yb_r4
    call utl_tmg_start(137,'--------YbArraysCopy')
    !$OMP PARALLEL DO PRIVATE (localObsIndex, bodyIndex, memberIndex2)
    do localObsIndex = 1, numLocalObs
      bodyIndex = localBodyIndices(localObsIndex)
      do memberIndex2 = 1, nEnsGain
        YbTinvRCopy_pert(localObsIndex,memberIndex2) = YbTinvR_pert(memberIndex2,localObsIndex)
      end do
    end do
    !$OMP END PARALLEL DO
    if (eob_simObsAssim) then
      !$OMP PARALLEL DO PRIVATE (localObsIndex, bodyIndex, memberIndex2)
      do localObsIndex = 1, numLocalObs
        bodyIndex = localBodyIndices(localObsIndex)
        do memberIndex2 = 1, nEnsGain
          YbTinvRCopy_mean(localObsIndex,memberIndex2) = YbTinvR_mean(memberIndex2,localObsIndex)             
        end do
      end do
      !$OMP END PARALLEL DO
    end if
    call utl_tmg_stop(137)

    call utl_tmg_start(138,'--------YbTinvRYb1')
    call enkf_calcYbTinvRYb(nEnsGain, nEnsGain, maxNumLocalObs, numLocalObs, &
                            YbTinvRYb_pert, YbTinvRCopy_pert, &
                            ensObsGain_mpiglobal, localBodyIndices, &
                            YbTinvRYb_mean, YbTinvRCopy_mean)
    call utl_tmg_stop(138)

    ! computing YbTinvRYb that uses modulated and original ensembles for perturbation update
    if ( trim(algorithm) == 'CVLETKF-ME' .or. &
         trim(algorithm) == 'LETKF-Gain-ME' ) then

      call utl_tmg_start(139,'--------YbTinvRYb2')
      call enkf_calcYbTinvRYb(nEnsGain, nEns, maxNumLocalObs, numLocalObs, &
                              YbTinvRYb_mod, YbTinvRCopy_pert, &
                              ensObs_mpiglobal, localBodyIndices)
      call utl_tmg_stop(139)

    end if !CVLETKF-ME or LETKF-GAIN-ME

    call utl_tmg_stop(136)

    ! Rest of the computation of local weights for this grid point
    if (numLocalObs > 0) then

      if (trim(algorithm) == 'LETKF') then
        !
        ! Weight calculation for standard LETKF algorithm
        !

        ! Add second term of PaInv
        PaInv_pert(:,:) = YbTinvRYb_pert(:,:)
        do memberIndex = 1, nEns
          PaInv_pert(memberIndex,memberIndex) = PaInv_pert(memberIndex,memberIndex) + real(nEns - 1,8)
        end do
        if (eob_simObsAssim) then
          PaInv_mean(:,:) = YbTinvRYb_mean(:,:)            
          do memberIndex = 1, nEns
            PaInv_mean(memberIndex,memberIndex) = PaInv_mean(memberIndex,memberIndex) + real(nEns - 1,8)
          end do
        end if

        ! Compute Pa and sqrt(Pa) matrices from PaInv
        Pa_pert(:,:) = PaInv_pert(:,:)
        call utl_tmg_start(135,'------EigenDecomp')
        call utl_matInverse(Pa_pert, nEns, inverseSqrt_opt=PaSqrt_pert)
        if (eob_simObsAssim) then
          Pa_mean(:,:) = PaInv_mean(:,:)
          call utl_matInverse(Pa_mean, nEns)
        end if
        call utl_tmg_stop(135)

        ! Compute ensemble mean local weights as Pa * YbTinvR * (obs - meanYb)
        weightsTemp(:) = 0.0d0
        do localObsIndex = 1, numLocalObs
          bodyIndex = localBodyIndices(localObsIndex)
          do memberIndex = 1, nEns
            weightsTemp(memberIndex) = weightsTemp(memberIndex) +   &
                                       YbTinvR_mean(memberIndex,localObsIndex) *  &
                                       ( ensObs_mpiglobal%obsValue(bodyIndex) - &
                                         ensObs_mpiglobal%meanYb(bodyIndex) )
          end do
        end do

        weightsMeanLatLon(:,1) = 0.0d0
        do memberIndex2 = 1, nEns
          do memberIndex1 = 1, nEns
            weightsMeanLatLon(memberIndex1,1) =  &
                 weightsMeanLatLon(memberIndex1,1) +  &
                 Pa_mean(memberIndex1,memberIndex2)*weightsTemp(memberIndex2)
          end do
        end do

        ! Compute ensemble perturbation weights: [(Nens-1)^1/2*PaSqrt]
        weightsMembersLatLon(:,:) = sqrt(real(nEns - 1,8)) * PaSqrt_pert(:,:)

      else if (trim(algorithm) == 'LETKF-Gain') then
        !
        ! Weight calculation for standard LETKF algorithm
        !

        ! Compute eigenValues/Vectors of Yb^T R^-1 Yb = E * Lambda * E^T
        call utl_tmg_start(135,'------EigenDecomp')
        tolerance = 1.0D-50
        call utl_eigenDecomp(YbTinvRYb_pert, eigenValues_pert, eigenVectors_pert, tolerance, matrixRank)
        if (eob_simObsAssim) then
          call utl_eigenDecomp(YbTinvRYb_mean, eigenValues_mean, eigenVectors_mean, tolerance, matrixRank)
        end if
        call utl_tmg_stop(135)

        ! Compute ensemble mean local weights as E * (Lambda + (Nens-1)*I)^-1 * E^T * YbTinvR * (obs - meanYb)
        weightsTemp(:) = 0.0d0
        do localObsIndex = 1, numLocalObs
          bodyIndex = localBodyIndices(localObsIndex)
          do memberIndex = 1, nEns
            weightsTemp(memberIndex) = weightsTemp(memberIndex) +   &
                                       YbTinvR_mean(memberIndex,localObsIndex) *  &
                                       ( ensObs_mpiglobal%obsValue(bodyIndex) - &
                                         ensObs_mpiglobal%meanYb(bodyIndex) )
          end do
        end do
        weightsTemp2(:) = 0.0d0
        do memberIndex2 = 1, matrixRank
          do memberIndex1 = 1, nEns
            weightsTemp2(memberIndex2) = weightsTemp2(memberIndex2) +   &
                                         eigenVectors_mean(memberIndex1,memberIndex2) *  &
                                         weightsTemp(memberIndex1)
          end do
        end do
        do memberIndex = 1, matrixRank
          weightsTemp2(memberIndex) = weightsTemp2(memberIndex) *  &
                                      1.0D0/(eigenValues_mean(memberIndex) + real(nEns - 1,8))
        end do
        weightsMeanLatLon(:,1) = 0.0d0
        do memberIndex2 = 1, matrixRank
          do memberIndex1 = 1, nEns
            weightsMeanLatLon(memberIndex1,1) =  &
                 weightsMeanLatLon(memberIndex1,1) +   &
                 eigenVectors_mean(memberIndex1,memberIndex2) *  &
                 weightsTemp2(memberIndex2)
          end do
        end do

        ! Compute ensemble perturbation weights: 
        ! Wa = [ - (Nens-1)^1/2 * E *
        !        {(Nens-1)^-1/2*I - (Lambda + (Nens-1)*I)^-1/2} * Lambda^-1 *
        !        E^T * YbTinvRYb ]
        ! Loop over members within the current sub-ensemble being updated
        do memberIndex = 1, nEns

          ! E^T * YbTinvRYb
          weightsTemp(:) = 0.0d0
          do memberIndex2 = 1, matrixRank
            do memberIndex1 = 1, nEns
              weightsTemp(memberIndex2) = weightsTemp(memberIndex2) +  &
                                          eigenVectors_pert(memberIndex1,memberIndex2) *  &
                                          YbTinvRYb_pert(memberIndex1,memberIndex)
            end do
          end do

          ! {(Nens-1)^-1/2*I - (Lambda + (Nens-1)*I)^-1/2} Lambda^-1 * previous_result
          
          do memberIndex1 = 1, matrixRank
            weightsTemp(memberIndex1) = weightsTemp(memberIndex1) *  &
                                        ( 1.0D0/sqrt(real(nEns - 1,8)) -   &
                                          1.0D0/sqrt(eigenValues_pert(memberIndex1) +  &
                                                      real(nEns - 1,8)) )
            weightsTemp(memberIndex1) = weightsTemp(memberIndex1) /  &
                                        eigenValues_pert(memberIndex1)
          end do

          ! E * previous_result
          weightsMembersLatLon(:,memberIndex) = 0.0d0
          do memberIndex2 = 1, matrixRank
            do memberIndex1 = 1, nEns
              weightsMembersLatLon(memberIndex1,memberIndex) =   &
                   weightsMembersLatLon(memberIndex1,memberIndex) +   &
                   eigenVectors_pert(memberIndex1,memberIndex2) *  &
                   weightsTemp(memberIndex2)
            end do
          end do

          ! -1 * (Nens-1)^1/2 * previous_result
          weightsMembersLatLon(:,memberIndex) =  &
               -1.0D0 * sqrt(real(nEns - 1,8)) *  &
               weightsMembersLatLon(:,memberIndex)

          ! I + previous_result
          weightsMembersLatLon(memberIndex,memberIndex) =  &
               1.0D0 + weightsMembersLatLon(memberIndex,memberIndex)

        end do

        ! Remove the weights mean computed over the columns
        do memberIndex = 1, nEns
          weightsMembersLatLon(memberIndex,:) =  &
               weightsMembersLatLon(memberIndex,:) - &
               sum(weightsMembersLatLon(memberIndex,:))/real(nEns,8)
        end do

      else if (trim(algorithm) == 'LETKF-Gain-ME') then
        !
        ! Weight calculation for standard LETKF algorithm with modulated ensemble
        !

        ! Compute eigenValues/Vectors of Yb^T R^-1 Yb = E * Lambda * E^T
        call utl_tmg_start(135,'------EigenDecomp')
        tolerance = 1.0D-50
        call utl_eigenDecomp(YbTinvRYb_pert, eigenValues_pert, eigenVectors_pert, tolerance, matrixRank)
        if (eob_simObsAssim) then
          call utl_eigenDecomp(YbTinvRYb_mean, eigenValues_mean, eigenVectors_mean, tolerance, matrixRank)
        end if
        call utl_tmg_stop(135)

        ! Compute ensemble mean local weights as E * (Lambda + (Nens-1)*I)^-1 * E^T * YbTinvR * (obs - meanYb)
        weightsTemp(:) = 0.0d0
        do localObsIndex = 1, numLocalObs
          bodyIndex = localBodyIndices(localObsIndex)
          do memberIndex = 1, nEnsGain
            weightsTemp(memberIndex) = weightsTemp(memberIndex) +   &
                                       YbTinvR_mean(memberIndex,localObsIndex) *  &
                                       (ensObs_mpiglobal%obsValue(bodyIndex) - &
                                        ensObs_mpiglobal%meanYb(bodyIndex))
          end do
        end do
        weightsTemp2(:) = 0.0d0
        do memberIndex2 = 1, matrixRank
          do memberIndex1 = 1, nEnsGain
            weightsTemp2(memberIndex2) = weightsTemp2(memberIndex2) +   &
                                         eigenVectors_mean(memberIndex1,memberIndex2) *  &
                                         weightsTemp(memberIndex1)
          end do
        end do
        do memberIndex = 1, matrixRank
          weightsTemp2(memberIndex) = weightsTemp2(memberIndex) *  &
                                      1.0D0/(eigenValues_mean(memberIndex) + real(nEnsGain - 1,8))
        end do
        weightsMeanLatLon(:,1) = 0.0d0
        do memberIndex2 = 1, matrixRank
          do memberIndex1 = 1, nEnsGain
            weightsMeanLatLon(memberIndex1,1) =  &
                 weightsMeanLatLon(memberIndex1,1) +   &
                 eigenVectors_mean(memberIndex1,memberIndex2) *  &
                 weightsTemp2(memberIndex2)
          end do
        end do

        ! Compute ensemble perturbation weights: 
        ! Wa = [ - (Nens-1)^1/2 * E *
        !        {(Nens-1)^-1/2*I - (Lambda + (Nens-1)*I)^-1/2} * Lambda^-1 *
        !        E^T * YbTinvRYb_mod ]
        ! Loop over members within the current sub-ensemble being updated
        do memberIndex = 1, nEns

          ! E^T * YbTinvRYb_mod
          weightsTemp(:) = 0.0d0
          do memberIndex2 = 1, matrixRank
            do memberIndex1 = 1, nEnsGain
              weightsTemp(memberIndex2) = weightsTemp(memberIndex2) +  &
                                          eigenVectors_pert(memberIndex1,memberIndex2) *  &
                                          YbTinvRYb_mod(memberIndex1,memberIndex)
            end do
          end do

          ! {(Nens-1)^-1/2*I - (Lambda + (Nens-1)*I)^-1/2} Lambda^-1 * previous_result

          do memberIndex1 = 1, matrixRank
            weightsTemp(memberIndex1) = weightsTemp(memberIndex1) *  &
                                        ( 1.0D0/sqrt(real(nEnsGain - 1,8)) -   &
                                          1.0D0/sqrt(eigenValues_pert(memberIndex1) +  &
                                                      real(nEnsGain - 1,8)) )
            weightsTemp(memberIndex1) = weightsTemp(memberIndex1) /  &
                                        eigenValues_pert(memberIndex1)
          end do

          ! E * previous_result
          weightsMembersLatLon(:,memberIndex) = 0.0d0
          do memberIndex2 = 1, matrixRank
            do memberIndex1 = 1, nEnsGain
              weightsMembersLatLon(memberIndex1,memberIndex) =   &
                   weightsMembersLatLon(memberIndex1,memberIndex) +   &
                   eigenVectors_pert(memberIndex1,memberIndex2) *  &
                   weightsTemp(memberIndex2)
            end do
          end do

          ! -1 * (Nens-1)^1/2 * previous_result
          weightsMembersLatLon(:,memberIndex) =  &
               -1.0D0 * sqrt(real(nEnsGain - 1,8)) *  &
               weightsMembersLatLon(:,memberIndex)

        end do

        ! Remove the weights mean computed over the columns
        do memberIndex = 1, nEnsGain
          weightsMembersLatLon(memberIndex,:) =  &
               weightsMembersLatLon(memberIndex,:) - &
               sum(weightsMembersLatLon(memberIndex,:))/real(nEns,8)
        end do

      else if (trim(algorithm) == 'CVLETKF') then
        !
        ! Weight calculation for cross-validation LETKF algorithm
        !

        ! Compute eigenValues/Vectors of Yb^T R^-1 Yb = E * Lambda * E^T
        call utl_tmg_start(135,'------EigenDecomp')
        tolerance = 1.0D-50
        call utl_eigenDecomp(YbTinvRYb_mean, eigenValues_mean, eigenVectors_mean, tolerance, matrixRank)
        call utl_tmg_stop(135)

        ! Compute ensemble mean local weights as E * (Lambda + (Nens-1)*I)^-1 * E^T * YbTinvR * (obs - meanYb)
        weightsTemp(:) = 0.0d0
        do localObsIndex = 1, numLocalObs
          bodyIndex = localBodyIndices(localObsIndex)
          do memberIndex = 1, nEns
            weightsTemp(memberIndex) = weightsTemp(memberIndex) +   &
                                       YbTinvR_mean(memberIndex,localObsIndex) *  &
                                       ( ensObs_mpiglobal%obsValue(bodyIndex) - &
                                         ensObs_mpiglobal%meanYb(bodyIndex) )
          end do
        end do
        weightsTemp2(:) = 0.0d0
        do memberIndex2 = 1, matrixRank
          do memberIndex1 = 1, nEns
            weightsTemp2(memberIndex2) = weightsTemp2(memberIndex2) +   &
                                         eigenVectors_mean(memberIndex1,memberIndex2) *  &
                                         weightsTemp(memberIndex1)
          end do
        end do
        do memberIndex = 1, matrixRank
          weightsTemp2(memberIndex) = weightsTemp2(memberIndex) *  &
                                      1.0D0/(eigenValues_mean(memberIndex) + real(nEns - 1,8))
        end do
        weightsMeanLatLon(:,1) = 0.0d0
        do memberIndex2 = 1, matrixRank
          do memberIndex1 = 1, nEns
            weightsMeanLatLon(memberIndex1,1) =  &
                 weightsMeanLatLon(memberIndex1,1) +   &
                 eigenVectors_mean(memberIndex1,memberIndex2) *  &
                 weightsTemp2(memberIndex2)
          end do
        end do

        ! Compute ensemble perturbation weights: 
        ! Wa = [ I - (Nens-1)^1/2 * E * 
        !        {(Nens-1)^-1/2*I - (Lambda + (Nens-1)*I)^-1/2} * Lambda^-1 *
        !        E^T * YbTinvRYb ]
        ! Loop over sub-ensembles

        !$OMP PARALLEL DO PRIVATE(subEnsIndex, memberIndexCV, memberIndexCV1, memberIndexCV2, &
        !$OMP                     memberIndex, memberIndex1, memberIndex2, weightsTemp, tolerance, &
        !$OMP                     YbTinvRYb_CV, eigenValues_CV, eigenVectors_CV, matrixRank)
        do subEnsIndex = 1, numSubEns

          ! Use complement (independent) ens to get eigenValues/Vectors of Yb^T R^-1 Yb = E*Lambda*E^T
          call utl_tmg_start(135,'------EigenDecomp')
          do memberIndexCV2 = 1, nEnsIndependentPerSubEns
            memberIndex2 = memberIndexSubEnsComp(memberIndexCV2, subEnsIndex)
            do memberIndexCV1 = 1, nEnsIndependentPerSubEns
              memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
              YbTinvRYb_CV(memberIndexCV1,memberIndexCV2) = YbTinvRYb_pert(memberIndex1,memberIndex2)
            end do
          end do
          tolerance = 1.0D-50
          call utl_eigenDecomp(YbTinvRYb_CV, eigenValues_CV, eigenVectors_CV, tolerance, matrixRank)
          call utl_tmg_stop(135)

          ! Loop over members within the current sub-ensemble being updated
          do memberIndexCV = 1, nEnsPerSubEns

            ! This is index of member being updated
            memberIndex = memberIndexSubEns(memberIndexCV, subEnsIndex)

            ! E^T * YbTinvRYb
            weightsTemp(:) = 0.0d0
            do memberIndex2 = 1, matrixRank
              do memberIndexCV1 = 1, nEnsIndependentPerSubEns
                memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
                weightsTemp(memberIndex2) = weightsTemp(memberIndex2) +  &
                                            eigenVectors_CV(memberIndexCV1,memberIndex2) *  &
                                            YbTinvRYb_pert(memberIndex1,memberIndex)
              end do
            end do

            ! {(Nens-1)^-1/2*I - (Lambda + (Nens-1)*I)^-1/2} Lambda^-1 * previous_result

            do memberIndex1 = 1, matrixRank
              weightsTemp(memberIndex1) = weightsTemp(memberIndex1) *  &
                                          ( 1.0D0/sqrt(real(nEnsIndependentPerSubEns - 1,8)) -   &
                                            1.0D0/sqrt(eigenValues_CV(memberIndex1) +  &
                                                      real(nEnsIndependentPerSubEns - 1,8)) )
              weightsTemp(memberIndex1) = weightsTemp(memberIndex1) /  &
                                          eigenValues_CV(memberIndex1)
            end do

            ! E * previous_result
            weightsMembersLatLon(:,memberIndex) = 0.0d0
            do memberIndex2 = 1, matrixRank
              do memberIndexCV1 = 1, nEnsIndependentPerSubEns
                memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
                weightsMembersLatLon(memberIndex1,memberIndex) =   &
                     weightsMembersLatLon(memberIndex1,memberIndex) +   &
                     eigenVectors_CV(memberIndexCV1,memberIndex2) *  &
                     weightsTemp(memberIndex2)
              end do
            end do

            ! -1 * (Nens-1)^1/2 * previous_result
            weightsMembersLatLon(:,memberIndex) =  &
                 -1.0D0 * sqrt(real(nEnsIndependentPerSubEns - 1,8)) *  &
                 weightsMembersLatLon(:,memberIndex)

            ! I + previous_result
            weightsMembersLatLon(memberIndex,memberIndex) =  &
                 1.0D0 + weightsMembersLatLon(memberIndex,memberIndex)

          end do ! memberIndexCV
        end do ! subEnsIndex
        !$OMP END PARALLEL DO

        ! Remove the weights mean computed over the columns
        do memberIndex = 1, nEns
          weightsMembersLatLon(memberIndex,:) =  &
               weightsMembersLatLon(memberIndex,:) - &
               sum(weightsMembersLatLon(memberIndex,:))/real(nEns,8)
        end do

      else if (trim(algorithm) == 'CVLETKF-ME') then
        !
        ! Weight calculation for cross-validation LETKF algorithm
        !

        ! Compute eigenValues/Vectors of Yb^T R^-1 Yb = E * Lambda * E^T
        call utl_tmg_start(135,'------EigenDecomp')
        tolerance = 1.0D-50
        call utl_eigenDecomp(YbTinvRYb_mean, eigenValues_mean, eigenVectors_mean, tolerance, matrixRank)
        call utl_tmg_stop(135)
        !if (matrixRank < (nEns-1)) then
        !  write(*,*) 'YbTinvRYb is rank deficient =', matrixRank, nEns, numLocalObs
        !end if

        ! Compute ensemble mean local weights as E * (Lambda + (Nens-1)*I)^-1 * E^T * YbTinvR * (obs - meanYb)
        weightsTemp(:) = 0.0d0
        do localObsIndex = 1, numLocalObs
          bodyIndex = localBodyIndices(localObsIndex)
          do memberIndex = 1, nEnsGain
            weightsTemp(memberIndex) = weightsTemp(memberIndex) +   &
                                       YbTinvR_mean(memberIndex,localObsIndex) *  &
                                       ( ensObs_mpiglobal%obsValue(bodyIndex) - &
                                         ensObs_mpiglobal%meanYb(bodyIndex) )
          end do
        end do
        weightsTemp2(:) = 0.0d0
        do memberIndex2 = 1, matrixRank
          do memberIndex1 = 1, nEnsGain
            weightsTemp2(memberIndex2) = weightsTemp2(memberIndex2) +   &
                                         eigenVectors_mean(memberIndex1,memberIndex2) *  &
                                         weightsTemp(memberIndex1)
          end do
        end do
        do memberIndex = 1, matrixRank
          weightsTemp2(memberIndex) = weightsTemp2(memberIndex) *  &
                                      1.0D0/(eigenValues_mean(memberIndex) + real(nEnsGain - 1,8))
        end do
        weightsMeanLatLon(:,1) = 0.0d0
        do memberIndex2 = 1, matrixRank
          do memberIndex1 = 1, nEnsGain
            weightsMeanLatLon(memberIndex1,1) =  &
                 weightsMeanLatLon(memberIndex1,1) +   &
                 eigenVectors_mean(memberIndex1,memberIndex2) *  &
                 weightsTemp2(memberIndex2)
          end do
        end do

        ! Compute ensemble perturbation weights: 
        ! Wa = [ - (Nens-1)^1/2 * E *
        !        {(Nens-1)^-1/2*I - (Lambda + (Nens-1)*I)^-1/2} * Lambda^-1 *
        !        E^T * YbTinvRYb_mod ]
        ! Loop over sub-ensembles
        !$OMP PARALLEL DO PRIVATE(subEnsIndex, memberIndexCV, memberIndexCV1, memberIndexCV2, &
        !$OMP                     memberIndex, memberIndex1, memberIndex2, weightsTemp, tolerance, &
        !$OMP                     YbTinvRYb_CV, eigenValues_CV, eigenVectors_CV, matrixRank)
        do subEnsIndex = 1, numSubEns

          ! Use complement (independent) ens to get eigenValues/Vectors of Yb^T R^-1 Yb = E*Lambda*E^T
          call utl_tmg_start(135,'------EigenDecomp')
          do memberIndexCV2 = 1, nEnsIndependentPerSubEns
            memberIndex2 = memberIndexSubEnsComp(memberIndexCV2, subEnsIndex)
            do memberIndexCV1 = 1, nEnsIndependentPerSubEns
              memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
              YbTinvRYb_CV(memberIndexCV1,memberIndexCV2) = YbTinvRYb_pert(memberIndex1,memberIndex2)
            end do
          end do
          tolerance = 1.0D-50
          call utl_eigenDecomp(YbTinvRYb_CV, eigenValues_CV, eigenVectors_CV, tolerance, matrixRank)
          call utl_tmg_stop(135)

          ! Loop over members within the current sub-ensemble being updated
          do memberIndexCV = 1, nEnsPerSubEns

            ! This is index of member being updated
            memberIndex = memberIndexSubEns(memberIndexCV, subEnsIndex)

            ! E^T * YbTinvRYb
            weightsTemp(:) = 0.0d0
            do memberIndex2 = 1, matrixRank
              do memberIndexCV1 = 1, nEnsIndependentPerSubEns
                memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
                weightsTemp(memberIndex2) = weightsTemp(memberIndex2) +  &
                                            eigenVectors_CV(memberIndexCV1,memberIndex2) *  &
                                            YbTinvRYb_mod(memberIndex1,memberIndex)
              end do
            end do

            ! {(Nens-1)^-1/2*I - (Lambda + (Nens-1)*I)^-1/2} Lambda^-1 * previous_result

            do memberIndex1 = 1, matrixRank
              weightsTemp(memberIndex1) = weightsTemp(memberIndex1) *  &
                                          ( 1.0D0/sqrt(real(nEnsIndependentPerSubEns - 1,8)) -   &
                                            1.0D0/sqrt(eigenValues_CV(memberIndex1) +  &
                                                       real(nEnsIndependentPerSubEns - 1,8)) )
              weightsTemp(memberIndex1) = weightsTemp(memberIndex1) /  &
                                          eigenValues_CV(memberIndex1)
            end do

            ! E * previous_result
            weightsMembersLatLon(:,memberIndex) = 0.0d0
            do memberIndex2 = 1, matrixRank
              do memberIndexCV1 = 1, nEnsIndependentPerSubEns
                memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
                weightsMembersLatLon(memberIndex1,memberIndex) =   &
                     weightsMembersLatLon(memberIndex1,memberIndex) +   &
                     eigenVectors_CV(memberIndexCV1,memberIndex2) *  &
                     weightsTemp(memberIndex2)
              end do
            end do

            ! -1 * (Nens-1)^1/2 * previous_result
            weightsMembersLatLon(:,memberIndex) =  &
                 -1.0D0 * sqrt(real(nEnsIndependentPerSubEns - 1,8)) *  &
                 weightsMembersLatLon(:,memberIndex)

          end do ! memberIndexCV
        end do ! subEnsIndex
        !$OMP END PARALLEL DO

        ! Remove the weights mean computed over the columns
        do memberIndex = 1, nEnsGain
          weightsMembersLatLon(memberIndex,:) =  &
               weightsMembersLatLon(memberIndex,:) - &
               sum(weightsMembersLatLon(memberIndex,:))/real(nEns,8)
        end do

      else if (trim(algorithm) == 'CVLETKF-PERTOBS') then
        !
        ! Weight calculation for perturbed-obs cross-validation LETKF algorithm
        !

        ! Compute eigenValues/Vectors of Yb^T R^-1 Yb = E * Lambda * E^T
        call utl_tmg_start(135,'------EigenDecomp')
        tolerance = 1.0D-50
        call utl_eigenDecomp(YbTinvRYb_mean, eigenValues_mean, eigenVectors_mean, tolerance, matrixRank)
        call utl_tmg_stop(135)

        ! Compute ensemble mean local weights as E * (Lambda + (Nens-1)*I)^-1 * E^T * YbTinvR * (obs - meanYb)
        weightsTemp(:) = 0.0d0
        do localObsIndex = 1, numLocalObs
          bodyIndex = localBodyIndices(localObsIndex)
          do memberIndex = 1, nEns
            weightsTemp(memberIndex) = weightsTemp(memberIndex) +   &
                                       YbTinvR_mean(memberIndex,localObsIndex) *  &
                                       ( ensObs_mpiglobal%obsValue(bodyIndex) - &
                                         ensObs_mpiglobal%meanYb(bodyIndex) )
          end do
        end do
        weightsTemp2(:) = 0.0d0
        do memberIndex2 = 1, matrixRank
          do memberIndex1 = 1, nEns
            weightsTemp2(memberIndex2) = weightsTemp2(memberIndex2) +   &
                                         eigenVectors_mean(memberIndex1,memberIndex2) *  &
                                         weightsTemp(memberIndex1)
          end do
        end do
        do memberIndex = 1, matrixRank
          weightsTemp2(memberIndex) = weightsTemp2(memberIndex) *  &
                                      1.0D0/(eigenValues_mean(memberIndex) + real(nEns - 1,8))
        end do
        weightsMeanLatLon(:,1) = 0.0d0
        do memberIndex2 = 1, matrixRank
          do memberIndex1 = 1, nEns
            weightsMeanLatLon(memberIndex1,1) =  &
                 weightsMeanLatLon(memberIndex1,1) +   &
                 eigenVectors_mean(memberIndex1,memberIndex2) *  &
                 weightsTemp2(memberIndex2)
          end do
        end do

        ! Compute ensemble perturbation weights using mean increment weights 
        ! formula, but with subset of members: 
        ! wa_i = I_i + E * (Lambda + (Nens-1)*I)^-1 * E^T * YbTinvR * (obs + randpert_i - Yb_i)
        ! Wa   = wa_i - mean_over_i(wa_i) 
        !
        ! Loop over sub-ensembles
        !$OMP PARALLEL DO PRIVATE(subEnsIndex, memberIndexCV, memberIndexCV1, memberIndexCV2, &
        !$OMP                     memberIndex, memberIndex1, memberIndex2, weightsTemp, tolerance, &
        !$OMP                     YbTinvRYb_CV, eigenValues_CV, eigenVectors_CV, matrixRank, &
        !$OMP                     weightsTemp2, localObsIndex, bodyIndex)
        do subEnsIndex = 1, numSubEns

          ! Use complement (independent) ens to get eigenValues/Vectors of Yb^T R^-1 Yb = E*Lambda*E^T
          call utl_tmg_start(135,'------EigenDecomp')
          do memberIndexCV2 = 1, nEnsIndependentPerSubEns
            memberIndex2 = memberIndexSubEnsComp(memberIndexCV2, subEnsIndex)
            do memberIndexCV1 = 1, nEnsIndependentPerSubEns
              memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
              YbTinvRYb_CV(memberIndexCV1,memberIndexCV2) = YbTinvRYb_pert(memberIndex1,memberIndex2)
            end do
          end do
          tolerance = 1.0D-50
          call utl_eigenDecomp(YbTinvRYb_CV, eigenValues_CV, eigenVectors_CV, tolerance, matrixRank)
          call utl_tmg_stop(135)

          ! Loop over members within the current sub-ensemble being updated
          do memberIndexCV = 1, nEnsPerSubEns

            ! This is index of member being updated (i'th member)
            memberIndex = memberIndexSubEns(memberIndexCV, subEnsIndex)

            ! YbTinvRYb * (obsValue + randPert_i - Yb_i)
            weightsTemp(:) = 0.0d0
            do localObsIndex = 1, numLocalObs
              bodyIndex = localBodyIndices(localObsIndex)
              do memberIndexCV1 = 1, nEnsIndependentPerSubEns
                memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
                weightsTemp(memberIndexCV1) =  & 
                     weightsTemp(memberIndexCV1) +   &
                     YbTinvR_pert(memberIndex1,localObsIndex) *  &
                     ( ensObs_mpiglobal%obsValue(bodyIndex) +  &
                       ensObs_mpiglobal%randPert_r4(memberIndex,bodyIndex) -  &
                       ( ensObs_mpiglobal%meanYb(bodyIndex) +  &
                         ensObs_mpiglobal%Yb_r4(memberIndex,bodyIndex) ) )
              end do
            end do

            ! E^T * previous_result
            weightsTemp2(:) = 0.0d0
            do memberIndex2 = 1, matrixRank
              do memberIndex1 = 1, nEnsIndependentPerSubEns
                weightsTemp2(memberIndex2) = weightsTemp2(memberIndex2) +   &
                                             eigenVectors_CV(memberIndex1,memberIndex2) *  &
                                             weightsTemp(memberIndex1)
              end do
            end do

            ! [lambda + (N_indep-1)*I]^-1 * previous_result
            do memberIndex1 = 1, matrixRank
              weightsTemp2(memberIndex1) =  &
                   weightsTemp2(memberIndex1) *  &
                   1.0D0/(eigenValues_CV(memberIndex1) + real(nEnsIndependentPerSubEns - 1,8))
            end do

            ! E * previous_result
            weightsMembersLatLon(:,memberIndex) = 0.0d0
            do memberIndex2 = 1, matrixRank
              do memberIndexCV1 = 1, nEnsIndependentPerSubEns
                memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
                weightsMembersLatLon(memberIndex1,memberIndex) =  &
                     weightsMembersLatLon(memberIndex1,memberIndex) +   &
                     eigenVectors_CV(memberIndexCV1,memberIndex2) *  &
                     weightsTemp2(memberIndex2)
              end do
            end do

            ! I + previous_result
            weightsMembersLatLon(memberIndex,memberIndex) =  &
                 1.0D0 + weightsMembersLatLon(memberIndex,memberIndex)

          end do ! memberIndexCV
        end do ! subEnsIndex
        !$OMP END PARALLEL DO

        ! Remove the weights mean computed over the columns
        do memberIndex = 1, nEns
          weightsMembersLatLon(memberIndex,:) =  &
               weightsMembersLatLon(memberIndex,:) - &
               sum(weightsMembersLatLon(memberIndex,:))/real(nEns,8)
        end do

      else

        call utl_abort('UNKNOWN LETKF ALGORITHM')

      end if

    else

      ! no obs near this grid point, mean weights zero, member weights identity
      weightsMeanLatLon(:,1) = 0.0d0
      weightsMembersLatLon(:,:) = 0.0d0
      do memberIndex = 1, nEns
        if ( useModulatedEns ) then
          do eigenVectorColumnIndex = 1, numRetainedEigen 
            memberIndexInModEns = (eigenVectorColumnIndex - 1) * nEns + memberIndex
            weightsMembersLatLon(memberIndexInModEns,memberIndex) = 1.0d0
          end do
        else
          weightsMembersLatLon(memberIndex,memberIndex) = 1.0d0
        end if
      end do

    end if ! numLocalObs > 0

    call utl_tmg_stop(134)

    firstCall = .false.

  end subroutine enkf_LETKFcomputeWeights

  !----------------------------------------------------------------------
  ! enkf_calcYbTinvRYb (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_calcYbTinvRYb(nEns1, nEns2, maxNumLocalObs, numLocalObs, &
                                YbTinvRYb_pert, YbTinvRCopy_pert, &
                                ensObs_mpiglobal, localBodyIndices,  &
                                YbTinvRYb_mean, YbTinvRCopy_mean)
    !
    !:Purpose: Compute the background covariance in ensemble space
    !
    implicit none

    ! Arguments:
    integer, intent(in)            :: nEns1
    integer, intent(in)            :: nEns2
    integer, intent(in)            :: maxNumLocalObs
    integer, intent(in)            :: numLocalObs
    type(struct_eob), intent(in)   :: ensObs_mpiglobal
    integer, intent(in)            :: localBodyIndices(maxNumLocalObs)
    real(8), intent(out)           :: YbTinvRYb_pert(nEns1,nEns2)
    real(8), intent(in)            :: YbTinvRCopy_pert(maxNumLocalObs,nEns1)
    real(8), intent(out), optional :: YbTinvRYb_mean(nEns1,nEns2)
    real(8), intent(in),  optional :: YbTinvRCopy_mean(maxNumLocalObs,nEns1)

    ! Locals:
    integer :: threadIndex, omp_get_thread_num
    integer :: memberIndex1, memberIndex2, localObsIndex, bodyIndex
    real(4), allocatable :: YbCopy2_r4(:,:)
    real(4), allocatable :: YbCopy_r4(:,:)

    allocate(YbCopy_r4(maxNumLocalObs,nEns2))
    allocate(YbCopy2_r4(numLocalObs,mmpi_numThread))

    call utl_tmg_start(137,'--------YbArraysCopy')
    !$OMP PARALLEL DO PRIVATE (localObsIndex, bodyIndex, memberIndex2)
    do localObsIndex = 1, numLocalObs
      bodyIndex = localBodyIndices(localObsIndex)
      do memberIndex2 = 1, nEns2
        YbCopy_r4(localObsIndex,memberIndex2) = ensObs_mpiglobal%Yb_r4(memberIndex2,bodyIndex)
      end do
    end do
    !$OMP END PARALLEL DO
    call utl_tmg_stop(137)

    if (nEns1 == nEns2) then
      ! When nEns1 equals nEns2 we assume matrix is symmetric

      !$OMP PARALLEL PRIVATE (memberIndex1, memberIndex2, threadIndex)
      threadIndex = 1 + omp_get_thread_num()
      !$OMP DO
      do memberIndex2 = 1, nEns2
        YbCopy2_r4(:,threadIndex) = YbCopy_r4(1:numLocalObs,memberIndex2)
        do memberIndex1 = 1, memberIndex2 ! compute only upper triangle
          YbTinvRYb_pert(memberIndex1,memberIndex2) = 0.0d0
          YbTinvRYb_pert(memberIndex1,memberIndex2) =  &
               YbTinvRYb_pert(memberIndex1,memberIndex2) +  &
               sum(YbTinvRCopy_pert(1:numLocalObs,memberIndex1) * YbCopy2_r4(1:numLocalObs,threadIndex))             
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      ! copy upper triangle to lower triangle (symmetric matrix)

      !$OMP PARALLEL DO PRIVATE (memberIndex1, memberIndex2)
      do memberIndex2 = 1, nEns2
        do memberIndex1 = memberIndex2+1, nEns1
          YbTinvRYb_pert(memberIndex1,memberIndex2) =  &
               YbTinvRYb_pert(memberIndex2,memberIndex1)
        end do
      end do
      !$OMP END PARALLEL DO

      if (eob_simObsAssim .and. present(YbTinvRYb_mean)) then     
        YbTinvRYb_mean(:,:) = 0.0D0
        !$OMP PARALLEL DO PRIVATE (memberIndex1, memberIndex2)
        do memberIndex2 = 1, nEns2
          do memberIndex1 = 1, memberIndex2 ! compute only upper triangle
            YbTinvRYb_mean(memberIndex1,memberIndex2) =  &
                 YbTinvRYb_mean(memberIndex1,memberIndex2) +  &
                 sum(YbTinvRCopy_mean(1:numLocalObs,memberIndex1) * YbCopy_r4(1:numLocalObs,memberIndex2))              
          end do
        end do
        !$OMP END PARALLEL DO
        ! copy upper triangle to lower triangle (symmetric matrix)
        do memberIndex2 = 1, nEns2
          do memberIndex1 = memberIndex2+1, nEns1
            YbTinvRYb_mean(memberIndex1,memberIndex2) =  &
                 YbTinvRYb_mean(memberIndex2,memberIndex1)
          end do
        end do
      end if

    else
      ! When nEns1 different from nEns2 we cannot assume matrix is symmetric

      !$OMP PARALLEL PRIVATE (memberIndex1, memberIndex2, threadIndex)
      threadIndex = 1 + omp_get_thread_num()
      !$OMP DO
      do memberIndex2 = 1, nEns2
        YbCopy2_r4(:,threadIndex) = YbCopy_r4(1:numLocalObs,memberIndex2)
        do memberIndex1 = 1, nEns1
          YbTinvRYb_pert(memberIndex1,memberIndex2) = 0.0d0
          YbTinvRYb_pert(memberIndex1,memberIndex2) =  &
               YbTinvRYb_pert(memberIndex1,memberIndex2) +  &
               sum(YbTinvRCopy_pert(1:numLocalObs,memberIndex1) * YbCopy2_r4(1:numLocalObs,threadIndex))             
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      if (eob_simObsAssim .and. present(YbTinvRYb_mean)) then     
        YbTinvRYb_mean(:,:) = 0.0D0
        !$OMP PARALLEL DO PRIVATE (memberIndex1, memberIndex2)
        do memberIndex2 = 1, nEns2
          do memberIndex1 = 1, nEns1
            YbTinvRYb_mean(memberIndex1,memberIndex2) =  &
                 YbTinvRYb_mean(memberIndex1,memberIndex2) +  &
                 sum(YbTinvRCopy_mean(1:numLocalObs,memberIndex1) * YbCopy_r4(1:numLocalObs,memberIndex2))              
          end do
        end do
        !$OMP END PARALLEL DO
      end if

    end if

    deallocate(YbCopy_r4)
    deallocate(YbCopy2_r4)

  end subroutine enkf_calcYbTinvRYb

  !----------------------------------------------------------------------
  ! enkf_computeVertLocation (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_computeVertLocation(vertLocation_r4,stateVectorMeanTrl)
    !
    !:Purpose:  Compute extract global 3D vertical location field from supplied
    !           stateVector. Can be either logPressure or depth levels.
    !
    implicit none

    ! Arguments:
    real(4), allocatable, intent(inout) :: vertLocation_r4(:,:,:)
    type(struct_gsv),     intent(inout) :: stateVectorMeanTrl

    ! Locals:
    integer          :: nLev_M, nLev_depth, nLev_vertLocation, levIndex, nsize, ierr
    real(4), pointer :: vertLocation_ptr_r4(:,:,:)
    type(struct_gsv) :: stateVectorMeanTrlPressure
    type(struct_gsv) :: stateVectorMeanTrlPressure_1step

    write(*,*) 'enkf_computeVertLocation: starting'

    nLev_M = gsv_getNumLev(stateVectorMeanTrl, 'MM')
    nLev_depth = gsv_getNumLev(stateVectorMeanTrl, 'DP')
    if ( nLev_M > 0 .and. nLev_depth > 0 ) then
      call utl_abort('enkf_computeVertLocation: both momentum and depth levels exist.')
    else if ( nLev_M == 0 .and. nLev_depth == 0 ) then
      call utl_abort('enkf_computeVertLocation: neither momentum nor depth levels exist.')
    end if
    nLev_vertLocation = max(nLev_M, nLev_depth)

    allocate(vertLocation_r4(stateVectorMeanTrl%hco%ni, &
                             stateVectorMeanTrl%hco%nj, &
                             nLev_vertLocation))

    if ( nLev_M > 0 ) then ! log pressure for NWP fields

      ! Compute background ens mean 3D log pressure and make mpiglobal for vertical localization
      call gsv_allocate( stateVectorMeanTrlPressure, tim_nstepobsinc,  &
                         stateVectorMeanTrl%hco, stateVectorMeanTrl%vco, dateStamp_opt=tim_getDateStamp(),  &
                         mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                         dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0','P_M','P_T'/) )
      call gsv_zero(stateVectorMeanTrlPressure)
      call gsv_copy(stateVectorMeanTrl, stateVectorMeanTrlPressure, allowVarMismatch_opt=.true.)
      call gvt_transform(stateVectorMeanTrlPressure,'ZandP_nl')
      if (mmpi_myid == 0) then
        call gsv_allocate( stateVectorMeanTrlPressure_1step, 1,  &
                           stateVectorMeanTrl%hco, stateVectorMeanTrl%vco, dateStamp_opt=tim_getDateStamp(),  &
                           mpi_local_opt=.false., &
                           dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0','P_M','P_T'/) )
      end if
      call gsv_transposeTilesToStep(stateVectorMeanTrlPressure_1step, stateVectorMeanTrlPressure, (tim_nstepobsinc+1)/2)
      call gsv_deallocate(stateVectorMeanTrlPressure)
      if (mmpi_myid == 0) then
        call gsv_getField(stateVectorMeanTrlPressure_1step,vertLocation_ptr_r4,'P_M')
        vertLocation_r4(:,:,:) = log(vertLocation_ptr_r4(:,:,:))
      end if
      nsize = stateVectorMeanTrlPressure%ni * stateVectorMeanTrlPressure%nj * nLev_M
      call rpn_comm_bcast(vertLocation_r4, nsize, 'mpi_real4', 0, 'GRID', ierr)

    else if ( nLev_depth > 0 ) then ! depth for ocean fields

      ! fill in all horizontal grid points with the same profile of depth values
      do levIndex = 1, nLev_depth
        write(*,*) 'setting vertLocation for levIndex =', levIndex, &
                   ', depth = ', stateVectorMeanTrl%vco%depths(levIndex)
        vertLocation_r4(:,:,levIndex) = stateVectorMeanTrl%vco%depths(levIndex)
      end do

    end if

    write(*,*) 'enkf_computeVertLocation: finished'

  end subroutine enkf_computeVertLocation

  !----------------------------------------------------------------------
  ! enkf_LETKFsetupMpiDistribution (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_LETKFsetupMpiDistribution(numLatLonMpiGlobal, myNumLatLonRecv,  &
                                            myLatIndexesRecv, myLonIndexesRecv, &
                                            latIndexesSendMpiGlobal, lonIndexesSendMpiGlobal, &
                                            procIndexesSendMpiGlobal, &
                                            numProcsSendMpiGlobal, wInterpInfo)
    !
    ! :Purpose: Setup for distribution of grid points over mpi tasks by building a global
    !           list of which MPI tasks need the computed weights for which grid points.
    !           Also build local list for grid points required to recv.
    !
    implicit none

    ! Arguments:
    integer,                     intent(out) :: numLatLonMpiGlobal
    integer,                     intent(out) :: myNumLatLonRecv
    integer, allocatable,        intent(out) :: myLatIndexesRecv(:)
    integer, allocatable,        intent(out) :: myLonIndexesRecv(:)
    integer, allocatable,        intent(out) :: latIndexesSendMpiGlobal(:)
    integer, allocatable,        intent(out) :: lonIndexesSendMpiGlobal(:)
    integer, allocatable,        intent(out) :: procIndexesSendMpiGlobal(:,:)
    integer, allocatable,        intent(out) :: numProcsSendMpiGlobal(:)
    type(struct_enkfInterpInfo), intent(in)  :: wInterpInfo

    ! Locals:
    integer :: latIndex, lonIndex, procIndex, latLonIndex, myLatLonIndex, latLonIndexMpiGlobal
    integer :: myLonBeg, myLonEnd, myLatBeg, myLatEnd
    integer :: myLonBegHalo, myLonEndHalo, myLatBegHalo, myLatEndHalo
    integer :: numLatLonRecvMax, myNumLatLon, numLatLonMax, ierr
    integer, allocatable :: allLatIndexesRecv(:,:), allLonIndexesRecv(:,:)
    integer, allocatable :: allLatIndexesSend(:,:), allLonIndexesSend(:,:)
    integer, allocatable :: localLatIndexesSend(:), localLonIndexesSend(:)
    integer              :: allNumLatLonRecv(mmpi_nprocs), allNumLatLon(mmpi_nprocs)

    myLonBeg = wInterpInfo%myLonBeg
    myLonEnd = wInterpInfo%myLonEnd
    myLatBeg = wInterpInfo%myLatBeg
    myLatEnd = wInterpInfo%myLatEnd

    myLonBegHalo = wInterpInfo%myLonBegHalo
    myLonEndHalo = wInterpInfo%myLonEndHalo
    myLatBegHalo = wInterpInfo%myLatBegHalo
    myLatEndHalo = wInterpInfo%myLatEndHalo

    write(*,*) 'enkf_LETKFsetupMpiDistribution: starting'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! First, count the number of grid points where weights needed locally (for recv-ing)
    ! Note, this includes the "halo" needed for interpolation
    myNumLatLonRecv = 0
    do latIndex = myLatBegHalo, myLatEndHalo
      LON_LOOP1: do lonIndex = myLonBegHalo, myLonEndHalo
        ! If this lat-lon is to be interpolated, then skip calculation
        if (wInterpInfo%numIndexes(lonIndex,latIndex) > 0) cycle LON_LOOP1
        myNumLatLonRecv = myNumLatLonRecv + 1
      end do LON_LOOP1
    end do

    ! Communicate to all mpi tasks
    call rpn_comm_allgather(myNumLatLonRecv, 1, "mpi_integer",  &
                            allNumLatLonRecv, 1,"mpi_integer", "GRID", ierr)
    numLatLonRecvMax = maxval(allNumLatLonRecv)
    write(*,*) 'enkf_LETKFsetupMpiDistribution: allNumLatLonRecv =', allNumLatLonRecv(:)
    write(*,*) 'enkf_LETKFsetupMpiDistribution: numLatLonRecvSum =', sum(allNumLatLonRecv)
    write(*,*) 'enkf_LETKFsetupMpiDistribution: numLatLonRecvMax =', numLatLonRecvMax

    ! Now create a list of grid point indexes where weights needed locally (for recv-ing)
    allocate(myLatIndexesRecv(numLatLonRecvMax))
    allocate(myLonIndexesRecv(numLatLonRecvMax))
    myLatIndexesRecv(:) = -1
    myLonIndexesRecv(:) = -1
    myNumLatLonRecv = 0
    do latIndex = myLatBegHalo, myLatEndHalo
      LON_LOOP2: do lonIndex = myLonBegHalo, myLonEndHalo
        ! If this lat-lon is to be interpolated, then skip calculation
        if (wInterpInfo%numIndexes(lonIndex,latIndex) > 0) cycle LON_LOOP2
        myNumLatLonRecv = myNumLatLonRecv + 1

        myLatIndexesRecv(myNumLatLonRecv) = latIndex
        myLonIndexesRecv(myNumLatLonRecv) = lonIndex
      end do LON_LOOP2
    end do

    ! Communicate to all mpi tasks this list of grid point lat-lon indexes
    allocate(allLatIndexesRecv(numLatLonRecvMax, mmpi_nprocs))
    allocate(allLonIndexesRecv(numLatLonRecvMax, mmpi_nprocs))
    call rpn_comm_allgather(myLatIndexesRecv, numLatLonRecvMax, "mpi_integer",  &
                            allLatIndexesRecv, numLatLonRecvMax, "mpi_integer",  &
                            "GRID", ierr)
    call rpn_comm_allgather(myLonIndexesRecv, numLatLonRecvMax, "mpi_integer",  &
                            allLonIndexesRecv, numLatLonRecvMax, "mpi_integer",  &
                            "GRID", ierr)

    ! Now count number of local grid points without the halo
    myNumLatLon = 0
    do latIndex = myLatBeg, myLatEnd
      LON_LOOP3: do lonIndex = myLonBeg, myLonEnd
        ! If this lat-lon is to be interpolated, then skip calculation
        if (wInterpInfo%numIndexes(lonIndex,latIndex) > 0) cycle LON_LOOP3
        myNumLatLon = myNumLatLon + 1
      end do LON_LOOP3
    end do

    ! Communicate to all mpi tasks
    call rpn_comm_allreduce(myNumLatLon, numLatLonMpiGlobal, &
                            1,"mpi_integer","mpi_sum","GRID",ierr)
    call rpn_comm_allgather(myNumLatLon,  1, "mpi_integer",  &
                            allNumLatLon, 1, "mpi_integer",  &
                            "GRID", ierr)
    numLatLonMax = maxval(allNumLatLon)

    ! Build global lists of lat-lon indexes and list of mpi tasks where each needs to be sent
    allocate(allLatIndexesSend(numLatLonMax,mmpi_nprocs))
    allocate(allLonIndexesSend(numLatLonMax,mmpi_nprocs))
    allocate(latIndexesSendMpiGlobal(numLatLonMpiGlobal))
    allocate(lonIndexesSendMpiGlobal(numLatLonMpiGlobal))
    allocate(localLatIndexesSend(numLatLonMax))
    allocate(localLonIndexesSend(numLatLonMax))
    latIndexesSendMpiGlobal(:) = -1
    lonIndexesSendMpiGlobal(:) = -1
    localLatIndexesSend(:) = -1
    localLonIndexesSend(:) = -1

    myNumLatLon = 0
    do latIndex = myLatBeg, myLatEnd
      LON_LOOP4: do lonIndex = myLonBeg, myLonEnd
        ! If this lat-lon is to be interpolated, then skip calculation
        if (wInterpInfo%numIndexes(lonIndex,latIndex) > 0) cycle LON_LOOP4
        myNumLatLon = myNumLatLon + 1

        localLatIndexesSend(myNumLatLon) = latIndex
        localLonIndexesSend(myNumLatLon) = lonIndex
      end do LON_LOOP4
    end do

    call rpn_comm_allgather(localLatIndexesSend, numLatLonMax, "mpi_integer",  &
                            allLatIndexesSend,   numLatLonMax, "mpi_integer",  &
                            "GRID", ierr)
    call rpn_comm_allgather(localLonIndexesSend, numLatLonMax, "mpi_integer",  &
                            allLonIndexesSend,   numLatLonMax, "mpi_integer",  &
                            "GRID", ierr)

    ! Reorganize into single dimension list
    latLonIndexMpiGlobal = 0
    do procIndex = 1, mmpi_nprocs
      do myLatLonIndex = 1, allNumLatLon(procIndex)
        latLonIndexMpiGlobal = latLonIndexMpiGlobal + 1

        latIndexesSendMpiGlobal(latLonIndexMpiGlobal) = allLatIndexesSend(myLatLonIndex, procIndex)
        lonIndexesSendMpiGlobal(latLonIndexMpiGlobal) = allLonIndexesSend(myLatLonIndex, procIndex)
      end do
    end do

    
    ! Figure out which mpi tasks I will need to send my results to
    allocate(procIndexesSendMpiGlobal(numLatLonMpiGlobal,10))
    allocate(numProcsSendMpiGlobal(numLatLonMpiGlobal))
    procIndexesSendMpiGlobal(:,:) = -1
    numProcsSendMpiGlobal(:) = 0
    do latLonIndex = 1, numLatLonMpiGlobal
      do procIndex = 1, mmpi_nprocs
        if ( any( (latIndexesSendMpiGlobal(latLonIndex) == allLatIndexesRecv(1:allNumLatLonRecv(procIndex), procIndex)) .and.  &
                  (lonIndexesSendMpiGlobal(latLonIndex) == allLonIndexesRecv(1:allNumLatLonRecv(procIndex), procIndex)) ) ) then
          numProcsSendMpiGlobal(latLonIndex) = numProcsSendMpiGlobal(latLonIndex) + 1
          procIndexesSendMpiGlobal(latLonIndex,numProcsSendMpiGlobal(latLonIndex)) = procIndex
        end if
      end do
    end do

    write(*,*) 'enkf_LETKFsetupMpiDistribution: lat/lon/proc indexes I need to receive:'
    do latLonIndex = 1, myNumLatLonRecv
      write(*,*) myLatIndexesRecv(latLonIndex), myLonIndexesRecv(latLonIndex)
    end do

    if (mmpi_myid == 0) then
      write(*,*) 'enkf_LETKFsetupMpiDistribution: number of lat/lon indexes globally =', numLatLonMpiGlobal
      write(*,*) 'enkf_LETKFsetupMpiDistribution: global list of indexMpiGlobal/lat/lon/proc indexes:'
      do latLonIndexMpiGlobal = 1, numLatLonMpiGlobal
        write(*,*) latLonIndexMpiGlobal,  &
                   latIndexesSendMpiGlobal(latLonIndexMpiGlobal),  &
                   lonIndexesSendMpiGlobal(latLonIndexMpiGlobal),  &
                   procIndexesSendMpiGlobal(latLonIndexMpiGlobal,1:numProcsSendMpiGlobal(latLonIndexMpiGlobal))
      end do
    end if

    write(*,*) 'enkf_LETKFsetupMpiDistribution: done'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine enkf_LETKFsetupMpiDistribution

  !----------------------------------------------------------------------
  ! enkf_LETKFgetMpiGlobalTags (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_LETKFgetMpiGlobalTags(latLonTagMpiGlobal,myLatIndexesRecv,myLonIndexesRecv)
    implicit none

    ! Arguments:
    integer, intent(out) :: latLonTagMpiGlobal(:,:)
    integer, intent(in)  :: myLatIndexesRecv(:)
    integer, intent(in)  :: myLonIndexesRecv(:)

    ! Locals:
    integer :: ierr, ni, nj, lonIndex, latIndex
    integer :: countTags, myNumLatLonRecv, numLatLonRecvMax
    integer, allocatable :: allNumLatLonRecv(:)
    integer, allocatable :: allLatIndexesRecv(:,:), allLonIndexesRecv(:,:)

    write(*,*) 'enkf_LETKFgetMpiGlobalTags: Starting'

    ni = size(latLonTagMpiGlobal,1)
    nj = size(latLonTagMpiGlobal,2)

    myNumLatLonRecv = size(myLatIndexesRecv)
    allocate(allNumLatLonRecv(mmpi_nprocs))
    call rpn_comm_allgather(myNumLatLonRecv, 1, "mpi_integer",  &
                            allNumLatLonRecv, 1,"mpi_integer", "GRID", ierr)
    numLatLonRecvMax = maxval(allNumLatLonRecv)

    ! Communicate to all mpi tasks this list of grid point lat-lon indexes
    allocate(allLatIndexesRecv(numLatLonRecvMax, mmpi_nprocs))
    allocate(allLonIndexesRecv(numLatLonRecvMax, mmpi_nprocs))
    call rpn_comm_allgather(myLatIndexesRecv, numLatLonRecvMax, "mpi_integer",  &
                            allLatIndexesRecv, numLatLonRecvMax, "mpi_integer",  &
                            "GRID", ierr)
    call rpn_comm_allgather(myLonIndexesRecv, numLatLonRecvMax, "mpi_integer",  &
                            allLonIndexesRecv, numLatLonRecvMax, "mpi_integer",  &
                            "GRID", ierr)

    latLonTagMpiGlobal(:,:) = 0
    !$OMP PARALLEL DO PRIVATE(latIndex, lonIndex)
    do lonIndex = 1, ni
      do latIndex = 1, nj
        if (any(lonIndex == allLonIndexesRecv(:,:) .and. latIndex == allLatIndexesRecv(:,:))) then
          latLonTagMpiGlobal(lonIndex,latIndex) = 1
        end if
      end do
    end do
    !$OMP END PARALLEL DO

    countTags = 0
    do lonIndex = 1, ni
      do latIndex = 1, nj
        if (latLonTagMpiGlobal(lonIndex,latIndex) == 1) then
          countTags = countTags + 1
          latLonTagMpiGlobal(lonIndex,latIndex) = countTags
        end if
      end do
    end do
    write(*,*) 'number of Recv grid points found = ', maxval(latLonTagMpiGlobal(:,:))

    if (maxval(latLonTagmpiGlobal(:,:)) > (mmpi_maxTagValue - 2)) then
      write(*,*) 'maximum allowable tag value = ', mmpi_maxTagValue - 2
      call utl_abort('enkf_LETKFgetMpiGlobalTags: mpi tag values exceeded max allowable value')
    end if
    
    write(*,*) 'enkf_LETKFgetMpiGlobalTags: Finished'

  end subroutine enkf_LETKFgetMpiGlobalTags

  !----------------------------------------------------------------------
  ! enkf_latLonAlreadyFound (private function)
  !----------------------------------------------------------------------
  function enkf_latLonAlreadyFound(allLatIndexesRecv, allLonIndexesRecv, latLonIndex, procIndex) result(found)
    implicit none

    ! Arguments:
    integer, intent(in) :: allLatIndexesRecv(:,:)
    integer, intent(in) :: allLonIndexesRecv(:,:)
    integer, intent(in) :: latLonIndex
    integer, intent(in) :: procIndex
    ! Result:
    logical :: found

    ! Locals:
    integer :: latLonIndex2, procIndex2, numLatLonRecvMax

    numLatLonRecvMax = size(allLatIndexesRecv, 1)

    ! check on all previous mpi tasks if this lat/lon has already been encountered
    found = .false.
    do procIndex2 = 1, procIndex-1
      WEIGHTS1LEV_LOOP: do latLonIndex2 = 1, numLatLonRecvMax
        if (allLatIndexesRecv(latLonIndex2, procIndex2) < 0) cycle WEIGHTS1LEV_LOOP
        if ( (allLatIndexesRecv(latLonIndex, procIndex) ==  &
              allLatIndexesRecv(latLonIndex2, procIndex2)) .and.  &
             (allLonIndexesRecv(latLonIndex, procIndex) ==  &
              allLonIndexesRecv(latLonIndex2, procIndex2)) ) then
          found = .true.
          exit WEIGHTS1LEV_LOOP
        end if
      end do WEIGHTS1LEV_LOOP
    end do

  end function enkf_latLonAlreadyFound

  !--------------------------------------------------------------------------
  ! enkf_setupInterpInfo
  !--------------------------------------------------------------------------
  subroutine enkf_setupInterpInfo(wInterpInfo, hco, weightLatLonStep,  &
                                  myLonBeg,myLonEnd,myLatBeg,myLatEnd)
    !
    ! :Purpose: Setup the weights and lat/lon indices needed to bilinearly
    !           interpolate the LETKF weights from a coarse grid to the full
    !           resolution grid. The coarseness of the grid is specified by
    !           the weightLatLonStep argument.
    !
    implicit none

    ! Arguments:
    type(struct_enkfInterpInfo), intent(out) :: wInterpInfo
    type(struct_hco),            intent(in)  :: hco
    integer,                     intent(in)  :: weightLatLonStep
    integer,                     intent(in)  :: myLonBeg
    integer,                     intent(in)  :: myLonEnd
    integer,                     intent(in)  :: myLatBeg
    integer,                     intent(in)  :: myLatEnd

    ! Locals:
    integer :: lonIndex, latIndex, ni, nj
    integer :: myLonBegHalo, myLonEndHalo, myLatBegHalo, myLatEndHalo
    real(8) :: interpWeightLon, interpWeightLat
    logical :: includesYinYangBndry

    ni = hco%ni
    nj = hco%nj

    myLonBegHalo = 1 + weightLatLonStep * floor(real(myLonBeg - 1)/real(weightLatLonStep))
    myLonEndHalo = min(ni, 1 + weightLatLonStep * ceiling(real(myLonEnd - 1)/real(weightLatLonStep)))
    myLatBegHalo = 1 + weightLatLonStep * floor(real(myLatBeg - 1)/real(weightLatLonStep))
    myLatEndHalo = min(nj, 1 + weightLatLonStep * ceiling(real(myLatEnd - 1)/real(weightLatLonStep)))
    write(*,*) 'enkf_setupInterpInfo: myLonBeg/End, myLatBeg/End (original)  = ',  &
               myLonBeg, myLonEnd, myLatBeg, myLatEnd
    write(*,*) 'enkf_setupInterpInfo: myLonBeg/End, myLatBeg/End (with Halo) = ',  &
               myLonBegHalo, myLonEndHalo, myLatBegHalo, myLatEndHalo
    write(*,*) 'enkf_setupInterpInfo: myLonCount, myLatCount (with Halo) = ', &
               myLonEndHalo-myLonBegHalo+1, myLatEndHalo-myLatBegHalo+1
    write(*,*) 'enkf_setupInterpInfo: number of local gridpts where weights computed = ',  &
               ( 1 + ceiling(real(myLonEndHalo - myLonBegHalo) / real(weightLatLonStep)) ) *  &
               ( 1 + ceiling(real(myLatEndHalo - myLatBegHalo) / real(weightLatLonStep)) )
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    wInterpInfo%latLonStep   = weightLatLonStep
    wInterpInfo%myLonBegHalo = myLonBegHalo
    wInterpInfo%myLonEndHalo = myLonEndHalo
    wInterpInfo%myLatBegHalo = myLatBegHalo
    wInterpInfo%myLatEndHalo = myLatEndHalo
    wInterpInfo%myLonBeg = myLonBeg
    wInterpInfo%myLonEnd = myLonEnd
    wInterpInfo%myLatBeg = myLatBeg
    wInterpInfo%myLatEnd = myLatEnd

    allocate(wInterpInfo%numIndexes(myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
    if (weightLatLonStep > 1) then
      ! Figure out if this tile straddles Yin-Yang boundary
      if (hco%grtyp == 'U' .and. myLatBegHalo <= nj/2 .and. myLatEndHalo >= ((nj/2)+1)) then
        includesYinYangBndry = .true.
      else
        includesYinYangBndry = .false.
      end if
      allocate(wInterpInfo%lonIndexes(4,myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
      allocate(wInterpInfo%latIndexes(4,myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
      allocate(wInterpInfo%interpWeights(4,myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
      wInterpInfo%lonIndexes(:,:,:) = 0
      wInterpInfo%latIndexes(:,:,:) = 0
      wInterpInfo%interpWeights(:,:,:) = 0.0D0
      ! Determine which lat-lon are interpolated (wInterpInfo%numIndexes>0)
      wInterpInfo%numIndexes(:,:) = 4
      do latIndex = myLatBegHalo, myLatEndHalo, weightLatLonStep
        do lonIndex = myLonBegHalo, myLonEndHalo, weightLatLonStep
          wInterpInfo%numIndexes(lonIndex,latIndex) = 0
        end do
      end do
      ! Ensure weights are computed along edge of domain
      if (myLonEndHalo == ni .and.  &
          myLatEndHalo == nj) then
        wInterpInfo%numIndexes(myLonEndHalo,myLatEndHalo) = 0
      end if
      if (myLonEndHalo == ni) then
        do latIndex = myLatBegHalo, myLatEndHalo, weightLatLonStep
          wInterpInfo%numIndexes(myLonEndHalo,latIndex) = 0
        end do
        ! Ensure weights are computed along both sides of Yin-Yang boundary
        if (includesYinYangBndry) then
          wInterpInfo%numIndexes(ni,nj/2) = 0
          wInterpInfo%numIndexes(ni,(nj/2)+1) = 0
          write(*,*) 'enkf_setupInterpInfo: Yin-Yang boundary (lon,lat1,lat2) =',  &
                     ni, nj/2, (nj/2)+1
        end if
      end if
      if (myLatEndHalo == nj) then
        do lonIndex = myLonBegHalo, myLonEndHalo, weightLatLonStep
          wInterpInfo%numIndexes(lonIndex,myLatEndHalo) = 0
        end do
      end if
      ! Ensure weights are computed along both sides of Yin-Yang boundary
      if (includesYinYangBndry) then
        do lonIndex = myLonBegHalo, myLonEndHalo, weightLatLonStep
          wInterpInfo%numIndexes(lonIndex,nj/2) = 0
          wInterpInfo%numIndexes(lonIndex,(nj/2)+1) = 0
          write(*,*) 'enkf_setupInterpInfo: Yin-Yang boundary (lon,lat1,lat2) =',  &
                     lonIndex, nj/2, (nj/2)+1
        end do
      end if
      ! For lon-only interpolation
      do latIndex = myLatBegHalo, myLatEndHalo
        if (wInterpInfo%numIndexes(myLonBegHalo,latIndex) > 0) cycle
        do lonIndex = myLonBegHalo, myLonEndHalo
          if (wInterpInfo%numIndexes(lonIndex,latIndex) == 0) cycle
          ! Find nearest grid point with a value towards left
          wInterpInfo%numIndexes(lonIndex,latIndex) = 2
          wInterpInfo%lonIndexes(1,lonIndex,latIndex) = myLonBegHalo +  &
               weightLatLonStep * floor(real(lonIndex - myLonBegHalo)/real(weightLatLonStep)) 
          wInterpInfo%lonIndexes(2,lonIndex,latIndex) = min(ni,  &
               wInterpInfo%lonIndexes(1,lonIndex,latIndex) + weightLatLonStep)
          wInterpInfo%latIndexes(1,lonIndex,latIndex) = latIndex
          wInterpInfo%latIndexes(2,lonIndex,latIndex) = latIndex
          wInterpInfo%interpWeights(1,lonIndex,latIndex) =   &
               real(wInterpInfo%lonIndexes(2,lonIndex,latIndex) - lonIndex, 8)/real(weightLatLonStep, 8)
          wInterpInfo%interpWeights(2,lonIndex,latIndex) = 1.0D0 -  &
               wInterpInfo%interpWeights(1,lonIndex,latIndex)
        end do
      end do
      ! For lat-only interpolation
      do latIndex = myLatBegHalo, myLatEndHalo
        do lonIndex = myLonBegHalo, myLonEndHalo, weightLatLonStep
          if (wInterpInfo%numIndexes(lonIndex,myLatBegHalo) > 0) cycle
          if (wInterpInfo%numIndexes(lonIndex,latIndex) == 0) cycle
          ! Find nearest grid point with a value towards bottom
          wInterpInfo%numIndexes(lonIndex,latIndex) = 2
          wInterpInfo%lonIndexes(1,lonIndex,latIndex) = lonIndex
          wInterpInfo%lonIndexes(2,lonIndex,latIndex) = lonIndex
          wInterpInfo%latIndexes(1,lonIndex,latIndex) = myLatBegHalo +  &
               weightLatLonStep * floor(real(latIndex - myLatBegHalo)/real(weightLatLonStep)) 
          wInterpInfo%latIndexes(2,lonIndex,latIndex) = min(nj,  &
               wInterpInfo%latIndexes(1,lonIndex,latIndex) + weightLatLonStep)
          ! Ensure we do not interpolate values across Yin-Yang boundary
          if (includesYinYangBndry) then
            if (latIndex <= nj/2) then
              wInterpInfo%latIndexes(2,lonIndex,latIndex) = min(nj/2, wInterpInfo%latIndexes(2,lonIndex,latIndex))
            else if(latIndex >= (nj/2)+1) then
              wInterpInfo%latIndexes(1,lonIndex,latIndex) = max((nj/2)+1, wInterpInfo%latIndexes(1,lonIndex,latIndex))
            end if
          end if
          wInterpInfo%interpWeights(1,lonIndex,latIndex) =  &
               real(wInterpInfo%latIndexes(2,lonIndex,latIndex) - latIndex, 8)/real(weightLatLonStep, 8)
          wInterpInfo%interpWeights(2,lonIndex,latIndex) = 1.0D0 -  &
               wInterpInfo%interpWeights(1,lonIndex,latIndex)
        end do
      end do
      ! For interior points needing 2D interpolation
      do latIndex = myLatBegHalo, myLatEndHalo
        do lonIndex = myLonBegHalo, myLonEndHalo
          if (wInterpInfo%numIndexes(lonIndex,latIndex) == 0) cycle ! no interpolation
          if (wInterpInfo%lonIndexes(1,lonIndex,latIndex) /= 0) cycle ! already set up
          wInterpInfo%numIndexes(lonIndex,latIndex) = 4
          ! 1. bottom-left indexes
          wInterpInfo%lonIndexes(1,lonIndex,latIndex) = myLonBegHalo +  &
               weightLatLonStep * floor(real(lonIndex - myLonBegHalo)/real(weightLatLonStep)) 
          wInterpInfo%latIndexes(1,lonIndex,latIndex) = myLatBegHalo +  &
               weightLatLonStep * floor(real(latIndex - myLatBegHalo)/real(weightLatLonStep)) 
          ! 2. bottom-right indexes
          wInterpInfo%lonIndexes(2,lonIndex,latIndex) = min(ni,  &
               wInterpInfo%lonIndexes(1,lonIndex,latIndex) + weightLatLonStep)
          wInterpInfo%latIndexes(2,lonIndex,latIndex) = wInterpInfo%latIndexes(1,lonIndex,latIndex)
          ! 3. top-left indexes
          wInterpInfo%lonIndexes(3,lonIndex,latIndex) = wInterpInfo%lonIndexes(1,lonIndex,latIndex)
          wInterpInfo%latIndexes(3,lonIndex,latIndex) = min(nj,  &
               wInterpInfo%latIndexes(1,lonIndex,latIndex) + weightLatLonStep)
          ! 4. top-right indexes
          wInterpInfo%lonIndexes(4,lonIndex,latIndex) = wInterpInfo%lonIndexes(2,lonIndex,latIndex)
          wInterpInfo%latIndexes(4,lonIndex,latIndex) = wInterpInfo%latIndexes(3,lonIndex,latIndex)
          ! Ensure we do not interpolate values across Yin-Yang boundary
          if (includesYinYangBndry) then
            if (latIndex <= nj/2) then
              wInterpInfo%latIndexes(3,lonIndex,latIndex) = min(nj/2, wInterpInfo%latIndexes(3,lonIndex,latIndex))
              wInterpInfo%latIndexes(4,lonIndex,latIndex) = min(nj/2, wInterpInfo%latIndexes(4,lonIndex,latIndex))
            else if(latIndex >= (nj/2)+1) then
              wInterpInfo%latIndexes(1,lonIndex,latIndex) = max((nj/2)+1, wInterpInfo%latIndexes(1,lonIndex,latIndex))
              wInterpInfo%latIndexes(2,lonIndex,latIndex) = max((nj/2)+1, wInterpInfo%latIndexes(2,lonIndex,latIndex))
            end if
          end if
          ! one-dimensional weights in lon and lat directions
          interpWeightLon = real(wInterpInfo%lonIndexes(4,lonIndex,latIndex) - lonIndex, 8) /  &
                            real(weightLatLonStep, 8)
          interpWeightLat = real(wInterpInfo%latIndexes(4,lonIndex,latIndex) - latIndex, 8) /  &
                            real(weightLatLonStep, 8)
          ! four interpolation weights
          wInterpInfo%interpWeights(1,lonIndex,latIndex) = interpWeightLon * interpWeightLat
          wInterpInfo%interpWeights(2,lonIndex,latIndex) = (1.0D0 - interpWeightLon) * interpWeightLat
          wInterpInfo%interpWeights(3,lonIndex,latIndex) = interpWeightLon * (1.0D0 - interpWeightLat)
          wInterpInfo%interpWeights(4,lonIndex,latIndex) = (1.0D0 - interpWeightLon) * (1.0D0 - interpWeightLat)
        end do
      end do
    else
      ! no interpolation, all weights are computed
      wInterpInfo%numIndexes(:,:) = 0
    end if
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine enkf_setupInterpInfo

  !--------------------------------------------------------------------------
  ! enkf_interpWeights
  !--------------------------------------------------------------------------
  subroutine enkf_interpWeights(wInterpInfo, weights)
    !
    ! :Purpose: Perform the bilinear interpolation of the weights
    !           using the precalculated interpolation info.
    !
    implicit none

    ! Arguments:
    type(struct_enkfInterpInfo), intent(in)  :: wInterpInfo
    real(8),                     intent(out) :: weights(1:,1:,wInterpInfo%myLonBegHalo:,wInterpInfo%myLatBegHalo:)

    ! Locals:
    integer :: myLonBegHalo, myLonEndHalo, myLatBegHalo, myLatEndHalo
    integer :: lonIndex, latIndex, memberIndex1, memberIndex2, interpIndex
    integer :: interpLonIndex, interpLatIndex, numMembers1, numMembers2
    integer :: totalCount(mmpi_numthread)
    integer, external :: omp_get_thread_num
    logical, save :: firstCall = .true.

    myLonBegHalo = wInterpInfo%myLonBegHalo
    myLonEndHalo = wInterpInfo%myLonEndHalo
    myLatBegHalo = wInterpInfo%myLatBegHalo
    myLatEndHalo = wInterpInfo%myLatEndHalo
    numMembers1 = size(weights,1)
    numMembers2 = size(weights,2)
    totalCount(:) = 0

    !$OMP PARALLEL DO PRIVATE(latIndex, lonIndex, interpIndex, interpLatIndex, interpLonIndex, memberIndex1, memberIndex2)
    do latIndex = wInterpInfo%myLatBeg, wInterpInfo%myLatEnd
      do lonIndex = wInterpInfo%myLonBeg, wInterpInfo%myLonEnd
        if (wInterpInfo%numIndexes(lonIndex,latIndex) <= 0) cycle
        weights(:,:,lonIndex,latIndex) = 0.0D0
        if (wInterpInfo%lonIndexes(1,lonIndex,latIndex) == 0) cycle

        totalCount(omp_get_thread_num()+1) = totalCount(omp_get_thread_num()+1) + wInterpInfo%numIndexes(lonIndex,latIndex)

        do interpIndex = 1, wInterpInfo%numIndexes(lonIndex,latIndex)
          interpLonIndex = wInterpInfo%lonIndexes(interpIndex,lonIndex,latIndex)
          interpLatIndex = wInterpInfo%latIndexes(interpIndex,lonIndex,latIndex)

          do memberIndex2 = 1, numMembers2
            do memberIndex1 = 1, numMembers1
              weights(memberIndex1,memberIndex2,lonIndex,latIndex) =  &
                   weights(memberIndex1,memberIndex2,lonIndex,latIndex) + &
                   wInterpInfo%interpWeights(interpIndex,lonIndex,latIndex) *  &
                   weights(memberIndex1,memberIndex2,interpLonIndex,interpLatIndex)
            end do
          end do

        end do ! interpIndex

      end do ! lonIndex
    end do ! latIndex
    !$OMP END PARALLEL DO

    if (firstCall) write(*,*) 'enkf_interpWeights: totalCount = ', totalCount(:)
    firstCall = .false.

  end subroutine enkf_interpWeights

  !--------------------------------------------------------------------------
  ! enkf_modifyAMSUBobsError
  !--------------------------------------------------------------------------
  subroutine enkf_modifyAMSUBobsError(obsSpaceData)
    implicit none

    ! Arguments:
    type(struct_obs), target, intent(inout) :: obsSpaceData

    ! Locals:
    real(pre_obsReal), parameter :: AMSUB_trop_oer = 1.0 ! assumed value for AMSU-B obs error in tropics
    integer            :: headerIndex, bodyIndex, bodyIndexBeg, bodyIndexEnd, codeType
    real(pre_obsReal)  :: lat_obs

    ! for AMSUB observations set the observation error std dev equal to 1.0
    ! in the larger tropical area where the spread-skill correlation suggests 
    ! that the data are accurate (.i.e |lat|<40. ). Otherwise don't reduce the
    ! observational error.
    do headerIndex = 1, obs_numheader(obsSpaceData)
      lat_obs = obs_headElem_r(obsSpaceData, obs_lat, headerIndex)
      codeType = obs_headElem_i(obsSpaceData, obs_ity, headerIndex)
      lat_obs = lat_obs * MPC_DEGREES_PER_RADIAN_R8
      if ( abs(lat_obs) < 40. .and. (codeType == codtyp_get_codtyp('amsub') .or.  &
                                     codeType == codtyp_get_codtyp('mhs') .or.  &
                                     codeType == codtyp_get_codtyp('mwhs2')) ) then
        bodyIndexBeg = obs_headElem_i(obsSpaceData, obs_rln, headerIndex)
        bodyIndexEnd = obs_headElem_i(obsSpaceData, obs_nlv, headerIndex) + bodyIndexBeg - 1
        do bodyIndex = bodyIndexBeg, bodyIndexEnd
          call obs_bodySet_r(obsSpaceData, obs_oer, bodyIndex, AMSUB_trop_oer)
        end do
      end if
    end do

  end subroutine enkf_modifyAMSUBobsError

  !--------------------------------------------------------------------------
  ! enkf_rejectHighLatIR
  !--------------------------------------------------------------------------
  subroutine enkf_rejectHighLatIR(obsSpaceData)
    implicit none

    ! Arguments:
    type(struct_obs), target, intent(inout) :: obsSpaceData

    ! Locals:
    integer           :: headerIndex, bodyIndex, bodyIndexBeg, bodyIndexEnd, codeType
    real(pre_obsReal) :: lat_obs

    ! reject all HIR radiance observation in arctic and antarctic (.i.e |lat|>60. )
    do headerIndex = 1, obs_numheader(obsSpaceData)
      lat_obs = obs_headElem_r(obsSpaceData, obs_lat, headerIndex)
      codeType = obs_headElem_i(obsSpaceData, obs_ity, headerIndex)
      lat_obs = lat_obs * MPC_DEGREES_PER_RADIAN_R8
      if ( abs(lat_obs) > 60. .and. tvs_isIdBurpHyperSpectral(codeType) ) then
        write(*,*) 'enkf_rejectHighLatIR: !!!!!!!!--------WARNING--------!!!!!!!!'
        write(*,*) 'enkf_rejectHighLatIR: This HIR radiance profile was rejected because |lat|>60.'
        write(*,*) 'enkf_rejectHighLatIR: latidude= ', lat_obs, 'codtyp= ', codeType
        bodyIndexBeg = obs_headElem_i(obsSpaceData, obs_rln, headerIndex)
        bodyIndexEnd = obs_headElem_i(obsSpaceData, obs_nlv, headerIndex) + bodyIndexBeg - 1
        do bodyIndex = bodyIndexBeg, bodyIndexEnd
          call obs_bodySet_i(obsSpaceData, obs_ass, bodyIndex, obs_notAssimilated)
          ! also set the 'rejected by selection process' flag (bit 11)
          call obs_bodySet_i( obsSpaceData, obs_flg, bodyIndex,  &
                              ibset( obs_bodyElem_i( obsSpaceData, obs_flg, bodyIndex ), 11) )
        end do
      end if
    end do

  end subroutine enkf_rejectHighLatIR

  !--------------------------------------------------------------------------
  ! enkf_getModulatedState
  !--------------------------------------------------------------------------
  subroutine enkf_getModulatedState( stateVector_in, stateVectorMeanTrl, &
                                     vLocalizeLengthScale, numRetainedEigen, nEns, &
                                     eigenVectorColumnIndex, stateVector_out, &
                                     beSilent )
    !
    !:Purpose: Compute vertical localization matrix, and the corresponding
    !          eigenvectors/eigenvalues, to obtain modulated stateVector.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)    :: stateVector_in
    type(struct_gsv), intent(in)    :: stateVectorMeanTrl
    real(8),          intent(in)    :: vLocalizeLengthScale
    integer,          intent(in)    :: numRetainedEigen
    integer,          intent(in)    :: nEns
    integer,          intent(in)    :: eigenVectorColumnIndex
    type(struct_gsv), intent(inout) :: stateVector_out
    logical,          intent(in)    :: beSilent

    ! Locals:
    real(4)          :: modulationFactor_r4
    real(4), pointer :: field_out_r4(:,:,:,:)
    integer :: nLev, nlev_out, levIndex, latIndex, lonIndex
    integer :: lon1, lon2, lat1, lat2
    integer :: varIndex, stepIndex, eigenVectorLevelIndex
    character(len=4) :: varName

    call utl_tmg_start(130,'--getModulatedState')

    if ( .not. beSilent ) write(*,*) 'enkf_getModulatedState: START'

    if ( stateVector_in%dataKind /= 4 ) then
      call utl_abort('enkf_getModulatedState: only dataKind=4 is implemented')
    end if

    nLev = stateVector_in%vco%nLev_M
    if ( vLocalizeLengthScale <= 0.0d0 .or. nLev <= 1 ) then
      call utl_abort('enkf_getModulatedState: no vertical localization')
    end if

    ! Compute perturbation by subtracting ensMean
    call gsv_copy(stateVector_in, stateVector_out, beSilent_opt=beSilent)
    call gsv_add(stateVectorMeanTrl, stateVector_out, scaleFactor_opt=-1.0d0)

    lon1 = stateVector_out%myLonBeg
    lon2 = stateVector_out%myLonEnd
    lat1 = stateVector_out%myLatBeg
    lat2 = stateVector_out%myLatEnd

    ! Compute modulated member perturbation from original member perturbation:
    !   v'_k = (Nens*nLamda/(Nens - 1))^1/2 * Lambda^1/2 * E * x'_k
    step_loop: do stepIndex = 1, stateVector_out%numStep
      var_loop: do varIndex = 1, vnl_numvarmax
        varName = vnl_varNameList(varIndex)
        if ( .not. gsv_varExist(stateVector_out,varName) ) cycle var_loop

        nlev_out  = stateVector_out%varNumLev(varIndex)

        call gsv_getField(statevector_out,field_out_r4,varName)

        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2
            do levIndex = 1, nlev_out
              if ( nlev_out == 1 ) then
                eigenVectorLevelIndex = nLev
              else
                eigenVectorLevelIndex = levIndex
              end if

              call getModulationFactor( stateVector_in%vco, eigenVectorLevelIndex, &
                                        eigenVectorColumnIndex, numRetainedEigen, &
                                        nEns, vLocalizeLengthScale, &
                                        modulationFactor_r4, beSilent_opt=beSilent )

              field_out_r4(lonIndex,latIndex,levIndex,stepIndex) = &
                                 field_out_r4(lonIndex,latIndex,levIndex,stepIndex) * &
                                 modulationFactor_r4
            end do
          end do
        end do

      end do var_loop
    end do step_loop

    ! Now add to ensMean to get modulated member
    ! v_k = v'_k + v_mean
    call gsv_add(stateVectorMeanTrl, stateVector_out)

    if ( .not. beSilent ) write(*,*) 'enkf_getModulatedState: END'

    call utl_tmg_stop(130)

  end subroutine enkf_getModulatedState

  !--------------------------------------------------------------------------
  ! enkf_setupModulationFactor
  !--------------------------------------------------------------------------
  subroutine enkf_setupModulationFactor(vco, numRetainedEigen, nEns, vLocalizeLengthScale, &
                                        beSilent)
    !
    !:Purpose: setup modulationFactorArray by calling getModulationFactor for first time. 
    !
    implicit none

    ! Arguments:
    type(struct_vco), pointer, intent(in) :: vco
    integer,                   intent(in) :: numRetainedEigen
    integer,                   intent(in) :: nEns
    real(8),                   intent(in) :: vLocalizeLengthScale
    logical,                   intent(in) :: beSilent

    ! Locals:
    integer :: eigenVectorColumnIndex
    integer :: eigenVectorLevelIndex
    real(4) :: modulationFactor_r4

    eigenVectorColumnIndex = 1
    eigenVectorLevelIndex = 1
    call getModulationFactor(vco, eigenVectorLevelIndex, &
                             eigenVectorColumnIndex, numRetainedEigen, &
                             nEns, vLocalizeLengthScale, &
                             modulationFactor_r4, beSilent_opt=beSilent)
     
  end subroutine enkf_setupModulationFactor

  !--------------------------------------------------------------------------
  ! getModulationFactor
  !--------------------------------------------------------------------------
  subroutine getModulationFactor( vco, eigenVectorLevelIndex, &
                                  eigenVectorColumnIndex, numRetainedEigen, &
                                  nEns, vLocalizeLengthScale, &
                                  modulationFactor_r4, beSilent_opt )
    !
    !:Purpose: compute modulation factor needed to multiply ensemble
    !          perturbation to get the modulated perturbation:
    !          (Nens*nLambda/(Nens - 1))^1/2 * Lambda^1/2
    !
    implicit none

    ! Arguments:
    type(struct_vco), pointer, intent(in)  :: vco
    integer,                   intent(in)  :: eigenVectorLevelIndex
    integer,                   intent(in)  :: eigenVectorColumnIndex
    integer,                   intent(in)  :: numRetainedEigen
    integer,                   intent(in)  :: nEns
    real(8),                   intent(in)  :: vLocalizeLengthScale
    real(4),                   intent(out) :: modulationFactor_r4
    logical, optional,         intent(in)  :: beSilent_opt

    ! Locals:
    integer             :: levIndex1, levIndex2, eigenIndex
    integer             :: nLev, nLev_M, nLev_depth, matrixRank
    real(8)             :: zr, zcorr, pSurfRef
    real(8)             :: tolerance
    real(8), pointer    :: pressureProfile(:)
    real(8), allocatable, save :: eigenValues(:)
    real(8), allocatable, save :: eigenVectors(:,:)
    real(8), allocatable, save :: verticalLocalizationMat(:,:)
    real(8), allocatable, save :: verticalLocalizationMatLowRank(:,:)
    real(4), allocatable, save :: modulationFactorArray_r4(:,:)
    logical :: beSilent

    logical, save :: firstCall = .true.

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .false.
    end if

    ! Compute vertical localization matrix and its eigenValues/Vectors on first call
    if ( firstCall ) then
      firstCall = .false.
      if ( mmpi_myid == 0 .and. .not. beSilent ) then
        write(*,*) 'getModulationFactor: computing eigenValues/Vectors'
      end if

      nLev_M = vco%nLev_M
      nLev_depth = vco%nlev_depth
      nLev = max(nLev_M,nLev_depth)

      allocate(eigenValues(nLev))
      allocate(eigenVectors(nLev,nLev))
      allocate(verticalLocalizationMat(nLev,nLev))
      allocate(verticalLocalizationMatLowRank(nLev,nLev))
      allocate(modulationFactorArray_r4(numRetainedEigen,nLev))
      verticalLocalizationMatLowRank(:,:) = 0.0d0

      pSurfRef = 101000.D0
      call czp_fetch1DLevels(vco, pSurfRef, profM_opt=pressureProfile)

      call lfn_Setup(LocFunctionWanted='FifthOrder')

      ! Calculate 5'th order function
      do levIndex1 = 1, nLev
        do levIndex2 = 1, nLev
          zr = abs(log(pressureProfile(levIndex2)) - log(pressureProfile(levIndex1)))
          zcorr = lfn_response(zr,vLocalizeLengthScale)
          verticalLocalizationMat(levIndex1,levIndex2) = zcorr
        end do
      end do

      ! Compute eigenValues/Vectors of vertical localization matrix
      tolerance = 1.0D-50
      call utl_eigenDecomp(verticalLocalizationMat, eigenValues, eigenVectors, &
                           tolerance, matrixRank)
      if ( matrixRank < numRetainedEigen ) then
        write(*,*) 'matrixRank=', matrixRank
        call utl_abort('getModulationFactor: verticalLocalizationMat is rank deficient=')
      end if

      ! Compute low-ranked vertical localization matrix
      do levIndex1 = 1, nLev
        do levIndex2 = 1, nLev
          do eigenIndex = 1, numRetainedEigen
            verticalLocalizationMatLowRank(levIndex1,levIndex2) = verticalLocalizationMatLowRank(levIndex1,levIndex2) + & 
                                                                  eigenVectors(levIndex1,eigenIndex) * &
                                                                  eigenVectors(levIndex2,eigenIndex) * &
                                                                  eigenValues(eigenIndex)
          end do
        end do
      end do

      ! now compute the 2D modulationFactor array
      do levIndex1 = 1, nLev
        do eigenIndex = 1, numRetainedEigen
          modulationFactorArray_r4(eigenIndex,levIndex1) = real( &
                        1 / sqrt(verticalLocalizationMatLowRank(levIndex1,levIndex1)) * &
                        eigenVectors(levIndex1,eigenIndex) * &
                        eigenValues(eigenIndex) ** 0.5 * &
                        (nEns * numRetainedEigen / (nEns - 1)) ** 0.5,4)
        end do
      end do

      if ( mmpi_myid == 0 .and. .not. beSilent ) then
        do levIndex1 = 1, numRetainedEigen
          write(*,*) 'getModulationFactor: eigen mode=', levIndex1, ', eigenVectors=', eigenVectors(:,levIndex1)
        end do
        write(*,*) 'getModulationFactor: eigenValues=', eigenValues(1:numRetainedEigen)

        do levIndex1 = 1, nLev
          write(*,*) 'getModulationFactor: verticalLocalizationMat for lev ', levIndex1, '=', verticalLocalizationMat(levIndex1,:)
          write(*,*) 'getModulationFactor: verticalLocalizationMatLowRank for lev ', levIndex1, '=', verticalLocalizationMatLowRank(levIndex1,:)
        end do
      end if
    end if

    modulationFactor_r4 = modulationFactorArray_r4(eigenVectorColumnIndex,eigenVectorLevelIndex)
  
  end subroutine getModulationFactor

  !--------------------------------------------------------------------------
  ! enkf_getLocalizationRadius
  !--------------------------------------------------------------------------
  subroutine enkf_getLocalizationRadius(hLocalize, hLocalizePressure, &
                                        anlVertLocation, hLinearLoc, hLoc)
    !
    !:Purpose: get the localization radius, interpolated or not, at a given pressure
    !
    implicit none

    ! Arguments:
    real(8), intent(in)  :: hLocalize(:)         ! the list of localization radii (m)
    real(8), intent(in)  :: hLocalizePressure(:) ! the pressures where the radius changes (log(P))
    real(8), intent(in)  :: anlVertLocation      ! the gridpoint vertical coordinate in log(P)
    logical, intent(in)  :: hLinearLoc           ! apply linear vertical interpolation for the localization radius
    real(8), intent(out) :: hLoc                 ! the gridpoint localization radius

    ! Locals:
    integer              :: hLocIndex
    logical              :: hLocNotAssigned

    ! radius is constant
    if ( all(hLocalize(:) == hLocalize(1)) ) then
      hLoc = hLocalize(1)

    ! radius varies vertically, and is linearly interpolated with log(P)
    else if (hLinearLoc) then
      hLocNotAssigned = .true.
      do hLocIndex = 1, count(hLocalize > 0.0d0) - 1
        ! piece-wise linear interpolation
        if ( anlVertLocation >= hLocalizePressure(hLocIndex)   .and. &
             anlVertLocation <= hLocalizePressure(hLocIndex+1) .and. &
             hLocNotAssigned ) then
          hLoc = hLocalize(hLocIndex) + &
                 (anlVertLocation - hLocalizePressure(hLocIndex)) * &
                 (hLocalize(hLocIndex+1) - hLocalize(hLocIndex))  / &
                 (hLocalizePressure(hLocIndex+1) - hLocalizePressure(hLocIndex))
          hLocNotAssigned = .false.
        ! constant radius value near the top of the atmosphere
        else if ( anlVertLocation <= hLocalizePressure(hLocIndex)   .and. &
                  hLocNotAssigned ) then
          hLoc = hLocalize(1)
        ! constant radius value near the bottom of the atmosphere
        else if ( anlVertLocation >= hLocalizePressure(hLocIndex+1) .and. &
                  hLocNotAssigned ) then
          hLoc = hLocalize(count(hLocalize > 0.0d0))
        end if
      end do

    ! radius varies vertically, but is not interpolated
    else
      hLocIndex = 1 + count(anlVertLocation > hLocalizePressure(:)) - count(hLocalizePressure <= 0.0d0)
      hLoc = hLocalize(hLocIndex)

    end if

  end subroutine enkf_getLocalizationRadius

end module enkf_mod

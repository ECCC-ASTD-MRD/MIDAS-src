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

module enkf_mod
  ! MODULE enkf_mod (prefix='enkf' category='1. High-level functionality')
  !
  ! :Purpose: Various routines that are useful for implementing
  !           an EnKF in MIDAS, including the LETKF.
  !
  use mpi, only : mpi_statuses_ignore ! this is the mpi library module
  use midasMpi_mod
  use utilities_mod
  use mathPhysConstants_mod
  use timeCoord_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use ensembleStateVector_mod
  use gridStateVector_mod
  use obsSpaceData_mod
  use tovs_nl_mod
  use ensembleObservations_mod
  use gridVariableTransforms_mod
  use localizationFunction_mod
  use varNameList_mod
  use codePrecision_mod
  use codTyp_mod
  implicit none
  save
  private

  ! public types
  public :: struct_enkfInterpInfo

  ! public procedures
  public :: enkf_setupInterpInfo, enkf_LETKFanalyses, enkf_modifyAMSUBobsError
  public :: enkf_rejectHighLatIR, enkf_getModulatedState

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
  end type struct_enkfInterpInfo

  integer, external :: get_max_rss
  logical :: debug

contains

  !----------------------------------------------------------------------
  ! enkf_LETKFanalyses
  !----------------------------------------------------------------------
  subroutine enkf_LETKFanalyses(algorithm, numSubEns, randomShuffleSubEns,  &
                                ensembleAnl, ensembleTrl, &
                                ensObs_mpiglobal, ensObsGain_mpiglobal, &
                                stateVectorMeanAnl, &
                                wInterpInfo, maxNumLocalObs,  &
                                hLocalize, hLocalizePressure, vLocalize,  &
                                mpiDistribution, numRetainedEigen)
    ! :Purpose: Local subroutine containing the code for computing
    !           the LETKF analyses for all ensemble members, ensemble
    !           mean.
    implicit none

    ! Arguments
    character(len=*)            :: algorithm
    integer                     :: numSubEns
    logical                     :: randomShuffleSubEns
    type(struct_ens), pointer   :: ensembleTrl
    type(struct_ens)            :: ensembleAnl
    type(struct_eob), target    :: ensObs_mpiglobal
    type(struct_eob)            :: ensObsGain_mpiglobal
    type(struct_gsv)            :: stateVectorMeanAnl
    type(struct_enkfInterpInfo) :: wInterpInfo
    integer                     :: maxNumLocalObs
    real(8)                     :: hLocalize(:)
    real(8)                     :: hLocalizePressure(:)
    real(8)                     :: vLocalize
    character(len=*)            :: mpiDistribution
    integer                     :: numRetainedEigen

    ! Locals
    integer :: nEns, nEnsPerSubEns, nEnsIndependentPerSubEns
    integer :: nEnsPerSubEns_mod, nEnsIndependentPerSubEns_mod
    integer :: nLev_M, nLev_depth, nLev_weights
    integer :: memberIndex, memberIndex1, memberIndex2, ierr, matrixRank
    integer :: memberIndexCV, memberIndexCV1, memberIndexCV2
    integer :: procIndex, procIndexSend, hLocIndex
    integer :: latIndex, lonIndex, stepIndex, varLevIndex, levIndex, levIndex2
    integer :: bodyIndex, localObsIndex, numLocalObs, numLocalObsFound
    integer :: countMaxExceeded, maxCountMaxExceeded, numGridPointWeights
    integer :: myNumLatLonRecv, myNumLatLonSend
    integer :: latLonIndex, subEnsIndex, subEnsIndex2
    integer :: sendTag, recvTag, nsize, numRecv, numSend
    integer :: myLonBeg, myLonEnd, myLatBeg, myLatEnd, numVarLev
    integer :: myLonBegHalo, myLonEndHalo, myLatBegHalo, myLatEndHalo
    integer :: imode, dateStamp, timePrint, datePrint, randomSeed, newDate
    integer :: nEnsUsed, eigenVectorColumnIndex, eigenVectorLevelIndex
    integer :: memberIndexInModEns, retainedEigenIndex, nLev
    real(8) :: anlLat, anlLon, anlVertLocation, distance, tolerance, localization
    real(8) :: modulationFactor

    integer, allocatable :: localBodyIndices(:)
    integer, allocatable :: myLatIndexesRecv(:), myLonIndexesRecv(:)
    integer, allocatable :: myLatIndexesSend(:), myLonIndexesSend(:)
    integer, allocatable :: myNumProcIndexesSend(:)
    integer, allocatable :: myProcIndexesRecv(:), myProcIndexesSend(:,:)
    integer, allocatable :: requestIdRecv(:), requestIdSend(:)
    integer, allocatable :: memberIndexSubEns(:,:), memberIndexSubEnsComp(:,:)
    integer, allocatable :: randomMemberIndexArray(:), latLonTagMpiGlobal(:,:)
    integer, allocatable :: memberIndexSubEns_mod(:,:), memberIndexSubEnsComp_mod(:,:)

    real(8), allocatable :: distances(:)
    real(8), allocatable :: PaInv(:,:), PaSqrt(:,:), Pa(:,:), YbTinvR(:,:), YbTinvRYb(:,:)
    real(8), allocatable :: YbTinvRYb_CV(:,:), YbTinvRYb_CV_mod(:,:), YbTinvRYb_mod(:,:)
    real(8), allocatable :: eigenValues(:), eigenVectors(:,:)
    real(8), allocatable :: eigenValues_CV(:), eigenVectors_CV(:,:)
    real(8), allocatable :: eigenValues_CV_mod(:), eigenVectors_CV_mod(:,:)
    real(8), allocatable :: weightsTemp(:), weightsTemp2(:)
    real(8), allocatable :: weightsMembers(:,:,:,:), weightsMembersLatLon(:,:,:)
    real(8), allocatable :: weightsMean(:,:,:,:), weightsMeanLatLon(:,:,:)
    real(8), allocatable :: memberAnlPert(:)
    real(4), allocatable :: vertLocation_r4(:,:,:)

    real(4), pointer     :: meanTrl_ptr_r4(:,:,:,:), meanAnl_ptr_r4(:,:,:,:), meanInc_ptr_r4(:,:,:,:)
    real(4), pointer     :: memberTrl_ptr_r4(:,:,:,:), memberAnl_ptr_r4(:,:,:,:)
    real(4)              :: pert_r4

    character(len=4)     :: varLevel
    character(len=2)     :: varKind

    type(struct_hco), pointer :: hco_ens
    type(struct_vco), pointer :: vco_ens
    type(struct_gsv)          :: stateVectorMeanInc
    type(struct_gsv)          :: stateVectorMeanTrl

    logical :: hLocalizeIsConstant, firstTime = .true.
    logical :: printBeforeComm = .true.
    logical :: printAfterComm = .true.
    integer :: latIndexPrint, lonIndexPrint

    call utl_tmg_start(100,'--LETKFanalysis')

    write(*,*) 'enkf_LETKFanalyses: starting'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    ! Set things up for the redistribution of work across mpi tasks
    !
    call enkf_LETKFsetupMpiDistribution(myNumLatLonRecv, myNumLatLonSend,   &
                                        myLatIndexesRecv, myLonIndexesRecv,   &
                                        myLatIndexesSend, myLonIndexesSend,   &
                                        myProcIndexesRecv, myProcIndexesSend, &
                                        myNumProcIndexesSend, mpiDistribution, wInterpInfo)
    allocate(requestIdSend(3*myNumLatLonSend*maxval(myNumProcIndexesSend)))
    allocate(requestIdRecv(3*myNumLatLonRecv))

    nEns       = ens_getNumMembers(ensembleAnl)
    if ( numRetainedEigen > 0 ) then
      nEnsUsed   = nEns * numRetainedEigen
    else
      nEnsUsed   = nEns
    end if
    nLev_M     = ens_getNumLev(ensembleAnl, 'MM')
    nLev_depth = ens_getNumLev(ensembleAnl, 'DP')
    nLev_weights = max(nLev_M,nLev_depth)
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

    !
    ! Compute gridded 3D ensemble weights
    !
    allocate(localBodyIndices(maxNumLocalObs))
    allocate(distances(maxNumLocalObs))
    allocate(YbTinvR(nEnsUsed,maxNumLocalObs))
    allocate(YbTinvRYb(nEnsUsed,nEnsUsed))
    if ( trim(algorithm) == 'CVLETKF-ME' ) allocate(YbTinvRYb_mod(nEnsUsed,nEns))
    allocate(eigenValues(nEnsUsed))
    allocate(eigenVectors(nEnsUsed,nEnsUsed))
    allocate(PaInv(nEnsUsed,nEnsUsed))
    allocate(PaSqrt(nEnsUsed,nEnsUsed))
    allocate(Pa(nEnsUsed,nEnsUsed))
    allocate(memberAnlPert(nEns))
    allocate(weightsTemp(nEnsUsed))
    allocate(weightsTemp2(nEnsUsed))
    weightsTemp(:) = 0.0d0
    weightsTemp2(:) = 0.0d0
    ! Weights for mean analysis
    allocate(weightsMean(nEnsUsed,1,myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
    weightsMean(:,:,:,:) = 0.0d0
    allocate(weightsMeanLatLon(nEnsUsed,1,myNumLatLonSend))
    weightsMeanLatLon(:,:,:) = 0.0d0
    ! Weights for member analyses
    allocate(weightsMembers(nEnsUsed,nEns,myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
    weightsMembers(:,:,:,:) = 0.0d0
    allocate(weightsMembersLatLon(nEnsUsed,nEns,myNumLatLonSend))
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
      nEnsIndependentPerSubEns = nEns - nEnsPerSubEns
      allocate(YbTinvRYb_CV(nEnsIndependentPerSubEns,nEnsIndependentPerSubEns))
      allocate(eigenValues_CV(nEnsIndependentPerSubEns))
      allocate(eigenVectors_CV(nEnsIndependentPerSubEns,nEnsIndependentPerSubEns))
      allocate(memberIndexSubEns(nEnsPerSubEns,numSubEns))
      allocate(memberIndexSubEnsComp(nEnsIndependentPerSubEns,numSubEns))
      if ( numRetainedEigen > 0 ) then
        nEnsPerSubEns_mod = nEnsPerSubEns * numRetainedEigen
        nEnsIndependentPerSubEns_mod = nEnsUsed - nEnsPerSubEns_mod
        allocate(YbTinvRYb_CV_mod(nEnsIndependentPerSubEns_mod,nEnsIndependentPerSubEns_mod))
        allocate(eigenValues_CV_mod(nEnsIndependentPerSubEns_mod))
        allocate(eigenVectors_CV_mod(nEnsIndependentPerSubEns_mod,nEnsIndependentPerSubEns_mod))
        allocate(memberIndexSubEns_mod(nEnsPerSubEns_mod,numSubEns))
        allocate(memberIndexSubEnsComp_mod(nEnsIndependentPerSubEns_mod,numSubEns))
      end if
      if (.not.randomShuffleSubEns) then
        ! form subensembles with contiguous sequential groups of members
        do subEnsIndex = 1, numSubEns
          do memberIndex = 1, nEnsPerSubEns
            memberIndexSubEns(memberIndex,subEnsIndex) =  &
                 (subEnsIndex-1)*nEnsPerSubEns + memberIndex
          end do
        end do
        if ( numRetainedEigen > 0 ) then
          do subEnsIndex = 1, numSubEns
            memberIndex2 = 0
            do memberIndex = 1, nEnsPerSubEns
              do eigenVectorColumnIndex = 1, numRetainedEigen
                !memberIndex2 = memberIndex2 + memberIndex
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
        if ( numRetainedEigen > 0 ) then
          call utl_abort('enkf_LETKFanalyses: randomShuffleSubEns not implemented for algorithm:' // &
                          trim(algorithm))
        end if

        ! compute random seed from the date for randomly forming subensembles
        imode = -3 ! stamp to printable date and time: YYYYMMDD, HHMMSShh
        dateStamp = tim_getDateStamp()
        ierr = newdate(dateStamp, datePrint, timePrint, imode)
        timePrint = timePrint/1000000
        datePrint =  datePrint*100 + timePrint
        ! Remove the century, keeping 2 digits of the year
        randomSeed = datePrint - 100000000*(datePrint/100000000)
        allocate(randomMemberIndexArray(nEns))
        do memberIndex = 1, nEns
          randomMemberIndexArray(memberIndex) = memberIndex
        end do
        call utl_randomOrderInt(randomMemberIndexArray,randomSeed)
        write(*,*) 'enkf_LETKFanalyses: seed for random shuffle of sub ens = ', randomSeed
        write(*,*) 'enkf_LETKFanalyses: randomOrder = ', randomMemberIndexArray(:)
        do subEnsIndex = 1, numSubEns
          do memberIndex = 1, nEnsPerSubEns
            memberIndexSubEns(memberIndex,subEnsIndex) =  &
                 randomMemberIndexArray((subEnsIndex-1)*nEnsPerSubEns + memberIndex)
          end do
        end do
      end if

      do subEnsIndex = 1, numSubEns
        memberIndex = 1
        do subEnsIndex2 = 1, numSubEns
          if (subEnsIndex2 == subEnsIndex) cycle

          memberIndexSubEnsComp(memberIndex:memberIndex+nEnsPerSubEns-1,subEnsIndex) =  &
            memberIndexSubEns(:,subEnsIndex2)

          memberIndex = memberIndex + nEnsPerSubEns
        end do

        if ( numRetainedEigen > 0 ) then
          memberIndex2 = 1
          do subEnsIndex2 = 1, numSubEns
            if (subEnsIndex2 == subEnsIndex) cycle

              memberIndexSubEnsComp_mod(memberIndex2:memberIndex2+nEnsPerSubEns_mod-1,subEnsIndex) =  &
                memberIndexSubEns_mod(:,subEnsIndex2)

            memberIndex2 = memberIndex2 + nEnsPerSubEns_mod
          end do
        end if

      end do

      write(*,*) 'nEns, numSubEns, nEnsPerSubEns, nEnsIndependentPerSubEns = ',  &
                 nEns, numSubEns, nEnsPerSubEns, nEnsIndependentPerSubEns
      do subEnsIndex = 1, numSubEns
        write(*,*) 'memberIndexSubEns = '
        write(*,*) memberIndexSubEns(:,subEnsIndex)
        if ( numRetainedEigen > 0 ) then
          write(*,*) 'memberIndexSubEns_mod = '
          write(*,*) memberIndexSubEns_mod(:,subEnsIndex)
        end if

        write(*,*) 'memberIndexSubEnsComp = '
        write(*,*) memberIndexSubEnsComp(:,subEnsIndex)
        if ( numRetainedEigen > 0 ) then
          write(*,*) 'memberIndexSubEnsComp_mod = '
          write(*,*) memberIndexSubEnsComp_mod(:,subEnsIndex)
        end if

      end do

    end if ! if CVLETKF(-PERTOBS)(-ME) algorithm

    call lfn_Setup(LocFunctionWanted='FifthOrder')

    ! compute 3D field of vertical location needed for localization
    if (vLocalize > 0.0d0) then
      call enkf_computeVertLocation(vertLocation_r4,stateVectorMeanTrl)
    end if

    call utl_tmg_start(108,'----Barr')
    call rpn_comm_barrier('GRID',ierr)
    call utl_tmg_stop(108)

    ! get mpi global list of tags used for mpi send/recv
    call utl_tmg_start(109, '----GetGlobalTags')
    allocate(latLonTagMpiGlobal(stateVectorMeanAnl%ni,stateVectorMeanAnl%nj))
    call enkf_LETKFgetMpiGlobalTags(latLonTagMpiGlobal,myLatIndexesRecv,myLonIndexesRecv)
    call utl_tmg_stop(109)

    ! Compute the weights for ensemble mean and members
    countMaxExceeded = 0
    maxCountMaxExceeded = 0
    numGridPointWeights = 0
    LEV_LOOP: do levIndex = 1, nLev_weights
      write(*,*) 'computing ensemble updates for vertical level = ', levIndex

      !
      ! First post all recv instructions for communication of weights
      !
      call utl_tmg_start(101,'----CommWeights')
      numSend = 0
      numRecv = 0
      do latLonIndex = 1, myNumLatLonRecv
        latIndex = myLatIndexesRecv(latLonIndex)
        lonIndex = myLonIndexesRecv(latLonIndex)
        procIndex = myProcIndexesRecv(latLonIndex)
        recvTag = latLonTagMpiGlobal(lonIndex,latIndex)

        nsize = nEnsUsed
        numRecv = numRecv + 1
        call mpi_irecv( weightsMean(:,1,lonIndex,latIndex),  &
                        nsize, mmpi_datyp_real8, procIndex-1, recvTag,  &
                        mmpi_comm_grid, requestIdRecv(numRecv), ierr )
        nsize = nEnsUsed * nEns
        numRecv = numRecv + 1
        recvTag = recvTag + maxval(latLonTagMpiGlobal(:,:))
        call mpi_irecv( weightsMembers(:,:,lonIndex,latIndex),  &
                        nsize, mmpi_datyp_real8, procIndex-1, recvTag,  &
                        mmpi_comm_grid, requestIdRecv(numRecv), ierr )
      end do
      call utl_tmg_stop(101)

      LATLON_LOOP: do latLonIndex = 1, myNumLatLonSend
        latIndex = myLatIndexesSend(latLonIndex)
        lonIndex = myLonIndexesSend(latLonIndex)

        numGridPointWeights = numGridPointWeights + 1

        ! lat-lon of the grid point for which we are doing the analysis
        anlLat = hco_ens%lat2d_4(lonIndex,latIndex)
        anlLon = hco_ens%lon2d_4(lonIndex,latIndex)
        hLocalizeIsConstant = all(hLocalize(:) == hLocalize(1))
        if (vLocalize > 0.0d0 .or. .not.hLocalizeIsConstant) then
          anlVertLocation = vertLocation_r4(lonIndex,latIndex,levIndex)
        end if

        ! Find which horizontal localization value to use for this analysis level
        if (hLocalizeIsConstant) then
          hLocIndex = 1
        else
          hLocIndex = 1 + count(anlVertLocation > hLocalizePressure(:))
        end if

        ! Get list of nearby observations and distances to gridpoint
        call utl_tmg_start(102,'----GetLocalBodyIndices')
        numLocalObs = eob_getLocalBodyIndices(ensObs_mpiglobal, localBodyIndices,     &
                                              distances, anlLat, anlLon, anlVertLocation,  &
                                              hLocalize(hLocIndex), vLocalize, numLocalObsFound)
        if (numLocalObsFound > maxNumLocalObs) then
          countMaxExceeded = countMaxExceeded + 1
          maxCountMaxExceeded = max(maxCountMaxExceeded, numLocalObsFound)
        end if
        call utl_tmg_stop(102)

        call utl_tmg_start(103,'----CalculateWeights')

        ! Extract initial quantities YbTinvR and first term of PaInv (YbTinvR*Yb)
        do localObsIndex = 1, numLocalObs
          bodyIndex = localBodyIndices(localObsIndex)

          ! Compute value of localization function
          ! Horizontal
          localization = lfn_Response(distances(localObsIndex),hLocalize(hLocIndex))
          ! Vertical when NOT using modulated ensembles - use pressures at the grid point (not obs) location
          if ( vLocalize > 0.0d0 .and. numRetainedEigen == 0 ) then
            distance = abs( anlVertLocation - ensObs_mpiglobal%vertLocation(bodyIndex) )
            localization = localization * lfn_Response(distance,vLocalize)
          end if

          do memberIndex = 1, nEnsUsed
            YbTinvR(memberIndex,localObsIndex) =  &
                 ensObsGain_mpiglobal%Yb_r4(memberIndex, bodyIndex) * &
                 localization * ensObsGain_mpiglobal%obsErrInv(bodyIndex)
          end do
        end do ! localObsIndex

        call utl_tmg_start(105,'------CalcYbTinvRYb')
        YbTinvRYb(:,:) = 0.0D0
        do localObsIndex = 1, numLocalObs
          bodyIndex = localBodyIndices(localObsIndex)
          !$OMP PARALLEL DO PRIVATE (memberIndex1, memberIndex2)
          do memberIndex2 = 1, nEnsUsed
            do memberIndex1 = 1, nEnsUsed
              YbTinvRYb(memberIndex1,memberIndex2) =  &
                   YbTinvRYb(memberIndex1,memberIndex2) +  &
                   YbTinvR(memberIndex1,localObsIndex) * ensObsGain_mpiglobal%Yb_r4(memberIndex2, bodyIndex)
            end do
          end do
          !$OMP END PARALLEL DO

          ! computing YbTinvRYb that uses modulated and original ensembles for perturbation update
          if ( trim(algorithm) == 'CVLETKF-ME' ) then
            if (localObsIndex == 1) YbTinvRYb_mod(:,:) = 0.0D0
            !$OMP PARALLEL DO PRIVATE (memberIndex1, memberIndex2)
            do memberIndex2 = 1, nEns
              do memberIndex1 = 1, nEnsUsed
                YbTinvRYb_mod(memberIndex1,memberIndex2) =  &
                     YbTinvRYb_mod(memberIndex1,memberIndex2) +  &
                     YbTinvR(memberIndex1,localObsIndex) * ensObsGain_mpiglobal%Yb_r4(memberIndex2, bodyIndex)
              end do
            end do
            !$OMP END PARALLEL DO
          end if

        end do ! localObsIndex
        call utl_tmg_stop(105)

        ! Rest of the computation of local weights for this grid point
        if (numLocalObs > 0) then

          if (trim(algorithm) == 'LETKF') then
            !
            ! Weight calculation for standard LETKF algorithm
            !

            ! Add second term of PaInv
            PaInv(:,:) = YbTinvRYb(:,:)
            do memberIndex = 1, nEns
              PaInv(memberIndex,memberIndex) = PaInv(memberIndex,memberIndex) + real(nEns - 1,8)
            end do

            ! Compute Pa and sqrt(Pa) matrices from PaInv
            Pa(:,:) = PaInv(:,:)
            call utl_tmg_start(104,'------EigenDecomp')
            call utl_matInverse(Pa, nEns, inverseSqrt_opt=PaSqrt)
            call utl_tmg_stop(104)

            ! Compute ensemble mean local weights as Pa * YbTinvR * (obs - meanYb)
            weightsTemp(:) = 0.0d0
            do localObsIndex = 1, numLocalObs
              bodyIndex = localBodyIndices(localObsIndex)
              do memberIndex = 1, nEns
                weightsTemp(memberIndex) = weightsTemp(memberIndex) +   &
                                           YbTinvR(memberIndex,localObsIndex) *  &
                                           ( ensObs_mpiglobal%obsValue(bodyIndex) - &
                                             ensObs_mpiglobal%meanYb(bodyIndex) )
              end do
            end do

            weightsMeanLatLon(:,1,latLonIndex) = 0.0d0
            do memberIndex2 = 1, nEns
              do memberIndex1 = 1, nEns
                weightsMeanLatLon(memberIndex1,1,latLonIndex) =  &
                     weightsMeanLatLon(memberIndex1,1,latLonIndex) +  &
                     Pa(memberIndex1,memberIndex2)*weightsTemp(memberIndex2)
              end do
            end do

            ! Compute ensemble perturbation weights: [(Nens-1)^1/2*PaSqrt]
            weightsMembersLatLon(:,:,latLonIndex) = sqrt(real(nEns - 1,8)) * PaSqrt(:,:)

          else if (trim(algorithm) == 'LETKF-ME') then
            !
            ! Weight calculation for standard LETKF algorithm with modulated ensemble
            !

            ! Compute eigenValues/Vectors of Yb^T R^-1 Yb = E * Lambda * E^T
            call tmg_start(90,'LETKF-eigenDecomp')
            tolerance = 1.0D-50
            call utl_eigenDecomp(YbTinvRYb, eigenValues, eigenVectors, tolerance, matrixRank)
            call tmg_stop(90)

            ! Compute ensemble mean local weights as E * (Lambda + (Nens-1)*I)^-1 * E^T * YbTinvR * (obs - meanYb)
            weightsTemp(:) = 0.0d0
            do localObsIndex = 1, numLocalObs
              bodyIndex = localBodyIndices(localObsIndex)
              do memberIndex = 1, nEnsUsed
                weightsTemp(memberIndex) = weightsTemp(memberIndex) +   &
                                           YbTinvR(memberIndex,localObsIndex) *  &
                                           ( ensObs_mpiglobal%obsValue(bodyIndex) - &
                                             ensObs_mpiglobal%meanYb(bodyIndex) )
              end do
            end do
            weightsTemp2(:) = 0.0d0
            do memberIndex2 = 1, matrixRank
              do memberIndex1 = 1, nEnsUsed
                weightsTemp2(memberIndex2) = weightsTemp2(memberIndex2) +   &
                                             eigenVectors(memberIndex1,memberIndex2) *  &
                                             weightsTemp(memberIndex1)
              end do
            end do
            do memberIndex = 1, matrixRank
              weightsTemp2(memberIndex) = weightsTemp2(memberIndex) *  &
                                          1.0D0/(eigenValues(memberIndex) + real(nEnsUsed - 1,8))
            end do
            weightsMeanLatLon(:,1,latLonIndex) = 0.0d0
            do memberIndex2 = 1, matrixRank
              do memberIndex1 = 1, nEnsUsed
                weightsMeanLatLon(memberIndex1,1,latLonIndex) =  &
                     weightsMeanLatLon(memberIndex1,1,latLonIndex) +   &
                     eigenVectors(memberIndex1,memberIndex2) *  &
                     weightsTemp2(memberIndex2)
              end do
            end do

            ! Compute ensemble perturbation weights: [(Nens-1)^1/2*PaSqrt]
            !weightsMembersLatLon(:,:,latLonIndex) = sqrt(real(nEnsUsed - 1,8)) * PaSqrt(:,:)
          else if (trim(algorithm) == 'CVLETKF') then
            !
            ! Weight calculation for cross-validation LETKF algorithm
            !

            ! Compute eigenValues/Vectors of Yb^T R^-1 Yb = E * Lambda * E^T
            call utl_tmg_start(104,'------EigenDecomp')
            tolerance = 1.0D-50
            call utl_eigenDecomp(YbTinvRYb, eigenValues, eigenVectors, tolerance, matrixRank)
            call utl_tmg_stop(104)
            !if (matrixRank < (nEns-1)) then
            !  write(*,*) 'YbTinvRYb is rank deficient =', matrixRank, nEns, numLocalObs
            !end if

            ! Compute ensemble mean local weights as E * (Lambda + (Nens-1)*I)^-1 * E^T * YbTinvR * (obs - meanYb)
            weightsTemp(:) = 0.0d0
            do localObsIndex = 1, numLocalObs
              bodyIndex = localBodyIndices(localObsIndex)
              do memberIndex = 1, nEns
                weightsTemp(memberIndex) = weightsTemp(memberIndex) +   &
                                           YbTinvR(memberIndex,localObsIndex) *  &
                                           ( ensObs_mpiglobal%obsValue(bodyIndex) - &
                                             ensObs_mpiglobal%meanYb(bodyIndex) )
              end do
            end do
            weightsTemp2(:) = 0.0d0
            do memberIndex2 = 1, matrixRank
              do memberIndex1 = 1, nEns
                weightsTemp2(memberIndex2) = weightsTemp2(memberIndex2) +   &
                                             eigenVectors(memberIndex1,memberIndex2) *  &
                                             weightsTemp(memberIndex1)
              end do
            end do
            do memberIndex = 1, matrixRank
              weightsTemp2(memberIndex) = weightsTemp2(memberIndex) *  &
                                          1.0D0/(eigenValues(memberIndex) + real(nEns - 1,8))
            end do
            weightsMeanLatLon(:,1,latLonIndex) = 0.0d0
            do memberIndex2 = 1, matrixRank
              do memberIndex1 = 1, nEns
                weightsMeanLatLon(memberIndex1,1,latLonIndex) =  &
                     weightsMeanLatLon(memberIndex1,1,latLonIndex) +   &
                     eigenVectors(memberIndex1,memberIndex2) *  &
                     weightsTemp2(memberIndex2)
              end do
            end do

            ! Compute ensemble perturbation weights: 
            ! Wa = [ I - (Nens-1)^1/2 * E * 
            !        {(Nens-1)^-1/2*I - (Lambda + (Nens-1)*I)^-1/2} * Lambda^-1 *
            !        E^T * YbTinvRYb ]
            ! Loop over sub-ensembles
            do subEnsIndex = 1, numSubEns

              ! Use complement (independent) ens to get eigenValues/Vectors of Yb^T R^-1 Yb = E*Lambda*E^T
              call utl_tmg_start(104,'------EigenDecomp')
              do memberIndexCV2 = 1, nEnsIndependentPerSubEns
                memberIndex2 = memberIndexSubEnsComp(memberIndexCV2, subEnsIndex)
                do memberIndexCV1 = 1, nEnsIndependentPerSubEns
                  memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
                  YbTinvRYb_CV(memberIndexCV1,memberIndexCV2) = YbTinvRYb(memberIndex1,memberIndex2)
                end do
              end do
              tolerance = 1.0D-50
              call utl_eigenDecomp(YbTinvRYb_CV, eigenValues_CV, eigenVectors_CV, tolerance, matrixRank)
              call utl_tmg_stop(104)

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
                                                YbTinvRYb(memberIndex1,memberIndex)
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
                weightsMembersLatLon(:,memberIndex,latLonIndex) = 0.0d0
                do memberIndex2 = 1, matrixRank
                  do memberIndexCV1 = 1, nEnsIndependentPerSubEns
                    memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
                    weightsMembersLatLon(memberIndex1,memberIndex,latLonIndex) =   &
                         weightsMembersLatLon(memberIndex1,memberIndex,latLonIndex) +   &
                         eigenVectors_CV(memberIndexCV1,memberIndex2) *  &
                         weightsTemp(memberIndex2)
                  end do
                end do

                ! -1 * (Nens-1)^1/2 * previous_result
                weightsMembersLatLon(:,memberIndex,latLonIndex) =  &
                     -1.0D0 * sqrt(real(nEnsIndependentPerSubEns - 1,8)) *  &
                     weightsMembersLatLon(:,memberIndex,latLonIndex)

                ! I + previous_result
                weightsMembersLatLon(memberIndex,memberIndex,latLonIndex) =  &
                     1.0D0 + weightsMembersLatLon(memberIndex,memberIndex,latLonIndex)

              end do ! memberIndexCV
            end do ! subEnsIndex

            ! Remove the weights mean computed over the columns
            do memberIndex = 1, nEns
              weightsMembersLatLon(memberIndex,:,latLonIndex) =  &
                   weightsMembersLatLon(memberIndex,:,latLonIndex) - &
                   sum(weightsMembersLatLon(memberIndex,:,latLonIndex))/real(nEns,8)
            end do

          else if (trim(algorithm) == 'CVLETKF-ME') then
            !
            ! Weight calculation for cross-validation LETKF algorithm
            !

            ! Compute eigenValues/Vectors of Yb^T R^-1 Yb = E * Lambda * E^T
            call tmg_start(90,'LETKF-eigenDecomp')
            tolerance = 1.0D-50
            call utl_eigenDecomp(YbTinvRYb, eigenValues, eigenVectors, tolerance, matrixRank)
            call tmg_stop(90)
            !if (matrixRank < (nEns-1)) then
            !  write(*,*) 'YbTinvRYb is rank deficient =', matrixRank, nEns, numLocalObs
            !end if

            ! Compute ensemble mean local weights as E * (Lambda + (Nens-1)*I)^-1 * E^T * YbTinvR * (obs - meanYb)
            weightsTemp(:) = 0.0d0
            do localObsIndex = 1, numLocalObs
              bodyIndex = localBodyIndices(localObsIndex)
              do memberIndex = 1, nEnsUsed
                weightsTemp(memberIndex) = weightsTemp(memberIndex) +   &
                                           YbTinvR(memberIndex,localObsIndex) *  &
                                           ( ensObs_mpiglobal%obsValue(bodyIndex) - &
                                             ensObs_mpiglobal%meanYb(bodyIndex) )
              end do
            end do
            weightsTemp2(:) = 0.0d0
            do memberIndex2 = 1, matrixRank
              do memberIndex1 = 1, nEnsUsed
                weightsTemp2(memberIndex2) = weightsTemp2(memberIndex2) +   &
                                             eigenVectors(memberIndex1,memberIndex2) *  &
                                             weightsTemp(memberIndex1)
              end do
            end do
            do memberIndex = 1, matrixRank
              weightsTemp2(memberIndex) = weightsTemp2(memberIndex) *  &
                                          1.0D0/(eigenValues(memberIndex) + real(nEnsUsed - 1,8))
            end do
            weightsMeanLatLon(:,1,latLonIndex) = 0.0d0
            do memberIndex2 = 1, matrixRank
              do memberIndex1 = 1, nEnsUsed
                weightsMeanLatLon(memberIndex1,1,latLonIndex) =  &
                     weightsMeanLatLon(memberIndex1,1,latLonIndex) +   &
                     eigenVectors(memberIndex1,memberIndex2) *  &
                     weightsTemp2(memberIndex2)
              end do
            end do

            ! Compute ensemble perturbation weights: 
            ! Wa = [ - (Nens-1)^1/2 * E *
            !        {(Nens-1)^-1/2*I - (Lambda + (Nens-1)*I)^-1/2} * Lambda^-1 *
            !        E^T * YbTinvRYb_mod ]
            ! Loop over sub-ensembles
            do subEnsIndex = 1, numSubEns

              ! Use complement (independent) ens to get eigenValues/Vectors of Yb^T R^-1 Yb = E*Lambda*E^T
              call tmg_start(90,'LETKF-eigenDecomp')
              do memberIndexCV2 = 1, nEnsIndependentPerSubEns_mod
                memberIndex2 = memberIndexSubEnsComp_mod(memberIndexCV2, subEnsIndex)
                do memberIndexCV1 = 1, nEnsIndependentPerSubEns_mod
                  memberIndex1 = memberIndexSubEnsComp_mod(memberIndexCV1, subEnsIndex)
                  YbTinvRYb_CV_mod(memberIndexCV1,memberIndexCV2) = YbTinvRYb(memberIndex1,memberIndex2)
                end do
              end do
              tolerance = 1.0D-50
              call utl_eigenDecomp(YbTinvRYb_CV_mod, eigenValues_CV_mod, eigenVectors_CV_mod, tolerance, matrixRank)
              call tmg_stop(90)

              ! Loop over members within the current sub-ensemble being updated
              do memberIndexCV = 1, nEnsPerSubEns

                ! This is index of member being updated
                memberIndex = memberIndexSubEns(memberIndexCV, subEnsIndex)

                ! E^T * YbTinvRYb
                weightsTemp(:) = 0.0d0
                do memberIndex2 = 1, matrixRank
                  do memberIndexCV1 = 1, nEnsIndependentPerSubEns_mod
                    memberIndex1 = memberIndexSubEnsComp_mod(memberIndexCV1, subEnsIndex)
                    weightsTemp(memberIndex2) = weightsTemp(memberIndex2) +  &
                                                eigenVectors_CV_mod(memberIndexCV1,memberIndex2) *  &
                                                YbTinvRYb_mod(memberIndex1,memberIndex)
                  end do
                end do

                ! {(Nens-1)^-1/2*I - (Lambda + (Nens-1)*I)^-1/2} Lambda^-1 * previous_result

                do memberIndex1 = 1, matrixRank
                  weightsTemp(memberIndex1) = weightsTemp(memberIndex1) *  &
                                              ( 1.0D0/sqrt(real(nEnsIndependentPerSubEns_mod - 1,8)) -   &
                                                1.0D0/sqrt(eigenValues_CV_mod(memberIndex1) +  &
                                                           real(nEnsIndependentPerSubEns_mod - 1,8)) )
                  weightsTemp(memberIndex1) = weightsTemp(memberIndex1) /  &
                                              eigenValues_CV_mod(memberIndex1)
                end do

                ! E * previous_result
                weightsMembersLatLon(:,memberIndex,latLonIndex) = 0.0d0
                do memberIndex2 = 1, matrixRank
                  do memberIndexCV1 = 1, nEnsIndependentPerSubEns_mod
                    memberIndex1 = memberIndexSubEnsComp_mod(memberIndexCV1, subEnsIndex)
                    weightsMembersLatLon(memberIndex1,memberIndex,latLonIndex) =   &
                         weightsMembersLatLon(memberIndex1,memberIndex,latLonIndex) +   &
                         eigenVectors_CV_mod(memberIndexCV1,memberIndex2) *  &
                         weightsTemp(memberIndex2)
                  end do
                end do

                ! -1 * (Nens-1)^1/2 * previous_result
                weightsMembersLatLon(:,memberIndex,latLonIndex) =  &
                     -1.0D0 * sqrt(real(nEnsIndependentPerSubEns_mod - 1,8)) *  &
                     weightsMembersLatLon(:,memberIndex,latLonIndex)

              end do ! memberIndexCV
            end do ! subEnsIndex

            ! Remove the weights mean computed over the columns
            do memberIndex = 1, nEnsUsed
              weightsMembersLatLon(memberIndex,:,latLonIndex) =  &
                   weightsMembersLatLon(memberIndex,:,latLonIndex) - &
                   sum(weightsMembersLatLon(memberIndex,:,latLonIndex))/real(nEns,8)
            end do

          else if (trim(algorithm) == 'CVLETKF-PERTOBS') then
            !
            ! Weight calculation for perturbed-obs cross-validation LETKF algorithm
            !

            ! Compute eigenValues/Vectors of Yb^T R^-1 Yb = E * Lambda * E^T
            call utl_tmg_start(104,'------EigenDecomp')
            tolerance = 1.0D-50
            call utl_eigenDecomp(YbTinvRYb, eigenValues, eigenVectors, tolerance, matrixRank)
            call utl_tmg_stop(104)
            !if (matrixRank < (nEns-1)) then
            !  write(*,*) 'YbTinvRYb is rank deficient =', matrixRank, nEns, numLocalObs
            !end if

            ! Compute ensemble mean local weights as E * (Lambda + (Nens-1)*I)^-1 * E^T * YbTinvR * (obs - meanYb)
            weightsTemp(:) = 0.0d0
            do localObsIndex = 1, numLocalObs
              bodyIndex = localBodyIndices(localObsIndex)
              do memberIndex = 1, nEns
                weightsTemp(memberIndex) = weightsTemp(memberIndex) +   &
                                           YbTinvR(memberIndex,localObsIndex) *  &
                                           ( ensObs_mpiglobal%obsValue(bodyIndex) - &
                                             ensObs_mpiglobal%meanYb(bodyIndex) )
              end do
            end do
            weightsTemp2(:) = 0.0d0
            do memberIndex2 = 1, matrixRank
              do memberIndex1 = 1, nEns
                weightsTemp2(memberIndex2) = weightsTemp2(memberIndex2) +   &
                                             eigenVectors(memberIndex1,memberIndex2) *  &
                                             weightsTemp(memberIndex1)
              end do
            end do
            do memberIndex = 1, matrixRank
              weightsTemp2(memberIndex) = weightsTemp2(memberIndex) *  &
                                          1.0D0/(eigenValues(memberIndex) + real(nEns - 1,8))
            end do
            weightsMeanLatLon(:,1,latLonIndex) = 0.0d0
            do memberIndex2 = 1, matrixRank
              do memberIndex1 = 1, nEns
                weightsMeanLatLon(memberIndex1,1,latLonIndex) =  &
                     weightsMeanLatLon(memberIndex1,1,latLonIndex) +   &
                     eigenVectors(memberIndex1,memberIndex2) *  &
                     weightsTemp2(memberIndex2)
              end do
            end do

            ! Compute ensemble perturbation weights using mean increment weights 
            ! formula, but with subset of members: 
            ! wa_i = I_i + E * (Lambda + (Nens-1)*I)^-1 * E^T * YbTinvR * (obs + randpert_i - Yb_i)
            ! Wa   = wa_i - mean_over_i(wa_i) 
            !
            ! Loop over sub-ensembles
            do subEnsIndex = 1, numSubEns

              ! Use complement (independent) ens to get eigenValues/Vectors of Yb^T R^-1 Yb = E*Lambda*E^T
              call utl_tmg_start(104,'------EigenDecomp')
              do memberIndexCV2 = 1, nEnsIndependentPerSubEns
                memberIndex2 = memberIndexSubEnsComp(memberIndexCV2, subEnsIndex)
                do memberIndexCV1 = 1, nEnsIndependentPerSubEns
                  memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
                  YbTinvRYb_CV(memberIndexCV1,memberIndexCV2) = YbTinvRYb(memberIndex1,memberIndex2)
                end do
              end do
              tolerance = 1.0D-50
              call utl_eigenDecomp(YbTinvRYb_CV, eigenValues_CV, eigenVectors_CV, tolerance, matrixRank)
              call utl_tmg_stop(104)

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
                         YbTinvR(memberIndex1,localObsIndex) *  &
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
                weightsMembersLatLon(:,memberIndex,latLonIndex) = 0.0d0
                do memberIndex2 = 1, matrixRank
                  do memberIndexCV1 = 1, nEnsIndependentPerSubEns
                    memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
                    weightsMembersLatLon(memberIndex1,memberIndex,latLonIndex) =  &
                         weightsMembersLatLon(memberIndex1,memberIndex,latLonIndex) +   &
                         eigenVectors_CV(memberIndexCV1,memberIndex2) *  &
                         weightsTemp2(memberIndex2)
                  end do
                end do

                ! I + previous_result
                weightsMembersLatLon(memberIndex,memberIndex,latLonIndex) =  &
                     1.0D0 + weightsMembersLatLon(memberIndex,memberIndex,latLonIndex)

              end do ! memberIndexCV
            end do ! subEnsIndex

            ! Remove the weights mean computed over the columns
            do memberIndex = 1, nEns
              weightsMembersLatLon(memberIndex,:,latLonIndex) =  &
                   weightsMembersLatLon(memberIndex,:,latLonIndex) - &
                   sum(weightsMembersLatLon(memberIndex,:,latLonIndex))/real(nEns,8)
            end do

          else

            call utl_abort('UNKNOWN LETKF ALGORITHM')

          end if

!latIndex = myLatIndexesSend(latLonIndex)
!lonIndex = myLonIndexesSend(latLonIndex)
!latIndexPrint = 0
!lonIndexPrint = 0 
!if ( printBeforeComm .and. levIndex == 47 ) then
!  printBeforeComm = .false. 
!  latIndexPrint = latIndex
!  lonIndexPrint = lonIndex
!  write(*,*) 'maziar: before comm, numLocalObs=', numLocalObs, ', latIndexPrint=', latIndexPrint, &
!             ', lonIndexPrint=', lonIndexPrint
!  do memberIndex2 = 1, nEns
!    do memberIndex1 = 1, nEns
!      do eigenVectorColumnIndex = 1, numRetainedEigen
!        memberIndexInModEns = (eigenVectorColumnIndex - 1) * nEns + memberIndex1
!        write(*,*) 'maziar: before comm, memberIndex2=', memberIndex2, ', memberIndex1=', memberIndex1, &
!        ', eigenVectorColumnIndex=', eigenVectorColumnIndex, &
!        ', memberIndexInModEns=', memberIndexInModEns, &
!        ', weightsMembersLatLon=', weightsMembersLatLon(memberIndexInModEns,memberIndex2,latLonIndex)
!      end do
!    end do
!  end do
!end if

        else

          !! no obs near this grid point, mean weights zero, member weights identity
          !weightsMeanLatLon(:,1,latLonIndex) = 0.0d0
          !weightsMembersLatLon(:,:,latLonIndex) = 0.0d0
          !do memberIndex = 1, nEns
          !  if ( numRetainedEigen > 0 ) then
          !    do eigenVectorColumnIndex = 1, numRetainedEigen 
          !      memberIndexInModEns = (eigenVectorColumnIndex - 1) * nEns + memberIndex
          !      weightsMembersLatLon(memberIndexInModEns,memberIndex,latLonIndex) = 1.0d0
          !    end do
          !  else
          !    weightsMembersLatLon(memberIndex,memberIndex,latLonIndex) = 1.0d0
          !  end if
          !end do

        end if ! numLocalObs > 0

        call utl_tmg_stop(103)

        !
        ! Now post all send instructions (each lat-lon may be sent to multiple tasks)
        !
        call utl_tmg_start(101,'----CommWeights')
        latIndex = myLatIndexesSend(latLonIndex)
        lonIndex = myLonIndexesSend(latLonIndex)
        do procIndex = 1, myNumProcIndexesSend(latLonIndex)
          sendTag = latLonTagMpiGlobal(lonIndex,latIndex)
          procIndexSend = myProcIndexesSend(latLonIndex, procIndex)

          nsize = nEnsUsed
          numSend = numSend + 1
          call mpi_isend( weightsMeanLatLon(:,1,latLonIndex),  &
                          nsize, mmpi_datyp_real8, procIndexSend-1, sendTag,  &
                          mmpi_comm_grid, requestIdSend(numSend), ierr )
          nsize = nEnsUsed * nEns
          numSend = numSend + 1
          sendTag = sendTag + maxval(latLonTagMpiGlobal(:,:))
          call mpi_isend( weightsMembersLatLon(:,:,latLonIndex),  &
                          nsize, mmpi_datyp_real8, procIndexSend-1, sendTag,  &
                          mmpi_comm_grid, requestIdSend(numSend), ierr )
        end do
        call utl_tmg_stop(101)

      end do LATLON_LOOP

      !
      ! Wait for communiations to finish before continuing
      !
      call utl_tmg_start(101,'----CommWeights')
      if (firstTime) write(*,*) 'numSend/Recv = ', numSend, numRecv
      firstTime = .false.

      if ( numRecv > 0 ) then
        call mpi_waitAll(numRecv, requestIdRecv(1:numRecv), MPI_STATUSES_IGNORE, ierr)
      end if

      if ( numSend > 0 ) then
        call mpi_waitAll(numSend, requestIdSend(1:numSend), MPI_STATUSES_IGNORE, ierr)
      end if

      call utl_tmg_stop(101)

      !
      ! Interpolate weights from coarse to full resolution
      !
      call utl_tmg_start(106,'----InterpolateWeights')
      !if (wInterpInfo%latLonStep > 1) then
      !  call enkf_interpWeights(wInterpInfo, weightsMean)
      !  call enkf_interpWeights(wInterpInfo, weightsMembers)
      !end if
      call utl_tmg_stop(106)

      call utl_tmg_start(107,'----ApplyWeights')

      if ( levIndex == 47 .and. numRetainedEigen > 0 ) then
      !if ( levIndex == 47 .and. numRetainedEigen > 0 .and. mpi_myid == 0 ) then
      !if ( numRetainedEigen > 0 ) then
        do latIndex = myLatBegHalo, myLatEndHalo
          do lonIndex = myLonBegHalo, myLonEndHalo
            ! If this lat-lon is to be interpolated, then skip calculation
            if (wInterpInfo%numIndexes(lonIndex,latIndex) > 0) cycle

            anlLat = hco_ens%lat2d_4(lonIndex,latIndex) * MPC_DEGREES_PER_RADIAN_R8 
            anlLon = hco_ens%lon2d_4(lonIndex,latIndex) * MPC_DEGREES_PER_RADIAN_R8
            if ( abs(anlLat+43.0) < 1 .and. abs(anlLon-75.0) < 1 ) then
              if ( latIndex == myLatBegHalo .and. lonIndex == myLonBegHalo ) &
                write(*,*) 'maziar: max(anlLat)=', maxval(hco_ens%lat2d_4) * MPC_DEGREES_PER_RADIAN_R8, &
                ', min(anlLat)=', minval(hco_ens%lat2d_4) * MPC_DEGREES_PER_RADIAN_R8, &
                ', max(anlLon)=', maxval(hco_ens%lon2d_4) * MPC_DEGREES_PER_RADIAN_R8, &
                ', min(anlLon)=', minval(hco_ens%lon2d_4) * MPC_DEGREES_PER_RADIAN_R8

              write(*,*) 'maziar: anlLat=', anlLat, ', anlLon=', anlLon
              do memberIndex2 = 1, nEns
                do memberIndex1 = 1, nEns
                  do eigenVectorColumnIndex = 1, numRetainedEigen
                    memberIndexInModEns = (eigenVectorColumnIndex - 1) * nEns + memberIndex1

                    if ( memberIndex1 == 1 ) then
                      write(*,*) 'maziar: weightsMean anlLat/anlLon, memberIndex2=', memberIndex2, &
                      ', eigenVectorColumnIndex=', eigenVectorColumnIndex, &
                      ', memberIndexInModEns=', memberIndexInModEns, &
                      ', weightsMean=', weightsMean(memberIndexInModEns,1,lonIndex,latIndex)
                    end if

                    write(*,*) 'maziar: weightsMember anlLat/anlLon, memberIndex2=', memberIndex2, ', memberIndex1=', memberIndex1, &
                    ', eigenVectorColumnIndex=', eigenVectorColumnIndex, &
                    ', memberIndexInModEns=', memberIndexInModEns, &
                    ', weightsMembers=', weightsMembers(memberIndexInModEns,memberIndex2,lonIndex,latIndex)
                  end do
                end do
              end do
            end if

          end do 
        end do
      end if

      !
      ! Apply the weights to compute the ensemble mean and members
      !
      call gsv_getField(stateVectorMeanInc,meanInc_ptr_r4)
      call gsv_getField(stateVectorMeanTrl,meanTrl_ptr_r4)
      call gsv_getField(stateVectorMeanAnl,meanAnl_ptr_r4)

      !$OMP PARALLEL DO PRIVATE(latIndex, lonIndex, varLevIndex, varLevel, varKind, levIndex2, memberTrl_ptr_r4, memberAnl_ptr_r4), &
      !$OMP PRIVATE(memberAnlPert, stepIndex, memberIndex, memberIndex2, memberIndex1)
      do latIndex = myLatBeg, myLatEnd
        LON_LOOP5: do lonIndex = myLonBeg, myLonEnd

          ! skip this grid point if all weights zero (no nearby obs)
          if (all(weightsMean(:,1,lonIndex,latIndex) == 0.0d0)) cycle LON_LOOP5

          ! Compute the ensemble mean increment and analysis
          do varLevIndex = 1, numVarLev
            ! Only treat varLevIndex values that correspond with current levIndex
            varLevel = vnl_varLevelFromVarname(gsv_getVarNameFromK(stateVectorMeanInc,varLevIndex))
            if (varLevel == 'SF'   .or. varLevel == 'SFMM' .or. &
                varLevel == 'SFTH' .or. varLevel == 'SS') then
              varKind = vnl_varKindFromVarname(gsv_getVarNameFromK(stateVectorMeanInc,varLevIndex))
              if (varKind == 'OC') then
                levIndex2 = 1
              else
                levIndex2 = nLev_weights
              end if
            else if (varLevel == 'MM' .or. varLevel == 'TH' .or. varLevel == 'DP') then
              levIndex2 = gsv_getLevFromK(stateVectorMeanInc,varLevIndex)
            else if (varLevel == 'OT') then
              ! Most (all?) variables using the 'other' coordinate are surface
              levIndex2 = nLev_weights
            else
              write(*,*) 'varLevel = ', varLevel
              call utl_abort('enkf_LETKFanalyses: unknown varLevel')
            end if
            if (levIndex2 /= levIndex) cycle
            memberTrl_ptr_r4 => ens_getOneLev_r4(ensembleTrl,varLevIndex)
            nLev = gsv_getNumLevFromVarName( stateVectorMeanInc, &
                                             gsv_getVarNameFromK(stateVectorMeanInc,varLevIndex) )
            do stepIndex = 1, tim_nstepobsinc
              ! mean increment
              if ( numRetainedEigen > 0 ) then
                do memberIndex = 1, nEns
                  do eigenVectorColumnIndex = 1, numRetainedEigen

                    ! maziar: does this cover all cases?
                    if ( nLev == 1 ) then
                      eigenVectorLevelIndex = nLev
                    else
                      eigenVectorLevelIndex = levIndex
                    end if

                    call getModulationFactor( stateVectorMeanInc%vco, eigenVectorLevelIndex, &
                                              eigenVectorColumnIndex, numRetainedEigen, &
                                              nEns, vLocalize, &
                                              modulationFactor )
                    pert_r4 = memberTrl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) -  &
                                        meanTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex)

                    !if ( latIndex == myLatBeg .and. lonIndex == myLonBeg .and. debug ) then
                    !  write(*,*) 'maziar: MEAN INCR, varName=', gsv_getVarNameFromK(stateVectorMeanInc,varLevIndex),', levIndex=', levIndex, ', stepIndex=', stepIndex
                    !  write(*,*) 'maziar: modulationFactor=', modulationFactor
                    !  write(*,*) 'maziar: mean trial=', meanTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex)
                    !  write(*,*) 'maziar: original member pert=', pert_r4
                    !end if
!if ( latIndex == myLatBeg .and. lonIndex == myLonBeg .and. stepIndex == 1 .and. &
!    memberIndex == 1 .and. trim(gsv_getVarNameFromK(stateVectorMeanInc,varLevIndex)) == 'TT' ) then
!  write(*,*) 'maziar: levIndex=', levIndex, &
!            ', eigenVectorColumnIndex=', eigenVectorColumnIndex, &
!            ', modulationFactor=', modulationFactor
!end if


                    pert_r4 = pert_r4 * real(modulationFactor,4)

                    memberIndexInModEns = (eigenVectorColumnIndex - 1) * nEns + memberIndex
                    meanInc_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) =  &
                         meanInc_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) +  &
                         weightsMean(memberIndexInModEns,1,lonIndex,latIndex) * pert_r4

                    !if ( latIndex == myLatBeg .and. lonIndex == myLonBeg .and. debug ) then
                    !  write(*,*) 'maziar: modulated member pert=', pert_r4
                    !  write(*,*) 'maziar: memberIndexInModEns=',  memberIndexInModEns
                    !  write(*,*) 'maziar: mean increment=',  meanInc_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex)
                    !end if

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
          do varLevIndex = 1, numVarLev
            ! Only treat varLevIndex values that correspond with current levIndex
            varLevel = vnl_varLevelFromVarname(gsv_getVarNameFromK(stateVectorMeanInc,varLevIndex))
            if (varLevel == 'SF'   .or. varLevel == 'SFMM' .or. &
                varLevel == 'SFTH' .or. varLevel == 'SS') then
              varKind = vnl_varKindFromVarname(gsv_getVarNameFromK(stateVectorMeanInc,varLevIndex))
              if (varKind == 'OC') then
                levIndex2 = 1
              else
                levIndex2 = nLev_weights
              end if
            else if (varLevel == 'MM' .or. varLevel == 'TH' .or. varLevel == 'DP') then
              levIndex2 = gsv_getLevFromK(stateVectorMeanInc,varLevIndex)
            else if (varLevel == 'OT') then
              ! Most (all?) variables using the 'other' coordinate are surface
              levIndex2 = nLev_weights
            else
              write(*,*) 'varLevel = ', varLevel
              call utl_abort('enkf_LETKFanalyses: unknown varLevel')
            end if
            if (levIndex2 /= levIndex) cycle
            memberTrl_ptr_r4 => ens_getOneLev_r4(ensembleTrl,varLevIndex)
            memberAnl_ptr_r4 => ens_getOneLev_r4(ensembleAnl,varLevIndex)
            nLev = gsv_getNumLevFromVarName( stateVectorMeanInc, &
                                             gsv_getVarNameFromK(stateVectorMeanInc,varLevIndex) )
            do stepIndex = 1, tim_nstepobsinc
              ! Compute analysis member perturbation
              memberAnlPert(:) = 0.0d0

              anlLat = hco_ens%lat2d_4(lonIndex,latIndex) * MPC_DEGREES_PER_RADIAN_R8 
              anlLon = hco_ens%lon2d_4(lonIndex,latIndex) * MPC_DEGREES_PER_RADIAN_R8

              if ( numRetainedEigen > 0 ) then
                do memberIndex2 = 1, nEns
                  do memberIndex1 = 1, nEns
                    do eigenVectorColumnIndex = 1, numRetainedEigen

                      ! Index of the modulated ensemble member corresponding to original
                      ! ensemble member index (memberIndex1) and eigenVectorColumnIndex.
                      memberIndexInModEns = (eigenVectorColumnIndex - 1) * nEns + memberIndex1

                      ! maziar: does this cover all cases?
                      if ( nLev == 1 ) then
                        eigenVectorLevelIndex = nLev
                      else
                        eigenVectorLevelIndex = levIndex
                      end if

                      call getModulationFactor( stateVectorMeanInc%vco, eigenVectorLevelIndex, &
                                                eigenVectorColumnIndex, numRetainedEigen, &
                                                nEns, vLocalize, &
                                                modulationFactor )

                      ! Compute background ensemble perturbations for the modulated ensemble (Xb_Mod)
                      pert_r4 = memberTrl_ptr_r4(memberIndex1,stepIndex,lonIndex,latIndex) -  &
                                          meanTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex)
                      pert_r4 = pert_r4 * real(modulationFactor,4)
                      
                      if ( abs(anlLat+43.0) < 1 .and. abs(anlLon-75.0) < 1 .and. &
                        trim(gsv_getVarNameFromK(stateVectorMeanInc,varLevIndex)) == 'TT' ) then
                        write(*,*) 'maziar: after comm, levIndex=', levIndex, &
                        ', memberIndex2=', memberIndex2, ', memberIndex1=', memberIndex1, &
                        ', eigenVectorColumnIndex=', eigenVectorColumnIndex, &
                        ', memberIndexInModEns=', memberIndexInModEns, &
                        ', modulationFactor=', modulationFactor, &
                        ', pert=', pert_r4, &
                        ', weightsMembers=', weightsMembers(memberIndexInModEns,memberIndex2,lonIndex,latIndex)
                      end if

                      ! sum Xb_Mod * Wa over all modulated ensembles to get member perturbations for
                      !   original ensemble (memberIndex2)
                      memberAnlPert(memberIndex2) = memberAnlPert(memberIndex2) + &
                           weightsMembers(memberIndexInModEns,memberIndex2,lonIndex,latIndex) *  pert_r4
                    end do ! eigenVectorIndex
                  end do ! memberIndex1

                  ! Compute final member perturbations by removing background original ensemble perturbations
                  memberAnlPert(memberIndex2) = (memberTrl_ptr_r4(memberIndex2,stepIndex,lonIndex,latIndex) -  &
                                                 meanTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex)) + &
                                                 memberAnlPert(memberIndex2)

                end do ! memberIndex2


                !if ( latIndex == myLatBeg .and. lonIndex == myLonBeg .and. debug ) then
                !  write(*,*) 'maziar: Anal pert, varName=', gsv_getVarNameFromK(stateVectorMeanInc,varLevIndex),', levIndex=', levIndex, ', stepIndex=', stepIndex
                !  write(*,*) 'maziar: after remove mean weight sum(memberAnlPert)=', sum(memberAnlPert(:))
                !end if

              else
                do memberIndex2 = 1, nEns
                  do memberIndex1 = 1, nEns
                    if ( abs(anlLat+43.0) < 1 .and. abs(anlLon-75.0) < 1 .and. &
                      trim(gsv_getVarNameFromK(stateVectorMeanInc,varLevIndex)) == 'TT' ) then
                      write(*,*) 'maziar: after comm, levIndex=', levIndex, &
                      ', memberIndex2=', memberIndex2, ', memberIndex1=', memberIndex1, &
                      ', pert=', memberTrl_ptr_r4(memberIndex1,stepIndex,lonIndex,latIndex) -  &
                      meanTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex), &
                      ', weightsMembers=', weightsMembers(memberIndex1,memberIndex2,lonIndex,latIndex)
                    end if
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

        end do LON_LOOP5
      end do
      !$OMP END PARALLEL DO

      call utl_tmg_stop(107)

    end do LEV_LOOP

    if (countMaxExceeded > 0) then
      write(*,*) 'enkf_LETKFanalyses: WARNING: Found more local obs than specified max number at ', &
                 real(100*countMaxExceeded)/real(numGridPointWeights), '% of grid points.'
      write(*,*) '                      Maximum number found was ', maxCountMaxExceeded,  &
                 ' which is greater than specified number ', maxNumLocalObs
      write(*,*) '                      Therefore will keep closest obs only.'
    end if

    call utl_tmg_start(108,'----Barr')
    call rpn_comm_barrier('GRID',ierr)
    call utl_tmg_stop(108)

    call gsv_deallocate(stateVectorMeanInc)
    call gsv_deallocate(stateVectorMeanTrl)

    write(*,*) 'enkf_LETKFanalyses: done'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    call utl_tmg_stop(100)

  end subroutine enkf_LETKFanalyses

  !----------------------------------------------------------------------
  ! enkf_computeVertLocation (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_computeVertLocation(vertLocation_r4,stateVectorMeanTrl)
    !
    !:Purpose:  Compute extract global 3D vertical location field from supplied
    !           stateVector. Can be either logPressure or depth levels.
    !
    implicit none

    ! Arguments
    real(4), allocatable, intent(inout) :: vertLocation_r4(:,:,:)
    type(struct_gsv),     intent(inout) :: stateVectorMeanTrl

    ! Locals
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
  ! enkf_setupMpiDistribution (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_LETKFsetupMpiDistribution(myNumLatLonRecv, myNumLatLonSend, &
                                            myLatIndexesRecv, myLonIndexesRecv, &
                                            myLatIndexesSend, myLonIndexesSend, &
                                            myProcIndexesRecv, myProcIndexesSend, &
                                            myNumProcIndexesSend, mpiDistribution, wInterpInfo)
    !
    ! :Purpose: Setup for distribution of grid points over mpi tasks.
    !
    implicit none

    ! Arguments
    integer              :: myNumLatLonRecv, myNumLatLonSend
    integer, allocatable :: myLatIndexesRecv(:), myLonIndexesRecv(:)
    integer, allocatable :: myLatIndexesSend(:), myLonIndexesSend(:)
    integer, allocatable :: myProcIndexesRecv(:), myProcIndexesSend(:,:)
    integer, allocatable :: myNumProcIndexesSend(:)
    character(len=*)     :: mpiDistribution
    type(struct_enkfInterpInfo) :: wInterpInfo

    ! Locals
    integer :: latIndex, lonIndex, procIndex, procIndexSend, latLonIndex
    integer :: myLonBegHalo, myLonEndHalo, myLatBegHalo, myLatEndHalo
    integer :: numLatLonRecvMax, numLatLonTotalUnique, ierr
    integer, allocatable :: allLatIndexesRecv(:,:), allLonIndexesRecv(:,:)
    integer, allocatable :: allLatIndexesSend(:,:), allLonIndexesSend(:,:)
    integer, allocatable :: allNumLatLonRecv(:), allNumLatLonSend(:)

    myLonBegHalo = wInterpInfo%myLonBegHalo
    myLonEndHalo = wInterpInfo%myLonEndHalo
    myLatBegHalo = wInterpInfo%myLatBegHalo
    myLatEndHalo = wInterpInfo%myLatEndHalo

    write(*,*) 'enkf_LETKFsetupMpiDistribution: starting'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if (trim(mpiDistribution) == 'TILES') then

      ! First, determine number of grid points needed locally (for recv-ing)
      myNumLatLonRecv = 0
      do latIndex = myLatBegHalo, myLatEndHalo
        LON_LOOP0: do lonIndex = myLonBegHalo, myLonEndHalo
          ! If this lat-lon is to be interpolated, then skip calculation
          if (wInterpInfo%numIndexes(lonIndex,latIndex) > 0) cycle LON_LOOP0
          myNumLatLonRecv = myNumLatLonRecv + 1
        end do LON_LOOP0
      end do

      write(*,*) 'enkf_LETKFsetupMpiDistribution: myNumLatLonRecv =', myNumLatLonRecv

      ! Determine list of grid point indexes where weights needed locally (for recv-ing)
      allocate(myLatIndexesRecv(myNumLatLonRecv))
      allocate(myLonIndexesRecv(myNumLatLonRecv))
      allocate(myProcIndexesRecv(myNumLatLonRecv))
      myNumLatLonRecv = 0
      do latIndex = myLatBegHalo, myLatEndHalo
        LON_LOOP1: do lonIndex = myLonBegHalo, myLonEndHalo
          ! If this lat-lon is to be interpolated, then skip calculation
          if (wInterpInfo%numIndexes(lonIndex,latIndex) > 0) cycle LON_LOOP1
          myNumLatLonRecv = myNumLatLonRecv + 1

          myLatIndexesRecv(myNumLatLonRecv) = latIndex
          myLonIndexesRecv(myNumLatLonRecv) = lonIndex
          myProcIndexesRecv(myNumLatLonRecv) = mmpi_myid+1
        end do LON_LOOP1
      end do

      ! No communication, so send info equals recv info
      myNumLatLonSend = myNumLatLonRecv
      allocate(myLatIndexesSend(myNumLatLonSend))
      allocate(myLonIndexesSend(myNumLatLonSend))
      allocate(myProcIndexesSend(myNumLatLonSend,1))
      allocate(myNumProcIndexesSend(myNumLatLonSend))

      myLatIndexesSend(:) = myLatIndexesRecv(:)
      myLonIndexesSend(:) = myLonIndexesRecv(:)
      myProcIndexesSend(:,1) = myProcIndexesRecv(:)
      myNumProcIndexesSend(:) = 1

    else if (trim(mpiDistribution) == 'ROUNDROBIN') then

      ! First, determine number of grid points needed locally (for recv-ing)
      myNumLatLonRecv = 0
      do latIndex = myLatBegHalo, myLatEndHalo
        LON_LOOP2: do lonIndex = myLonBegHalo, myLonEndHalo
          ! If this lat-lon is to be interpolated, then skip calculation
          if (wInterpInfo%numIndexes(lonIndex,latIndex) > 0) cycle LON_LOOP2
          myNumLatLonRecv = myNumLatLonRecv + 1
        end do LON_LOOP2
      end do

      ! Communicate to all mpi tasks
      allocate(allNumLatLonRecv(mmpi_nprocs))
      call rpn_comm_allgather(myNumLatLonRecv, 1, "mpi_integer",  &
                              allNumLatLonRecv, 1,"mpi_integer", "GRID", ierr)
      numLatLonRecvMax = maxval(allNumLatLonRecv)
      write(*,*) 'enkf_LETKFsetupMpiDistribution: allNumLatLonRecv =', allNumLatLonRecv(:)
      write(*,*) 'enkf_LETKFsetupMpiDistribution: numLatLonRecvSum =', sum(allNumLatLonRecv)
      write(*,*) 'enkf_LETKFsetupMpiDistribution: numLatLonRecvMax =', numLatLonRecvMax

      ! Determine list of grid point indexes where weights needed locally (for recv-ing)
      allocate(myLatIndexesRecv(numLatLonRecvMax))
      allocate(myLonIndexesRecv(numLatLonRecvMax))
      allocate(myProcIndexesRecv(numLatLonRecvMax))
      myLatIndexesRecv(:) = -1
      myLonIndexesRecv(:) = -1
      myProcIndexesRecv(:) = -1
      myNumLatLonRecv = 0
      do latIndex = myLatBegHalo, myLatEndHalo
        LON_LOOP3: do lonIndex = myLonBegHalo, myLonEndHalo
          ! If this lat-lon is to be interpolated, then skip calculation
          if (wInterpInfo%numIndexes(lonIndex,latIndex) > 0) cycle LON_LOOP3
          myNumLatLonRecv = myNumLatLonRecv + 1

          myLatIndexesRecv(myNumLatLonRecv) = latIndex
          myLonIndexesRecv(myNumLatLonRecv) = lonIndex
        end do LON_LOOP3
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

      ! From these lat-lon lists, create unique master list of all grid points where weights computed
      ! and assign to mpi tasks for doing the calculation and for send-ing
      allocate(myLatIndexesSend(numLatLonRecvMax))
      allocate(myLonIndexesSend(numLatLonRecvMax))
      myLatIndexesSend(:) = -1
      myLonIndexesSend(:) = -1
      numLatLonTotalUnique = 0
      myNumLatLonSend = 0
      do procIndex = 1, mmpi_nprocs
        WEIGHTS1LEV_LOOP: do latLonIndex = 1, allNumLatLonRecv(procIndex)
          if (enkf_latLonAlreadyFound(allLatIndexesRecv, allLonIndexesRecv, latLonIndex, procIndex)) &
               cycle WEIGHTS1LEV_LOOP
          ! Count the total number of weights
          numLatLonTotalUnique = numLatLonTotalUnique + 1

          ! Round-robin distribution of master list across mpi tasks
          procIndexSend = 1 + mod(numLatLonTotalUnique-1, mmpi_nprocs)

          ! Store the lat-lon indexes of the weights I am responsible for
          if (procIndexSend == (mmpi_myid+1)) then
            myNumLatLonSend = myNumLatLonSend + 1
            myLatIndexesSend(myNumLatLonSend) =  &
                 allLatIndexesRecv(latLonIndex, procIndex)
            myLonIndexesSend(myNumLatLonSend) =  &
                 allLonIndexesRecv(latLonIndex, procIndex)
          end if
        end do WEIGHTS1LEV_LOOP
      end do
      write(*,*) 'enkf_LETKFsetupMpiDistribution: number of lat/lon points where weights to be computed =',  &
                 numLatLonTotalUnique

      ! Communicate to all mpi tasks this list of grid point lat-lon indexes
      allocate(allNumLatLonSend(mmpi_nprocs))
      call rpn_comm_allgather(myNumLatLonSend, 1, "mpi_integer",  &
                              allNumLatLonSend, 1,"mpi_integer", "GRID", ierr)
      allocate(allLatIndexesSend(numLatLonRecvMax, mmpi_nprocs))
      allocate(allLonIndexesSend(numLatLonRecvMax, mmpi_nprocs))
      call rpn_comm_allgather(myLatIndexesSend, numLatLonRecvMax, "mpi_integer",  &
                              allLatIndexesSend, numLatLonRecvMax, "mpi_integer",  &
                              "GRID", ierr)
      call rpn_comm_allgather(myLonIndexesSend, numLatLonRecvMax, "mpi_integer",  &
                              allLonIndexesSend, numLatLonRecvMax, "mpi_integer",  &
                              "GRID", ierr)

      ! Figure out which mpi tasks I will need to send my results to
      allocate(myProcIndexesSend(myNumLatLonSend,mmpi_nprocs))
      allocate(myNumProcIndexesSend(myNumLatLonSend))
      myProcIndexesSend(:,:) = -1
      myNumProcIndexesSend(:) = 0
      do latLonIndex = 1, myNumLatLonSend
        do procIndex = 1, mmpi_nprocs
          if ( any( (myLatIndexesSend(latLonIndex) == allLatIndexesRecv(1:allNumLatLonRecv(procIndex), procIndex)) .and.  &
                    (myLonIndexesSend(latLonIndex) == allLonIndexesRecv(1:allNumLatLonRecv(procIndex), procIndex)) ) ) then
            myNumProcIndexesSend(latLonIndex) = myNumProcIndexesSend(latLonIndex) + 1
            myProcIndexesSend(latLonIndex,myNumProcIndexesSend(latLonIndex)) = procIndex
          end if
        end do
      end do

      ! Figure out which mpi tasks I will receive the results from
      do latLonIndex = 1, myNumLatLonRecv
        do procIndex = 1, mmpi_nprocs
          if ( any( (myLatIndexesRecv(latLonIndex) == allLatIndexesSend(1:allNumLatLonSend(procIndex), procIndex)) .and.  &
                    (myLonIndexesRecv(latLonIndex) == allLonIndexesSend(1:allNumLatLonSend(procIndex), procIndex)) ) ) then
            myProcIndexesRecv(latLonIndex) = procIndex
          end if
        end do
      end do

    else
      call utl_abort('enkf_LETKFsetupMpiDistribution: unknown MPI distribution selected')
    end if

    write(*,*) 'enkf_LETKFsetupMpiDistribution: lat/lon/proc indexes I need to receive:'
    do latLonIndex = 1, myNumLatLonRecv
      write(*,*) myLatIndexesRecv(latLonIndex), myLonIndexesRecv(latLonIndex),  &
                 myProcIndexesRecv(latLonIndex)
    end do

    write(*,*) 'enkf_LETKFsetupMpiDistribution: number of lat/lon indexes I am responsible for =', myNumLatLonSend
    write(*,*) 'enkf_LETKFsetupMpiDistribution: the lat/lon/proc indexes I am responsible for:'
    do latLonIndex = 1, myNumLatLonSend
      write(*,*) myLatIndexesSend(latLonIndex), myLonIndexesSend(latLonIndex),  &
                 myProcIndexesSend(latLonIndex,1:myNumProcIndexesSend(latLonIndex))
    end do

    write(*,*) 'enkf_LETKFsetupMpiDistribution: done'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine enkf_LETKFsetupMpiDistribution

  !----------------------------------------------------------------------
  ! enkf_LETKFgetMpiGlobalTags (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_LETKFgetMpiGlobalTags(latLonTagMpiGlobal,myLatIndexesRecv,myLonIndexesRecv)
    implicit none
    ! Arguments:
    integer :: latLonTagMpiGlobal(:,:)
    integer :: myLatIndexesRecv(:)
    integer :: myLonIndexesRecv(:)

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
    
    write(*,*) 'enkf_LETKFgetMpiGlobalTags: Finished'

  end subroutine enkf_LETKFgetMpiGlobalTags

  !----------------------------------------------------------------------
  ! enkf_latLonAlreadyFound (private function)
  !----------------------------------------------------------------------
  function enkf_latLonAlreadyFound(allLatIndexesRecv, allLonIndexesRecv, latLonIndex, procIndex) result(found)
    implicit none
    ! Arguments:
    integer :: allLatIndexesRecv(:,:), allLonIndexesRecv(:,:)
    integer :: latLonIndex, procIndex
    logical :: found

    ! Locals:
    integer :: latLonIndex2, procIndex2, numLatLonRecvMax

    numLatLonRecvMax = size(allLatIndexesRecv, 1)

    ! check on all previous mpi tasks if this lat/lon has already been encountered
    found = .false.
    do procIndex2 = 1, procIndex-1
      WEIGHTS1LEV_LOOP2: do latLonIndex2 = 1, numLatLonRecvMax
        if (allLatIndexesRecv(latLonIndex2, procIndex2) < 0) cycle WEIGHTS1LEV_LOOP2
        if ( (allLatIndexesRecv(latLonIndex, procIndex) ==  &
              allLatIndexesRecv(latLonIndex2, procIndex2)) .and.  &
             (allLonIndexesRecv(latLonIndex, procIndex) ==  &
              allLonIndexesRecv(latLonIndex2, procIndex2)) ) then
          found = .true.
          exit WEIGHTS1LEV_LOOP2
        end if
      end do WEIGHTS1LEV_LOOP2
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

    ! Arguments
    type(struct_enkfInterpInfo) :: wInterpInfo
    type(struct_hco) :: hco
    integer :: weightLatLonStep
    integer :: myLonBeg
    integer :: myLonEnd
    integer :: myLatBeg
    integer :: myLatEnd

    ! Locals
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

    ! Arguments
    type(struct_enkfInterpInfo) :: wInterpInfo
    real(8) :: weights(1:,1:,wInterpInfo%myLonBegHalo:,wInterpInfo%myLatBegHalo:)

    ! Locals
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

    !$OMP PARALLEL DO PRIVATE(latIndex, lonIndex, interpLatIndex, interpLonIndex, memberIndex1, memberIndex2)
    do latIndex = myLatBegHalo, myLatEndHalo
      do lonIndex = myLonBegHalo, myLonEndHalo
        if (wInterpInfo%numIndexes(lonIndex,latIndex) > 0) then

          ! Interpolation for ensemble member perturbation weight fields
          weights(:,:,lonIndex,latIndex) = 0.0D0
          if (wInterpInfo%lonIndexes(1,lonIndex,latIndex) == 0) cycle ! temporary until all interpolation setup completed
          do interpIndex = 1, wInterpInfo%numIndexes(lonIndex,latIndex)
            interpLonIndex = wInterpInfo%lonIndexes(interpIndex,lonIndex,latIndex)
            interpLatIndex = wInterpInfo%latIndexes(interpIndex,lonIndex,latIndex)

            totalCount(omp_get_thread_num()+1) = totalCount(omp_get_thread_num()+1) + 1
            do memberIndex2 = 1, numMembers2
              do memberIndex1 = 1, numMembers1
                weights(memberIndex1,memberIndex2,lonIndex,latIndex) =  &
                     weights(memberIndex1,memberIndex2,lonIndex,latIndex) + &
                     wInterpInfo%interpWeights(interpIndex,lonIndex,latIndex) *  &
                     weights(memberIndex1,memberIndex2,interpLonIndex,interpLatIndex)
              end do
            end do

          end do ! interpIndex
        end if ! numIndexes > 0
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

    ! arguments:
    type(struct_obs), target  :: obsSpaceData

    ! locals:
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

    ! arguments:
    type(struct_obs), target  :: obsSpaceData

    ! locals:
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
                                     eigenVectorColumnIndex, &
                                     stateVector_out, debug_opt )
    !
    !:Purpose: Compute vertical localization matrix, and the corresponding
    !          eigenvectors/eigenvalues, to obtain modulated stateVector.
    !
    implicit none

    ! Aguments:
    type(struct_gsv), intent(in) :: stateVector_in
    type(struct_gsv), intent(in) :: stateVectorMeanTrl
    real(8), intent(in) :: vLocalizeLengthScale
    integer, intent(in) :: numRetainedEigen
    integer, intent(in) :: nEns
    integer, intent(in) :: eigenVectorColumnIndex
    type(struct_gsv), intent(inout) :: stateVector_out
    logical, optional   :: debug_opt

    ! Locals:
    real(8)          :: modulationFactor
    real(4), pointer :: field_out_r4(:,:,:,:)
    real(4), pointer :: field_mean_r4(:,:,:,:)

    integer :: nLev, nlev_out, levIndex, latIndex, lonIndex
    integer :: lon1, lon2, lat1, lat2
    integer :: varIndex, stepIndex, eigenVectorLevelIndex

    character(len=4) :: varName

    write(*,*) 'enkf_getModulatedState: START'

    if ( stateVector_in%dataKind /= 4 ) then
      call utl_abort('enkf_getModulatedState: only dataKind=4 is implemented')
    end if

    nLev = stateVector_in%vco%nLev_M
    if ( vLocalizeLengthScale <= 0.0d0 .or. nLev <= 1 ) then
      call utl_abort('enkf_getModulatedState: no vertical localization')
    end if

    if ( present(debug_opt) ) then
      debug = debug_opt
    else
      debug = .false.
    end if

    ! Compute perturbation by subtracting ensMean
    call gsv_copy(stateVector_in, stateVector_out)
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
        if ( varName /= 'Z_T ' .and. varName /= 'Z_M ' .and. varName /= 'P_T ' .and. varName /= 'P_M ' ) call gsv_getField(stateVectorMeanTrl,field_mean_r4,varName)

        !!$OMP PARALLEL DO PRIVATE(latIndex,lonIndex,levIndex)
        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2
            do levIndex = 1, nlev_out

              if ( latIndex == lat1 .and. lonIndex == lon1 .and. debug ) then
                write(*,*) 'maziar: lat1=', lat1, ', lon1=', lon1
                write(*,*) 'maziar: varName=', varName,', levIndex=', levIndex, ', stepIndex=', stepIndex
                if ( varName /= 'Z_T ' .and. varName /= 'Z_M ' .and. varName /= 'P_T ' .and. varName /= 'P_M ' ) write(*,*) 'maziar: mean field=', field_mean_r4(lonIndex,latIndex,levIndex,stepIndex)
                write(*,*) 'maziar: original pert field=', field_out_r4(lonIndex,latIndex,levIndex,stepIndex)
              end if

              ! maziar: does this cover all cases?
              if ( nlev_out == 1 ) then
                eigenVectorLevelIndex = nLev
              else
                eigenVectorLevelIndex = levIndex
              end if

              call getModulationFactor( stateVector_in%vco, eigenVectorLevelIndex, &
                                        eigenVectorColumnIndex, numRetainedEigen, &
                                        nEns, vLocalizeLengthScale, &
                                        modulationFactor )

              if ( latIndex == lat1 .and. lonIndex == lon1 .and. debug ) then
                write(*,*) 'maziar: modulationFactor=', modulationFactor
              end if

              field_out_r4(lonIndex,latIndex,levIndex,stepIndex) = &
                                 field_out_r4(lonIndex,latIndex,levIndex,stepIndex) * &
                                 real(modulationFactor,4)

              if ( latIndex == lat1 .and. lonIndex == lon1 .and. debug ) then
                write(*,*) 'maziar: modulated pert field=', field_out_r4(lonIndex,latIndex,levIndex,stepIndex)
              end if

            end do
          end do
        end do
        !!$OMP END PARALLEL DO

      end do var_loop
    end do step_loop

    ! Now add to ensMean to get modulated member
    ! v_k = v'_k + v_mean
    call gsv_add(stateVectorMeanTrl, stateVector_out)

    write(*,*) 'enkf_getModulatedState: END'

  end subroutine enkf_getModulatedState

  !--------------------------------------------------------------------------
  ! getModulationFactor
  !--------------------------------------------------------------------------
  subroutine getModulationFactor( vco, eigenVectorLevelIndex, &
                                  eigenVectorColumnIndex, numRetainedEigen, &
                                  nEns, vLocalizeLengthScale, &
                                  modulationFactor )
    !
    !:Purpose: compute modulation factor needed to multiply ensemble
    !          perturbation to get the modulated perturbation:
    !          (Nens*nLambda/(Nens - 1))^1/2 * Lambda^1/2
    !
    implicit none

    ! Aguments:
    type(struct_vco), pointer, intent(in) :: vco
    integer, intent(in) :: eigenVectorLevelIndex
    integer, intent(in) :: eigenVectorColumnIndex
    integer, intent(in) :: numRetainedEigen
    integer, intent(in) :: nEns
    real(8), intent(in) :: vLocalizeLengthScale
    real(8), intent(out) :: modulationFactor

    ! Locals:
    integer             :: levIndex1, levIndex2, eigenIndex, status
    integer             :: nLev, matrixRank
    real(8)             :: zr, zcorr, pSurfRef
    real(8)             :: tolerance
    real(8), pointer    :: pressureProfile(:)
    real(8), allocatable, save :: eigenValues(:)
    real(8), allocatable, save :: eigenVectors(:,:)
    real(8), allocatable, save :: verticalLocalizationMat(:,:)
    real(8), allocatable, save :: verticalLocalizationMatLowRank(:,:)

    logical, save :: firstCall = .true.

    ! Compute vertical localization matrix and its eigenValues/Vectors on first call
    if ( firstCall ) then
      firstCall = .false.
      write(*,*) 'getModulationFactor: computing eigenValues/Vectors'

      nLev = vco%nLev_M

      allocate(eigenValues(nLev))
      allocate(eigenVectors(nLev,nLev))
      allocate(verticalLocalizationMat(nLev,nLev))
      allocate(verticalLocalizationMatLowRank(nLev,nLev))
      verticalLocalizationMatLowRank(:,:) = 0.0d0

      pSurfRef = 101000.D0
      nullify(pressureProfile)
      status = vgd_levels( vco%vgrid, &
                           ip1_list=vco%ip1_M, &
                           levels=pressureProfile, &
                           sfc_field=pSurfRef, &
                           in_log=.false. )
      if ( status /= VGD_OK ) then
        call utl_abort('getModulationFactor: ERROR from vgd_levels')
      end if

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

      !eigenValues(2:nLev) = 0.0d0
if (mpi_myid == 0) then
  do levIndex1 = 1, numRetainedEigen
    write(*,*) 'maziar: eigen mode=', levIndex1, ', eigenVectors=', eigenVectors(:,levIndex1)
  end do
  write(*,*) 'maziar: eigenValues=', eigenValues(1:numRetainedEigen)

  do levIndex1 = 1, nLev
    write(*,*) 'maziar: verticalLocalizationMat for lev ', levIndex1, '=', verticalLocalizationMat(levIndex1,:)
    write(*,*) 'maziar: verticalLocalizationMatLowRank for lev ', levIndex1, '=', verticalLocalizationMatLowRank(levIndex1,:)
  end do
end if
    end if

    modulationFactor = 1 / sqrt(verticalLocalizationMatLowRank(eigenVectorLevelIndex,eigenVectorLevelIndex)) * &
                       eigenVectors(eigenVectorLevelIndex,eigenVectorColumnIndex) * &
                       eigenValues(eigenVectorColumnIndex) ** 0.5 * &
                       (nEns * numRetainedEigen / (nEns - 1)) ** 0.5

  end subroutine getModulationFactor

end module enkf_mod

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
  use mpi_mod
  use utilities_mod
  use mathPhysConstants_mod
  use columnData_mod
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
  public :: enkf_rejectHighLatIR

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

contains

  !----------------------------------------------------------------------
  ! enkf_LETKFanalyses
  !----------------------------------------------------------------------
  subroutine enkf_LETKFanalyses(algorithm, numSubEns,  &
                                ensembleAnl, ensembleTrl, ensObs_mpiglobal,  &
                                stateVectorMeanAnl, &
                                wInterpInfo, maxNumLocalObs,  &
                                hLocalize, hLocalizePressure, vLocalize,  &
                                alphaRTPP, mpiDistribution)
    ! :Purpose: Local subroutine containing the code for computing
    !           the LETKF analyses for all ensemble members, ensemble
    !           mean.
    implicit none

    ! Arguments
    character(len=*)            :: algorithm
    integer                     :: numSubEns
    type(struct_ens), pointer   :: ensembleTrl
    type(struct_ens)            :: ensembleAnl
    type(struct_eob)            :: ensObs_mpiglobal
    type(struct_gsv)            :: stateVectorMeanAnl
    type(struct_enkfInterpInfo) :: wInterpInfo
    integer                     :: maxNumLocalObs
    real(8)                     :: hLocalize(:)
    real(8)                     :: hLocalizePressure(:)
    real(8)                     :: vLocalize
    real(8)                     :: alphaRTPP
    character(len=*)            :: mpiDistribution

    ! Locals
    integer :: nEns, nEnsPerSubEns, nEnsIndependentPerSubEns, nLev_M, ierr, matrixRank
    integer :: memberIndex, memberIndex1, memberIndex2
    integer :: memberIndexCV, memberIndexCV1, memberIndexCV2
    integer :: procIndex, procIndexSend, hLocIndex
    integer :: latIndex, lonIndex, stepIndex, varLevIndex, levIndex, levIndex2
    integer :: bodyIndex, localObsIndex, numLocalObs, numLocalObsFound
    integer :: countMaxExceeded, maxCountMaxExceeded, numGridPointWeights
    integer :: myNumLatLonRecv, myNumLatLonSend, numLatLonRecvMax
    integer :: numLatLonTotalUnique, latLonIndex, subEnsIndex, subEnsIndex2
    integer :: sendTag, recvTag, nsize, numRecv, numSend
    integer :: myLonBeg, myLonEnd, myLatBeg, myLatEnd, numVarLev
    integer :: myLonBegHalo, myLonEndHalo, myLatBegHalo, myLatEndHalo
    real(8) :: anlLat, anlLon, anlLogPres, distance, tolerance, localization

    integer, allocatable :: localBodyIndices(:)
    integer, allocatable :: myLatIndexesRecv(:), myLonIndexesRecv(:)
    integer, allocatable :: myLatIndexesSend(:), myLonIndexesSend(:)
    integer, allocatable :: myNumProcIndexesSend(:)
    integer, allocatable :: myProcIndexesRecv(:), myProcIndexesSend(:,:)
    integer, allocatable :: requestIdRecv(:), requestIdSend(:)
    integer, allocatable :: memberIndexSubEns(:,:), memberIndexSubEnsComp(:,:)

    real(8), allocatable :: distances(:)
    real(8), allocatable :: PaInv(:,:), PaSqrt(:,:), Pa(:,:), YbTinvR(:,:), YbTinvRYb(:,:)
    real(8), allocatable :: YbTinvRYb_CV(:,:)
    real(8), allocatable :: eigenValues(:), eigenVectors(:,:)
    real(8), allocatable :: eigenValues_CV(:), eigenVectors_CV(:,:)
    real(8), allocatable :: weightsTemp(:), weightsTemp2(:)
    real(8), allocatable :: weightsMembers(:,:,:,:), weightsMembersLatLon(:,:,:)
    real(8), allocatable :: weightsMean(:,:,:,:), weightsMeanLatLon(:,:,:)
    real(8), allocatable :: memberAnlPert(:)
    real(4), allocatable :: logPres_M_r4(:,:,:)

    real(4), pointer     :: meanTrl_ptr_r4(:,:,:,:), meanAnl_ptr_r4(:,:,:,:), meanInc_ptr_r4(:,:,:,:)
    real(4), pointer     :: memberTrl_ptr_r4(:,:,:,:), memberAnl_ptr_r4(:,:,:,:)

    character(len=4)     :: varLevel

    type(struct_hco), pointer :: hco_ens
    type(struct_vco), pointer :: vco_ens
    type(struct_gsv)          :: stateVectorMeanInc
    type(struct_gsv)          :: stateVectorMeanTrl

    logical :: firstTime = .true.

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

    nEns = ens_getNumMembers(ensembleAnl)
    nLev_M = ens_getNumLev(ensembleAnl, 'MM')
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
    allocate(YbTinvR(nEns,maxNumLocalObs))
    allocate(YbTinvRYb(nEns,nEns))
    allocate(eigenValues(nEns))
    allocate(eigenVectors(nEns,nEns))
    allocate(PaInv(nEns,nEns))
    allocate(PaSqrt(nEns,nEns))
    allocate(Pa(nEns,nEns))
    allocate(memberAnlPert(nEns))
    allocate(weightsTemp(nEns))
    allocate(weightsTemp2(nEns))
    weightsTemp(:) = 0.0d0
    weightsTemp2(:) = 0.0d0
    ! Weights for mean analysis
    allocate(weightsMean(nEns,1,myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
    weightsMean(:,:,:,:) = 0.0d0
    allocate(weightsMeanLatLon(nEns,1,myNumLatLonSend))
    weightsMeanLatLon(:,:,:) = 0.0d0
    ! Weights for member analyses
    allocate(weightsMembers(nEns,nEns,myLonBegHalo:myLonEndHalo,myLatBegHalo:myLatEndHalo))
    weightsMembers(:,:,:,:) = 0.0d0
    allocate(weightsMembersLatLon(nEns,nEns,myNumLatLonSend))
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

    ! Quantities needed for CVLETKF and CVLETKF-PERTOBS
    if (trim(algorithm) == 'CVLETKF' .or. trim(algorithm) == 'CVLETKF-PERTOBS') then
      nEnsPerSubEns = nEns / numSubEns
      if ( (nEnsPerSubEns * numSubEns) /= nEns ) then
        call utl_abort('enkf_LETKFanalyses: ensemble size not divisible by numSubEnsembles')
      end if
      if (numSubEns <= 1) then
        call utl_abort('enkf_LETKFanalyses: for CVLETKF(-PERTOBS) algorithm, numSubEns must be greater than 1')
      end if
      nEnsIndependentPerSubEns = nEns - nEnsPerSubEns
      allocate(YbTinvRYb_CV(nEnsIndependentPerSubEns,nEnsIndependentPerSubEns))
      allocate(eigenValues_CV(nEnsIndependentPerSubEns))
      allocate(eigenVectors_CV(nEnsIndependentPerSubEns,nEnsIndependentPerSubEns))
      allocate(memberIndexSubEns(nEnsPerSubEns,numSubEns))
      allocate(memberIndexSubEnsComp(nEnsIndependentPerSubEns,numSubEns))
      do subEnsIndex = 1, numSubEns
        do memberIndex = 1, nEnsPerSubEns
          memberIndexSubEns(memberIndex,subEnsIndex) =  &
               (subEnsIndex-1)*nEnsPerSubEns + memberIndex
        end do
      end do
      do subEnsIndex = 1, numSubEns
        memberIndex = 1
        do subEnsIndex2 = 1, numSubEns
          if (subEnsIndex2 == subEnsIndex) cycle

          memberIndexSubEnsComp(memberIndex:memberIndex+nEnsPerSubEns-1,subEnsIndex) =  &
            memberIndexSubEns(:,subEnsIndex2)
          memberIndex = memberIndex + nEnsPerSubEns
        end do
      end do

      write(*,*) 'nEns, numSubEns, nEnsPerSubEns, nEnsIndependentPerSubEns = ',  &
                 nEns, numSubEns, nEnsPerSubEns, nEnsIndependentPerSubEns
      do subEnsIndex = 1, numSubEns
        write(*,*) 'memberIndexSubEns = '
        write(*,*) memberIndexSubEns(:,subEnsIndex)
        write(*,*) 'memberIndexSubEnsComp = '
        write(*,*) memberIndexSubEnsComp(:,subEnsIndex)
      end do
    end if ! if CVLETKF(-PERTOBS) algorithm

    call lfn_Setup(LocFunctionWanted='FifthOrder')

    ! compute 3D field of log(pressure) needed for localization
    call enkf_computeLogPresM(logPres_M_r4,stateVectorMeanTrl)

    ! Compute the weights for ensemble mean and members
    countMaxExceeded = 0
    maxCountMaxExceeded = 0
    numGridPointWeights = 0
    LEV_LOOP: do levIndex = 1, nLev_M
      write(*,*) 'computing ensemble updates for vertical level = ', levIndex

      !
      ! First post all recv instructions for communication of weights
      !
      call tmg_start(103,'LETKF-commWeights')
      numSend = 0
      numRecv = 0
      do latLonIndex = 1, myNumLatLonRecv
        latIndex = myLatIndexesRecv(latLonIndex)
        lonIndex = myLonIndexesRecv(latLonIndex)
        procIndex = myProcIndexesRecv(latLonIndex)
        recvTag = (latIndex-1)*stateVectorMeanAnl%ni + lonIndex

        nsize = nEns
        numRecv = numRecv + 1
        call mpi_irecv( weightsMean(:,1,lonIndex,latIndex),  &
                        nsize, mpi_datyp_real8, procIndex-1, recvTag,  &
                        mpi_comm_grid, requestIdRecv(numRecv), ierr )
        nsize = nEns*nEns
        numRecv = numRecv + 1
        recvTag = recvTag + stateVectorMeanAnl%ni*stateVectorMeanAnl%nj
        call mpi_irecv( weightsMembers(:,:,lonIndex,latIndex),  &
                        nsize, mpi_datyp_real8, procIndex-1, recvTag,  &
                        mpi_comm_grid, requestIdRecv(numRecv), ierr )
      end do
      call tmg_stop(103)

      LATLON_LOOP: do latLonIndex = 1, myNumLatLonSend
        latIndex = myLatIndexesSend(latLonIndex)
        lonIndex = myLonIndexesSend(latLonIndex)

        numGridPointWeights = numGridPointWeights + 1

        ! lat-lon of the grid point for which we are doing the analysis
        anlLat = hco_ens%lat2d_4(lonIndex,latIndex)
        anlLon = hco_ens%lon2d_4(lonIndex,latIndex)
        anlLogPres = logPres_M_r4(lonIndex,latIndex,levIndex)

        ! Find which horizontal localization value to use for this analysis level
        hLocIndex = 1 + count(anlLogPres > hLocalizePressure(:))

        ! Get list of nearby observations and distances to gridpoint
        call tmg_start(9,'LETKF-getLocalBodyIndices')
        numLocalObs = eob_getLocalBodyIndices(ensObs_mpiglobal, localBodyIndices,     &
                                              distances, anlLat, anlLon, anlLogPres,  &
                                              hLocalize(hLocIndex), vLocalize, numLocalObsFound)
        if (numLocalObsFound > maxNumLocalObs) then
          countMaxExceeded = countMaxExceeded + 1
          maxCountMaxExceeded = max(maxCountMaxExceeded, numLocalObsFound)
        end if
        call tmg_stop(9)

        call tmg_start(91,'LETKF-calcWeights')

        ! Extract initial quantities YbTinvR and first term of PaInv (YbTinvR*Yb)
        do localObsIndex = 1, numLocalObs
          bodyIndex = localBodyIndices(localObsIndex)

          ! Compute value of localization function
          call tmg_start(18,'LETKF-locFunction')
          ! Horizontal
          localization = lfn_Response(distances(localObsIndex),hLocalize(hLocIndex))
          ! Vertical - use pressures at the grid point (not obs) location
          if (vLocalize > 0) then
            distance = abs( anlLogPres - ensObs_mpiglobal%logPres(bodyIndex) )
            localization = localization * lfn_Response(distance,vLocalize)
          end if
          call tmg_stop(18)

          do memberIndex = 1, nEns
            YbTinvR(memberIndex,localObsIndex) =  &
                 ensObs_mpiglobal%Yb_r4(memberIndex, bodyIndex) * &
                 localization * ensObs_mpiglobal%obsErrInv(bodyIndex)
          end do

          if (localObsIndex == 1) YbTinvRYb(:,:) = 0.0D0
          !$OMP PARALLEL DO PRIVATE (memberIndex1, memberIndex2)
          do memberIndex2 = 1, nEns
            do memberIndex1 = 1, nEns
              YbTinvRYb(memberIndex1,memberIndex2) =  &
                   YbTinvRYb(memberIndex1,memberIndex2) +  &
                   YbTinvR(memberIndex1,localObsIndex) * ensObs_mpiglobal%Yb_r4(memberIndex2, bodyIndex)
            end do
          end do
          !$OMP END PARALLEL DO

        end do ! localObsIndex

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
            call tmg_start(90,'LETKF-eigenDecomp')
            call utl_matInverse(Pa, nEns, inverseSqrt_opt=PaSqrt)
            call tmg_stop(90)

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

            ! Compute ensemble perturbation weights: (1-alphaRTPP)*[(Nens-1)^1/2*PaSqrt]+alphaRTPP*I
            weightsMembersLatLon(:,:,latLonIndex) = (1.0d0 - alphaRTPP) * sqrt(real(nEns - 1,8)) * PaSqrt(:,:)
            do memberIndex = 1, nEns
              weightsMembersLatLon(memberIndex,memberIndex,latLonIndex) = alphaRTPP +  &
                   weightsMembersLatLon(memberIndex,memberIndex,latLonIndex)
            end do

          else if (trim(algorithm) == 'CVLETKF') then
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
            ! Wa = (1-alphaRTPP) * 
            !      [ I - (Nens-1)^1/2 * E * 
            !        {(Nens-1)^-1/2*I - (Lambda + (Nens-1)*I)^-1/2} * Lambda^-1 *
            !        E^T * YbTinvRYb ]
            !      + alphaRTPP*I
            ! Loop over sub-ensembles
            do subEnsIndex = 1, numSubEns

              ! Use complement (independent) ens to get eigenValues/Vectors of Yb^T R^-1 Yb = E*Lambda*E^T
              call tmg_start(90,'LETKF-eigenDecomp')
              do memberIndexCV2 = 1, nEnsIndependentPerSubEns
                memberIndex2 = memberIndexSubEnsComp(memberIndexCV2, subEnsIndex)
                do memberIndexCV1 = 1, nEnsIndependentPerSubEns
                  memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
                  YbTinvRYb_CV(memberIndexCV1,memberIndexCV2) = YbTinvRYb(memberIndex1,memberIndex2)
                end do
              end do
              tolerance = 1.0D-50
              call utl_eigenDecomp(YbTinvRYb_CV, eigenValues_CV, eigenVectors_CV, tolerance, matrixRank)
              call tmg_stop(90)

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

            ! Apply RTPP
            weightsMembersLatLon(:,:,latLonIndex) =  &
                 (1.0d0 - alphaRTPP) * weightsMembersLatLon(:,:,latLonIndex)
            do memberIndex = 1, nEns
              weightsMembersLatLon(memberIndex,memberIndex,latLonIndex) = alphaRTPP +  &
                   weightsMembersLatLon(memberIndex,memberIndex,latLonIndex)
            end do

          else if (trim(algorithm) == 'CVLETKF-PERTOBS') then
            !
            ! Weight calculation for perturbed-obs cross-validation LETKF algorithm
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
            ! Wa   = (1-alphaRTPP) * Wa
            !
            ! Loop over sub-ensembles
            do subEnsIndex = 1, numSubEns

              ! Use complement (independent) ens to get eigenValues/Vectors of Yb^T R^-1 Yb = E*Lambda*E^T
              call tmg_start(90,'LETKF-eigenDecomp')
              do memberIndexCV2 = 1, nEnsIndependentPerSubEns
                memberIndex2 = memberIndexSubEnsComp(memberIndexCV2, subEnsIndex)
                do memberIndexCV1 = 1, nEnsIndependentPerSubEns
                  memberIndex1 = memberIndexSubEnsComp(memberIndexCV1, subEnsIndex)
                  YbTinvRYb_CV(memberIndexCV1,memberIndexCV2) = YbTinvRYb(memberIndex1,memberIndex2)
                end do
              end do
              tolerance = 1.0D-50
              call utl_eigenDecomp(YbTinvRYb_CV, eigenValues_CV, eigenVectors_CV, tolerance, matrixRank)
              call tmg_stop(90)

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

            ! Apply RTPP
            weightsMembersLatLon(:,:,latLonIndex) =  &
                 (1.0d0 - alphaRTPP) * weightsMembersLatLon(:,:,latLonIndex)
            do memberIndex = 1, nEns
              weightsMembersLatLon(memberIndex,memberIndex,latLonIndex) = alphaRTPP +  &
                   weightsMembersLatLon(memberIndex,memberIndex,latLonIndex)
            end do

          else

            call utl_abort('UNKNOWN LETKF ALGORITHM')

          end if

        else
          ! no observations near this grid point, set weights to zero
          weightsMeanLatLon(:,1,latLonIndex) = 0.0d0
          weightsMembersLatLon(:,:,latLonIndex) = 0.0d0
        end if ! numLocalObs > 0

        call tmg_stop(91)

        !
        ! Now post all send instructions (each lat-lon may be sent to multiple tasks)
        !
        call tmg_start(103,'LETKF-commWeights')
        latIndex = myLatIndexesSend(latLonIndex)
        lonIndex = myLonIndexesSend(latLonIndex)
        do procIndex = 1, myNumProcIndexesSend(latLonIndex)
          sendTag = (latIndex-1)*stateVectorMeanAnl%ni + lonIndex
          procIndexSend = myProcIndexesSend(latLonIndex, procIndex)

          nsize = nEns
          numSend = numSend + 1
          call mpi_isend( weightsMeanLatLon(:,1,latLonIndex),  &
                          nsize, mpi_datyp_real8, procIndexSend-1, sendTag,  &
                          mpi_comm_grid, requestIdSend(numSend), ierr )
          nsize = nEns*nEns
          numSend = numSend + 1
          sendTag = sendTag + stateVectorMeanAnl%ni*stateVectorMeanAnl%nj
          call mpi_isend( weightsMembersLatLon(:,:,latLonIndex),  &
                          nsize, mpi_datyp_real8, procIndexSend-1, sendTag,  &
                          mpi_comm_grid, requestIdSend(numSend), ierr )
        end do
        call tmg_stop(103)

      end do LATLON_LOOP

      !
      ! Wait for communiations to finish before continuing
      !
      call tmg_start(103,'LETKF-commWeights')
      if (firstTime) write(*,*) 'numSend/Recv = ', numSend, numRecv
      firstTime = .false.

      if ( numRecv > 0 ) then
        call mpi_waitAll(numRecv, requestIdRecv(1:numRecv), MPI_STATUSES_IGNORE, ierr)
      end if

      if ( numSend > 0 ) then
        call mpi_waitAll(numSend, requestIdSend(1:numSend), MPI_STATUSES_IGNORE, ierr)
      end if

      call tmg_stop(103)

      !
      ! Interpolate weights from coarse to full resolution
      !
      call tmg_start(92,'LETKF-interpolateWeights')
      if (wInterpInfo%latLonStep > 1) then
        call enkf_interpWeights(wInterpInfo, weightsMean)
        call enkf_interpWeights(wInterpInfo, weightsMembers)
      end if
      call tmg_stop(92)

      call tmg_start(100,'LETKF-applyWeights')

      !
      ! Apply the weights to compute the ensemble mean and members
      !
      meanInc_ptr_r4 => gsv_getField_r4(stateVectorMeanInc)
      meanTrl_ptr_r4 => gsv_getField_r4(stateVectorMeanTrl)
      meanAnl_ptr_r4 => gsv_getField_r4(stateVectorMeanAnl)
      do latIndex = myLatBeg, myLatEnd
        LON_LOOP5: do lonIndex = myLonBeg, myLonEnd

          ! skip this grid point if all weights zero (no nearby obs)
          if (all(weightsMean(:,1,lonIndex,latIndex) == 0.0d0)) cycle LON_LOOP5

          ! Compute the ensemble mean increment and analysis
          do varLevIndex = 1, numVarLev
            ! Only treat varLevIndex values that correspond with current levIndex
            varLevel = vnl_varLevelFromVarname(gsv_getVarNameFromK(stateVectorMeanInc,varLevIndex))
            if (varLevel == 'SF' .or. varLevel == 'SFMM' .or. varLevel == 'SFTH') then
              levIndex2 = nLev_M
            else
              levIndex2 = gsv_getLevFromK(stateVectorMeanInc,varLevIndex)
            end if
            if (levIndex2 /= levIndex) cycle
            memberTrl_ptr_r4 => ens_getOneLev_r4(ensembleTrl,varLevIndex)
            do stepIndex = 1, tim_nstepobsinc
              ! mean increment
              do memberIndex = 1, nEns
                meanInc_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) =  &
                     meanInc_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) +  &
                     weightsMean(memberIndex,1,lonIndex,latIndex) *  &
                     (memberTrl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) -  &
                      meanTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex))
              end do ! memberIndex
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
            if (varLevel == 'SF' .or. varLevel == 'SFMM' .or. varLevel == 'SFTH') then
              levIndex2 = nLev_M
            else
              levIndex2 = gsv_getLevFromK(stateVectorMeanInc,varLevIndex)
            end if
            if (levIndex2 /= levIndex) cycle
            memberTrl_ptr_r4 => ens_getOneLev_r4(ensembleTrl,varLevIndex)
            memberAnl_ptr_r4 => ens_getOneLev_r4(ensembleAnl,varLevIndex)
            do stepIndex = 1, tim_nstepobsinc
              ! Compute analysis member perturbation
              memberAnlPert(:) = 0.0d0
              do memberIndex2 = 1, nEns
                do memberIndex1 = 1, nEns
                  memberAnlPert(memberIndex2) = memberAnlPert(memberIndex2) + &
                       weightsMembers(memberIndex1,memberIndex2,lonIndex,latIndex) *  &
                       (memberTrl_ptr_r4(memberIndex1,stepIndex,lonIndex,latIndex) -  &
                       meanTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex))
                end do ! memberIndex1
              end do ! memberIndex2
              ! Add analysis member perturbation to mean analysis
              memberAnl_ptr_r4(:,stepIndex,lonIndex,latIndex) =  &
                   meanAnl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) + memberAnlPert(:)
            end do ! stepIndex
          end do ! varLevIndex

        end do LON_LOOP5
      end do

      call tmg_stop(100)

    end do LEV_LOOP

    if (countMaxExceeded > 0) then
      write(*,*) 'enkf_LETKFanalyses: WARNING: Found more local obs than specified max number at ', &
                 real(100*countMaxExceeded)/real(numGridPointWeights), '% of grid points.'
      write(*,*) '                      Maximum number found was ', maxCountMaxExceeded,  &
                 ' which is greater than specified number ', maxNumLocalObs
      write(*,*) '                      Therefore will keep closest obs only.'
    end if

    call tmg_start(19,'LETKF-barr')
    call rpn_comm_barrier('GRID',ierr)
    call tmg_stop(19)

    call gsv_deallocate(stateVectorMeanInc)
    call gsv_deallocate(stateVectorMeanTrl)

    write(*,*) 'enkf_LETKFanalyses: done'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine enkf_LETKFanalyses

  !----------------------------------------------------------------------
  ! enkf_computeLogPresM (private subroutine)
  !----------------------------------------------------------------------
  subroutine enkf_computeLogPresM(logPres_M_r4,stateVectorMeanTrl)
    !
    !:Purpose:  Compute extract global 3D log pressure field from supplied
    !           stateVector.
    !
    implicit none

    ! Arguments
    real(4), allocatable :: logPres_M_r4(:,:,:)
    type(struct_gsv) :: stateVectorMeanTrl

    ! Locals
    integer          :: nLev_M, nsize, ierr
    real(4), pointer :: logPres_M_ptr_r4(:,:,:)
    type(struct_gsv) :: stateVectorMeanTrlPressure
    type(struct_gsv) :: stateVectorMeanTrlPressure_1step

    nLev_M = gsv_getNumLev(stateVectorMeanTrl, 'MM')

    ! Compute background ens mean 3D log pressure and make mpiglobal for vertical localization
    call gsv_allocate( stateVectorMeanTrlPressure, tim_nstepobsinc,  &
                       stateVectorMeanTrl%hco, stateVectorMeanTrl%vco, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0','P_M','P_T'/) )
    call gsv_zero(stateVectorMeanTrlPressure)
    call gsv_copy(stateVectorMeanTrl, stateVectorMeanTrlPressure, allowMismatch_opt=.true.)
    call gvt_transform(stateVectorMeanTrlPressure,'PsfcToP_nl')
    if (mpi_myid == 0) then
      call gsv_allocate( stateVectorMeanTrlPressure_1step, 1,  &
                         stateVectorMeanTrl%hco, stateVectorMeanTrl%vco, dateStamp_opt=tim_getDateStamp(),  &
                         mpi_local_opt=.false., &
                         dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0','P_M','P_T'/) )
    end if
    call gsv_transposeTilesToStep(stateVectorMeanTrlPressure_1step, stateVectorMeanTrlPressure, (tim_nstepobsinc+1)/2)
    call gsv_deallocate(stateVectorMeanTrlPressure)
    allocate(logPres_M_r4(stateVectorMeanTrlPressure%ni, stateVectorMeanTrlPressure%nj, nLev_M))
    if (mpi_myid == 0) then
      logPres_M_ptr_r4 => gsv_getField3d_r4(stateVectorMeanTrlPressure_1step,'P_M')
      logPres_M_r4(:,:,:) = log(logPres_M_ptr_r4(:,:,:))
    end if
    nsize = stateVectorMeanTrlPressure%ni * stateVectorMeanTrlPressure%nj * nLev_M
    call rpn_comm_bcast(logPres_M_r4, nsize, 'mpi_real4', 0, 'GRID', ierr)

  end subroutine enkf_computeLogPresM

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
          myProcIndexesRecv(myNumLatLonRecv) = mpi_myid+1
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
      allocate(allNumLatLonRecv(mpi_nprocs))
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
      allocate(allLatIndexesRecv(numLatLonRecvMax, mpi_nprocs))
      allocate(allLonIndexesRecv(numLatLonRecvMax, mpi_nprocs))
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
      do procIndex = 1, mpi_nprocs
        WEIGHTS1LEV_LOOP: do latLonIndex = 1, allNumLatLonRecv(procIndex)
          if (enkf_latLonAlreadyFound(allLatIndexesRecv, allLonIndexesRecv, latLonIndex, procIndex)) &
               cycle WEIGHTS1LEV_LOOP
          ! Count the total number of weights
          numLatLonTotalUnique = numLatLonTotalUnique + 1

          ! Round-robin distribution of master list across mpi tasks
          procIndexSend = 1 + mod(numLatLonTotalUnique-1, mpi_nprocs)

          ! Store the lat-lon indexes of the weights I am responsible for
          if (procIndexSend == (mpi_myid+1)) then
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
      allocate(allNumLatLonSend(mpi_nprocs))
      call rpn_comm_allgather(myNumLatLonSend, 1, "mpi_integer",  &
                              allNumLatLonSend, 1,"mpi_integer", "GRID", ierr)
      allocate(allLatIndexesSend(numLatLonRecvMax, mpi_nprocs))
      allocate(allLonIndexesSend(numLatLonRecvMax, mpi_nprocs))
      call rpn_comm_allgather(myLatIndexesSend, numLatLonRecvMax, "mpi_integer",  &
                              allLatIndexesSend, numLatLonRecvMax, "mpi_integer",  &
                              "GRID", ierr)
      call rpn_comm_allgather(myLonIndexesSend, numLatLonRecvMax, "mpi_integer",  &
                              allLonIndexesSend, numLatLonRecvMax, "mpi_integer",  &
                              "GRID", ierr)

      ! Figure out which mpi tasks I will need to send my results to
      allocate(myProcIndexesSend(myNumLatLonSend,mpi_nprocs))
      allocate(myNumProcIndexesSend(myNumLatLonSend))
      myProcIndexesSend(:,:) = -1
      myNumProcIndexesSend(:) = 0
      do latLonIndex = 1, myNumLatLonSend
        do procIndex = 1, mpi_nprocs
          if ( any( (myLatIndexesSend(latLonIndex) == allLatIndexesRecv(1:allNumLatLonRecv(procIndex), procIndex)) .and.  &
                    (myLonIndexesSend(latLonIndex) == allLonIndexesRecv(1:allNumLatLonRecv(procIndex), procIndex)) ) ) then
            myNumProcIndexesSend(latLonIndex) = myNumProcIndexesSend(latLonIndex) + 1
            myProcIndexesSend(latLonIndex,myNumProcIndexesSend(latLonIndex)) = procIndex
          end if
        end do
      end do

      ! Figure out which mpi tasks I will receive the results from
      do latLonIndex = 1, myNumLatLonRecv
        do procIndex = 1, mpi_nprocs
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

    myLonBegHalo = wInterpInfo%myLonBegHalo
    myLonEndHalo = wInterpInfo%myLonEndHalo
    myLatBegHalo = wInterpInfo%myLatBegHalo
    myLatEndHalo = wInterpInfo%myLatEndHalo
    numMembers1 = size(weights,1)
    numMembers2 = size(weights,2)

    do latIndex = myLatBegHalo, myLatEndHalo
      do lonIndex = myLonBegHalo, myLonEndHalo
        if (wInterpInfo%numIndexes(lonIndex,latIndex) > 0) then

          ! Interpolation for ensemble member perturbation weight fields
          weights(:,:,lonIndex,latIndex) = 0.0D0
          if (wInterpInfo%lonIndexes(1,lonIndex,latIndex) == 0) cycle ! temporary until all interpolation setup completed
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
        end if ! numIndexes > 0
      end do ! lonIndex
    end do ! latIndex

  end subroutine enkf_interpWeights

  !--------------------------------------------------------------------------
  ! enkf_modifyAMSUBobsError
  !--------------------------------------------------------------------------
  subroutine enkf_modifyAMSUBobsError(obsSpaceData)
    implicit none

    ! arguments:
    type(struct_obs), target  :: obsSpaceData

    ! locals:
    real(obs_real), parameter :: AMSUB_trop_oer = 1.0 ! assumed value for AMSU-B obs error in tropics
    integer            :: headerIndex, bodyIndex, bodyIndexBeg, bodyIndexEnd, codeType
    real(obs_real)     :: lat_obs

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
    integer        :: headerIndex, bodyIndex, bodyIndexBeg, bodyIndexEnd, codeType
    real(obs_real) :: lat_obs

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

end module enkf_mod
